!
! S3COM
! Copyright (c) 2022, University of Lille
!
! Redistribution and use in source and binary forms, with or without modification,
! are permitted provided that the following conditions are met:
!  1. Redistributions of source code must retain the above copyright notice,
!     this list of conditions and the following disclaimer.
!  2. Redistributions in binary form must reproduce the above copyright notice,
!     this list of conditions and the following disclaimer in the documentation
!     and/or other materials provided with the distribution.
!  3. Neither the name of the copyright holder nor the names of its contributors
!     may be used to endorse or promote products derived from this software without
!     specific prior written permission.

! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
! EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
! OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
! SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! History
! Jan 2022 - O. Sourdeval - Original version
!

module mod_ret_run
   
   use s3com_types,     only: wp, wpi, type_nml, type_model, type_s3com, type_rttov_opt, type_cld
   use s3com_config,    only: nstates
   use mod_ret_utils,   only: ret_update, ret_update_mat, ret_cov_mat, ret_J, J_meas, J_state, ret_K, ret_Sx, ret_Xip1, idx_liq
   use mod_rttov_utils, only: find_ret_idx_rttov
   use mod_utils_math,  only: inverse
   
   implicit none
   
   private
   public :: ret_run
   
contains
   
   !> @brief Run the retrieval algorithm using the optimal estimation (OE) method
   !! @details This subroutine runs the OE to retrieve the liquid water path and the cloud top droplet number concentration of 
   !! single-layer homogeneous liquid clouds.
   !! @param[in] rttov_atm_ret   rttov atmospheric structure
   !! @param[in] rttov_opt       rttov options structure
   !! @param[inout] s3com        s3com retrieval structure
   !! @param[in] cld             cloud structure
   subroutine ret_run(rttov_atm_ret, rttov_opt, s3com, cld)
      
      ! Input
      type(type_model), intent(in)  :: rttov_atm_ret
      type(type_rttov_opt), intent(in) :: rttov_opt
      type(type_cld), intent(in) :: cld
      
      ! Input/output
      type(type_s3com), intent(inout)  :: s3com
      
      ! Internal
      integer(wpi) :: iter, ipoint, imeas
      integer(wpi), dimension(:), allocatable :: idx_ret
      
      idx_ret = find_ret_idx_rttov(s3com)
      
      s3com%ret%flag_conv_test(idx_ret) = .false.
      
      iter = 1
      
      s3com%ret%gamma(idx_ret) = 10.0_wp
      
      ! Run the OE method using the Levenberg-Marquardt method
      do while (any(.not. s3com%ret%flag_conv_test(idx_ret)))
         
         write(6,*) "Iteration:", iter
         
         call ret_update_mat(rttov_atm_ret, rttov_opt, s3com, cld)
         
         ! Xi
         do ipoint = 1, size(idx_ret)
            
            call inverse(s3com%ret%Se(idx_ret(ipoint),:,:),s3com%ret%Se_inv(idx_ret(ipoint),:,:),s3com%ret%nmeas)
            call inverse(s3com%ret%Sa(idx_ret(ipoint),:,:),s3com%ret%Sa_inv(idx_ret(ipoint),:,:),s3com%ret%nstates)
            
            s3com%ret%Kt(idx_ret(ipoint),:,:) = transpose(s3com%ret%K(idx_ret(ipoint),:,:))
            
            s3com%ret%Ji(idx_ret(ipoint)) = ret_J(s3com,idx_ret(ipoint))
            
            write(6,*) "Before update:"
            write(6,*) " - Cost function =", s3com%ret%Ji(idx_ret(ipoint))
            if (iter .eq. 1) then
               write(6,*) " - LWP =", s3com%ret%Xi(idx_ret(ipoint),1), "+/-", sqrt(s3com%ret%Sa(idx_ret(ipoint),1,1))
               write(6,*) " - CDNC_top =", s3com%ret%Xi(idx_ret(ipoint),2)*1.0E-06_wp, "+/-", &
               sqrt(s3com%ret%Sa(idx_ret(ipoint),2,2))*1.0E-06_wp
            else
               write(6,*) " - LWP =", s3com%ret%Xi(idx_ret(ipoint),1), "+/-", sqrt(s3com%ret%Sx(idx_ret(ipoint),1,1))
               write(6,*) " - CDNC_top =", s3com%ret%Xi(idx_ret(ipoint),2)*1.0E-06_wp, "+/-", &
               sqrt(s3com%ret%Sx(idx_ret(ipoint),2,2))*1.0E-06_wp
            endif
            write(6,*)
            
         enddo
         
         call ret_update(rttov_atm_ret, rttov_opt, size(idx_ret), idx_ret, s3com, cld)
         
         ! Xi+1
         do ipoint = 1, size(idx_ret)
            write(6,*) "After update:"
            write(6,*) " - Cost function =", s3com%ret%Jip1(idx_ret(ipoint))
            write(6,*) " - LWP =", s3com%ret%Xip1(idx_ret(ipoint),1), "+/-", sqrt(s3com%ret%Sx(idx_ret(ipoint),1,1))
            write(6,*) " - CDNC =", s3com%ret%Xip1(idx_ret(ipoint),2)*1.0E-06_wp, "+/-", &
            sqrt(s3com%ret%Sx(idx_ret(ipoint),2,2))*1.0E-06_wp
            write(6,*)
         enddo
         
         ! Ajust the inverse of step size for the next iteration
         do ipoint = 1, size(idx_ret)
            if (s3com%ret%Jip1(idx_ret(ipoint)) / s3com%ret%Ji(idx_ret(ipoint)) .gt. 1.0_wp) then
               write(6,*) "Cost function has increased!"
               s3com%ret%gamma(idx_ret(ipoint)) = s3com%ret%gamma(idx_ret(ipoint)) * 10.0_wp
            else
               s3com%ret%gamma(idx_ret(ipoint)) = s3com%ret%gamma(idx_ret(ipoint)) / 2.0_wp
            endif
         enddo
         
         !write(*,*) "gamma =", s3com%ret%gamma
         !write(*,*) "s3com%ret%Ji_meas(idx_ret)", s3com%ret%Ji_meas(idx_ret)
         !write(*,*) "s3com%ret%Ji_state(idx_ret)", s3com%ret%Ji_state(idx_ret)
         
         ! Convergence criteria
         ! J_meas < 0.6
         s3com%ret%flag_conv_test(idx_ret) = (s3com%ret%Ji_meas(idx_ret) .lt. (s3com%ret%nmeas/10.0_wp) &
         ! J_state < 0.6
         .and. s3com%ret%Ji_state(idx_ret) .lt. (s3com%ret%nstates/10.0_wp) & 
         ! Additional criterion to stop the OE when the cost function stagnates 
         .or. abs(s3com%ret%Jip1(idx_ret) - s3com%ret%Ji(idx_ret)) .lt. 1.0E-04_wp)
         
         write(6,*) "s3com%ret%flag_conv_test(idx_ret) =", s3com%ret%flag_conv_test(idx_ret)
         
         do ipoint = 1, size(idx_ret)
            write(6,*) "Y =", (s3com%ret%Y(ipoint,imeas), imeas=1,s3com%ret%nmeas)
            write(6,*) "F =", (s3com%ret%F(ipoint,imeas), imeas=1,s3com%ret%nmeas)
         enddo
         
         write(6,*)
         
         s3com%ret%n_iter(idx_ret) = iter
         
         if (iter .gt. 50) exit
         
         iter = iter + 1
         
      enddo
      
      write(*,*) "Number of iterations =", s3com%ret%n_iter(idx_ret)
      write(*,*) "Ended!"
      
   end subroutine ret_run
   
end module mod_ret_run

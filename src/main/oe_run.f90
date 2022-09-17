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

MODULE MOD_OE_RUN

  use s3com_types, only: wp, type_rttov_atm, type_rttov_opt, type_s3com
  USE mod_oe_utils, only:get_jacobian, get_covmat, update_oemat, oe_cost, oe_Sx, oe_Xip1, oe_cost_meas, idx_ice, update_oe
  USE mod_rttov_utils, only: idx_rttov 
  use mod_utils_math, only: inverse
  
  IMPLICIT NONE

CONTAINS

  SUBROUTINE oe_run(oe,rttov_atm_oe,rttov_opt)

    ! Input variables
    TYPE(type_rttov_opt), INTENT(IN) :: rttov_opt

    ! Input/Output variables
    TYPE(type_rttov_atm), INTENT(INOUT) :: rttov_atm_oe    
    TYPE(type_s3com), INTENT(INOUT) :: oe

    ! Internal variables
    TYPE(type_s3com) :: oe_ip1
    
    INTEGER(KIND=4) :: ipoint, iloop, iloop1, iloop_g
    INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: idx_iwp, idx_oe
    
    LOGICAL, DIMENSION(:), ALLOCATABLE :: flag_rttov_save
    
    
    
    idx_oe = idx_rttov(oe)
    allocate(flag_rttov_save(size(oe%flag_rttov)))

    
    oe%g_meas = 1.0D5
    iloop=1; ipoint=1

    oe%flag_testconv(:) = .TRUE.     

    write(*,*) size(idx_oe)

    oe%flag_testconv(idx_oe) = .FALSE.

    oe%iwp(idx_oe) = oe%Xa(idx_oe,1)

    DO WHILE(ANY(.NOT.oe%flag_testconv))
       write(*,*) "iloop:", iloop        

       flag_rttov_save = oe%flag_rttov

       call update_oemat(rttov_atm_oe,rttov_opt,oe)

       oe%stepsize(idx_oe) = 10.0D0

       DO ipoint=1,size(idx_oe)
          call inverse(oe%Se(idx_oe(ipoint),:,:),oe%Se_i(idx_oe(ipoint),:,:),oe%nmeas)
          call inverse(oe%Sa(idx_oe(ipoint),:,:),oe%Sa_i(idx_oe(ipoint),:,:),oe%nstates)

          oe%Xi(idx_oe(ipoint),:) = (/oe%iwp(idx_oe(ipoint))/)
          oe%Kt(idx_oe(ipoint),:,:) = TRANSPOSE(oe%K(idx_oe(ipoint),:,:))
          oe%Gi(idx_oe(ipoint)) = oe_cost(oe,idx_oe(ipoint))

          ! write(*,*) "Before update: point=", idx_oe(ipoint)," cost=", &
          !      oe%Gi(idx_oe(ipoint)), "iwp=", oe%Xi(idx_oe(ipoint),1)
       END DO

       oe%i_stepsize(idx_oe(:)) = 0

       CALL update_oe(oe,oe_ip1,size(idx_oe),idx_oe,rttov_atm_oe,rttov_opt)

       ! DO ipoint=1,size(idx_oe)
       !    write(*,*) "After update: point=", idx_oe(ipoint)," cost=", &
       !         oe_ip1%Gip1(idx_oe(ipoint)), "iwp=", oe_ip1%Xip1(idx_oe(ipoint),1)
       ! END DO

       ! Only run an update over pixels that won't work, reset the original flag after that
       ! oe%flag_rttov(idx_oe) = oe_ip1%Gip1(idx_oe)/oe_ip1%Gi(idx_oe).GT.1.0D0; idx_oe = idx_rttov(oe)

       ! iloop_g = 0
       ! DO WHILE(ANY(oe%flag_rttov).AND.(iloop_g.LT.10))
       !    write(*,*) "cost has increased!"
       !    oe%stepsize(idx_oe) = oe%stepsize(idx_oe) * 5.0D0           
       !    CALL update_oe(oe,oe_ip1,size(idx_oe),idx_oe,rttov_atm_oe,rttov_opt)

       !    oe%flag_rttov = oe_ip1%Gip1/oe_ip1%Gi.GT.1.0D0; idx_oe = idx_rttov(oe)
       !    iloop_g = iloop_g + 1
       ! END DO
       ! oe%flag_rttov = flag_rttov_save; idx_oe = idx_rttov(oe)

       DO ipoint=1,size(idx_oe)
          IF(oe_ip1%i_stepsize(idx_oe(ipoint)).GT.0) &
               oe_ip1%stepsize(idx_oe(ipoint)) = oe_ip1%stepsize(idx_oe(ipoint)) / (2.0D0*oe_ip1%i_stepsize(idx_oe(ipoint)))
          oe_ip1%stepsize(idx_oe(ipoint)) = oe_ip1%stepsize(idx_oe(ipoint)) / 5.0D0
       END DO

       oe = oe_ip1
       DO ipoint=1,size(idx_oe)        
          oe%Xi(idx_oe(ipoint),:) = oe%Xip1(idx_oe(ipoint),:)
          oe%g_meas(idx_oe(ipoint)) = oe_cost_meas(oe,idx_oe(ipoint))
       END DO

       oe%flag_testconv(idx_oe) = oe%g_meas(idx_oe).LT.(oe%nmeas/1.0)
       oe%flag_rttov(idx_oe) = .NOT.oe%flag_testconv(idx_oe); idx_oe = idx_rttov(oe)
       oe%n_iter(idx_oe) = iloop

       if(iloop.GT.15) exit

       iloop = iloop + 1

    end do

    deallocate(flag_rttov_save)
    write(*,*) 'finished'


  END SUBROUTINE oe_run
    
END MODULE MOD_OE_RUN

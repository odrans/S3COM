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

!> Useful tools used for and by RTTOV
module mod_rttov_utils
   
   use s3com_types,    only: wp, wpi, type_s3com
   use rttov_const,    only: errorstatus_success
   use rttov_unix_env, only: rttov_exit
   
   implicit none
   
   private
   public :: find_idx_rttov, find_ret_idx_rttov, check_rttov_status, get_rttov_model
   
contains
   
   ! ============================================================================================================================
   !> @brief Return the number of profiles to be processed in RTTOV
   !! @param[in] s3com        s3com data structure
   !! @param[out] idx_rttov   number of profiles to be processed
   function find_idx_rttov(s3com) result(idx_rttov)
      
      ! Input
      type(type_s3com), intent(in) :: s3com
      
      ! Output
      integer(wpi), dimension(:), allocatable :: idx_rttov
      
      ! Internal
      integer(wpi), dimension(s3com%npoints) :: idx_all
      integer(wpi) :: idx, i
      
      idx_all = 0
      idx = 1
      
      do i = 1, s3com%npoints
         if (s3com%flag_rttov(i)) then
            idx_all(idx) = i
            idx = idx + 1
         endif
      enddo
      
      idx = idx - 1
      
      allocate(idx_rttov(idx)); idx_rttov(1:idx) = idx_all(1:idx)
      
      return
      
   end function find_idx_rttov
   ! ============================================================================================================================
   
   ! ============================================================================================================================
   !> @brief Return the number of profiles to be processed in RTTOV for the inversion of liquid cloud properties
   !! @param[in] s3com            s3com data structure
   !! @param[out] ret_idx_rttov   number of profiles to be processed
   function find_ret_idx_rttov(s3com) result(ret_idx_rttov)
      
      ! Input
      type(type_s3com), intent(in) :: s3com
      
      ! Output
      integer(wpi), dimension(:), allocatable :: ret_idx_rttov
      
      ! Internal
      integer(wpi), dimension(s3com%ret%npoints) :: idx_all
      integer(wpi) :: idx, i
      
      idx_all = 0
      idx = 1
      
      do i = 1, s3com%ret%npoints
         if (s3com%ret%flag_rttov(i)) then
            idx_all(idx) = i
            idx = idx + 1
         endif
      enddo
      
      idx = idx - 1
      
      allocate(ret_idx_rttov(idx)); ret_idx_rttov(1:idx) = idx_all(1:idx)
      
      return
      
   end function find_ret_idx_rttov
   ! ============================================================================================================================
   
   ! ============================================================================================================================
   !> @brief Return error status and location
   !! @param[in] status     error status
   !! @param[in] location   location of the error
   subroutine check_rttov_status(status, location)
      
      ! Input
      integer, intent(in) :: status
      character(len=*), intent(in) :: location
      
      if (status /= errorstatus_success) then
         write(6,*) location
         call rttov_exit(status)
      endif
      
      return
      
   end subroutine check_rttov_status
   ! ============================================================================================================================
   
   ! ============================================================================================================================
   !> @brief Return the name of the RTTOV model to use (direct, K or TL)
   !! @param[in] s3com    s3com data structure
   !! @param[out] model   name of the RTTOV model
   function get_rttov_model(s3com) result(model)
      
      ! Input
      type(type_s3com), intent(in) :: s3com
      
      ! Output
      character(len=8) :: model
      
      if(s3com%jac%do_jacobian_calc) then
         model = "jacobian"
      else if(s3com%k_tl%do_k_tl_calc) then
         model = "TL"
      else
         model = "direct"
      endif
      
      return
      
   end function get_rttov_model
   ! ============================================================================================================================
   
end module mod_rttov_utils

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

module mod_ret_cloud_model
   
   use s3com_types,  only: type_model, type_s3com
   use s3com_config, only: wp, lwp_lay_threshold
   
   implicit none
   
   private
   public :: init_cloud_z, init_cloud_prof
   
contains
   
   ! ----------------------------------------------------------------------------------------------------------------------------
   !> @brief
   !> @detail This subroutine extracts the cloud base and top altitude for a liquid cloud, and calculate the liquid water path of
   ! the cloud
   ! @param[in] rttov_atm data structure
   ! @param[inout] s3com RET retrieval data structure
   subroutine init_cloud_z(rttov_atm, s3com)
      
      ! Inputs/outputs
      type(type_model), intent(in) :: rttov_atm
      type(type_s3com), intent(inout) :: s3com
      
      ! Internal
      integer(kind=4) :: npoints, nlayers, ipoint, ilayer
      
      ! Get the number of points and layers from rttov_atm
      npoints = rttov_atm%npoints
      nlayers = rttov_atm%nlayers
      
      ! Loop over each points
      do ipoint = 1, npoints
         
         ! Loop over each layer
         do ilayer = 1, nlayers
            
            ! If the liquid water path of the current layer is above the threshold value,
            ! set the index of this layer as the top of the liquid cloud layer and exit
            if (rttov_atm%lwc(ipoint,ilayer)*rttov_atm%dz(ipoint,ilayer) .gt. lwp_lay_threshold) then
               s3com%ret%ztop_liq_idx(ipoint) = ilayer
               exit
            endif
            
         enddo
         
         ! Loop over each layer in reverse order
         do ilayer = nlayers, 1, -1
            
            ! If the liquid water path of the current layer is above the threshold value,
            ! set the index of this layer as the base of the liquid cloud layer and exit
            if (rttov_atm%lwc(ipoint,ilayer)*rttov_atm%dz(ipoint,ilayer) .gt. lwp_lay_threshold) then
               s3com%ret%zbase_liq_idx(ipoint) = ilayer
               exit
            endif
            
         enddo
         
         ! If both the top and base of the liquid cloud layer are found,
         ! calculate the corresponding heights and liquid water path of the cloud
         if (s3com%ret%ztop_liq_idx(ipoint) .ne. 0 .and. s3com%ret%zbase_liq_idx(ipoint) .ne. 0) then
            s3com%ret%ztop_liq(ipoint)  = rttov_atm%z(ipoint,s3com%ret%ztop_liq_idx(ipoint))
            s3com%ret%zbase_liq(ipoint) = rttov_atm%z(ipoint,s3com%ret%zbase_liq_idx(ipoint))
            s3com%ret%lwp(ipoint)       = sum(rttov_atm%lwc(ipoint,1:nlayers)*rttov_atm%dz(ipoint,1:nlayers))
         endif
         
         ! If the liquid water path of the cloud is under the threshold value, this is not a liquid cloud
         if (s3com%ret%lwp(ipoint) .lt. 1e-4) s3com%ret%lwp(ipoint) = 0._wp
         
      enddo
      
   end subroutine init_cloud_z
   ! ----------------------------------------------------------------------------------------------------------------------------
   
   ! ----------------------------------------------------------------------------------------------------------------------------
   subroutine init_cloud_prof(rttov_atm, s3com)
      
      use mod_rttov_utils, only: find_ret_idx_rttov
      
      type(type_model), intent(inout) :: rttov_atm
      type(type_s3com), target, intent(inout) :: s3com
      integer(kind=4), dimension(:), allocatable :: idx_ret
      integer(kind=4) :: npoints, nlayers, idx, ipoint
      
      npoints = rttov_atm%npoints
      nlayers = rttov_atm%nlayers
      
      idx_ret = find_ret_idx_rttov(s3com)
      
      do ipoint = 1, size(idx_ret)
         
         idx = idx_ret(ipoint)
         
         if (s3com%ret%lwp(idx) .gt. 1e-4) then
            
            s3com%ret%lwc(idx,s3com%ret%ztop_liq_idx(idx):s3com%ret%zbase_liq_idx(idx)) = s3com%ret%lwp(idx) / &
            sum(rttov_atm%dz(idx,s3com%ret%ztop_liq_idx(idx):s3com%ret%zbase_liq_idx(idx)))
            
            s3com%ret%clc(idx,s3com%ret%ztop_liq_idx(idx):s3com%ret%zbase_liq_idx(idx)) = 1.0_wp
            
            s3com%ret%cod(idx) = s3com%ret%Xi(idx,1)
            
            s3com%ret%reff(idx,s3com%ret%ztop_liq_idx(idx):s3com%ret%zbase_liq_idx(idx)) = s3com%ret%Xi(idx,2)
            
         endif
         
      enddo
      
      rttov_atm%lwc  = s3com%ret%lwc
      rttov_atm%clc  = s3com%ret%clc
      rttov_atm%reff = s3com%ret%reff
      rttov_atm%cod  = s3com%ret%cod
      
   end subroutine init_cloud_prof
   ! ----------------------------------------------------------------------------------------------------------------------------
   
end module mod_ret_cloud_model

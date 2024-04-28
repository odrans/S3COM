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
   
   use s3com_types,  only: type_model, type_s3com, type_cld
   use s3com_config, only: wp, wpi, iwp_threshold, lwp_lay_threshold, cw, k, rholiq, pi
   
   implicit none
   
   private
   public :: init_cloud_z, init_cloud_prof
   
contains
   
   ! ============================================================================================================================
   !> @brief Calculate the cloud top and base altitudes and the liquid water path of a single-layer homogeneous liquid cloud
   !! @detail This subroutine determines the type of the cloud scene, extracts the cloud top and base altitudes of a single-layer
   !! homogeneous liquid cloud and calculates its liquid water path.
   !! @param[in] rttov_atm   rttov atmospheric structure
   !! @param[inout] s3com    s3com retrieval structure
   subroutine init_cloud_z(rttov_atm, s3com)
      
      ! Input/output
      type(type_model), intent(in) :: rttov_atm
      type(type_s3com), intent(inout) :: s3com
      
      ! Internal
      real(wp), dimension(rttov_atm%npoints,rttov_atm%nlayers) :: diff_liq
      integer(wpi) :: npoints, nlayers, ipoint, ilayer
      
      ! Get the number of points and layers from rttov_atm
      npoints = rttov_atm%npoints
      nlayers = rttov_atm%nlayers
      
      diff_liq = 0.0_wp
      
      ! Loop over each point
      do ipoint = 1, npoints
         
         ! Determine the type of the cloud scene (no cloud, single-layer cloud, multi-layer clouds)
         ! ----------------------------------------------------------------------------------------------------------------------
         ! Calculate the ice water path of the atmospheric column
         s3com%ret%iwp(ipoint) = 0.0_wp
         
         do ilayer = 1, nlayers
            s3com%ret%iwp(ipoint) = s3com%ret%iwp(ipoint) + rttov_atm%iwc(ipoint,ilayer) * rttov_atm%dz(ipoint,ilayer)
         enddo
         
         ! Calculate the difference of the cloud liquid water content between two successive layers
         do ilayer = 1, nlayers-1
            diff_liq(ipoint,ilayer) = rttov_atm%lwc(ipoint,ilayer+1) - rttov_atm%lwc(ipoint,ilayer)
         enddo
         
         ! Count the number of slope changes in the vertical profile of the cloud liquid water content
         ! If n_slope_liq = O --> no liquid cloud or cloud-free
         ! If n_slope_liq = 1 --> homogeneous single-layer liquid cloud
         ! If n_slope_liq > 1 --> heterogeneous single-layer liquid cloud and/or multi-layer conditions
         s3com%ret%n_slope_liq(ipoint) = 0
         
         do ilayer = 1, nlayers-2
            if (diff_liq(ipoint,ilayer) .gt. 0._wp .and. diff_liq(ipoint,ilayer+1) .lt. 0._wp) then
               s3com%ret%n_slope_liq(ipoint) = s3com%ret%n_slope_liq(ipoint) + 1
            endif
         enddo
         ! ----------------------------------------------------------------------------------------------------------------------
         
         ! Extract the cloud top and base and calculate the liquid water path of a single-layer homogeneous liquid cloud only
         ! ----------------------------------------------------------------------------------------------------------------------
         if (s3com%ret%n_slope_liq(ipoint) .eq. 1 .and. s3com%ret%iwp(ipoint) .lt. iwp_threshold) then
            
            ! Loop over each layer
            do ilayer = 1, nlayers
               
               ! If the liquid water path of the current layer is above the threshold value,
               ! set the index of the layer as the top of the liquid cloud layer and exit
               if (rttov_atm%lwc(ipoint,ilayer)*rttov_atm%dz(ipoint,ilayer) .gt. lwp_lay_threshold) then
                  s3com%ret%ztop_liq_idx(ipoint) = ilayer
                  exit
               endif
               
            enddo
            
            ! Loop over each layer in reverse order
            do ilayer = nlayers, 1, -1
               
               ! If the liquid water path of the current layer is above the threshold value,
               ! set the index of the layer as the base of the liquid cloud layer and exit
               if (rttov_atm%lwc(ipoint,ilayer)*rttov_atm%dz(ipoint,ilayer) .gt. lwp_lay_threshold) then
                  s3com%ret%zbase_liq_idx(ipoint) = ilayer+1
                  exit
               endif
               
            enddo
            
            ! If both top and base of the liquid cloud layer are found,
            ! calculate the corresponding heights and liquid water path of the cloud
            if (s3com%ret%ztop_liq_idx(ipoint) .gt. 0 .and. s3com%ret%zbase_liq_idx(ipoint) .gt. 0) then
               s3com%ret%ztop_liq(ipoint)  = rttov_atm%z(ipoint,s3com%ret%ztop_liq_idx(ipoint))
               s3com%ret%zbase_liq(ipoint) = rttov_atm%z(ipoint,s3com%ret%zbase_liq_idx(ipoint))
               
               s3com%ret%lwp(ipoint) = 0.0_wp
               do ilayer = s3com%ret%zbase_liq_idx(ipoint), s3com%ret%ztop_liq_idx(ipoint), -1
                  s3com%ret%lwp(ipoint) = s3com%ret%lwp(ipoint) + rttov_atm%lwc(ipoint,ilayer) * rttov_atm%dz(ipoint,ilayer)
               enddo
            endif
            
            !write(6,*) "IWP =", s3com%ret%iwp(ipoint), "kg/m2"
            !write(6,*) "LWP =", s3com%ret%lwp(ipoint), "kg/m2"
            !write(6,*) "Z_top =", s3com%ret%ztop_liq_idx(ipoint), "; ", "Z_base =", s3com%ret%zbase_liq_idx(ipoint)
            
         endif
         ! ----------------------------------------------------------------------------------------------------------------------
         !write(6,*) "ipoint:", ipoint, "n_slope:", s3com%ret%n_slope_liq
      enddo
      
   end subroutine init_cloud_z
   ! ============================================================================================================================
   
   ! ============================================================================================================================
   !> @brief Define the microphysics of the cloud model
   !! @detail This subroutine defines the microphysics of the cloud model. Currently, two cloud models can be used for the 
   !! retrieval algorithm:
   !! - the adiabatic cloud model;
   !! - the homogeneous cloud model where all properties are vertically homogeneous (deprecated for now).
   !! @param[inout] rttov_atm   rttov atmospheric structure
   !! @param[inout] s3com       s3com retrieval structure
   subroutine init_cloud_prof(rttov_atm, s3com)
      
      use mod_rttov_utils, only: find_ret_idx_rttov
      
      ! Input/output
      type(type_model), intent(inout) :: rttov_atm
      type(type_s3com), intent(inout) :: s3com
      
      ! Internal
      integer(wpi), dimension(:), allocatable :: idx_ret
      integer(wpi) :: idx, ipoint, ilayer
      logical :: adiab_cld_model ! If true, uses the adiabatic cloud model
      
      ! Get the index of single-layer homogeneous liquid cloud pixels
      idx_ret = find_ret_idx_rttov(s3com)
      
      adiab_cld_model = .true.
      
      if (adiab_cld_model) then
      
      ! Adiabatic cloud model
      ! -------------------------------------------------------------------------------------------------------------------------
      ! Loop over each single-layer homogeneous liquid cloud pixel
      do ipoint = 1, size(idx_ret)
         
         idx = idx_ret(ipoint)
         
         ! Adiabatic liquid water path
         s3com%ret%lwp_ad(idx) = s3com%ret%Xi(idx,1)
         !write(6,*) "LWP_ad =", s3com%ret%lwp_ad(idx)
         
         ! Cloud geometric height
         s3com%ret%H(idx) = s3com%ret%ztop_liq(idx) - s3com%ret%zbase_liq(idx)
         
         ! Adiabacity factor
         s3com%ret%fad(idx) = (2.0_wp*s3com%ret%lwp_ad(idx))/(cw*s3com%ret%H(idx)**2.0_wp)
         !write(6,*) "f_ad =", s3com%ret%fad(idx)
         
         ! Cloud droplet number concentration at cloud top - CDNC_top
         s3com%ret%cdnc_ad(idx,s3com%ret%ztop_liq_idx(idx)) = s3com%ret%Xi(idx,2)
         
         ! Loop over each layer between cloud top and cloud base
         do ilayer = s3com%ret%ztop_liq_idx(idx), s3com%ret%zbase_liq_idx(idx)
            
            ! Adiabatic liquid water content
            s3com%ret%lwc_ad(idx,ilayer) = s3com%ret%fad(idx)*cw*(rttov_atm%z(idx,ilayer) - &
                                           rttov_atm%z(idx,s3com%ret%zbase_liq_idx(idx)))
            
         enddo
         
         ! Calculate corresponding liquid water path
         s3com%ret%lwp(idx) = 0.0_wp
         
         do ilayer = s3com%ret%ztop_liq_idx(idx), s3com%ret%zbase_liq_idx(idx)
            s3com%ret%lwp(idx) = s3com%ret%lwp(idx) + s3com%ret%lwc_ad(idx,ilayer)*rttov_atm%dz(idx,ilayer)
         enddo
         
         !write(6,*) "LWP_calc =", s3com%ret%lwp(idx)
         
         ! Correct the adiabatic liquid water content to get the right liquid water path
         do ilayer = s3com%ret%ztop_liq_idx(idx), s3com%ret%zbase_liq_idx(idx)
            s3com%ret%lwc_corr(idx,ilayer) = s3com%ret%lwc_ad(idx,ilayer)*(s3com%ret%lwp_ad(idx)/s3com%ret%lwp(idx))
         enddo
         
         ! Calculate corresponding liquid water path
         s3com%ret%lwp(idx) = 0.0_wp
         
         do ilayer = s3com%ret%ztop_liq_idx(idx), s3com%ret%zbase_liq_idx(idx)
            s3com%ret%lwp(idx) = s3com%ret%lwp(idx) + s3com%ret%lwc_corr(idx,ilayer)*rttov_atm%dz(idx,ilayer)
         enddo
         
         !write(6,*) "LWP_calc_corr =", s3com%ret%lwp(idx)
         
         ! Loop over each layer between cloud top and cloud base
         do ilayer = s3com%ret%ztop_liq_idx(idx), s3com%ret%zbase_liq_idx(idx)
            
            ! Cloud fraction
            s3com%ret%clc(idx,ilayer) = 1.0_wp
            
            ! Adiabatic cloud droplet number concentration
            s3com%ret%cdnc_ad(idx,ilayer) = s3com%ret%cdnc_ad(idx,s3com%ret%ztop_liq_idx(idx))
            
            ! Adiabatic cloud droplet effective radius
            s3com%ret%re_ad(idx,ilayer) = ((3.0_wp*s3com%ret%lwc_corr(idx,ilayer)) / &
            (4.0_wp*pi*k*rholiq*s3com%ret%cdnc_ad(idx,s3com%ret%ztop_liq_idx(idx))))**(1.0_wp/3.0_wp)
            s3com%ret%re_ad(idx,ilayer) = s3com%ret%re_ad(idx,ilayer) * 1.0E+06_wp !< conversion from m to um (default input in RTTOV)
            
            !write(6,*) ilayer, rttov_atm%z(idx,ilayer), s3com%ret%cdnc_ad(idx,ilayer), s3com%ret%lwc_corr(idx,ilayer), &
            !           s3com%ret%re_ad(idx,ilayer)
            
         enddo
         
      enddo
      ! -------------------------------------------------------------------------------------------------------------------------
      
      rttov_atm%clc  = s3com%ret%clc
      rttov_atm%lwc  = s3com%ret%lwc_corr
      rttov_atm%reff = s3com%ret%re_ad
      rttov_atm%cdnc = s3com%ret%cdnc_ad
      rttov_atm%lwp  = s3com%ret%lwp_ad
      
      else
      
      ! Homogeneous vertical cloud model (deprecated for now)
      ! -------------------------------------------------------------------------------------------------------------------------
      ! Loop over each pixel containing liquid clouds
      do ipoint = 1, size(idx_ret)
         
         idx = idx_ret(ipoint)
         
         ! Liquid water path
         s3com%ret%lwp_hom(idx) = s3com%ret%Xi(idx,1)
         
         ! Cloud droplet number concentration at cloud top
         s3com%ret%cdnc_hom(idx,s3com%ret%ztop_liq_idx(idx)) = s3com%ret%Xi(idx,2)
         
         ! Loop over each layer between cloud base and cloud top
         do ilayer = s3com%ret%zbase_liq_idx(idx), s3com%ret%ztop_liq_idx(idx), -1
            
            ! Cloud fraction
            s3com%ret%clc(idx,ilayer) = 1.0_wp
            
            ! Cloud droplet number concentration
            s3com%ret%cdnc_hom(idx,ilayer) = s3com%ret%cdnc_hom(idx,s3com%ret%ztop_liq_idx(idx))
            
            ! Liquid water content
            s3com%ret%lwc_hom(idx,ilayer) = s3com%ret%lwp_hom(idx) / &
            sum(rttov_atm%dz(idx,s3com%ret%ztop_liq_idx(idx):s3com%ret%zbase_liq_idx(idx)))
            
            ! Cloud droplet effective radius
            s3com%ret%re_hom(idx,ilayer) = s3com%ret%re_hom(idx,s3com%ret%ztop_liq_idx(idx))
            
         enddo
         
      enddo
      
      rttov_atm%clc  = s3com%ret%clc
      rttov_atm%lwc  = s3com%ret%lwc_hom
      rttov_atm%reff = s3com%ret%re_hom
      rttov_atm%cdnc = s3com%ret%cdnc_hom
      rttov_atm%lwp  = s3com%ret%lwp_hom
      
      endif
      ! -------------------------------------------------------------------------------------------------------------------------
      
   end subroutine init_cloud_prof
   ! ============================================================================================================================
   
end module mod_ret_cloud_model

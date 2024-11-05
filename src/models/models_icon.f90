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

!> @brief Ensemble of subroutines useful to load and process ICON model simulations
module mod_icon
   
   use s3com_types,     only: wp, wpi, type_icon
   use mod_io_icon,     only: icon_read
   use mod_io_utils,    only: extract_coordinates
   use s3com_config,    only: rd, rv, Q_ext, pi, rholiq, cw, k, lwp_threshold, lwp_lay_threshold, iwp_threshold
   use mod_cld_phys_sb, only: re_sb
   
   implicit none
   
   private
   public :: icon_load, icon_free
   
contains
   
   ! ============================================================================================================================
   !> @brief Load ICON model simulations
   !! @details This subroutine loads ICON model simulations from a NetCDF file and processes them to obtain the required 
   !! variables for S3COM.
   !! @param[in] fname   input NetCDF file name containing ICON simulations
   !! @param[out] icon   ICON data structure
   subroutine icon_load(fname, icon)
      
      ! Input
      character(len=256), intent(in) :: fname
      
      ! Input/output
      type(type_icon), intent(out) :: icon
      
      ! Internal
      integer(wpi) :: nlayers, npoints
      
      ! Extract the number of vertical layers and grid points in the input files
      call extract_coordinates(fname, nlayers, npoints)
      
      ! Initialize the ICON array
      call icon_init(npoints, nlayers, icon)
      
      ! Read input netcdf file containing ICON outputs
      call icon_read(fname, icon)
      
      ! Post-process ICON data
      call icon_process(icon)
      
   end subroutine icon_load
   ! ============================================================================================================================
   
   ! ============================================================================================================================
   !> @brief Initialize ICON data structure
   !! @details This subroutine initializes the ICON data structure by allocating the required arrays.
   !! All arrays are initialized to zero.
   !! @param[in] npoints   number of grid points
   !! @param[in] nlayers   number of vertical layers
   !! @param[inout] icon   ICON data structure
   subroutine icon_init(npoints, nlayers, icon)
      
      ! Input
      type(type_icon) :: icon
      integer(wpi), intent(in) :: npoints, nlayers
      
      ! Internal
      integer(wpi) :: nlevels
      
      ! Number of atmospheric levels
      nlevels = nlayers + 1
      
      ! Get the number of pixels, atmospheric layers and atmospheric levels
      icon%npoints = npoints
      icon%nlayers = nlayers
      icon%nlevels = nlevels
      
      ! 1D variables
      allocate(icon%height(nlayers), source = 0)
      allocate(icon%height_2(nlevels), source = 0)
      
      ! 2D variables (integer)
      allocate(icon%ztop_liq_idx(npoints), source = 0)
      allocate(icon%zbase_liq_idx(npoints), source = 0)
      
      ! 2D variables (real)
      allocate(icon%lon(npoints), source = 0._wp)
      allocate(icon%lat(npoints), source = 0._wp)
      allocate(icon%lon_orig(npoints), source = 0._wp)
      allocate(icon%lat_orig(npoints), source = 0._wp)
      allocate(icon%point_orig(npoints), source = 0)
      allocate(icon%topography(npoints), source = 0._wp)
      allocate(icon%landmask(npoints), source = 0._wp)
      allocate(icon%ps(npoints), source = 0._wp)
      allocate(icon%ts(npoints), source = 0._wp)
      allocate(icon%t_2m(npoints), source = 0._wp)
      allocate(icon%q_2m(npoints), source = 0._wp)
      allocate(icon%u_10m(npoints), source = 0._wp)
      allocate(icon%v_10m(npoints), source = 0._wp)
      allocate(icon%lwp(npoints), source = 0._wp)
      allocate(icon%iwp(npoints), source = 0._wp)
      allocate(icon%cod(npoints), source = 0._wp)
      allocate(icon%ztop_liq(npoints), source = 0._wp)
      allocate(icon%zbase_liq(npoints), source = 0._wp)
      allocate(icon%reff_top(npoints), source = 0._wp)
      allocate(icon%cdnc_top(npoints), source = 0._wp)
      allocate(icon%lwp_stratocumulus(npoints), source = 0._wp)
      allocate(icon%lwp_stratocumulus_filter(npoints), source = 0._wp)
      allocate(icon%cod_stratocumulus(npoints), source = 0._wp)
      
      ! 3D variables on atmospheric levels
      allocate(icon%z_ifc(npoints, nlevels), source = 0._wp)
      allocate(icon%p_ifc(npoints, nlevels), source = 0._wp)
      allocate(icon%t_ifc(npoints, nlevels), source = 0._wp)
      allocate(icon%q_ifc(npoints, nlevels), source = 0._wp)
      
      ! 3D variables in atmospheric layers
      allocate(icon%z(npoints, nlayers), source = 0._wp)
      allocate(icon%p(npoints, nlayers), source = 0._wp)
      allocate(icon%t(npoints, nlayers), source = 0._wp)
      allocate(icon%q(npoints, nlayers), source = 0._wp)
      allocate(icon%clc(npoints, nlayers), source = 0._wp)
      allocate(icon%clw(npoints, nlayers), source = 0._wp)
      allocate(icon%cli(npoints, nlayers), source = 0._wp)
      allocate(icon%qnc(npoints, nlayers), source = 0._wp)
      allocate(icon%qr(npoints, nlayers), source = 0._wp)
      allocate(icon%qs(npoints, nlayers), source = 0._wp)
      allocate(icon%dz(npoints, nlayers), source = 0._wp)
      allocate(icon%rho(npoints, nlayers), source = 0._wp)
      allocate(icon%tv(npoints, nlayers), source = 0._wp)
      allocate(icon%lwc(npoints, nlayers), source = 0._wp)
      allocate(icon%iwc(npoints, nlayers), source = 0._wp)
      allocate(icon%cdnc(npoints, nlayers), source = 0._wp)
      allocate(icon%reff(npoints, nlayers), source = 0._wp)
      allocate(icon%beta_ext(npoints, nlayers), source = 0._wp)
      
   end subroutine icon_init
   ! ============================================================================================================================
   
   ! ============================================================================================================================
   !> @brief Process ICON data structure
   !! @details This subroutine processes the ICON data structure to obtain the required variables for S3COM.
   !! @note Current processing steps:
   !! - The layer depth is computed as the difference between two successive levels enclosing the layer.
   !! - The moist air density is computed from the pressure and virtual temperature. Used to convert cloud properties from 
   !! @units{kg kg-1} to @units{kg m-3}.
   !! - The effective radius is computed from the water content and number concentration consistently with the Seifert and 
   !! Beheng (2005) scheme used in ICON.
   !! - All quantities are in SI units except for the effective radius which is in @units{um}.
   !! @param[inout] icon   ICON data structure
   subroutine icon_process(icon)
      
      ! Input/output
      type(type_icon), intent(inout) :: icon
      
      ! Internal
      integer(wpi) :: i, j, nlevels, nlayers, npoints
      integer(wpi), dimension(icon%npoints) :: n_slope_liq
      real(wp), dimension(icon%npoints,icon%nlayers) :: diff_liq
      
      nlevels = icon%nlevels
      nlayers = icon%nlayers
      npoints = icon%npoints
      
      ! Atmospheric structure
      ! -------------------------------------------------------------------------------------------------------------------------
      ! Layer depth
      do i = 1, npoints
         do j = 1, nlayers
            icon%dz(i,j) = icon%z_ifc(i,j) - icon%z_ifc(i,j+1)
         enddo
      enddo
      ! -------------------------------------------------------------------------------------------------------------------------
      
      ! Atmospheric moisture
      ! -------------------------------------------------------------------------------------------------------------------------
      ! Virtual temperature @units{K}
      icon%tv = icon%t * (1 + 0.608_wp * icon%q)
      
      ! Moist air density @units{kg m-3}
      icon%rho = icon%p / (rd * icon%tv)
      
      ! Liquid water content @units{kg m-3}
      icon%lwc = icon%clw * icon%rho
      
      ! Ice water content @units{kg m-3}
      icon%iwc = icon%cli * icon%rho
      
      ! Liquid water path @units{kg m-2}
      do i = 1, npoints
         icon%lwp(i) = 0.0_wp
         do j = 1, nlayers
            icon%lwp(i) = icon%lwp(i) + icon%lwc(i,j) * icon%dz(i,j)
         enddo
      enddo
      
      ! Ice water path @units{kg m-2}
      do i = 1, npoints
         icon%iwp(i) = 0.0_wp
         do j = 1, nlayers
            icon%iwp(i) = icon%iwp(i) + icon%iwc(i,j) * icon%dz(i,j)
         enddo
      enddo
      
      ! Saturation vapour pressure of liquid water (Pa); Murphy and Koop (2005)
      !icon%es_w = exp(54.842763 - 6763.22 / icon%t - 4.210 * log(icon%t) + 0.000367 * icon%t + &
      !            tanh(0.0415 * (icon%t - 218.8)) * (53.878 - 1331.22 / icon%t - 9.44523 * log(icon%t) + 0.014025 * icon%t))
      
      ! Saturation vapour pressure of ice water (Pa); Murphy and Koop (2005)
      !icon%es_i = exp(9.550426 - 5723.265 / icon%t + 3.53068 * log(icon%t) - 0.00728332 * icon%t)
      ! -------------------------------------------------------------------------------------------------------------------------
      
      ! Cloud properties
      ! -------------------------------------------------------------------------------------------------------------------------
      ! Cloud droplet number concentration @units{\# m-3}
      icon%cdnc = icon%qnc * icon%rho
      
      ! Cloud droplet effective radius @units{um}
      do i = 1, npoints
         do j = 1, nlayers
            if(icon%cdnc(i,j) .gt. 0.0_wp) then
               icon%reff(i,j) = re_sb(icon%lwc(i,j), icon%cdnc(i,j), 1)
            endif
            !write(*,*) j, ";", icon%z(i,j), ";", icon%lwc(i,j), ";", icon%reff(i,j)*1.0E+06_wp, ";", icon%cdnc(i,j)*1.0E-06_wp
         enddo
      enddo
      
      icon%reff = icon%reff * 1.0E+06_wp !< Conversion from @units{m} to @units{um} (default input in RTTOV)
      
      ! Cloud droplet extinction coefficient @units{m-1}
      do i = 1, npoints
         do j = 1, nlayers
            if(icon%reff(i,j) .gt. 0.0_wp) then
               icon%beta_ext(i,j) = (3.0_wp/4.0_wp)*(Q_ext/rholiq)*(icon%lwc(i,j)/(icon%reff(i,j)*1.0E-06_wp))
            endif
         enddo
      enddo
      
      ! Cloud optical depth @units{-}
      do i = 1, npoints
         icon%cod(i) = 0.0_wp
         do j = 1, nlayers
            icon%cod(i) = icon%cod(i) + icon%beta_ext(i,j) * icon%dz(i,j)
         enddo
      enddo
      ! -------------------------------------------------------------------------------------------------------------------------
      
      ! Cloud top properties of single-layer homogeneous liquid clouds only
      ! -------------------------------------------------------------------------------------------------------------------------
      ! Cloud droplet effective radius (um) and number concentration (m-3)
      do i = 1, npoints
         
         ! Calculate the difference of the liquid water content between two successive layers
         do j = 1, nlayers-1
            diff_liq(i,j) = icon%lwc(i,j+1) - icon%lwc(i,j)
         enddo
         
         n_slope_liq(i) = 0
         
         ! Count the number of slope changes in the liquid water content vertical profile
         ! If n_slope_liq = O --> cloud-free pixel
         ! If n_slope_liq = 1 --> single-layer homogeneous liquid cloud
         ! If n_slope_liq > 1 --> single-layer heterogeneous liquid cloud and/or multi-layer conditions
         do j = 1, nlayers-2
            if (diff_liq(i,j) .gt. 0.0_wp .and. diff_liq(i,j+1) .lt. 0.0_wp) then
               n_slope_liq(i) = n_slope_liq(i) + 1
            endif
         enddo
         
         ! For single-layer homogeneous liquid clouds only
         if (n_slope_liq(i) .eq. 1) then
            
            ! Loop over each layer
            do j = 1, nlayers
               
               ! Determination of the cloud top index
               if (icon%lwc(i,j)*icon%dz(i,j) .gt. lwp_lay_threshold) then
                  icon%ztop_liq_idx(i) = j
                  !write(6,*) icon%ztop_liq_idx(i)
                  exit
               endif
               
            enddo
            
            ! Loop over each layer in reverse order
            do j = nlayers, 1, -1
               
               ! Determination of the cloud base index
               if (icon%lwc(i,j)*icon%dz(i,j) .gt. lwp_lay_threshold) then
                  icon%zbase_liq_idx(i) = j+1
                  !write(6,*) icon%zbase_liq_idx(i)
                  exit
               endif
               
            enddo
            
            ! Calculation of the cloud liquid water path and cloud optical depth
            if (icon%ztop_liq_idx(i) .gt. 0 .and. icon%zbase_liq_idx(i) .gt. 0) then
               
               !icon%ztop_liq(i)  = icon%z(i,icon%ztop_liq_idx(i))
               !icon%zbase_liq(i) = icon%z(i,icon%zbase_liq_idx(i))
               
               icon%lwp_stratocumulus(i) = 0.0_wp
               icon%cod_stratocumulus(i) = 0.0_wp
               
               ! Loop between cloud base and cloud top
               do j = icon%zbase_liq_idx(i), icon%ztop_liq_idx(i), -1
                  icon%lwp_stratocumulus(i) = icon%lwp_stratocumulus(i) + icon%lwc(i,j) * icon%dz(i,j)
                  icon%cod_stratocumulus(i) = icon%cod_stratocumulus(i) + icon%beta_ext(i,j) * icon%dz(i,j)
               enddo
               
            endif
            
         endif
         
         ! For optically thick single-layer homogeneous liquid clouds (no multi-layer conditions with other liquid or ice clouds)
         if (icon%lwp_stratocumulus(i) .gt. lwp_threshold .and. & !< Distinction between a cloud-free pixel and a liquid cloud pixel
            n_slope_liq(i) .eq. 1 .and.                         & !< Single-layer homogeneous liquid cloud
            icon%cod_stratocumulus(i) .gt. 5.0_wp .and.         & !< Optically thick liquid cloud
            icon%iwp(i) .lt. iwp_threshold) then                  !< No ice cloud present
            
            icon%lwp_stratocumulus_filter(i) = icon%lwp_stratocumulus(i)
            
            ! Determination of the cloud top properties
            icon%cdnc_top(i) = icon%cdnc(i,icon%ztop_liq_idx(i))
            icon%reff_top(i) = icon%reff(i,icon%ztop_liq_idx(i))
            
            !write(*,*)
            !write(6,*) "lon =", icon%lon(i), "lat =", icon%lat(i), "LWP_obs (kg/m2) =", icon%lwp_stratocumulus_filter(i), &
            !           "Re_ztop_obs (um) =", icon%reff_top(i), "CDNC_ztop_obs (m-3) =", icon%cdnc_top(i)
            !write(*,*) icon%lon(i), icon%lat(i), icon%lwp_stratocumulus_filter(i), icon%reff_top(i), icon%cdnc_top(i)
            
         endif
         
      enddo
      ! -------------------------------------------------------------------------------------------------------------------------
      
   end subroutine icon_process
   ! ============================================================================================================================
   
   ! ============================================================================================================================
   !> @brief Free ICON data structure
   !! @details This subroutine frees the ICON data structure by deallocating the required arrays.
   !! @param[inout] icon   ICON data structure
   subroutine icon_free(icon)
      
      type(type_icon), intent(inout) :: icon
      
      ! 1D variables
      deallocate(icon%height, icon%height_2)
      
      ! 2D variables
      deallocate(icon%ztop_liq_idx, icon%zbase_liq_idx)
      deallocate(icon%lon, icon%lat, icon%lon_orig, icon%lat_orig, icon%point_orig, icon%topography, icon%landmask, icon%ps, icon%ts, icon%t_2m, &
                 icon%q_2m, icon%u_10m, icon%v_10m, icon%cod, icon%lwp, icon%iwp, icon%ztop_liq, icon%zbase_liq, icon%reff_top, &
                 icon%cdnc_top, icon%lwp_stratocumulus, icon%lwp_stratocumulus_filter, icon%cod_stratocumulus)
      
      ! 3D variables on atmospheric levels
      deallocate(icon%z_ifc, icon%p_ifc, icon%t_ifc, icon%q_ifc)
      
      ! 3D variables in atmospheric layers
      deallocate(icon%z, icon%p, icon%t, icon%q, icon%clc, icon%clw, icon%cli, icon%qnc, icon%qr, icon%qs, icon%dz, icon%rho, &
                 icon%tv, icon%lwc, icon%iwc, icon%cdnc, icon%reff, icon%beta_ext)
      
   end subroutine icon_free
   ! ============================================================================================================================
   
end module mod_icon

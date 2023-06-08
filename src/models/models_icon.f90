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
   
   use s3com_types,     only: wp, type_icon
   use mod_io_icon,     only: icon_read
   use mod_io_utils,    only: extract_coordinates
   use s3com_config,    only: rd, rv, epsilon, mu, nu, a, b, Q_ext, rholiq, lwp_lay_threshold
   use mod_cld_phys_sb, only: re_sb
   
   implicit none
   
   private
   public :: icon_load, icon_free
   
contains
   
   ! ----------------------------------------------------------------------------------------------------------------------------
   !> @brief Load ICON model simulations
   !> @details This subroutine loads ICON model simulations from a NetCDF file and processes them to obtain the required variables for S3COM.
   !> @param[in] fname Input NetCDF file name containing ICON simulations
   !> @param[out] icon ICON data structure
   subroutine icon_load(fname, icon)
      
      ! Inputs
      character(len=256), intent(in) :: fname
      
      ! Inputs/Outputs
      type(type_icon), intent(out) :: icon
      
      ! Internal
      integer(kind=4) :: nlayers, npoints
      
      ! Extract the number of vertical layers and grid points in the input files
      call extract_coordinates(fname, nlayers, npoints)
      
      ! Initialize the ICON array
      call icon_init(npoints, nlayers, icon)
      
      ! Read input netcdf file containing ICON outputs
      call icon_read(fname, icon)
      
      ! Post-process ICON data
      call icon_process(icon)
      
   end subroutine icon_load
   ! ----------------------------------------------------------------------------------------------------------------------------
   
   ! ----------------------------------------------------------------------------------------------------------------------------
   !> @brief Process ICON data structure
   !> @details This subroutine processes the ICON data structure to obtain the required variables for S3COM.
   !> @note Current processing steps:
   ! - The layer depth is approximated to twice the difference between the interface and the layer center.
   ! - The moist air density is computed from the pressure and virtual temperature. Used to convert cloud properties from @units{kg/kg} to @units{kg/m3}.
   ! - The effective radius is computed from the water content and number concentration consistently with the Seifert and Beheng (2005) scheme used in ICON
   ! - All quantities in SI units except for the effective radius which is in @units{um}.
   !> @todo Why not compute dz from ifc(i+1) - ifc(i)?
   !> @param[inout] icon ICON data structure
   subroutine icon_process(icon)
      
      ! Inputs/outputs
      type(type_icon), intent(inout) :: icon
      
      ! Internal
      integer(kind=4) :: i, j, nlevels, nlayers, npoints
      
      nlevels = icon%nlevels
      nlayers = icon%nlayers
      npoints = icon%npoints
      
      ! Atmospheric structure
      ! -------------------------------------------------------------------------------------------------------------------------
      ! Layer depth (approximate!! why not ifc(i+1) - ifc(i)?)
      icon%dz(:,1:nlayers) = abs(icon%z_ifc(:,1:nlayers) - icon%z(:,1:nlayers)) * 2._wp
      
      ! Atmospheric moisture
      ! -------------------------------------------------------------------------------------------------------------------------
      ! Virtual temperature
      icon%tv = icon%t * (1 + 0.608 * icon%q)
      
      ! Moist air density (kg/m3)
      icon%rho = icon%p / (rd * icon%tv)
      
      ! Saturation vapour pressure of liquid water (Pa); Murphy and Koop (2005)
      !icon%es_w = exp(54.842763 - 6763.22 / icon%t - 4.210 * log(icon%t) + 0.000367 * icon%t + &
      !            tanh(0.0415 * (icon%t - 218.8)) * (53.878 - 1331.22 / icon%t - 9.44523 * log(icon%t) + 0.014025 * icon%t))
      
      ! Saturation vapour pressure of ice water (Pa); Murphy and Koop (2005)
      !icon%es_i = exp(9.550426 - 5723.265 / icon%t + 3.53068 * log(icon%t) - 0.00728332 * icon%t)
      
      ! Cloud properties
      ! -------------------------------------------------------------------------------------------------------------------------
      ! Cloud liquid water content (kg/m3)
      icon%lwc = icon%clw * icon%rho
      icon%lwc(:,:) = 0._wp            !for retrievals
      icon%lwc(:,111:122) = 2._wp * 1e-4 !for retrievals
      
      ! Cloud ice water content (kg/m3)
      icon%iwc = icon%cli * icon%rho
      
      ! Liquid water path (kg/m2)
      icon%lwp = sum(icon%lwc * icon%dz)
      
      ! Ice water path (kg/m2)
      icon%iwp = sum(icon%iwc * icon%dz)
      
      ! Cloud droplet number concentration (particules/m3)
      icon%cdnc = icon%qnc * icon%rho
      icon%cdnc(:,:) = 0._wp               !for retrievals
      icon%cdnc(:,111:122) = 16.586_wp * 1e6 !for retrievals
      
      ! Cloud fraction
      icon%clc(:,:) = 0._wp       !for retrievals
      icon%clc(:,111:122) = 1._wp !for retrievals
      
      ! Cloud liquid water effective radius (m)
      do i = 1, npoints
         do j = 1, nlayers
            if(icon%cdnc(i,j) .gt. 0._wp) then
               icon%reff(i,j) = re_sb(icon%lwc(i,j), icon%cdnc(i,j), 1)
            endif
         enddo
      enddo
      
      icon%reff = icon%reff * 1e6 !! m to um (default input in RTTOV)
      
      ! Cloud droplet extinction coefficient (m-1)
      do i = 1, npoints
         do j = 1, nlayers
            if(icon%reff(i,j) .gt. 0) then
               icon%beta_ext(i,j) = (3._wp/4._wp)*(Q_ext/rholiq)*(icon%lwc(i,j)/(icon%reff(i,j)*1e-6))
            endif
         enddo
      enddo
      
      ! Cloud optical depth
      do i = 1, npoints
         icon%cod(i) = 0._wp
         do j = 1, nlayers
            icon%cod(i) = icon%cod(i) + icon%beta_ext(i,j)*icon%dz(i,j)
         enddo
      enddo
      
      ! Cloud droplet effective radius (um) and number concentration (m-3) at cloud top
      do i = 1, npoints
         do j = 1, nlayers
            if (icon%lwc(i,j)*icon%dz(i,j) .gt. lwp_lay_threshold .and. icon%cod(i) .gt. 5._wp) then
               icon%ztop_liq_idx(i) = j
               icon%reff_top(i) = icon%reff(i,icon%ztop_liq_idx(i))
               icon%cdnc_top(i) = icon%cdnc(i,icon%ztop_liq_idx(i))
               exit
            endif
         enddo
      enddo
      
   end subroutine icon_process
   ! ----------------------------------------------------------------------------------------------------------------------------
   
   ! ----------------------------------------------------------------------------------------------------------------------------
   !> @brief Initialize ICON data structure
   !> @details This subroutine initializes the ICON data structure by allocating the required arrays.
   !> All arrays are initialized to zero.
   !> @todo Not all arrays are initialized to zero. Redo following the NWPSAF procedure.
   !> @param[in] npoints Number of grid points
   !> @param[in] nlayers Number of vertical layers
   !> @param[inout] icon ICON data structure
   subroutine icon_init(npoints, nlayers, icon)
      
      ! Inputs
      type(type_icon) :: icon
      integer(kind=4), intent(in) :: npoints, nlayers
      
      ! Internal
      integer(kind=4) :: nlevels
      
      nlevels = nlayers + 1
      
      icon%npoints = npoints
      icon%nlayers = nlayers
      icon%nlevels = nlevels
      
      allocate(icon%height(nlayers),       source = 0)
      allocate(icon%height_2(nlevels),     source = 0)
      allocate(icon%ztop_liq_idx(npoints), source = 0)
      
      ! 2D variables
      allocate(icon%lon(npoints),        source = 0._wp)
      allocate(icon%lat(npoints),        source = 0._wp)
      allocate(icon%lon_orig(npoints),   source = 0._wp)
      allocate(icon%lat_orig(npoints),   source = 0._wp)
      allocate(icon%topography(npoints), source = 0._wp)
      allocate(icon%landmask(npoints),   source = 0._wp)
      allocate(icon%ps(npoints),         source = 0._wp)
      allocate(icon%ts(npoints),         source = 0._wp)
      allocate(icon%t_2m(npoints),       source = 0._wp)
      allocate(icon%q_2m(npoints),       source = 0._wp)
      allocate(icon%u_10m(npoints),      source = 0._wp)
      allocate(icon%v_10m(npoints),      source = 0._wp)
      allocate(icon%lwp(npoints),        source = 0._wp)
      allocate(icon%iwp(npoints),        source = 0._wp)
      allocate(icon%cod(npoints),        source = 0._wp)
      allocate(icon%reff_top(npoints),   source = 0._wp)
      allocate(icon%cdnc_top(npoints),   source = 0._wp)
      
      ! 3D variables on atmospheric levels
      allocate(icon%z_ifc(npoints, nlevels), source = 0._wp)
      allocate(icon%p_ifc(npoints, nlevels), source = 0._wp)
      allocate(icon%t_ifc(npoints, nlevels), source = 0._wp)
      allocate(icon%q_ifc(npoints, nlevels), source = 0._wp)
      
      ! 3D variables in atmospheric layers
      allocate(icon%z(npoints, nlayers),        source = 0._wp)
      allocate(icon%p(npoints, nlayers),        source = 0._wp)
      allocate(icon%t(npoints, nlayers),        source = 0._wp)
      allocate(icon%q(npoints, nlayers),        source = 0._wp)
      allocate(icon%clc(npoints, nlayers),      source = 0._wp)
      allocate(icon%clw(npoints, nlayers),      source = 0._wp)
      allocate(icon%cli(npoints, nlayers),      source = 0._wp)
      allocate(icon%qnc(npoints, nlayers),      source = 0._wp)
      allocate(icon%qr(npoints, nlayers),       source = 0._wp)
      allocate(icon%qs(npoints, nlayers),       source = 0._wp)
      allocate(icon%dz(npoints, nlayers),       source = 0._wp)
      allocate(icon%rho(npoints, nlayers),      source = 0._wp)
      allocate(icon%tv(npoints, nlayers),       source = 0._wp)
      allocate(icon%lwc(npoints, nlayers),      source = 0._wp)
      allocate(icon%iwc(npoints, nlayers),      source = 0._wp)
      allocate(icon%cdnc(npoints, nlayers),     source = 0._wp)
      allocate(icon%reff(npoints, nlayers),     source = 0._wp)
      allocate(icon%beta_ext(npoints, nlayers), source = 0._wp)
      
   end subroutine icon_init
   ! ----------------------------------------------------------------------------------------------------------------------------
   
   ! ----------------------------------------------------------------------------------------------------------------------------
   !> @brief Free ICON data structure
   !> @details This subroutine frees the ICON data structure by deallocating the required arrays.
   !> @param[inout] icon ICON data structure
   subroutine icon_free(icon)
      
      type(type_icon), intent(inout) :: icon
      
      deallocate(icon%height, icon%height_2, icon%lon, icon%lat, icon%lon_orig, icon%lat_orig, &
                 icon%topography, icon%landmask, icon%ps, icon%ts, icon%t_2m, icon%q_2m, icon%u_10m, icon%v_10m, icon%cod, &
                 icon%p, icon%z, icon%z_ifc, icon%p_ifc, icon%t_ifc, icon%q_ifc, &
                 icon%t, icon%q, icon%clc, icon%clw, icon%cli, icon%qnc, icon%qr, icon%qs, icon%dz, &
                 icon%rho, icon%tv, icon%lwc, icon%iwc, icon%lwp, icon%iwp, icon%cdnc, icon%reff, icon%beta_ext, &
                 icon%ztop_liq_idx, icon%reff_top, icon%cdnc_top)
      
   end subroutine icon_free
   ! ----------------------------------------------------------------------------------------------------------------------------
   
end module mod_icon

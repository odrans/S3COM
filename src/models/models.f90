
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

!> @brief Ensemble of subroutines loading model data
module mod_models
   
   use s3com_types,    only: wp, type_icon, type_model, type_nml, type_nwpsaf
   use mod_icon,       only: icon_load, icon_free
   use mod_nwpsaf,     only: nwpsaf_load, nwpsaf_free
   use mod_utils_phys, only: solar_angles
   use mod_io_verbose, only: verbose_model
   
   implicit none
   
   private
   public :: models_load, models_free
   
contains
   
   ! ----------------------------------------------------------------------------------------------------------------------------
   !> @brief Load model simulations to be used in S3COM
   !> @details This subroutine loads model data from their respective netCDF files and stores them consistently in the `model` structure.
   !> This subroutine is currently set up for the ICON and NWPSAF models.
   !> @param[in] nml Namelist data, containing the name of the model and path to its simulation outputs
   !> @param[out] model Model data structure containing all necessary variables for S3COM
   subroutine models_load(nml, model)
      
      ! Inputs
      type(type_nml), intent(in) :: nml
      
      ! Output variables
      type(type_model), intent(out) :: model
      
      ! Internal
      type(type_icon) :: icon
      type(type_nwpsaf) :: nwpsaf
      integer(kind=4) :: nlayers, npoints
      
      ! Load input data
      select case (nml%model_name)
         case ("ICON")
            call icon_load(nml%fname_in, icon)
            npoints = icon%npoints
            nlayers = icon%nlayers
         
         case ("NWPSAF")
            call nwpsaf_load(nml%fname_in, nwpsaf)
            npoints = nwpsaf%npoints
            nlayers = nwpsaf%nlayers
      end select
      
      ! Initialize model array
      call models_init(model, npoints, nlayers)
      
      ! Set up model structure with input data
      select case (nml%model_name)
         case ("ICON")
            call models_setup_icon(model, icon)
         case ("NWPSAF")
            call models_setup_nwpsaf(model, nwpsaf)
      end select
      
      ! Set up solar angles corresponding to model data
      call models_setup_solar(model)
      
      ! Clear input data
      select case (nml%model_name)
         case ("ICON")
            call icon_free(icon)
         case ("NWPSAF")
            call nwpsaf_free(nwpsaf)
      end select
      
      call verbose_model(nml, model)
      
   end subroutine models_load
   ! ----------------------------------------------------------------------------------------------------------------------------
   
   ! ----------------------------------------------------------------------------------------------------------------------------
   !> @brief Initialize the `model` data structure
   !> @details This subroutine initializes the `model` structure by allocating the required arrays.
   !> All arrays are initialized to zero.
   !> @param[in] npoints Number of points in the model simulation
   !> @param[in] nlayers Number of layers in the model simulation
   !> @param[out] model Model data structure containing all necessary variables for S3COM
   subroutine models_init(model, npoints, nlayers)
      
      ! Inputs
      integer(kind = 4), intent(in) :: npoints, nlayers
      
      ! Outputs
      type(type_model), intent(out) :: model
      
      ! Internal
      integer(kind = 4) :: nlevels
      
      nlevels = nlayers + 1
      
      model%npoints   =  npoints
      model%nlayers   =  nlayers
      model%nlevels   =  nlevels
      
      model%date = (/0, 0, 0/)
      model%time = (/0, 0, 0/)
      
      allocate(model%height(nlayers), source = 0)
      allocate(model%height_2(nlevels), source = 0)
      allocate(model%point(npoints), source = 0)
      
      ! 2D fields
      allocate(model%lat(npoints), source = 0._wp)
      allocate(model%lon(npoints), source = 0._wp)
      allocate(model%lat_orig(npoints), source = 0._wp)
      allocate(model%lon_orig(npoints), source = 0._wp)
      allocate(model%topography(npoints), source = 0._wp)
      allocate(model%landmask(npoints), source = 0._wp)
      allocate(model%ps(npoints), source = 0._wp)
      allocate(model%ts(npoints), source = 0._wp)
      allocate(model%t_2m(npoints), source = 0._wp)
      allocate(model%q_2m(npoints), source = 0._wp)
      allocate(model%u_10m(npoints), source = 0._wp)
      allocate(model%v_10m(npoints), source = 0._wp)
      allocate(model%iwp(npoints), source = 0._wp)
      allocate(model%lwp(npoints), source = 0._wp)
      allocate(model%cod(npoints), source = 0._wp)
      allocate(model%reff_top(npoints), source = 0._wp)
      allocate(model%cdnc_top(npoints), source = 0._wp)
      allocate(model%sunzenangle(npoints), source = 0._wp)
      allocate(model%sunazangle(npoints), source = 0._wp)
      
      ! 3D fields at atmospheric levels
      allocate(model%zh(npoints, nlevels), source = 0._wp)
      allocate(model%p(npoints, nlevels), source = 0._wp)
      allocate(model%t(npoints, nlevels), source = 0._wp)
      allocate(model%q(npoints, nlevels), source = 0._wp)
      allocate(model%co2(npoints, nlevels), source = 0._wp)
      allocate(model%ch4(npoints, nlevels), source = 0._wp)
      allocate(model%n2o(npoints, nlevels), source = 0._wp)
      allocate(model%s2o(npoints, nlevels), source = 0._wp)
      allocate(model%co(npoints, nlevels), source = 0._wp)
      allocate(model%o3(npoints, nlevels), source = 0._wp)
      
      ! 3D fields in atmospheric layers
      allocate(model%z(npoints, nlayers), source = 0._wp)
      allocate(model%dz(npoints, nlayers), source = 0._wp)
      allocate(model%clc(npoints, nlayers), source = 0._wp)
      allocate(model%iwc(npoints, nlayers), source = 0._wp)
      allocate(model%lwc(npoints, nlayers), source = 0._wp)
      allocate(model%reff(npoints, nlayers), source = 0._wp)
      allocate(model%cdnc(npoints, nlayers), source = 0._wp)
      allocate(model%beta_ext(npoints, nlayers), source = 0._wp)
      
   end subroutine models_init
   ! ----------------------------------------------------------------------------------------------------------------------------
   
   ! ----------------------------------------------------------------------------------------------------------------------------
   !> @brief Free the `model` data structure
   !> @details This subroutine frees the `model` structure by deallocating the arrays.
   !> @param[inout] model Model data structure containing all necessary variables for S3COM
   subroutine models_free(model)
      
      ! Inputs/outputs
      type(type_model), intent(inout) :: model
      
      deallocate(model%height, model%height_2, model%point)
      
      ! 2D fields
      deallocate(model%lat, model%lon, model%lat_orig, model%lon_orig, model%topography, model%landmask, model%ps, model%ts, &
                 model%t_2m, model%q_2m, model%u_10m, model%v_10m, model%iwp, model%lwp, model%cod, model%cdnc_top, &
                 model%reff_top, model%sunzenangle, model%sunazangle)
      
      ! 3D fields at atmospheric levels
      deallocate(model%zh, model%p, model%t, model%q, model%co2, model%ch4, model%n2o, model%s2o, model%co, model%o3)
      
      ! 3D fields in atmospheric layers
      deallocate(model%z, model%dz, model%clc, model%reff, model%cdnc, model%iwc, model%lwc, model%beta_ext)
      
   end subroutine models_free
   ! ----------------------------------------------------------------------------------------------------------------------------
   
   ! ----------------------------------------------------------------------------------------------------------------------------
   !> @brief Set up the `model` data structure by importing the ICON data
   !> @param[inout] model Model data structure containing all necessary variables for S3COM
   !> @param[in] icon structure containing the ICON data
   subroutine models_setup_icon(model, icon)
      
      ! Inputs
      type(type_icon), intent(in) :: icon
      
      ! Inputs/outputs
      type(type_model), intent(inout) :: model
      
      ! Internal
      real(wp) :: hour, minute, sec
      character(len=8) :: icon_date
      
      ! Extract information from the ICON time (%Y%m%d.%f)
      ! -----------------------------------------------------------------------------
      hour = real(icon%time - int(icon%time), wp) * 24._wp; model%time(1) = int(hour)
      minute = (hour - int(hour)) * 60; model%time(2) = int(minute)
      sec = (minute - int(minute)) * 60; model%time(3) = int(sec)
      
      write(icon_date, "(I8)") int(icon%time)
      read(icon_date(1:4), "(I4)") model%date(3)
      read(icon_date(5:6), "(I2)") model%date(2)
      read(icon_date(7:8), "(I2)") model%date(1)
      ! -----------------------------------------------------------------------------
      
      model%nlat     = icon%nlat
      model%nlon     = icon%nlon
      model%mode     = icon%mode
      
      model%npoints  = icon%npoints
      model%nlevels  = icon%nlevels
      model%nlayers  = icon%nlayers
      
      model%height   = icon%height
      model%height_2 = icon%height_2
      
      ! 2D fields
      model%lat        = icon%lat
      model%lon        = icon%lon
      model%lat_orig   = icon%lat_orig
      model%lon_orig   = icon%lon_orig
      model%topography = icon%topography
      model%landmask   = icon%landmask
      model%ps         = icon%ps
      model%ts         = icon%ts
      model%t_2m       = icon%t_2m
      model%q_2m       = icon%q_2m
      model%u_10m      = icon%u_10m
      model%v_10m      = icon%v_10m
      model%lwp        = icon%lwp
      model%iwp        = icon%iwp
      model%cod        = icon%cod
      model%cdnc_top   = icon%cdnc_top
      model%reff_top   = icon%reff_top
      
      ! 3D fields on atmospheric levels
      model%p  = icon%p_ifc
      model%t  = icon%t_ifc
      model%q  = icon%q_ifc
      model%zh = icon%z_ifc
      
      ! 3D fields in atmospheric layers
      model%z        = icon%z
      model%dz       = icon%dz
      model%clc      = icon%clc
      model%iwc      = icon%iwc
      model%lwc      = icon%lwc
      model%cdnc     = icon%cdnc
      model%reff     = icon%reff
      model%beta_ext = icon%beta_ext
      
   end subroutine models_setup_icon
   ! ----------------------------------------------------------------------------------------------------------------------------
   
   ! ----------------------------------------------------------------------------------------------------------------------------
   !> @brief Set up the `model` data structure by importing the NWPSAF data
   !> @warning The time for the NWPSAF data is set to 12:00:00
   !> @param[inout] model Model data structure containing all necessary variables for S3COM
   !> @param[in] nwpsaf structure containing the NWPSAF data
   subroutine models_setup_nwpsaf(model, nwpsaf)
      
      ! Inputs
      type(type_nwpsaf), intent(in) :: nwpsaf
      
      ! Outputs
      type(type_model), intent(inout) :: model
      
      ! Internal
      real(wp) :: hour, minute, sec
      ! character(len=8) :: nwpsaf_date
      
      ! Extract information from the NWPSAF time (%Y%m%d.%f)
      ! ----------------------------------------------------------------------------------------------------------
      hour = 12
      minute = 0
      sec = 0
      
      model%date(1) = nwpsaf%day(1)
      model%date(2) = nwpsaf%month(1)
      model%date(3) = nwpsaf%year(1)
      if(model%npoints > 1) write(*,*) "Warning: The date from the first pixel is currently used (to be improved)"
      ! ----------------------------------------------------------------------------------------------------------
      
      model%nlat       = nwpsaf%nlat
      model%nlon       = nwpsaf%nlon
      model%mode       = nwpsaf%mode
      
      model%npoints    = nwpsaf%npoints
      model%nlevels    = nwpsaf%nlevels
      model%nlayers    = nwpsaf%nlayers
      
      model%point      = nwpsaf%point
      
      model%lat_orig   = nwpsaf%lat_orig
      model%lon_orig   = nwpsaf%lon_orig
      
      model%height     = nwpsaf%height
      model%height_2   = nwpsaf%height_2
      
      model%lat        = nwpsaf%lat
      model%lon        = nwpsaf%lon
      model%topography = nwpsaf%elevation
      model%u_10m      = nwpsaf%u10
      model%v_10m      = nwpsaf%v10
      model%ts         = nwpsaf%tsurf
      model%ps         = nwpsaf%psurf
      model%q_2m       = nwpsaf%q2m
      model%t_2m       = nwpsaf%t2m
      model%landmask   = nwpsaf%lsm
      
      model%p          = nwpsaf%paph
      model%z          = nwpsaf%altitudeh
      model%dz         = nwpsaf%dz
      model%t          = nwpsaf%temph
      model%q          = nwpsaf%humh
      model%clc        = nwpsaf%cc
      model%iwc        = nwpsaf%iwc
      model%lwc        = nwpsaf%lwc
      model%reff       = nwpsaf%reff
      model%cdnc       = nwpsaf%cdnc
      
   end subroutine models_setup_nwpsaf
   ! ----------------------------------------------------------------------------------------------------------------------------
   
   ! ----------------------------------------------------------------------------------------------------------------------------
   !> @brief Compute the solar angles corresponding to model simualtions
   !> @details The solar zenith and azimuth angles are computed using the `solar_angles` function using
   !> the latitude, longitude, date and time of the model simulation.
   !> @todo RTTOV v13.2 has a function to compute the solar angles. This function should be used instead.
   subroutine models_setup_solar(model)
      
      ! Inputs/outputs
      type(type_model), intent(inout) :: model
      
      ! Internal
      real(wp) :: lon
      integer(kind = 4) :: i
      
      do i = 1, model%npoints
         
         ! This function requires longitude to be (-180, 180)
         lon = model%lon(i)
         
         if(lon > 180) lon = model%lon(i) - 360
         
         call solar_angles(model%lat(i), lon, model%date, model%time, model%sunzenangle(i), model%sunazangle(i))
         
      enddo
      
   end subroutine models_setup_solar
   ! ----------------------------------------------------------------------------------------------------------------------------
   
end module mod_models

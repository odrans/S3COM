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

module s3com_types

  implicit none

  private
  public :: wp, dp
  public :: type_s3com, type_nwpsaf, type_model, type_rttov_opt, type_nml, type_icon
  public :: type_cld, type_cld_mie

  !!Few kind definitions for variables
  integer, parameter :: sp = selected_real_kind(6, 37)     !< Kind for single precision reals
  integer, parameter :: dp = selected_real_kind(12, 307)   !< Kind for double precision reals
  integer, parameter :: wp = sp                            !< Kind for working precision reals
  integer, parameter :: wpi = selected_int_kind(4)         !< Kind for working precision integers

  !> @brief Namelist variables
  !> @details These variables are directly read from the namelist file that is provided as argument to the S3COM executable
  type type_nml
     character(len=256) ::      &
          path_s3com,           &    !< Path to S3COM directory
          path_rttov,           &    !< Path to RTTOV directory
          fname_in,             &    !< Name of the input model file (e.g. ICON or NWPSAF)
          path_out,             &    !< Path to the repository containing where the output files will be written
          suffix_out,           &    !< Suffix added to the output filenames
          model_name                 !< Name of the physical model. Currently supported: ICON, NWPSAF
     integer(kind=4) ::         &
          npoints_it,           &    !< Number of subset data points (chunks) to process in each iteration (only relevant to optimize memory usage)
          platform,             &    !< Platform ID for RTTOV
          satellite,            &    !< Satellite ID for RTTOV
          instrument,           &    !< Instrument ID for RTTOV
          nchannels,            &    !< Number of instrument channels to simulate
          ir_scatt_model,       &    !< Scattering model for IR source term: 1=DOM; 2=Chou-scaling
          vis_scatt_model,      &    !< Scattering model for solar source term: 1=DOM; 2=single-scattering; 3=MFASIS
          dom_nstreams,         &    !< Number of streams for DOM scattering model
          dom_nmoments,         &    !< Number of moments for discrete ordinate method
          rttov_nthreads,       &    !< Number of threads for RTTOV calculations
          gas_unit,             &    !< Gas units for atmospheric profiles
          ice_scheme,           &    !< Scheme used for ice cloud optical properties: 1=Baum; 2=Baran 2014; 3=Baran 2018. Not relevant if `user_cld_opt_param` is true.
          clw_scheme                 !< Scheme used for liquid cloud optical properties: 1=OPAC; 2=Deff. Not relevant if `user_cld_opt_param` is true.
     integer(kind=4), dimension(:), allocatable :: &
          channel_list               !< List of satellite channels to simulate (should be of dimension nchannels)
     logical :: &
          user_cld_opt_param,   &    !< If true, users can supply their own cloud optical properties (`ice_scheme` and `clw_scheme` are then not used)
          flag_retrievals,      &    !< Flag to specify if retrievals should be performed
          flag_output_atm,      &    !< Flag to specify if atmospheric profiles should be output
          flag_output_jac,      &    !< Flag to specify if Jacobian matrices should be output
          flag_output_k_tl,     &    !< Flag to specify if K matrices should be output
          do_jacobian_calc,     &    !< Flag to specify if Jacobian matrices should be calculated
          do_k_tl_calc,         &    !< Flag to specify if K matrices should be calculated
          do_opdep_calc,        &    !< If false, disables the RTTOV gas optical depth calculation (default = true)
          dom_rayleigh,         &    !< If true, Rayleigh scattering is included in the DOM model
          mmr_cldaer,           &    !< Cloud and gas units: true => kg/kg (cld+aer); false => g/m3 (cld), cm-3 (aer)
          ozone_data,           &    !< Set to true when supplying a profile of ozone, false to use climatology from RTTOV
          add_refrac,           &    !< If true accounts for atmospheric refraction
          add_clouds,           &    !< If true, clouds are included in the RTTOV model
          add_aerosols               !< If true, aerosols are included in the RTTOV model
  end type type_nml


  !> @brief Model outputs from NWPSAF simulations
  !! @details They are either read form the NetCDF file or calculated from these model output
  type type_nwpsaf
     integer(kind=4) ::        &
          npoints,             &     !< Total number of grid points
          nlevels,             &     !< Number of vertical levels
          nlayers,             &     !< Number of vertical layers (typically, nlevels - 1)
          nlat,                &     !< Number of latitude points in the grid
          nlon,                &     !< Number of longitude points in the grid
          mode                       !< Model grid type (1: track, 2: lon-lat, 3: lat-lon)
     integer(kind=4), dimension(:), allocatable :: &
          height,              &     !< Index of vertical layers
          height_2                   !< Index of vertical levels
     integer(kind=4), dimension(:), allocatable :: &
          point,               &     !< Index of grid points
          day,                 &     !< Day of the simulation
          month,               &     !< Month of the simulation
          year                       !< Year of the simulation
     real(wp), dimension(:), allocatable :: &
          lon,                 &     !< Longitude @units{degrees East}
          lat,                 &     !< Latitude @units{degrees North}
          lon_orig,            &     !< Longitude that won't be regridded @units{degrees East}
          lat_orig,            &     !< Latitude that won't be regridded @units{degrees North}
          elevation,           &     !< Surface height @units{m}
          lsm,                 &     !< Land/sea mask (0/1)
          psurf,               &     !< Surface pressure @units{Pa}
          tsurf,               &     !< Skin temperature @units{K}
          t2m,                 &     !<  2-m temperature @units{K}
          q2m,                 &     !<  2-m specific humidity @units{kg/kg}
          u10,                 &     !< U-component of 10-m wind @units{m/s}
          v10                        !< V-component of 10-m wind @units{m/s}
     real(wp), dimension(:,:), allocatable :: &
          p,                  &      !< Layer pressure @units{Pa}
          p_ifc,              &      !< Pressure at half-level center @units{Pa}
          altitude,           &      !< Layer height @units{m}
          altitudeh,          &      !< Height at half-levels center @units{m}
          t,                  &      !< Temperature @units{K}
          t_ifc,              &      !< Temperature at half-levels center @units{K}
          q,                  &      !< Specific humidity @units{kg/kg}
          q_ifc,              &      !< Specific humidity at half level center @units{kg/kg}
          clc,                &      !< Total cloud fraction (0-1)
          clw,                &      !< Specific cloud water content @units{kg/kg}
          cli,                &      !< Specific cloud ice content @units{kg/kg}
          qnc,                &      !< Cloud droplet number concentration @units{\# m^-3}
          qr,                 &      !< Rain mixing ratio @units{kg/kg}
          qs,                 &      !< Snow mixing ratio @units{kg/kg}
          dz,                 &      !< Layer thickness @units{m}
          rho,                &      !< Air density used for liquid clouds @units{kg m^-3}
          tv,                 &      !< Virtual temperature @units{K}
          lwc,                &      !< Liquid water content @units{kg m^-3}
          iwc,                &      !< Ice water content @units{kg m^-3}
          cdnc,               &      !< Cloud droplet number concentration @units{\# m^-3}
          Reff                       !< Cloud liquid water effective radius @units{m}
  end type type_nwpsaf

  !> @brief Model outputs from ICON simulations
  !! @details They are either read form the NetCDF file or calculated from these model output
  type type_icon
     integer(kind=4) :: &
          npoints,             &     !< Total number of grid points
          nlevels,             &     !< Number of vertical levels
          nlayers,             &     !< Number of vertical layers (typically, nlevels - 1)
          nlat,                &     !< Number of latitude points in the grid
          nlon,                &     !< Number of longitude points in the grid
          mode                       !< Model grid type (1: track, 2: lon-lat, 3: lat-lon)
     real(dp) ::               &
          time                       !< Time of the simulation (format is \%Y\%m\%d.\%f)
     integer(kind=4), dimension(:), allocatable :: &
          height,              &     !< Index of vertical layers
          height_2                   !< Index of vertical levels
     real(wp), dimension(:), allocatable :: &
          lon,                 &     !< Longitude @units{degrees East}
          lat,                 &     !< Latitude @units{degrees North}
          lon_orig,            &     !< Longitude that won't be regridded @units{degrees East}
          lat_orig,            &     !< Latitude that won't be regridded @units{degrees North}
          topography,          &     !< Surface height @units{m}
          landmask,            &     !< Land/sea mask (0/1)
          ps,                  &     !< Surface pressure @units{Pa}
          ts,                  &     !< Skin temperature @units{K}
          t_2m,                &     !< 2-m temperature @units{K}
          q_2m,                &     !< 2-m specific humidity @units{kg/kg}
          u_10m,               &     !< U-component of 10-m wind @units{m/s}
          v_10m                      !< V-component of 10-m wind @units{m/s}
     real(wp), dimension(:,:), allocatable :: &
          p,                  &      !< Layer pressure @units{Pa}
          p_ifc,              &      !< Pressure at half-level center @units{Pa}
          z,                  &      !< Layer height @units{m}
          z_ifc,              &      !< Height at half-levels center @units{m}
          t,                  &      !< Temperature @units{K}
          t_ifc,              &      !< Temperature at half-levels center @units{K}
          q,                  &      !< Specific humidity @units{kg/kg}
          q_ifc,              &      !< Specific humidity at half level center @units{kg/kg}
          clc,                &      !< Total cloud fraction (0-1)
          clw,                &      !< Specific cloud water content @units{kg/kg}
          cli,                &      !< Specific cloud ice content @units{kg/kg}
          qnc,                &      !< Cloud droplet number concentration @units{\# m^-3}
          qr,                 &      !< Rain mixing ratio @units{kg/kg}
          qs,                 &      !< Snow mixing ratio @units{kg/kg}
          dz,                 &      !< Layer thickness @units{m}
          rho,                &      !< Air density used for liquid clouds @units{kg m^-3}
          tv,                 &      !< Virtual temperature @units{K}
          lwc,                &      !< Liquid water content @units{kg m^-3}
          iwc,                &      !< Ice water content @units{kg m^-3}
          cdnc,               &      !< Cloud droplet number concentration @units{\# m^-3}
          Reff                       !< Cloud liquid water effective radius @units{m}
  end type type_icon

  !! Type containing variables stored for model outputs
  type type_model
     integer(kind=4) :: &
          nlevels, npoints, nlayers, nlat, nlon, mode, &!Dimensions
          idx_start,     & ! Starting index for subset profile
          idx_end          ! Ending index for subset profile
     integer(kind=4), dimension(:), allocatable :: &
          height, height_2
     integer(kind=4), dimension(3) :: &
          time, &   ! day, month, year
          date      ! hour, minute, second
     integer(kind=4), dimension(:), allocatable :: &
          point
     real(wp), dimension(:), allocatable :: &
          lon_orig,                           & !Longitude that won't be regridded (degrees east)
          lat_orig,                           & !Latitude  that won't be regridded (degress north)
          lon,                                & !Longitude (degrees east)
          lat,                                & !Latitude (degress north)
          topography,                         & !Surface height
          landmask,                           & !Land/sea mask (0/1)
          ps,                                 & !Surface pressure (Pa)
          ts,                                 & !Skin temperature (K)
          t_2m,                               & !2m temperature (K)
          q_2m,                               & !2m specific water vapor content (kg/kg)
          u_10m,                              & !U-component of wind (m/s)
          v_10m,                              & !V-component of wind (m/s)
          sunzenangle,                        & !Solar zenith angle
          sunazangle                            !Solar azimuth angle
     real(wp), dimension(:,:), allocatable :: &
          o3,                                 & !ozone
          co2,                                & !Carbon dioxide
          ch4,                                & !Methane
          n2o,                                & !n2o
          s2o,                                & !s2o
          co,                                 & !Carbon monoxide
          p,                                  & !Model pressure levels (pa)
          z,                                  & !Model level height (m)
          zh,                                 & !Model level height at half-levels (m)
          dz,                                 & !Layer thickness (m)
          t,                                  & !Temperature (K)
          q,                                  & !Specific humidity (kg/kg)
          clc,                                & !Total cloud fraction (0-1)
          qnc,                                & !Cloud droplet number concentration (particules/kg)
          qr,                                 & !Rain mixing ratio (kg/kg)
          qs,                                 & !Snow mixing ratio (kg/kg)
          lwc,                                & !Liquid water content (kg/m3)
          iwc,                                & !Ice water content (kg/m3)
          cdnc,                               & !Cloud droplet number concentration (1/m3)
          Reff
  end type type_model
  
  type type_rttov_opt
     integer ::       &
          dosolar,    &
          nchannels,  &
          nthreads,   &
          platform,   &
          satellite,  &
          instrument, &
          month,      &
          gas_units,   &
          ice_scheme, &
          clw_scheme, &
          ir_scatt_model, &
          vis_scatt_model, &
          dom_nstreams, &
          dom_nmoments
     logical :: &
          mmr_cldaer, &
          ozone_data, &
          co2_data  , &
          n2o_data  , &
          ch4_data  , &
          co_data   , &
          so2_data, &
          add_clouds, &
          add_aerosols, &
          add_refrac, &
          do_opdep_calc, &
          dom_rayleigh, &
          user_cld_opt_param
     integer, dimension(:), allocatable :: &
          channel_list
     character(len = 32) :: &
          platform_name, &
          inst_name, &
          sat_name
     real(wp) :: &
          zenangle, azangle
  end type type_rttov_opt

  type type_rttov_cld_opt_param
     integer :: &
          nmom
     integer, dimension(:), allocatable :: &
          mom
     real(wp), dimension(:), allocatable :: &
          phangle
     real(wp), dimension(:,:), allocatable :: &
          abs, &
          sca, &
          bpr
     real(wp), dimension(:,:,:), allocatable :: &
          legcoef, &
          pha
  end type type_rttov_cld_opt_param

  !!Type containing variables used by S3COM for retrievals
  type type_s3com_rad
     real(kind=wp), dimension(:), allocatable ::   &
          wavelength
     real(kind=wp), dimension(:,:), allocatable :: &
          y,                                       &
          f,                                       &
          f_ref_total,                             &
          f_ref_clear,                             &
          f_bt_total,                              &
          f_bt_clear,                              &
          f_rad_total,                             &
          f_rad_clear
  end type type_s3com_rad

  type type_s3com_jac
     real(kind=wp), dimension(:), allocatable ::   &
          wavelength
     real(kind=wp), dimension(:,:,:), allocatable :: &
          p,                                         &
          t,                                         &
          cfrac,                                     &
          clwde
     logical :: do_jacobian_calc
  end type type_s3com_jac

  type type_s3com_k_tl
     real(kind=wp), dimension(:), allocatable ::   &
          wavelength
     real(kind=wp), dimension(:,:,:), allocatable :: &
          t
     logical :: do_k_tl_calc
  end type type_s3com_k_tl
    
  type type_s3com_atm
     real(kind=wp), dimension(:,:), allocatable :: &
          t,                                       &
          z,                                       &
          clc,                                     &
          cdnc,                                    &
          reff,                                    &
          lwc
  end type type_s3com_atm

  type type_s3com_opt
     type(type_rttov_opt) :: rttov
  end type type_s3com_opt

  type type_s3com
     integer(kind=4), dimension(3) :: &
          time, &   ! day, month, year
          date      ! hour, minute, second
     integer(kind=4) :: &
          npoints,      &
          nlevels,      &
          nlayers,      &
          nmeas,        &
          nstates,      &
          idx_start,    &
          idx_end
     logical, dimension(:), allocatable :: &
          flag_rttov
     type(type_nml) :: nml
     type(type_s3com_rad) :: rad
     type(type_s3com_atm) :: atm
     type(type_s3com_jac) :: jac
     type(type_s3com_k_tl) :: k_tl
     type(type_s3com_opt) :: opt
  end type type_s3com

  type type_cld_mie
     integer(kind=4) :: &
          nang, &
          nchan, &
          nrad, &
          nmom
     character(len=128) :: &
          fn_mie,          &
          fn_legcoef
     integer(kind=4), dimension(:), allocatable :: &
          chan_id, mom
     real(kind=wp), dimension(:), allocatable :: &
          chan_wl, &
          radius, &
          angle
     real(kind=wp), dimension(:,:), allocatable :: &
          Cext, &
          Csca, &
          Cabs, &
          w0
     real(kind=wp), dimension(:,:,:), allocatable :: &
          pha, &
          legcoef
  end type type_cld_mie

  type type_cld
     type(type_cld_mie) :: mie
  end type type_cld

  !!Type containing variables used by S3COM for retrievals
  type type_s3com_obsolete
     integer(kind=4) :: &
          nstates,        &
          nmeas,          &
          npoints
     integer(kind=4), dimension(:), allocatable :: &
          ztop_ice_idx,                              &
          zbase_ice_idx,                             &
          ztop_liquid_idx,                           &
          zbase_liquid_idx,                          &
          n_iter, i_stepsize
     real(kind=wp), dimension(:), allocatable :: &
          g, gi, gip1, g_meas, stepsize
     logical, dimension(:), allocatable :: &
          flag_rttov, flag_testconv
     real(kind=wp), dimension(:,:), allocatable :: &
          Y,                                         &
          y_refl_total,                              &
          y_refl_clear,                              &
          y_bt_total,                                &
          y_bt_clear,                                &
          y_rad_total,                               &
          y_rad_clear,                               &
          y_rad_cloudy,                              &
          F,                                         &
          f_refl_total,                              &
          f_refl_clear,                              &
          f_bt_total,                                &
          f_bt_clear,                                &
          f_rad_total,                               &
          f_rad_clear,                               &
          f_rad_cloudy,                              &
          Xi,                                        &
          Xip1,                                      &
          Xa,                                        &
          brdf,                                      &
          emissivity, &
          t, &
          z, &
          clc, &
          reff, &
          cdnc
     real(kind=wp), dimension(:,:,:), allocatable :: &
          K,                                           &
          Kt,                                          &
          Sy,                                          &
          Sf,                                          &
          Se,                                          &
          Se_i,                                        &
          Sx,                                          &
          Sx_i,                                        &
          Sa,                                          &
          Sa_i
     real(kind=wp), dimension(:,:), allocatable :: &
          iwc,                                       & !Output measurement vector
          lwc,                                       &
          cla
     real(kind=wp), dimension(:), allocatable :: &
          ztop_ice,                                &
          zbase_ice,                               &
          iwp,                                     &
          iwp_model,                               &
          ztop_liquid,                             &
          zbase_liquid,                            &
          lwp,                                     &
          lwp_model,                               &
          cdnc_ret,                                &
          cdnc_model
  end type type_s3com_obsolete

contains


end module s3com_types

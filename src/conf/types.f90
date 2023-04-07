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
  integer, parameter :: wpi = selected_int_kind(8)         !< Kind for working precision integers

  !> @brief Structure containing namelist variables
  !> @details These variables are directly read from the namelist file that is provided as argument to the S3COM executable
  type type_nml
     character(len=256) ::      &
          path_s3com,           &    !< Path to S3COM directory
          path_rttov,           &    !< Path to RTTOV directory
          fname_in,             &    !< Name of the input model file (e.g. ICON or NWPSAF)
          path_out,             &    !< Path to the repository containing where the output files will be written
          suffix_out,           &    !< Suffix added to the output filenames
          model_name                 !< Name of the physical model. Currently supported: ICON, NWPSAF
     integer(wpi) ::            &
          npoints_it,           &    !< Number of subset data points (chunks) to process in each iteration (only relevant to optimize memory usage)
          platform,             &    !< Platform ID for RTTOV
          satellite,            &    !< Satellite ID for RTTOV
          instrument,           &    !< Instrument ID for RTTOV
          nchannels,            &    !< Number of instrument channels to simulate
          ir_scatt_model,       &    !< Scattering model for IR source term: 1=DOM; 2=Chou-scaling
          vis_scatt_model,      &    !< Scattering model for solar source term: 1=DOM; 2=single-scattering; 3=MFASIS
          dom_nstreams,         &    !< Number of streams for DOM scattering model
          dom_nmoments,         &    !< Number of moments for discrete ordinate method
          rttov_nthreads,       &    !< Number of threads for RTTOV calculations (if openMP is enabled)
          gas_unit,             &    !< Gas units for atmospheric profiles (1: kg/kg; 2: ppmv over moist air; 0: ppmv over dry air)
          ice_scheme,           &    !< Scheme used for ice cloud optical properties: 1=Baum; 2=Baran 2014; 3=Baran 2018. Not relevant if `user_cld_opt_param` is true.
          clw_scheme                 !< Scheme used for liquid cloud optical properties: 1=OPAC; 2=Deff. Not relevant if `user_cld_opt_param` is true.
     integer(wpi), dimension(:), allocatable :: &
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
          mmr_cldaer,           &    !< Cloud and aerosol units: true => kg/kg (cld+aer); false => g/m3 (cld), cm-3 (aer)
          ozone_data,           &    !< Set to true when supplying a profile of ozone, false to use climatology from RTTOV
          add_refrac,           &    !< If true accounts for atmospheric refraction
          add_clouds,           &    !< If true, clouds are included in the RTTOV model
          add_aerosols               !< If true, aerosols are included in the RTTOV model
  end type type_nml


  !> @brief Structure for model outputs from NWPSAF simulations
  !! @details They are either read form the NetCDF file (@input flag) or calculated from these model output
  !! @note This time is missing for NWP-SAF simulations, later set to 12:00:00 UTC for all data points
  type type_nwpsaf
     integer(wpi) ::        &
          npoints,             &     !< Total number of grid points
          nlevels,             &     !< Number of vertical levels
          nlayers,             &     !< Number of vertical layers (typically, nlevels - 1)
          nlat,                &     !< Number of latitude points in the grid
          nlon,                &     !< Number of longitude points in the grid
          mode                       !< Model grid type (1: track, 2: lon-lat, 3: lat-lon)
     integer(wpi), dimension(:), allocatable :: &
          height,              &     !< Index of vertical layers
          height_2                   !< Index of vertical levels
     integer(wpi), dimension(:), allocatable :: &
          point,               &     !< Index of grid points @input
          day,                 &     !< Day of the simulation @input
          month,               &     !< Month of the simulation @input
          year                       !< Year of the simulation @input
     real(wp), dimension(:), allocatable :: &
          lon,                 &     !< Longitude @units{degrees East}
          lat,                 &     !< Latitude @units{degrees North}
          lon_orig,            &     !< Longitude that won't be regridded @units{degrees East}
          lat_orig,            &     !< Latitude that won't be regridded @units{degrees North}
          elevation,           &     !< Surface height @units{m} @input
          lsm,                 &     !< Land/sea mask (0/1) @input
          psurf,               &     !< Surface pressure @units{Pa} @input
          tsurf,               &     !< Skin temperature @units{K} @input
          t2m,                 &     !<  2-m temperature @units{K} @input
          q2m,                 &     !<  2-m specific humidity @units{kg/kg}
          u10,                 &     !< Zonal 10-m wind @units{m/s} @input
          v10                        !< Meridional 10-m wind @units{m/s} @input
     real(wp), dimension(:,:), allocatable :: &
          pap,                 &      !< Layer pressure @units{Pa} @input
          paph,                &      !< Pressure at half-level center @units{Pa} @input
          altitude,            &      !< Layer height @units{m} @input
          altitudeh,           &      !< Height at half-levels center @units{m} @input
          temp,                &      !< Air temperature @units{K} @input
          temph,               &      !< Air temperature at half-levels center @units{K} @input
          hum,                 &      !< Specific humidity @units{kg/kg} @input
          humh,                &      !< Specific humidity at half level center @units{kg/kg} @input
          cc,                  &      !< Total cloud fraction (0-1) @input
          dz,                  &      !< Geometrical thickness of layer @units{m}
          rho,                 &      !< Air density used for liquid clouds @units{kg m^-3}
          tv,                  &      !< Virtual temperature @units{K}
          lwc,                 &      !< Liquid water content @units{kg m^-3} @input
          iwc,                 &      !< Ice water content @units{kg m^-3} @input
          cdnc,                &      !< Cloud droplet number concentration @units{\# m^-3}
          Reff                        !< Cloud liquid water effective radius @units{m}
  end type type_nwpsaf

  !> @brief Structure for model outputs from ICON simulations
  !! @details They are either read form the NetCDF file (@input flag) or calculated from these model output
  type type_icon
     integer(wpi) :: &
          npoints,             &     !< Total number of grid points
          nlevels,             &     !< Number of vertical levels
          nlayers,             &     !< Number of vertical layers (typically, nlevels - 1)
          nlat,                &     !< Number of latitude points in the grid
          nlon,                &     !< Number of longitude points in the grid
          mode                       !< Model grid type (1: track, 2: lon-lat, 3: lat-lon)
     real(dp) ::               &
          time                       !< Time of the simulation (format is \%Y\%m\%d.\%f UTC) @input
     integer(wpi), dimension(:), allocatable :: &
          height,              &     !< Index of vertical layers @input
          height_2                   !< Index of vertical levels @input
     real(wp), dimension(:), allocatable :: &
          lon,                 &     !< Longitude @units{degrees East} @input
          lat,                 &     !< Latitude @units{degrees North} @input
          lon_orig,            &     !< Longitude that won't be regridded @units{degrees East}
          lat_orig,            &     !< Latitude that won't be regridded @units{degrees North}
          topography,          &     !< Surface height @units{m} @input
          landmask,            &     !< Land/sea mask (0/1) @input
          ps,                  &     !< Surface pressure @units{Pa} @input
          ts,                  &     !< Skin temperature @units{K} @input
          t_2m,                &     !< 2-m temperature @units{K} @input
          q_2m,                &     !< 2-m specific humidity @units{kg/kg} @input
          u_10m,               &     !< Zonal 10-m wind @units{m/s} @input
          v_10m                      !< Meridional 10-m wind @units{m/s} @input
     real(wp), dimension(:,:), allocatable :: &
          p,                  &      !< Layer pressure @units{Pa} @input
          p_ifc,              &      !< Pressure at half-level center @units{Pa} @input
          z,                  &      !< Layer height @units{m} @input
          z_ifc,              &      !< Height at half-levels center @units{m} @input
          t,                  &      !< Temperature @units{K} @input
          t_ifc,              &      !< Temperature at half-levels center @units{K} @input
          q,                  &      !< Specific humidity @units{kg/kg} @input
          q_ifc,              &      !< Specific humidity at half level center @units{kg/kg} @input
          clc,                &      !< Total cloud fraction (0-1)  @input
          clw,                &      !< Specific cloud water content @units{kg/kg} @input
          cli,                &      !< Specific cloud ice content @units{kg/kg} @input
          qnc,                &      !< Cloud droplet number concentration @units{\# m^-3} @input
          qr,                 &      !< Rain mixing ratio @units{kg/kg} @input
          qs,                 &      !< Snow mixing ratio @units{kg/kg} @input
          dz,                 &      !< Layer thickness @units{m}
          rho,                &      !< Air density used for liquid clouds @units{kg m^-3}
          tv,                 &      !< Virtual temperature @units{K}
          lwc,                &      !< Liquid water content @units{kg m^-3}
          iwc,                &      !< Ice water content @units{kg m^-3}
          cdnc,               &      !< Cloud droplet number concentration @units{\# m^-3}
          Reff                       !< Cloud liquid water effective radius @units{m}
  end type type_icon

  !> @brief General structure for all atmospheric data from model outputs
  !! @details It is used to store model outputs consistently and completed with other variables.
  !! This is the structure that is used in the main program and passed for radiation calculations.
  type type_model
     integer(wpi) ::           &
          npoints,             &     !< Total number of grid points
          nlevels,             &     !< Number of vertical levels
          nlayers,             &     !< Number of vertical layers (typically, nlevels - 1)
          nlat,                &     !< Number of latitude points in the grid
          nlon,                &     !< Number of longitude points in the grid
          mode,                &     !< Model grid type (1: track, 2: lon-lat, 3: lat-lon)
          idx_start,           &     !< Starting point index for the subset grid
          idx_end                    !< Ending point index for the subset grid
     integer(wpi), dimension(:), allocatable :: &
          height,              &     !< Index of vertical layers
          height_2                   !< Index of vertical levels
     integer(wpi), dimension(3) :: &
          time,                &     !< Time of the day, UTC @units{/hour, minute, second/}
          date                       !< Day of the year @units{/day, month, year/}
     integer(wpi), dimension(:), allocatable :: &
          point
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
          u_10m,               &     !< Zonal 10-m wind @units{m/s}
          v_10m,               &     !< Meridional 10-m wind @units{m/s}
          sunzenangle,         &     !< Solar zenith angle @units{degrees}
          sunazangle                 !< Solar azimuth angle @units{degrees}
     real(wp), dimension(:,:), allocatable :: &
          o3,                  &     !< Ozone concentrations on model levels for user-defined gas profiles, see namelist configuration @units{gas_unit}
          co2,                 &     !< CO2 concentrations, similarly to o3
          ch4,                 &     !< CH4 concentrations, similarly to o3
          n2o,                 &     !< N2O concentrations, similarly to o3
          s2o,                 &     !< S2O concentrations, similarly to o3
          co,                  &     !< CO concentrations, similarly to o3
          p,                   &     !< Pressure on model levels @units{Pa}
          z,                   &     !< Altitude on model levels @units{m}
          dz,                  &     !< Geometrical thickness of model layer @units{m}
          t,                   &     !< Temperature on model levels @units{K}
          q,                   &     !< Specific humidity on model levels @units{kg/kg}
          clc,                 &     !< Total cloud fraction in model layer (0-1)
          lwc,                 &     !< Liquid water content in model layer @units{kg/m3}
          iwc,                 &     !< Ice water content in model layer @units{kg/m3}
          cdnc,                &     !< Cloud droplet number concentration in model layer @units{\# m^-3}
          reff                       !< Liquid water cloud effective radius in model layer @units{m}
  end type type_model

  !> @brief Structure containing main RTTOV options
  !! @details Data from namelist are indicated with @input. Others are set in `rttov_setup_opt`
  type type_rttov_opt
     integer(wpi) ::           &
          dosolar,             &     !< Activation of solar radiation calculations (0: false, 1: true)
          nchannels,           &     !< Number of instrument channels to simulate @input
          nthreads,            &     !< Number of threads for RTTOV calculations (if openMP is available) @input
          platform,            &     !< Platform ID for RTTOV @input
          satellite,           &     !< Satellite ID for RTTOV @input
          instrument,          &     !< Instrument ID for RTTOV @input
          month,               &     !< Month of the year, used to load surface brdf and emissivity data
          gas_units,           &     !< Gas units for atmospheric profiles (1: kg/kg; 2: ppmv over moist air; 0: ppmv over dry air)
          ice_scheme,          &     !< Scheme used for ice cloud optical properties: 1=Baum; 2=Baran 2014; 3=Baran 2018. Not relevant if `user_cld_opt_param` is true. @input
          clw_scheme,          &     !< Scheme used for liquid cloud optical properties: 1=OPAC; 2=Deff @input
          ir_scatt_model,      &     !< Scattering model for IR source term: 1=DOM; 2=Chou-scaling @input
          vis_scatt_model,     &     !< Scattering model for solar source term: 1=DOM; 2=single-scattering; 3=MFASIS @input
          dom_nstreams,        &     !< Number of streams for discrete ordinate method @input
          dom_nmoments               !< Number of moments for discrete ordinate method @input
     logical ::                &
          mmr_cldaer,          &     !< Cloud and aerosol units: true => kg/kg (cld+aer); false => g/m3 (cld), cm-3 (aer) @input
          ozone_data,          &     !< True when supplying a profile of ozone, false to use RTTOV climatologies @input
          co2_data,            &     !< True when supplying a profile of CO2, false to use RTTOV climatologies
          n2o_data,            &     !< True when supplying a profile of N2O, false to use RTTOV climatologies
          ch4_data,            &     !< True when supplying a profile of CH4, false to use RTTOV climatologies
          co_data,             &     !< True when supplying a profile of CO, false to use RTTOV climatologies
          so2_data,            &     !< True when supplying a profile of SO2, false to use RTTOV climatologies
          add_clouds,          &     !< True to add clouds to the RTTOV calculations @input
          add_aerosols,        &     !< True to add aerosols to the RTTOV calculations @input
          add_refrac,          &     !< True to add atmospheric refraction @input
          do_opdep_calc,       &     !< If false, disables the RTTOV gas optical depth calculation @input
          dom_rayleigh,        &     !< If false, disables the RTTOV Rayleigh scattering calculation @input
          user_cld_opt_param         !< If true, users can supply their own cloud optical properties (`ice_scheme` and `clw_scheme` are then not used) @input
     integer, dimension(:), allocatable :: &
          channel_list               !< List of channels to simulate @input
     character(len = 32) ::    &
          platform_name,       &     !< Platform name
          inst_name,           &     !< Instrument name
          sat_name                   !< Satellite name
     real(wp) :: &
          zenangle,            &     !< Satellite zenith angle @units{degrees}
          azangle                    !< Satellite azimuth angle @units{degrees}
  end type type_rttov_opt

  !> @brief Structure containing cloud optical properties from Mie calculations
  type type_cld_mie
     integer(wpi) ::           &
          nang,                &     !< Number of scattering angles for which Mie calculations were performed
          nchan,               &     !< Number of instrumental channels
          nrad,                &     !< Number of available effective radii on which Mie calculations were performed
          nmom                       !< Number of coefficients for Legendre polynomial decomposition
     character(len=128) ::     &
          fn_mie,              &     !< Name of the file containing the Mie optical properties for a given instrument
          fn_legcoef                 !< Name of the file containing the corresponding Legendre polynomial coefficients
     integer(wpi), dimension(:), allocatable :: &
          chan_id,             &     !< ID of the instrument channel in RTTOV
          mom                        !< Moments number of the Legendre polynomial decomposition
     real(kind=wp), dimension(:), allocatable :: &
          chan_wl,             &     !< Wavelength of the instrument channel @units{um}
          radius,              &     !< Effective radius of the spherical cloud particles @units{um}
          angle                      !< Scattering angle @units{degrees}
     real(kind=wp), dimension(:,:), allocatable :: &
          Cext,                &     !< Extinction cross-section coefficient @units{um^2}
          Csca,                &     !< Scattering cross-section coefficient @units{um^2}
          Cabs,                &     !< Absorption cross-section coefficient @units{um^2}
          w0                         !< Single-scattering albedo
     real(kind=wp), dimension(:,:,:), allocatable :: &
          pha,                 &     !< Azimuthally-averaged phase function
          legcoef                    !< Coefficients of the Legendre polynomial decomposition of the phase function
  end type type_cld_mie

  !> @brief Structure containing all cloud optical properties
  type type_cld
     type(type_cld_mie) :: mie
  end type type_cld



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


  ! Type containing variables used by S3COM for retrievals
  ! type type_s3com_obsolete
  !    integer(kind=4) :: &
  !         nstates,        &
  !         nmeas,          &
  !         npoints
  !    integer(kind=4), dimension(:), allocatable :: &
  !         ztop_ice_idx,                              &
  !         zbase_ice_idx,                             &
  !         ztop_liquid_idx,                           &
  !         zbase_liquid_idx,                          &
  !         n_iter, i_stepsize
  !    real(kind=wp), dimension(:), allocatable :: &
  !         g, gi, gip1, g_meas, stepsize
  !    logical, dimension(:), allocatable :: &
  !         flag_rttov, flag_testconv
  !    real(kind=wp), dimension(:,:), allocatable :: &
  !         Y,                                         &
  !         y_refl_total,                              &
  !         y_refl_clear,                              &
  !         y_bt_total,                                &
  !         y_bt_clear,                                &
  !         y_rad_total,                               &
  !         y_rad_clear,                               &
  !         y_rad_cloudy,                              &
  !         F,                                         &
  !         f_refl_total,                              &
  !         f_refl_clear,                              &
  !         f_bt_total,                                &
  !         f_bt_clear,                                &
  !         f_rad_total,                               &
  !         f_rad_clear,                               &
  !         f_rad_cloudy,                              &
  !         Xi,                                        &
  !         Xip1,                                      &
  !         Xa,                                        &
  !         brdf,                                      &
  !         emissivity, &
  !         t, &
  !         z, &
  !         clc, &
  !         reff, &
  !         cdnc
  !    real(kind=wp), dimension(:,:,:), allocatable :: &
  !         K,                                           &
  !         Kt,                                          &
  !         Sy,                                          &
  !         Sf,                                          &
  !         Se,                                          &
  !         Se_i,                                        &
  !         Sx,                                          &
  !         Sx_i,                                        &
  !         Sa,                                          &
  !         Sa_i
  !    real(kind=wp), dimension(:,:), allocatable :: &
  !         iwc,                                       & !Output measurement vector
  !         lwc,                                       &
  !         cla
  !    real(kind=wp), dimension(:), allocatable :: &
  !         ztop_ice,                                &
  !         zbase_ice,                               &
  !         iwp,                                     &
  !         iwp_model,                               &
  !         ztop_liquid,                             &
  !         zbase_liquid,                            &
  !         lwp,                                     &
  !         lwp_model,                               &
  !         cdnc_ret,                                &
  !         cdnc_model
  ! end type type_s3com_obsolete

contains


end module s3com_types

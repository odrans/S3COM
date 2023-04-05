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

  !!Few kind definitions for variables
  integer, parameter :: sp = selected_real_kind(6, 37)
  integer, parameter :: dp = selected_real_kind(12, 307)
  integer, parameter :: wp = sp

  ! Namelist
  type type_nml
     character(len=256) :: &
          path_rttov, &
          fname_in, &
          path_out, &
          suffix_out, &
          model_name
     integer(kind=4) :: &
          npoints_it, &
          month, &
          platform, &
          satellite, &
          instrument, &
          nchannels, &
          ir_scatt_model, &
          vis_scatt_model, &
          dom_nstreams, &
          rttov_nthreads, &
          gas_unit, &
          ice_scheme, &
          clw_scheme
     integer(kind=4), dimension(:), allocatable :: &
          channel_list
     logical :: &
          flag_retrievals,  &
          flag_output_atm,  &
          flag_output_jac,  &
          flag_output_k_tl, &
          do_jacobian_calc, &
          do_k_tl_calc,     &
          do_opdep_calc,    &
          dom_rayleigh,     &
          mmr_cldaer,       &
          ozone_data,       &
          add_refrac,       &
          add_clouds,       &
          add_aerosols,     &
          switchrad
  end type type_nml

  !!Type containing variables from NWP-SAF simulations
  type type_nwpsaf
     integer(kind=4) :: &
          npoints, &
          nlevels, &
          nlayers, &
          nlat, &
          nlon, &
          mode
     real(dp) :: &
          time
     integer(kind=4), dimension(:), allocatable :: &
          height, height_2 ! height index
     integer(kind=4), dimension(:), allocatable :: &
          point,                            &
          day,                              &
          month,                            &
          year
     real(wp), dimension(:), allocatable :: &
          lon,                              & !Longitude (degrees east)
          lat,                              & !Latitude (degress north)
          lon_orig,                         & !Longitude that won't be regridded (degrees east)
          lat_orig,                         & !Latitude  that won't be regridded (degress north)
          elevation,                        & !Surface height
          lsm,                              & !Land/sea mask (0/1)
          psurf,                            & !Surface pressure (Pa)
          tsurf,                            & !Skin temperature (K)
          t2m,                              & ! 2-m temperature (K)
          q2m,                              & ! 2-m specific humidity
          u10,                              & !U-component of wind (m/s)
          v10                                 !V-component of wind (m/s)
     real(wp), dimension(:,:), allocatable :: &
          p,                                    & !Layer pressure (Pa)
          p_ifc,                                & !Pressure at half-level center (Pa)
          z,                                    & !Layer height (m)
          z_ifc,                                & !Height at half-levels center (m)
          t,                                    & !Temperature (K)
          t_ifc,                                & !Temperature at half-levels center
          q,                                    & !Specific humidity (kg/kg)
          q_ifc,                                & !Specific humidity at half level center (kg/kg)
          clc,                                  & !Total cloud fraction (0-1)
          clw,                                  & !Specific cloud water content (kg/kg)
          cli,                                  & !Specific cloud ice content (kg/kg)
          qnc,                                  & !Cloud droplet number concentration (particules/kg)
          qr,                                   & !Rain mixing ratio (kg/kg)
          qs,                                   & !Snow mixing ratio (kg/kg)
          dz,                                   & !Layer thickness (m)
          rho,                                  & !Air density used for liquid clouds (kg/m3)
          tv,                                   & !Virtual temperature (K)
          lwc,                                  & !Liquid water content (kg/m3)
          iwc,                                  & !Ice water content (kg/m3)
          cdnc,                                 & !Cloud droplet number concentration (1/m3)
          Reff                                    !Cloud liquid water effective radius (m)
  end type type_nwpsaf

  !!Type containing variables from ICON simulations
  type type_icon
     integer(kind=4) :: &
          nlevels, &
          npoints, &
          nlayers, &
          nlat, &
          nlon, &
          mode
     real(dp) :: &
          time
     integer(kind=4), dimension(:), allocatable :: &
          height,                                  &
          height_2,                                & ! height index
          ztop_liq_idx
     real(wp), dimension(:), allocatable :: &
          lon_orig,                         & !Longitude that won't be regridded (degrees east)
          lat_orig,                         & !Latitude  that won't be regridded (degress north)
          lon,                              & !Longitude (degrees east)
          lat,                              & !Latitude (degress north)
          topography,                       & !Surface height
          landmask,                         & !Land/sea mask (0/1)
          ps,                               & !Surface pressure (Pa)
          ts,                               & !Skin temperature (K)
          t_2m,                             & !2m temperature (K)
          q_2m,                             & !2m specific water vapor content (kg/kg)
          u_10m,                            & !U-component of wind (m/s)
          v_10m,                            & !V-component of wind (m/s)
          lwp,                              & !Liquid water path (g/m2)
          iwp,                              & !Ice water path (g/m2)
          cod,                              & !Cloud optical depth (-)
          reff_top                            !Cloud droplet effective radius at cloud top (µm)
     real(wp), dimension(:,:), allocatable :: &
          p,                                    & !Layer pressure (Pa)
          p_ifc,                                & !Pressure at half-level center (Pa)
          z,                                    & !Layer height (m)
          z_ifc,                                & !Height at half-levels center (m)
          t,                                    & !Temperature (K)
          t_ifc,                                & !Temperature at half-levels center
          q,                                    & !Specific humidity (kg/kg)
          q_ifc,                                & !Specific humidity at half level center (kg/kg)
          clc,                                  & !Total cloud fraction (0-1)
          clw,                                  & !Specific cloud water content (kg/kg)
          cli,                                  & !Specific cloud ice content (kg/kg)
          qnc,                                  & !Cloud droplet number concentration (kg-1)
          qr,                                   & !Rain mixing ratio (kg/kg)
          qs,                                   & !Snow mixing ratio (kg/kg)
          dz,                                   & !Layer thickness (m)
          rho,                                  & !Air density used for liquid clouds (kg/m3)
          tv,                                   & !Virtual temperature (K)
          es_w,                                 & !Saturation vapour pressure of liquid water (Pa)
          es_i,                                 & !Saturation vapour pressure of ice water (Pa)
          lwc,                                  & !Liquid water content (kg/m3)
          iwc,                                  & !Ice water content (kg/m3)
          cdnc,                                 & !Cloud droplet number concentration (cm-3)
          reff,                                 & !Cloud droplet effective radius (µm)
          beta_ext                                !Cloud droplet extinction coefficient (m-1)
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
          lwp,                                & !Liquid water path (g/m2)
          iwp,                                & !Ice water path (g/m2)
          cod,                                & !Cloud optical depth (-)
          reff_top,                           & !Cloud droplet effective radius at cloud top (um)
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
          qnc,                                & !Cloud droplet number concentration (kg-1)
          qr,                                 & !Rain mixing ratio (kg/kg)
          qs,                                 & !Snow mixing ratio (kg/kg)
          lwc,                                & !Liquid water content (kg/m3)
          iwc,                                & !Ice water content (kg/m3)
          cdnc,                               & !Cloud droplet number concentration (cm-3)
          reff,                               & !Cloud droplet effective radius (µm)
          beta_ext                              !Cloud droplet extinction coefficient (m-1)
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
          dom_nstreams
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
          switchrad,  &
          do_opdep_calc, &
          dom_rayleigh
     integer, dimension(:), allocatable :: &
          channel_list
     character(len = 32) :: &
          platform_name, &
          inst_name
     real(wp) :: &
          zenangle, azangle
  end type type_rttov_opt

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
     real(kind=wp), dimension(:,:), allocatable ::   &
          lwp
     real(kind=wp), dimension(:,:,:), allocatable :: &
          p,                                         &
          t,                                         &
          cfrac,                                     &
          clwde,                                     &
          cdnc
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
     real(kind=wp), dimension(:), allocatable :: &
        cod
     real(kind=wp), dimension(:,:), allocatable :: &
          z,                                       &
          dz,                                      &
          lwc,                                     &
          cdnc,                                    &
          reff,                                    &
          beta_ext
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



end module s3com_types

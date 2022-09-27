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

MODULE s3com_types

  IMPLICIT NONE

  !!Few kind definitions for variables
  INTEGER, PARAMETER :: sp = selected_real_kind(6, 37)
  INTEGER, PARAMETER :: dp = selected_real_kind(12, 307)
  INTEGER, PARAMETER :: wp = sp

  ! Namelist
  TYPE type_nml
     CHARACTER(LEN = 256) :: &
          path_rttov, &
          fname_in, &
          path_out, &
          suffix_out
     INTEGER(KIND=4) :: &
          npoints_it, &
          month, &
          platform, &
          satellite, &
          instrument, &
          nchannels, &
          ir_scatt_model, &
          vis_scatt_model, &
          dom_nstreams, &
          rttov_nthreads
     INTEGER(KIND = 4), DIMENSION(:), ALLOCATABLE :: &
          channel_list
     LOGICAL :: &
          flag_retrievals, &
          addrefrac, &
          dom_rayleigh, &
          flag_output_atm
  END TYPE type_nml

  !!Type containing variables used by S3COM for retrievals
  TYPE type_s3com
     INTEGER(KIND=4) :: &
          nstates,        &
          nmeas,          &
          npoints
     INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: &
          ztop_ice_idx,                              &
          zbase_ice_idx,                             &
          ztop_liquid_idx,                           &
          zbase_liquid_idx,                          &
          n_iter, i_stepsize
     REAL(KIND=wp), DIMENSION(:), ALLOCATABLE :: &
          g, gi, gip1, g_meas, stepsize
     LOGICAL, DIMENSION(:), ALLOCATABLE :: &
          flag_rttov, flag_testconv
     REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE :: &
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

     REAL(KIND=wp), DIMENSION(:,:,:), ALLOCATABLE :: &
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
     REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE :: &
          iwc,                                       & !Output measurement vector
          lwc,                                       &
          cla
     REAL(KIND=wp), DIMENSION(:), ALLOCATABLE :: &
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
  END TYPE type_s3com

  !!Type containing variables from ICON simulations
  TYPE type_icon
     INTEGER(KIND=4) :: &
          Nlevels, Npoints, Nlat, Nlon, mode !Dimensions
     INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: &
          height
     REAL(wp) ::  &
          co2,      & !Carbon dioxide
          ch4,      & !Methane
          n2o,      & !n2o
          co          !Carbon monoxide
     REAL(wp), DIMENSION(:), ALLOCATABLE :: &
          lon_orig,                           & !Longitude that won't be regridded (degrees east)
          lat_orig,                           & !Latitude  that won't be regridded (degress north)
          lon,                                & !Longitude (degrees east)
          lat,                                & !Latitude (degress north)
          orography,                          & !Surface height
          landmask,                           & !Land/sea mask (0/1)
          psfc,                               & !Surface pressure (Pa)
          skt,                                & !Skin temperature (K)
          t2m,                                & !2m temperature (K)
          q2m,                                & !2m specific water vapor content (kg/kg)
          u_wind,                             & !U-component of wind (m/s)
          v_wind                                !V-component of wind (m/s)
     REAL(wp), DIMENSION(:,:), ALLOCATABLE :: &
          p,                                    & !Model pressure levels (pa)
          z,                                    & !Model level height (m)
          zh,                                   & !Model level height at half-levels (m)
          dz,                                   & !Layer thickness (m)
          t,                                    & !Temperature (K)
          sh,                                   & !Specific humidity (kg/kg)
          tca,                                  & !Total cloud fraction (0-1)
          clw,                                  & !Specific cloud water content (kg/kg)
          cli,                                  & !Specific cloud ice content (kg/kg)
          qnc,                                  & !Cloud droplet number concentration (particules/kg)
          lwc,                                  & !Liquid water content (kg/m3)
          iwc,                                  & !Ice water content (kg/m3)
          t_celcius,                            & !Temperature (degrees Celcius)
          e_sat,                                & !Saturation vapour pressure of water (Pa)
          ssh,                                  & !Saturation specific humidity (kg/kg)
          wv,                                   & !Humidity ratio (0-1)
          rh,                                   & !Relative humidity (%)
          pv,                                   & !Pressure of water vapour (Pa)
          pd,                                   & !Pressure of dry air (Pa)
          rho,                                  & !Air density used for liquid clouds (kg/m3)
          rho_atm,                              & !Air density used for ice clouds (kg/m3)
          cdnc,                                 & !Cloud droplet number concentration (1/m3)
          Reff,                                 & !Cloud liquid water effective radius (m)
          Deff,                                 & !Cloud liquid water effective diameter (m)
          beta_ext,                             & !Cloud extinction coefficient (1/m)
          tv,                                   & !Virtual temperature (K)
          dz_cod,                               & !Layer thickness (m)
          cod                                     !Cloud optical depth (unitless)
  END TYPE type_icon

  !!Type containing variables from RTTOV simulations
  TYPE type_rttov_atm
     INTEGER, POINTER :: &
          idx_start,       & !Index of starting ICON point
          idx_end,         & !Index of ending ICON point
          nPoints,         & !Number of profiles to simulate
          nLevels            !Number of levels
     REAL(wp), POINTER :: &
          co2,              & !Carbon dioxide
          ch4,              & !Methane
          n2o,              & !n2o
          co                  !Carbon monoxide
     REAL(wp), DIMENSION(:), POINTER :: &
          h_surf,                         & !Surface height
          u_surf,                         & !U component of surface wind
          v_surf,                         & !V component of surface wind
          t_skin,                         & !Surface skin temperature
          p_surf,                         & !Surface pressure
          t2m,                            & !2 m Temperature
          q2m,                            & !2 m Specific humidity
          lsmask,                         & !land-sea mask
          lat,                            & !Latitude
          lon                               !Longitude
     REAL(wp), DIMENSION(:,:), POINTER :: &
          z,                              & !Height at model levels
          dz,                             & !Layer depth at model levels
          p,                              & !Pressure at model levels
          t,                              & !Temperature
          q,                              & !Specific humidity
          o3,                             & !Ozone
          reff,                           &
          cdnc

     !!These fields below are needed ONLY for the RTTOV all-sky brightness temperature
     REAL(wp), DIMENSION(:,:), POINTER :: &
          tca,                              & !Cloud fraction
          iwc,                              & !Cloud ice water content
          lwc                                 !Cloud liquid water content
  END TYPE type_rttov_atm

  TYPE type_rttov_opt
     INTEGER ::       &
          dosolar,    &
          nchannels,  &
          nthreads,   &
          manthreads, &
          platform,   &
          satellite,  &
          instrument, &
          month         !Month (needed for surface emissivity calculation)
     INTEGER, DIMENSION(:), ALLOCATABLE :: &
          channel_list
     REAL(wp) :: &
          zenangle, azangle, sunzenangle, sunazangle
  END TYPE type_rttov_opt

END MODULE s3com_types

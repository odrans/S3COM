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

MODULE MOD_RTTOV

  USE s3com_types,        ONLY: wp, type_rttov_atm, type_rttov_opt, type_s3com
  USE mod_rttov_utils, ONLY: idx_rttov

  !!rttov_const contains useful RTTOV constants
  USE rttov_const, ONLY:    &
       errorstatus_success,   &
       errorstatus_fatal,     &
       platform_name,         &
       inst_name,             &
       surftype_sea,          &
       surftype_land,         &
       watertype_fresh_water, &
       watertype_ocean_water, &
       sensor_id_mw,          &
       sensor_id_po,          &
       wcl_id_stco,           &
       wcl_id_stma

  USE rttov_types, ONLY: &
       rttov_options,      &
       rttov_coefs,        &
       rttov_profile,      &
       rttov_transmission, &
       rttov_radiance,     &
       rttov_chanprof,     &
       rttov_emissivity,   &
       rttov_reflectance,  &
       rttov_opt_param

  !!The rttov_emis_atlas_data type must be imported separately
  USE mod_rttov_emis_atlas, ONLY: &
       rttov_emis_atlas_data,       &
       atlas_type_ir, atlas_type_mw

  !!The rttov_brdf_atlas_data type must be imported separately
  USE mod_rttov_brdf_atlas, ONLY: rttov_brdf_atlas_data

  USE rttov_unix_env, ONLY: rttov_exit

  USE parkind1, ONLY: jpim, jprb, jplm

  IMPLICIT NONE

#include "rttov_direct.interface"
#include "rttov_parallel_direct.interface"
#include "rttov_read_coefs.interface"
#include "rttov_dealloc_coefs.interface"
#include "rttov_alloc_direct.interface"
#include "rttov_user_options_checkinput.interface"
#include "rttov_print_opts.interface"
#include "rttov_print_profile.interface"
#include "rttov_skipcommentline.interface"

  !!Use emissivity atlas
#include "rttov_setup_emis_atlas.interface"
#include "rttov_get_emis.interface"
#include "rttov_deallocate_emis_atlas.interface"

  !!Use BRDF atlas
#include "rttov_setup_brdf_atlas.interface"
#include "rttov_get_brdf.interface"
#include "rttov_deallocate_brdf_atlas.interface"

  INTEGER(KIND=jpim), PARAMETER :: iup   = 20 !Unit for input profile file
  INTEGER(KIND=jpim), PARAMETER :: ioout = 21 !Unit for output

  !!==========================================================================================================================!!
  !! RTTOV variables/structures                                                                                               !!
  !!==========================================================================================================================!!

  TYPE(rttov_options)              :: opts                     ! Options structure
  TYPE(rttov_coefs)                :: coefs                    ! Coefficients structure
  TYPE(rttov_chanprof),    POINTER :: chanprof(:)    => NULL() ! Input channel/profile list
  LOGICAL(KIND=jplm),      POINTER :: calcemis(:)    => NULL() ! Flag to indicate calculation of emissivity within RTTOV
  TYPE(rttov_emissivity),  POINTER :: emissivity(:)  => NULL() ! Input/output surface emissivity
  LOGICAL(KIND=jplm),      POINTER :: calcrefl(:)    => NULL() ! Flag to indicate calculation of BRDF within RTTOV
  TYPE(rttov_reflectance), POINTER :: reflectance(:) => NULL() ! Input/output surface BRDF
  TYPE(rttov_profile),     POINTER :: profiles(:)    => NULL() ! Input profiles
  TYPE(rttov_transmission)         :: transmission             ! Output transmittances
  TYPE(rttov_radiance)             :: radiance                 ! Output radiances
  TYPE(rttov_opt_param)            :: cld_opt_param            ! Input cloud optical parameters
  TYPE(rttov_emis_atlas_data)      :: emis_atlas               ! Data structure for emissivity atlas
  TYPE(rttov_brdf_atlas_data)      :: brdf_atlas               ! Data structure for BRDF atlas

  INTEGER(KIND=jpim)               :: errorstatus              ! Return error status of RTTOV subroutine calls

  INTEGER(KIND=jpim) :: atlas_type
  INTEGER(KIND=jpim) :: alloc_status
  CHARACTER(LEN=11)  :: NameOfRoutine = 'example_fwd'

  !!==========================================================================================================================!!
  !! Variables for input                                                                                                      !!
  !!==========================================================================================================================!!

  CHARACTER(LEN=256) :: coef_filename
  CHARACTER(LEN=256) :: prof_filename, cld_coef_filename
  INTEGER(KIND=jpim) :: nthreads
  INTEGER(KIND=jpim) :: dosolar
  INTEGER(KIND=jpim) :: nlevels
  INTEGER(KIND=jpim) :: nprof
  INTEGER(KIND=jpim) :: nchannels
  INTEGER(KIND=jpim) :: nchanprof
  INTEGER(KIND=jpim), ALLOCATABLE :: channel_list(:)
  REAL(KIND=jprb)    :: trans_out(10)

  ! Loop variables
  INTEGER(KIND=jpim) :: j, jch
  INTEGER(KIND=jpim) :: np, nch
  INTEGER(KIND=jpim) :: ilev, nprint
  INTEGER(KIND=jpim) :: iprof, joff, ichan
  INTEGER            :: ios, imonth

  ! Initialization parameters
  INTEGER ::   &
       platform, & !RTTOV platform
       sensor,   & !RTTOV instrument
       satellite

CONTAINS

  !!==========================================================================================================================!!
  !! SUBROUTINE rttov_column                                                                                                  !!
  !!==========================================================================================================================!!

  SUBROUTINE run_rttov(rttov_atm,rttov_opt,oe,dealloc)

    USE s3com_config, ONLY: rd

    !!Inputs variables
    TYPE(type_rttov_atm), INTENT(IN) :: rttov_atm
    TYPE(type_rttov_opt), INTENT(IN) :: rttov_opt

    LOGICAL, INTENT(IN) :: dealloc !Flag to determine whether to deallocate RTTOV types

    !!Inout/Outputs variables
    TYPE(type_s3com), INTENT(INOUT) :: oe

    !!Local variables
    INTEGER, DIMENSION(:), ALLOCATABLE :: list_points
    INTEGER                            :: errorstatus, alloc_status, idx_prof, ilevel

    errorstatus = 0_jpim

    nthreads = 36 !35

    list_points = idx_rttov(oe)
    nprof = size(list_points); nlevels = rttov_atm%nlevels

    oe%f_refl_total(list_points,:) = 0._wp
    oe%f_refl_clear(list_points,:) = 0._wp
    oe%f_bt_total(list_points,:)   = 0._wp
    oe%f_bt_clear(list_points,:)   = 0._wp
    oe%f_rad_total(list_points,:)  = 0._wp
    oe%f_rad_clear(list_points,:)  = 0._wp
    oe%f_rad_cloudy(list_points,:) = 0._wp
    oe%brdf(list_points,:)         = 0._wp
    oe%emissivity(list_points,:)   = 0._wp

    IF (nprof .EQ. 0) RETURN

    !!--------------------------------------------------------------------------------------------------------------------!!
    !! 3. Allocate RTTOV input and output structures
    !!--------------------------------------------------------------------------------------------------------------------!!

    !Determine the total number of radiances to simulate (nchanprof).
    !In this example we simulate all specified channels for each profile, but
    !in general one can simulate a different number of channels for each profile.

    nchanprof = nchannels * nprof

    !!Allocate structures for rttov_direct

    CALL rttov_alloc_direct(      &
         errorstatus,             &
         1_jpim,                  & !1 => allocate
         nprof,                   &
         nchanprof,               &
         nlevels,                 &
         chanprof,                &
         opts,                    &
         profiles,                &
         coefs,                   &
         transmission,            &
         radiance,                &
         calcemis=calcemis,       &
         emissivity=emissivity,   &
         calcrefl=calcrefl,       &
         reflectance=reflectance, &
         init=.TRUE._jplm)

    IF (errorstatus /= errorstatus_success) THEN
       WRITE(*,*) 'allocation error for rttov_direct structures'
       CALL rttov_exit(errorstatus)
    ENDIF

    !!--------------------------------------------------------------------------------------------------------------------!!
    !! 4. Build the list of profile/channel indices in chanprof                                                           !!
    !!--------------------------------------------------------------------------------------------------------------------!!

    nch = 0_jpim

    DO j = 1, nprof
       DO jch = 1, nchannels
          nch = nch + 1_jpim
          chanprof(nch)%prof = j
          chanprof(nch)%chan = channel_list(jch)
       ENDDO
    ENDDO

    !!--------------------------------------------------------------------------------------------------------------------!!
    !! 5. Read profile data                                                                                               !!
    !!--------------------------------------------------------------------------------------------------------------------!!

    !!Gas units for profiles
    profiles(:)%gas_units = 1 ! Units: kg/kg

    !!Loop over all profiles and read data for each one
    DO iprof = 1, nprof

       idx_prof = list_points(iprof)

       !!Pressure (hPa), temperature (K) and water vapour (kg/kg)
       profiles(iprof)%p(:) = rttov_atm%p(idx_prof,:)*1E-2
       profiles(iprof)%t(:) = rttov_atm%t(idx_prof,:)
       profiles(iprof)%q(:) = rttov_atm%q(idx_prof,:)

       !!Air variables in 2 meters
       profiles(iprof)%s2m%t = rttov_atm%t2m(idx_prof)         !2m temperature (K)
       profiles(iprof)%s2m%q = rttov_atm%q2m(idx_prof)         !2m water vapour (kg/kg)
       profiles(iprof)%s2m%p = rttov_atm%p_surf(idx_prof)*1E-2 !Surface pressure (hPa)
       profiles(iprof)%s2m%u = rttov_atm%u_surf(idx_prof)      !10m zonal wind (m/s)
       profiles(iprof)%s2m%v = rttov_atm%v_surf(idx_prof)      !10m meridional wind (m/s)
       profiles(iprof)%s2m%wfetc = 100000                      !Used typical value given in documentation

       !!Skin variables
       profiles(iprof)%skin%t = rttov_atm%t_skin(idx_prof)
       profiles(iprof)%skin%salinity = 0.0                        !tmp, use other typical value
       profiles(iprof)%skin%fastem = (/3.0, 5.0, 15.0, 0.1, 0.3/) !Typical RTTOV default for land

       !!Surface type and water type
       IF (rttov_atm%lsmask(iprof) < 0.5) THEN
          profiles(iprof)%skin%surftype = surftype_sea
       ELSE
          profiles(iprof)%skin%surftype = surftype_land
       ENDIF

       profiles(iprof)%skin%watertype = watertype_fresh_water !tmp, adapt this to truth,
       !fresh more likely for ICON-DE simulations

       !!Elevation (km), latitude (deg) and longitude (deg)
       profiles(iprof)%elevation = rttov_atm%h_surf(idx_prof)*1E-3
       profiles(iprof)%latitude  = rttov_atm%lat(idx_prof)
       profiles(iprof)%longitude = rttov_atm%lon(idx_prof)

       !!Satellite and solar angles (deg)
       profiles(iprof)%zenangle    = rttov_opt%zenangle
       profiles(iprof)%azangle     = rttov_opt%azangle
       profiles(iprof)%sunzenangle = rttov_opt%sunzenangle
       profiles(iprof)%sunazangle  = rttov_opt%sunazangle

       profiles(iprof)%mmr_cldaer = .FALSE. !Logical flag to set cloud and aerosol
       !Units: true => kg/kg (cld+aer); false => g/m3 (cld), cm-3 (aer)
       !(must be the same for all profiles)

       !!Cloud variables for simple cloud scheme, set cfraction to 0. to turn this off (VIS/IR only)
       profiles(iprof)%cfrac = rttov_atm%tca(idx_prof,:)

       !! Used by OPAC
       profiles(iprof)%cloud(1,:) = rttov_atm%lwc(idx_prof,:)*1E3 !(kg/m3)

       !!Ice cloud input profiles
       profiles(iprof)%ice_scheme = 3 !Cloud ice water scheme: 1=Baum; 2=Baran 2014; 3=Baran 2018
       profiles(iprof)%cloud(6,:) = rttov_atm%iwc(idx_prof,:)*1E3 !(kg/m3)

       !!Liquid cloud input profiles
       profiles(:)%clw_scheme = 2 !Cloud liquid water scheme: 1=OPAC; 2=“Deff”
       profiles(iprof)%clwde(:) = rttov_atm%reff(idx_prof,:)*2.0 ! Need the diameter

    ENDDO

    !!--------------------------------------------------------------------------------------------------------------------!!
    !! 6. Specify surface emissivity and reflectance                                                                      !!
    !!--------------------------------------------------------------------------------------------------------------------!!

    !In this example we have no values for input emissivity
    !emissivity(:)%emis_in = 0._jprb

    !!Use emissivity atlas
    CALL rttov_get_emis( &
         errorstatus,    &
         opts,           &
         chanprof,       &
         profiles,       &
         coefs,          &
         emis_atlas,     &
         emissivity(:)%emis_in)

    IF (errorstatus /= errorstatus_success) THEN
       WRITE(*,*) 'error reading emissivity atlas'
       CALL rttov_exit(errorstatus)
    ENDIF

    !!Calculate emissivity within RTTOV where the atlas emissivity value is zero or less
    calcemis(:) = (emissivity(:)%emis_in <= 0._jprb)

    IF (opts%rt_ir%addsolar) THEN
       !In this example we have no values for input reflectances
       !reflectance(:)%refl_in = 0._jprb

       !!Use BRDF atlas
       CALL rttov_get_brdf(         &
            errorstatus,            &
            opts,                   &
            chanprof,               &
            profiles,               &
            coefs,                  &
            brdf_atlas,             &
            reflectance(:)%refl_in)

       IF (errorstatus /= errorstatus_success) THEN
          WRITE(*,*) 'error reading BRDF atlas'
          CALL rttov_exit(errorstatus)
       ENDIF

       !!Calculate BRDF within RTTOV where the atlas BRDF value is zero or less
       calcrefl(:) = (reflectance(:)%refl_in <= 0._jprb)
    ENDIF

    !!Use the RTTOV emissivity and BRDF calculations over sea surfaces
    DO j = 1, SIZE(chanprof)
       IF (profiles(chanprof(j)%prof)%skin%surftype == surftype_sea) THEN
          calcemis(j) = .TRUE.
          calcrefl(j) = .TRUE.
       ENDIF
    ENDDO

    !!Use default cloud top BRDF for simple cloud in VIS/NIR channels
    reflectance(:)%refl_cloud_top = 0._jprb

    !!--------------------------------------------------------------------------------------------------------------------!!
    !! 7. Call RTTOV forward model                                                                                        !!
    !!--------------------------------------------------------------------------------------------------------------------!!

    IF (nthreads <= 1) THEN
       CALL rttov_direct(              &
            errorstatus,               & !out   error flag
            chanprof,                  & !in    channel and profile index structure
            opts,                      & !in    options structure
            profiles,                  & !in    profile array
            coefs,                     & !in    coefficients structure
            transmission,              & !inout computed transmittances
            radiance,                  & !inout computed radiances
            calcemis    = calcemis,    & !in    flag for internal emissivity calcs
            emissivity  = emissivity,  & !inout input/output emissivities per channel
            calcrefl    = calcrefl,    & !in    flag for internal BRDF calcs
            reflectance = reflectance)   !inout input/output BRDFs per channel
    ELSE
       CALL rttov_parallel_direct(     &
            errorstatus,               & !out   error flag
            chanprof,                  & !in    channel and profile index structure
            opts,                      & !in    options structure
            profiles,                  & !in    profile array
            coefs,                     & !in    coefficients structure
            transmission,              & !inout computed transmittances
            radiance,                  & !inout computed radiances
            calcemis    = calcemis,    & !in    flag for internal emissivity calcs
            emissivity  = emissivity,  & !inout input/output emissivities per channel
            calcrefl    = calcrefl,    & !in    flag for internal BRDF calcs
            reflectance = reflectance, & !inout input/output BRDFs per channel
            nthreads    = nthreads)      !in    number of threads to use
    ENDIF

    IF (errorstatus /= errorstatus_success) THEN
       WRITE (*,*) 'rttov_direct error'
       CALL rttov_exit(errorstatus)
    ENDIF

    !!Output the results
    DO iprof = 1, nprof

       idx_prof = list_points(iprof)
       joff = (iprof-1_jpim) * nchannels
       ichan = 1

       DO j = 1+joff, nchannels+joff

          oe%f_refl_total(idx_prof,ichan) = radiance%refl(j)
          oe%f_refl_clear(idx_prof,ichan) = radiance%refl_clear(j)
          oe%f_bt_total(idx_prof,ichan)   = radiance%bt(j)
          oe%f_bt_clear(idx_prof,ichan)   = radiance%bt_clear(j)
          oe%brdf(idx_prof,ichan)         = reflectance(j)%refl_out
          oe%emissivity(idx_prof,ichan)   = emissivity(j)%emis_out
          oe%f_rad_total(idx_prof,ichan)  = radiance%total(j)*coefs%coef%ff_cwn(chanprof(j)%chan)**2*1E-7 !(W/m2/sr/um)
          oe%f_rad_clear(idx_prof,ichan)  = radiance%clear(j)*coefs%coef%ff_cwn(chanprof(j)%chan)**2*1E-7 !(W/m2/sr/um)
          oe%f_rad_cloudy(idx_prof,ichan) = radiance%cloudy(j)

          ichan = ichan + 1

       ENDDO

    ENDDO

    !!--------------------------------------------------------------------------------------------------------------------!!
    !! 8. Deallocate all RTTOV arrays and structures                                                                      !!
    !!--------------------------------------------------------------------------------------------------------------------!!

    !!Deallocate structures for rttov_direct
    CALL rttov_alloc_direct(    &
         errorstatus,           &
         0_jpim,                & !0 => deallocate
         nprof,                 &
         nchanprof,             &
         nlevels,               &
         chanprof,              &
         opts,                  &
         profiles,              &
         coefs,                 &
         transmission,          &
         radiance,              &
         calcemis=calcemis,     &
         emissivity=emissivity, &
         calcrefl=calcrefl,     &
         reflectance=reflectance)

    IF (errorstatus /= errorstatus_success) THEN
       WRITE(*,*) 'deallocation error for rttov_direct structures'
       CALL rttov_exit(errorstatus)
    ENDIF

    IF (dealloc) THEN
       CALL rttov_dealloc_coefs(errorstatus, coefs)

       IF (errorstatus /= errorstatus_success) THEN
          WRITE(*,*) 'coefs deallocation error'
       ENDIF

       CALL rttov_deallocate_emis_atlas(emis_atlas) !Deallocate emissivity atlas

       IF (opts%rt_ir%addsolar) THEN
          CALL rttov_deallocate_brdf_atlas(brdf_atlas) !Deallocate BRDF atlas
       ENDIF
    ENDIF

  END SUBROUTINE run_rttov

END MODULE MOD_RTTOV

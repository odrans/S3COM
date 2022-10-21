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

module MOD_RTTOV

  use s3com_types,        only: wp, type_model, type_rttov_opt, type_s3com, type_s3com_new
  use mod_rttov_utils, only: idx_rttov

  !!rttov_const contains useful RTTOV constants
  use rttov_const, only:    &
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

  use rttov_types, only: &
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
  use mod_rttov_emis_atlas, only: &
       rttov_emis_atlas_data,       &
       atlas_type_ir, atlas_type_mw

  !!The rttov_brdf_atlas_data type must be imported separately
  use mod_rttov_brdf_atlas, only: rttov_brdf_atlas_data

  use rttov_unix_env, only: rttov_exit

  use parkind1, only: jpim, jprb, jplm

  implicit none

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

  integer(KIND=jpim), parameter :: iup   = 20 !Unit for input profile file
  integer(KIND=jpim), parameter :: ioout = 21 !Unit for output

  !!==========================================================================================================================!!
  !! RTTOV variables/structures                                                                                               !!
  !!==========================================================================================================================!!

  type(rttov_options)              :: opts                     ! Options structure
  type(rttov_coefs)                :: coefs                    ! Coefficients structure
  type(rttov_chanprof),    pointer :: chanprof(:)    => null() ! Input channel/profile list
  logical(KIND=jplm),      pointer :: calcemis(:)    => null() ! Flag to indicate calculation of emissivity within RTTOV
  type(rttov_emissivity),  pointer :: emissivity(:)  => null() ! Input/output surface emissivity
  logical(KIND=jplm),      pointer :: calcrefl(:)    => null() ! Flag to indicate calculation of BRDF within RTTOV
  type(rttov_reflectance), pointer :: reflectance(:) => null() ! Input/output surface BRDF
  type(rttov_profile),     pointer :: profiles(:)    => null() ! Input profiles
  type(rttov_transmission)         :: transmission             ! Output transmittances
  type(rttov_radiance)             :: radiance                 ! Output radiances
  type(rttov_opt_param)            :: cld_opt_param            ! Input cloud optical parameters
  type(rttov_emis_atlas_data)      :: emis_atlas               ! Data structure for emissivity atlas
  type(rttov_brdf_atlas_data)      :: brdf_atlas               ! Data structure for BRDF atlas

  integer(KIND=jpim)               :: errorstatus              ! Return error status of RTTOV subroutine calls

  integer(KIND=jpim) :: atlas_type

  character(LEN=11)  :: NameOfRoutine = 'example_fwd'

  !!==========================================================================================================================!!
  !! Variables for input                                                                                                      !!
  !!==========================================================================================================================!!

  character(LEN=256) :: coef_filename
  character(LEN=256) :: prof_filename, cld_coef_filename
  integer(KIND=jpim) :: nthreads
  integer(KIND=jpim) :: dosolar
  integer(KIND=jpim) :: nlevels
  integer(KIND=jpim) :: nprof
  integer(KIND=jpim) :: nchannels
  integer(KIND=jpim) :: nchanprof
  integer(KIND=jpim), allocatable :: channel_list(:)
  real(KIND=jprb)    :: trans_out(10)

  ! Loop variables
  integer(KIND=jpim) :: j, jch
  integer(KIND=jpim) :: np, nch
  integer(KIND=jpim) :: ilev, nprint
  integer(KIND=jpim) :: iprof, joff, ichan
  integer            :: ios, imonth

  ! Initialization parameters
  integer ::   &
       platform, & !RTTOV platform
       sensor,   & !RTTOV instrument
       satellite

contains

  !!==========================================================================================================================!!
  !! SUBROUTINE rttov_column                                                                                                  !!
  !!==========================================================================================================================!!

  subroutine run_rttov(rttov_atm,rttov_opt,oe,dealloc)

    use s3com_config, only: rd

    !!Inputs variables
    type(type_model), intent(IN) :: rttov_atm
    type(type_rttov_opt), intent(IN) :: rttov_opt

    logical, intent(IN) :: dealloc !Flag to determine whether to deallocate RTTOV types

    !!Inout/Outputs variables
    type(type_s3com_new), intent(INOUT) :: oe

    !!Local variables
    integer, dimension(:), allocatable :: list_points
    integer                            :: errorstatus, idx_prof

    errorstatus = 0_jpim

    nthreads = rttov_opt%nthreads

    list_points = idx_rttov(oe)
    nprof = size(list_points); nlevels = rttov_atm%nlevels

    if (nprof .eq. 0) return

    !!--------------------------------------------------------------------------------------------------------------------!!
    !! 3. Allocate RTTOV input and output structures
    !!--------------------------------------------------------------------------------------------------------------------!!

    !Determine the total number of radiances to simulate (nchanprof).
    !In this example we simulate all specified channels for each profile, but
    !in general one can simulate a different number of channels for each profile.

    nchanprof = nchannels * nprof

    !!Allocate structures for rttov_direct

    call rttov_alloc_direct(      &
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
         init=.true._jplm)

    if (errorstatus /= errorstatus_success) then
       write(*,*) 'allocation error for rttov_direct structures'
       call rttov_exit(errorstatus)
    endif

    !!--------------------------------------------------------------------------------------------------------------------!!
    !! 4. Build the list of profile/channel indices in chanprof                                                           !!
    !!--------------------------------------------------------------------------------------------------------------------!!

    nch = 0_jpim

    do j = 1, nprof
       do jch = 1, nchannels
          nch = nch + 1_jpim
          chanprof(nch)%prof = j
          chanprof(nch)%chan = channel_list(jch)
       enddo
    enddo

    !!--------------------------------------------------------------------------------------------------------------------!!
    !! 5. Read profile data                                                                                               !!
    !!--------------------------------------------------------------------------------------------------------------------!!

    !!Gas units for profiles
    profiles(:)%gas_units = 1 ! Units: kg/kg

    !!Loop over all profiles and read data for each one
    do iprof = 1, nprof

       idx_prof = list_points(iprof)

       !!Pressure (hPa), temperature (K) and water vapour (kg/kg)
       profiles(iprof)%p(:) = rttov_atm%p(idx_prof,:)*1E-2
       profiles(iprof)%t(:) = rttov_atm%t(idx_prof,:)
       profiles(iprof)%q(:) = rttov_atm%q(idx_prof,:)


       !!Air variables in 2 meters
       profiles(iprof)%s2m%t = rttov_atm%t_2m(idx_prof)        !2m temperature (K)
       profiles(iprof)%s2m%q = rttov_atm%q_2m(idx_prof)        !2m water vapour (kg/kg)
       profiles(iprof)%s2m%p = rttov_atm%ps(idx_prof)*1E-2     !Surface pressure (hPa)
       profiles(iprof)%s2m%u = rttov_atm%u_10m(idx_prof)       !10m zonal wind (m/s)
       profiles(iprof)%s2m%v = rttov_atm%v_10m(idx_prof)       !10m meridional wind (m/s)
       profiles(iprof)%s2m%wfetc = 100000                      !Used typical value given in documentation

       !!Skin variables
       profiles(iprof)%skin%t = rttov_atm%ts(idx_prof)
       profiles(iprof)%skin%salinity = 0.0                        !tmp, use other typical value
       profiles(iprof)%skin%fastem = (/3.0, 5.0, 15.0, 0.1, 0.3/) !Typical RTTOV default for land

       !!Surface type and water type
       if (rttov_atm%landmask(iprof) < 0.5) then
          profiles(iprof)%skin%surftype = surftype_sea
       else
          profiles(iprof)%skin%surftype = surftype_land
       endif

       profiles(iprof)%skin%watertype = watertype_fresh_water !tmp, adapt this to truth,
       !fresh more likely for ICON-DE simulations

       !!Elevation (km), latitude (deg) and longitude (deg)
       profiles(iprof)%elevation = rttov_atm%topography(idx_prof)*1E-3
       profiles(iprof)%latitude  = rttov_atm%lat(idx_prof)
       profiles(iprof)%longitude = rttov_atm%lon(idx_prof)

       !!Satellite and solar angles (deg)
       profiles(iprof)%zenangle    = rttov_opt%zenangle
       profiles(iprof)%azangle     = rttov_opt%azangle
       profiles(iprof)%sunzenangle = rttov_atm%sunzenangle(idx_prof)
       profiles(iprof)%sunazangle  = rttov_atm%sunazangle(idx_prof)

       profiles(iprof)%mmr_cldaer = .false. !Logical flag to set cloud and aerosol
       !Units: true => kg/kg (cld+aer); false => g/m3 (cld), cm-3 (aer)

       !!Cloud variables for simple cloud scheme, set cfraction to 0. to turn this off (VIS/IR only)
       profiles(iprof)%cfrac = rttov_atm%clc(idx_prof,:)

       !! Used by OPAC
       profiles(iprof)%cloud(1,:) = rttov_atm%lwc(idx_prof,:)*1E3 !(kg/m3)

       !!Ice cloud input profiles
       profiles(iprof)%ice_scheme = 3 !Cloud ice water scheme: 1=Baum; 2=Baran 2014; 3=Baran 2018
       profiles(iprof)%cloud(6,:) = rttov_atm%iwc(idx_prof,:)*1E3 !(kg/m3)

       !!Liquid cloud input profiles
       profiles(:)%clw_scheme = 2 !Cloud liquid water scheme: 1=OPAC; 2=“Deff”
       profiles(iprof)%clwde(:) = rttov_atm%reff(idx_prof,:)*2.0 ! Need the diameter

    enddo

    !!--------------------------------------------------------------------------------------------------------------------!!
    !! 6. Specify surface emissivity and reflectance                                                                      !!
    !!--------------------------------------------------------------------------------------------------------------------!!

    !!Use emissivity atlas
    call rttov_get_emis( &
         errorstatus,    &
         opts,           &
         chanprof,       &
         profiles,       &
         coefs,          &
         emis_atlas,     &
         emissivity(:)%emis_in)

    if (errorstatus /= errorstatus_success) then
       write(*,*) 'error reading emissivity atlas'
       call rttov_exit(errorstatus)
    endif

    !!Calculate emissivity within RTTOV where the atlas emissivity value is zero or less
    calcemis(:) = (emissivity(:)%emis_in <= 0._jprb)

    if (opts%rt_ir%addsolar) then

       !!Use BRDF atlas
       call rttov_get_brdf(         &
            errorstatus,            &
            opts,                   &
            chanprof,               &
            profiles,               &
            coefs,                  &
            brdf_atlas,             &
            reflectance(:)%refl_in)

       if (errorstatus /= errorstatus_success) then
          write(*,*) 'error reading BRDF atlas'
          call rttov_exit(errorstatus)
       endif

       !!Calculate BRDF within RTTOV where the atlas BRDF value is zero or less
       calcrefl(:) = (reflectance(:)%refl_in <= 0._jprb)
    endif

    !!Use the RTTOV emissivity and BRDF calculations over sea surfaces
    do j = 1, size(chanprof)
       if (profiles(chanprof(j)%prof)%skin%surftype == surftype_sea) then
          calcemis(j) = .true.
          calcrefl(j) = .true.
       endif
    enddo

    !!Use default cloud top BRDF for simple cloud in VIS/NIR channels
    reflectance(:)%refl_cloud_top = 0._jprb

    !!--------------------------------------------------------------------------------------------------------------------!!
    !! 7. Call RTTOV forward model                                                                                        !!
    !!--------------------------------------------------------------------------------------------------------------------!!

    if (nthreads <= 1) then
       call rttov_direct(              &
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
    else
       call rttov_parallel_direct(     &
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
    endif

    if (errorstatus /= errorstatus_success) then
       write (*,*) 'rttov_direct error'
       call rttov_exit(errorstatus)
    endif

    !!Output the results
    do iprof = 1, nprof

       idx_prof = list_points(iprof)
       joff = (iprof-1_jpim) * nchannels
       ichan = 1

       do j = 1+joff, nchannels+joff

          oe%rad%f_ref_total(idx_prof,ichan)   = radiance%refl(j)
          oe%rad%f_ref_clear(idx_prof,ichan)   = radiance%refl_clear(j)
          oe%rad%f_bt_total(idx_prof,ichan)    = radiance%bt(j)
          oe%rad%f_bt_clear(idx_prof,ichan)    = radiance%bt_clear(j)
          oe%rad%f_rad_total(idx_prof,ichan)   = radiance%total(j)*coefs%coef%ff_cwn(chanprof(j)%chan)**2*1E-7 !(W/m2/sr/um)
          oe%rad%f_rad_clear(idx_prof,ichan)   = radiance%clear(j)*coefs%coef%ff_cwn(chanprof(j)%chan)**2*1E-7 !(W/m2/sr/um)

          ! oe%brdf(idx_prof,ichan)         = reflectance(j)%refl_out
          ! oe%emissivity(idx_prof,ichan)   = emissivity(j)%emis_out

          ichan = ichan + 1

       enddo

    enddo

    !!--------------------------------------------------------------------------------------------------------------------!!
    !! 8. Deallocate all RTTOV arrays and structures                                                                      !!
    !!--------------------------------------------------------------------------------------------------------------------!!

    !!Deallocate structures for rttov_direct
    call rttov_alloc_direct(    &
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

    if (errorstatus /= errorstatus_success) then
       write(*,*) 'deallocation error for rttov_direct structures'
       call rttov_exit(errorstatus)
    endif

    if (dealloc) then
       call rttov_dealloc_coefs(errorstatus, coefs)

       if (errorstatus /= errorstatus_success) then
          write(*,*) 'coefs deallocation error'
       endif

       call rttov_deallocate_emis_atlas(emis_atlas) !Deallocate emissivity atlas

       if (opts%rt_ir%addsolar) then
          call rttov_deallocate_brdf_atlas(brdf_atlas) !Deallocate BRDF atlas
       endif
    endif

  end subroutine run_rttov

end module MOD_RTTOV

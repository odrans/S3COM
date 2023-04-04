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

module mod_rttov

  use s3com_types,        only: wp, type_model, type_rttov_opt, type_s3com, type_cld
  use mod_rttov_utils,    only: find_idx_rttov, check_rttov_status
    
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

  use parkind1, only: jpim, jprb, jplm
  
  use rttov_unix_env, only: rttov_exit

  implicit none

  private
  public :: run_rttov, opts, coefs, emis_atlas, brdf_atlas

#include "rttov_direct.interface"
#include "rttov_tl.interface"
#include "rttov_k.interface"

#include "rttov_parallel_direct.interface"
#include "rttov_parallel_tl.interface"
#include "rttov_parallel_k.interface"

#include "rttov_alloc_direct.interface"
#include "rttov_alloc_tl.interface"
#include "rttov_alloc_k.interface"

#include "rttov_init_opt_param.interface"
#include "rttov_read_coefs.interface"
#include "rttov_dealloc_coefs.interface"
#include "rttov_init_emis_refl.interface"
#include "rttov_init_prof.interface"
#include "rttov_init_rad.interface"
#include "rttov_init_transmission.interface"
#include "rttov_user_options_checkinput.interface"
#include "rttov_print_opts.interface"
#include "rttov_print_profile.interface"
#include "rttov_skipcommentline.interface"

! Use emissivity atlas
#include "rttov_setup_emis_atlas.interface"
#include "rttov_get_emis.interface"
#include "rttov_deallocate_emis_atlas.interface"

! Use BRDF atlas
#include "rttov_setup_brdf_atlas.interface"
#include "rttov_get_brdf.interface"
#include "rttov_deallocate_brdf_atlas.interface"

  type(rttov_options)              :: opts                       ! Options structure
  type(rttov_coefs)                :: coefs                      ! Coefficients structure
  type(rttov_chanprof),    pointer :: chanprof(:)      => null() ! Input channel/profile list
  logical(kind=jplm),      pointer :: calcemis(:)      => null() ! Flag to indicate calculation of emissivity within RTTOV
  type(rttov_emissivity),  pointer :: emissivity(:)    => null() ! Input/output surface emissivity
  type(rttov_emissivity),  pointer :: emissivity_k(:)  => null() ! Emissivity Jacobians
  type(rttov_emissivity),  pointer :: emissivity_tl(:)  => null() ! Input/output surface emissivity perturbations
  logical(kind=jplm),      pointer :: calcrefl(:)      => null() ! Flag to indicate calculation of BRDF within RTTOV
  type(rttov_reflectance), pointer :: reflectance(:)   => null() ! Input/output surface BRDF
  type(rttov_reflectance), pointer :: reflectance_k(:) => null() ! Reflectance Jacobians
  type(rttov_reflectance), pointer :: reflectance_tl(:) => null() ! Input/output surface BRDF perturbations
  type(rttov_profile),     pointer :: profiles(:)      => null() ! Input profiles
  type(rttov_profile),     pointer :: profiles_k(:)    => null() ! Output Jacobians
  type(rttov_profile),     pointer :: profiles_tl(:)    => null() ! Input atmospheric profile and surface variable perturbations
  type(rttov_transmission)         :: transmission               ! Output transmittances
  type(rttov_transmission)         :: transmission_k             ! Transmittance Jacobians
  type(rttov_transmission)         :: transmission_tl             ! Output transmittance perturbations
  type(rttov_radiance)             :: radiance                   ! Output radiances
  type(rttov_radiance)             :: radiance_k                 ! Radiance Jacobians
  type(rttov_radiance)             :: radiance_tl                 ! Output radiance, BT and BRF perturbations
  type(rttov_opt_param)            :: cld_opt_param              ! Input cloud optical parameters
  type(rttov_emis_atlas_data)      :: emis_atlas                 ! Data structure for emissivity atlas
  type(rttov_brdf_atlas_data)      :: brdf_atlas                 ! Data structure for BRDF atlas


  integer(kind=jpim) :: nthreads
  integer(kind=jpim) :: nlevels, nlayers
  integer(kind=jpim) :: nprof
  integer(kind=jpim) :: nchannels
  integer(kind=jpim) :: nchanprof
  integer(kind=jpim) :: max_mom, nleg

  real(kind=jprb), parameter :: tl_perturbation = -0.01_jprb

  ! Loop variables
  integer(kind=jpim) :: j, jch
  integer(kind=jpim) :: nch
  integer(kind=jpim) :: ilev, ilay
  integer(kind=jpim) :: iprof, joff, ichan, ichanprof, n_true
  integer(kind=jpim) :: idx_reff

contains

  subroutine run_rttov(rttov_atm, rttov_opt, s3com, cld, dealloc)

    !!Inputs variables
    type(type_model), intent(in) :: rttov_atm
    type(type_rttov_opt), intent(in) :: rttov_opt
    type(type_cld), intent(in) :: cld

    logical, intent(IN) :: dealloc !Flag to determine whether to deallocate RTTOV types

    !!Inout/Outputs variables
    type(type_s3com), intent(inout) :: s3com

    !!Local variables
    integer, dimension(:), allocatable :: list_points
    integer                            :: errorstatus, idx_prof

    errorstatus = 0_jpim

    !> Find how any profiles are to be processed (using flag_rttov)
    n_true = count(s3com%flag_rttov); allocate(list_points(n_true))
    list_points = find_idx_rttov(s3com)

    !> Set up a few useful dimensions
    nthreads = rttov_opt%nthreads
    max_mom = cld%mie%nmom
    nleg = max_mom + 1
    nchannels = rttov_opt%nchannels
    nlevels = rttov_atm%nlevels
    nlayers = rttov_atm%nlayers
    nprof = size(list_points)
    nchanprof = nchannels * nprof

    if (nprof == 0) return

    !> Allocate RTTOV input and output structures depending on the model that will be called
    !> Current options are direct, jacobian (K) and tangent linear (TL)
    !> Note that K and TL also include the direct model calls
    !> ----------------------------------------------------------------------------------------------------
    if(s3com%jac%do_jacobian_calc) then     ! Jacobien
       call rttov_alloc_k(                   &
            errorstatus,                     &
            1_jpim,                          & !1 => allocate
            nprof,                           &
            nchanprof,                       &
            nlevels,                         &
            chanprof,                        &
            opts,                            &
            profiles,                        &
            profiles_k,                      &
            coefs,                           &
            transmission,                    &
            transmission_k,                  &
            radiance,                        &
            radiance_k,                      &
            calcemis      = calcemis,        &
            emissivity    = emissivity,      &
            emissivity_k  = emissivity_k,    &
            calcrefl      = calcrefl,        &
            reflectance   = reflectance,     &
            reflectance_k = reflectance_k,   &
            init          = .true._jplm)
    else if(s3com%k_tl%do_k_tl_calc) then    ! Tangent linear
       call rttov_alloc_tl(                  &
            errorstatus,                     &
            1_jpim,                          & !1 => allocate
            nprof,                           &
            nchanprof,                       &
            nlevels,                         &
            chanprof,                        &
            opts,                            &
            profiles,                        &
            profiles_tl,                     &
            coefs,                           &
            transmission,                    &
            transmission_tl,                 &
            radiance,                        &
            radiance_tl,                     &
            calcemis       = calcemis,       &
            emissivity     = emissivity,     &
            emissivity_tl  = emissivity_tl,  &
            calcrefl       = calcrefl,       &
            reflectance    = reflectance,    &
            reflectance_tl = reflectance_tl, &
            init           = .true._jplm)
    else ! Direct model
       call rttov_alloc_direct(               &
            errorstatus,                      &
            1_jpim,                           & !1 => allocate
            nprof,                            &
            nchanprof,                        &
            nlevels,                          &
            chanprof,                         &
            opts,                             &
            profiles,                         &
            coefs,                            &
            transmission,                     &
            radiance,                         &
            calcemis      = calcemis,         &
            emissivity    = emissivity,       &
            calcrefl      = calcrefl,         &
            reflectance   = reflectance,      &
            cld_maxnmom = max_mom, &
            cld_nphangle = cld%mie%nang,      &
            cld_opt_param = cld_opt_param,    &
            init          = .true._jplm)
    endif
    call check_rttov_status(errorstatus, "rttov_alloc")
    ! ----------------------------------------------------------------------------------------------------


    ! Build the list of profile/channel indices in chanprof
    ! ----------------------------------------------------------------------------------------------------
    nch = 0_jpim

    do j = 1, nprof
       do jch = 1, nchannels
          nch = nch + 1_jpim
          chanprof(nch)%prof = j
          chanprof(nch)%chan = rttov_opt%channel_list(jch)
       enddo
    enddo
    ! ----------------------------------------------------------------------------------------------------


    !!--------------------------------------------------------------------------------------------------------------------!!
    !! 5. Read profile data                                                                                               !!
    !!--------------------------------------------------------------------------------------------------------------------!!

    !!Gas units for profiles
    profiles(1:nprof)%gas_units = rttov_opt%gas_units

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

       profiles(iprof)%skin%watertype = watertype_fresh_water !tmp, adapt this

       !!Elevation (km), latitude (deg) and longitude (deg)
       profiles(iprof)%elevation = rttov_atm%topography(idx_prof)*1E-3
       profiles(iprof)%latitude  = rttov_atm%lat(idx_prof)
       profiles(iprof)%longitude = rttov_atm%lon(idx_prof)

       !!Satellite and solar angles (deg)
       profiles(iprof)%zenangle    = rttov_opt%zenangle
       profiles(iprof)%azangle     = rttov_opt%azangle
       profiles(iprof)%sunzenangle = rttov_atm%sunzenangle(idx_prof)
       profiles(iprof)%sunazangle  = rttov_atm%sunazangle(idx_prof)

       profiles(iprof)%mmr_cldaer = rttov_opt%mmr_cldaer !Logical flag to set cloud and aerosol
       !Units: true => kg/kg (cld+aer); false => g/m3 (cld), cm-3 (aer)

       ! Cloud variables for simple cloud scheme, set cfraction to 0. to turn this off (VIS/IR only)
       profiles(iprof)%cfrac(:) = rttov_atm%clc(idx_prof,:)

       ! call rttov_print_profile (profiles(iprof))
       ! stop
    end do

    if(opts%rt_ir%user_cld_opt_param) then
       cld_opt_param%nmom = cld%mie%nmom
       cld_opt_param%phangle(:) = cld%mie%angle(:)
       ichanprof = 0_jpim
       do iprof = 1, nprof
          idx_prof = list_points(iprof)
          do ichan = 1, nchannels
             ichanprof = ichanprof + 1_jpim
             do ilay = 1, rttov_atm%nlayers
                if(rttov_atm%reff(idx_prof,ilay) > 0) then
                   idx_reff = minloc(abs(cld%mie%radius(:) - rttov_atm%reff(idx_prof,ilay)), 1)
                   cld_opt_param%abs(ilay, ichanprof) = cld%mie%Cabs(idx_reff, ichan) * 1E-12 * rttov_atm%cdnc(idx_prof,ilay) * 1E3 ! converting the absorption coefficient from cm2 to km-1
                   cld_opt_param%sca(ilay, ichanprof) = cld%mie%Csca(idx_reff, ichan) * 1E-12 * rttov_atm%cdnc(idx_prof,ilay) * 1E3 ! converting the scattering coefficient from cm2 to km-1
                   cld_opt_param%legcoef(1:nleg, ilay, ichanprof) = cld%mie%legcoef(1:nleg, idx_reff, ichan)
                   cld_opt_param%pha(:, ilay, ichanprof) = cld%mie%pha(:, idx_reff, ichan)
                end if
             end do
          enddo
       enddo

       ! If doing solar calculations pre-calculate some phase angle data for scattering calculations
       if (opts % rt_ir % addsolar) then
          call rttov_init_opt_param(errorstatus, opts, cld_opt_param)
          call check_rttov_status(errorstatus, "rttov_init_opt_param")
       end if

    else

       do iprof = 1, nprof
          idx_prof = list_points(iprof)

          ! Used by OPAC
          profiles(iprof)%cloud(1,:) = rttov_atm%lwc(idx_prof,:)*1E3 !(kg/m3)

          ! Ice cloud input profiles
          profiles(iprof)%ice_scheme = rttov_opt%ice_scheme !Cloud ice water scheme: 1=Baum; 2=Baran 2014; 3=Baran 2018
          profiles(iprof)%cloud(6,:) = rttov_atm%iwc(idx_prof,:)*1E3 !(kg/m3 -> g/m3)

          ! Liquid cloud input profiles
          profiles(:)%clw_scheme = rttov_opt%clw_scheme !Cloud liquid water scheme: 1=OPAC; 2=“Deff”
          profiles(iprof)%clwde(:) = rttov_atm%reff(idx_prof,:)*2.0 ! Need the diameter

       end do
    end if


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
    call check_rttov_status(errorstatus, "rttov_get_emis")

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
       call check_rttov_status(errorstatus, "error reading BRDF atlas")

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

    ! Let RTTOV provide diffuse surface reflectances
    reflectance(:) % diffuse_refl_in = 0._jprb

    !!--------------------------------------------------------------------------------------------------------------------------!!
    !! 7. Call RTTOV K, TL or direct model                                                                                      !!
    !!--------------------------------------------------------------------------------------------------------------------------!!
    if(s3com%jac%do_jacobian_calc) then     !! Jacobian model

       ! Initialise RTTOV Jacobian structures to zero
       call rttov_init_prof(profiles_k(:))
       call rttov_init_rad(radiance_k)
       call rttov_init_transmission(transmission_k)
       call rttov_init_emis_refl(emissivity_k, reflectance_k)

       ! Set input perturbation in radiance_k:
       radiance_k%total(:) = 1._jprb
       radiance_k%bt(:) = 1._jprb

       if (nthreads <= 1) then
          call rttov_k(                       &
               errorstatus,                   & ! out   error flag
               chanprof,                      & ! in    channel and profile index structure
               opts,                          & ! in    options structure
               profiles,                      & ! in    profile array
               profiles_k,                    & ! inout Jacobian array
               coefs,                         & ! in    coefficients structure
               transmission,                  & ! inout computed transmittances
               transmission_k,                & ! inout transmittance Jacobians
               radiance,                      & ! inout computed radiances
               radiance_k,                    & ! inout input radiance/BT perturbation
               calcemis      = calcemis,      & ! in    flag for internal emissivity calcs
               emissivity    = emissivity,    & ! inout input/output emissivities per channel
               emissivity_k  = emissivity_k,  & ! inout emissivity Jacobians
               calcrefl      = calcrefl,      & ! in    flag for internal BRDF calcs
               reflectance   = reflectance,   & ! inout input/output BRDFs per channel
               reflectance_k = reflectance_k)   ! inout BRDF Jacobians
       else
          call rttov_parallel_k(              &
               errorstatus,                   & ! out   error flag
               chanprof,                      & ! in    channel and profile index structure
               opts,                          & ! in     options structure
               profiles,                      & ! in    profile array
               profiles_k,                    & ! inout Jacobian array
               coefs,                         & ! in    coefficients structure
               transmission,                  & ! inout computed transmittances
               transmission_k,                & ! inout transmittance Jacobians
               radiance,                      & ! inout computed radiances
               radiance_k,                    & ! inout input radiance/BT perturbation
               calcemis      = calcemis,      & ! in    flag for internal emissivity calcs
               emissivity    = emissivity,    & ! inout input/output emissivities per channel
               emissivity_k  = emissivity_k,  & ! inout emissivity Jacobians
               calcrefl      = calcrefl,      & ! in    flag for internal BRDF calcs
               reflectance   = reflectance,   & ! inout input/output BRDFs per channel
               reflectance_k = reflectance_k, & ! inout BRDF Jacobians
               nthreads      = nthreads)        ! in    number of threads to use
       endif
       call check_rttov_status(errorstatus, "rttov_k")
    else if (s3com%k_tl%do_k_tl_calc) then !! TL model
       call rttov_init_prof(profiles_tl(:))
       call rttov_init_rad(radiance_tl)
       call rttov_init_transmission(transmission_tl)
       call rttov_init_emis_refl(emissivity_tl, reflectance_tl)

       do iprof = 1, nprof
          idx_prof = list_points(iprof)
          joff = (iprof-1_jpim) * nchannels
          ichan = 1
          do j = 1+joff, nchannels+joff
             do ilev = 1, nlevels
                profiles_tl(iprof)%t(ilev) = tl_perturbation * profiles(iprof)%t(ilev)

                !! Call the tangent linear model
                if (nthreads <= 1) then
                   call rttov_tl(                        &
                        errorstatus,                     &
                        chanprof,                        &
                        opts,                            &
                        profiles,                        &
                        profiles_tl,                     &
                        coefs,                           &
                        transmission,                    &
                        transmission_tl,                 &
                        radiance,                        &
                        radiance_tl,                     &
                        calcemis       = calcemis,       &
                        emissivity     = emissivity,     &
                        emissivity_tl  = emissivity_tl,  &
                        calcrefl       = calcrefl,       &
                        reflectance    = reflectance,    &
                        reflectance_tl = reflectance_tl)
                else
                   call rttov_parallel_tl(               &
                        errorstatus,                     &
                        chanprof,                        &
                        opts,                            &
                        profiles,                        &
                        profiles_tl,                     &
                        coefs,                           &
                        transmission,                    &
                        transmission_tl,                 &
                        radiance,                        &
                        radiance_tl,                     &
                        calcemis       = calcemis,       &
                        emissivity     = emissivity,     &
                        emissivity_tl  = emissivity_tl,  &
                        calcrefl       = calcrefl,       &
                        reflectance    = reflectance,    &
                        reflectance_tl = reflectance_tl, &
                        nthreads       = nthreads)
                endif
                call check_rttov_status(errorstatus, 'rttov_tl')

                s3com%k_tl%t(idx_prof,ilev,ichan) = real(radiance_tl%total(j)*coefs%coef%ff_cwn(chanprof(j)%chan)**2*1E-7, wp) / &
                     real(profiles_tl(iprof)%t(ilev), wp)
                write(6,*) s3com%k_tl%t(idx_prof,ilev,ichan)
                stop
                profiles_tl(iprof)%t(ilev) = 0._jprb
             enddo
             ichan = ichan + 1
          enddo
       enddo
    else  ! Direct model
       if (nthreads <= 1) then
          call rttov_direct(              &
               errorstatus,               & ! out   error flag
               chanprof,                  & ! in    channel and profile index structure
               opts,                      & ! in    options structure
               profiles,                  & ! in    profile array
               coefs,                     & ! in    coefficients structure
               transmission,              & ! inout computed transmittances
               radiance,                  & ! inout computed radiances
               calcemis    = calcemis,    & ! in    flag for internal emissivity calcs
               emissivity  = emissivity,  & ! inout input/output emissivities per channel
               calcrefl    = calcrefl,    & ! in    flag for internal BRDF calcs
               cld_opt_param = cld_opt_param, &! in    cloud optical parameters
               reflectance = reflectance)  ! inout input/output BRDFs per channel
       else
          call rttov_parallel_direct(     &
               errorstatus,               & ! out   error flag
               chanprof,                  & ! in    channel and profile index structure
               opts,                      & ! in    options structure
               profiles,                  & ! in    profile array
               coefs,                     & ! in    coefficients structure
               transmission,              & ! inout computed transmittances
               radiance,                  & ! inout computed radiances
               calcemis    = calcemis,    & ! in    flag for internal emissivity calcs
               emissivity  = emissivity,  & ! inout input/output emissivities per channel
               calcrefl    = calcrefl,    & ! in    flag for internal BRDF calcs
               reflectance = reflectance, & ! inout input/output BRDFs per channel
               cld_opt_param = cld_opt_param, &! in    cloud optical parameters
               nthreads    = nthreads)      ! in    number of threads to use
       endif
       call check_rttov_status(errorstatus, 'rttov_direct')
    endif

    do iprof = 1, nprof
       idx_prof = list_points(iprof)
       joff = (iprof-1_jpim) * nchannels
       ichan = 1

       do j = 1+joff, nchannels+joff
          ! Radiances, BTs and BRFs
          s3com%rad%f_rad_total(idx_prof,ichan) = real(radiance%total(j) * coefs%coef%ff_cwn(chanprof(j)%chan)**2._wp * 1E-7, wp) !(W/m2/sr/um)
          s3com%rad%f_rad_clear(idx_prof,ichan) = real(radiance%clear(j)*coefs%coef%ff_cwn(chanprof(j)%chan)**2._wp * 1E-7, wp) !(W/m2/sr/um)
          s3com%rad%f_bt_total(idx_prof,ichan)  = real(radiance%bt(j), wp)
          s3com%rad%f_bt_clear(idx_prof,ichan)  = real(radiance%bt_clear(j), wp)
          s3com%rad%f_ref_total(idx_prof,ichan) = real(radiance%refl(j), wp)
          s3com%rad%f_ref_clear(idx_prof,ichan) = real(radiance%refl_clear(j), wp)

          if(s3com%jac%do_jacobian_calc) then
             do ilev = 1, profiles_k(ichan)%nlevels
                s3com%jac%t(idx_prof,ilev,ichan) = real(profiles_k(ichan)%t(ilev)*coefs%coef%ff_cwn(chanprof(j)%chan)**2*1E-7, wp)
             enddo
          end if

          ichan = ichan + 1
       enddo
    enddo



    !> Deallocate all RTTOV arrays and structures
    !> ----------------------------------------------------------------------------------------------------
    if(s3com%jac%do_jacobian_calc) then ! K model
       call rttov_alloc_k(                &
            errorstatus,                  &
            0_jpim,                       & !0 => deallocate
            nprof,                        &
            nchanprof,                    &
            nlevels,                      &
            chanprof,                     &
            opts,                         &
            profiles,                     &
            profiles_k,                   &
            coefs,                        &
            transmission,                 &
            transmission_k,               &
            radiance,                     &
            radiance_k,                   &
            calcemis      = calcemis,     &
            emissivity    = emissivity,   &
            emissivity_k  = emissivity_k, &
            calcrefl      = calcrefl,     &
            reflectance   = reflectance,  &
            reflectance_k = reflectance_k)
    else if (s3com%k_tl%do_k_tl_calc) then  !! TL model
       call rttov_alloc_tl(                 &
            errorstatus,                    &
            0_jpim,                         & !0 => deallocate
            nprof,                          &
            nchanprof,                      &
            nlevels,                        &
            chanprof,                       &
            opts,                           &
            profiles,                       &
            profiles_tl,                    &
            coefs,                          &
            transmission,                   &
            transmission_tl,                &
            radiance,                       &
            radiance_tl,                    &
            calcemis       = calcemis,      &
            emissivity     = emissivity,    &
            emissivity_tl  = emissivity_tl, &
            calcrefl       = calcrefl,      &
            reflectance    = reflectance,   &
            reflectance_tl = reflectance_tl)
    else !! Direct
       call rttov_alloc_direct(         &
            errorstatus,                &
            0_jpim,                     & !0 => deallocate
            nprof,                      &
            nchanprof,                  &
            nlevels,                    &
            chanprof,                   &
            opts,                       &
            profiles,                   &
            coefs,                      &
            transmission,               &
            radiance,                   &
            calcemis      = calcemis,   &
            emissivity    = emissivity, &
            calcrefl      = calcrefl,   &
            reflectance   = reflectance, &
            cld_maxnmom = max_mom, &
            cld_nphangle = cld%mie%nang,      &
            cld_opt_param = cld_opt_param)
    endif
    call check_rttov_status(errorstatus, 'rttov_dealloc')
    !> ----------------------------------------------------------------------------------------------------

    if (dealloc) then
       call rttov_dealloc_coefs(errorstatus, coefs)
       call check_rttov_status(errorstatus, 'rttov_dealloc_coefs')

       call rttov_deallocate_emis_atlas(emis_atlas) !Deallocate emissivity atlas

       if (opts%rt_ir%addsolar) then
          call rttov_deallocate_brdf_atlas(brdf_atlas) !Deallocate BRDF atlas
       endif
    endif

    deallocate(list_points)

   end subroutine run_rttov

end module mod_rttov

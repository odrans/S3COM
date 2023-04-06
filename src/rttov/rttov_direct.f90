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

! Initialize a few structures for RTTOV (opts, coefs, emis_atlas, brdf_atlas)
module mod_rttov_direct

  use s3com_types,  only: wp, type_cld
  use mod_rttov_utils, only: check_rttov_status
  use mod_rttov_opts, only: opts
  use mod_rttov_coefs, only: coefs
  use mod_rttov_atlas, only: emis_atlas, brdf_atlas
  use rttov_types, only: &
       rttov_profile,      &!
       rttov_transmission, &!
       rttov_radiance,     &!
       rttov_chanprof,     &!
       rttov_emissivity,   &!
       rttov_reflectance,  &!
       rttov_opt_param!


  use parkind1, only: jpim, jprb, jplm

  implicit none

  private
  public :: rttov_direct_init, rttov_direct_run, rttov_direct_free, &
       profiles, chanprof, transmission, &
       radiance, calcemis, &
       calcrefl, reflectance, cld_opt_param, emissivity


#include "rttov_direct.interface"
#include "rttov_parallel_direct.interface"
#include "rttov_alloc_direct.interface"

  type(rttov_chanprof),    pointer :: chanprof(:)      => null() ! Input channel/profile list
  type(rttov_profile),     pointer :: profiles(:)      => null() ! Input profiles
  type(rttov_transmission)         :: transmission               ! Output transmittances
  type(rttov_radiance)             :: radiance                   ! Output radiances
  logical(kind=jplm),      pointer :: calcemis(:)      => null() ! Flag to indicate calculation of emissivity within RTTOV
  logical(kind=jplm),      pointer :: calcrefl(:)      => null() ! Flag to indicate calculation of BRDF within RTTOV
  type(rttov_reflectance), pointer :: reflectance(:)   => null()
  type(rttov_opt_param)            :: cld_opt_param              ! Input cloud optical parameters
  type(rttov_emissivity),  pointer :: emissivity(:)    => null() ! Input/output surface emissivity

contains

  ! Initialize the emissivity and BRDF atlas
  subroutine rttov_direct_init(nprof, nchanprof, nlevels, cld)

    integer, intent(in) :: nprof, nchanprof, nlevels
    type(type_cld), intent(in) :: cld

    integer :: errorstatus

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
            cld_maxnmom   = cld%mie%nmom,     &
            cld_nphangle  = cld%mie%nang,     &
            cld_opt_param = cld_opt_param,    &
            init          = .true._jplm)

    call check_rttov_status(errorstatus, "rttov_alloc_direct")

  end subroutine rttov_direct_init


  ! Run the RTTOV direct model
  subroutine rttov_direct_run(nthreads)

    integer, intent(in) :: nthreads

    integer :: errorstatus

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

    call check_rttov_status(errorstatus, "rttov_direct")

  end subroutine rttov_direct_run

  subroutine rttov_direct_free(nprof, nchanprof, nlevels, cld)

    integer, intent(in) :: nprof, nchanprof, nlevels
    type(type_cld), intent(in) :: cld

    integer :: errorstatus

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
         cld_maxnmom = cld%mie%nmom, &
         cld_nphangle = cld%mie%nang,      &
         cld_opt_param = cld_opt_param)

    call check_rttov_status(errorstatus, "rttov_alloc_direct")

  end subroutine rttov_direct_free


end module mod_rttov_direct

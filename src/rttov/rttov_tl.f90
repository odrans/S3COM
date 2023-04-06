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
module mod_rttov_tl

  use s3com_types,  only: wp, type_cld
  use mod_rttov_utils, only: check_rttov_status
  use mod_rttov_opts, only: opts
  use mod_rttov_coefs, only: coefs
  use mod_rttov_atlas, only: emis_atlas, brdf_atlas
  use mod_rttov_direct, only: profiles, chanprof, transmission, &
       radiance, calcemis, &
       calcrefl, reflectance, cld_opt_param, emissivity
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
  public :: rttov_tl_init, rttov_tl_run, rttov_tl_free, &
       emissivity_tl, reflectance_tl, profiles_tl, transmission_tl, &
       radiance_tl

#include "rttov_tl.interface"
#include "rttov_parallel_tl.interface"
#include "rttov_alloc_tl.interface"

  type(rttov_emissivity),  pointer :: emissivity_tl(:)  => null() ! Emissivity Jacobians
  type(rttov_reflectance), pointer :: reflectance_tl(:) => null() ! Reflectance Jacobians
  type(rttov_profile),     pointer :: profiles_tl(:)    => null() ! Output Jacobians
  type(rttov_transmission)         :: transmission_tl             ! Transmittance Jacobians
  type(rttov_radiance)             :: radiance_tl                 ! Radiance Jacobians

contains

  ! Initialize the emissivity and BRDF atlas
  subroutine rttov_tl_init(nprof, nchanprof, nlevels, cld)

    integer, intent(in) :: nprof, nchanprof, nlevels
    type(type_cld), intent(in) :: cld

    integer :: errorstatus

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

    call check_rttov_status(errorstatus, "rttov_alloc_tl")

  end subroutine rttov_tl_init

  ! Run the RTTOV direct model
  subroutine rttov_tl_run(nthreads)

    integer, intent(in) :: nthreads

    integer :: errorstatus

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

    call check_rttov_status(errorstatus, "rttov_tl")

  end subroutine rttov_tl_run

  subroutine rttov_tl_free(nprof, nchanprof, nlevels, cld)

    integer, intent(in) :: nprof, nchanprof, nlevels
    type(type_cld), intent(in) :: cld

    integer :: errorstatus

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

       call check_rttov_status(errorstatus, "rttov_alloc_tl")

  end subroutine rttov_tl_free


end module mod_rttov_tl

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

!> Initialize the RTTOV TL model
module mod_rttov_tl
   
   use s3com_types,      only: wp, type_cld
   use mod_rttov_utils,  only: check_rttov_status
   use mod_rttov_opts,   only: opts
   use mod_rttov_coefs,  only: coefs
   use mod_rttov_atlas,  only: emis_atlas, brdf_atlas
   use mod_rttov_direct, only: profiles, chanprof, transmission, radiance, calcemis, calcrefl, reflectance, cld_opt_param, &
                               emissivity
   use rttov_types,      only: rttov_profile, rttov_transmission, rttov_radiance, rttov_chanprof, rttov_emissivity, &
                               rttov_reflectance, rttov_opt_param
   use parkind1,         only: jpim, jprb, jplm
   
   implicit none
   
   private
   public :: rttov_tl_init, rttov_tl_run, rttov_tl_free, emissivity_tl, reflectance_tl, profiles_tl, transmission_tl, &
             radiance_tl, cld_opt_param_tl
   
#include "rttov_tl.interface"
#include "rttov_parallel_tl.interface"
#include "rttov_alloc_tl.interface"
   
   type(rttov_emissivity),  pointer :: emissivity_tl(:)  => null() !< inout   surface emissivities TL
   type(rttov_reflectance), pointer :: reflectance_tl(:) => null() !< inout   reflectances TL
   type(rttov_profile),     pointer :: profiles_tl(:)    => null() !< inout   profiles TL
   type(rttov_transmission)         :: transmission_tl             !< inout   transmittances TL
   type(rttov_radiance)             :: radiance_tl                 !< inout   radiances TL
   type(rttov_opt_param)            :: cld_opt_param_tl            !< in      cloud optical parameters TL
   
contains
   
   ! ============================================================================================================================
   !> @brief Initialize the RTTOV arrays and structrures for the tangent linear (TL) model
   !! @param[in] nprof       number of profiles per call to RTTOV
   !! @param[in] nchanprof   number of channels simulated per call to RTTOV
   !! @param[in] nlevels     number of profile levels
   !! @param[in] cld         cloud data structure
   subroutine rttov_tl_init(nprof, nchanprof, nlevels, cld)
      
      ! Input
      integer, intent(in) :: nprof, nchanprof, nlevels
      type(type_cld), intent(in) :: cld
      
      ! Internal
      integer :: errorstatus
      
      ! Allocate the RTTOV arrays and structures for the TL model
      call rttov_alloc_tl(                      &
           errorstatus,                         & !< out     error flag
           1_jpim,                              & !< in      1 => allocate
           nprof,                               & !< in      number of profiles per call to RTTOV
           nchanprof,                           & !< in      number of channels simulated per call to RTTOV
           nlevels,                             & !< in      number of profile levels
           chanprof,                            & !< in      channel and profile index structure
           opts,                                & !< in      options structure
           profiles,                            & !< in      profile array
           profiles_tl,                         & !< inout   TL array
           coefs,                               & !< in      coefficients structure
           transmission,                        & !< inout   computed transmittances
           transmission_tl,                     & !< inout   transmittance TL
           radiance,                            & !< inout   computed radiances
           radiance_tl,                         & !< inout   radiance TL
           calcemis         = calcemis,         & !< in      flag for internal emissivity calculations
           emissivity       = emissivity,       & !< inout   surface emissivities per channel
           emissivity_tl    = emissivity_tl,    & !< inout   surface emissivities TL
           calcrefl         = calcrefl,         & !< in      flag for internal BRDF calculations
           reflectance      = reflectance,      & !< inout   BRDFs per channel
           reflectance_tl   = reflectance_tl,   & !< inout   BRDFs TL
           cld_maxnmom      = cld%mie%nmom,     & !< in      maximum number of coefficients in Legendre expansion of cloud phase functions
           cld_nphangle     = cld%mie%nang,     & !< in      number of phase angles over which cloud phase functions are to be defined
           cld_opt_param    = cld_opt_param,    & !< in      cloud optical parameters
           cld_opt_param_tl = cld_opt_param_tl, & !< inout   cloud optical parameters TL
           init             = .true._jplm)        !< in      initialize the newly allocated structures
      
      call check_rttov_status(errorstatus, "rttov_alloc_tl")
      
   end subroutine rttov_tl_init
   ! ============================================================================================================================
   
   ! ============================================================================================================================
   !> @brief Run the RTTOV TL model
   !! @param[in] nthreads   number of threads
   subroutine rttov_tl_run(nthreads)
      
      ! Input
      integer, intent(in) :: nthreads
      
      ! Internal
      integer :: errorstatus
      
      ! Call the RTTOV TL model
      if (nthreads <= 1) then
         call rttov_tl(                        &
              errorstatus,                     & !< out     error flag
              chanprof,                        & !< in      channel and profile index structure
              opts,                            & !< in      options structure
              profiles,                        & !< in      profile array
              profiles_tl,                     & !< inout   TL array
              coefs,                           & !< in      coefficients structure
              transmission,                    & !< inout   computed transmittances
              transmission_tl,                 & !< inout   transmittance TL
              radiance,                        & !< inout   computed radiances
              radiance_tl,                     & !< inout   radiance TL
              calcemis       = calcemis,       & !< in      flag for internal emissivity calculations
              emissivity     = emissivity,     & !< inout   surface emissivities per channel
              emissivity_tl  = emissivity_tl,  & !< inout   surface emissivities TL
              calcrefl       = calcrefl,       & !< in      flag for internal BRDF calculations
              reflectance    = reflectance,    & !< inout   BRDFs per channel
              reflectance_tl = reflectance_tl)   !< inout   BRDFs TL
      else
         call rttov_parallel_tl(               &
              errorstatus,                     & !< out     error flag
              chanprof,                        & !< in      channel and profile index structure
              opts,                            & !< in      options structure
              profiles,                        & !< in      profile array
              profiles_tl,                     & !< inout   TL array
              coefs,                           & !< in      coefficients structure
              transmission,                    & !< inout   computed transmittances
              transmission_tl,                 & !< inout   transmittance TL
              radiance,                        & !< inout   computed radiances
              radiance_tl,                     & !< inout   radiance TL
              calcemis       = calcemis,       & !< in      flag for internal emissivity calculations
              emissivity     = emissivity,     & !< inout   surface emissivities per channel
              emissivity_tl  = emissivity_tl,  & !< inout   surface emissivities TL
              calcrefl       = calcrefl,       & !< in      flag for internal BRDF calculations
              reflectance    = reflectance,    & !< inout   BRDFs per channel
              reflectance_tl = reflectance_tl, & !< inout   BRDFs TL
              nthreads       = nthreads)         !< in      number of threads to use
      endif
      
      call check_rttov_status(errorstatus, "rttov_tl")
      
   end subroutine rttov_tl_run
   ! ============================================================================================================================
   
   ! ============================================================================================================================
   !> @brief Free the RTTOV arrays and structures for the TL model
   !! @param[in] nprof       number of profiles per call to RTTOV
   !! @param[in] nchanprof   number of channels simulated per call to RTTOV
   !! @param[in] nlevels     number of profile levels
   !! @param[in] cld         cloud data structure
   subroutine rttov_tl_free(nprof, nchanprof, nlevels, cld)
      
      ! Input
      integer, intent(in) :: nprof, nchanprof, nlevels
      type(type_cld), intent(in) :: cld
      
      ! Internal
      integer :: errorstatus
      
      ! Deallocate the RTTOV arrays and structures for the TL model
      call rttov_alloc_tl(                      &
           errorstatus,                         & !< out     error flag
           0_jpim,                              & !< in      0 => deallocate
           nprof,                               & !< in      number of profiles per call to RTTOV
           nchanprof,                           & !< in      number of channels simulated per call to RTTOV
           nlevels,                             & !< in      number of profile levels
           chanprof,                            & !< in      channel and profile index structure
           opts,                                & !< in      options structure
           profiles,                            & !< in      profile array
           profiles_tl,                         & !< inout   TL array
           coefs,                               & !< in      coefficients structure
           transmission,                        & !< inout   computed transmittances
           transmission_tl,                     & !< inout   transmittance TL
           radiance,                            & !< inout   computed radiances
           radiance_tl,                         & !< inout   radiance TL
           calcemis         = calcemis,         & !< in      flag for internal emissivity calculations
           emissivity       = emissivity,       & !< inout   surface emissivities per channel
           emissivity_tl    = emissivity_tl,    & !< inout   surface emissivities TL
           calcrefl         = calcrefl,         & !< in      flag for internal BRDF calculations
           reflectance      = reflectance,      & !< inout   BRDFs per channel
           reflectance_tl   = reflectance_tl,   & !< inout   BRDFs TL
           cld_maxnmom      = cld%mie%nmom,     & !< in      maximum number of coefficients in Legendre expansion of cloud phase functions
           cld_nphangle     = cld%mie%nang,     & !< in      number of phase angles over which cloud phase functions are to be defined
           cld_opt_param    = cld_opt_param,    & !< in      cloud optical parameters
           cld_opt_param_tl = cld_opt_param_tl)   !< inout   cloud optical parameters TL
           
      call check_rttov_status(errorstatus, "rttov_alloc_tl")
      
   end subroutine rttov_tl_free
   ! ============================================================================================================================
   
end module mod_rttov_tl

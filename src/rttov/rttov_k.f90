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

!> Initialize the RTTOV K model
module mod_rttov_k
   
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
   public :: rttov_k_init, rttov_k_run, rttov_k_free, emissivity_k, reflectance_k, profiles_k, transmission_k, radiance_k, &
             cld_opt_param_k
   
#include "rttov_k.interface"
#include "rttov_parallel_k.interface"
#include "rttov_alloc_k.interface"
   
   type(rttov_emissivity),  pointer :: emissivity_k(:)  => null() !< inout   surface emissivities Jacobians
   type(rttov_reflectance), pointer :: reflectance_k(:) => null() !< inout   reflectance Jacobians
   type(rttov_profile),     pointer :: profiles_k(:)    => null() !< inout   profiles Jacobians
   type(rttov_transmission)         :: transmission_k             !< inout   transmittances Jacobians
   type(rttov_radiance)             :: radiance_k                 !< inout   radiances Jacobians
   type(rttov_opt_param)            :: cld_opt_param_k            !< in      cloud optical parameters Jacobians
   
contains
   
   ! ============================================================================================================================
   !> @brief Initialize the RTTOV arrays and structrures for the Jacobian (K) model
   !! @param[in] nprof       number of profiles per call to RTTOV
   !! @param[in] nchanprof   number of channels simulated per call to RTTOV
   !! @param[in] nlevels     number of profile levels
   !! @param[in] cld         cloud data structure
   subroutine rttov_k_init(nprof, nchanprof, nlevels, cld)
      
      ! Input
      integer, intent(in) :: nprof, nchanprof, nlevels
      type(type_cld), intent(in) :: cld
      
      ! Internal
      integer :: errorstatus
      
      ! Allocate the RTTOV arrays and structures for the K model
      call rttov_alloc_k(                     &
           errorstatus,                       & !< out     error flag
           1_jpim,                            & !< in      1 => allocate
           nprof,                             & !< in      number of profiles per call to RTTOV
           nchanprof,                         & !< in      number of channels simulated per call to RTTOV
           nlevels,                           & !< in      number of profile levels
           chanprof,                          & !< in      channel and profile index structure
           opts,                              & !< in      options structure
           profiles,                          & !< in      profile array
           profiles_k,                        & !< inout   Jacobian array
           coefs,                             & !< in      coefficients structure
           transmission,                      & !< inout   computed transmittances
           transmission_k,                    & !< inout   transmittance Jacobians
           radiance,                          & !< inout   computed radiances
           radiance_k,                        & !< inout   input radiance/BT perturbation
           calcemis        = calcemis,        & !< in      flag for internal emissivity calculations
           emissivity      = emissivity,      & !< inout   surface emissivities per channel
           emissivity_k    = emissivity_k,    & !< inout   surface emissivities Jacobians
           calcrefl        = calcrefl,        & !< in      flag for internal BRDF calculations
           reflectance     = reflectance,     & !< inout   BRDFs per channel
           reflectance_k   = reflectance_k,   & !< inout   BRDFs Jacobians
           cld_maxnmom     = cld%mie%nmom,    & !< in      maximum number of coefficients in Legendre expansion of cloud phase functions
           cld_nphangle    = cld%mie%nang,    & !< in      number of phase angles over which cloud phase functions are to be defined
           cld_opt_param   = cld_opt_param,   & !< in      cloud optical parameters
           cld_opt_param_k = cld_opt_param_k, & !< inout   cloud optical parameters Jacobians
           init            = .true._jplm)       !< in      initialize the newly allocated structures
      
      call check_rttov_status(errorstatus, "rttov_alloc_k")
      
   end subroutine rttov_k_init
   ! ============================================================================================================================
   
   ! ============================================================================================================================
   !> @brief Run the RTTOV K model
   !! @param[in] nthreads   number of threads
   subroutine rttov_k_run(nthreads)
      
      ! Input
      integer, intent(in) :: nthreads
      
      ! Internal
      integer :: errorstatus
      
      ! Call the RTTOV K model
      if (nthreads <= 1) then
         call rttov_k(                           &
              errorstatus,                       & !< out     error flag
              chanprof,                          & !< in      channel and profile index structure
              opts,                              & !< in      options structure
              profiles,                          & !< in      profile array
              profiles_k,                        & !< inout   Jacobian array
              coefs,                             & !< in      coefficients structure
              transmission,                      & !< inout   computed transmittances
              transmission_k,                    & !< inout   transmittance Jacobians
              radiance,                          & !< inout   computed radiances
              radiance_k,                        & !< inout   input radiance/BT perturbation
              calcemis        = calcemis,        & !< in      flag for internal emissivity calculations
              emissivity      = emissivity,      & !< inout   surface emissivities per channel
              emissivity_k    = emissivity_k,    & !< inout   surface emissivities Jacobians
              calcrefl        = calcrefl,        & !< in      flag for internal BRDF calculations
              reflectance     = reflectance,     & !< inout   BRDFs per channel
              reflectance_k   = reflectance_k,   & !< inout   BRDF Jacobians
              cld_opt_param   = cld_opt_param,   & !< in      cloud optical properties
              cld_opt_param_k = cld_opt_param_k)   !< inout   cloud optical properties Jacobians
      else
         call rttov_parallel_k(                  &
              errorstatus,                       & !< out     error flag
              chanprof,                          & !< in      channel and profile index structure
              opts,                              & !< in      options structure
              profiles,                          & !< in      profile array
              profiles_k,                        & !< inout   Jacobian array
              coefs,                             & !< in      coefficients structure
              transmission,                      & !< inout   computed transmittances
              transmission_k,                    & !< inout   transmittance Jacobians
              radiance,                          & !< inout   computed radiances
              radiance_k,                        & !< inout   input radiance/BT perturbation
              calcemis        = calcemis,        & !< in      flag for internal emissivity calculations
              emissivity      = emissivity,      & !< inout   surface emissivities per channel
              emissivity_k    = emissivity_k,    & !< inout   surface emissivities Jacobians
              calcrefl        = calcrefl,        & !< in      flag for internal BRDF calculations
              reflectance     = reflectance,     & !< inout   BRDFs per channel
              reflectance_k   = reflectance_k,   & !< inout   BRDF Jacobians
              cld_opt_param   = cld_opt_param,   & !< in      cloud optical properties
              cld_opt_param_k = cld_opt_param_k, & !< inout   cloud optical properties Jacobians
              nthreads        = nthreads)          !< in      number of threads to use
      endif
      
      call check_rttov_status(errorstatus, "rttov_k")
      
   end subroutine rttov_k_run
   ! ============================================================================================================================
   
   ! ============================================================================================================================
   !> @brief Free the RTTOV arrays and structures for the K model
   !! @param[in] nprof       number of profiles per call to RTTOV
   !! @param[in] nchanprof   number of channels simulated per call to RTTOV
   !! @param[in] nlevels     number of profile levels
   !! @param[in] cld         cloud data structure
   subroutine rttov_k_free(nprof, nchanprof, nlevels, cld)
      
      ! Input
      integer, intent(in) :: nprof, nchanprof, nlevels
      type(type_cld), intent(in) :: cld
      
      ! Internal
      integer :: errorstatus
      
      ! Deallocate the RTTOV arrays and structures for the K model
      call rttov_alloc_k(                     &
           errorstatus,                       & !< out     error flag
           0_jpim,                            & !< in      0 => deallocate
           nprof,                             & !< in      number of profiles per call to RTTOV
           nchanprof,                         & !< in      number of channels simulated per call to RTTOV
           nlevels,                           & !< in      number of profile levels
           chanprof,                          & !< in      channel and profile index structure
           opts,                              & !< in      options structure
           profiles,                          & !< in      profile array
           profiles_k,                        & !< inout   Jacobian array
           coefs,                             & !< in      coefficients structure
           transmission,                      & !< inout   computed transmittances
           transmission_k,                    & !< inout   transmittance Jacobians
           radiance,                          & !< inout   computed radiances
           radiance_k,                        & !< inout   input radiance/BT perturbation
           calcemis        = calcemis,        & !< in      flag for internal emissivity calculations
           emissivity      = emissivity,      & !< inout   surface emissivities per channel
           emissivity_k    = emissivity_k,    & !< inout   surface emissivities Jacobians
           calcrefl        = calcrefl,        & !< in      flag for internal BRDF calculations
           reflectance     = reflectance,     & !< inout   BRDFs per channel
           reflectance_k   = reflectance_k,   & !< inout   BRDF Jacobians
           cld_maxnmom     = cld%mie%nmom,    & !< in      maximum number of coefficients in Legendre expansion of cloud phase functions
           cld_nphangle    = cld%mie%nang,    & !< in      number of phase angles over which cloud phase functions are to be defined
           cld_opt_param   = cld_opt_param,   & !< in      cloud optical parameters
           cld_opt_param_k = cld_opt_param_k)   !< inout   cloud optical parameters Jacobians
      
      call check_rttov_status(errorstatus, "rttov_alloc_k")
      
   end subroutine rttov_k_free
   ! ============================================================================================================================
   
end module mod_rttov_k

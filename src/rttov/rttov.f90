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

!> Run RTTOV
module mod_rttov
   
   use s3com_types,      only: wp, type_model, type_rttov_opt, type_s3com, type_cld
   use mod_rttov_utils,  only: find_idx_rttov, find_ret_idx_rttov, check_rttov_status, get_rttov_model
   use mod_rttov_opts,   only: opts
   use mod_rttov_coefs,  only: coefs
   use mod_rttov_atlas,  only: emis_atlas, brdf_atlas
   use mod_rttov_direct, only: profiles, chanprof, transmission, radiance, calcemis, calcrefl, reflectance, cld_opt_param, &
                               emissivity, rttov_direct_init, rttov_direct_run, rttov_direct_free
   use mod_rttov_k,      only: rttov_k_init, rttov_k_run, rttov_k_free, emissivity_k, reflectance_k, profiles_k, &
                               transmission_k, radiance_k, cld_opt_param_k
   use mod_rttov_tl,     only: rttov_tl_init, rttov_tl_run, rttov_tl_free, emissivity_tl, reflectance_tl, profiles_tl, &
                               transmission_tl, radiance_tl
   
   ! rttov_const contains useful RTTOV constants
   use rttov_const, only:     &
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
   
   use rttov_types, only:  &
       rttov_options,      &
       rttov_coefs,        &
       rttov_profile,      &
       rttov_transmission, &
       rttov_radiance,     &
       rttov_chanprof,     &
       rttov_emissivity,   &
       rttov_reflectance,  &
       rttov_opt_param
   
   ! The rttov_emis_atlas_data type must be imported separately
   use mod_rttov_emis_atlas, only: rttov_emis_atlas_data, atlas_type_ir, atlas_type_mw

   ! The rttov_brdf_atlas_data type must be imported separately
   use mod_rttov_brdf_atlas, only: rttov_brdf_atlas_data
   
   use parkind1, only: jpim, jprb, jplm
   use rttov_unix_env, only: rttov_exit
   
   implicit none
   
   private
   public :: run_rttov

#include "rttov_init_opt_param.interface"
#include "rttov_read_coefs.interface"
#include "rttov_init_emis_refl.interface"
#include "rttov_init_prof.interface"
#include "rttov_init_rad.interface"
#include "rttov_init_transmission.interface"
#include "rttov_user_options_checkinput.interface"
#include "rttov_print_opts.interface"
#include "rttov_print_profile.interface"
#include "rttov_skipcommentline.interface"

! Use emissivity atlas
#include "rttov_get_emis.interface"

! Use BRDF atlas
#include "rttov_get_brdf.interface"
   
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
   
   ! ============================================================================================================================
   !> @brief Run RTTOV
   !! @param[in] rttov_atm   RTTOV atmospheric structure
   !! @param[in] rttov_opt   RTTOV options structure
   !! @param[in] cld         cloud structure
   !! @param[inout] s3com    s3com structure
   subroutine run_rttov(rttov_atm, rttov_opt, s3com, cld)
      
      ! Input
      type(type_model), intent(in) :: rttov_atm
      type(type_rttov_opt), intent(in) :: rttov_opt
      type(type_cld), intent(in) :: cld
      
      ! Input/output
      type(type_s3com), intent(inout) :: s3com
      
      ! Internal
      integer, dimension(:), allocatable :: list_points
      integer :: errorstatus, idx_prof
      character(len=8) :: model_rttov
      
      errorstatus = 0_jpim

      ! Find how many profiles are to be processed
      if (s3com%nml%flag_retrievals) then
         n_true = count(s3com%ret%flag_rttov); allocate(list_points(n_true))
         list_points = find_ret_idx_rttov(s3com) !< Inversion of liquid cloud properties only
      else
         n_true = count(s3com%flag_rttov); allocate(list_points(n_true))
         list_points = find_idx_rttov(s3com)     !< In general
      endif

      ! Set up a few useful dimensions
      nthreads  = rttov_opt%nthreads
      max_mom   = cld%mie%nmom
      nleg      = max_mom + 1
      nchannels = rttov_opt%nchannels
      nlevels   = rttov_atm%nlevels
      nlayers   = rttov_atm%nlayers
      nprof     = size(list_points)
      nchanprof = nchannels * nprof
      
      if (nprof == 0) return

      ! Allocate RTTOV input and output structures depending on the model that will be called
      ! Current options are direct, jacobian (K) and tangent linear (TL). K and TL also include the direct call.
      ! Note that there are no outputs to these functions, all arrays are allocated internally.
      ! -------------------------------------------------------------------------------------------------------------------------
      model_rttov = get_rttov_model(s3com)

      select case (model_rttov)
         case ("direct")
            call rttov_direct_init(nprof, nchanprof, nlevels, cld)
         case ("jacobian")
            call rttov_k_init(nprof, nchanprof, nlevels, cld)
         case ("TL")
            call rttov_tl_init(nprof, nchanprof, nlevels, cld)
      end select
      ! -------------------------------------------------------------------------------------------------------------------------
      
      ! Build the list of profile/channel indices in chanprof
      ! -------------------------------------------------------------------------------------------------------------------------
      nch = 0_jpim
      do j = 1, nprof
         do jch = 1, nchannels
            nch = nch + 1_jpim
            chanprof(nch)%prof = j
            chanprof(nch)%chan = rttov_opt%channel_list(jch)
         enddo
      enddo
      ! -------------------------------------------------------------------------------------------------------------------------

      ! Read profile data
      ! -------------------------------------------------------------------------------------------------------------------------
      ! Gas units for profiles
      profiles(1:nprof)%gas_units = rttov_opt%gas_units

      ! Loop over all profiles and read data for each one
      do iprof = 1, nprof
         
         idx_prof = list_points(iprof)
         
         ! Pressure @units{hPa}, temperature @units{K} and water vapour @units{kg kg-1}
         profiles(iprof)%p(:) = rttov_atm%p(idx_prof,:)*1.0E-02_wp
         profiles(iprof)%t(:) = rttov_atm%t(idx_prof,:)
         profiles(iprof)%q(:) = rttov_atm%q(idx_prof,:)
         
         ! Air variables in 2 meters
         profiles(iprof)%s2m%p = rttov_atm%ps(idx_prof)*1.0E-02_wp !< Surface pressure      @units{hPa}
         profiles(iprof)%s2m%t = rttov_atm%t_2m(idx_prof)          !< 2m temperature        @units{K}
         profiles(iprof)%s2m%q = rttov_atm%q_2m(idx_prof)          !< 2m water vapour       @units{kg kg-1}
         profiles(iprof)%s2m%u = rttov_atm%u_10m(idx_prof)         !< 10m zonal wind        @units{m s-1}
         profiles(iprof)%s2m%v = rttov_atm%v_10m(idx_prof)         !< 10m meridional wind   @units{m s-1}
         profiles(iprof)%s2m%wfetc = 100000                        !< Used typical value given in documentation
         
         ! Skin variables
         profiles(iprof)%skin%t = rttov_atm%ts(idx_prof)
         profiles(iprof)%skin%salinity = 0.0                        !< Temporary, use other typical value
         profiles(iprof)%skin%fastem = (/3.0, 5.0, 15.0, 0.1, 0.3/) !< Typical RTTOV default for land
         
         ! Surface type and water type
         if (rttov_atm%landmask(iprof) < 0.5) then
            profiles(iprof)%skin%surftype = surftype_sea
         else
            profiles(iprof)%skin%surftype = surftype_land
         endif
         
         profiles(iprof)%skin%watertype = watertype_fresh_water !< Temporary, adapt this
         
         ! Elevation @units{km}, latitude @units{degrees North} and longitude @units{degrees East}
         profiles(iprof)%elevation = rttov_atm%topography(idx_prof)*1.0E-03_wp
         profiles(iprof)%latitude  = rttov_atm%lat(idx_prof)
         profiles(iprof)%longitude = rttov_atm%lon(idx_prof)
         
         ! Satellite and solar viewing angles @units{°}
         profiles(iprof)%zenangle    = rttov_opt%zenangle
         profiles(iprof)%azangle     = rttov_opt%azangle
         profiles(iprof)%sunzenangle = rttov_atm%sunzenangle(idx_prof)
         profiles(iprof)%sunazangle  = rttov_atm%sunazangle(idx_prof)

         profiles(iprof)%mmr_cldaer = rttov_opt%mmr_cldaer !< Logical flag to set cloud and aerosol
                                                           !< Units: true => kg/kg (cld+aer); false => g/m3 (cld), cm-3 (aer)

         ! Cloud variables for simple cloud scheme, set cfraction to 0. to turn this off (VIS/IR only)
         profiles(iprof)%cfrac(:) = rttov_atm%clc(idx_prof,:)
         
         ! call rttov_print_profile (profiles(iprof))
         
      enddo

      if (opts%rt_ir%user_cld_opt_param) then

         ! Cloud variables for user defined cloud scheme
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
                     ! Converting the absorption coefficient from @units{cm^2} to @units{km-1}
                     cld_opt_param%abs(ilay,ichanprof) = cld%mie%Cabs(idx_reff,ichan)*1.0E-12_wp * &
                          rttov_atm%cdnc(idx_prof,ilay)*1.0E+03_wp
                     ! Converting the scattering coefficient from @units{cm^2} to @units{km-1}
                     cld_opt_param%sca(ilay,ichanprof) = cld%mie%Csca(idx_reff,ichan)*1.0E-12_wp * &
                          rttov_atm%cdnc(idx_prof,ilay)*1.0E+03_wp
                     cld_opt_param%legcoef(1:nleg,ilay,ichanprof) = cld%mie%legcoef(1:nleg,idx_reff,ichan)
                     cld_opt_param%pha(:,ilay,ichanprof) = cld%mie%pha(:,idx_reff,ichan)

                  endif

               enddo

            enddo
            !call rttov_print_profile (profiles(iprof))
         enddo

         if (opts%rt_ir%addsolar) then
            call rttov_init_opt_param(errorstatus, opts, cld_opt_param)
            call check_rttov_status(errorstatus, "rttov_init_opt_param")
         endif
         
      else

         ! Cloud variables for parameterized cloud scheme
         do iprof = 1, nprof
            
            idx_prof = list_points(iprof)
            
            ! Used by OPAC
            ! Converting the liquid water content from @units{kg m-3} to @units{g m-3}
            profiles(iprof)%cloud(1,:) = rttov_atm%lwc(idx_prof,:) * 1.0E+03_wp 
            
            ! Ice cloud input profiles
            ! Cloud ice water scheme: 1=Baum; 2=Baran 2014; 3=Baran 2018
            profiles(iprof)%ice_scheme = rttov_opt%ice_scheme
            
            ! Converting the ice water content from @units{kg m-3} to @units{g m-3}
            profiles(iprof)%cloud(6,:) = rttov_atm%iwc(idx_prof,:) * 1.0E+03_wp
            
            ! Liquid cloud input profiles
            profiles(:)%clw_scheme = rttov_opt%clw_scheme                  !< Cloud liquid water scheme: 1=OPAC; 2=“Deff”
            profiles(iprof)%clwde(:) = rttov_atm%reff(idx_prof,:) * 2.0_wp !< RTTOV needs the effective diameter
            
         enddo
         
      endif
      ! -------------------------------------------------------------------------------------------------------------------------

      ! Specify surface emissivity and reflectance
      ! -------------------------------------------------------------------------------------------------------------------------
      ! Use emissivity atlas
      call rttov_get_emis( &
           errorstatus,    &
           opts,           &
           chanprof,       &
           profiles,       &
           coefs,          &
           emis_atlas,     &
           emissivity(:)%emis_in)
      call check_rttov_status(errorstatus, "rttov_get_emis")
      
      ! Calculate emissivity within RTTOV where the atlas emissivity value is zero or less
      calcemis(:) = (emissivity(:)%emis_in <= 0._jprb)
      
      if (opts%rt_ir%addsolar) then
         
         ! Use BRDF atlas
         call rttov_get_brdf(         &
              errorstatus,            &
              opts,                   &
              chanprof,               &
              profiles,               &
              coefs,                  &
              brdf_atlas,             &
              reflectance(:)%refl_in)
         call check_rttov_status(errorstatus, "error reading BRDF atlas")
         
         ! Calculate BRDF within RTTOV where the atlas BRDF value is zero or less
         calcrefl(:) = (reflectance(:)%refl_in <= 0._jprb)
         
      endif
      
      ! Use the RTTOV emissivity and BRDF calculations over sea surfaces
      do j = 1, size(chanprof)
         if (profiles(chanprof(j)%prof)%skin%surftype == surftype_sea) then
            calcemis(j) = .true.
            calcrefl(j) = .true.
         endif
      enddo
      
      ! Use default cloud top BRDF for simple cloud in VIS/NIR channels
      reflectance(:)%refl_cloud_top = 0._jprb
      
      ! Let RTTOV provide diffuse surface reflectances
      reflectance(:)%diffuse_refl_in = 0._jprb
      ! -------------------------------------------------------------------------------------------------------------------------

      ! Call RTTOV direct, K or TL model
      ! -------------------------------------------------------------------------------------------------------------------------
      select case (model_rttov)
         case ("direct")
            call rttov_direct_run(nthreads)
         case ("jacobian")
            ! Initialize RTTOV Jacobian structures to zero
            call rttov_init_prof(profiles_k(:))
            call rttov_init_rad(radiance_k)
            call rttov_init_transmission(transmission_k)
            call rttov_init_emis_refl(emissivity_k, reflectance_k)
            
            ! Set input perturbation in radiance_k:
            radiance_k%total(:) = 1._jprb
            radiance_k%bt(:) = 1._jprb
            
            call rttov_k_run(nthreads)
         case ("TL")
            write(6, *) "Warning: TL mode not well tested"
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
                     
                     call rttov_tl_run(nthreads)
                     
                     s3com%k_tl%t(idx_prof,ichan,ilev) = &
                     real((radiance_tl%total(j)*coefs%coef%ff_cwn(chanprof(j)%chan)**2.0_wp*1.0E-07_wp / &
                     profiles_tl(iprof)%t(ilev)), wp)
                     !write(6,*) s3com%k_tl%t(idx_prof,ichan,ilev)
                     profiles_tl(iprof)%t(ilev) = 0._jprb
                     
                  enddo
                  
               ichan = ichan + 1
               
            enddo
            
         enddo
         
      end select
      ! -------------------------------------------------------------------------------------------------------------------------
      
      ! -------------------------------------------------------------------------------------------------------------------------
      do iprof = 1, nprof
         
         idx_prof = list_points(iprof)
         joff = (iprof-1_jpim) * nchannels
         ichan = 1
         
         do j = 1+joff, nchannels+joff
            
            ! Radiances, brightness temperatures, reflectances, surface albedo and surface emissivity
            ! Converting radiances in all-skies from @units{mW m-2 sr-1 cm-1} to @units{W m-2 sr-1 um-1}
            s3com%rad%f_rad_total(idx_prof,ichan) = &
            real(radiance%total(j)*coefs%coef%ff_cwn(chanprof(j)%chan)**2.0_wp*1.0E-07_wp, wp)
            ! Converting radiances in clear-sky from @units{mW m-2 sr-1 cm-1} to @units{W m-2 sr-1 um-1}
            s3com%rad%f_rad_clear(idx_prof,ichan) = &
            real(radiance%clear(j)*coefs%coef%ff_cwn(chanprof(j)%chan)**2.0_wp*1.0E-07_wp, wp)
            s3com%rad%f_bt_total(idx_prof,ichan)  = real(radiance%bt(j), wp)
            s3com%rad%f_bt_clear(idx_prof,ichan)  = real(radiance%bt_clear(j), wp)
            s3com%rad%f_ref_total(idx_prof,ichan) = real(radiance%refl(j), wp)
            s3com%rad%f_ref_clear(idx_prof,ichan) = real(radiance%refl_clear(j), wp)
            s3com%rad%brdf(idx_prof,ichan)        = real(reflectance(j)%refl_out, wp)
            s3com%rad%emiss(idx_prof,ichan)       = real(emissivity(j)%emis_out, wp)
            
            if (s3com%jac%do_jacobian_calc) then
               
               ! Surface Jacobians
               s3com%jac%brdf(idx_prof,ichan)  = real(reflectance_k(j)%refl_out, wp)
               s3com%jac%emiss(idx_prof,ichan) = real(emissivity_k(j)%emis_out, wp)
               
               ! Atmospheric profile Jacobians
               do ilev = 1, profiles_k(ichan)%nlevels
                  s3com%jac%t(idx_prof,ichan,ilev) = &
                  real(profiles_k(j)%t(ilev)*coefs%coef%ff_cwn(chanprof(j)%chan)**2.0_wp*1.0E-07_wp, wp)
                  s3com%jac%q(idx_prof,ichan,ilev) = &
                  real(profiles_k(j)%q(ilev)*coefs%coef%ff_cwn(chanprof(j)%chan)**2.0_wp*1.0E-07_wp, wp)
               enddo
               
               ! Cloud profile Jacobians
               do ilay = 1, profiles_k(ichan)%nlayers
                  s3com%jac%cfrac(idx_prof,ichan,ilay) = &
                  real(profiles_k(ichan)%cfrac(ilay)*coefs%coef%ff_cwn(chanprof(j)%chan)**2.0_wp*1.0E-07_wp, wp)
               enddo
               
               if (.not. opts%rt_ir%user_cld_opt_param .and. rttov_opt%clw_scheme .eq. 2) then !< Deff scheme only
                  do ilay = 1, profiles_k(ichan)%nlayers
                     s3com%jac%clwde(idx_prof,ichan,ilay) = &
                     real(profiles_k(ichan)%clwde(ilay)*coefs%coef%ff_cwn(chanprof(j)%chan)**2.0_wp*1.0E-07_wp, wp)
                  enddo
               endif
               
            endif
            
            ichan = ichan + 1
            
         enddo
         
         ! Forward model simulation vector (for retrievals)
         s3com%ret%F = s3com%rad%f_rad_total
         
      enddo
      ! -------------------------------------------------------------------------------------------------------------------------
      
      ! Deallocate all RTTOV arrays and structures
      ! -------------------------------------------------------------------------------------------------------------------------
      select case (model_rttov)
         case ("direct")
            call rttov_direct_free(nprof, nchanprof, nlevels, cld)
         case ("jacobian")
            call rttov_k_free(nprof, nchanprof, nlevels, cld)
         case ("tl")
            call rttov_tl_free(nprof, nchanprof, nlevels, cld)
      end select
      
      deallocate(list_points)
      ! -------------------------------------------------------------------------------------------------------------------------
      
   end subroutine run_rttov
   ! ============================================================================================================================
   
end module mod_rttov

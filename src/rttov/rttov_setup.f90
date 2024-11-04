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

!> Set up the RTTOV options and atmospheric structures
module mod_rttov_setup
   
   use s3com_types, only: wp, wpi, type_model
   
   implicit none
   
   private
   public :: rttov_setup_opt, rttov_setup_atm
   
contains
   
   ! ============================================================================================================================
   !> @brief Set up the RTTOV options structure
   !! @param[in] nml          namelist structure
   !! @param[in] zenangle     satellite zenith angle
   !! @param[in] azangle      satellite azimuth angle
   !! @param[out] rttov_opt   RTTOV options structure
   !! param[inout] s3com      s3com data structure
   subroutine rttov_setup_opt(nml, zenangle, azangle, rttov_opt, s3com)
      
      use s3com_types, only: type_rttov_opt, type_nml, type_s3com
      
      ! Input
      type(type_nml), intent(in) :: nml
      real(wp), intent(in) :: zenangle, azangle
      
      ! Input/output
      type(type_s3com), intent(inout) :: s3com
      
      ! Output
      type(type_rttov_opt), intent(out) :: rttov_opt
      
      rttov_opt%platform   = nml%platform
      rttov_opt%satellite  = nml%satellite
      rttov_opt%instrument = nml%instrument
      rttov_opt%dosolar    = 1
      rttov_opt%nchannels  = nml%nchannels
      rttov_opt%nthreads   = nml%rttov_nthreads
      
      rttov_opt%ir_scatt_model  = nml%ir_scatt_model
      rttov_opt%vis_scatt_model = nml%vis_scatt_model
      
      rttov_opt%do_opdep_calc = nml%do_opdep_calc
      rttov_opt%dom_nstreams  = nml%dom_nstreams
      rttov_opt%dom_nmoments  = nml%dom_nmoments
      rttov_opt%dom_rayleigh  = nml%dom_rayleigh
      
      rttov_opt%add_refrac   = nml%add_refrac
      rttov_opt%add_clouds   = nml%add_clouds
      rttov_opt%add_aerosols = nml%add_aerosols
      
      rttov_opt%gas_units  = 1 !< @units{kg kg-1}
      rttov_opt%mmr_cldaer = nml%mmr_cldaer
      
      rttov_opt%user_cld_opt_param = nml%user_cld_opt_param
      rttov_opt%ice_scheme         = nml%ice_scheme
      rttov_opt%clw_scheme         = nml%clw_scheme
      
      rttov_opt%ozone_data = nml%ozone_data
      rttov_opt%co2_data   = .false.
      rttov_opt%n2o_data   = .false.
      rttov_opt%ch4_data   = .false.
      rttov_opt%co_data    = .false.
      rttov_opt%so2_data   = .false.
      
      allocate(rttov_opt%channel_list(rttov_opt%nchannels))
      
      rttov_opt%channel_list = nml%channel_list
      rttov_opt%month        = s3com%date(2)
      rttov_opt%zenangle     = zenangle
      rttov_opt%azangle      = azangle
      
      s3com%opt%rttov = rttov_opt
      
   end subroutine rttov_setup_opt
   ! ============================================================================================================================
   
   ! ============================================================================================================================
   !> @brief Set up the RTTOV atmospheric structure
   !! @param[in] idx_start    index of starting point
   !! @param[in] idx_end      index of ending point
   !! @param[in] model        model data structure
   !! @param[out] rttov_atm   RTTOV atmospheric structure
   subroutine rttov_setup_atm(idx_start, idx_end, model, rttov_atm)
      
      ! Input
      integer, intent(in) :: idx_start, idx_end
      type(type_model), intent(in) :: model
      
      ! Output
      type(type_model), intent(out) :: rttov_atm
      
      ! Internal
      integer(wpi) :: nidx
      
      nidx = idx_end - idx_start + 1
      
      rttov_atm%npoints   = nidx
      rttov_atm%idx_start = idx_start
      rttov_atm%idx_end   = idx_end
      rttov_atm%nlevels   = model%nlevels
      rttov_atm%nlayers   = model%nlayers
      
      rttov_atm%lat         = model%lat(idx_start:idx_end)
      rttov_atm%lon         = model%lon(idx_start:idx_end)
      rttov_atm%landmask    = model%landmask(idx_start:idx_end)
      rttov_atm%topography  = model%topography(idx_start:idx_end)
      rttov_atm%ps          = model%ps(idx_start:idx_end)
      rttov_atm%ts          = model%ts(idx_start:idx_end)
      rttov_atm%t_2m        = model%t_2m(idx_start:idx_end)
      rttov_atm%q_2m        = model%q_2m(idx_start:idx_end)
      rttov_atm%u_10m       = model%u_10m(idx_start:idx_end)
      rttov_atm%v_10m       = model%v_10m(idx_start:idx_end)
      rttov_atm%iwp         = model%iwp(idx_start:idx_end)
      rttov_atm%lwp         = model%lwp(idx_start:idx_end)
      rttov_atm%cod         = model%cod(idx_start:idx_end)
      rttov_atm%sunzenangle = model%sunzenangle(idx_start:idx_end)
      rttov_atm%sunazangle  = model%sunazangle(idx_start:idx_end)
      
      rttov_atm%z    = model%z(idx_start:idx_end,:)
      rttov_atm%p    = model%p(idx_start:idx_end,:)    
      rttov_atm%t    = model%t(idx_start:idx_end,:)
      rttov_atm%q    = model%q(idx_start:idx_end,:)
      rttov_atm%clc  = model%clc(idx_start:idx_end,:)
      rttov_atm%iwc  = model%iwc(idx_start:idx_end,:)
      rttov_atm%lwc  = model%lwc(idx_start:idx_end,:)
      rttov_atm%reff = model%reff(idx_start:idx_end,:)
      rttov_atm%cdnc = model%cdnc(idx_start:idx_end,:)
      rttov_atm%dz   = model%dz(idx_start:idx_end,:)
      
   end subroutine rttov_setup_atm
   ! ============================================================================================================================
   
end module mod_rttov_setup

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

!> Inspired by https://cerfacs.fr/coop/fortran-namelist-workedex
module mod_io_namelist
   
   use, intrinsic :: iso_fortran_env, only: stderr => error_unit
   use s3com_types,    only: wpi, type_nml
   use mod_io_verbose, only: verbose_namelist
   
   implicit none
   
   private
   public :: namelist_load
   
contains
   
   ! ============================================================================================================================
   !> @brief Load namelist structure
   type(type_nml) function namelist_load() result(nml)
      
      character(len=256) :: fname_nml
      
      ! Set the namelist file
      if (command_argument_count() .ne. 1) then
         write(*, *) "Namelist not provided"
         stop
      else
         call get_command_argument(1, fname_nml)
      endif
      
      call read_namelist(fname_nml, nml)
      
      call verbose_namelist(fname_nml, nml)
      
   end function namelist_load
   ! ============================================================================================================================
   
   ! ============================================================================================================================
   !> @brief Import namelist data from the namelist file
   !! @param[in] file_path   path of the namelist file
   !! @param[out] nml        namelist structure
   subroutine read_namelist(file_path, nml)
      
      ! Input
      character(len=*), intent(in) :: file_path
      
      ! Output
      type(type_nml), intent(out) :: nml
      
      ! Internal
      integer :: file_unit, iostat, i
      
      ! Namelist variables
      character(len=256) :: fname_in, path_rttov, path_out, suffix_out, model_name, path_s3com
      logical :: flag_retrievals, flag_output_atm, flag_output_jac, flag_output_k_tl, do_jacobian_calc, do_k_tl_calc, &
                 do_opdep_calc, add_refrac, dom_rayleigh, mmr_cldaer, ozone_data, add_aerosols, add_clouds, user_cld_opt_param
      integer(wpi) :: npoints_it, nchannels, platform, satellite, instrument, ir_scatt_model, vis_scatt_model, dom_nstreams, &
                      rttov_nthreads, ice_scheme, clw_scheme, dom_nmoments
      integer(wpi), dimension(:), allocatable :: channel_list
      integer(wpi), dimension(2) :: channel_seq
      
      ! Namelist definition
      ! -------------------------------------------------------------------------------------------------------------------------
      namelist /general/   &
         path_s3com,       &
         path_out,         &
         suffix_out,       &
         fname_in,         &
         flag_retrievals,  &
         flag_output_atm,  &
         flag_output_jac,  &
         flag_output_k_tl, &
         npoints_it,       &
         nchannels,        &
         model_name
      
      namelist /rttov/ &
         platform,     &
         satellite,    &
         instrument,   &
         channel_list, &
         channel_seq
      
      namelist /rttov_init/  &
         path_rttov,         &
         do_jacobian_calc,   &
         do_opdep_calc,      &
         ir_scatt_model,     &
         vis_scatt_model,    &
         dom_nstreams,       &
         dom_nmoments,       &
         dom_rayleigh,       &
         rttov_nthreads,     &
         user_cld_opt_param, &
         ice_scheme,         &
         clw_scheme,         &
         mmr_cldaer,         &
         ozone_data,         &
         add_aerosols,       &
         add_clouds,         &
         add_refrac
      ! -------------------------------------------------------------------------------------------------------------------------
      
      ! Open the namelist file
      call open_namelist(file_path, file_unit, iostat)
      
      ! Read namelist data
      read (nml=general, iostat=iostat, unit=file_unit)
      allocate(channel_list(nchannels))
      read (nml=rttov, iostat=iostat, unit=file_unit)
      read (nml=rttov_init, iostat=iostat, unit=file_unit)
      
      ! Close the namelist file
      call close_namelist(file_path, file_unit, iostat)
      
      ! Define the channel list
      if (channel_seq(1) > 0) channel_list = (/(i, i=channel_seq(1), channel_seq(2), 1)/)
      
      ! TL model is not used by default (not well tested)
      do_k_tl_calc = .false.
      
      ! Namelist structure
      ! -------------------------------------------------------------------------------------------------------------------------
      nml%path_rttov         = path_rttov
      nml%path_s3com         = path_s3com
      nml%path_out           = path_out
      nml%suffix_out         = suffix_out
      nml%fname_in           = fname_in
      nml%npoints_it         = npoints_it
      nml%nchannels          = nchannels
      nml%model_name         = model_name
      nml%channel_list       = channel_list
      nml%platform           = platform
      nml%satellite          = satellite
      nml%instrument         = instrument
      nml%do_jacobian_calc   = do_jacobian_calc
      nml%do_k_tl_calc       = do_k_tl_calc
      nml%do_opdep_calc      = do_opdep_calc
      nml%ir_scatt_model     = ir_scatt_model
      nml%vis_scatt_model    = vis_scatt_model
      nml%dom_rayleigh       = dom_rayleigh
      nml%dom_nstreams       = dom_nstreams
      nml%dom_nmoments       = dom_nmoments
      nml%rttov_nthreads     = rttov_nthreads
      nml%flag_retrievals    = flag_retrievals
      nml%flag_output_atm    = flag_output_atm
      nml%flag_output_jac    = flag_output_jac
      nml%flag_output_k_tl   = flag_output_k_tl
      nml%ice_scheme         = ice_scheme
      nml%clw_scheme         = clw_scheme
      nml%mmr_cldaer         = mmr_cldaer
      nml%ozone_data         = ozone_data
      nml%add_aerosols       = add_aerosols
      nml%add_clouds         = add_clouds
      nml%add_refrac         = add_refrac
      nml%user_cld_opt_param = user_cld_opt_param
      ! -------------------------------------------------------------------------------------------------------------------------
      
   end subroutine read_namelist
   ! ============================================================================================================================
   
   ! ============================================================================================================================
   !> @brief Open the namelist file
   !! @param[in] file_path    path of the namelist file
   !! @param[out] file_unit   unit number
   !! @param[out] iostat      status flag
   subroutine open_namelist(file_path, file_unit, iostat)
      
      ! Input
      character(len=*), intent(in) :: file_path
      
      ! Output
      integer, intent(out) :: file_unit, iostat
      
      inquire (file=file_path, iostat=iostat)
      
      if (iostat /= 0) then
         write (stderr, '(3a)') 'Error: file "', trim(file_path), '" not found!'
      endif
      
      open (action='read', file=file_path, iostat=iostat, newunit=file_unit)
      
   end subroutine open_namelist
   ! ============================================================================================================================
   
   ! ============================================================================================================================
   !> @brief Close the namelist file
   !! @param[in] file_path    path of the namelist file
   !! @param[out] file_unit   unit number
   !! @param[out] iostat      status flag
   subroutine close_namelist(file_path, file_unit, iostat)
      
      ! Input
      character(len=*), intent(in) :: file_path
      
      ! Output
      integer, intent(in) :: file_unit, iostat
      
      ! Internal
      character(len=1000) :: line
      
      if (iostat /= 0) then
         write (stderr, '(2a)') 'Error reading file :"', trim(file_path)
         write (stderr, '(a, i0)') 'iostat was:"', iostat
         backspace(file_unit)
         read(file_unit,fmt='(A)') line
         write(stderr,'(A)') 'Invalid line : '//trim(line)
      endif
      
      close (file_unit)
      
   end subroutine close_namelist
   ! ============================================================================================================================
   
end module mod_io_namelist

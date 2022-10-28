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

module mod_io_namelist
  !! Inspired by https://cerfacs.fr/coop/fortran-namelist-workedex

  use, intrinsic :: iso_fortran_env, only: stderr => error_unit
  use s3com_types,         only: type_nml

  implicit none

  private
  public :: namelist_load

contains

  type(type_nml) function namelist_load() result(nml)

    character(len=256) :: fname_nml

    ! Set the namelist file
    if(command_argument_count().ne.1) then
       write(*,*) "Namelist not provided"
       stop
    else
       call get_command_argument(1, fname_nml)
       write(*,*) "Namelist file: ", trim(fname_nml)
    endif

    call read_namelist(fname_nml, nml)

  end function namelist_load

  subroutine read_namelist(file_path, nml)

    character(len=*),  intent(in)  :: file_path
    integer                        :: file_unit, iostat

    ! Namelist variables
    character(len=256) :: fname_in, path_rttov, path_out, suffix_out
    logical :: flag_retrievals, flag_output_atm, do_opdep_calc, addrefrac, dom_rayleigh

    integer(kind = 4) :: month, npoints_it, nchannels, platform, satellite, instrument, &
         ir_scatt_model, vis_scatt_model, dom_nstreams, rttov_nthreads
    integer(kind = 4), dimension(:), allocatable :: channel_list

    type(type_nml), intent(out)        :: nml

    ! Namelist definition===============================
    namelist /general/ &
         fname_in, &
         path_out, &
         suffix_out, &
         month, &
         flag_retrievals, &
         flag_output_atm, &
         npoints_it, &
         nchannels

    namelist /rttov_init/ &
         path_rttov, &
         do_opdep_calc, &
         addrefrac, &
         ir_scatt_model, &
         vis_scatt_model, &
         dom_nstreams, &
         dom_rayleigh, &
         rttov_nthreads

    namelist /rttov/ &
         channel_list, &
         platform, &
         satellite, &
         instrument
    ! Namelist definition===============================

    call open_namelist(file_path, file_unit, iostat)
    if (iostat /= 0) then
       !! write here what to do if opening failed"
       return
    end if

    read (nml=general, iostat=iostat, unit=file_unit)
    allocate(channel_list(nchannels))

    read (nml=rttov_init, iostat=iostat, unit=file_unit)
    read (nml=rttov, iostat=iostat, unit=file_unit)

    call close_namelist(file_path, file_unit, iostat)
    if (iostat /= 0) then
       !! write here what to do if reading failed"
       return
    end if

    nml%path_rttov = path_rttov
    nml%path_out = path_out
    nml%suffix_out = suffix_out
    nml%fname_in = fname_in
    nml%flag_retrievals = flag_retrievals
    nml%month = month
    nml%npoints_it = npoints_it
    nml%nchannels = nchannels
    nml%channel_list = channel_list
    nml%platform = platform
    nml%satellite = satellite
    nml%instrument = instrument
    nml%do_opdep_calc = do_opdep_calc
    nml%addrefrac = addrefrac
    nml%ir_scatt_model = ir_scatt_model
    nml%vis_scatt_model = vis_scatt_model
    nml%dom_rayleigh = dom_rayleigh
    nml%dom_nstreams = dom_nstreams
    nml%rttov_nthreads = rttov_nthreads
    nml%flag_output_atm = flag_output_atm

  end subroutine read_namelist


  !! Namelist helpers

  subroutine open_namelist(file_path, file_unit, iostat)
    !! Check whether file exists, with consitent error message
    !! return the file unit
    character(len=*),  intent(in)  :: file_path
    integer,  intent(out) :: file_unit, iostat

    inquire (file=file_path, iostat=iostat)
    if (iostat /= 0) then
       write (stderr, '(3a)') 'Error: file "', trim(file_path), '" not found!'
    end if
    open (action='read', file=file_path, iostat=iostat, newunit=file_unit)
  end subroutine open_namelist

  subroutine close_namelist(file_path, file_unit, iostat)
    !! Check the reading was OK
    !! return error line IF not
    !! close the unit
    character(len=*),  intent(in)  :: file_path
    character(len=1000) :: line
    integer,  intent(in) :: file_unit, iostat

    if (iostat /= 0) then
       write (stderr, '(2a)') 'Error reading file :"', trim(file_path)
       write (stderr, '(a, i0)') 'iostat was:"', iostat
       backspace(file_unit)
       read(file_unit,fmt='(A)') line
       write(stderr,'(A)') &
            'Invalid line : '//trim(line)
    end if
    close (file_unit)
  end subroutine close_namelist

end module mod_io_namelist

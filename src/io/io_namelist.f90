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

MODULE mod_io_namelist
  !! Inspired from https://cerfacs.fr/coop/fortran-namelist-workedex

USE, INTRINSIC :: iso_fortran_env, ONLY: stderr => error_unit

CONTAINS

  subroutine read_namelist(file_path, nml)
    !! Read some parmeters,  Here we use a namelist
    !! but if you were to change the storage format (TOML,or home-made),
    !! this signature would not change

    USE s3com_types,         ONLY: type_nml

    character(len=*),  intent(in)  :: file_path
    integer                        :: file_unit, iostat

    ! Namelist variables
    character(len=256) :: fname_out, fname_in
    logical :: flag_retrievals
    integer(kind = 4) :: month, npoints_it

    TYPE(type_nml), intent(out)        :: nml

    ! Namelist definition===============================
    namelist /ORDER/ &
         fname_out, &
         fname_in, &
         month, &
         flag_retrievals, &
         npoints_it

    fname_out = "undefined"
    fname_in = "undefined"
    month = -999
    npoints_it = 1
    flag_retrievals = .FALSE.
    ! Namelist definition===============================

    call open_namelist(file_path, file_unit, iostat)
    if (iostat /= 0) then
       !! write here what to do if opening failed"
       return
    end if

    read (nml=ORDER, iostat=iostat, unit=file_unit)
    call close_namelist(file_path, file_unit, iostat)
    if (iostat /= 0) then
       !! write here what to do if reading failed"
       return
    end if

    nml%fname_out = fname_out
    nml%fname_in = fname_in
    nml%flag_retrievals = flag_retrievals
    nml%month = month
    nml%npoints_it = npoints_it

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

END MODULE mod_io_namelist

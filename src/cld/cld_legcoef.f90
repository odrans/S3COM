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

!> Module dealing with computing and loading the coefficients of the Legendre polynomials
module mod_cld_legcoef

  use netcdf
  use rttov_const, only :errorstatus_success
  use parkind1, ONLY : jprb
  use s3com_types, only: wp, type_cld_mie, type_cld
  use mod_io_utils, only: map_ll_to_point, check_netcdf_status

  implicit none

#include "rttov_legcoef_calc.interface"

  private
  public :: cld_legcoef_load

contains

  subroutine cld_legcoef_load(cld_mie)

    type(type_cld_mie) :: cld_mie

    logical :: legcoef_exist

    cld_mie%fn_legcoef = fn_legcoef(cld_mie%fn_mie, cld_mie%nmom)

    inquire(file=cld_mie%fn_legcoef, exist=legcoef_exist)

    if(legcoef_exist) then
       call cld_legcoef_read(cld_mie)
    else
       write(*,"(A)") "Legendre coefficients are not available, computing them (might take a while)"
       call cld_legcoef_compute(cld_mie)
       call cld_legcoef_write(cld_mie)
    end if

  end subroutine cld_legcoef_load
  
  ! Compute Legendre coefficients from phase function
  subroutine cld_legcoef_compute(cld_mie)

    type(type_cld_mie) :: cld_mie

    integer :: wvl_id, size_id, status, nleg
    real(jprb), allocatable, dimension(:) :: pha, angle, legcoef

    nleg = cld_mie%nmom + 1
    allocate(cld_mie%mom(nleg))

    allocate(cld_mie%legcoef(nleg, cld_mie%nrad, cld_mie%nchan))

    allocate(pha(cld_mie%nang))
    allocate(angle(cld_mie%nang))
    allocate(legcoef(nleg))

    angle = cld_mie%angle

    do wvl_id = 1, cld_mie%nchan
       write(*,*) wvl_id
       do size_id = 1, cld_mie%nrad
          pha(:) = cld_mie%pha(:, size_id, wvl_id)
          call rttov_legcoef_calc(status, pha, angle, cld_mie%nmom, legcoef)
          cld_mie%legcoef(:, size_id, wvl_id) = real(legcoef, wp)
          if (status /= errorstatus_success) write(*,*) "rttov_legcoef_calc failed"
       end do
    end do

    deallocate(pha, angle, legcoef)

  end subroutine cld_legcoef_compute

  ! Read the Legendre data
  subroutine cld_legcoef_read(cld_mie)

    integer :: ncid
    type(type_cld_mie) :: cld_mie

    integer :: varid
    integer :: status, nleg

    nleg = cld_mie%nmom + 1

    ! Open the file
    status = nf90_open(cld_mie%fn_legcoef, nf90_nowrite, ncid)
    call check_netcdf_status(status, "nf90_open")

    allocate(cld_mie%legcoef(nleg, cld_mie%nrad, cld_mie%nchan))
    status = nf90_inq_varid(ncid, "legcoef", varid)
    call check_netcdf_status(status, "nf90_inq_varid")
    status = nf90_get_var(ncid, varid, cld_mie%legcoef)

    status = nf90_close(ncid)
    call check_netcdf_status(status, "nf90_close")

  end subroutine cld_legcoef_read


  ! Write the Legendre data to a netCDF file
  subroutine cld_legcoef_write(cld_mie)

    integer :: ncid
    type(type_cld_mie), intent(in) :: cld_mie

    integer :: varid_legcoef, dimid_nmom, dimid_size, dimid_chan, dimids_legcoef(3)
    integer :: status, nleg

    nleg = cld_mie%nmom + 1

    status = nf90_create(cld_mie%fn_legcoef, nf90_netcdf4, ncid)
    call check_netcdf_status(status, "nf90_inq_dimid")

    status = nf90_def_dim(ncid, "mom", nleg, dimid_nmom)
    call check_netcdf_status(status, "nf90_def_dim")

    status = nf90_def_dim(ncid, "size", cld_mie%nrad, dimid_size)
    call check_netcdf_status(status, "nf90_def_dim")

    status = nf90_def_dim(ncid, "chan", cld_mie%nchan, dimid_chan)
    call check_netcdf_status(status, "nf90_def_dim")

    dimids_legcoef = (/ dimid_nmom, dimid_size, dimid_chan /)
    status = nf90_def_var(ncid, "legcoef", nf90_real, dimids_legcoef, varid_legcoef)
    call check_netcdf_status(status, "nf90_def_var")
    status = nf90_put_var(ncid, varid_legcoef, cld_mie%legcoef)
    call check_netcdf_status(status, "nf90_put_var")

    status = nf90_close(ncid)
    call check_netcdf_status(status, "nf90_close")

  end subroutine cld_legcoef_write

  ! Output the name of the netcdf file containing the legendre coeficients
  function fn_legcoef(fn, nmom) result(fn_new)

    character(len=*), intent(in) :: fn
    integer, intent(in) :: nmom

    character(len=:), allocatable :: fn_new
    character(len=10) :: nmom_str
    integer :: len_fn, pos_nc

    ! Convert nmom to a string
    write(nmom_str, "(I0)") nmom

    ! Get the length of the fn
    len_fn = len_trim(fn)

    ! Find the position of ".nc" in the fn
    pos_nc = index(fn, ".nc")

    if (pos_nc > 0) then
       ! Modify the fn by inserting nmom just before ".nc"
       fn_new = fn(1:pos_nc-1) // "_legcoef_" // trim(nmom_str) // ".nc"
    else
       ! ".nc" extension not found, append nmom and ".nc" to the fn
       fn_new = trim(fn) // "_legcoef_" // trim(nmom_str) // ".nc"
    end if

  end function fn_legcoef




end module mod_cld_legcoef

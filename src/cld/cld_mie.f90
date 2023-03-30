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

module mod_cld_mie

  use netcdf
  use rttov_const, only :errorstatus_success
  use parkind1, ONLY : jprb
  use s3com_types, only: type_cld_mie, type_cld

  implicit none

#include "rttov_legcoef_calc.interface"

  private
  public :: cld_mie_load

contains

  subroutine check_netcdf_status(status, location)
    integer, intent(in) :: status
    character(len=*), intent(in) :: location

    if (status /= nf90_noerr) then
       write(*,*) "Error at ", location, ": ", trim(nf90_strerror(status))
       stop
    end if
  end subroutine check_netcdf_status

  ! Load the Mie scattering optical properties
  subroutine cld_mie_load()

    type(type_cld) :: cld
    type(type_cld_mie) :: cld_mie

    call cld_mie_read(cld_mie)

    call cld_mie_legcoef(cld_mie)

    cld%mie = cld_mie

  end subroutine cld_mie_load

  ! Read the Mie scattering data from a netCDF file
  subroutine cld_mie_read(cld_mie)

    integer :: ncid
    type(type_cld_mie) :: cld_mie

    integer :: varid
    integer :: status

    ! Open the file
    status = nf90_open("/home/b/b380333/mie_scat_liqcld_eos_2_modis.nc", nf90_nowrite, ncid)
    call check_netcdf_status(status, "nf90_open")

    ! Read the number of wavelengths
    status = nf90_inq_varid(ncid, "chan", varid)
    call check_netcdf_status(status, "nf90_inq_varid")
    status = nf90_inquire_dimension(ncid, varid, len=cld_mie%nchan)
    call check_netcdf_status(status, "nf90_inquire_dimension")

    ! Read the number of sizes
    status = nf90_inq_varid(ncid, "size", varid)
    call check_netcdf_status(status, "nf90_inq_varid")
    status = nf90_inquire_dimension(ncid, varid, len=cld_mie%nrad)
    call check_netcdf_status(status, "nf90_inquire_dimension")

    ! Read the number of angles
    status = nf90_inq_varid(ncid, "angle", varid)
    call check_netcdf_status(status, "nf90_inq_varid")
    status = nf90_inquire_dimension(ncid, varid, len=cld_mie%nang)
    call check_netcdf_status(status, "nf90_inquire_dimension")

    ! Read the wavelengths
    allocate(cld_mie%chan_id(cld_mie%nchan))
    status = nf90_inq_varid(ncid, "chan", varid)
    call check_netcdf_status(status, "nf90_inq_varid")
    status = nf90_get_var(ncid, varid, cld_mie%chan_id)
    call check_netcdf_status(status, "nf90_get_var")

    ! Read the sizes
    allocate(cld_mie%radius(cld_mie%nrad))
    status = nf90_inq_varid(ncid, "size", varid)
    call check_netcdf_status(status, "nf90_inq_varid")
    status = nf90_get_var(ncid, varid, cld_mie%radius)
    call check_netcdf_status(status, "nf90_get_var")

    ! Read the angles
    allocate(cld_mie%angle(cld_mie%nang))
    status = nf90_inq_varid(ncid, "angle", varid)
    call check_netcdf_status(status, "nf90_inq_varid")
    status = nf90_get_var(ncid, varid, cld_mie%angle)
    call check_netcdf_status(status, "nf90_get_var")

    ! Read the extinction coefficient
    allocate(cld_mie%Cext(cld_mie%nrad, cld_mie%nchan))
    status = nf90_inq_varid(ncid, "Cext", varid)
    call check_netcdf_status(status, "nf90_inq_varid")
    status = nf90_get_var(ncid, varid, cld_mie%Cext)
    call check_netcdf_status(status, "nf90_get_var")

    ! Read the phase function
    allocate(cld_mie%pha(cld_mie%nang, cld_mie%nrad, cld_mie%nchan))
    status = nf90_inq_varid(ncid, "phase_function", varid)
    call check_netcdf_status(status, "nf90_inq_varid")
    status = nf90_get_var(ncid, varid, cld_mie%pha)

    ! Close the file
    status = nf90_close(ncid)
    call check_netcdf_status(status, "nf90_close")

  end subroutine cld_mie_read

  ! Compute Legendre coefficients from phase function
  subroutine cld_mie_legcoef(cld_mie)

    type(type_cld_mie) :: cld_mie

    integer :: wvl_id, size_id, status, nleg
    real(jprb), allocatable, dimension(:) :: pha, angle, legcoef

    cld_mie%nmom = 128

    nleg = cld_mie%nmom + 1

    allocate(cld_mie%legcoef(nleg, cld_mie%nrad, cld_mie%nchan))

    allocate(pha(cld_mie%nang))
    allocate(angle(cld_mie%nang))
    allocate(legcoef(nleg))

    angle = cld_mie%angle

    do wvl_id = 1, cld_mie%nchan
       do size_id = 1, cld_mie%nrad
         pha(:) = cld_mie%pha(:, size_id, wvl_id)
         call rttov_legcoef_calc(status, pha, angle, cld_mie%nmom, legcoef)
         !! cld_mie%legcoef(:, size_id, wvl_id) = legcoef
         if (status /= errorstatus_success) write(*,*) "rttov_legcoef_calc failed"
      end do
   end do

   deallocate(pha, angle, legcoef)

 end subroutine cld_mie_legcoef


end module mod_cld_mie

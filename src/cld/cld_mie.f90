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
  use s3com_types, only: type_cld_mie, type_cld, type_nml
  use mod_io_utils, only: check_netcdf_status
  use mod_cld_legcoef, only: cld_legcoef_load

  implicit none

  private
  public :: cld_mie_load

contains

  ! Load the Mie scattering optical properties
  subroutine cld_mie_load(nml, cld)

    type(type_nml), intent(in) :: nml

    type(type_cld), intent(out) :: cld
    type(type_cld_mie) :: cld_mie

    cld_mie%nmom = nml%dom_nmoments

    cld_mie%fn_mie = "/home/b/b380333/mie_scat_liqcld_eos_2_modis.nc"

    call cld_mie_read(cld_mie)

    call cld_legcoef_load(cld_mie)

    cld%mie = cld_mie

  end subroutine cld_mie_load

  ! Read the Mie scattering data from a netCDF file
  subroutine cld_mie_read(cld_mie)

    integer :: ncid
    type(type_cld_mie) :: cld_mie

    integer :: varid
    integer :: status

    ! Open the file
    status = nf90_open(cld_mie%fn_mie, nf90_nowrite, ncid)
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

    ! Read the single scattering albedo
    allocate(cld_mie%w0(cld_mie%nrad, cld_mie%nchan))
    status = nf90_inq_varid(ncid, "w0", varid)
    call check_netcdf_status(status, "nf90_inq_varid")
    status = nf90_get_var(ncid, varid, cld_mie%w0)
    call check_netcdf_status(status, "nf90_get_var")

    ! Read the phase function
    allocate(cld_mie%pha(cld_mie%nang, cld_mie%nrad, cld_mie%nchan))
    status = nf90_inq_varid(ncid, "phase_function", varid)
    call check_netcdf_status(status, "nf90_inq_varid")
    status = nf90_get_var(ncid, varid, cld_mie%pha)

    ! Close the file
    status = nf90_close(ncid)
    call check_netcdf_status(status, "nf90_close")

    allocate(cld_mie%Csca(cld_mie%nrad, cld_mie%nchan))
    cld_mie%Csca = cld_mie%Cext * cld_mie%w0

    allocate(cld_mie%Cabs(cld_mie%nrad, cld_mie%nchan))
    cld_mie%Cabs = cld_mie%Cext * (1. - cld_mie%w0)

  end subroutine cld_mie_read




end module mod_cld_mie

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

! Description:
!> @file
!!   Allocate and load the user-defined Mie cloud properties
!!
module mod_cld_mie

  use netcdf
  use rttov_const, only :errorstatus_success
  use parkind1, ONLY : jprb
  use s3com_types, only: type_cld_mie, type_cld, type_s3com
  use mod_io_utils, only: check_netcdf_status
  use mod_cld_legcoef, only: cld_legcoef_load

  implicit none

  private
  public :: cld_mie_load

contains

  !> @brief
  !!   General call to subroutines needed to load optical properties for liquid clouds from user-defined files
  !!
  !! @details
  !!
  !!   This initializes and sets the `cld%mie` structure, which contains user-defined liquid cloud optical properties needed by RTTOV.
  !!
  !!   These properties are read from netCDF files, currently expected to be located in $PATH/data/opt_prop. Current properties are
  !!   generated from a Mie code. They are loaded for a given instrument. S3COM stops if the property files are not found.
  !!
  !!   These stored Mie optical properties later required by RTTOV are:
  !!   - the absorption cross-section (in um^2)
  !!   - the scattering cross-section (in um^2)
  !!   - the phase function for defined angles
  !!   - the legendre coefficients of the phase function
  !!
  !!   Important: Note that all cloud properties are defined for a normalized droplet size distribution n(r)!
  !!   This means that the integral of n(r) is here 1, and converting the absorption and scattering cross-sections to the
  !!   total absorption and scattering coefficients (typically in km^-1) requires multiplying by the droplet number concentration.
  !!
  !! @param[out]    cld                   cld structure (currently only contains cld%mie)
  !! @param[in]     s3com                 s3com structure
  subroutine cld_mie_load(s3com, cld)

    type(type_s3com), intent(in) :: s3com
    type(type_cld), intent(out) :: cld

    type(type_cld_mie) :: cld_mie
    character(len = 128) :: dir_mie, inst_id, fn_mie
    logical :: file_exists

    dir_mie = trim(s3com%nml%path_s3com)//"/data/cld_optprop"
    inst_id = trim(s3com%opt%rttov%platform_name)//trim(s3com%opt%rttov%sat_name)//trim(s3com%opt%rttov%inst_name)
    fn_mie = trim(dir_mie)//"/mie_scat_liqcld_"//trim(inst_id)//".nc"

    inquire(file=fn_mie, exist=file_exists)
    if (.not. file_exists) then
       write(6,*) fn_mie, "does not exist. Check or use RTTOV parameterizations for cloud properties"
    end if

    ! Read the optical properties from cld_mie%fn_mie
    call cld_mie_read(fn_mie, cld_mie)

    ! Load the corresponding legendre polynomial coeficients (on nmom moments)
    cld_mie%nmom = s3com%nml%dom_nmoments
    call cld_legcoef_load(cld_mie)

    cld%mie = cld_mie

  end subroutine cld_mie_load


  !> @brief
  !!   Initializes the cld_mie structure and loads the Mie scattering data from a netCDF file
  !!
  !! @details
  !!
  !!   This reads Mie optical properties from a user-defined NetCDF file and stores them in the `cld_mie` structure.
  !!   The following variables are set in `cld_mie`:
  !!   - nchan: number of wavelengths for the given instrument
  !!   - nrad: number of radius on which Mie properties are defined
  !!   - nang: number of angles on which the phase function is computed
  !!   - chan_id: ID of the instrument channel in RTTOV
  !!   - radius: effective radii (in um)
  !!   - angle: phase function angles (in degrees)
  !!   - Csca: the scattering cross-section coefficient (in um^2). Computed as Cext * w0, with
  !!   Cext the extinction cross-section and w0 the single-scattering albedo
  !!   - Cabs: the absorption cross-section coefficient (in um^2). Computed as Cext * (1-w0).
  !!   - pha: the phase function
  !!
  !! @param[in]     fn_mie                name of the netCDF file containing the Mie scattering data
  !! @param[out]    cld_mie               cld_mie structure
  subroutine cld_mie_read(fn_mie, cld_mie)

    character(len = 128), intent(in) :: fn_mie
    type(type_cld_mie), intent(out) :: cld_mie

    integer :: ncid, varid, status

    cld_mie%fn_mie = fn_mie

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

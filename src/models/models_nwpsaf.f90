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

!> @brief Ensemble of subroutines useful to load and process NWPSAF model simulations
!! @warning Due to NWPSAF being based on one-moment microphysics, only the cloud water content is provided.
!! Therefore, only ice clouds are currently simulated in S3COM for NWPSAF (using the Baran 2018 scheme).
!! @todo Add support for liquid clouds by finding how to reconstruct the droplet size distribution from LWC.
module mod_nwpsaf

  use s3com_types,         only: wp, type_nwpsaf
  use mod_io_utils,         only: extract_coordinates
  use mod_io_nwpsaf,       only: nwpsaf_read

  use s3com_config,        only: rd, rv, epsilon, mu, nu, a, b, Q_ext, rholiq

  implicit none

  private
  public :: nwpsaf_load, nwpsaf_free

contains

  !> @brief Load NWPSAF model simulations
  !! @details This subroutine loads NWPSAF model simulations from a NetCDF file and processes them to obtain the required variables for S3COM.
  !! @param[in] fname Input NetCDF file name containing NWPSAF simulations
  !! @param[inout] nwpsaf NWPSAF data structure
  subroutine nwpsaf_load(fname, nwpsaf)

    ! Inputs
    character(len=256), intent(in) :: fname

    ! Inputs/Outputs
    type(type_nwpsaf), intent(inout) :: nwpsaf

    ! Internal
    integer(kind=4) :: nlayers, npoints

    ! Extract the number of vertical layers and grid points in the input files
    call extract_coordinates(fname, nlayers, npoints)

    ! Initialize the NWPSAF array
    call nwpsaf_init(npoints, nlayers, nwpsaf)

    ! Read input netcdf file containing NWPSAF outputs
    call nwpsaf_read(fname, nwpsaf)

    ! Post-process NWPSAF data
    call nwpsaf_process(nwpsaf)

  end subroutine nwpsaf_load

  !> @brief Process NWPSAF data structure
  !! @details This subroutine processes the NWPSAF data structure to obtain the required variables for S3COM.
  !! @note Current approximations:
  !! - NWPSAF doesn't include 2-m specific humidity, so it is set to the specific humidity at surface level
  !! - NWPSAF doesn't include cloud effective radius, set to 0 for now. Not used for ice clouds anyway.
  !! @param[inout] nwpsaf NWPSAF data structure
  subroutine nwpsaf_process(nwpsaf)

    type(type_nwpsaf), intent(inout) :: nwpsaf

    ! Internal
    integer(kind=4) :: nlevels, nlayers, npoints

    nlevels = nwpsaf%nlevels
    nlayers = nwpsaf%nlayers
    npoints = nwpsaf%npoints

    nwpsaf%q2m(1:npoints) = nwpsaf%hum(1:npoints, nlayers) ! NWPSAF doesn't include 2-m specific humidity
    nwpsaf%paph(:,1) = 1E-2 ! RTTOV requires values strictly greater than 0

    ! NWPSAF doesn't include cloud effective radius, set to 0 for now. Not used for ice clouds anyway.
    nwpsaf%reff = 0._wp

  end subroutine nwpsaf_process

  !> @brief Initialize NWPSAF data structure
  !! @details This subroutine initializes the NWPSAF data structure by allocating the required arrays.
  !! All arrays are initialized to zero.
  !! @param[in] npoints Number of grid points
  !! @param[in] nlayers Number of vertical layers
  !! @param[inout] nwpsaf NWPSAF data structure
  subroutine nwpsaf_init(npoints, nlayers, nwpsaf)

    type(type_nwpsaf), intent(inout) :: nwpsaf
    integer(kind=4), intent(in) :: npoints, nlayers
    integer(kind=4) :: nlevels

    nlevels = nlayers + 1

    nwpsaf%npoints = npoints
    nwpsaf%nlayers = nlayers
    nwpsaf%nlevels = nlevels

    allocate(nwpsaf%height(nlayers), source = 0)
    allocate(nwpsaf%height_2(nlevels), source = 0)

    !! 2D variables
    allocate(nwpsaf%lon(npoints), source = 0._wp)
    allocate(nwpsaf%lat(npoints), source = 0._wp)
    allocate(nwpsaf%lat_orig(npoints), source = 0._wp)
    allocate(nwpsaf%lon_orig(npoints), source = 0._wp)
    allocate(nwpsaf%elevation(npoints), source = 0._wp)
    allocate(nwpsaf%lsm(npoints), source = 0._wp)
    allocate(nwpsaf%psurf(npoints), source = 0._wp)
    allocate(nwpsaf%tsurf(npoints), source = 0._wp)
    allocate(nwpsaf%t2m(npoints), source = 0._wp)
    allocate(nwpsaf%q2m(npoints), source = 0._wp)
    allocate(nwpsaf%u10(npoints), source = 0._wp)
    allocate(nwpsaf%v10(npoints), source = 0._wp)

    allocate(nwpsaf%point(npoints), source = 0)
    allocate(nwpsaf%month(npoints), source = 0)
    allocate(nwpsaf%year(npoints), source = 0)
    allocate(nwpsaf%day(npoints), source = 0)

    !! 3D variables on atmospheric levels
    allocate(nwpsaf%altitudeh(npoints, nlevels), source = 0._wp)
    allocate(nwpsaf%temph(npoints, nlevels), source = 0._wp)
    allocate(nwpsaf%paph(npoints, nlevels), source = 0._wp)
    allocate(nwpsaf%humh(npoints, nlevels), source = 0._wp)

    !! 3D variables in atmospheric layers
    allocate(nwpsaf%altitude(npoints, nlayers), source = 0._wp)
    allocate(nwpsaf%temp(npoints, nlayers), source = 0._wp)
    allocate(nwpsaf%pap(npoints, nlayers), source = 0._wp)
    allocate(nwpsaf%hum(npoints, nlayers), source = 0._wp)
    allocate(nwpsaf%cc(npoints, nlayers), source = 0._wp)
    allocate(nwpsaf%dz(npoints, nlayers), source = 0._wp)
    allocate(nwpsaf%rho(npoints, nlayers), source = 0._wp)
    allocate(nwpsaf%tv(npoints, nlayers), source = 0._wp)
    allocate(nwpsaf%lwc(npoints, nlayers), source = 0._wp)
    allocate(nwpsaf%iwc(npoints, nlayers), source = 0._wp)
    allocate(nwpsaf%cdnc(npoints, nlayers), source = 0._wp)
    allocate(nwpsaf%reff(npoints, nlayers), source = 0._wp)

  end subroutine nwpsaf_init

  !> @brief Free NWPSAF data structure
  !! @details This subroutine frees the NWPSAF data structure by deallocating the required arrays.
  !! @param[inout] nwpsaf NWPSAF data structure
  subroutine nwpsaf_free(nwpsaf)

    type(type_nwpsaf), intent(inout) :: nwpsaf

    deallocate(nwpsaf%height, nwpsaf%height_2, nwpsaf%lon, nwpsaf%lat, nwpsaf%lon_orig, nwpsaf%lat_orig, &
         nwpsaf%elevation, nwpsaf%lsm, nwpsaf%psurf, nwpsaf%tsurf, nwpsaf%t2m, nwpsaf%q2m, nwpsaf%u10, nwpsaf%v10, &
         nwpsaf%pap, nwpsaf%altitude, nwpsaf%altitudeh, nwpsaf%paph, nwpsaf%temph, nwpsaf%humh, &
         nwpsaf%temp, nwpsaf%hum, nwpsaf%cc, nwpsaf%dz, &
         nwpsaf%rho, nwpsaf%tv, nwpsaf%lwc, nwpsaf%iwc, nwpsaf%cdnc, nwpsaf%reff, &
         nwpsaf%day, nwpsaf%month, nwpsaf%year, nwpsaf%point)

  end subroutine nwpsaf_free

end module mod_nwpsaf

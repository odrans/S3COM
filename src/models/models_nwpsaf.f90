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

module mod_nwpsaf

  use s3com_types,         only: wp, type_nwpsaf
  use mod_io_icon,         only: extract_coordinates
  use mod_io_nwpsaf,       only: nwpsaf_read

  use s3com_config,        only: rd, rv, epsilon, mu, nu, a, b, Q_ext, rholiq

  implicit none

  private
  public :: nwpsaf_load, nwpsaf_clear

contains

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


  subroutine nwpsaf_process(nwpsaf)

    type(type_nwpsaf), intent(inout) :: nwpsaf

    ! Internal
    integer(kind=4) :: i, j, nlevels, nlayers, npoints

    nlevels = nwpsaf%nlevels
    nlayers = nwpsaf%nlayers
    npoints = nwpsaf%npoints

    nwpsaf%q2m(1:npoints) = nwpsaf%q(1:npoints, nlayers) ! NWPSAF doesn't include 2-m specific humidity
    nwpsaf%p_ifc(:,1) = 1E-2 ! RTTOV requires values strictly greater than 0

    !   ! Atmospheric structure
    !   ! ----------------------------------------------------------------------------------------------------

    !   ! Layer depth (approximate!! why not ifc(i+1) - ifc(i)?)
    !   nwpsaf%dz(:,1:nlayers) = abs(nwpsaf%z_ifc(:,1:nlayers) - nwpsaf%z(:,1:nlayers)) * 2._wp

    !   ! ----------------------------------------------------------------------------------------------------

    !   ! Atmospheric moisture
    !   ! ----------------------------------------------------------------------------------------------------
    !   ! Virtual temperature
    !   nwpsaf%tv = nwpsaf%t * (1 + 0.608 * nwpsaf%q)

    !   ! Moist air density (kg/m3)
    !   nwpsaf%rho = nwpsaf%p / (rd * nwpsaf%tv)

    !   ! ! Saturation vapour pressure of liquid water (Pa); Murphy and Koop (2005)
    !   ! nwpsaf%es_w = EXP(54.842763 - 6763.22 / nwpsaf%t - 4.210 * LOG(nwpsaf%t) + 0.000367 * nwpsaf%t + TANH(0.0415 * (nwpsaf%t - 218.8)) * &
    !   !      (53.878 - 1331.22 / nwpsaf%t - 9.44523 * LOG(nwpsaf%t) + 0.014025 * nwpsaf%t) )

    !   ! ! Saturation vapour pressure of ice water (Pa); Murphy and Koop (2005)
    !   ! nwpsaf%es_i = EXP(9.550426 - 5723.265 / nwpsaf%t + 3.53068 * LOG(nwpsaf%t) - 0.00728332 * nwpsaf%t)

    !   ! ----------------------------------------------------------------------------------------------------


    !   ! Cloud properties
    !   ! ----------------------------------------------------------------------------------------------------
    !   ! Cloud liquid water content (kg/m3)
    !   nwpsaf%lwc = nwpsaf%clw * nwpsaf%rho

    !   ! Cloud ice water content (kg/m3)
    !   nwpsaf%iwc = nwpsaf%cli * nwpsaf%rho

    !   ! Cloud droplet number concentration (particules/m3)
    !   nwpsaf%cdnc = nwpsaf%qnc * nwpsaf%rho

    !   !Cloud liquid water effective radius (m)
    !   do i = 1, Npoints
    !      do j = 1, Nlayers
    !         if(nwpsaf%cdnc(i,j) .gt. 0) then
    !            nwpsaf%Reff(i,j) = (a/2._wp)*(gamma((nu+1._wp+3._wp*b)/mu)/gamma((nu+1._wp+2._wp*b)/mu))*(nwpsaf%lwc(i,j) / &
    !            nwpsaf%cdnc(i,j))**b*(gamma((nu+1._wp)/mu)/gamma((nu+2._wp)/mu))**b
    !         end if
    !      end do
    !   end do
    !   nwpsaf%Reff = nwpsaf%Reff * 1E6 !! m to um (default input in RTTOV)

    !   !!Cloud extinction coefficient in m-1
    !   ! nwpsaf%beta_ext = (3._wp/4._wp)*(Q_ext/rholiq)*(nwpsaf%lwc/nwpsaf%Reff)

    !   ! !!Layer thickness
    !   ! nwpsaf%dz_cod(:,1:nlevels-1) = ((rd*nwpsaf%tv)/grav)*(log(nwpsaf%p(:,1:nlevels+1)/nwpsaf%p(:,1:nlevels)))

    !   ! !!Cloud optical depth
    !   ! nwpsaf%cod(:,1:nlevels-1) = ((nwpsaf%beta_ext(:,1:nlevels)+nwpsaf%beta_ext(:,1:nlevels+1))/2._wp)*nwpsaf%dz_cod(:,1:nlevels-1)
    !   ! ----------------------------------------------------------------------------------------------------

  end subroutine nwpsaf_process

  subroutine nwpsaf_init(npoints, nlayers, nwpsaf)

    type(type_nwpsaf) :: nwpsaf
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
    allocate(nwpsaf%lat,nwpsaf%lat_orig, nwpsaf%lon_orig, nwpsaf%elevation, nwpsaf%lsm, &
         nwpsaf%psurf, nwpsaf%tsurf, nwpsaf%t2m, nwpsaf%q2m, nwpsaf%u10, nwpsaf%v10, &
         mold = nwpsaf%lon)

    allocate(nwpsaf%point(npoints), source = 0)
    allocate(nwpsaf%month, nwpsaf%year, nwpsaf%day, &
         mold = nwpsaf%point)

    !! 3D variables on atmospheric levels
    allocate(nwpsaf%z_ifc(npoints, nlevels), source = 0._wp)
    allocate(nwpsaf%p_ifc, nwpsaf%t_ifc, nwpsaf%q_ifc, &
         mold = nwpsaf%z_ifc)

    !! 3D variables in atmospheric layers
    allocate(nwpsaf%z(npoints, nlayers), source = 0._wp)
    allocate(nwpsaf%p, nwpsaf%t, nwpsaf%q, nwpsaf%clc, nwpsaf%clw, nwpsaf%cli, &
         nwpsaf%qnc, nwpsaf%qr, nwpsaf%qs, nwpsaf%dz, nwpsaf%rho, nwpsaf%tv, nwpsaf%lwc, &
         nwpsaf%iwc, nwpsaf%cdnc, nwpsaf%reff, &
         mold = nwpsaf%z)

  end subroutine nwpsaf_init


  subroutine nwpsaf_clear(nwpsaf)

    type(type_nwpsaf), intent(inout) :: nwpsaf

    deallocate(nwpsaf%height, nwpsaf%height_2, nwpsaf%lon, nwpsaf%lat, nwpsaf%lon_orig, nwpsaf%lat_orig, &
         nwpsaf%elevation, nwpsaf%lsm, nwpsaf%psurf, nwpsaf%tsurf, nwpsaf%t2m, nwpsaf%q2m, nwpsaf%u10, nwpsaf%v10, &
         nwpsaf%p, nwpsaf%z, nwpsaf%z_ifc, nwpsaf%p_ifc, nwpsaf%t_ifc, nwpsaf%q_ifc, &
         nwpsaf%t, nwpsaf%q, nwpsaf%clc, nwpsaf%clw, nwpsaf%cli, nwpsaf%qnc, nwpsaf%qr, nwpsaf%qs, nwpsaf%dz, &
         nwpsaf%rho, nwpsaf%tv, nwpsaf%lwc, nwpsaf%iwc, nwpsaf%cdnc, nwpsaf%Reff, &
         nwpsaf%day, nwpsaf%month, nwpsaf%year, nwpsaf%point)

  end subroutine nwpsaf_clear

end module mod_nwpsaf

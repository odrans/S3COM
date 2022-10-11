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

MODULE MOD_ICON

  USE s3com_types,         ONLY: wp, type_icon
  USE mod_read_icon,       ONLY: icon_read, extract_coordinates

  USE s3com_config,        ONLY: rd, rv, epsilon, mu, nu, a, b, Q_ext, rholiq

  IMPLICIT NONE

CONTAINS

  SUBROUTINE icon_load(fname, icon)

    ! Inputs
    CHARACTER(LEN=256), INTENT(IN) :: fname

    ! Inputs/Outputs
    TYPE(type_icon), INTENT(INOUT) :: icon

    ! Internal
    INTEGER(KIND=4) :: nlevels, npoints

    write(6, "(1X, A, 1X, A)", advance='no') "Reading ICON inputs:", trim(fname)

    ! Extract the number of vertical levels and grid points in the input files
    CALL extract_coordinates(fname, nlevels, npoints)

    ! Initialize the ICON array
    CALL icon_init(npoints, nlevels, icon)

    ! Read input netcdf file containing ICON outputs
    CALL icon_read(fname, icon)

    ! Post-process ICON data
    CALL icon_process(icon)

    write(6, "(A)") "... Done"

  END SUBROUTINE icon_load


  SUBROUTINE icon_process(icon)

    TYPE(type_icon), INTENT(INOUT) :: icon

    ! Internal
    INTEGER(KIND=4) :: i, j, nlevels, npoints


    nlevels = icon%nlevels
    npoints = icon%npoints


    ! Atmospehric structure
    ! ----------------------------------------------------------------------------------------------------

    ! Layer depth
    icon%dz(:,1:nlevels) = abs(icon%zh(:,1:nlevels) - icon%z(:,1:nlevels)) * 2._wp

    ! ----------------------------------------------------------------------------------------------------

    ! Atmospehric moisture
    ! ----------------------------------------------------------------------------------------------------
    ! Virtual temperature
    icon%tv = icon%t * (1 + 0.608 * icon%q)

    ! Moist air density (kg/m3)
    icon%rho = icon%p / (rd * icon%tv)

    ! ! Saturation vapour pressure of liquid water (Pa); Murphy and Koop (2005)
    ! icon%es_w = EXP(54.842763 - 6763.22 / icon%t - 4.210 * LOG(icon%t) + 0.000367 * icon%t + TANH(0.0415 * (icon%t - 218.8)) * &
    !      (53.878 - 1331.22 / icon%t - 9.44523 * LOG(icon%t) + 0.014025 * icon%t) )

    ! ! Saturation vapour pressure of ice water (Pa); Murphy and Koop (2005)
    ! icon%es_i = EXP(9.550426 - 5723.265 / icon%t + 3.53068 * LOG(icon%t) - 0.00728332 * icon%t)

    ! ----------------------------------------------------------------------------------------------------


    ! Cloud properties
    ! ----------------------------------------------------------------------------------------------------
    ! Cloud liquid water content (kg/m3)
    icon%lwc = icon%clw * icon%rho

    ! Cloud ice water content (kg/m3)
    icon%iwc = icon%cli * icon%rho

    ! Cloud droplet number concentration (particules/m3)
    icon%cdnc = icon%qnc * icon%rho

    !Cloud liquid water effective radius (m)
    DO i = 1, Npoints
       DO j = 1, Nlevels
          IF(icon%cdnc(i,j) .GT. 0) THEN
             icon%Reff(i,j) = (a/2._wp)*(gamma((nu+1._wp+3._wp*b)/mu)/gamma((nu+1._wp+2._wp*b)/mu))*(icon%lwc(i,j) / icon%cdnc(i,j))**b*&
                  (gamma((nu+1._wp)/mu)/gamma((nu+2._wp)/mu))**b
          END IF
       END DO
    END DO
    icon%Reff = icon%Reff * 1E6 !! m to um (default input in RTTOV)

    !!Cloud extinction coefficient in m-1
    ! icon%beta_ext = (3._wp/4._wp)*(Q_ext/rholiq)*(icon%lwc/icon%Reff)

    ! !!Layer thickness
    ! icon%dz_cod(:,1:nlevels-1) = ((rd*icon%tv)/grav)*(log(icon%p(:,1:nlevels+1)/icon%p(:,1:nlevels)))

    ! !!Cloud optical depth
    ! icon%cod(:,1:nlevels-1) = ((icon%beta_ext(:,1:nlevels)+icon%beta_ext(:,1:nlevels+1))/2._wp)*icon%dz_cod(:,1:nlevels-1)
    ! ----------------------------------------------------------------------------------------------------



  END SUBROUTINE icon_process

  SUBROUTINE icon_init(npoints, nlevels, y)

    TYPE(type_icon) :: y
    INTEGER(KIND=4), INTENT(IN) :: npoints, nlevels

    y%nlevels = nlevels; y%npoints = npoints

    ALLOCATE(y%height(nlevels)); y%height = 0._wp
    ALLOCATE(y%lon(npoints)); y%lon = 0._wp
    ALLOCATE(y%lat(npoints)); y%lat = 0._wp
    ALLOCATE(y%lon_orig(npoints)); y%lon_orig = 0._wp
    ALLOCATE(y%lat_orig(npoints)); y%lat_orig = 0._wp
    ALLOCATE(y%orography(npoints)); y%orography = 0._wp
    ALLOCATE(y%landmask(npoints)); y%landmask = 0._wp
    ALLOCATE(y%psfc(npoints)); y%psfc = 0._wp
    ALLOCATE(y%skt(npoints)); y%skt = 0._wp
    ALLOCATE(y%t2m(npoints)); y%t2m = 0._wp
    ALLOCATE(y%q2m(npoints)); y%q2m = 0._wp
    ALLOCATE(y%u_wind(npoints)); y%u_wind = 0._wp
    ALLOCATE(y%v_wind(npoints)); y%v_wind = 0._wp
    ALLOCATE(y%p(npoints, nlevels)); y%p = 0._wp
    ALLOCATE(y%z(npoints, nlevels)); y%z = 0._wp
    ALLOCATE(y%zh(npoints, nlevels+1)); y%zh = 0._wp
    ALLOCATE(y%t(npoints, nlevels)); y%t = 0._wp
    ALLOCATE(y%q(npoints, nlevels)); y%q = 0._wp
    ALLOCATE(y%tca(npoints, nlevels)); y%tca = 0._wp
    ALLOCATE(y%clw(npoints, nlevels)); y%clw = 0._wp
    ALLOCATE(y%cli(npoints, nlevels)); y%cli = 0._wp
    ALLOCATE(y%qnc(npoints, nlevels)); y%qnc = 0._wp
    ALLOCATE(y%qr(npoints, nlevels)); y%qr = 0._wp
    ALLOCATE(y%qs(npoints, nlevels)); y%qs = 0._wp
    ALLOCATE(y%dz(npoints, nlevels)); y%dz = 0._wp
    ALLOCATE(y%es_w(npoints, nlevels)); y%es_w = 0._wp
    ALLOCATE(y%es_i(npoints, nlevels)); y%es_i = 0._wp
    ALLOCATE(y%rho(npoints, nlevels)); y%rho = 0._wp
    ALLOCATE(y%tv(npoints, nlevels)); y%tv = 0._wp
    ALLOCATE(y%lwc(npoints, nlevels)); y%lwc = 0._wp
    ALLOCATE(y%iwc(npoints, nlevels)); y%iwc = 0._wp
    ALLOCATE(y%cdnc(npoints, nlevels)); y%cdnc = 0._wp
    ALLOCATE(y%Reff(npoints, nlevels)); y%Reff = 0._wp
    ALLOCATE(y%beta_ext(npoints, nlevels)); y%beta_ext = 0._wp
    ALLOCATE(y%dz_cod(npoints, nlevels)); y%dz_cod = 0._wp
    ALLOCATE(y%cod(npoints, nlevels)); y%cod = 0._wp

  END SUBROUTINE icon_init


  SUBROUTINE icon_clear(y)

    TYPE(type_icon), INTENT(INOUT) :: y

    DEALLOCATE(y%height, y%lon, y%lat, y%lon_orig, y%lat_orig, &
         y%orography, y%landmask, y%psfc, y%skt, y%t2m, y%q2m, y%u_wind, y%v_wind, &
         y%p, y%z, y%zh, y%t, y%q, y%tca, y%clw, y%cli, y%qnc, y%qr, y%qs, y%dz, &
         y%es_w, y%es_i, y%rho, y%tv, &
         y%lwc, y%iwc, y%cdnc, y%Reff, y%beta_ext, y%dz_cod, y%cod)

  END SUBROUTINE icon_clear




END MODULE MOD_ICON

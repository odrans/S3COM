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

    ! Extract the number of vertical levels and grid points in the input files
    CALL extract_coordinates(fname, nlevels, npoints)

    ! Initialize the ICON array
    CALL icon_init(npoints, nlevels, icon)

    ! Read input netcdf file containing ICON outputs
    CALL icon_read(fname, icon)

    ! Post-process ICON data
    CALL icon_process(icon)

  END SUBROUTINE icon_load


  SUBROUTINE icon_process(icon)

    TYPE(type_icon), INTENT(INOUT) :: icon

    ! Internal
    INTEGER(KIND=4) :: i, j, nlevels, npoints

    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: &
         t_celcius, &
         e_sat, &
         e, &
         q_sat, &
         rho, &
         r, &
         pd, &
         rh, &
         tv

    nlevels = icon%nlevels
    npoints = icon%npoints

    ALLOCATE(&
         rh(npoints, nlevels), &
         e_sat(npoints, nlevels), &
         e(npoints, nlevels), &
         pd(npoints, nlevels), &
         r(npoints, nlevels), &
         q_sat(npoints, nlevels), &
         rho(npoints, nlevels), &
         t_celcius(npoints, nlevels), &
         tv(npoints, nlevels))


    ! Atmospehric structure
    ! ----------------------------------------------------------------------------------------------------

    ! Layer depth
    icon%dz(:,1:nlevels) = abs(icon%zh(:,1:nlevels) - icon%z(:,1:nlevels)) * 2._wp

    ! Air temperature in Celcius
    t_celcius = icon%t - 273.15_wp

    ! ----------------------------------------------------------------------------------------------------


    ! Atmospehric moisture
    ! ----------------------------------------------------------------------------------------------------
    ! Saturation vapour pressure of water (Pa); Monteith and Unsworth (2008) Tetens' formula for T >= 0 °C
    e_sat = 610.78_wp * exp((17.27_wp * t_celcius) / (t_celcius + 237.3_wp))

    ! Saturation specific humidity (kg/kg)
    q_sat = epsilon * ( e_sat / (icon%p - (1._wp - epsilon) * e_sat) )

    ! Water vapor mixing ratio (0-1)
    r = icon%sh / (1._wp - icon%sh)

    ! Relative humidity (%)
    rh = (r / q_sat) * 100._wp

    ! Pressure of water vapor (Pa)
    e = e_sat * ( rh / 100._wp )

    ! Pressure of dry air (Pa)
    pd = icon%p - e

    ! Moist air density (kg/m3)
    rho = ( pd / ( rd * icon%t ) )+( e / (rv * icon%t) )

    ! Virtual temperature
    tv = icon%t / (1._wp - (1._wp - epsilon) * (e / icon%p))
    ! ----------------------------------------------------------------------------------------------------


    ! Cloud properties
    ! ----------------------------------------------------------------------------------------------------
    ! Cloud liquid water content (kg/m3)
    icon%lwc = icon%clw * rho

    ! Cloud ice water content (kg/m3)
    icon%iwc = icon%cli * rho

    ! Cloud droplet number concentration (particules/m3)
    icon%cdnc = icon%qnc * rho

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
    icon%beta_ext = (3._wp/4._wp)*(Q_ext/rholiq)*(icon%lwc/icon%Reff)

    ! !!Layer thickness
    ! icon%dz_cod(:,1:nlevels-1) = ((rd*icon%tv)/grav)*(log(icon%p(:,1:nlevels+1)/icon%p(:,1:nlevels)))

    ! !!Cloud optical depth
    ! icon%cod(:,1:nlevels-1) = ((icon%beta_ext(:,1:nlevels)+icon%beta_ext(:,1:nlevels+1))/2._wp)*icon%dz_cod(:,1:nlevels-1)
    ! ----------------------------------------------------------------------------------------------------



    stop


  END SUBROUTINE icon_process

  SUBROUTINE icon_init(npoints, nlevels, y)

    TYPE(type_icon) :: y
    INTEGER(KIND=4), INTENT(IN) :: npoints, nlevels

    ALLOCATE(y%lon(Npoints)); y%lon = 0._wp
    ALLOCATE(y%lat(Npoints)); y%lat = 0._wp
    ALLOCATE(y%orography(Npoints)); y%orography = 0._wp
    ALLOCATE(y%landmask(Npoints)); y%landmask = 0._wp
    ALLOCATE(y%psfc(Npoints)); y%psfc = 0._wp
    ALLOCATE(y%skt(Npoints)); y%skt = 0._wp
    ALLOCATE(y%t2m(Npoints)); y%t2m = 0._wp
    ALLOCATE(y%q2m(Npoints)); y%q2m = 0._wp
    ALLOCATE(y%u_wind(Npoints)); y%u_wind = 0._wp
    ALLOCATE(y%v_wind(Npoints)); y%v_wind = 0._wp
    ALLOCATE(y%z(Npoints,Nlevels)); y%z = 0._wp
    ALLOCATE(y%zh(Npoints,Nlevels+1)); y%zh = 0._wp
    ALLOCATE(y%p(Npoints,Nlevels)); y%p = 0._wp
    ALLOCATE(y%t(Npoints,Nlevels)); y%t = 0._wp
    ALLOCATE(y%sh(Npoints,Nlevels)); y%sh = 0._wp
    ALLOCATE(y%tca(Npoints,Nlevels)); y%tca = 0._wp
    ALLOCATE(y%clw(Npoints,Nlevels)); y%clw = 0._wp
    ALLOCATE(y%cli(Npoints,Nlevels)); y%cli = 0._wp
    ALLOCATE(y%qnc(Npoints,Nlevels)); y%qnc = 0._wp
    ALLOCATE(y%iwc(Npoints,Nlevels)); y%iwc = 0._wp
    ALLOCATE(y%dz(Npoints,Nlevels)); y%dz = 0._wp
    ALLOCATE(y%t_celcius(Npoints,Nlevels)); y%t_celcius = 0._wp
    ALLOCATE(y%e_sat(Npoints,Nlevels)); y%e_sat = 0._wp
    ALLOCATE(y%ssh(Npoints,Nlevels)); y%ssh = 0._wp
    ALLOCATE(y%wv(Npoints,Nlevels)); y%wv = 0._wp
    ALLOCATE(y%rh(Npoints,Nlevels)); y%rh = 0._wp
    ALLOCATE(y%pv(Npoints,Nlevels)); y%pv = 0._wp
    ALLOCATE(y%pd(Npoints,Nlevels)); y%pd = 0._wp
    ALLOCATE(y%lwc(Npoints,Nlevels)); y%lwc = 0._wp
    ALLOCATE(y%cdnc(Npoints,Nlevels)); y%cdnc = 0._wp
    ALLOCATE(y%Reff(Npoints,Nlevels)); y%Reff = 0._wp
    ALLOCATE(y%Deff(Npoints,Nlevels)); y%Deff = 0._wp
    ALLOCATE(y%beta_ext(Npoints,Nlevels)); y%beta_ext = 0._wp
    ALLOCATE(y%tv(Npoints,Nlevels)); y%tv = 0._wp
    ALLOCATE(y%dz_cod(Npoints,Nlevels)); y%dz_cod = 0._wp
    ALLOCATE(y%cod(Npoints,Nlevels)); y%cod = 0._wp
    ALLOCATE(y%height(Nlevels)); y%height = 0._wp

    y%nlevels = Nlevels; y%npoints = Npoints

  END SUBROUTINE icon_init


  SUBROUTINE icon_clear(y)

    TYPE(type_icon), INTENT(INOUT) :: y

    DEALLOCATE(&
         y%lon, y%lat, y%orography, y%landmask, y%psfc, &
         y%skt, y%t2m, y%q2m, y%u_wind, y%v_wind, y%z, &
         y%zh, y%p, y%t, y%sh, y%tca, y%clw, y%cli,y%qnc, &
         y%rho_atm, y%iwc, y%dz, y%t_celcius, y%e_sat, y%ssh, &
         y%wv, y%rh, y%pv, y%pd, y%rho, y%lwc, y%cdnc, y%Reff, &
         y%Deff, y%beta_ext,  y%tv, y%dz_cod, y%cod, y%height)

  END SUBROUTINE icon_clear




END MODULE MOD_ICON

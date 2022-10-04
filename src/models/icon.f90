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

    nlevels = icon%nlevels
    npoints = icon%npoints

    !!Air density (kg/m3)
    icon%rho_atm = icon%p/(rd*icon%t)
    !!Cloud ice water content (kg/m3)
    icon%iwc = icon%cli*icon%rho_atm
    !!Layer depth
    icon%dz(:,1:nlevels) = abs(icon%zh(:,1:nlevels)-icon%z(:,1:nlevels))*2._wp

    !!Temperature (degrees Celsius)
    icon%t_celcius = icon%t-273.15_wp
    !!Saturation vapour pressure of water (Pa) - Monteith and Unsworth (2008) Tetens' formula for T >= 0 °C
    icon%e_sat = 610.78_wp*exp((17.27_wp*icon%t_celcius)/(icon%t_celcius+237.3_wp))
    !!Saturation specific humidity (kg/kg)
    icon%ssh = epsilon*(icon%e_sat/(icon%p-(1._wp-epsilon)*icon%e_sat))
    !!Water vapor mixing ratio (0-1)
    icon%wv = icon%sh/(1._wp-icon%sh)
    !!Relative humidity (%)
    icon%rh = (icon%wv/icon%ssh)*100.
    !!Pressure of water vapor (Pa)
    icon%pv = icon%e_sat*(icon%rh/100.)
    !!Pressure of dry air (Pa)
    icon%pd = icon%p-icon%pv
    !!Air density (kg/m3)
    icon%rho = (icon%pd/(rd*icon%t))+(icon%pv/(rv*icon%t))
    !!Cloud liquid water content (kg/m3)
    icon%lwc = icon%clw*icon%rho
    !!Cloud droplet number concentration (particules/m3)
    icon%cdnc = icon%rho * icon%qnc
    !!Cloud liquid water effective radius (m)
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
    !!Virtual temperature
    icon%tv = icon%t/(1._wp-(1._wp-epsilon)*(icon%pv/icon%p))
    !!Layer thickness
    icon%dz_cod(:,1:nlevels-1) = ((rd*icon%tv)/grav)*(log(icon%p(:,1:nlevels+1)/icon%p(:,1:nlevels)))
    !!Cloud optical depth
    icon%cod(:,1:nlevels-1) = ((icon%beta_ext(:,1:nlevels)+icon%beta_ext(:,1:nlevels+1))/2._wp)*icon%dz_cod(:,1:nlevels-1)


  END SUBROUTINE icon_process

  SUBROUTINE icon_init(npoints, nlevels, y)

    TYPE(type_icon) :: y
    INTEGER(KIND=4), INTENT(IN) :: npoints, nlevels

    ALLOCATE(y%lon(Npoints),         &
         y%lat(Npoints),               &
         y%orography(Npoints),         &
         y%landmask(Npoints),          &
         y%psfc(Npoints),              &
         y%skt(Npoints),               &
         y%t2m(Npoints),               &
         y%q2m(Npoints),               &
         y%u_wind(Npoints),            &
         y%v_wind(Npoints),            &
         y%z(Npoints,Nlevels),         &
         y%zh(Npoints,Nlevels+1),      &
         y%p(Npoints,Nlevels),         &
         y%t(Npoints,Nlevels),         &
         y%sh(Npoints,Nlevels),        &
         y%tca(Npoints,Nlevels),       &
         y%clw(Npoints,Nlevels),       &
         y%cli(Npoints,Nlevels),       &
         y%qnc(Npoints,Nlevels),       &
         y%rho_atm(Npoints,Nlevels),   &
         y%iwc(Npoints,Nlevels),       &
         y%dz(Npoints,Nlevels),        &
         y%t_celcius(Npoints,Nlevels), &
         y%e_sat(Npoints,Nlevels),     &
         y%ssh(Npoints,Nlevels),       &
         y%wv(Npoints,Nlevels),        &
         y%rh(Npoints,Nlevels),        &
         y%pv(Npoints,Nlevels),        &
         y%pd(Npoints,Nlevels),        &
         y%rho(Npoints,Nlevels),       &
         y%lwc(Npoints,Nlevels),       &
         y%cdnc(Npoints,Nlevels),      &
         y%Reff(Npoints,Nlevels),      &
         y%Deff(Npoints,Nlevels),      &
         y%beta_ext(Npoints,Nlevels),  &
         y%tv(Npoints,Nlevels),        &
         y%dz_cod(Npoints,Nlevels),    &
         y%cod(Npoints,Nlevels),      &
         y%height(Nlevels))

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

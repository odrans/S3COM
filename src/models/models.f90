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

MODULE MOD_MODELS

  USE s3com_types,         ONLY: wp, type_icon, type_model
  USE mod_icon,            ONLY: icon_load, icon_clear
  USE mod_utils_math,      ONLY: solar_angles

  IMPLICIT NONE

CONTAINS

  SUBROUTINE models_load(fname, model)

    ! Inputs
    CHARACTER(LEN=256), INTENT(IN) :: fname

    ! Output variables
    TYPE(type_model), INTENT(OUT) :: model

    ! Internal
    INTEGER(KIND=4) :: nlevels, npoints
    TYPE(type_icon) :: icon

    ! Load input data
    CALL icon_load(fname, icon)

    ! Coordinates
    npoints = icon%npoints
    nlevels = icon%nlevels

    ! Initialize model array
    CALL models_init(model, npoints, nlevels)

    ! Set up part of the model with ICON data
    CALL models_setup_icon(model, icon)

    ! Set up solar angles corresponding to model data
    CALL models_setup_solar(model)

    ! Clear input data
    CALL icon_clear(icon)

    write(*,*) "Model input loaded"

  END SUBROUTINE models_load


  SUBROUTINE models_init(model, npoints, nlevels)

    ! Input variables
    INTEGER(KIND = 4), INTENT(IN) :: npoints, nlevels

    ! Output variables
    TYPE(type_model), INTENT(OUT)   :: model

    ! Internal variables

    model%nPoints   =  npoints
    model%nLevels   =  nlevels

    !! 2D fields
    ALLOCATE(model%lat(npoints)); model%lat = 0._wp
    ALLOCATE(model%lon(npoints)); model%lon = 0._wp
    ALLOCATE(model%lat_orig(npoints)); model%lat_orig = 0._wp
    ALLOCATE(model%lon_orig(npoints)); model%lon_orig = 0._wp

    ALLOCATE(model%orography(npoints)); model%orography = 0._wp
    ALLOCATE(model%u_wind(npoints)); model%u_wind = 0._wp
    ALLOCATE(model%v_wind(npoints)); model%v_wind = 0._wp
    ALLOCATE(model%skt(npoints)); model%skt = 0._wp
    ALLOCATE(model%psfc(npoints)); model%psfc = 0._wp
    ALLOCATE(model%q2m(npoints)); model%q2m = 0._wp
    ALLOCATE(model%t2m(npoints)); model%t2m = 0._wp
    ALLOCATE(model%landmask(npoints)); model%landmask = 0._wp
    ALLOCATE(model%sunzenangle(npoints)); model%sunzenangle = 0._wp
    ALLOCATE(model%sunazangle(npoints)); model%sunazangle = 0._wp

    !! 3D fields
    ALLOCATE(model%co2(npoints, nlevels)); model%co2 = 0._wp
    ALLOCATE(model%ch4(npoints, nlevels)); model%ch4 = 0._wp
    ALLOCATE(model%n2o(npoints, nlevels)); model%n2o = 0._wp
    ALLOCATE(model%s2o(npoints, nlevels)); model%s2o = 0._wp
    ALLOCATE(model%co(npoints, nlevels)); model%co = 0._wp
    ALLOCATE(model%p(npoints, nlevels)); model%p = 0._wp
    ALLOCATE(model%z(npoints, nlevels)); model%z = 0._wp
    ALLOCATE(model%dz(npoints, nlevels)); model%dz = 0._wp
    ALLOCATE(model%t(npoints, nlevels)); model%t = 0._wp
    ALLOCATE(model%sh(npoints, nlevels)); model%sh = 0._wp
    ALLOCATE(model%tca(npoints, nlevels)); model%tca = 0._wp
    ALLOCATE(model%reff(npoints, nlevels)); model%reff = 0._wp
    ALLOCATE(model%cdnc(npoints, nlevels)); model%cdnc = 0._wp
    ALLOCATE(model%iwc(npoints, nlevels)); model%iwc = 0._wp
    ALLOCATE(model%lwc(npoints, nlevels)); model%lwc = 0._wp


  END SUBROUTINE models_init


  SUBROUTINE models_setup_icon(model, icon)

    ! Input variables
    TYPE(type_icon), INTENT(IN)    :: icon

    ! Output variables
    TYPE(type_model), INTENT(INOUT)   :: model

    model%Nlat = icon%nlat
    model%Nlon = icon%nlon
    model%mode = icon%mode

    model%nPoints   =  icon%nPoints
    model%nLevels   =  icon%nLevels

    model%lat       =  icon%lat
    model%lon       =  icon%lon
    model%orography =  icon%orography
    model%u_wind    =  icon%u_wind
    model%v_wind    =  icon%v_wind
    model%skt       =  icon%skt
    model%psfc      =  icon%psfc
    model%q2m       =  icon%q2m
    model%t2m       =  icon%t2m
    model%landmask  =  icon%landmask

    model%p         =  icon%p
    model%z         =  icon%z
    model%dz        =  icon%dz
    model%t         =  icon%t
    model%sh        =  icon%q
    model%tca       =  icon%tca
    model%iwc       =  icon%iwc
    model%lwc       =  icon%lwc
    model%reff      =  icon%reff
    model%cdnc      =  icon%cdnc

  END SUBROUTINE models_setup_icon

  SUBROUTINE models_setup_solar(model)

    TYPE(type_model), INTENT(INOUT)   :: model

    INTEGER(KIND = 4), DIMENSION(3) :: date, time
    INTEGER(KIND = 4) :: i

    !! This needs to be replaced by real date
    date = (/02, 05, 2013/)
    time = (/12, 00, 00/)

    DO i = 1, model%npoints

       CALL solar_angles(model%lat(i), model%lon(i), date, time, &
            model%sunzenangle(i), model%sunazangle(i))

    END DO

  END SUBROUTINE models_setup_solar

END MODULE MOD_MODELS

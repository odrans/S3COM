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

module mod_models

  use s3com_types,         only: wp, type_icon, type_model, type_nml
  use mod_icon,            only: icon_load, icon_clear
  use mod_utils_math,      only: solar_angles

  implicit none

  private
  public :: models_load, models_init, models_deinit, models_setup_icon, models_setup_solar

contains

  subroutine models_load(nml, model)

    ! Inputs
    type(type_nml), intent(IN)     :: nml

    ! Output variables
    type(type_model), intent(OUT) :: model

    ! Internal
    integer(kind=4) :: nlayers, npoints
    type(type_icon) :: icon

    ! Load input data
    call icon_load(nml%fname_in, icon)

    ! Coordinates
    npoints = icon%npoints
    nlayers = icon%nlayers

    ! Initialize model array
    call models_init(model, npoints, nlayers)

    ! Set up part of the model with ICON data
    call models_setup_icon(model, icon)

    ! Set up solar angles corresponding to model data
    call models_setup_solar(model)

    ! Clear input data
    call icon_clear(icon)

    write(*,*) "Model input loaded"

  end subroutine models_load


  subroutine models_init(model, npoints, nlayers)

    ! Input variables
    integer(kind = 4), intent(IN) :: npoints, nlayers

    ! Output variables
    type(type_model), intent(OUT)   :: model

    ! Internal variables
    integer(kind = 4) :: nlevels

    nlevels = nlayers + 1

    model%npoints   =  npoints
    model%nlayers   =  nlayers
    model%nlevels   =  nlevels

    model%date = (/0, 0, 0/)
    model%time = (/0, 0, 0/)

    !! 2D fields
    allocate(model%lat(npoints), source = 0._wp)
    allocate(model%lon, model%lat_orig, model%lon_orig, &
         model%topography, model%u_10m, model%v_10m, &
         model%ts, model%ps, model%q_2m, model%t_2m, &
         model%landmask, model%sunzenangle, model%sunazangle, &
         mold = model%lat)

    !! 3D fields at atmospheric levels
    allocate(model%p(npoints, nlevels), source = 0._wp)
    allocate(model%t, model%q, model%co2, model%ch4, &
         model%n2o, model%s2o, model%co, &
         mold = model%p)

    !! 3D fields in atmospheric layers
    allocate(model%z(npoints, nlayers), source = 0._wp)
    allocate(model%dz, model%clc, model%reff, model%cdnc, &
         model%iwc, model%lwc, &
         mold = model%z)

  end subroutine models_init

  subroutine models_deinit(model)

    ! Output variables
    type(type_model), intent(inout)   :: model

    deallocate(model%lat, model%lon, model%lat_orig, model%lon_orig, &
         model%topography, model%u_10m, model%v_10m, &
         model%ts, model%ps, model%q_2m, model%t_2m, &
         model%landmask, model%sunzenangle, model%sunazangle, &
         model%p, model%t, model%q, model%co2, model%ch4, &
         model%n2o, model%s2o, model%co, &
         model%z, model%dz, model%clc, model%reff, model%cdnc, &
         model%iwc, model%lwc)

  end subroutine models_deinit

  subroutine models_setup_icon(model, icon)

    ! Input variables
    type(type_icon), intent(in)    :: icon

    ! Output variables
    type(type_model), intent(inout)   :: model

    real(wp) :: hour, minute, sec
    character(len=8) :: icon_date

    ! Extract information from the ICON time (%Y%m%d.%f)
    ! ----------------------------------------------------------------------------------------------------
    hour = (icon%time - int(icon%time)) * 24; model%time(1) = int(hour)
    minute = (hour - int(hour)) * 60; model%time(2) = int(minute)
    sec = (minute - int(minute)) * 60; model%time(3) = int(sec)

    write(icon_date, "(I8)") int(icon%time)
    read(icon_date(1:4), "(I)") model%date(3)
    read(icon_date(5:6), "(I)") model%date(2)
    read(icon_date(7:8), "(I)") model%date(1)
    ! ----------------------------------------------------------------------------------------------------

    model%nlat = icon%nlat
    model%nlon = icon%nlon
    model%mode = icon%mode

    model%npoints   =  icon%npoints
    model%nlevels   =  icon%nlevels
    model%nlayers   =  icon%nlayers

    model%lat       =  icon%lat
    model%lon       =  icon%lon
    model%topography =  icon%topography
    model%u_10m     =  icon%u_10m
    model%v_10m     =  icon%v_10m
    model%ts        =  icon%ts
    model%ps        =  icon%ps
    model%q_2m      =  icon%q_2m
    model%t_2m      =  icon%t_2m
    model%landmask  =  icon%landmask

    model%p         =  icon%p_ifc
    model%z         =  icon%z_ifc
    model%dz        =  icon%dz
    model%t         =  icon%t_ifc
    model%q         =  icon%q_ifc
    model%clc       =  icon%clc
    model%iwc       =  icon%iwc
    model%lwc       =  icon%lwc
    model%reff      =  icon%reff
    model%cdnc      =  icon%cdnc

  end subroutine models_setup_icon

  subroutine models_setup_solar(model)

    type(type_model), intent(inout)   :: model

    integer(kind = 4) :: i

    do i = 1, model%npoints

       call solar_angles(model%lat(i), model%lon(i), model%date, model%time, &
            model%sunzenangle(i), model%sunazangle(i))

    end do

  end subroutine models_setup_solar

end module mod_models

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

MODULE s3com_types

  IMPLICIT NONE

  integer, parameter :: sp = selected_real_kind(6, 37)
  integer, parameter :: dp = selected_real_kind(12, 307)
  integer, parameter :: wp = sp

  TYPE type_s3com
     REAL(kind=wp), DIMENSION(:,:) :: &
          OE_Y          ! Measurement vector from the optimal estimation method
     INTEGER(kind = 4) :: &
          nchannels
  END TYPE type_s3com

  type type_icon
     integer(kind=4) :: &
          nlevels, npoints, nhydro, nlat, nlon, mode! Dimensions
     real(wp) :: &
          maemsfc_lw

     real(wp),dimension(:),allocatable:: &
          lon_orig,  & ! Longitude that won't be regridded (deg)
          lat_orig,  & ! Latitude  that won't be regridded (deg)
          lon,       & ! Longitude (deg)
          lat,       & ! Latitude (deg)
          skt,       & ! Skin temperature (K)
          psfc,      & ! Surface pressure (Pa)
          landmask,  & ! Land/sea mask (0/1)
          seaice,    & ! sea ice mask (0/1)
          orography, & ! surface height (I think?)
          u_wind,    & ! U-component of wind (m/s)
          v_wind,    & ! V-component of wind (m/s)
          sunlit       ! Sunlit flag
     real(wp),dimension(:,:),allocatable :: &
          p,         & ! Model pressure levels (pa)
          p_lev,     & ! Model pressure levels (pa)          
          ph,        & ! Moddel pressure @ half levels (pa)
          zlev,      & ! Model level height (m)
          zlev_half, & ! Model level height @ half-levels (m)
          T,         & ! Temperature (K)
          T_lev,     & ! Temperature (K)          
          sh,        & ! Specific humidity (kg/kg)
          sh_lev,    & ! Specific humidity (kg/kg)          
          rh,        & ! Relative humidity (1)
          tca,       & ! Total cloud fraction (1)
          cca,       & ! Convective cloud fraction (1)
          dtau_s,    & ! 0.67micron optical depth (stratiform cloud) (1)
          dtau_c,    & ! 0.67micron optical depth (convective cloud) (1)
          dem_s,     & ! 11micron emissivity (stratiform cloud) 
          dem_c,     & ! 11microm emissivity (convective cloud)
          lwc,       & ! Liquid water content (kg/m3)
          iwc          ! Ice water content (kg/m3)

     real(wp),dimension(:,:,:),allocatable :: &
          Reff         ! Subcolumn effective radius
  end type type_icon

  type type_rttov_atm
     integer, pointer :: &
          npoints,      & ! Number of profiles to simulate
          nlevels,      & ! Number of levels
          nlayers         ! Number of layers
     real(wp), dimension(:), pointer :: &
          lat,          & ! Latitude
          lon,          & ! Longitude
          t_skin,       & ! Surface skin temperature
          h_surf,       & ! Surface height
          p_surf,       & ! Surface pressure
          u_surf,       & ! U component of surface wind
          v_surf,       & ! V component of surface wind
          t2m,          & ! 2-m Temperature
          q2m,          & ! 2-m Specific humidity
          lsmask        & ! land-sea mask
     real(wp), dimension(:,:), pointer :: &
          p,            & ! Pressure @ model levels
          z,            & ! Height @ model levels
          t,            & ! Temperature
          q               ! Specific humidity
     real(wp), dimension(:,:), pointer :: &
          tca,          & ! Cloud fraction
          iwc,          & ! ice water content
          lwc,          & ! liquid water content
          reff,         & ! droplet effective radius
          cdnc          & ! cloud droplet number concentration
  end type type_rttov_atm

  type type_rttov_opt
     
     integer ::     &
          dosolar,          &
          nchannels,        &
          platform,         &
          satellite,        &
          instrument,       &
          month                ! Month (needed for surface emissivity calculation)
     integer,dimension(:),allocatable :: &
          channel_list
     real(wp) :: &
          zenangle,azangle,sunzenangle,sunazangle
  end type type_rttov_opt
  

  
END MODULE s3com_types

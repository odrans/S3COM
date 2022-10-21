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

module MOD_WRITE_OUTPUT

  use s3com_types,  only: type_icon, wp, type_nml, type_model, type_s3com
  use netcdf
  use mod_read_icon, only: map_point_to_ll

  implicit none

contains

  ! Write all S3COM outputs
  subroutine write_output(s3com, model, nml)

    ! Input variables
    type(type_model), intent(IN) :: model
    type(type_s3com), intent(IN) :: s3com
    type(type_nml), intent(IN) :: nml

    call write_output_rad(s3com, model, nml)

    ! IF(nml%flag_output_atm) THEN
    !    CALL write_output_atm(icon, oe, nml, atm)
    ! END IF

    ! IF(nml%flag_retrievals) THEN
    !    CALL write_output_ret(icon, oe, nml)
    ! END IF

  end subroutine write_output


  ! Write radiation outputs, mainly satellite measurements simulated by RTTOV
  subroutine write_output_rad(s3com, model, nml)

    ! Input variables
    type(type_model),       intent(IN) :: model
    type(type_s3com),      intent(IN) :: s3com
    type(type_nml),        intent(IN) :: nml

    ! Local variables
    real(KIND=wp), dimension(model%nlon, model%nlat, nml%nchannels) :: &
         gridded_f_ref_total,  &
         gridded_f_ref_clear,  &
         gridded_f_bt_total,   &
         gridded_f_bt_clear,   &
         gridded_f_rad_total,  &
         gridded_f_rad_clear

    integer(KIND=4) :: ncid, errst

    integer(KIND=4) ::      &
         varid_lon,         &
         varid_lat,         &
         varid_chan,        &
         varid_ref_total,   &
         varid_ref_clear,   &
         varid_bt_total,    &
         varid_bt_clear,    &
         varid_rad_total,   &
         varid_rad_clear

    integer(KIND=4) ::      &
         dimid_lon,         &
         dimid_lat,         &
         dimid_chan,        &
         dimid_latlon(2),   &
         dimid_latlonchan(3)

    character(LEN = 256) :: fn_out_rad, suffix

    suffix = trim(nml%suffix_out)
    if(trim(suffix) .ne. "") suffix = "_"//trim(suffix)//"_"

    fn_out_rad = trim(nml%path_out)//"S3COM"//trim(suffix)//"_rad.nc"

    call map_point_to_ll(model%Nlon, model%Nlat, model%mode, x2=s3com%rad%f_ref_total, y3=gridded_f_ref_total)
    call map_point_to_ll(model%Nlon, model%Nlat, model%mode, x2=s3com%rad%f_ref_clear, y3=gridded_f_ref_clear)
    call map_point_to_ll(model%Nlon, model%Nlat, model%mode, x2=s3com%rad%f_bt_total,  y3=gridded_f_bt_total)
    call map_point_to_ll(model%Nlon, model%Nlat, model%mode, x2=s3com%rad%f_bt_clear,  y3=gridded_f_bt_clear)
    call map_point_to_ll(model%Nlon, model%Nlat, model%mode, x2=s3com%rad%f_rad_total, y3=gridded_f_rad_total)
    call map_point_to_ll(model%Nlon, model%Nlat, model%mode, x2=s3com%rad%f_rad_clear, y3=gridded_f_rad_clear)

    errst = nf90_create(fn_out_rad, NF90_CLOBBER, ncid)

    errst = nf90_def_dim(ncid, "lon", model%nlon, dimid_lon)
    errst = nf90_def_dim(ncid, "lat", model%nlat, dimid_lat)
    errst = nf90_def_dim(ncid, "chan", nml%nchannels, dimid_chan)

    dimid_latlon     = (/dimid_lon, dimid_lat/)
    dimid_latlonchan = (/dimid_lon, dimid_lat, dimid_chan/)

    errst = nf90_def_var(ncid, "lon",             NF90_REAL, dimid_lon,        varid_lon)
    errst = nf90_def_var(ncid, "lat",             NF90_REAL, dimid_lat,        varid_lat)
    errst = nf90_def_var(ncid, "chan",            NF90_REAL, dimid_chan,       varid_chan)
    errst = nf90_def_var(ncid, "ref_total",       NF90_REAL, dimid_latlonchan, varid_ref_total)
    errst = nf90_def_var(ncid, "ref_clear",       NF90_REAL, dimid_latlonchan, varid_ref_clear)
    errst = nf90_def_var(ncid, "BT_total",        NF90_REAL, dimid_latlonchan, varid_bt_total)
    errst = nf90_def_var(ncid, "BT_clear",        NF90_REAL, dimid_latlonchan, varid_bt_clear)
    errst = nf90_def_var(ncid, "rad_total",       NF90_REAL, dimid_latlonchan, varid_rad_total)
    errst = nf90_def_var(ncid, "rad_clear",       NF90_REAL, dimid_latlonchan, varid_rad_clear)

    errst = nf90_put_att(ncid, varid_lon,        "units", "degrees_east")
    errst = nf90_put_att(ncid, varid_lat,        "units", "degrees_north")
    errst = nf90_put_att(ncid, varid_chan,       "units", "um")
    errst = nf90_put_att(ncid, varid_ref_total,  "units", "-")
    errst = nf90_put_att(ncid, varid_ref_clear,  "units", "-")
    errst = nf90_put_att(ncid, varid_bt_total,   "units", "K")
    errst = nf90_put_att(ncid, varid_bt_clear,   "units", "K")
    errst = nf90_put_att(ncid, varid_rad_total,  "units", "W/m2/sr/um")
    errst = nf90_put_att(ncid, varid_rad_clear,  "units", "W/m2/sr/um")

    errst = nf90_put_att(ncid, varid_lon,        "description", "Longitude")
    errst = nf90_put_att(ncid, varid_lat,        "description", "Latitude")
    errst = nf90_put_att(ncid, varid_chan,       "description", "Central wavelength of instrument channel")
    errst = nf90_put_att(ncid, varid_ref_total,  "description", "All-sky upwelling reflectance at TOA")
    errst = nf90_put_att(ncid, varid_ref_clear,  "description", "Clear-sky upwelling reflectance at TOA")
    errst = nf90_put_att(ncid, varid_bt_total,   "description", "All-sky upwelling brightness temperature at TOA")
    errst = nf90_put_att(ncid, varid_bt_clear,   "description", "Clear-sky upwelling brightness temperature at TOA")
    errst = nf90_put_att(ncid, varid_rad_total,  "description", "All-sky upwelling radiance at TOA")
    errst = nf90_put_att(ncid, varid_rad_clear,  "description", "Clear-sky upwelling radiance at TOA")

    errst = nf90_enddef(ncid)

    errst = nf90_put_var(ncid, varid_lon,        model%lon_orig)
    errst = nf90_put_var(ncid, varid_lat,        model%lat_orig)
    errst = nf90_put_var(ncid, varid_chan,       s3com%rad%wavelength)
    errst = nf90_put_var(ncid, varid_ref_total,  gridded_f_ref_total)
    errst = nf90_put_var(ncid, varid_ref_clear,  gridded_f_ref_clear)
    errst = nf90_put_var(ncid, varid_bt_total,   gridded_f_bt_total)
    errst = nf90_put_var(ncid, varid_bt_clear,   gridded_f_bt_clear)
    errst = nf90_put_var(ncid, varid_rad_total,  gridded_f_rad_total)
    errst = nf90_put_var(ncid, varid_rad_clear,  gridded_f_rad_clear)

    errst = nf90_close(ncid)

  end subroutine write_output_rad

  ! Write atmospheric outputs
  subroutine write_output_atm(icon, oe, nml, atm)

    ! Input variables
    type(type_icon),       intent(IN) :: icon
    type(type_s3com),      intent(IN) :: oe
    type(type_nml),        intent(IN) :: nml
    type(type_s3com),      intent(IN) :: atm

    real(KIND=wp), dimension(icon%Nlon, icon%Nlat, icon%Nlevels) :: &
         gridded_atm_t, &
         gridded_atm_z, &
         gridded_atm_clc, &
         gridded_atm_cdnc, &
         gridded_atm_reff, &
         gridded_atm_lwc

    integer(KIND=4) :: ncid, errst
    integer(KIND=4) :: &
         varid_lon, &
         varid_lat, &
         varid_lev, &
         varid_atm_t, &
         varid_atm_z, &
         varid_atm_clc, &
         varid_atm_cdnc, &
         varid_atm_reff, &
         varid_atm_lwc

    integer(KIND=4) :: dimid_lon, dimid_lat, dimid_lev, dimid_latlonlev(3)

    character(LEN = 256) :: fn_out_atm, suffix

    suffix = trim(nml%suffix_out)
    if(trim(suffix) .ne. "") suffix = "_"//trim(suffix)//"_"

    fn_out_atm = trim(nml%path_out)//"S3COM"//trim(suffix)//"_atm.nc"

    call map_point_to_ll(icon%Nlon, icon%Nlat, icon%mode, x2=atm%atm%t,     y3=gridded_atm_t)
    call map_point_to_ll(icon%Nlon, icon%Nlat, icon%mode, x2=atm%atm%z,     y3=gridded_atm_z)
    call map_point_to_ll(icon%Nlon, icon%Nlat, icon%mode, x2=atm%atm%clc,   y3=gridded_atm_clc)
    call map_point_to_ll(icon%Nlon, icon%Nlat, icon%mode, x2=atm%atm%cdnc,  y3=gridded_atm_cdnc)
    call map_point_to_ll(icon%Nlon, icon%Nlat, icon%mode, x2=atm%atm%reff,  y3=gridded_atm_reff)
    call map_point_to_ll(icon%Nlon, icon%Nlat, icon%mode, x2=atm%atm%lwc,   y3=gridded_atm_lwc)

    errst = nf90_create(fn_out_atm, NF90_CLOBBER, ncid)

    errst = nf90_def_dim(ncid, "Longitude", icon%Nlon,       dimid_lon)
    errst = nf90_def_dim(ncid, "Latitude",  icon%Nlat,       dimid_lat)
    errst = nf90_def_dim(ncid, "height",    icon%Nlevels,    dimid_lev)

    dimid_latlonlev = (/dimid_lon, dimid_lat, dimid_lev/)

    errst = nf90_def_var(ncid, "Longitude",       NF90_REAL, dimid_lon,        varid_lon)
    errst = nf90_def_var(ncid, "Latitude",        NF90_REAL, dimid_lat,        varid_lat)
    errst = nf90_def_var(ncid, "height",          NF90_REAL, dimid_lev,        varid_lev)
    errst = nf90_def_var(ncid, "ta",              NF90_REAL, dimid_latlonlev,  varid_atm_t)
    errst = nf90_def_var(ncid, "z",               NF90_REAL, dimid_latlonlev,  varid_atm_z)
    errst = nf90_def_var(ncid, "clc",             NF90_REAL, dimid_latlonlev,  varid_atm_clc)
    errst = nf90_def_var(ncid, "cdnc",            NF90_REAL, dimid_latlonlev,  varid_atm_cdnc)
    errst = nf90_def_var(ncid, "lwc",             NF90_REAL, dimid_latlonlev,  varid_atm_lwc)
    errst = nf90_def_var(ncid, "reff",            NF90_REAL, dimid_latlonlev,  varid_atm_reff)

    errst = nf90_def_var_fill(ncid, varid_atm_reff, 0, 0)
    errst = nf90_def_var_fill(ncid, varid_atm_cdnc, 0, 0)
    errst = nf90_def_var_fill(ncid, varid_atm_lwc, 0, 0)

    errst = nf90_put_att(ncid, varid_lon,        "units", "degrees_east")
    errst = nf90_put_att(ncid, varid_lat,        "units", "degrees_north")
    errst = nf90_put_att(ncid, varid_lev,        "units", "")
    errst = nf90_put_att(ncid, varid_atm_z,      "units", "m")
    errst = nf90_put_att(ncid, varid_atm_t,      "units", "K")
    errst = nf90_put_att(ncid, varid_atm_clc,    "units", "K")
    errst = nf90_put_att(ncid, varid_atm_cdnc,   "units", "m-3")
    errst = nf90_put_att(ncid, varid_atm_lwc,    "units", "kg m-3")
    errst = nf90_put_att(ncid, varid_atm_reff,   "units", "um")

    errst = nf90_enddef(ncid)

    errst = nf90_put_var(ncid, varid_lon,        icon%lon_orig)
    errst = nf90_put_var(ncid, varid_lat,        icon%lat_orig)
    errst = nf90_put_var(ncid, varid_lev,        icon%height)
    errst = nf90_put_var(ncid, varid_atm_t,      gridded_atm_t)
    errst = nf90_put_var(ncid, varid_atm_z,      gridded_atm_z)
    errst = nf90_put_var(ncid, varid_atm_clc,    gridded_atm_clc)
    errst = nf90_put_var(ncid, varid_atm_cdnc,   gridded_atm_cdnc)
    errst = nf90_put_var(ncid, varid_atm_lwc ,   gridded_atm_lwc)
    errst = nf90_put_var(ncid, varid_atm_reff,   gridded_atm_reff)

    errst = nf90_close(ncid)

  end subroutine write_output_atm


! Write retrieval outputs
  ! subroutine write_output_ret(icon, oe, nml)

  !   ! Input variables
  !   type(type_icon),       intent(IN) :: icon
  !   type(type_s3com),      intent(IN) :: oe
  !   type(type_nml),        intent(IN) :: nml

  !   real(KIND=wp), dimension(icon%Nlon, icon%Nlat) :: &
  !        gridded_iwp_model, &
  !        gridded_iwp_ret, &
  !        gridded_g

  !   integer(KIND=4) :: ncid, errst
  !   integer(KIND=4) :: &
  !        varid_lon, &
  !        varid_lat, &
  !        varid_iwp_ret, &
  !        varid_iwp_mod, &
  !        varid_g

  !   integer(KIND=4) :: dimid_lon, dimid_lat, dimid_chan, dimid_latlon(2)

  !   character(LEN = 256) :: fn_out_ret, suffix

  !   suffix = trim(nml%suffix_out)
  !   if(trim(suffix) .ne. "") suffix = "_"//trim(suffix)//"_"

  !   fn_out_ret = trim(nml%path_out)//"S3COM"//trim(suffix)//"_ret.nc"

  !   call map_point_to_ll(icon%Nlon, icon%Nlat, icon%mode, x1=oe%Xip1(:,1),    y2=gridded_iwp_ret)
  !   call map_point_to_ll(icon%Nlon, icon%Nlat, icon%mode, x1=oe%iwp_model(:), y2=gridded_iwp_model)
  !   call map_point_to_ll(icon%Nlon, icon%Nlat, icon%mode, x1=oe%gip1(:),      y2=gridded_g)

  !   errst = nf90_create(fn_out_ret, NF90_CLOBBER, ncid)

  !   errst = nf90_def_dim(ncid, "Longitude", icon%Nlon,       dimid_lon)
  !   errst = nf90_def_dim(ncid, "Latitude",  icon%Nlat,       dimid_lat)

  !   dimid_latlon     = (/dimid_lon, dimid_lat/)

  !   errst = nf90_def_var(ncid, "Longitude",       NF90_REAL, dimid_lon,        varid_lon)
  !   errst = nf90_def_var(ncid, "Latitude",        NF90_REAL, dimid_lat,        varid_lat)
  !   errst = nf90_def_var(ncid, "iwp_ret",         NF90_REAL, dimid_latlon,     varid_iwp_ret)
  !   errst = nf90_def_var(ncid, "iwp_model",       NF90_REAL, dimid_latlon,     varid_iwp_mod)

  !   errst = nf90_put_att(ncid, varid_lon,        "units", "degrees_east")
  !   errst = nf90_put_att(ncid, varid_lat,        "units", "degrees_north")
  !   errst = nf90_put_att(ncid, varid_iwp_ret,    "units", "kg/m2")
  !   errst = nf90_put_att(ncid, varid_iwp_mod,    "units", "kg/m2")

  !   errst = nf90_enddef(ncid)

  !   errst = nf90_put_var(ncid, varid_lon,        icon%lon_orig)
  !   errst = nf90_put_var(ncid, varid_lat,        icon%lat_orig)
  !   errst = nf90_put_var(ncid, varid_iwp_ret,    gridded_iwp_ret)
  !   errst = nf90_put_var(ncid, varid_iwp_mod,    gridded_iwp_model)

  !   errst = nf90_close(ncid)

  ! end subroutine write_output_ret




end module MOD_WRITE_OUTPUT

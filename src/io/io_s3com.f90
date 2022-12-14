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

module mod_io_s3com
   
   use s3com_types,  only: type_icon, wp, type_nml, type_model, type_s3com
   use mod_io_utils, only: map_point_to_ll
   use netcdf
   
   implicit none
   
   private
   public :: write_output, write_output_rad, write_output_atm
   
contains
   
   !!Write all S3COM outputs
   subroutine write_output(s3com, model, nml)
      
      !!Input variables
      type(type_s3com), intent(in) :: s3com
      type(type_model), intent(in) :: model
      type(type_nml),   intent(in) :: nml
      
      call write_output_rad(s3com, model, nml)
      
      if (s3com%nml%flag_output_jac) then
         call write_output_jac(s3com, model, nml)
         write(*,*) "outputting jacobian results"
      endif
      
      !if (nml%flag_output_atm) then
      !   call write_output_atm(icon, oe, nml, atm)
      !endif
      
      !if (nml%flag_retrievals) then
      !   call write_output_ret(icon, oe, nml)
      !endif
      
   end subroutine write_output
   
   !!Write radiation outputs, mainly satellite measurements simulated by RTTOV
   subroutine write_output_rad(s3com, model, nml)
      
      !!Input variables
      type(type_s3com), intent(in) :: s3com
      type(type_model), intent(in) :: model
      type(type_nml),   intent(in) :: nml
      
      !!Local variables
      real(kind=wp), dimension(model%nlon, model%nlat, nml%nchannels) :: &
         gridded_f_ref_total, &
         gridded_f_ref_clear, &
         gridded_f_bt_total,  &
         gridded_f_bt_clear,  &
         gridded_f_rad_total, &
         gridded_f_rad_clear
         
      integer(kind=4) :: ncid, errst
      
      integer(kind=4) ::  &
         varid_lon,       &
         varid_lat,       &
         varid_chan,      &
         varid_ref_total, &
         varid_ref_clear, &
         varid_bt_total,  &
         varid_bt_clear,  &
         varid_rad_total, &
         varid_rad_clear
         
      integer(kind=4) ::  &
         dimid_lon,       &
         dimid_lat,       &
         dimid_chan,      &
         dimid_latlon(2), &
         dimid_latlonchan(3)
         
      character(LEN = 256) :: fn_out_rad, suffix, attr_instrument
      
      suffix = trim(nml%suffix_out)
      if(trim(suffix) .ne. "") suffix = "_"//trim(suffix)
      
      fn_out_rad = trim(nml%path_out)//"S3COM"//trim(suffix)//"_rad.nc"
      
      call map_point_to_ll(model%nlon, model%nlat, model%mode, x2=s3com%rad%f_ref_total, y3=gridded_f_ref_total)
      call map_point_to_ll(model%nlon, model%nlat, model%mode, x2=s3com%rad%f_ref_clear, y3=gridded_f_ref_clear)
      call map_point_to_ll(model%nlon, model%nlat, model%mode, x2=s3com%rad%f_bt_total,  y3=gridded_f_bt_total)
      call map_point_to_ll(model%nlon, model%nlat, model%mode, x2=s3com%rad%f_bt_clear,  y3=gridded_f_bt_clear)
      call map_point_to_ll(model%nlon, model%nlat, model%mode, x2=s3com%rad%f_rad_total, y3=gridded_f_rad_total)
      call map_point_to_ll(model%nlon, model%nlat, model%mode, x2=s3com%rad%f_rad_clear, y3=gridded_f_rad_clear)
      
      errst = nf90_create(fn_out_rad, NF90_CLOBBER, ncid)
      
      errst = nf90_def_dim(ncid, "lon", model%nlon, dimid_lon)
      errst = nf90_def_dim(ncid, "lat", model%nlat, dimid_lat)
      errst = nf90_def_dim(ncid, "chan", nml%nchannels, dimid_chan)
      
      dimid_latlon     = (/dimid_lon, dimid_lat/)
      dimid_latlonchan = (/dimid_lon, dimid_lat, dimid_chan/)
      
      errst = nf90_def_var(ncid, "lon",       NF90_REAL, dimid_lon,        varid_lon)
      errst = nf90_def_var(ncid, "lat",       NF90_REAL, dimid_lat,        varid_lat)
      errst = nf90_def_var(ncid, "chan",      NF90_REAL, dimid_chan,       varid_chan)
      errst = nf90_def_var(ncid, "ref_total", NF90_REAL, dimid_latlonchan, varid_ref_total)
      errst = nf90_def_var(ncid, "ref_clear", NF90_REAL, dimid_latlonchan, varid_ref_clear)
      errst = nf90_def_var(ncid, "bt_total",  NF90_REAL, dimid_latlonchan, varid_bt_total)
      errst = nf90_def_var(ncid, "bt_clear",  NF90_REAL, dimid_latlonchan, varid_bt_clear)
      errst = nf90_def_var(ncid, "rad_total", NF90_REAL, dimid_latlonchan, varid_rad_total)
      errst = nf90_def_var(ncid, "rad_clear", NF90_REAL, dimid_latlonchan, varid_rad_clear)
      
      errst = nf90_put_att(ncid, varid_lon,       "standard_name", "longitude")
      errst = nf90_put_att(ncid, varid_lat,       "standard_name", "latitude")
      errst = nf90_put_att(ncid, varid_chan,      "standard_name", "chan")
      errst = nf90_put_att(ncid, varid_ref_total, "standard_name", "ref_total")
      errst = nf90_put_att(ncid, varid_ref_clear, "standard_name", "ref_clear")
      errst = nf90_put_att(ncid, varid_bt_total,  "standard_name", "bt_total")
      errst = nf90_put_att(ncid, varid_bt_clear,  "standard_name", "bt_clear")
      errst = nf90_put_att(ncid, varid_rad_total, "standard_name", "rad_total")
      errst = nf90_put_att(ncid, varid_rad_clear, "standard_name", "rad_clear")
      
      errst = nf90_put_att(ncid, varid_lon,       "long_name", "Longitude")
      errst = nf90_put_att(ncid, varid_lat,       "long_name", "Latitude")
      errst = nf90_put_att(ncid, varid_chan,      "long_name", "Central wavelength of instrument channel")
      errst = nf90_put_att(ncid, varid_ref_total, "long_name", "All-sky upwelling reflectance at TOA")
      errst = nf90_put_att(ncid, varid_ref_clear, "long_name", "Clear-sky upwelling reflectance at TOA")
      errst = nf90_put_att(ncid, varid_bt_total,  "long_name", "All-sky upwelling brightness temperature at TOA")
      errst = nf90_put_att(ncid, varid_bt_clear,  "long_name", "Clear-sky upwelling brightness temperature at TOA")
      errst = nf90_put_att(ncid, varid_rad_total, "long_name", "All-sky upwelling radiance at TOA")
      errst = nf90_put_att(ncid, varid_rad_clear, "long_name", "Clear-sky upwelling radiance at TOA")
      
      errst = nf90_put_att(ncid, varid_lon,       "units", "degrees_east")
      errst = nf90_put_att(ncid, varid_lat,       "units", "degrees_north")
      errst = nf90_put_att(ncid, varid_chan,      "units", "um")
      errst = nf90_put_att(ncid, varid_ref_total, "units", "-")
      errst = nf90_put_att(ncid, varid_ref_clear, "units", "-")
      errst = nf90_put_att(ncid, varid_bt_total,  "units", "K")
      errst = nf90_put_att(ncid, varid_bt_clear,  "units", "K")
      errst = nf90_put_att(ncid, varid_rad_total, "units", "W/m2/sr/um")
      errst = nf90_put_att(ncid, varid_rad_clear, "units", "W/m2/sr/um")
      
      attr_instrument = trim(s3com%opt%rttov%inst_name)//"/"//trim(s3com%opt%rttov%platform_name)
      
      errst = nf90_put_att(ncid, nf90_global, 'RTTOV_instrument', trim(attr_instrument))
      
      errst = nf90_enddef(ncid)
      
      errst = nf90_put_var(ncid, varid_lon,       model%lon_orig)
      errst = nf90_put_var(ncid, varid_lat,       model%lat_orig)
      errst = nf90_put_var(ncid, varid_chan,      s3com%rad%wavelength)
      errst = nf90_put_var(ncid, varid_ref_total, gridded_f_ref_total)
      errst = nf90_put_var(ncid, varid_ref_clear, gridded_f_ref_clear)
      errst = nf90_put_var(ncid, varid_bt_total,  gridded_f_bt_total)
      errst = nf90_put_var(ncid, varid_bt_clear,  gridded_f_bt_clear)
      errst = nf90_put_var(ncid, varid_rad_total, gridded_f_rad_total)
      errst = nf90_put_var(ncid, varid_rad_clear, gridded_f_rad_clear)
      
      errst = nf90_close(ncid)
      
   end subroutine write_output_rad
   
   !!Write Jacobian outputs calculated by RTTOV
   subroutine write_output_jac(s3com, model, nml)
      
      !!Input variables
      type(type_s3com), intent(in) :: s3com
      type(type_model), intent(in) :: model
      type(type_nml),   intent(in) :: nml
      
      !!Local variables
      real(kind=wp), dimension(model%nlon, model%nlat, model%nlevels, nml%nchannels) :: &
         gridded_jac_p, &
         gridded_jac_t
      
      real(kind=wp), dimension(model%nlon, model%nlat, model%nlayers, nml%nchannels) :: &
         gridded_jac_clc, &
         gridded_jac_deff
         
      integer(kind=4) :: ncid, errst
      
      integer(kind=4) :: &
         varid_lon,      &
         varid_lat,      &
         varid_lev,      &
         varid_lev_2,    &
         varid_chan,     &
         varid_jac_p,    &
         varid_jac_t,    &
         varid_jac_clc,  &
         varid_jac_deff
         
      integer(kind=4) ::      &
         dimid_lon,           &
         dimid_lat,           &
         dimid_lev,           &
         dimid_lev_2,         &
         dimid_chan,          &
         dimid_latlon(2),     &
         dimid_latlonchan(3), &
         dimid_latlonlevchan(4), &
         dimid_latlonlev2chan(4)
         
      character(LEN = 256) :: fn_out_jac, suffix, attr_instrument
      
      suffix = trim(nml%suffix_out)
      if(trim(suffix) .ne. "") suffix = "_"//trim(suffix)
      
      fn_out_jac = trim(nml%path_out)//"S3COM"//trim(suffix)//"_jac.nc"
      
      call map_point_to_ll(model%nlon, model%nlat, model%mode, x3=s3com%jac%p, y4=gridded_jac_p)
      write(*,*) "gridded_jac_p OK"
      call map_point_to_ll(model%nlon, model%nlat, model%mode, x3=s3com%jac%t, y4=gridded_jac_t)
      write(*,*) "gridded_jac_t OK"
      call map_point_to_ll(model%nlon, model%nlat, model%mode, x3=s3com%jac%cfrac, y4=gridded_jac_clc)
      write(*,*) "gridded_jac_cfrac OK"
      call map_point_to_ll(model%nlon, model%nlat, model%mode, x3=s3com%jac%clwde, y4=gridded_jac_deff)
      write(*,*) "gridded_jac_deff OK"
      
      errst = nf90_create(fn_out_jac, NF90_CLOBBER, ncid)
      
      errst = nf90_def_dim(ncid, "lon",   model%nlon,    dimid_lon)
      errst = nf90_def_dim(ncid, "lat",   model%nlat,    dimid_lat)
      errst = nf90_def_dim(ncid, "lev",   model%nlayers, dimid_lev)
      errst = nf90_def_dim(ncid, "lev_2", model%nlevels, dimid_lev_2)
      errst = nf90_def_dim(ncid, "chan",  nml%nchannels, dimid_chan)
      
      !dimid_latlon        = (/dimid_lon, dimid_lat/)
      !dimid_latlonchan    = (/dimid_lon, dimid_lat, dimid_chan/)
      dimid_latlonlevchan  = (/dimid_lon, dimid_lat, dimid_lev, dimid_chan/)
      dimid_latlonlev2chan = (/dimid_lon, dimid_lat, dimid_lev_2, dimid_chan/)
      
      errst = nf90_def_var(ncid, "lon",      NF90_REAL, dimid_lon,            varid_lon)
      errst = nf90_def_var(ncid, "lat",      NF90_REAL, dimid_lat,            varid_lat)
      errst = nf90_def_var(ncid, "chan",     NF90_REAL, dimid_chan,           varid_chan)
      errst = nf90_def_var(ncid, "height",   NF90_REAL, dimid_lev,            varid_lev)
      errst = nf90_def_var(ncid, "height_2", NF90_REAL, dimid_lev_2,          varid_lev_2)
      errst = nf90_def_var(ncid, "jac_p",    NF90_REAL, dimid_latlonlev2chan, varid_jac_p)
      errst = nf90_def_var(ncid, "jac_t",    NF90_REAL, dimid_latlonlev2chan, varid_jac_t)
      errst = nf90_def_var(ncid, "jac_clc",  NF90_REAL, dimid_latlonlevchan,  varid_jac_clc)
      errst = nf90_def_var(ncid, "jac_reff", NF90_REAL, dimid_latlonlevchan,  varid_jac_deff)
      
      errst = nf90_put_att(ncid, varid_lon,      "standard_name", "longitude")
      errst = nf90_put_att(ncid, varid_lat,      "standard_name", "latitude")
      errst = nf90_put_att(ncid, varid_lev,      "standard_name", "altitude")
      errst = nf90_put_att(ncid, varid_lev_2,    "standard_name", "altitude_2")
      errst = nf90_put_att(ncid, varid_chan,     "standard_name", "chan")
      errst = nf90_put_att(ncid, varid_jac_p,    "standard_name", "p_jac")
      errst = nf90_put_att(ncid, varid_jac_t,    "standard_name", "T_jac")
      errst = nf90_put_att(ncid, varid_jac_clc,  "standard_name", "clc_jac")
      errst = nf90_put_att(ncid, varid_jac_deff, "standard_name", "deff_jac")
      
      errst = nf90_put_att(ncid, varid_lon,      "long_name", "Longitude")
      errst = nf90_put_att(ncid, varid_lat,      "long_name", "Latitude")
      errst = nf90_put_att(ncid, varid_lev,      "long_name", "Altitude")
      errst = nf90_put_att(ncid, varid_lev_2,    "long_name", "Altitude_2")
      errst = nf90_put_att(ncid, varid_chan,     "long_name", "Central wavelength of instrument channel")
      errst = nf90_put_att(ncid, varid_jac_p,    "long_name", "Pressure Jacobian")
      errst = nf90_put_att(ncid, varid_jac_t,    "long_name", "Temperature Jacobian")
      errst = nf90_put_att(ncid, varid_jac_clc,  "long_name", "Cloud cover Jacobian")
      errst = nf90_put_att(ncid, varid_jac_deff, "long_name", "Cloud droplet effective diameter Jacobian")
      
      errst = nf90_put_att(ncid, varid_lon,      "units", "degrees_east")
      errst = nf90_put_att(ncid, varid_lat,      "units", "degrees_north")
      errst = nf90_put_att(ncid, varid_lev,      "units", "")
      errst = nf90_put_att(ncid, varid_lev_2,    "units", "")
      errst = nf90_put_att(ncid, varid_chan,     "units", "um")
      errst = nf90_put_att(ncid, varid_jac_p,    "units", "(W/m2/sr/um)/hPa")
      errst = nf90_put_att(ncid, varid_jac_t,    "units", "(W/m2/sr/um)/K")
      errst = nf90_put_att(ncid, varid_jac_clc,  "units", "(W/m2/sr/um)")
      errst = nf90_put_att(ncid, varid_jac_deff, "units", "(W/m2/sr/um)/um")
      
      attr_instrument = trim(s3com%opt%rttov%inst_name)//"/"//trim(s3com%opt%rttov%platform_name)
      
      errst = nf90_put_att(ncid, nf90_global, 'RTTOV_instrument', trim(attr_instrument))
      
      errst = nf90_enddef(ncid)
      
      errst = nf90_put_var(ncid, varid_lon,      model%lon_orig)
      errst = nf90_put_var(ncid, varid_lat,      model%lat_orig)
      errst = nf90_put_var(ncid, varid_lev,      model%height)
      errst = nf90_put_var(ncid, varid_lev_2,    model%height_2)
      errst = nf90_put_var(ncid, varid_chan,     s3com%rad%wavelength)
      errst = nf90_put_var(ncid, varid_jac_p,    gridded_jac_p)
      errst = nf90_put_var(ncid, varid_jac_t,    gridded_jac_t)
      errst = nf90_put_var(ncid, varid_jac_clc,  gridded_jac_clc)
      errst = nf90_put_var(ncid, varid_jac_deff, gridded_jac_deff)
      
      errst = nf90_close(ncid)
      
   end subroutine write_output_jac
   
  ! Write atmospheric outputs
  subroutine write_output_atm(icon, oe, nml, atm)

    ! Input variables
    type(type_icon),       intent(in) :: icon
    type(type_s3com),      intent(in) :: oe
    type(type_nml),        intent(in) :: nml
    type(type_s3com),      intent(in) :: atm

    real(kind=wp), dimension(icon%Nlon, icon%Nlat, icon%Nlevels) :: &
         gridded_atm_t, &
         gridded_atm_z, &
         gridded_atm_clc, &
         gridded_atm_cdnc, &
         gridded_atm_reff, &
         gridded_atm_lwc

    integer(kind=4) :: ncid, errst
    integer(kind=4) :: &
         varid_lon, &
         varid_lat, &
         varid_lev, &
         varid_atm_t, &
         varid_atm_z, &
         varid_atm_clc, &
         varid_atm_cdnc, &
         varid_atm_reff, &
         varid_atm_lwc

    integer(kind=4) :: dimid_lon, dimid_lat, dimid_lev, dimid_latlonlev(3)

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
  !   type(type_icon),       intent(in) :: icon
  !   type(type_s3com),      intent(in) :: oe
  !   type(type_nml),        intent(in) :: nml

  !   real(kind=wp), dimension(icon%Nlon, icon%Nlat) :: &
  !        gridded_iwp_model, &
  !        gridded_iwp_ret, &
  !        gridded_g

  !   integer(kind=4) :: ncid, errst
  !   integer(kind=4) :: &
  !        varid_lon, &
  !        varid_lat, &
  !        varid_iwp_ret, &
  !        varid_iwp_mod, &
  !        varid_g

  !   integer(kind=4) :: dimid_lon, dimid_lat, dimid_chan, dimid_latlon(2)

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

end module mod_io_s3com

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
   
   use s3com_types,  only: wp, wpi, type_s3com, type_model, type_nml
   use mod_io_utils, only: map_point_to_ll
   use netcdf
   
   implicit none
   
   private
   public :: write_output, write_output_rad, write_output_k, write_output_k_tl, write_output_atm, write_output_ret
   
contains
   
   ! ============================================================================================================================
   !> @brief Write all S3COM outputs
   !! @param[in] s3com   s3com structure
   !! @param[in] model   model structure
   !! @param[in] nml     namelist structure
   subroutine write_output(s3com, model, nml)
      
      ! Input
      type(type_s3com), intent(in) :: s3com
      type(type_model), intent(in) :: model
      type(type_nml),   intent(in) :: nml
      
      ! Write radiation data outputs (by default)
      call write_output_rad(s3com, model, nml)
      
      ! Optional: write jacobian data outputs from the K model
      if (s3com%nml%flag_output_jac) then
         call write_output_k(s3com, model, nml)
      endif
      
      ! Optional: write jacobian data outputs from the TL model
      if (s3com%nml%flag_output_k_tl) then
         call write_output_k_tl(s3com, model, nml)
      endif
      
      ! Optional: write atmospheric data outputs
      if (s3com%nml%flag_output_atm) then
         call write_output_atm(s3com, model, nml)
      endif
      
      ! Optional: write retrieval data outputs
      if (s3com%nml%flag_retrievals) then
         call write_output_ret(s3com, model, nml)
      endif
      
   end subroutine write_output
   ! ============================================================================================================================
   
   ! ============================================================================================================================
   !> @brief Write radiation outputs, mainly satellite measurements simulated by RTTOV
   !! @details This subroutine outputs radiation data simulated by RTTOV at the top of atmosphere (TOA)
   !! @param[in] s3com   s3com structure
   !! @param[in] model   model structure
   !! @param[in] nml     namelist structure
   subroutine write_output_rad(s3com, model, nml)
      
      ! Input
      type(type_s3com), intent(in) :: s3com
      type(type_model), intent(in) :: model
      type(type_nml), intent(in) :: nml
      
      ! Internal
      real(wp), dimension(model%nlon, model%nlat, nml%nchannels) :: &
         gridded_f_ref_total, &
         gridded_f_ref_clear, &
         gridded_f_bt_total,  &
         gridded_f_bt_clear,  &
         gridded_f_rad_total, &
         gridded_f_rad_clear
      
      integer(wpi) :: ncid, errst
      
      integer(wpi) ::     &
         varid_pnt,       &
         varid_lon,       &
         varid_lat,       &
         varid_chan,      &
         varid_ref_total, &
         varid_ref_clear, &
         varid_bt_total,  &
         varid_bt_clear,  &
         varid_rad_total, &
         varid_rad_clear
      
      integer(wpi) ::       &
         dimid_lon,         &
         dimid_lat,         &
         dimid_chan,        &
         dimid_pnt,         &
         dimid_latlon(2),   &
         dimid_pntchan(2),  &
         dimid_latlonchan(3)
      
      character(len=256) :: fn_out_rad, suffix, attr_instrument
      
      suffix = trim(nml%suffix_out)
      if (trim(suffix) .ne. "") suffix = "_"//trim(suffix)
      
      fn_out_rad = trim(nml%path_out)//"S3COM"//trim(suffix)//"_rad.nc"
      
      if (model%mode > 1) then
         call map_point_to_ll(model%nlon, model%nlat, model%mode, x2=s3com%rad%f_ref_total, y3=gridded_f_ref_total)
         call map_point_to_ll(model%nlon, model%nlat, model%mode, x2=s3com%rad%f_ref_clear, y3=gridded_f_ref_clear)
         call map_point_to_ll(model%nlon, model%nlat, model%mode, x2=s3com%rad%f_bt_total,  y3=gridded_f_bt_total)
         call map_point_to_ll(model%nlon, model%nlat, model%mode, x2=s3com%rad%f_bt_clear,  y3=gridded_f_bt_clear)
         call map_point_to_ll(model%nlon, model%nlat, model%mode, x2=s3com%rad%f_rad_total, y3=gridded_f_rad_total)
         call map_point_to_ll(model%nlon, model%nlat, model%mode, x2=s3com%rad%f_rad_clear, y3=gridded_f_rad_clear)
      endif
      
      errst = nf90_create(fn_out_rad, NF90_CLOBBER, ncid)
      errst = nf90_def_dim(ncid, "chan", nml%nchannels, dimid_chan)
      
      if (model%mode > 1) then
         errst = nf90_def_dim(ncid, "lon", model%nlon, dimid_lon)
         errst = nf90_def_dim(ncid, "lat", model%nlat, dimid_lat)

         dimid_latlon     = (/dimid_lon, dimid_lat/)
         dimid_latlonchan = (/dimid_lon, dimid_lat, dimid_chan/)
      else
         errst = nf90_def_dim(ncid, "point", model%npoints, dimid_pnt)
         dimid_pntchan = (/dimid_pnt, dimid_chan/)
      endif
      
      if (model%mode > 1) then
         errst = nf90_def_var(ncid, "lon",       NF90_REAL, dimid_lon,        varid_lon)
         errst = nf90_def_var(ncid, "lat",       NF90_REAL, dimid_lat,        varid_lat)
         errst = nf90_def_var(ncid, "chan",      NF90_REAL, dimid_chan,       varid_chan)
         errst = nf90_def_var(ncid, "ref_total", NF90_REAL, dimid_latlonchan, varid_ref_total)
         errst = nf90_def_var(ncid, "ref_clear", NF90_REAL, dimid_latlonchan, varid_ref_clear)
         errst = nf90_def_var(ncid, "bt_total",  NF90_REAL, dimid_latlonchan, varid_bt_total)
         errst = nf90_def_var(ncid, "bt_clear",  NF90_REAL, dimid_latlonchan, varid_bt_clear)
         errst = nf90_def_var(ncid, "rad_total", NF90_REAL, dimid_latlonchan, varid_rad_total)
         errst = nf90_def_var(ncid, "rad_clear", NF90_REAL, dimid_latlonchan, varid_rad_clear)
      else
         errst = nf90_def_var(ncid, "point",     NF90_INT,  dimid_pnt,     varid_pnt)
         errst = nf90_def_var(ncid, "lon",       NF90_REAL, dimid_lon,     varid_lon)
         errst = nf90_def_var(ncid, "lat",       NF90_REAL, dimid_lat,     varid_lat)
         errst = nf90_def_var(ncid, "lon",       NF90_REAL, dimid_pnt,     varid_lon)
         errst = nf90_def_var(ncid, "lat",       NF90_REAL, dimid_pnt,     varid_lat)
         errst = nf90_def_var(ncid, "chan",      NF90_REAL, dimid_chan,    varid_chan)
         errst = nf90_def_var(ncid, "ref_total", NF90_REAL, dimid_pntchan, varid_ref_total)
         errst = nf90_def_var(ncid, "ref_clear", NF90_REAL, dimid_pntchan, varid_ref_clear)
         errst = nf90_def_var(ncid, "bt_total",  NF90_REAL, dimid_pntchan, varid_bt_total)
         errst = nf90_def_var(ncid, "bt_clear",  NF90_REAL, dimid_pntchan, varid_bt_clear)
         errst = nf90_def_var(ncid, "rad_total", NF90_REAL, dimid_pntchan, varid_rad_total)
         errst = nf90_def_var(ncid, "rad_clear", NF90_REAL, dimid_pntchan, varid_rad_clear)
      endif
      
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

      errst = nf90_put_var(ncid, varid_lon,  model%lon_orig)
      errst = nf90_put_var(ncid, varid_lat,  model%lat_orig)
      errst = nf90_put_var(ncid, varid_chan, s3com%rad%wavelength)

      if (model%mode > 1) then
         errst = nf90_put_var(ncid, varid_ref_total, gridded_f_ref_total)
         errst = nf90_put_var(ncid, varid_ref_clear, gridded_f_ref_clear)
         errst = nf90_put_var(ncid, varid_bt_total,  gridded_f_bt_total)
         errst = nf90_put_var(ncid, varid_bt_clear,  gridded_f_bt_clear)
         errst = nf90_put_var(ncid, varid_rad_total, gridded_f_rad_total)
         errst = nf90_put_var(ncid, varid_rad_clear, gridded_f_rad_clear)
      else
         errst = nf90_put_var(ncid, varid_pnt, model%point_orig)
         errst = nf90_put_var(ncid, varid_ref_total, s3com%rad%f_ref_total)
         errst = nf90_put_var(ncid, varid_ref_clear, s3com%rad%f_ref_clear)
         errst = nf90_put_var(ncid, varid_bt_total,  s3com%rad%f_bt_total)
         errst = nf90_put_var(ncid, varid_bt_clear,  s3com%rad%f_bt_clear)
         errst = nf90_put_var(ncid, varid_rad_total, s3com%rad%f_rad_total)
         errst = nf90_put_var(ncid, varid_rad_clear, s3com%rad%f_rad_clear)
      endif
      
      errst = nf90_close(ncid)
      
   end subroutine write_output_rad
   ! ============================================================================================================================
   
   ! ============================================================================================================================
   !> @brief Write Jacobian outputs calculated by RTTOV using the K model
   !! @details This subroutine outputs Jacobian data calculated by the RTTOV K model
   !! @param[in] s3com   s3com structure
   !! @param[in] model   model structure
   !! @param[in] nml     namelist structure
   subroutine write_output_k(s3com, model, nml)
      
      ! Input
      type(type_s3com), intent(in) :: s3com
      type(type_model), intent(in) :: model
      type(type_nml), intent(in) :: nml
      
      ! Internal
      real(wp), dimension(model%nlon, model%nlat, nml%nchannels, model%nlevels) :: &
         gridded_k_t,                                                              &
         gridded_k_q
      
      real(wp), dimension(model%nlon, model%nlat, nml%nchannels, model%nlayers) :: &
         gridded_k_clc,                                                            &
         gridded_k_deff
      
      integer(wpi) :: ncid, errst
      
      integer(wpi) :: &
         varid_lon,   &
         varid_lat,   &
         varid_lay,   &
         varid_lev,   &
         varid_chan,  &
         varid_k_t,   &
         varid_k_q,   &
         varid_k_clc, &
         varid_k_deff
      
      integer(wpi) ::            &
         dimid_lon,              &
         dimid_lat,              &
         dimid_lay,              &
         dimid_lev,              &
         dimid_chan,             &
         dimid_latlon(2),        &
         dimid_latlonchan(3),    &
         dimid_latlonlaychan(4), &
         dimid_latlonlevchan(4)
      
      character(len=256) :: fn_out_k, suffix, attr_instrument
      
      suffix = trim(nml%suffix_out)
      if (trim(suffix) .ne. "") suffix = "_"//trim(suffix)
      
      fn_out_k = trim(nml%path_out)//"S3COM"//trim(suffix)//"_k.nc"
      
      call map_point_to_ll(model%nlon, model%nlat, model%mode, x3=s3com%jac%t,     y4=gridded_k_t)
      call map_point_to_ll(model%nlon, model%nlat, model%mode, x3=s3com%jac%q,     y4=gridded_k_q)
      call map_point_to_ll(model%nlon, model%nlat, model%mode, x3=s3com%jac%cfrac, y4=gridded_k_clc)
      call map_point_to_ll(model%nlon, model%nlat, model%mode, x3=s3com%jac%clwde, y4=gridded_k_deff)
      
      errst = nf90_create(fn_out_k, NF90_CLOBBER, ncid)
      
      errst = nf90_def_dim(ncid, "lon",   model%nlon,    dimid_lon)
      errst = nf90_def_dim(ncid, "lat",   model%nlat,    dimid_lat)
      errst = nf90_def_dim(ncid, "layer", model%nlayers, dimid_lay)
      errst = nf90_def_dim(ncid, "level", model%nlevels, dimid_lev)
      errst = nf90_def_dim(ncid, "chan",  nml%nchannels, dimid_chan)
      
      dimid_latlon        = (/dimid_lon, dimid_lat/)
      dimid_latlonchan    = (/dimid_lon, dimid_lat, dimid_chan/)
      dimid_latlonlaychan = (/dimid_lon, dimid_lat, dimid_chan, dimid_lay/)
      dimid_latlonlevchan = (/dimid_lon, dimid_lat, dimid_chan, dimid_lev/)
      
      errst = nf90_def_var(ncid, "lon",    NF90_REAL, dimid_lon,           varid_lon)
      errst = nf90_def_var(ncid, "lat",    NF90_REAL, dimid_lat,           varid_lat)
      errst = nf90_def_var(ncid, "chan",   NF90_REAL, dimid_chan,          varid_chan)
      errst = nf90_def_var(ncid, "lay",    NF90_REAL, dimid_lay,           varid_lay)
      errst = nf90_def_var(ncid, "lev",    NF90_REAL, dimid_lev,           varid_lev)
      errst = nf90_def_var(ncid, "k_t",    NF90_REAL, dimid_latlonlevchan, varid_k_t)
      errst = nf90_def_var(ncid, "k_q",    NF90_REAL, dimid_latlonlevchan, varid_k_q)
      errst = nf90_def_var(ncid, "k_clc",  NF90_REAL, dimid_latlonlaychan, varid_k_clc)
      errst = nf90_def_var(ncid, "k_deff", NF90_REAL, dimid_latlonlaychan, varid_k_deff)
      
      errst = nf90_put_att(ncid, varid_lon,    "standard_name", "longitude")
      errst = nf90_put_att(ncid, varid_lat,    "standard_name", "latitude")
      errst = nf90_put_att(ncid, varid_lay,    "standard_name", "layer")
      errst = nf90_put_att(ncid, varid_lev,    "standard_name", "level")
      errst = nf90_put_att(ncid, varid_chan,   "standard_name", "chan")
      errst = nf90_put_att(ncid, varid_k_t,    "standard_name", "k_t")
      errst = nf90_put_att(ncid, varid_k_q,    "standard_name", "k_q")
      errst = nf90_put_att(ncid, varid_k_clc,  "standard_name", "k_clc")
      errst = nf90_put_att(ncid, varid_k_deff, "standard_name", "k_deff")
      
      errst = nf90_put_att(ncid, varid_lon,    "long_name", "Longitude")
      errst = nf90_put_att(ncid, varid_lat,    "long_name", "Latitude")
      errst = nf90_put_att(ncid, varid_lay,    "long_name", "Layer index")
      errst = nf90_put_att(ncid, varid_lev,    "long_name", "Level index")
      errst = nf90_put_att(ncid, varid_chan,   "long_name", "Central wavelength of instrument channel")
      errst = nf90_put_att(ncid, varid_k_t,    "long_name", "Temperature Jacobian")
      errst = nf90_put_att(ncid, varid_k_q,    "long_name", "Humidity Jacobian")
      errst = nf90_put_att(ncid, varid_k_clc,  "long_name", "Cloud cover Jacobian")
      errst = nf90_put_att(ncid, varid_k_deff, "long_name", "Cloud droplet effective diameter Jacobian")
      
      errst = nf90_put_att(ncid, varid_lon,    "units", "degrees_east")
      errst = nf90_put_att(ncid, varid_lat,    "units", "degrees_north")
      errst = nf90_put_att(ncid, varid_lay,    "units", "")
      errst = nf90_put_att(ncid, varid_lev,    "units", "")
      errst = nf90_put_att(ncid, varid_chan,   "units", "um")
      errst = nf90_put_att(ncid, varid_k_t,    "units", "(W/m2/sr/um)/K")
      errst = nf90_put_att(ncid, varid_k_q,    "units", "(W/m2/sr/um)/(kg/kg)")
      errst = nf90_put_att(ncid, varid_k_clc,  "units", "(W/m2/sr/um)")
      errst = nf90_put_att(ncid, varid_k_deff, "units", "(W/m2/sr/um)/um")
      
      attr_instrument = trim(s3com%opt%rttov%inst_name)//"/"//trim(s3com%opt%rttov%platform_name)
      
      errst = nf90_put_att(ncid, nf90_global, 'RTTOV_instrument', trim(attr_instrument))
      
      errst = nf90_enddef(ncid)
      
      errst = nf90_put_var(ncid, varid_lon,    model%lon_orig)
      errst = nf90_put_var(ncid, varid_lat,    model%lat_orig)
      errst = nf90_put_var(ncid, varid_lay,    model%height)
      errst = nf90_put_var(ncid, varid_lev,    model%height_2)
      errst = nf90_put_var(ncid, varid_chan,   s3com%rad%wavelength)
      errst = nf90_put_var(ncid, varid_k_t,    gridded_k_t)
      errst = nf90_put_var(ncid, varid_k_q,    gridded_k_q)
      errst = nf90_put_var(ncid, varid_k_clc,  gridded_k_clc)
      errst = nf90_put_var(ncid, varid_k_deff, gridded_k_deff)
      
      errst = nf90_close(ncid)
      
   end subroutine write_output_k
   ! ============================================================================================================================
   
   ! ============================================================================================================================
   !> @brief Write jacobian outputs calculated by RTTOV using the TL model
   !! @details This subroutine outputs Jacobian data calculated by the RTTOV TL model
   !! @param[in] s3com   s3com structure
   !! @param[in] model   model structure
   !! @param[in] nml     namelist structure
   subroutine write_output_k_tl(s3com, model, nml)
      
      ! Input
      type(type_s3com), intent(in) :: s3com
      type(type_model), intent(in) :: model
      type(type_nml), intent(in) :: nml
      
      ! Internal
      real(wp), dimension(model%nlon, model%nlat, nml%nchannels, model%nlevels) :: &
         gridded_k_tl_t
      
      integer(wpi) :: ncid, errst
      
      integer(wpi) ::     &
         varid_k_tl_lon,  &
         varid_k_tl_lat,  &
         varid_k_tl_lev,  &
         varid_k_tl_chan, &
         varid_k_tl_t
      
      integer(wpi) ::         &
         dimid_lon,           &
         dimid_lat,           &
         dimid_lev,           &
         dimid_chan,          &
         dimid_latlon(2),     &
         dimid_latlonchan(3), &
         dimid_latlonlevchan(4)
      
      character(len=256) :: fn_out_k_tl, suffix, attr_instrument
      
      suffix = trim(nml%suffix_out)
      if (trim(suffix) .ne. "") suffix = "_"//trim(suffix)
      
      fn_out_k_tl = trim(nml%path_out)//"S3COM"//trim(suffix)//"_k_tl.nc"
      
      call map_point_to_ll(model%nlon, model%nlat, model%mode, x3=s3com%k_tl%t, y4=gridded_k_tl_t)
      
      errst = nf90_create(fn_out_k_tl, NF90_CLOBBER, ncid)
      
      errst = nf90_def_dim(ncid, "lon",   model%nlon,    dimid_lon)
      errst = nf90_def_dim(ncid, "lat",   model%nlat,    dimid_lat)
      errst = nf90_def_dim(ncid, "level", model%nlevels, dimid_lev)
      errst = nf90_def_dim(ncid, "chan",  nml%nchannels, dimid_chan)
      
      dimid_latlon        = (/dimid_lon, dimid_lat/)
      dimid_latlonchan    = (/dimid_lon, dimid_lat, dimid_chan/)
      dimid_latlonlevchan = (/dimid_lon, dimid_lat, dimid_chan, dimid_lev/)
      
      errst = nf90_def_var(ncid, "lon",    NF90_REAL, dimid_lon,           varid_k_tl_lon)
      errst = nf90_def_var(ncid, "lat",    NF90_REAL, dimid_lat,           varid_k_tl_lat)
      errst = nf90_def_var(ncid, "lev",    NF90_REAL, dimid_lev,           varid_k_tl_lev)
      errst = nf90_def_var(ncid, "chan",   NF90_REAL, dimid_chan,          varid_k_tl_chan)
      errst = nf90_def_var(ncid, "k_tl_t", NF90_REAL, dimid_latlonlevchan, varid_k_tl_t)
      
      errst = nf90_put_att(ncid, varid_k_tl_lon,  "standard_name", "longitude")
      errst = nf90_put_att(ncid, varid_k_tl_lat,  "standard_name", "latitude")
      errst = nf90_put_att(ncid, varid_k_tl_lev,  "standard_name", "level")
      errst = nf90_put_att(ncid, varid_k_tl_chan, "standard_name", "chan")
      errst = nf90_put_att(ncid, varid_k_tl_t,    "standard_name", "k_tl_t")
      
      errst = nf90_put_att(ncid, varid_k_tl_lon,  "long_name", "Longitude")
      errst = nf90_put_att(ncid, varid_k_tl_lat,  "long_name", "Latitude")
      errst = nf90_put_att(ncid, varid_k_tl_lev,  "long_name", "Level index")
      errst = nf90_put_att(ncid, varid_k_tl_chan, "long_name", "Central wavelength of instrument channel")
      errst = nf90_put_att(ncid, varid_k_tl_t,    "long_name", "Temperature Jacobian")
      
      errst = nf90_put_att(ncid, varid_k_tl_lon,  "units", "degrees_east")
      errst = nf90_put_att(ncid, varid_k_tl_lat,  "units", "degrees_north")
      errst = nf90_put_att(ncid, varid_k_tl_lev,  "units", "")
      errst = nf90_put_att(ncid, varid_k_tl_chan, "units", "um")
      errst = nf90_put_att(ncid, varid_k_tl_t,    "units", "W/m2/sr/um/K")
      
      attr_instrument = trim(s3com%opt%rttov%inst_name)//"/"//trim(s3com%opt%rttov%platform_name)
      
      errst = nf90_put_att(ncid, nf90_global, 'RTTOV_instrument', trim(attr_instrument))
      
      errst = nf90_enddef(ncid)
      
      errst = nf90_put_var(ncid, varid_k_tl_lon,  model%lon_orig)
      errst = nf90_put_var(ncid, varid_k_tl_lat,  model%lat_orig)
      errst = nf90_put_var(ncid, varid_k_tl_lev,  model%height_2)
      errst = nf90_put_var(ncid, varid_k_tl_chan, s3com%rad%wavelength)
      errst = nf90_put_var(ncid, varid_k_tl_t,    gridded_k_tl_t)
      
      errst = nf90_close(ncid)
      
   end subroutine write_output_k_tl
   ! ============================================================================================================================
   
   ! ============================================================================================================================
   !> @brief Write atmospheric outputs calculated from ICON-LEM
   !! @details This subroutine outputs atmospheric data calculated from ICON-LEM
   !! @param[in] s3com   s3com structure
   !! @param[in] model   model structure
   !! @param[in] nml     namelist structure
   subroutine write_output_atm(s3com, model, nml)
      
      ! Input
      type(type_s3com), intent(in) :: s3com
      type(type_model), intent(in) :: model
      type(type_nml), intent(in) :: nml
      
      ! Internal
      real(wp), dimension(model%nlon, model%nlat) :: &
         gridded_atm_lwp, gridded_atm_re_ztop, gridded_atm_cdnc_ztop
      
      integer(wpi) :: ncid, errst
      
      integer(wpi) ::       &
         varid_atm_lon,     &
         varid_atm_lat,     &
         varid_atm_lwp,     &
         varid_atm_re_ztop, &
         varid_atm_cdnc_ztop
      
      integer(wpi) ::    &
         dimid_lon,      &
         dimid_lat,      &
         dimid_latlon(2)
      
      character(len=256) :: fn_out_atm, suffix, attr_instrument
      
      suffix = trim(nml%suffix_out)
      if (trim(suffix) .ne. "") suffix = "_"//trim(suffix)
      
      fn_out_atm = trim(nml%path_out)//"S3COM"//trim(suffix)//"_atm.nc"
      
      call map_point_to_ll(model%nlon, model%nlat, model%mode, x1=model%lwp_stratocumulus_filter, y2=gridded_atm_lwp)
      call map_point_to_ll(model%nlon, model%nlat, model%mode, x1=model%reff_top,                 y2=gridded_atm_re_ztop)
      call map_point_to_ll(model%nlon, model%nlat, model%mode, x1=model%cdnc_top,                 y2=gridded_atm_cdnc_ztop)
      
      errst = nf90_create(fn_out_atm, NF90_CLOBBER, ncid)
      
      errst = nf90_def_dim(ncid, "lon", model%nlon, dimid_lon)
      errst = nf90_def_dim(ncid, "lat", model%nlat, dimid_lat)
      
      dimid_latlon = (/dimid_lon, dimid_lat/)
      
      errst = nf90_def_var(ncid, "lon",       NF90_REAL, dimid_lon,    varid_atm_lon)
      errst = nf90_def_var(ncid, "lat",       NF90_REAL, dimid_lat,    varid_atm_lat)
      errst = nf90_def_var(ncid, "lwp",       NF90_REAL, dimid_latlon, varid_atm_lwp)
      errst = nf90_def_var(ncid, "re_ztop",   NF90_REAL, dimid_latlon, varid_atm_re_ztop)
      errst = nf90_def_var(ncid, "cdnc_ztop", NF90_REAL, dimid_latlon, varid_atm_cdnc_ztop)
      
      errst = nf90_put_att(ncid, varid_atm_lon,       "standard_name", "longitude")
      errst = nf90_put_att(ncid, varid_atm_lat,       "standard_name", "latitude")
      errst = nf90_put_att(ncid, varid_atm_lwp,       "standard_name", "lwp")
      errst = nf90_put_att(ncid, varid_atm_re_ztop,   "standard_name", "re_ztop")
      errst = nf90_put_att(ncid, varid_atm_cdnc_ztop, "standard_name", "cdnc_ztop")
      
      errst = nf90_put_att(ncid, varid_atm_lon,       "long_name", "Longitude")
      errst = nf90_put_att(ncid, varid_atm_lat,       "long_name", "Latitude")
      errst = nf90_put_att(ncid, varid_atm_lwp,       "long_name", "Liquid water path")
      errst = nf90_put_att(ncid, varid_atm_re_ztop,   "long_name", "Cloud droplet effective radius at cloud top")
      errst = nf90_put_att(ncid, varid_atm_cdnc_ztop, "long_name", "Cloud droplet number concentration at cloud top")
      
      errst = nf90_put_att(ncid, varid_atm_lon,       "units", "degrees_east")
      errst = nf90_put_att(ncid, varid_atm_lat,       "units", "degrees_north")
      errst = nf90_put_att(ncid, varid_atm_lwp,       "units", "kg/m2")
      errst = nf90_put_att(ncid, varid_atm_re_ztop,   "units", "um")
      errst = nf90_put_att(ncid, varid_atm_cdnc_ztop, "units", "m-3")
      
      attr_instrument = trim(s3com%opt%rttov%inst_name)//"/"//trim(s3com%opt%rttov%platform_name)
      
      errst = nf90_put_att(ncid, nf90_global, 'RTTOV_instrument', trim(attr_instrument))
      
      errst = nf90_enddef(ncid)
      
      errst = nf90_put_var(ncid, varid_atm_lon,       model%lon_orig)
      errst = nf90_put_var(ncid, varid_atm_lat,       model%lat_orig)
      errst = nf90_put_var(ncid, varid_atm_lwp,       gridded_atm_lwp)
      errst = nf90_put_var(ncid, varid_atm_re_ztop,   gridded_atm_re_ztop)
      errst = nf90_put_var(ncid, varid_atm_cdnc_ztop, gridded_atm_cdnc_ztop)
      
      errst = nf90_close(ncid)
      
   end subroutine write_output_atm
   ! ============================================================================================================================
   
   ! ============================================================================================================================
   !> @brief Write retrieval outputs simulated by S3COM
   !! @details This subroutines outputs the observed and retrieved cloud properties simulated by S3COM
   !! @param[in] s3com   s3com structure
   !! @param[in] model   model structure
   !! @param[in] nml     namelist structure
   subroutine write_output_ret(s3com, model, nml)
      
      ! Input
      type(type_s3com), intent(in) :: s3com
      type(type_model), intent(in) :: model
      type(type_nml), intent(in) :: nml
      
      ! Internal
      real(wp), dimension(model%nlon, model%nlat) :: &
         gridded_lwp_model,                          &
         gridded_cdnc_top_model,                     &
         gridded_lwp_ret,                            &
         gridded_cdnc_top_ret,                       &
         gridded_J
      
      integer(wpi) :: ncid, errst
      
      integer(wpi) ::        &
         varid_lon,          &
         varid_lat,          &
         varid_lwp_mod,      &
         varid_cdnc_top_mod, &
         varid_lwp_ret,      &
         varid_cdnc_top_ret, &
         varid_J
      
      integer(wpi) :: dimid_lon, dimid_lat, dimid_latlon(2)
      
      character(len=256) :: fn_out_ret, suffix, attr_instrument
      
      suffix = trim(nml%suffix_out)
      if (trim(suffix) .ne. "") suffix = "_"//trim(suffix)
      
      fn_out_ret = trim(nml%path_out)//"S3COM"//trim(suffix)//"_ret.nc"
      
      call map_point_to_ll(model%nlon, model%nlat, model%mode, x1=model%lwp_stratocumulus_filter(:), y2=gridded_lwp_model)
      call map_point_to_ll(model%nlon, model%nlat, model%mode, x1=s3com%ret%X(:,1),                  y2=gridded_lwp_ret)
      call map_point_to_ll(model%nlon, model%nlat, model%mode, x1=model%cdnc_top(:),                 y2=gridded_cdnc_top_model)
      call map_point_to_ll(model%nlon, model%nlat, model%mode, x1=s3com%ret%X(:,2),                  y2=gridded_cdnc_top_ret)
      call map_point_to_ll(model%nlon, model%nlat, model%mode, x1=s3com%ret%J(:),                    y2=gridded_J)
      
      errst = nf90_create(fn_out_ret, nf90_clobber, ncid)
      
      errst = nf90_def_dim(ncid, "lon", model%nlon, dimid_lon)
      errst = nf90_def_dim(ncid, "lat", model%nlat, dimid_lat)
      
      dimid_latlon = (/dimid_lon, dimid_lat/)
      
      errst = nf90_def_var(ncid, "lon",          nf90_real, dimid_lon,    varid_lon)
      errst = nf90_def_var(ncid, "lat",          nf90_real, dimid_lat,    varid_lat)
      errst = nf90_def_var(ncid, "lwp_mod",      nf90_real, dimid_latlon, varid_lwp_mod)
      errst = nf90_def_var(ncid, "lwp_ret",      nf90_real, dimid_latlon, varid_lwp_ret)
      errst = nf90_def_var(ncid, "cdnc_top_mod", nf90_real, dimid_latlon, varid_cdnc_top_mod)
      errst = nf90_def_var(ncid, "cdnc_top_ret", nf90_real, dimid_latlon, varid_cdnc_top_ret)
      errst = nf90_def_var(ncid, "J",            nf90_real, dimid_latlon, varid_J)
      
      errst = nf90_put_att(ncid, varid_lon,          "standard_name", "longitude")
      errst = nf90_put_att(ncid, varid_lat,          "standard_name", "latitude")
      errst = nf90_put_att(ncid, varid_lwp_mod,      "standard_name", "lwp_mod")
      errst = nf90_put_att(ncid, varid_lwp_ret,      "standard_name", "lwp_ret")
      errst = nf90_put_att(ncid, varid_cdnc_top_mod, "standard_name", "cdnc_top_mod")
      errst = nf90_put_att(ncid, varid_cdnc_top_ret, "standard_name", "cdnc_top_ret")
      errst = nf90_put_att(ncid, varid_J,            "standard_name", "J")
      
      errst = nf90_put_att(ncid, varid_lon,          "long_name", "Longitude")
      errst = nf90_put_att(ncid, varid_lat,          "long_name", "Latitude")
      errst = nf90_put_att(ncid, varid_lwp_mod,      "long_name", "Modelled cloud liquid water path")
      errst = nf90_put_att(ncid, varid_lwp_ret,      "long_name", "Retrieved cloud liquid water path")
      errst = nf90_put_att(ncid, varid_cdnc_top_mod, "long_name", "Modelled cloud droplet number concentration at cloud top")
      errst = nf90_put_att(ncid, varid_cdnc_top_ret, "long_name", "Retrieved cloud droplet number concentration at cloud top")
      errst = nf90_put_att(ncid, varid_J,            "long_name", "Cost function")
      
      errst = nf90_put_att(ncid, varid_lon,          "units", "degrees_east")
      errst = nf90_put_att(ncid, varid_lat,          "units", "degrees_north")
      errst = nf90_put_att(ncid, varid_lwp_mod,      "units", "kg/m2")
      errst = nf90_put_att(ncid, varid_lwp_ret,      "units", "kg/m2")
      errst = nf90_put_att(ncid, varid_cdnc_top_mod, "units", "m-3")
      errst = nf90_put_att(ncid, varid_cdnc_top_ret, "units", "m-3")
      errst = nf90_put_att(ncid, varid_J,            "units", "-")
      
      attr_instrument = trim(s3com%opt%rttov%inst_name)//"/"//trim(s3com%opt%rttov%platform_name)
      
      errst = nf90_put_att(ncid, nf90_global, 'RTTOV_instrument', trim(attr_instrument))
      
      errst = nf90_enddef(ncid)
      
      errst = nf90_put_var(ncid, varid_lon,          model%lon_orig)
      errst = nf90_put_var(ncid, varid_lat,          model%lat_orig)
      errst = nf90_put_var(ncid, varid_lwp_mod,      gridded_lwp_model)
      errst = nf90_put_var(ncid, varid_lwp_ret,      gridded_lwp_ret)
      errst = nf90_put_var(ncid, varid_cdnc_top_mod, gridded_cdnc_top_model)
      errst = nf90_put_var(ncid, varid_cdnc_top_ret, gridded_cdnc_top_ret)
      errst = nf90_put_var(ncid, varid_J,            gridded_J)
      
      errst = nf90_close(ncid)
      
   end subroutine write_output_ret
   ! ============================================================================================================================
   
end module mod_io_s3com

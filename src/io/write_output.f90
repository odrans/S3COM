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

MODULE MOD_WRITE_OUTPUT
   
   USE s3com_types,  ONLY: type_s3com, type_icon, wp, type_nml
   USE netcdf
   
   IMPLICIT NONE
   
   CONTAINS
      
      SUBROUTINE write_output(icon, oe, nml)
         
         USE mod_read_icon, ONLY: map_point_to_ll
         
         ! Input variables
         TYPE(type_icon),    INTENT(IN) :: icon
         TYPE(type_s3com),      INTENT(IN) :: oe
         TYPE(type_nml),      INTENT(IN) :: nml

         ! Local variables
         REAL(KIND=wp), DIMENSION(icon%Nlon, icon%Nlat, nml%nchannels) :: &
            gridded_y_refl_total, gridded_y_refl_clear,                     &
            gridded_y_bt_total, gridded_y_bt_clear,                         &
            gridded_y_rad_total, gridded_y_rad_clear, gridded_y_rad_cloudy, &
            gridded_brdf, gridded_emiss
         REAL(KIND=wp), DIMENSION(icon%Nlon, icon%Nlat) :: gridded_iwp_model, gridded_iwp_ret, gridded_g
         
         INTEGER(KIND=4) :: ichannel, ncid, errst
         INTEGER(KIND=4) ::                     &
            varid_lon, varid_lat, varid_chan,   &
            varid_refl_total, varid_refl_clear, &
            varid_bt_total, varid_bt_clear,     &
            varid_brdf, varid_emiss,            &
            varid_iwp_ret, varid_iwp_mod,       &
            varid_g,                            &
            varid_rad_total, varid_rad_clear, varid_rad_cloudy
         INTEGER(KIND=4) :: dimid_lon, dimid_lat, dimid_chan, dimid_latlon(2), dimid_latlonchan(3)
         
         !DO ichannel=1,rttov_nchannels
         !   call map_point_to_ll(icon%Nlon,icon%Nlat,icon%mode,x1=oe%y(:,ichannel),y2=gridded_y(:,:,ichannel))
         !ENDDO
         
         CALL map_point_to_ll(icon%Nlon, icon%Nlat, icon%mode, x2=oe%y_refl_total, y3=gridded_y_refl_total)
         CALL map_point_to_ll(icon%Nlon, icon%Nlat, icon%mode, x2=oe%y_refl_clear, y3=gridded_y_refl_clear)
         CALL map_point_to_ll(icon%Nlon, icon%Nlat, icon%mode, x2=oe%y_bt_total,   y3=gridded_y_bt_total)
         CALL map_point_to_ll(icon%Nlon, icon%Nlat, icon%mode, x2=oe%y_bt_clear,   y3=gridded_y_bt_clear)
         CALL map_point_to_ll(icon%Nlon, icon%Nlat, icon%mode, x2=oe%y_rad_total,  y3=gridded_y_rad_total)
         CALL map_point_to_ll(icon%Nlon, icon%Nlat, icon%mode, x2=oe%y_rad_clear,  y3=gridded_y_rad_clear)
         CALL map_point_to_ll(icon%Nlon, icon%Nlat, icon%mode, x2=oe%y_rad_cloudy, y3=gridded_y_rad_cloudy)
         CALL map_point_to_ll(icon%Nlon, icon%Nlat, icon%mode, x2=oe%brdf,         y3=gridded_brdf)
         CALL map_point_to_ll(icon%Nlon, icon%Nlat, icon%mode, x2=oe%emissivity,   y3=gridded_emiss)
         CALL map_point_to_ll(icon%Nlon, icon%Nlat, icon%mode, x1=oe%Xip1(:,1),    y2=gridded_iwp_ret)
         CALL map_point_to_ll(icon%Nlon, icon%Nlat, icon%mode, x1=oe%iwp_model(:), y2=gridded_iwp_model)
         CALL map_point_to_ll(icon%Nlon, icon%Nlat, icon%mode, x1=oe%gip1(:),      y2=gridded_g)   
         
         errst = nf90_create(nml%fname_out,NF90_CLOBBER,ncid)
         
         errst = nf90_def_dim(ncid, "Longitude", icon%Nlon,       dimid_lon)
         errst = nf90_def_dim(ncid, "Latitude",  icon%Nlat,       dimid_lat)
         errst = nf90_def_dim(ncid, "Channel",   nml%nchannels, dimid_chan)
         
         dimid_latlon     = (/dimid_lon, dimid_lat/)
         dimid_latlonchan = (/dimid_lon, dimid_lat, dimid_chan/)
         
         errst = nf90_def_var(ncid, "Longitude",       NF90_REAL, dimid_lon,        varid_lon)
         errst = nf90_def_var(ncid, "Latitude",        NF90_REAL, dimid_lat,        varid_lat)
         errst = nf90_def_var(ncid, "Channel",         NF90_REAL, dimid_chan,       varid_chan)
         errst = nf90_def_var(ncid, "BRF_total",       NF90_REAL, dimid_latlonchan, varid_refl_total)
         errst = nf90_def_var(ncid, "BRF_clear",       NF90_REAL, dimid_latlonchan, varid_refl_clear)
         errst = nf90_def_var(ncid, "BT_total",        NF90_REAL, dimid_latlonchan, varid_bt_total)
         errst = nf90_def_var(ncid, "BT_clear",        NF90_REAL, dimid_latlonchan, varid_bt_clear)
         errst = nf90_def_var(ncid, "Radiance_total",  NF90_REAL, dimid_latlonchan, varid_rad_total)
         errst = nf90_def_var(ncid, "Radiance_clear",  NF90_REAL, dimid_latlonchan, varid_rad_clear)
         errst = nf90_def_var(ncid, "Radiance_cloudy", NF90_REAL, dimid_latlonchan, varid_rad_cloudy)
         errst = nf90_def_var(ncid, "BRDF",            NF90_REAL, dimid_latlonchan, varid_brdf)
         errst = nf90_def_var(ncid, "Emissivity",      NF90_REAL, dimid_latlonchan, varid_emiss)
         errst = nf90_def_var(ncid, "iwp_ret",         NF90_REAL, dimid_latlon,     varid_iwp_ret)
         errst = nf90_def_var(ncid, "iwp_model",       NF90_REAL, dimid_latlon,     varid_iwp_mod)

         errst = nf90_put_att(ncid, varid_lon,        "units", "degrees_east")
         errst = nf90_put_att(ncid, varid_lat,        "units", "degrees_north")
         errst = nf90_put_att(ncid, varid_chan,       "units", "")
         errst = nf90_put_att(ncid, varid_refl_total, "units", "")
         errst = nf90_put_att(ncid, varid_refl_clear, "units", "")
         errst = nf90_put_att(ncid, varid_bt_total,   "units", "K")
         errst = nf90_put_att(ncid, varid_bt_clear,   "units", "K")
         errst = nf90_put_att(ncid, varid_rad_total,  "units", "W/m2/sr/um")
         errst = nf90_put_att(ncid, varid_rad_clear,  "units", "W/m2/sr/um")
         errst = nf90_put_att(ncid, varid_rad_cloudy, "units", "W/m2/sr/um")
         errst = nf90_put_att(ncid, varid_brdf,       "units", "sr-1")
         errst = nf90_put_att(ncid, varid_emiss,      "units", "")
         errst = nf90_put_att(ncid, varid_iwp_ret,    "units", "kg/m2")
         errst = nf90_put_att(ncid, varid_iwp_mod,    "units", "kg/m2")

         errst = nf90_enddef(ncid)
         
         errst = nf90_put_var(ncid, varid_lon,        icon%lon_orig)
         errst = nf90_put_var(ncid, varid_lat,        icon%lat_orig)
         errst = nf90_put_var(ncid, varid_chan,       nml%channel_list)
         errst = nf90_put_var(ncid, varid_refl_total, gridded_y_refl_total(:,:,:))
         errst = nf90_put_var(ncid, varid_refl_clear, gridded_y_refl_clear(:,:,:))
         errst = nf90_put_var(ncid, varid_bt_total,   gridded_y_bt_total(:,:,:))
         errst = nf90_put_var(ncid, varid_bt_clear,   gridded_y_bt_clear(:,:,:))
         errst = nf90_put_var(ncid, varid_rad_total,  gridded_y_rad_total(:,:,:))
         errst = nf90_put_var(ncid, varid_rad_clear,  gridded_y_rad_clear(:,:,:))
         errst = nf90_put_var(ncid, varid_rad_cloudy, gridded_y_rad_cloudy(:,:,:))
         errst = nf90_put_var(ncid, varid_brdf,       gridded_brdf(:,:,:))
         errst = nf90_put_var(ncid, varid_emiss,      gridded_emiss(:,:,:))
         errst = nf90_put_var(ncid, varid_iwp_ret,    gridded_iwp_ret(:,:))
         errst = nf90_put_var(ncid, varid_iwp_mod,    gridded_iwp_model(:,:))

         errst = nf90_close(ncid)
         
         write(*,*) "Done writting"
         
      END SUBROUTINE write_output

END MODULE MOD_WRITE_OUTPUT

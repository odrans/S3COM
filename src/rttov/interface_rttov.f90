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


MODULE MOD_RTTOV_INTERFACE

  USE s3com_types,  ONLY: wp, type_rttov_opt, type_nml

  USE rttov_const, ONLY: &
       surftype_sea,      &
       surftype_land,     &
       sensor_id_mw,      &
       sensor_id_po,      &
       inst_name,         &
       platform_name

  USE mod_rttov, ONLY: platform, satellite, sensor, nChannels,       &
       opts, errorstatus_success, rttov_exit, coefs,                  &
       emis_atlas, brdf_atlas, atlas_type, dosolar, imonth, &
       dosolar, channel_list

  !!The rttov_emis_atlas_data type must be imported separately
  USE mod_rttov_emis_atlas, ONLY:  &
       rttov_emis_atlas_data,       &
       atlas_type_ir, atlas_type_mw

  !!The rttov_brdf_atlas_data type must be imported separately
  USE mod_rttov_brdf_atlas, ONLY : rttov_brdf_atlas_data

  IMPLICIT NONE

#include "rttov_direct.interface"
#include "rttov_parallel_direct.interface"
#include "rttov_read_coefs.interface"
#include "rttov_dealloc_coefs.interface"
#include "rttov_alloc_direct.interface"
#include "rttov_user_options_checkinput.interface"
#include "rttov_print_opts.interface"
#include "rttov_print_profile.interface"
#include "rttov_skipcommentline.interface"

  !!Use emissivity atlas
#include "rttov_setup_emis_atlas.interface"
#include "rttov_get_emis.interface"
#include "rttov_deallocate_emis_atlas.interface"

  !!Use BRDF atlas
#include "rttov_setup_brdf_atlas.interface"
#include "rttov_get_brdf.interface"
#include "rttov_deallocate_brdf_atlas.interface"

CONTAINS

  SUBROUTINE RTTOV_INIT(rttov_opt, nml)

    TYPE(type_rttov_opt), INTENT(in) :: rttov_opt
    TYPE(type_nml), intent(IN) :: nml

    !!Local variables
    CHARACTER(len=256) :: coef_filename, cld_coef_filename, sat, path_emis_atlas, path_brdf_atlas, path_rttov_2
    INTEGER :: errorstatus

    imonth    = rttov_opt%month
    dosolar   = rttov_opt%dosolar
    nChannels = rttov_opt%nchannels

    ALLOCATE(channel_list(nchannels))
    channel_list(1:nchannels)=rttov_opt%channel_list

    IF (rttov_opt%satellite .NE. 0) THEN
       WRITE(sat,*) rttov_opt%satellite
       sat="_"//trim(adjustl(sat))//"_"
    END IF

    coef_filename = trim(nml%path_rttov)//"/rtcoef_rttov13/rttov13pred54L/rtcoef_"//&
         trim(platform_name(rttov_opt%platform))//trim(sat)//trim(inst_name(rttov_opt%instrument))//"_o3.dat"

    cld_coef_filename = trim(nml%path_rttov)//"/rtcoef_rttov13/cldaer_visir/sccldcoef_"//&
         trim(platform_name(rttov_opt%platform))//trim(sat)//trim(inst_name(rttov_opt%instrument))//".dat"

    path_emis_atlas = trim(nml%path_rttov)//'/emis_data'
    path_brdf_atlas = trim(nml%path_rttov)//'/brdf_data'

    !!-----------------------------------------------------------------------------------------------------------------------!!
    !! 1. Initialise RTTOV options structure                                                                                 !!
    !!-----------------------------------------------------------------------------------------------------------------------!!

    IF (dosolar == 1) THEN
       opts%rt_ir%addsolar = .TRUE.              !Solar radiation included
    ELSE
       opts%rt_ir%addsolar = .FALSE.             !Solar radiation not included (default = false)
    ENDIF

    opts%interpolation%addinterp       = .TRUE.  !If true input profiles may be supplied on user-defined levels, and internal
    !interpolation is used (default = false)
    opts%interpolation%interp_mode     = 1       !Interpolation method
    !opts%interpolation%reg_limit_extrap = .TRUE.

    opts%rt_all%addrefrac              =  nml%addrefrac  !If true RTTOV calculations accounts for atmospheric refraction (default = true)
    opts%rt_ir%addaerosl               = .FALSE. !If true accounts for scattering due to aerosols (default = false)
    opts%rt_ir%addclouds               = .TRUE.  !If true accounts for scattering due to clouds (default = false)

    opts%rt_ir%ir_scatt_model          = nml%ir_scatt_model      !Scattering model for emission source term:
    !1 => DOM; 2 => Chou-scaling
    opts%rt_ir%vis_scatt_model         = nml%vis_scatt_model       !Scattering model for solar source term:
    !1 => DOM; 2 => single-scattering; 3 => MFASIS
    opts%rt_ir%dom_nstreams            = nml%dom_nstreams       !Number of streams for Discrete Ordinates (DOM)
    opts%rt_ir%dom_rayleigh            = nml%dom_rayleigh !Enables Rayleigh multiple-scattering in solar DOM simulations

    opts%rt_all%ozone_data             = .FALSE. !Set the relevant flag to .TRUE. when supplying a profile of the given
    opts%rt_all%co2_data               = .FALSE. !trace gas (ensure the coefficient file supports the gas)
    opts%rt_all%n2o_data               = .FALSE.
    opts%rt_all%ch4_data               = .FALSE.
    opts%rt_all%co_data                = .FALSE.
    opts%rt_all%so2_data               = .FALSE.

    opts%rt_mw%clw_data                = .FALSE.

    opts%config%verbose                = .FALSE.  !If false only messages for fatal errors are output (default = true)
    opts%config%do_checkinput          = .TRUE. !If true checks whether input profiles are within both absolute and regression
    !limits (default = true)

    opts%rt_all%switchrad              = .TRUE.

    !!-----------------------------------------------------------------------------------------------------------------------!!
    !! 2. Read coefficients                                                                                                  !!
    !!-----------------------------------------------------------------------------------------------------------------------!!

    !!Read optical depth and cloud coefficient files together
    CALL rttov_read_coefs(errorstatus, coefs, opts, file_coef=coef_filename, file_sccld=cld_coef_filename)
    IF (errorstatus /= errorstatus_success) THEN
       WRITE(*,*) 'fatal error reading coefficients'
       CALL rttov_exit(errorstatus)
    ENDIF

    !!Ensure input number of channels is not higher than number stored in coefficient file
    IF (nchannels > coefs%coef%fmv_chn) THEN
       nchannels = coefs%coef%fmv_chn
    ENDIF

    !!Ensure the options and coefficients are consistent
    CALL rttov_user_options_checkinput(errorstatus, opts, coefs)
    IF (errorstatus /= errorstatus_success) THEN
       WRITE(*,*) 'error in rttov options'
       CALL rttov_exit(errorstatus)
    ENDIF

    !!Initialise the RTTOV emissivity atlas
    !This loads the default IR/MW atlases: use the atlas_id argument to select alternative atlases
    IF (coefs%coef%id_sensor == sensor_id_mw .OR. coefs%coef%id_sensor == sensor_id_po) THEN
       atlas_type = atlas_type_mw !MW atlas
    ELSE
       atlas_type = atlas_type_ir !IR atlas
    ENDIF

    CALL rttov_setup_emis_atlas(                               &
         errorstatus,                                          &
         opts,                                                 &
         imonth,                                               &
         atlas_type,                                           & !Selects MW (1) or IR (2)
         emis_atlas,                                           &
         path = path_emis_atlas, & !Default path to atlas data
         coefs = coefs)                                          !This is mandatory for the CNRM MW atlas, ignored by TELSEM2;
    !if supplied for IR atlases they are initialised for this
    !sensor and this makes the atlas much faster to access

    IF (errorstatus /= errorstatus_success) THEN
       WRITE(*,*) 'error initialising emissivity atlas'
       CALL rttov_exit(errorstatus)
    ENDIF

    IF (opts%rt_ir%addsolar) THEN
       !!Initialise the RTTOV BRDF atlas
       CALL rttov_setup_brdf_atlas(                              &
            errorstatus,                                         &
            opts,                                                &
            imonth,                                              &
            brdf_atlas,                                          &
            path = path_brdf_atlas, & !Default path to atlas data
            coefs = coefs)                                         !If supplied the BRDF atlas is initialised for this sensor
       !and this makes the atlas much faster to access
       !
       IF (errorstatus /= errorstatus_success) THEN
          WRITE(*,*) 'error initialising BRDF atlas'
          CALL rttov_exit(errorstatus)
       ENDIF

    ENDIF

  END SUBROUTINE RTTOV_INIT

END MODULE  MOD_RTTOV_INTERFACE

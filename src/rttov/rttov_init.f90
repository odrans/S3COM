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

! Initialize a few structures for RTTOV (opts, coefs, emis_atlas, brdf_atlas)
module mod_rttov_interface

  use s3com_types,  only: wp, type_rttov_opt, type_nml, type_s3com

  use rttov_const, only: surftype_sea, surftype_land, sensor_id_mw, &
       sensor_id_po, inst_name, platform_name, errorstatus_success

  use mod_rttov, only: opts, coefs, emis_atlas, brdf_atlas

  !!The rttov_emis_atlas_data type must be imported separately
  use mod_rttov_emis_atlas, only: rttov_emis_atlas_data, atlas_type_ir, atlas_type_mw

  !!The rttov_brdf_atlas_data type must be imported separately
  use mod_rttov_brdf_atlas, only : rttov_brdf_atlas_data

  use rttov_unix_env, only: rttov_exit

  implicit none

  private
  public :: rttov_init

#include "rttov_direct.interface"
#include "rttov_parallel_direct.interface"
#include "rttov_alloc_direct.interface"
#include "rttov_k.interface"
#include "rttov_parallel_k.interface"
#include "rttov_alloc_k.interface"
#include "rttov_read_coefs.interface"
#include "rttov_dealloc_coefs.interface"
#include "rttov_init_emis_refl.interface"
#include "rttov_init_prof.interface"
#include "rttov_init_rad.interface"
#include "rttov_init_transmission.interface"
#include "rttov_user_options_checkinput.interface"
#include "rttov_print_opts.interface"
#include "rttov_print_profile.interface"
#include "rttov_skipcommentline.interface"

  ! Use emissivity atlas
#include "rttov_setup_emis_atlas.interface"
#include "rttov_get_emis.interface"
#include "rttov_deallocate_emis_atlas.interface"

  ! Use BRDF atlas
#include "rttov_setup_brdf_atlas.interface"
#include "rttov_get_brdf.interface"
#include "rttov_deallocate_brdf_atlas.interface"

contains

  subroutine rttov_init(rttov_opt, s3com)

    type(type_rttov_opt), intent(in) :: rttov_opt

    type(type_s3com), intent(inout) :: s3com

    !!Local variables
    character(len=256) :: coef_filename, cld_coef_filename, sat, path_emis_atlas, path_brdf_atlas, path_rttov_2
    integer(kind=4) :: errorstatus, imonth, nchannels, atlas_type
    
    imonth    = rttov_opt%month
    nChannels = rttov_opt%nchannels

    if (rttov_opt%satellite .ne. 0) then
       write(sat,*) rttov_opt%satellite
       sat="_"//trim(adjustl(sat))//"_"
    end if

    coef_filename = trim(s3com%nml%path_rttov)//"/rtcoef_rttov13/rttov13pred54L/rtcoef_"//&
         trim(platform_name(rttov_opt%platform))//trim(sat)//trim(inst_name(rttov_opt%instrument))//"_o3.dat"

    cld_coef_filename = trim(s3com%nml%path_rttov)//"/rtcoef_rttov13/cldaer_visir/sccldcoef_"//&
         trim(platform_name(rttov_opt%platform))//trim(sat)//trim(inst_name(rttov_opt%instrument))//".dat"

    path_emis_atlas = trim(s3com%nml%path_rttov)//'/emis_data'
    path_brdf_atlas = trim(s3com%nml%path_rttov)//'/brdf_data'

    s3com%opt%rttov%platform_name = platform_name(rttov_opt%platform)
    s3com%opt%rttov%inst_name = inst_name(rttov_opt%instrument)

    !!-----------------------------------------------------------------------------------------------------------------------!!
    !! 1. Initialise RTTOV options structure                                                                                 !!
    !!-----------------------------------------------------------------------------------------------------------------------!!

    if (rttov_opt%dosolar == 1) then
       opts%rt_ir%addsolar = .true.              !Solar radiation included
    else
       opts%rt_ir%addsolar = .false.             !Solar radiation not included (default = false)
    endif

    opts%interpolation%addinterp       = .true.  !If true input profiles may be supplied on user-defined levels, and internal
    !interpolation is used (default = false)
    opts%interpolation%interp_mode     = 4       !Interpolation method
    opts%interpolation%lgradp = .false.

    opts%dev%do_opdep_calc             = s3com%nml%do_opdep_calc !If false disables the RTTOV gas optical depth calculation (default = true)
    opts%rt_all%addrefrac              = s3com%nml%addrefrac  !If true RTTOV calculations accounts for atmospheric refraction (default = true)
    opts%rt_ir%addaerosl               = .false. !If true accounts for scattering due to aerosols (default = false)
    opts%rt_ir%addclouds               = .true.  !If true accounts for scattering due to clouds (default = false)

    opts%rt_ir%ir_scatt_model          = s3com%nml%ir_scatt_model      !Scattering model for emission source term:
    !1 => DOM; 2 => Chou-scaling
    opts%rt_ir%vis_scatt_model         = s3com%nml%vis_scatt_model       !Scattering model for solar source term:
    !1 => DOM; 2 => single-scattering; 3 => MFASIS
    opts%rt_ir%dom_nstreams            = s3com%nml%dom_nstreams       !Number of streams for Discrete Ordinates (DOM)
    opts%rt_ir%dom_rayleigh            = s3com%nml%dom_rayleigh !Enables Rayleigh multiple-scattering in solar DOM simulations

    opts%rt_all%ozone_data             = .false. !Set the relevant flag to .true. when supplying a profile of the given
    opts%rt_all%co2_data               = .false. !trace gas (ensure the coefficient file supports the gas)
    opts%rt_all%n2o_data               = .false.
    opts%rt_all%ch4_data               = .false.
    opts%rt_all%co_data                = .false.
    opts%rt_all%so2_data               = .false.

    opts%rt_mw%clw_data                = .false.

    opts%config%verbose                = .false.  !If false only messages for fatal errors are output (default = true)
    opts%config%do_checkinput          = .true. !If true checks whether input profiles are within both absolute and regression
    !limits (default = true)

    !opts%rt_all%switchrad              = .true. ! Input K perturbation in BT
    opts%rt_all%switchrad              = .false. ! Input K perturbation in radiance
    
    !!-----------------------------------------------------------------------------------------------------------------------!!
    !! 2. Read coefficients                                                                                                  !!
    !!-----------------------------------------------------------------------------------------------------------------------!!

    ! Read optical depth and cloud coefficient files together
    call rttov_read_coefs(errorstatus, coefs, opts, file_coef=coef_filename, file_sccld=cld_coef_filename)
    if (errorstatus /= errorstatus_success) then
       write(*,*) 'fatal error reading coefficients'
       call rttov_exit(errorstatus)
    endif

    ! Ensure input number of channels is not higher than number stored in coefficient file
    if (nchannels > coefs%coef%fmv_chn) then
       nchannels = coefs%coef%fmv_chn
    endif

    ! Ensure the options and coefficients are consistent
    call rttov_user_options_checkinput(errorstatus, opts, coefs)
    if (errorstatus /= errorstatus_success) then
       write(*,*) 'error in rttov options'
       call rttov_exit(errorstatus)
    endif

    ! Initialise the RTTOV emissivity atlas
    ! This loads the default IR/MW atlases: use the atlas_id argument to select alternative atlases
    if (coefs%coef%id_sensor == sensor_id_mw .or. coefs%coef%id_sensor == sensor_id_po) then
       atlas_type = atlas_type_mw !MW atlas
    else
       atlas_type = atlas_type_ir !IR atlas
    endif

    call rttov_setup_emis_atlas(                               &
         errorstatus,                                          &
         opts,                                                 &
         imonth,                                               &
         atlas_type,                                           & !Selects MW (1) or IR (2)
         emis_atlas,                                           &
         path = path_emis_atlas, & !Default path to atlas data
         coefs = coefs)                                          !This is mandatory for the CNRM MW atlas, ignored by TELSEM2;

    ! If supplied for IR atlases they are initialised for this
    ! sensor and this makes the atlas much faster to access

    if (errorstatus /= errorstatus_success) then
       write(*,*) 'error initialising emissivity atlas'
       call rttov_exit(errorstatus)
    endif

    if (opts%rt_ir%addsolar) then
       !!Initialise the RTTOV BRDF atlas
       call rttov_setup_brdf_atlas(                              &
            errorstatus,                                         &
            opts,                                                &
            imonth,                                              &
            brdf_atlas,                                          &
            path = path_brdf_atlas, & !Default path to atlas data
            coefs = coefs)                                         !If supplied the BRDF atlas is initialised for this sensor
       !and this makes the atlas much faster to access
       !
       if (errorstatus /= errorstatus_success) then
          write(*,*) 'error initialising BRDF atlas'
          call rttov_exit(errorstatus)
       endif

    endif

    s3com%rad%wavelength = 10000._wp / coefs%coef%ff_cwn(rttov_opt%channel_list(:))
    
  end subroutine rttov_init

end module mod_rttov_interface

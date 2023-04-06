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

module mod_rttov_opts

  use rttov_const, only: sensor_id_po, inst_name, platform_name
  use s3com_types,  only: wp, type_rttov_opt, type_nml, type_s3com
  use mod_rttov_utils, only: check_rttov_status
  use rttov_types, only: rttov_options

  implicit none

  private
  public :: rttov_opts_init, opts

#include "rttov_print_opts.interface"

  type(rttov_options)              :: opts                       ! Options structure

contains

  !> 1. Initialize the RTTOV options structure
  !> 2. The platform and instrument are set in the s3com structure
  subroutine rttov_opts_init(rttov_opt, s3com)

    type(type_rttov_opt), intent(in) :: rttov_opt

    type(type_s3com), intent(inout) :: s3com

    s3com%opt%rttov%platform_name = platform_name(rttov_opt%platform)
    s3com%opt%rttov%inst_name = inst_name(rttov_opt%instrument)

    opts%interpolation%addinterp       = .true.                       ! If true input profiles may be supplied on user-defined levels, and internal
    opts%interpolation%interp_mode     = 1                            ! Pressure-level interpolation method
    opts%interpolation%lgradp          = .false.                      ! If true, the pressure gradient is calculated from the input pressure levels

    opts%dev%do_opdep_calc             = rttov_opt%do_opdep_calc      ! If false, disables the RTTOV gas optical depth calculation (default = true)

    if (rttov_opt%dosolar == 1) then
       opts%rt_ir%addsolar = .true.                                   ! Solar radiation included (default = false)
    else
       opts%rt_ir%addsolar = .false.                                  ! Solar radiation not included (default = false)
    endif
    opts%rt_ir%addaerosl               = rttov_opt%add_aerosols       ! If true, aerosols are included in the RTTOV model
    opts%rt_ir%addclouds               = rttov_opt%add_clouds         ! If true, clouds are included in the RTTOV model
    opts%rt_ir%ir_scatt_model          = rttov_opt%ir_scatt_model     ! Scattering model for IR source term: 1 => DOM; 2 => Chou-scaling
    opts%rt_ir%vis_scatt_model         = rttov_opt%vis_scatt_model    ! Scattering model for solar source term: 1: DOM; 2: single-scattering; 3: MFASIS
    opts%rt_ir%dom_nstreams            = rttov_opt%dom_nstreams       ! Number of streams for DOM scattering model
    opts%rt_ir%dom_rayleigh            = rttov_opt%dom_rayleigh       ! If true, Rayleigh scattering is included in the DOM model
    opts%rt_ir%user_cld_opt_param      = rttov_opt%user_cld_opt_param ! If true, the user can supply cloud optical parameters

    opts%rt_all%ozone_data             = rttov_opt%ozone_data         ! Set the relevant flag to .true. when supplying a profile of the given
    opts%rt_all%co2_data               = rttov_opt%co2_data           ! gas. The profile must be supplied in the input profile structure
    opts%rt_all%n2o_data               = rttov_opt%n2o_data           ! If false, the gas is assumed to be constant and the value
    opts%rt_all%ch4_data               = rttov_opt%ch4_data           ! supplied in the RTTOV options structure is used.
    opts%rt_all%co_data                = rttov_opt%co_data
    opts%rt_all%so2_data               = rttov_opt%so2_data
    opts%rt_all%addrefrac              = rttov_opt%add_refrac         ! If true accounts for atmospheric refraction (default = false)
    opts%rt_all%switchrad              = .false.                      ! Input K perturbation in radiance (false) or BT (true)

    opts%rt_mw%clw_data                = .false.

    opts%config%apply_reg_limits       = .true.                       ! If true the regression limits are set as hard limits and replace the value. It will remove the warning in verbose.
    opts%config%verbose                = .false.                      ! If false only messages for fatal errors are output (default = true)
    opts%config%do_checkinput          = .false.                      ! If true checks whether input profiles are within both absolute and regression

    ! call rttov_print_opts(opts)

  end subroutine rttov_opts_init

end module mod_rttov_opts

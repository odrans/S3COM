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

module mod_rttov_coefs

  use s3com_types,  only: wp, type_rttov_opt, type_nml, type_s3com
  use mod_rttov_utils, only: check_rttov_status
  use rttov_const, only: inst_name, platform_name
  use mod_rttov_opts, only: opts
  use rttov_types, only: rttov_coefs

  use mod_io_verbose, only: verbose_rttov

  implicit none

  private
  public :: rttov_coefs_init, rttov_coefs_deinit, coefs

#include "rttov_dealloc_coefs.interface"
#include "rttov_read_coefs.interface"
#include "rttov_user_options_checkinput.interface"
#include "rttov_print_opts.interface"

  type(rttov_coefs)                :: coefs                      ! Coefficients structure

contains

  subroutine rttov_coefs_init(rttov_opt, s3com)

    type(type_rttov_opt), intent(in) :: rttov_opt

    type(type_s3com), intent(inout) :: s3com

    !!Local variables
    character(len=256) :: coef_filename, cld_coef_filename, sat, file_format
    integer(kind=4) :: errorstatus, nchannels

    nchannels = rttov_opt%nchannels

    if (rttov_opt%satellite .ne. 0) then
       write(sat,*) rttov_opt%satellite
       sat="_"//trim(adjustl(sat))//"_"
    end if
    s3com%opt%rttov%sat_name = sat

    file_format = ".dat"
    if(trim(inst_name(rttov_opt%instrument)) == "iasi") file_format = ".H5"

    coef_filename = trim(s3com%nml%path_rttov)//"/rtcoef_rttov13/rttov13pred54L/rtcoef_"//&
         trim(platform_name(rttov_opt%platform))//trim(sat)//trim(inst_name(rttov_opt%instrument))//"_o3"//trim(file_format)

    cld_coef_filename = trim(s3com%nml%path_rttov)//"/rtcoef_rttov13/cldaer_visir/sccldcoef_"//&
         trim(platform_name(rttov_opt%platform))//trim(sat)//trim(inst_name(rttov_opt%instrument))//trim(file_format)

    write(*,*) "coef_filename:", trim(cld_coef_filename)
    s3com%opt%rttov%platform_name = platform_name(rttov_opt%platform)
    s3com%opt%rttov%inst_name = inst_name(rttov_opt%instrument)

    ! Read optical depth and cloud coefficient files together
    call rttov_read_coefs(errorstatus, coefs, opts, file_coef=coef_filename, file_sccld=cld_coef_filename)
    call check_rttov_status(errorstatus, 'rttov_read_coefs')

    ! Ensure input number of channels is not higher than number stored in coefficient file
    if (nchannels > coefs%coef%fmv_chn) nchannels = coefs%coef%fmv_chn

    ! Ensure the options and coefficients are consistent
    call rttov_user_options_checkinput(errorstatus, opts, coefs)
    call check_rttov_status(errorstatus, 'rttov_user_options_checkinput')

    s3com%rad%wavelength = 10000._wp / real(coefs%coef%ff_cwn(rttov_opt%channel_list(:)), wp)

    call verbose_rttov(s3com, rttov_opt, opts)

    ! call rttov_print_info (coefs)

  end subroutine rttov_coefs_init

  subroutine rttov_coefs_deinit()

    integer(kind=4) :: errorstatus

    call rttov_dealloc_coefs(errorstatus, coefs)
    call check_rttov_status(errorstatus, 'rttov_dealloc_coefs')

  end subroutine rttov_coefs_deinit

end module mod_rttov_coefs

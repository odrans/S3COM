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

module mod_rttov_setup

  use s3com_types, only: wp

  implicit none

contains

  subroutine rttov_setup_opt(model, nml, zenangle, azangle, rttov_opt)

    use s3com_types,  only: type_rttov_opt, type_nml, type_model
    use s3com_config, only: RTTOV_DOSOLAR

    type(type_nml), intent(in) :: nml
    type(type_model), intent(in) :: model
    real(wp), intent(in) :: zenangle, azangle

    type(type_rttov_opt), intent(out) :: rttov_opt

    rttov_opt%platform   = nml%platform
    rttov_opt%satellite  = nml%satellite
    rttov_opt%instrument = nml%instrument
    rttov_opt%dosolar    = RTTOV_DOSOLAR
    rttov_opt%nchannels  = nml%nchannels
    rttov_opt%nthreads = nml%rttov_nthreads

    allocate(rttov_opt%channel_list(rttov_opt%nchannels))

    rttov_opt%channel_list = nml%channel_list
    rttov_opt%month        = model%date(2)
    rttov_opt%zenangle     = zenangle
    rttov_opt%azangle      = azangle

  end subroutine rttov_setup_opt

  subroutine rttov_setup_atm(idx_start, idx_end, model, rttov_atm)

    use s3com_types, only: type_model

    integer, intent(in)            :: idx_start, idx_end
    type(type_model), intent(in)   :: model

    type(type_model), intent(out)  :: rttov_atm

    integer :: nidx

    nidx = idx_end - idx_start + 1

    rttov_atm%npoints     = nidx
    rttov_atm%idx_start   = idx_start
    rttov_atm%idx_end     = idx_end
    rttov_atm%nlevels     = model%nlevels
    rttov_atm%nlayers     = model%nlayers

    rttov_atm%lat         = model%lat(idx_start:idx_end)
    rttov_atm%lon         = model%lon(idx_start:idx_end)
    rttov_atm%landmask    = model%landmask(idx_start:idx_end)
    rttov_atm%topography  = model%topography(idx_start:idx_end)
    rttov_atm%ps          = model%ps(idx_start:idx_end)
    rttov_atm%ts          = model%ts(idx_start:idx_end)
    rttov_atm%t_2m        = model%t_2m(idx_start:idx_end)
    rttov_atm%q_2m        = model%q_2m(idx_start:idx_end)
    rttov_atm%u_10m       = model%u_10m(idx_start:idx_end)
    rttov_atm%v_10m       = model%v_10m(idx_start:idx_end)
    rttov_atm%sunzenangle = model%sunzenangle(idx_start:idx_end)
    rttov_atm%sunazangle  = model%sunazangle(idx_start:idx_end)

    rttov_atm%p           = model%p(idx_start:idx_end,:)
    rttov_atm%z           = model%z(idx_start:idx_end,:)
    rttov_atm%t           = model%t(idx_start:idx_end,:)
    rttov_atm%q           = model%q(idx_start:idx_end,:)
    rttov_atm%clc         = model%clc(idx_start:idx_end,:)
    rttov_atm%iwc         = model%iwc(idx_start:idx_end,:)
    rttov_atm%lwc         = model%lwc(idx_start:idx_end,:)
    rttov_atm%reff        = model%reff(idx_start:idx_end,:)

  end subroutine rttov_setup_atm

end module mod_rttov_setup

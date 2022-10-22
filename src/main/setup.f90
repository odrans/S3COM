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

module mod_s3com_setup

  use s3com_types,  only: wp, type_nml, type_icon, type_model, type_s3com

  use s3com_config, only: nstates, apriori_iwp

  implicit none

  public :: s3com_init, s3com_subset, s3com_update

contains

  subroutine s3com_init(nml, model, s3com)

    type(type_nml), intent(IN)        :: nml
    type(type_model), intent(IN)      :: model

    type(type_s3com), intent(OUT)     :: s3com

    integer(KIND = 4) :: npoints, nlayers, nlevels, nstates, nmeas

    npoints = model%npoints
    nlevels = model%nlevels
    nlayers = model%nlayers
    nstates = 1
    nmeas = nml%nchannels

    s3com%npoints = npoints
    s3com%nlayers = nlayers
    s3com%nlevels = nlevels
    s3com%nstates = nstates
    s3com%nmeas = nmeas

    allocate(s3com%rad%wavelength(nmeas), source = 0._wp)

    allocate(s3com%rad%y(npoints, nmeas), source = 0._wp)
    allocate(s3com%rad%f, s3com%rad%f_ref_total, s3com%rad%f_ref_clear, s3com%rad%f_bt_total, &
         s3com%rad%f_bt_clear, s3com%rad%f_rad_total, s3com%rad%f_rad_clear, &
         mold = s3com%rad%y)


    s3com%nml = nml

    return

  end subroutine s3com_init

  subroutine s3com_subset(s3com, rttov_atm, oe)

    type(type_s3com), intent(IN)     :: s3com
    type(type_model), intent(IN)     :: rttov_atm

    type(type_s3com), intent(OUT) :: oe

    integer(KIND=4) :: idx_start, idx_end, npoints

    oe%idx_start = rttov_atm%idx_start
    oe%idx_end = rttov_atm%idx_end
    oe%npoints = oe%idx_end - oe%idx_start + 1

    allocate(oe%flag_rttov(oe%npoints)); oe%flag_rttov = .true.

    oe%npoints = s3com%npoints
    oe%nmeas = s3com%nmeas

    allocate(oe%rad%f_ref_total(oe%npoints, oe%nmeas), source = 0._wp)
    allocate(oe%rad%f_ref_clear, oe%rad%f_bt_total, &
         oe%rad%f_bt_clear, oe%rad%f_rad_total, oe%rad%f_rad_clear, &
         mold = oe%rad%f_ref_total)

    return

  end subroutine s3com_subset

  subroutine s3com_update(s3com, oe)

    type(type_s3com), intent(INOUT) :: s3com

    type(type_s3com), intent(INOUT) :: oe

    ! Radiation data
    s3com%rad%f_ref_total(oe%idx_start:oe%idx_end, 1:oe%nmeas) = oe%rad%f_ref_total(1:oe%npoints, 1:oe%nmeas)
    s3com%rad%f_ref_clear(oe%idx_start:oe%idx_end, 1:oe%nmeas) = oe%rad%f_ref_clear(1:oe%npoints, 1:oe%nmeas)
    s3com%rad%f_bt_total(oe%idx_start:oe%idx_end, 1:oe%nmeas) = oe%rad%f_bt_total(1:oe%npoints, 1:oe%nmeas)
    s3com%rad%f_bt_clear(oe%idx_start:oe%idx_end, 1:oe%nmeas) = oe%rad%f_bt_clear(1:oe%npoints, 1:oe%nmeas)
    s3com%rad%f_rad_total(oe%idx_start:oe%idx_end, 1:oe%nmeas) = oe%rad%f_rad_total(1:oe%npoints, 1:oe%nmeas)
    s3com%rad%f_rad_clear(oe%idx_start:oe%idx_end, 1:oe%nmeas) = oe%rad%f_rad_clear(1:oe%npoints, 1:oe%nmeas)

    deallocate(oe%rad%f_ref_total, oe%rad%f_ref_clear, oe%rad%f_bt_total, oe%rad%f_bt_clear, &
         oe%rad%f_rad_total, oe%rad%f_rad_clear)

    ! General
    deallocate(oe%flag_rttov)

    return

  end subroutine s3com_update

end module mod_s3com_setup

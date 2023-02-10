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

program s3com_main

  use s3com_types,         only: wp, type_rttov_opt, type_nml, type_model, type_s3com
  use mod_io_namelist,     only: namelist_load
  use mod_rttov_interface, only: rttov_init
  use mod_rttov_setup,     only: rttov_setup_opt, rttov_setup_atm
  use mod_rttov,           only: run_rttov
  use mod_s3com_setup,     only: s3com_init, s3com_subset, s3com_update
  use mod_models,          only: models_load, models_deinit
  use mod_utils_math,      only: n_chunks
  use mod_io_s3com,        only: write_output
  use mod_rttov_utils,     only: idx_rttov
  !  use mod_oe_utils,        only: idx_ice
  !  use mod_model_cloud,     only: init_zcloud, init_cloudprof
  !  use mod_oe_run,          only: oe_run

  implicit none

  type(type_model)              :: model, rttov_atm  
  type(type_rttov_opt)          :: rttov_opt
  type(type_s3com)              :: s3com, s3com_chunk
  type(type_nml)                :: nml

  real(kind=wp) :: zenangle, azangle

  integer(kind=4), dimension(:), allocatable :: idx_iwp, idx_oe
  integer(kind=4) :: nchunks, idx_start, idx_end, ichunk
  integer(kind=4) :: nlayers, npoints, npoints_it

  logical :: flag_oe

  ! Read namelist file
  nml = namelist_load()

  ! Temporary: setting the viewing satellite angles
  zenangle = 0._wp; azangle = 0._wp

  ! Load atmospheric data from selected models (`model` created)
  call models_load(nml, model)

  ! Initialise the s3com structure
  call s3com_init(nml, model, s3com)
  npoints = s3com%npoints
  nlayers = s3com%nlayers

  ! Setup the RTTOV optics (`rttov_opt` created)
  call rttov_setup_opt(nml, zenangle, azangle, rttov_opt, s3com)

  ! Initialize RTTOV (loading surface data and optical properties)
  call rttov_init(rttov_opt, s3com)

  ! Loop over data points, by chunks
  nChunks = n_chunks(npoints, nml%npoints_it)
  loop_chunks : do iChunk = 1, nChunks

     write(6,*) ichunk, "/", nchunks

     ! Compute starting and end points for the chunk
     if (nChunks .EQ. 1) then
        idx_start = 1; idx_end = npoints
     else
        idx_start = (iChunk-1)* nml%npoints_it + 1
        idx_end = iChunk * nml%npoints_it
        if (idx_end .gt. npoints) idx_end = npoints
     end if

     ! Subset the atmosphere for RTTOV
     call rttov_setup_atm(idx_start, idx_end, model, rttov_atm)

     ! Subset the relevant part of the s3com structure
     call s3com_subset(idx_start, idx_end, s3com, s3com_chunk)

     if (nml%flag_retrievals) then

        ! call rttov_setup_atm(idx_start, idx_end, model, rttov_atm_oe)

        ! ! Extract the cloud position for ice and liquid phase and modify rttov_atm
        ! call init_zcloud(rttov_atm_oe,atm_oe)

        ! ! Set conditions to use RTTOV (now: need an ice cloud)
        ! idx_iwp = idx_ice(atm_oe,.TRUE.)
        ! atm_oe%flag_rttov(idx_iwp) = .TRUE.

        ! idx_oe = idx_rttov(atm_oe)
        ! IF (SIZE(idx_oe) .EQ. 0) THEN
        !    CYCLE
        ! ENDIF

        ! ! Update the cloud profiles in rttov_atm_oe using parameters from atm_oe
        ! call init_cloudprof(rttov_atm_oe, atm_oe)

        ! ! Run rttov on the selected atmosphere
        ! call run_rttov(rttov_atm_oe, rttov_opt, atm_oe, dealloc = .FALSE.)

     else

        s3com_chunk%flag_rttov(:) = .TRUE.
        call run_rttov(rttov_atm, rttov_opt, s3com_chunk, dealloc = .FALSE.)

     end if

     call s3com_update(s3com, s3com_chunk)

  end do loop_chunnks

  ! Write output file
  call write_output(s3com, model, nml)

  ! write(*,*)
  ! write(*,*) s3com%rad%wavelength(1), s3com%rad%wavelength(5)
  ! write(*,*) abs(0.741737 - s3com%rad%f_ref_total(1,1)), abs(0.2682695 - s3com%rad%f_ref_total(1,5))

  ! Deallocate arrays
  call models_deinit(model)

end program s3com_main

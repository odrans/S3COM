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

PROGRAM S3COM_MAIN

  USE s3com_types,         ONLY: wp, type_rttov_atm, type_rttov_opt, type_nml, type_model, type_s3com_new, type_s3com_new_ss
  USE mod_io_namelist,     ONLY: namelist_load
  USE mod_rttov_interface, ONLY: rttov_init
  USE mod_rttov_setup,     ONLY: rttov_setup_opt, rttov_setup_atm
  USE mod_rttov,           ONLY: run_rttov
  USE mod_atm_init,        ONLY: atm_init, atm_update, s3com_init, s3com_subset
  USE mod_models,          ONLY: models_load
  USE mod_utils_math,      ONLY: n_chunks
!  USE mod_model_cloud,     ONLY: init_zcloud, init_cloudprof
  USE mod_write_output,    ONLY: write_output
!  USE mod_oe_utils,        ONLY: idx_ice
  USE mod_rttov_utils,     ONLY: idx_rttov
!  USE mod_oe_run,          ONLY: oe_run

  IMPLICIT NONE

  TYPE(type_model), TARGET      :: model
  TYPE(type_rttov_atm)          :: rttov_atm
  TYPE(type_rttov_opt)          :: rttov_opt
  TYPE(type_s3com_new), TARGET  :: s3com
  TYPE(type_nml)                :: nml
  TYPE(type_s3com_new_ss)       :: oe

  REAL(KIND=wp) :: zenangle, azangle

  INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: idx_iwp, idx_oe
  INTEGER(KIND=4) :: nchunks, idx_start, idx_end, ichunk
  INTEGER(KIND=4) :: nlayers, npoints, npoints_it

  LOGICAL :: flag_oe

  ! Read namelist file (`nml` created)
  CALL namelist_load(nml)

  ! Temporary: setting the viewing angles
  zenangle = 0._wp; azangle = 0._wp       !Viewing satellite angles

  ! Load atmospheric data from selected models (`model` created)
  CALL models_load(nml, model)
  npoints = model%npoints
  nlayers = model%nlayers

  CALL s3com_init(nml, model, s3com)

  ! Setup the RTTOV optics (`rttov_opt` created)
  CALL rttov_setup_opt(model, nml, zenangle, azangle, rttov_opt)

  ! Initialize RTTOV (loading surface data and optical properties)
  CALL rttov_init(rttov_opt, nml, s3com)

  ! Loop over data points, by chunks
  nChunks = n_chunks(npoints, nml%npoints_it)
  DO iChunk = 1, nChunks

     WRITE(6,*) ichunk, "/", nchunks

     ! Compute starting and end points for the chunk
     IF (nChunks .EQ. 1) THEN
        idx_start = 1; idx_end = npoints
     ELSE
        idx_start = (iChunk-1)*npoints_it+1; idx_end = iChunk*npoints_it
        IF (idx_end .GT. npoints) idx_end=npoints
     END IF

     ! Subset the atmosphere (`rttov_atm` created as pointer to `model`)
     CALL rttov_setup_atm(idx_start, idx_end, model, rttov_atm)

     ! Point to the corresponding parts of the s3com structure
     CALL s3com_subset(s3com, rttov_atm, oe)

     IF (nml%flag_retrievals) THEN

        ! CALL rttov_setup_atm(idx_start, idx_end, model, rttov_atm_oe)

        ! ! Extract the cloud position for ice and liquid phase and modify rttov_atm
        ! CALL init_zcloud(rttov_atm_oe,atm_oe)

        ! ! Set conditions to use RTTOV (now: need an ice cloud)
        ! idx_iwp = idx_ice(atm_oe,.TRUE.)
        ! atm_oe%flag_rttov(idx_iwp) = .TRUE.

        ! idx_oe = idx_rttov(atm_oe)
        ! IF (SIZE(idx_oe) .EQ. 0) THEN
        !    CYCLE
        ! ENDIF

        ! ! Update the cloud profiles in rttov_atm_oe using parameters from atm_oe
        ! CALL init_cloudprof(rttov_atm_oe, atm_oe)

        ! ! Run rttov on the selected atmosphere
        ! CALL run_rttov(rttov_atm_oe, rttov_opt, atm_oe, dealloc = .FALSE.)

     ELSE

        oe%flag_rttov(:) = .TRUE.
        CALL run_rttov(rttov_atm, rttov_opt, oe, dealloc = .FALSE.)

     ENDIF

  ENDDO

  ! Write output file
  CALL write_output(s3com, model, nml)

END PROGRAM S3COM_MAIN

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

PROGRAM S3COM

  USE s3com_types,         ONLY: wp, type_rttov_atm, type_rttov_opt, type_s3com, type_nml, type_model
  USE mod_io_namelist,     ONLY: namelist_load
  USE mod_rttov_interface, ONLY: rttov_init
  USE mod_rttov_setup,     ONLY: rttov_setup_opt, rttov_setup_atm
  USE mod_rttov,           ONLY: run_rttov
  USE mod_atm_init,        ONLY: atm_init, atm_update
  USE mod_models,          ONLY: models_load
  USE mod_model_cloud,     ONLY: init_zcloud, init_cloudprof
  USE mod_write_output,    ONLY: write_output
  USE mod_oe_utils,        ONLY: idx_ice
  USE mod_rttov_utils,     ONLY: idx_rttov
  USE mod_oe_run,          ONLY: oe_run

  IMPLICIT NONE

  TYPE(type_model), TARGET         :: model
  TYPE(type_rttov_atm)    :: rttov_atm, rttov_atm_oe
  TYPE(type_rttov_opt)    :: rttov_opt
  TYPE(type_s3com)        :: atm_out, atm_oe, atm_oe_ip1, oe_tmp
  TYPE(type_nml)          :: nml

  CHARACTER(LEN = 32) :: fname_nml

  REAL(KIND=wp) :: zenangle, azangle

  INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: idx_iwp, idx_oe
  INTEGER(KIND=4) :: i, loc, j
  INTEGER(KIND=4) :: nchunks, idx_start, idx_end, nPtsPerIt, ichunk
  INTEGER(KIND=4) :: nidx
  INTEGER(KIND=4) :: Nlevels, Npoints, npoints_it

  LOGICAL :: flag_oe, dealloc_rttov

  ! Load namelist
  CALL namelist_load(nml)

  ! Temporary: setting the viewing angles
  zenangle = 0._wp; azangle = 0._wp       !Viewing satellite angles

  ! Load selected model inputs
  CALL models_load(nml%fname_in, model)
  npoints = model%npoints
  nlevels = model%nlevels

  ! Setup the RTTOV optics
  CALL rttov_setup_opt(model, zenangle, azangle, rttov_opt, nml)

  ! Initialize RTTOV (load data)
  CALL rttov_init(rttov_opt, nml)

  ! Initialize `atm_out`: the full atmospheric model to be outputed by S3COM
  flag_oe = .FALSE.
  CALL atm_init(1, npoints, nlevels, atm_out, flag_oe, nml)

  npoints_it = nml%npoints_it
  nChunks = nPoints/nPoints_it
  IF (MOD(npoints,npoints_it)/=0) nchunks = nchunks + 1
  IF (nPoints .EQ. nPoints_it) nChunks = 1

  flag_oe = .TRUE.

  DO iChunk = 1, nChunks

     WRITE(6,*) ichunk, "/", nchunks

     IF (nChunks .EQ. 1) THEN
        idx_start = 1; idx_end = nPoints
     ELSE
        idx_start = (iChunk-1)*nPoints_it+1; idx_end = iChunk*nPoints_it
        IF (idx_end .GT. nPoints) idx_end=nPoints
     END IF

     ! Subset the atmosphere based on the model inputs
     CALL rttov_setup_atm(idx_start, idx_end, model, rttov_atm)

     ! Initialize `atm_oe`: a subset atmospheric model used within the optimal estimation framework
     CALL atm_init(rttov_atm%idx_start, rttov_atm%idx_end, nlevels, atm_oe, flag_oe, nml)

     IF (nml%flag_retrievals) THEN

        CALL rttov_setup_atm(idx_start, idx_end, model, rttov_atm_oe)

        ! Extract the cloud position for ice and liquid phase and modify rttov_atm
        CALL init_zcloud(rttov_atm_oe,atm_oe)

        ! Set conditions to use RTTOV (now: need an ice cloud)
        idx_iwp = idx_ice(atm_oe,.TRUE.)
        atm_oe%flag_rttov(idx_iwp) = .TRUE.

        idx_oe = idx_rttov(atm_oe)
        IF (SIZE(idx_oe) .EQ. 0) THEN
           CYCLE
        ENDIF

        ! Update the cloud profiles in rttov_atm_oe using parameters from atm_oe
        CALL init_cloudprof(rttov_atm_oe, atm_oe)

        ! Run rttov on the selected atmosphere
        CALL run_rttov(rttov_atm_oe, rttov_opt, atm_oe, dealloc = .FALSE.)

     ELSE

        atm_oe%flag_rttov(:) = .TRUE.
        CALL run_rttov(rttov_atm, rttov_opt, atm_oe, dealloc = .FALSE.)

     ENDIF

     atm_oe%y_refl_total = atm_oe%f_refl_total; atm_oe%y_refl_clear = atm_oe%f_refl_clear
     atm_oe%y_bt_total = atm_oe%f_bt_total; atm_oe%y_bt_clear = atm_oe%f_bt_clear
     atm_oe%y_rad_total = atm_oe%f_rad_total; atm_oe%y_rad_clear = atm_oe%f_rad_clear; atm_oe%y_rad_cloudy = atm_oe%f_rad_cloudy

     IF(nml%flag_output_atm) THEN
        atm_oe%t = rttov_atm%t
        atm_oe%z = rttov_atm%z
        atm_oe%clc = rttov_atm%clc
        atm_oe%cdnc = rttov_atm%cdnc
        atm_oe%lwc = rttov_atm%lwc
        atm_oe%reff = rttov_atm%reff
     END IF

     IF (nml%flag_retrievals) THEN
        atm_oe%iwp_model = atm_oe%iwp

        CALL oe_run(atm_oe, rttov_atm_oe, rttov_opt)
     ENDIF

     ! Update the final atmosperic model
     CALL atm_update(idx_start, idx_end, atm_oe, atm_out)

  ENDDO

  ! Write output file
  ! CALL write_output(model, atm_out, nml, atm_out)

END PROGRAM S3COM

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

  USE s3com_types,         ONLY: wp, type_rttov_atm, type_rttov_opt, type_icon, type_s3com, type_nml
  USE mod_read_icon,       ONLY: read_icon, construct_icon, map_point_to_ll, extract_coordinates
  USE mod_io_namelist,     ONLY: read_namelist
  USE mod_rttov_interface, ONLY: rttov_init
  USE mod_rttov_setup,     ONLY: rttov_setup_opt, rttov_setup_atm
  USE mod_rttov,           ONLY: run_rttov
  USE mod_setup_atm,       ONLY: setup_atm, update_atm
  USE mod_model_cloud,     ONLY: init_zcloud, init_cloudprof
  USE mod_write_output,    ONLY: write_output
  USE mod_oe_utils,        ONLY: idx_ice
  USE mod_rttov_utils,     ONLY: idx_rttov
  USE mod_oe_run,          ONLY: oe_run

  IMPLICIT NONE

  TYPE(type_icon), TARGET :: icon
  TYPE(type_rttov_atm)    :: rttov_atm, rttov_atm_oe
  TYPE(type_rttov_opt)    :: rttov_opt
  TYPE(type_s3com)        :: atm, oe, oe_ip1, oe_tmp
  TYPE(type_nml)          :: nml

  REAL(KIND=wp) :: zenangle, azangle, sunzenangle, sunazangle

  INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: idx_iwp, idx_oe
  INTEGER(KIND=4) :: i, loc, j
  INTEGER(KIND=4) :: nchunks, idx_start, idx_end, nPtsPerIt, ichunk
  INTEGER(KIND=4) :: nidx
  INTEGER(KIND=4) :: Nlevels, Npoints, npoints_it

  LOGICAL :: flag_oe, dealloc_rttov

  ! Read from file.
  call read_namelist('config.nml', nml)
  npoints_it = nml%npoints_it

  ! Extract the number of vertical levels and grid points in the input files
  call extract_coordinates(nml%fname_in, Nlevels, Npoints)

  ! Temporary: setting the viewing and solar angles
  zenangle = 0._wp; azangle = 0._wp       !Viewing satellite angles
  sunzenangle = 0._wp; sunazangle = 0._wp !Viewing solar angles

  !zenangle = 30._wp; azangle = 0._wp
  !sunzenangle = 37._wp; sunazangle = 0._wp

  ! Construct the icon pointer
  CALL construct_icon(npoints, nlevels, icon)

  ! Read input netcdf file containing ICON outputs
  CALL read_icon(nml%fname_in, icon)

  ! Setup the RTTOV optics
  CALL rttov_setup_opt(zenangle, azangle, sunzenangle, sunazangle, nml%month, rttov_opt)

  ! Initialize RTTOV (load data)
  CALL rttov_init(rttov_opt)

  ! Setup the overall atmospheric model used for radiative transfer
  flag_oe = .FALSE.
  CALL setup_atm(1, icon%npoints, icon%nlevels, atm, flag_oe)

  nChunks = icon%nPoints/nPoints_it
  IF (MOD(icon%npoints,npoints_it)/=0) nchunks = nchunks + 1
  IF (icon%nPoints .EQ. nPoints_it) nChunks = 1

  flag_oe = .TRUE.

  DO iChunk = 1, nChunks

     WRITE(6,*) ichunk, "/", nchunks

     IF (nChunks .EQ. 1) THEN
        idx_start = 1; idx_end = nPoints
     ELSE
        idx_start = (iChunk-1)*nPoints_it+1; idx_end = iChunk*nPoints_it
        IF (idx_end .GT. nPoints) idx_end=nPoints
     END IF

     ! Subset the atmosphere based on the ICON input file
     CALL rttov_setup_atm(idx_start, idx_end, icon, rttov_atm)
     CALL rttov_setup_atm(idx_start, idx_end, icon, rttov_atm_oe)

     ! Setup the oe variables
     CALL setup_atm(rttov_atm%idx_start, rttov_atm%idx_end, icon%nlevels, oe, flag_oe)

     IF (nml%flag_retrievals) THEN

        ! Extract the cloud position for ice and liquid phase and modify rttov_atm
        CALL init_zcloud(rttov_atm_oe,oe)

        ! Set conditions to use RTTOV (now: need an ice cloud)
        idx_iwp = idx_ice(oe,.TRUE.)
        oe%flag_rttov(idx_iwp) = .TRUE.

        idx_oe = idx_rttov(oe)
        IF (SIZE(idx_oe) .EQ. 0) THEN
           CYCLE
        ENDIF

        ! Update the cloud profiles in rttov_atm_oe using parameters from oe
        CALL init_cloudprof(rttov_atm_oe, oe)

        ! Run rttov on the selected atmosphere
        CALL run_rttov(rttov_atm_oe, rttov_opt, oe, dealloc = .FALSE.)

     ELSE

        oe%flag_rttov(:) = .TRUE.
        CALL run_rttov(rttov_atm, rttov_opt, oe, dealloc = .FALSE.)

     ENDIF

     oe%y_refl_total = oe%f_refl_total; oe%y_refl_clear = oe%f_refl_clear
     oe%y_bt_total = oe%f_bt_total; oe%y_bt_clear = oe%f_bt_clear
     oe%y_rad_total = oe%f_rad_total; oe%y_rad_clear = oe%f_rad_clear; oe%y_rad_cloudy = oe%f_rad_cloudy
     oe%iwp_model = oe%iwp

     IF (nml%flag_retrievals) THEN
        CALL oe_run(oe, rttov_atm_oe, rttov_opt)
     ENDIF

     ! Update the final atmosperic model
     CALL update_atm(idx_start, idx_end, oe, atm)

  ENDDO

  ! Write output file
  CALL write_output(nml%fname_out, icon, atm)

END PROGRAM S3COM

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

MODULE mod_model_cloud
   
   USE s3com_types,  ONLY: type_rttov_atm, type_s3com
   USE s3com_config, ONLY: wp, iwp_lay_threshold, lwp_lay_threshold
   
   IMPLICIT NONE
   
   CONTAINS
      
      SUBROUTINE init_zcloud(rttov,oe)
         
         TYPE(type_rttov_atm), INTENT(IN) :: rttov
         TYPE(type_s3com), INTENT(INOUT)     :: oe
         
         INTEGER(KIND=4) :: ipoint, ilevel
         
         ! oe%ztop_ice(:) = 0._wp; oe%zbase_ice(:) = 0._wp
         ! oe%iwp (:) = 0._wp
         ! oe%ztop_ice_idx(:) = 0; oe%zbase_ice_idx(:) = 0

         ! DO ipoint = 1, rttov%npoints

         !    DO ilevel = 1, rttov%nlevels

         !       IF(rttov%iwc(ipoint,ilevel)*rttov%dz(ipoint,ilevel)*1E3 .GT. iwp_lay_threshold) THEN
         !          oe%ztop_ice_idx(ipoint) = ilevel
         !          EXIT
         !       ENDIF

         !    ENDDO

         !    DO ilevel = rttov%nlevels, 1, -1

         !       IF(rttov%iwc(ipoint,ilevel)*rttov%dz(ipoint,ilevel)*1E3 .GT. iwp_lay_threshold) then
         !          oe%zbase_ice_idx(ipoint) = ilevel
         !          EXIT
         !       ENDIF

         !    ENDDO

         !    IF(oe%ztop_ice_idx(ipoint) /= 0 .AND. oe%zbase_ice_idx(ipoint) /= 0) THEN
         !       oe%ztop_ice(ipoint)  = rttov%z(ipoint,oe%ztop_ice_idx(ipoint))
         !       oe%zbase_ice(ipoint) = rttov%z(ipoint,oe%zbase_ice_idx(ipoint) + 1) ! 1 is added to take the top altitude of the lower layer
         !       oe%iwp(ipoint) = SUM(rttov%iwc(ipoint,:) * rttov%dz(ipoint,:))
         !       !oe%lwp(ipoint) = SUM(rttov%lwc(ipoint,:) * rttov%dz(ipoint,:))
         !    ENDIF

         !    IF(oe%iwp(ipoint) .LT. 1E-4) THEN
         !       oe%iwp(ipoint) = 0._wp
         !    ENDIF

         ! ENDDO

      END SUBROUTINE init_zcloud
      
      SUBROUTINE init_cloudprof(rttov,oe)
         
         USE mod_rttov_utils, ONLY: idx_rttov
         
         TYPE(type_rttov_atm),  INTENT(INOUT) :: rttov
         TYPE(type_s3com), TARGET, INTENT(INOUT) :: oe
         
         INTEGER(KIND=4) :: ipoint, ilevel, idx
         INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: idx_oe
         
         ! idx_oe = idx_rttov(oe)

         ! oe%cla = rttov%tca

         ! DO ipoint = 1, size(idx_oe)

         !    idx = idx_oe(ipoint)
         !    oe%iwc(idx,:) = 0._wp

         !    IF(oe%iwp(idx) .GT. 1E-4) THEN
         !       oe%iwc(idx,oe%ztop_ice_idx(idx):oe%zbase_ice_idx(idx)) = oe%iwp(idx)/&
         !       &SUM(rttov%dz(idx,oe%ztop_ice_idx(idx):oe%zbase_ice_idx(idx)))
         !       oe%cla(idx,oe%ztop_ice_idx(idx):oe%zbase_ice_idx(idx)) = 1.0_wp
         !    ENDIF

         ! ENDDO

         ! rttov%iwc => oe%iwc
         ! rttov%tca => oe%cla
         
      END SUBROUTINE init_cloudprof
      
END MODULE mod_model_cloud

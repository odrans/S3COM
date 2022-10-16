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

MODULE mod_rttov_utils
   
   USE s3com_types, ONLY: wp, type_s3com_new_ss
   
   IMPLICIT NONE
   
   CONTAINS
   
      FUNCTION idx_rttov(oe)
      
         TYPE(type_s3com_new_ss), INTENT(IN) :: oe
         INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: idx_rttov
         INTEGER(KIND=4), DIMENSION(oe%npoints) :: idx_all
         INTEGER(KIND=4) :: idx, i
      
         idx_all = 0
         idx = 1
      
         DO i = 1, oe%npoints
            IF(oe%flag_rttov(i)) THEN
               idx_all(idx) = i
               idx = idx + 1
            END IF
         END DO
      
         idx = idx-1
      
         ALLOCATE(idx_rttov(idx)); idx_rttov(1:idx) = idx_all(1:idx)
      
         RETURN
      
      END FUNCTION idx_rttov
      
END MODULE mod_rttov_utils

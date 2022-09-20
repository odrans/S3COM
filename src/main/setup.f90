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

MODULE mod_setup_atm
   
   USE s3com_types,  ONLY: type_s3com, wp, type_nml
   USE s3com_config, ONLY: nstates, apriori_iwp
   
   IMPLICIT NONE
   
   CONTAINS
      
      SUBROUTINE update_atm(idx_start,idx_end,oe,s3com)
         
         ! Input variables
         INTEGER(kind=4), INTENT(IN) :: idx_start, idx_end
         TYPE(type_s3com), INTENT(IN) :: oe
         
         ! Output variables
         TYPE(type_s3com), INTENT(INOUT) :: s3com
         
         s3com%y_refl_total(idx_start:idx_end,:) = oe%y_refl_total(:,:)
         s3com%y_refl_clear(idx_start:idx_end,:) = oe%y_refl_clear(:,:)
         s3com%y_bt_total(idx_start:idx_end,:)   = oe%y_bt_total(:,:)
         s3com%y_bt_clear(idx_start:idx_end,:)   = oe%y_bt_clear(:,:)
         s3com%y_rad_total(idx_start:idx_end,:)  = oe%y_rad_total(:,:)
         s3com%y_rad_clear(idx_start:idx_end,:)  = oe%y_rad_clear(:,:)
         s3com%y_rad_cloudy(idx_start:idx_end,:) = oe%y_rad_cloudy(:,:)
         s3com%brdf(idx_start:idx_end,:)         = oe%brdf(:,:)
         s3com%emissivity(idx_start:idx_end,:)   = oe%emissivity(:,:)
         s3com%Xip1(idx_start:idx_end,:)         = oe%Xip1(:,:)
         s3com%iwp_model(idx_start:idx_end)      = oe%iwp_model(:)
         s3com%Gip1(idx_start:idx_end)           = oe%Gip1(:)
         
      END SUBROUTINE update_atm
      
      SUBROUTINE setup_atm(idx_start, idx_end, nlevels, y, flag_oe, nml)

        INTEGER(kind=4) :: nchannels

        ! Input variables
        INTEGER(kind=4), INTENT(IN) :: nlevels, idx_start, idx_end
        LOGICAL :: flag_oe

        ! Output variables
        TYPE(type_s3com), INTENT(OUT)   :: y
        TYPE(type_nml), INTENT(IN)      :: nml

        ! Internal variables
        INTEGER(KIND=4) :: npoints, ipoint

         nchannels = nml%nchannels

         npoints = idx_end - idx_start + 1
         
         y%npoints = npoints
         y%nstates = nstates
         y%nmeas = nml%nchannels

         ALLOCATE(y%y_rad_total(npoints,nchannels)); y%y_rad_total = 0._wp
         ALLOCATE(y%f_rad_total(npoints,nchannels)); y%f_rad_total = 0._wp
         ALLOCATE(y%y_rad_clear(npoints,nchannels)); y%y_rad_clear = 0._wp
         ALLOCATE(y%f_rad_clear(npoints,nchannels)); y%f_rad_clear = 0._wp
         ALLOCATE(y%y_rad_cloudy(npoints,nchannels)); y%y_rad_cloudy = 0._wp
         ALLOCATE(y%f_rad_cloudy(npoints,nchannels)); y%f_rad_cloudy = 0._wp
         ALLOCATE(y%y_refl_total(npoints,nchannels)); y%y_refl_total = 0._wp
         ALLOCATE(y%f_refl_total(npoints,nchannels)); y%f_refl_total = 0._wp
         ALLOCATE(y%y_refl_clear(npoints,nchannels)); y%y_refl_clear = 0._wp
         ALLOCATE(y%f_refl_clear(npoints,nchannels)); y%f_refl_clear = 0._wp
         ALLOCATE(y%y_bt_total(npoints,nchannels)); y%y_bt_total = 0._wp
         ALLOCATE(y%f_bt_total(npoints,nchannels)); y%f_bt_total = 0._wp
         ALLOCATE(y%y_bt_clear(npoints,nchannels)); y%y_bt_clear = 0._wp
         ALLOCATE(y%f_bt_clear(npoints,nchannels)); y%f_bt_clear = 0._wp
         ALLOCATE(y%Xa(npoints,nstates)); y%Xa = 0._wp
         ALLOCATE(y%g(npoints)); y%g = 0._wp
         ALLOCATE(y%g_meas(npoints)); y%g_meas = 0._wp
         ALLOCATE(y%ztop_ice(npoints)); y%ztop_ice = 0._wp
         ALLOCATE(y%zbase_ice(npoints)); y%zbase_ice = 0._wp
         ALLOCATE(y%ztop_ice_idx(npoints)); y%ztop_ice_idx = 0
         ALLOCATE(y%zbase_ice_idx(npoints)); y%zbase_ice_idx = 0
         ALLOCATE(y%iwp(npoints)); y%iwp = 0._wp
         ALLOCATE(y%brdf(npoints,nchannels)); y%brdf = 0._wp
         ALLOCATE(y%emissivity(npoints,nchannels)); y%emissivity = 0._wp
         ALLOCATE(y%iwp_model(npoints)); y%iwp_model = 0._wp
         ALLOCATE(y%xip1(npoints,nstates)); y%xip1 = 0._wp
         ALLOCATE(y%gip1(npoints)); y%gip1 = 0._wp
         
         IF(flag_oe) THEN
            ALLOCATE(y%iwc(npoints,nlevels)); y%iwc = 0._wp
            ALLOCATE(y%cla(npoints,nlevels)); y%cla = 0._wp
            ALLOCATE(y%xi(npoints,nstates)); y%xi = 0._wp
            ALLOCATE(y%gi(npoints)); y%gi = 0._wp
            ALLOCATE(y%stepsize(npoints)); y%stepsize = 0._wp
            ALLOCATE(y%i_stepsize(npoints)); y%i_stepsize = 0
            ALLOCATE(y%n_iter(npoints)); y%n_iter = 0
            ALLOCATE(y%k(npoints,nstates,nchannels)); y%k = 0._wp
            ALLOCATE(y%kt(npoints,nstates,nchannels)); y%kt = 0._wp
            ALLOCATE(y%sy(npoints,nchannels,nchannels)); y%sy = 0._wp
            ALLOCATE(y%sf(npoints,nchannels,nchannels)); y%sf = 0._wp
            ALLOCATE(y%se(npoints,nchannels,nchannels)); y%se = 0._wp
            ALLOCATE(y%sx(npoints,nstates,nstates)); y%sx = 0._wp
            ALLOCATE(y%sa(npoints,nstates,nstates)); y%sa = 0._wp
            ALLOCATE(y%se_i(npoints,nchannels,nchannels)); y%se_i = 0._wp
            ALLOCATE(y%sx_i(npoints,nstates,nstates)); y%sx_i = 0._wp
            ALLOCATE(y%sa_i(npoints,nstates,nstates)); y%sa_i = 0._wp
            ALLOCATE(y%flag_rttov(npoints)); y%flag_rttov = .FALSE.
            ALLOCATE(y%flag_testconv(npoints)); y%flag_testconv = .FALSE.
         END IF
         
         y%Xa(:,1) = apriori_iwp
         
      END SUBROUTINE setup_atm
      
END MODULE mod_setup_atm

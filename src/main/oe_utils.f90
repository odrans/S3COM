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

MODULE MOD_OE_UTILS
   
   USE s3com_types,        ONLY: wp, type_rttov_atm, type_rttov_opt, type_s3com
   USE mod_model_cloud, ONLY: init_cloudprof
   USE mod_rttov,       ONLY: run_rttov
   USE s3com_config,       ONLY: nstates, apriori_iwp_error
   USE mod_utils_math,  ONLY: inverse
   
   IMPLICIT NONE
   
   CONTAINS
      
      FUNCTION idx_ice(oe,flag_ice)
      
         TYPE(type_s3com), INTENT(IN) :: oe
         REAL(KIND=wp), DIMENSION(:), ALLOCATABLE :: idx_ice
         REAL(KIND=wp), DIMENSION(oe%npoints) :: idx_all
         INTEGER(KIND=4) :: idx, i
         LOGICAL, INTENT(IN) :: flag_ice
         LOGICAL :: test
         
         idx_all = 0
         idx=1
         
         DO i = 1, oe%npoints
            
            IF(flag_ice) THEN
               test = oe%iwp(i) .GT. 1E-4
            ELSE
               test = oe%iwp(i) .LT. 1E-4
            ENDIF
            
            IF(test) THEN
               idx_all(idx) = i
               idx = idx + 1
            ENDIF
            
         ENDDO
         
         idx=idx-1
         
         ALLOCATE(idx_ice(idx)); idx_ice = idx_all(1:idx)
         
         RETURN
         
      END FUNCTION idx_ice
      
      REAL(KIND=wp) FUNCTION oe_cost(oe,ipoint)
         
         TYPE(type_s3com), INTENT(IN) :: oe
         INTEGER(KIND=4), INTENT(IN) :: ipoint
         
         oe_cost = oe_cost_meas(oe,ipoint) + oe_cost_state(oe,ipoint)
         
         RETURN
         
      END FUNCTION oe_cost
      
      REAL(KIND=wp) FUNCTION oe_cost_meas(oe,ipoint)
         
         TYPE(type_s3com), INTENT(IN) :: oe
         INTEGER(KIND=4), INTENT(IN) :: ipoint
         
         oe_cost_meas = DOT_PRODUCT(oe%Y(ipoint,:) - oe%F(ipoint,:), MATMUL(oe%Se_i(ipoint,:,:),oe%Y(ipoint,:) - oe%F(ipoint,:)))
         
         RETURN
         
      END FUNCTION oe_cost_meas
      
      REAL(KIND=wp) FUNCTION oe_cost_state(oe,ipoint)
         
         TYPE(type_s3com), INTENT(IN) :: oe
         INTEGER(KIND=4), INTENT(IN) :: ipoint
         
         oe_cost_state = DOT_PRODUCT(oe%Xi(ipoint,:) - oe%Xa(ipoint,:),MATMUL(oe%Sa_i(ipoint,:,:),oe%Xi(ipoint,:) - oe%Xa(ipoint,:)))
         
         RETURN
         
      END FUNCTION oe_cost_state
      
      FUNCTION oe_Sx(oe,ipoint)
         
         REAL(KIND=wp), DIMENSION(nstates,nstates) :: oe_Sx, Sx_i
         TYPE(type_s3com), INTENT(IN) :: oe
         INTEGER(KIND=4), INTENT(IN) :: ipoint
         
         Sx_i(:,:) = (1.0_wp + oe%stepsize(ipoint)) * oe%Sa_i(ipoint,:,:) + MATMUL(oe%Kt(ipoint,:,:),MATMUL(oe%Se_i(ipoint,:,:),oe%K(ipoint,:,:)))
         CALL inverse(Sx_i, oe_Sx, nstates)
    
    RETURN
  END FUNCTION oe_Sx

  FUNCTION oe_Xip1(oe,ipoint)

    REAL(KIND=wp), DIMENSION(nstates) :: oe_Xip1
    
    TYPE(type_s3com), INTENT(IN) :: oe
    INTEGER(KIND=4), INTENT(IN) :: ipoint
    
    oe_Xip1(:) = oe%Xi(ipoint,:) + MATMUL(oe%Sx(ipoint,:,:),(MATMUL(oe%Kt(ipoint,:,:),MATMUL(oe%Se_i(ipoint,:,:),(oe%Y(ipoint,:)-oe%F(ipoint,:)))) &
         - MATMUL(oe%Sa_i(ipoint,:,:),(oe%Xi(ipoint,:)-oe%Xa(ipoint,:)))))
    
    RETURN
  END FUNCTION oe_Xip1
  

  

  SUBROUTINE get_jacobian(rttov_atm,rttov_opt,oe)

    USE mod_rttov_utils, only: idx_rttov

    ! Parameters
    REAL(kind=wp), PARAMETER :: diwp=1E-4

    ! Input variables
    TYPE(type_rttov_opt), INTENT(IN) :: rttov_opt
    TYPE(type_rttov_atm), INTENT(IN) :: rttov_atm
    TYPE(type_s3com), INTENT(INOUT) :: oe

    ! Local variables
    TYPE(type_s3com) :: oe_pert
    TYPE(type_rttov_atm) :: rttov_atm_pert

    REAL(kind=wp), DIMENSION(rttov_atm%nPoints,rttov_opt%nchannels) :: Fp1, Fm1
    INTEGER(KIND=4) :: ipoint
    INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: idx_oe
    
    oe_pert = oe; rttov_atm_pert = rttov_atm

    idx_oe = idx_rttov(oe)
    
    oe_pert%iwp(idx_oe) = oe%iwp(idx_oe) + diwp
    call init_cloudprof(rttov_atm_pert,oe_pert)
    call run_rttov(rttov_atm_pert,rttov_opt,oe_pert,dealloc=.FALSE.)
    Fp1 = oe_pert%F

    oe_pert%iwp(idx_oe) = oe%iwp(idx_oe) - diwp    
    call init_cloudprof(rttov_atm_pert,oe_pert)
    call run_rttov(rttov_atm_pert,rttov_opt,oe_pert,dealloc=.FALSE.)
    Fm1 = oe_pert%F

    DO ipoint=1,oe_pert%npoints
       if(oe%flag_rttov(ipoint)) oe%K(ipoint,1,:) = (Fp1(ipoint,:) - Fm1(ipoint,:))/(2.0_wp*diwp)
    END DO

  END SUBROUTINE get_jacobian


  SUBROUTINE get_covmat(oe)

    ! Input variables
    TYPE(type_s3com), INTENT(INOUT) :: oe

    INTEGER(KIND=4) :: &
         ichannel, istate

    DO ichannel=1, oe%nmeas
       oe%Sy(:,ichannel,ichannel) = (oe%y(:,ichannel) * 1E-2)**2
    END DO

    DO istate=1,nstates
       oe%Sa(:,istate,istate) = ( oe%Xa(:,istate) * apriori_iwp_error )**2
    END DO

    oe%Se = oe%Sy

  END SUBROUTINE get_covmat


  SUBROUTINE update_oemat(rttov_atm,rttov_opt,oe)

    ! Input variables
    TYPE(type_rttov_opt), INTENT(IN) :: rttov_opt
    TYPE(type_rttov_atm), INTENT(IN) :: rttov_atm
    TYPE(type_s3com), INTENT(INOUT) :: oe

    ! Local variables
    TYPE(type_rttov_atm) :: rttov_atm_new

    rttov_atm_new = rttov_atm

    call init_cloudprof(rttov_atm_new,oe)
    call run_rttov(rttov_atm_new,rttov_opt,oe,.FALSE.)
    call get_jacobian(rttov_atm_new,rttov_opt,oe)
    call get_covmat(oe)
    
  END SUBROUTINE update_oemat


  SUBROUTINE update_oe(oe_in,oe_ip1,nidx,idx_oe,rttov_atm,rttov_opt)

    ! Input variables
    INTEGER(KIND=4), INTENT(IN) :: nidx
    INTEGER(KIND=4), DIMENSION(nidx), INTENT(IN) :: idx_oe
    TYPE(type_rttov_opt), INTENT(IN) :: rttov_opt

    TYPE(type_s3com), INTENT(IN) :: oe_in
    TYPE(type_rttov_atm), INTENT(INOUT) :: rttov_atm
    
    ! Local variables
    TYPE(type_s3com), INTENT(OUT) :: oe_ip1
    TYPE(type_s3com) :: oe
    INTEGER(KIND=4) :: ipoint

    oe = oe_in
    
    ! First run of optimal estimation, update Xi
    DO ipoint=1,size(idx_oe)
       oe%Sx(idx_oe(ipoint),:,:) = oe_Sx(oe,idx_oe(ipoint))
       oe%Xip1(idx_oe(ipoint),:) = oe_Xip1(oe,idx_oe(ipoint))
    END DO

    ! Check if Xi is physically consistent, else change stepsize and update again
    DO ipoint=1,size(idx_oe)
       DO WHILE(oe%Xip1(idx_oe(ipoint),1).LT.1.0E-4)
          oe%stepsize(idx_oe(ipoint)) = oe%stepsize(idx_oe(ipoint)) * 2.0D0
          oe%i_stepsize(idx_oe(ipoint)) = oe%i_stepsize(idx_oe(ipoint)) + 1

          oe%Sx(idx_oe(ipoint),:,:) = oe_Sx(oe,idx_oe(ipoint))
          oe%Xip1(idx_oe(ipoint),:) = oe_Xip1(oe,idx_oe(ipoint))
       END DO
    END DO

    oe_ip1=oe
    oe_ip1%iwp(idx_oe) = oe%Xip1(idx_oe,1)
    call update_oemat(rttov_atm,rttov_opt,oe_ip1)

    DO ipoint=1,size(idx_oe)
       oe_ip1%Gip1(idx_oe(ipoint)) = oe_cost(oe_ip1,idx_oe(ipoint))
    END DO
    oe%Gip1 = oe_ip1%Gip1
    
  END SUBROUTINE update_oe
  


  ! SUBROUTINE copy_oe(oe_in,oe_out,nidx,idx_oe)

  !   TYPE(type_s3com), INTENT(IN) :: oe_in
  !   TYPE(type_s3com), INTENT(OUT) :: oe_out
  !   INTEGER(KIND=4), INTENT(IN) :: nidx
  !   INTEGER(KIND=4), DIMENSION(nidx), INTENT(IN) :: idx_oe

  !   TYPE(type_s3com) :: oe_tmp

  !   oe_out = oe_in
    
    

  !   s3com%y(npoints,nchannels),s3com%f(npoints,nchannels),s3com%y_clear(npoints,nchannels),s3com%f_clear(npoints,nchannels)
  !   s3com%x(npoints,nstates),s3com%xa(npoints,nstates), s3com%g(npoints), s3com%g_meas(npoints)
  !   s3com%ztop_ice(npoints),s3com%zbase_ice(npoints),s3com%ztop_ice_idx(npoints),s3com%zbase_ice_idx(npoints)
  !   s3com%ztop_liq(npoints),s3com%zbase_liq(npoints),s3com%ztop_liq_idx(npoints),s3com%zbase_liq_idx(npoints)
  !   s3com%iwp(npoints),s3com%lwp(npoints),s3com%brdf(npoints,nchannels)
  !   s3com%xi(npoints,nstates),s3com%xip1(npoints,nstates), s3com%gi(npoints),s3com%gip1(npoints),s3com%stepsize(npoints)
  !   s3com%iloop_stepsize(npoints), s3com%iloop_iter(npoints)
  !   s3com%k(npoints,nstates,nchannels), s3com%kt(npoints,nstates,nchannels)
  !   s3com%sy(npoints,nchannels,nchannels), s3com%sf(npoints,nchannels,nchannels)
  !   s3com%se(npoints,nchannels,nchannels), s3com%sx(npoints,nstates,nstates), s3com%sa(npoints,nstates,nstates)
  !   s3com%se_i(npoints,nchannels,nchannels), s3com%sx_i(npoints,nstates,nstates), s3com%sa_i(npoints,nstates,nstates)
  !   s3com%idx_orig(npoints),s3com%flag_rttov(npoints),s3com%flag_stepsize(npoints)
    
   
  ! END SUBROUTINE copy_oe
    

END MODULE MOD_OE_UTILS

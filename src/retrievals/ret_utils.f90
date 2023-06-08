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

module mod_ret_utils
   
   use s3com_types,         only: wp, type_model, type_rttov_opt, type_s3com, type_cld
   use s3com_config,        only: nstates, apriori_cod_error, apriori_reff_error
   use mod_ret_cloud_model, only: init_cloud_prof
   use mod_rttov_utils,     only: find_ret_idx_rttov
   use mod_rttov,           only: run_rttov
   use mod_utils_math,      only: inverse
   
   implicit none
   
   public
   
contains
   
   ! ----------------------------------------------------------------------------------------------------------------------------
   !> @brief Initialize retrieval data structure
   !> @details This subroutine initializes the retrieval data structure by allocating the required arrays. All arrays are 
   !> initialized to zero.
   !> @param[inout] s3com retrieval data structure
   subroutine ret_init(s3com)
      
      ! Inputs/outputs
      type(type_s3com), intent(inout) :: s3com
      
      ! Internal
      integer(kind=4) :: npoints, nlevels, nlayers, nmeas, nstates, nbparams
      
      s3com%ret%npoints  = s3com%npoints
      s3com%ret%nlevels  = s3com%nlevels
      s3com%ret%nlayers  = s3com%nlayers
      s3com%ret%nmeas    = s3com%nmeas
      s3com%ret%nstates  = s3com%nstates
      
      npoints  = s3com%ret%npoints
      nlevels  = s3com%ret%nlevels
      nlayers  = s3com%ret%nlayers
      nmeas    = s3com%ret%nmeas
      nstates  = s3com%ret%nstates
      nbparams = 3*(s3com%ret%nlevels + 1) ! number of non-retrieved parameters of the forward model
      
      ! Cloud properties
      allocate(s3com%ret%ztop_liq_idx(npoints), source = 0)
      allocate(s3com%ret%zbase_liq_idx(npoints), source = 0)
      allocate(s3com%ret%ztop_liq(npoints), source = 0._wp)
      allocate(s3com%ret%zbase_liq(npoints), source = 0._wp)
      allocate(s3com%ret%lwp(npoints), source = 0._wp)
      allocate(s3com%ret%cod(npoints), source = 0._wp)
      allocate(s3com%ret%lwc(npoints, nlayers), source = 0._wp)
      allocate(s3com%ret%clc(npoints, nlayers), source = 0._wp)
      allocate(s3com%ret%reff(npoints, nlayers), source = 0._wp)
      allocate(s3com%ret%cdnc(npoints, nlayers), source = 0._wp)
      
      ! Optimal estimation method
      !allocate(s3com%ret%Y(npoints, nmeas), source = 0._wp)
      !allocate(s3com%ret%F(npoints, nmeas), source = 0._wp)
      allocate(s3com%ret%n_iter(npoints), source = 0)
      allocate(s3com%ret%Ji(npoints), source = 0._wp)
      allocate(s3com%ret%Jip1(npoints), source = 0._wp)
      allocate(s3com%ret%Ji_meas(npoints), source = 0._wp)
      allocate(s3com%ret%Xa(npoints, nstates), source = 0._wp)
      allocate(s3com%ret%Xi(npoints, nstates), source = 0._wp)
      allocate(s3com%ret%Xip1(npoints, nstates), source = 0._wp)
      
      ! Error variance-covariance matrices
      allocate(s3com%ret%Se(npoints, nmeas, nmeas), source = 0._wp)
      allocate(s3com%ret%Se_inv(npoints, nmeas, nmeas), source = 0._wp)
      allocate(s3com%ret%Sy(npoints, nmeas, nmeas), source = 0._wp)
      allocate(s3com%ret%Sf(npoints, nmeas, nmeas), source = 0._wp)
      allocate(s3com%ret%Sa(npoints, nstates, nstates), source = 0._wp)
      allocate(s3com%ret%Sa_inv(npoints, nstates, nstates), source = 0._wp)
      allocate(s3com%ret%Sx(npoints, nstates, nstates), source = 0._wp)
      allocate(s3com%ret%Sx_inv(npoints, nstates, nstates), source = 0._wp)
      allocate(s3com%ret%K(npoints, nmeas, nstates), source = 0._wp)
      allocate(s3com%ret%Kt(npoints, nstates, nmeas), source = 0._wp)
      
      ! Non retrieved parameters of the forward model
      allocate(s3com%ret%error_t(npoints, nlevels), source = 0._wp)
      allocate(s3com%ret%error_q(npoints, nlevels), source = 0._wp)
      allocate(s3com%ret%error_clc(npoints, nlayers), source = 0._wp)
      allocate(s3com%ret%error_surf_brdf(npoints, nmeas), source = 0._wp)
      allocate(s3com%ret%error_surf_emiss(npoints, nmeas), source = 0._wp)
      allocate(s3com%ret%error_cld_top(npoints), source = 0._wp)
      allocate(s3com%ret%error_cld_base(npoints), source = 0._wp)
      
      allocate(s3com%ret%sigma_t(npoints, nlevels), source = 0._wp)
      allocate(s3com%ret%sigma_q(npoints, nlevels), source = 0._wp)
      allocate(s3com%ret%sigma_clc(npoints, nlayers), source = 0._wp)
      allocate(s3com%ret%sigma_surf_brdf(npoints), source = 0._wp)
      allocate(s3com%ret%sigma_surf_emiss(npoints), source = 0._wp)
      allocate(s3com%ret%sigma_cld_top(npoints), source = 0._wp)
      allocate(s3com%ret%sigma_cld_base(npoints), source = 0._wp)
      
      allocate(s3com%ret%Kb(npoints, nmeas, nbparams), source = 0._wp)
      allocate(s3com%ret%Kbt(npoints, nbparams, nmeas), source = 0._wp)
      allocate(s3com%ret%Sb(npoints, nbparams, nbparams), source = 0._wp)
      
      ! Flags
      allocate(s3com%ret%flag_rttov(npoints)); s3com%ret%flag_rttov = .true.
      allocate(s3com%ret%flag_conv_test(npoints)); s3com%ret%flag_conv_test = .true.
      
   end subroutine ret_init
   ! ----------------------------------------------------------------------------------------------------------------------------
   
   ! ----------------------------------------------------------------------------------------------------------------------------
   !> @brief Determination of the cloud phase (liquid or not) and accounting of the number of liquid clouds
   !> @param[in] s3com retrieval data structure
   !> @param[in] flag_liq
   function idx_liq(s3com, flag_liq)
      
      ! Inputs
      type(type_s3com), intent(in) :: s3com
      logical, intent(in) :: flag_liq
      
      ! Internal
      real(kind=wp), dimension(:), allocatable :: idx_liq
      real(kind=wp), dimension(s3com%ret%npoints) :: idx_all
      integer(kind=4) :: idx, i
      logical :: test
      
      idx_all = 0
      idx = 1
      
      do i = 1, s3com%ret%npoints
         
         if(flag_liq) then
            test = s3com%ret%lwp(i) .gt. 1e-4
         else
            test = s3com%ret%lwp(i) .lt. 1e-4
         endif
         
         if(test) then
            idx_all(idx) = i
            idx = idx + 1
         endif
         
      enddo
      
      idx = idx - 1
      
      allocate(idx_liq(idx)); idx_liq = idx_all(1:idx)
      
      return
      
   end function idx_liq
   ! ----------------------------------------------------------------------------------------------------------------------------
   
   ! ----------------------------------------------------------------------------------------------------------------------------
   !> @brief Cost function: sum of the contributions from a priori and from measurements
   !> @param[in] s3com retrieval data structure
   !> @param[in] ipoint
   real(kind=wp) function ret_J(s3com, ipoint)
      
      ! Inputs
      type(type_s3com), intent(in) :: s3com
      integer(kind=4), intent(in) :: ipoint
      
      ret_J = J_state(s3com,ipoint) + J_meas(s3com,ipoint)
      
      return
      
   end function ret_J
   ! ----------------------------------------------------------------------------------------------------------------------------
   
   ! ----------------------------------------------------------------------------------------------------------------------------
   !> @brief Cost function: contribution from a priori
   !> @param[in] s3com retrieval data structure
   !> @param[in] ipoint
   real(kind=wp) function J_state(s3com, ipoint)
      
      ! Inputs
      type(type_s3com), intent(in) :: s3com
      integer(kind=4), intent(in) :: ipoint
      
      J_state = dot_product(s3com%ret%Xi(ipoint,:)-s3com%ret%Xa(ipoint,:), &
                matmul(s3com%ret%Sa_inv(ipoint,:,:),s3com%ret%Xi(ipoint,:)-s3com%ret%Xa(ipoint,:)))
      
      return
      
   end function J_state
   ! ----------------------------------------------------------------------------------------------------------------------------
   
   ! ----------------------------------------------------------------------------------------------------------------------------
   !> @brief Cost function: contribution from measurements
   !> @param[in] s3com retrieval data structure
   !> @param[in] ipoint
   real(kind=wp) function J_meas(s3com, ipoint)
      
      ! Inputs
      type(type_s3com), intent(in) :: s3com
      integer(kind=4), intent(in) :: ipoint
      
      J_meas = dot_product(s3com%ret%Y(ipoint,:)-s3com%ret%F(ipoint,:), &
               matmul(s3com%ret%Se_inv(ipoint,:,:),s3com%ret%Y(ipoint,:)-s3com%ret%F(ipoint,:)))
      
      return
      
   end function J_meas
   ! ----------------------------------------------------------------------------------------------------------------------------
   
   ! ----------------------------------------------------------------------------------------------------------------------------
   !> @brief State vector best estimate error variance-covariance matrix
   !> @param[in] s3com retrieval data structure
   !> @param[in] ipoint
   function ret_Sx(s3com, ipoint)
      
      ! Inputs
      type(type_s3com), intent(in) :: s3com
      integer(kind=4), intent(in) :: ipoint
      
      ! Internal
      real(kind=wp), dimension(nstates,nstates) :: ret_Sx, Sx_inv
      
      Sx_inv(:,:) = s3com%ret%Sa_inv(ipoint,:,:) + &
                    matmul(s3com%ret%Kt(ipoint,:,:),matmul(s3com%ret%Se_inv(ipoint,:,:),s3com%ret%K(ipoint,:,:)))
      
      call inverse(Sx_inv, ret_Sx, nstates)
      
      return
      
   end function ret_Sx
   ! ----------------------------------------------------------------------------------------------------------------------------
   
   ! ----------------------------------------------------------------------------------------------------------------------------
   !> @brief State vector best estimate following the Gauss-Newton method
   !> @param[in] s3com retrieval data structure
   !> @param[in] ipoint
   function ret_Xip1(s3com, ipoint)
      
      ! Inputs
      type(type_s3com), intent(in) :: s3com
      integer(kind=4), intent(in) :: ipoint
      
      ! Internal
      real(kind=wp), dimension(nstates) :: ret_Xip1
      
      ret_Xip1(:) = s3com%ret%Xi(ipoint,:) + &
                    matmul(s3com%ret%Sx(ipoint,:,:),(matmul(s3com%ret%Kt(ipoint,:,:), &
                    matmul(s3com%ret%Se_inv(ipoint,:,:),(s3com%ret%Y(ipoint,:)-s3com%ret%F(ipoint,:)))) - &
                    matmul(s3com%ret%Sa_inv(ipoint,:,:),(s3com%ret%Xi(ipoint,:)-s3com%ret%Xa(ipoint,:)))))
      
      return
      
   end function ret_Xip1
   ! ----------------------------------------------------------------------------------------------------------------------------
   
   ! ----------------------------------------------------------------------------------------------------------------------------
   !> @brief Jacobian matrix K
   !> @details This subroutine calculates the sensitivity of each observation to each state to be retrieved.
   !> @param[inout] rttov_atm
   !> @param[in] rttov_opt
   !> @param[inout] s3com retrieval data structure
   subroutine ret_K(rttov_atm, rttov_opt, s3com, cld)
      
      ! Inputs
      type(type_rttov_opt), intent(in) :: rttov_opt
      type(type_cld), intent(in) :: cld
      
      ! Inputs/outputs
      type(type_model), intent(inout) :: rttov_atm
      type(type_s3com), intent(inout) :: s3com
      
      ! Internal
      type(type_model) :: rttov_atm_pert
      real(kind=wp), dimension(rttov_atm%npoints, rttov_opt%nchannels) :: Fp1, Fm1
      integer(kind=4), dimension(:), allocatable :: idx_ret
      integer(kind=4) :: ipoint, idx, imeas!, ilayer
      
      ! Local parameters
      real(kind=wp), parameter :: delta_pert = 0.01_wp ! input perturbation of 1 %
      
      rttov_atm_pert = rttov_atm
      
      idx_ret = find_ret_idx_rttov(s3com)
      
      ! Xi(1) = COD
      ! -----------
      do ipoint = 1, size(idx_ret)
         
         idx = idx_ret(ipoint)
         
         Fp1 = 0.0_wp; Fm1 = 0.0_wp
         
         rttov_atm_pert%cod(idx) = rttov_atm%cod(idx) * (1 + delta_pert)
         call run_rttov(rttov_atm_pert, rttov_opt, s3com, cld)
         Fp1 = s3com%rad%f_rad_total
         
         rttov_atm_pert%cod(idx) = rttov_atm%cod(idx) * (1 - delta_pert)
         call run_rttov(rttov_atm_pert, rttov_opt, s3com, cld)
         Fm1 = s3com%rad%f_rad_total
         
         do imeas = 1, rttov_opt%nchannels
            if (s3com%ret%flag_rttov(idx)) then
               s3com%ret%K(idx,imeas,1) = (Fp1(idx,imeas)-Fm1(idx,imeas))/(2.0_wp*rttov_atm%cod(idx)*delta_pert)
            endif
         enddo
         
      enddo
      ! -----------
      
      ! Xi(2) = reff
      ! ------------
!      do ipoint = 1, size(idx_ret)
!         
!         idx = idx_ret(ipoint)
!         
!         do ilayer = 1, rttov_atm_pert%nlayers
!            
!            Fp1 = 0.0_wp; Fm1 = 0.0_wp
!            
!            rttov_atm_pert%reff(idx,ilayer) = rttov_atm%reff(idx,ilayer) * (1 + delta_pert)
!            call run_rttov(rttov_atm_pert, rttov_opt, s3com, dealloc=.false.)
!            Fp1 = s3com%rad%f_rad_total
!            
!            rttov_atm_pert%reff(idx,ilayer) = rttov_atm%reff(idx,ilayer) * (1 - delta_pert)
!            call run_rttov(rttov_atm_pert, rttov_opt, s3com, dealloc=.false.)
!            Fm1 = s3com%rad%f_rad_total
!            
!            do imeas = 1, rttov_opt%nchannels
!               if (s3com%ret%flag_rttov(idx)) then
!                  if (rttov_atm%reff(idx,ilayer) .gt. 0.0_wp) then
!                     s3com%ret%K(idx,imeas,ilayer+1) = &
!                     (Fp1(idx,imeas)-Fm1(idx,imeas))/(2.0_wp*rttov_atm%reff(idx,ilayer)*delta_pert)
!                  endif
!               endif
!            enddo
!            
!         enddo
!         
!      enddo
      
      do ipoint = 1, size(idx_ret)
         
         idx = idx_ret(ipoint)
         
         Fp1 = 0.0_wp; Fm1 = 0.0_wp
         
         rttov_atm_pert%reff(idx,s3com%ret%ztop_liq_idx(idx):s3com%ret%zbase_liq_idx(idx)) = &
         rttov_atm%reff(idx,s3com%ret%ztop_liq_idx(idx):s3com%ret%zbase_liq_idx(idx)) + 0.01_wp
         
         call run_rttov(rttov_atm_pert, rttov_opt, s3com, cld)
         
         Fp1 = s3com%rad%f_rad_total
         
         rttov_atm_pert%reff(idx,s3com%ret%ztop_liq_idx(idx):s3com%ret%zbase_liq_idx(idx)) = &
         rttov_atm%reff(idx,s3com%ret%ztop_liq_idx(idx):s3com%ret%zbase_liq_idx(idx)) - 0.01_wp
         
         call run_rttov(rttov_atm_pert, rttov_opt, s3com, cld)
         
         Fm1 = s3com%rad%f_rad_total
         
         do imeas = 1, rttov_opt%nchannels
            if (s3com%ret%flag_rttov(idx)) then
               s3com%ret%K(idx,imeas,2) = (Fp1(idx,imeas)-Fm1(idx,imeas))/(2.0_wp*0.01_wp)
            endif
         enddo
         
      enddo
      ! ------------
      
   end subroutine ret_K
   ! ----------------------------------------------------------------------------------------------------------------------------
   
   ! ----------------------------------------------------------------------------------------------------------------------------
   !> @brief Forward model error variance-covariance matrix calculation
   !> @details This subroutine calculates the sensitivity of each observation to each non-retrieved forward model parameter 
   !> (Jacobian matrix Kb), the non-retrieved forward model parameters error variance-covariance matrix (Sb) and the forward 
   !> model error variance-covariance matrix (Sf).
   !> @param[inout] rttov_atm
   !> @param[in] rttov_opt
   !> @param[inout] s3com retrieval/jacobian/radiance data structure
   subroutine ret_Sf(rttov_atm, rttov_opt, s3com)
      
      ! Inputs
      type(type_rttov_opt), intent(in) :: rttov_opt
      
      ! Inputs/outputs
      type(type_model), intent(inout) :: rttov_atm
      type(type_s3com), intent(inout) :: s3com
      
      ! Internal
      integer(kind=4), dimension(:), allocatable :: idx_ret
      integer(kind=4) :: npoints, nlevels, nlayers, nmeas
      integer(kind=4) :: ipoint, idx, imeas, ilayer, ilevel
      integer(kind=4) :: idx_t, idx_q, idx_clc, idx_surf_alb, idx_surf_emiss, idx_cld_top, idx_cld_base
      
      npoints = rttov_atm%npoints
      nlevels = rttov_atm%nlevels
      nlayers = rttov_atm%nlayers
      nmeas   = rttov_opt%nchannels
      
      idx_ret = find_ret_idx_rttov(s3com)
      
      do ipoint = 1, size(idx_ret)
         
         idx = idx_ret(ipoint)
         
         do imeas = 1, rttov_opt%nchannels
            
            if (s3com%ret%flag_rttov(idx)) then
               
               ! Jacobian matrix Kb
               s3com%ret%Kb(idx,imeas,1:nlevels)               = s3com%jac%t(idx,imeas,1:nlevels)     ! temperature profile
               s3com%ret%Kb(idx,imeas,nlevels+1:2*nlevels)     = s3com%jac%q(idx,imeas,1:nlevels)     ! humidity profile
               s3com%ret%Kb(idx,imeas,2*nlevels+1:3*nlevels-1) = s3com%jac%cfrac(idx,imeas,1:nlayers) ! cloud fraction profile
               s3com%ret%Kb(idx,imeas,3*nlevels:3*nlevels)     = s3com%jac%brdf(idx,imeas)            ! surface brdf/albedo
               s3com%ret%Kb(idx,imeas,3*nlevels+1:3*nlevels+1) = s3com%jac%emiss(idx,imeas)           ! surface emissivity
               s3com%ret%Kb(idx,imeas,3*nlevels+2:3*nlevels+2) = 0.0_wp                               ! cloud top altitude
               s3com%ret%Kb(idx,imeas,3*nlevels+3:3*nlevels+3) = 0.0_wp                               ! cloud base altitude
               
            endif
            
         enddo
         
         s3com%ret%Kbt(idx,:,:) = transpose(s3com%ret%Kb(idx,:,:))
         
      enddo
      
      do ipoint = 1, size(idx_ret)
         
         idx = idx_ret(ipoint)
         
         if (s3com%ret%flag_rttov(idx)) then
            
            ! Attribution of errors/uncertainties
            s3com%ret%error_t(idx,1:nlevels) = 1.0_wp / rttov_atm%t(idx,1:nlevels) ! temperature: error of 1 K per level
            s3com%ret%error_q(idx,1:nlevels) = 0.1_wp                              ! humidity: error of 10 % per level
            s3com%ret%error_clc(idx,1:nlayers) = 0.01_wp                           ! cloud fraction: error of 1 % per layer
            s3com%ret%error_surf_brdf(idx,1:nmeas) = 0.01_wp                       ! surface brdf/albedo: error of 1 % per channel
            s3com%ret%error_surf_emiss(idx,1:nmeas) = 0.01_wp                      ! surface emissivity: error of 1 % per channel
            s3com%ret%error_cld_top(idx) = 100.0_wp / s3com%ret%ztop_liq(idx)      ! cloud top altitude: error of 100 m
            s3com%ret%error_cld_base(idx) = 100.0_wp / s3com%ret%zbase_liq(idx)    ! cloud base altitude: error of 100 m
            ! -------------------------------------------------------------------------------------------------------------------
            
            ! Variances
            ! -------------------------------------------------------------------------------------------------------------------
            ! Temperature profile
            s3com%ret%sigma_t(idx,1:nlevels) = rttov_atm%t(idx,1:nlevels)*s3com%ret%error_t(idx,1:nlevels)
            
            ! Humidity profile
            s3com%ret%sigma_q(idx,1:nlevels) = rttov_atm%q(idx,1:nlevels)*s3com%ret%error_q(idx,1:nlevels)
            
            ! Cloud fraction profile
            s3com%ret%sigma_clc(idx,1:nlayers) = rttov_atm%clc(idx,1:nlayers)*s3com%ret%error_clc(idx,1:nlayers)
            
            ! Surface brdf/albedo
            s3com%ret%sigma_surf_brdf(idx) = sum(s3com%rad%brdf(idx,1:nmeas)*s3com%ret%error_surf_brdf(idx,1:nmeas))
            
            ! Surface emissivity
            s3com%ret%sigma_surf_emiss(idx) = sum(s3com%rad%emiss(idx,1:nmeas)*s3com%ret%error_surf_emiss(idx,1:nmeas))
            
            !sigma_brdf(idx,1:rttov_opt%nchannels) = 0.0_wp; sigma_emiss(idx,1:rttov_opt%nchannels) = 0.0_wp
            
            !do imeas = 1, rttov_opt%nchannels
            !   sigma_brdf(idx) = sigma_brdf(idx) + s3com%rad%brdf(idx,imeas)*error_surf_brdf
            !   sigma_emiss(idx) = sigma_emiss(idx) + s3com%rad%emissivity(idx,imeas)*error_emiss
            !enddo
            
            ! Cloud top altitude
            s3com%ret%sigma_cld_top(idx) = s3com%ret%ztop_liq(idx)*s3com%ret%error_cld_top(idx)
            
            ! Cloud base altitude
            s3com%ret%sigma_cld_base(idx) = s3com%ret%zbase_liq(idx)*s3com%ret%error_cld_base(idx)
            ! -------------------------------------------------------------------------------------------------------------------
            
            ! Non-retrieved forward model parameters error variance-covariance matrix (Sb is diagonal)
            ! -------------------------------------------------------------------------------------------------------------------
            ! Temperature profile
            idx_t = 0
            do ilevel = 1, nlevels
               s3com%ret%Sb(idx,idx_t+ilevel,idx_t+ilevel) = s3com%ret%sigma_t(idx,ilevel)**2
            enddo
            
            ! Humidity profile
            idx_q = nlevels
            do ilevel = 1, nlevels
               s3com%ret%Sb(idx,idx_q+ilevel,idx_q+ilevel) = s3com%ret%sigma_q(idx,ilevel)**2
            enddo
            
            ! Cloud fraction profile
            idx_clc = 2*nlevels
            do ilayer = 1, nlayers
               s3com%ret%Sb(idx,idx_clc+ilayer,idx_clc+ilayer) = s3com%ret%sigma_clc(idx,ilayer)**2
            enddo
            
            ! Surface brdf/albedo
            idx_surf_alb = 3*nlevels
            s3com%ret%Sb(idx,idx_surf_alb,idx_surf_alb) = s3com%ret%sigma_surf_brdf(idx)**2
            
            ! Surface emissivity
            idx_surf_emiss = 3*nlevels+1
            s3com%ret%Sb(idx,idx_surf_emiss,idx_surf_emiss) = s3com%ret%sigma_surf_emiss(idx)**2
            
            ! Cloud top altitude
            idx_cld_top = 3*nlevels+2
            s3com%ret%Sb(idx,idx_cld_top,idx_cld_top) = s3com%ret%sigma_cld_top(idx)**2
            
            ! Cloud base altitude
            idx_cld_base = 3*nlevels+3
            s3com%ret%Sb(idx,idx_cld_base,idx_cld_base) = s3com%ret%sigma_cld_base(idx)**2
            ! -------------------------------------------------------------------------------------------------------------------
            
         endif
         
      enddo
      
      do ipoint = 1, size(idx_ret)
         
         idx = idx_ret(ipoint)
         
         if (s3com%ret%flag_rttov(idx)) then
            ! Forward model error variance-covariance matrix (Sf is non diagonal)
            s3com%ret%Sf(idx,:,:) = matmul(s3com%ret%Kb(idx,:,:),matmul(s3com%ret%Sb(idx,:,:),s3com%ret%Kbt(idx,:,:)))
         endif
         
      enddo
      write(*,*) "s3com%ret%Sf", s3com%ret%Sf
   end subroutine ret_Sf
   ! ----------------------------------------------------------------------------------------------------------------------------
   
   ! ----------------------------------------------------------------------------------------------------------------------------
   !> @brief Calculation of the error variance-covariance matrices
   !> @param[inout] s3com retrieval data structure
   subroutine ret_cov_mat(s3com)
      
      ! Inputs/outputs
      type(type_s3com), intent(inout) :: s3com
      
      ! Internal
      integer(kind=4) :: ipoint, ichannel!, istate
      
      do ipoint = 1, s3com%ret%npoints
         
         ! A priori error variance-covariance matrix (Sa is diagonal, non-diagonal terms are considered to be null)
         s3com%ret%Sa(ipoint,1,1) = (s3com%ret%Xa(ipoint,1)*apriori_cod_error)**2
         s3com%ret%Sa(ipoint,2,2) = (s3com%ret%Xa(ipoint,2)*apriori_reff_error)**2
         
         do ichannel = 1, s3com%ret%nmeas
            
            ! Measurement error variance-covariance matrix (Sy is diagonal, non-diagonal terms are considered to be null)
            ! Uncertainty of 5 % in MODIS radiances for the reflective solar bands (from 0.41 to 2.3 um)
            ! Xiong et al., 2018: Updates of Moderate Resolution Imaging Spectroradiometer on-orbit calibration uncertainty 
            ! assessments
            s3com%ret%Sy(ipoint,ichannel,ichannel) = (s3com%ret%Y(ipoint,ichannel)*0.05)**2
            
         enddo
         
      enddo
      
      ! Total error variance-covariance matrix in the measurement space (Se is non diagonal)
      s3com%ret%Se = s3com%ret%Sy + s3com%ret%Sf
      
   end subroutine ret_cov_mat
   ! ----------------------------------------------------------------------------------------------------------------------------
   
   ! ----------------------------------------------------------------------------------------------------------------------------
   !> @brief Update the error variance-covariance matrices
   !> @param[in] rttov_atm
   !> @param[in] rttov_opt
   !> @param[inout] s3com retrieval data structure
   subroutine ret_update_mat(rttov_atm, rttov_opt, s3com, cld)
      
      ! Inputs
      type(type_model), intent(in) :: rttov_atm
      type(type_rttov_opt), intent(in) :: rttov_opt
      type(type_cld), intent(in) :: cld
      
      ! Inputs/outputs
      type(type_s3com), intent(inout) :: s3com
      
      ! Internal
      type(type_model) :: rttov_atm_new
      
      rttov_atm_new = rttov_atm
      
      call init_cloud_prof(rttov_atm_new, s3com)
      call run_rttov(rttov_atm_new, rttov_opt, s3com, cld)
      call ret_K(rttov_atm_new, rttov_opt, s3com, cld)
      call ret_Sf(rttov_atm_new, rttov_opt, s3com)
      call ret_cov_mat(s3com)
      
   end subroutine ret_update_mat
   ! ----------------------------------------------------------------------------------------------------------------------------
   
   ! ----------------------------------------------------------------------------------------------------------------------------
   !> @brief Calculation of the state vector best estimate and its associated cost function
   !> @details This subroutine calculates the state vector best estimate and the corresponding cost function.
   !> @param[in] rttov_atm
   !> @param[in] rttov_opt
   !> @param[in] nidx
   !> @param[in] idx_ret
   !> @param[inout] s3com retrieval data structure
   subroutine ret_update(rttov_atm, rttov_opt, nidx, idx_ret, s3com, cld)
      
      ! Inputs
      type(type_model), intent(in) :: rttov_atm
      type(type_rttov_opt), intent(in) :: rttov_opt
      type(type_cld), intent(in) :: cld
      integer(kind=4), intent(in) :: nidx
      integer(kind=4), dimension(nidx), intent(in) :: idx_ret
      
      ! Inputs/outputs
      type(type_s3com), intent(inout) :: s3com
      
      ! Internal
      integer(kind=4) :: ipoint
      
      ! First run of optimal estimation
      do ipoint = 1, size(idx_ret)
         s3com%ret%Sx(idx_ret(ipoint),:,:) = ret_Sx(s3com,idx_ret(ipoint))
         s3com%ret%Xip1(idx_ret(ipoint),:) = ret_Xip1(s3com,idx_ret(ipoint))
      enddo
      
      ! Update Xi
      do ipoint = 1, size(idx_ret)
         !if (s3com%ret%Xip1(idx_ret(ipoint),2) .lt. 1._wp) s3com%ret%Xip1(idx_ret(ipoint),2) = 1._wp   ! min Deff is 2 um
         !if (s3com%ret%Xip1(idx_ret(ipoint),2) .gt. 25._wp) s3com%ret%Xip1(idx_ret(ipoint),2) = 25._wp ! max Deff is 52 um
         s3com%ret%Xi(idx_ret(ipoint),:) = s3com%ret%Xip1(idx_ret(ipoint),:)
      enddo
      
      ! Update error variance-covariance matrices
      call ret_update_mat(rttov_atm, rttov_opt, s3com, cld)
      
      ! Update cost function
      do ipoint = 1, size(idx_ret)
         
         call inverse(s3com%ret%Se(idx_ret(ipoint),:,:),s3com%ret%Se_inv(idx_ret(ipoint),:,:),s3com%ret%nmeas)
         call inverse(s3com%ret%Sa(idx_ret(ipoint),:,:),s3com%ret%Sa_inv(idx_ret(ipoint),:,:),s3com%ret%nstates)
         
         s3com%ret%Kt(idx_ret(ipoint),:,:) = transpose(s3com%ret%K(idx_ret(ipoint),:,:))
         
         s3com%ret%Jip1(idx_ret(ipoint)) = ret_J(s3com,idx_ret(ipoint))
         
      enddo
      
   end subroutine ret_update
   ! ----------------------------------------------------------------------------------------------------------------------------
   
   ! ----------------------------------------------------------------------------------------------------------------------------
   !> @brief Free retrieval data structure
   !> @details This subroutine frees the retrieval data structure by deallocating the required arrays.
   !> @param[inout] s3com retrieval data structure
   subroutine ret_free(s3com)
      
      type(type_s3com), intent(inout) :: s3com
      
      ! Cloud properties
      deallocate(s3com%ret%ztop_liq_idx, s3com%ret%zbase_liq_idx, s3com%ret%ztop_liq, s3com%ret%zbase_liq)
      deallocate(s3com%ret%lwp, s3com%ret%cod, s3com%ret%lwc, s3com%ret%clc, s3com%ret%reff, s3com%ret%cdnc)
      
      ! Optimal estimation method
      deallocate(s3com%ret%Y, s3com%ret%F)
      deallocate(s3com%ret%n_iter, s3com%ret%Ji, s3com%ret%Jip1, s3com%ret%Ji_meas, s3com%ret%Xa, s3com%ret%Xi, s3com%ret%Xip1)
      
      ! Error variance-covariance matrices
      deallocate(s3com%ret%Se, s3com%ret%Se_inv, s3com%ret%Sy, s3com%ret%Sf, s3com%ret%Sa, s3com%ret%Sa_inv, s3com%ret%Sx, &
                 s3com%ret%Sx_inv, s3com%ret%K, s3com%ret%Kt)
      
      ! Non-retrieved parameters of the forward model
      deallocate(s3com%ret%error_t, s3com%ret%error_q, s3com%ret%error_clc, s3com%ret%error_surf_brdf, &
                 s3com%ret%error_surf_emiss, s3com%ret%error_cld_top, s3com%ret%error_cld_base)
      deallocate(s3com%ret%sigma_t, s3com%ret%sigma_q, s3com%ret%sigma_clc, s3com%ret%sigma_surf_brdf, &
                 s3com%ret%sigma_surf_emiss, s3com%ret%sigma_cld_top, s3com%ret%sigma_cld_base)
      deallocate(s3com%ret%Kb, s3com%ret%Kbt, s3com%ret%Sb)
      
      ! Flags
      deallocate(s3com%ret%flag_rttov)
      deallocate(s3com%ret%flag_conv_test)
      
   end subroutine ret_free
   ! ----------------------------------------------------------------------------------------------------------------------------
   
end module mod_ret_utils

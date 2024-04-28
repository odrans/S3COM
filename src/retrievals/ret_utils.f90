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
   
   use s3com_types,         only: wp, wpi, type_model, type_rttov_opt, type_s3com, type_cld
   use s3com_config,        only: nstates, lwp_apriori_error, cdnc_apriori_error, iwp_threshold, lwp_threshold
   use mod_ret_cloud_model, only: init_cloud_prof
   use mod_rttov_utils,     only: find_ret_idx_rttov
   use mod_rttov,           only: run_rttov
   use mod_utils_math,      only: inverse
   
   implicit none
   
   public
   
contains
   
   ! ============================================================================================================================
   !> @brief Initialize retrieval data structure
   !! @details This subroutine initializes the retrieval data structure by allocating the required arrays. All arrays are 
   !! initialized to zero.
   !! @param[inout] s3com   s3com retrieval structure
   subroutine ret_init(s3com)
      
      ! Input/output
      type(type_s3com), intent(inout) :: s3com
      
      ! Internal
      integer(wpi) :: npoints, nlevels, nlayers, nmeas, nstates, nbparams
      
      s3com%ret%npoints = s3com%npoints
      s3com%ret%nlevels = s3com%nlevels
      s3com%ret%nlayers = s3com%nlayers
      s3com%ret%nmeas   = s3com%nmeas
      s3com%ret%nstates = s3com%nstates
      
      npoints  = s3com%ret%npoints
      nlevels  = s3com%ret%nlevels
      nlayers  = s3com%ret%nlayers
      nmeas    = s3com%ret%nmeas
      nstates  = s3com%ret%nstates
      nbparams = 3*(s3com%ret%nlevels + 1) !< Number of non-retrieved parameters of the forward model
      
      ! Cloud properties
      allocate(s3com%ret%ztop_liq_idx(npoints), source = 0)
      allocate(s3com%ret%zbase_liq_idx(npoints), source = 0)
      allocate(s3com%ret%n_slope_liq(npoints), source = 0)
      allocate(s3com%ret%ztop_liq(npoints), source = 0._wp)
      allocate(s3com%ret%zbase_liq(npoints), source = 0._wp)
      allocate(s3com%ret%H(npoints), source = 0._wp)
      allocate(s3com%ret%fad(npoints), source = 0._wp)
      allocate(s3com%ret%lwp(npoints), source = 0._wp)
      allocate(s3com%ret%lwp_ad(npoints), source = 0._wp)
      allocate(s3com%ret%lwp_hom(npoints), source = 0._wp)
      allocate(s3com%ret%iwp(npoints), source = 0._wp)
      allocate(s3com%ret%cod(npoints), source = 0._wp)
      allocate(s3com%ret%clc(npoints, nlayers), source = 0._wp)
      allocate(s3com%ret%lwc_ad(npoints, nlayers), source = 0._wp)
      allocate(s3com%ret%lwc_corr(npoints, nlayers), source = 0._wp)
      allocate(s3com%ret%lwc_hom(npoints, nlayers), source = 0._wp)
      allocate(s3com%ret%re_ad(npoints, nlayers), source = 0._wp)
      allocate(s3com%ret%re_hom(npoints, nlayers), source = 0._wp)
      allocate(s3com%ret%cdnc_ad(npoints, nlayers), source = 0._wp)
      allocate(s3com%ret%cdnc_hom(npoints, nlayers), source = 0._wp)
      
      ! Optimal estimation method
      allocate(s3com%ret%Y(npoints, nmeas), source = 0._wp)
      allocate(s3com%ret%F(npoints, nmeas), source = 0._wp)
      allocate(s3com%ret%n_iter(npoints), source = 0)
      allocate(s3com%ret%gamma(npoints), source = 0._wp)
      allocate(s3com%ret%J(npoints), source = 0._wp)
      allocate(s3com%ret%Ji(npoints), source = 0._wp)
      allocate(s3com%ret%Jip1(npoints), source = 0._wp)
      allocate(s3com%ret%Ji_meas(npoints), source = 0._wp)
      allocate(s3com%ret%Ji_state(npoints), source = 0._wp)
      allocate(s3com%ret%X(npoints, nstates), source = 0._wp)
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
      
      ! Non-retrieved parameters of the forward model
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
      allocate(s3com%ret%flag_rttov(npoints)); s3com%ret%flag_rttov = .false.         !< Flag to use RTTOV
      allocate(s3com%ret%flag_conv_test(npoints)); s3com%ret%flag_conv_test = .false. !< Flag to stop the OE
      
   end subroutine ret_init
   ! ============================================================================================================================
   
   ! ============================================================================================================================
   !> @brief Count the number of single-layer homogeneous liquid clouds
   !! @param[in] s3com      s3com retrieval structure
   !! @param[in] flag_liq   flag
   function idx_liq(s3com, flag_liq)
      
      ! Input
      type(type_s3com), intent(in) :: s3com
      logical, intent(in) :: flag_liq
      
      ! Internal
      integer(wpi), dimension(:), allocatable :: idx_liq
      integer(wpi), dimension(s3com%ret%npoints) :: idx_all
      integer(wpi) :: idx, i
      logical :: test
      
      idx_all = 0; idx = 1
      
      do i = 1, s3com%ret%npoints
         
         test = .false.
         
         if (flag_liq) then
            test = s3com%ret%n_slope_liq(i) .eq. 1 .and. s3com%ret%lwp(i) .gt. lwp_threshold .and. &
                   s3com%ret%iwp(i) .lt. iwp_threshold
         endif
         
         if (test) then
            idx_all(idx) = i
            idx = idx + 1
         endif
         
      enddo
      
      idx = idx - 1
      
      allocate(idx_liq(idx), source=0); idx_liq = idx_all(1:idx)
      
      return
      
   end function idx_liq
   ! ============================================================================================================================
   
   ! ============================================================================================================================
   !> @brief Count the number of single-layer homogeneous liquid clouds
   !! @param[in] s3com      s3com retrieval structure
   !! @param[in] flag_liq   flag
   subroutine find_ret_idx_cld_liq(s3com, flag_liq)
      
      ! Input
      logical, intent(in) :: flag_liq
      type(type_s3com), intent(inout) :: s3com
      
      ! Input/output
      
      
      ! Internal
      integer(wpi) :: ipoint
      !integer(wpi), intent(in) :: ipoint
      logical :: test
      !logical, dimension(s3com%ret%npoints) :: flag_rttov
      logical, dimension(:), allocatable :: flag_rttov
      
      allocate(flag_rttov(s3com%ret%npoints), source = .false.)
      
      write(6,*) "flag_liq:", flag_liq
      
      do ipoint = 1, s3com%ret%npoints
         !write(6,*) "ipoint:", ipoint
         test = .false.
         
         if (flag_liq) then
            test = s3com%ret%n_slope_liq(ipoint) .eq. 1 .and. s3com%ret%lwp(ipoint) .gt. lwp_threshold .and. &
                   s3com%ret%iwp(ipoint) .lt. iwp_threshold
         endif
         
         if (test) then
            flag_rttov(ipoint) = .true.
         else
            flag_rttov(ipoint) = .false.
         endif
         
         write(6,*) s3com%ret%n_slope_liq(ipoint), s3com%ret%lwp(ipoint), s3com%ret%iwp(ipoint), flag_rttov(ipoint)
         
         !s3com%ret%flag_rttov(i) = flag_rttov(i)
         !write(6,*) "s3com%ret%flag_rttov:", s3com%ret%flag_rttov(i)
      enddo
      
      s3com%ret%flag_rttov = flag_rttov
      
      do ipoint = 1, s3com%ret%npoints
         write(6,*) "s3com%ret%flag_rttov:", s3com%ret%flag_rttov(ipoint)
      enddo
      !return
      
   end subroutine find_ret_idx_cld_liq
   ! ============================================================================================================================
   
   ! ============================================================================================================================
   !> @brief Count the number of single-layer homogeneous liquid clouds
   !! @param[in] s3com      s3com retrieval structure
   !! @param[in] flag_liq   flag
   function find_ret_idx_liq(s3com, flag_liq, ipoint) result(flag_rttov)
      
      ! Input
      logical, intent(in) :: flag_liq
      type(type_s3com), intent(in) :: s3com
      
      ! Input/output
      
      
      ! Internal
      integer(wpi), intent(in) :: ipoint
      logical :: test
      !logical, dimension(s3com%ret%npoints) :: flag_rttov
      logical, dimension(:), allocatable :: flag_rttov
      
      allocate(flag_rttov(s3com%ret%npoints), source = .false.)
      write(6,*) "flag_liq:", flag_liq
      !do ipoint = 1, s3com%ret%npoints
         write(6,*) "ipoint:", ipoint
         test = .false.
         
         if (flag_liq) then
            test = s3com%ret%n_slope_liq(ipoint) .eq. 1 .and. s3com%ret%lwp(ipoint) .gt. lwp_threshold .and. &
                   s3com%ret%iwp(ipoint) .lt. iwp_threshold
         endif
         write(6,*) s3com%ret%n_slope_liq(ipoint), s3com%ret%lwp(ipoint), lwp_threshold, s3com%ret%iwp(ipoint), iwp_threshold
         if (test) then
            flag_rttov(ipoint) = .true.
         else
            flag_rttov(ipoint) = .false.
         endif
         
         write(6,*) "flag_rttov:", flag_rttov(ipoint)
         
         !s3com%ret%flag_rttov(i) = flag_rttov(i)
         !write(6,*) "s3com%ret%flag_rttov:", s3com%ret%flag_rttov(i)
      !enddo
      
      !s3com%ret%flag_rttov = flag_rttov
      !write(6,*) "s3com%ret%flag_rttov:", s3com%ret%flag_rttov(i)
      !return
      
   end function find_ret_idx_liq
   ! ============================================================================================================================
   
   ! ============================================================================================================================
   !> @brief Cost function: sum of the contributions from a priori and from measurements
   !! @param[in] s3com    s3com retrieval structure
   !! @param[in] ipoint   index of point i
   real(wp) function ret_J(s3com, ipoint)
      
      ! Input
      type(type_s3com), intent(in) :: s3com
      integer(wpi), intent(in) :: ipoint
      
      ret_J = J_state(s3com,ipoint) + J_meas(s3com,ipoint)
      
      return
      
   end function ret_J
   ! ============================================================================================================================
   
   ! ============================================================================================================================
   !> @brief Cost function: contribution from a priori
   !! @param[in] s3com    s3com retrieval structure
   !! @param[in] ipoint   index of point i
   real(wp) function J_state(s3com, ipoint)
      
      ! Input
      type(type_s3com), intent(in) :: s3com
      integer(wpi), intent(in) :: ipoint
      
      J_state = dot_product(s3com%ret%Xi(ipoint,:)-s3com%ret%Xa(ipoint,:), &
                matmul(s3com%ret%Sa_inv(ipoint,:,:),s3com%ret%Xi(ipoint,:)-s3com%ret%Xa(ipoint,:)))
      
      return
      
   end function J_state
   ! ============================================================================================================================
   
   ! ============================================================================================================================
   !> @brief Cost function: contribution from measurements
   !! @param[in] s3com    s3com retrieval structure
   !! @param[in] ipoint   index of point i
   real(wp) function J_meas(s3com, ipoint)
      
      ! Input
      type(type_s3com), intent(in) :: s3com
      integer(wpi), intent(in) :: ipoint
      
      J_meas = dot_product(s3com%ret%Y(ipoint,:)-s3com%ret%F(ipoint,:), &
               matmul(s3com%ret%Se_inv(ipoint,:,:),s3com%ret%Y(ipoint,:)-s3com%ret%F(ipoint,:)))
      
      return
      
   end function J_meas
   ! ============================================================================================================================
   
   ! ============================================================================================================================
   !> @brief State vector best estimate error variance-covariance matrix
   !! @param[in] s3com    s3com retrieval structure
   !! @param[in] ipoint   index of point i
   function ret_Sx(s3com, ipoint)
      
      ! Input
      type(type_s3com), intent(in) :: s3com
      integer(wpi), intent(in) :: ipoint
      
      ! Internal
      real(wp), dimension(nstates,nstates) :: ret_Sx, Sx_inv
      
      ! Levenberg-Marquardt method
      Sx_inv(:,:) = (1.0_wp + s3com%ret%gamma(ipoint))*s3com%ret%Sa_inv(ipoint,:,:) + &
                    matmul(s3com%ret%Kt(ipoint,:,:),matmul(s3com%ret%Se_inv(ipoint,:,:),s3com%ret%K(ipoint,:,:)))
      
      call inverse(Sx_inv, ret_Sx, nstates)
      
      return
      
   end function ret_Sx
   ! ============================================================================================================================
   
   ! ============================================================================================================================
   !> @brief State vector best estimate following the Levenberg-Marquardt method
   !! @param[in] s3com    s3com retrieval structure
   !! @param[in] ipoint   index of point i
   function ret_Xip1(s3com, ipoint)
      
      ! Input
      type(type_s3com), intent(in) :: s3com
      integer(wpi), intent(in) :: ipoint
      
      ! Internal
      real(wp), dimension(nstates) :: ret_Xip1
      
      ret_Xip1(:) = s3com%ret%Xi(ipoint,:) + &
                    matmul(s3com%ret%Sx(ipoint,:,:),(matmul(s3com%ret%Kt(ipoint,:,:), &
                    matmul(s3com%ret%Se_inv(ipoint,:,:),(s3com%ret%Y(ipoint,:)-s3com%ret%F(ipoint,:)))) - &
                    matmul(s3com%ret%Sa_inv(ipoint,:,:),(s3com%ret%Xi(ipoint,:)-s3com%ret%Xa(ipoint,:)))))
      
      return
      
   end function ret_Xip1
   ! ============================================================================================================================
   
   ! ============================================================================================================================
   !> @brief Jacobian matrix K
   !! @details This subroutine calculates the sensitivity of each observation to each state to be retrieved.
   !! @param[inout] rttov_atm   rttov atmospheric structure
   !! @param[in] rttov_opt      rttov options structure
   !! @param[inout] s3com       s3com structure
   !! @param[in] cld            cloud structure
   subroutine ret_K(rttov_atm, rttov_opt, s3com, cld)
      
      ! Input
      type(type_rttov_opt), intent(in) :: rttov_opt
      type(type_cld), intent(in) :: cld
      
      ! Input/output
      type(type_model), intent(in) :: rttov_atm
      type(type_s3com), intent(inout) :: s3com
      
      ! Internal
      type(type_s3com) :: s3com_pert
      type(type_model) :: rttov_atm_pert
      real(wp), dimension(rttov_atm%npoints, rttov_opt%nchannels) :: Fp1dx, Fp2dx, Fm1dx, Fm2dx
      integer(wpi), dimension(:), allocatable :: idx_ret
      integer(wpi) :: ipoint, idx, imeas
      real(wp), parameter ::           &
         lwp_delta_pert  = 1.0E-01_wp, & !< Input perturbation of 10 %
         cdnc_delta_pert = 1.0E-01_wp    !< Input perturbation of 10 %
      
      s3com_pert = s3com; rttov_atm_pert = rttov_atm
      
      idx_ret = find_ret_idx_rttov(s3com)
      
      ! Jacobian of the cloud liquid water path (LWP)
      ! -------------------------------------------------------------------------------------------------------------------------
      do ipoint = 1, size(idx_ret)
         
         idx = idx_ret(ipoint)
         
         Fp1dx = 0.0_wp; Fm1dx = 0.0_wp; Fp2dx = 0.0_wp; Fm2dx = 0.0_wp
         
         ! Positive perturbation
         ! ----------------------------------------------------------------------------------------------------------------------
         s3com_pert%ret%Xi(idx,1) = s3com%ret%Xi(idx,1) * (1 + lwp_delta_pert)
         call init_cloud_prof(rttov_atm_pert, s3com_pert)
         call run_rttov(rttov_atm_pert, rttov_opt, s3com_pert, cld)
         Fp1dx = s3com_pert%ret%F
         !write(*,*) "LWP_k=", "Fp1dx=", Fp1dx
         
         s3com_pert%ret%Xi(idx,1) = s3com%ret%Xi(idx,1) * (1 + 2.0_wp*lwp_delta_pert)
         call init_cloud_prof(rttov_atm_pert, s3com_pert)
         call run_rttov(rttov_atm_pert, rttov_opt, s3com_pert, cld)
         Fp2dx = s3com_pert%ret%F
         !write(*,*) "LWP_k=", "Fp2dx=", Fp2dx
         ! ----------------------------------------------------------------------------------------------------------------------
         
         ! Negative perturbation
         ! ----------------------------------------------------------------------------------------------------------------------
         s3com_pert%ret%Xi(idx,1) = s3com%ret%Xi(idx,1) * (1 - lwp_delta_pert)
         call init_cloud_prof(rttov_atm_pert, s3com_pert)
         call run_rttov(rttov_atm_pert, rttov_opt, s3com_pert, cld)
         Fm1dx = s3com_pert%ret%F
         !write(*,*) "LWP_k=", "Fm1dx=", Fm1dx
         
         s3com_pert%ret%Xi(idx,1) = s3com%ret%Xi(idx,1) * (1 - 2.0_wp*lwp_delta_pert)
         call init_cloud_prof(rttov_atm_pert, s3com_pert)
         call run_rttov(rttov_atm_pert, rttov_opt, s3com_pert, cld)
         Fm2dx = s3com_pert%ret%F
         !write(*,*) "LWP_k=", "Fm2dx=", Fm2dx
         ! ----------------------------------------------------------------------------------------------------------------------
         
         ! Jacobian
         ! ----------------------------------------------------------------------------------------------------------------------
         do imeas = 1, rttov_opt%nchannels
            ! Order 2
            !write(*,*) "s3com%ret%Xi(idx,1)=", s3com%ret%Xi(idx,1)
            s3com%ret%K(idx,imeas,1) = (Fp1dx(idx,imeas)-Fm1dx(idx,imeas)) / (2.0_wp*s3com%ret%Xi(idx,1)*lwp_delta_pert)
            ! Order 3
            !s3com%ret%K(idx,imeas,1) = (-Fp2dx(idx,imeas) + 9.0_wp*Fp1dx(idx,imeas) - 9.0_wp*Fm1dx(idx,imeas) + Fm2dx(idx,imeas)) / &
            !(24.0_wp*s3com%ret%Xi(idx,1)*lwp_delta_pert)
            ! Order 4
            !s3com%ret%K(idx,imeas,1) = (-Fp2dx(idx,imeas) + 8.0_wp*Fp1dx(idx,imeas) - 8.0_wp*Fm1dx(idx,imeas) + Fm2dx(idx,imeas)) / &
            !(12.0_wp*s3com%ret%Xi(idx,1)*lwp_delta_pert)
            !write(*,*) "K1 =", s3com%ret%K(idx,imeas,1)
         enddo
         ! ----------------------------------------------------------------------------------------------------------------------
         
      enddo
      
      ! Reset Xi to its initial value
      s3com_pert%ret%Xi(idx,1) = s3com%ret%Xi(idx,1)
      ! -------------------------------------------------------------------------------------------------------------------------
      
      ! Jacobian of the cloud drouplet number concentration (CDNC)
      ! -------------------------------------------------------------------------------------------------------------------------
      do ipoint = 1, size(idx_ret)
         
         idx = idx_ret(ipoint)
         
         Fp1dx = 0.0_wp; Fm1dx = 0.0_wp; Fp2dx = 0.0_wp; Fm2dx = 0.0_wp
         
         ! Positive perturbation
         ! ----------------------------------------------------------------------------------------------------------------------
         s3com_pert%ret%Xi(idx,2) = s3com%ret%Xi(idx,2) * (1 + cdnc_delta_pert)
         !s3com_pert%ret%Xi(idx,2) = s3com%ret%Xi(idx,2) + cdnc_delta_pert
         call init_cloud_prof(rttov_atm_pert, s3com_pert)
         call run_rttov(rttov_atm_pert, rttov_opt, s3com_pert, cld)
         Fp1dx = s3com_pert%ret%F
         !write(*,*) "CNDC_k=", "Fp1dx=", Fp1dx
         
         s3com_pert%ret%Xi(idx,2) = s3com%ret%Xi(idx,2) * (1 + 2.0_wp*cdnc_delta_pert)
         !s3com_pert%ret%Xi(idx,2) = s3com%ret%Xi(idx,2) + 2.0_wp*cdnc_delta_pert
         call init_cloud_prof(rttov_atm_pert, s3com_pert)
         call run_rttov(rttov_atm_pert, rttov_opt, s3com_pert, cld)
         Fp2dx = s3com_pert%ret%F
         !write(*,*) "CNDC_k=", "Fp2dx=", Fp2dx
         ! ----------------------------------------------------------------------------------------------------------------------
         
         ! Negative perturbation
         ! ----------------------------------------------------------------------------------------------------------------------
         s3com_pert%ret%Xi(idx,2) = s3com%ret%Xi(idx,2) * (1 - cdnc_delta_pert)
         !s3com_pert%ret%Xi(idx,2) = s3com%ret%Xi(idx,2) - cdnc_delta_pert
         call init_cloud_prof(rttov_atm_pert, s3com_pert)
         call run_rttov(rttov_atm_pert, rttov_opt, s3com_pert, cld)
         Fm1dx = s3com_pert%ret%F
         !write(*,*) "CNDC_k=", "Fm1dx=", Fm1dx
         
         s3com_pert%ret%Xi(idx,2) = s3com%ret%Xi(idx,2) * (1 - 2.0_wp*cdnc_delta_pert)
         !s3com_pert%ret%Xi(idx,2) = s3com%ret%Xi(idx,2) - 2.0_wp*cdnc_delta_pert
         call init_cloud_prof(rttov_atm_pert, s3com_pert)
         call run_rttov(rttov_atm_pert, rttov_opt, s3com_pert, cld)
         Fm2dx = s3com_pert%ret%F
         !write(*,*) "CNDC_k=", "Fm2dx=", Fm2dx
         ! ----------------------------------------------------------------------------------------------------------------------
         
         ! Jacobian
         ! ----------------------------------------------------------------------------------------------------------------------
         do imeas = 1, rttov_opt%nchannels
            ! Order 2
            s3com%ret%K(idx,imeas,2) = (Fp1dx(idx,imeas)-Fm1dx(idx,imeas))/(2.0_wp*s3com%ret%Xi(idx,2)*cdnc_delta_pert)
            ! Order 4
            !s3com%ret%K(idx,imeas,2) = (-Fp2dx(idx,imeas) + 8.0_wp*Fp1dx(idx,imeas) - 8.0_wp*Fm1dx(idx,imeas) + Fm2dx(idx,imeas)) / &
            !(12.0_wp*s3com%ret%Xi(idx,2)*cdnc_delta_pert)
            !write(*,*) "K2 =", s3com%ret%K(idx,imeas,2)
         enddo
         ! ----------------------------------------------------------------------------------------------------------------------
         
      enddo
      
      ! Reset Xi to its initial value
      s3com_pert%ret%Xi(idx,2) = s3com%ret%Xi(idx,2)
      ! -------------------------------------------------------------------------------------------------------------------------
      
   end subroutine ret_K
   ! ============================================================================================================================
   
   ! ============================================================================================================================
   !> @brief Forward model error variance-covariance matrix calculation
   !! @details This subroutine calculates the sensitivity of each observation to each non-retrieved forward model parameter 
   !! (Jacobian matrix Kb), the non-retrieved forward model parameters error variance-covariance matrix (Sb) and the forward 
   !! model error variance-covariance matrix (Sf).
   !! @param[inout] rttov_atm   rttov atmospheric structure
   !! @param[in] rttov_opt      rttov options structure
   !! @param[inout] s3com       s3com structure
   subroutine ret_Sf(rttov_atm, rttov_opt, s3com)
      
      ! Input
      type(type_rttov_opt), intent(in) :: rttov_opt
      
      ! Input/output
      type(type_model), intent(inout) :: rttov_atm
      type(type_s3com), intent(inout) :: s3com
      
      ! Internal
      integer(wpi), dimension(:), allocatable :: idx_ret
      integer(wpi) :: npoints, nlevels, nlayers, nmeas
      integer(wpi) :: ipoint, idx, imeas, ilayer, ilevel
      integer(wpi) :: idx_t, idx_q, idx_clc, idx_surf_alb, idx_surf_emiss, idx_cld_top, idx_cld_base
      
      npoints = rttov_atm%npoints
      nlevels = rttov_atm%nlevels
      nlayers = rttov_atm%nlayers
      nmeas   = rttov_opt%nchannels
      
      idx_ret = find_ret_idx_rttov(s3com)
      
      ! Jacobian matrix Kb
      ! -------------------------------------------------------------------------------------------------------------------------
      do ipoint = 1, size(idx_ret)
         
         idx = idx_ret(ipoint)
         
         do imeas = 1, rttov_opt%nchannels
            
            !if (s3com%ret%flag_rttov(idx)) then
               
               s3com%ret%Kb(idx,imeas,1:nlevels)               = s3com%jac%t(idx,imeas,1:nlevels)     !< Temperature profile
               s3com%ret%Kb(idx,imeas,nlevels+1:2*nlevels)     = s3com%jac%q(idx,imeas,1:nlevels)     !< Humidity profile
               s3com%ret%Kb(idx,imeas,2*nlevels+1:3*nlevels-1) = s3com%jac%cfrac(idx,imeas,1:nlayers) !< Cloud fraction profile
               s3com%ret%Kb(idx,imeas,3*nlevels:3*nlevels)     = s3com%jac%brdf(idx,imeas)            !< Surface brdf/albedo
               s3com%ret%Kb(idx,imeas,3*nlevels+1:3*nlevels+1) = s3com%jac%emiss(idx,imeas)           !< Surface emissivity
               s3com%ret%Kb(idx,imeas,3*nlevels+2:3*nlevels+2) = 0.0_wp                               !< Cloud top altitude
               s3com%ret%Kb(idx,imeas,3*nlevels+3:3*nlevels+3) = 0.0_wp                               !< Cloud base altitude
               
            !endif
            
         enddo
         
         s3com%ret%Kbt(idx,:,:) = transpose(s3com%ret%Kb(idx,:,:))
         
      enddo
      ! -------------------------------------------------------------------------------------------------------------------------
      
      ! Non-retrieved forward model parameters error variance-covariance matrix (Sb is diagonal)
      ! -------------------------------------------------------------------------------------------------------------------------
      do ipoint = 1, size(idx_ret)
         
         idx = idx_ret(ipoint)
         
         !if (s3com%ret%flag_rttov(idx)) then
            
            ! Attribution of errors/uncertainties
            ! -------------------------------------------------------------------------------------------------------------------
            s3com%ret%error_t(idx,1:nlevels) = 1.0_wp / rttov_atm%t(idx,1:nlevels) !< Temperature: error of 1 K per level
            s3com%ret%error_q(idx,1:nlevels) = 0.1_wp                              !< Humidity: error of 10 % per level
            s3com%ret%error_clc(idx,1:nlayers) = 0.01_wp                           !< Cloud fraction: error of 1 % per layer
            s3com%ret%error_surf_brdf(idx,1:nmeas) = 0.01_wp                       !< Surface brdf/albedo: error of 1 % per channel
            s3com%ret%error_surf_emiss(idx,1:nmeas) = 0.01_wp                      !< Surface emissivity: error of 1 % per channel
            s3com%ret%error_cld_top(idx) = 100.0_wp / s3com%ret%ztop_liq(idx)      !< Cloud top altitude: error of 100 m
            s3com%ret%error_cld_base(idx) = 100.0_wp / s3com%ret%zbase_liq(idx)    !< Cloud base altitude: error of 100 m
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
            
            ! Cloud top altitude
            s3com%ret%sigma_cld_top(idx) = s3com%ret%ztop_liq(idx)*s3com%ret%error_cld_top(idx)
            
            ! Cloud base altitude
            s3com%ret%sigma_cld_base(idx) = s3com%ret%zbase_liq(idx)*s3com%ret%error_cld_base(idx)
            ! -------------------------------------------------------------------------------------------------------------------
            
            ! Covariances
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
            
         !endif
         
      enddo
      ! -------------------------------------------------------------------------------------------------------------------------
      
      ! Forward model error variance-covariance matrix (Sf is non diagonal)
      ! -------------------------------------------------------------------------------------------------------------------------
      do ipoint = 1, size(idx_ret)
         
         idx = idx_ret(ipoint)
         
         !if (s3com%ret%flag_rttov(idx)) then
            s3com%ret%Sf(idx,:,:) = matmul(s3com%ret%Kb(idx,:,:),matmul(s3com%ret%Sb(idx,:,:),s3com%ret%Kbt(idx,:,:)))
         !endif
         
      enddo
      ! -------------------------------------------------------------------------------------------------------------------------
      
   end subroutine ret_Sf
   ! ============================================================================================================================
   
   ! ============================================================================================================================
   !> @brief Calculation of the error variance-covariance matrices
   !! @details This subroutine calculates the a priori error variance-covariance matrix Sa, the measurement error 
   !! variance-covariance matrix Sy, and the total error variance-covariance matrix in the measurement space Se.
   !! @param[inout] s3com retrieval data structure
   subroutine ret_cov_mat(s3com)
      
      ! Input/output
      type(type_s3com), intent(inout) :: s3com
      
      ! Internal
      integer(wpi) :: ipoint, ichannel
      
      ! Loop over each point
      do ipoint = 1, s3com%ret%npoints
         
         ! A priori error variance-covariance matrix (Sa is diagonal, non-diagonal terms are considered to be null)
         ! ----------------------------------------------------------------------------------------------------------------------
         s3com%ret%Sa(ipoint,1,1) = (s3com%ret%Xa(ipoint,1)*lwp_apriori_error)**2.0_wp
         s3com%ret%Sa(ipoint,2,2) = (s3com%ret%Xa(ipoint,2)*cdnc_apriori_error)**2.0_wp
         ! ----------------------------------------------------------------------------------------------------------------------
         
         ! Measurement error variance-covariance matrix (Sy is diagonal, non-diagonal terms are considered to be null)
         ! ----------------------------------------------------------------------------------------------------------------------
         ! Loop over each channel
         do ichannel = 1, s3com%ret%nmeas
            
            ! Uncertainties are taken from Xiong et al., 2018: "Updates of Moderate Resolution Imaging Spectroradiometer on-orbit
            ! calibration uncertainty assessments"
            
            ! Uncertainty of 5 % in MODIS radiances for the reflective solar bands (1-19, 26)
            if (s3com%nml%channel_list(ichannel) .ge. 1 .and. s3com%nml%channel_list(ichannel) .le. 19 .or. &
               s3com%nml%channel_list(ichannel) .eq. 26) then
               s3com%ret%Sy(ipoint,ichannel,ichannel) = (s3com%ret%Y(ipoint,ichannel)*0.05_wp)**2.0_wp
            endif
            
            ! Uncertainty of 1 % in MODIS radiances for the thermal emissive bands (20-25, 27-36)
            if (s3com%nml%channel_list(ichannel) .ge. 20 .and. s3com%nml%channel_list(ichannel) .le. 25 .or. &
               s3com%nml%channel_list(ichannel) .ge. 27 .and. s3com%nml%channel_list(ichannel) .le. 36) then
               s3com%ret%Sy(ipoint,ichannel,ichannel) = (s3com%ret%Y(ipoint,ichannel)*0.01_wp)**2.0_wp
            endif
            
            ! Uncertainty of 0.75 % in MODIS radiance for the thermal emissive band 20
            if (s3com%nml%channel_list(ichannel) .eq. 20) then
               s3com%ret%Sy(ipoint,ichannel,ichannel) = (s3com%ret%Y(ipoint,ichannel)*0.0075_wp)**2.0_wp
            endif
            
            ! Uncertainty of 0.50 % in MODIS radiance for the thermal emissive band 31
            if (s3com%nml%channel_list(ichannel) .eq. 31) then
               s3com%ret%Sy(ipoint,ichannel,ichannel) = (s3com%ret%Y(ipoint,ichannel)*0.0050_wp)**2.0_wp
            endif
            
         enddo
         ! ----------------------------------------------------------------------------------------------------------------------
         
      enddo
      
      ! Total error variance-covariance matrix in the measurement space (Se is non diagonal)
      ! -------------------------------------------------------------------------------------------------------------------------
      s3com%ret%Se = s3com%ret%Sy + s3com%ret%Sf
      ! -------------------------------------------------------------------------------------------------------------------------
      
   end subroutine ret_cov_mat
   ! ============================================================================================================================
   
   ! ============================================================================================================================
   !> @brief Update the error variance-covariance matrices
   !! @param[in] rttov_atm   rttov atmospheric structure
   !! @param[in] rttov_opt   rttov optical structure
   !! @param[inout] s3com    s3com retrieval structure
   !! @param[in] cld         cloud structure
   subroutine ret_update_mat(rttov_atm, rttov_opt, s3com, cld)
      
      ! Input
      type(type_model), intent(in) :: rttov_atm
      type(type_rttov_opt), intent(in) :: rttov_opt
      type(type_cld), intent(in) :: cld
      
      ! Input/output
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
   ! ============================================================================================================================
   
   ! ============================================================================================================================
   !> @brief Calculation of the state vector best estimate and the cost function
   !! @details This subroutine calculates the state vector best estimate and the cost function.
   !! @param[in] rttov_atm   rttov atmospheric structure
   !! @param[in] rttov_opt   rttov options structure
   !! @param[in] nidx        number of 
   !! @param[in] idx_ret     index
   !! @param[in] cld         cloud structure
   !! @param[inout] s3com    s3com retrieval structure
   subroutine ret_update(rttov_atm, rttov_opt, nidx, idx_ret, s3com, cld)
      
      ! Input
      type(type_model), intent(in) :: rttov_atm
      type(type_rttov_opt), intent(in) :: rttov_opt
      type(type_cld), intent(in) :: cld
      integer(wpi), intent(in) :: nidx
      integer(wpi), dimension(nidx), intent(in) :: idx_ret
      
      ! Input/output
      type(type_s3com), intent(inout) :: s3com
      
      ! Internal
      integer(wpi) :: ipoint
      
      ! First run of optimal estimation
      do ipoint = 1, size(idx_ret)
         s3com%ret%Sx(idx_ret(ipoint),:,:) = ret_Sx(s3com,idx_ret(ipoint))
         s3com%ret%Xip1(idx_ret(ipoint),:) = ret_Xip1(s3com,idx_ret(ipoint))
      enddo
      
      ! Ensure the physical consistency of Xi, else change stepsize and update again
      do ipoint = 1, size(idx_ret)
         do while (s3com%ret%Xip1(idx_ret(ipoint),1) .le. 1.0E-04_wp .or. & !< Lower limit of LWP for a warm stratocumulus (1E-4 kg m-2)
            s3com%ret%Xip1(idx_ret(ipoint),2) .lt. 0.0_wp)                  !< Lower limit of CDNC for a warm stratocumulus (0 cm-3)
            
            s3com%ret%gamma(idx_ret(ipoint)) = s3com%ret%gamma(idx_ret(ipoint)) * 5.0_wp
            
            s3com%ret%Sx(idx_ret(ipoint),:,:) = ret_Sx(s3com,idx_ret(ipoint))
            s3com%ret%Xip1(idx_ret(ipoint),:) = ret_Xip1(s3com,idx_ret(ipoint))
         enddo
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
         s3com%ret%Ji_meas(idx_ret(ipoint)) = J_meas(s3com,idx_ret(ipoint))
         s3com%ret%Ji_state(idx_ret(ipoint)) = J_state(s3com,idx_ret(ipoint))
         
      enddo
      
   end subroutine ret_update
   ! ============================================================================================================================
   
   ! ============================================================================================================================
   !> @brief Free retrieval data structure
   !! @details This subroutine frees the retrieval data structure by deallocating the required arrays.
   !! @param[inout] s3com   s3com retrieval structure
   subroutine ret_free(s3com)
      
      ! Input/output
      type(type_s3com), intent(inout) :: s3com
      
      ! Cloud properties
      deallocate(s3com%ret%ztop_liq_idx, s3com%ret%zbase_liq_idx, s3com%ret%n_slope_liq)
      deallocate(s3com%ret%ztop_liq, s3com%ret%zbase_liq, s3com%ret%H, s3com%ret%fad, s3com%ret%lwp, s3com%ret%lwp_ad, &
                 s3com%ret%lwp_hom, s3com%ret%iwp, s3com%ret%cod, s3com%ret%clc, s3com%ret%lwc_ad, s3com%ret%lwc_corr, &
                 s3com%ret%lwc_hom, s3com%ret%re_ad, s3com%ret%re_hom, s3com%ret%cdnc_ad, s3com%ret%cdnc_hom)
      
      ! Optimal estimation method
      deallocate(s3com%ret%Y, s3com%ret%F, s3com%ret%n_iter, s3com%ret%gamma, s3com%ret%J, s3com%ret%Ji, s3com%ret%Jip1, &
                 s3com%ret%Ji_meas, s3com%ret%Ji_state, s3com%ret%X, s3com%ret%Xa, s3com%ret%Xi, s3com%ret%Xip1)
      
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
      deallocate(s3com%ret%flag_rttov, s3com%ret%flag_conv_test)
      
   end subroutine ret_free
   ! ============================================================================================================================
   
end module mod_ret_utils

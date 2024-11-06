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

module mod_s3com_setup
   
   use s3com_types,  only: wp, wpi, type_nml, type_icon, type_model, type_s3com
   use s3com_config, only: nstates
   
   implicit none
   
   public :: s3com_init, s3com_subset, s3com_update
   
contains
   
   ! ============================================================================================================================
   !> @brief Initialize the s3com structure
   !! @param[in] nml      namelist structure
   !! @param[in] model    model structure
   !! @param[out] s3com   s3com structure
   subroutine s3com_init(nml, model, s3com)
      
      ! Input
      type(type_nml), intent(in) :: nml
      type(type_model), intent(in) :: model
      
      ! Output
      type(type_s3com), intent(out) :: s3com
      
      ! Internal
      integer(wpi) :: npoints, nlayers, nlevels, nmeas
      
      npoints = model%npoints
      nlevels = model%nlevels
      nlayers = model%nlayers
      nmeas   = nml%nchannels
      
      s3com%npoints  = npoints
      s3com%nlevels  = nlevels
      s3com%nlayers  = nlayers
      s3com%nstates  = nstates
      s3com%nmeas    = nmeas
      
      s3com%date = model%date
      s3com%time = model%time
      
      allocate(s3com%rad%wavelength(nmeas), source = 0._wp)
      
      ! Radiation data
      ! -------------------------------------------------------------------------------------------------------------------------
      allocate(s3com%rad%f_rad_total(npoints, nmeas), source = 0._wp)
      allocate(s3com%rad%f_rad_clear(npoints, nmeas), source = 0._wp)
      allocate(s3com%rad%f_ref_total(npoints, nmeas), source = 0._wp)
      allocate(s3com%rad%f_ref_clear(npoints, nmeas), source = 0._wp)
      allocate(s3com%rad%f_bt_total(npoints, nmeas), source = 0._wp)
      allocate(s3com%rad%f_bt_clear(npoints, nmeas), source = 0._wp)
      allocate(s3com%rad%emiss(npoints, nmeas), source = 0._wp)
      allocate(s3com%rad%brdf(npoints, nmeas), source = 0._wp)
      ! -------------------------------------------------------------------------------------------------------------------------
      
      ! Jacobian data (K model)
      ! -------------------------------------------------------------------------------------------------------------------------
      allocate(s3com%jac%emiss(npoints, nmeas), source = 0._wp)
      allocate(s3com%jac%brdf(npoints, nmeas), source = 0._wp)
      
      allocate(s3com%jac%p(npoints, nmeas, nlevels), source = 0._wp)
      allocate(s3com%jac%t(npoints, nmeas, nlevels), source = 0._wp)
      allocate(s3com%jac%q(npoints, nmeas, nlevels), source = 0._wp)
      
      allocate(s3com%jac%cfrac(npoints, nmeas, nlayers), source = 0._wp)
      allocate(s3com%jac%clwde(npoints, nmeas, nlayers), source = 0._wp)
      ! -------------------------------------------------------------------------------------------------------------------------
      
      ! Jacobian data (TL model)
      ! -------------------------------------------------------------------------------------------------------------------------
      allocate(s3com%k_tl%t(npoints, nmeas, nlevels), source = 0._wp)
      allocate(s3com%k_tl%q(npoints, nmeas, nlevels), source = 0._wp)
      ! -------------------------------------------------------------------------------------------------------------------------
      
      s3com%nml = nml
      
   end subroutine s3com_init
   ! ============================================================================================================================
   
   ! ============================================================================================================================
   !> @brief Subset the relevant part of the s3com structure
   !! @param[in] idx_start      starting index
   !! @param[in] idx_end        ending index
   !! @param[in] s3com          s3com structure
   !! @param[out] s3com_chunk   s3com structure
   subroutine s3com_subset(idx_start, idx_end, s3com, s3com_chunk)
      
      ! Input
      type(type_s3com), intent(in) :: s3com
      integer(wpi), intent(in) :: idx_start, idx_end
      
      ! Output
      type(type_s3com), intent(out) :: s3com_chunk

      s3com_chunk%idx_start = idx_start
      s3com_chunk%idx_end   = idx_end
      s3com_chunk%npoints   = idx_end - idx_start + 1

      allocate(s3com_chunk%flag_rttov(s3com_chunk%npoints)); s3com_chunk%flag_rttov = .true.

      s3com_chunk%nlevels = s3com%nlevels
      s3com_chunk%nlayers = s3com%nlayers
      s3com_chunk%nmeas   = s3com%nmeas
      s3com_chunk%nstates = s3com%nstates

      s3com_chunk%nml%flag_retrievals = s3com%nml%flag_retrievals

      ! Radiation data
      ! -------------------------------------------------------------------------------------------------------------------------
      allocate(s3com_chunk%rad%f_rad_total(s3com_chunk%npoints, s3com_chunk%nmeas), source = 0._wp)
      allocate(s3com_chunk%rad%f_rad_clear(s3com_chunk%npoints, s3com_chunk%nmeas), source = 0._wp)
      allocate(s3com_chunk%rad%f_ref_total(s3com_chunk%npoints, s3com_chunk%nmeas), source = 0._wp)
      allocate(s3com_chunk%rad%f_ref_clear(s3com_chunk%npoints, s3com_chunk%nmeas), source = 0._wp)
      allocate(s3com_chunk%rad%f_bt_total(s3com_chunk%npoints, s3com_chunk%nmeas), source = 0._wp)
      allocate(s3com_chunk%rad%f_bt_clear(s3com_chunk%npoints, s3com_chunk%nmeas), source = 0._wp)
      allocate(s3com_chunk%rad%emiss(s3com_chunk%npoints, s3com_chunk%nmeas), source = 0._wp)
      allocate(s3com_chunk%rad%brdf(s3com_chunk%npoints, s3com_chunk%nmeas), source = 0._wp)
      ! -------------------------------------------------------------------------------------------------------------------------
      
      ! Jacobian data (K model)
      ! -------------------------------------------------------------------------------------------------------------------------
      s3com_chunk%jac%do_jacobian_calc = s3com%nml%do_jacobian_calc
      
      allocate(s3com_chunk%jac%emiss(s3com_chunk%npoints, s3com_chunk%nmeas), source = 0._wp)
      allocate(s3com_chunk%jac%brdf(s3com_chunk%npoints, s3com_chunk%nmeas), source = 0._wp)
      
      allocate(s3com_chunk%jac%p(s3com_chunk%npoints, s3com_chunk%nmeas, s3com_chunk%nlevels), source = 0._wp)
      allocate(s3com_chunk%jac%t(s3com_chunk%npoints, s3com_chunk%nmeas, s3com_chunk%nlevels), source = 0._wp)
      allocate(s3com_chunk%jac%q(s3com_chunk%npoints, s3com_chunk%nmeas, s3com_chunk%nlevels), source = 0._wp)
      
      allocate(s3com_chunk%jac%cfrac(s3com_chunk%npoints, s3com_chunk%nmeas, s3com_chunk%nlayers), source = 0._wp)
      allocate(s3com_chunk%jac%clwde(s3com_chunk%npoints, s3com_chunk%nmeas, s3com_chunk%nlayers), source = 0._wp)
      ! -------------------------------------------------------------------------------------------------------------------------
      
      ! Jacobian data (TL model)
      ! -------------------------------------------------------------------------------------------------------------------------
      s3com_chunk%k_tl%do_k_tl_calc = s3com%nml%do_k_tl_calc
      
      allocate(s3com_chunk%k_tl%t(s3com_chunk%npoints, s3com_chunk%nmeas, s3com_chunk%nlevels), source = 0._wp)
      allocate(s3com_chunk%k_tl%q(s3com_chunk%npoints, s3com_chunk%nmeas, s3com_chunk%nlevels), source = 0._wp)
      ! -------------------------------------------------------------------------------------------------------------------------
      
   end subroutine s3com_subset
   ! ============================================================================================================================
   
   ! ============================================================================================================================
   !> @brief Update
   !! @param[inout] s3com         s3com structure
   !! @param[inout] s3com_chunk   s3com structure
   subroutine s3com_update(s3com, s3com_chunk)
      
      ! Input/output
      type(type_s3com), intent(inout) :: s3com
      type(type_s3com), intent(inout) :: s3com_chunk
      
      ! Internal
      integer(wpi) :: idx_start, idx_end, npoints, nlevels, nlayers, nmeas, nstates
      
      idx_start = s3com_chunk%idx_start
      idx_end   = s3com_chunk%idx_end
      npoints   = s3com_chunk%npoints
      nlevels   = s3com_chunk%nlevels
      nlayers   = s3com_chunk%nlayers
      nmeas     = s3com_chunk%nmeas
      nstates   = s3com_chunk%nstates
      
      ! Radiation data
      ! -------------------------------------------------------------------------------------------------------------------------
      s3com%rad%f_rad_total(idx_start:idx_end, 1:nmeas) = s3com_chunk%rad%f_rad_total(1:npoints, 1:nmeas)
      s3com%rad%f_rad_clear(idx_start:idx_end, 1:nmeas) = s3com_chunk%rad%f_rad_clear(1:npoints, 1:nmeas)
      s3com%rad%f_ref_total(idx_start:idx_end, 1:nmeas) = s3com_chunk%rad%f_ref_total(1:npoints, 1:nmeas)
      s3com%rad%f_ref_clear(idx_start:idx_end, 1:nmeas) = s3com_chunk%rad%f_ref_clear(1:npoints, 1:nmeas)
      s3com%rad%f_bt_total(idx_start:idx_end, 1:nmeas)  = s3com_chunk%rad%f_bt_total(1:npoints, 1:nmeas)
      s3com%rad%f_bt_clear(idx_start:idx_end, 1:nmeas)  = s3com_chunk%rad%f_bt_clear(1:npoints, 1:nmeas)
      s3com%rad%emiss(idx_start:idx_end, 1:nmeas)       = s3com_chunk%rad%emiss(1:npoints, 1:nmeas)
      s3com%rad%brdf(idx_start:idx_end, 1:nmeas)        = s3com_chunk%rad%brdf(1:npoints, 1:nmeas)
      
      deallocate(s3com_chunk%rad%f_rad_total, s3com_chunk%rad%f_rad_clear, s3com_chunk%rad%f_ref_total, &
                 s3com_chunk%rad%f_ref_clear, s3com_chunk%rad%f_bt_total, s3com_chunk%rad%f_bt_clear, s3com_chunk%rad%emiss, &
                 s3com_chunk%rad%brdf)
      ! -------------------------------------------------------------------------------------------------------------------------
      
      ! Jacobian data (K model)
      ! -------------------------------------------------------------------------------------------------------------------------
      s3com%jac%brdf(idx_start:idx_end, 1:nmeas)  = s3com_chunk%jac%brdf(1:npoints, 1:nmeas)
      s3com%jac%emiss(idx_start:idx_end, 1:nmeas) = s3com_chunk%jac%emiss(1:npoints, 1:nmeas)
      
      s3com%jac%p(idx_start:idx_end, 1:nmeas, 1:nlevels) = s3com_chunk%jac%p(1:npoints, 1:nmeas, 1:nlevels)
      s3com%jac%t(idx_start:idx_end, 1:nmeas, 1:nlevels) = s3com_chunk%jac%t(1:npoints, 1:nmeas, 1:nlevels)
      s3com%jac%q(idx_start:idx_end, 1:nmeas, 1:nlevels) = s3com_chunk%jac%q(1:npoints, 1:nmeas, 1:nlevels)
      
      s3com%jac%cfrac(idx_start:idx_end, 1:nmeas, 1:nlayers) = s3com_chunk%jac%cfrac(1:npoints, 1:nmeas, 1:nlayers)
      s3com%jac%clwde(idx_start:idx_end, 1:nmeas, 1:nlayers) = s3com_chunk%jac%clwde(1:npoints, 1:nmeas, 1:nlayers)
      
      s3com%jac%do_jacobian_calc = s3com_chunk%jac%do_jacobian_calc
      
      deallocate(s3com_chunk%jac%brdf, s3com_chunk%jac%emiss, s3com_chunk%jac%p, s3com_chunk%jac%t, s3com_chunk%jac%q, &
                 s3com_chunk%jac%cfrac, s3com_chunk%jac%clwde)
      ! -------------------------------------------------------------------------------------------------------------------------
      
      ! Jacobian data (TL model)
      ! -------------------------------------------------------------------------------------------------------------------------
      s3com%k_tl%t(idx_start:idx_end, 1:nmeas, 1:nlevels) = s3com_chunk%k_tl%t(1:npoints, 1:nmeas, 1:nlevels)
      s3com%k_tl%q(idx_start:idx_end, 1:nmeas, 1:nlevels) = s3com_chunk%k_tl%q(1:npoints, 1:nmeas, 1:nlevels)
      
      s3com%k_tl%do_k_tl_calc = s3com_chunk%k_tl%do_k_tl_calc
      
      deallocate(s3com_chunk%k_tl%t, s3com_chunk%k_tl%q)
      ! -------------------------------------------------------------------------------------------------------------------------
      
      ! General
      ! -------------------------------------------------------------------------------------------------------------------------
      deallocate(s3com_chunk%flag_rttov)
      ! -------------------------------------------------------------------------------------------------------------------------
      
   end subroutine s3com_update
   ! ============================================================================================================================
   
end module mod_s3com_setup

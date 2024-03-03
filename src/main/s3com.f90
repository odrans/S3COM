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

program s3com_main
   
   use s3com_types,         only: wp, type_rttov_opt, type_nml, type_model, type_s3com, type_cld
   use s3com_config,        only: nstates, apriori_cod, apriori_reff
   use mod_io_namelist,     only: namelist_load
   use mod_rttov_atlas,     only: rttov_atlas_init, rttov_atlas_deinit
   use mod_rttov_opts,      only: rttov_opts_init
   use mod_rttov_coefs,     only: rttov_coefs_init, rttov_coefs_deinit
   use mod_rttov_setup,     only: rttov_setup_opt, rttov_setup_atm
   use mod_rttov,           only: run_rttov
   use mod_s3com_setup,     only: s3com_init, s3com_subset, s3com_update
   use mod_models,          only: models_load, models_free
   use mod_utils_math,      only: n_chunks
   use mod_io_s3com,        only: write_output
   use mod_rttov_utils,     only: find_idx_rttov, find_ret_idx_rttov
   use mod_cld_mie,         only: cld_mie_load
   use mod_ret_cloud_model, only: init_cloud_z, init_cloud_prof
   use mod_ret_utils,       only: ret_init, ret_free, idx_liq
   use mod_ret_run,         only: ret_run
   
   implicit none
   
   type(type_model)     :: model, rttov_atm, rttov_atm_ret
   type(type_rttov_opt) :: rttov_opt
   type(type_s3com)     :: s3com, s3com_chunk
   type(type_nml)       :: nml
   type(type_cld)       :: cld
   
   real(kind=wp) :: zenangle, azangle
   
   integer(kind=4), dimension(:), allocatable :: idx_lwp, idx_ret
   integer(kind=4) :: nchunks, idx_start, idx_end, ichunk
   integer(kind=4) :: nlayers, npoints
   
   ! Read namelist file
   nml = namelist_load()
   
   ! Temporary: setting the viewing satellite angles
   zenangle = 0._wp; azangle = 0._wp
   
   ! Load atmospheric data from selected models (`model` created)
   call models_load(nml, model)
   
   ! Initialise the s3com structure
   call s3com_init(nml, model, s3com)
   npoints = s3com%npoints
   nlayers = s3com%nlayers
   
   ! Setup the RTTOV options (`rttov_opt` created)
   call rttov_setup_opt(nml, zenangle, azangle, rttov_opt, s3com)
   
   ! Initialize RTTOV (loading surface data and optical properties)
   ! Structures are created and allocated internally, no outputs here
   call rttov_opts_init(rttov_opt, s3com)
   call rttov_coefs_init(rttov_opt, s3com)
   call rttov_atlas_init(rttov_opt, s3com)
   
   ! Load cloud optical properties
   if (rttov_opt%user_cld_opt_param) then
      call cld_mie_load(s3com, cld)
   endif
   
   ! Loop over data points, by chunks
   nChunks = n_chunks(npoints, nml%npoints_it)
   
   loop_chunks : do iChunk = 1, nChunks
      
      write(6,*) ichunk, "/", nchunks
      
      ! Compute starting and end points for the chunk
      if (nChunks .EQ. 1) then
         idx_start = 1; idx_end = npoints
      else
         idx_start = (iChunk-1)* nml%npoints_it + 1
         idx_end = iChunk * nml%npoints_it
         if (idx_end .gt. npoints) idx_end = npoints
      end if
      
      ! Subset the atmosphere for RTTOV
      call rttov_setup_atm(idx_start, idx_end, model, rttov_atm)
      
      ! Subset the relevant part of the s3com structure
      call s3com_subset(idx_start, idx_end, s3com, s3com_chunk)
      
      ! Run RTTOV
      s3com_chunk%flag_rttov(:) = .true.
      call run_rttov(rttov_atm, rttov_opt, s3com_chunk, cld)
      
      !
      call s3com_update(s3com, s3com_chunk)
      
      ! Calculate the radiance corresponding to the "true synthetic" cloud
      s3com%ret%Y = s3com%rad%f_rad_total
      
      ! Retrievals
      if (nml%flag_retrievals) then
         
         ! Initialize the retrieval arrays
         call ret_init(s3com)
         
         ! Subset the atmosphere for RTTOV
         call rttov_setup_atm(idx_start, idx_end, model, rttov_atm_ret)
         
         ! Extract cloud top and base altitudes for the liquid phase only
         call init_cloud_z(rttov_atm_ret, s3com)
         
         ! Set conditions to use RTTOV for the liquid clouds only
         idx_lwp = idx_liq(s3com, .true.)
         s3com%ret%flag_rttov(idx_lwp) = .true.
         
         ! Calculate the corresponding index, i.e., the number of pixels with liquid clouds only
         idx_ret = find_ret_idx_rttov(s3com)
         
         if (size(idx_ret) .eq. 0) then
            cycle
         endif
         
         ! Define the apriori values
         s3com%ret%Xa(:,1) = apriori_cod  ! Xa(1): COD
         s3com%ret%Xa(:,2) = apriori_reff ! Xa(2): reff - Deff maximum is 52 um when using "Deff scheme"
         
         ! Set Xi = Xa
         s3com%ret%Xi(:,:) = s3com%ret%Xa(:,:)
         
         ! Update the cloud profiles in rttov_atm_ret using parameters from S3COM
         call init_cloud_prof(rttov_atm_ret, s3com)
         
         ! Run RTTOV
         s3com%ret%flag_ret = nml%flag_retrievals
         call run_rttov(rttov_atm_ret, rttov_opt, s3com, cld)
         
         ! Perform the retrievals
         call ret_run(rttov_atm_ret, rttov_opt, s3com, cld)
         
         ! Free the retrieval arrays
         call ret_free(s3com)
         
      endif
      
   enddo loop_chunks
   
   ! Write output file
   call write_output(s3com, model, nml)
   
   ! Deallocate arrays
   call models_free(model)
   
   call rttov_coefs_deinit()
   call rttov_atlas_deinit()
   
end program s3com_main

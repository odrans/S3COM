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

!> Initialize RTTOV atlas structures
module mod_rttov_atlas
   
   use s3com_types,          only: wp, wpi, type_rttov_opt, type_nml, type_s3com
   use mod_rttov_utils,      only: check_rttov_status
   use rttov_const,          only: sensor_id_mw, sensor_id_po
   use mod_rttov_opts,       only: opts
   use mod_rttov_coefs,      only: coefs
   use mod_rttov_emis_atlas, only: rttov_emis_atlas_data, atlas_type_ir, atlas_type_mw
   use mod_rttov_brdf_atlas, only: rttov_brdf_atlas_data
   
   implicit none
   
   private
   public :: rttov_atlas_init, rttov_atlas_deinit, emis_atlas, brdf_atlas

! Use emissivity atlas
#include "rttov_setup_emis_atlas.interface"
#include "rttov_get_emis.interface"
#include "rttov_deallocate_emis_atlas.interface"

! Use BRDF atlas
#include "rttov_setup_brdf_atlas.interface"
#include "rttov_get_brdf.interface"
#include "rttov_deallocate_brdf_atlas.interface"

   type(rttov_emis_atlas_data) :: emis_atlas !< Data structure for emissivity atlas
   type(rttov_brdf_atlas_data) :: brdf_atlas !< Data structure for BRDF atlas
   
contains
   
   ! ============================================================================================================================
   !> @brief Initialize the emissivity and BRDF atlas
   !! @param[in] rttov_opt   rttov options structure
   !! @param[inout] s3com    s3com data structure
   subroutine rttov_atlas_init(rttov_opt, s3com)
      
      ! Input
      type(type_rttov_opt), intent(in) :: rttov_opt
      
      ! Input/output
      type(type_s3com), intent(inout) :: s3com
      
      ! Internal
      character(len=256) :: path_emis_atlas, path_brdf_atlas
      integer(wpi) :: errorstatus, imonth, atlas_type
      
      imonth = rttov_opt%month
      
      path_emis_atlas = trim(s3com%nml%path_rttov)//'/emis_data'
      path_brdf_atlas = trim(s3com%nml%path_rttov)//'/brdf_data'
      
      ! Initialize the RTTOV emissivity atlas
      ! This loads the default IR/MW atlases: use the atlas_id argument to select alternative atlases
      if (coefs%coef%id_sensor == sensor_id_mw .or. coefs%coef%id_sensor == sensor_id_po) then
         atlas_type = atlas_type_mw !< MW atlas
      else
         atlas_type = atlas_type_ir !< IR atlas
      endif
      
      call rttov_setup_emis_atlas(  &
           errorstatus,             &
           opts,                    &
           imonth,                  &
           atlas_type,              & !< Selects MW (1) or IR (2)
           emis_atlas,              &
           path  = path_emis_atlas, &
           coefs = coefs)
      
      call check_rttov_status(errorstatus, 'rttov_setup_emis_atlas')
      
      ! Initialise the RTTOV BRDF atlas
      if (opts%rt_ir%addsolar) then
         call rttov_setup_brdf_atlas(  &
              errorstatus,             &
              opts,                    &
              imonth,                  &
              brdf_atlas,              &
              path  = path_brdf_atlas, &
              coefs = coefs)
         
         call check_rttov_status(errorstatus, 'rttov_setup_brdf_atlas')
      endif
      
   end subroutine rttov_atlas_init
   ! ============================================================================================================================
   
   ! ============================================================================================================================
   !> @brief Deallocate the emissivity and BRDF atlas
   subroutine rttov_atlas_deinit()
      
      ! Deallocate emissivity atlas
      call rttov_deallocate_emis_atlas(emis_atlas)
      
      ! Deallocate BRDF atlas
      if (opts%rt_ir%addsolar) call rttov_deallocate_brdf_atlas(brdf_atlas)
      
   end subroutine rttov_atlas_deinit
   ! ============================================================================================================================
   
end module mod_rttov_atlas

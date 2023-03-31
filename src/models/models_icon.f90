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

module mod_icon
   
   use s3com_types,  only: wp, type_icon
   use mod_io_icon,  only: icon_read
   use mod_io_utils, only: extract_coordinates
   use s3com_config, only: rd, rv, epsilon, mu, nu, a, b, Q_ext, rholiq, lwp_lay_threshold
   
   implicit none
   
   private
   public :: icon_load, icon_clear
   
   contains
      
      subroutine icon_load(fname, icon)
         
         !! Inputs
         character(len=256), intent(in) :: fname
         
         !! Inputs/Outputs
         type(type_icon), intent(inout) :: icon
         
         !! Internal
         integer(kind=4) :: nlayers, npoints
         
         !! Extract the number of vertical layers and grid points in the input files
         call extract_coordinates(fname, nlayers, npoints)
         
         !! Initialize the ICON array
         call icon_init(npoints, nlayers, icon)
         
         !! Read input netcdf file containing ICON outputs
         call icon_read(fname, icon)
         
         !! Post-process ICON data
         call icon_process(icon)
         
      end subroutine icon_load
      
      subroutine icon_process(icon)
         
         type(type_icon), intent(inout) :: icon
         
         !! Internal
         integer(kind=4) :: i, j, nlevels, nlayers, npoints
         
         nlevels = icon%nlevels
         nlayers = icon%nlayers
         npoints = icon%npoints
         
         !! Atmospheric structure
         !! ---------------------------------------------------------------------------------------------------------------------
         !! Layer depth (approximate!! why not ifc(i+1) - ifc(i)?)
         icon%dz(:,1:nlayers) = abs(icon%z_ifc(:,1:nlayers) - icon%z(:,1:nlayers)) * 2._wp
         
         !! Atmospheric moisture
         !! ----------------------------------------------------------------------------------------------------
         !! Virtual temperature
         icon%tv = icon%t * (1 + 0.608 * icon%q)
         
         !! Moist air density (kg/m3)
         icon%rho = icon%p / (rd * icon%tv)
         
         !! Saturation vapour pressure of liquid water (Pa); Murphy and Koop (2005)
         !icon%es_w = exp(54.842763 - 6763.22 / icon%t - 4.210 * log(icon%t) + 0.000367 * icon%t + &
         !            tanh(0.0415 * (icon%t - 218.8)) * (53.878 - 1331.22 / icon%t - 9.44523 * LOG(icon%t) + 0.014025 * icon%t))
         
         !! Saturation vapour pressure of ice water (Pa); Murphy and Koop (2005)
         !icon%es_i = exp(9.550426 - 5723.265 / icon%t + 3.53068 * log(icon%t) - 0.00728332 * icon%t)
         
         !! Cloud properties
         !! ---------------------------------------------------------------------------------------------------------------------
         !! Cloud liquid water content (kg/m3)
         icon%lwc = icon%clw * icon%rho
         
         !! Cloud ice water content (kg/m3)
         icon%iwc = icon%cli * icon%rho
         
         !! Cloud liquid droplet number concentration (particules/m3)
         icon%cdnc = icon%qnc * icon%rho
         
         !! Cloud liquid water effective radius (m)
         do i = 1, npoints
            do j = 1, nlayers
               if(icon%cdnc(i,j) .gt. 0) then
                  icon%reff(i,j) = (a/2._wp)*(gamma((nu+1._wp+3._wp*b)/mu)/gamma((nu+1._wp+2._wp*b)/mu))*(icon%lwc(i,j) / &
                                   icon%cdnc(i,j))**b*(gamma((nu+1._wp)/mu)/gamma((nu+2._wp)/mu))**b
               endif
            enddo
         enddo
         
         icon%reff = icon%reff * 1E6 !m to um (default input in RTTOV)
         
         !! Cloud liquid droplet extinction coefficient (m-1)
         do i = 1, npoints
            do j = 1, nlayers
               if(icon%reff(i,j) .gt. 0) then
                  icon%beta_ext(i,j) = (3._wp/4._wp)*(Q_ext/rholiq)*(icon%lwc(i,j)/(icon%reff(i,j)*1E-6))
               endif
            enddo
         enddo
         
         !! Cloud optical depth
         do i = 1, npoints
            icon%cod(i) = 0._wp
            do j = 1, nlayers
               icon%cod(i) = icon%cod(i) + icon%beta_ext(i,j)*icon%dz(i,j)
            enddo
         enddo
         
         !! Cloud liquid water effective radius at cloud top (um)
         do i = 1, npoints
            do j = 1, nlayers
               if (icon%lwc(i,j)*1E3*icon%dz(i,j) .gt. lwp_lay_threshold .and. icon%cod(i) .gt. 5._wp) then
                  icon%ztop_liq_idx(i) = j
                  icon%reff_top(i) = icon%reff(i,icon%ztop_liq_idx(i))
                  !write(*,*) icon%ztop_liq_idx(i), icon%reff_top(i)
                  exit
               endif
            enddo
         enddo
         
      end subroutine icon_process
      
      subroutine icon_init(npoints, nlayers, icon)
         
         type(type_icon) :: icon
         integer(kind=4), intent(in) :: npoints, nlayers
         integer(kind=4) :: nlevels
         
         nlevels = nlayers + 1
         
         icon%npoints = npoints
         icon%nlayers = nlayers
         icon%nlevels = nlevels
         
         allocate(icon%height(nlayers), source = 0)
         allocate(icon%height_2(nlevels), source = 0)
         allocate(icon%ztop_liq_idx(npoints), source = 0)
         
         !! 2D variables
         allocate(icon%lon(npoints), source = 0._wp)
         allocate(icon%lat, icon%lon_orig, icon%lat_orig, icon%topography, icon%landmask, &
                  icon%ps, icon%ts, icon%t_2m, icon%q_2m, icon%u_10m, icon%v_10m, icon%cod, icon%reff_top, &
                  mold = icon%lon)
         
         !! 3D variables on atmospheric levels
         allocate(icon%z_ifc(npoints, nlevels), source = 0._wp)
         allocate(icon%p_ifc, icon%t_ifc, icon%q_ifc, &
                  mold = icon%z_ifc)
         
         !! 3D variables in atmospheric layers
         allocate(icon%z(npoints, nlayers), source = 0._wp)
         allocate(icon%p, icon%t, icon%q, icon%clc, icon%clw, icon%cli, &
                  icon%qnc, icon%qr, icon%qs, icon%dz, icon%rho, icon%tv, icon%lwc, &
                  icon%iwc, icon%cdnc, icon%reff, icon%beta_ext, &
                  mold = icon%z)
         
      end subroutine icon_init
      
      subroutine icon_clear(icon)
         
         type(type_icon), intent(inout) :: icon
         
         deallocate(icon%height, icon%height_2, icon%lon, icon%lat, icon%lon_orig, icon%lat_orig, &
                    icon%topography, icon%landmask, icon%ps, icon%ts, icon%t_2m, icon%q_2m, icon%u_10m, icon%v_10m, icon%cod, &
                    icon%p, icon%z, icon%z_ifc, icon%p_ifc, icon%t_ifc, icon%q_ifc, &
                    icon%t, icon%q, icon%clc, icon%clw, icon%cli, icon%qnc, icon%qr, icon%qs, icon%dz, &
                    icon%rho, icon%tv, icon%lwc, icon%iwc, icon%cdnc, icon%reff, icon%beta_ext, icon%ztop_liq_idx, icon%reff_top)
         
      end subroutine icon_clear
      
end module mod_icon

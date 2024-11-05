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

module mod_io_icon
   
   use netcdf
   use s3com_types,    only: wp, dp, wpi, type_icon
   use mod_utils_fort, only: s3com_error
   use mod_io_utils,   only: map_ll_to_point, check_netcdf_status
   
   implicit none
   
   private
   public :: icon_read
   
contains
   
   ! ============================================================================================================================
   !> @brief Import ICON simulation data from the ICON NetCDF file
   !! @param[in] fname     ICON NetCDF file name
   !! @param[inout] icon   ICON data structure
   subroutine icon_read(fname, icon)
      
      ! Parameters
      character(len=64), parameter :: routine_name = 'READ_ICON'
      integer, parameter :: NMAX_DIM = 5
      
      ! Input
      character(len=256), intent(in) :: fname
      
      ! Input/output
      type(type_icon), intent(inout) :: icon
      
      ! Internal
      character(len=256) :: errmsg
      character(len=256) :: dimname(NMAX_DIM), vname
      
      integer(wpi) :: idim, dimsize(NMAX_DIM), vdimid(NMAX_DIM)
      integer(wpi) :: ncid, ndims, nvars, ngatts, recdim, errst, vid, vrank
      integer(wpi) :: dim1, dim2, dim3, dim4, dim5, i, j, k
      integer(wpi) :: npoints, nlevels, nlayers
      integer, dimension(:), allocatable :: plon, plat, point

      real(wp), dimension(:), allocatable :: lat, lon, ll
      real(wp), allocatable :: x1(:), x2(:,:), x3(:,:,:), x4(:,:,:,:), x5(:,:,:,:,:) !< temporary arrays
      real(dp), dimension(1) :: x0
      
      logical :: Llat, Llon, Lpoint
      
      ! Get the number of points, atmospheric levels and layers from the ICON data structure
      npoints = icon%npoints; nlevels = icon%nlevels; nlayers = icon%nlayers

      ! Open the ICON NetCDF file
      errst = nf90_open(fname, nf90_nowrite, ncid)
      call check_netcdf_status(errst, 'nf90_open')
      
      ! Get the number of dimensions, variables, global attributes, and the unlimitted dimension ID
      errst = nf90_inquire(ncid, ndims, nvars, ngatts, recdim)
      call check_netcdf_status(errst, 'nf90_inquire')
      
      Llat = .false.; Llon = .false.; Lpoint = .false.

      ! Get the dimension IDs and names
      ! -------------------------------------------------------------------------------------------------------------------------
      do idim = 1,ndims

         errst = nf90_Inquire_Dimension(ncid, idim, name=dimname(idim), len=dimsize(idim))
         call check_netcdf_status(errst, 'nf90_Inquire_Dimension')

         if (trim(dimname(idim)) .eq. 'point') then
            Lpoint = .true.
            if (npoints /= dimsize(idim)) then
               errmsg = 'Number of points selected is greater than in input file '//trim(fname)
               call s3com_error(routine_name,errmsg)
            endif
            icon%nlon = npoints
            icon%nlat = npoints
         endif
         
         if (trim(dimname(idim)) .eq. 'lon') then
            Llon = .true.; icon%nlon = dimsize(idim)
         endif
         
         if (trim(dimname(idim)) .eq. 'lat') then
            Llat = .true.; icon%nlat = dimsize(idim)
         endif
         
      enddo

      allocate(point(icon%npoints), lon(icon%nlon), lat(icon%nlat))
      ! --------------------------------------------------------------------------------------------------------------------------

      ! Extract coordinates
      ! -------------------------------------------------------------------------------------------------------------------------
      if (Lpoint) then ! 1D mode
         icon%nlon = npoints
         icon%nlat = npoints
         icon%mode = 1
      else if (Llon .and. Llat) then ! 2D mode
         if (npoints /= icon%nlon * icon%nlat) then
            errmsg = 'Number of points selected is different from indicated in input file '//trim(fname)
            call s3com_error(routine_name,errmsg)
         endif
         icon%mode = 2 ! Don't know yet if (lon, lat) or (lat, lon) at this point
      else
         errmsg = trim(fname)//' file contains wrong dimensions'
         call s3com_error(routine_name,errmsg)
      endif

      if (Lpoint) then
         errst = nf90_inq_varid(ncid, 'point', vid)
         call check_netcdf_status(errst, 'nf90_inq_varid')
         errst = nf90_get_var(ncid, vid, point, start = (/1/), count = (/icon%npoints/))
         call check_netcdf_status(errst, 'nf90_get_var')
      end if

      errst = nf90_inq_varid(ncid, 'lon', vid)
      call check_netcdf_status(errst, 'nf90_inq_varid')
      errst = nf90_get_var(ncid, vid, lon, start = (/1/), count = (/icon%nlon/))
      call check_netcdf_status(errst, 'nf90_get_var')
      
      errst = nf90_inq_varid(ncid, 'lat', vid)
      call check_netcdf_status(errst, 'nf90_inq_varid')
      errst = nf90_get_var(ncid, vid, lat, start = (/1/), count = (/icon%nlat/))
      call check_netcdf_status(errst, 'nf90_get_var')
      
      errst = nf90_inq_varid(ncid, 'height', vid)
      call check_netcdf_status(errst, 'nf90_inq_varid')
      errst = nf90_get_var(ncid, vid, icon%height, start = (/1/), count = (/icon%nlayers/))
      call check_netcdf_status(errst, 'nf90_get_var')
      
      errst = nf90_inq_varid(ncid, 'height_2', vid)
      call check_netcdf_status(errst, 'nf90_inq_varid')
      errst = nf90_get_var(ncid, vid, icon%height_2, start = (/1/), count = (/icon%nlevels/))
      call check_netcdf_status(errst, 'nf90_get_var')

      icon%lon_orig(1:npoints) = lon
      icon%lat_orig(1:npoints) = lat
      icon%point_orig(1:npoints) = point
      ! -------------------------------------------------------------------------------------------------------------------------
      
      ! Extract all variables
      ! -------------------------------------------------------------------------------------------------------------------------
      do vid = 1, nvars
         
         vdimid = 0
         errst = nf90_Inquire_Variable(ncid, vid, name=vname, ndims=vrank, dimids=vdimid)
         call check_netcdf_status(errst, 'nf90_Inquire_Variable')
         
         ! Read in into temporary array of correct shape
         select case (vrank)
            case (0)
               errst = nf90_get_var(ncid, vid, x0, start=(/1/), count=(/1/))
            case (1)
               dim1 = dimsize(vdimid(1))
               allocate(x1(dim1))
               errst = nf90_get_var(ncid, vid, x1, start=(/1/), count=(/dim1/))
            case (2)
               dim1 = dimsize(vdimid(1))
               dim2 = dimsize(vdimid(2))
               allocate(x2(dim1,dim2))
               errst = nf90_get_var(ncid, vid, x2, start=(/1,1/), count=(/dim1,dim2/))
            case (3)
               dim1 = dimsize(vdimid(1))
               dim2 = dimsize(vdimid(2))
               dim3 = dimsize(vdimid(3))
            allocate(x3(dim1,dim2,dim3))
            errst = nf90_get_var(ncid, vid, x3, start=(/1,1,1/), count=(/dim1,dim2,dim3/))
            if ((icon%mode == 2) .or. (icon%mode == 3)) then
               if ((dim1 == icon%nlon) .and. (dim2 == icon%nlat)) then
                  icon%mode = 2
               else if ((dim1 == icon%nlat) .and. (dim2 == icon%nlon)) then
                  icon%mode = 3
               else
                  errmsg = 'Wrong mode for variable '//trim(vname)
                  call s3com_error(routine_name,errmsg)
               endif
            endif
            case (4)
               dim1 = dimsize(vdimid(1))
               dim2 = dimsize(vdimid(2))
               dim3 = dimsize(vdimid(3))
               dim4 = dimsize(vdimid(4))
               allocate(x4(dim1,dim2,dim3,dim4))
               errst = nf90_get_var(ncid, vid, x4, start=(/1,1,1,1/), count=(/dim1,dim2,dim3,dim4/))
            case (5)
               dim1 = dimsize(vdimid(1))
               dim2 = dimsize(vdimid(2))
               dim3 = dimsize(vdimid(3))
               dim4 = dimsize(vdimid(4))
               dim5 = dimsize(vdimid(5))
               allocate(x5(dim1,dim2,dim3,dim4,dim5))
               errst = nf90_get_var(ncid, vid, x5, start=(/1,1,1,1,1/), count=(/dim1,dim2,dim3,dim4,dim5/))
         end select
         
         call check_netcdf_status(errst, 'nf90_get_var')
         
         ! Map to the right input argument
         select case (trim(vname))
            
            ! 1D variable
            case('time')                                                           !< Time
               icon%time = x0(1)
            ! 2D variables
            case ('topography_c')                                                  !< Orography
               if (Lpoint) then
                  icon%topography(1:npoints) = x1(1:npoints)
               else
                  call map_ll_to_point(dim1,dim2,npoints,x2=x2,y1=icon%topography)
               endif
            case ('FR_LAND')                                                       !< Land use class fraction
               if (Lpoint) then
                  icon%landmask(1:npoints) = x1(1:npoints)
               else
                  call map_ll_to_point(dim1,dim2,npoints,x2=x2,y1=icon%landmask)
               endif
            case ('ps')                                                            !< Surface pressure
               if (Lpoint) then
                  icon%ps(1:npoints) = x1(1:npoints)
               else
                  call map_ll_to_point(dim1,dim2,npoints,x2=x2,y1=icon%ps)
               endif
            case ('t_s')                                                           !< Skin temperature
               if (Lpoint) then
                  icon%ts(1:npoints) = x1(1:npoints)
               else
                  call map_ll_to_point(dim1,dim2,npoints,x2=x2,y1=icon%ts)
               endif
            case ('tas')                                                           !< Temperature in 2m
               if (Lpoint) then
                  icon%t_2m(1:npoints) = x1(1:npoints)
               else
                  call map_ll_to_point(dim1,dim2,npoints,x2=x2,y1=icon%t_2m)
               endif
            case ('huss')                                                          !< Specific water vapour content in 2m
               if (Lpoint) then
                  icon%q_2m(1:npoints) = x1(1:npoints)
               else
                  call map_ll_to_point(dim1,dim2,npoints,x2=x2,y1=icon%q_2m)
               endif
            case ('u_10m')                                                         !< Zonal wind in 10 m
               if (Lpoint) then
                  icon%u_10m(1:npoints) = x1(1:npoints)
               else
                  call map_ll_to_point(dim1,dim2,npoints,x2=x2,y1=icon%u_10m)
               endif
            case ('v_10m')                                                         !< Meridional wind in 10 m
               if (Lpoint) then
                  icon%v_10m(1:npoints) = x1(1:npoints)
               else
                  call map_ll_to_point(dim1,dim2,npoints,x2=x2,y1=icon%v_10m)
               endif
            
            ! 3D variables
            case ('z_mc')                                                          !< Geometric height at full level center
               if (Lpoint) then
                  icon%z(1:npoints,:) = x2(1:npoints,1:nlayers)
               else
                  call map_ll_to_point(dim1,dim2,npoints,x3=x3,y2=icon%z)
               endif
            case ('z_ifc')                                                         !< Geometric height at half level center
               if (Lpoint) then
                  icon%z_ifc(1:npoints,1:nlevels) = x2(1:npoints,1:nlevels)
               else
                  call map_ll_to_point(dim1,dim2,npoints,x3=x3,y2=icon%z_ifc)
               endif
            case ('ta_ifc')                                                        !< Temperature at half level center
               if (Lpoint) then
                  icon%t_ifc(1:npoints,1:nlevels) = x2(1:npoints,1:nlevels)
               else
                  call map_ll_to_point(dim1,dim2,npoints,x3=x3,y2=icon%t_ifc)
               endif
            case ('pres_ifc')                                                      !< Pressures at half level center
               if (Lpoint) then
                  icon%p_ifc(1:npoints,1:nlevels) = x2(1:npoints,1:nlevels)
               else
                  call map_ll_to_point(dim1,dim2,npoints,x3=x3,y2=icon%p_ifc)
               endif
            case ('hus_ifc')                                                       !< Specific humidity at half level center
               if (Lpoint) then
                  icon%q_ifc(1:npoints,1:nlevels) = x2(1:npoints,1:nlevels)
               else
                  call map_ll_to_point(dim1,dim2,npoints,x3=x3,y2=icon%q_ifc)
               endif
            case ('pres')                                                          !< Air pressure at full level center
               if (Lpoint) then
                  icon%p(1:npoints,:) = x2(1:npoints,1:nlayers)
               else
                  call map_ll_to_point(dim1,dim2,npoints,x3=x3,y2=icon%p)
               endif
            case ('ta')                                                            !< Air temperature at full level center
               if (Lpoint) then
                  icon%t(1:npoints,:) = x2(1:npoints,1:nlayers)
               else
                  call map_ll_to_point(dim1,dim2,npoints,x3=x3,y2=icon%T)
               endif
            case ('hus')                                                           !< Specific humidity at full level center
               if (Lpoint) then
                  icon%q(1:npoints,:) = x2(1:npoints,1:nlayers)
               else
                  call map_ll_to_point(dim1,dim2,npoints,x3=x3,y2=icon%q)
               endif
            case ('clc')                                                           !< Cloud cover at full level center
               if (Lpoint) then
                  icon%clc(1:npoints,1:nlayers) = x2(1:npoints,1:nlayers)
               else
                  call map_ll_to_point(dim1,dim2,npoints,x3=x3,y2=icon%clc)
               endif
            case ('clw')                                                           !< Specific cloud water content at full level center
               if (Lpoint) then
                  icon%clw(1:npoints,:) = x2(1:npoints,1:nlayers)
               else
                  call map_ll_to_point(dim1,dim2,npoints,x3=x3,y2=icon%clw)
               endif
            case ('cli')                                                           !< Specific cloud ice content at full level center
               if (Lpoint) then
                  icon%cli(1:npoints,:) = x2(1:npoints,1:nlayers)
               else
                  call map_ll_to_point(dim1,dim2,npoints,x3=x3,y2=icon%cli)
               endif
            case ('qnc')                                                           !< Cloud droplet number concentration at full level center
               if (Lpoint) then
                  icon%qnc(1:npoints,:) = x2(1:npoints,1:nlayers)
               else
                  call map_ll_to_point(dim1,dim2,npoints,x3=x3,y2=icon%qnc)
               endif
            case ('qr')                                                            !< Rain mixing ratio at full level center
               if (Lpoint) then
                  icon%qr(1:npoints,:) = x2(1:npoints,1:nlayers)
               else
                  call map_ll_to_point(dim1,dim2,npoints,x3=x3,y2=icon%qr)
               endif
            case ('qs')                                                            !< Snow mixing ratio at full level center
               if (Lpoint) then
                  icon%qs(1:npoints,:) = x2(1:npoints,1:nlayers)
               else
                  call map_ll_to_point(dim1,dim2,npoints,x3=x3,y2=icon%qs)
               endif
            
         end select
         
         ! Free memory
         if (vrank == 1) deallocate(x1)
         if (vrank == 2) deallocate(x2)
         if (vrank == 3) deallocate(x3)
         if (vrank == 4) deallocate(x4)
         if (vrank == 5) deallocate(x5)
         
      enddo
      ! -------------------------------------------------------------------------------------------------------------------------
      
      ! Fill in the lat/lon vectors with the right values for 2D modes
      ! -------------------------------------------------------------------------------------------------------------------------
      ! This might be helpful if the inputs are 2D (gridded) and you want outputs in 1D mode
      allocate(plon(npoints),plat(npoints))
      allocate(ll(npoints))

      if (icon%mode == 1) then !< (point)
         ! lon and lat are already in the right format
          icon%lon(1:npoints) = lon(1:npoints)
          icon%lat(1:npoints) = lat(1:npoints)
       else if (icon%mode == 2) then !< (lon, lat)
         ll = lat
         do j = 1, dim2
            do i = 1, dim1
               k = (j - 1) * dim1 + i
               plon(k) = i
               plat(k) = j
            enddo
         enddo
         icon%lon(1:npoints) = lon(plon(1:npoints))
         icon%lat(1:npoints) = ll(plat(1:npoints))
      else if (icon%mode == 3) then !< (lat, lon)
         ll = lon
         do j = 1, dim2
            do i = 1, dim1
               k = (j - 1) * dim1 + i
               lon(k) = ll(j)
               lat(k) = lat(i)
            enddo
         enddo
         icon%lon(1:npoints) = ll(plon(1:npoints))
         icon%lat(1:npoints) = lat(plat(1:npoints))
      endif
      
      deallocate(plon,plat)
      deallocate(ll)
      
      errst = nf90_close(ncid)
      call check_netcdf_status(errst, 'nf90_close')
      ! -------------------------------------------------------------------------------------------------------------------------
      
   end subroutine icon_read
   ! ============================================================================================================================
   
end module mod_io_icon

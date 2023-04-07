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

module mod_io_nwpsaf

  use netcdf
  use s3com_types, only: wp, dp, type_nwpsaf
  use mod_utils_fort, only: s3com_error
  use mod_io_utils, only: map_ll_to_point, check_netcdf_status

  implicit none

  private
  public :: nwpsaf_read

contains

  !> @brief Import model simulation data from the NWPSAF NetCDF file
  !! @param[in] fname NetCDF file name
  !! @param[inout] nwpsaf NWPSAF data structure
  subroutine nwpsaf_read(fname, nwpsaf)

    !!Parameters
    character(len=64), parameter :: routine_name = 'READ_NWPSAF'
    integer, parameter :: NMAX_DIM = 5

    !!Inputs
    character(len=256), intent(in) :: fname

    !!Inputs/Outputs
    type(type_nwpsaf), intent(inout) :: nwpsaf

    !!Local variables
    character(len=256) :: errmsg
    character(len=256) :: vname, dimname(nmax_dim)

    integer(kind=4)                    :: dimsize(NMAX_DIM), vdimid(NMAX_DIM), idim
    integer(kind=4)                    :: ncid, ndims, nvars, ngatts, recdim, status, vid, vrank
    integer(kind=4)                    :: dim1, dim2, dim3, dim4, dim5, i
    integer(kind=4)                    :: npoints, nlevels, nlayers

    real(wp), dimension(:), allocatable :: lat, lon
    real(wp), allocatable               :: x1(:), x2(:,:), x3(:,:,:), x4(:,:,:,:), x5(:,:,:,:,:) ! Temporary arrays
    real(dp), dimension(1) :: x0

    logical :: Lpoint

    npoints = nwpsaf%npoints; nlevels = nwpsaf%nlevels; nlayers = nwpsaf%nlayers

    ! Open the file
    status = nf90_open(fname, nf90_nowrite, ncid)
    call check_netcdf_status(status, "nf90_open")

    ! Get the number of dimensions, variables, global attributes, and the unlimitted dimension ID
    status = nf90_inquire(ncid, ndims, nvars, ngatts, recdim)
    call check_netcdf_status(status, "nf90_inquire")

    ! Get the dimension IDs and sizes
    dimname(:) = ""
    dimsize(:) = 0
    do idim = 1,ndims
       status = nf90_Inquire_Dimension(ncid, idim, name=dimname(idim), len=dimsize(idim))
    enddo

    Lpoint = .true.

    nwpsaf%nlon = npoints
    nwpsaf%nlat = npoints
    nwpsaf%mode = 1

    allocate(lon(nwpsaf%nlon), lat(nwpsaf%nlat))

    ! Load all dimensions
    status = nf90_inq_varid(ncid, 'lon', vid)
    call check_netcdf_status(status, "nf90_inq_varid")
    status = nf90_get_var(ncid, vid, lon, start = (/1/), count = (/nwpsaf%nlon/))
    call check_netcdf_status(status, "nf90_get_var")

    status = nf90_inq_varid(ncid, 'lat', vid)
    call check_netcdf_status(status, "nf90_inq_varid")
    status = nf90_get_var(ncid, vid, lat, start = (/1/), count = (/nwpsaf%nlat/))
    call check_netcdf_status(status, "nf90_get_var")

    status = nf90_inq_varid(ncid, 'height', vid)
    call check_netcdf_status(status, "nf90_inq_varid")
    status = nf90_get_var(ncid, vid, nwpsaf%height, start = (/1/), count = (/nwpsaf%nlayers/))
    call check_netcdf_status(status, "nf90_get_var")

    status = nf90_inq_varid(ncid, 'height_2', vid)
    call check_netcdf_status(status, "nf90_inq_varid")
    status = nf90_get_var(ncid, vid, nwpsaf%height_2, start = (/1/), count = (/nwpsaf%nlevels/))
    call check_netcdf_status(status, "nf90_get_var")

    nwpsaf%lon = lon; nwpsaf%lat = lat
    nwpsaf%lon_orig = lon; nwpsaf%lat_orig = lat

    ! Convert longitudes to 0-360
    do i = 1, npoints
       if(nwpsaf%lon(i) < 0) nwpsaf%lon(i) = 360 + nwpsaf%lon(i)
    end do
    !!========================================================================================================================!!

    !!========================================================================================================================!!
    !! Extract all variables                                                                                                  !!
    !!========================================================================================================================!!
    do vid = 1, nvars

       vdimid = 0
       status = nf90_Inquire_Variable(ncid, vid, name=vname, ndims=vrank, dimids=vdimid)
       call check_netcdf_status(status, "nf90_Inquire_Variable")

       !! Read in into temporary array of correct shape
       select case (vrank)
        case (0)
          status = nf90_get_var(ncid, vid, x0, start=(/1/), count=(/1/))
       case (1)
          dim1 = dimsize(vdimid(1))
          allocate(x1(dim1))
          status = nf90_get_var(ncid, vid, x1, start=(/1/), count=(/dim1/))
       case (2)
          dim1 = dimsize(vdimid(1))
          dim2 = dimsize(vdimid(2))
          allocate(x2(dim1,dim2))
          status = nf90_get_var(ncid, vid, x2, start=(/1,1/), count=(/dim1,dim2/))
       case (3)
          dim1 = dimsize(vdimid(1))
          dim2 = dimsize(vdimid(2))
          dim3 = dimsize(vdimid(3))
          allocate(x3(dim1,dim2,dim3))
          status = nf90_get_var(ncid, vid, x3, start=(/1,1,1/), count=(/dim1,dim2,dim3/))
          if ((nwpsaf%mode == 2) .or. (nwpsaf%mode == 3)) then
             if ((dim1 == nwpsaf%nlon) .and. (dim2 == nwpsaf%nlat)) then
                nwpsaf%mode = 2
             else if ((dim1 == nwpsaf%nlat) .and. (dim2 == nwpsaf%nlon)) then
                nwpsaf%mode = 3
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
          status = nf90_get_var(ncid, vid, x4, start=(/1,1,1,1/), count=(/dim1,dim2,dim3,dim4/))
       case (5)
          dim1 = dimsize(vdimid(1))
          dim2 = dimsize(vdimid(2))
          dim3 = dimsize(vdimid(3))
          dim4 = dimsize(vdimid(4))
          dim5 = dimsize(vdimid(5))
          allocate(x5(dim1,dim2,dim3,dim4,dim5))
          status = nf90_get_var(ncid, vid, x5, start=(/1,1,1,1,1/), count=(/dim1,dim2,dim3,dim4,dim5/))
       end select
       call check_netcdf_status(status, "nf90_get_var")

       ! Map to the right input argument
       select case (trim(vname))

          ! 2D variables
       case ('point') ! point
          nwpsaf%point(1:npoints) = int(x1(1:npoints))
       case ('elevation') ! Orography
          nwpsaf%elevation(1:npoints) = x1(1:npoints)
       case ('lsm') ! Land-sea mask
          nwpsaf%lsm(1:npoints) = x1(1:npoints)
       case ('psurf') ! Surface pressure
          nwpsaf%psurf(1:npoints) = x1(1:npoints)
       case ('tsurf') ! Skin temperature
          nwpsaf%tsurf(1:npoints) = x1(1:npoints)
       case ('t2m') ! Temperature in 2m
          nwpsaf%t2m(1:npoints) = x1(1:npoints)
       case ('u10') ! Zonal wind in 10 m
          nwpsaf%u10(1:npoints) = x1(1:npoints)
       case ('v10') ! Meridional wind in 10 m
          nwpsaf%v10(1:npoints) = x1(1:npoints)
       case ('day') ! Meridional wind in 10 m
          nwpsaf%day(1:npoints) = int(x1(1:npoints))
       case ('month') ! Meridional wind in 10 m
          nwpsaf%month(1:npoints) = int(x1(1:npoints))
       case ('year') ! Meridional wind in 10 m
          nwpsaf%year(1:npoints) = int(x1(1:npoints))

          ! 3D variables
       case ('altitude') !Geometric height at full level center
          nwpsaf%altitude(1:npoints,:) = x2(1:npoints,1:nlayers)
       case ('altitudeh') !Geometric height at half level center
          nwpsaf%altitudeh(1:npoints,1:nlevels) = x2(1:npoints,1:nlevels)
       case ('temp') !Air temperature
          nwpsaf%temp(1:npoints,:) = x2(1:npoints,1:nlayers)
       case ('temph') !Temperature at half level center
          nwpsaf%temph(1:npoints,1:nlevels) = x2(1:npoints,1:nlevels)
       case ('hum') !Specific humidity
          nwpsaf%hum(1:npoints,:) = x2(1:npoints,1:nlayers)
       case ('humh') !Specific humidity at half level center
          nwpsaf%humh(1:npoints,1:nlevels) = x2(1:npoints,1:nlevels)
       case ('pap') !Air pressure
          nwpsaf%pap(1:npoints,:) = x2(1:npoints,1:nlayers)
       case ('paph') !Pressures at half level center
          nwpsaf%paph(1:npoints,1:nlevels) = x2(1:npoints,1:nlevels)
       case ('cc') !Cloud cover
          nwpsaf%cc(1:npoints,1:nlayers) = x2(1:npoints,1:nlayers)
       case ('lwc') !Specific cloud water content
          nwpsaf%lwc(1:npoints,:) = x2(1:npoints,1:nlayers)
       case ('iwc') !Specific cloud ice content
          nwpsaf%iwc(1:npoints,:) = x2(1:npoints,1:nlayers)
       end select

       ! Free memory
       if(vrank == 1) deallocate(x1)
       if(vrank == 2) deallocate(x2)
       if(vrank == 3) deallocate(x3)
       if(vrank == 4) deallocate(x4)
       if(vrank == 5) deallocate(x5)

    enddo

    !!=======================================================================================================================!!

    status = nf90_close(ncid)
    call check_netcdf_status(status, "nf90_close")

  end subroutine nwpsaf_read

end module mod_io_nwpsaf

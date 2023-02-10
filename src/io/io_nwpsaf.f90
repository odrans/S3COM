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
  use mod_io_utils, only: map_ll_to_point

  implicit none

  private
  public :: nwpsaf_read

  !!Types to be used as arrays of pointers
  type var1d
     character(LEN=16) :: name
     character(LEN=16) :: units
     integer :: dimsid(3)
     integer :: dimssz(2)
     integer :: vid
     logical :: lout
     real(wp), pointer, dimension(:) :: pntr
  end type var1d

  type var2d
     character(LEN=16) :: name
     character(LEN=16) :: units
     integer :: dimsid(4)
     integer :: dimssz(3)
     integer :: vid
     logical :: lout
     real(wp), pointer, dimension(:,:) :: pntr
  end type var2d

  type var3d
     character(LEN=16) :: name
     character(LEN=16) :: units
     integer :: dimsid(5)
     integer :: dimssz(4)
     integer :: vid
     logical :: lout
     real(wp), pointer, dimension(:,:,:) :: pntr
  end type var3d

contains

  subroutine nwpsaf_read(fname, nwpsaf)

    !!Parameters
    character(LEN=64), parameter :: routine_name = 'READ_NWPSAF'
    integer, parameter :: NMAX_DIM = 5

    !!Inputs
    character(LEN=256), intent(in) :: fname

    !!Inputs/Outputs
    type(type_nwpsaf), intent(inout) :: nwpsaf

    !!Local variables
    character(LEN=256) :: errmsg, straux
    character(LEN=256) :: dimname(NMAX_DIM), vname

    integer(KinD=4)                    :: idim, dimsize(NMAX_DIM), vdimid(NMAX_DIM)
    integer(KinD=4)                    :: ncid, ndims, nvars, ngatts, recdim, errst, vid, vrank
    integer(KinD=4)                    :: dim1, dim2, dim3, dim4, dim5, i, j, k
    integer(KinD=4)                    :: npoints, nlevels, nlayers
    integer, dimension(:), allocatable :: plon, plat

    real(wp), dimension(:), allocatable :: lat, lon, ll, height
    real(wp), allocatable               :: x1(:), x2(:,:), x3(:,:,:), x4(:,:,:,:), x5(:,:,:,:,:) ! Temporary arrays
    real(dp), dimension(1) :: x0

    logical :: Llat, Llon, Lpoint

    npoints = nwpsaf%npoints; nlevels = nwpsaf%nlevels; nlayers = nwpsaf%nlayers

    !!========================================================================================================================!!
    !! Checking the opening of the NWPSAF input NetCDF file                                                                     !!
    !!========================================================================================================================!!

    errst = nf90_open(fname, nf90_nowrite, ncid)
    if (errst/=0)  then
       errmsg = "Couldn't open "//trim(fname)
       call s3com_error(routine_name,errmsg)
    endif

    errst = nf90_inquire(ncid, ndims, nvars, ngatts, recdim)
    !!========================================================================================================================!!

    !!========================================================================================================================!!
    !! Checking the dimensions (track or lat-lon)                                                                             !!
    !!========================================================================================================================!!
    Lpoint = .true.
    !!========================================================================================================================!!

    !!========================================================================================================================!!
    !! Extract coordinates                                                                                                    !!
    !!========================================================================================================================!!
    nwpsaf%nlon = npoints
    nwpsaf%nlat = npoints
    nwpsaf%mode = 1

    allocate(lon(nwpsaf%nlon), lat(nwpsaf%nlat))

    errst = nf90_inq_varid(ncid, 'lon', vid)
    errst = nf90_get_var(ncid, vid, lon, start = (/1/), count = (/nwpsaf%nlon/))

    errst = nf90_inq_varid(ncid, 'lat', vid)
    errst = nf90_get_var(ncid, vid, lat, start = (/1/), count = (/nwpsaf%nlat/))

    errst = nf90_inq_varid(ncid, 'height', vid)
    errst = nf90_get_var(ncid, vid, nwpsaf%height, start = (/1/), count = (/nwpsaf%nlayers/))

    errst = nf90_inq_varid(ncid, 'height_2', vid)
    errst = nf90_get_var(ncid, vid, nwpsaf%height_2, start = (/1/), count = (/nwpsaf%nlevels/))

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
       errst = nf90_Inquire_Variable(ncid, vid, NAME=vname, ndims=vrank, dimids=vdimid)
       if (errst /= 0) then
          write(straux, *) vid
          errmsg = 'Error in nf90_Inquire_Variable, vid '//trim(straux)
          call s3com_error(routine_name,errmsg,errcode=errst)
       endif

       !!Read in into temporary array of correct shape
       if (vrank == 0) then
          errst = nf90_get_var(ncid, vid, x0, start=(/1/), count=(/1/))
       endif
       if (vrank == 1) then
          dim1 = dimsize(vdimid(1))
          allocate(x1(dim1))
          errst = nf90_get_var(ncid, vid, x1, start=(/1/), count=(/dim1/))
       endif
       if (vrank == 2) then
          dim1 = dimsize(vdimid(1))
          dim2 = dimsize(vdimid(2))
          allocate(x2(dim1,dim2))
          errst = nf90_get_var(ncid, vid, x2, start=(/1,1/), count=(/dim1,dim2/))
       endif
       if (vrank == 3) then
          dim1 = dimsize(vdimid(1))
          dim2 = dimsize(vdimid(2))
          dim3 = dimsize(vdimid(3))
          allocate(x3(dim1,dim2,dim3))
          errst = nf90_get_var(ncid, vid, x3, start=(/1,1,1/), count=(/dim1,dim2,dim3/))
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
       endif
       if (vrank == 4) then
          dim1 = dimsize(vdimid(1))
          dim2 = dimsize(vdimid(2))
          dim3 = dimsize(vdimid(3))
          dim4 = dimsize(vdimid(4))
          allocate(x4(dim1,dim2,dim3,dim4))
          errst = nf90_get_var(ncid, vid, x4, start=(/1,1,1,1/), count=(/dim1,dim2,dim3,dim4/))
       endif
       if (vrank == 5) then
          dim1 = dimsize(vdimid(1))
          dim2 = dimsize(vdimid(2))
          dim3 = dimsize(vdimid(3))
          dim4 = dimsize(vdimid(4))
          dim5 = dimsize(vdimid(5))
          allocate(x5(dim1,dim2,dim3,dim4,dim5))
          errst = nf90_get_var(ncid, vid, x5, start=(/1,1,1,1,1/), count=(/dim1,dim2,dim3,dim4,dim5/))
       endif
       if (errst /= 0) then
          write(straux, *)  vid
          errmsg = 'Error in nf90_get_var, vid '//trim(straux)
          call s3com_error(routine_name,errmsg,errcode=errst)
       endif

       ! Map to the right input argument
       select case (trim(vname))

          ! 2D variables
       case ('point') ! point
          nwpsaf%point(1:npoints) = x1(1:npoints)
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
          nwpsaf%day(1:npoints) = x1(1:npoints)
       case ('month') ! Meridional wind in 10 m
          nwpsaf%month(1:npoints) = x1(1:npoints)
       case ('year') ! Meridional wind in 10 m
          nwpsaf%year(1:npoints) = x1(1:npoints)

      ! 3D variables
       case ('altitude') !Geometric height at full level center
             nwpsaf%z(1:npoints,:) = x2(1:npoints,1:nlayers)
       case ('altitudeh') !Geometric height at half level center
             nwpsaf%z_ifc(1:npoints,1:nlevels) = x2(1:npoints,1:nlevels)
       case ('temph') !Temperature at half level center
             nwpsaf%t_ifc(1:npoints,1:nlevels) = x2(1:npoints,1:nlevels)
       case ('paph') !Pressures at half level center
             nwpsaf%p_ifc(1:npoints,1:nlevels) = x2(1:npoints,1:nlevels)
       case ('humh') !Specific humidity at half level center
             nwpsaf%q_ifc(1:npoints,1:nlevels) = x2(1:npoints,1:nlevels)
       case ('pap') !Air pressure
             nwpsaf%p(1:npoints,:) = x2(1:npoints,1:nlayers)
       case ('temp') !Air temperature
             nwpsaf%t(1:npoints,:) = x2(1:npoints,1:nlayers)
       case ('hum') !Specific humidity
             nwpsaf%q(1:npoints,:) = x2(1:npoints,1:nlayers)
       case ('cc') !Cloud cover
          nwpsaf%clc(1:npoints,1:nlayers) = x2(1:npoints,1:nlayers)
       case ('clw') !Specific cloud water content
          nwpsaf%clw(1:npoints,:) = x2(1:npoints,1:nlayers)
       case ('cli') !Specific cloud ice content
          nwpsaf%cli(1:npoints,:) = x2(1:npoints,1:nlayers)
       end select

       ! Free memory
       if(vrank == 1) deallocate(x1)
       if(vrank == 2) deallocate(x2)
       if(vrank == 3) deallocate(x3)
       if(vrank == 4) deallocate(x4)
       if(vrank == 5) deallocate(x5)

    enddo

    !!=======================================================================================================================!!

    errst = nf90_close(ncid)
    if (errst /= 0) then
       errmsg = "Error in nf90_close"
       call s3com_error(routine_name,errmsg,errcode=errst)
    endif


  end subroutine nwpsaf_read

end module mod_io_nwpsaf

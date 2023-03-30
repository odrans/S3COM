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
  use s3com_types, only: wp, dp, type_icon
  use mod_utils_fort, only: s3com_error
  use mod_io_utils, only: map_ll_to_point

  implicit none

  private
  public :: icon_read

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

  subroutine icon_read(fname, icon)

    !!Parameters
    character(LEN=64), parameter :: routine_name = 'READ_ICON'
    integer, parameter :: NMAX_DIM = 5

    !!Inputs
    character(LEN=256), intent(in) :: fname

    !!Inputs/Outputs
    type(type_icon), intent(inout) :: icon

    !!Local variables
    character(LEN=256) :: errmsg, straux
    character(LEN=256) :: dimname(NMAX_DIM), vname

    integer(KinD=4)                    :: idim, dimsize(NMAX_DIM), vdimid(NMAX_DIM)
    integer(KinD=4)                    :: ncid, ndims, nvars, ngatts, recdim, errst, vid, vrank
    integer(KinD=4)                    :: dim1, dim2, dim3, dim4, dim5, i, j, k
    integer(KinD=4)                    :: npoints, nlevels, nlayers
    integer, dimension(:), allocatable :: plon, plat

    real(wp), dimension(:), allocatable :: lat, lon, ll
    real(wp), allocatable               :: x1(:), x2(:,:), x3(:,:,:), x4(:,:,:,:), x5(:,:,:,:,:) ! Temporary arrays
    real(dp), dimension(1) :: x0

    logical :: Llat, Llon, Lpoint

    npoints = icon%npoints; nlevels = icon%nlevels; nlayers = icon%nlayers

    !!========================================================================================================================!!
    !! Checking the opening of the ICON input NetCDF file                                                                     !!
    !!========================================================================================================================!!

    errst = nf90_open(fname, nf90_nowrite, ncid)
    if (errst/=0)  then
       errmsg = "Couldn't open "//trim(fname)
       call s3com_error(routine_name,errmsg)
    endif

    !!========================================================================================================================!!

    !!========================================================================================================================!!
    !! Checking the dimensions (track or lat-lon)                                                                             !!
    !!========================================================================================================================!!

    errst = nf90_inquire(ncid, ndims, nvars, ngatts, recdim)
    if (errst /= 0) then
       errmsg = "Error in nf90_inquire"
       call s3com_error(routine_name, errmsg, errcode=errst)
    endif

    Llat = .false.; Llon = .false.; Lpoint = .false.

    do idim = 1,ndims
       errst = nf90_Inquire_Dimension(ncid, idim, NAME=dimname(idim), LEN=dimsize(idim))
       if (errst /= 0) then
          write(straux, *) idim
          errmsg = "Error in nf90_Inquire_Dimension, idim: "//trim(straux)
          call s3com_error(routine_name, errmsg)
       endif

       ! IF ((trim(dimname(idim)) .EQ. 'height') .AND. (nlevels > dimsize(idim))) THEN
       !    errmsg = 'Number of levels selected is greater than in input file '//trim(fname)
       !    CALL s3com_error(routine_name, errmsg)
       ! ENDIF

       if (trim(dimname(idim)) .eq. 'point') then
          Lpoint = .true.
          if (npoints /= dimsize(idim)) then
             errmsg = 'Number of points selected is greater than in input file '//trim(fname)
             call s3com_error(routine_name,errmsg)
          endif
       endif

       if (trim(dimname(idim)) .eq. 'lon') then
          Llon = .true.; icon%nlon = dimsize(idim)
       endif

       if (trim(dimname(idim)) .eq. 'lat') then
          Llat = .true.; icon%nlat = dimsize(idim)
       endif
    enddo

    allocate(lon(icon%nlon), lat(icon%nlat))

    !!========================================================================================================================!!

    !!========================================================================================================================!!
    !! Extract coordinates                                                                                                    !!
    !!========================================================================================================================!!

    if (Llon .and. Llat) then ! 2D mode
       if (npoints /= icon%nlon*icon%nlat) then
          errmsg = 'Number of points selected is different from indicated in input file '//trim(fname)
          call s3com_error(routine_name,errmsg)
       endif
       lon = -1.0E30; lat = -1.0E30
       icon%mode = 2 ! Don't know yet if (lon,lat) or (lat,lon) at this point
    else if (Lpoint) then ! 1D mode
       icon%nlon = npoints
       icon%nlat = npoints
       icon%mode = 1
    else
       errmsg = trim(fname)//' file contains wrong dimensions'
       call s3com_error(routine_name,errmsg)
    endif

    errst = nf90_inq_varid(ncid, 'lon', vid)
    if (errst /= 0) then
       errmsg = "Error in nf90_inq_varid, var: lon"
       call s3com_error(routine_name,errmsg,errcode=errst)
    endif

    errst = nf90_get_var(ncid, vid, lon, start = (/1/), count = (/icon%nlon/))
    if (errst /= 0) then
       errmsg = "Error in nf90_get_var, var: lon"
       call s3com_error(routine_name,errmsg,errcode=errst)
    endif

    errst = nf90_inq_varid(ncid, 'lat', vid)
    if (errst /= 0) then
       errmsg = "Error in nf90_inq_varid, var: lat"
       call s3com_error(routine_name,errmsg,errcode=errst)
    endif

    errst = nf90_get_var(ncid, vid, lat, start = (/1/), count = (/icon%nlat/))
    if (errst /= 0) then
       errmsg = "Error in nf90_get_var, var: lat"
       call s3com_error(routine_name,errmsg,errcode=errst)
    endif

    errst = nf90_inq_varid(ncid, 'height', vid)
    if (errst /= 0) then
       errmsg = "Error in nf90_inq_varid, var: height"
       call s3com_error(routine_name,errmsg,errcode=errst)
    endif
    
    errst = nf90_get_var(ncid, vid, icon%height, start = (/1/), count = (/icon%nlayers/))
    if (errst /= 0) then
       errmsg = "Error in nf90_get_var, var: height"
       call s3com_error(routine_name,errmsg,errcode=errst)
    endif

    errst = nf90_inq_varid(ncid, 'height_2', vid)
    if (errst /= 0) then
       errmsg = "Error in nf90_inq_varid, var: height_2"
       call s3com_error(routine_name,errmsg,errcode=errst)
    endif
    
    errst = nf90_get_var(ncid, vid, icon%height_2, start = (/1/), count = (/icon%nlevels/))
    if (errst /= 0) then
       errmsg = "Error in nf90_get_var, var: height_2"
       call s3com_error(routine_name,errmsg,errcode=errst)
    endif
    
    icon%lon_orig = lon; icon%lat_orig = lat

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

       !!Map to the right input argument
       select case (trim(vname))

          !! 1D variable
       case("time") ! Time
          icon%time = x0(1)
          !!2D variables
       case ('topography_c') !Orography
          if (Lpoint) then
             icon%topography(1:npoints) = x1(1:npoints)
          else
             call map_ll_to_point(dim1,dim2,npoints,x2=x2,y1=icon%topography)
          endif
       case ('FR_LAND') !Land use class fraction
          if (Lpoint) then
             icon%landmask(1:npoints) = x1(1:npoints)
          else
             call map_ll_to_point(dim1,dim2,npoints,x2=x2,y1=icon%landmask)
          endif
       case ('ps') !Surface pressure
          if (Lpoint) then
             icon%ps(1:npoints) = x1(1:npoints)
          else
             call map_ll_to_point(dim1,dim2,npoints,x2=x2,y1=icon%ps)
          endif
       case ('t_s') !Skin temperature
          if (Lpoint) then
             icon%ts(1:npoints) = x1(1:npoints)
          else
             call map_ll_to_point(dim1,dim2,npoints,x2=x2,y1=icon%ts)
          endif
       case ('tas') !Temperature in 2m
          if (Lpoint) then
             icon%t_2m(1:npoints) = x1(1:npoints)
          else
             call map_ll_to_point(dim1,dim2,npoints,x2=x2,y1=icon%t_2m)
          endif
       case ('huss') !Specific water vapour content in 2m
          if (Lpoint) then
             icon%q_2m(1:npoints) = x1(1:npoints)
          else
             call map_ll_to_point(dim1,dim2,npoints,x2=x2,y1=icon%q_2m)
          endif
       case ('u_10m') !Zonal wind in 10 m
          if (Lpoint) then
             icon%u_10m(1:npoints) = x1(1:npoints)
          else
             call map_ll_to_point(dim1,dim2,npoints,x2=x2,y1=icon%u_10m)
          endif
       case ('v_10m') !Meridional wind in 10 m
          if (Lpoint) then
             icon%v_10m(1:npoints) = x1(1:npoints)
          else
             call map_ll_to_point(dim1,dim2,npoints,x2=x2,y1=icon%v_10m)
          endif

          !!3D variables
       case ('z_mc') !Geometric height at full level center
          if (Lpoint) then
             icon%z(1:npoints,:) = x2(1:npoints,1:nlayers)
          else
             call map_ll_to_point(dim1,dim2,npoints,x3=x3,y2=icon%z)
          endif
       case ('z_ifc') !Geometric height at half level center
          if (Lpoint) then
             icon%z_ifc(1:npoints,1:nlevels) = x2(1:npoints,1:nlevels)
          else
             call map_ll_to_point(dim1,dim2,npoints,x3=x3,y2=icon%z_ifc)
          endif
       case ('ta_ifc') !Temperature at half level center
          if (Lpoint) then
             icon%t_ifc(1:npoints,1:nlevels) = x2(1:npoints,1:nlevels)
          else
             call map_ll_to_point(dim1,dim2,npoints,x3=x3,y2=icon%t_ifc)
          endif
       case ('pres_ifc') !Pressures at half level center
          if (Lpoint) then
             icon%p_ifc(1:npoints,1:nlevels) = x2(1:npoints,1:nlevels)
          else
             call map_ll_to_point(dim1,dim2,npoints,x3=x3,y2=icon%p_ifc)
          endif
       case ('hus_ifc') !Specific humidity at half level center
          if (Lpoint) then
             icon%q_ifc(1:npoints,1:nlevels) = x2(1:npoints,1:nlevels)
          else
             call map_ll_to_point(dim1,dim2,npoints,x3=x3,y2=icon%q_ifc)
          endif
       case ('pres') !Air pressure
          if (Lpoint) then
             icon%p(1:npoints,:) = x2(1:npoints,1:nlayers)
          else
             call map_ll_to_point(dim1,dim2,npoints,x3=x3,y2=icon%p)
          endif
       case ('ta') !Air temperature
          if (Lpoint) then
             icon%T(1:npoints,:) = x2(1:npoints,1:nlayers)
          else
             call map_ll_to_point(dim1,dim2,npoints,x3=x3,y2=icon%T)
          endif
       case ('hus') !Specific humidity
          if (Lpoint) then
             icon%q(1:npoints,:) = x2(1:npoints,1:nlayers)
          else
             call map_ll_to_point(dim1,dim2,npoints,x3=x3,y2=icon%q)
          endif
       case ('clc') !Cloud cover
          if (Lpoint) then
             icon%clc(1:npoints,1:nlayers) = x2(1:npoints,1:nlayers)
          else
             call map_ll_to_point(dim1,dim2,npoints,x3=x3,y2=icon%clc)
          endif
       case ('clw') !Specific cloud water content
          if (Lpoint) then
             icon%clw(1:npoints,:) = x2(1:npoints,1:nlayers)
          else
             call map_ll_to_point(dim1,dim2,npoints,x3=x3,y2=icon%clw)
          endif
       case ('cli') !Specific cloud ice content
          if (Lpoint) then
             icon%cli(1:npoints,:) = x2(1:npoints,1:nlayers)
          else
             call map_ll_to_point(dim1,dim2,npoints,x3=x3,y2=icon%cli)
          endif
       case ('qnc') !Cloud droplet number concentration
          if (Lpoint) then
             icon%qnc(1:npoints,:) = x2(1:npoints,1:nlayers)
          else
             call map_ll_to_point(dim1,dim2,npoints,x3=x3,y2=icon%qnc)
          endif
       case ('qr') !Rain mixing ratio
          if (Lpoint) then
             icon%qr(1:npoints,:) = x2(1:npoints,1:nlayers)
          else
             call map_ll_to_point(dim1,dim2,npoints,x3=x3,y2=icon%qr)
          endif
       case ('qs') !Snow mixing ratio
          if (Lpoint) then
             icon%qs(1:npoints,:) = x2(1:npoints,1:nlayers)
          else
             call map_ll_to_point(dim1,dim2,npoints,x3=x3,y2=icon%qs)
          endif

       end select

       !! Free memory
       if(vrank == 1) deallocate(x1)
       if(vrank == 2) deallocate(x2)
       if(vrank == 3) deallocate(x3)
       if(vrank == 4) deallocate(x4)
       if(vrank == 5) deallocate(x5)

    enddo

    !!=======================================================================================================================!!

    !!=======================================================================================================================!!
    !! Fill in the lat/lon vectors with the right values for 2D modes                                                        !!
    !!=======================================================================================================================!!

    !This might be helpful if the inputs are 2D (gridded) and you want outputs in 1D mode
    allocate(plon(npoints),plat(npoints))
    allocate(ll(npoints))

    if(icon%mode == 2) then !(lon,lat)
       ll = lat
       do j = 1, dim2
          do i = 1, dim1
             k = (j-1)*dim1 + i
             plon(k) = i
             plat(k) = j
          enddo
       enddo
       icon%lon(1:npoints) = lon(plon(1:npoints))
       icon%lat(1:npoints) = ll(plat(1:npoints))
    else if (icon%mode == 3) then !(lat,lon)
       ll = lon
       do j = 1, dim2
          do i = 1, dim1
             k = (j-1)*dim1 + i
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
    if (errst /= 0) then
       errmsg = "Error in nf90_close"
       call s3com_error(routine_name,errmsg,errcode=errst)
    endif


  end subroutine icon_read

end module mod_io_icon

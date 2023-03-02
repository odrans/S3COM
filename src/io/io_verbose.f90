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

module mod_io_verbose

  use s3com_types, only: type_nml, type_model, type_s3com, type_rttov_opt
  use rttov_types, only: rttov_options

  implicit none

  character(len=*), parameter :: hyphens = "-------------------------------------------------------------------"
  character(len=*), parameter :: s3com_version = "v0.9.0-beta"
  character(len=*), parameter :: rttov_version = "v13.1"

  private
  public :: verbose_namelist, verbose_model, verbose_rttov

contains

  subroutine verbose_namelist(fname_nml, nml)

    type(type_nml), intent(in) :: nml

    character(len = 256), intent(in) :: fname_nml

    character(len = 256) :: fn_out_rad

    fn_out_rad = trim(nml%path_out)//"S3COM"//trim(nml%suffix_out)//"_rad.nc"

    write(*,*)
    write(*,'(A)') hyphens
    write(*,'(1X, 3A)') "General setup: S3COM (", s3com_version, ")"
    write(*,'(2X, A, 1X, A)') "- Namelist file:", trim(fname_nml)
    write(*,'(2X, A)') "- S3COM output: "
    write(*,'(4X, A, 1X, A, 1X, 3A)') "| radiances:", write_bool(.TRUE.), "(", trim(fn_out_rad), ")"
    write(*,'(4X, A, 1X, A)') "| retrievals: ", write_bool(nml%flag_retrievals)
    write(*,'(4X, A, 1X, A)') "| jacobians: ", write_bool(nml%flag_output_jac)
    write(*,'(4X, A, 1X, A)') "| atmospheric data: ", write_bool(nml%flag_output_atm)
    write(*,*) hyphens

  end subroutine verbose_namelist

  subroutine verbose_model(nml, model)

    type(type_nml), intent(in)   :: nml
    type(type_model), intent(in)   :: model

    character(len=7), dimension(3), parameter :: mode_desc = (/"track  ", "lon-lat", "lat-lon"/)

    write(*,*)
    write(*,'(A)') hyphens
    write(*,"(1X, A, 1X, A)") "Physical model: ", trim(nml%model_name)
    write(*,"(2X, A, 1X, A)") "- Input file:", trim(nml%fname_in)
    write(*,'(2X, A)') "- Grid description:"
    write(*,"(4X, A, 1X, A)") "| type:", trim(mode_desc(model%mode))
    if(model%mode == 1) write(*,"(4X, A, 1X, I5)") "| npoints:", model%npoints
    if(model%mode > 1) then
      write(*,"(4X, A, 1X, I5, &
         2(1X, A, 1X, I3), A)") "| npoints:", model%npoints, &
         "(nlat:", model%nlat,", nlon:", model%nlon, ")"
    end if



    write(*,"(4X, A, 1X, I3)") "| nlayers:", model%nlayers
    write(*,'(A)') hyphens

  end subroutine verbose_model


  subroutine verbose_rttov(s3com, rttov_opt, opts)

    type(type_rttov_opt), intent(in) :: rttov_opt
    type(type_s3com), intent(in) :: s3com
    type(rttov_options), intent(in) :: opts

    character(len = 12), dimension(2), parameter :: ir_scatt_model = (/"DOM         ", "Chou-scaling"/)
    character(len=17), dimension(3), parameter :: vis_scatt_model = (/"DOM              ", &
    "Single-scattering", "MFASIS           "/)
    character(len = 19), dimension(2), parameter :: gas_units = (/"kg/kg over moist air ", "ppmv over moist air  "/)
    character(len = 20), dimension(2), parameter :: cld_units = (/"kg/kg over moist air", "g m-3               "/)
    character(len = 20), dimension(2), parameter :: aer_units = (/"kg/kg over moist air", "cm-3                "/)
    character(len = 3), dimension(6), parameter :: gas_rttov = (/"O3 ", "CO2", "N2O", "CH4", "CO ", "SO2"/)
    character(len = 10), dimension(3), parameter :: ice_scheme = (/"Baum      ", "Baran 2014", "Baran 2018"/)
    character(len = 4), dimension(2), parameter :: clw_scheme = (/"OPAC", "Deff"/)

    logical, dimension(6) :: gas_rttov_used

    integer(4) :: nstream, i, mmr_cldaer
    character(len=32) :: ir_model, vis_model, nstreams_char, testchar, testchar2
    logical :: flag_gases

    flag_gases = any([opts%dev%do_opdep_calc, opts%rt_all%addrefrac])

    write(nstreams_char, '(I3)') rttov_opt%dom_nstreams

    mmr_cldaer = 1
    if(.NOT. rttov_opt%mmr_cldaer) mmr_cldaer = 2

    gas_rttov_used = [rttov_opt%ozone_data, rttov_opt%co2_data, rttov_opt%n2o_data, &
         rttov_opt%ch4_data, rttov_opt%co_data, rttov_opt%so2_data]

    ir_model = ir_scatt_model(rttov_opt%ir_scatt_model)
    if(rttov_opt%ir_scatt_model == 1) ir_model = trim(ir_model)//" (nstream: "//trim(adjustl(nstreams_char))//")"

    vis_model = vis_scatt_model(rttov_opt%vis_scatt_model)
    if(rttov_opt%vis_scatt_model == 1) vis_model = trim(vis_model)//" (nstream: "//trim(adjustl(nstreams_char))//")"

    write(*,*)
    write(*,'(A)') hyphens
    write(*,'(1X, A)') "Radiative model: RTTOV ("//rttov_version//")"
    write(*,'(2X, A, 1X, A)') "- Path:", trim(s3com%nml%path_rttov)
    write(*,'(2X, A, 1X, I3)') "- # parallel threads:", rttov_opt%nthreads
    write(*,'(2X, A)') "- Instrument:"
    write(*,'(4X, A, 1X, A)') "| name/platform:", trim(s3com%opt%rttov%inst_name)//"/"//trim(s3com%opt%rttov%platform_name)
    write(*,'(4X, A, 1X, I4)') "| nchannels:", rttov_opt%nchannels
    write(*,'(4X, A, 2(1X, F6.3, 1X, A))') "| wavelength range:", &
    minval(s3com%rad%wavelength), "-", maxval(s3com%rad%wavelength), "um"
    write(*,'(2X, A)') "- Radiative transfer:"
    write(*,'(4X, A, 1X, A)') "| IR model:", ir_model
    write(*,'(4X, A, 1X, A)') "| VIS model:", vis_model
    write(*,'(2X, A)') "- Atmosphere:"
    write(*,'(4X, 4(A))') "| gases: ", write_bool(flag_gases)
    if(flag_gases) then
       write(*,'(6X, 4(A))') "| absorption: ", write_bool(opts%dev%do_opdep_calc), &
                             "; refraction: ", write_bool(opts%rt_all%addrefrac)
       write(*,'(6X, A, $)') "| from RTTOV: "
       do i = 1, size(gas_rttov)
          if (.NOT.gas_rttov_used(i)) then
             write(*,'(2A,$)') trim(gas_rttov(i)), " "
          end if
       end do
       write(*,'(1X, A, $)') "; from user: H2O "
       do i = 1, size(gas_rttov)
          if (gas_rttov_used(i)) then
             write(*,'(2A,$)') trim(gas_rttov(i)), " "
          end if
       end do
       write(*,'(A)')
       write(*,'(6X, 4(A))') "| input units: ", trim(gas_units(rttov_opt%gas_units))
    end if
    write(*,'(4X, A, 1X, A, 1X, A, 1X, I3)') "| clouds:", write_bool(opts%rt_ir%addclouds)
    if(opts%rt_ir%addclouds) then
       write(*,'(6X, 5(A))') "| schemes: ", "ice: ", trim(ice_scheme(rttov_opt%ice_scheme)), &
                             " ; liquid: ", trim(clw_scheme(rttov_opt%clw_scheme))
       write(*,'(6X, 2(A))') "| input units: ", trim(cld_units(mmr_cldaer))
    end if
    write(*,'(4X, A, 1X, A)') "| aerosols:", write_bool(opts%rt_ir%addaerosl)
    if(opts%rt_ir%addaerosl) then
       write(*,'(6X, 4(A))') "| input units: ", trim(aer_units(mmr_cldaer))
    end if
    write(*,'(A)') hyphens

  end subroutine verbose_rttov


  character(len=5) function write_bool(flag)
    logical, intent(in) :: flag
    if (flag) then
       write_bool = "TRUE"
    else
       write_bool = "FALSE"
    end if
  end function write_bool


end module mod_io_verbose

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

module mod_cld_phys_sb

  use s3com_types, only: wp

  implicit none

  private
  public :: re_sb

contains

  !> Compute the effective radius of a cloud particle consistent with the Seifert-Beheng microphysics scheme for
  !> a given hydrometeor class
  !> Inputs:
  !> - L: liquid water content [kg m-3]
  !> - N: number concentration [m-3]
  !> - hydro_class: 1: droplet; 2: ice crystal; 3: rain; 4: snow crystal; 5: graupel; 6: hail
  !> Returns:
  !> - re: effective radius [m]
  function re_sb(L, N, hydro_class) result(re)

    integer(kind = 4) :: hydro_class

    ! Cloud droplet parameters; hydro_class = 1
    real(wp), parameter :: nu_drop = 1.0_wp
    real(wp), parameter :: mu_drop = 1.0_wp
    real(wp), parameter :: a_drop = 1.24E-1
    real(wp), parameter :: b_drop = 1.0_wp / 3.0_wp

    ! Cloud ice crystal parameters; hydro_class = 2
    real(wp), parameter :: nu_ice = 0.0_wp
    real(wp), parameter :: mu_ice = 1.0_wp / 3.0_wp
    real(wp), parameter :: a_ice = 0.835_wp
    real(wp), parameter :: b_ice = 0.390_wp

    ! Rain parameters; hydro_class = 3
    real(wp), parameter :: nu_rain = 0.0_wp
    real(wp), parameter :: mu_rain = 1.0_wp / 3.0_wp
    real(wp), parameter :: a_rain = 1.24E-1
    real(wp), parameter :: b_rain = 0.333_wp

    ! Snow parameters; hydro_class = 4
    real(wp), parameter :: nu_snow = 0.0_wp
    real(wp), parameter :: mu_snow = 0.5_wp
    real(wp), parameter :: a_snow = 5.13_wp
    real(wp), parameter :: b_snow = 0.5_wp

    ! Graupel parameters; hydro_class = 5
    real(wp), parameter :: nu_graupel = 1.0_wp
    real(wp), parameter :: mu_graupel = 1.0_wp / 3.0_wp
    real(wp), parameter :: a_graupel = 0.142_wp
    real(wp), parameter :: b_graupel = 0.314_wp

    ! Hail parameters; hydro_class = 6
    real(wp), parameter :: nu_hail = 1.0_wp
    real(wp), parameter :: mu_hail = 1.0_wp / 3.0_wp
    real(wp), parameter :: a_hail = 0.1366_wp
    real(wp), parameter :: b_hail = 0.333_wp

    real(wp), intent(in) :: L
    real(wp), intent(in) :: N

    real(wp) :: re

    select case(hydro_class)
    case(1)
       re = re_sb_all(L, N, a_drop, b_drop, nu_drop, mu_drop)
    case(2)
       re = re_sb_all(L, N, a_ice, b_ice, nu_ice, mu_ice)
    case(3)
       re = re_sb_all(L, N, a_rain, b_rain, nu_rain, mu_rain)
    case(4)
       re = re_sb_all(L, N, a_snow, b_snow, nu_snow, mu_snow)
    case(5)
       re = re_sb_all(L, N, a_graupel, b_graupel, nu_graupel, mu_graupel)
    case(6)
       re = re_sb_all(L, N, a_hail, b_hail, nu_hail, mu_hail)
    end select

  end function re_sb

  !> Compute the effective radius of a cloud particle consistent with the Seifert-Beheng microphysics scheme
  !> Valid equation for all hydrometeor classes. Based on Seifert and Beheng (2005) doi:10.1007/s00703-005-0112-4
  function re_sb_all(L, N, a, b, nu, mu) result(re)

    real(wp) :: re

    real(wp), intent(in) :: a
    real(wp), intent(in) :: b
    real(wp), intent(in) :: nu
    real(wp), intent(in) :: mu
    real(wp), intent(in) :: L
    real(wp), intent(in) :: N

    re = a / 2.0_wp * gamma( (3._wp * b + nu + 1) / mu ) / gamma((2._wp * b + nu + 1) / mu ) * (L / N * gamma((nu + 1) / mu) / gamma((nu + 2) / mu) )**b

  end function re_sb_all


end module mod_cld_phys_sb

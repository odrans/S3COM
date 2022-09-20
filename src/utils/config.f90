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

MODULE S3COM_CONFIG

   USE s3com_types, ONLY: wp
   
   IMPLICIT NONE
   
   REAL(KIND=wp), PARAMETER ::  &
      apriori_iwp       = 1E-1, & ! in kg.m-2
      apriori_iwp_error = 10.0    ! in %
   
   !!Parallel options
   INTEGER(kind=4) :: &
      RTTOV_NTHREADS = 36 !35
   
   !!S3COM options
   REAL(KIND=wp), PARAMETER ::  &
      iwp_lay_threshold = 1E-4, &
      lwp_lay_threshold = 1E-4
   
   INTEGER(KIND=4), PARAMETER :: &
      nstates = 1
   
   !!RTTOV optics
   INTEGER(kind=4), PARAMETER :: &
      RTTOV_DOSOLAR    = 1

   !!Mixing ratios of trace gases
   REAL(wp) ::               &
      mr_co2 = 5.241e-04_wp, &
      mr_ch4 = 9.139e-07_wp, &
      mr_n2o = 4.665e-07_wp, &
      mr_co  = 2.098e-07_wp
   
   REAL(wp), PARAMETER :: &
      tmelt  = 273.15_wp, & !Melting temperature of ice/snow (K)
      rhoice = 917._wp,   & !Density of ice (kg/m3)
      rholiq = 1000._wp     !Density of liquid water (kg/m3)
   
   !!Molecular weights (g/mol)
   REAL(wp), PARAMETER ::  &
      amw   = 18.01534_wp, & !Water
      amd   = 28.9644_wp,  & !Dry air
      amO3  = 47.9983_wp,  & !Ozone
      amCO2 = 44.0096_wp,  & !Cabone dioxide
      amCH4 = 16.0426_wp,  & !Methane
      amN2O = 44.0129_wp,  & !Nitrogen dioxide
      amCO  = 28.0102_wp     !Carbone monoxide
   
   !!WMO/SI value
   REAL(wp), PARAMETER :: &
      avo  = 6.023E23_wp, & !Avogadro constant used by ISCCP simulator (1/mol)
      grav = 9.806650_wp    !Av. gravitational acceleration used by ISCCP simulator (m/s2)
   
   !!Thermodynamic constants for the dry and moist atmosphere
   REAL(wp), PARAMETER ::    &
      rd      = 287.058_wp,  & !Gas constant for dry air (J/K/Kg)
      cpd     = 1004.64_wp,  & !Specific heat at constant pressure for dry air (J/K/Kg)
      rv      = 461.51_wp,   & !Gas constant for water vapor (J/K/Kg)
      cpv     = 1869.46_wp,  & !Specific heat at constant pressure for water vapor (J/K/Kg)
      km      = 1.38e-23_wp, & !Boltzmann constant (J/K)
      epsilon = 0.622_wp
   
   !!Constants for the droplets size distribution
   REAL(wp), PARAMETER ::    &
      a  = 0.124_wp,         &
      b  = 1._wp/3._wp,      &
      mu = 1._wp,            &
      nu = 1._wp
   
   !!Optics
   REAL(wp), PARAMETER :: &
      Q_ext = 2._wp
   
END MODULE S3COM_CONFIG

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

module s3com_config
   
   use s3com_types, only: wp, wpi
   
   implicit none
   
   public
   
   ! ============================================================================================================================
   ! Retrieval options
   ! -----------------
   real(wp), parameter ::               &
      lwp_apriori        = 1.08E-01_wp, &      !< Liquid water path a priori value                             @units{kg m-2}
      cdnc_apriori       = 1.22E+08_wp, &      !< Cloud droplet number concentration a priori value            @units{\# m-3}
      !lwp_apriori        = 6.82E-02_wp, &      !< Liquid water path a priori value                             @units{kg m-2}
      !cdnc_apriori       = 7.48E+07_wp, &      !< Cloud droplet number concentration a priori value            @units{\# m-3}
      lwp_apriori_error  = 5.0_wp,      &      !< Liquid water path a priori error (500 %)                     @units{-}
      cdnc_apriori_error = 5.0_wp,      &      !< Cloud droplet number concentration a priori error (500 %)    @units{-}
      cw                 = 1.9E-06_wp,  &      !< Condensation rate of water. Default is 1.9E-06.              @units{kg m-4}
      k                  = 0.80_wp             !< Droplet size distribution shape parameter. Default is 0.8.   @units{-}
   
   integer(wpi), parameter ::           &
      nstates = 2                              !< Size of the state vector X
   ! ============================================================================================================================
   
   ! ============================================================================================================================
   ! S3COM options
   ! -------------
   real(wp), parameter ::               &
      lwp_threshold     = 1.0E-04_wp,   &      !< Threshold of the cloud liquid water path                     @units{kg m-2}
      iwp_threshold     = 1.0E-04_wp,   &      !< Threshold of the cloud ice water path                        @units{kg m-2}
      lwp_lay_threshold = 1.0E-04_wp,   &      !< Threshold of the cloud layer liquid water path               @units{kg m-2}
      iwp_lay_threshold = 1.0E-04_wp           !< Threshold of the cloud layer ice water path                  @units{kg m-2}
   ! ============================================================================================================================
   
   ! ============================================================================================================================
   ! RTTOV options
   ! -------------
   integer(wpi), parameter ::           &
      rttov_dosolar = 1
   ! ============================================================================================================================
   
   ! ============================================================================================================================
   ! Thermodynamic constants for water phases
   ! ----------------------------------------
   real(wp), parameter ::               &
      tmelt  = 273.15_wp,               &      !< Melting temperature of ice/snow                              @units{K}
      rhoice = 917.0_wp,                &      !< Density of ice                                               @units{kg m-3}
      rholiq = 1000.0_wp                       !< Density of liquid water                                      @units{kg m-3}
   ! ============================================================================================================================
   
   ! ============================================================================================================================
   ! Mixing ratios of trace gases
   ! ----------------------------
   real(wp) ::                          &
      mr_co2 = 5.241E-04_wp,            &      !< Carbon dioxide
      mr_ch4 = 9.139E-07_wp,            &      !< Methane
      mr_n2o = 4.665E-07_wp,            &      !< Nitrogen dioxide
      mr_co  = 2.098E-07_wp                    !< Carbon monoxide
   ! ============================================================================================================================
   
   ! ============================================================================================================================
   ! Molecular weights
   ! -----------------
   real(wp), parameter ::               &
      amw   = 18.01534_wp,              &      !< Water                                                        @units{g mol-1}
      amd   = 28.9644_wp,               &      !< Dry air                                                      @units{g mol-1}
      amO3  = 47.9983_wp,               &      !< Ozone                                                        @units{g mol-1}
      amCO2 = 44.0096_wp,               &      !< Cabone dioxide                                               @units{g mol-1}
      amCH4 = 16.0426_wp,               &      !< Methane                                                      @units{g mol-1}
      amN2O = 44.0129_wp,               &      !< Nitrogen dioxide                                             @units{g mol-1}
      amCO  = 28.0102_wp                       !< Carbone monoxide                                             @units{g mol-1}
   ! ============================================================================================================================
   
   ! ============================================================================================================================
   ! WMO/SI value
   ! ------------
   real(wp), parameter ::               &
      avo  = 6.023E+23_wp,              &      !< Avogadro constant used by ISCCP simulator                    @units{mol-1}
      grav = 9.806650_wp                       !< Av. gravitational acceleration used by ISCCP simulator       @units{m s-2}
   ! ============================================================================================================================
   
   ! ============================================================================================================================
   ! Thermodynamic constants for the dry and moist atmosphere
   ! --------------------------------------------------------
   real(wp), parameter ::               &
      rd  = 287.058_wp,                 &      !< Gas constant for dry air                                     @units{J K-1 kg-1}
      cpd = 1004.64_wp,                 &      !< Specific heat at constant pressure for dry air               @units{J K-1 kg-1}
      rv  = 461.51_wp,                  &      !< Gas constant for water vapor                                 @units{J K-1 kg-1}
      cpv = 1869.46_wp,                 &      !< Specific heat at constant pressure for water vapor           @units{J K-1 kg-1}
      km  = 1.38E-23_wp                        !< Boltzmann constant                                           @units{J K-1}
   ! ============================================================================================================================
   
   ! ============================================================================================================================
   ! Optics
   ! ------
   real(wp), parameter ::               &
      Q_ext = 2.0_wp,                   &      !< Extinction efficiency factor                                 @units{-}
      pi = 4.0_wp * atan(1.0_wp)               !< Numerical definition of PI                                   @units{-}
   ! ============================================================================================================================
   
end module s3com_config

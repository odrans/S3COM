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

MODULE MOD_RTTOV_SETUP
   
   USE s3com_types, ONLY: wp
   
   IMPLICIT NONE
   
   CONTAINS
   
      SUBROUTINE RTTOV_SETUP_OPT(zenangle,azangle,sunzenangle,sunazangle,rttov_opt, nml)
         
         USE s3com_types,  ONLY: type_rttov_opt, type_nml
         USE s3com_config, ONLY:                                 &
             RTTOV_DOSOLAR
      
            !!Input variables
            REAL(wp), INTENT(IN) :: zenangle, azangle, sunzenangle, sunazangle
            
            !!Output variables
            TYPE(type_rttov_opt), INTENT(OUT) :: rttov_opt
            TYPE(type_nml), INTENT(IN) :: nml

            rttov_opt%platform   = nml%platform
            rttov_opt%satellite  = nml%satellite
            rttov_opt%instrument = nml%instrument
            rttov_opt%dosolar    = RTTOV_DOSOLAR
            rttov_opt%nchannels  = nml%nchannels
            
            ALLOCATE(rttov_opt%channel_list(rttov_opt%nchannels))
            
            rttov_opt%channel_list = nml%channel_list
            rttov_opt%month        = nml%month
            rttov_opt%zenangle     = zenangle
            rttov_opt%azangle      = azangle
            rttov_opt%sunzenangle  = sunzenangle
            rttov_opt%sunazangle   = sunazangle
            
      END SUBROUTINE RTTOV_SETUP_OPT
      
      SUBROUTINE rttov_setup_atm(idx_start,idx_end,icon,rttov_atm)
         
         USE s3com_types, ONLY: type_rttov_atm, type_icon
         
         !!Input variables
         INTEGER, TARGET, INTENT(IN)         :: idx_start, idx_end
         TYPE(type_icon), TARGET, INTENT(IN) :: icon
         
         !!Output variables
         TYPE(type_rttov_atm)  :: rttov_atm
         INTEGER, SAVE, TARGET :: nidx
         
         nidx = idx_end - idx_start + 1
         
         rttov_atm%nPoints   => nidx
         rttov_atm%nLevels   => icon%nLevels
         rttov_atm%idx_start => idx_start
         rttov_atm%idx_end   => idx_end
         rttov_atm%co2       => icon%co2
         rttov_atm%ch4       => icon%ch4
         rttov_atm%n2o       => icon%n2o
         rttov_atm%co        => icon%co
         rttov_atm%h_surf    => icon%orography(idx_start:idx_end)
         rttov_atm%u_surf    => icon%u_wind(idx_start:idx_end)
         rttov_atm%v_surf    => icon%v_wind(idx_start:idx_end)
         rttov_atm%t_skin    => icon%skt(idx_start:idx_end)
         rttov_atm%p_surf    => icon%psfc(idx_start:idx_end)
         rttov_atm%q2m       => icon%q2m(idx_start:idx_end)
         rttov_atm%t2m       => icon%t2m(idx_start:idx_end)
         rttov_atm%lsmask    => icon%landmask(idx_start:idx_end)
         rttov_atm%lat       => icon%lat(idx_start:idx_end)
         rttov_atm%lon       => icon%lon(idx_start:idx_end)
         rttov_atm%p         => icon%p(idx_start:idx_end,:)
         rttov_atm%z         => icon%z(idx_start:idx_end,:)
         rttov_atm%dz        => icon%dz(idx_start:idx_end,:)
         rttov_atm%t         => icon%t(idx_start:idx_end,:)
         rttov_atm%q         => icon%sh(idx_start:idx_end,:)
         rttov_atm%tca       => icon%tca(idx_start:idx_end,:)
         rttov_atm%iwc       => icon%iwc(idx_start:idx_end,:)
         rttov_atm%lwc       => icon%lwc(idx_start:idx_end,:)
         rttov_atm%reff      => icon%reff(idx_start:idx_end,:)
         rttov_atm%cdnc      => icon%cdnc(idx_start:idx_end,:)
         
      END SUBROUTINE rttov_setup_atm

END MODULE MOD_RTTOV_SETUP

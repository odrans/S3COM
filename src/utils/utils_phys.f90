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

module mod_utils_phys
   
   use s3com_types,    only: wp, wpi
   use mod_utils_math, only: day_number
   
   implicit none
   
   private
   public :: solar_angles, sunae
   
contains
   
   ! ============================================================================================================================
   !> @brief Calculate the solar viewing angles at a specific location and time
   !! @param[in] lat            latitude @units{degrees North}
   !! @param[in] lon            longitude @units{degrees East}
   !! @param[in] date           date of the simulation
   !! @param[in] time           time of the simulation
   !! @param[out] sunzenangle   solar zenith angle @units{°}
   !! @param[out] sunazangle    solar azimuth angle @units{°}
   subroutine solar_angles(lat, lon, date, time, sunzenangle, sunazangle)
      
      ! Input
      real(wp), intent(in) :: lat, lon
      integer(wpi), dimension(3), intent(in) :: date, time
      
      ! Output
      real(wp), intent(out) :: sunzenangle, sunazangle
      
      ! Internal
      real(wp) :: elevation, dec, soldst, hour
      integer(wpi) :: julian
      
      ! Compute the julian date (day of year)
      call day_number(date(1), date(2), date(3), julian)
      
      ! Compute the fractional hour
      hour = time(1) + time(2) / 60._wp + time(3) / 3600._wp
      
      ! Get the angles
      call sunae(date(3), julian, hour, lat, lon, sunazangle, elevation, dec, soldst)
      
      ! Convert elevation into solar zenith angle
      sunzenangle = 90 - elevation
      
      return
      
   end subroutine solar_angles
   ! ============================================================================================================================
   
   ! ============================================================================================================================
   !> @brief Calculate the local solar azimuth and elevation angles, and the distance to and angle subtended by the Sun, at a 
   !! specific location and time using approximate formulas in The Astronomical Almanac. See details below.
   subroutine sunae(YEAR, DAY, HOUR, LAT, LONG, AZ, EL, SOLDIA, SOLDST)
      
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! RCS version control information:
      ! $Header: sunae.f,v 1.3 96/05/30 09:30:15 wiscombe Exp $
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      !     Calculates the local solar azimuth and elevation angles, and
      !     the distance to and angle subtended by the Sun, at a specific
      !     location and time using approximate formulas in The Astronomical
      !     Almanac.  Accuracy of angles is 0.01 deg or better (the angular
      !     width of the Sun is about 0.5 deg, so 0.01 deg is more than
      !     sufficient for most applications).
      
      !     Unlike many GCM (and other) sun angle routines, this
      !     one gives slightly different sun angles depending on
      !     the year.  The difference is usually down in the 4th
      !     significant digit but can slowly creep up to the 3rd
      !     significant digit after several decades to a century.
      
      !     A refraction correction appropriate for the "US Standard
      !     Atmosphere" is added, so that the returned sun position is
      !     the APPARENT one.  The correction is below 0.1 deg for solar
      !     elevations above 9 deg.  To remove it, comment out the code
      !     block where variable REFRAC is set, and replace it with
      !     REFRAC = 0.0.
      
      !   References:
      
      !     Michalsky, J., 1988: The Astronomical Almanac's algorithm for
      !        approximate solar position (1950-2050), Solar Energy 40,
      !        227-235 (but the version of this program in the Appendix
      !        contains errors and should not be used)
      
      !     The Astronomical Almanac, U.S. Gov't Printing Office, Washington,
      !        D.C. (published every year): the formulas used from the 1995
      !        version are as follows:
      !        p. A12: approximation to sunrise/set times
      !        p. B61: solar elevation ("altitude") and azimuth
      !        p. B62: refraction correction
      !        p. C24: mean longitude, mean anomaly, eclipti! longitude,
      !                obliquity of ecliptic, right ascension, declination,
      !                Earth-Sun distance, angular diameter of Sun
      !        p. L2:  Greenwich mean sidereal time (ignoring T^2, T^3 terms)
      
      !   Authors:  Dr. Joe Michalsky (joe@asrc.albany.edu)
      !             Dr. Lee Harrison (lee@asrc.albany.edu)
      !             Atmospheric Sciences Research Center
      !             State University of New York
      !             Albany, New York
      
      !   Modified by:  Dr. Warren Wiscombe (wiscombe@climate.gsfc.nasa.gov)
      !                 NASA Goddard Space Flight Center
      !                 Code 913
      !                 Greenbelt, MD 20771
      
      !   WARNING:  Do not use this routine outside the year range
      !             1950-2050.  The approximations become increasingly
      !             worse, and the calculation of Julian date becomes
      !             more involved.
      
      !   Input:
      
      !      YEAR     year (INTEGER; range 1950 to 2050)
      
      !      DAY      day of year at LAT-LONG location (INTEGER; range 1-366)
      
      !      HOUR     hour of DAY [GMT or UT] (REAL; range -13.0 to 36.0)
      !               = (local hour) + (time zone number)
      !                 + (Daylight Savings Time correction; -1 or 0)
      !               where (local hour) range is 0 to 24,
      !               (time zone number) range is -12 to +12, and
      !               (Daylight Time correction) is -1 if on Daylight Time
      !               (summer half of year), 0 otherwise;
      !               Example: 8:30 am Eastern Daylight Time would be
      
      !                           HOUR = 8.5 + 5 - 1 = 12.5
      
      !      LAT      latitude [degrees]
      !               (REAL; range -90.0 to 90.0; north is positive)
      
      !      LONG     longitude [degrees]
      !               (REAL; range -180.0 to 180.0; east is positive)
      
      !   Output:
      
      !      AZ       solar azimuth angle (measured east from north,
      !               0 to 360 degs)
      
      !      EL       solar elevation angle [-90 to 90 degs];
      !               solar zenith angle = 90 - EL
      
      !      SOLDIA   solar diameter [degs]
      
      !      SOLDST   distance to sun [Astronomical Units, AU]
      !               (1 AU = mean Earth-sun distance = 1.49597871E+11 m
      !                in IAU 1976 System of Astronomical Constants)
      
      !   Local Variables:
      
      !     DEC       Declination (radians)
      
      !     ECLONG    Ecliptic longitude (radians)
      
      !     GMST      Greenwich mean sidereal time (hours)
      
      !     HA        Hour angle (radians, -pi to pi)
      
      !     JD        Modified Julian date (number of days, including
      !               fractions thereof, from Julian year J2000.0);
      !               actual Julian date is JD + 2451545.0
      
      !     LMST      Local mean sidereal time (radians)
      
      !     MNANOM    Mean anomaly (radians, normalized to 0 to 2*pi)
      
      !     MNLONG    Mean longitude of Sun, corrected for aberration
      !               (deg; normalized to 0-360)
      
      !     OBLQEC    Obliquity of the ecliptic (radians)
      
      !     RA        Right ascension  (radians)
      
      !     REFRAC    Refraction correction for US Standard Atmosphere (degs)
      
      ! --------------------------------------------------------------------
      !   Uses double precision for safety and because Julian dates can
      !   have a large number of digits in their full form (but in practice
      !   this version seems to work fine in single precision if you only
      !   need about 3 significant digits and aren't doing precise climate
      !   change or solar tracking work).
      ! --------------------------------------------------------------------
      
      !   Why does this routine require time input as Greenwich Mean Time
      !   (GMT; also called Universal Time, UT) rather than "local time"?
      !   Because "local time" (e.g. Eastern Standard Time) can be off by
      !   up to half an hour from the actual local time (called "local mean
      !   solar time").  For society's convenience, "local time" is held
      !   constant across each of 24 time zones (each technically 15 longitude
      !   degrees wide although some are distorted, again for societal
      !   convenience).  Local mean solar time varies continuously around a
      !   longitude belt;  it is not a step function with 24 steps.
      !   Thus it is far simpler to calculate local mean solar time from GMT,
      !   by adding 4 min for each degree of longitude the location is
      !   east of the Greenwich meridian or subtracting 4 min for each degree
      !   west of it.
      
      ! --------------------------------------------------------------------
      
      !   TIME
      !
      !   The measurement of time has become a complicated topic.  A few
      !   basic facts are:
      !
      !   (1) The Gregorian calendar was introduced in 1582 to replace
      !   Julian calendar; in it, every year divisible by four is a leap
      !   year just as in the Julian calendar except for centurial years
      !   which must be exactly divisible by 400 to be leap years.  Thus
      !   2000 is a leap year, but not 1900 or 2100.
      
      !   (2) The Julian day begins at Greenwich noon whereas the calendar
      !   day begins at the preceding midnight;  and Julian years begin on
      !   "Jan 0" which is really Greenwich noon on Dec 31.  True Julian
      !   dates are a continous count of day numbers beginning with JD 0 on
      !   1 Jan 4713 B.C.  The term "Julian date" is widely misused and few
      !   people understand it; it is best avoided unless you want to study
      !   the Astronomical Almanac and learn to use it correctly.
      
      !   (3) Universal Time (UT), the basis of civil timekeeping, is
      !   defined by a formula relating UT to GMST (Greenwich mean sidereal
      !   time).  UTC (Coordinated Universal Time) is the time scale
      !   distributed by most broadcast time services.  UT, UTC, and other
      !   related time measures are within a few sec of each other and are
      !   frequently used interchangeably.
      
      !   (4) Beginning in 1984, the "standard epoch" of the astronomical
      !   coordinate system is Jan 1, 2000, 12 hr TDB (Julian date
      !   2,451,545.0, denoted J2000.0).  The fact that this routine uses
      !   1949 as a point of reference is merely for numerical convenience.
      ! --------------------------------------------------------------------
      
      !     .. Scalar Arguments ..
      integer(wpi) YEAR, DAY
      real(wp) AZ, EL, HOUR, LAT, LONG, SOLDIA, SOLDST
      
      !     .. Local Scalars ..
      logical PASS1
      integer(wpi) DELTA, LEAP
      real(wp) DEC, DEN, ECLONG, GMST, HA, JD, LMST, MNANOM, MNLONG, NUM, OBLQEC, PI, RA, RPD, REFRAC, TIME, TWOPI
      
      !     .. Intrinsic Functions ..
      intrinsic AINT, ASIN, ATAN, COS, MOD, SIN, TAN
      
      !     .. Data statements ..
      save PASS1, PI, TWOPI, RPD
      data PASS1 /.true./
      !     ..
      
      if (YEAR .lt. 1950   .or. YEAR .gt. 2050)  stop 'SUNAE--bad input variable YEAR'
      if (DAY  .lt. 1      .or. DAY  .gt. 366)   stop 'SUNAE--bad input variable DAY'
      if (HOUR .lt. -13.0  .or. HOUR .gt. 36.0)  stop 'SUNAE--bad input variable HOUR'
      if (LAT  .lt. -90.0  .or. LAT  .gt. 90.0)  stop 'SUNAE--bad input variable LAT'
      if (LONG .lt. -180.0 .or. LONG .gt. 180.0) stop 'SUNAE--bad input variable LONG'
      
      if (PASS1) then
         PI    = 2. * asin( 1.0 )
         TWOPI = 2. * PI
         RPD   = PI / 180.
         PASS1 = .false.
      endif
      
      !                    ** current Julian date (actually add 2,400,000
      !                    ** for true JD);  LEAP = leap days since 1949;
      !                    ** 32916.5 is midnite 0 jan 1949 minus 2.4e6
      
      DELTA = YEAR - 1949
      LEAP  = DELTA / 4
      JD    = 32916.5 + (DELTA * 365 + LEAP + DAY) + HOUR / 24.
      
      !                    ** last yr of century not leap yr unless divisible
      !                    ** by 400 (not executed for the allowed YEAR range,
      !                    ** but left in so our successors can adapt this for
      !                    ** the following 100 years)
      
      if (mod(YEAR, 100) .eq. 0 .and. mod(YEAR, 400) .ne. 0) JD = JD - 1.
      
      !                     ** ecliptic coordinates
      !                     ** 51545.0 + 2.4e6 = noon 1 jan 2000
      
      TIME  = JD - 51545.0
      
      !                    ** force mean longitude between 0 and 360 degs
      
      MNLONG = 280.460 + 0.9856474 * TIME
      MNLONG = mod(MNLONG, 360.)
      if (MNLONG .lt. 0.) MNLONG = MNLONG + 360.
      
      !                    ** mean anomaly in radians between 0 and 2*pi
      
      MNANOM = 357.528 + 0.9856003 * TIME
      MNANOM = mod(MNANOM, 360.)
      if (MNANOM .lt. 0.) MNANOM = MNANOM + 360.
      
      MNANOM = MNANOM*RPD
      
      !                    ** ecliptic longitude and obliquity
      !                    ** of ecliptic in radians
      
      ECLONG = MNLONG + 1.915 * sin(MNANOM) + 0.020 * sin(2.*MNANOM)
      ECLONG = mod(ECLONG, 360.)
      if (ECLONG .lt. 0.) ECLONG = ECLONG + 360.
      
      OBLQEC = 23.439 - 0.0000004 * TIME
      ECLONG = ECLONG * RPD
      OBLQEC = OBLQEC * RPD
      
      !                    ** right ascension
      
      NUM = cos(OBLQEC) * sin(ECLONG)
      DEN = cos(ECLONG)
      RA  = atan(NUM / DEN)
      
      !                    ** Force right ascension between 0 and 2*pi
      
      if (DEN .lt. 0.0) then
         RA = RA + PI
      else if (NUM .lt. 0.0) then
         RA = RA + TWOPI
      endif
      
      !                    ** declination
      
      DEC = asin(sin(OBLQEC) * sin(ECLONG))
      
      !                    ** Greenwich mean sidereal time in hours
      
      GMST = 6.697375 + 0.0657098242 * TIME + HOUR
      
      !                    ** Hour not changed to sidereal time since
      !                    ** 'time' includes the fractional day
      
      GMST = mod(GMST, 24.)
      if (GMST .lt. 0.) GMST = GMST + 24.
      
      !                    ** local mean sidereal time in radians
      
      LMST = GMST + LONG / 15.
      LMST = mod(LMST, 24.)
      if (LMST .lt. 0.) LMST = LMST + 24.
      
      LMST = LMST * 15. * RPD
      
      !                    ** hour angle in radians between -pi and pi
      
      HA = LMST - RA
      
      if (HA .lt. - PI) HA = HA + TWOPI
      if (HA .gt. PI) HA = HA - TWOPI
      
      !                    ** solar azimuth and elevation
      
      EL = asin(sin(DEC) * sin(LAT*RPD) + cos(DEC) * cos(LAT*RPD) * cos(HA))
      
      AZ = asin(-cos(DEC) * sin(HA) / cos(EL))
      
      !                    ** Put azimuth between 0 and 2*pi radians
      
      if (sin(DEC) - sin(EL) * sin(LAT*RPD) .ge. 0.) then
         if (sin(AZ) .lt. 0.) AZ = AZ + TWOPI
      else
         AZ = PI - AZ
      endif
      
      !                     ** Convert elevation and azimuth to degrees
      EL = EL / RPD
      AZ = AZ / RPD
      
      !  ======== Refraction correction for U.S. Standard Atmos. ==========
      !      (assumes elevation in degs) (3.51823=1013.25 mb/288 K)
      
      if (EL .ge. 19.225) then
         
         REFRAC = 0.00452 * 3.51823 / tan(EL*RPD)
         
      else if (EL .gt. -0.766 .and. EL .lt. 19.225) then
         
         REFRAC = 3.51823 * (0.1594 + EL * (0.0196 + 0.00002 * EL)) / (1. + EL * (0.505 + 0.0845 * EL))
         
      else if (EL .le. -0.766) then
         
         REFRAC = 0.0
         
      endif
      
      EL = EL + REFRAC
      ! ===================================================================
      
      !                   ** distance to sun in A.U. & diameter in degs
      
      SOLDST = 1.00014 - 0.01671 * cos(MNANOM) - 0.00014 * cos(2.*MNANOM)
      SOLDIA = 0.5332 / SOLDST
      
      if (EL .lt. -90.0 .or. EL .gt. 90.0)  stop 'SUNAE--output argument EL out of range'
      if (AZ .lt. 0.0   .or. AZ .gt. 360.0) stop 'SUNAE--output argument AZ out of range'
      
      return
      
   end subroutine sunae
   ! ============================================================================================================================
   
end module mod_utils_phys

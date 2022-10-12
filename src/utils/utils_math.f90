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

MODULE mod_utils_math

use s3com_types, only: wp
  
IMPLICIT NONE

CONTAINS

  subroutine inverse(a,c,n)
    !============================================================
    ! Inverse matrix
    ! Method: Based on Doolittle LU factorization for Ax=b
    ! Alex G. December 2009
    !-----------------------------------------------------------
    ! input ...
    ! a(n,n) - array of coefficients for matrix A
    ! n      - dimension
    ! output ...
    ! c(n,n) - inverse matrix of A
    ! comments ...
    ! the original matrix a(n,n) will be destroyed 
    ! during the calculation
    !===========================================================
    implicit none 
    integer n
    real(wp) :: a(n,n), c(n,n)
    real(wp) :: L(n,n), U(n,n), b(n), d(n), x(n)
    real(wp) :: coeff
    integer i, j, k

    ! step 0: initialization for matrices L and U and b
    ! Fortran 90/95 aloows such operations on matrices
    L=0.0
    U=0.0
    b=0.0

    ! step 1: forward elimination
    do k=1, n-1
       do i=k+1,n
          coeff=a(i,k)/a(k,k)
          L(i,k) = coeff
          do j=k+1,n
             a(i,j) = a(i,j)-coeff*a(k,j)
          end do
       end do
    end do

    ! Step 2: prepare L and U matrices 
    ! L matrix is a matrix of the elimination coefficient
    ! + the diagonal elements are 1.0
    do i=1,n
       L(i,i) = 1.0
    end do
    ! U matrix is the upper triangular part of A
    do j=1,n
       do i=1,j
          U(i,j) = a(i,j)
       end do
    end do

    ! Step 3: compute columns of the inverse matrix C
    do k=1,n
       b(k)=1.0
       d(1) = b(1)
       ! Step 3a: Solve Ld=b using the forward substitution
       do i=2,n
          d(i)=b(i)
          do j=1,i-1
             d(i) = d(i) - L(i,j)*d(j)
          end do
       end do
       ! Step 3b: Solve Ux=d using the back substitution
       x(n)=d(n)/U(n,n)
       do i = n-1,1,-1
          x(i) = d(i)
          do j=n,i+1,-1
             x(i)=x(i)-U(i,j)*x(j)
          end do
          x(i) = x(i)/u(i,i)
       end do
       ! Step 3c: fill the solutions x(n) into column k of C
       do i=1,n
          c(i,k) = x(i)
       end do
       b(k)=0.0
    end do
  end subroutine inverse


  subroutine day_number(day, month, year, julian)

    INTEGER(KIND = 4), INTENT(IN) :: day, month, year

    INTEGER(KIND = 4), INTENT(OUT) :: julian

    if (month.le.2) then
       julian = 31 * (month - 1) + day
       return
    endif
    if (month.gt.8) then
       julian = 31 * (month - 1) - ((month - 2) / 2) - 2 + day
    else
       julian= 31 * (month-1)-((month-1)/2)-2+day
    endif
    if(year.ne.0 .and. mod(year,4).eq.0) julian = julian + 1

    return
  end subroutine day_number


  SUBROUTINE solar_angles(lat, lon, date, time, sunzenangle, sunazangle)

    REAL(KIND = wp), INTENT(IN) :: lat, lon
    INTEGER(KIND = 4), DIMENSION(3), INTENT(IN) :: date, time

    REAL(KIND = wp), INTENT(OUT) :: sunzenangle, sunazangle

    REAL(KIND = wp) :: elevation, dec, soldst, hour
    INTEGER(KIND = 4) :: julian

    ! Compute the julian date (day of year)
    call day_number(date(1), date(2), date(3), julian)

    ! Compute the fractional hour
    hour = time(1) + time(2) / 60._wp + time(3) / 3600._wp

    ! Get the angles
    call sunae(date(3), julian, hour, lat, lon, sunazangle, elevation, dec, soldst)

    ! Convert elevation into solar zenith angle
    sunzenangle = 90 - elevation

    RETURN

  END SUBROUTINE solar_angles



END MODULE mod_utils_math

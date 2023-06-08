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

module mod_utils_math
   
   use s3com_types, only: wp
   
   implicit none
   
   private
   public :: inverse, n_chunks, day_number
   
contains
   
   ! ----------------------------------------------------------------------------------------------------------------------------
   subroutine inverse(a, c, n)
      
      ! ======================================================
      ! Inverse matrix
      ! Method: Based on Doolittle LU factorization for Ax = b
      ! Alex G. December 2009
      ! ------------------------------------------------------
      ! input ...
      ! a(n,n) - array of coefficients for matrix A
      ! n      - dimension
      ! output ...
      ! c(n,n) - inverse matrix of A
      ! comments ...
      ! the original matrix a(n,n) will be destroyed
      ! during the calculation
      ! ======================================================
      
      implicit none
      
      integer n
      real(wp) :: a(n,n), c(n,n)
      real(wp) :: L(n,n), U(n,n), b(n), d(n), x(n)
      real(wp) :: coeff
      integer i, j, k
      
      ! step 0: initialization for matrices L and U and b
      ! Fortran 90/95 aloows such operations on matrices
      L = 0.0_wp
      U = 0.0_wp
      b = 0.0_wp
      
      ! Step 1: forward elimination
      do k = 1, n-1
         do i = k+1, n
            coeff=a(i,k)/a(k,k)
            L(i,k) = coeff
            do j = k+1, n
               a(i,j) = a(i,j) - coeff*a(k,j)
            enddo
         enddo
      enddo
      
      ! Step 2: prepare L and U matrices
      ! L matrix is a matrix of the elimination coefficient
      ! + the diagonal elements are 1.0
      do i = 1, n
         L(i,i) = 1.0_wp
      enddo
      
      ! U matrix is the upper triangular part of A
      do j = 1, n
         do i = 1, j
            U(i,j) = a(i,j)
         enddo
      enddo
      
      ! Step 3: compute columns of the inverse matrix C
      do k = 1, n
         b(k) = 1.0_wp
         d(1) = b(1)
         
         ! Step 3a: Solve Ld=b using the forward substitution
         do i = 2, n
            d(i) = b(i)
            do j = 1, i-1
               d(i) = d(i) - L(i,j)*d(j)
            enddo
         enddo
         
         ! Step 3b: Solve Ux=d using the back substitution
         x(n) = d(n) / U(n,n)
         do i = n-1, 1, -1
            x(i) = d(i)
            do j = n, i+1, -1
               x(i) = x(i) - U(i,j)*x(j)
            enddo
            x(i) = x(i)/u(i,i)
         enddo
         
         ! Step 3c: fill the solutions x(n) into column k of C
         do i = 1, n
            c(i,k) = x(i)
         enddo
         
         b(k) = 0.0_wp
      enddo
      
   end subroutine inverse
   ! ----------------------------------------------------------------------------------------------------------------------------
   
   ! ----------------------------------------------------------------------------------------------------------------------------
   function n_chunks(npoints, npoints_it) result(nchunks)
      
      ! Inputs
      integer(kind=4), intent(in) :: npoints, npoints_it
      
      ! Internal
      integer(kind=4) :: nchunks
      
      nChunks = npoints / npoints_it
      if (mod(npoints,npoints_it)/=0) nchunks = nchunks + 1
      if (npoints .eq. npoints_it) nChunks = 1
      
   end function n_chunks
   ! ----------------------------------------------------------------------------------------------------------------------------
   
   ! ----------------------------------------------------------------------------------------------------------------------------
   subroutine day_number(day, month, year, julian)
      
      ! Inputs
      integer(kind=4), intent(in) :: day, month, year
      
      ! Outputs
      integer(kind = 4), intent(out) :: julian
      
      if (month .le. 2) then
         julian = 31 * (month - 1) + day
         return
      endif
      
      if (month .gt. 8) then
         julian = 31 * (month - 1) - ((month - 2) / 2) - 2 + day
      else
         julian = 31 * (month - 1) - ((month - 1) / 2) - 2 + day
      endif
      
      if(year.ne.0 .and. mod(year,4).eq.0) julian = julian + 1
      
      return
      
   end subroutine day_number
   ! ----------------------------------------------------------------------------------------------------------------------------
   
end module mod_utils_math

!
!   Copyright 1974-2020 George Delic, HiPERISM Consulting LLC
!
!   subroutine SIN
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU General Public License for more details.
!      https://www.gnu.org/licenses/gpl-3.0.en.html
!      https://www.gnu.org/licenses/gpl-3.0.txt
!   A copy is included in this distribution
!
! PURPOSE: compute the Bessel function value by recurrence
!
! AUTHOR:  George Delic, Ph.D., for this source code
!          Reference equation Eq.(31) in
!          G. Delic,
!          The legendre series and a quadrature formula for its coefficients,
!          Journal of Computational Physics, vol. 14 (1974), pp. 254-268.
!
! ORIGIN:  Source code developed by AUTHOR
!
! LANGUAGE:  FOTRAN IV and Fortran 77 (ANSI X3J3) X3.9-1978
!            with minor additions from Fortran 90/95
!
! This code was developed on the Telefunken TR-440 with Fortran IV in 1974
! which provided a word length with some 24 decimals in double precision:
!     IMPLICIT REAL*8(A-H,O-Z)
! However, current commodity processors need to be forced to use quadruple
! precision using IMPLICIT REAL*16(A-H,O-Z) and constants defined as:
!     1.0Q+00, etc.
! In addition the compiler must be forced to avoid optimizations that
! would truncate precision. In the current example, the Intel compiler is
! invoked with:
!
! -mieee-fp -init=arrays,zero -double-size 128 -real-size 128 -integer-size 32
!
      SUBROUTINE SIN
! precision using IMPLICIT REAL*16(A-H,O-Z) and constants defined as:
      IMPLICIT REAL*16(A-H,O-Z)
      DIMENSION BES(71),FAC(71),BESLOG(2)
      COMMON /C/ KMIN,KMAX
      COMMON /E/ XIN(26)

      REAL (SELECTED_REAL_KIND (32, 70) ) BES, FAC, BESLOG, XIN, XREEL
      REAL (SELECTED_REAL_KIND (32, 70) ) E, FAC1, XFN, FAC2, FAC3, FAC4
      REAL (SELECTED_REAL_KIND (32, 70) ) S, SNORM, SNORM1, A, A2, A3
      REAL (SELECTED_REAL_KIND (32, 70) ) FAC5, SUM, TERM, TEST

      DATA XREEL/1.0Q+00/
      DATA NLMAX/71/
      E=DEXP(1.0Q0)
      FAC1=0.5Q+00*(DLOG10(E)-DLOG10(2.0Q+00))
C     FAC(1)=0.0Q0
C     FAC(2)=0.0Q0
C     XFN=1.0Q0
C     DO 6 IFA=3,NLMAX
C     XFN=XFN+1.0Q0
C   6 FAC(IFA)=FAC(IFA-1)+DLOG(XFN)
! read in upper limit IONUM
      READ (5,31)IONUM
      WRITE(6,31)IONUM
! for the range 1 =< IO =< IONUM read factorial NF!
      DO 6 IO=1,IONUM
      READ (5,31)NF,FAC(IO)
    6 WRITE(6,31)NF,FAC(IO)
      NLMIN=NLMAX
! set argunment to 1.0
      FAC2=E*XREEL
      FAC2=DLOG10(FAC2)
      IN=NLMIN-1
      DO 13 I=IN,NLMIN
      N=I-1
      FAC3=I
      FAC4=N+I
      FAC5=N
   13 BESLOG(2-NLMIN+I)=FAC1-FAC3*DLOG10(FAC4)+FAC5*FAC2
! set first two highest orders in the recurrence
! following Eq.(38)
      BES(NLMIN)=10.0Q0**BESLOG(2)
      BES(NLMIN-1)=10.0Q0**BESLOG(1)
      S=NLMIN+NLMIN-1
      SNORM=S*BES(NLMIN)*BES(NLMIN)
      S=NLMIN+NLMIN-3
      SNORM=SNORM+S*BES(NLMIN-1)*BES(NLMIN-1)
! start accumulation for sigmna  in SNORM
! then perform downward recurrence
      DO 16 IO=2,IN
      JO=NLMIN-IO
      A=JO+JO+1
      A2=JO+JO-1
      BES(JO)=A*BES(JO+1)/XREEL-BES(JO+2)
   16 SNORM=SNORM+A2*BES(JO)*BES(JO)
! after recurrence compute normalization SNORM from recurrence
! and compare to exact value in SNORM1
      SNORM=1.0Q0/DSQRT(SNORM)
      SNORM1=DSIN(XREEL)/BES(1)
      WRITE(6,30)SNORM
      WRITE(6,30)SNORM1
! normalize recurrence values from recurrence to obtain value
! of Bessel Function with argument "1.0" for use in Eq.(37)
      DO 17 J=1,NLMIN
   17 BES(J)=SNORM*BES(J)
      A3=DSIN(XREEL)/XREEL
      A3=BES(1)/A3
! write value of the Bessel Function of lowest order 
! from both methods (they should be identical)
! then write value of the Bessel Function for all orders
! required in Eq.(37)
      WRITE(6,30)XREEL,A3
      WRITE(6,31)((IO,BES(IO)),IO=1,NLMAX)
! Accumulate sum in Eq.(37) and test for magnitude of last term
      DO 7 IO=KMIN,KMAX
      SUM=0.0Q0
      DO 8 IK=1,IONUM
      TERM=BES(IO+IK-1)/FAC(IK)
      TEST=TERM/BES(IO)
      IF(TEST.LT.0.10Q-22)GO TO 34
    8 SUM=SUM+TERM
   34 XIN(IO)=SUM+SUM
! check magnitude of terms XIN(:), TERM,  TEST values then return
      WRITE(6,1)IO
      WRITE(6,31)IK,XIN(IO),IK,TERM,IK,TEST
    7 CONTINUE
      RETURN
    2 WRITE(6,3)
      STOP
    1 FORMAT(1X ,I3)
    3 FORMAT(//,21H SHOULD NOT EXCEED 20)
   30 FORMAT(1X ,4D31.24)
   31 FORMAT(4(1X ,I2,D31.24))
   32 FORMAT(1X ,I4)
      END

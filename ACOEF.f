!
!   Copyright 1974-2020 George Delic, HiPERISM Consulting LLC
!
!   subroutine ACOEF
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
! PURPOSE: compute the quadrature weights A(k,j)
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
!
      SUBROUTINE ACOEF
! precision using IMPLICIT REAL*16(A-H,O-Z) and constants defined as:
      IMPLICIT REAL*16(A-H,O-Z)
      DIMENSION P(26)

      REAL (SELECTED_REAL_KIND (32, 70) ) P, ACOF, ZEROS, FZEROS, ZEROF
      REAL (SELECTED_REAL_KIND (32, 70) ) EINS, PAK, PZ, XIT, YIT, YIT1
      REAL (SELECTED_REAL_KIND (32, 70) ) XK

      COMMON /A/ ACOF(5850),IARTH(25),B(350),IBRTH(25),ZEROS(1326),IZRTH
     1(51),FZEROS(1326)
      COMMON /C/ KMIN,KMAX
      EINS=1.0Q0
      ISUM=0
      DO 1 IK=KMIN,KMAX
      K=IK-1
      ISUM=ISUM+K*IK
      IARTH(K)=ISUM-IK*IK
      XK=K
      N=K+1
      DO 2 I=1,K
      JIK=I*IK
      P(1)=1.0Q0
      IF(IK.GE.3)GO TO 8
      P(2)=ZEROS(1)
      GO TO 6
    8 P(2)=ZEROS(IZRTH(K)+I)
      DO 3 IT=3,K
      LIT=IT-1
      YIT=LIT
      XIT= EINS/YIT
      LIT1=LIT-1
      YIT1=LIT1
    3 P(IT)=((YIT+YIT1)*P(2)*P(LIT)-YIT1*P(LIT1))*XIT
    6 PAK=(1.0Q0-P(2)**2)/(XK*P(K))
      PZ=P(2)
      DO 4 IN=1,N
      INK=JIK+IN
      P(1)=1.0Q0
      P(2)=ZEROS(IZRTH(N+K)+IN)
      DO 5 IT=3,IK
      LIT=IT-1
      YIT=LIT
      XIT= EINS/YIT
      LIT1=LIT-1
      YIT1=LIT1
    5 P(IT)=((YIT+YIT1)*P(2)*P(LIT)-YIT1*P(LIT1))*XIT
      ZEROF=P(2)-PZ
      IF(ZEROF.NE.0.0Q+00)GO TO 7
      IF(P(IK).NE.0.0Q+00)GO TO 9
      ACOF(IARTH(K)+INK)=0.0Q+00
      GO TO 4
    7 ACOF(IARTH(K)+INK)=P(IK)*PAK/ZEROF
      GO TO 4
    9 ACOF(IARTH(K)+INK)=0.0Q+00
      WRITE(6,10)K,IN,I
   10 FORMAT(10X,4I4)
    4 CONTINUE
    2 CONTINUE
    1 CONTINUE
      RETURN
      END

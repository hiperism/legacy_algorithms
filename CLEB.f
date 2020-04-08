!
!   Copyright 1974-2020 George Delic, HiPERISM Consulting LLC
!
!   subroutine CLEB 
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
! PURPOSE: compute the Clebsch-Gordon coefficient C(ijk,000)
!
! AUTHOR:  George Delic, Ph.D., for this modification                
!          Reference equation Eq.(24) in
!          G. Delic,
!          The legendre series and a quadrature formula for its coefficients,
!          Journal of Computational Physics, vol. 14 (1974), pp. 254-268.
!
! ORIGIN:  based in part on source code provided by B.A. Robson
!          (Australina National University, Canberra, Australia) and
!          T. Tamura, ORNL Technical report ORNL-2877,
!          https://www.osti.gov/biblio/4593484
!
! LANGUAGE:  FOTRAN IV and Fortran 77 (ANSI X3J3) X3.9-1978
!            with minor additions from Fortran 90/95
!
!
! This code was developed on the Telefunken TR-440 with Fortran IV in 1973
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
! Coefficients that are printed out correspond to the following equations
! in the accompanying report IKDA73-8.pdf
!
      SUBROUTINE CLEB
! precision using IMPLICIT REAL*16(A-H,O-Z) and constants defined as:
      IMPLICIT REAL*16(A-H,O-Z)
      COMMON /B/ G(500),RAC,J1,J2,J3,M1,M2,M3
      DIMENSION I(11)

      REAL (SELECTED_REAL_KIND (32, 70) ) G, RAC
      REAL (SELECTED_REAL_KIND (32, 70) ) C, A, B, H, D, E, F, S, Q, T 

      EQUIVALENCE (I(1),I1),(I(2),I2),(I(3),I3),(I(4),I4),(I(5),I5),(I(6
     1),I6),(I(7),I7),(I(8),I8),(I(9),I9),(I(10),I10),(I(11),I11)
      RAC=0.0Q0
C  TEST M1+M2=M3
      IF(M1+M2-M3)300,40,300
C  TEST TABLE SIZE
   40 I(10)=(J1+J2+J3)/2+2
      N=I(10)
      I(11)=J3+2
      IF(I(10)-500)70,70,50
   50 WRITE (6,60) I(10),IM,J1,J2,J3,M1,M2,M3
      GO TO 300
   60 FORMAT(11H0TABLE SIZE ,6I5)
   70 I(1)=J1+J2-J3
      I(2)=J2+J3-J1
      I(3)=J3+J1-J2
      I(4)=J1-M1
      I(5)=J1+M1
      I(6)=J2-M2
      I(7)=J2+M2
      I(8)=J3-M3
      I(9)=J3+M3
C  CHECK I(J)=EVEN,TRIANGULAR INEQUALITY,M LESS J,FIND NO OF TERMS
      DO 110 J=1,9
      K=I(J)/2
      IF(I(J)-K-K)300,80,300
   80 IF(K)300,90,90
   90 IF(K-N)100,110,110
  100 N=K
  110 I(J)=K+1
      IF(M3)115,400,115
  115 IL=0
      LA=I1-I5
      LB=I1-I6
      IF(IL-LA)120,130,130
  120 IL=LA
  130 IF(IL-LB)140,145,145
  140 IL=LB
C  FORM COEFFICIENT OF SUM
  145 C=(G(I11)-G(I11-1)+G(I1)+G(I2)+G(I3)-G(I10)+G(I4)+G(I5)+G(I6)+G(I7
     1)+G(I8)+G(I9))/2.0Q0
      N1=I1-IL
      N2=I4-IL
      N3=I7-IL
      N4=IL+1
      N5=IL-LA+1
      N6=IL-LB+1
      C=C-G(N1)-G(N2)-G(N3)-G(N4)-G(N5)-G(N6)
      C=DEXP(C)
      IF(MOD(IL,2).EQ.0)GO TO 500
      C=-C
  500 IF(N)300,150,160
  150 RAC=C
      GO TO 300
C  FORM SUM
  160 A=N1-1
      B=N2-1
      H=N3-1
      D=N4
      E=N5
      F=N6
      S=1
      Q=N-1
      DO 170 J=1,N
      T=(A-Q)/(D+Q)*(B-Q)/(E+Q)*(H-Q)/(F+Q)
      S=1.0-S*T
      Q=Q-1.0
  170 CONTINUE
      RAC=C*S
C     return to caller here
  300 RETURN
C  SPECIAL FORMULA FOR M3=0 AND M1=0 OR 1/2
  400 IF(MOD(I10,2).EQ.0)GO TO 420
  410 K=1
      GO TO 430
  420 K=0
  430 IF(M1)115,440,460
  440 L=0
      IF(K)300,480,300
  460 IF(M1-1)115,470,115
  470 L=1
  480 X=L
      M=I3+(I1+K+1)/2-L
      N4=I10/2+K
      N5=I4+I5
      N6=I6+I7
      N1=(I1+1-K)/2
      N2=(I2+1+K-L)/2
      N3=(I3+1+K-L)/2
      RAC=DEXP((G(I11)-G(I11-1)+G(I1)+G(I2)+G(I3)-G(I10))/2.0Q0+G(N4)
     1-G(N1)-G(N2)-G(N3)+X*(G(3)-(G(N5)-G(N5-1)+G(N6)-G(N6-1))/2.0Q0))
      IF(MOD(M,2).EQ.0)GO TO 300
      RAC=-RAC
      GO TO 300
      END

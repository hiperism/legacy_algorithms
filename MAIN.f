!
!   Copyright 1974-2020 George Delic, HiPERISM Consulting LLC
!
!   program MAIN 
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
! PURPOSE: compute the quadrature weights and abscissas
!          using the zeros x(k) of the Legendre Polynomials
!          P(K), for 1 =< K =< 25
!          then apply to test example
!
! AUTHOR:  George Delic, Ph.D., for this source code
!          Reference equations Eq.(30,36) in
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
      PROGRAM MAIN
! precision using IMPLICIT REAL*16(A-H,O-Z) and constants defined as:
      IMPLICIT REAL*16(A-H,O-Z)
      DIMENSION GK(25),GKP(25)

      REAL (SELECTED_REAL_KIND (32, 70) ) GK, GKP, ACOF, B, ZEROS, XIN
      REAL (SELECTED_REAL_KIND (32, 70) ) FZEROS, XF1, SUM, SUM1, S, TS
      REAL (SELECTED_REAL_KIND (32, 70) ) TB, SUI, SU, DERM, TERM, SUM2
      REAL (SELECTED_REAL_KIND (32, 70) ) XT, YT, SAM, XF2

      COMMON /A/ ACOF(5850),IARTH(25),B(350),IBRTH(25),ZEROS(1326),IZRTH
     1(51),FZEROS(1326)
      COMMON /C/ KMIN,KMAX
      COMMON /E/XIN(26)
! read in K range - check the values are odd, exit if they are even
      READ (5,12) KMIN,KMAX
      WRITE(6,12) KMIN,KMAX
      IF(MOD(KMIN,2).EQ.0)GO TO 4
      IF(MOD(KMAX,2).EQ.0)GO TO 4
      NM=KMAX+KMAX+1
      KMIN=KMIN+1
      KMAX=KMAX+1
! compute terms in Eq.(37) with downward recurrence of Bessel Function
      CALL SIN
! If previously calculated then read in the weigths B and A
!     READ(60)IZRTH,IARTH,IBRTH,ZEROS,B,ACOF
!     GO TO 1001
! Otherwise compute weights B and A, and abscissas for K range
      DO 6 IIN=KMIN,KMAX
      IRIN=IIN-1
      IBRTH(IRIN)=IRIN*IRIN+IRIN
      IBRTH(IRIN)=IBRTH(IRIN)/2-1
      READ (5,52)IRINO
      IF(IRINO.NE.IRIN)GO TO 4
      IPAR=1
      B(IBRTH(IRIN)+IRIN+1)=0.0Q0
      IF(MOD(IRIN,2).NE.0)GO TO 7
      IPAR=IPAR+1
      IRIN=IRIN+1
!
! read in B coefficients
!
    7 READ(5,5)((IC,ACOF(IO)),IO=1,IRIN)
      DO 8 ION=1,IRIN
      GO TO (9,10),IPAR
C   9 B(IBRTH(IIN-1)+ION)=-ACOF(ION)
    9 B(IBRTH(IIN-1)+ION)=+ACOF(ION)
      GO TO 8
   10 B(IBRTH(IIN-1)+ION)=ACOF(ION)
    8 CONTINUE
    6 CONTINUE
!
! read in range of K
!
      READ (5,51)KMI9,KMA9
      WRITE(6,51)KMI9,KMA9
      IO1=2*KMA9+1
      IO2=KMA9
      DO 53 IO=1,IO1
      IZRTH(IO)=IO*(IO-1)
      IZRTH(IO)=IZRTH(IO)/2
      READ(5,52)KIN
!p    WRITE(6,52)KIN
      IF(KIN.NE.IO)GO TO 4
      IKO=(IO-1)/2+1
      IF(MOD(IO,2).EQ.0)GO TO 54
      IF(IO.EQ.1)GO TO 54
      IKO=IKO-1
      ZEROS(IZRTH(IO)+IKO+1)=0.0Q0
!
! read in zeros of Legendre Polynomial
!
   54 READ(5,5)((JJ,ZEROS(IZRTH(IO)+II)),II=1,IKO)
!
! optional print out
!
!p    WRITE(6,5)((JJ,ZEROS(IZRTH(IO)+II)),II=1,IKO)
      DO 55 II=1,IKO
   55 ZEROS(IZRTH(IO)+IO-II+1)=-ZEROS(IZRTH(IO)+II)
   53 CONTINUE
!
! call to compute A(k,j) coefficients in Eq. 31
!
      CALL ACOEF
 1001 CONTINUE
C     XF1=SQRT(2.0Q0)
C     XF1=XF1*XF1*XF1
      DO 1 K=1,NM
      DO 1 IN=1,K
      FZEROS(IZRTH(K)+IN)=DEXP(ZEROS(IZRTH(K)+IN))
C     FZEROS(IZRTH(K)+IN)=1.0D0/DSQRT(1.0D0-ZEROS(IZRTH(K)+IN))
    1 CONTINUE
!
! form sum in Eq.(30)
!
      DO 2 JO=KMIN,KMAX
      SUM=0.0Q0
      SUM1=0.0Q0
      SUM2=0.0Q0
      K=JO-1
      GKP(K)=0.0Q0
      WRITE(6,81)K
      N=K+K+1
      IPARK=1
      IF(MOD(K,2).NE.0)GO TO 20
      IPARK=IPARK+1
   20 CONTINUE
      GO TO (203,204),IPARK
  204 S=0.0Q0
      DO 216 J=1,K
      INK=K+1+J*JO
  216 S=S+ACOF(IARTH(K)+INK)*FZEROS(IZRTH(K)+J)
! form first term in Eq.(30)
      GKP(K)=GKP(K)+B(IBRTH(K)+K+1)*(FZEROS(IZRTH(N)+K+1)-S)
  203 CONTINUE
      N2=N+1
      DO 13 IN=1,K
!
! print out B (Eq. 32) coefficients
!
      WRITE(6,75)IN,B(IBRTH(K)+IN)
      S=0.0Q0
      GO TO (201,202),IPARK
  201 DO 214 J=1,K
      INK=IN+J*JO
      TS=FZEROS(IZRTH(K)+JO-J)-FZEROS(IZRTH(K)+J)
  214 S=S+ACOF(IARTH(K)+INK)*TS
      TB=FZEROS(IZRTH(N)+N2-IN)-FZEROS(IZRTH(N)+IN)-S
      GKP(K)=GKP(K)-B(IBRTH(K)+IN)*TB
      GO TO 200
  202 DO 215 J=1,K
      INK=IN+J*JO
      TS=FZEROS(IZRTH(K)+JO-J)+FZEROS(IZRTH(K)+J)
  215 S=S+ACOF(IARTH(K)+INK)*TS
      TB=FZEROS(IZRTH(N)+N2-IN)+FZEROS(IZRTH(N)+IN)-S
      GKP(K)=GKP(K)+B(IBRTH(K)+IN)*TB
  200 CONTINUE
! form second term in Eq.(30)
      SUI=0.0Q0
      SU=0.0Q0
      DO 14 J=1,K
      INK=IN+J*JO
!
! print A (Eq. 33) coefficients 
!
      WRITE(6,77)IN,J,ACOF(IARTH(K)+INK)
      SUI=SUI+ACOF(IARTH(K)+INK)*FZEROS(IZRTH(K)+JO-J)
   14 SU=SU+ACOF(IARTH(K)+INK)*FZEROS(IZRTH(K)+J)
      DERM=FZEROS(IZRTH(N)+IN)-SU
      TERM=  +B(IBRTH(K)+IN)*(FZEROS(IZRTH(N)+IN)-SU)

! uncomment for diagnostic print out
!     WRITE(6,3)IN,FZEROS(IZRTH(N)+IN),SU,DERM,TERM

      IF(TERM)60,61,62
   60 SUM1=SUM1+TERM
      GO TO 61
   62 SUM2=SUM2+TERM
   61 GO TO (18,19),IPARK
   18 TERM=  -B(IBRTH(K)+IN)*(FZEROS(IZRTH(N)+N2-IN)-SUI)
      GO TO 66
   19 TERM=  +B(IBRTH(K)+IN)*(FZEROS(IZRTH(N)+N2-IN)-SUI)
   66 DERM=FZEROS(IZRTH(N)+N2-IN)-SUI

! uncomment for diagnostic print out
!     WRITE(6,3)IN,FZEROS(IZRTH(N)+N2-IN),SUI,DERM,TERM

      IF(TERM)63,13,65
   63 SUM1=SUM1+TERM
      GO TO 13
   65 SUM2=SUM2+TERM
   13 CONTINUE
      GO TO (23,22),IPARK
   22 JIK=K+1
      WRITE(6,75)JIK,B(IBRTH(K)+JIK)
      SU=0.0Q0
      DO 21 J=1,K
      INK=JIK+J*JO
      WRITE(6,77)JIK,J,ACOF(IARTH(K)+INK)
   21 SU=SU+ACOF(IARTH(K)+INK)*FZEROS(IZRTH(K)+J)
      SUM=SUM+B(IBRTH(K)+K+1)*(FZEROS(IZRTH(N)+K+1)-SU)
      DERM=FZEROS(IZRTH(K)+K+1)-SU
      TERM=B(IBRTH(K)+K+1)*DERM

! uncomment for diagnostic print out
!     WRITE(6,3)IN,FZEROS(IZRTH(N)+K+1),SU,DERM,TERM

   23 SUM=SUM+SUM1+SUM2
      XT=SUM/XIN(JO)
      YT=DABS(XIN(JO)-SUM)
C     DO 15 J=1,K
C     SAM=0.0Q0
C     DO 16 IN=1,N
C     KJ=K*J
C  16 SAM=SAM+B(IBRTH(K)-IN)*ACOF(IARTH(K)-KJ)
C  15 GK(J)=SAM
C     WRITE(6,17)((JR,GK(JR)),JR=1,K)
C     WRITE(6,3)K,XIN(JO),SUM,XT,YT
C     WRITE(7,103)K,XIN(JO)
C     WRITE(7,103)K,SUM
C     WRITE(7,103)K,XT
C     WRITE(7,103)K,YT
! compare quadrature sum of Eq.(30) and exact result in Eq.(37)
!                            GKP(:)                 XIN(:)
!                            vvvvvvv                vvvvvv
! write cummulative terms in Eq.(25) and compare to XIN(K)
!                 K Eq.(37) SUM ratio
!     WRITE(6,103)K,XIN(JO),SUM,XT
    2 CONTINUE
! ratio       XT =        GKP(:) / XIN(:)
! difference  YT =  DABS( GKP(:) - XIN(:) )
      DO 217 IO=KMIN,KMAX
      K=IO-1
      XT=GKP(K)/XIN(IO)
      YT=DABS(XIN(IO)-GKP(K))
!
! print results of Table II
!
      WRITE(6,3) K,XIN(IO),GKP(K),XT,YT
  217 CONTINUE
!     STOP
      READ (5,12,END=4)KMIN,KMAX
      WRITE(6,12)KMIN,KMAX
      KM=KMAX+KMAX+1
      WRITE(6,112)
      WRITE(6,79)KMIN,KM
      DO 71 IO=1,KM
      WRITE(6,80)IO
!p    WRITE(7,80)IO
      IKO=(IO-1)/2+1
      IF(MOD(IO,2).EQ.0)GO TO 70
      IF(IO.EQ.1)GO TO 70
   70 CONTINUE
!p    WRITE(7,72)((II,ZEROS(IZRTH(IO)+II)),II=1,IKO)
   71 CONTINUE
      WRITE(6,112)
!
! print out B (Eq. 32) and A (Eq. 33) coefficients 
!
      DO 74 KK=KMIN,KMAX
      KO=KK
      IPARK=1
      WRITE(6,81) KK
!p    WRITE(7,81) KK
      IF(MOD(KK,2).NE.0)GO TO 73
      KO=KO+1
      IPARK=IPARK+1
   73 CONTINUE
      DO 74 IN=1,KO
!     WRITE(6,75) IN,B(IBRTH(KK)+IN)
      WRITE(7,75) IN,B(IBRTH(KK)+IN)
      DO 76 J=1,KK
      INK=IN+J*(KK+1)
!     WRITE(6,77)IN,J,ACOF(IARTH(KK)+INK)
      WRITE(7,77)IN,J,ACOF(IARTH(KK)+INK)
   76 CONTINUE
   74 CONTINUE
C     DO 2 JO=KMIN,KMAX
C     XIN(JO)=XF1/XF2
C     K=JO-1
C     XF2=K+K+1
C     XT=DABS(GK(JO)/XIN(JO))
C     WRITE(6,3)K,XIN(JO),GK(JO),XT
C   2 CONTINUE
    4 STOP
    3 FORMAT(1X ,I3,4D32.24)
  103 FORMAT(1X ,I3,4D28.20)
    5 FORMAT(2(1X,I4,D31.24))
   12 FORMAT(1X ,2I3)
   17 FORMAT(5(1X ,I4,D16.8))
   51 FORMAT(1X ,2I5)
   52 FORMAT(1X ,I4)
   72 FORMAT(10X,I3,D27.20)
   75 FORMAT(10X,4H  B(,I2,1H),2H =,D27.20)
   77 FORMAT(10X,2HA(,I2,1H,,I2,2H)=,D27.20)
!  78 FORMAT(1H0)
   78 FORMAT("\")
   79 FORMAT(10X,22HZEROS OF P(K) FOR K = ,I2,4H TO ,I2)
   80 FORMAT(10X,4HK = ,I2)
   81 FORMAT(/,15X,8HG ( K = ,I2,2H ))
  111 FORMAT(//)
  112 FORMAT(1X)
      END

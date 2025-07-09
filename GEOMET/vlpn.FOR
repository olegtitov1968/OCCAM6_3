      SUBROUTINE LAGINT (X,Y,N,XINT,YOUT)

C************************************************************************
C
C     SUBROUTINE LAGINT (X,Y,N,XINT,YOUT)
C
C*************************************************************
C
C ------- RECOMMENDED BY IERS !!!
C
C     THIS SUBROUTINE PERFORMS LAGRANGIAN INTERPOLATION
C     WITHIN A SET OF (X,Y) PAIRS TO GIVE THE Y
C     VALUE CORRESPONDING TO XINT.  THIS PROGRAM USES A
C     WINDOW OF 4 DATA POINTS TO PERFORM THE INTERPOLATION.
C     IF THE WINDOW SIZE NEEDS TO BE CHANGED, THIS CAN BE
C     DONE BY CHANGING THE INDICES IN THE DO LOOPS FOR
C     VARIABLES M AND J.
C
C     PARAMETERS ARE :
C     X     - ARRAY OF VALUES OF THE INDEPENDENT VARIABLE
C     Y     - ARRAY OF FUNCTION VALUES CORRESPONDING TO X
C     N     - NUMBER OF POINTS
C     XINT  - THE X-VALUE FOR WHICH ESTIMATE OF Y IS DESIRED
C     YOUT  - THE Y VALUE RETURNED TO CALLER
C
C     Distributed by Daniel Gambis (1997)
C
c   Last Revision: 12 April 2021 (OT)
C
C   Copyright (C) Oleg Titov, Geoscience Australia, 2001-2021
C
C   This program is free software: you can redistribute it and/or modify it
C   under the terms of the GNU General Public License as published by the
C   Free Software Foundation, either version 3 of the License, or (at your
C   option) any later version.
C
C   This program is distributed in the hope that it will be useful, but
C   WITHOUT ANY WARRANTY; without even the implied WARRANTY of 
C   MERCHANTABILITY of FITNESS FOR A PARTICULAR PURPOSE.
C   See the GNU General Public License for more details.
C
C   You should have received a copy of the GNU General Public License
C   along with the program. If not, see <http://www.gnu.org/license/>.
C
C

C
C
C*************************************************************
      REAL*8 X(N),Y(N),XINT,YOUT,TERM
      INTEGER N,I,J
C
      YOUT = 0.0D0
      DO 5 I = 1,N-1
        IF ( XINT .GE. X(I) .AND. XINT .LT. X(I+1) ) K = I
    5 CONTINUE
      IF ( K .LT. 2 ) K = 2
      IF ( K .GT. N-2 ) K = N-2
      DO 20 M = K-1,K+2
        TERM = Y(M)
        DO 10 J = K-1,K+2
          IF ( M .NE. J ) THEN
            TERM = TERM * (XINT - X(J))/(X(M) - X(J))
          END IF
   10   CONTINUE
        YOUT = YOUT + TERM
   20 CONTINUE
      RETURN
      END

************************************************************************

      SUBROUTINE EPHEM

C*************************************************************
C
C   SUBROUTINE EPHEM
C
C     THIS SUBROUTINE READS THE TABULAR VALUES FOR THE EARTH
C   ORIENTATION PARAMETERS AND EPHEMERIS FOR EARTH SUN AND MOON
C   FROM THE FILE "EPHEM.DAT"
C
C     IF EOP ARE OF 1-DAY SERIES, IT REMOVES THE SHORT PERIOD
C   VARIATIONS CORRECTION FROM TABULAR DATA OF UT1-UTC
C
C     IF "EPHEM.DAT" DOES NOT EXIST IT FLAGS AN ERROR MESSAGE
C
C   REVISION 1991 SEPTEMBER 30 (N.Z.L.)
C   REVISION 1993 APRIL 22  READING NEW EPH. FORMAT (N.Z.L.)
C   REVISION 1993 NOV. 17  INCREASE DIMENSIONS (N.Z.L.)
C   REVISION 1997 APRIL 23; (OT) READING NEW EOP1 FORMAT
C   REVISION August 4, 1998, (OT). Increasing of a priori EOP accuracy
C   REVISION January 31, 2001 (OT).
C   REVISION MARCH 11, 2001 (OT).   COMMON-BLOCK /LOG/ deleted
C   REVISION September 24, 2001 (OT).   The Sun geocentric velocities
C                 have been added to /GEOC/
C   REVISION June 19, 2002 (OT).  COMMON-BLOCK /H/ canceled
C   LAST REVISION June 04, 2010 (OT).  Time scale changed
C
C*************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION MJD
      CHARACTER LINE*250,BODY*5
      COMMON /BARIC/XC(30),YC(30),ZC(30),XV(30),YV(30),ZV(30),
     *              BCTIM(30),IPB
      COMMON /GEOC/ XG(30),YG(30),ZG(30),VXG(30),VYG(30),VZG(30),
     *              BTIM(30),IPG
      COMMON /MOON/ TMOON(30),RAM(6,30),DECM(6,30),PARM(6,30),
     *              XMON(30),YMON(30),ZMON(30),IPM

      COMMON  /G/   UTTJ, XP(30), YP(30), MJD(30), UT1KOR(30),
     *              NPAR, IP, JQA, IDELT

      COMMON  /NUT/  de(30), dp(30)
      COMMON /PHYS/ C, FL, AH, AU, J20
      COMMON /MATH/ PI, TWOPI, ROG, ROS, ROZ, TOL

      OPEN (16,FILE='EPHEM.DAT',STATUS='UNKNOWN',ERR=100)
      REWIND 16
      IP  = 0
      IPM = 0
      IPG = 0
      IPB = 0

1     CONTINUE

      READ (16,201,END=2) LINE
      READ (LINE,'(A5)') BODY
201   FORMAT (A250)

C----------------------------------------------------------------
C  CHECK IF THE LINE CORRESPONDS TO EARTH ROTATION PARAMETERS
C  ONE-DAY SERIES WITH SHORT PERIOD UT1 VARIATIONS REMOVED (EOP1)
C  AND READ JULIAN DATE, X WOBBLE, Y WOBBLE AND UT1-UTC
C----------------------------------------------------------------

      IF (BODY.EQ.'EOP1 ') THEN
         IP = IP+1
         READ (LINE,211)
     *        MJD(IP),XP(IP),YP(IP),UT1KOR(IP),idelt,dp(ip),de(ip)

c211      FORMAT (6X,F7.1,1x,2(F8.6),1x,F10.7,i4,f8.5,f8.5)

c
c    Change of format (OT) - 04.01.2008
c

211      FORMAT (6X,F7.1,1x,2(F8.6),1x,F10.7,i4,f9.6,f9.6)

C----------------------------------------------------------------
C MJD IS SUPPOSED TO BE IN TDT, WHICH IS NOT VERY CLEAR. PROBABLY
C IT IS REALLY IN UTC, BUT IT MUST BE CHECKED. CALC 7.4 ASSUMES IT
C TO BE CT AND THEN INTERPOLATE ALSO IN CT. THE REAL OCCAM VERSION
C SHOUL TAKE IT AS UTC AND INTERPOLATE AS UTC, AT LEAST UNTIL THIS
C POINT IS CROSS CHECKED PROPERLY.
C FOR THE TEST, OCCAM WILL FOLLOW CALC STRATEGY FOR A WHILE.
C----------------------------------------------------------------

c+  04.06.2010

         DTDT = MJD(IP)+(32.184D0 + IDELT)/86400.D0
c         DTDT = MJD(IP)
c-  04.06.2010

         D2000 = (DTDT - J20 - 0.5D0) / 36525.D0

c
c    Subtraction of short frequencies < 35d before interpolation
c

c+  04.06.2010
c
c         CALL CUT1VR (D2000,DELUT1)
         CALL CUT1VS (D2000,DELUT1)

c-  04.06.2010

         UT1KOR(IP) = UT1KOR(IP) - DELUT1

         GO TO 1
      ENDIF


C----------------------------------------------------------------
C   CHECK IF BODY IS MOON AND GET X,Y,Z COORDINATES AT J2000
C  (PEP EPHEMERIS DATA)
C----------------------------------------------------------------

      IF (BODY.EQ.'MOOND') THEN
         IPM=IPM+1
         READ (LINE,415) TMOON(IPM),XMON(IPM),YMON(IPM),ZMON(IPM)
415      FORMAT (6X,F7.1,3D24.16)
         GO TO 1
      ENDIF

C----------------------------------------------------------------
C   CHECK IF BODY IS EARTH AND GET VALUES FOR EACH DAY OF
C  EARTH BARICENTRIC COORDINATES (PEP EPHEMERIS)
C----------------------------------------------------------------

      IF (BODY.EQ.'EARTD') THEN
         IPB = IPB+1
         READ (LINE,400) BCTIM(IPB),XC(IPB),YC(IPB),ZC(IPB),
     *    XV(IPB),YV(IPB),ZV(IPB)
400      FORMAT (6X,F7.1,6D24.16)

         GO TO 1
      ENDIF

C----------------------------------------------------------------
C   CHECK IF BODY IS SUN AND GET GEOCENTRIC COORDINATES FOR
C  EACH DAY  (PEP EPHEMERIS)
C----------------------------------------------------------------

      IF (BODY.EQ.'SUND ') THEN
         IPG = IPG+1
         READ (LINE,440) BTIM(IPG),XG(IPG),YG(IPG),ZG(IPG),
     *    VXG(IPG), VYG(IPG), VZG(IPG)

440      FORMAT (6X,F7.1,6D24.16)

         GO TO 1
      ENDIF

2     CONTINUE


      CLOSE (16)
      RETURN

100   write (*,*) ' ***** ERROR in file EPHEM.DAT ***** '

      RETURN
      END


************************************************************************
C     CORRECTIONS OF THE HOUR ANGLE
***********************************************************************

      SUBROUTINE CUT1VR (D2000, DELUT1)
C
C-----------------------------------------------------------------------
C CORRECTION FOR PERIODIC UT1-VARIATIONS
C REFERENCES: YODER ET AL. 1981
C THEORETICAL VALUE K/C=0.944 USED (=CK)
C (HAS ALREADY BEEN CONSIDERED IN AMP !!!!! )
C
C PROGRAMMER: H. SCHUH, 1982, N. ZARRAOA, 1990
C
C THE COMPUTATION OF ANGULAR ELEMENTS HAS BEEN MODIFIED TO FULFILL THE
C IERS RECOMMENDATION (1989). (NZL)
C THE ALGORITHM IMPLEMENTED IN THE DERIVATION OF NUTATION HAS BEEN
C USED HERE. (NZL, NOV 90)
C-----------------------------------------------------------------------
C  LAST REVISION  1990 NOVEMBER 30   (N.Z.L.)
C
      IMPLICIT DOUBLE PRECISION  (A-H,O-Z)
C
      COMMON /PHYS/ C, FL, AH, AU, J20
      COMMON /MATH/ PI, TWOPI, ROG, ROS, ROZ, TOL

      DIMENSION  AMARG(5,41), AMP(41), B(5,5), A(5)
C
      DATA AMARG /
     A 1., 0., 2., 2., 2.,  2., 0., 2., 0., 1.,  2., 0., 2., 0., 2.,
     B 0., 0., 2., 2., 1.,  0., 0., 2., 2., 2.,  1., 0., 2., 0., 0.,
     C 1., 0., 2., 0., 1.,  1., 0., 2., 0., 2.,  3., 0., 0., 0., 0.,
     D-1., 0., 2., 2., 1., -1., 0., 2., 2., 2.,  1., 0., 0., 2., 0.,
     E 2., 0., 2.,-2., 2.,  0., 1., 2., 0., 2.,  0., 0., 2., 0., 0.,
     F 0., 0., 2., 0., 1.,  0., 0., 2., 0., 2.,  2., 0., 0., 0.,-1.,
     G 2., 0., 0., 0., 0.,  2., 0., 0., 0., 1.,  0.,-1., 2., 0., 2.,
     H 0., 0., 0., 2.,-1.,  0., 0., 0., 2., 0.,  0., 0., 0., 2., 1.,
     I 0.,-1., 0., 2., 0.,  1., 0., 2.,-2., 1.,  1., 0., 2.,-2., 2.,
     J 1., 1., 0., 0., 0., -1., 0., 2., 0., 0., -1., 0., 2., 0., 1.,
     K-1., 0., 2., 0., 2.,  1., 0., 0., 0.,-1.,  1., 0., 0., 0., 0.,
     L 1., 0., 0., 0., 1.,  0., 0., 0., 1., 0.,  1.,-1., 0., 0., 0.,
     M-1., 0., 0., 2.,-1., -1., 0., 0., 2., 0., -1., 0., 0., 2., 1.,
     N 1., 0.,-2., 2.,-1., -1.,-1., 0., 2., 0.  /
C
      DATA AMP /
     A -0.02D0, -0.04D0, -0.10D0, -0.05D0, -0.12D0, -0.04D0, -0.41D0,
     B -0.99D0, -0.02D0, -0.08D0, -0.20D0, -0.08D0,  0.02D0,  0.03D0,
     C -0.30D0, -3.21D0, -7.76D0,  0.02D0, -0.34D0,  0.02D0, -0.02D0,
     D  0.05D0, -0.73D0, -0.05D0, -0.05D0,  0.05D0,  0.10D0,  0.04D0,
     E  0.05D0,  0.18D0,  0.44D0,  0.53D0, -8.26D0,  0.54D0,  0.05D0,
     F -0.06D0,  0.12D0, -1.82D0,  0.13D0,  0.02D0, -0.09D0 /
C
      DATA B /
     A  485866.733D0, 1325.D0,  715922.633D0,  31.310D0,  0.064D0,
     B 1287099.804D0,   99.D0, 1292581.224D0,  -0.577D0, -0.012D0,
     C  335778.877D0, 1342.D0,  295263.137D0, -13.257D0,  0.011D0,
     D 1072261.307D0, 1236.D0, 1105601.328D0,  -6.891D0,  0.019D0,
     E  450160.280D0,   -5.D0, -482890.539D0,   7.455D0,  0.008D0 /

      SER = 1.296D6

      DCY  = D2000
      DCY2 = DCY * DCY
      DCY3 = DCY * DCY2

C----------------------------------
C   COMPUTE ANGULAR ELEMENTS
C----------------------------------

      DO I = 1, 5
         A(I) = ((B(5,I)*DCY+B(4,I))*DCY+B(3,I)+B(2,I)*SER)*DCY+B(1,I)
         A(I) = A(I)/ROS
         A(I) = DMOD(A(I),TWOPI)
      end do

      DELUT1 = 0.D0

      EL = A(1)
      ELSTRI = A(2)
      F = A(3)
      D = A(4)
      O = A(5)
C--------------------------------------------------------------
C  COMPUTE THE SUM OF THE CONTRIBUTION OF EACH COMPONENT
C--------------------------------------------------------------

      DO  I = 1, 41
         ARG =   AMARG(1,I) * EL + AMARG(2,I) * ELSTRI + AMARG(3,I) * F
     *         + AMARG(4,I) * D  + AMARG(5,I) * O
         DELUT1 = DELUT1 + AMP(I) * DSIN(ARG)
      end do

C CONVERSION TO MSEC

      DELUT1 = DELUT1 * 1.D-1

C--------------------------------------------------------
C CONVERT DELUT1 FROM MSEC TO SEC
C--------------------------------------------------------

      DELUT1 = (DELUT1 * 1.D-3)

      RETURN
      END

************************************************************************

      SUBROUTINE CUT1VS (D2000, DELUT1)
C
C-----------------------------------------------------------------------
C CORRECTION FOR PERIODIC UT1-VARIATIONS : IERS STANDARDS 1992
C REFERENCES: YODER ET AL. 1981
C             BROSCH ET AL. 1989
C
C THIS ROUTINE IS SIMILAR TO THE PREVIOUS ONE, BUT COLLECTS BOTH
C THE SOLID EARTH AND OCEANIC TERMS.
C
C *** PENDING ***
C CHECK IF THIS ROUTINE CAN BE USED ALSO WHEN USNO UT1-UTC IS USED
C CHECK IF THE COMPONENTS OF THE NEXT ROUTINE ARE INCLUDED OR NOT
C ***************
C
C PROGRAMMER: N. ZARRAOA, 1993
C
C-----------------------------------------------------------------------
C REVISION  1993 MARCH 5   (N.Z.L.)
C
C REVISION June 19, 2002 (OT).  COMMON-BLOCK /H/ canceled
C LAST REVISION  2006 MAY  17  (OT)
C
      IMPLICIT DOUBLE PRECISION  (A-H,O-Z)
C
      COMMON /PHYS/ C, FL, AH, AU, J20
      COMMON /MATH/ PI, TWOPI, ROG, ROS, ROZ, TOL

      DIMENSION  AMARG(5,62), AMPSIN(62), AMPCOS(62), B(5,5), A(5)
      DIMENSION  AMARG1(5,32),AMARG2(5,30)
      EQUIVALENCE (AMARG1(1,1),AMARG(1,1)),(AMARG2(1,1),AMARG(1,33))
C
      DATA AMARG1 /
     A 1., 0., 2., 2., 2.,  2., 0., 2., 0., 1.,  2., 0., 2., 0., 2.,
     B 0., 0., 2., 2., 1.,  0., 0., 2., 2., 2.,  1., 0., 2., 0., 0.,
     C 1., 0., 2., 0., 1.,  1., 0., 2., 0., 2.,  3., 0., 0., 0., 0.,
     D-1., 0., 2., 2., 1., -1., 0., 2., 2., 2.,  1., 0., 0., 2., 0.,
     E 2., 0., 2.,-2., 2.,  0., 1., 2., 0., 2.,  0., 0., 2., 0., 0.,
     F 0., 0., 2., 0., 1.,  0., 0., 2., 0., 2.,  2., 0., 0., 0.,-1.,
     G 2., 0., 0., 0., 0.,  2., 0., 0., 0., 1.,  0.,-1., 2., 0., 2.,
     H 0., 0., 0., 2.,-1.,  0., 0., 0., 2., 0.,  0., 0., 0., 2., 1.,
     I 0.,-1., 0., 2., 0.,  1., 0., 2.,-2., 1.,  1., 0., 2.,-2., 2.,
     J 1., 1., 0., 0., 0., -1., 0., 2., 0., 0., -1., 0., 2., 0., 1.,
     K-1., 0., 2., 0., 2.,  1., 0., 0., 0.,-1./

      DATA AMARG2 /                              1., 0., 0., 0., 0.,
     L 1., 0., 0., 0., 1.,  0., 0., 0., 1., 0.,  1.,-1., 0., 0., 0.,
     M-1., 0., 0., 2.,-1., -1., 0., 0., 2., 0., -1., 0., 0., 2., 1.,
     N 1., 0.,-2., 2.,-1., -1.,-1., 0., 2., 0.,  0., 2., 2.,-2., 2.,
     O 0., 1., 2.,-1., 1.,  0., 1., 2.,-1., 2.,  0., 0., 2.,-1., 0.,
     P 0., 0., 2.,-2., 1.,  0., 0., 2.,-2., 2.,  0., 2., 0., 0., 0.,
     Q 2., 0., 0.,-2.,-1.,  2., 0., 0.,-2., 0.,  2., 0., 0.,-2., 1.,
     R 0.,-1., 2.,-2., 1.,  0., 1., 0., 0.,-1.,  0.,-1., 2.,-1., 2.,
     S 0., 1., 0., 0., 0.,  0., 1., 0., 0., 1.,  1., 0., 0.,-1., 0.,
     T 2., 0.,-2., 0., 0., -2., 0., 2., 0., 1., -1., 1., 0., 1., 0.,
     U 0., 0., 0., 0., 2.,  0., 0., 0., 0., 1. /

C
C  DATA AMPSIN AND AMPCOS FROM IERS CONVENSTIONS 2003
C
      DATA AMPSIN /
     A -0.02D0, -0.04D0, -0.10D0, -0.05D0, -0.12D0, -0.04D0, -0.41D0,
     B -1.00D0, -0.02D0, -0.08D0, -0.20D0, -0.08D0,  0.02D0,  0.03D0,
     C -0.30D0, -3.22D0, -7.79D0,  0.02D0, -0.34D0,  0.02D0, -0.02D0,
     D  0.05D0, -0.74D0, -0.05D0, -0.05D0,  0.05D0,  0.10D0,  0.04D0,
     E  0.05D0,  0.18D0,  0.44D0,  0.54D0, -8.33D0,  0.55D0,  0.05D0,
     F -0.06D0,  0.12D0, -1.84D0,  0.13D0,  0.02D0, -0.09D0, -0.06D0,
     G  0.03D0, -1.91D0,  0.26D0,  1.18D0,-49.06D0, -0.20D0,  0.05D0,
     H -0.56D0,  0.04D0, -0.05D0,  0.09D0,  0.82D0,-15.65D0, -0.14D0,
     I  0.03D0, -0.14D0,  0.43D0, -0.04D0,  8.20D0, -1689.54D0 /

      DATA AMPCOS /
     A  0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,
     B  0.01D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,
     C  0.00D0,  0.02D0,  0.05D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,
     D  0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,
     E  0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.06D0,  0.00D0,  0.00D0,
     F  0.00D0,  0.00D0,  0.01D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,
     G  0.00D0,  0.02D0,  0.00D0, -0.01D0,  0.43D0,  0.00D0,  0.00D0,
     H  0.01D0,  0.00D0,  0.00D0,  0.00D0, -0.01D0,  0.15D0,  0.00D0,
     I  0.00D0,  0.00D0, -0.01D0,  0.00D0,  0.11D0, -25.04D0 /

C
      DATA B /
     A  485866.733D0, 1325.D0,  715922.633D0,  31.310D0,  0.064D0,
     B 1287099.804D0,   99.D0, 1292581.224D0,  -0.577D0, -0.012D0,
     C  335778.877D0, 1342.D0,  295263.137D0, -13.257D0,  0.011D0,
     D 1072261.307D0, 1236.D0, 1105601.328D0,  -6.891D0,  0.019D0,
     E  450160.280D0,   -5.D0, -482890.539D0,   7.455D0,  0.008D0 /

      SER = 1.296D6

      DCY  = D2000
      DCY2 = DCY * DCY
      DCY3 = DCY * DCY2

C----------------------------------
C   COMPUTE ANGULAR ELEMENTS
C----------------------------------

      DO  I = 1, 5
         A(I) = ((B(5,I)*DCY+B(4,I))*DCY+B(3,I)+B(2,I)*SER)*DCY+B(1,I)
         A(I) = A(I)/ROS
         A(I) = DMOD(A(I),TWOPI)
      END DO

      DELUT1 = 0.D0

      EL = A(1)
      ELSTRI = A(2)
      F = A(3)
      D = A(4)
      O = A(5)

C--------------------------------------------------------------
C  COMPUTE THE SUM OF THE CONTRIBUTION OF EACH COMPONENT
C--------------------------------------------------------------

      DO  I = 1, 62
         ARG =   AMARG(1,I) * EL + AMARG(2,I) * ELSTRI + AMARG(3,I) * F
     *         + AMARG(4,I) * D  + AMARG(5,I) * O
         DELUT1 = DELUT1 + AMPSIN(I) * DSIN(ARG) + AMPCOS(I)*DCOS(ARG)
      END DO

C CONVERSION TO MSEC

      DELUT1 = DELUT1 * 1.D-1

C--------------------------------------------------------
C CONVERT DELUT1 FROM MSEC TO SEC
C--------------------------------------------------------

      DELUT1 = (DELUT1 * 1.D-3)

      RETURN
      END

************************************************************************
C    END OF THE CORRECTIONS TO HOUR ANGLE
************************************************************************

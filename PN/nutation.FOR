*********************************************************************

      SUBROUTINE MP2000 (ALPHA, DELTA, D50IN, EQE, m_eph, idatj,
     *  utj, RQU, JMAP, sn, OFFN, dx, dy, 
     *  RAAP, DEAP, PREC, PR, VR, IDELT, NUTMOD)
C
C-----------------------------------------------------------------------
C PROGRAMMERS: S. HEITZ, U. GREPEL, H. SCHUH, J. CAMPBELL, FEB. 84
C              N. ZARRAOA, OCT. 91
C-----------------------------------------------------------------------
C
C   APPLY PRECESION AND NUTATION TO 2000.0 SOURCE COORDINATES AND
C  TRANSFORM THEM TO APPARENT DATE COORDINATES
C
C  REVISION   1991 OCTOBER 25  (N.Z.L.) (CORRECTION IF NO EPHEMERIS
C                  AVAILABLE FOR SUN COORDINATES)
C  REVISION   1993 FEBRUARY 13 (NZL) INCLUSION OF GEODESIC NUTATION
C  REVISION   1993 NOVEMBER 17 (NZL) INCREASED DIMENSIONS
C  REVISION   1999 FEBRUARY 10 (OT)  HERRING MODEL FOR NUTATION
C  REVISION   2000 DECEMBER 09 (OT)
C  Revision:  2001  MARCH, 11 (OT)   COMMON-BLOCK /LOG/ deleted
C                                        MHB 2000 model added
C  Revision:  2001  September, 24 (OT)   Geocentric velocities of the SUN
C                  are from EPHEM.DAT, changes in /GEOC/ also
C  Revision:  2002  January, 11 (OT)  Geodetic nutation is applicable
C                  only for IAU 1980 model (remark by Z.Malkin)
C  Revision:  2005  December, 13 (OT)  Transition to the mean system of date
C                  Herring's nutation model is off
C  Revision:  2006  January, 04 (OT)  Direct call of DELE405 is added.
C                                      Options for ephemerides approach is added
C                                      Jupiter is added to /BAR/
C  Last revision: 2007 February, 08 (OT) Second term of aberrational effect is updated
C-----------------------------------------------------------------------
C 1. DECLARATIONS
C-----------------------------------------------------------------------

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER  RA(6)

      DIMENSION A(7), RW(6), B(5,5), DVELB(3), DCORB(3), XS(3), VS(3),
     1           RQU(3), DCG(3), DVG(3), EP(3,3), XLUN(3), R2(3), XM(3)
      DIMENSION XW(9,106),X1(9,20),X2(9,20),X3(9,20),X4(9,20),X5(9,20),
     *          X6(9,6), sn(3,3), OFFN(3), E(3,3), PR(3), PREC(3,3)
      DIMENSION VR(3), VR2(3), D2(3), PN(3,3)
      DIMENSION rrde(6), rrds(6), rrdm(6), rrdjup(6)              ! 04.01.2006

      real*8 dpsi_fcn, deps_fcn

C  STATEMENT SUBSTITUTED BY A DO-LOOP ASSIGNEMENT

c
c  Change on 3.11.2005
c
      COMMON  /BAR/ DCORBN(3,3),DVELBN(3,3),DCORG(3,3),
     *              DVELG(3,3), CMOON(3,3), xjup(3,3)            ! 13.07.2009
      COMMON /PHYS/ C, FL, AH, AU, J20
      COMMON /MATH/ PI, TWOPI, ROG, ROS, ROZ, TOL


C   COMMON WITH BARYCENTRIC DATA FROM ASTRONOMICAL ALMANAC (2000.0)

      COMMON /BARIC/ XBRC(30),YBRC(30),ZBRC(30),XBRV(30),YBRV(30),
     1  ZBRV(30), BCTIM(30),IP1
      COMMON /GEOC/ XG(30),YG(30),ZG(30),VXG(30),VYG(30),VZG(30),
     * BTIM(30),IP2

      COMMON  /NUT/  de(30), dp(30)

     

C-----------------------------------------------------------------------
C 2. DATA STATEMENTS
C-----------------------------------------------------------------------

c      open (2, file='nutation.dat', status='unknown')
C
C 2.1  AXES OF ROTATION FOR THE TRANSFORMATION BETWEEM MEAN EQUINOX OF
C      EPOCH AND TRUE EQINOX OF DATE
C-----------------------------------------------------------------------

      DATA RA / 3,2,3,1,3,1 /

C
C***********************************************************************
C
C               1980 IAU THEORY OF NUTATION (WAHR THEORY)
C           TABLE OF MULTIPLES OF ARGUMENTS AND COEFFICIENTS
C
C                   MULTIPLE OF            LONGITUDE        OBLIQUITY
C              L    L'   F    D  OMEGA   COEFF. OF SIN    COEFF. OF COS
C**********************************************************************

      DATA X1/0.d0,0.d0,0.d0,0.d0,1.d0,-171996.d0,-174.2d0,92025.d0,
     / 8.9d0, 0.d0, 0.d0, 2.d0,-2.d0, 2.d0,-13187.d0,-1.6d0,5736.d0,
     /-3.1d0, 0.d0, 0.d0, 2.d0, 0.d0, 2.d0,-2274.d0, -0.2d0, 977.d0,
     /-0.5d0, 0.d0, 0.d0, 0.d0, 0.d0, 2.d0, 2062.d0, 0.2d0,-895.d0,
     / 0.5d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1426.d0,-3.4d0, 54.d0,
     /-0.1d0, 1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 712.d0, 0.1d0, -7.d0, 0.d0,
     /      0.d0, 1.d0, 2.d0,-2.d0, 2.d0,-517.d0, 1.2d0,224.d0, -0.6d0,
     /      0.d0, 0.d0, 2.d0, 0.d0, 1.d0,-386.d0,-0.4d0,200.d0,  0.0d0,
     /      1.d0, 0.d0, 2.d0, 0.d0, 2.d0,-301.d0, 0.0d0,129.d0, -0.1d0,
     /      0.d0,-1.d0, 2.d0,-2.d0, 2.d0, 217.d0,-0.5d0,-95.d0,  0.3d0,
     /      1.d0, 0.d0, 0.d0,-2.d0, 0.d0,-158.d0, 0.0d0, -1.d0,  0.0d0,
     /      0.d0, 0.d0, 2.d0,-2.d0, 1.d0, 129.d0, 0.1d0,-70.d0,  0.0d0,
     /     -1.d0, 0.d0, 2.d0, 0.d0, 2.d0, 123.d0, 0.0d0,-53.d0,  0.0d0,
     /      1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 63.d0, 0.1d0, -33.d0,  0.0d0,
     /      0.d0, 0.d0, 0.d0, 2.d0, 0.d0, 63.d0, 0.0d0, -2.d0,  0.0d0,
     /     -1.d0, 0.d0, 2.d0, 2.d0, 2.d0,-59.d0, 0.0d0, 26.d0,  0.0d0,
     /     -1.d0, 0.d0, 0.d0, 0.d0, 1.d0,-58.d0,-0.1d0, 32.d0,  0.0d0,
     /      1.d0, 0.d0, 2.d0, 0.d0, 1.d0,-51.d0, 0.0d0, 27.d0,  0.0d0,
     /      2.d0, 0.d0, 0.d0,-2.d0, 0.d0, 48.d0, 0.0d0,  1.d0,  0.0d0,
     /     -2.d0, 0.d0, 2.d0, 0.d0, 1.d0, 46.d0, 0.0d0,-24.d0,  0.0d0/
      DATA X2/ 0.d0,0.d0,2.d0,2.d0, 2.d0,-38.d0, 0.0d0, 16.d0,  0.0d0,
     /      2.d0, 0.d0, 2.d0, 0.d0, 2.d0,-31.d0, 0.0d0, 13.d0,  0.0d0,
     /      2.d0, 0.d0, 0.d0, 0.d0, 0.d0, 29.d0, 0.0d0, -1.d0,  0.0d0,
     /      1.d0, 0.d0, 2.d0,-2.d0, 2.d0, 29.d0, 0.0d0,-12.d0,  0.0d0,
     /      0.d0, 0.d0, 2.d0, 0.d0, 0.d0, 26.d0, 0.0d0, -1.d0,  0.0d0,
     /      0.d0, 0.d0, 2.d0,-2.d0, 0.d0,-22.d0, 0.0d0,  0.d0,  0.0d0,
     /     -1.d0, 0.d0, 2.d0, 0.d0, 1.d0, 21.d0, 0.0d0,-10.d0,  0.0d0,
     /      0.d0, 2.d0, 0.d0, 0.d0, 0.d0, 17.d0,-0.1d0,  0.d0,  0.0d0,
     /      0.d0, 2.d0, 2.d0,-2.d0, 2.d0,-16.d0, 0.1d0,  7.d0,  0.0d0,
     /     -1.d0, 0.d0, 0.d0, 2.d0, 1.d0, 16.d0, 0.0d0, -8.d0,  0.0d0,
     /      0.d0, 1.d0, 0.d0, 0.d0, 1.d0,-15.d0, 0.0d0,  9.d0,  0.0d0,
     /      1.d0, 0.d0, 0.d0,-2.d0, 1.d0,-13.d0, 0.0d0,  7.d0,  0.0d0,
     /      0.d0,-1.d0, 0.d0, 0.d0, 1.d0,-12.d0, 0.0d0,  6.d0,  0.0d0,
     /      2.d0, 0.d0,-2.d0, 0.d0, 0.d0, 11.d0, 0.0d0,  0.d0,  0.0d0,
     /     -1.d0, 0.d0, 2.d0, 2.d0, 1.d0,-10.d0, 0.0d0,  5.d0,  0.0d0,
     /      1.d0, 0.d0, 2.d0, 2.d0, 2.d0, -8.d0, 0.0d0,  3.d0,  0.0d0,
     /      0.d0,-1.d0, 2.d0, 0.d0, 2.d0, -7.d0, 0.0d0,  3.d0,  0.0d0,
     /      0.d0, 0.d0, 2.d0, 2.d0, 1.d0, -7.d0, 0.0d0,  3.d0,  0.0d0,
     /      1.d0, 1.d0, 0.d0,-2.d0, 0.d0, -7.d0, 0.0d0,  0.d0,  0.0d0,
     /      0.d0, 1.d0, 2.d0, 0.d0, 2.d0,  7.d0, 0.0d0, -3.d0,  0.0d0/
      DATA X3/-2.d0,0.d0,0.d0,2.d0, 1.d0, -6.d0, 0.0d0,  3.d0,  0.0d0,
     /      0.d0, 0.d0, 0.d0, 2.d0, 1.d0, -6.d0, 0.0d0,  3.d0,  0.0d0,
     /      2.d0, 0.d0, 2.d0,-2.d0, 2.d0,  6.d0, 0.0d0, -3.d0,  0.0d0,
     /      1.d0, 0.d0, 0.d0, 2.d0, 0.d0,  6.d0, 0.0d0,  0.d0,  0.0d0,
     /      1.d0, 0.d0, 2.d0,-2.d0, 1.d0,  6.d0, 0.0d0, -3.d0,  0.0d0,
     /      0.d0, 0.d0, 0.d0,-2.d0, 1.d0, -5.d0, 0.0d0,  3.d0,  0.0d0,
     /      0.d0,-1.d0, 2.d0,-2.d0, 1.d0, -5.d0, 0.0d0,  3.d0,  0.0d0,
     /      2.d0, 0.d0, 2.d0, 0.d0, 1.d0, -5.d0, 0.0d0,  3.d0,  0.0d0,
     /      1.d0,-1.d0, 0.d0, 0.d0, 0.d0,  5.d0, 0.0d0,  0.d0,  0.0d0,
     /      1.d0, 0.d0, 0.d0,-1.d0, 0.d0, -4.d0, 0.0d0,  0.d0,  0.0d0,
     /      0.d0, 0.d0, 0.d0, 1.d0, 0.d0, -4.d0, 0.0d0,  0.d0,  0.0d0,
     /      0.d0, 1.d0, 0.d0,-2.d0, 0.d0, -4.d0, 0.0d0,  0.d0,  0.0d0,
     /      1.d0, 0.d0,-2.d0, 0.d0, 0.d0,  4.d0, 0.0d0,  0.d0,  0.0d0,
     /      2.d0, 0.d0, 0.d0,-2.d0, 1.d0,  4.d0, 0.0d0, -2.d0,  0.0d0,
     /      0.d0, 1.d0, 2.d0,-2.d0, 1.d0,  4.d0, 0.0d0, -2.d0,  0.0d0,
     /      1.d0, 1.d0, 0.d0, 0.d0, 0.d0, -3.d0, 0.0d0,  0.d0,  0.0d0,
     /      1.d0,-1.d0, 0.d0,-1.d0, 0.d0, -3.d0, 0.0d0,  0.d0,  0.0d0,
     /     -1.d0,-1.d0, 2.d0, 2.d0, 2.d0, -3.d0, 0.0d0,  1.d0,  0.0d0,
     /      0.d0,-1.d0, 2.d0, 2.d0, 2.d0, -3.d0, 0.0d0,  1.d0,  0.0d0,
     /      1.d0,-1.d0, 2.d0, 0.d0, 2.d0, -3.d0, 0.0d0,  1.d0,  0.0d0/
      DATA X4/ 3.d0,0.d0,2.d0,0.d0, 2.d0, -3.d0, 0.0d0,  1.d0,  0.0d0,
     /     -2.d0, 0.d0, 2.d0, 0.d0, 2.d0, -3.d0, 0.0d0,  1.d0,  0.0d0,
     /      1.d0, 0.d0, 2.d0, 0.d0, 0.d0,  3.d0, 0.0d0,  0.d0,  0.0d0,
     /     -1.d0, 0.d0, 2.d0, 4.d0, 2.d0, -2.d0, 0.0d0,  1.d0,  0.0d0,
     /      1.d0, 0.d0, 0.d0, 0.d0, 2.d0, -2.d0, 0.0d0,  1.d0,  0.0d0,
     /     -1.d0, 0.d0, 2.d0,-2.d0, 1.d0, -2.d0, 0.0d0,  1.d0,  0.0d0,
     /      0.d0,-2.d0, 2.d0,-2.d0, 1.d0, -2.d0, 0.0d0,  1.d0,  0.0d0,
     /     -2.d0, 0.d0, 0.d0, 0.d0, 1.d0, -2.d0, 0.0d0,  1.d0,  0.0d0,
     /      2.d0, 0.d0, 0.d0, 0.d0, 1.d0,  2.d0, 0.0d0, -1.d0,  0.0d0,
     /      3.d0, 0.d0, 0.d0, 0.d0, 0.d0,  2.d0, 0.0d0,  0.d0,  0.0d0,
     /      1.d0, 1.d0, 2.d0, 0.d0, 2.d0,  2.d0, 0.0d0, -1.d0,  0.0d0,
     /      0.d0, 0.d0, 2.d0, 1.d0, 2.d0,  2.d0, 0.0d0, -1.d0,  0.0d0,
     /      1.d0, 0.d0, 0.d0, 2.d0, 1.d0, -1.d0, 0.0d0,  0.d0,  0.0d0,
     /      1.d0, 0.d0, 2.d0, 2.d0, 1.d0, -1.d0, 0.0d0,  1.d0,  0.0d0,
     /      1.d0, 1.d0, 0.d0,-2.d0, 1.d0, -1.d0, 0.0d0,  0.d0,  0.0d0,
     /      0.d0, 1.d0, 0.d0, 2.d0, 0.d0, -1.d0, 0.0d0,  0.d0,  0.0d0,
     /      0.d0, 1.d0, 2.d0,-2.d0, 0.d0, -1.d0, 0.0d0,  0.d0,  0.0d0,
     /      0.d0, 1.d0,-2.d0, 2.d0, 0.d0, -1.d0, 0.0d0,  0.d0,  0.0d0,
     /      1.d0, 0.d0,-2.d0, 2.d0, 0.d0, -1.d0, 0.0d0,  0.d0,  0.0d0,
     /      1.d0, 0.d0,-2.d0,-2.d0, 0.d0, -1.d0, 0.0d0,  0.d0,  0.0d0/
      DATA X5/ 1.d0,0.d0,2.d0,-2.d0,0.d0, -1.d0, 0.0d0,  0.d0,  0.0d0,
     /      1.d0, 0.d0, 0.d0,-4.d0, 0.d0, -1.d0, 0.0d0,  0.d0,  0.0d0,
     /      2.d0, 0.d0, 0.d0,-4.d0, 0.d0, -1.d0, 0.0d0,  0.d0,  0.0d0,
     /      0.d0, 0.d0, 2.d0, 4.d0, 2.d0, -1.d0, 0.0d0,  0.d0,  0.0d0,
     /      0.d0, 0.d0, 2.d0,-1.d0, 2.d0, -1.d0, 0.0d0,  0.d0,  0.0d0,
     /     -2.d0, 0.d0, 2.d0, 4.d0, 2.d0, -1.d0, 0.0d0,  1.d0,  0.0d0,
     /      2.d0, 0.d0, 2.d0, 2.d0, 2.d0, -1.d0, 0.0d0,  0.d0,  0.0d0,
     /      0.d0,-1.d0, 2.d0, 0.d0, 1.d0, -1.d0, 0.0d0,  0.d0,  0.0d0,
     /      0.d0, 0.d0,-2.d0, 0.d0, 1.d0, -1.d0, 0.0d0,  0.d0,  0.0d0,
     /      0.d0, 0.d0, 4.d0,-2.d0, 2.d0,  1.d0, 0.0d0,  0.d0,  0.0d0,
     /      0.d0, 1.d0, 0.d0, 0.d0, 2.d0,  1.d0, 0.0d0,  0.d0,  0.0d0,
     /      1.d0, 1.d0, 2.d0,-2.d0, 2.d0,  1.d0, 0.0d0, -1.d0,  0.0d0,
     /      3.d0, 0.d0, 2.d0,-2.d0, 2.d0,  1.d0, 0.0d0,  0.d0,  0.0d0,
     /     -2.d0, 0.d0, 2.d0, 2.d0, 2.d0,  1.d0, 0.0d0, -1.d0,  0.0d0,
     /     -1.d0, 0.d0, 0.d0, 0.d0, 2.d0,  1.d0, 0.0d0, -1.d0,  0.0d0,
     /      0.d0, 0.d0,-2.d0, 2.d0, 1.d0,  1.d0, 0.0d0,  0.d0,  0.0d0,
     /      0.d0, 1.d0, 2.d0, 0.d0, 1.d0,  1.d0, 0.0d0,  0.d0,  0.0d0,
     /     -1.d0, 0.d0, 4.d0, 0.d0, 2.d0,  1.d0, 0.0d0,  0.d0,  0.0d0,
     /      2.d0, 1.d0, 0.d0,-2.d0, 0.d0,  1.d0, 0.0d0,  0.d0,  0.0d0,
     /      2.d0, 0.d0, 0.d0, 2.d0, 0.d0,  1.d0, 0.0d0,  0.d0,  0.0d0/
      DATA X6/ 2.d0,0.d0,2.d0,-2.d0,1.d0,  1.d0, 0.0d0, -1.d0,  0.0d0,
     /      2.d0, 0.d0,-2.d0, 0.d0, 1.d0,  1.d0, 0.0d0,  0.d0,  0.0d0,
     /      1.d0,-1.d0, 0.d0,-2.d0, 0.d0,  1.d0, 0.0d0,  0.d0,  0.0d0,
     /     -1.d0, 0.d0, 0.d0, 1.d0, 1.d0,  1.d0, 0.0d0,  0.d0,  0.0d0,
     /     -1.d0,-1.d0, 0.d0, 2.d0, 1.d0,  1.d0, 0.0d0,  0.d0,  0.0d0,
     /      0.d0, 1.d0, 0.d0, 1.d0, 0.d0,  1.d0, 0.0d0,  0.d0,  0.0d0/

C-----------------------------------------------------------------------
C 2.5  COEFFICIENTS OF SERIES EXPANSION FOR BASIC ARGUMENTS
C
C     VALUES GIVEN BY IERS TECHNICAL NOTE NO. 3, 1989 IN ARCSECS
C-----------------------------------------------------------------------

      DATA B /
     A  485866.733D0, 1325.D0,  715922.633D0,  31.310D0,  0.064D0,
     B 1287099.804D0,   99.D0, 1292581.224D0,  -0.577D0, -0.012D0,
     C  335778.877D0, 1342.D0,  295263.137D0, -13.257D0,  0.011D0,
     D 1072261.307D0, 1236.D0, 1105601.328D0,  -6.891D0,  0.019D0,
     E  450160.280D0,   -5.D0, -482890.539D0,   7.455D0,  0.008D0 /

C-----------------------------------------------------------------------
C 2.6 OTHER CONSTANTS
C-----------------------------------------------------------------------

      SER = 1.296D6
      SOR = 1.D0/ROS

C-----------------------------------------------------------------------
C  ASSIGN THE NUTATION COEFICCIENTS IN A DO-LOOP WAY INSTEAD OF THE
C  FORMER EQUIVALENCE STATEMENT
C-----------------------------------------------------------------------

      DO i=1,9
         DO j=1,20
            XW (i,j)    = X1(i,j)
            XW (i,20+j) = X2(i,j)
            XW (i,40+j) = X3(i,j)
            XW (i,60+j) = X4(i,j)
            XW (i,80+j) = X5(i,j)
            IF (j.le.6) XW(i,100+j) = X6(i,j)
         ENDDO
      ENDDO
C-----------------------------------------------------------------------
C 3. JULIAN CENTURIES SINCE BASIC EPOCH
C-----------------------------------------------------------------------

      D50  = D50IN

C----------------------------------------------------------------------
C  ADDITION OF EXTRA SECS TO TRANSFORM UT IN ET OR TDT
C  IN 1984 : ET - UT = 54 SEC
C  FROM 1984 ON: TDT-UT = 32.184 + INTEGER NO. OF SECONDS
C  READ FROM POLE TABLES.
C-----------------------------------------------------------------------

      D50  = D50 + (32.184D0 + IDELT)/ 86400.D0


c      D50K = D50KIN / 36525.D0
      DCY  = D50    / 36525.D0
      DCY2 = DCY * DCY
      DCY3 = DCY * DCY2
      DCY4 = DCY2 * DCY2
      DCY5 = DCY2 * DCY3

C-----------------------------------------------------------------------
C 4. NO PROPER MOTION OF QUASARS
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C 5. COMPUTATION OF UNIT VECORS IN THE DIRECTION OF THE OBJECT
C    IN THE CATALOGUE SYSTEM
C-----------------------------------------------------------------------

      SID = DSIN (DELTA)
      COD = DCOS (DELTA)
      SIA = DSIN (ALPHA)
      COA = DCOS (ALPHA)

      RQU(1) = COD * COA
      RQU(2) = COD * SIA
      RQU(3) = SID

C-----------------------------------------------------------------------
C 6. COMPUTATION OF NEWCOMB PARAMETERS
C    (IAU-VALUES OF 1976; LIESKE 1977, P. 15)
C-----------------------------------------------------------------------

      ZETA0 = ((-1.39D-4 * D50K + 1.39656D0) * D50K + 2306.2181D0) *DCY
     1       + (-3.45D-4 * D50K + 0.30188D0) * DCY2
     2       + 0.017998D0 * DCY3

      ZET   = ((-1.39D-4 * D50K + 1.39656D0) * D50K + 2306.2181D0) *DCY
     1       + ( 0.66D-4 * D50K + 1.09468D0) * DCY2
     2       + 0.018203D0 * DCY3

      THETA = ((-2.17D-4 * D50K - 0.85330D0) * D50K + 2004.3109D0) *DCY
     1       + (-2.17D-4 * D50K - 0.42665D0) * DCY2
     2       - 0.041833D0 * DCY3

C-----------------------------------------------------------------------
C 6. COMPUTATION OF NEWCOMB PARAMETERS
C    (IAU-VALUES OF 2000; CAPITAINE et al, 2002, p.572)
C-----------------------------------------------------------------------
c
c      ZETA0 =  2306.0809506d0*DCY + 0.3019015d0*DCY2 + 0.0179663d0*DCY3
c     *      -  0.0000327d0*DCY4 - 0.0000002d0*DCY5 + 2.5976176d0
c
c      ZET =  2306.0803226d0*DCY + 1.0947790d0*DCY2 + 0.0182273d0*DCY3
c     *      +  0.0000470d0*DCY4 - 0.0000003d0*DCY5 - 2.5976176d0
c
c      THETA =  2004.1917476d0*DCY - 0.4269353d0*DCY2 - 0.0418251d0*DCY3
c     *      -  0.0000601d0*DCY4 - 0.0000001d0*DCY5


      RW(1) =  ZETA0 / ROS        !  Changed on 8.11.2005
      RW(2) =  THETA / ROS
      RW(3) =  ZET   / ROS        !  Changed on 8.11.2005

C-----------------------------------------------------------------------
C 7. SERIES EXPANSION FOR NUTATION
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C 7.1 COMPUTING NUMERICAL VALUES FOR THE BASIC ARGUMENTS A(1), ..., A(5)
C-----------------------------------------------------------------------

      DO I = 1, 5
         A(I) = ((B(5,I)*DCY+B(4,I))*DCY+B(3,I)+B(2,I)*SER)*DCY+B(1,I)
         A(I) = A(I)*SOR
         A(I) = DMOD(A(I),TWOPI)
      END DO

C-----------------------------------------------------------------------
C 7.2 EXCENTRICITY OF ECLIPTIC (EXPLANATORY SUPPLEMENT, P. 98)
C-----------------------------------------------------------------------

      A(6) = (- 9.445D-17 * D50 - 1.1513D-9) * D50 + 0.016709114d0

C-----------------------------------------------------------------------
C 7.3 MEAN OBLIQUITY OF ECLIPTIC (LIESKE, P. 15)
C-----------------------------------------------------------------------

      A(7) = ((0.001813D0 * D50K - 0.00059D0) * D50K - 46.8150D0) *D50K
     1     + ((0.005439D0 * D50K - 0.00117D0) * D50K - 46.8150D0) *DCY
     2     +  (0.005439D0 * D50K - 0.00059D0)                     *DCY2
     3     +   0.001813D0                                         *DCY3
     4     +  84381.448D0
      A(7) = A(7) / ROS         !  epsilon      !  Mean obliquity, radian

C-----------------------------------------------------------------------
C 7.3 MEAN OBLIQUITY OF ECLIPTIC (IAU-VALUES OF 2000; CAPITAINE et al, 2002, p.572)
C-----------------------------------------------------------------------
C

c      A(7) = 84381.448d0 - 46.84024d0*DCY -
c     * 0.00059d0*DCY2 + 0.001813*DCY3
c
c
c      A(7) = A(7) / ROS         !  epsilon      !  Mean obliquity


C-----------------------------------------------------------------------
C 7.4 INITIALIZATION OF ANGLES (NUTATION)
C-----------------------------------------------------------------------

      RW(4) =   A(7)      !        epsilon      !  Mean obliquity
      RW(5) =   0.D0
      RW(6) = - A(7)      !      - epsilon      !  Mean obliquity  with negative sign

C-----------------------------------------------------------------------
C 7.6 SERIES EXPANSION
C-----------------------------------------------------------------------

      if (nutmod.eq.1) then           !  IAU 1980 model
         DO I = 1, 106
            ARG = 0.D0
            DO J = 1, 5
               ARG = ARG + XW(J,I) * A(J)
            END DO
c
c  24.10.2005  Change to transform from the true system of the date to
c              the mean system of the date
c
            RW(5) = RW(5) + (XW(6,I)+XW(7,I)*DCY)*1.D-4 * DSIN(ARG)/ROS
            RW(4) = RW(4) + (XW(8,I)+XW(9,I)*DCY)*1.D-4 * DCOS(ARG)/ROS

c
         end do
      end if

c      print *, ' ut_h = ', 2400000.5d0+IDATJ+UTJ/TWOPI
c      print *, ' djdeph = ', D50 + 2451545.D0

      if (nutmod.eq.2) then               !  IAU2000a model (IERS 2003)


c
c  TDT to be used instead of UT here !
c
         	
c         ut = IDATJ+UTJ/TWOPI           !
c         ut_h = 2400000.5d0 + ut        !

         DJDEPH = D50 + 2451545.D0      !  Change on 1-May-2009

c
c  IAU2000a model replaces the old MHB2000 one
c
c         
         call iau2000a ( djdeph, dpsi_ls, deps_ls,
     .                    dpsi_plan, deps_plan,
ccc     .                    CX, CY,                   !  19-July-2019 (OT)
     *                    dpsi_fcn, deps_fcn,
     .                    dpsi_prec, deps_prec,
     .                    dpsi_tot , deps_tot )




         dp1 = dpsi_tot * 1.d-03   !  mas -> arcsec
         de1 = deps_tot * 1.d-03   !  mas -> arcsec


c
c  24.10.2005  Change to turansform from the true system of the date to
c              the mean system of the date
c
                                                              
         rw(5) = rw(5) + dp1/ros           !  true longitude, radian
         rw(4) = rw(4) + de1/ros           !  true obliquity, radian


c  *********** Part for alternative reduction for FCN,
c              used in April - July 2019, to be checked
c   commented on 19 July 2019 (OT)
c
c  FCN
c

ccc         denom = (RW(5)*DCOS(A(7)))**2 + 1.d0

         
ccc         CX = 1.d-3 * CX     !  muas -> mas
ccc        CY = 1.d-3 * CY     !  muas -> mas
                                                                   
ccc         dpsi_fcn = -(CX - RW(5)*DCOS(A(7)) * CY)/denom   ! mas
ccc         deps_fcn = -(CX*RW(5) * DCOS(A(7)) + CY)/denom   ! mas


ccc         rw(5) = rw(5) + 1.d-3*dpsi_fcn/ros            !  from mas to radian
ccc         rw(4) = rw(4) + 1.d-3*deps_fcn/ros            !  from mas to radian

******* end of temporarily used part (OT)


      end if
C-----------------------------------------------------------------------
C 7.65 GEODESIC NUTATION
C   FUKUSHIMA (1990) FROM IERS STANDARDS (1992)
C
C   CORRECTION TO THE NUTATION IN LONGITUDE TERM
C-----------------------------------------------------------------------


      if (nutmod.eq.1) then

         DPSIG = (-0.000153d0 * DSIN(A(2)) -
     *          0.000002D0 * DSIN (2.d0*A(2))) / ROS
         RW(5) = RW(5) + DPSIG                                  !   psi

      end if

C-----------------------------------------------------------------------
C 7.7 EQUATION OF EQUINOXES
C   STANDARD TREATMENT IS TO MULTIPLY THE NUTATION IN LONGITUDE BY THE
C     COSINE OF THE MEAN OBLIQUITY
C
C-----------------------------------------------------------------------


c
c    The true longitude * cosine of the mean obliquity
c
       EQE =   RW (5) * DCOS (A(7))  !  for the mean system of the date



c
c  24.10.2005  Change to transform from the true system of the date to
c              the mean system of the date
c
      OFFN(1) = RW(5)           !        psi           !  true longitude
      OFFN(2) = RW(4)           !        epsilon       !  true obliquity
      OFFN(3) = A(7)            !        epsilon       !  mean obliquity


      

      eps_0 = pi * (84381.448d0 / 3600.d0 * 180.d0)    ! e_0, rad  

      dx = dp1 * dsin(a(7)) + de1 * offn(1) * dcos(eps_0)
      dy = de1 - dp1 * offn(1) * dcos(eps_0) * dsin(a(7))
  
        


c
c  Matrix for nutation 24.10.05    N=R1(-eps)*R3(dpsi)*R1(eps+deps)
c

      sn(1,1)= dcos(rw(5))
      sn(2,1)= -dcos(a(7))*dsin(rw(5))
      sn(3,1)= -dsin(a(7))*dsin(rw(5))
      sn(1,2)= dcos(rw(4))*dsin(rw(5))
      sn(2,2)= dcos(a(7))*dcos(rw(4))*dcos(rw(5))+dsin(a(7))*dsin(rw(4))
      sn(3,2)= dsin(a(7))*dcos(rw(4))*dcos(rw(5))-dcos(a(7))*dsin(rw(4))
      sn(1,3)= dsin(rw(4))*dsin(rw(5))
      sn(2,3)= dcos(a(7))*dsin(rw(4))*dcos(rw(5))-dsin(a(7))*dcos(rw(4))
      sn(3,3)= dsin(a(7))*dsin(rw(4))*dcos(rw(5))+dcos(a(7))*dcos(rw(4))

c      print *, sn(2,3),  de1/ros
c      print *, sn(3,2), -de1/ros

c
c  Matrix for precession 24.10.05
c

      PREC(1,1)= dcos(rw(1))*dcos(rw(2))*dcos(rw(3))-
     *   dsin(rw(1))*dsin(rw(3))
      PREC(2,1)=-dsin(rw(1))*dcos(rw(2))*dcos(rw(3))-
     *   dcos(rw(1))*dsin(rw(3))
      PREC(3,1)= -dsin(rw(2))*dcos(rw(3))
      PREC(1,2)= dcos(rw(1))*dcos(rw(2))*dsin(rw(3))+
     *   dsin(rw(1))*dcos(rw(3))
      PREC(2,2)=-dsin(rw(1))*dcos(rw(2))*dsin(rw(3))+
     *   dcos(rw(1))*dcos(rw(3))
      PREC(3,2)= -dsin(rw(2))*dsin(rw(3))
      PREC(1,3)= dcos(rw(1))*dsin(rw(2))
      PREC(2,3)= -dsin(rw(1))*dsin(rw(2))
      PREC(3,3)= dcos(rw(2))

      DO I = 1, 3
         
          Do j=1,3

             PN(I,j) = 0.D0        !  P*N
          
             DO  K = 1, 3

                PN(i,j) =  PN(i,j) + prec(i,k) * sn(k,j)        !  P*N

             ENDDO

          end do

      END DO


cc      xx = rw(5) * dsin(a(7)) + rw(2) + (rw(3)+ rw(1))*de1/ros
cc      yy = rw(2)*rw(1) -de1/ros
cc      zz = - rw(5)*dcos(a(7)) - (rw(1) + rw(3))


c      print *, pn(1,1), 1.d0 - xx**2/2.d0 - zz**2/2.d0
c      print *, pn(2,1), zz
c      print *, pn(2,2), 1.d0 - yy**2/2.d0 - zz**2/2.d0
c      print *, pn(3,1), -xx + yy*zz
c      print *, pn(3,2), xx*zz + yy
c      print *, pn(1,2), xx*yy-zz
c      print *, pn(1,3), xx
c      print *, pn(2,3), -yy
c      print *, pn(3,3), 1.d0 - xx**2/2.d0 - yy**2/2.d0

c      print *, pn(1,1), 1.d0-2.d0*rw(1)**2 - rw(2)**2/2.d0 
c     * - rw(2)*rw(5) * a(7)
c      print *, pn(2,1), - 2.d0*rw(1) - rw(5)
c      print *, pn(2,2), 
c     *   1.d0 - 2.d0*rw(1)**2 - 2.d0*rw(1)*rw(5)
c      print *, pn(3,1), -rw(2) - rw(5)*a(7)
c      print *, pn(3,2), -rw(2)*rw(5)-de1/ros - rw(1)*rw(2)
c      print *, pn(1,2), rw(5) + 2.d0*rw(1) + de1*rw(2)/ros
c      print *, pn(1,3), rw(2) + rw(5)*a(7) + 2.d0*rw(1)*de1/ros
c      print *, pn(2,3), de1/ros  - rw(1)*rw(2) 
c      print *, pn(3,3), 1.d0 - rw(2)**2/2.d0 


C SET TDT TIME FOR COMPUTATION OF EPHEMERIS

      DJDEPH = D50 + 51544.5D0

c
c  -)  insert DE405/LE405 ephemerides for Sun, Moon, Earth, Jupiter
c
      if (m_eph.eq.1) then

         call DELE_405 (djdeph + 2400000.5d0,
     *    rrde, rrds, rrdm, rrdjup)

         do i=1,3
            dvelb(i)=rrde(i+3)          !  The Earth barycentric velocities
            dcorb(i)=rrde(i)            !  The Earth barycentric coordinates
            DCG(i)  = rrds(i)           !  The Sun geocentric coordinates
            DVG(i)  = rrds(i+3)         !  The Sun geocentric velocities
            XLUN(i) = rrdm(i)           !  The Moon geocentric coordinates
            Xjup(i,jmap) = rrdjup(i)    !  The Jupiter geocentric coordinates
         end do

      else

         CALL BARINT (DJDEPH,DCORB,DVELB)

C**************************************************************
C     COMPUTE THE GEOCENTRIC COORDINATES AND VELOCITIES OF THE SUN
C**************************************************************

         CALL SUNINT (DJDEPH,DCG,DVG)         !  Interpolation for velocities
                                               !  was added in September, 2001
c         print *, djdeph, dcg(1) 

C-----------------------------------------------------------------
C   COMPUTE THE GEOCENTRIC COORDINATES OF THE MOON IF
C   PEP EPHEMERIS ARE IN USE
C------------------------------------------------------------------

         CALL MOONINT (DJDEPH,XLUN)

         do i=1,3
            Xjup(i,jmap) = 0.d0     !  Jupiter is not calculated for this option
         end do

      end if

      DO  I = 1, 3
         DVELBN(I,JMAP) = DVELB(i)      !  The Earth barycentric velocities
         DCORBN(I,JMAP) = DCORB(i)      !  The Earth barycentric coordinates
      END DO


c
c     Calculation of the vectors s*P and Vearth*P (P  - precession matrix)
c     to be used in NUTPAR1 code
c

      DO I = 1, 3
         PR(I) = 0.D0        !  s*P
         VR(I) = 0.D0        !  Vearth*P
         XS(I) = 0.D0
         VS(I) = 0.D0
         XM(I) = 0.D0
         DO  J = 1, 3
            PR(I) =  PR(I) + RQU(J) * PREC(J,I)        !  s*P
            VR(I) =  VR(I) + DVELB(J) * PREC(J,I)      !  Vearth*P
            XS(I) =  XS(I) + DCG(J) * PREC(J,I)        !  Xsun*P
            VS(I) =  VS(I) + DVG(J) * PREC(J,I)        !  Vsun*P
            XM(I) =  XM(I) + XLUN(J) * PREC(J,I)       !  Xmoon*P
         ENDDO
      END DO

      DO I = 1, 3

         R2(I) = 0.D0
         VR2(I) = 0.D0
         DCORG(I,JMAP) = 0.d0           !  The Sun geocentric coordinates
         DVELG(I,JMAP) = 0.d0           !  The Sun geocentric velocities
         CMOON(I,JMAP) = 0.d0           !  The Moon geocentric coordinates

         DO  J = 1, 3

            R2(I) =  R2(I) + PR(J) * SN(J,I)       !  s*P*N        apparent positions
            VR2(I) = VR2(I) + VR(J) * SN(J,I)      !  Vearth*P*N   apparent velocities

            DCORG(I,JMAP) = DCORG (I,JMAP) + XS(J) * SN(J,I)  !  Xsun*P*N   apparent positions for tides
            DVELG(I,JMAP) = DVELG (I,JMAP) + VS(J) * SN(J,I)  !  Vsun*P*N
            CMOON(I,JMAP) = CMOON (I,JMAP) + XM(J) * SN(J,I)  !  Xmoon*P*N

         ENDDO

      END DO



C-----------------------------------------------------------------------
C 9. ANNUAL ABERRATION
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     ABERRATION (P. STUMPFF (1979), GL. 34 - 36)
C     ANNUAL ABERRATION MUST BE APPLIED JUST FOR VISIBLE POSITIONS
C     FOR ZENITH DISTANCE AND AZIMUTH
C-----------------------------------------------------------------------

c
c   Change on 03.05.2006; the coordinates corrected for P*N to be used here
c

      VELR    = - (R2(1) * VR2(1) + R2(2) * VR2(2) + R2(3) * VR2(3))
      VE2    = DSQRT(VR2(1)**2 + VR2(2)**2 + VR2(3)**2)

c
c   Change on 03.05.2006; the coordinates corrected for P*N to be used here
c
      DO I = 1, 3

c
c        Second term of the aberrational effect is updated on 08.02.2007
c

         d2(i) = 0.5d0 * (2.d0 * velr**2 * r2(i) - ve2**2 * r2(i) -
     *                    velr*vr2(i))/(c**2)

         R2(I) = R2(I) + (VR2(I) + R2(I) * VELR)/C + d2(i)

      END DO

C-----------------------------------------------------------------------
C 10. SPHERICAL POLAR COORDINATES
C-----------------------------------------------------------------------


      RAAP = DATAN2(R2(2),R2(1))
      IF (R2(2) .LT. 0) RAAP =  TWOPI + RAAP
      DEAP = DATAN2(R2(3),DSQRT(R2(1)*R2(1) + R2(2)*R2(2)))


      RETURN
      END




C************************************************************************
C
C
      SUBROUTINE NUTPAR1 (rqu,HSG,sn,PREC,PR,VR,OFFN,
     * DKDPSI,DKDEPS,DKDPSI_v, DKDEPS_v)
C
C-----------------------------------------------------------------------
C THIS SUBROUTINE COMPUTES THE PARTIAL DERIVATIVES OF THE SOURCE POSITION
C WITH RESPECT TO EPSILON AND PSI
C
C PROGRAMMER:  N. ZARRAOA, NOV 12, 1991
C
C  REVISION  1991 NOVEMBER 12   (N.Z.L.)
C  REVISION  2001 APRIL 17   (OT)  Remark made by Johannes Boehm
C  LAST REVISION  2002 June 19,   (OT)  HSG correction added
C
C-----------------------------------------------------------------------

      IMPLICIT DOUBLE PRECISION  (A-H,O-Z)
C
      DIMENSION PSI1(3),PSI2(3), PSI1_V(3),PSI2_V(3)
      DIMENSION DNDPSI(3,3),DNDEPS(3,3),SP(3,3),OFFN(3), rqu(3),
     * DSDPSI(3,3),DSDEPS(3,3), PR(3), VR(3)
      DIMENSION sn(3,3), PREC(3,3), S2(3,3), DKDPSI(3), DKDEPS(3)
      DIMENSION DKDPSI_v(3), DKDEPS_v(3)

      COMMON /PHYS/ C, FL, AH, AU, J20
      COMMON /MATH/ PI, TWOPI, ROG, ROS, ROZ, TOL

CC-----------------------------------------------------------------------
CC   COMPUTE THE PARTIAL DERIVATIVE OF THE NUTATION MATRIX WITH
CC  RESPECT TO EPSILON AND PSI
CC-----------------------------------------------------------------------
C

      DNDPSI(1,1) = -DSIN(OFFN(1))
      DNDPSI(1,2) =  DCOS(OFFN(1)) * DCOS(OFFN(2))
      DNDPSI(1,3) =  DCOS(OFFN(1)) * DSIN(OFFN(2))
      DNDPSI(2,1) = -DCOS(OFFN(1)) * DCOS(OFFN(3))
      DNDPSI(2,2) = -DSIN(OFFN(1)) * DCOS(OFFN(2)) * DCOS(OFFN(3))
      DNDPSI(2,3) = -DSIN(OFFN(1)) * DSIN(OFFN(2)) * DCOS(OFFN(3))
      DNDPSI(3,1) = -DCOS(OFFN(1)) * DSIN(OFFN(3))
c
c  OT 1-May-2009
c
c      DNDPSI(3,2) = -DSIN(OFFN(1)) * DSIN(OFFN(2)) * DCOS(OFFN(3))
c
c
      DNDPSI(3,2) = -DSIN(OFFN(1)) * DCOS(OFFN(2)) * DSIN(OFFN(3))   !  1-05-2009
      DNDPSI(3,3) = -DSIN(OFFN(1)) * DSIN(OFFN(2)) * DSIN(OFFN(3))

      DNDEPS(1,1) =  0.D0
      DNDEPS(1,2) = -DSIN(OFFN(1)) * DSIN(OFFN(2))
      DNDEPS(1,3) =  DSIN(OFFN(1)) * DCOS(OFFN(2))
      DNDEPS(2,1) = 0.D0
      DNDEPS(2,2) = -DCOS(OFFN(1)) * DCOS(OFFN(3)) * DSIN(OFFN(2)) +
     *               DSIN(OFFN(3)) * DCOS(OFFN(2))
      DNDEPS(2,3) =  DCOS(OFFN(1)) * DCOS(OFFN(3)) * DCOS(OFFN(2)) +
     *               DSIN(OFFN(3)) * DSIN(OFFN(2))
      DNDEPS(3,1) = 0.D0
      DNDEPS(3,2) = -DCOS(OFFN(1)) * DSIN(OFFN(3)) * DSIN(OFFN(2)) -
     *               DCOS(OFFN(3)) * DCOS(OFFN(2))
      DNDEPS(3,3) =  DCOS(OFFN(1)) * DSIN(OFFN(3)) * DCOS(OFFN(2)) -
     *               DCOS(OFFN(3)) * DSIN(OFFN(2))

CC-----------------------------------------------------------------------
CC   COMPUTE THE MATRIX FOR THE DIURNAL SPIN AND ITS DERIVATIVES WITH
CC  RESPECT TO EPSILON AND PSI
CC-----------------------------------------------------------------------

      SP(1,1) =  DCOS (HSG)
      SP(1,2) = -DSIN (HSG)
      SP(1,3) = 0.D0
      SP(2,1) =  DSIN (HSG)
      SP(2,2) =  DCOS (HSG)
      SP(2,3) = 0.D0
      SP(3,1) = 0.D0
      SP(3,2) = 0.D0
      SP(3,3) = 1.D0

      S2(1,1) = -DSIN (HSG) * DCOS(OFFN(3))
      S2(1,2) = -DCOS (HSG) * DCOS(OFFN(3))
      S2(1,3) = 0.D0
      S2(2,1) = DCOS (HSG) * DCOS(OFFN(3))
      S2(2,2) = -DSIN (HSG) * DCOS(OFFN(3))
      S2(2,3) = 0.D0
      S2(3,1) = 0.D0
      S2(3,2) = 0.D0
      S2(3,3) = 0.D0


C-----------------------------------------------------------------------
C  COMPUTE THE PARTIAL DERIVATIVES BY MATRIX MULTIPLICATION:
C    DKDPSI0 = s*P*DNDPSI*R3(-H) + s*P*N*dR3(-H)/dH*cos(eta)
C    DKDEPS0 = s*DNDEPS*R3(-H)
C-----------------------------------------------------------------------

      CALL ABCX (PR,1,3,0,'N',DNDPSI,3,3,0,'N',SP,3,3,0,'N',
     *           PSI1,NRX,NCX,0,IER)

      CALL ABCX (PR,1,3,0,'N',sn,3,3,0,'N',S2,3,3,0,'N',
     *           PSI2,NRX,NCX,0,IER)

      CALL ABCX (PR,1,3,0,'N',DNDEPS,3,3,0,'N',SP,3,3,0,'N',
     *           DKDEPS,NRX,NCX,0,IER)


C-----------------------------------------------------------------------
C  COMPUTE THE PARTIAL DERIVATIVES BY MATRIX MULTIPLICATION:
C    DKDPSI0 = Ve*P*DNDPSI*R3(-H) + Ve*P*N*dR3(-H)/dH*cos(eta)
C    DKDEPS0 = Ve*P*DNDEPS*R3(-H)
C
C  Aberration term at the nutation partials!
C-----------------------------------------------------------------------


      CALL ABCX (VR,1,3,0,'N',DNDPSI,3,3,0,'N',SP,3,3,0,'N',
     *           PSI1_v,NRX,NCX,0,IER)

      CALL ABCX (VR,1,3,0,'N',sn,3,3,0,'N',S2,3,3,0,'N',
     *           PSI2_v,NRX,NCX,0,IER)

      CALL ABCX (VR,1,3,0,'N',DNDEPS,3,3,0,'N',SP,3,3,0,'N',
     *           DKDEPS_v,NRX,NCX,0,IER)


      DO I=1,3
         DKDPSI(I)   =  PSI1(I) + PSI2(I)
         DKDPSI_v(I) = (PSI1_v(I) + PSI2_v(I))/c
         DKDEPS_v(I) =  DKDEPS_v(I)/c
      ENDDO

      RETURN
      END

************************************************************************
C    END OF THE NUTATION PARTIALS PART
************************************************************************

      SUBROUTINE BARINT (DJDEPH,DCORB,DVELB)
C
C   THIS ROUTINE PERFORMS INTERPOLATION FROM A GIVEN TABLE UP TO
C  A GIVEN PRECISION. IT IS USED TO INTERPOLATE BARICENTRIC
C  COORDINATES AND VELOCITIES OF EARTH. IT ALSO COMPUTES THE
C  BARICENTRIC ACCELERATION OF EARTH.
C
C  REVISION NESTOR ZARRAOA, OCTOBER 1990
C  REVISION  1993 NOV. 17 (NZL) INCREASED DIMENSIONS
C  LAST REVISION  2001 MARCH 11 (OT)
C  LAST REVISION  2001 SEPTEMBER 24 (OT), Interpolation LAGINT for
C                 coordinates, velocities and accelerations
C
C*****************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION DCORB(3),DVELB(3)
      COMMON /BARIC/ XBRC(30),YBRC(30),ZBRC(30),XBRV(30),YBRV(30),
     1 ZBRV(30), BCTIM(30),IP1
C
C   GET INTERPOLATED VALUES FOR COORDINATES
C

      CALL LAGINT (BCTIM, XBRC, IP1, DJDEPH, DCORB(1) )
      CALL LAGINT (BCTIM, YBRC, IP1, DJDEPH, DCORB(2) )
      CALL LAGINT (BCTIM, ZBRC, IP1, DJDEPH, DCORB(3) )

C
C   GET INTERPOLATED VALUES FOR VELOCITIES (ON DATE)
C

      CALL LAGINT (BCTIM, XBRV, IP1, DJDEPH, DVELB(1) )
      CALL LAGINT (BCTIM, YBRV, IP1, DJDEPH, DVELB(2) )
      CALL LAGINT (BCTIM, ZBRV, IP1, DJDEPH, DVELB(3) )


      RETURN
      END

***********************************************************************

      SUBROUTINE SUNINT (DJDEPH,DCORG,DVELG)

C----------------------------------------------------------------------
C
C   THIS ROUTINE PERFORMS INTERPOLATION FROM A GIVEN TABLE UP TO
C  A GIVEN PRECISION. IT IS USED TO INTERPOLATE GEOCENTRIC SUN
C  COORDINATES AND VELOCITIES.
C
C  REVISION NESTOR ZARRAOA, 4 DECEMBER 1990
C  REVISION  1993 NOV. 17  (NZL)  INCREASED DIMENSIONS
C  REVISION  2001 MARCH 11 (OT)
C  LAST REVISION  2001 SEPTEMBER 24 (OT), Interpolation LAGINT for
C                 coordinates and velocities, also changes in /GEOC/
C
C-----------------------------------------------------------------------

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION DCORG(3),DVELG(3)
      COMMON /GEOC/ XG(30),YG(30),ZG(30),VXG(30),VYG(30),VZG(30),
     * BTIM(30),IP2

      DO I=1,3
         DCORG(I)=0.D0
         DVELG(I)=0.d0
      END DO

C
C   GET INTERPOLATED VALUES FOR COORDINATES
C

      CALL LAGINT (BTIM, XG, IP2, DJDEPH, DCORG(1) )
      CALL LAGINT (BTIM, YG, IP2, DJDEPH, DCORG(2) )
      CALL LAGINT (BTIM, ZG, IP2, DJDEPH, DCORG(3) )

C
C   GET INTERPOLATED VALUES FOR VELOCITIES
C
      CALL LAGINT (BTIM, VXG, IP2, DJDEPH, DVELG(1) )
      CALL LAGINT (BTIM, VYG, IP2, DJDEPH, DVELG(2) )
      CALL LAGINT (BTIM, VZG, IP2, DJDEPH, DVELG(3) )


      RETURN
      END

**********************************************************************

      SUBROUTINE MOONINT (DJDEPH,XLUN)

C------------------------------------------------------------------
C
C   THIS ROUTINE PERFORMS INTERPOLATION FROM A GIVEN TABLE UP TO
C  A GIVEN PRECISION.
C
C  REVISION NESTOR ZARRAOA, 4 DECEMBER 1990
C  REVISION  1993 NOV. 17  (NZL)  INCREASED DIMENSIONS
C  REVISION  2001 MARCH 11 (OT)
C  LAST REVISION  2001 SEPTEMBER 24 (OT), Interpolation LAGINT for
C                 coordinates.
C
C------------------------------------------------------------------

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XLUN(3)
      COMMON /MOON/ TMOON(30),RAM(6,30),DECM(6,30),PARM(6,30),
     *              XMON(30), YMON(30), ZMON(30), IPM

      DO I=1,3
         XLUN(I)=0.D0
      END DO

C
C   GET INTERPOLATED VALUES FOR COORDINATES
C

      CALL LAGINT (TMOON, XMON, IPM, DJDEPH, XLUN(1) )
      CALL LAGINT (TMOON, YMON, IPM, DJDEPH, XLUN(2) )
      CALL LAGINT (TMOON, ZMON, IPM, DJDEPH, XLUN(3) )

C   THAT'S ALL FOLKS

      RETURN
      END



      SUBROUTINE PARAM
C
C   THIS ROUTINE INITIALIZES:
C
C    - THE MAIN PHYSICAL PARAMETERS USED IN OCCAM PACKAGE
C       . C, VELOCITY OF LIGHT
C       . FL, EARTH'S FLATTENING FACTOR
C       . AH, EARTH'S EQUATORIAL RADIUS
C       . AU, ASTRONOMICAL UNIT
C       . GMS, GRAVITATIONAL CONSTANT FOR THE SUN (G*M)
C       . GMM, GRAVITATIONAL CONSTANT FOR THE MOON (G*M)
C       . GME, GRAVITATIONAL CONSTANT FOR THE EARTH (G*M)
C       . mu, relationship Mass of Moon/Mass of Earth
C
C   All Love numbers for elastic model !!!
C
C       . H02, LOVE NUMBER h0                             0.6026
C       . L02, LOVE NUMBER l0                             0.0831
C       . H22, LOVE NUMBER h2  (latitude dependence)     -0.0006
C       . L22, LOVE NUMBER l2  (latitude dependence)      0.0002
C       . H3,  LOVE NUMBER h3                             0.292
C       . L3,  LOVE NUMBER l3                             0.015
C       . HI_1, LOVE NUMBER hi for diurnal tides         -0.0025
C       . LI_1, LOVE NUMBER li for diurnal tides         -0.0007
C       . HI_2, LOVE NUMBER hi for semidiurnal tides     -0.0022
C       . LI_2, LOVE NUMBER li for semidiurnal tides     -0.0007
C       . L1_1, LOVE NUMBER l1 for diurnal tides          0.0012
C       . L1_2, LOVE NUMBER l1 for semidiurnal tides      0.0024
C
C    - MATHEMATICAL PARAMETERS
C       . PI (3.14159...)
C       . TWOPI, 2*PI
C
C    - CONVERSION FACTORS
C       . ROG, FROM RADIANS TO DEGREES
C       . ROS, FROM RADIANS TO ARC SECONDS
C       . ROZ, FROM RADIANS TO TIME SECONDS
C
C    - OTHER CONSTANTS
C       . TOL, TOLERANCE (O VALUE FOR COMPUTER CALCULUS)
C       . J20, MODIFIED JULIAN DATE FOR DAY JD 2000.0
C
C   WHERE NOT ESPECIFICLY STATED, IERS RECOMMENDED VALUES ARE USED
C
C   PROGRAMMER: N. ZARRAOA
C
C        REVISION  1990 NOVEMBER     (NZL). IERS-89 STANDARDS
C        REVISION  1991 NOVEMBER     (NZL). TESTS WITH CALC
C        REVISION  1993 FEBRUARY 11  (NZL). IERS-92 STANDARDS.
C        REVISION  1996 APRIL    14  (OT).  IERS-96 STANDARDS.
C        REVISION  2001 FEBRUARY 07  (OT).  IERS-2000 STANDARDS.
C   LAST REVISION  2005 DECEMBER 22  (OT).  IERS-2003 STANDARDS.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION L02, L22, L3, LI_1, LI_2, L1_1, L1_2, mu

      COMMON /PHYS/ C, FL1, AH, AU, J20
      COMMON /MATH/ PI, TWOPI, ROG, ROS, ROZ, TOL
      COMMON /GRAV/ GMS, GMM, GME, GM_Mer, GM_Ven, GM_Mar, GM_Jup,
     * GM_Sat, GM_Ura, GM_Nep
      COMMON /LOVE/ H02, L02, H22, L22, H3, L3, HI_1, LI_1, HI_2, LI_2,
     *              L1_1, L1_2

C   VALUES FROM IERS (2003)

       PI = 3.1415926535897932D0
       C  = 2.99792458D8                !  m/sec       !  06-12-2012
       FL1 = 298.25642D0                !  Flattening factor
       TOL  = 1.D-14
       J20 = 51544
       AH  = 6378136.6D0                !  Earth radius
       CINV = 4.990047838061D+2         !  Astronomical Unit in seconds


       GMS = 0.13271242076D21           !  G*Mass of Sun
       mu  = 0.0123000383d0
       GM_Mer = GMS/6023600.d0          !  G*Mass of Mercury
       GM_Ven = GMS/408523.71d0         !  G*Mass of Venus
       GME = 3.986004418D14             !  G*Mass of Earth
       GM_Mar = GMS/3098708.d0          !  G*Mass of Mars
       GM_Jup = GMS/1047.3486d0         !  G*Mass of Jupiter
       GM_Sat = GMS/3497.898d0          !  G*Mass of Saturn
       GM_Ura = GMS/22902.98d0          !  G*Mass of Uranus
       GM_Nep = GMS/19412.24d0          !  G*Mass of Neptune


       H02 = 0.6026D0
       L02 = 0.0831D0
       H22 = -0.0006D0
       L22 = 0.0002D0
       HI_1 = -0.0025D0
       LI_1 = -0.0007D0
       H3  = 0.292D0
       L3  = 0.015D0
       HI_2 = -0.0022D0
       LI_2 = -0.0007D0
       L1_1 = 0.0012D0
       L1_2 = 0.0024D0

      GMM = mu * GME

      TWOPI = 2.D0 * PI

c    1296000 - arcseconds in 360 degrees

      ROG   = 360.D0  / TWOPI
      ROS   = 1296.D3 / TWOPI !  arcseconds in radian
      ROZ   = 8.64D4  / TWOPI !  seconds in radian

      AU = CINV * C           !  Astronomical unit in meters

      RETURN
      END

*******************************************************************

      SUBROUTINE ARC (U, IT, IMT, AST, NAT)
C------------------------------------------------------------------
C   SUBROUTINE ARC TRANSFORMS ANGLES IN RADIANS TO DEGREES (OR HOURS)
C  MINUTES AND SECONDS
C
C   INPUT VARIABLES
C     - U.- ANGLE IN RADIANS
C     - NAT.- 0 = CONVERSION TO DEGREES
C             1 = CONVERSION TO HOURS
C
C   OUTPUT VARIABLES
C     - IT.- INTEGER NUMBER OF DEGREES (HOURS)
C     - IMS.- INTEGER NUMBER OF MINUTES
C     - AST.- SECONDS
C
C   EXTERNAL
C     - COMMON H.- PHYSICAL CONSTANTS AND CONVERSION PARAMETERS
C
C  LATEST REVISION  1990 NOVEMBER 19  (N.Z.L.)
C-------------------------------------------------------------------

      IMPLICIT DOUBLE PRECISION  (A-H,O-Z)

      COMMON /MATH/ PI, TWOPI, ROG, ROS, ROZ, TOL

      IF (NAT .EQ. 0)  THEN
         T = U * ROG
      ELSE
         T = U * 24.D0 / TWOPI
      ENDIF
C
      IT  = T
      AT  = IT
      T   = T - AT
      T   = T * 60.D0
      T   = DABS(T)
      IMT = T
      AMT = IMT
      T   = T - AMT
      AST = T * 60.D0
C
      RETURN
      END

*******************************************************************

      SUBROUTINE TRANSF (X, Y, Z, B, L, ELH, R, IMODE)
C-----------------------------------------------------------------------
C TRANSFORMATION FROM GEOCENTRIC X, Y, Z TO ELLIPSOIDAL COORDINATES
C B, L, H  AND  THE OTHER WAY AROUND
C
C PROGRAMMERS: J. CAMPBELL, H. SCHUH
C
C   INPUT VARIABLES
C     - AH.- MAJOR SEMIAXIS OF ELLIPSOID
C     - FL1.- EARTH'S FLATTENING FACTOR
C     - IMODE .- IF 0 computes latitude, longitude and height from x,y,z
C              - IF 1 computes x, y, z from latitude, longitude and height
C
C   INPUT/OUTPUT VARIABLES
C     - X,Y,Z.- GEOCENTRIC CARTESIAN COORDINATES OF STATION
C     - B.- ELLIPSOIDAL LATITUDE
C     - L.- ELLIPSOIDAL LONGITUDE (REAL*8)
C     - H.- ELLIPSOIDAL HEIGHT
C     - R.- GEOCENTRIC DISTANCE
C
C   EXTERNAL
C     - COMMON PHYS.- PHYSICAL CONSTANTS
C     - COMMON MATH.- CONVERSION PARAMETERS
C
C        REVISION  1990 NOVEMBER 28  (N.Z.L.)
C        REVISION  1993 FEBRUARY 11  (NZL). REMOVE UNNECESSARY VARIABLES
C   LAST REVISION  1993 DECEMBER  9  (NZL). TWO-WAYS TRANSFORMATION
C-----------------------------------------------------------------------

      IMPLICIT DOUBLE PRECISION  (A-H,O-Z)
      DOUBLE PRECISION  L, LS, N

      COMMON /PHYS/ C, FL1, AH, AU, J20
      COMMON /MATH/ PI, TWOPI, ROG, ROS, ROZ, TOL

      FL  = 1.D0 / FL1
      BH  = AH * (1.D0 - FL)
      AH2 = AH * AH
      BH2 = BH * BH
      CC  = AH2 / BH
      EQ  = 2.D0 * FL - FL * FL
      ES2 = (AH2 - BH2) / BH2

C================================================================
C PART 1.  TRANSFORM FROM X, Y, Z  TO GEOD. LATIT., LONGIT., HEIGHT
C================================================================

      IF (IMODE.EQ.0) THEN

         XY=X*X+Y*Y

         if (xy.lt.0.d0) then
            print *, ' xy in transf less than zero'
            print *, ' x = ', ' y = ', ' xy = ', xy
         end if

         PP  = DSQRT(XY)

C-----------------------------------------------------------------------
C APPROXIMATE VALUE FOR N
C-----------------------------------------------------------------------

         N   = AH
         K   = 0
         ELH = 0.D0
C
    4 CONTINUE
C
         K = K + 1
C
         PE  = 1.D0 - EQ * (N / (N + ELH))
         TAB = Z / (PP * PE)
         B   = DATAN(TAB)
         CB  = DCOS(B)
         CB2 = CB * CB
         csb=1.D0 + ES2 * CB2
         if (csb.lt.0.d0) print *, ' csb in transf less than zero'

         N   = CC / DSQRT(CSB)
         ELH = (PP / CB) - N

C-----------------------------------------------------------------------
C 5 ITERATIONS
C-----------------------------------------------------------------------

         IF (K .NE. 4)   GOTO 4


         L = DATAN2(Y,X)
         IF (Y .LT. 0)  L = L + TWOPI
         R = DSQRT(X**2 + Y**2 + Z**2)

         CALL ARC (B, IB1, IB2, B3, 0)
         CALL ARC (L, IL1, IL2, LS, 0)

C=================================================================
C PART 2.  TRANSFORM FROM GEOD. LATIT., LONGIT., HEIGHT TO X, Y, Z
C=================================================================

      ELSE IF (IMODE.EQ.1) THEN

         N = AH2 / DSQRT (AH2 * DCOS (B)**2 + BH2 * DSIN(B)**2)

         X = (N + ELH) * DCOS (B) * DCOS (L)
         Y = (N + ELH) * DCOS (B) * DSIN (L)
         Z = (N*BH2/AH2 + ELH) * DSIN (B)

      ENDIF

      RETURN
      END

C *******************************************************************
C     SUBRUTINA JULDAT
C *******************************************************************
C
C USAGE : call juldat (iy,id,ih,im,s,idj,utj,ipaso)
C
C PROGRAMMER: N. ZARRAOA
C
C ARGUMENTS: IY    - Year
C            ID    - Day of year
C            IH    - Hour
C            IM    - Minutes
C            S     - Seconds
C            IDJ   - Modified Julian Day (integer)
C            UTJ   - Fraction of Julian Day (rad)
C            IPASO - If ipaso.eq.0  then   ydhms --------> Julian date
C                    if ipaso.ne.0  then  julian date ---> ydhms
C
C REQUIRED ROUTINES : None
C
C LASTEST REVISION : September 20, 1989
C LAST REVISION: 12, January, 2001 (Oleg Titov)  Problem 2001 has been fixed
C
C***************************************************************************

      SUBROUTINE JULDAT (IY,ID,IH,IM,S,IDJ,UTJ,IPASO)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      COMMON /PHYS/ C, FL1, AH, AU, J20
      COMMON /MATH/ PI, TWOPI, ROG, ROS, ROZ, TOL

C
C    INITIALIZING  ( J20 = MOD. JULIAN DATE FOR 2000 JAN. 01.0 )
C
C   PASS FROM YEAR + DOY TO JULIAN DATE
C

      IF (IPASO.EQ.0) THEN

        IF (IY.LE.2000) THEN
           IDIF=(IY-2000)/4
           IDIF=IDIF+(IY-2000)*365
           JDYEAR=J20+IDIF
        ELSE
           IDIF=(IY-2001)/4 + 1   !   It is following from 366 days in 2000
           IDIF=IDIF+(IY-2000)*365
           JDYEAR=J20+IDIF
        END IF

C   GET TIME IN RADIANS

        UTJ = (dble(IH)/24.D0+dble(IM)/1440.D0+S/86400.D0)*TWOPI
        IDJ=JDYEAR+(ID-1)

      ENDIF

C   PASS FROM JULIAN DATE TO YEAR + DOY

      IF (IPASO.NE.0) THEN

        IC0 = IDJ-J20
        IC=IC0/1461

C   CHANGE THE NUMBER OF 4 YEARS CYCLES IF DATE BEFORE 2000

        IF (IC0.LT.0) IC=IC-1

        IR=IC0-IC*1461
        IY = 2000+IC*4

        IF (IR.LT.366) ID = IR+1

        IF (IR.GE.366) THEN
            IR = IR-366
            IY1 = IR/365
            IY = IY + 1 + IY1
            ID = IR - IY1*365 + 1
        ENDIF

        HORAS=(UTJ/TWOPI)*24.d0
        IH=HORAS+0.0001d0
        PMIN=(HORAS-IH)*60.d0
        IM=PMIN+0.0001d0
        S=(PMIN-IM)*60.d0

      ENDIF

      RETURN
      END

C *******************************************************************
C     FUNCTION DOTPR
C *******************************************************************
C
C USAGE : DOTPR(A,B)
C
C DOT PRODUCT OF 3-DIM VECTORS A AND B
C

      DOUBLE PRECISION FUNCTION DOTPR(A,B)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(3),B(3)
      DOTPR = A(1)*B(1) + A(2)*B(2) + A(3)*B(3)
      RETURN
      END


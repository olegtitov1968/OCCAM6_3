      SUBROUTINE CPOSIT (is,ISTCAT,HSG,
     *    R_after_PN, V_after_PN, XLUN_after_PN,
     *    rrot, CVECT, pres)

C
C    SUBROUTINE CPOSIT
C
C   THIS SUBROUTINE COMPUTE CORRECTIONS DUE TO DIFFERENT FENOMENA THAT
C  AFFECTS VLBI OBSERVABLES OR MODELING, AS EARTH DEFORMATIONS
C
C  REVISION  1990 NOVEMBER 19   (N.Z.L.)
C  REVISION  1993 MARCH 3 (N.Z.L.) INCLUSION OF IERS-92 MODELS
C  REVISION  1993 DEC. 9 (NZL) PETROV COMMENTS
C  REVISION  1997 MAR. 30 (OT)  Atmospheric loading
C  REVISION  1997 APR. 09 (OT)  Mapping function of A.Niell
C  REVISION  2000 JULY 14 (VT) Mattews model for the solid tides,
C                              (OT) Secular polar motion for polar tide
C  REVISION  2000 DEC 24  (VT) Schernecks Atmo loading,
C                                   refined thermal deformation
C  REVISION 2001 January 24  (OT)  Gravity parameters removed
C  REVISION 2001 MARCH, 24  (OT)  COMMON /BAR1/ instead of /BAR/
C  REVISION 2001 APRIL, O6 (OT) Sub. atmosph is not used more
c  revision 2002 July 29 jboehm, new discrete models are included
C  REVISION 2006 MARCH, 24 (OT) Correction for OCCAM 6.2
C  REVISION 2006 SEPTEMBER, 04 (OT) Common blocks /BAR1/, /MOON/, /CORR/
C                                        deleted
C  REVISION 2010 June, 04 (OT) Martix R4 added
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
C************************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'OCCAM_N.FI'   ! To call 'nstations'

      LOGICAL OPTTD, OPTOC,OPTAT,OPTPT,OPTLOV,NOOCL(nstations)
c+ jboehm
      logical optsnow, optsoil, optoc1
c- jboehm
      DOUBLE PRECISION NSOCL
      DIMENSION CTIDE(3),COCE(3),CATM(3),CPOL(3),CVECT(3),rrot(3,3)

      dimension V_after_PN(3), R_after_PN(3), XLUN_after_PN(3)

      character*8 sta
      character*(8*nstations) atmsta
      dimension atmmjd(nstations,30)
      dimension atmx(nstations,30),atmy(nstations,30),atmz(nstations,30)
c+ jboehm
      character*(8*nstations) ocesta
      dimension ocemjd(nstations,30)
      dimension ocex(nstations,30),ocey(nstations,30),ocez(nstations,30)
c
      character*(8*nstations) snosta
      dimension snomjd(nstations,30)
      dimension snox(nstations,30),snoy(nstations,30),snoz(nstations,30)
c
      character*80 soista
      dimension soimjd(nstations,30)
      dimension soix(nstations,30),soiy(nstations,30),soiz(nstations,30)
c+ jboehm
      dimension csno(3),csoi(3),coce1(3)
c- jboehm
      COMMON /COCPOS/ GLONG,PHI,ELHGT, RX, RY, RZ
      COMMON /DATEb / IDATJ,UTJ,XJ,YJ,IDELT
      COMMON /OPT/ OPTTD,OPTOC,OPTAT,OPTPT,OPTLOV
c+ jboehm
     .             ,optsnow,optsoil,optoc1
c- jboehm
      COMMON /OCLOAD/ NOOCL,ROCL(nstations,2,noclotds),
     * NSOCL(nstations,2,noclotds), EWOCL(nstations,2,noclotds)
      COMMON /PHYS/ C, FL, AH, AU, J20
      COMMON /MATH/ PI, TWOPI, ROG, ROS, ROZ, TOL
      COMMON /atscher/ istatm,atmsta,atmmjd,atmx,atmy,atmz
c+ jboehm
      COMMON /ocean/   istoce,ocesta,ocemjd,ocex,ocey,ocez
      COMMON /snow/    istsno,snosta,snomjd,snox,snoy,snoz
      COMMON /soil/    istsoi,soista,soimjd,soix,soiy,soiz
c- jboehm



C  GET TIME FROM 1900.0 IN JULIAN CENTURIES (T00)
C  CHECK WITH THE MODEL IF IT IS JD MODIFIED OR NOT

      DO I=1,3
         CTIDE(I)=0.D0
         COCE(I)=0.D0
         CATM(I)=0.D0
         CPOL(I)=0.D0
c+ jboehm
         csno(i)=0.d0
         csoi(i)=0.d0
         coce1(i)=0.d0
c- jboehm
         CVECT(I)=0.D0
      end do

      DATJ = IDATJ + UTJ/TWOPI
      T00 = (DATJ - 15019.5D0) / 36525.D0

C
C  CALL CORRECTIONS DUE TO EARTH DEFORMATIONS. THESE CORRECTIONS
C  MEAN CHANGES IN THE COORDINATES OF THE STATIONS APPLIED TO THE
C  GEOMETRIC MODEL FOR THE THEORETICAL DELAY
C
C
C     CORRECTION DUE TO EARTH TIDES.
C     The matthews-model is introduced


      CALL matthew (datj,HSG,RX,RY,RZ,
     *    R_after_PN, V_after_PN, XLUN_after_PN,
     *    rrot,CTIDE(1),CTIDE(2),CTIDE(3))

C
C     CORRECTION DUE TO OCEAN LOADING
C
      LOADFC=1
      IF (.NOT.NOOCL(ISTCAT)) LOADFC=0

C    DAYF WILL DEPEND ON THE FORMAT FOR THE TIME IN THE ATM. LOADING
C  DATA. THIS IS  JUST A PROVISIONAL VALUE

      DAYF = IDATJ + UTJ /TWOPI
      CALL JULDAT (IY,ID,IH,IM,S,IDATJ,UTJ,1)
      IFAC = 365
      IF (MOD(IY,4).EQ.0) IFAC=366
      YEAR = dble(IY)+(dble(ID)+dble(IH)/23.93447d0+dble(IM)/1440.d0+
     * dble(S)/86400.d0)/(dble(IFAC)+0.2422d0)


      IF (OPTOC.AND.LOADFC.EQ.0) THEN
         TIMDT = IDATJ + UTJ/TWOPI + (32.184+IDELT)/86400.d0
         CALL COCLD (TIMDT,RX,RY,RZ,ISTCAT,COCE(1),COCE(2),COCE(3))
      ENDIF



C
C    CORRECTION DUE TO ATMOSPHERICAL LOADING
C     atmospheric loading scherneck

      IF (OPTAT) THEN

        
         read(20,rec=istcat) jp,sta
         call atmo_sch(datj,istatm,atmsta,atmmjd,
     *                atmx,atmy,atmz,
     *                sta,catm(1),catm(2),catm(3))

c         print *, catm(1), catm(2), catm(3)

      ENDIF

c+ jboehm
C
C    CORRECTION DUE TO OCEAN LOADING (DISCRETE)
C

      IF (OPToc1) THEN
         read(20,rec=istcat) jp,sta
         call atmo_sch(datj,istoce,ocesta,ocemjd,
     *                ocex,ocey,ocez,
     *                sta,coce1(1),coce1(2),coce1(3))
      ENDIF

C
C    CORRECTION DUE TO SNOW LOADING
C

      IF (OPTsnow) THEN
         read(20,rec=istcat) jp,sta
         call atmo_sch(datj,istsno,snosta,snomjd,
     *                snox,snoy,snoz,
     *                sta,csno(1),csno(2),csno(3))
      ENDIF

C
C    CORRECTION DUE TO SOIL MOISTURE
C

      IF (OPTsoil) THEN
         read(20,rec=istcat) jp,sta
         call atmo_sch(datj,istsoi,soista,soimjd,
     *                soix,soiy,soiz,
     *                sta,csoi(1),csoi(2),csoi(3))
      ENDIF
c- jboehm

C
C     CORRECTION DUE TO POLE TIDE
C

      IF (OPTPT) THEN

         CALL CPOLTD (GLONG, PHI, YEAR, RX, RY, RZ, XJ, YJ,
     *                   CPOL(1), CPOL(2), CPOL(3))
      ENDIF
C
C    COMPUTE TOTAL CORRECTION TO STATION COORDINATES AND RETURN
C
      DO I=1,3

c
c     Remove atmosphere loading at all due to strange problem with auxiliar.dat
c     on 02-12-2011
c

ccc       CVECT(i) = CVECT(i)

         CVECT(I) = CVECT(I)+CTIDE(I)+COCE(I)+         CPOL(I)

c         print *, catm(1), catm(2), catm(3)

c         CVECT(I) = CVECT(I)+CTIDE(I)+COCE(I)+CATM(I)+CPOL(I)
c+ jboehm
     .   + csno(i) + csoi(i) + coce1(i)
c- jboehm

c        write (13,800) is, cvect(i), ctide(i), coce(i),
c     *   catm(i), cpol(i)
c 800     format (i3,5(2x,f18.6))

      end do

      RETURN
      END


***********************************************************************
C    BEGINNING OF THE PACKAGE FOR THE OCEAN LOADING CORRECTION
***********************************************************************

      SUBROUTINE COCLD (TIMDT,X, Y, Z, IST, CX, CY, CZ)

C***********************************************************************
C
C OCEAN LOADING SITE DISPLACEMENT
C REF.: MERIT STANDARDS, APRIL 81, SECOND DRAFT, APPENDICES 7,11
C       AND C. C. GOADS, J. GEOPHYS. RES., 85, P. 2679 - 2683, MAY 1980
C
C PROGRAMMER: H. SCHUH  (ORIGINAL VERSION CALLED COMPHD) MODIFIED BY
C N. ZARRAOA (1989) IN ORDER TO CORRECT A SINGLE STATION
C
C COMPUTES THE STATION HEIGHT DISPLACEMENT FOR OCEAN LOADING
C (FOR 9 OCEAN TIDES)
C
C  LATEST REVISION  1990 NOVEMBER 19   (N.Z.L.)
C
C*********************************************************************

      IMPLICIT DOUBLE PRECISION  (A-H,O-Z)
      LOGICAL  NOOCL(20)
      DOUBLE PRECISION NSOCL
      DIMENSION  ANGLE(11)
C
      COMMON /OCLOAD/ NOOCL,ROCL(20,2,11),NSOCL(20,2,11),EWOCL(20,2,11)
      COMMON /PHYS/ C, FL, AH, AU, J20
      COMMON /MATH/ PI, TWOPI, ROG, ROS, ROZ, TOL

C
C-----------------------------------------------------------------------
C ORDER OF THE OCEAN TIDE QUANTITIES IN OCLOAD(I,J,K)
C (AMPLITUDE(I=1) AND PHASE(I=2))
C
C J = 2 - M2
C J = 3 - S2
C J = 4 - N2
C J = 5 - K2
C J = 6 - K1
C J = 7 - O1
C J = 8 - P1
C J = 9 - Q1
C J =10 - MF
C

C
      DELTAR = 0.D0
      DELTNS = 0.D0
      DELTEW = 0.D0

C
C-----------------------------------------------------------------------
C P - GEOCENTRIC LATITUDE
C G - GEOCENTRIC LONGITUDE
C-----------------------------------------------------------------------

      P = DATAN(Z / (DSQRT(X**2 + Y**2)))
      G = DATAN2(Y,X)
      IF (Y .LT. 0.D0)  G = G + TWOPI
      R0 = DSQRT(X**2 + Y**2 + Z**2)
      CF = DCOS(P)
      SF = DSIN(P)
      CL = DCOS(G)
      SL = DSIN(G)
C
      CALL ANG (TIMDT, ANGLE)


      DO  J = 1, 11
C-----------------------------------------------------------------------
C COMPUTATION OF HEIGHT DISPLACEMENT (CM)
C-----------------------------------------------------------------------
           DELTAR = ROCL(IST,1,J) * DCOS(ANGLE(J)-
     1             (ROCL(IST,2,J)/ROG)) + DELTAR
           DELTNS = NSOCL(IST,1,J)* DCOS(ANGLE(J)-
     1             (NSOCL(IST,2,J)/ROG)) + DELTNS
           DELTEW = EWOCL(IST,1,J)* DCOS(ANGLE(J)-
     1             (EWOCL(IST,2,J)/ROG)) + DELTEW

      END DO

C*******************************************************************
C   CHANGE THE SIGN OF THE HORIZONTAL CONTRIBUTIONS WHEN SCHERNECK
C   VALUES ARE USED FOR AMPLITUDE AND PHASE.
C   THIS STRATEGY FOLLOWS CALC V7.0         NZL 11/10/90
C
C   CHANGE OF SIGN HAS BEEN ELIMINATED. IN THIS ROUTINE, OCCAM TAKES
C   N-S AND E-W AS POSITIVE TOWARDS SOUTH AND WEST RESPECTIVELY,
C   WHICH IS IN ACCORDANCE WITH THE IERS STANDARDS (1992) AND OPPOSITE
C   TO THE CALC'S APPROACH.                 NZL  5/MARCH/93
C*******************************************************************

c      DELTAR = DELTAR / 100.D0    ! from cm -> m
c      DELTNS = DELTNS / 100.D0    !
c      DELTEW = DELTEW / 100.D0    !

C-----------------------------------------------------------------------
C CX,CY,CZ, CORRECTION TO X,Y,Z COMPONENTS OF STATION VECTOR
C-----------------------------------------------------------------------

      CX = DELTAR*CF*CL + DELTNS*SF*CL + DELTEW*SL       !  m
      CY = DELTAR*CF*SL + DELTNS*SF*SL - DELTEW*CL       !  m
      CZ = DELTAR*SF    - DELTNS*CF                         !  m


      RETURN
      END
C
C***********************************************************************
C
      SUBROUTINE ANG (timdt, ANGLE)
C
C-----------------------------------------------------------------------
C PROGRAMMER: C. GOAD, H. SCHUH
C
C REF.: MERIT STANDARDS, APPENDIX 11
C
C COMPUTES THE ANGULAR ARGUMENT WHICH DEPENDS ON TIME
C FOR 9 TIDAL ARGUMENT CALCULATIONS
C
C DATE - JULIAN DATE WITH LAST FIVE FIGURES. EX. 44040.45366
C
C ANGLE - ANGULAR ARGUMENT FOR SCHWIDERSKI COMPUTATION
C
C   C A U T I O N
C   =============
C
C SCHWIDERSKI MODIFIES THE ANGULAR ARGUMENTS OF THE DIURNAL TERMS
C BY +/- 90 DEGREES. THEREFORE HIS DIURNAL PHASES CANNOT BE USED
C WITH THE STANDARD DOODSEN OR CARTWRIGHT CONVENTIONS
C-----------------------------------------------------------------------
C  LATEST REVISION  1990 NOVEMBER 19   (N.Z.L.)
C
      IMPLICIT DOUBLE PRECISION  (A-H,O-Z)
      DIMENSION  ANGFAC(4,11), ANGLE(11), SPEED(11)
C
      COMMON /PHYS/ C, FL, AH, AU, J20
      COMMON /MATH/ PI, TWOPI, ROG, ROS, ROZ, TOL

C-----------------------------------------------------------------------
C SPEED OF ALL TERMS IN RADIANS PER SEC
C-----------------------------------------------------------------------
      EQUIVALENCE  (SPEED(1),SIGM2), (SPEED(2),SIGS2), (SPEED(3),SIGN2),
     1             (SPEED(4),SIGK2), (SPEED(5),SIGK1), (SPEED(6),SIGO1),
     2             (SPEED(7),SIGP1), (SPEED(8),SIGQ1), (SPEED(9),SIGMF),
     3             (SPEED(10),SIGMM), (SPEED(11),SIGSSA)
C
      DATA  SIGM2  / 1.40519D-4 /
      DATA  SIGS2  / 1.45444D-4 /
      DATA  SIGN2  / 1.37880D-4 /
      DATA  SIGK2  / 1.45842D-4 /
      DATA  SIGK1  /  .72921D-4 /
      DATA  SIGO1  /  .67598D-4 /
      DATA  SIGP1  /  .72523D-4 /
      DATA  SIGQ1  /  .64959D-4 /
      DATA  SIGMF  / .053234D-4 /
      DATA  SIGMM  / .026392D-4 /
      DATA  SIGSSA / .003982D-4 /
      DATA  ANGFAC / 2.D0, -2.D0, 2*0.D0, 4*0.D0,
     3               2.D0, -3.D0, 1.D0, 0.D0, 2.D0, 3*0.D0,
     5               1.D0, 2*0.D0, 0.25D0, 1.D0, -2.D0,0.D0 ,-0.25D0,
     7               -1.D0, 2*0.D0, -0.25D0, 1.D0, -3.D0, 1.D0, -0.25D0,
     9               0.D0, 2.D0, 2*0.D0, 0.D0, 1.D0, -1.D0, 0.D0,
     A               2.D0, 3*0.D0 /

C-----------------------------------------------------------------------
C ID - DAY OF YEAR
C-----------------------------------------------------------------------

      ID = timdt

C-----------------------------------------------------------------------
C FDAY - FRACTIONAL PART OF DAY IN SECONDS
C-----------------------------------------------------------------------
      FDAY  = (timdt - ID) * 86400.D0

C  NO. OF DAYS FROM 1975.0

      ICAPD = ID-42412
      CAPT  = (27392.500528D0 + 1.000000035D0 * ICAPD) / 36525.D0

C-----------------------------------------------------------------------
C H0 - MEAN LONGITUDE OF SUN AT BEGINNING OF DAY
C-----------------------------------------------------------------------
      H0 = (279.69668D0 + (36000.768930485D0 + 3.03D-4 * CAPT) * CAPT)
     *     / ROG

C-----------------------------------------------------------------------
C S0 - MEAN LONGITUDE OF MOON AT BEGINNING OF DAY
C-----------------------------------------------------------------------

      S0 = (((1.9D-6 * CAPT - .001133D0) * CAPT + 481267.88314137D0)
     *     * CAPT + 270.434358D0) / ROG

C-----------------------------------------------------------------------
C P0 - MEAN LONGITUDE OF LUNAR PERIGEE AT BEGINNING OF DAY
C-----------------------------------------------------------------------

      P0 = (((-1.2D-5 * CAPT - .010325D0) * CAPT + 4069.0340329577D0)
     *     * CAPT + 334.329653D0) / ROG

      DO  K = 1, 11
         ANGLE(K) = SPEED(K) * FDAY + ANGFAC(1,K) * H0 + ANGFAC(2,K)
     *              * S0 + ANGFAC(3,K) * P0 + ANGFAC(4,K) * TWOPI
         ANGLE(K) = DMOD(ANGLE(K),TWOPI)
         IF (ANGLE(K) .LT. 0.D0)  ANGLE(K) = ANGLE(K) + TWOPI
      END DO

      RETURN
      END

************************************************************************
C    END OF THE PACKAGE FOR OCEAN LOADING CORRECTION
************************************************************************

************************************************************************
C    ROUTINE FOR POLE TIDE CORRECTION
************************************************************************

      SUBROUTINE CPOLTD (GLONG, PHI, YEAR, RX, RY, RZ,
     *   XJ, YJ, DX, DY, DZ)
C
C-----------------------------------------------------------------------
C THIS SUBROUTINE COMPUTES THE DISPLACEMENTS OF THE VLBI-STATIONS
C CAUSED BY POLAR MOTION (POLE TIDE)
C DPHI CHANGED TO NEGATIVE SIGN (ACCORDING J.WAHR, 1985)
C
C REF.: CHUCK YODER, PRIV. COMM., MAR. 84
C       JOHN WAHR, 1985
C
C
C PROGRAMMER: H. SCHUH, NOV. 1984; N. ZARRAOA
C-----------------------------------------------------------------------
C       REVISION  1990 DECEMBER 5   (N.Z.L.)
C       REVISION  1993 MARCH 5   (NZL) DISPLACEMENTS REWRITTEN
C       REVISION  1993 DECEMBER 9 (NZL) PETROV COMMENTS
C       REVISION  1997 April   14 (OT) IERS - 96
C       REVISION  2001 January 24 (OT) Mean pole secular motion
C       REVISION  2003 June 13 (OT) Another signs for S(r), S(phi), S(lam)!!
C                         and new mean pole secular motion (IERS-2000)
C       REVISION  2003 June 16 (OT) New approach for Dx, Dy, Dz calculation
C  LAST REVISION  2004 September (OT) ERROR! Co-latitude instead of latitude!!!
c
c             ATTENTION !!!
c
c  IERS 1996 Conventions recommends to make use 'average' Love numbers
c  for pole tide, namely, h=0.6027, l=0.0836
c
C

      IMPLICIT DOUBLE PRECISION  (A-H,O-Z)
      DOUBLE PRECISION  L2, H2

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc TESMER 02.11.04 also add the bernese mean pole for pole tide cc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      CHARACTER YESLIN*1
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc TESMER 02.11.04 also add the bernese mean pole for pole tide cc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      COMMON /PHYS/ C, FL, AH, AU, J20
      COMMON /MATH/ PI, TWOPI, ROG, ROS, ROZ, TOL

      EARTHR = AH
      H2=0.6027d0
      L2=0.0836d0
      R0 = DSQRT(RX**2 + RY**2 + RZ**2)

C-----------------------------------------------------------------------
C P - GEOCENTRIC LATITUDE;   G - GEOCENTRIC LONGITUDE
C-----------------------------------------------------------------------

      CG   = DCOS(GLONG)
      SG   = DSIN(GLONG)
      p    = pi/2.d0 - phi

C--------------------------------------------------------------
C     CALCULATION OF MEAN POLE COORDINATES
C--------------------------------------------------------------

        CALL MEAN_POLE (YEAR, XJ1, YJ1)    !  X,Y - mas (units)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc TESMER 02.11.04 also add the bernese mean pole for pole tide cc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      xj =   xj*ros - xj1/1000.d0
      yj = -(yj*ros - yj1/1000.d0)

C---------------------------------------------------------
C   COMPUTE VERTICAL DISPLACEMENT BEING POSITIVE UPWARDS
C---------------------------------------------------------

c      DR = -32.d0 * DSIN (2.d0*P) * (XJ*CG + YJ*SG)
      DR = -33.d0 * DSIN (2.d0*P) * (XJ*CG + YJ*SG)  !  Changed on 19.09.2006
      dr =  dr/1000.d0            ! mm -> m         16.06.2003

C-------------------------------------------------------------------
C   COMPUTE HORIZONTAL DISPLACEMENT ALONG THE NORTH SOUTH DIRECTION
C  BEING POSITIVE TOWARDS THE NORTH
C-------------------------------------------------------------------

      DPHI = -9.d0 * DCOS(2.d0*P) * (XJ*CG + YJ*SG)

C-----------------------------------------------------------------------
C CONVERT DPHI FROM (M) TO RAD
C-----------------------------------------------------------------------
      DPHI = DPHI/1000.d0        ! mm -> m        16.06.2003

C-------------------------------------------------------------------
C   COMPUTE HORIZONTAL DISPLACEMENT ALONG THE EAST WEST DIRECTION
C  BEING POSITIVE TOWARDS THE EAST
C-------------------------------------------------------------------

      DLAM = 9.d0 * DCOS(P) * (XJ*SG - YJ*CG)
      dlam = dlam/1000.d0       ! mm -> m         16.06.2003

C-----------------------------------------------------------------------
C   P - COLATITUDE !!
C-----------------------------------------------------------------------


      DX =  cg*dcos(p)*dphi - sg*dlam + cg*dsin(p)*dr      !  m
      DY =  sg*dcos(p)*dphi + cg*dlam + sg*dsin(p)*dr      !  m
      DZ =    -dsin(p)*dphi           +    dcos(p)*dr      !  m

      RETURN
      END


C************************************************************************
C


      SUBROUTINE MEAN_POLE (YEAR, XJ, YJ)

C
C        PREPARED   14, JULY, 2000  (OLEG TITOV)
C        UPDATED    20, JANUARY, 2001  (OLEG TITOV) IERS CONV. 2000
C   LAST UPDATED    13, JUNE, 2003  (OLEG TITOV) NEW IERS CONV. 2000
C
C*********************************************************************

      IMPLICIT DOUBLE PRECISION  (A-H,O-Z)

c
c     Recommended by Harald
c

c      X0 = -0.7d0       !  mas
c      Y0 = 179.4d0      !  mas
c      X1 = 0.862d0      !  mas/year
c      Y1 = 3.217d0      !  mas/year
c      YEAR0 = 1950.d0

C     The values are recommended by IERS Conv. 2000 (McCarthy, 2002)

      X0 = 54.d0       !  mas
      Y0 = 357.d0      !  mas
      X1 = 0.83d0      !  mas/year
      Y1 = 3.95d0      !  mas/year
      YEAR0 = 2000.d0

      XJ = X0 + X1*(YEAR - YEAR0)
      YJ = Y0 + Y1*(YEAR - YEAR0)

      RETURN
      END


***********************************************************************
C    END OF THE PACKAGE FOR THE EARTH TIDES CORRECTION
***********************************************************************



      SUBROUTINE INTPT1 (FNAME,FLAT,FLON,UN)

CCCCCCCCC+CCCCCCCCC+CCCCCCCCC+CCCCCCCCC+CCCCCCCCC+CCCCCCCCC+CCCCCCCCC+C-CCCCCCC+
c
c  This subroutine is based on the program 'intpt.f' by R.H. Rapp that was
c  created to interpolate geoidal heights for EGM96. Further information and
c  the original source code can be found at:
c
c  http://164.214.2.59/GandG/wgs-84/egm96.html
c
c  The subroutine intpt1 can be used to interpolate values in the
c  z200-grid (geopotential heights of the 200 hPa pressure level) and the
c  smfw3-grid that are used for the isobaric mapping functions by A. Niell.
c
c  Attention: Some parameters might have to be changed in the source code.
c             They are marked with 'jboehm'
c
c  input parameters
c  ----------------
c  fname:  name of the ascii file, e.g. 'h94012.h00' or 'w94012.h00'
c          for the format description see below
c  flat:   latitude  in radians
c  flon:   longitude in radians
c
c  output parameters
c  -----------------
c  un:     interpolated value
c
c  jboehm, 2002 May 8
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


c ... original description is following
C
C GRID FILE STRUCTURE:
C
C  This FORTRAN program will interpolate values from a WORLD-WIDE grid file.  
C  The grid file must have a header as record 1 with the following format:
C
C       south    north    west    east    spacing N-S   spacing W-E
C
C  The grid data are in subsequent records arranged from North to South, West 
C  to East, i.e. record 2 would be the northern most latitude band of data.
C                                                                         
C-------------------------------------------------------------------------------
C
C  Several interpolation schemes are possible using this program. By setting
C  the parameter IWINDO to a particular value, one may use Bi-linear interpo-
C  lation or a spline interpolation with a window of size of IWINDO x IWINDO.
C  Check the subroutine "INTERP" for more information.                    
C                                                                      
C  This program is set up to extend the longitude values by a minimum of 2 posts 
C  on the positive side of the east longitude limit and a minimum of 2 posts to 
C  the negative side of the west longitude.          
C
C-------------------------------------------------------------------------------
C
C  Change the following variables in the parameter statements below to match
C  the grid file information and desired interpolation method:
C
C                IWINDO = Size of window of interpolation and method. 
C                             IWINDO=0 gives bilinear interpolation. 
C                             Check subroutine "INTERP" below for more 
C                             information.
C
C                ****NOTE: If IWINDO=0 be sure to set NBDR=4 *****
C                             In general, NBDR=2*IWINDO
C
C
C                NLAT = Number of values in latitude.
C                NLON = Number of values in longitude. (Leave the NBDR 
C                            parameter in the NLON parameter declaration.
C                            This controls the extension of the grid in 
C                            longitude.)
C                DLAT = Spacing in latitude.
C                DLON = Spacing in longitude.
C
C-------------------------------------------------------------------------------
C
C FILE INFORMATION:
C
C  INPUT DATA GRID: UNIT = 1
C
C    For structure see above.
C
C
C-------------------------------------------------------------------------------
C
C  This program was provided by Professor Richard H. Rapp of The Ohio State
C    University, December, 1996.  It was put into its present form by the National
C    Imagery and Mapping Agency (NIMA), December, 1996.  
c
C
CCCCCCCCC+CCCCCCCCC+CCCCCCCCC+CCCCCCCCC+CCCCCCCCC+CCCCCCCCC+CCCCCCCCC+C-CCCCCCC+

      implicit real*8(a-h,o-z)
      real*8 south,north,west,east,dphi,dlam

c+  jboehm, mind that the following lines might have to be changed
      PARAMETER(IWINDO=4)
      PARAMETER(NBDR=2*IWINDO)
      PARAMETER(NLAT=91,NLON=(145+NBDR),DLAT=2D0,DLON=2.5D0)
c-  jboehm

      DIMENSION A1(NLON)
      REAL*4 H (nlat, nlon)
      CHARACTER NAM*41
      REAL*8 XMIN,XMAX,RDX
      DATA XMIN, XMAX, RDX, NAM/-50.,50., 2.,
     .     'GEOID HTS'/
      character*10 fname

c  coordinates from radians to degrees
      RHO   = 57.29577951D0
      flat1 = flat*rho
      flon1 = flon*rho

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCC                                                           CCCCCCCCCCC
CCCCCCCCCC                 INPUT/OUPUT FILE NAMES                    CCCCCCCCCCC
CCCCCCCCCC            MAKE CHANGES TO FILE NAMES BELOW               CCCCCCCCCCC
CCCCCCCCCC                                                           CCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c data grid file

c+  jboehm, mind to set the path, e.g.
      open(1,file='../imf/'//fname,status='old')
c      open(1,file=fname,status='old')
c-  jboehm

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c read header of input grid file.
      read(1,*)south,north,west,east,dphi,dlam
c     write(*,*)'input file header:'
c     write(*,998)south,north,west,east,dphi,dlam,
c    &  int((north-south)/dphi+1),int((east-west)/dlam+1)
c     write(*,*)'Southwest point in grid storage array:'
c     write(*,998)south,west-nbdr/2*dlam
998   format(6f8.2,2i6)


C read input grid file and store in array H
      DO I=1,NLAT
         READ(1,*) (A1(K),K=1,nlon-nbdr)
         DO J=1,nlon-nbdr
            H(NLAT-I+1,J+nbdr/2)=A1(J)
         ENDDO
         DO K=1,nbdr/2
            H(NLAT-I+1,K)=A1(nlon-nbdr-nbdr/2+K)
3           H(NLAT-I+1,k+nlon-nbdr/2)=A1(K)
         ENDDO
      ENDDO

c  jboehm
      close(1)

c call interpolation subroutine
      CALL INTERP(IWINDO,12.D0,H,south,west-nbdr/2*dlam,DLAT,DLON,
     .            NLAT,NLON,NLAT,NLON,FLAT1,FLON1,UN)

      RETURN
      END

      SUBROUTINE INTERP(IWINDO,DMIN,H,PHIS,DLAW,DDFI,DDLA,NPHI,NDLA,
     .                  IPDIM,ILDIM,PHI,DLA,VALINT)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE FOR INTERPOLATION OF VALUES FROM A STANDARD DTM-GRID  C
C     TO INDIVIDUAL STATION LOCATIONS.                                 C
C                                                                      C
C                                                                      C
C     INPUT PARAMETERS...                                              C
C     ===================                                              C
C     IWINDO...    A SPLINE WINDOW OF SIZE 'IWINDO' X 'IWINDO' WILL BE C
C                  USED AROUND EACH STATION. IF 'IWINDO' IS 0 OR 1,    C
C                  BILINEAR INTERPOLATION WILL BE USED.                C
C     DMIN...      MINIMUM ACCEPTABLE DISTANCE FROM THE GRID EDGE IN   C
C                  KM (USEFUL FOR FFT GRIDS).                          C
C     H...         2D DATA ARRAY (ELEMENT (1,1) IN SW CORNER).         C
C     PHIS,DLAW... LATITUDE AND LONGITUDE OF SW GRID POINT.            C
C     DDFI,DDLA... GRID SPACING IN LATITUDE AND LONGITUDE DIRECTION.   C
C     NPHI,NDLA... NUMBER OF GRID POINTS IN LATITUDE AND LONGITUDE     C
C                  DIRECTION.                                          C
C     IPDIM,ILDIM..DIMENSIONS OF 2D DATA ARRAY 'H' AS DECLARED IN THE  C
C                  CALLING PROGRAM.                                    C
C     PHI,DLA...   LATITUDE AND LONGITUDE OF INTERPOLATION POINT.      C
C                                                                      C
C                                                                      C
C     OUTPUT PARAMETERS...                                             C
C     ====================                                             C
C     VALINT...    INTERPOLATED VALUE.                                 C
C                                                                      C
C                                                                      C
C     EXECUTION TIME ON CDC 990 IS...                                  C
C     ===============================                                  C
C     +------------------+-------------------+-------------------+     C
C     I  INTERPOLATION   I  OPT=LOW          I  OPT=HIGH         I     C
C     I------------------I-------------------I-------------------I     C
C     I  BILINEAR        I  1.44 MSEC/STAT.  I  1.44 MSEC/STAT.  I     C
C     I  3 X 3 SPLINE    I  1.53 MSEC/STAT.  I  1.51 MSEC/STAT.  I     C
C     I  5 X 5 SPLINE    I  1.70 MSEC/STAT.  I  1.67 MSEC/STAT.  I     C
C     I  7 X 7 SPLINE    I  2.02 MSEC/STAT.  I  1.74 MSEC/STAT.  I     C
C     I  9 X 9 SPLINE    I  2.31 MSEC/STAT.  I  2.00 MSEC/STAT.  I     C
C     +------------------+-------------------+-------------------+     C
C                                                                      C
C                                                                      C
C     PROGRAM CREATION BY...   H. DENKER          MAY 30, 1987         C
C                              H. DENKER          MARCH 13, 1989       C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      PARAMETER(IPA1=20)
      IMPLICIT REAL*8(A-H,O-Z)
      LOGICAL LODD
      REAL*4 H(IPDIM,ILDIM)
      DIMENSION A(IPA1),R(IPA1),Q(IPA1),HC(IPA1)
      IDIM1=IPA1
      RHO=57.29577951D0
      REARTH=6371000.D0
      IF(IWINDO.LT.2) IWINDO=2
      IF(IWINDO.GT.IDIM1) IWINDO=IDIM1
      ILIM=DMIN*1000.*RHO/(REARTH*DDFI)
      JLIM=DMIN*1000.*RHO/(REARTH*DDLA*COS((PHIS+DDFI*NPHI/2.)/RHO))
      LODD=(IWINDO/2)*2.NE.IWINDO
      RI=(PHI-PHIS)/DDFI
      RJ=(DLA-DLAW)/DDLA
      IF(LODD) THEN
        I0=RI-0.5
        J0=RJ-0.5
      ELSE
        I0=RI
        J0=RJ
      ENDIF
      I0=I0-IWINDO/2+1
      J0=J0-IWINDO/2+1
      II=I0+IWINDO-1
      JJ=J0+IWINDO-1
      IF(I0.LT.0 .OR. II.GE.NPHI .OR. J0.LT.0 .OR. JJ.GE.NDLA) THEN
        WRITE(6,7008) PHI,DLA
        VALINT=999999.
        RETURN
      ELSEIF(I0.LT.ILIM .OR. II.GT.NPHI-ILIM .OR. J0.LT.JLIM .OR.
     .  JJ.GT.NDLA-JLIM) THEN
c        IF(NPOINT.LE.ILIST) WRITE(6,7009) PHI,DLA
        VALINT=999999.
        RETURN
      ENDIF
7008  FORMAT(' ',2F10.6,' STATION TOO NEAR GRID BOUNDARY  - NO INT.'
     .,' POSSIBLE|')
7009  FORMAT(' ',2F10.6,' STATION OUTSIDE ACCEPTABLE AREA - NO INT.'
     .,' PERFORMED|')
      IF(IWINDO.GT.2) THEN
        DO 110 I=1,IWINDO
          DO 111 J=1,IWINDO
            A(J)=H(I0+I,J0+J)
111       CONTINUE
          CALL INITSP(A,IWINDO,R,Q)
          HC(I)=SPLINE(RJ-J0+1.,A,IWINDO,R)
110     CONTINUE
        CALL INITSP(HC,IWINDO,R,Q)
        VALINT=SPLINE(RI-I0+1.,HC,IWINDO,R)
      ELSE
        VALINT=BILIN(RI+1.,RJ+1.,H,NPHI,NDLA,IPDIM,ILDIM)
      ENDIF
      RETURN
      END

      FUNCTION BILIN(RI,RJ,A,IMAX,JMAX,IADIM,JADIM)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C                           B I L I N                                  C
C                                                                      C
C  INTERPOLATES VALUES IN AN ARRAY A USING BILINEAR                    C
C  (PARABOLIC HYPERBOLOID) INTERPOLATION.                              C
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
C  PARAMETERS:                                                         C
C                                                                      C
C  BILIN...       INTERPOLATED VALUE                                   C
C                                                                      C
C  RI, RJ...      INTERPOLATION ARGUMENT, (1,1) IN LOWER LEFT CORNER,  C
C                 (IMAX, JMAX) IN UPPER RIGHT.                         C
C                                                                      C
C  A...           INTEGER*2 ARRAY WITH ARGUMENTS                       C
C                                                                      C
C  IMAX, JMAX...  NUMBER OF POINTS IN GRID                             C
C                                                                      C
C  IADIM, JADIM...DECLARED DIMENSIONS OF 'A'                           C
C                                                                      C
C  OUTSIDE AREA COVERED BY 'A' THE FUNCTION RETURNS THE VALUE OF       C
C  THE NEAREST BOUNDARY POINT.                                         C
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
C  PROGRAMMER:                                                         C
C  RENE FORSBERG, JULY 1983                                            C
C                                                                      C
C  MODIFICATIONS BY:                                                   C
C  HEINER DENKER, 07/01/1987                                           C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IMPLICIT REAL*8 (A-H, O-Z)
      DIMENSION A(IADIM, JADIM)
      REAL*4 A
      IN = IFRAC(RI)
      IE = IFRAC(RJ)
      RN = RI - IN
      RE = RJ - IE
      IF (IN.LT.1) THEN
        IN = 1
        RN = 0.0
      ELSEIF (IN.GE.IMAX) THEN
        IN = IMAX-1
        RN = 1.0
      ENDIF
      IF (IE.LT.1) THEN
        IE = 1
        RE = 0.0
      ELSEIF (IE.GE.JMAX) THEN
        IE = JMAX-1
        RE = 1.0
      ENDIF
      RNM1=1.-RN
      REM1=1.-RE
      BILIN = RNM1*REM1*A(IN,IE) +
     .RN*REM1*A(IN+1,IE) + RNM1*RE*A(IN,IE+1) +
     .RN*RE*A(IN+1,IE+1)
      RETURN
      END

      SUBROUTINE INITSP(Y, N, R, Q)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C                      I N I T S P                                     C
C                                                                      C
C  INITIALIZATION PROCEDURE FOR FAST 1-DIMENSIONAL EQUIDISTANT         C
C  SPLINE INTERPOLATION, WITH FREE BOUNDARY END CONDITIONS             C
C  REFERENCE: JOSEF STOER: EINFUHRUNG IN DIE NUMERISCHE MATHEMATIK     C
C  I, SPRINGER 1972, PAGE 82 AND 86.                                   C
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
C  PARAMETERS (REAL):                                                  C
C                                                                      C
C  Y...   GIVEN VALUES, Y(1), ..., Y(N)                                C
C                                                                      C
C  R...   SPLINE MOMENTS (1 ... N), TO BE USED BY FUNCTION 'SPLINE'    C
C                                                                      C
C  Q...   WORK-ARRAY, DECLARED AT LEAST 1:N                            C
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
C  RENE FORSBERG, JULY 1983                                            C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Y(N), R(N), Q(N)
      Q(1) = 0.0
      R(1) = 0.0
      DO 11 K = 2, N-1
        P = Q(K-1)/2+2
        Q(K) = -0.5/P
        R(K) = (3*(Y(K+1)-2*Y(K)+Y(K-1)) - R(K-1)/2)/P
   11 CONTINUE
      R(N) = 0.0
      DO 12 K = N-1, 2, -1
        R(K) = Q(K)*R(K+1)+R(K)
   12 CONTINUE
      RETURN
      END

      FUNCTION SPLINE(X, Y, N, R)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C                          S P L I N E                                 C
C                                                                      C
C  FAST ONE-DIMENSIONAL EQUIDISTANT SPLINE INTERPOLATION FUNCTION.     C
C  REFERENCE: JOSEF STOER: EINFUHRUNG IN DIE NUMERISCHE MATHEMATIK     C
C  I, SPRINGER 1972, PAGE 81.                                          C
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
C  PARAMETERS:                                                         C
C                                                                      C
C  X...  INTERPOLATION ARGUMENT (REAL), X = 1 FIRST DATA-POINT,        C
C        X = N LAST DATA-POINT. OUTSIDE THE RANGE LINEAR EXTRA-        C
C        POLATION IS USED.                                             C
C                                                                      C
C  Y...  REAL*8 ARRAY, 1 .. N : DATA VALUES                            C
C                                                                      C
C  R...  DO: SPLINE MOMENTS CALCULATED BY SUBROUTINE 'INITSP'          C
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
C  PROGRAMMER:                                                         C
C  RENE FORSBERG, JUNE 1983                                            C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IMPLICIT REAL*8 (A-H, P-Z)
      DIMENSION Y(N), R(N)
      IF (X.LT.1) THEN
        SPLINE = Y(1) + (X-1)*(Y(2)-Y(1)-R(2)/6)
      ELSEIF (X.GT.N) THEN
        SPLINE = Y(N) + (X-N)*(Y(N)-Y(N-1)+R(N-1)/6)
      ELSE
        J = IFRAC(X)
        XX = X - J
        SPLINE = Y(J) +
     .           XX * ((Y(J+1)-Y(J)-R(J)/3-R(J+1)/6) +
     .           XX * (R(J)/2 +
     .           XX * (R(J+1)-R(J))/6))
      ENDIF
      RETURN
      END

      FUNCTION IFRAC(R)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C                    F U N C T I O N   I F R A C                       C
C                  ===============================                     C
C                                                                      C
C  SUBROUTINE GIVING TRUE INTEGER PART OF A REAL E.G.                  C
C                                                                      C
C    FOR   1. = R < 2.   IS    IFRAC = 1                               C
C    FOR   0. = R < 1.   IS    IFRAC = 0                               C
C    FOR  -1. = R < 0.   IS    IFRAC =-1                               C
C    FOR  -2. = R <-1.   IS    IFRAC =-2                               C
C                                                                      C
C  RF, JUNE 1983                                                       C
C  HD, JANUARY 1987                                                    C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IMPLICIT REAL*8(A-H,O-Z)
      IFRAC=R
      IF (R.GE.0) RETURN
      IF (R.EQ.IFRAC) RETURN
      IFRAC = IFRAC - 1
      RETURN
      END


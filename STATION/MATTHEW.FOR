      SUBROUTINE matthew (tim,HSG,X1,X2,X3,
     *    R_after_PN, V_after_PN, XLUN_after_PN,
     *    rrot,ds1,ds2,ds3)
C
C-----------------------------------------------------------------------
c       matthews - tides modell is called here
c
c       gives to matthews:      xsit (1-3)  [m]    station coord.
c                               xsol (1-3)  [m]    geoz. sun ephems
c                               xlun (1-3)  [m]    geoz. monn ephems
c                               iyr
c                               imonth
c                               iday
c                               fhr
c
c       gets from matthews:     dsit (1-3)  [m]    station displacements
c
c       Tesmer 9.5.00
C       Revision 24.09.2001 (OT)    Small cosmetic changes
C       Last revision 04.06.2010 (OT)    Matrix R4 added
C
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION  (A-H,O-Z)
      DIMENSION XSOL(3), XLUN(3), XSUN(3), XSIT(3), DSIT(3), rrot(3,3)
      DIMENSION XSOL1(3), XLUN1(3), dsit1(3)
      dimension V_after_PN(3), R_after_PN(3), XLUN_after_PN(3)
C
c      COMMON /BAR1/ RBAR0(3),VBAR0(3)
c      COMMON /MOON/ XMOON(3)
      common /phys/ c

      XSIT(1) = X1
      XSIT(2) = X2
      XSIT(3) = X3

C
C***********************************************************************
C   GET APPARENT RECTANGULAR GEOCENTRIC COORDINATES OF THE SUN
C   SUBSTRACTING THE LIGHT TIME FROM SUN TO EARTH  (NZL 4/12/90)
C   ROTATE AN ANGLE HSG TO PUT THE SUN ON TERRESTRIAL FRAME
C***********************************************************************

C
C   Calculation of RBAR0 has been changed in STATION.FOR in September, 2001
C   The geocentric velocities are from EPHEM.DAT instead of interpolation
C   from coordinates
C

      dsol = 0.d0

	do i=1,3
         dsol=dsol + R_after_PN(i)**2
c	   write (13,*) r_after_pn(i), v_after_pn(i), xlun_after_pn(i)
      end do

	dsol = dsqrt(dsol)

      do i=1,3

        xsun(i) = R_after_PN(i) - dsol/c * V_after_PN(i)

      end do




c  rotate by inv(R3*R2*R1) to have in TRF - related frame

      do i=1,3
 	   xsol(i) = 0.d0
         xlun(i) = 0.d0

c         write (13,*) (rrot(k,i), k=1,3)

	  do k=1,3

            xsol(i) = xsol(i) + rrot(k,i)*xsun(k)
            xlun(i) = xlun(i) + rrot(k,i)*xlun_after_pn(k)

	  end do

        

c        write (13,*) 'after ROT', xsol(i), xlun(i)

	end do

c      print *, xsun(1), xsun(2), xsun(3)
c      print *, xsol(1), xsol(2), xsol(3)



c      XSOL1(1) = XSUN(1)*DCOS(HSG) + XSUN(2)*DSIN(HSG)
c      XSOL1(2) = -XSUN(1)*DSIN(HSG) + XSUN(2)*DCOS(HSG)
c      XSOL1(3) = XSUN(3)

c      XLUN1(1) = XMOON(1)*DCOS(HSG) + XMOON(2)*DSIN(HSG)
c      XLUN1(2) = -XMOON(1)*DSIN(HSG) + XMOON(2)*DCOS(HSG)
c      XLUN1(3) = XMOON(3)

C  ROTATE BY HSG

      XSOL1(1) = XSUN(1)*DCOS(HSG) + XSUN(2)*DSIN(HSG)
      XSOL1(2) = -XSUN(1)*DSIN(HSG) + XSUN(2)*DCOS(HSG)
      XSOL1(3) = XSUN(3)

      XLUN1(1) = Xlun_after_pn(1)*DCOS(HSG) + Xlun_after_pn(2)*DSIN(HSG)
      XLUN1(2) =-XLUN_after_PN(1)*DSIN(HSG) + XLUN_after_PN(2)*DCOS(HSG)
      XLUN1(3) = XLUN_after_PN(3)

c the routine TIDE needs time in IYear, IMonth, IDay, fraction of day [hours]

      TIMJUL=TIM + 2400000.5D0

      call datjul (TIMJUL,iYEAR,iMONTH,DAY)

      iDAY=int(DAY)
      FHR=(DAY-iDAY)*24.D0

c      write (13,*) xsit(1), xsit(2), xsit(3)


c TIDE is the matthews model, coded by dehant et al.

      call TIDE (XSIT,iYEAR,iMONTH,iDAY,FHR,XSOL,XLUN,DSIT)
      call TIDE (XSIT,iYEAR,iMONTH,iDAY,FHR,XSOL1,XLUN1,DSIT1)

c      print *,'dd = ',dsit(1)-dsit1(1),dsit(2)-dsit1(2),dsit(3)-dsit1(3)

      ds1=dsit(1)
      ds2=dsit(2)
      ds3=dsit(3)

      return
      end

C****************************************************************************************
      subroutine datjul(jd,year,month,day)

ccc   Meeus, Jean: Astronomische Algorithmen, Barth, 1992, S.76f.       ccc
ccc   Umwandlung des Julianischen Tages in buergerliches Datum nach     ccc
ccc   dem gregorianischen oder julianischen (vor dem 15.10.1582)        ccc
ccc   Kalender.                                                         ccc

      implicit none
      integer year,month,z,a,b,c,d,e,alpha
      double precision jd,day,f

      z = jd + 0.5D0
      f = jd + 0.5D0 - z
      if (z.lt.2299161) then
        a = z
      else
        alpha = (z-1867216.25D0)/36524.25D0
        a = z + 1 + alpha - int(alpha/4.D0)
      endif
      b = a + 1524
      c = (b-122.1D0)/365.25D0
      d = 365.25*c
      e = (b-d)/30.6001D0
      day = b - d - int(30.6001*e) + f
      if (e.lt.14) then
        month = e-1
      else
        month = e-13
      endif
      if (month.gt.2) then
        year = c - 4716
      else
        year = c - 4715
      endif
      return
      end


C****************************************************************************************

       SUBROUTINE TIDE(XSTA,iYR,iMONTH,iDAY,FHR,XSUN,XMON,DXTIDE)
C
C PURPOSE    :  COMPUTATION OF TIDAL CORRECTIONS OF STATION DISPLACEMENTS
C               CAUSED BY LUNAR AND SOLAR GRAVITATIONAL ATTRACTION
C               (SEE IERS STANDARDS 1996)
C
C         STEP 1 (HERE GENERAL DEGREE 2 AND 3 CORRECTIONS +
C                CALL ST1IDIU + CALL ST1ISEM + CALL ST1L1)
C         + STEP 2 (CALL STEP2DIU + CALL STEP2LON + CALL STEP2IDIU)
C
C IT HAS BEEN DECIDED THAT THE STEP 3 NON-CORRECTION FOR PERMANENT TIDE
C WOULD NOT BE APPLIED IN ORDER TO AVOID JUMP IN THE REFERENCE FRAME
C (THIS STEP 3 MUST ADDED IN ORDER TO GET THE NON-TIDAL STATION POSITION
C AND TO BE CONFORMED WITH THE IAG RESOLUTION.)
C
C    INPUT :  XSTA(I),I=1,2,3: GEOCENTRIC POSITION OF THE STATION
C             XSUN(I),I=1,2,3: GEOC. POSITION OF THE SUN
C             XMON(I),I=1,2,3: GEOC. POSITION OF THE MOON
C             IYR : YEAR
C             IMONTH : MONTH
C             IDAY : DAY
C             FHR=hr+zmin/60.+sec/3600. : HR IN THE DAY
C   OUTPUT :  DXTIDE(I),I=1,2,3: DISPLACEMENT VECTOR
C
C SUBROUTINES CALLED  :  SPROD
C                        ST1IDIU
C                        ST1ISEM
C                        ST1L1
C                        STEP2DIU
C                        STEP2LON
C                        STEP2IDIU
C                        STEP2ILON
C
C AUTHOR IERS 1996 :  V. DEHANT, S. MATHEWS AND J. GIPSON
C    (TEST BETWEEN TWO SUBROUTINES)
C
C AUTHOR IERS 2000 :  V. DEHANT AND S. MATHEWS
C    (TEST IN THE BERNESE PROGRAM BY C. BRUYNINX)
C
C CREATED    :  96/03/23              LAST MODIFIED :  00/05/17 14:10
C
C
C*************************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      COMMON /PHYS/ C, FL1, AH, AU, J20
      COMMON /GRAV/ GMS, GMM, GME
C
      DOUBLE PRECISION XSTA(3),XSUN(3),XMON(3),DXTIDE(3),XCORSTA(3)
      DOUBLE PRECISION H20,L20,H3,L3,H2,L2
      DOUBLE PRECISION mass_ratio_sun,mass_ratio_moon
C
C NOMINAL SECOND DEGREE AND THIRD DEGREE LOVE NUMBERS AND SHIDA NUMBERS
C ---------------------------------------------------------------------

      DATA H20/0.6078D0/,L20/0.0847D0/,H3/0.292D0/,L3/0.015D0/

C
C SCALAR PRODUCT OF STATION VECTOR WITH SUN/MOON VECTOR
C -----------------------------------------------------


      scs = xsun(1)*xsta(1)+xsun(2)*xsta(2)+xsun(3)*xsta(3)
      scm = xmon(1)*xsta(1)+xmon(2)*xsta(2)+xmon(3)*xsta(3)

      rsta = dsqrt(xsta(1)*xsta(1) + xsta(2)*xsta(2) + xsta(3)*xsta(3))
      rsun = dsqrt(xsun(1)*xsun(1) + xsun(2)*xsun(2) + xsun(3)*xsun(3))
      rmon = dsqrt(xmon(1)*xmon(1) + xmon(2)*xmon(2) + xmon(3)*xmon(3))

cc      CALL SPROD(XSTA,XSUN,SCS,RSTA,RSUN)
cc      CALL SPROD(XSTA,XMON,SCM,RSTA,RMON)
      SCSUN=SCS/RSTA/RSUN
      SCMON=SCM/RSTA/RMON


C
C COMPUTATION OF NEW H2 AND L2
C ----------------------------

      COSPHI=DSQRT(XSTA(1)**2+XSTA(2)**2)/RSTA
      H2=H20-0.0006d0*(1.d0-1.5d0*COSPHI**2)
      L2=L20+0.0002d0*(1.d0-1.5d0*COSPHI**2)

C
C P2-TERM
C -------

      P2SUN=3.d0*(H2/2.d0-L2)*SCSUN**2-H2/2.d0
      P2MON=3.d0*(H2/2.d0-L2)*SCMON**2-H2/2.d0

C
C P3-TERM
C -------

      P3SUN=2.5d0*(H3-3.d0*L3)*SCSUN**3+1.5d0*(L3-H3)*SCSUN
      P3MON=2.5d0*(H3-3.d0*L3)*SCMON**3+1.5d0*(L3-H3)*SCMON

C
C TERM IN DIRECTION OF SUN/MOON VECTOR
C ------------------------------------

      X2SUN=3.d0*L2*SCSUN
      X2MON=3.d0*L2*SCMON
      X3SUN=3.d0*L3/2.d0*(5.d0*SCSUN**2-1.d0)
      X3MON=3.d0*L3/2.d0*(5.d0*SCMON**2-1.d0)

C
C FACTORS FOR SUN/MOON
C --------------------

      MASS_RATIO_SUN= GMS/GME                        ! Changed on 04.01.2006
      MASS_RATIO_MOON = GMM/GME                      ! Changed on 04.01.2006
      RE = ah                                        ! Changed on 04.01.2006

      FAC2SUN=MASS_RATIO_SUN*RE*(RE/RSUN)**3
      FAC2MON=MASS_RATIO_MOON*RE*(RE/RMON)**3
      FAC3SUN=FAC2SUN*(RE/RSUN)
      FAC3MON=FAC2MON*(RE/RMON)

C
C TOTAL DISPLACEMENT
C ------------------

      DO 10 I=1,3

        DXTIDE(I)=FAC2SUN*( X2SUN*XSUN(I)/RSUN + P2SUN*XSTA(I)/RSTA ) +
     1            FAC2MON*( X2MON*XMON(I)/RMON + P2MON*XSTA(I)/RSTA ) +
     2            FAC3SUN*( X3SUN*XSUN(I)/RSUN + P3SUN*XSTA(I)/RSTA ) +
     3            FAC3MON*( X3MON*XMON(I)/RMON + P3MON*XSTA(I)/RSTA )
c         write (13,*) 'dx1', dxtide(i)
10    CONTINUE

      call zero_vec8(xcorsta)

C
C CORRECTIONS FOR THE OUT-OF-PHASE PART OF LOVE NUMBERS (PART H_2^(0)I
C            AND L_2^(0)I )

C FIRST, FOR THE DIURNAL BAND

      CALL ST1IDIU(XSTA,XSUN,XMON,FAC2SUN,FAC2MON,XCORSTA)

      DO 11 I=1,3
         DXTIDE(I)=DXTIDE(I)+XCORSTA(I)
11    CONTINUE

C
C SECOND, FOR THE SEMI-DIURNAL BAND
C

      CALL ST1ISEM(XSTA,XSUN,XMON,FAC2SUN,FAC2MON,XCORSTA)

      DO 12 I=1,3
         DXTIDE(I)=DXTIDE(I)+XCORSTA(I)
12    CONTINUE

C
C CORRECTIONS FOR THE LATITUDE DEPENDENCE OF LOVE NUMBERS (PART L^(1) )
C

      CALL ST1L1(XSTA,XSUN,XMON,FAC2SUN,FAC2MON,XCORSTA)

      DO 13 I=1,3
         DXTIDE(I)=DXTIDE(I)+XCORSTA(I)
13    CONTINUE

C
C CONSIDER CORRECTIONS FOR STEP 2
C

C CORRECTIONS FOR THE DIURNAL BAND:

C
C  FIRST, WE NEED TO KNOW THE DATE CONVERTED IN JULIAN CENTURIES
c

c   1) CALL THE SUBROUTINE COMPUTING THE JULIAN DATE
c
c      EXPRESSION OF THE HOURS, MINUTES AND SECONDES IN FRACTION OF DAY
c

      djuld=fjldy(iyr,imonth,iday,fhr)


C convert to centuries.

      T=(DJULD-2451545.d0)/36525.d0
c
c   2) CALL THE SUBROUTINE COMPUTING THE CORRECTION OF UTC TIME
c

      CALL DUTC(IYR,IMONTH,IDAY,DTT)

      fhr=fhr+dtt/3600.

C
C  SECOND, WE CAN CALL THE SUBROUTINE STEP2DIU, FOR THE DIURNAL BAND CORRECTIONS,
C   (in-phase and out-of-phase frequency dependence):
C

      CALL STEP2DIU(XSTA,FHR,T,XCORSTA)

      DO 14 I=1,3
         DXTIDE(I)=DXTIDE(I)+XCORSTA(I)
14    CONTINUE

C
C  CORRECTIONS FOR THE LONG-PERIOD BAND,
C   (in-phase and out-of-phase frequency dependence):
C

      call STEP2LON(XSTA,FHR,T,XCORSTA)

      DO 15 I=1,3
         DXTIDE(I)=DXTIDE(I)+XCORSTA(I)

c         write (13,*) 'dx2', dxtide(i)

15    CONTINUE



C
C CONSIDER CORRECTIONS FOR STEP 3
C

C UNCORRECT FOR THE PERMANENT TIDE

C
C      PI=3.141592654
C      SINPHI=XSTA(3)/RSTA

C      COSPHI=dsqrt(XSTA(1)**2+XSTA(2)**2)/RSTA
C      COSLA=XSTA(1)/COSPHI/RSTA
C      SINLA=XSTA(2)/COSPHI/RSTA
C      DR=-DSQRT(5./4./PI)*H2*0.31460*(3./2.*SINPHI**2-0.5)
C      DN=-DSQRT(5./4./PI)*L2*0.31460*3.*COSPHI*SINPHI
C      DXTIDE(1)=DXTIDE(1)-DR*COSLA*COSPHI+DN*COSLA*SINPHI
C      DXTIDE(2)=DXTIDE(2)-DR*SINLA*COSPHI+DN*SINLA*SINPHI
C      DXTIDE(3)=DXTIDE(3)-DR*SINPHI      -DN*COSPHI
C

      RETURN
      END

C  *************************************************************
c
c    COMPUTING OF THE UTC TIME
c
      subroutine dutc(iyear,imonth,iday,dut)

c  Made from the subroutine sent by Dennis McCarthy

      implicit none
      INTEGER iyear,imonth,iday
      INTEGER i,j
      DOUBLE PRECISION dut,an
      double precision utc(2,150)

C     DIFFERENCES UTC-ET
C
C
       data ((utc(i,j), i= 1, 2), j= 1, 98) /

     1 1955.5, 31.59,
     2 1956.5, 32.06,
     3 1957.5, 31.81999,
     4 1958.5, 32.68999,
     5 1959.5, 33.04999,
     6 1960.5, 33.15999,
     7 1961.5, 33.59,
     8 1962.0, 34.032,
     9 1962.5, 34.235,
     x 1963.0, 34.441,
     1 1963.5, 34.644,
     2 1964.0, 34.95,
     3 1964.5, 35.28599,
     4 1965.0, 35.72499,
     5 1965.5, 36.15999,
     6 1966.0, 36.498,
     7 1966.5, 36.96799,
     8 1967.0, 37.444,
     9 1967.5, 37.91299,
     * 1968.0, 38.38999,
     1 1968.25, 38.526,
     2 1968.5, 38.75999,
     3 1968.75, 39.0,
     4 1969.0, 39.23799,
     5 1969.25, 39.472,
     6 1969.5, 39.707,
     7 1969.75, 39.94599,
     8 1970.0, 40.185,
     9 1970.25, 40.41999,
     x 1970.5, 40.65399,
     1 1970.75, 40.89199,
     2 1971.0, 41.131,
     3 1971.085, 41.21099,
     4 1971.162, 41.284,
     5 1971.247, 41.36399,
     6 1971.329, 41.442,
     7 1971.414, 41.52199,
     8 1971.496, 41.59999,
     9 1971.581, 41.68,
     * 1971.666, 41.761,
     1 1971.748, 41.838,
     2 1971.833, 41.91899,
     3 1971.915, 41.99599,
     4 1971.99999, 42.18399,
     5 1972.0, 42.18399,
     6 1972.49999, 42.18399,
     7 1972.5, 43.18399,
     8 1972.99999, 43.18399,
     9 1973.0, 44.18399,
     x 1973.99999, 44.18399,
     1 1974.0, 45.18399,
     2 1974.99999, 45.18399,
     3 1975.0, 46.18399,
     4 1975.99999, 46.18399,
     5 1976.0, 47.18399,
     6 1976.99999, 47.18399,
     7 1977.0, 48.18399,
     8 1977.99999, 48.18399,
     9 1978.0, 49.18399,
     * 1978.99999, 49.18399,
     1 1979.0, 50.18399,
     2 1979.99999, 50.18399,
     3 1980.0, 51.18399,
     4 1981.49999, 51.18399,
     5 1981.5, 52.18399,
     6 1982.49999, 52.18399,
     7 1982.5, 53.18399,
     8 1983.49999, 53.18399,
     9 1983.5, 54.18399,
     x 1985.49999, 54.18399,
     1 1985.5, 55.18399,
     2 1987.99999, 55.18399,
     3 1988.0, 56.18399,
     4 1989.99999, 56.18399,
     5 1990.0, 57.18399,
     6 1990.99999, 57.18399,
     7 1991.0, 58.18399,
     8 1992.49999, 58.18399,
     9 1992.5, 59.18399,
     * 1993.49999, 59.18399,
     1 1993.5, 60.18399,
     2 1994.49999, 60.18399,
     3 1994.5, 61.18399,
     4 1995.99999, 61.18399,
     5 1996.0, 62.18399,
     6 1997.49999, 62.18399,
     7 1997.5, 63.18399,
     8 1998.99999, 63.18399,
     9 1999.0, 64.18399,
     * 2005.99999, 64.18399,
     1 2006.0, 65.18399,
     2 2008.99999, 65.18399,
     3 2009.0, 66.18399,
     4 2012.49999, 66.18399,
     5 2012.5, 67.18399,
     6 2015.49999, 67.18399,
     7 2016.99999, 68.18399,
     8 2017.0,  68.18399 /

      an=iyear+imonth/12d0+iday/365.25d0

      if (an.lt.1955.5) then
       print *,'THE LEAP SECONDS ARE NOT GIVEN BEFORE 1955. '
       PRINT *,'THE CORRECTIONS COMPUTED FOR THE PRIOR YEARS ARE '
       print *,'COMPUTED WITH UTC-ET=31.59 2 SECOND.'
       dut=31.59
       goto 510
      endif

      if (an.gt.2017.0d0) then

       dut=69.184
       goto 510

      endif

      do 5 i=1,97
        if (an.gt.utc(1,i).AND.an.gt.UTC(1,i+1)) then
          dut=utc(2,i)
          goto 510
        endif
   5  continue

 510  return
      end

C********************************************************************

      DOUBLE PRECISION function fjldy(iyr,imon,iday,rhour)

C compute julian day given year,month,day, and real hours
C Uses algorthim of Jean Meeus in "Astronomical Algorithms"
C
C Written by  John Gipson
C    April 8, 1996
C This uses the astronomical convention in which the year which
C historians call -1 BC is represented as 0.

C Should work for all years > -4716.  Checked on the following dates:
C
C   16. April   -1410, 0.00 H is FJLDY = 1206160.5D0    C
C   31. January -1100, 0.00 H is FJLDY = 1319312.5D0    C
C   24. January -0800, 0.00 H is FJLDY = 1428880.5D0    C
C   17. January -0500, 0.00 H is FJLDY = 1538448.5D0    C
C   10. January -0200, 0.00 H is FJLDY = 1648016.5D0    C
C   03. January   100, 0.00 H is FJLDY = 1757584.5D0    C
C   29. February  400, 0.00 H is FJLDY = 1867216.5D0    C
C   20. December  699, 0.00 H is FJLDY = 1976720.5D0    C
C   15. February 1000, 0.00 H is FJLDY = 2086352.5D0    C
C   08. February 1300, 0.00 H is FJLDY = 2195920.5D0    C
C   11. February 1600, 0.00 H is FJLDY = 2305488.5D0    C
C   06. February 1900, 0.00 H is FJLDY = 2415056.5D0    C
C   01. January  1988, 0.00 H is FJLDY = 2447161.5D0    C
C   01. February 1988, 0.00 H is FJLDY = 2447192.5D0    C
C   29. February 1988, 0.00 H is FJLDY = 2447220.5D0    C
C   01. March    1988, 0.00 H is FJLDY = 2447221.5D0    C
C   01. February 2200, 0.00 H is FJLDY = 2524624.5D0    C
C   27. January  2500, 0.00 H is FJLDY = 2634192.5D0    C
C   23. January  2800, 0.00 H is FJLDY = 2743760.5D0    C
C   22. December 3002, 0.00 H is FJLDY = 2817872.5D0    C

      implicit none
      INTEGER iyr,imon,iday
      DOUBLE PRECISION rhour
      INTEGER jyr,jmon

      IF(imon .GT. 2) then
         jmon=imon
         jyr=iyr
      else
         jmon=imon+12
         jyr=iyr-1
      endif

C Factor of .1 in months and years insures that roundoff errors will not be a
C problem. For example 30.6*5=183.  Some machines will calculate this as
C 182.9999..., and truncate to 182.

      fjldy=INT(365.25d0*float(jyr+4716)+.1)
     >    +INT(30.6*float(jmon+1)+0.1) +iday-1524.5+rhour/24.

C
C  Check to see if gregorian or julian date. These handle leap years differently,
C  as well as have an offset. Switchover date was Oct 4 1582. Next day was Oct 15 1582.
C

      IF(iyr*10000.+imon*100.+iday .GT. 15821015.)
     >    fjldy=fjldy+float(2-INT(jyr/100)+INT(INT(jyr/100)/4))

      return
      end

C*********************************************************************************************
C
C
      SUBROUTINE SPROD(X,Y,SCAL,R1,R2)
C
C  COMPUTATION OF THE SCALAR-PRODUCT OF TWO VECTORS AND THEIR NORMS
C
C  INPUT :  X(I),I=1,2,3: COMPONENTS OF VECTOR X
C           Y(I),I=1,2,3: COMPONENTS OF VECTOR Y
C  OUTPUT :  SCAL: SCALAR PRODUCT OF X AND Y
C            R1,R2  : LENGTHS OF THE TWO VECTORS X AND Y
C

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION X(3),Y(3)
      R1=DSQRT(X(1)**2+X(2)**2+X(3)**2)
      R2=DSQRT(Y(1)**2+Y(2)**2+Y(3)**2)
      SCAL=X(1)*Y(1)+X(2)*Y(2)+X(3)*Y(3)
      RETURN
      END

C
C-------------------------------------------------------------------------
C

      SUBROUTINE ST1L1(XSTA,XSUN,XMON,FAC2SUN,FAC2MON,XCORSTA)

C
C THIS SUBROUTINE GIVES THE CORRECTIONS INDUCED BY THE LATITUDE DEPENDENCE
C GIVEN BY L^(1) IN MAHTEWS ET AL (1991)
C
C       INPUT : XSTA,XSUN,XMON,FAC3SUN,FAC3MON
C       OUTPUT : XCORSTA
C

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

C compute the norm of a vector

      DOUBLE PRECISION norm8
      DIMENSION XSTA(3),XSUN(3),XMON(3),XCORSTA(3)
      DOUBLE PRECISION L1,L1D,L1SD

      DATA L1D/0.0012d0/,L1SD/0.0024d0/
      RSTA=norm8(xsta)
      SINPHI=XSTA(3)/RSTA
      COSPHI=DSQRT(XSTA(1)**2+XSTA(2)**2)/RSTA
      SINLA=XSTA(2)/COSPHI/RSTA
      COSLA=XSTA(1)/COSPHI/RSTA
      RMON=norm8(XMON)
      rsun=norm8(xsun)
C
C FOR THE DIURNAL BAND
C

      L1=L1D
      DNSUN=-L1*SINPHI**2*FAC2SUN*XSUN(3)*(XSUN(1)*COSLA+XSUN(2)*SINLA)
     &            /RSUN**2
      DNMON=-L1*SINPHI**2*FAC2MON*XMON(3)*(XMON(1)*COSLA+XMON(2)*SINLA)
     &            /RMON**2
      DESUN=L1*SINPHI*(COSPHI**2-SINPHI**2)*FAC2SUN*XSUN(3)*
     1 (XSUN(1)*SINLA-XSUN(2)*COSLA)/RSUN**2
      DEMON=L1*SINPHI*(COSPHI**2-SINPHI**2)*FAC2MON*XMON(3)*
     1 (XMON(1)*SINLA-XMON(2)*COSLA)/RMON**2
      DE=3.d0*(DESUN+DEMON)
      DN=3.d0*(DNSUN+DNMON)
      XCORSTA(1)=-DE*SINLA-DN*SINPHI*COSLA
      XCORSTA(2)=DE*COSLA-DN*SINPHI*SINLA
      XCORSTA(3)=DN*COSPHI

C
C FOR THE SEMI-DIURNAL BAND
C

      L1=L1SD
      COSTWOLA=COSLA**2-SINLA**2
      SINTWOLA=2.d0*COSLA*SINLA
      DNSUN=-L1/2.d0*SINPHI*COSPHI*FAC2SUN*((XSUN(1)**2-XSUN(2)**2)*
     1 COSTWOLA+2.d0*XSUN(1)*XSUN(2)*SINTWOLA)/RSUN**2
      DNMON=-L1/2.d0*SINPHI*COSPHI*FAC2MON*((XMON(1)**2-XMON(2)**2)*
     1 COSTWOLA+2.d0*XMON(1)*XMON(2)*SINTWOLA)/RMON**2
      DESUN=-L1/2.d0*SINPHI**2*COSPHI*FAC2SUN*((XSUN(1)**2-XSUN(2)**2)*
     1 SINTWOLA-2.d0*XSUN(1)*XSUN(2)*COSTWOLA)/RSUN**2
      DEMON=-L1/2.d0*SINPHI**2*COSPHI*FAC2MON*((XMON(1)**2-XMON(2)**2)*
     1 SINTWOLA-2.d0*XMON(1)*XMON(2)*COSTWOLA)/RMON**2
      DE=3.d0*(DESUN+DEMON)
      DN=3.d0*(DNSUN+DNMON)

      XCORSTA(1)=XCORSTA(1)-DE*SINLA-DN*SINPHI*COSLA
      XCORSTA(2)=XCORSTA(2)+DE*COSLA-DN*SINPHI*SINLA
      XCORSTA(3)=XCORSTA(3)+DN*COSPHI

      RETURN
      END

C*************************************************************************

C     Last change:  VD   17 May 00   1:20 pm
C  THESE ARE THE SUBROUTINES FOR THE STEP2 OF THE TIDAL CORRECTIONS.
C  THEY ARE CALLED TO ACCOUNT FOR THE FREQUENCY DEPENDENCE
C  OF THE LOVE NUMBERS.
C

      SUBROUTINE STEP2DIU(XSTA,FHR,T,XCORSTA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DOUBLE PRECISION XSTA(3),XCORSTA(3),DATDI(9,31)
      DOUBLE PRECISION deg2rad

      DATA deg2rad/0.0174532925d0/
      DATA ((DATDI(i,j),i=1,9),j=1,31)/

     * -3., 0., 2., 0., 0.,-0.01,-0.01, 0.0 , 0.0,
     * -3., 2., 0., 0., 0.,-0.01,-0.01, 0.0 , 0.0,
     * -2., 0., 1.,-1., 0.,-0.02,-0.01, 0.0 , 0.0,
     * -2., 0., 1., 0., 0.,-0.08,-0.05, 0.01,-0.02,
     * -2., 2.,-1., 0., 0.,-0.02,-0.01, 0.0 , 0.0,
     * -1., 0., 0.,-1., 0.,-0.10,-0.05, 0.0 ,-0.02,
     * -1., 0., 0., 0., 0.,-0.51,-0.26,-0.02,-0.12,
     * -1., 2., 0., 0., 0., 0.01, 0.0 , 0.0 , 0.0,
     *  0.,-2., 1., 0., 0., 0.01, 0.0 , 0.0 , 0.0,
     *  0., 0.,-1., 0., 0., 0.02, 0.01, 0.0 , 0.0,
     *  0., 0., 1., 0., 0., 0.06, 0.02, 0.0 , 0.01,
     *  0., 0., 1., 1., 0., 0.01, 0.0 , 0.0 , 0.0,
     *  0., 2.,-1., 0., 0., 0.01, 0.0 , 0.0 , 0.0,
     *  1.,-3., 0., 0., 1.,-0.06, 0.0 , 0.0 , 0.0,
     *  1.,-2., 0., 1., 0., 0.01, 0.0 , 0.0 , 0.0,
     *  1.,-2., 0., 0., 0.,-1.23,-0.05, 0.06,-0.06,
     *  1.,-1., 0., 0.,-1., 0.02, 0.0 , 0.0 , 0.0,
     *  1.,-1., 0., 0., 1., 0.04, 0.0 , 0.0 , 0.0,
     *  1., 0., 0.,-1., 0.,-0.22, 0.01, 0.01, 0.0,
     *  1., 0., 0., 0., 0.,12.02,-0.45,-0.66, 0.17,
     *  1., 0., 0., 1., 0., 1.73,-0.07,-0.10, 0.02,
     *  1., 0., 0., 2., 0.,-0.04, 0.0 , 0.0 , 0.0,
     *  1., 1., 0., 0.,-1.,-0.50, 0.0 , 0.03, 0.0,
     *  1., 1., 0., 0., 1., 0.01, 0.0 , 0.0 , 0.0,
     *  0., 1., 0., 1.,-1.,-0.01, 0.0 , 0.0 , 0.0,
     *  1., 2.,-2., 0., 0.,-0.01, 0.0 , 0.0 , 0.0,
     *  1., 2., 0., 0., 0.,-0.12, 0.01, 0.01, 0.0,
     *  2.,-2., 1., 0., 0.,-0.01, 0.0 , 0.0 , 0.0,
     *  2., 0.,-1., 0., 0.,-0.02, 0.02, 0.0 , 0.01,
     *  3., 0., 0., 0., 0., 0.0 , 0.01, 0.0 , 0.01,
     *  3., 0., 0., 1., 0., 0.0 , 0.01, 0.0 , 0.0/

      S=218.31664563D0+481267.88194D0*T-0.0014663889D0*T**2
     1 +0.00000185139D0*T**3
      TAU=fhr*15.D0+280.4606184D0+36000.7700536D0*T+0.00038793D0*T**2
     1 -0.0000000258D0*T**3-S

      PR=1.396971278d0*T+0.000308889d0*T**2+0.000000021d0*T**3
     1 +0.000000007d0*T**4

      S=S+PR

      H=280.46645D0+36000.7697489D0*T+0.00030322222D0*T**2
     1 +0.000000020d0*T**3-0.00000000654d0*T**4

      P=83.35324312D0+4069.01363525D0*T-0.01032172222D0*T**2
     1 -0.0000124991D0*T**3+0.00000005263D0*T**4

      ZNS=234.95544499D0 +1934.13626197D0*T-0.00207561111D0*T**2
     1 -0.00000213944D0*T**3+0.00000001650D0*T**4

      PS=282.93734098D0+1.71945766667D0*T+0.00045688889D0*T**2
     1 -0.00000001778D0*T**3-0.00000000334D0*T**4

C Reduce angles to between 0 and 360.

      s=  dmod(s,360.d0)
      tau=dmod(tau,360.d0)
      h=  dmod(h,360.d0)
      p=  dmod(p,360.d0)
      zns=dmod(zns,360.d0)
      ps=dmod(ps,360.d0)

c      WRITE(2,'(6f10.3)') tau,s,h,p,zns,ps

      RSTA=DSQRT(XSTA(1)**2+XSTA(2)**2+XSTA(3)**2)
      SINPHI=XSTA(3)/RSTA
      COSPHI=DSQRT(XSTA(1)**2+XSTA(2)**2)/RSTA

      COSLA=XSTA(1)/COSPHI/RSTA
      SINLA=XSTA(2)/COSPHI/RSTA
      ZLA = DATAN2(XSTA(2),XSTA(1))

      DO 99 I=1,3
         XCORSTA(I)=0.d0
99    continue

      DO 98 J=1,31
         THETAF=(TAU+DATDI(1,J)*S+DATDI(2,J)*H+DATDI(3,J)*P+
     1   DATDI(4,J)*ZNS+DATDI(5,J)*PS)*deg2rad
         DR=DATDI(6,J)*2.d0*SINPHI*COSPHI*SIN(THETAF+ZLA)+
     1   DATDI(7,J)*2.d0*SINPHI*COSPHI*COS(THETAF+ZLA)
         DN=DATDI(8,J)*(COSPHI**2-SINPHI**2)*DSIN(THETAF+ZLA)+
     1   DATDI(9,J)*(COSPHI**2-SINPHI**2)*DCOS(THETAF+ZLA)
         DE=DATDI(8,J)*SINPHI*DCOS(THETAF+ZLA)+
     1   DATDI(9,J)*SINPHI*DSIN(THETAF+ZLA)
         XCORSTA(1)=XCORSTA(1)+DR*COSLA*COSPHI-DE*SINLA
     1   -DN*SINPHI*COSLA
         XCORSTA(2)=XCORSTA(2)+DR*SINLA*COSPHI+DE*COSLA
     1   -DN*SINPHI*SINLA
         XCORSTA(3)=XCORSTA(3)+DR*SINPHI+DN*COSPHI

98    CONTINUE

      DO 97 I=1,3
         XCORSTA(I)=XCORSTA(I)/1000.d0
97    CONTINUE

      RETURN
      END

C
C  *************************************************************
C

      SUBROUTINE STEP2LON(XSTA,FHR,T,XCORSTA)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DOUBLE PRECISION XSTA(3),XCORSTA(3),DATDI(9,5)
      DOUBLE PRECISION deg2rad

      DATA deg2rad/0.0174532925d0/

      DATA ((DATDI(i,j),i=1,9),j=1,5)/

     *   0, 0, 0, 1, 0,   0.47, 0.23, 0.16, 0.07,
     *   0, 2, 0, 0, 0,  -0.20,-0.12,-0.11,-0.05,
     *   1, 0,-1, 0, 0,  -0.11,-0.08,-0.09,-0.04,
     *   2, 0, 0, 0, 0,  -0.13,-0.11,-0.15,-0.07,
     *   2, 0, 0, 1, 0,  -0.05,-0.05,-0.06,-0.03/

C

      S=218.31664563D0+481267.88194D0*T-0.0014663889D0*T**2
     1 +0.00000185139D0*T**3

      PR=1.396971278d0*T+0.000308889d0*T**2+0.000000021d0*T**3
     1 +0.000000007d0*T**4

      S=S+PR

      H=280.46645D0+36000.7697489D0*T+0.00030322222D0*T**2
     1 +0.000000020d0*T**3-0.00000000654d0*T**4

      P=83.35324312D0+4069.01363525D0*T-0.01032172222D0*T**2
     1 -0.0000124991D0*T**3+0.00000005263D0*T**4

      ZNS=234.95544499D0 +1934.13626197D0*T-0.00207561111D0*T**2
     1 -0.00000213944D0*T**3+0.00000001650D0*T**4

      PS=282.93734098D0+1.71945766667D0*T+0.00045688889D0*T**2
     1 -0.00000001778D0*T**3-0.00000000334D0*T**4

      RSTA=DSQRT(XSTA(1)**2+XSTA(2)**2+XSTA(3)**2)
      SINPHI=XSTA(3)/RSTA
      COSPHI=DSQRT(XSTA(1)**2+XSTA(2)**2)/RSTA
      COSLA=XSTA(1)/COSPHI/RSTA
      SINLA=XSTA(2)/COSPHI/RSTA

C reduce angles to between 0 and 360.

      s=  dmod(s,360.d0)
      tau=dmod(tau,360.d0)
      h=  dmod(h,360.d0)
      p=  dmod(p,360.d0)
      zns=dmod(zns,360.d0)
      ps=dmod(ps,360.d0)

      dr_tot=0.
      dn_tot=0.

      DO 99 I=1,3
         XCORSTA(I)=0.
99    continue

      DO 98 J=1,5

         THETAF=(DATDI(1,J)*S+DATDI(2,J)*H+DATDI(3,J)*P+
     1   DATDI(4,J)*ZNS+DATDI(5,J)*PS)*DEG2RAD
         DR=DATDI(6,J)*(3.*SINPHI**2-1.)/2.*DCOS(THETAF)+
     1   DATDI(8,J)*(3.*SINPHI**2-1.)/2.*DSIN(THETAF)
         DN=DATDI(7,J)*(COSPHI*SINPHI*2.d0)*DCOS(THETAF)+
     1   DATDI(9,J)*(COSPHI*SINPHI*2.d0)*DSIN(THETAF)
         DE=0.d0
         dr_tot=dr_tot+dr
         dn_tot=dn_tot+dn
         XCORSTA(1)=XCORSTA(1)+DR*COSLA*COSPHI-DE*SINLA
     1   -DN*SINPHI*COSLA
         XCORSTA(2)=XCORSTA(2)+DR*SINLA*COSPHI+DE*COSLA
     1   -DN*SINPHI*SINLA
          XCORSTA(3)=XCORSTA(3)+DR*SINPHI+DN*COSPHI
98    CONTINUE

      DO 97 I=1,3
         XCORSTA(I)=XCORSTA(I)/1000.d0
97    CONTINUE

      RETURN

      END

C**************************************************************************************************
C

      SUBROUTINE STEP2IDIU(XSTA,FHR,T,XCORSTA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DOUBLE PRECISION XSTA(3),XCORSTA(3),DATDI(9,15)
      DOUBLE PRECISION deg2rad
      DATA deg2rad/0.0174532925d0/
      DATA ((DATDI(i,j),i=1,9),j=1,15)/
     * -3., 0., 2., 0., 0.,-0.01,-0.01, 0.0 , 0.0,
     * -3., 2., 0., 0., 0.,-0.01,-0.01, 0.0 , 0.0,
     * -2., 0., 1.,-1., 0.,-0.02,-0.01, 0.0 , 0.0,
     * -2., 0., 1., 0., 0.,-0.09,-0.05, 0.01,-0.02,
     * -2., 2.,-1., 0., 0.,-0.02,-0.01, 0.0 , 0.0,
     * -1., 0., 0.,-1., 0.,-0.10,-0.05, 0.01,-0.02,
     * -1., 0., 0., 0., 0.,-0.53,-0.26, 0.06,-0.12,
     *  0., 0., 1., 0., 0., 0.06, 0.02,-0.01, 0.01,
     *  1.,-2., 0., 0., 0.,-1.17,-0.05, 0.08,-0.06,
     *  1., 0., 0.,-1., 0.,-0.21, 0.01, 0.01, 0.0,
     *  1., 0., 0., 0., 0.,11.29,-0.44,-0.70, 0.17,
     *  1., 0., 0., 1., 0., 1.63,-0.07,-0.10, 0.02,
     *  1., 2., 0., 0., 0.,-0.11, 0.01, 0.01, 0.0,
     *  2., 0.,-1., 0., 0.,-0.02, 0.02, 0.0 , 0.01,
     *  3., 0., 0., 0., 0., 0.0 , 0.01, 0.0 , 0.01/

      S=218.31664563D0+481267.88194D0*T-0.0014663889D0*T**2
     1 +0.00000185139D0*T**3

      TAU=fhr*15.D0+280.4606184D0+36000.7700536D0*T+0.00038793D0*T**2
     1 -0.0000000258D0*T**3-S

      PR=1.396971278d0*T+0.000308889d0*T**2+0.000000021d0*T**3
     1 +0.000000007d0*T**4

      S=S+PR
      H=280.46645D0+36000.7697489D0*T+0.00030322222D0*T**2
     1 +0.000000020d0*T**3-0.00000000654d0*T**4
      P=83.35324312D0+4069.01363525D0*T-0.01032172222D0*T**2
     1 -0.0000124991D0*T**3+0.00000005263D0*T**4
      ZNS=234.95544499D0 +1934.13626197D0*T-0.00207561111D0*T**2
     1 -0.00000213944D0*T**3+0.00000001650D0*T**4
      PS=282.93734098D0+1.71945766667D0*T+0.00045688889D0*T**2
     1 -0.00000001778D0*T**3-0.00000000334D0*T**4

C Reduce angles to between 0 and 360.

      s=  dmod(s,360.d0)
      tau=dmod(tau,360.d0)
      h=  dmod(h,360.d0)
      p=  dmod(p,360.d0)
      zns=dmod(zns,360.d0)
      ps=dmod(ps,360.d0)

c      WRITE(2,'(6f10.3)') tau,s,h,p,zns,ps

      RSTA=DSQRT(XSTA(1)**2+XSTA(2)**2+XSTA(3)**2)
      SINPHI=XSTA(3)/RSTA
      COSPHI=DSQRT(XSTA(1)**2+XSTA(2)**2)/RSTA

      COSLA=XSTA(1)/COSPHI/RSTA
      SINLA=XSTA(2)/COSPHI/RSTA
      ZLA = DATAN2(XSTA(2),XSTA(1))

      DO 99 I=1,3
         XCORSTA(I)=0.d0
99    continue

      DO 98 J=1,11
         THETAF=(TAU+DATDI(1,J)*S+DATDI(2,J)*H+DATDI(3,J)*P+
     1   DATDI(4,J)*ZNS+DATDI(5,J)*PS)*deg2rad
         DR=DATDI(7,J)*2.*SINPHI*COSPHI*DCOS(THETAF+ZLA)
         DN=DATDI(9,J)*(COSPHI**2-SINPHI**2)*DCOS(THETAF+ZLA)
         DE=-DATDI(9,J)*SINPHI*DSIN(THETAF+ZLA)
         XCORSTA(1)=XCORSTA(1)+DR*COSLA*COSPHI-DE*SINLA
     1   -DN*SINPHI*COSLA
         XCORSTA(2)=XCORSTA(2)+DR*SINLA*COSPHI+DE*COSLA
     1   -DN*SINPHI*SINLA
         XCORSTA(3)=XCORSTA(3)+DR*SINPHI+DN*COSPHI
98    CONTINUE

      DO 97 I=1,3
         XCORSTA(I)=XCORSTA(I)/1000.d0
97    CONTINUE

      RETURN
      END

C
C  *************************************************************
C
C-------------------------------------------------------------------------

      SUBROUTINE ST1IDIU(XSTA,XSUN,XMON,FAC2SUN,FAC2MON,XCORSTA)

C
C THIS SUBROUTINE GIVES THE OUT-OF-PHASE CORRECTIONS INDUCED BY
C MANTLE INELASTICITY IN THE DIURNAL BAND
C
C       INPUT : XSTA,XSUN,XMON,FAC2SUN,FAC2MON
C      OUTPUT : XCORSTA
C

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION norm8
      DIMENSION XSTA(3),XSUN(3),XMON(3),XCORSTA(3)
      DATA DHI/-0.0025d0/,DLI/-0.0007d0/
      RSTA=NORM8(XSTA)
      SINPHI=XSTA(3)/RSTA
      COSPHI=DSQRT(XSTA(1)**2+XSTA(2)**2)/RSTA
      COS2PHI=COSPHI**2-SINPHI**2
      SINLA=XSTA(2)/COSPHI/RSTA
      COSLA=XSTA(1)/COSPHI/RSTA
      RMON=NORM8(XMON)
      RSUN=NORM8(XSUN)
      DRSUN=-3.d0*DHI*SINPHI*COSPHI*FAC2SUN*XSUN(3)*(XSUN(1)*
     1            SINLA-XSUN(2)*COSLA)/RSUN**2
      DRMON=-3.d0*DHI*SINPHI*COSPHI*FAC2MON*XMON(3)*(XMON(1)*
     1            SINLA-XMON(2)*COSLA)/RMON**2
      DNSUN=-3.d0*DLI*COS2PHI*FAC2SUN*XSUN(3)*(XSUN(1)*SINLA-
     1            XSUN(2)*COSLA)/RSUN**2
      DNMON=-3.d0*DLI*COS2PHI*FAC2MON*XMON(3)*(XMON(1)*SINLA-
     1            XMON(2)*COSLA)/RMON**2
      DESUN=-3.d0*DLI*SINPHI*FAC2SUN*XSUN(3)*
     1            (XSUN(1)*COSLA+XSUN(2)*SINLA)/RSUN**2
      DEMON=-3.d0*DLI*SINPHI*FAC2MON*XMON(3)*
     1            (XMON(1)*COSLA+XMON(2)*SINLA)/RMON**2

      DR=DRSUN+DRMON
      DN=DNSUN+DNMON
      DE=DESUN+DEMON

      XCORSTA(1)=DR*COSLA*COSPHI-DE*SINLA-DN*SINPHI*COSLA
      XCORSTA(2)=DR*SINLA*COSPHI+DE*COSLA-DN*SINPHI*SINLA
      XCORSTA(3)=DR*SINPHI+DN*COSPHI

      RETURN
      END
C
C-------------------------------------------------------------------------

      SUBROUTINE ST1ISEM(XSTA,XSUN,XMON,FAC2SUN,FAC2MON,XCORSTA)

C
C THIS SUBROUTINE GIVES THE OUT-OF-PHASE CORRECTIONS INDUCED BY
C MANTLE INELASTICITY IN THE DIURNAL BAND
C
C       INPUT : XSTA,XSUN,XMON,FAC2SUN,FAC2MON
C      OUTPUT : XCORSTA
C

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION norm8
      DIMENSION XSTA(3),XSUN(3),XMON(3),XCORSTA(3)
      DATA DHI/-0.0022/,DLI/-0.0007/
      RSTA=NORM8(XSTA)
      SINPHI=XSTA(3)/RSTA
      COSPHI=DSQRT(XSTA(1)**2+XSTA(2)**2)/RSTA
      SINLA=XSTA(2)/COSPHI/RSTA
      COSLA=XSTA(1)/COSPHI/RSTA
      COSTWOLA=COSLA**2-SINLA**2
      SINTWOLA=2.d0*COSLA*SINLA
      RMON=NORM8(XMON)
      RSUN=NORM8(XSUN)

      DRSUN=-0.75d0*DHI*COSPHI**2*FAC2SUN*((XSUN(1)**2-XSUN(2)**2)*
     1 SINTWOLA-2.d0*XSUN(1)*XSUN(2)*COSTWOLA)/RSUN**2
      DRMON=-0.75d0*DHI*COSPHI**2*FAC2MON*((XMON(1)**2-XMON(2)**2)*
     1 SINTWOLA-2.d0*XMON(1)*XMON(2)*COSTWOLA)/RMON**2
      DNSUN=1.5d0*DLI*SINPHI*COSPHI*FAC2SUN*((XSUN(1)**2-XSUN(2)**2)*
     1 SINTWOLA-2.d0*XSUN(1)*XSUN(2)*COSTWOLA)/RSUN**2
      DNMON=1.5d0*DLI*SINPHI*COSPHI*FAC2MON*((XMON(1)**2-XMON(2)**2)*
     1 SINTWOLA-2.d0*XMON(1)*XMON(2)*COSTWOLA)/RMON**2
      DESUN=-1.5d0*DLI*COSPHI*FAC2SUN*((XSUN(1)**2-XSUN(2)**2)*
     1 COSTWOLA+2.d0*XSUN(1)*XSUN(2)*SINTWOLA)/RSUN**2
      DEMON=-1.5d0*DLI*COSPHI*FAC2MON*((XMON(1)**2-XMON(2)**2)*
     1 COSTWOLA+2.d0*XMON(1)*XMON(2)*SINTWOLA)/RMON**2

      DR=DRSUN+DRMON
      DN=DNSUN+DNMON
      DE=DESUN+DEMON

      XCORSTA(1)=DR*COSLA*COSPHI-DE*SINLA-DN*SINPHI*COSLA
      XCORSTA(2)=DR*SINLA*COSPHI+DE*COSLA-DN*SINPHI*SINLA
      XCORSTA(3)=DR*SINPHI+DN*COSPHI

      RETURN
      END

C
C*********************************************************************

      DOUBLE PRECISION function norm8(a)
      DOUBLE PRECISION a(3)
      norm8=dSQRT(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))
      return
      end

C
C*********************************************************************

      subroutine zero_vec8(v)
      DOUBLE PRECISION v(3)
      v(1)=0.d0
      v(2)=0.d0
      v(3)=0.d0
      return
      end



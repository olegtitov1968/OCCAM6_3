

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
C
C   PN IS PART OF THE PACKAGE OCCAM 
C
C   THIS PROGRAM COMPUTES THE INFORMATION RELATED TO TIME OF OBSERVATION
C  LIKE PRECESION AND NUTATION, EARTH ROTATION PARAMETERS, GREENWICH
C  SIDEREAL TIME, ETC.
C   IT ALSO COMPUTES THE PARTIAL DERIVATIVES OF THE SOURCE POSITION
C  WITH RESPECT TO THE NUTATION ANGLES EPSILON AND PSI (DKDEPS, DKDPSI)
C
C   IT STORES ALL THIS INFORMATION ON FILE SORTIM FOR FURTHER USE
C
C   INCORPORATION OF THE NON-ROTATING ORIGIN EXPRESSION IS PENDING
C
C   REVISION  1991 OCTOBER 19  (N.Z.L.)
C   Revision: 1997 April, 23 (OT). DELE403 does not have
C                  accelerations of the Earth. But due to old format
C                  of reading the corresponding values will be equal to 0.
C   Revision: 1997 June, 25 (OT). Variations in Greenwich Sidereal
C                  Time due to ascending node of the lunar orbit
C                  (Aoki, Kinoshita; 1983).
C   Revision: 1999 February,10 (OT). HERRING MODEL FOR NUTATION
C   Revision: 2000  December,24 (Tesmer)
C                  correction for HSG (~l. 445) is clear now.
C                  (to have C04 in the ICRF and ITRF system if C04 series
C                   computed by IERS after 97.0 in use)
C                  (see IERS Gazettes 03, 08 and 16).
C   Revision: 2001  January,21 (OT)   UT1INT has been deleted
C   Revision: 2001  January,24 (OT)   Eanes model has been implemented
C                  for subdiurnal EOP variations instead of Gipson one
C   Revision: 2001  MARCH, 11 (OT)   COMMON-BLOCK /LOG/ deleted
C                                         MBH 2000 nutation model added
C   REVISION 2002 June 19, (OT) HSG correction added to subroutine NUTPAR
C                                   Common-block /H/ deleted
C   Revision:  2005  December, 13 (OT)  Transition to the mean system of date
C                  Herring's nutation model is off
C   Revision:  2006  January, 04 (OT)  Direct call of DELE405 is added.
C                                      Options for ephemerides approach is added
C                                      Jupiter is added to /BAR/
C   Revision:  2006 February, 22 (OT)  Linear interpolation is added
C
C   Revision:  2006  January, 04 (OT)  DELE405 declared as logical
C   REVISION 2006 October, 03 (O.T.)  Jupiter is written into SORTIM
C

C   LAST REVISION 2021, May 01 (O.T.)
C
C   Copyright (C) 2001-2021 Oleg Titov, Geoscience Australia
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
C**************************************************************************

      IMPLICIT DOUBLE PRECISION  (A-H, O-Z) 
      CHARACTER  SOURCE*8, QQUIT*1, C1*1, C2*1
      DOUBLE PRECISION MJD
      LOGICAL WRT,DOSHRT,UT1S, dtidal, linear, dele405
      DIMENSION RQU(3),RQP1(3),RQM1(3), PRB(3), PRA(3), PR(3), VR(3)
      DIMENSION XJ(3),YJ(3),DKDEPS(3),DKDPSI(3), VRA(3), VRB(3)
      DIMENSION DKDPSI_v(3), DKDEPS_v(3)
c
c   24.10.2005 (OT) - matrices for rotation due to nutation
c
      DIMENSION SNB(3,3), SNA(3,3), SN(3,3)
      DIMENSION PRECB(3,3), PRECA(3,3), PREC(3,3)
      DIMENSION OFFNB(3), OFFNA(3), OFFN(3), pn(3,3)
      DIMENSION pm(2)  !  (OT) 01.09.2006

C  2004 May 12 (jboehm)
      dimension xp_linear(2), yp_linear(2), ut1kor_linear(2)

      COMMON /BAR/ RBAR(3,3),VBAR(3,3),DCORG(3,3),DVELG(3,3),CMOON(3,3)
     * ,XJUP(3,3)                                          ! 04.01.2006
      COMMON /G/  UTTJ, XP(30), YP(30), MJD(30), UT1KOR(30),
     1              NPAR, IP, JQA, IDELT
      COMMON /BARIC/ XBRC(30),YBRC(30),ZBRC(30),XBRV(30),YBRV(30),
     1 ZBRV(30), BCTIM(30),IP1
      COMMON  /STA/ R1(3), R2(3), B(3), V2(3), A2(3)
      COMMON /MOON/ TMOON(30),RAM(6,30),DECM(6,30),PARM(6,30),
     *              XMON(30), YMON(30), ZMON(30), IPM
      COMMON /PHYS/ C, FL, AH, AU, J20
      COMMON /MATH/ PI, TWOPI, ROG, ROS, ROZ, TOL
c
c  Change  on 3.11.2005
c
c      COMMON /PNMAT/ E(3,3)

C------------------------
C   GET PARAMETERS
C------------------------

      CALL PARAM

C   INITIALIZATION

      DOSHRT=.TRUE.
      DTIDAL=.TRUE.
      UT1S = .FALSE.
      WRT=.TRUE.
      DELE405 = .TRUE.

      open (11, file='precession.dat',access='append', status='unknown')
      open (12, file = 'solution.dat', status='unknown')

C
C  CHECK IF DEFAULT MODELS FOR THE SHORT PERIOD UT1 VARIATIONS,
C  AND ANNUAL ABERRATION MUST BE KEPT OR CHANGED
C

524   WRITE (*,521)
521   FORMAT (//' DEFAULT MODELS WILL BE APPLIED FOR:',
     * /'  A.- SHORT PERIOD UT1 VARIATIONS '
     * /'  B.- EFFECT OF OCEANIC TIDES ON UT1 AND POLAR MOTION'//
     * ' DO YOU WANT TO CHANGE THIS ? (Y,(N),Q TO QUIT)')
      READ (*,458) QQUIT
458   FORMAT(A1)

      IF (QQUIT.EQ.'Q'.OR.QQUIT.EQ.'q') STOP
      IF (QQUIT.EQ.'y'.OR.QQUIT.EQ.'Y') THEN

         WRITE (*,522)
522      FORMAT (/' ENTER WHICH EFFECT YOU WANT TO NEGLECT'/' - - ')

         READ (*,'(A1,1X,A1,1X)') C1,C2

         WRITE (*,523) C1,C2
523      FORMAT (' NEGLECTING ',4(A1,1X))
         PRINT *, ' ARE YOU SURE ?  ((Y),N,Q)'

         READ (*,458) QQUIT

         IF (QQUIT.EQ.'Q'.OR.QQUIT.EQ.'q') STOP
         IF (QQUIT.EQ.'N'.OR.QQUIT.EQ.'n') GO TO 524

         IF (C1.EQ.'A'.OR.C2.EQ.'A') then
            DOSHRT = .FALSE.
            write (12,496)
         end if
         IF (C1.EQ.'a'.OR.C2.EQ.'a') then
            DOSHRT = .FALSE.
            write (12,496)
         end if
         IF (C1.EQ.'B'.OR.C2.EQ.'B') then
            DTIDAL = .FALSE.
            write (12,497)
         end if
         IF (C1.EQ.'b'.OR.C2.EQ.'b') then
            DTIDAL = .FALSE.
            write (12,497)
         end if
      else if (QQUIT.EQ.'N'.OR.QQUIT.EQ.'n') then
         write (12,495)
      ENDIF

 495  format (1x,'Short Period UT1 variations and oceanic tides are inc
     *luded')
 496  format (1x,'Short Period UT1 variations are NOT included')
 497  format (1x,'Oceanic tides are NOT included')

C------------------------------------------------------------------------
C  READING EARTH ORIENTATION PARAMETERS AND READING EPHEMERIS FOR EARTH,
C  SUN AND MOON (TABLE EXTRACTED FROM IERS C 04)
C------------------------------------------------------------------------

      CALL EPHEM

C  OPEN AND READ SORTIM

      CALL UNFOLD (19,'SORTIM')
      read (19,rec=1) nrec

      IF (DTIDAL) THEN

          MTIDAL = 2

          WRITE (*,621)
621       FORMAT (//' CHOOSE MODEL TO BE APPLIED FOR: ' /
     * '  A.- RAY MODEL - IERS 96 CONVENTIONS                ' /
     * '  B.- EANES MODEL - IERS 2000 CONVENTIONS - (ON DEFAULT)' //
     *   ' SELECT YOUR PREFERENCE  (A, (B), Q TO QUIT)')
          READ (*,458) QQUIT
          IF (QQUIT.EQ.'A'.OR.QQUIT.EQ.'a') then
             MTIDAL = 1
             write (12,500)
          else if (QQUIT.eq.'B'.OR.QQUIT.EQ.'b') then
             write (12,501)
          end if
          IF (QQUIT.EQ.'Q'.OR.QQUIT.EQ.'q') STOP

      END IF
 
 500  format (1x, ' Model of subdiurnal EOP variations by Ray ')
 501  format (1x, ' Model of subdiurnal EOP variations by Eanes ')

C
C  CHECK IF YOU WANT TO APPLY FOR NUTATION MODEL
C
      NUT = 1

      WRITE (*,622)
622   FORMAT (/' DO YOU WANT TO APPLY FOR                      '/
     * '           A - IAU 1980 NUT. MODEL?      (ON DEFAULT)  '/
     * '           B - MBH_2000 MODEL?                         '//
     *   ' SELECT YOUR PREFERENCE  ((A), B, C, Q TO QUIT)')
      READ (*,458) QQUIT
      IF (QQUIT.EQ.'B'.OR.QQUIT.EQ.'b') then
         NUT = 2
         write (12,503)
      else  IF (QQUIT.EQ.'A'.OR.QQUIT.EQ.'a') then
         write (12,504)
      end if
      IF (QQUIT.EQ.'Q'.OR.QQUIT.EQ.'q') STOP

 503  format (1x, ' MBH (2000) nutation model ')
 504  format (1x, ' IAU 1980 nutation model ')
      write (12,*)


      IF (DOSHRT) THEN
         WRITE (*,531)
531      FORMAT (//' FOR THE TIDAL UT1 VARIATIONS YOU CAN',
     *   ' CHOOSE BETWEEN:'//
     * '  A.- IERS-STANDARDS (1989) (UT1-UT1R SHORT PERIOD < 35 DAYS)'/
     * '  B.- IERS-CONV. (1996) (UT1-UT1S SHORT AND LONG PERIODS)'//
     *   ' SELECT YOUR PREFERENCE  ((A), B, Q TO QUIT)')
         READ (*,458) QQUIT

         IF (QQUIT.EQ.'Q'.OR.QQUIT.EQ.'q') STOP
         IF (QQUIT.EQ.'B'.OR.QQUIT.EQ.'b') UT1S = .TRUE.
      ENDIF

c
c   Selection of the approach for the ephemerides DELE405
c   calculation  (04.01.2006)
c

      M_EPH = 1

          WRITE (*,623)
623       FORMAT (//' CHOOSE APPROACH TO BE APPLIED FOR EPHEMERIDES: ' /
     * '  A.- INTERPOLATION BY 9 POINTS                           ' /
     * '  B.- DIRECT CALCULATIONS (REQUIRES DELE405 EPHEMERIS!! ON DEFAU
     *LT)' //
     * ' SELECT YOUR PREFERENCE  (A, (B), Q TO QUIT)')
          READ (*,458) QQUIT
          IF (QQUIT.EQ.'A'.OR.QQUIT.EQ.'a') then
             M_EPH = 2
             write (12,700)
          else if (QQUIT.eq.'B'.OR.QQUIT.EQ.'b') then
             write (12,701)
          end if
          IF (QQUIT.EQ.'Q'.OR.QQUIT.EQ.'q') STOP

 700  format (1x,'Ephemerides are calculated by 9-point interpolation ')
 701  format (1x,'Ephemerides are calculated directly from the DELE405')


C  2004 May 12 (jboehm)
C
C  ask the user for the kind of interpolation
      linear = .FALSE.
      WRITE (*,132)
132   FORMAT (/' WHICH KIND OF INTERPOLATION FOR EOPS?         '/
     * '           A - LAGRANGE      (ON DEFAULT)              '/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCC TESMER 20.08.04 clearer way to find the a priori EOP lines CCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     * '           B - LINEAR INTERPOLATION                    '/
     * '           B - LINEAR INTERPOLATION                    '//
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCC TESMER 20.08.04 clearer way to find the a priori EOP lines CCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     *   ' SELECT YOUR PREFERENCE  ((A), B, Q TO QUIT)')
      READ (*,458) QQUIT
      IF (QQUIT.EQ.'B'.OR.QQUIT.EQ.'b') linear = .TRUE.
      IF (QQUIT.EQ.'Q'.OR.QQUIT.EQ.'q') STOP
C
C  get julian date for first and last observation

      read (19,rec=1) j, source, ra2000, de2000, Idatj, UTJ
      tmjd1 = Idatj + UTJ/TWOPI

      read (19,rec=j) j, source, ra2000, de2000, Idatj, UTJ
      tmjd2 = Idatj + UTJ/TWOPI

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCC TESMER 20.08.04 clearer way to find the a priori EOP lines CCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  ask the user how to compute the rate
      print*, ' '
      WRITE (*,133)
133   FORMAT (/' HOW TO COMPUTE THE A PRIORI EOP-LINE?     '/
     * '         A - "SOMEHOW" A LAGRANGE WAY          '/
     * '         B - OFFSET: A PRIORI VALUE NEAREST TO MID OF SESSION,'/
     * '             RATE:   WEIGHTED MEAN OF THE TWO A PRIORI RATES'/
     * '                     (FROM MID VALUE TO VALUES BEFORE/AFTER)'/
     * '                     -> this is standard for LSM-approach <-'/
     * '         C - OFFSET: WEIGHTED MEAN OF THREE NEAREST A PRIORIS,'/
     * '             RATE:   WEIGHTED MEAN OF THE TWO A PRIORI RATES'/
     * '                     (FROM MID VALUE TO VALUES BEFORE/AFTER)'//
     * ' SELECT YOUR PREFERENCE  ((A), B, C, Q TO QUIT)')
      READ (*,458) QQUIT
      IF (QQUIT.EQ.'Q'.OR.QQUIT.EQ.'q') STOP

      IF (QQUIT.EQ.'B'.OR.QQUIT.EQ.'b'.OR.
     *    QQUIT.EQ.'C'.OR.QQUIT.EQ.'c') THEN


        if (int(tmjd1).eq.int(tmjd2)) then
c IN CASE OF SESSION IS COMPLETELY BETWEEN TWO 0h VALUES, TAKE IT EASY
c IN CASE OF SESSION IS COMPLETELY BETWEEN TWO 0h VALUES, TAKE IT EASY
c IN CASE OF SESSION IS COMPLETELY BETWEEN TWO 0h VALUES, TAKE IT EASY
          print*, 'session is completely between two 0h values'
          print*, '=> linear interpolation is done directly using'
          print*, '   the two enclosing a priori 0h values !!!'

          tmjd_less=dble(int(tmjd1))
          tmjd_more=dble(int(tmjd1))+1.d0

c searching for the time of offset
          midmjd=nint((tmjd1+tmjd2)/2.d0)
          tmjd_offs=dble(midmjd)

          do iii=1,IP
            if(MJD(iii).eq.tmjd_less) then
              xp_less=XP(iii)
              yp_less=YP(iii)
              ut1kor_less=UT1KOR(iii)
            end if
            if(MJD(iii).eq.tmjd_more) then
              xp_more=XP(iii)
              yp_more=YP(iii)
              ut1kor_more=UT1KOR(iii)
            end if
            if(MJD(iii).eq.tmjd_offs) then
              xp_offs=XP(iii)
              yp_offs=YP(iii)
              ut1kor_offs=UT1KOR(iii)
            end if
          end do

          xrat=(xp_more-xp_less)
          yrat=(yp_more-yp_less)
          ut1korrat=(ut1kor_more-ut1kor_less)

c values are exactly a line between a prioris
          xp_linear(1)=xp_offs+xrat*(tmjd1-tmjd_offs)
          xp_linear(2)=xp_offs+xrat*(tmjd2-tmjd_offs)
          yp_linear(1)=yp_offs+yrat*(tmjd1-tmjd_offs)
          yp_linear(2)=yp_offs+yrat*(tmjd2-tmjd_offs)

          ut1kor_linear(1)=ut1kor_offs+ut1korrat*(tmjd1-tmjd_offs)
          ut1kor_linear(2)=ut1kor_offs+ut1korrat*(tmjd2-tmjd_offs)


        else

c IN CASE OF SESSION IS NOT COMPLETELY BETWEEN TWO 0h VALUES !!!
c IN CASE OF SESSION IS NOT COMPLETELY BETWEEN TWO 0h VALUES !!!
c IN CASE OF SESSION IS NOT COMPLETELY BETWEEN TWO 0h VALUES !!!

c get the three 0h apriori values before, during and after session
          midmjd=nint((tmjd1+tmjd2)/2.d0)
          tmjd_mids=dble(midmjd)
          tmjd_less=tmjd_mids-1.d0
          tmjd_more=tmjd_mids+1.d0

          do iii=1,IP
            if(MJD(iii).eq.tmjd_less) then
              xp_less=XP(iii)
              yp_less=YP(iii)
              ut1kor_less=UT1KOR(iii)
            end if
            if(MJD(iii).eq.tmjd_mids) then
              xp_mids=XP(iii)
              yp_mids=YP(iii)
              ut1kor_mids=UT1KOR(iii)
            end if
            if(MJD(iii).eq.tmjd_more) then
              xp_more=XP(iii)
              yp_more=YP(iii)
              ut1kor_more=UT1KOR(iii)
            end if
          end do

c the rate is the weighted mean of the two rates,
c to have the determined rate in the unit [1/days]
c instead of the unit [1/duration of session] !!!
          weight_less=abs(tmjd_mids-tmjd1)
          weight_more=abs(tmjd2-tmjd_mids)
          wfac=weight_less+weight_more
          weight_less=weight_less/wfac
          weight_more=weight_more/wfac

          xval_less=(xp_mids-xp_less)*weight_less
          xval_more=(xp_more-xp_mids)*weight_more
          xrat=xval_less+xval_more

          yval_less=(yp_mids-yp_less)*weight_less
          yval_more=(yp_more-yp_mids)*weight_more
          yrat=yval_less+yval_more

          ut1korval_less=(ut1kor_mids-ut1kor_less)*weight_less
          ut1korval_more=(ut1kor_more-ut1kor_mids)*weight_more
          ut1korrat=ut1korval_less+ut1korval_more

c the offset of the a priori line is defined by the 0h value in the mid
C CASE B CASE B CASE B CASE B CASE B CASE B CASE B CASE B CASE B CASE B
          xp_linear(1)=xp_mids+xrat*(tmjd1-tmjd_mids)
          xp_linear(2)=xp_mids+xrat*(tmjd2-tmjd_mids)
          yp_linear(1)=yp_mids+yrat*(tmjd1-tmjd_mids)
          yp_linear(2)=yp_mids+yrat*(tmjd2-tmjd_mids)
          ut1kor_linear(1)=ut1kor_mids+ut1korrat*(tmjd1-tmjd_mids)
          ut1kor_linear(2)=ut1kor_mids+ut1korrat*(tmjd2-tmjd_mids)


          IF (QQUIT.EQ.'C'.OR.QQUIT.EQ.'c') THEN

c the offset is the weighted mean of the three values,
c to have the determined rate in the unit [days]
c instead of the unit [duration of session] !!!
C CASE C CASE C CASE C CASE C CASE C CASE C CASE C CASE C CASE C CASE C
            tmjdavg=(tmjd1+tmjd2)/2.d0
            weight_less=1.d0/abs(tmjd_less-tmjdavg)
            weight_mids=1.d0/abs(tmjd_mids-tmjdavg)
            weight_more=1.d0/abs(tmjd_more-tmjdavg)

            wfac=weight_less+weight_mids+weight_more
            weight_less=weight_less/wfac
            weight_mids=weight_mids/wfac
            weight_more=weight_more/wfac

            xp_offs=xp_less*weight_less
     *             +xp_mids*weight_mids
     *             +xp_more*weight_more
            yp_offs=yp_less*weight_less
     *             +yp_mids*weight_mids
     *             +yp_more*weight_more
            ut1kor_offs=ut1kor_less*weight_less
     *             +ut1kor_mids*weight_mids
     *             +ut1kor_more*weight_more

            xp_linear(1)=xp_offs+xrat*(tmjd1-tmjd_mids)
            xp_linear(2)=xp_offs+xrat*(tmjd2-tmjd_mids)
            yp_linear(1)=yp_offs+yrat*(tmjd1-tmjd_mids)
            yp_linear(2)=yp_offs+yrat*(tmjd2-tmjd_mids)
            ut1kor_linear(1)=ut1kor_offs+ut1korrat*(tmjd1-tmjd_mids)
            ut1kor_linear(2)=ut1kor_offs+ut1korrat*(tmjd2-tmjd_mids)

          end if

        end if

      else

C CASE A CASE A CASE A CASE A CASE A CASE A CASE A CASE A CASE A CASE A
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCC TESMER 20.08.04 clearer way to find the a priori EOP lines CCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  get the reference values for the linear interpolation
        CALL LAGINT (  MJD, UT1KOR, IP, tmjd1, ut1kor_linear(1)   )
        CALL LAGINT (  MJD, UT1KOR, IP, tmjd2, ut1kor_linear(2)   )
        CALL LAGINT (  MJD,     XP, IP, tmjd1, xp_linear(1)   )
        CALL LAGINT (  MJD,     XP, IP, tmjd2, xp_linear(2)   )
        CALL LAGINT (  MJD,     YP, IP, tmjd1, yp_linear(1)   )
        CALL LAGINT (  MJD,     YP, IP, tmjd2, yp_linear(2)   )
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCC TESMER 20.08.04 clearer way to find the a priori EOP lines CCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      end if
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCC TESMER 20.08.04 clearer way to find the a priori EOP lines CCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


C  end as of 2004 May 12  (jboehm)


C
C     Begin to calculate
C
      do 1 i=1,nrec

c         if (mod(i,100).eq.0) write (*,671) i,nrec
c671      format ('+ ***** PROCESSING RECORD NO. ',I4,'   OF ',I4)

         read (19,rec=i) j,source,ra2000, de2000, Idatj,UTJ, d00K

         JDJ=IDATJ
C


         ut = IDATJ + UTJ/TWOPI
         UTp = UT + (1.D0/86400.D0)
         UTm = UT - (1.D0/86400.D0)




C-----------------------------------------------------------------------
C SOURCE SYSTEM 2000 IN USE CALL MP2000
C   GET TIME IN JULIAN DAYS SINCE JD 2451545.0 OR MJD 51544.5
C   THIS TIME IN JULIAN CENTURIES IS USED FOR THE NUTATION
C   ARGUMENTS IN MP2000
C-----------------------------------------------------------------------

         D2000  = ut  - 51544.5D0
         D2000p = utp - 51544.5D0
         D2000m = utm - 51544.5D0
         D00K  = 0.D0

C-----------------------------------------------------------------------
C   *** CALLING  MP2000 ***
C
C   PRECESION AND NUTATION OF SOURCE COORDINATES TO DATE COORDINATES
C  IS PERFORMED FOR ONE SECOND AFTER, ONE BEFORE AND ON TIME
C-----------------------------------------------------------------------


c
c   Update on 13-Oct-2009 (OT), d00k deleted
c

         jm1 = 1
         jm2 = 2
         jm3 = 3 

         CALL MP2000 (ra2000, de2000, D2000p, EQEP1, m_eph, idatj,
     * utj,RQP1, jm1, SNB, OFFNB, dxb, dyb, 
     * RAB, DEB, PRECB, PRB, VRB, IDELT, NUT, xm1, ym1, zm1)

         CALL MP2000 (ra2000, de2000, D2000m, EQEM1, m_eph, idatj,
     * utj,RQM1, jm2, SNA, OFFNA, dxa, dya, 
     *  RAA, DEA, PRECA, PRA, VRA, IDELT, NUT, xp1, yp1, zp1)

         CALL MP2000 (ra2000, de2000, D2000 , EQE0,  m_eph, idatj,
     * utj,RQU ,jm3, SN, OFFN, dx, dy, 
     * RAAP, DEAP, PREC, PR, VR, IDELT, NUT, xx, yy, zz)


C----------------------------------------------
C COMPUTE THE GREENWICH SIDERIAL TIME GST0
C
C   EXPRESSION FROM IERS NUMERICAL STANDARDS, IERS TN3
C  THE PREVIOUSLY CODED EXPRESIONS (IN RADIANS) PRODUCED
C  A DIFFERENCE OF 2.D-10 WHICH RESULTED IN 2 PSEC DIFF.
C  IN THE RESIDUALS
C----------------------------------------------

         DAYMOD = IDATJ
         T      = (DAYMOD - 51544.5D0) / 36525.D0              !  centuries
         T2     =  T * T
         T3     = T2 * T
         HAG    =  24110.54841D0
ccc         HBG    =  8640184.80000D0   
         HBG    =   8640184.812866D0
         HCG    =   0.093104D0
         


         HDG    =  -6.2D-6



         HSGA   = (HDG * T3 + HCG * T2 + HBG * T + HAG)/86400.D0

         HSGA   = DMOD(HSGA,1.D0) + 1.D0
         HSGA   = HSGA * TWOPI
         GST0   = HSGA

c
c    HSGA is converted from radians to hours
c
         CALL ARC (HSGA, IHSH, IHSM, HSS, 1)

C-----------------------------------------------------------------------
C COMPUTE: CSID - RATIO OF SIDEREAL TIME TO SOLAR TIME
C          OMEGA - ROTATIONAL ANGULAR VELOCITY OF THE EARTH (MOYER,1981
C-----------------------------------------------------------------------

   
 
         CSID= -5.9D-15*T2 + 5.9006D-11*T + 1.002737909350795D0  ! Aoki (1982)
         


         OMEGA = CSID / ROZ



C-----------------------------------------------------------------------
C EFFECT OF CHANGE OF NUTATION IN OMEGA
C-----------------------------------------------------------------------

         DEQEDT = (EQEP1 - EQEM1) / 2.D0

         if (jmap.eq.3) OMEGA = OMEGA + DEQEDT

         utjP = utj + (1.D0 / 86400.D0) * TWOPI
         utjM = utj - (1.D0 / 86400.D0) * TWOPI
         utj0 = utj

         

         DO 154  JMAP = 1, 3

            IF (JMAP .EQ. 1)  utj= utjP
            IF (JMAP .EQ. 2)  utj= utjM
            IF (JMAP .EQ. 3)  utj= utj0

C--------------------------------------------------------------------
C **** THE TIME TAG FOR THE INTERPOLATED VALUES OF ASTRONOMICAL AND
C      EOP DATA IS NOT CLEAR.
C   PRESENTLY TDT IS USED IN CALC 7.4
C   HOWEVER, IT LOOKS MORE APPROPRIATE TO USE UTC FOR THE EOP TABLES
C   AND TDT FOR EPHEMERIS
C---------------------------------------------------------------------
C  APPLY CORRECTION TO PASS FROM UTC TO TDT
C---------------------------------------------------------------------

            DJDM   = IDATJ + UTJ/TWOPI

cc            DJDM   = IDATJ + UTJ/TWOPI + (32.184+IDELT)/86400.D0

            DJDMP1 = DJDM + 1.D0 / 86400.D0
            DJDMM1 = DJDM - 1.D0 / 86400.D0

C-----------------------------------------------------------------------
C IERS UT1-UTC VALUES ARE THE 1-DAY SERIES.
C SHORT PERIOD TERMS IN UT1 HAVE BEEN SUBSTRACTED BEFORE INTERPOLATION
C AND MUST BE RESTORED AFTER
C-----------------------------------------------------------------------

            IF (MJD(1).LE.IDATJ.AND.MJD(IP).GT.IDATJ) GO TO 371

C IF AN ERROR IS DETECTED IN EOP TABLES FLAG A MESSAGE

            WRITE (*,*) ' ***** ERROR FOUND IN THE EOP TABLES ***** '
            WRITE (*,*) ' HIT ENTER TO CONTINUE, CTRL+C TO QUIT'
            READ (*,'(I1)') IDUMMY

371         CONTINUE

c
c    EOP Interpolation (linear or langrangian)
c
C  2004 May 12 (jboehm)
C
c    LINEAR interpolation scheme

            if (linear) then

               dut1p1 = ut1kor_linear(1) +
     .                  (ut1kor_linear(2)-ut1kor_linear(1))/
     .                  (tmjd2           -tmjd1)*
     .                  (djdmp1          -tmjd1)
               dut1m1 = ut1kor_linear(1) +
     .                  (ut1kor_linear(2)-ut1kor_linear(1))/
     .                  (tmjd2           -tmjd1)*
     .                  (djdmm1          -tmjd1)
               dut1   = ut1kor_linear(1) +
     .                  (ut1kor_linear(2)-ut1kor_linear(1))/
     .                  (tmjd2           -tmjd1)*
     .                  (djdm            -tmjd1)
               xj(jmap)= xp_linear(1) +
     .                  (xp_linear(2)    -xp_linear(1))/
     .                  (tmjd2           -tmjd1)*
     .                  (djdm            -tmjd1)
               yj(jmap)= yp_linear(1) +
     .                  (yp_linear(2)    -yp_linear(1))/
     .                  (tmjd2           -tmjd1)*
     .                  (djdm            -tmjd1)
            else

c    Lagrangian interpolation scheme (recommended for Ray's model
c    application) by IERS (short-period terms were removed before)
c    (in VLPN.FOR )

               CALL LAGINT (  MJD, UT1KOR, IP, DJDMP1, DUT1P1   )     
               CALL LAGINT (  MJD, UT1KOR, IP, DJDMM1, DUT1M1   )
               CALL LAGINT (  MJD, UT1KOR, IP,   DJDM,   DUT1   )
               CALL LAGINT (  MJD,     XP, IP,   DJDM, XJ(JMAP) )
               CALL LAGINT (  MJD,     YP, IP,   DJDM, YJ(JMAP) )
            end if
C  end as of 2004 May 12 (jboehm)

C
C APPLY SHORT TERM VARIATIONS CORRECTIONS IF REQUESTED
C

            IF (DOSHRT) THEN

C   GET TDT TO COMPUTE THE CORRECTION (CENTURIES SINCE J2000.0)

        DTDT = IDATJ + UTJ/TWOPI + (32.184D0 + IDELT + 1.d0)/86400.D0

              D2000 = (DTDT - J20 - 0.5D0) / 36525.D0
c
c    Adding of short frequencies < 35d after interpolation
c
              IF (.NOT.UT1S) CALL CUT1VR (D2000,DELUT1)
              IF (UT1S) CALL CUT1VS (D2000,DELUT1)

              DUT1P1 = DUT1P1 + DELUT1

              

        DTDT = IDATJ + UTJ/TWOPI + (32.184D0 + IDELT - 1.d0)/86400.D0

              D2000 = (DTDT - J20 - 0.5D0) / 36525.D0
              IF (.NOT.UT1S) CALL CUT1VR (D2000,DELUT1)
              IF (UT1S) CALL CUT1VS (D2000,DELUT1)

              DUT1M1 = DUT1M1 + DELUT1

              DTDT = IDATJ + UTJ/TWOPI + (32.184D0 + IDELT)/86400.D0
              D2000 = (DTDT - J20 - 0.5D0) / 36525.D0
              IF (.NOT.UT1S) CALL CUT1VR (D2000,DELUT1)
              IF (UT1S) CALL CUT1VS (D2000,DELUT1)

              DUT1 = DUT1 + DELUT1

          ENDIF

C   COMPUTE FINAL OMEGA TAKING INTO ACCOUNT THE DERIVATIVES OF THE UT1
C  CORRECTION WITH TIME




      if (i.eq.1.and.jmap.eq.3) then
  
         write (11, 778) idatj + utj/twopi, xx, yy, zz, dx, dy
 778     format (f12.5,5(2x,f13.7))

          DO ii = 1, 3
                   
             DO  J = 1, 3

                 PN (ii,J) = 0.D0

                 do k=1,3
      
                    PN(ii,J) = PN(ii,J) + PN(ii,K) * SN(K,J)  !  P*N   apparent positions for tides

                 end do

             end do  

          ENDDO


      end if  
           

c           if (jmap.eq.3) OMEGA = OMEGA * (1.d0 + (DUT1P1-DUT1M1)/2.d0)
 



c
c      Intraday variations of Omega ; Added by OT on 14 January 2019
c


           CALL OMEGA_TIDE (Dateb, domega)

c           print *, omega, domega

           OMEGA = OMEGA + domega

c           print *, omega
c
c - EnD - 14.01.2019
c

          
           DUT11  = DUT1

           DUT1OM = DUT11 * OMEGA      !  (ut1-utc)*r   Nestor's line

                
           GST01  = GST0

           IF (JMAP .EQ. 1)  EQE = EQEP1
           IF (JMAP .EQ. 2)  EQE = EQEM1
           IF (JMAP .EQ. 3)  EQE = EQE0

C-----------------------------------------------------------------------
C GMST - GREENWICH MEAN SIDEREAL TIME
C-----------------------------------------------------------------------

         
           GMST = GST01 + utj*csid + DUT1OM      !         Radians

c
c          Variarions of +/- 1.d-7 may change Omega in a reasonable way
c

                   

           HSG = GMST + EQE 

C-----------------------------------------------------------------------
C COMPUTE GREENWICH SIDEREAL TIME HSG
C-----------------------------------------------------------------------

C
C  New expression; recommended by IERS after 01.01.1997.
C

           DATEb = IDATJ + UTJ/TWOPI + (32.184+IDELT)/86400.d0

           T = (DATEb - 51544.5D0)/36525.0D0

           T2 = T*T

           AS_NODE_MOON = -0.00005939d0*T2*T2 + 0.007702d0*T2*T
     *     + 7.4722d0*T2 - 6962890.2665d0*T + 450160.398036d0
           AS_NODE_MOON = DMOD(AS_NODE_MOON,1296000d0)
           ARG = DMOD(PI * AS_NODE_MOON/648000.0D0,6.28318530718D0)

           SAS=dsin(arg)
           SAS2=dsin(2.d0*arg)

           

          


c     THE FOLLOWING CORRECTION IS AFTER CONVENTIONS 96, WHICH
c     HAS TO BE INCLUDED FOR IERS EOP COMPUTED AFTER 1997.0

           cor = (0.00264d0*sas+0.000063d0*sas2)/(15.d0*ROZ)    !  0.000176/roz
         
           
           IF (JMAP .EQ. 1)  THEN
             
              HSGP1 = HSG + cor 

           ELSE IF (JMAP .EQ. 2)  THEN
              
              HSGM1 = HSG + cor

           END IF
           
           if (jmap.eq.3) HSG = HSG + cor

c          
C  APPLY DUIRNAL AND SEMIDUIRNAL TIDE CORRECTIONS
C  RAY   (1996) MODEL IS RECOMMENDED BY IERS CONVENTIONS (1996)
C  EANES (2000) MODEL IS RECOMMENDED BY IERS CONVENTIONS (2000)
C

           IF (DTIDAL) THEN

              IF (MTIDAL.EQ.1) CALL RAY   ( DATEb, CORX, CORY, CORUT1 )
              IF (MTIDAL.EQ.2) CALL eanes ( DATEb, CORX, CORY, CORUT1 ) ! 24.01.2001

              call pmsdnut (dateb, pm)

              CORUT1 = CORUT1 * OMEGA

              IF (JMAP .EQ. 1)  THEN
             
                 HSGP1 = HSGP1 + corut1

              ELSE IF (JMAP .EQ. 2)  THEN
              
                 HSGM1 = HSGM1 + corut1

              END IF
           
              if (jmap.eq.3) HSG = HSG + corut1

                            
              XJ(JMAP) = XJ(JMAP) + CORX
              YJ(JMAP) = YJ(JMAP) + CORY

           ENDIF

C-----------------------------------------------------------------------
C   COMPUTE THE PARTIAL DERIVATIVES OF THE SOURCE POSITION WITH RESPECT
C  TO EPSILON AND PSI, NUTATION OFFSETS
C    NEGLECT THE DELAY-RATE DERIVATIVES BY THE MOMENT
C-----------------------------------------------------------------------

          IF (JMAP.EQ.3) then

               CALL NUTPAR1 (rqu,HSG,SN,PREC,PR,VR,OFFN,
     *         DKDPSI,DKDEPS,DKDPSI_v,DKDEPS_v)

          end if


C   COMPUTE GREENWICH HOUR ANGLE H

c
c   RAAP, DEAP - coordinates of the source corrected for aberration
c   It is necessary for further calculations of mapping function, etc.
c


c           H   = HSG - RAAP
c           H   = DMOD(H,TWOPI)

c           IF (H .LT. 0.D0)  H = H + TWOPI

           IF (JMAP .EQ. 1)  THEN

              HP1   = HSGP1 - RAB

           ELSE IF (JMAP .EQ. 2)  THEN

              HM1 = HSGM1 - RAA

           END IF

           if (jmap.eq.3) H = HSG - RAAP

           HP1 = DMOD(HP1,TWOPI)
           HM1 = DMOD(HM1,TWOPI)
           H = DMOD(H,TWOPI)

           if (HP1.LT.0.d0) HP1 = HP1 + TWOPI
           if (HM1.LT.0.d0) HM1 = HM1 + TWOPI
           if (H.LT.0.d0)   H = H + TWOPI


           YJ(JMAP) = YJ(JMAP) / ROS
           XJ(JMAP) = XJ(JMAP) / ROS

154      CONTINUE

         IDATJ=JDJ
         DELT = dble(IDELT)

         if (i.eq.1) print *, xj(3), yj(3)       

C---------------------------------------------------------------------
C   WRITE RESULTS TO SORTIM
C---------------------------------------------------------------------

c
c   Changes in the list of parameter on 3.11.2005 by OT
c
         WRITE (19,rec=i) j,source,ra2000, de2000, Idatj,UTJ, DELT,nut,
     +   RAB, DEB, RQP1, RAA, DEA, RQM1, RAAP, DEAP, RQU,
     +   RBAR, VBAR,
     +   HSGP1,HP1,HSGM1,HM1,HSG,H,OMEGA,XJ,YJ,
     +   CMOON,DCORG,DVELG,
     +   DKDPSI,DKDEPS,DKDPSI_v,DKDEPS_v,   !  More changes on OPENING.FOR  on 14.11.2005
     +   SNB,SNA,SN,PRECB, PRECA, PREC,PR, xx, yy, zz


1     continue


C  write the results into an ascii file
      open (31,file='eop_apriori.txt')
C  specifiy nutation model
      if (nut.eq.1) write(31,'(a23)') 'Nutation:      IAU 1980'
c
c ###########  change by OT on 22.02.2006 (Herring's model is off)
c
      if (nut.eq.2) write(31,'(a23)') 'Nutation:      MBH_2000'
c      if (nut.eq.2) write(31,'(a22)') 'Nutation:      Herring'
c      if (nut.eq.3) write(31,'(a23)') 'Nutation:      MBH_2000'
c ###########  end of change

C  ocean tides
      if (dtidal) then
        write(31,'(a18)')    'Ocean tides:   yes'
       else
        write(31,'(a17)')    'Ocean tides:   no'
      end if
      if (mtidal.eq.1) then
        write(31,'(a24)') '               Ray Model'
       else
        write(31,'(a26)') '               Eanes Model'
      end if
C  ut1 variations
      if (doshrt) then
        write(31,'(a18)')    'UT1 variation: yes'
       else
        write(31,'(a17)')    'UT1 variation: no'
      end if
      if (ut1s) then
        write(31,'(a22)') '               UT1-UTS'
       else
        write(31,'(a22)') '               UT1-UTR'
      end if
C  interpolation of xpol,ypol,dut1
      if (linear) then
        write(31,'(a25)') 'linear interpolation: yes'
       else
        write(31,'(a24)') 'linear interpolation: no'
      end if
      write(31,'(a3,3x,2(1x,f10.4))') 'mjd',tmjd1,tmjd2
      write(31,'(a2,4x,2(1x,f10.6))') 'xp', (xp_linear(i),i=1,2)
      write(31,'(a2,4x,2(1x,f10.6))') 'yp', (yp_linear(i),i=1,2)



c
c Start:  Updated by Oleg Titov on 15.05.2008
c
c

c      print *, tmjd1, tmjd2

      DLIN1 = (tmjd1 - J20 - 0.5D0) / 36525.D0
      DLIN2 = (tmjd2 - J20 - 0.5D0) / 36525.D0

c      print *, dlin1, dlin2

      CALL CUT1VR (DLIN1,DELUT_LIN1)
      CALL CUT1VR (DLIN2,DELUT_LIN2)

c      print *, delut_lin1, delut_lin2
c      print *, ut1kor_linear(1), ut1kor_linear(2)

      ut1kor_linear(1) =ut1kor_linear(1) + DELUT_LIN1
      ut1kor_linear(2) =ut1kor_linear(2) + DELUT_LIN2

c      print *, ut1kor_linear(1), ut1kor_linear(2)

c Finish:  Updated by Oleg Titov on 15.05.2008
c
c     File eop_apriori.txt is read by KVLBI.FOR
c
      write(31,'(a3,3x,2(1x,f10.7))') 'ut1',(ut1kor_linear(i),i=1,2)
      close (31)

      STOP
      end


      PROGRAM STATIO
 
C--------------------------------------------------------------------- 
C
C    IT COMPUTES CORRECTIONS DUE TO STATION POSITION AND TIME, LIKE
C  EARTH TIDES, OCEAN LOADING, POLE TIDE, VARIATION OF HOUR ANGLES
C
C    IT ALSO PROVIDES THE PARTIAL DERIVATIVES WITH RESPECT TO STATIONS
C  AND SOURCES
C
C    ALL THE INFORMATION IS STORED ON FILE STATIM FOR FURTHER USE
C
C REVISION 1991 OCTOBER 25 (N.Z.L.)
C REVISION 1993 MAY 7  (NZL)  PUT PARTIALS RESPECT TO EOP
C REVISION 1993 NOVEMBER 15  (NZL)  ADD VELOCITY MODEL
C REVISION 1997 MARCH    30  (OT)   Atmospheric loading
C REVISION 1997 APRIL    09  (OT)   Mapping function of A.Niell
C REVISION 1999 JANUARY  15  (OT)   The number of stations is corrected
C REVISION 1999 APRIL    06  (OT)   Ocean loading reading has been changed
C REVISION 2000 JULY     14  (OT)   Polar tide has been checked
C REVISION 2000 NOVEMBER 20  (OT)   Troposphere gradients have been added
C REVISION 2000 NOVEMBER 20  (VT)   Atmospheric loading computed by
C                                          Scherneck added,
C                                          Antenna deformation refined
C REVISION 2001 JANUARY  31  (OT)   COMMON /zenith/ removed;
C                                          ZWET removed from STATIM
C REVISION 2001  MARCH  24  (OT)    COMMON /BAR1/ instead of /BAR/
C REVISION 2001 APRIL 06 (OT)  Remarks by Volker Tesmer have been applied for
C REVISION 2001 September, 24 (OT)  Geocentiric velocities of the Sun
C       are from EPHEM.DAT rather than calculationg from geocentric coordinates
C REVISION 2002 MAY, 04 (OT) Advanced formula for thermal deformation
C                                 has been applied
C REVISION 2002 JUNE 12, (OT)  'NUT' is read in SORTIM (19)
C LAST REVISION 2002 JUNE 19, (OT)  Rotation R3 moved from PARTIALS.FOR
c revision 2002 July 29 jboehm, new discrete models are included:
c                                - isobaric mapping function
c                                - snow loading
c                                - soil moisture
c                                - ocean loading (additional)
c resivion 2003 Feb. 17, raytraced mapping function was added
c REVISION 2003 May 19 jboehm: iavail, iavailr and iavailq
c                               were changed to integer arrays
c REVISION 2004 January 27 (OT): Gilcreek postseismic motion
C REVISION 2004 August  19 (O.T.)  Corrected bug for the Lena Skurihina
C                     correction term for thermal deformation
C REVISION 2006 March 02 (O.T.)  Correction for the OCCAM 6.2
C REVISION 2006 March 22 (O.T.)  Correction for the Lena Skurihina's
C               term is made as option; numerous changes due to implementation of
C               VMF1 and apriori gradients; changes due to transition to the
C               mean system of the Equinox
C LAST REVISION 2006 October, 03 (O.T.)  Jupiter is added
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
C---------------------------------------------------------------------

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INCLUDE 'OCCAM_N.FI'

      CHARACTER STAN*8,QUASAR*8,AXTYP*4,STAT*8,CYES*1,STNAME*8,IVETO*1
c+ vtesmer
      CHARACTER*1 YESMAR, YESTRC, YESAZD, YESAGR
c- vtesmer
      character*(8*nstations) atmsta, ocesta, snosta, soista
      character stc*4

      LOGICAL OPTTD,OPTOC,OPTAT,OPTPT,OPTLOV,NOOCL(nstations),
     *        nothe(nstations)
c+ jboehm
      logical optsnow, optsoil, optoc1
c- jboehm
c+ otitov
      logical opttd1
c- otitov
      DOUBLE PRECISION NSOCL

c atmosphere loading
      dimension atmmjd(nstations,30)
      dimension atmx(nstations,30),atmy(nstations,30),atmz(nstations,30)
      dimension tmjddec(30),dxatm(30),dyatm(30),dzatm(30)
c+ jboehm discrete ocean loading
      dimension ocemjd(nstations,30)
      dimension ocex(nstations,30),ocey(nstations,30),ocez(nstations,30)
      dimension dxoce(30),dyoce(30),dzoce(30)
c snow loading
      dimension snomjd(nstations,30)
      dimension snox(nstations,30),snoy(nstations,30),snoz(nstations,30)
      dimension dxsno(30),dysno(30),dzsno(30)
c soil moisture
      dimension soimjd(nstations,30)
      dimension soix(nstations,30),soiy(nstations,30),soiz(nstations,30)
      dimension dxsoi(30),dysoi(30),dzsoi(30)
c- jboehm

      DIMENSION RQU(3),RQP1(3),RQM1(3),VBAR(3,3),RBAR(3,3),ABAR(3,3)
      DIMENSION HSG(3),H(3),HSTA(3),HA(3),Z(3),HA2(3)
      DIMENSION PCOORD(5),PCOORR(5),DE(3)
      DIMENSION EROT(3,3),EDROT(3,3,2),XPOLE(3),YPOLE(3)
      DIMENSION CORZ(3), SB(3,3), SA(3,3), SS(3,3)
      DIMENSION CRX(3),CRY(3),CRZ(3),CGLONG(3),CPHI(3),CELHGT(3),CR0(3)
      DIMENSION CGLONG2(3),CPHI2(3),CELHGT2(3),CR02(3)
      DIMENSION CRX2(3),CRY2(3),CRZ2(3),CRX_C(3),CRY_C(3),CRZ_C(3),pp(3)
      DIMENSION CRX1(3),CRY1(3),CRZ1(3)

      DIMENSION crx4(3), cry4(3), crz4(3)

      DIMENSION crx5(3), cry5(3), crz5(3)


      DIMENSION crx6(3), cry6(3), crz6(3)


      DIMENSION ppolar_rate(3), p4(3), DSDUT_V(3), pv(3)
      

c+ otitov 04.06.2010
      DIMENSION v2(3), w2x(3), w3(3), w2(3), r4(3,3), r5(3,3), r5_1(3,3)
      dimension qe(3,3), rr3(3,3), trf(3), qd(3,3), qf(3), w4(3)
      dimension V_after_PN(3), R_after_PN(3), XLUN_after_PN(3)

      Dimension A2(3), a4(3)   ! Geocentric acceleration
      Dimension Vel_p(3), vel_m(3) 

c- otitov

c+ jboehm iveto(4) to iveto(7)
c+ otitov 22.03.2006    iveto(7) to iveto(8)

      DIMENSION IVETO(8),CAZ(3),SAZ(3)
      DIMENSION AZIMUTH(3)

c- otitov
c- jboehm
      dimension a_grad(2)
      DIMENSION CMOON(3,3),DCORG(3,3), DVELG(3,3), xjup(3,3)   !  Added on 03.10.2006
      DIMENSION SNB(3,3), SNA(3,3), SN(3,3), CR(3), CRA(3), CRB(3)
      DIMENSION CR_1(3), CRA_1(3), CRB_1(3)
      DIMENSION PREC(3,3), PRECB(3,3), PRECA(3,3), PR(3), hh(3,3)

      DIMENSION dsdut_t(3), dsdxw_t(3), dsdyw_t(3)
 
c+ jboehm for imf
      dimension tmjdimf(8),dz200(8),dsmfw3(8),dazim(8),delev(8),
     .          adz200(nstations,8),adsmfw3(nstations,8),
     .          adazim(nstations,8),adelev(nstations,8),
     .          tmjdimf1(8),iavail(nstations)
c- jboehm
c+ jboehm for rmf
      dimension tmjdrmf(8),dah(8),dbh(8),dch(8),daw(8),dbw(8),dcw(8),
     .          adah(nstations,8),adbh(nstations,8),adch(nstations,8),
     .          adaw(nstations,8),adbw(nstations,8),adcw(nstations,8),
     .          tmjdrmf1(8),iavailr(nstations)
c- jboehm
c+ jboehm for qmf
      dimension tmjdqmf(8),dahq(8),dawq(8),
     .          adahq(nstations,8),adawq(nstations,8),
     .          tmjdqmf1(8),iavailq(nstations)
c- jboehm
c+ jboehm for vmf1, 2005-Oct-28
      dimension tmjdvmf(8),dahv(8),dawv(8),
     .          adahv(nstations,8),adawv(nstations,8),
     .          tmjdvmf1(8),iavailv(nstations)
     .          ,press1(8), tempe1(8)                                   !  24-Mar-2014  
     .          ,temperat(nstations,8), pressure(nstations,8)           !  24-Mar-2014 

      dimension hde(8), wde(8), hd1(nstations,8), wd1(nstations,8)   !  7-Nov-2019

c- jboehm
c+ jboehm 2005-Oct-24 for gra (a priori gradients davis 93)
      dimension tmjdgra(8),dgrnh(8),dgreh(8),dgrnw(8),dgrew(8),
     .          adgrnh(nstations,8),adgreh(nstations,8),
     .          adgrnw(nstations,8),adgrew(nstations,8),
     .          tmjdgra1(8),iavailg(nstations)
C- jboehm 2005-Oct-24 for gra



      COMMON /CORR/ CTIDE(3),cmniel_w(2),cmniel_d(2),AXKT(2),therm_d
      COMMON /PART/ PCOORD,PCOORR,PNUTAT(2),PPOLAR(3,3)
      COMMON /PDATA/ SD,CD,SIDGRE,DKDPSI(3),DKDEPS(3),DSDXW(3,3),
     * DSDYW(3,3),DSDUT(3,3), DKDPSI_v(3),DKDEPS_v(3)
      COMMON /DATEb / IDATJ,UTJ,XJM,YJM,IDELT
      COMMON /COCPOS/ GLONG,PHI,ELHGT,RX,RY,RZ
      COMMON /COCDEL/ CGLONG,CPHI,CELHGT,CRX,CRY,CRZ,HA,Z,DE
      COMMON /COCDEL2/ CGLONG2,CPHI2,CELHGT2,HA2

      COMMON /PARM/ OMEGA
      COMMON /METDAT/ TEMP,PRES,HUMD,IFAC
      COMMON /OPT/ opttd,OPTOC,OPTAT,OPTPT,OPTLOV
c+ jboehm
     .             ,optsnow,optsoil,optoc1
c- jboehm
c+ otitov  22.03.2006
     .             ,opttd1
c- otitov  22.03.2006
      COMMON /OCLOAD/ NOOCL,ROCL(nstations,2,noclotds),
     *           NSOCL(nstations,2,noclotds),EWOCL(nstations,2,noclotds)
      COMMON /THERM/ ga(nstations),hd(nstations),gf(nstations),
     *     hf(nstations),hp(nstations),hs(nstations),
     *     hv(nstations),ifo(nstations),ref_temp(nstations)

cc      COMMON /BAR1/ RBAR0(3),VBAR0(3)    !  Usable for tidal  effects in MATTHEW.FOR
cc      COMMON /MOON/ XLUN(3)        !  Usable for tidal  effects in MATTHEW.FOR

      COMMON /PHYS/ C, FL, AH, AU, J20
      COMMON /MATH/ PI, TWOPI, ROG, ROS, ROZ, TOL
      COMMON /atscher/ istatm,atmsta,atmmjd,atmx,atmy,atmz
c+ jboehm
      COMMON /ocean/   istoce,ocesta,ocemjd,ocex,ocey,ocez
      COMMON /snow/    istsno,snosta,snomjd,snox,snoy,snoz
      COMMON /soil/    istsoi,soista,soimjd,soix,soiy,soiz
      common /corrimf/ tmjdimf,adz200,adsmfw3,adazim,adelev,iavail,
     .                 azimuth
      common /corrrmf/ tmjdrmf,adah,adbh,adch,adaw,adbw,adcw,iavailr
      common /corrqmf/ tmjdqmf,adahq,adawq,iavailq
      common /corrout/ z200,smfw3,azim,elev,dimfh,dimfw,gimfh,grn,gre,
     .                 drmfh,drmfw,dqmfh,dqmfw
c- jboehm
C+ jboehm, 2005-Sep-5
      common /corrgmf/ gmfh,gmfw
C- jboehm, 2005-Sep-5
C+ jboehm, 2005-Oct-24
      common /corrgra/ tmjdgra,adgrnh,adgreh,adgrnw,adgrew,iavailg,
     .                 grnh,greh,grnw,grew
C- jboehm, 2005-Oct-24
c+ jboehm, 2005-Oct-28

      common /corrvmf/ tmjdvmf,adahv,adawv,hd1,wd1,iavailv,vmf1h,vmf1w
     * ,vmf1hd,vmf1wd,pressure,temperat    !   24-Mar-2014; 5-Nov-2019 (OT)

c- jboehm, 2005-Oct-28

C-------------------
C  GET PARAMETERS
C-------------------

      CALL PARAM

      open (12, access='append', file='solution.dat', status='unknown')
      open (13, file='statio1.dat', status='unknown')
c      open (15, file='r5.dat', status='unknown')
      open (102, file='station.dat', status='unknown')

      opttd = .TRUE.
      OPTOC = .TRUE.
      OPTAT = .TRUE.
      OPTPT = .TRUE.
      OPTLOV=.FALSE.
c+ jboehm
      OPTSNOW=.TRUE.
      OPTSOIL=.TRUE.
      optoc1 =.TRUE.
c- jboehm
c+ otitov   22.03.2006
      opttd1=.FALSE.
c- otitov

C--------------------------------------
C     INITIALIZATION OF VARIABLES
C--------------------------------------

      DO I=1,nstations
         NOOCL(I)=.TRUE.
         nothe(I)=.TRUE.
         DO J=1,2
            DO K=1,noclotds
               ROCL(I,J,K)=0.D0
               NSOCL(I,J,K)=0.D0
               EWOCL(I,J,K)=0.D0
            END DO
         end do
      end do

C
C    CHECK OPTIONS SELECTED TO ADD CORRECTIONS TO STATION POSITION
C

601   WRITE (*,500)
500   FORMAT ('  BY DEFAULT ALL CORRECTIONS WILL BE ADDED TO THE DATA',
     *' BASES'/' DO YOU WANT TO CHANGE THIS ?    (Y / (N))     HIT "Q"',
     *' TO QUIT')
      READ (*,600,ERR=601) CYES
600   FORMAT (A1)

      IF (CYES.EQ.'Q'.OR.CYES.EQ.'q') STOP 'PROGRAM TERMINATED BY USER'

      IF (CYES.EQ.'Y'.OR.CYES.EQ.'y') THEN

c+ jboehm
c+   otitov 22.03.2006  - Extra thermal deformation option is added
602      continue
         WRITE (*,505)
505      FORMAT (' THIS CORRECTIONS ARE APPLIED...'/
     *' A. ANTENNA DEF.   B. OCEAN LOAD.   C. ATMOS. LOAD.   D. POLE ',
     *'TIDE   E. ANTENNA DEF EXTRA.  F. SNOW LOADING   G. SOIL MOISTURE
     * H. OCEAN LOADING 1'/
     *' SELECT THE CORRECTIONS TO NEGLECT (or ADD for the "E" '/
     *'- - - - - - - -')
c- jboehm

         READ (*,605) IVETO
605      FORMAT (8(A1,1X))             !  Changed by Johannes

         WRITE (*,510) IVETO
c+ jboehm 4 to 8
510      FORMAT (' THIS CORRECTIONS WILL BE NEGLECTED (added)'/3X,7(A1,1
     * X)/ ' DO YOU WANT TO CHANGE THIS ?   (Y / (N))  HIT "Q" TO QUIT')
c- jboehm

         READ (*,600) CYES

         IF (CYES.EQ.'Q'.OR.CYES.EQ.'q') STOP 'PROGRAM TERMINATED'

         IF (CYES.EQ.'Y'.OR.CYES.EQ.'y') GO TO 602

C    CHANGE FLAG VALUES TO THE CORRECTIONS

         DO ID=1,8  !  Changed by Oleg on 22.03.2006
            PRINT *, IVETO(ID)
            IF (IVETO(ID).EQ.'A'.OR.IVETO(ID).EQ.'a') then
               opttd = .FALSE.
               write (12,491)
            end if
            IF (IVETO(ID).EQ.'B'.OR.IVETO(ID).EQ.'b') then
               OPTOC = .FALSE.
               write (12,492)
            end if
            IF (IVETO(ID).EQ.'C'.OR.IVETO(ID).EQ.'c') then
               OPTAT = .FALSE.
               write (12,493)
            end if
            IF (IVETO(ID).EQ.'D'.OR.IVETO(ID).EQ.'d') then
               OPTPT = .FALSE.
               write (12,494)
            end if
            IF (IVETO(ID).EQ.'E'.OR.IVETO(ID).EQ.'e') then
               opttd1 = .true.
               write (12,495)
            end if
            IF (IVETO(ID).EQ.'F'.OR.IVETO(ID).EQ.'f') then
               OPTSNOW = .FALSE.
               write (12,496)
            end if
            IF (IVETO(ID).EQ.'G'.OR.IVETO(ID).EQ.'g') then
               OPTSOIL = .FALSE.
               write (12,497)
            end if
            IF (IVETO(ID).EQ.'H'.OR.IVETO(ID).EQ.'h') then
               OPTOC1 = .FALSE.
               write (12,498)
            end if
c-   otitov 22.03.2006  - Extra thermal deformation option is added
         end do
      else
        write (12,490)
      end if
      write (12,*)

 490  format (1x,'Seven corrections are implemented, no extra thermal')
 491  format (1x,'Antenna thermal deformations are NOT included')
 492  format (1x,'Ocean loading is NOT included')
 493  format (1x,'Atmosphere loading is NOT included')
 494  format (1x,'Pole tide correction is NOT included')
c+   otitov 22.03.2006  - Extra thermal deformation option is added
 495  format (1x,'Extra antenna thermal deformations is added')
 496  format (1x,'Snow loading is NOT included')
 497  format (1x,'Soil moisture is NOT included')
 498  format (1x,'Ocean loading 1 is NOT included')
c-   otitov 22.03.2006  - Extra thermal deformation option is added

C
C     SELECT THE PARTIAL DERIVATIVES TO BE INCLUDED ON THE DATA SET
C    FORGET ABOUT LOVE NUMBERS DERIVATIVES BY NOW


C      WRITE (*,100)
C100   FORMAT (//' BY DEFAULT WE WILL EXCLUDE THE PARTIAL DERIVATIVES'/
C     1' WITH RESPECT TO LOVE NUMBERS FROM THE DATA SET'//
C     2' DO YOU WANT TO INCLUDE THEM ?   (Y/(N))    HIT "Q" TO QUIT'/)

C      READ (*,200) CYES
C200   FORMAT (A1)

C      IF (CYES.EQ.'Q'.OR.CYES.EQ.'q') STOP  'PROGRAM TERMINATED BY USER'
C      IF (CYES.EQ.'Y'.OR.CYES.EQ.'y')     OPTLOV=.TRUE.

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc TESMER 02.11.04 also add the bernese way for hydrostatic ZD ccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C    do the a priori hydrostatic ZD the bernese way ?
      WRITE (*,105)
105   FORMAT (/' DO YOU WANT TO : '/
     1' - REFINE APRIORI HYDROSTATIC ZD WITH MARINI (DEFAULT) => "Y"'/
     2' - OR NOT (WHICH IS THE BERNESE WAY) => "N"'/
     3' - or quit "Q" '/)
      YESMAR='Y'
      READ (*,106) YESMAR
106   FORMAT (A1)
      IF (YESMAR.EQ.'Q'.OR.YESMAR.EQ.'q') STOP  'PROGRAM TERMINATED '
      IF (YESMAR.EQ.'N'.OR.YESMAR.EQ.'n') print*,'U USE BERNESE WAY!'
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc TESMER 02.11.04 also add the bernese way for hydrostatic ZD ccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc TESMER 23.02.05 provide trop refraction correction as option cc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      WRITE (*,115)
115   FORMAT (/' DO YOU WANT TO : '/
     1' - correct ZD due to tropospheric refraction (DEFAULT) => "Y"'/
     2' - OR NOT => "N"'/
     3' - or quit "Q" '/)
      YESTRC='Y'
      READ (*,116) YESTRC
116   FORMAT (A1)
      IF (YESTRC.EQ.'Q'.OR.YESTRC.EQ.'q') STOP  'PROGRAM TERMINATED '
      IF (YESTRC.EQ.'N'.OR.YESTRC.EQ.'n') print*,'refr. not corr.  !'
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc TESMER 23.02.05 provide trop refraction correction as option cc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc TESMER 23.02.05 NOW USE const. hydr. ZD OF BERNESE apriori ccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      WRITE (*,117)
117   FORMAT (/' DO YOU WANT TO : '/
     1' - CALCULATE APR. ZD FROM MEAS. SURF. PRESS. (DEFAULT) => "Y"'/
     2' - OR USE CONSTANT ZD FROM "hydr_zd_apriori.txt" => "N"'/
     3' - or quit "Q" '/)
      YESAZD='Y'
      READ (*,118) YESAZD
118   FORMAT (A1)
      IF (YESAZD.EQ.'Q'.OR.YESAZD.EQ.'q') STOP  'PROGRAM TERMINATED '
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc TESMER 26.07.05 READ APRIORI ZD VALUES FROM FILE           ccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      IF (YESAZD.EQ.'N'.OR.YESAZD.EQ.'n') then
        print*,'CONSTANT APR ZD  !'
        OPEN (44,FILE='hydr_zd_apriori.txt',STATUS='OLD',ERR=683)
      endif
684   goto 685
683   YESAZD='Y'
      print*,' '
      print*,'File "hydr_zd_apriori.txt" does not exist!'
      print*,'ignoring this option, computation done as always'
      print*,' '
685   continue
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc TESMER 26.07.05 READ APRIORI ZD VALUES FROM FILE           ccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc TESMER 23.02.05 NOW USE const. hydr. ZD OF BERNESE apriori ccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc TESMER 23.02.05 NOW USE apriori gradients from GSFC        ccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      WRITE (*,119)
c+ 2005-Oct-28, jboehm
119   FORMAT (/' DO YOU WANT TO : '/
     1' - USE 0.0 AS APRIORI GRADIENTS (DEFAULT) => "Y"'/
     2' - USE "gradients_apriori.txt" AS APRIORI GRADIENTS => "N"'/
     3' - USE GRADIENTS FROM GRA_ALL => "G"'/
     4' - or quit "Q" '/)
      YESAGR='Y'
      READ (*,120) YESAGR
      IF (YESAGR.EQ.'') YESAGR='Y'
c- 2005-Oct-28, jboehm
120   FORMAT (A1)
      IF (YESAGR.EQ.'Q'.OR.YESAGR.EQ.'q') STOP  'PROGRAM TERMINATED '
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc TESMER 26.07.05 READ APRIORI GRADIENT VALUES FROM FILE     ccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IF (YESAGR.EQ.'N'.OR.YESAGR.EQ.'n') then
        print*,'APR GRA applied  !'
        OPEN (43,FILE='gradients_apriori.txt',STATUS='OLD',ERR=693)
      endif
c+ jboehm, 2005-Oct-28
      IF (YESAGR.EQ.'G'.OR.YESAGR.EQ.'g') then
        print*,'apriori gradients from NWM applied  !'
      endif
c- jboehm, 2005-Oct-28
694   goto 695
693   YESAGR='Y'
      print*,' '
      print*,'File "gradients_apriori.txt" does not exist!'
      print*,'values are set to zero'
      print*,' '
695   continue
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc TESMER 26.07.05 READ APRIORI GRADIENT VALUES FROM FILE     ccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc TESMER 23.02.05 NOW USE apriori gradients from GSFC        ccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      write (*,289)
289   FORMAT (//' ***** PROGRAM RUNNING *****       Please Wait'//)

C   OPENING FILES

      CALL UNFOLD (17,'BASTIM')
      CALL UNFOLD (18,'STATIM')
      CALL UNFOLD (19,'SORTIM')
      CALL UNFOLD (20,'STACAT')

C     PART No. 1.   COMPUTATION OF ELLIPSOIDAL COORDINATES FOR ALL THE
C  STATIONS INVOLVED IN THIS RUNNING AND GENERATE "STACAT", STATIONS
C  CATALOGUE

      IF (OPTOC) OPEN (11,FILE='LOADING.OCE',STATUS='OLD',ERR=653)
      IF (OPTTD) OPEN (21,FILE='THERMAL.DEF',STATUS='OLD',ERR=673)

      READ (20,REC=1) IREC

C
C PART No. 1, COMPUTE CORRECTIONS
C

      CRVX = 0.d0
      CRVY = 0.d0
      CRVZ = 0.d0 
    
      do i = 1,3

         v2(i) = 0.d0

      end do
 

      DO 1 I=1,IREC
         READ (20,REC=I) JST,STAN,RX,RY,RZ,EPOCH,VX,VY,VZ,AXTYP,OFFS,
     *    stc,nst

         CALL TRANSF (RX,RY,RZ,PHI,GLONG,ELHGT,R0,0)

C      COMPUTE SIN AND COSINE OF PHI, AND THE LENGTH OF THE EQUATORIAL
C  COMPONENT OF THE POSITION VECTOR

         RXY = DSQRT(RX * RX + RY * RY)

         WRITE (20,REC=I) JST,STAN,RX,RY,RZ,EPOCH,VX,VY,VZ,AXTYP,
     *                    OFFS,stc,nst,PHI,GLONG,ELHGT,R0,RXY

C   GET OCEAN-LOADING AMPLITUDES AND PHASES FOR THE STATIONS

         LOADFC=0
         IF (OPTOC) THEN
            REWIND 11
292         CONTINUE
            READ (11,109,END=991) STNAME
109         FORMAT (A8)
            IF (STAN.EQ.STNAME) then

               NOOCL(I)=.FALSE.
               READ (11,165,END=1) (ROCL(I,1,J),J=1,noclotds)
               READ (11,165,END=1) (EWOCL(I,1,J),J=1,noclotds)
               READ (11,165,END=1) (NSOCL(I,1,J),J=1,noclotds)
               READ (11,166,END=1) (ROCL(I,2,J),J=1,noclotds)
               READ (11,166,END=1) (EWOCL(I,2,J),J=1,noclotds)
               READ (11,166,END=1) (NSOCL(I,2,J),J=1,noclotds)
165            FORMAT (1x,11(F7.5))
166            FORMAT (1x,11(F7.1))

            end if

            IF (NOOCL(I)) GO TO 292

C IF NO OCEAN LOADING, SEND A MESSAGE

991         continue
            IF (NOOCL(I)) THEN
               WRITE (*,*) ' NO OCEAN LOADING DATA FOR STATION ',STAN
               WRITE (*,*) '         - HIT ANY KEY TO CONTINUE - '
               READ (*,'(A1)') CYES
            ENDIF
         ENDIF

1     CONTINUE               !      end of first cycle

C
C PART No. 2, COMPUTE CORRECTIONS
C

c set stationcounter for atmospheric loading to zero

      istatm=0
      do iatm=1,8*nstations
        atmsta(iatm:iatm)=' '
      end do
c+ jboehm
      istoce=0
      do ioce=1,8*nstations
        ocesta(ioce:ioce)=' '
      end do
c
      istsno=0
      do isno=1,8*nstations
        snosta(isno:isno)=' '
      end do
c
      istsoi=0
      do isoi=1,8*nstations
        soista(isoi:isoi)=' '
      end do
c- jboehm

      READ (20,REC=1) IREC

     
      do 2 i=1,irec

         READ (20, REC=I) JP,STAN

c+ jboehm
         call read_imf (stan,tmjdimf1,dz200,dsmfw3,dazim,delev,iav_imf)
         if (iav_imf.eq.1) then        !   18.01.2005 OT
           iavail(i) = 1
           do it = 1,8
             tmjdimf(it) = tmjdimf1(it)
           end do
         end if
         do it = 1,8
           adz200(i,it)  = dz200(it)
           adsmfw3(i,it) = dsmfw3(it)
           adazim(i,it)  = dazim(it)
           adelev(i,it)  = delev(it)
         end do
c- jboehm

c+ jboehm for raytraced mapping function
         call read_rmf (stan,tmjdrmf1,dah,dbh,dch,daw,dbw,dcw,iav_rmf)

         if (iav_rmf.eq.1) then       !   18.01.2005 OT
           iavailr(i) = 1
           do it = 1,8
             tmjdrmf(it) = tmjdrmf1(it)
           end do
         end if
         do it = 1,8
           adah(i,it)  = dah(it)
           adbh(i,it)  = dbh(it)
           adch(i,it)  = dch(it)
           adaw(i,it)  = daw(it)
           adbw(i,it)  = dbw(it)
           adcw(i,it)  = dcw(it)
         end do
c- jboehm

c+ jboehm for quasi-raytraced mapping function
         call read_vmf (stan,tmjdqmf1,dahq,dawq,iav_vmf)
         if (iav_vmf.eq.1) then       !   18.01.2005 OT
           iavailq(i) = 1
           do it = 1,8
             tmjdqmf(it) = tmjdqmf1(it)
           end do
         end if
         do it = 1,8
           adahq(i,it)  = dahq(it)
           adawq(i,it)  = dawq(it)
         end do
c- jboehm

c+ 2005-Oct-28, jboehm for Vienna Mapping Function 1 VMF1

c        call read_vm1 (stan,tmjdvmf1,dahv,dawv,iav_vm1) 

c
c    hydrostatic delay (hd) and wet delay (wd) added  7-Nov-2019
c
      call read_vm1 (stan,tmjdvmf1,dahv,dawv,hde,wde,
     *   press1,tempe1,iav_vm1)   !   24-Mar-2014; pressure and temperature are from VM1 file 

         

         if (iav_vm1.eq.1) then

           iavailv(i) = 1
           do it = 1,8

             tmjdvmf(it) = tmjdvmf1(it)

           end do

         end if

         do it = 1,8

           adahv(i,it)  = dahv(it)
           adawv(i,it)  = dawv(it)

           hd1(i,it) = hde(it)   !  9-Nov-2019
           wd1(i,it) = wde(it)   !  9-Nov-2019
                         
           pressure(i,it)  = press1(it)                         !    24-Mar-2014
           temperat(i,it)  = tempe1(it)                         !    24-Mar-2014

         end do

c- 2005-Oct-28, jboehm

c+ jboehm 2005-Oct-24 for davis gradients

         call read_gra (stan,tmjdgra1,dgrnh,dgreh,dgrnw,dgrew,iav_grad)

         if (iav_grad.eq.1) then
           iavailg(i) = 1
           do it = 1,8
             tmjdgra(it) = tmjdgra1(it)
           end do
         end if
         do it = 1,8
           adgrnh(i,it)  = dgrnh(it)
           adgreh(i,it)  = dgreh(it)
           adgrnw(i,it)  = dgrnw(it)
           adgrew(i,it)  = dgrew(it)
         end do
c- jboehm 2005-Oct-24

C ATMOSPHERIC LOADING

         if (optat) then

           iista=i

c VERY NEW SUBROUTINE FOR ATMLOADING FROM SCHERNECK
c+ jboehm
           call atm_sch('ATM',stan,iista,tmjddec,dxatm,dyatm,dzatm)
c- jboehm

c if corrections found, give them on fields
c
c =>
c
c istatm 	      nstations with corrections
c atmsta(1:8*istatm)  names of stations with corr in row
c atmmjd(istatm,j)    j times of corrections of stations
c atmx(istatm,j)      j dx corrections of stations
c atmy(istatm,j)      j dy corrections of stations
c atmz(istatm,j)      j dz corrections of stations

          if (tmjddec(1).ne.0) then
            istatm=istatm+1
            atmsta(((istatm*8)-7):(istatm*8))=stan(1:8)
            do j=1,30
              atmmjd(istatm,j)=tmjddec(j)
              atmx(istatm,j)=dxatm(j)
              atmy(istatm,j)=dyatm(j)
              atmz(istatm,j)=dzatm(j)
            end do
          end if

        end if

c+ jboehm
cc
C OCEAN LOADING (DISCRETE)
cc

        if (optoc1) then

          iista=i

          call atm_sch('OCE',stan,iista,tmjddec,dxoce,dyoce,dzoce)

          if (tmjddec(1).ne.0) then
            istoce=istoce+1
            ocesta(((istoce*8)-7):(istoce*8))=stan(1:8)
            do j=1,30
              ocemjd(istoce,j)=tmjddec(j)
              ocex(istoce,j)=dxoce(j)
              ocey(istoce,j)=dyoce(j)
              ocez(istoce,j)=dzoce(j)
            end do
          end if

        end if

cc
C SNOW LOADING
cc

        if (optsnow) then

          iista=i

          call atm_sch('SNO',stan,iista,tmjddec,dxsno,dysno,dzsno)

          if (tmjddec(1).ne.0) then
            istsno=istsno+1
            snosta(((istsno*8)-7):(istsno*8))=stan(1:8)
            do j=1,30
              snomjd(istsno,j)=tmjddec(j)
              snox(istsno,j)=dxsno(j)
              snoy(istsno,j)=dysno(j)
              snoz(istsno,j)=dzsno(j)
            end do
          end if

        end if

cc
C SOIL MOISTURE LOADING
cc

        if (optsoil) then

          iista=i

c VERY NEW SUBROUTINE FOR ATMLOADING FROM SCHERNECK
          call atm_sch('SOI',stan,iista,tmjddec,dxsoi,dysoi,dzsoi)

          if (tmjddec(1).ne.0) then
            istsoi=istsoi+1
            soista(((istsoi*8)-7):(istsoi*8))=stan(1:8)
            do j=1,30
              soimjd(istsoi,j)=tmjddec(j)
              soix(istsoi,j)=dxsoi(j)
              soiy(istsoi,j)=dysoi(j)
              soiz(istsoi,j)=dzsoi(j)
            end do
          end if

        end if

c ATNENNA DEFORMATION
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c read also info if prime or secondary focus
         if (opttd) then
            rewind 21
792         CONTINUE
            READ (21,112,END=995) STNAME
112         FORMAT (A8)
            IF (STAN.NE.STNAME) go to 792
            nothe(i)=.FALSE.
            read (21,121,end=32) gf(i),hf(i),ga(i),
     *                           hp(i),hv(i),hs(i),hd(i),
     *                           ifo(i),ref_temp(i)
121      format (f9.5,f6.1,f11.6,4(f6.1),i6,f6.1)

C
C IF NO thermal antenna deformation info, SEND A MESSAGE
C

995         IF (nothe(i)) THEN
      WRITE (*,*) ' ***** NO THERMAL DEF. COEFFICIENT FOR STATION ',STAN
               WRITE (*,*) '         - HIT ANY KEY TO CONTINUE - '
               READ (*,'(A1)') CYES
               gf(i)=0.d0
               hf(i)=0.d0
               ga(i)=0.d0
               hp(i)=0.d0
               hv(i)=0.d0
               hs(i)=0.d0
               hd(i)=0.d0
               ref_temp(i)=999.d0
            ENDIF

C IF NO reference temperature, SEND A MESSAGE
            if (ref_temp(i).eq.999.d0) then
              write(*,*) 'no reference temperature for ',stan
              write(*,996)
            end if
 996     format ('thermal deformation will be zero'//)

         ENDIF
 32      continue

 2    continue             !              end of 2nd cycle

C
C PART No. 3, COMPUTE CORRECTIONS AND PARTIAL DERIVATIVES
C

c+ jboehm's nice toy
c      open (11,file='INFO/statcorr.txt')
c      write (11,'(a8,a12,7a27)') 'STATION','MJD','SOLID EARTH TIDE',
c     .          'OCEAN LOADING','ATMOSPHERE LOADING','POLE TIDE',
c     .          'SNOW LOADING','SOIL MOISTURE','OCEAN LOADING 1'
c      write (11,'(a24)') '# geocentric coordinates'
c      open (12,file='imf_info.dat')
c- jboehm


      do i = 1,3

         crx(i) = 0.d0
         cry(i) = 0.d0
         crz(i) = 0.d0

         crx4(i) = 0.d0
         cry4(i) = 0.d0
         crz4(i) = 0.d0

      end do
         

      READ (18,REC=1) IREC1

      DO 3 I=1,IREC1


         READ (18,rec=I) JSTA,ISTCAT,ISOUR,TEMP,PRES,HUMD,IFAC,CABL
C
C OPEN FILE SORTIM    List of parameters changed on 3.11.2005 !!
C
         READ (19,rec=iSOUR) jSOR,quasar,ra2000, de2000,Idatj,UTJ0,DELT,
     +                   nut,RAB, DE(1), RQP1,
     +                       RAA, DE(2), RQM1,
     +                       RAAP, DE(3), RQU, RBAR, VBAR,
     +   HSG(1),H(1),HSG(2),H(2),HSG(3),H(3),OMEGA,XPOLE,YPOLE,
     +   CMOON, DCORG, DVELG,
     +   DKDPSI, DKDEPS, DKDPSI_v, DKDEPS_v,
     +   SNB, SNA, SN, PRECB, PRECA, PREC, PR


         IDELT = DELT + .0001

         READ (20,REC=ISTCAT) JP,STAT,RX,RY,RZ,EPOCH,VX,VY,VZ,AXTYP,
     *                        OFFS,stc,nst,PHI,GLONG, ELHGT,R0,RXY

         IY = EPOCH
         ID = (EPOCH-IY)*365.25
         IH = 0
         IM = 0
         S  = 0.d0

C    ADD VELOCITY MODEL

         CALL JULDAT (IY,ID,IH,IM,S,IDJEP,UTJEP,0)

         TS = (IDATJ+UTJ0/twopi-IDJEP)/365.25d0



c+ OT
c
c  Gilcreek postseismic motion
c
         i_earthquake=52581    !   Earthquake epoch

         if (idatj.ge.i_earthquake) then

            TG = (IDATJ-UTJEP-52581.d0)/365.25d0


            if (STAT.EQ.'GILCREEK') then

               tau=0.25d0

               du = 0.008d0
               dn = -0.0478d0 - 0.0237d0*(1.d0-dexp(-tg/tau))
               dea = 0.0226d0 + 0.0086d0*(1.d0-dexp(-tg/tau))

               offx= -dsin(phi)*dcos(glong)*dn - dsin(glong)*dea +
     *  dcos(phi)*dcos(glong)*du
               offy= -dsin(phi)*dsin(glong)*dn  +dcos(glong)*dea +
     *  dcos(phi)*dsin(glong)*du
               offz=  dcos(phi)*dn + dsin(phi)*du

c            print *, tg, phi, glong
c            print *, rx, offx
c            print *, ry, offy
c            print *, rz, offz


               RX = RX + offx
               RY = RY + offy
               RZ = RZ + offz


            end if

         end if
c-  OT
         RX = RX + VX*TS
         RY = RY + VY*TS
         RZ = RZ + VZ*TS

         if (i.le.5.or.(i.ge.200.and.i.le.206)) write(102,1021) i,STAT, 
     *          RX, RY, RZ, VX, VY, VZ,
     *          VX*TS, VY*TS, VZ*TS, TS, epoch
1021     format (i5,2x,a8, 11(2x,f16.6))
           

         CALL TRANSF (RX,RY,RZ,PHI,GLONG,ELHGT,R0,0)

C    NOW THREE CALLS TO CORRECTIONS ROUTINES IN ORDER TO COMPUTE THEM
C  ALSO ONE SECOND AFTER AND BEFORE OF OBSERVATION TIME.

         DO 154 IJ=1,3

            IF (IJ.EQ.1) UTJ = UTJ0 + (1.D0/86400.D0)*TWOPI
            IF (IJ.EQ.2) UTJ = UTJ0 - (1.D0/86400.D0)*TWOPI
            IF (IJ.EQ.3) UTJ = UTJ0

C---------------------------------------------------------------
C   GET POLAR WOBBLE
C---------------------------------------------------------------

            XJM = XPOLE(IJ)
            YJM = YPOLE(IJ)

            DO  IK = 1,3

               V_after_PN(IK) = DVELG(IK,3)      !   Geocentric velocities of the Sun
               R_after_PN(IK) = DCORG(IK,3)      !   Geocentric positions of the Sun
               XLUN_after_PN(IK) = CMOON(IK,3)       !   Geocentric positions of the Moon

            END DO

c + OT    04.06.2010

C----------------------------------------------------------------
C     APPLY POLE COORDINATES CORRECTION AND COMPUTE DERIVATIVES
C  OF COORDINATES WITH RESPECT TO X, Y WOBBLE and UT1
C----------------------------------------------------------------

    

            CALL POLROT(XPOLE(IJ),YPOLE(IJ),EROT,EDROT,IJ) !  Matrix EROT=R1(y)*R2(x)

 
             rr3(1,1) = DCOS(HSG(IJ))
             rr3(2,1) = DSIN(HSG(IJ))
             rr3(3,1) = 0.d0
             rr3(1,2) = -DSIN(HSG(IJ))
             rr3(2,2) = DCOS(HSG(IJ))
             rr3(3,2) = 0.d0
             rr3(1,3) = 0.d0
             rr3(2,3) = 0.d0
             rr3(3,3) = 1.d0
 
             do ii=1,3
                do j=1,3
                   r4(ii,j) = 0.d0
                   do k=1,3
                      r4(ii,j ) = r4(ii,j) + rr3(ii,k) * erot(k,j)    ! R3*R1*R2
                   end do
                end do
             end do
 
          

c - OT    04.06.2010

C----------------------------------------------------------------------
C    CALL ROUTINE CPOSIT TO COMPUTE CORRECTIONS TO STATION COORDINATES
C  DUE TO EARTH TIDES, OCEAN & ATMOSPHERIC LOADING AND POLE TIDE
c+ jboehm
c  and due to snow loading and soil moisture
c- jboehm
C    APPLY CORRECTION TO STATION POSITION
C
C---------------------------------------------------------------------

c+  OT 04.06.2010
c
c            CALL CPOSIT (ij,ISTCAT,HSG(IJ),pres)  ! Requires geocentic positions of the Sun and Moon

            CALL CPOSIT (ij,ISTCAT,HSG(IJ),
     *       R_after_PN, V_after_PN, XLUN_after_PN,
     *       r4,CTIDE,pres)  ! Requires geocentic positions of the Sun and Moon
c
c-  OT 04.06.2010

          
             CRX0 = RX + CTIDE(1)
             CRY0 = RY + CTIDE(2)
             CRZ0 = RZ + CTIDE(3)

             V2(1) = - OMEGA * CRY0
             V2(2) =   OMEGA * CRX0
             V2(3) =   0.D0

             trf(1) = crx0
             trf(2) = cry0 
             trf(3) = crz0  


c
c    Rotation R1(y)*R2(x)*TRS
c



            CRX(IJ) = EROT(1,1)*CRX0 + EROT(1,2)*CRY0 + EROT(1,3)*CRZ0
            CRY(IJ) = EROT(2,1)*CRX0 + EROT(2,2)*CRY0 + EROT(2,3)*CRZ0
            CRZ(IJ) = EROT(3,1)*CRX0 + EROT(3,2)*CRY0 + EROT(3,3)*CRZ0


c
c    Rotation R1(y)*R2(x)*TRS, but without tidal correction
c



            CRX4(IJ) = EROT(1,1)* RX + EROT(1,2)* RY + EROT(1,3)* RZ
            CRY4(IJ) = EROT(2,1)* RX + EROT(2,2)* RY + EROT(2,3)* RZ
            CRZ4(IJ) = EROT(3,1)* RX + EROT(3,2)* RY + EROT(3,3)* RZ


c             CRX4(IJ) = RX
c             CRY4(IJ) = RY
c             CRZ4(IJ) = RZ



c 
c           Transformation of velocities Rotation R1(y)*R2(x)* V2
c
           if (ij.eq.3) then

              CRVX = EROT(1,1)*V2(1) + EROT(1,2)*V2(2) + EROT(1,3)*V2(3)
              CRVY = EROT(2,1)*V2(1) + EROT(2,2)*V2(2) + EROT(2,3)*V2(3)
              CRVZ = EROT(3,1)*V2(1) + EROT(3,2)*V2(2) + EROT(3,3)*V2(3)

           end if  


c
c   Spherical coordinates for true system of Equinox
c

            CALL TRANSF (CRX(IJ),CRY(IJ),CRZ(IJ),CPHI2(ij),CGLONG2(ij),
     *                   CELHGT2(ij),CR02(ij),0)

            

            IF (IJ.EQ.1) THEN   !  One second before

               DPDXW1 = EDROT(1,1,1)*CRX0 +                       !
     *                    EDROT(1,2,1)*CRY0 + EDROT(1,3,1)*CRZ0   !  R1(y)*dR2(X)/dx*TRS
               DPDXW2 = EDROT(2,1,1)*CRX0 +
     *                    EDROT(2,2,1)*CRY0 + EDROT(2,3,1)*CRZ0
               DSDXW(1,3) = EDROT(3,1,1)*CRX0 +
     *                    EDROT(3,2,1)*CRY0 + EDROT(3,3,1)*CRZ0
               DPDYW1 = EDROT(1,1,2)*CRX0 +                       !
     *                    EDROT(1,2,2)*CRY0 + EDROT(1,3,2)*CRZ0   !  dR1(y)/dy*R2(X)*TRS
               DPDYW2 = EDROT(2,1,2)*CRX0 +
     *                    EDROT(2,2,2)*CRY0 + EDROT(2,3,2)*CRZ0
               DSDYW(1,3) = EDROT(3,1,2)*CRX0 +
     *                    EDROT(3,2,2)*CRY0 + EDROT(3,3,2)*CRZ0

               DSDXW(1,1) = DPDXW1*DCOS(HSG(1)) - DPDXW2*DSIN(HSG(1))
               DSDXW(1,2) = DPDXW1*DSIN(HSG(1)) + DPDXW2*DCOS(HSG(1))

               DSDYW(1,1) = DPDYW1*DCOS(HSG(1)) - DPDYW2*DSIN(HSG(1))
               DSDYW(1,2) = DPDYW1*DSIN(HSG(1)) + DPDYW2*DCOS(HSG(1))

               DSDUT(1,1) = -CRX(IJ)*DSIN(HSG(1)) - CRY(IJ)*DCOS(HSG(1))
               DSDUT(1,2) = CRX(IJ)*DCOS(HSG(1)) - CRY(IJ)*DSIN(HSG(1))



C  PUT RATIO BETWEEN SIDEREAL TIME AND UT1 IN THE PARTIALS

               CSID  = 1.002737909350795D0


               DSDUT(1,1) = DSDUT(1,1) * CSID
               DSDUT(1,2) = DSDUT(1,2) * CSID
               DSDUT(1,3) = 0.D0

            ENDIF

            IF (IJ.EQ.2) THEN  !  One second after

               DPDXW1 = EDROT(1,1,1)*CRX0 +
     *                    EDROT(1,2,1)*CRY0 + EDROT(1,3,1)*CRZ0
               DPDXW2 = EDROT(2,1,1)*CRX0 +
     *                    EDROT(2,2,1)*CRY0 + EDROT(2,3,1)*CRZ0
               DSDXW(2,3) = EDROT(3,1,1)*CRX0 +
     *                    EDROT(3,2,1)*CRY0 + EDROT(3,3,1)*CRZ0
               DPDYW1 = EDROT(1,1,2)*CRX0 +
     *                    EDROT(1,2,2)*CRY0 + EDROT(1,3,2)*CRZ0
               DPDYW2 = EDROT(2,1,2)*CRX0 +
     *                    EDROT(2,2,2)*CRY0 + EDROT(2,3,2)*CRZ0
               DSDYW(2,3) = EDROT(3,1,2)*CRX0 +
     *                    EDROT(3,2,2)*CRY0 + EDROT(3,3,2)*CRZ0

               DSDXW(2,1) = DPDXW1*DCOS(HSG(2)) - DPDXW2*DSIN(HSG(2))
               DSDXW(2,2) = DPDXW1*DSIN(HSG(2)) + DPDXW2*DCOS(HSG(2))

               DSDYW(2,1) = DPDYW1*DCOS(HSG(2)) - DPDYW2*DSIN(HSG(2))
               DSDYW(2,2) = DPDYW1*DSIN(HSG(2)) + DPDYW2*DCOS(HSG(2))

               DSDUT(2,1) = -CRX(IJ)*DSIN(HSG(2)) - CRY(IJ)*DCOS(HSG(2))
               DSDUT(2,2) = CRX(IJ)*DCOS(HSG(2)) - CRY(IJ)*DSIN(HSG(2))

C  PUT RATIO BETWEEN SIDEREAL TIME AND UT1 IN THE PARTIALS

               CSID  = 1.002737909350795D0

               DSDUT(2,1) = DSDUT(2,1) * CSID
               DSDUT(2,2) = DSDUT(2,2) * CSID
               DSDUT(2,3) = 0.D0

            ENDIF

            IF (IJ.EQ.3) THEN   ! in time

               DPDXW1 = EDROT(1,1,1)*CRX0 +                        !
     *                    EDROT(1,2,1)*CRY0 + EDROT(1,3,1)*CRZ0    !  R1(y)*dR2(X)/dx
               DPDXW2 = EDROT(2,1,1)*CRX0 +
     *                    EDROT(2,2,1)*CRY0 + EDROT(2,3,1)*CRZ0
               DSDXW(3,3) = EDROT(3,1,1)*CRX0 +
     *                    EDROT(3,2,1)*CRY0 + EDROT(3,3,1)*CRZ0
               DPDYW1 = EDROT(1,1,2)*CRX0 +                        !
     *                    EDROT(1,2,2)*CRY0 + EDROT(1,3,2)*CRZ0    !  dR1(y)/dy*R2(X)
               DPDYW2 = EDROT(2,1,2)*CRX0 +
     *                    EDROT(2,2,2)*CRY0 + EDROT(2,3,2)*CRZ0
               DSDYW(3,3) = EDROT(3,1,2)*CRX0 +
     *                    EDROT(3,2,2)*CRY0 + EDROT(3,3,2)*CRZ0

c
c             R3*(-H)*R1(y)*dR2(x)/dx * TRS
c
               DSDXW(3,1) = DPDXW1*DCOS(HSG(3)) - DPDXW2*DSIN(HSG(3))
               DSDXW(3,2) = DPDXW1*DSIN(HSG(3)) + DPDXW2*DCOS(HSG(3))

c
c             R3*(-H)*dR1(y)/dy*R2(x) * TRS
c
               DSDYW(3,1) = DPDYW1*DCOS(HSG(3)) - DPDYW2*DSIN(HSG(3))
               DSDYW(3,2) = DPDYW1*DSIN(HSG(3)) + DPDYW2*DCOS(HSG(3))

c
c             dR3*(-H)/dH*R1(y)*R2(x) * TRS
c
               DSDUT(3,1) = -CRX(IJ)*DSIN(HSG(3)) - CRY(IJ)*DCOS(HSG(3))
               DSDUT(3,2) = CRX(IJ)*DCOS(HSG(3)) - CRY(IJ)*DSIN(HSG(3))

C  PUT RATIO BETWEEN SIDEREAL TIME AND UT1 IN THE PARTIALS

               CSID  = 1.002737909350795D0

               DSDUT(3,1) = DSDUT(3,1) * CSID
               DSDUT(3,2) = DSDUT(3,2) * CSID
               DSDUT(3,3) = 0.D0


c for delay rate

c
c             R3*(-dH/dt)*R1(y)*dR2(x)/dx * TRS
c
               DSDXW_t(1) = -DPDXW1*DSIN(HSG(3)) - DPDXW2*DCOS(HSG(3))
               DSDXW_t(2) = DPDXW1*DCOS(HSG(3)) - DPDXW2*DSIN(HSG(3))
               DSDXW_t(3) = 0.d0

c
c             R3*(-dH/dt)*dR1(y)/dy*R2(x) * TRS
c
               DSDYW_t(1) = -DPDYW1*DSIN(HSG(3)) - DPDYW2*DCOS(HSG(3))
               DSDYW_t(2) = DPDYW1*DCOS(HSG(3)) - DPDYW2*DSIN(HSG(3))
               DSDYW_t(3) = 0.d0

c
c             dR3*(-dH/dt)/dH*R1(y)*R2(x) * TRS
c
               DSDUT_t(1) = -CRX(IJ)*DCOS(HSG(3)) + CRY(IJ)*DSIN(HSG(3))
               DSDUT_t(2) = -CRX(IJ)*DSIN(HSG(3)) - CRY(IJ)*DCOS(HSG(3))
               DSDUT_t(3) = 0.d0 

c
c             To estimate UT1-UTC with delay rate
c
c             dR3*(-H)/dH*R1(y)*R2(x) * V2
c

               DSDUT_V(1) = -CRVX*DSIN(HSG(3)) - CRVY*DCOS(HSG(3))
               DSDUT_V(2) = CRVX*DCOS(HSG(3)) - CRVY*DSIN(HSG(3))

C  PUT RATIO BETWEEN SIDEREAL TIME AND UT1 IN THE PARTIALS

               CSID  = 1.002737909350795D0

               DSDUT_V(1) = DSDUT_V(1) * CSID
               DSDUT_V(2) = DSDUT_V(2) * CSID
               DSDUT_V(3) = 0.D0

            ENDIF


c
c    Rotation R3(-H)*R1(y)*R2(x)*TRS
c

            
            CRX2(1) = CRX(IJ)*DCOS(HSG(IJ)) - CRY(IJ)*DSIN(HSG(IJ))
            CRX2(2) = CRX(IJ)*DSIN(HSG(IJ)) + CRY(IJ)*DCOS(HSG(IJ))
            CRX2(3) = CRZ(IJ)

c
c    Rotation R3(-H)*R1(y)*R2(x)*TRS, but without tidal effect 
c



            CRX5(1) = CRX4(IJ)*DCOS(HSG(IJ)) - CRY4(IJ)*DSIN(HSG(IJ))
            CRX5(2) = CRX4(IJ)*DSIN(HSG(IJ)) + CRY4(IJ)*DCOS(HSG(IJ))
            CRX5(3) = CRZ4(IJ)
         
            

c            CRX6(ij) = CRX4(IJ)*DCOS(HSG(IJ)) - CRY4(IJ)*DSIN(HSG(IJ))
c            CRY6(ij) = CRX4(IJ)*DSIN(HSG(IJ)) + CRY4(IJ)*DCOS(HSG(IJ))
c            CRZ6(ij) = CRZ4(IJ)

c
c    CPHI2, CGLONG2, CR02 - not corrected for precession and nutation !!!
c    to be used for for apparent positions, not for reductions and partials
c

c
c    Rotation P*N*R3(-H)*R1(y)*R2(x)*TRS
c
c    First of all S = P * N - first step
c
            do ii=1,3

               do j=1,3

                  SB(ii,j)=0.d0
                  SA(ii,j)=0.d0
                  SS(ii,j)=0.d0

                  do k=1,3

                     SB(ii,j) = SB(ii,j) + PRECB(ii,k)*SNB(k,j)  ! P*N
                     SA(ii,j) = SA(ii,j) + PRECA(ii,k)*SNA(k,j)  ! P*N
                     SS(ii,j) = SS(ii,j) + PREC(ii,k)*SN(k,j)    ! P*N

                  end do
               end do
            end do


c
c    Rotation P*N*R3(-H)*R1(y)*R2(x)*TRS
c
c    and P*N*R3(-H)*R1(y)*R2(x)*V2 for ij = 3
c
c    Then     S*R3(-H)*R1(y)*R2(x)*TRS  - second step
c
            do ii=1,3

               
               if (ij.eq.1) then

                  pp(ii) = 0.d0
                  p4(ii) = 0.d0

c                  if (i.eq.1) print *, ij, crx2(1), crx2(2)

                  do j=1,3

                     pp(ii) = pp(ii) + sb(ii,j)*crx2(j)
                     p4(ii) = p4(ii) + sb(ii,j)*crx5(j)

                  end do

                  crx1(1) = pp(1)
                  cry1(1) = pp(2)
                  crz1(1) = pp(3)

                  crx6(1) = p4(1)
                  cry6(1) = p4(2)
                  crz6(1) = p4(3)


               end if

               if (ij.eq.2) then

                  pp(ii) = 0.d0
                  p4(ii) = 0.d0

                  do j=1,3

                     pp(ii) = pp(ii) + sa(ii,j)*crx2(j)
                     p4(ii) = p4(ii) + sa(ii,j)*crx5(j)


                  end do

                  crx1(2) = pp(1)
                  cry1(2) = pp(2)
                  crz1(2) = pp(3)


                  crx6(2) = p4(1)
                  cry6(2) = p4(2)
                  crz6(2) = p4(3)

               end if

               if (ij.eq.3) then

                  pp(ii) = 0.d0
                  p4(ii) = 0.d0

                  pv(ii) = 0.d0 

c                  if (i.eq.1) print *, ij, crx2(1), crx2(2)

                  do j=1,3

                     pp(ii) = pp(ii) + ss(ii,j)*crx2(j)
                     p4(ii) = p4(ii) + ss(ii,j)*crx5(j)

c
c To estimate UT1-UTC with delay rate
c
c for velocity V2  
C                       P*N * dR3*(-H)/dH*R1(y)*R2(x) * V2
c


                     pv(ii) = pv(ii) + ss(ii,j) * dsdut_v(j)

                  end do

                  crx1(3) = pp(1)
                  cry1(3) = pp(2)
                  crz1(3) = pp(3)

                  crx6(3) = p4(1)
                  cry6(3) = p4(2)
                  crz6(3) = p4(3)


               end if

            end do

            CALL TRANSF (CRX1(IJ),CRY1(IJ),CRZ1(IJ),CPHI(IJ),CGLONG(IJ),
     *                   CELHGT(IJ),CR0(IJ),0)


c
c    CPHI, CGLONG, CR0 - corrected for precession and nutation
c    to be used for reductions and partials, not for apparent positions
c

            SD = DSIN(DE(ij))
            CD = DCOS(DE(ij))


C  GET LOCAL HOUR ANGLE OF SOURCE AND SIDEREAL TIME


            HSTA(IJ) = HSG(IJ) + CGLONG(IJ)
            HA(IJ)   = H(IJ)   + CGLONG(IJ)
            HA2(IJ)   = H(IJ)   + CGLONG2(IJ)  !  for apparent positions


            IF (HA(IJ) .LT. 0.D0)  HA(IJ) = HA(IJ) + TWOPI
            IF (HA2(IJ) .LT. 0.D0)  HA2(IJ) = HA2(IJ) + TWOPI  !  for apparent positions

            CHA=DCOS(HA(IJ))
            SHA=DSIN(HA(IJ))

            CHA2=DCOS(HA2(IJ))       !  for apparent positions
            SHA2=DSIN(HA2(IJ))       !  for apparent positions


C-----------------------------------------------------------------------
C COMPUTE COSZ, SINZ AND THE ELEVATIONS
C The parameters are to be corrected for aberration - done in NUTPAR.for
C
C-----------------------------------------------------------------------


            SLT = DSIN(CPHI2(IJ))
            CLT = DCOS(CPHI2(IJ))

            CZ1 = SLT*SD + CLT*CD*CHA2
            SZ1 = DSQRT(1.D0 - CZ1 * CZ1)

            Z(IJ) = DATAN2(SZ1,CZ1)            !       Zenith distance

            ELEV = PI/2.d0-Z(IJ)             !       Elevation

C   COMPUTE CORRECTION TO ZENITH DISTANCE DUE TO TROPOSPHERIC REFRACTION
C  CORZ = REFRACT. FACTOR / TANGENT OF ELEVATION
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc TESMER 23.02.05 provide trop refraction correction as option cc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


            if (YESTRC.EQ.'N'.OR.YESTRC.EQ.'n') then

               corz(ij) = 0.d0
               if (i.eq.1) print*, 'corz is 0 now !!!'

            else

               CORZ(IJ) = 313.D-6 / DTAN(PI/2.D0 - Z(IJ))
              


               IF (Z(IJ).GT.PI/2.D0) CORZ(IJ) = -CORZ(IJ)

               if(i.eq.1.and.ij.eq.3) write (12,590)
 590           format (' Correction for tropospheric refraction is
     * implemented')

            end if
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc TESMER 23.02.05 provide trop refraction correction as option cc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

            CAZ(IJ) = (SD-CZ1*SLT)/(SZ1*CLT)
            SAZ(IJ) = -CD*SHA2/SZ1

            AZIMUTH(IJ) = DATAN2(SAZ(IJ),CAZ(IJ))       !   Azimuth


c
c  For troposphere gradient partial derivatives
c
            if (ij.eq.3) then

               a_grad(1) = dcos(azimuth(ij)) / dtan(elev)   !  North-South
               a_grad(2) = dsin(azimuth(ij)) / dtan(elev)   !  East-West


            end if


154      CONTINUE

     
        
c         w4(1) = (crx6(1) - crx6(2))/2.d0 
c         w4(2) = (cry6(1) - cry6(2))/2.d0
c         w4(3) = (crz6(1) - crz6(2))/2.d0 

      
         w4(1) = (crx1(1) - crx1(2))/2.d0 
         w4(2) = (cry1(1) - cry1(2))/2.d0
         w4(3) = (crz1(1) - crz1(2))/2.d0 




         vel_p(1) = crx1(1) - crx1(3)  ! velocity +0.5 sec
         vel_p(2) = cry1(1) - cry1(3)  
         vel_p(3) = crz1(1) - crz1(3)  



         vel_m(1) = crx1(3) - crx1(2)  ! velocity -0.5 sec
         vel_m(2) = cry1(3) - cry1(2)  
         vel_m(3) = crz1(3) - crz1(2)  


         a4(1) = vel_p(1) - vel_m(1)  ! acceleration
         a4(2) = vel_p(2) - vel_m(2)  
         a4(3) = vel_p(3) - vel_m(3)  


ccc         IF (IJ.EQ.3) THEN   ! in time
c               
c               DSDXW(3,3) = EDROT(3,1,1)*w4(1) +
c     *                    EDROT(3,2,1)*w4(2) + EDROT(3,3,1)*w4(3)
c               DSDYW(3,3) = EDROT(3,1,2)*w4(1) +
c     *                    EDROT(3,2,2)*w4(2) + EDROT(3,3,2)*w4(3)
c
c               DSDUT(3,3) = 0.D0
c
ccc         ENDIF


    
C----------------------------------------------------------------------
C   SIDGRE IS THE GREENWICH SIDEREAL TIME AT THE MOMENT
C   SIDGRE_dot is its partial
C-----------------------------------------------------------------------

         SIDGRE = HSG(3)
         SIDGRE_dot = (HSG(2) - HSG(1))/2.d0

C----------------------------------------------------------
C    CALL FOR THE CORRECTIONS DUE TO TROPOSPHERICAL DELAY
C----------------------------------------------------------

c+ jboehm also istcat, isour

         CALL CDELAY (ZDRY, istcat,isour, 
     *   iav_imf, iav_rmf, iav_vmf,
     *   iav_vm1, iav_grad, YESMAR,YESAZD, YESAGR)


c
c  PRES in CDELAY to be fixed for consistency to VMF1 data
c

c- jboehm

C----------------------------------------------------------
C    CALL FOR THE CORRECTIONS DUE TO THERMAL DEFORMATIONS OF ANTENNAS
C    (Haas, et al., 1998; Skurihina, 2000)
C----------------------------------------------------------

         IF (opttd) THEN

c+   otitov 22.03.2006
            CALL thermdef (STAT,AXTYP,OFFS,opttd1,CAZ,SAZ,Z,
     *       temp,corz,therm_d,istcat)      !  Change on 04.09.2006
c-   otitov 22.03.2006

         ENDIF

C----------------------------------------------------------
C    CALL FOR THE CORRECTIONS DUE TO AXIS OFFSET
C----------------------------------------------------------

         CALL AXIS (STAT,AXTYP,OFFS,CAZ,SAZ,CORZ,CPHI2,Z,DE,axkt)



         CALL FEED (STAT, AXTYP, CORZ, HSTA, CPHI2,Z,DE, FEED_R)

c         print *, FEED_R

c          axkt(2)  = axkt(2) + FEED_R * 1.d9 !  feed from rad to nrad

C----------------------------------------------------------
C   CALL FOR THE PARTIAL DERIVATIVES
C----------------------------------------------------------

c
c   Change on 3.11.2005; then 02.02.2007 (HSTA array added)
c

c+ OT   04.06.2010

         CALL PARTIA (RQU,HSTA,RBAR,VBAR,CRX1,CRY1,CRZ1,SN,PREC,SS,
     *    RA2000,DE2000,EROT,STAT,SIDGRE_dot,dsdxw_t,dsdyw_t,dsdut_t,
     *    ppolar_rate)  ! 29-Oct-2019

          
c         CALL PARTIA (RQU,HSTA,VBAR,RX,RY,RZ,SN,PREC,SS,
c     *    RA2000,DE2000,EROT,STAT)
c- OT   04.06.2010

c         write (13,*) 'ppolar = ', ppolar_rate(3), ' pv = ', pv(3)


c
c        for calculation of the partial derivative of delay rate on UT1-UTC
c

cc         print *, pv(1), pv(2), pv(3)

         dr_ut1 = -dotpr(rqu,pv)/c



C----------------------------------------------------------
C   WRITE ON FILE STATIM THE MODIFIED COORDINATES, THE STATION
C  ANGLES, MAPPING FUNCTIONS AND PARTIAL DERIVATIVES
C----------------------------------------------------------

      WRITE (18,REC=I) JSTA,ISTCAT,ISOUR,TEMP,PRES,HUMD,IFAC,CABL,
c
c  Change on 4.11.2005. Station coordinates are in the mean system of equinox
c
     *                 CRX1,CRY1,CRZ1,RXX,RYY,RZZ,         !  01.09.2006
     *                 CGLONG,CPHI,CELHGT,
     *      Z(1),HA(1),HSTA(1),Z(2),HA(2),HSTA(2),Z(3),HA(3),HSTA(3),
     *                 ZDRY,                   
     *                 cmniel_w,cmniel_d,
     *                 w4(1), w4(2), w4(3),           !  08.03.2017 OT - geocentric velocity of station
     *                 a4(1), a4(2), a4(3),           !  13.12.2018 OT - geocentric velocity of station 
     *            AXKT,
     *            PCOORD,PCOORR,PNUTAT,PPOLAR,
     *            ppolar_rate,                        !  03.01.2020
     *             therm_d,azimuth, a_grad,
     *            dr_ut1,                             !  01.01.2021  - OT, partial derivative of delay rate on UT1-UTC
c+ jboehm
     *            z200,smfw3,azim,elev,
     *            dimfh,dimfw,   !  dry and wet IMF maping functions
     *            gimfh,         !  dry gradient IMF mapping function
     *            grn,gre,       !  apriori north and east gradients
     *            iavail(istcat),
     *            drmfh,drmfw,    !  dry and wet ray traced mapping function
     *            iavailr(istcat),
     *            dqmfh,dqmfw,    !  dry and wet Vienna mapping function
     *            iavailq(istcat)

c- jboehm
C+ jboehm, 2005 Sep 5
     *            ,gmfh,gmfw            !  Global Mapping Function GMF
     *            ,grnh,greh,grnw,grew  ! apriori gradients
     *            ,iavailv(istcat)      ! VMF1
     *            ,vmf1h,vmf1w
     *            ,vmf1hd,vmf1wd      !  VMF1 for delay rate 5-Nov-2019 (OT)

c- jboehm, 2005-Oct-28
     


3     CONTINUE                 !          end of third cycle


c+ jboehm
c      close (11)
c- jboehm

      GO TO 655

653   WRITE (*,654)
654   FORMAT(///' ***** ERROR OPENING FILE LOADING.OCE *****'//
     *' PLEASE BE SURE THAT THE FILE EXISTS AND IS IN THIS DIRECTORY'/
     *' BEFORE RUNNING THIS PROGRAM AGAIN'/' HIT ENTER TO QUIT'//////)
      READ (*,'(I1)') IDUMMY

673   WRITE (*,674)
674   FORMAT(//' ***** ERROR OPENING FILE THERMAL.DEF *****'//
     *' PLEASE BE SURE THAT THE FILE EXISTS AND IS IN THIS DIRECTORY'/
     *' BEFORE RUNNING THIS PROGRAM AGAIN'/' HIT ENTER TO QUIT'//////)
      READ (*,'(I1)') IDUMMY
                                            
655   CONTINUE

      CLOSE (13)
      CLOSE (17)
      CLOSE (18)
      CLOSE (19)
      CLOSE (20)
      STOP
      END


C************************************************************************

c+  otitov   22.03.2006
      SUBROUTINE thermdef (STAT,AXTYP,OFFS,opttd1,CAZ,SAZ,Z,
     *   TEMP,corz,therm_d,ist)       !  change on 04.09.2006 - OT
c      SUBROUTINE thermdef (STAT,AXTYP,OFFS,CAZ,SAZ,TEMP,corz,ist)
c-  otitov   22.03.2006

C**********************************************************************
C
C   THIS ROUTINE COMPUTES THE CORRECTIONS TO THE DELAY AND RATE
C  DUE TO THERMAL DEFORMATION OF THE VLBI ANTENNA
C
C   REVISION 1999 AUGUST  4   (O.T.)
C   Revision TESMER 2000 August 9, in case of HOBART-ANTENNA X-Y2,
C					don't stop
C   Revision     tesmer 2000.11.20 cope with differnet focus
C   REVISION 2001 MAY    04 (O.T.)  X-YN ,X-YE and RICHMOND types are added
C   REVISION 2004 August  19 (O.T.)  Corrected bug for the Lena Skurihina
C                     correction term
C   REVISION 2006 March 22 (O.T.)  Correction for the Lena Skurihina's
C                                    term is made as option
C   REVISION 2006 August 28 (O.T.)  /COCDEL/ is off;
C                            Array Z is in SUBROUTINE line
C   REVISION 2006 SEPTEMBER 04 (O.T.) COMMON-BLOCK /CORR/ deleted
C                      IK=3 deleted; therm_d in subroutine line
C   LAST REVISION 2006 SEPTEMBER 18 (O.T.) Sign was changed to 'PLUS'
c   followed by the corresponding change in 'KVLBI.FOR' and 'COLL.FOR'
c
c*******************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INCLUDE 'OCCAM_N.FI'

      CHARACTER AXTYP*4,STAT*8,STT*1
      logical opttd1
      DIMENSION Z(3),HA(3),DE(3),CAZ(3),SAZ(3),CORZ(3)

      COMMON /THERM/ ga(nstations),hd(nstations),gf(nstations),
     * hf(nstations),hp(nstations),hs(nstations),hv(nstations),
     * ifo(nstations),ref_temp(nstations)
      COMMON /PHYS/ C, FL, AH, AU, J20

c
c     For delay calculation only IK=3 is necessary
c     For delay rate calculation only IK=1,2 is necessary
c

c      ik = 3  ! Inserted on 3.11.2005  deleted on 04.09.2006


      tempin=temp

c cope with reference temperature
C      temp0=20.d0                  !  reference temperature


      if (ref_temp(ist).ne.999.d0) then
        temp0=ref_temp(ist)
      else
        temp0=temp
      endif

c if no temperature is measured => temp is ref_temp
c => correction is zero
      if (temp.eq.-999.) then
         temp=temp0
      end if

c COPE WITH DIFFERENT FOCUS OF ANTENNAS
      focfac=1.8d0
c if antenna is prime focus (ifo=10) => focfac=0.9
      if(ifo(ist).eq.10) focfac=0.9d0
c if antenna is secondary focus (ifo=20) => focfac=1.8
      if(ifo(ist).eq.20) focfac=1.8d0

      SZ = DSIN (Z(3)-CORZ(3))        !    sin(zenith)
      CZ = DCOS (Z(3)-CORZ(3))        !    cos(zenith)
      CD = DCOS (DE(3))                !    cos(declination)
      SD = DSIN (DE(3))                !    sin(declination)

C
C  TREATING THE CASE OF AN AZIMUTH-ELEVATION MOUNTING
C----------------------------------------------------------

      IF (AXTYP.EQ.'AZEL')  THEN

         THERM_D =  ( gf(ist)*(temp-temp0)*hf(ist)*cz +
     *     ga(ist)*(temp-temp0)*(hp(ist)*cz+hv(ist)-focfac*hs(ist)) ) /C

c+  otitov   22.03.2006
           corr2 =  ga(ist)*(temp-temp0)*OFFS*sz /C        ! Elena Skurihina's term
c-  otitov   22.03.2006

C
C  TREATING THE CASE OF AN EQUATORIAL ANTENNA
C----------------------------------------------------------

         ELSE IF (AXTYP.EQ.'EQUA') THEN

             THERM_D =  ( gf(ist)*(temp-temp0)*hf(ist)*cz +
     *       ga(ist)*(temp-temp0)*(hp(ist)*cz+hv(ist)-focfac*hs(ist)+
     *       hd(ist)*cd) ) /C

c+  otitov   22.03.2006
             corr2 =  ga(ist)*(temp-temp0)*OFFS*cd /C        ! Elena Skurihina's term
c-  otitov   22.03.2006

C
C  TREATING THE CASE OF AN X/Y-fixed N-S axis ANTENNA
C----------------------------------------------------------


         ELSE IF (AXTYP.EQ.'X-Y1'.OR.AXTYP.EQ.'X-YN') THEN

            if ((1.d0-(sz*caz(3))**2).lt.0.d0) print *, 'XY1'
               THERM_D =  ( gf(ist)*(temp-temp0)*hf(ist)*cz +
     *                  ga(ist)*(temp-temp0)*(hp(ist)*cz+
     *                  hv(ist)-focfac*hs(ist)) ) /C
c+  otitov   22.03.2006
               corr2 =  ga(ist)*(temp-temp0)*OFFS
     *           *dsqrt(1.d0-(sz*caz(3))**2) /C    ! Elena Skurihina's term
c-  otitov   22.03.2006

C
C  TREATING THE CASE OF AN X/Y-fixed E-W axis ANTENNA
C----------------------------------------------------------

         ELSE IF (AXTYP.EQ.'X-Y2'.OR.AXTYP.EQ.'X-YE') THEN

            if ((1.d0-(sz*saz(3))**2).lt.0.d0) print *, 'XY2'
            THERM_D =  ( gf(ist)*(temp-temp0)*hf(ist)*cz +
     *                  ga(ist)*(temp-temp0)*(hp(ist)*cz+
     *                  hv(ist)-focfac*hs(ist)) ) /C
c+  otitov   22.03.2006
               corr2 =   ga(ist)*(temp-temp0)*OFFS*
     *         dsqrt(1.d0-(sz*saz(3))**2) /C    ! Elena Skurihina's term
c-  otitov   22.03.2006

C
C  TREATING other ANTENNA types
C----------------------------------------------------------

      ELSE

         WRITE (*,321) AXTYP,STAT
321   FORMAT (//' ***** WARNING *****  AXIS MODEL ',A4,' FOR STATION',
     *1X,A8,' NOT AVAILABLE'//' HIT "Q" TO QUIT, "C" TO CONTINUE'//)
            READ (*,'(A1)') STT
            IF (STT.NE.'C'.AND.STT.NE.'c') STOP

      ENDIF

C AXIS MODEL FOR RICHMOND   SHOULD BE APPLIED FOR !!!!
C-----------------------------
      IF (STAT.EQ.'RICHMOND')  THEN

         coeff = 1.D0 - (CZ*DSIN(.6817256D0)+
     *    SZ*DCOS(.6817256D0)*(CAZ(3)*DCOS(.0020944D0)-SAZ(3)*
     *      DSIN(.0020944D0)))**2

         if (coeff.lt.0.d0) print *, 'richmond'


         THERM_D =  ( gf(ist)*(temp-temp0)*hf(ist)*cz +
     *               ga(ist)*(temp-temp0)*(hp(ist)*cz+
     *               hv(ist)-focfac*hs(ist)) ) /C
c+  otitov   22.03.2006
            corr2 =  ga(ist)*(temp-temp0)*OFFS*DSQRT(coeff)/c      ! Elena Skurihina's term
c-  otitov   22.03.2006

      END IF

c+  otitov   22.03.2006

      if (opttd1) then
c        if (corr2.ne.0) print 77,stat,ga(ist),temp,temp0,sz,offs,c*corr2
c 77      format (a8,6(2x,f8.5))
         THERM_D = THERM_D - corr2   !  Last check on 19.09.2006 with therm = therm1 - therm2 later
      end if

c-  otitov   22.03.2006

      RETURN
      END



c      SUBROUTINE AXIS (STAT,AXTYP,OFFS,CAZ,SAZ,CORZ,azimuth)
      SUBROUTINE AXIS (STAT,AXTYP,OFFS,CAZ,SAZ,CORZ,CPHI,Z,DE,axkt)

C************************************************************************
C
C   THIS ROUTINE COMPUTES THE CORRECTIONS TO THE DELAY AND RATE
C  DUE TO AXIS OFFSET OF THE ANTENNA
C
C        REVISION 1990 NOVEMBRE 19 (N.Z.L.)
C        REVISION 2001 MARCH    04 (O.T.) X-YN and X-YE types are added
C   LAST REVISION 2006 SEPTEMBER 04 (O.T.) COMMON-BLOCK /CORR/ deleted
C
C************************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      CHARACTER AXTYP*4,STAT*8,STT*1
      DIMENSION Z(3),DE(3), CPHI(3)
      DIMENSION CAZ(3),SAZ(3),CORZ(3),AXKT(2)
      COMMON /PHYS/ C, FL, AH, AU, J20
      COMMON /MATH/ PI, TWOPI, ROG, ROS, ROZ, TOL

c
c     For delay calculation only IK=3 is necessary
c     For delay rate calculation only IK=1,2 is necessary
c



      DO IK=1,3

         SZ = DSIN (Z(IK)-CORZ(IK))
         CZ = DCOS (Z(IK)-CORZ(IK))
         CD = DCOS (DE(IK))
         SD = DSIN (DE(IK))

C  TREATING THE CASE OF AN AZIMUTH-ELEVATION MOUNTING
C----------------------------------------------------------

         IF (AXTYP.EQ.'AZEL')  THEN

               AXKT0 = - (OFFS * SZ) / C

C  TREATING THE CASE OF AN EQUATORIAL ANTENNA
C----------------------------------------------------------

         ELSE IF (AXTYP.EQ.'EQUA') THEN

C   COMPUTE ANGLE BETWEEN NORTH DIRECTION AND SOURCE DIRECTION
C  CORRECTED BY REFRACTION

            CSX = SZ*CAZ(IK)*DCOS(CPHI(IK)) + CZ*DSIN(CPHI(IK))
            snix  = dacos(csx)
            SNX = DSIN(SNIX)

            A0 = -OFFS * CD / C
            AXKT0 = -OFFS * SNX / C


C  TREATING THE CASE OF AN X-Y ANTENNA WITH FIXED AXIS NORTH-SOUTH
C-----------------------------------------------------------------

         ELSE IF (AXTYP.EQ.'X-Y1'.OR.AXTYP.EQ.'X-YN') THEN

            if ((1.d0-(sz*caz(ik))**2).lt.0.d0) print *, 'XY1'

            AXKT0 = -(OFFS * DSQRT(1.D0 - (SZ*CAZ(IK))**2)) / C

C  TREATING THE CASE OF AN X-Y ANTENNA WITH FIXED AXIS EAST-WEST
C---------------------------------------------------------------

         ELSE IF (AXTYP.EQ.'X-Y2'.OR.AXTYP.EQ.'X-YE') THEN

            if ((1.d0-(sz*saz(ik))**2).lt.0.d0) print *, 'XY2'

            AXKT0 = -(OFFS * DSQRT(1.D0 - (SZ*SAZ(IK))**2)) / C

         ELSE

            WRITE (*,321) AXTYP,STAT



321   FORMAT (//' ***** WARNING *****  AXIS MODEL ',A4,' FOR STATION',
     *1X,A8,' NOT AVAILABLE'//' HIT "Q" TO QUIT, "C" TO CONTINUE'//)
            READ (*,'(A1)') STT
            IF (STT.NE.'C'.AND.STT.NE.'c') STOP
         ENDIF

C AXIS MODEL FOR RICHMOND
C-----------------------------

         IF (STAT.EQ.'RICHMOND')  THEN

            AXKT0  = - (OFFS * DSQRT( 1.D0 - (CZ*DSIN(.6817256D0)+
     *      SZ*DCOS(.6817256D0)*(CAZ(IK)*DCOS(.0020944D0)-SAZ(IK)*
     *      DSIN(.0020944D0)))**2))/C

         ENDIF

         IF (IK.EQ.1) AXKTP1 = AXKT0       !   for delay rate
         IF (IK.EQ.2) AXKTM1 = AXKT0       !   for delay rate

      end do


      AXKT(1) = AXKT0         !      correction to the delay
    
      DAXKDT = (AXKTP1 - AXKTM1) / 2.D0
      AXKT(2)   = DAXKDT  * 1.D9        !  correction to the delay rate

      RETURN
      END

C************************************************************************

      SUBROUTINE POLROT (XP,YP,E,ED,JMAP)
C
C   THIS ROUTINE COMPUTES THE POLAR MOTION MATRIX TO CORRECT STATION
C  COORDINATES AND ITS DERIVATIVE RESPECT TO XP AND YP
C
C   THIS ROTATION IS COMPOSED OF A ROTATION WITH RESPECT AXIS 1 OF
C  THE ANGLE (+ Y WOBBLE) AND A ROTATION WITH RESPECT TO AXIS 2 OF
C  THE ANGLE (+ X WOBBLE)
C
C   REMEMBER THAT BIH SYSTEM XW,YW IS A LEFT HANDED SYSTEM AND THE
C  CARTESIAN COORDINATE SYSTEM IS A DIRECT ONE
C
C   LATEST REVISION 1990 NOVEMBRE 19  (N.Z.L.)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION E(3,3),ED(3,3,2)

      COX=DCOS(XP)
      SIX=DSIN(XP)
      COY=DCOS(YP)
      SIY=DSIN(YP)

C**********************************************************
C    CODE FOR ROTATION  R2(X) * R1(Y)
C**********************************************************

c      E(1,1) = COX
c      E(1,2) = SIX*SIY
c      E(1,3) = -SIX*COY
c      E(2,1) = 0.d0
c      E(2,2) = COY
c      E(2,3) = SIY
c      E(3,1) = SIX
c      E(3,2) = -COX*SIY
c      E(3,3) = COX*COY

C**********************************************************
C    CODE FOR ROTATION    R1(Y) * R2(X)
C**********************************************************
       E(1,1) = COX
       E(1,2) = 0.D0
       E(1,3) = -SIX
       E(2,1) = SIX*SIY
       E(2,2) = COY
       E(2,3) = SIY*COX
       E(3,1) = SIX*COY
       E(3,2) = -SIY
       E(3,3) = COX*COY

      IF (JMAP.EQ.3) THEN

C**********************************************************
C    CODE FOR PARTIAL RESPECT TO WX
C**********************************************************
         ED(1,1,1) = -SIX
         ED(1,2,1) = 0.d0
         ED(1,3,1) = -COX
         ED(2,1,1) = COX*SIY
         ED(2,2,1) = 0.d0
         ED(2,3,1) = -SIY*SIX
         ED(3,1,1) = COY*COX
         ED(3,2,1) = 0.d0
         ED(3,3,1) = -SIX*COY
C**********************************************************
C    CODE FOR PARTIAL RESPECT TO WY
C**********************************************************
         ED(1,1,2) = 0.D0
         ED(1,2,2) = 0.d0
         ED(1,3,2) = 0.d0
         ED(2,1,2) = SIX*COY
         ED(2,2,2) = -SIY
         ED(2,3,2) = COX*COY
         ED(3,1,2) = -SIY*SIX
         ED(3,2,2) = -COY
         ED(3,3,2) = -COX*SIY

      ENDIF
C**********************************************************

      RETURN
      END


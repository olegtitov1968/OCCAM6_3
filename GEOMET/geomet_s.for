      PROGRAM GEOMET     
 
C   PROGRAM GEOMET   
C  
C   THIS PROGRAM COMPUTES THE THEORETICAL DELAYS AND DELAY RATES FROM
C   - IERS-92 CONSENSUS MODEL - MODEL 1
C   - IAU-97 MODEL - MODEL 2
C
C     REVISION 1990 NOVEMBER 19  (N.Z.L.)
C     REVISION 1993 MARCH 6 (N.Z.L.) ADD CONSENSUS MODEL IERS-92
C     REVISION 1997 APRIL, 23 (O.T.) ADD CONSENSUS MODEL IERS-96
C     REVISION 2000 DECEMBER, 09 (OT)  Niell mapping function is used
C     REVISION 2001 JANUARY, 31 (OT)  COMMON /ZENITH/ removed
C                                       ZWET removed from STATIM
C     REVISION 2001 APRIL, 18 (OT)  New model by IAU-1997 has been added
C     REVISION 2001 AUGUST, 21 (OT)  On default model is A instead of B
C     REVISION 2001 September, 24 (OT)  COMMON-BLOCK /INF/ instead of /GEOC/
C     REVISION 2002 May, 21 (OT)  'Nut' is read in SORTIM
C     REVISION 2002 May, 21 (OT)  Cutting 'AMBIG' from BASTIM
C     REVISION 2004 May, 20 (OT)  Completion of the version 6.0
C     REVISION 2006 March, 24 (OT)  Completion of the version 6.2
C                IERS-96 model is off
C     REVISION  2006 JUNE, 27  (OT)  Partial derivatives for GAMMA added
C                     Check the reading of BASTIM  !!
C     REVISION  2006 SEPTEMBER, 19  (OT)  Earth gravity propagation effect
C                                           is used on default
C     REVISION 2006 October, 03 (O.T.)  Jupiter is added
C     REVISION 2006 October, 06 (O.T.)  Other planets were added;
C                                        Update of call QUASAR1
C     REVISION 2007 January, 07 (O.T.)  Gravitaional delay correction updated
C     REVISION 2007 July, 09 (O.T.)  Partials for acceleration vector added
C
C
C     LAST REVISION 2021, May 01 (O.T.)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IMPLICIT DOUBLE PRECISION  (A-H, O-Z)
      INCLUDE 'OCCAM_N.FI'
      DOUBLE PRECISION MJD
      CHARACTER quasar*8, CHR*1

c+  OT  19.09.2006
c
c      logical NOEAR
c
c-  OT  19.09.2006

      LOGICAL PPN,RATE,VARK, SCALEDX
      DIMENSION RX(2,3),RY(2,3),RZ(2,3),PHI(3),GLONG(3),ELH(3,3),
     1          rxy(2,3),ISTCAT(2)
      dimension hsta(2),ha(2),sha(2),cha(2),
     +RQU(3),RQP1(3),RQM1(3),HAR(2,3),HSTAR(2,3),Z(2,3)
      DIMENSION ISTA(2),TEMP(2),PRES(2),HUMD(2),IFAC(2)
      DIMENSION RMOON(3,3),DCORG(3,3), XJ(3), YJ(3)
      DIMENSION ZDRY(2), cmniel_d(2,2), cmniel_w(2)
      DIMENSION xe(3),ve(3),ae(3),xs(3),dk(3)
      DIMENSION DV2V1(3), acc_der(3), rb1(3), rb2(3), rb(3), bb12(3)
      DIMENSION rrdm(3), rrds(3), rrdmer(3), rrdven(3),
     *  rrdmar(3), rrdjup(3), rrdsat(3), rrdura(3), rrdnep(3)
	dimension rqu_gs(3)

      DIMENSION rrds_m(3), DVG_M(3)                           ! Solar positions and velocities for the middle of session 30-May-2014

      dimension dvg(3)  !  31-MAy-2013
 
c+ OT 28.03.2006
c      DIMENSION AXKT(2),PCOORD(5),PCOORR(5),PNUTAT(2),PPOLAR(3,3)
c      DIMENSION dimfh(2),dimfw(2),gimfh(2),dqmfh(2),dqmfw(2)
c      DIMENSION grn(2),gre(2),drmfh(2),drmfw(2),vmf1h(2),vmf1w(2)
c      DIMENSION agrad(2),gmfh(2),gmfw(2),grnh(2),greh(2),grnw(2),grew(2)
c      DIMENSION iavailr(nstations),iavailq(nstations),
c     *             iavailv(nstations),  iavail(nstations)
c- OT 28.03.2006


      COMMON  /H/     OMEGA, OMGA
      COMMON  /BAR/   RBAR(3,3), VBAR(3,3), ABAR(3,3)
      COMMON  /G/   UTTJ, XP(9), YP(9), MJD(9), UT1KOR(9), NPAR, IP, JQA
      COMMON  /STA/ R1(3), R2(3), B(3), V1(3), V2(3), A2(3)
      COMMON /PHYS/ C, FL, AH, AU, J20
      COMMON /MATH/ PI, TWOPI, ROG, ROS, ROZ, TOL
      COMMON /INF/ rrdm, rrdmer, rrdven, rrdmar, rrdjup,
     * rrdsat, rrdura, rrdnep, ZDRY, cmniel_d                         !  4-Oct-2006
      COMMON /GRAV/ GMS, GMM, GME, GM_Mer, GM_Ven, GM_Mar, GM_Jup,
     * GM_Sat, GM_Ura, GM_Nep

c
c   Just to read EPHEM.DAT
c


      COMMON /GEOC/ XG(30),YG(30),ZG(30),VXG(30),VYG(30),VZG(30),     !  31-May-2013  
     * BTIM(30),IP2

C------------------------
C    GET PARAMETERS
C------------------------

c      open (22,file='geomet_angle.dat',status='unknown')



      CALL PARAM

C   SELECT MODEL TO USE

      WRITE (*,418)
418   FORMAT (' SELECT THE RELATIVISTIC MODEL TO BE APPLIED '
     */ '       A.- IERS - 92 CONSENSUS MODEL'
     */ '       B.- IAU - 97 MODEL '/
     *' Enter your selection ((A)/B)       (Q to quit)'/)
      READ (*,419) CHR
419   FORMAT (A1)

      MODEL = 1
      IF (CHR.EQ.'Q'.OR.CHR.EQ.'q') STOP
      IF (CHR.EQ.'B'.OR.CHR.EQ.'b') MODEL = 2

c+  OT  19.09.2006
c
c
c      IF (MODEL.EQ.1) THEN
c
c        WRITE (*,618)
c618   FORMAT (' APPLY EARTH GRAVITY EFFECT ON LIGTH PROPAGATION ? ',
c     *' ((Y)/N)  (Q to quit)')
c        READ (*,419) CHR
c        NOEAR=.FALSE.
c
c        IF (CHR.EQ.'Q'.OR.CHR.EQ.'q') STOP
c        IF (CHR.EQ.'N'.OR.CHR.EQ.'n') NOEAR=.TRUE.
c
c      END IF
c
c      WRITE (*,420)
c420   FORMAT (//' *****  PROGRAM RUNNING  *****  Please Wait'/)

c-  OT  19.09.2006





      CALL UNFOLD (17,'BASTIM')
      CALL UNFOLD (18,'STATIM')
      CALL UNFOLD (19,'SORTIM')
      CALL UNFOLD (20,'STACAT')

C   READ ALL THE OBSERVATIONS IN ALL THE BASELINES


      RSUN = 0.d0
      TS_SUN2 = 0.d0
      DTS = 0.d0

      read (17,rec=1) nrec

      do ibasli=1,nrec

c         IF (MOD(IBASLI,100).EQ.0) WRITE (*,429) IBASLI,NREC
c429      FORMAT ('+ ***** PROCESSING RECORD NO. ',I4,' OF ',I4)

         READ (17,rec=ibasli) jBAS,ISOUR,ISTA(1),ISTA(2),D1,SD1,R01,SR1,
     *   DI1,SDI1,RI1,SRI1



C   READ THE SOURCE INFORMATION OF AN OBSERVATION - Change on 3.11.2005

         READ (19,rec=iSOUR) jSOR,quasar,ra2000,de2000,Idatj,UTJ,d00K,
     +                   nut,RAB, DEB, RQP1,
     +                       RAA, DEA, RQM1,
     +                       RAAP, DEAP, RQU, RBAR, VBAR,   !  call source position here!
     +   HSGP1,HP1,HSGM1,HM1,HSG,H,OMEGA,XJ,YJ,
     +   RMOON,DCORG

         TS = dble(IDATJ) + UTJ/twopi + 2400000.5d0
        
         if (ibasli.eq.1) Tstart = TS 
 
c
c     Direct call of dele405
c

 

c         call DELE405 (ts, rrds, rrdm, rrdmer, rrdven,
c     *    rrdmar, rrdjup, rrdsat, rrdura, rrdnep )

c
c      Call of the ephemeride files
c

          CALL EPHEM

          ts = ts -2400000.5d0 

          ts_m = TSTART - 2400000.5d0 +  0.5d0                 !    Moment for the middle of session 
 

          

C**************************************************************
C     COMPUTE THE GEOCENTRIC COORDINATES AND VELOCITIES OF THE SUN
C**************************************************************


         c_grav = 1.0d0*c
         epsilon = c/c_grav
         
         CALL SUNINT (ts ,rrds,DVG)         !  Interpolation for Jupiter  was added in 16 April, 2014


cccc       Lines for retardation
        
c         RESUN = DSQRT(RRDS(1)**2 + RRDS(2)**2 + RRDS(3)**2)
c         DTS = RESUN/C_grav
c         TS_SUN2 = TS - DTS/86400.d0
ccc         if (ibasli.le.10)  print *, RESUN, DTS, TS_SUN2 - TS
c         CALL SUNINT (ts_SUN2,rrds,DVG)         !  Interpolation for velocities
                                               !  was added in May, 2013




         CALL SUNINT (ts_m,rrds_m,DVG_m)         !  Interpolation for teh Sun positions to the middle of session    30 May, 2014
                                             

c         print *, dcorg(1,3), dcorg(2,3), dcorg(3,3)

c         rrds(1) = dcorg(1,3)
c         rrds(2) = dcorg(2,3)
c         rrds(3) = dcorg(3,3)

!          rrds is right for reduction, dcorg is wrong!

        


C-----------------------------------------------------------------
C   COMPUTE THE GEOCENTRIC COORDINATES OF THE MOON IF
C   PEP EPHEMERIS ARE IN USE
C------------------------------------------------------------------

         CALL MOONINT (ts,rrdm)



C   READ THE STATION COORDINATES FOR EACH TIME

         DO  K=1,2

            READ (18,REC=ISTA(K)) JSTA,ISTCAT(K),ISOUR2,TEMP(K),PRES(K),
     *      HUMD(K),IFAC(K),CAB,RX(K,1),RX(K,2),RX(K,3),RY(K,1),RY(K,2),
     *      RY(K,3),RZ(K,1),RZ(K,2),RZ(K,3),CRX,CRY,CRZ,GLONG,PHI,
     *      ELH(K,1),ELH(K,2),ELH(K,3),
     *      Z(K,1),HAR(K,1),HSTAR(K,1),Z(K,2),HAR(K,2),HSTAR(K,2),
     *      Z(K,3),HAR(K,3),HSTAR(K,3),ZDRY(K), cmniel_w,
     *      cmniel_d(k,1), cmniel_d(k,2)

C           RXY DEBE CALCULARSE EN STATIONS

       if (k.eq.1) CALL TRANSF (CRX,CRY,CRZ,PHI1,GLONG1,ELHGT1,R11,0)
       if (k.eq.2) CALL TRANSF (CRX,CRY,CRZ,PHI2,GLONG2,ELHGT2,R12,0)


            DO IJ0=1,3
               RXY(K,IJ0)=DSQRT(RX(K,IJ0)**2+RY(K,IJ0)**2)
            end do

c            print *, ' isour = ', isour, '   isour2 = ', isour2

            IF (ISOUR.NE.ISOUR2) THEN
               WRITE (*,933)
933            FORMAT (' *******  ERROR IN GEOMET, '/
     *' THE TIME CODE IS NOT THE SAME IN THE THREE DIRECT FILES USED'/
     *' HIT ENTER TO CONTINUE, CTRL+C TO QUIT'/)
               READ (*,'(I1)') IDUMMY
            ENDIF

         END DO



C
C  NUMERICAL CALCULATIONS
C
         IF (MODEL.EQ.1)  THEN

            DO JMAP=1,3

               DO I00 = 1,2
                  HSTA(I00) = HSTAR(I00,JMAP)
                  HA(I00)   = HAR(I00,JMAP)
                  SHA(I00) = DSIN(HA(I00))
                  CHA(I00) = DCOS(HA(I00))
               END DO

               HA1 = HSTAR(1,3)
               HA2 = HSTAR(2,3)

C-----------------------------------------------------------------------
C    COMPUTE GEOCENTRIC STATION AND VELOCITY VECTORS IN S4
C-----------------------------------------------------------------------

c
c + OT 04.06.2010
c

               CALL STAVEC (RX, RY, RZ, HA1, HA2, HSG, JMAP)


               DO i=1,3

                  XE(i)=Rbar(i,3)
                  VE(i)=Vbar(i,3)

               END DO

c
c - OT 04.06.2010
c

C-----------------------------------------------------------------------
C    COMPUTE MODEL DELAY
C-----------------------------------------------------------------------

c+  OT  19.09.2006
c
 
c   Scalar_t added on 17.12.2007
c

      IF (JMAP.EQ.1)   CALL VLBF (TS,QUASAR,RA2000,TP1,RQP1,JMAP,
     * TGravp1,RRDS,RRDS_M,
     * pgamm1, teta, cos_A, bb, Phi_a, AA, part_alpha)

      IF (JMAP.EQ.2)   CALL VLBF (TS,QUASAR,RA2000,TM1,RQM1,JMAP,
     * TGravm1,RRDS,RRDS_M,
     * pgamm2, teta, cos_A, bb, Phi_a, AA, part_alpha)

      IF (JMAP.EQ.3)   CALL VLBF (TS,QUASAR,RA2000,T0, RQU, JMAP, 
     * TGrav0,RRDS,RRDS_M,
     * pgamm, teta, cos_A, bb, Phi_a, AA, part_alpha)    !            TETA added on 22 May, 2014; cos_A on 28 Dec 2014, AA on 16, Jan, 2015, part_alpha on 3, Feb, 2015

     
c       IF (JMAP.EQ.1) print *, jmap, rqp1(1), rqp1(2), rqp1(3)
c       IF (JMAP.EQ.2) print *, jmap, rqm1(1), rqm1(2), rqm1(3)
c       IF (JMAP.EQ.3) print *, jmap, rqu(1), rqu(2), rqu(3)
c
c-  OT  19.09.2006

            END DO


         ELSE   ! model.eq.2        (IAU-1997)

            PPN   = .true.
            RATE  = .true.
            VARK  = .true.
            SCALEDX  = .true.        !  TT system is used

C
C Compute geocentric coordinates, velocities and accelerations
C   of the station (according to STAVEC subroutine from original
C   OCCAM distribution). The velocities and accelerations of the
C   stations must account for polar coordinates, tides, etc.
C

            HA1 = HSTAR(1,3)
            HA2 = HSTAR(2,3)

C-----------------------------------------------------------------------
C    COMPUTE GEOCENTRIC STATION AND VELOCITY VECTORS
C-----------------------------------------------------------------------




            CALL STAVEC (RX, RY, RZ, HA1, HA2, HSG, JMAP3)


c
c - OT 04.06.2010
c

C
C    Ephemerides of the Earth and the Sun.
C      It is not clear whether the way to compute the Earth's
C      acceleration at the moment of observation (BARINT subroutine)
C      implemented in OCCAM provides the accuracy needed to compute
C      delay rate observable analytically. In our opinion it is better
C      to use directly numerical differentiation of appropriate
C      Earth's velocity values.
C

            DO i=1,3

               XE(i)=Rbar(i,3)
               VE(i)=Vbar(i,3)
               AE(i)=(Vbar(i,1)-Vbar(i,2))*0.5d0

               XS(i)=XE(i)+DCORG(i,3)

            END DO

C
C    Numerical derivative of the direction to the source
C      (to account for the effect of precession and nutation
C       on the delay rate).
C

            DO i=1,3
               dk(i)=(RQP1(i)-RQM1(i))*0.5d0
            END DO

C
C    Calling QUASAR subroutine
C

            call quasar1 (t0,tgrav0,dtau,dtaugr,
     ,            PPN,RATE,VARK,SCALEDX,           !  05-Oct-2006 (OT)
     ,            RQU,DK,
     ,            R1,V1,R2,V2,A2,
     ,            XE,VE,AE,
     ,            XS,
     ,            GMS,GME,C,
     ,            pgamm)                           !  05-Oct-2006 (OT)



C
C  ADD BASIC PROPAGATION TERM AT STATION 1
C
            DO I=1,3
               DV2V1(I) = V2(I) - V1(I)
            ENDDO
            DATM = ZDRY(1)*cmniel_d(1,1) * DOTPR(DV2V1,rqu)/(C*C)

            t0 = t0 + datm          !   Added on 24.03.2006
C
C    Compute "artificial" values of time delay one second before and
C      one second after the time of observation to substitute
C      analytical value of delay rate into OCCAM
C
            tp1=t0+dtau
            tm1=t0-dtau
            tgravp1=tgrav0+dtaugr
            tgravm1=tgrav0-dtaugr

         END IF
c
c  Estimation of instantaneous acceleration vector    21-Jan-2016
c
    
        call acceler_der (rqu, acc_der)
cc       call acceler_der1 (rqu, acc_der)


c
c
c


         WRITE (17,rec=ibasli) JBAS,ISOUR,ISTA(1),ISTA(2),
     *   D1,SD1,R01,SR1,             !  Observed delay and delay rates with errors
     *   DI1,SDI1,RI1,SRI1,          !  Ionospheric delay and delay rates with errors
     *   TP1,T0,TM1,                 !  Theoretical delay for three moments
     *   TGRAVP1,TGRAV0,TGRAVM1               !  Theoretical delay rate for three moments
     *   ,pgamm, teta, cos_A, bb, Phi_a, AA                         !  28.12.2014, then on 16.01.2015
c     *   ,rrds_m(1), rrds_m(2), rrds_m(3)                           !  03 Feb 2015
     *   ,rrds(1), rrds(2), rrds(3)
     *   ,acc_der(1), acc_der(2), acc_der(3)



      END DO

      stop
      end

C***********************************************************************

                                                        
      SUBROUTINE STAVEC (RX, RY, RZ, HA1, HA2, HSG, JMAP)
C
C-----------------------------------------------------------------------
C GEOCENTRIC STATION VECTORS R1, R2, VELOCITY V2 AND ACCELERATION A2
C AT STATION N2
C
C PROGRAMMERS: J. CAMPBELL, H. SCHUH, FEB. 84
C-----------------------------------------------------------------------
C
C         REVISION  1990 NOVEMBER 19  (N.Z.L.)
C         REVISION  2007 JANUARY  04  (O.T.)
C  LATEST REVISION  2010 JUNE  04  (O.T.)
C
      IMPLICIT DOUBLE PRECISION  (A-H,O-Z)
C
      DIMENSION  RZ(2,3), RX(2,3), RY(2,3)
C
      COMMON  /H/   OMEGA, OMGA
      COMMON  /STA/ R1(3), R2(3), B(3), V1(3), V2(3), A2(3)

      OMEG2 = OMEGA * OMEGA


      SHA1 = DSIN(HA1)
      CHA1 = DCOS(HA1)
      SHA2 = DSIN(HA2)
      CHA2 = DCOS(HA2)

      R1(1) = RX(1,JMAP)
      R1(2) = RY(1,JMAP)
      R1(3) = RZ(1,JMAP)

      R2(1) = RX(2,JMAP)
      R2(2) = RY(2,JMAP)
      R2(3) = RZ(2,JMAP)

      DO I = 1, 3
c
c  The sign was changed to be compatible with the IERS Conventions;
c  It caused some change in the gravitational delay correction.
c


         B(I) = R2(I) - R1(I)           !

      END DO

      V2(1) = -OMEGA * RY(2,JMAP)
      V2(2) = OMEGA * RX(2,JMAP)
      V2(3) =   0.D0


      V1(1) =  -OMEGA * RY(1,JMAP)
      V1(2) =  OMEGA * RY(1,JMAP)
      V1(3) =   0.D0


      A2(1) = - OMEG2 * RX(2,JMAP)
      A2(2) =   OMEG2 * RY(2,JMAP)
      A2(3) =   0.D0

      RETURN
      END


C**********************************************************************

      SUBROUTINE VLBF (TS, SOURCE, RA2000, TAU, S, JMAP, 
     *                  TGRAV, RS, RS_M, p_gamma, 
     *                  Teta, COS_A, BB, Phi, AA, part_alpha)     ! 
 
C
C-----------------------------------------------------------------------
C PROGRAMMERS: J. CAMPBELL, H. SCHUH, FEB.84; N. ZARRAOA
C-----------------------------------------------------------------------
C
C   THIS ROUTINE APPLIES THE MODEL SELECTED,   IERS-92 MODEL
C
c   R1, R2 - VLBI station positions
c   R, V, A - Earth barycentric position, velocity, acceleration
c   RS - Sun position
c   S - quasar position
c
c
C    REVISION  1990 NOVEMBER 19  (N.Z.L.)
C    REVISION  1993 MARCH 4 (N.Z.L.) SIMPLIFICATION OF CODE ADD IERS-92 MODEL
C    REVISION  1993 DEC. 9  (NZL)  TGRAV CORRECTED FOR VEL. OF LIGHT
C    REVISION  1997 APRIL, 23  (OT)  IERS-96 MODEL
C    REVISION  2000 DECEMBER, 09  (OT)  Niell mapping function is used
C    REVISION  2006 MARCH, 24  (OT)  IERS-96 model is off
C    REVISION  2006 JUNE, 27  (OT)  Partial derivatives for GAMMA added
C    REVISION  2006 SEPTEMBER, 19  (OT)  Earth gravity propagation effect
C                                        is used on default
C    REVISION  2006 OCTOBER    04  (OT) Gravitational delay for all planets
C    REVISION 2007 January, 07 (OT)  Gravitaional delay correction updated
C LAST REVISION 2007 December, 17 (OT)  Partials for scalar_t added
C
C************************************************************************

      IMPLICIT DOUBLE PRECISION  (A-H,O-Z)

      CHARACTER*8 SOURCE, SL(26)

      DIMENSION S(3), RE(3), VE(3), AE(3), RS(3), R1MRS(3), R2MRS(3)
      DIMENSION DV2V1(3), AC(3), SA(3)
      DIMENSION RMer(3), RVen(3), RM(3), RMar(3), RJ(3), RSat(3),
     * RUra(3), RNep(3)
      DIMENSION RS_M(3)                                            !  30 May, 2014
      DIMENSION DRDA(3)                                            ! Delay rate, 25 Sep 2015
 

      COMMON  /BAR/ R(3,3), V(3,3), A(3,3)
      COMMON  /STA/ R1(3), R2(3), B(3), V1(3), V2(3), A2(3)
      COMMON /PHYS/ C, FL, AH, AU, J20
      COMMON /MATH/ PI, TWOPI, ROG, ROS, ROZ, TOL
      COMMON /GRAV/ GMS, GMM, GME, GM_Mer, GM_Ven, GM_Mar, GM_Jup,
     * GM_Sat, GM_Ura, GM_Nep
      COMMON /INF/ rm, rmer, rven, rmar, rj, rsat,
     * rura, rnep, ZDRY(2), cmniel_d(2,2)                   !  4-Oct-2006

C
C
      GSGE = GMS/GME

C
      CINV  = 1.D0 / C
      CINV2 = CINV * CINV
      CINV3 = CINV * CINV2
      CINV4 = CINV2 * CINV2

C-----------------------------------------------------------------------
C COMPUTE SCALAR PRODUCTS (S(I) = SOURCE VECTOR)
C-----------------------------------------------------------------------

      DO I=1,3
         RE(I) = R(I,JMAP)
         VE(I) = V(I,JMAP)
         AE(I) = A(I,JMAP)
      ENDDO

      VEXS = DOTPR (VE,S)
      V2XS = DOTPR (V2,S)
      BXVE = DOTPR (B,VE)
      RSUN = DSQRT (DOTPR (RS,RS))      !  Sun distance barycentric

      
      TAU = 0.d0
      TAU0 = 0.d0
      TAUEGR3 = 0.d0
      DCOSR1 = 0.d0
      DCOSR2 = 0.d0
      BETA = 0.d0
      T23 = 0.d0
      T32 = 0.d0

      p_gamma = 0.d0
      
C------------------------------------------------------------
C 1. BASIC DELAY
C-----------------------------------------------------------------------

      TAUG = DOTPR (B,S) * CINV         !   K*b0/c   -> [sec]


c

C 2. EFFECT OF RELATIVISTIC TIME DELAY (TAUR)
C
C    U = G*Msun/R(s-e)
C                                         2
C   [K*b0/c] * [ -(1+gamma)*U/c*c - abs(V) /2c*c - V*w2/c*c ]
C
C-----------------------------------------------------------------------

      TAUR = TAUG * CINV2 *                      ! Sign changed on 04.01.07
     *   (2.d0*GMS / RSUN + 0.5D0 * DOTPR (VE,VE) + DOTPR (V2,VE))


      TAUR2 = TAUG * CINV2 *                                          !                
     *   (0.5D0 * DOTPR (VE,VE) + DOTPR (V2,VE))    !   Coordinate terms is deleted

c
c    Coordinate term removed
c

      



C 3. C**3 - TERM
C-----------------------------------------------------------------------

      TAUC3 = 0.5D0 * CINV3 * VEXS * BXVE        ! Sign changed on 04.01.07


C 4. ABERRATION TERM
C-----------------------------------------------------------------------

      TTIME = BXVE * CINV2


c+ OT 3-Oct-2006

c
c  Geocentric positions of the Sun, Moon and planets are used
c  because barycentic position of the Geocenter is analytically substituted
c

    


      CALL GRAV_DELAY (R1,R2,S,V,B,Jmap,
     *         GMS, RS, RSun, RS_M, TGrav, TETA1 )      !   TETA added on 22 May, 2014, coordinates on 9 June,2014

c
c   TETA1 - for mean epoch here
c 

        CALL GRAV_DELAY2 (R1,R2,S,SOURCE,RA2000,V,B,Jmap,GMS, RS, RSun,
     * COS_A, COS_PSI, BB, phi, TETA, part_alpha, TGrav2, T23, T32,BETA)



c- OT 3-Oct-2006

C 6. CORRECTION DUE TO EARTH'S GRAVITY EFFECT ON
C    LIGHT PROPAGATION (TAUEGR)
C-----------------------------------------------------------------------

C  MODEL FROM IERS STANDARDS 1992

      TAUE = 2.d0 * GME * CINV3

      TAUE1 = DSQRT (DOTPR(R1,R1)) + DOTPR(R1,S)
      TAUE2 = DSQRT (DOTPR(R2,R2)) + DOTPR(R2,S)


      TAUEGR = TAUE * DLOG (TAUE1/TAUE2)

  
c      TAUEGR31 = 1.d0/(DSQRT(DOTPR(R1,R1)) * TAUE1)
c      TAUEGR32 = 1.d0/(DSQRT(DOTPR(R2,R2)) * TAUE2)

c      TAUEGR3 = 1.d12*(TAUEGR32 - TAUEGR31)/C      !   17 May 2016

c      dcosr1 = 180.d0*dacos(dotpr(r1,s)/dsqrt(dotpr(r1,r1)))/pi
c      dcosr2 = 180.d0*dacos(dotpr(r2,s)/dsqrt(dotpr(r2,r2)))/pi



C
C  CONSENSUS MODEL FROM IERS-1992 STANDARDS (TECH. NOTE 13)
C

      TINV = 1.D0 + CINV * (VEXS + V2XS)

c   TAU0 - in seconds!


ccc         TAU0 = (        TAUEGR - TAUG + TAUR2 - TAUC3 - TTIME ) / TINV

c      ELSE


         TAU0 = (TGRAV + TAUEGR - TAUG + TAUR - TAUC3 - TTIME )/TINV
        

c      END IF


c
c  Earth relativistic effect, coordinate term is added to the partial
c


c      p_gamma = 0.5d0 *( TGRAV + 2.d0*TAUG*CINV2*GMS/RSUN )/TINV      !  03 Feb 2015,  New design

      p_gamma = T32

c      p_gamma = 0.5d0 *Tauegr/TINV

c      if(jmap.eq.3) print*, BB, 2.d12*p_gamma, 1.d12*T32
     
       a_teta = (1.d0 - dcos(teta))/dsin(teta)
                        
c      p_gamma1 =  0.5d0*TGRAV/TINV         !  (20.07.2006 in Canberra)  Old design


C
C  ADD BASIC PROPAGATION TERM AT STATION 1
C

      DO I=1,3
         DV2V1(I) = V2(I) - V1(I)  !  meter/sec
      ENDDO

c
c  Geometric delay ( Chapter12, expr.11 )
c
c?
c?  Is it necessary to make selection for the dry mapping function here??
c?	

    
c
c   ZDRY in meters, therefore, one has to divide to CINV2  to get DATM in seconds
c

      DATM = CINV2 * ZDRY(1)*cmniel_d(1,1) * DOTPR(DV2V1,S)
      
      TAU = TAU0 + DATM


      RETURN
      END



      SUBROUTINE GRAV_DELAY (R1,R2,S,V,B,Jmap,GM_Body,RB,RBody,RB_M,                   
     * TGrav_B, TETA )                                    !   TETA added on 22 May, 2014


C  CORRECTION DUE TO RELATIVISTIC LIGHT BENDING (TGRAV) for Sun, Moon and planets
C-----------------------------------------------------------------------
C THIS IS THE ALGORITHM FROM IERS STANDARDS 1992
C
C  THE MODEL SHOULD BE COMPLETED WITH TERMS FOR THE REST OF PLANETS
C AND MOON. ALSO AN EXTRA TERM FOR OBSERVATIONS CLOSE TO THE SUN
C IS RECOMMENDED (SEE IERS-92 STANDARDS)
C
C
C   THIS INCLUDES A CORRECTION FOR THE TIME DELAY BETWEEN THE TIME THE
C  SIGNAL PASSED CLOSE TO THE SUN AND THE TIME IT ARRIVED TO EARTH
C
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      IMPLICIT DOUBLE PRECISION  (A-H,O-Z)

      COMMON /PHYS/ C, FL, AH, AU, J20

      DIMENSION S(3), R1(3), R2(3), V(3,3), RB(3), B(3), RB2(3)
      DIMENSION R1MRBody(3), R2MRBody(3)
      DIMENSION RMer(3), RVen(3), RM(3), RMar(3), RJ(3), RSat(3),
     * RUra(3), RNep(3)

      DIMENSION RB_M(3)

      CINV  = 1.D0 / C
      CINV2 = CINV * CINV
      CINV3 = CINV * CINV2

c       R1, R2 - geocentic positions of the first and second sites
c       RB - geocentric position of the gravity body
c       R1MRBody and R2MRBody - position of the first and second station wrt the gravity body

      RBod1 = 0.d0
      RBod2 = 0.d0 
      TRAV_B1 = 0.d0
      TRAV_B2 = 0.d0
      TGrav_B = 0.d0
      

      DO I=1,3

         
         R1MRBody(I) = R1(I) - RB(I)                                ! 6-Oct-06
c         R2MRBody(I) = R2(I) - (RB(I) + V(I,JMAP)*DOTPR(B,S)/C)    ! 6-Oct-06
         R2MRBody(I) = R2(I) - RB(I)                                ! 26-Jun-15
      


        RB2(i) = - RB_M(I)                            !  for mean epoch!

      ENDDO

      RBod1 = DSQRT (DOTPR(R1MRBody,R1MRBody))   ! Corrected distance from the body to site #1
      RBod2 = DSQRT (DOTPR(R2MRBody,R2MRBody))   ! Corrected distance from the body to site #2
      RBod2_M = DSQRT (DOTPR(RB_M,RB_M)) 

      
      RS1 = DOTPR(R1MRBody,S)         !   -R1 * cos(teta1)
      RS2 = DOTPR(R2MRBody,S)         !   -R2 * cos(teta2)


      TRAV_B1 = RBod1 + RS1           !  R1 - R1*cos(teta1)
      TRAV_B2 = RBod2 + RS2           !  R2 - R2*cos(teta2)

      TGrav_B = 2.d0 * GM_Body * CINV3  * DLOG (TRAV_B1/TRAV_B2)   !  seconds


c
c
c

      COS_TETA = -DOTPR(RB2,S)/RBod2_M                            ! 30 May 2014

      TETA = DACOS(COS_TETA)                                      ! 30 May 2014
  
      RETURN

      END




      SUBROUTINE GRAV_DELAY2 (R1,R2,S,SOURCE,RA2000,V,B,Jmap,GM_Body,RB,
     * RBody,COS_A, COS_PSI, BB, PHI, TETA, part_alpha, TGrav, T23, T34,
     * BETA)                                                        


CCORRECTION DUE TO RELATIVISTIC LIGHT BENDING (TGRAV) for Sun, Moon and planets
C-----------------------------------------------------------------------
C THIS IS THE ALGORITHM developed in 2013
C

C
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      IMPLICIT DOUBLE PRECISION  (A-H,O-Z)

      COMMON /PHYS/ C, FL, AH, AU, J20
      COMMON /MATH/ PI, TWOPI, ROG, ROS, ROZ, TOL

      DIMENSION S(3), R1(3), R2(3), V(3,3), RB(3), B(3), pole(3), B2(3)
      DIMENSION R2MRBody(3), R1MRBody(3)                                          !   25 Aug, 2015
      CHARACTER SOURCE*8
    
      

      CINV  = 1.D0 / C
      CINV2 = CINV * CINV
      CINV3 = CINV * CINV2

c       R1, R2 - geocentic positions of the first and second sites
c       RB - geocentric position of the gravity body
c       

      b1 = 0.d0       

      pole(1)=0.d0
      pole(2)=0.d0
      pole(3)=1.d0

      DA = 0.d0

      B1 = 0.d0
      T23 = 0.d0
      T32 = 0.d0
      T33 = 0.d0
      T30 = 0.d0
      T34 = 0.d0
      R2_Dist = 0.d0
      BETA = 0.d0
      phi = 0.d0

      do i=1,3

        b1 = b1 + (r2(i)-r1(i))**2

c        R2MRBody(I) = R2(I) - (RB(I) + V(I,JMAP)*DOTPR(B,S)/C)    ! 6-Oct-06
        R2MRBody(I) = R2(I) - RB(I)    ! 6-Oct-06
        R1MRBody(I) = R1(I) - RB(I)    ! 6-Oct-06

    
      end do

      b1 = dsqrt(b1)

      BB = DSQRT (DOTPR(B,B))    !  Baseline length
      RBod2 = DSQRT (DOTPR(R2MRBody,R2MRBody))   ! Corrected distance from the body to site #2
      RBod1 = DSQRT (DOTPR(R1MRBody,R1MRBody))   ! Corrected distance from the body to site #1

      COS_TETA = -DOTPR(R2MRBody,S)/RBod2         ! Angle between the body (wrt of the second station) and source (teta)
      COS_TETA_1 = -DOTPR(R1MRBody,S)/RBod1       ! Angle between the body (wrt of the first station) and source (teta)


c      R2E = DSQRT(DOTPR(R2,R2))  !   Geocentric distance to second station



c      COS_TETA = -DOTPR(pole,S)            ! Angle between the body and source (teta)

      RA_bas = Datan(b(2)/b(1))            ! Right ascension of the baseline vector , 16-Jan-2015
      DE_bas = Datan(b(3)/dsqrt(b(1)**2 + b(2)**2) )            ! Right ascension of the baseline vector , 03-Feb-2015

c      write (*,170) b(2), b(1), RA_BAS


      If (b(1).gt.0.d0.and.b(2).gt.0.d0) RA_bas = RA_bas 
      If (b(1).lt.0.d0.and.b(2).gt.0.d0) RA_bas = pi + RA_bas        ! from 90 to 180, 16-Jan-2015
      If (b(1).lt.0.d0.and.b(2).lt.0.d0) RA_bas = pi + RA_bas       ! from 180 to 270, 16-Jan-2015
      If (b(1).gt.0.d0.and.b(2).lt.0.d0) RA_bas = 2.d0*pi + RA_bas       ! from 270 to 360, 16-Jan-2015

     


      DA = RA_bas - RA2000                     ! difference in RA between baseline and source vectors    , 22-Oct-2015


c      RA_Sun = Datan(RB(2)/RB(1))            ! Right ascension of the baseline vector , 22-Oct-2015
c      If (rb(1).gt.0.d0.and.rb(2).gt.0.d0) RA_Sun = RA_Sun 
c      If (rb(1).lt.0.d0.and.rb(2).gt.0.d0) RA_Sun = pi + RA_Sun        ! from 90 to 180, 22-Oct-2015
c      If (rb(1).lt.0.d0.and.rb(2).lt.0.d0) RA_Sun = pi + RA_Sun       ! from 180 to 270, 22-Oct-2015
c      If (rb(1).gt.0.d0.and.rb(2).lt.0.d0) RA_Sun = 2.d0*pi + RA_Sun       ! from 270 to 360, 22-Oct-2015


c      DA_Sun = RA_Bas - RA_2000                 ! difference in RA between baseline and Sun vectors    , 22-Oct-2015


c      if(jmap.eq.3) write (*,170) b(1), b(2), ra_bas, RA_Sun, DA_Sun

c170   format (2(2x,f11.2),3(2x,f7.4))
    

      COS_PHI  = DOTPR(B,S)/BB                 ! Angle between the baseline vector and source (phi)


      COS_PSI  = DOTPR(R2MRBody,B)/(BB*RBod2)        ! Angle between the baseline vector and body (psi) as from station #2

c      COS_PSI  = DOTPR(pole,B)/(BB)        ! Angle between the baseline vector and body (psi)


      TETA = DACOS(COS_TETA)
      PHI = DACOS(COS_PHI)


      R2_Dist = DSQRT(DOTPR(R2,R2))         ! geocentric distance to station #2

      COS_BETA = -DOTPR(R2,S)/R2_Dist        ! Angle beta (geocentric distance and source)
      BETA = 180.d0 * DACOS (COS_BETA)/pi
      
c
c     SINUS of the Angles PHI and TETA
c

      SIN_TETA = DSIN(TETA)
      SIN_PHI = DSIN(PHI)     


      SIN_PSI = DSQRT (1.d0 - COS_PSI**2)

c      if(jmap.eq.3) write (*,*), DSIN(DA_Sun), SIN_PSI, SIN_PHI

c      if(jmap.eq.3) print *, '2', COS_PHI, COS_PHI**2 + SIN_PHI**2


      COS_A = - (COS_PSI + COS_PHI*COS_TETA)/(SIN_TETA * SIN_PHI)

 
      SIN_A =  DSIN(DA) * SIN_PSI /SIN_PHI


      Part_alpha =  BB * SIN_PHI * COS_A/RBod2       !  for alpha, meter
      

      SA = dsqrt(1.d0 - COS_A**2)   ! Sin_A

      CA2 = 1.d0 - 2.d0*COS_A**2
      SA2 = 2.d0* SA * COS_A
      SA4 = 2.d0*SA2 *CA2



      AA = Datan(SIN_A/COS_A)            ! Right ascension of the baseline vector 


c 
c     USE only for radio source structure calculations! Not for GR!   
c 
      If (COS_A.gt.0.d0.and.SIN_A.gt.0.d0) AA = AA 
      If (COS_A.lt.0.d0.and.SIN_A.gt.0.d0) AA = pi + AA        ! from 90 to 180, 16-Jan-2015
      If (COS_A.lt.0.d0.and.SIN_A.lt.0.d0) AA = pi + AA       ! from 180 to 270, 16-Jan-2015
      If (COS_A.gt.0.d0.and.SIN_A.lt.0.d0) AA = 2.d0*pi + AA       ! from 270 to 360, 16-Jan-2015


c
c   Redundant term (main geometry)
c                             


      T0 = -2.d0*GM_Body * cinv3 * bb * cos_phi/RBody 

     
c
c   First term (main deflection)
c


      T1 = 2.d0*GM_Body * cinv3 * bb * sin_phi * sin_teta *COS_A/
     *( RBody *(1.d0 - COS_TETA))

      a1 = 2.d0 * GM_Body * cinv2 * sin_teta/(RBody * (1.d0 - cos_teta))

c
c   Second term (minor deflection)
c


      T2 =  -(GM_Body * cinv3) * (bb * sin_phi * sin_teta *COS_A)**2/
     *( RBody *(1.d0 - COS_TETA))**2

   

c                                   
c   Third term (minor geometry)
c


      T3 = (GM_Body*cinv3)*bb**2 *(1.d0 - (cos_phi * cos_teta)**2)/
     *( RBody**2 *(1.d0 - COS_TETA))


      T32 = -4.d0*GM_Body*cinv3*COS_A**2*(bb*sin_phi/(RBody*teta))**2 

      T33 = 2.d0*GM_Body*cinv3*(BB*Sin_phi/(RBody*teta))**2

c      T30 = 2.d0*GM_Body*cinv3* CA2 *(BB*Sin_phi/(RBody*teta))**2  !  sin(phi)^2*cos(2A)

      GFF = 9.81d0     ! m/sec^2


      T34 =  ((BB*COS_PHI)**2 * 
     * COS_BETA/2.d0 + BB*COS_PHI*R2_Dist)*CINV3  ! sec**3/m

    
c      print *, BB, COS_PHI, COS_BETA, CINV3, T34


c      print *, '  redundant term = ', t0

c      if (180.d0*teta/pi.le.5.d0.and.jmap.eq.3) then

c         print *, 1.d12*(t2+t3 - (t32+t33)), COS_A, SA2

c      end if
                      	
c      TGrav = T0 + T1 + T2 + T3  !  seconds

       TGrav = T0 + T1 + T2 + T3 + T34

c      TGrav = T0 + T1    ! for estimation of sin^2(phi)*cos(2A) - BOTH

c      TGrav = T0 + T1 + T33  ! for estimation of 2.0 * sin^2(phi) * cos**2(A)


      T23 = T2 + T3 


c      TGrav = T0 + T1 + T32  ! for estimation of sin^2(phi)
	

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

c     print *, djdeph, dcorg(1)
c      print *, djdeph, dcorg(2)
c      print *, djdeph, dcorg(3)
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




      SUBROUTINE acceler_der (S,acc_der)

c
C     REVISION 2007 July, 09 (O.T.)  Advanced expression for partials
c
      IMPLICIT DOUBLE PRECISION  (A-H,O-Z)

      DIMENSION C2(3), C3(3), ACC_DER(3), VE(3), S(3)
      COMMON  /BAR/ R(3,3), V(3,3), A(3,3)
      COMMON  /STA/ R1(3), R2(3), B(3), V1(3), V2(3), A2(3)
      COMMON /PHYS/ C, FL, AH, AU, J20

c
c  Calculation of partials - to be used in COLL.FOR
C  Units: [sec]
c

      DO I=1,3
         VE(I) = V(I,3)
         c2(i) = 0.d0
         c3(I) = 0.d0
         acc_der(i) = 0.d0
      ENDDO


      BS = DOTPR (B,S)
      VVS = DOTPR (VE+V2,S)
      BVE = DOTPR (B,VE)

      C2(1) = ( BS * S(1) - B(1)) /c  ! ((s*b)*rqu(1)-b(1))/c   !  sec
      C2(2) = ( BS * S(2) - B(2)) /c  ! ((s*b)*rqu(2)-b(2))/c   !  sec
      C2(3) = ( BS * S(3) - B(3)) /c  ! ((s*b)*rqu(3)-b(3))/c   !  sec

c
c     [(s*b)*(VE(i)+W2(i)) + (s,VE+W2)*b(i) + (VE,b)*s(i) ] / c**2
c

      C3(1) = ( BS * (VE(1)+V2(1)) + VVS*B(1) + BVE*S(1)) /c**2    ! sec
      C3(2) = ( BS * (VE(2)+V2(2)) + VVS*B(2) + BVE*S(2)) /c**2
      C3(3) = ( BS * (VE(3)+V2(3)) + VVS*B(3) + BVE*S(3)) /c**2

c
c   Fourth term to be added !!!! Marked on 22.10.2008
c
c

      ACC_DER(1) = C2(1) + C3(1)
      ACC_DER(2) = C2(2) + C3(2)
      ACC_DER(3) = C2(3) + C3(3)

c      print *, bs/c, b(1)/c, (bs*s(1)-b(1))/c
c      print *, c2(1), c2(2), c2(3)
c      
      


      RETURN

      END





      SUBROUTINE acceler_der1 (S,c2)

c
C     REVISION 2016 June, 17 (O.T.)  Advanced expression for partials
c
      IMPLICIT DOUBLE PRECISION  (A-H,O-Z)

      DIMENSION C2(3), VE(3), S(3)
      COMMON  /BAR/ R(3,3), V(3,3), A(3,3)
      COMMON  /STA/ R1(3), R2(3), B(3), V1(3), V2(3), A2(3)
      COMMON /PHYS/ C, FL, AH, AU, J20

c
c  Calculation of partials - to be used in COLL.FOR
C  Units: [sec]
c

      BS = 0.d0 

      DO I=1,3

         VE(I) = V(I,3)
         c2(i) = 0.d0
         

      ENDDO

      BB = DSQRT(DOTPR(B,B))
      BS = DOTPR (B,S)
      BVE = DOTPR (B,VE)/c

      C2(1) =  - BS * VE(1)/(c*c)   ! (bs)*V(1)/c**2     sec
      C2(2) =  - BS * VE(2)/(c*c)   ! (bs)*V(2)/c**2     sec
      C2(3) =  - BS * VE(3)/(c*c)   ! (bs)*v(3)/c**2     sec


      RETURN

      END


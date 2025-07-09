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
C   LAST REVISION 2021, May 01 (O.T.)
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
      dimension rqu_gs(3), w2(2,3), a0(2,3)
 
      dimension ps(3), sp(3)

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
      COMMON  /STA/ R1(3), R2(3), B(3), V1(3), V2(3), A1(3), A2(3)
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

      open (22,file='geomet_angle.dat',status='unknown')



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

      nn1= 0
      s1 = 0.d0
      steta = 0.d0
      steta1 = 0.d0


      read (17,rec=1) nrec

      do ibasli=1,nrec

c         IF (MOD(IBASLI,100).EQ.0) WRITE (*,429) IBASLI,NREC
c429      FORMAT ('+ ***** PROCESSING RECORD NO. ',I4,' OF ',I4)

         READ (17,rec=ibasli) jBAS,ISOUR,ISTA(1),ISTA(2),D1,SD1,R01,SR1,
     *   DI1,SDI1,RI1,SRI1, time_start, time_length

         
C   READ THE SOURCE INFORMATION OF AN OBSERVATION - Change on 3.11.2005

         READ (19,rec=iSOUR) jSOR,quasar,ra2000,de2000,Idatj,UTJ,d00K,
     +                   nut,RAB, DEB, RQP1,
     +                       RAA, DEA, RQM1,
     +                       RAAP, DEAP, RQU, RBAR, VBAR,   !  call source position here!
     +   HSGP1,HP1,HSGM1,HM1,HSG,H,OMEGA,XJ,YJ,
     +   RMOON,DCORG
                    
         TS = dble(IDATJ) + UTJ/twopi + (32.184d0+d00k)/86400.d0 
     *  + 2400000.5d0

             
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
 
c          print *, ' ok '
          

C**************************************************************
C     COMPUTE THE GEOCENTRIC COORDINATES AND VELOCITIES OF THE SUN
C**************************************************************


         c_grav = 1.d0*c
         epsilon = c/c_grav
         
         CALL SUNINT (ts ,rrds,DVG)        

         RESUN = DSQRT(RRDS(1)**2 + RRDS(2)**2 + RRDS(3)**2)

c         if (ibasli.eq.1) print *, rrds(2)/rrds(1), rrds(3)/
c     * dsqrt(rrds(1)**2 + rrds(2)**2), RESUN


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
     *      Z(K,3),HAR(K,3),HSTAR(K,3),ZDRY(K), 
     *      cmniel_w,
     *      cmniel_d(k,1), cmniel_d(k,2)
     *      , w2(k,1), w2(k,2), w2(k,3)  !   geocentric velocity   08.03.2017 OT
     *      , a0(k,1), a0(k,2), a0(k,3)  !   geocentirc acceleration 13.12.2018 OT

C           RXY DEBE CALCULARSE EN STATIONS

       if (k.eq.1) CALL TRANSF (CRX,CRY,CRZ,PHI1,GLONG1,ELHGT1,R11,0)
       if (k.eq.2) CALL TRANSF (CRX,CRY,CRZ,PHI2,GLONG2,ELHGT2,R22,0)


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


            DO i=1,3

               XE(i)=Rbar(i,3)
               VE(i)=Vbar(i,3)

               AE(i)=(Vbar(i,1)-Vbar(i,2))*0.5d0

            END DO

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
      

             
c
c - OT 04.06.2010
c

C-----------------------------------------------------------------------
C    COMPUTE MODEL DELAY
C-----------------------------------------------------------------------


c        print *, ' OK 1'

        v2(1) = w2(2,1)
        v2(2) = w2(2,2)
        v2(3) = w2(2,3)

        v1(1) = w2(1,1)
        v1(2) = w2(1,2)
        v1(3) = w2(1,3)

        a2(1) = a0(2,1)
        a2(2) = a0(2,2)
        a2(3) = a0(2,3)

        a1(1) = a0(1,1)
        a1(2) = a0(1,2)
        a1(3) = a0(1,3)


        do i=1,3

           dv2v1(i) = v2(i) - v1(i)

        end do

c+  OT  19.09.2006
c
 
c   Scalar_t added on 17.12.2007
c

      IF (JMAP.EQ.1)   CALL VLBF (TS,QUASAR,RA2000,DE2000,TP1,RQP1,JMAP,
     * TGravp1,RRDS,RRDS_M, DTS, pgamm1, teta_m, ae, cos_A, bb, phi_a,
     * R11, R22, AA, part_alpha1,DIFF_A,RSUNp,frp,TGrav2p, Tau_Egr_p, 
     * TG1_P, fp, bwsp, parallax1, Freq_Gravp, Freq_Gravp1, Freq_Gravp2)

      IF (JMAP.EQ.2)   CALL VLBF (TS,QUASAR,RA2000,DE2000,TM1,RQM1,JMAP,
     * TGravm1,RRDS,RRDS_M, DTS, pgamm2, teta_p, ae, cos_A, bb, phi_a,
     * R11, R22, AA, part_alpha2,DIFF_A,RSUNq,frm,TGrav2m, Tau_Egr_m, 
     * TG1_M, fm, bwsm, parallax2, Freq_Gravm, Freq_Gravm1, Freq_Gravm2)


      IF (JMAP.EQ.3)   CALL VLBF (TS,QUASAR,RA2000,DE2000,T0, RQU, JMAP, 
     * TGrav0,RRDS,RRDS_M, DTS, pgamm, teta, ae, cos_A, bb, Phi_a, 
     * R11, R22, AA, part_alpha, DIFF_A, RSUN, freq0, TGrav2, Tau_Egr, 
     * TG1, ff,  bws, parallax, Freq_Grav, Freq_Grav1, Freq_Grav2)


c
c  Frequency is calculated for the Jmap = 3 only
c

      if (JMAP.eq.3) then 

c          nn1 = nn1 + 1

c          s1 = 2.d0*GMS*dsin(teta)*OMEGA/(C**2*RSUN*(1.d0-dcos(teta)))

c          steta = steta + s1

c          steta1 = steta/dble(nn1)

c          print *, s1, nn1, steta1

          T2 = (Tgrav2p  - Tgrav2m)/2.d0 
          T3 = (Tau_Egr_p - Tau_Egr_m)/2.d0  !  Earth gravitational delay

          TG2 = (TG1_P - TG1_M)/2.d0       !  Sun's light deflecton delay

          d_teta = (teta_p - teta_m)/2.d0

c
c   Important terms for general relativity propagation near the Sun
c   and the Earth
c


          freq0 = freq0 + Freq_Grav


c          freq0 = freq0 + TG2 

c          print *, TG2,  Freq_Grav, TG2 - Freq_Grav

          freq0 = freq0 + T3

         
     


      end if
 
c      r1st =  dsqrt(r1(1)**2 + r1(2)**2 + r1(3)**2)
c      r2st =   dsqrt(r2(1)**2 + r2(2)**2 + r2(3)**2)
c      print *, GME/(c**2*r1st),  GME/(c**2*r2st)


       
      if (ibasli.eq.1.and.jmap.eq.3) then


            print *, 'Julian date'   , ts
            print *, 'source vector ', rqu(1), rqu(2), rqu(3)
            print *, 'baseline (meter) ', b(1), b(2), b(3)
            print *, 'baseline length ', bb
            print *, ' omega ', omega

            print *, '1st station (meter) ', r1(1), r1(2), r1(3)
           
c            print *, r1st, GME/(r1st**2)


            print *, '2nd station (meter) ', r2(1), r2(2), r2(3)
        
c            print *, r2st, GME/(r2st**2)

            print *, 'barycentre vel vector (m/s) ', ve(1), ve(2), ve(3)
            print *, 'barycentre velocity (m/s) ', dsqrt(dotpr(Ve,Ve))  

                             
            print *, '2nd station vel vector (m/s)', v2(1), v2(2), v2(3)

            print *, '1st station velocity (m/s)', dsqrt(dotpr(v1,v1))
            print *, '2nd station velocity (m/s)', dsqrt(dotpr(v2,v2))

            print *, '2nd station acceler. vector ', a2(1), a2(2), a2(3)

            print *, 'partial picosecond ', 1.d12 * part_alpha


         
        end if

     
c       IF (JMAP.EQ.1) print *, jmap, rqp1(1), rqp1(2), rqp1(3)
c       IF (JMAP.EQ.2) print *, jmap, rqm1(1), rqm1(2), rqm1(3)
c       IF (JMAP.EQ.3) print *, jmap, rqu(1), rqu(2), rqu(3)
c
c-  OT  19.09.2006

            END DO        !  end of Jmap cycle


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




            CALL STAVEC (RX, RY, RZ, HA1, HA2, HSG, 3)


            BB = dsqrt(dotpr(b,b))

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

c            v2(1) = w2(2,1)
c            v2(2) = w2(2,2)
c            v2(3) = w2(2,3)

            v1(1) = w2(1,1)
            v1(2) = w2(1,2)
            v1(3) = w2(1,3)




            dtau = 0.d0
            dtaugr = 0.d0

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

c            
C    Compute "artificial" values of time delay one second before and
C      one second after the time of observation to substitute
C      analytical value of delay rate into OCCAM
C
            tp1=t0+dtau
            tm1=t0-dtau
            tgravp1=tgrav0+dtaugr
            tgravm1=tgrav0-dtaugr

         END IF


         DATM = ZDRY(1)*cmniel_d(1,1) * DOTPR(DV2V1,rqu)/(C*C)

         t0 = t0 + datm


c
c  Estimation of instantaneous acceleration vector    21-Jan-2016
c
    
        call acceler_der (rqu, acc_der)




        WRITE (17,rec=ibasli) JBAS,ISOUR,ISTA(1),ISTA(2),
     *  D1,SD1,R01,SR1,             !  Observed delay and delay rates with errors
     *  DI1,SDI1,RI1,SRI1,          !  Ionospheric delay and delay rates with errors
     *  TP1,T0,TM1,                 !  Theoretical delay for three moments
     *  TGravp1, TGrav0, TGravm1               !  Geometrical delay
     *  ,pgamm, teta, cos_A, bb         
     *  , Phi_a, AA                         !  28.12.2014, then on 16.01.2015
     *  , freq0, DIFF_A, RSun
     *  , ve(1), ve(2), ve(3)
     *  , parallax                                      !  parallax added on 16.11.2023
     *  , rbar(1,3), rbar(2,3), rbar(3,3)               !  Earth barycentric coordinates  13.05.2025
     
         

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
      COMMON  /STA/ R1(3), R2(3), B(3), V1(3), V2(3), A1(3), A2(3)

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

     

      RETURN
      END


C**********************************************************************

      SUBROUTINE VLBF (TS, SOURCE, RA2000, De2000, TAU, S, JMAP, 
     *    TGRAV, RS, RS_M, DTS, p_gamma, Teta, ae, COS_A, BB, phi,
     * R11, R22, AA, p_alpha, DIFF_A, RSun, freq0, TGrav2, Tauegr, 
     * TG1, FF, bws, parallax, Freq_Grav, Freq_Grav1, Freq_Grav2)
 
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
      DIMENSION DV2V1(3), AC(3), SA(3), DA2A1(3)
      DIMENSION RMer(3), RVen(3), RM(3), RMar(3), RJ(3), RSat(3),
     * RUra(3), RNep(3)
      DIMENSION RS_M(3)                                            !  30 May, 2014
      DIMENSION DRDA(3)                                            ! Delay rate, 25 Sep 2015
      dimension vv2(3)
      
      COMMON  /BAR/ R(3,3), V(3,3), A(3,3)
      COMMON  /STA/ R1(3), R2(3), B(3), V1(3), V2(3), A1(3), A2(3)
      COMMON /PHYS/ C, FL, AH, AU, J20
      COMMON /MATH/ PI, TWOPI, ROG, ROS, ROZ, TOL
      COMMON /GRAV/ GMS, GMM, GME, GM_Mer, GM_Ven, GM_Mar, GM_Jup,
     * GM_Sat, GM_Ura, GM_Nep
      COMMON /INF/ rm, rmer, rven, rmar, rj, rsat,
     * rura, rnep, ZDRY(2), cmniel_d(2,2)                   !  4-Oct-2006

      COMMON /H/ OMEGA, OMGA

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

         VV2(i) = Ve(i) + V2 (i)  !  Second station velocity, VV2 = V + w2
                                                                       

      ENDDO

     

      VEXS = DOTPR (VE,S)
      V2XS = DOTPR (V2,S)
      V1XS = DOTPR (V1,S)

      VV2EXS = DOTPR (VV2,S)

 
      BV2 = DOTPR (B,V2)

     

      BXVE = DOTPR (B,VE)
      BXV2 = DOTPR (B,VV2)

      RSUN = DSQRT (DOTPR (RS,RS))      !  Sun distance geocentric
      VSUN2  = DOTPR (VE,VE)            !  Earth barycentric velocity **2
      VEV2 = DOTPR (VV2, VV2)           !  Combined velocity VV2 **2
      W2 = dsqrt(dotpr(V2,V2))          !  Second station velocity w2

      VW2 = dotpr(VE,V2)                ! (V,w2)


c      BAS1 = DOTPR(B,B)                 !  Baseline length ** 2


      COS1 = VEXS/DSQRT(VSUN2)                      !  COS_TETA for velocity V
      COS_V2 = VV2EXS/DSQRT(VEV2)                   !  COS_TETA for velocity VV2

c      COS2 = BXVE/DSQRT(VSUN2 * BAS1)
c      COS_PH = DOTPR (B,S)/DSQRT(BAS1)

      SIN1 = dsqrt (1.d0 - COS1**2)                  !  SIN_teta for velocity   V 
      SIN_V2 = dsqrt (1.d0 - COS_V2**2)              !  SIN_teta for velocity   VV2

c      SIN_PH = dsqrt (1.d0 - COS_PH**2)
c      COS_D = - (COS2 + COS1 *COS_PHI)/(SIN1 * SIN_PH)

      


      
      TAU = 0.d0
      TAU0 = 0.d0
      TAUEGR3 = 0.d0
      DCOSR1 = 0.d0
      DCOSR2 = 0.d0
      BETA = 0.d0
      T23 = 0.d0
      T32 = 0.d0
      T33 = 0.d0

      TAUC31 = 0.d0


      cos_D1 = 0.d0
      d_ra = 0.d0
      d_de = 0.d0
      teta = 0.d0
      p_gamma = 0.d0
      p_alpah = 0.d0
      SIN_A = 0.d0
      der_gg = 0.d0  
    

      AEXS = DOTPR(AE,S)  
      A2XS = DOTPR(A2,S)  
      ASUN2 = DOTPR(AE,AE)  !  barycentric acceleration **2 

      Freq_b = 0.d0
      Freq0 = 0.d0
      Freq1 = 0.d0
      Freq2 = 0.d0
      Freq_Grav = 0.d0

 
      BV2 = DOTPR (B,V2)

      BAE = DOTPR(B,AE)   


      DO I=1,3

         DV2V1(I) = V2(I) - V1(I)  !  meter/sec
         DA2A1(i) = A2(I) - A1(I)  !  meter/sec**2

      ENDDO


      DVL = dsqrt(DOTPR(DV2V1,DV2V1))  !  Length of vector w2-w1

C------------------------------------------------------------
C 1. BASIC DELAY
C-----------------------------------------------------------------------

      TAUG = DOTPR (B,S) * CINV         !   K*b0/c   -> [sec]

c      if (jmap.eq.3.and.bb.le.3.d6) then 
c
c         write (22,*) 1.d12*(dotpr(r2,s)*v2xs-dotpr(r1,s)*v1xs)/(c*c),
c     *    1.d12*dotpr(r1,s)*v1xs/(c*c),
c     *    1.d12*dotpr(r2,s)*v2xs/(c*c),
c     *    1.d12*(dotpr(r2,s)*v2xs-dotpr(r1,s)*v1xs)/(c*c)    -  
c     *    1.d12*dotpr(r1,s)*(v2xs-v1xs)/(c*c),
c     *    1.d12*taug*v1xs/c,     1.d12*taug*v2xs/c,
c     *    1.d12*taug*(v1xs-v2xs)/c, taug
c  
c          write (22,*) ts, 1.d12*taug*v2xs/c, bb, source


c      end if

       bws = (bxve*v2xs + c*taug*vw2 -2.d0*c*taug*v2xs*vexs)/(c**3)      !  sec

          
c
c  Transformation
c

cc       DDT = (V2XS * DOTPR(R2,S) - V1XS * DOTPR(R1,S))*CINV2

c       print *, R22, R11, DDT



      bb = dsqrt(dotpr(b,b))


c      if (jmap.eq.3.and.(bb.ge.10.350d6.and.bb.le.10.360d6)) then
c  
c          write (22,*) ts, 1.d12 * bws, bb
c
c      end if  

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


      TAUR1 = TAUG * CINV2 *  (2.d0*GMS / RSUN)      !  second order term deleted 08-Feb-2017


      

      TAUR2 = TAUG * CINV2 *                                          !                
     *      (0.5D0 * DOTPR (VE,VE) + DOTPR (V2,VE))    !   Coordinate terms is deleted

      TAUR21 = TAUG *CINV2 * DOTPR(V2,VE)

cc      TAUR2 = TAUR2 - TAUR21 

      TAUR00 = -TAUG * CINV *                      
     *   (-2.d0*GMS *CINV / RSUN - 0.5D0 * DOTPR (VE,VE) *CINV 
     *    - VEXS + VEXS**2*CINV - V2XS )

      TAUR01 = TAUG * DOTPR (V2,VE)*CINV2

c      
c    Coordinate term removed
c

      TAUG_ABER = TAUG * V2XS * CINV 



C 3. C**3 - TERM
C-----------------------------------------------------------------------

      TAUC3 = 0.5D0 * CINV3 * VEXS * BXVE        ! Sign changed on 04.01.07
      



C 4. ABERRATION TERM
C-----------------------------------------------------------------------

      TTIME = BXVE * CINV2

      TTIME0 = - BXVE * CINV2 + BXVE * VEXS/2.d0 * CINV3

      TTIME01 = BXVE * V2XS * CINV3 - 2.d0 * TAUG*VEXS*V2XS*CINV2



cc      print *, taur01, ttime01


c+ OT 3-Oct-2006

c
c  Geocentric positions of the Sun, Moon and planets are used
c  because barycentic position of the Geocenter is analytically substituted
c

   
      CALL GRAV_DELAY (R1,R2,S,SOURCE,V,B,Jmap, DTS,
     *    GMS, RS, RSun, RS_M, TGrav, TETA1)      !   TETA added on 22 May, 2014, coordinates on 9 June,2014


    

c      CALL GRAV_DELAY3 (R1,R2,S,V,B,Jmap,
c     *         GMS, RS, RSun, RS_M, TGrav, TGrav1, T33, TETA1 )      !   TETA added on 22 May, 2014, coordinates on 9 June,2014


c
c   TETA1 - for mean epoch here
c 

        CALL GRAV_DELAY2 (R1,R2,S,SOURCE,RA2000,DE2000,
     *   V,B,Jmap,GMS, RS, RSun, COS_A, COS_PSI, BB, 
     *   phi, TETA, part_alpha, TGrav2, T23, TG1, AA, DIFF_A,
     *   parallax)

 


c- OT 3-Oct-2006

C 6. CORRECTION DUE TO EARTH'S GRAVITY EFFECT ON
C    LIGHT PROPAGATION (TAUEGR)
C-----------------------------------------------------------------------

C  MODEL FROM IERS STANDARDS 1992

      TAUE = 2.d0 * GME * CINV3

      TAUE1 = DSQRT (DOTPR(R1,R1)) + DOTPR(R1,S)
      TAUE2 = DSQRT (DOTPR(R2,R2)) + DOTPR(R2,S)


      TAUEGR = TAUE * DLOG (TAUE1/TAUE2)



C
C  MAIN ABERRATION TERM
C  CONSENSUS MODEL FROM IERS-1992 STANDARDS (TECH. NOTE 13)
C

      


      TINV = 1.D0 + CINV * (VEXS + V2XS )

    
    



c   TAU0 - in seconds!

                                                          
cc         TAU0 = (TAUEGR - TAUG + TAUR2 - TAUC3 - TTIME ) / TINV  !  No Solar general relativity effect

c      ELSE

          

      TAU0 = (TG1 +TAUEGR - TAUG + TAUR - TAUC3 - TTIME)/TINV  


c      print *, tau0, parallax

ccc      tau0 = tau0 + parallax  

c      print *, tau0
       

cc      TAU0 = (TGRAV + TAUEGR)/TINV - TAUG + TAUR00 + TTIME0 

cc     *  + TAUR01 + TTIME01      !  truncation of the second order w2 terms

c      print *, TAU0, TAU00, TAU0 - TAU00

c
c  Earth relativistic effect, coordinate term is added to the partial
c

c
c    Herewith p_gamma is a standard relativistic model
c

       p_gamma = 0.5d0 *(TGRAV + 2.d0*TAUG*CINV2*GMS/RSUN )/TINV      !  03 Feb 2015,  New design

c       p_gamma = TGrav2/2.d0

c
c    T_X - alternative gravitation effect, proportional to sin_teta
c

       p_alpha = -BB *DCOS(PHI) * VEXS * CINV2
        
    




C
C  ADD BASIC PROPAGATION TERM AT STATION 1
C

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



c
c    Frequency
c
 
    
C------------------------------------------------------------
C 1. BASIC FREQUENCY  -(W2-W1,s)/c
C-----------------------------------------------------------------------

      Freq_b = -DOTPR (DV2V1,S) /C         !   -(w2-w1,s)/c   -> [rad]

     
     
C
C 2. FActor for the basic frequency
C       
C    1 - (V,s)/c - (W2,s)/c + ((V,s)/c)**2 - 0.5*(V/c)**2 - (V,w2)/c**2
C                                         
C
C-----------------------------------------------------------------------


      FACT1 = 1.d0 - VEXS/C - V2XS/C + (VEXS/C)**2 - 0.5d0*VSUN2/(C*C)

c      FACT1 = FACT1  - GMS * CINV2/RSUN  !  redundant term


C
C 3. ACCELERATION TERMS  
C-----------------------------------------------------------------------
C

cc      Freq1 = TAUG*AEXS/C - BAE/(C*C)     !  for estimation of acceleration of second station

cc      Freq1 = TAUG * A2XS/C     ! for estimation of full barycentric acceleration
           
      Freq1 = TAUG*AEXS/C + TAUG*A2XS/C - BAE/(C*C)       !  full version


c
C 4. ABERRATION TERMS
C-----------------------------------------------------------------------
c
c   - (w2-w1,V)/c**2 + 0.5d0 * (w2-w1,V)(V,s)/2*c**3 - (w2-w1,s)(V,s)**2/c**3
c


      Freq2 =  -DOTPR (DV2V1,VE)/(C*C) 
     * + 0.5d0 * DOTPR(DV2V1,VE)*VEXS*CINV3
     * - 0.5d0 * DOTPR(DV2V1,S) * VEXS**2 * CINV3

c
c
c


c
c 5. TERMS FOR 0.1 frad
c-----------------------------------------------------------------------
c


      Freq5 = - 2.d0 * DOTPR(DV2V1,S) * VEXS * AEXS * CINV3
     *        + TAUG * DSQRT(VSUN2) * DSQRT(ASUN2)/(C*C)
     *        + 0.5d0 * BAE*VEXS/(C**3) 
     *        + 0.5d0 * BXVE * AEXS / (C**3)


      CALL GRAV_FREQ (R1,R2,S,VE, B, Jmap,GMS, RS, RSun, 
     *   DV2V1, Freq_Grav,  Freq_Grav1, Freq_Grav2)


    
      Freq0 = Freq_b * FACT1 + Freq1 + Freq2 + Freq5

c
c     Removal the major term if the full vector Omega is to be estimated
c

      Freq0 = Freq0 - Freq_b  ! Full OMEGA, if uncommented


      RETURN
      END



      SUBROUTINE GRAV_DELAY (R1,R2,S,SOURCE,V,B,Jmap,DTS,
     * GM_Body,RB,RBody,RB_M,                   
     * TGrav_B, TETA)                     !   TETA added on 22 May, 2014


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
      DIMENSION R1MRBody(3), R2MRBody(3),  R1M(3), R2M(3), VE(3)
      DIMENSION RMer(3), RVen(3), RM(3), RMar(3), RJ(3), RSat(3),
     * RUra(3), RNep(3)

      Character*8 SOURCE

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

      Rp1 = 0.d0
      Rp2 = 0.d0 
      TRAV_B1p = 0.d0
      TRAV_B2p = 0.d0
      TGrav_Bp = 0.d0  

      VS = 0.d0    
      


      DO I=1,3

         
         R1MRBody(I) = R1(I) - RB(I)                                ! 6-Oct-06

c         R2MRBody(I) = R2(I) - (RB(I) + V(I,JMAP)*DOTPR(B,S)/C)    ! 6-Oct-06

         R2MRBody(I) = R2(I) - RB(I)                                ! 26-Jun-15

      
         R1M(I) = R1MRbody(I) - DTS * V(I,JMAP)     ! retarded 
         R2M(I) = R2MRBody(I) - DTS * V(I,JMAP)    ! retarded


        RB2(i) = - RB_M(I)                            !  for mean epoch!

        VE(I) = V(I,JMAP)

 
c        VS = DOTPR(V(I,JMAP),S)                       ! for retardation
c        S(I) = S(I) - ( V(I,JMAP) - S(I)*VS)/C        ! for retardation

      ENDDO

      Rp1 = DSQRT (DOTPR(R1M,R1M))   ! Corrected distance from the body to site #1 retarded
      Rp2 = DSQRT (DOTPR(R2M,R2M))   ! Corrected distance from the body to site #2 retarded
      

      RBod1 = DSQRT (DOTPR(R1MRBody,R1MRBody))   ! Corrected distance from the body to site #1
      RBod2 = DSQRT (DOTPR(R2MRBody,R2MRBody))   ! Corrected distance from the body to site #2
      RBod2_M = DSQRT (DOTPR(RB_M,RB_M)) 


      RS1 = DOTPR(R1MRBody,S)         !   -R1 * cos(teta1)
      RS2 = DOTPR(R2MRBody,S)         !   -R2 * cos(teta2)

      RS1p = DOTPR(R1M,S)         !   -R1 * cos(teta1)        retarded
      RS2p = DOTPR(R2M,S)         !   -R2 * cos(teta2)        retarded


      TRAV_B1 = RBod1 + RS1           !  R1 - R1*cos(teta1)
      TRAV_B2 = RBod2 + RS2           !  R2 - R2*cos(teta2)

      TRAV_B1p = Rp1 + RS1p           !  R1 - R1*cos(teta1)   retarded
      TRAV_B2p = Rp2 + RS2p           !  R2 - R2*cos(teta2)   retarded 




      TGrav_B = 2.d0 * GM_Body * CINV3  * DLOG (TRAV_B1/TRAV_B2)   !  seconds

      TGrav_Bp = 2.d0 * GM_Body * CINV3  * DLOG (TRAV_B1p/TRAV_B2p)   !  seconds   retarded 

     
      BB = dsqrt(dotpr(B,B))


      COS_TETA = -DOTPR(RB2,S)/RBod2_M                            ! 30 May 2014

      TETA = DACOS(COS_TETA)                                      ! 30 May 2014
 
      RETURN

      END




      SUBROUTINE GRAV_DELAY2 (R1,R2,S,SOURCE,RA2000,De2000,
     * V,B,Jmap,GM_Body, RB,
     * RBody,COS_A, COS_PSI, BB, PHI, TETA, part_alpha, TGrav, T23, TG1,
     * AA, DIFF_A, parallax)
cc     * AA_BAS, DIFF_A, parallax)                                                       


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
      DIMENSION R2MRBody(3), R1MRBody(3), R2MRBodyV(3)
      DIMENSION Vector_N(3)                          !   18-Jun-2024
      dimension xx(9)                                !  06-Jul-2025 - for parallax
 
      CHARACTER SOURCE*8
    
      

      CINV  = 1.D0 / C
      CINV2 = CINV * CINV
      CINV3 = CINV * CINV2
      CINV4 = CINV2 * CINV2

c       R1, R2 - geocentic positions of the first and second sites
c       RB - geocentric position of the gravity body
c       

      b1 = 0.d0       

      pole(1)=0.d0
      pole(2)=0.d0
      pole(3)=1.d0

      DA = 0.d0
      DARG1 = 0.d0
      DARG2 = 0.d0

      BB = 0.d0

      B1 = 0.d0
      T23 = 0.d0
      T32 = 0.d0
      T33 = 0.d0
      T30 = 0.d0
      T34 = 0.d0
      R2_Dist = 0.d0
      BETA = 0.d0
      phi = 0.d0
      RBod2 = 0.d0  

      COS_A_BAS = 0.d0  
      COS_A_SUN = 0.d0  
      SIN_A_BAS = 0.d0  
      SIN_A_SUN = 0.d0  
      AA_SUN = 0.d0
      AA_BAS = 0.d0
      DIFF_A = 0.d0
      RA_2 = 0.d0
      DE_2 = 0.d0
      SIN_TETA = 0.d0  
      COS_TETA = 0.d0
      SIN_PHI = 0.d0
      TETA = 0.d0
      

      COS_A =0.d0
      COS_PHI = 0.d0
      COS_PSI = 0.d0
      SIN_PSI = 0.d0
      AA = 0.d0
      RA_bas = 0.d0
      DE_bas = 0.d0
      RA_Sun = 0.d0
      DE_SUn = 0.d0
      DA_Bas = 0.d0
      DA_Sun = 0.d0

      do i=1,3

        b1 = b1 + (r2(i)-r1(i))**2

c        R2MRBodyV(I) = R2(I) - (RB(I) + V(I,JMAP)*DOTPR(B,S)/C)    ! 6-Oct-06
      

        R2MRBody(I) = R2(I) - RB(I)    ! Vector changes its direction
        R1MRBody(I) = R1(I) - RB(I)    ! Vector changes its direction

        Vector_N(I) = 0.d0     

      end do

      b1 = dsqrt(b1)

      BB = DSQRT (DOTPR(B,B))    !  Baseline length

      VV = DSQRT (DOTPR(V,V))    !  Velocity vector length

      RBod2 = DSQRT (DOTPR(R2MRBody,R2MRBody))   ! Corrected distance from the body to site #2
      RBod1 = DSQRT (DOTPR(R1MRBody,R1MRBody))   ! Corrected distance from the body to site #1

      Vector_N(1) = R2MRBody(1)/RBod2
      Vector_N(2) = R2MRBody(2)/RBod2
      Vector_N(3) = R2MRBody(3)/RBod2


      COS_TETA = -DOTPR(R2MRBody,S)/RBod2         ! Angle between the body (wrt of the second station) and source (teta)  
      COS_TETA_1 = -DOTPR(R1MRBody,S)/RBod1       ! Angle between the body (wrt of the first station) and source (teta)



c      R2E = DSQRT(DOTPR(R2,R2))  !   Geocentric distance to second station



cc      COS_TETA = -DOTPR(pole,S)            ! Angle between the body and source (teta)   for source structure

      RA_bas = Datan(b(2)/b(1))                                 ! Right ascension of the baseline vector , 16-Jan-2015
      DE_bas = Datan(b(3)/dsqrt(b(1)**2 + b(2)**2) )            ! Declination of the baseline vector , 03-Feb-2015

c      write (*,170) b(2), b(1), RA_BAS


      If (b(1).gt.0.d0.and.b(2).gt.0.d0) RA_bas = RA_bas 
      If (b(1).lt.0.d0.and.b(2).gt.0.d0) RA_bas = pi + RA_bas        ! from 90 to 180, 16-Jan-2015
      If (b(1).lt.0.d0.and.b(2).lt.0.d0) RA_bas = pi + RA_bas       ! from 180 to 270, 16-Jan-2015
      If (b(1).gt.0.d0.and.b(2).lt.0.d0) RA_bas = 2.d0*pi + RA_bas       ! from 270 to 360, 16-Jan-2015

     
      RA_2 = datan(s(2)/s(1))
      DE_2 = datan(s(3)/dsqrt(s(1)**2 + s(2)**2))


      If (s(1).gt.0.d0.and.s(2).gt.0.d0) RA_2 = RA_2 
      If (s(1).lt.0.d0.and.s(2).gt.0.d0) RA_2 = pi + RA_2        ! from 90 to 180, 16-Jan-2015
      If (s(1).lt.0.d0.and.s(2).lt.0.d0) RA_2 = pi + RA_2       ! from 180 to 270, 16-Jan-2015
      If (s(1).gt.0.d0.and.s(2).lt.0.d0) RA_2 = 2.d0*pi + RA_2       ! from 270 to 360, 16-Jan-2015


      DA_Bas = RA_2 - RA_bas                     ! difference in RA between baseline and source vectors    , 22-Oct-2015 for source structure 


      RA_Sun = Datan(RB(2)/RB(1))                                     ! Right ascension of the Sun , 22-Oct-2015
      DE_Sun = Datan(RB(3)/dsqrt(RB(1)**2 + RB(2)**2))                ! Declination of the Sun , 20-Dec-2016

      If (rb(1).gt.0.d0.and.rb(2).gt.0.d0) RA_Sun = RA_Sun 
      If (rb(1).lt.0.d0.and.rb(2).gt.0.d0) RA_Sun = pi + RA_Sun        ! from 90 to 180, 22-Oct-2015
      If (rb(1).lt.0.d0.and.rb(2).lt.0.d0) RA_Sun = pi + RA_Sun       ! from 180 to 270, 22-Oct-2015
      If (rb(1).gt.0.d0.and.rb(2).lt.0.d0) RA_Sun = 2.d0*pi + RA_Sun       ! from 270 to 360, 22-Oct-2015

 
      DA_Sun = RA_2 - RA_Sun                 ! difference in RA between baseline and Sun vectors    , 22-Oct-2015


c      if(jmap.eq.3) write (*,170) b(1), b(2), ra_bas, RA_Sun, DA_Sun

c170   format (2(2x,f11.2),3(2x,f7.4))
    

      COS_PHI  = DOTPR(B,S)/BB                 ! Angle between the baseline vector and source (phi)

      COS_PHI_V = DOTPR(V,S)/VV        !  Angle between the velocity vector and source vector


      COS_PSI  = DOTPR(R2MRBody,B)/(BB*RBod2)        ! Angle between the baseline vector and body (psi) as from station #2
      
      COS_PSI_V = DOTPR(R2MRBody,V)/(VV*RBod2)      ! Angle between the velocity vector and body

cc      COS_PSI  = DOTPR(pole,B)/(BB)        ! Angle between the baseline vector and body (psi) for source structure


      TETA = DACOS(COS_TETA)               

c      TETA_V = DACOS(COS_TETA_V)


 
      PHI = DACOS(COS_PHI)

c      PHI_V = DACOS(COS_PHI_V)
      

      R2_Dist = DSQRT(DOTPR(R2,R2))         ! geocentric distance to station #2

      COS_BETA = -DOTPR(R2,S)/R2_Dist        ! Angle beta (geocentric distance and source)
      BETA = 180.d0 * DACOS (COS_BETA)/pi
      
c
c     SINUS of the Angles PHI and TETA
c

      SIN_TETA = DSIN(TETA)
      SIN_PHI = DSIN(PHI)     
c      SIN_PHI_V = DSIN(PHI_V)



      SIN_PSI = DSQRT (1.d0 - COS_PSI**2)

      COS_A = - (COS_PSI + COS_PHI*COS_TETA)/(SIN_TETA * SIN_PHI)           

c      COS_A_V = -(COS_PSI_V + COS_PHI_V*COS_TETA)/(SIN_TETA * SIN_PHI_V)
   
 
      SIN_A =  DSIN(DA_Bas) * SIN_PSI /SIN_PHI  !    


      SIN_A_BAS =  DSIN(DA_Bas) * DCOS(DE_BAS) /SIN_PHI  !   angle between direction to North Pole and baseline refered to the radio source 

      SIN_A_SUN =  DSIN(DA_Sun) * DCOS(DE_SUN) /SIN_TETA  !   angle between direction to North Pole and the Sun refered to the radio source

      COS_A_BAS =  (DSIN(DE_BAS) -DSIN(DE_2)*COS_PHI)/
     * (DCOS(DE_2) *SIN_PHI)                            !   angle between direction to North Pole and baseline refered to the radio source 

      COS_A_SUN =  (DSIN(DE_Sun) -DSIN(DE_2)*COS_TETA)/
     * (DCOS(DE_2)*SIN_TETA)                            !   angle between direction to North Pole and the Sun refered to the radio source





c      Part_alpha =  BB * SIN_PHI * COS_A/RBod2       !  for alpha, meter
      

      SA = dsqrt(1.d0 - COS_A**2)   ! Sin_A

      CA2 = 1.d0 - 2.d0*COS_A**2
      SA2 = 2.d0* SA * COS_A
      SA4 = 2.d0*SA2 *CA2



      AA = Datan(SIN_A/COS_A)            ! Right ascension of the baseline vector with respect to the Sun (for radio source structure)
     
     
      AA_SUN = DATAN(SIN_A_SUN/COS_A_SUN)   !   
      AA_BAS = DATAN(SIN_A_BAS/COS_A_BAS)


      If (cos_a_sun.gt.0.d0.and.sin_a_sun.gt.0.d0) AA_Sun = AA_Sun 
      If (cos_a_sun.lt.0.d0.and.sin_a_sun.gt.0.d0) AA_Sun = pi + AA_Sun        ! from 90 to 180, 22-Oct-2015
      If (cos_a_sun.lt.0.d0.and.sin_a_sun.lt.0.d0) AA_Sun = pi + AA_Sun       ! from 180 to 270, 22-Oct-2015
      If (cos_a_sun.gt.0.d0.and.sin_a_sun.lt.0.d0) AA_Sun=2.d0*pi+AA_Sun       ! from 270 to 360, 22-Oct-2015

c      if (jmap.eq.3.and.180.d0 * teta/pi.le.5.d0) then

c         print *, cos_a_bas, sin_a_bas, 180.d0 * AA_Bas/pi

c      end if


      If (cos_a_bas.gt.0.d0.and.sin_a_bas.gt.0.d0) AA_Bas = AA_Bas 
      If (cos_a_bas.lt.0.d0.and.sin_a_bas.gt.0.d0) AA_Bas = pi + AA_Bas        ! from 90 to 180, 22-Oct-2015
      If (cos_a_bas.lt.0.d0.and.sin_a_bas.lt.0.d0) AA_Bas = pi + AA_Bas       ! from 180 to 270, 22-Oct-2015
      If (cos_a_bas.gt.0.d0.and.sin_a_bas.lt.0.d0) AA_Bas=2.d0*pi+AA_Bas       ! from 270 to 360, 22-Oct-2015

c      if (jmap.eq.3.and.180.d0*teta/pi.le.5.d0) print*,180.d0*AA_Bas/pi

      DIFF_A = AA_SUN - AA_BAS               !   angle between direction to baseline and the Sun refered to the radio source

    
c      IF (jmap.eq.3.and. 180.d0 *TETA/PI. le. 5.d0) then

c         print *, DSIN(DIFF_A), SA, 180.d0 * AA_Bas/pi

c      END IF

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
     *( RBody *(1.d0 - COS_TETA))                                            !   change sign on 12-Mar-2024

      B0 = 6.d6     

   


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


c      RS1 = DOTPR(R1MRBody,S)         !   -R1 * cos(teta1)
c      RS2 = DOTPR(R2MRBody,S)         !   -R2 * cos(teta2)


      

c      TRAV_B1 = RBod1**2 - RS1**2    !  R1**2 - (R1*cos(teta1))**2 
c      TRAV_B2 = RBod2**2 - RS2**2    !  R2**2 - (R2*cos(teta2))**2 

c      TR21 = RBod1 - RS1              !  R1 + R1*cos(teta1) 
c      TR22 = RBod2 - RS2              !  R2 + R2*cos(teta2)


      parallax =  BB *cinv *sin_phi * sin_teta * cos_A !  [sec] 


c      parallaxN =(Dotpr(Vector_N,S)* Dotpr(B,S) - Dotpr(Vector_N,B))/c 
     

c
c     Parallax modeling
c

      xx(1) = -10.29   !   [muas]
      xx(2) =  11.97
      xx(3) =   4.13
      xx(4) = -14.32
      xx(5) =   1.18
      xx(6) =   4.85
      xx(7) =   1.10
      xx(8) =  -7.03
      xx(9) =  -4.09


      do i=1,9


        amp = xx(1) + xx(2) * dcos(de_2)*dcos(ra_2) +
     *        xx(3) * dcos(de_2)*dsin(ra_2) +
     *        xx(4) * dsin(de_2)+
     *        xx(5) * dcos(de_2)**2*dcos(2.d0*ra_2) +
     *        xx(6) * dcos(de_2)**2*dsin(2.d0*ra_2) +
     *        xx(7) * dsin(de_2)*dcos(de_2)*dcos(ra_2) +
     *        xx(8) * dsin(de_2)*dcos(de_2)*dsin(ra_2) +
     *        xx(9) * (3.d0*dsin(de_2)**2 - 1.d0)


      end do

      amp = amp/2.06265d11  !  [muasec -> rad]


      parallax_model =  amp * BB *cinv *sin_phi * sin_teta * cos_A !  [sec] 

      TG1 = T0 + T1 + parallax_model

c      print *, t0+t1, parallax_model, tg1
                           	
      TGrav = T0 + T1 + T2 + T3   !  seconds


      T23 = T2 + T3

      
 

      RETURN

      END



      SUBROUTINE GRAV_DELAY3 (R1,R2,S,V,B,Jmap,GM_Body,RB,RBody,RB_M,                   
     * TGrav_B0, TGrav_B1, T3, TETA )                                    !   TETA added on 22 May, 2014


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
      DIMENSION R2MRBody_x(3) 
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
      RBod2_X = 0.d0  
      TRAV_B1 = 0.d0
      TRAV_B2 = 0.d0
      TGrav_B = 0.d0
      TR21 = 0.d0
      TR22 = 0.d0  
      TR3 = 0.d0
      TR4 = 0.d0
      RS1 = 0.d0
      RS2 = 0.d0
      

      DO I=1,3

         
         R1MRBody(I) = R1(I) - RB(I)                                ! 6-Oct-06

         R2MRBody_X(I) = R2(I) - (RB(I) + V(I,JMAP)*DOTPR(B,S)/C)    ! 6-Oct-06

         R2MRBody(I) = R2(I) - RB(I)                                ! 26-Jun-15
      
         RB2(i) = - RB_M(I)                            !  for mean epoch!

      ENDDO

      RBod1 = DSQRT (DOTPR(R1MRBody,R1MRBody))   ! Corrected distance from the body to site #1
      RBod2 = DSQRT (DOTPR(R2MRBody,R2MRBody))   ! Corrected distance from the body to site #2
    
      RBod2_X = DSQRT (DOTPR(R2MRBody_X,R2MRBody_X))   ! Corrected distance from the body to site #2


      BB = DSQRT (DOTPR(B,B))                    ! Baseline lenght
      BS = DOTPR (B,S)                           ! B*cos(phi)  

      
      RS1 = DOTPR(R1MRBody,S)         !   -R1 * cos(teta1)
      RS2 = DOTPR(R2MRBody,S)         !   -R2 * cos(teta2)
      RS2_X = DOTPR(R2MRBody_X,S)     !   -R2-V * cos(teta2)

      CS1 = -RS1/Rbod1                !  cos(teta1)
      CS2 = -RS2/Rbod2                !  cos(teta2)
      CS2_X = -RS2_X/Rbod2                !  cos(teta2)


      RQ = 1.d26
      DC0 = 1.d50

      FF = RQ**2/DC0

      TRAV_B1 = RBod1**2  - RS1**2      !  R1**2 - (R1*cos(teta1))**2 
      TRAV_B2 = RBod2**2  - RS2**2      !  R2**2 - (R2*cos(teta2))**2 

      TR21 = RBod1 - RS1              !  R1 + R1*cos(teta1) 
      TR22 = RBod2 - RS2              !  R2 + R2*cos(teta2)

      TR3 = (Rbod1+RS1)/(Rbod2+RS2)

c     TR4 = (Rbod1+RS1)/(Rbod2_X+RS2_X)


      TR_01 = (TRAV_B1 + FF**2  ) * TR22
      TR_02 = (TRAV_B2 + FF**2  ) * TR21

c      TR_01 = (TRAV_B1 + FF**2  - FF**2 *RS1*RQ) * TR22
c      TR_02 = (TRAV_B2 + FF**2  - FF**2 *RS2*RQ) * TR21


      TR4 = (TRAV_B1 * TR22)/(TRAV_B2 * TR21)

      TR5 = TR_01/TR_02

      TR6 = FF**2*((TRAV_B2-TRAV_B1)/(TRAV_B2 *TR21 *(RBod2 + RS2)))

      TR8 = FF**2*((BB**2 - BS**2)/(TRAV_B2 *TR21 *(RBod2 + RS2)))   !  (b*sin(phi))**2

      TR7 = FF**2 * 2.d0* BB*RBod2/(TRAV_B2 *TR21 *(RBod2 + RS2))


c      print *, TR3, TR5, Tr6

      TGrav_B0 = 2.d0 * GM_Body * CINV3  * dlog(TR3)   !  seconds

      TGrav_B1 = 2.d0 *GM_Body *cinv3 * tr7

c      TGrav_B1 =  2.d0 * GM_Body * CINV3*  DLOG (TR6)   !  seconds
      
      T3 = 2.d0 * GM_Body * CINV3 * TR8
  
c      TGrav_B1 = TT5

c      print *, TGRAV_B0, TGRAV_B1, T3



      COS_TETA = -DOTPR(RB2,S)/RBod2_M                            ! 30 May 2014

      TETA = DACOS(COS_TETA)                                      ! 30 May 2014
  
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
      COMMON  /STA/ R1(3), R2(3), B(3), V1(3), V2(3), A1(3), A2(3)
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
      COMMON  /STA/ R1(3), R2(3), B(3), V1(3), V2(3), A1(3), A2(3)
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


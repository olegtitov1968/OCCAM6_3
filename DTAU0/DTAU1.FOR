      SUBROUTINE GETDAT (IYEAR,NOBS,name_no,i_no)

C************************************************************************
C
C     SUBROUTINE GETDAT
C
C     PERFORMS A READING OF NGS CARD
C
C    IT ALSO CALLS DIFFERENT SUBROUTINES WHICH CREATES FILES WITH
C   INFORMATION ABOUT CORRECTIONS APPLIED TO COORDINATES
C
C   REVISION  1990 NOVEMBER 19   (N.Z.L.)
C   REVISION  1998 JULY 21   (OT)
C   REVISION  1998 September 21   (OT)
C   REVISION  2000 July 12 (OT)    Baselines are not excluded
C   REVISION  2001 March 11 (OT)   COMMON /LOGIC/ deleted
C   REVISION   2001 APRIL 08 (OT)  Reading of the 1st NGS card was removed
C                                 Also CALL of subroutine SHORTB were removed
C   REVISION   2001 DECEMBER 18 (OT)  Remarks of Johannes Boehm
C   REVISION  2002  MAY 16  (OT) - Observations of sources out of
C                                       catalogue are deleted
C   REVISION  2003 JUNE 08 (OT) - Changes in OCCAM_N.FI and
C                                 OPENING.FOR (dictio) and NBASE description
C   LAST REVISION  2006  FEBRUARY 23  (OT) - All changes for OCCAM 6.2 version
C
C*************************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INCLUDE 'OCCAM_N.FI'

      CHARACTER BASFUL(nbases)*18, LINE*96, BIGSTA*6,BBAS*8
      CHARACTER BIGSOR*13
      character*8 FULSOR(nsources), FULSTA(nstations), name_no(400)  ! Change by Volker

      integer*4 NBASE(ntimmax)             !  09.06.2003
      

      COMMON /BSOR1/ BBAS(nobserv)
      COMMON /BASID0/ BASFUL,FULSOR,FULSTA, ista, isor, iobs
      COMMON /BASID1/ NDIM, nbas
      COMMON /STATUS/ ISTAT
      COMMON /LOCCAT/ CORSOR(nsources,2),CORSTA(nstations,3),
     *    VELSTA(nstations,3),EPOCH
      COMMON /BSOR/ BIGSOR(nscans,2),BIGSTA(nstat)

c      print *, ' ok '

      REWIND 4

C
C     OPEN DIRECT ACCESS FILES WHICH ARE THE INPUT FOR INDEPENDENT
C   PROGRAMS TO ADD CORRECTIONS TO THE OBSERVATIONS
C     THE FIRST FILE IS A DICTIONARY OF THE RECORDS IN EACH BASELINE
C     THERE ARE ONE FILE FOR THE SOURCES, ANOTHER FOR STATIONS AND
C   ONE MORE FOR BASELINES
C

      CALL UNFOLD (21,'DICTIO')
      CALL UNFOLD (19,'SORTIM')
      CALL UNFOLD (18,'STATIM')
      CALL UNFOLD (17,'BASTIM')

      OPEN (32,FILE='METEOR')

C
C     INITIALIZE THE DICTIONARY FILE

      DO K0=1,ntimmax
         NBASE(K0)=0
      END DO


      DO K2=1,nbas

         WRITE (21,REC=K2) K2,K2,K2,NBASE
         
      END DO

      print *, ' Program is reading the data file....'

C     INITIALIZE THE VARIABLES THAT KEEP THE STRUCTURE OF THE PREVIOUS
C  FILES

c
c   CHECK OF VARIABLES IN OCCAM_N.FI
c   NOBSERV USUALLY BE MORE THAN NSCANS AND NSTAT
c
      if (nobserv.lt.nscans) then
         print *, ' NOBSERV LESS THEN NSCANS  '
         RETURN
      end if
      if (nobserv.lt.nstat) then
         print *, ' NOBSERV LESS THEN NSTAT  '
         RETURN
      end if
C
C   NOBSERV USUALLY BE MORE THAN NSCANS AND NSTAT
C   THE FIRST 'DO' CYCLE IS A COMMON FOR ALL VARIABLES
C   TO OPTIMIZE COMPUTING TIME

      DO K=1,nobserv
         BBAS(K) = '        '
         IF (K.LE.NSCANS) THEN
            DO K1=1,2
               BIGSOR(K,K1) = '             '
            END DO
         END IF
         IF (K.LE.nstat) BIGSTA(K) = '      '
      END DO

C     READ THE NGS-DATA FILE

      NCT = 0

1     CONTINUE

      READ (4,201,END=999) LINE
201   FORMAT (A96)

      READ (LINE,202) ICOD0,ICOD1
202   FORMAT (74X,I4,I2)

C
C     LOOK FOR CARD  1, SITES, SOURCE AND TIME
C
      IF (ICOD1.EQ.1) THEN
         NCT=NCT+1
c         print *, ' nct = ', nct

cc         IF (MOD(NCT,500).EQ.0) WRITE (*,716) NCT,NOBS
cc716      FORMAT ('+ ***** WORKING ON RECORD  ',I6,' OF ',I6,' *****')

C
C     CALL FILLIN TO READ THE NGS DATA CARDS
C
c         print *, ' ok 00 '
         CALL FILLIN  (ICOD0,LINE,IYEAR,name_no,i_no)

      ENDIF

      GO TO 1

999   CLOSE (17)
      CLOSE (18)
      CLOSE (19)
      CLOSE (21)
      CLOSE (32)
      RETURN
      END

************************************************************************

      SUBROUTINE FILLIN (ICOD0,LINE,IYEAR,name_no,i_no)

C**********************************************************************
C
C     SUBROUTINE FILLIN (ICOD0,LINE,IYEAR)
C
C    THIS SUBROUTINE IS DESIGNED TO READ THE NGS FILE
C
C   91/10/15  CABLE CAL SIGN REVERSED IN ACCORDANCE WITH CALC 7.4  -NZL-
C
C   REVISION  1991 OCTOBER 15   (N.Z.L.)
C   REVISION  1999  SEPTEMBER 21   (OT)
C   REVISION  2001  APRIL 08  (OT)
C   REVISION  2002  MAY 16  (OT) - Observations of sources out of
C                                       catalogue are deleted
C   REVISION  2002  OCTOBER 16  (OT) - 1 sec added; write new sources
C                                           into file
C   REVISION  2003  JANUARY 23  (OT) - 1 sec extracted;
C   REVISION  2003 JUNE 08 (OT) - Changes in OCCAM_N.FI and OPENING.FOR (dictio)
C   REVISION  2004 OCTOBER 28 (OT) - Cable delay correction moved to
C                                         KVLBI.FOR
C   LAST REVISION  2008 NOVEMBER 14 (OT) - 3 sec extracted !!
C
C
C***********************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INCLUDE 'OCCAM_N.FI'

      CHARACTER BASFUL*18,LINE*96,CCH*1
      CHARACTER BASEL*18,SOUR*8,FULSOR*8,FULSTA*8,SOU*8
      CHARACTER BIGSOR*13,BIGSTA*6,BBAS*8,name_no(400)*8   ! Change by Volker

      INTEGER*4 NBASE(ntimmax),K2,KK2,KST1,KST2, i0
      

      COMMON /BSOR1/ BBAS(nobserv)
      COMMON /BASID0/ BASFUL(nbases),FULSOR(nsources),FULSTA(nstations),
     *ista, isor, iobs
      COMMON /BASID1/ NDIM, nbas
      COMMON /STATUS/ ISTAT
      COMMON /LOCCAT/ CORSOR(NSOURCES,2), CORSTA(NSTATIONS,3),
     * VELSTA (NSTATIONS,3),EPOCH
      COMMON /BSOR/ BIGSOR(nscans,2),BIGSTA(nstat)
      COMMON /MATH/ PI, TWOPI, ROG, ROS, ROZ, TOL

C     READ AGAIN THE CARD 1

      READ (LINE,203) BASEL,SOUR,IYEAR,IMON,IDAY,IHOU,IMIN,SEC
203   FORMAT (A18,2X,A8,1X,I4,1X,I2,1X,I2,1X,I2,1X,I2,1X,F14.10)




c
c    Observations of sources out of catalogue are deleted here
c


c      open (10, file = 'test1.dat', status='unknown')

c      write (10,*) sec


      k2= 0
      kst1 = 0
      kst2 = 0
      factor_index = 0.d0

       

      do i=1,i_no         ! The number of sources out of catalogue
         if(sour.eq.name_no(i)) then                 !  16.05.2002
           return                                    !
         end if
      end do

      CALL DDOY (IYEAR,IMON,IDAY,IHOU,IMIN,ITIM,IDOY)

C     INSERT THE DATA CARD IN THE DIRECT ACCESS FILES FOR FUTURE
C  CORRECTIONS


      CALL ACDIR (SOUR,IYEAR,IDOY,IHOU,IMIN,SEC,NREC0)
    
      CALL STADIR (BASEL,SOUR,IYEAR,IDOY,IHOU,IMIN,SEC,NREC1,NREC2)

      CALL BASDIR (BASEL,SOUR,IYEAR,IDOY,IHOU,IMIN,SEC,NREC3)

      IF (ISTAT.EQ.990) RETURN

      

C
C     INITIALIZE THE OBSERVED VALUES TO 0
C

      IFLAG=0
      DELAY=0.
      SIGDEL=0.
      RATE=0.
      SIGRAT=0.
      CORIDE=0.
      CORIRA=0.
      CORSGD=0.

      CORSGR=0.

      delay_cor = 3.d9         !   14.11.2008
      
C
C     LOOK FOR THE REST OF THE CARDS
C

2     CONTINUE

      READ (4,201,END=4) LINE
      READ (LINE,202) ICD0,ICD1
201   FORMAT (A96)
202   FORMAT (74X,I4,I2)

C     CHECK THE END OF AN OBSERVATION CARD

      IF (ICD0.NE.ICOD0) THEN
          IFLAG=1
          GO TO 3
      ENDIF

C     CARD  2, OBSERVED VALUES

      IF (ICD1.EQ.2) THEN
          READ (LINE,204) DELAY,SIGDEL,RATE,SIGRAT,IQC
204       FORMAT(F20.8,F10.5,F20.10,F10.5,I2)
c 204       format (f20.8,f10.6,f20.8,f10.6,i2)

c          print *, delay



            IF (DELAY.LT.-3.8d8.and.DELAY.GT.-4.2d8) then

               DELAY = DELAY + 4.0d8                        !   14.11.2008

            end if

            IF (DELAY.GT.3.8d8.and.DELAY.LT.4.2d8) then

               print *, delay, ' before '
 
               DELAY = DELAY - 4.0d8                        !   14.11.2008

               print *, delay

            end if


            IF (DELAY.LT.-9.d8.and.DELAY.GT.-1.1d9) then

               DELAY = DELAY + 1.d9                        !   14.11.2008

            end if
  
            IF (DELAY.GT.9.d8.and.DELAY.LT.1.1d9) then

               DELAY = DELAY - 1.d9                        !   14.11.2008

            end if
 
             IF (DELAY.LT.-1.9d9.and.DELAY.GT.-2.1d9) then
 
                DELAY = DELAY + 2.0d9                        !   14.11.2008
 
             end if
   
             IF (DELAY.GT.1.9d9.and.DELAY.LT.2.1d9) then
                                              
                DELAY = DELAY - 2.0d9                        !   14.11.2008
 
             end if

             IF (DELAY.LT.-4.9d9.and.DELAY.GT.-5.1d9) then
 
                DELAY = DELAY + 5.0d9                        !   14.11.2008
 
             end if
   
             IF (DELAY.GT.4.9d9.and.DELAY.LT.5.1d9) then
                                              
                DELAY = DELAY - 5.0d9                        !   14.11.2008
 
             end if


             IF (DELAY.LT.-2.9d9.and.DELAY.GT.-3.1d9) then
 
                DELAY = DELAY + 3.0d9                        !   14.11.2008
 
             end if
   
             IF (DELAY.GT.2.9d9.and.DELAY.LT.3.1d9) then
                                              
                DELAY = DELAY - 3.0d9                        !   14.11.2008
 
             end if


c          IF (DELAY.LE.-1.5d9) then
c
c                  DELAY=DELAY+delay_cor        !   14.11.2008
c
c          end if
c          IF (DELAY.GE.1.5d9)  then
c
c                DELAY=DELAY-delay_cor        !   14.11.2008
c
c          end if

C
C     CHECK THE QUALITY CODE OF THE OBSERVATION. IN THE FUTURE WE SHOULD
C  BE ABLE TO DISCERN THE MEANINGS OF THE DIFFERENT CODES TO MAKE A
C  BETTER SELECTION OF THE DATA DESIRED
C

          IF (IQC.NE.0) RETURN

      ENDIF

C      CARD  3, FRINGE AMPLITUDE AND PHASE



      IF (ICD1.EQ.3) CONTINUE

c          IF (ICD1.EQ.3) then


            

c             read (line,233) factor_index
c  233        format (67x,f3.1)

c           write (10,*) IQC, factor_index

c            if (factor_index.ne.0.d0) return

c           write (10,*) IQC

c            if (iqc.ne.0) return

c          end if

C      CARD  4, SYSTEM AND ANTENNA TEMPERATURE

      IF (ICD1.EQ.4) CONTINUE

C
C      CARD  5,
C      CABLE CAL CORRECTION TO BE MOVED TO ADJUSTMENT PART  !!!!!
C   !  28.10.2004 (OT)
C


      IF (ICD1.EQ.5) then
          READ (LINE,205) CAB1,CAB2
205       FORMAT (2F10.5)
      end if

C      CARD  6, METEOROLOGICAL DATA.

C     THIS VALUES MUST BE WRITTEN ON STATIM

      IF (ICD1.EQ.6) THEN
          READ (LINE,230) TDRY1,TDRY2,PRESS1,PRESS2,HUM1,HUM2,KW1,KW2
230       FORMAT (6F10.3,1X,I1,1X,I1)

      ENDIF

C      CARD  7, TIME INFORMATION

      IF (ICD1.EQ.7) THEN


c
c  Time1 - difference between the reference scan epoch and time of observation start (sec)
c  Time2 - lenght of scan (sec)
c

         READ (LINE,241) Time1, Time2

241      format (f10.5,f10.5)
        

      END IF 

C      CARD  8, IONOSPHERIC CORRECTION

c       print *, ' ok '

C     THIS VALUES MUST BE WRITTEN ON BASTIM




      IF (ICD1.EQ.8) THEN
          READ (LINE,206) DELION,SGDION,RATION,SGRION,IQCION
206       FORMAT (F20.10,F10.5,F20.10,F10.5,I3)
cc          print 206, DELION, SGDION, RATION, SGRION, IQCION
          IF (IQCION.NE.0) then
              RETURN
          else
             CORIDE = -1.D0*DELION
             CORIRA = -1.D0*RATION
             CORSGD=SGDION
             CORSGR=SGRION
          end if
      ENDIF

3     CONTINUE

C     IF IFLAG.NE.1 NOT ALL THE CARDS OF AN OBSERVATION HAVE BEEN
C  ALREADY READ, SO WE CONTINUE READING

      IF (IFLAG.NE.1) GO TO 2

4     CONTINUE

C     BACKSPACE IN ORDER TO PLACE THE FILE AT THE BEGINNING OF THE NEXT
C  OBSERVATION

      BACKSPACE 4

C
C     WRITE THE CABLE VALUES AND METEOROLOGICAL DATA ON STATIM
C

      READ (18,REC=NREC1) J001,ICAT1,ISO01,TD
ccc      print *, ' check ', J001, TD

      IF (TD.GT.998.D0) THEN
        WRITE (18,REC=NREC1) J001,ICAT1,ISO01,TDRY1,PRESS1,HUM1,KW1,CAB1
         READ (19,REC=ISO01) I,SOU,CRA,CDEC,IDT,UT0J
         T01 = IDT + UT0J/(2.D0*PI)
         WRITE (32,601) T01,TDRY1,PRESS1,HUM1,KW1,CAB1,ICAT1
      ENDIF

601   FORMAT (F11.5,2X,F6.2,2X,F7.2,2X,F6.2,2X,I1,2x,F8.5,2x,I2.2)

      READ (18,REC=NREC2) J002,ICAT2,ISO02,TD
      IF (TD.GT.998.D0) THEN
      WRITE (18,REC=NREC2) J002,ICAT2,ISO02,TDRY2,PRESS2,HUM2,KW2,CAB2
         READ (19,REC=ISO02) I,SOU,CRA,CDEC,IDT,UT0J
         T01 = IDT + UT0J/(2.D0*PI)
         WRITE (32,601) T01,TDRY2,PRESS2,HUM2,KW2,CAB2,ICAT2
      ENDIF

C
C     WRITE THE DATA VALUES WITH THE CORRECTIONS SELECTED ON BASTIM
C
C     DELAY IS IN MICROSECS AND RATE IS IN PICOSEC/SEC
C   SIGMAS IN NANOSEC AND PICOSEC/SEC



c     if (basel(1:8).eq.'KATH12M '.or.basel(11:18).eq.'KATH12M ')
c     *  then

c      DELAY = DELAY/1000.d0           !   CORCAB moved to adjustment part
c         print 222, t01, coride, delay

c        else

cc      if (sour.eq.'0235+164'.and.dabs(coride).ge.0.1d0) return

      DELAY = (DELAY + CORIDE)/1000.d0           !   CORCAB moved to adjustment part

ccc         delay = delay/1000.d0


c      end if

      RATE = RATE + CORIRA


      SIGDEL = DSQRT (SIGDEL**2 + CORSGD**2)
      SIGRAT = DSQRT (SIGRAT**2 + CORSGR**2)

C     WRITE ON BASTIM THE OBSERVED VALUES CORRECTED BY IONOSPHERE AND
C  CABLE (!!).  WRITE ALSO THE IONOSPHERIC CORRECTIONS TO BE ADDED TO
C  THE MODEL (!!)

      READ (17,REC=NREC3) JBAS,ISOU,IST1,IST2

      WRITE (17,REC=NREC3) JBAS,ISOU,IST1,IST2,DELAY,SIGDEL,RATE,SIGRAT,
     * DELION,SGDION,RATION,SGRION
     *,time1, time2                        !   28-May-2019 (OT) 

C
C      WRITE ON THE DICTIONARY FILE NREC3 ON THE CORRESPONDING RECORD
C
      DO I0=1,nbas
         IF (BASEL.EQ.BASFUL(I0)) then
              K2=I0
c              print *, k2, i0
        end if
      END DO

c      print *, ' nbas = ', nbas, k2


      READ (21,REC=K2) KK2,KST1,KST2,NBASE   ! K2 - order number of base

      IF (K2.EQ.1) KK2=nbas
      KST1=ICAT1
      KST2=ICAT2

           

c
c  NTIMMAX = 3000 if you change if you will have to change KVLBI_IE.FOR !!!
c
      DO I0=1,NTIMMAX   !  max = 3000 if you change if you will have to change KVLBI_IE.FOR
c         print *, ' i0 = ', i0

         IF(NBASE(I0).EQ.0) THEN
            NBASE(I0)=NREC3
            GO TO 12
         ENDIF
      END DO

 

C   IF THE PROGRAM COMES THIS WAY, THE LIMIT OF 500 OBSERVATION/BASELINE
C  HAS BEEN REACHED. SEND A WARNING MESSAGE

      WRITE (*,900) K2, NTIMMAX
900   FORMAT (///' ****** WARNING !!! *******'//
     *'THE BASELINE ',I3,' HAS MORE THAN ',I4,'DATA POINTS. ANY FURTHER'
     */'  OBSERVATION IN THESE BASELINE WILL BE IGNORED BY THE ESTIMATI'
     *'ON'/'  PROGRAMS'//' ENTER Q TO QUIT OR C TO CONTINUE',5X,
     *' (QUIT RECOMMENDED)'//)
      READ (*,400) CCH
400   FORMAT (A1)
      IF (CCH.NE.'C')  ISTAT=990
      IF (ISTAT.EQ.990) RETURN

C   IF NO PROBLEM ARISES WRITE THE UPDATED IDICTIO

      ambig=0.d0                    ! Change by Volker
12    continue

   
     
      WRITE (21,REC=K2) KK2,KST1,KST2,NBASE

c+ jboehm
      READ  (21,REC=1) KK2,KST1,KST2,NBASE,AMBIG
      WRITE (21,REC=1) nbas,KST1,KST2,NBASE,AMBIG
c- jboehm

      RETURN
      END

************************************************************************

      SUBROUTINE ACDIR (SOUR,IYEAR,IDOY,IHOUR,IMIN,SEC,NREC0)

C***************************************************************************
C
C     SUBROUTINE ACDIR(SOUR,IYEAR,IDOY,IHOU,IMIN,SEC,NREC0)
C
C    THIS SUBROUTINE GENERATES A DIRECT ACCESS FILE FOR FURTHER
C  CALCULATION OF PRECCESION AND NUTATION CORRECTIONS FOR THE SOURCES
C
C   REVISION  1990 NOVEMBER 19   (N.Z.L.)
C   REVISION  1998 JULY 21   (OT)
C   REVISION  1998 September 21   (OT)
C   REVISION  2001 MARCH  13   (OT)   Improvement of cycles
C   LAST REVISION  2003 JUNE 08 (OT) - Changes in OCCAM_N.FI and OPENING.FOR (dictio)
C
C***************************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INCLUDE 'OCCAM_N.FI'

      CHARACTER BASFUL*18
      CHARACTER SOUR*8,FULSOR*8,SOR0*8,FULSTA*8,CCH*1
      CHARACTER BIGSOR*13,BIGSTA*6,SOUR13*13,TIM13*13
      COMMON /BASID0/ BASFUL(NBASES),FULSOR(NSOURCES),FULSTA(NSTATIONS),
     * ISTA,ISOR, IOBS
      COMMON /BASID1/ NDIM, nbas
      COMMON /STATUS/ ISTAT
      COMMON /LOCCAT/ CORSOR(NSOURCES,2), CORSTA(NSTATIONS,3),
     *        VELSTA(NSTATIONS,3), EPOCH
      COMMON /BSOR/ BIGSOR(nscans,2),BIGSTA(nstat)

C
C     WRITE THE SOURCE NAME AND THE TIME (IN CHARACTER MODE) INTO TWO
C  VARIABLES TO COMPARE WITH A TABLE OF THE VALUES ALREADY WRITTEN
C
C    THIS TRICK IS DONE TO AVOID UNNECCESSARY REPETITIONS OF THE
C  COMPUTATIONS
C
C    THE TABLE (BIGSOR) IN THE COLUMN 2 WILL BE ALSO USED TO GENERATE
C  SIMILAR TABLES IN THE DIRECT ACCESS FILES FOR SOURCES AND BASELINES
C  IN ORDER TO DECREASE THE NUMBER OF BYTES USED IN THOSE TABLES
C  (OTHERWISE THEIR SIZE WOULD BE OUT OF RANGE)
C

      WRITE (SOUR13,100) SOUR
100   FORMAT (A8,'     ')
      ISEC = SEC + 1.D-7

      kcon = 0


      WRITE (TIM13,101) IYEAR,IDOY,IHOUR,IMIN,ISEC
101   FORMAT (I4.4,I3.3,I2.2,I2.2,I2.2)

C
C     COMPARE THE TWO VARIABLES WITH THE VALUES ALREADY WRITTEN ON
C  THE TABLE.
C     IF THEY ARE A NEW PAIR, WRITE THEM IN BIGSOR COLUMN ONE FOR
C  SOURCE NAME, COLUMN TWO FOR TIME, THEN WRITE  A NEW RECORD ON
C  THE DIRECT ACCESS FILE
C

      DO I=1,iobs

         IF (SOUR13.EQ.BIGSOR(I,1).AND.TIM13.EQ.BIGSOR(I,2)) THEN
            NREC0=I
            RETURN
         ENDIF

         IF (BIGSOR(I,1).EQ.'             ') THEN
            BIGSOR (I,1) = SOUR13
            BIGSOR (I,2) = TIM13
            NREC0=I
            DO J=1,isor
               IF (SOUR.EQ.FULSOR(J)) then
                  KCON=J
                  exit
               END IF
            end do

c            print 77, i, sour, kcon

c 77         format (i5,2x,a8,2x,i5)

            CALL JULDAT (IYEAR,IDOY,IHOUR,IMIN,SEC,IDATJ,UTJ,0)


            IF (I.GT.1) THEN
               READ (19,REC=1) NRC,SOR0,C1,C2,I3,C4,C5
               WRITE (19,REC=1) I,SOR0,C1,C2,I3,C4,C5
            ENDIF

            D0K = 0.D0
            WRITE (19,REC=I) I,SOUR,CORSOR(KCON,1),CORSOR(KCON,2),
     *       IDATJ,UTJ,D0K
            RETURN
         ENDIF

c         print *, ' i = ', i

      END DO
C
C    IF IT REACHES THIS POINT, AN ERROR MESSAGE MUST BE SENT TO SCREEN
C

      WRITE (*,900)
900   FORMAT (///' ****** WARNING !!! *******'//
     *'  THE PROGRAM MUST NOT REACH THIS POINT AT SUBROUTINE "ACDIR"'/
     *'  THERE ARE MORE THAN 5000 SCANS. YOU MUST CHANGE THIS PROGRAM'/
     *'  TO DEAL WITH THIS POSSIBILITY. NO FURTHER DATA WILL BE USABLE'
     *//' ENTER Q TO QUIT OR C TO CONTINUE',5X,' (QUIT RECOMMENDED)'//)
      READ (*,400) CCH
400   FORMAT (A1)
      IF (CCH.NE.'C')  ISTAT=990
      RETURN
      END


************************************************************************


      SUBROUTINE STADIR(BASEL,SOUR,IYEAR,IDOY,IHOU,IMIN,SEC,NREC1,NREC2)

C*************************************************************************
C
C      SUBROUTINE STADIR (BASEL,IYEAR,IDOY,IHOU,IMIN,SEC,NREC1,NREC2)
C
C     THIS SUBROUTINE GENERATES A DIRECT ACCESS FILE CONTAINING
C  INFORMATION INDEXED WITH THE PAIR STATION-TIME, IN ORDER TO
C  ADD CORRECTIONS DEPENDENT ON THE STATIONS WITH NO REPETITIVE
C  COMPUTATIONS
C
C   REVISION  1990 NOVEMBER 19   (N.Z.L.)
C   REVISION  1998 JULY 21   (OT)
C   REVISION  1998 September 21   (OT)
C   REVISION  2001 MARCH  13   (OT)   Improvement of cycles
C   LAST REVISION  2003 JUNE 08 (OT) - Changes in OCCAM_N.FI and OPENING.FOR (dictio)
C
C***************************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INCLUDE 'OCCAM_N.FI'

      CHARACTER BASFUL*18,BASEL*18,STAT(2)*8,FULSOR*8,FULSTA*8
      CHARACTER BIGSTA*6,BIGSOR*13,BBAS*8,STAT4(2)*6,TIM13*13,CCH*1
      CHARACTER*8 SOUR
      DIMENSION NREC(2),I01(2),ISTCAT(2)
      COMMON /BSOR1/ BBAS(nobserv)
      COMMON /BASID0/ BASFUL(NBASES),FULSOR(NSOURCES),FULSTA(NSTATIONS),
     * ISTA, ISOR, IOBS
      COMMON /BASID1/ NDIM, nbas
      COMMON /STATUS/ ISTAT
      COMMON /LOCCAT/ CORSOR(NSOURCES,2),CORSTA(NSTATIONS,3),
     * VELSTA(NSTATIONS,3),EPOCH
      COMMON /BSOR/ BIGSOR(nscans,2),BIGSTA(nstat)

C     READ THE STATIONS IN EACH OBSERVATION

      READ (BASEL,200) STAT

200   FORMAT (A8,2X,A8)

C     CHECK THE STATION NUMBER IN "FULSTA" AND THE TIME NUMBER IN
C  "BIGSOR(I,2)" TO GENERATE AN IDENTIFICATION CODE FOR THAT PAIR

      DO KP=1,2
        DO I=1,ista
           IF(STAT(KP).EQ.FULSTA(I)) THEN
              I01(KP)=I-1
              ISTCAT(KP)=I
              exit
           ENDIF
        END DO
      END DO

      ISEC = SEC + 1.D-7
      WRITE (TIM13,101) IYEAR,IDOY,IHOU,IMIN,ISEC
101   FORMAT (I4.4,I3.3,I2.2,I2.2,I2.2)

      DO I=1,iobs
         IF(TIM13.EQ.BIGSOR(I,2).AND.SOUR.EQ.BIGSOR(I,1)) THEN
            I02= I
            ISOUR=I
            exit
         ENDIF
      END DO

C     GET THE PAIR IDENTIFICATION CODE (CHARACTER*5)

      DO KP=1,2
         WRITE (STAT4(KP),100) I01(KP),I02
100      FORMAT (I2,I4.4)
         NREC(KP)=-1
      END DO

C
C     CHECK IF THE CODE HAS BEEN ALREADY WRITTEN IN THE DIRECT ACCESS
C  FILE.
C     IF IT HAS NOT, ADD THE CODE TO THE TABLE AND WRITE A NEW
C  RECORD ON THE FILE
C

cc      print *, ' nstat = ', nstat

      DO 5 KP=1,2

         DO I=1,nstat

            IF (STAT4(KP).EQ.BIGSTA(I)) THEN
                NREC(KP)=I
                GO TO 5
            ENDIF

            IF (BIGSTA(I).EQ.'      ') THEN
                BIGSTA (I) = STAT4(KP)
                NREC(KP)=I
                IF (I.GT.1) THEN

                   READ (18,REC=1) NRC,ISTC0,ISOU0,TD,PR,HU,KW,CA
                   WRITE (18,REC=1) I,ISTC0,ISOU0,TD,PR,HU,KW,CA

                ENDIF
                WRITE (18,REC=I) I,ISTCAT(KP),ISOUR,999.D0
                GO TO 5


            ENDIF

         END DO

C      THE PROGRAM SHOULD NEVER REACH THIS POINT

      WRITE (*,900)


900   FORMAT (///' ****** WARNING !!! *******'//
     *'  THE PROGRAM MUST NOT REACH THIS POINT AT SUBROUTINE "STADIR"'/
     *'  THERE ARE MORE THAN 7000 STATION OBSERVATIONS. YOU MUST CHANGE'
     */'  THIS PROGRAM TO DEAL WITH THIS POSSIBILITY. NO FURTHER DATA '/
     *'  WILL BE USABLE'//' ENTER Q TO QUIT OR C TO CONTINUE',5X,
     *' (QUIT RECOMMENDED)'//)

      READ (*,400) CCH
400   FORMAT (A1)
      IF (CCH.NE.'C')  ISTAT=990
      IF (ISTAT.NE.0) RETURN

5     CONTINUE

      NREC1=NREC(1)
      NREC2=NREC(2)

      RETURN
      END

************************************************************************

      SUBROUTINE BASDIR (BASEL,SOUR,IYEAR,IDOY,IHOU,IMIN,SEC,NREC4)

C***************************************************************************
C
C      SUBROUTINE BASDIR (BASEL,SOUR,IYEAR,IDOY,IHOU,IMIN,SEC,NREC4)
C
C     THIS SUBROUTINE GENERATES A DIRECT ACCESS FILE INDEXED WITH
C  RESPECT TO THE PAIR BASELINE-TIME, TO ADD BASELINE DEPENDENT
C  CORRECTIONS TO THE OBSERVATIONS. IT USES A TABLE WITH THE
C  IDENTIFICATION CODES OF THOSE PAIRS
C
C   REVISION  1990 NOVEMBER 19   (N.Z.L.)
C   REVISION  1998 JULY 21   (OT)
C   REVISION  1998 September 21   (OT)
C   REVISION  2001 MARCH  13   (OT)   Improvement of cycles
C   LAST REVISION  2003 JUNE 08 (OT) - Changes in OCCAM_N.FI and OPENING.FOR (dictio)
C
C************************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INCLUDE 'OCCAM_N.FI'

      CHARACTER BASFUL*18,BB*8,BBAS*8,STAT1*8,STAT2*8
      CHARACTER BASEL*18,FULSOR*8,FULSTA*8,CCH*1,SOUR*8
      CHARACTER*13 BIGSOR,BIGSTA*6,TIM13
      CHARACTER*6 CSTA1,CSTA2
      COMMON /BSOR1/ BBAS(nobserv)
      COMMON /BASID0/ BASFUL(NBASES),FULSOR(NSOURCES),FULSTA(NSTATIONS),
     * ISTA, ISOR, IOBS
      COMMON /BASID1/ NDIM, nbas
      COMMON /STATUS/ ISTAT
      COMMON /LOCCAT/ CORSOR(NSOURCES,2),CORSTA(NSTATIONS,3),
     * VELSTA(NSTATIONS,3),EPOCH
      COMMON /BSOR/ BIGSOR(nscans,2),BIGSTA(nstat)

C     CHECK THE BASELINE NUMBER IN "BASFUL" AND THE TIME EPOCH
C  NUMBER IN "BIGSOR(I,2)" IN ORDER TO GENERATE AN IDENTIFICATION
C  CODE FOR THAT PAIR

      ISEC = SEC + 1.D-7
      WRITE (TIM13,101) IYEAR,IDOY,IHOU,IMIN,ISEC
101   FORMAT (I4.4,I3.3,I2.2,I2.2,I2.2)

      DO I=1,nbas
         IF(BASEL.EQ.BASFUL(I)) THEN
            I01= I
            exit
         ENDIF
      end do

      DO I=1,iobs
         IF(TIM13.EQ.BIGSOR(I,2).AND.SOUR.EQ.BIGSOR(I,1)) THEN
            I02= I
            ISOUR=I
            exit
         ENDIF
      end do

C     CHECK THE RECORD CORRESPONDING TO EACH STATION AT THE EPOCH IN
C  "STATIM" FILE

      READ (BASEL,251) STAT1,STAT2
251   FORMAT (A8,2X,A8)

      DO I=1,ista
         IF (STAT1.EQ.FULSTA(I)) IST1=I-1
         IF (STAT2.EQ.FULSTA(I)) IST2=I-1
      END DO
      WRITE (CSTA1,151) IST1,I02
      WRITE (CSTA2,151) IST2,I02
151   FORMAT (I2,I4.4)

      DO I=1,nstat               !             nstat = 5000
         IF (CSTA1.EQ.BIGSTA(I)) ISTA1=I
         IF (CSTA2.EQ.BIGSTA(I)) ISTA2=I
      END DO

C     WRITE THE IDENTIFICATION CODE (CHARACTER*5)

      WRITE (BB,131) I01,I02
131   FORMAT (I4.4,I4.4)

C
C     CHECK IF THE PAIR HAS ALREADY BEEN WRITTEN
C     IF NO, ADD THE NEW CODE TO THE TABLE AND WRITE A NEW RECORD ON THE
C  DATA FILE
C

       DO I=1,nobserv
         IF (BB.EQ.BBAS(I)) THEN
            WRITE (*,341) BB
341         FORMAT (/' DUPLICATE OBSERVATION FOUND. CHECK LABEL.. ',A8/)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCC        TESMER 18.04.04       CCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            print*, '( observation no. ',i,' of ',nobserv,' )'
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCC        TESMER 18.04.04       CCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

            NREC4=I
            RETURN
         ENDIF
         IF (BBAS(I).EQ.'        ') THEN
            BBAS(I) = BB
            NREC4=I
            CALL JULDAT(IYEAR,IDOY,IHOU,IMIN,SEC,IDATJ,UTJ,0)
            IF (I.GT.1) THEN

               READ (17,REC=1) NRC,ISOU0,IST0,ISTA0,D0,SD0,R0,SR0,DI0,
     *                      SDI0,RI0,SRI0,time1, time2

               WRITE (17,REC=1) I,ISOU0,IST0,ISTA0,D0,SD0,R0,SR0,DI0,
     *                      SDI0,RI0,SRI0, time1, time2

            ENDIF

            WRITE (17,REC=I) I,ISOUR,ISTA1,ISTA2

            RETURN
         ENDIF
      END DO

C      THE PROGRAM SHOULD NEVER REACH THIS POINT

      WRITE (*,900)
900   FORMAT (///' ****** WARNING !!! *******'//
     *'  THE PROGRAM MUST NOT REACH THIS POINT AT SUBROUTINE "BASDIR"'/
     *'  THERE ARE MORE THAN 30000 DATA POINTS, PRESENT LIMIT OF OCCAM.'
     */'  HOW MANY DAYS WERE THEY OBSERVING ? MY GOD! 30000 POINTS !'/
     *'  YOU MUST CHANGE THIS PROGRAM TO DEAL WITH THIS (IF POSSIBLE!!)'
     */'  NO FURTHER DATA WILL BE USABLE (FURTHER DATA ?!!)'
     *//' ENTER Q TO QUIT OR C TO CONTINUE',5X,' (QUIT RECOMMENDED)'//)
      READ (*,400) CCH
400   FORMAT (A1)
      IF (CCH.NE.'C')  ISTAT=990

      RETURN
      END


*********************************************************************

      SUBROUTINE BUILDEPH (IYEAR)

C-------------------------------------------------------------------
C
C   THIS ROUTINE WILL EXTRACT FROM THE GLOBAL EPHEMERIS AND EOP
C   DATABASES, THE LINES CORRESPONDING TO THE EXPERIMENT AND
C   WRITE THEM TO THE FILE EPHEM.DAT FOR FURTHER USE
C
C   The files to be developed in accordance with information
C   about the leap seconds !!
C
C   REVISION  APRIL 21, 1993   NZL
C   REVISION  1993 NOV. 17  (NZL)  INCLUDE DATA FOR ALL EXPERIMENT
C			 EPH. BUILDING IS MOVED TO THE END OF PROGRAM
C   REVISION APRIL 23, 1997, (OT). New reading of ephemerides.
C   Changing of interpolation step
C   Virtual disk V is canceled
C   REVISION November, 1998, (OT). Increasing of a priori EOP accuracy
C   REVISION April, 13, 1999, (OT). Correction in IDELT
C   REVISION August, 03, 1999, (OT). Format for EPHEM.DAT file
C   REVISION   September, 22, 1999, (OT)
C   REVISION January, 16, 2000, (OT). Y2K problem has been solved
C   REVISION April, 16, 2001, (OT). The option for alternative part
C                                        has been removed
C   REVISION MAY, 26, 2004, (OT). Remark by Thomas Hobiger
C                                     (idelt=17,18,19 were added)
C   REVISION 23, February, 2006, (OT). Leap second on 1.1.2006
C   LAST REVISION 21, January, 2009, (OT). Leap second on 1.1.2009
C
C********************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER EPHNAM*59, EOPNAM*57,  YEAR*2, QASK*1
      CHARACTER*8 SOUR

C   INITIALIZE

C
C   READ INITIAL AND FINAL EPOCH IN DATA FILE
C

      CALL UNFOLD (19,'SORTIM')

      READ (19,REC=1) NT,SOUR,C1,C2,ID1,UT1
      READ (19,REC=NT) NT,SOUR,C1,C2,ID2,UT2

      IF (IYEAR.LT.2000) THEN
         IYEAR_1 = IYEAR-1900
      ELSE
         IYEAR_1 = IYEAR-2000
      END IF

      WRITE (YEAR,'(I2.2)') IYEAR_1

      OPEN (40, FILE='EPHEM.DAT')


C   WRITE DEFAULT PATH FOR THE EPHEMERIS DATA BASES

      WRITE (EPHNAM,100) IYEAR_1
100   FORMAT ('../EPHEM/EPHEM.',I2.2)
      WRITE (EOPNAM,101) IYEAR_1
101   FORMAT ('../EPHEM/EOP.',I2.2)

      OPEN (41,FILE=EPHNAM,STATUS='OLD',ERR=900)
      OPEN (42,FILE=EOPNAM,STATUS='OLD',ERR=901)


C  READ EPHEMERIS FILE AND EXTRACT THE 9 CLOSEST EPOCHS

      UT1 = UT1/(2.d0*3.1415926535897932d0)
      UT2 = UT2/(2.d0*3.1415926535897932d0)

      IF (UT1.GT.0.5) AJDINIT = ID1-1.D0
      IF (UT1.LE.0.5) AJDINIT = ID1-1.5D0

      IF (UT2.GT.0.5) AJDFIN = ID2+2.5D0
      IF (UT2.LE.0.5) AJDFIN = ID2+2.D0

1     CONTINUE

c
c     For DELE403 program we need only positions and velocities
c     for the Earth ephemerides (not acceleration).
c

      READ (41,*,END=910) AJDEP1,EX,EY,EZ,EXD,EYD,EZD
      READ (41,*,END=910) AJDEP2,SX,SY,SZ,SXD,SYD,SZD
      READ (41,*,END=910) AJDEP3,AMX,AMY,AMZ,AMXD,AMYD,AMZD

      if (ajdep1.ne.ajdep2) write (*,*) 'ERROR !!!   CHECK AN EPOCH ! '
      if (ajdep1.ne.ajdep3) write (*,*) 'ERROR !!!   CHECK AN EPOCH ! '

      ajdep=ajdep1

      IF (AJDEP.LT.AJDINIT) THEN
         GO TO 1
      ELSE IF (AJDEP.GT.AJDFIN) THEN
         GO TO 599
      ELSE

C  READ EPHEMERIS FOR EARTH, SUN AND MOON

      WRITE (40,120) AJDEP,EX,EY,EZ,EXD,EYD,EZD
120   FORMAT ('EARTD ',F7.1,6D24.16)

      WRITE (40,121) AJDEP,SX,SY,SZ,SXD,SYD,SZD
121   FORMAT ('SUND  ',F7.1,6D24.16)

      WRITE (40,122) AJDEP,AMX,AMY,AMZ,AMXD,AMYD,AMZD
122   FORMAT ('MOOND ',F7.1,6D24.16)

      GO TO 1

      ENDIF

C   END OF READING EPHEMERIS

C   START READING EOP VALUES. FIND 9 CLOSEST


599   continue

c
c     The operators are OK only for EOP1 series with 1-day intervals
c     For 5-day intervals the "4" should be changed on "20".
c

      AJDINIT = ID1 - 4.D0
      AJDFIN = ID2 + 4.D0

2     CONTINUE
      READ (42,124,END=920) AJDEP,XP,YP,UT1UTC,DP,DE

      IF (AJDEP.LT.AJDINIT) THEN
         GO TO 2
      ELSE IF (AJDEP.GT.AJDFIN) THEN
         GO TO 99

      ELSE

         if (ajdep.ge.43509.d0.and.ajdep.le.43873.d0)
     *            idelt = 17
         if (ajdep.ge.43874.d0.and.ajdep.le.44238.d0)
     *            idelt = 18
         if (ajdep.ge.44239.d0.and.ajdep.le.44785.d0)       !  1-January-1980
     *            idelt = 19
         if (ajdep.ge.44786.d0.and.ajdep.le.45150.d0)       !  1-July-1980
     *            idelt = 20
         if (ajdep.ge.45151.d0.and.ajdep.le.45515.d0)       !  1-July-1981
     *            idelt = 21
         if (ajdep.ge.45516.d0.and.ajdep.le.46246.d0)       !  1-July-1983
     *            idelt = 22
         if (ajdep.ge.46247.d0.and.ajdep.le.47160.d0)
     *            idelt = 23
         if (ajdep.ge.47161.d0.and.ajdep.le.47891.d0)       !  1-January-1988
     *            idelt = 24
         if (ajdep.ge.47892.d0.and.ajdep.le.48256.d0)       !  1-January-1990
     *            idelt = 25
         if (ajdep.ge.48257.d0.and.ajdep.le.48803.d0)       !  1-January-1991
     *            idelt = 26
         if (ajdep.ge.48804.d0.and.ajdep.le.49168.d0)       !  1-July-1992
     *            idelt = 27
         if (ajdep.ge.49169.d0.and.ajdep.le.49533.d0)       !  1-July-1993
     *            idelt = 28
         if (ajdep.ge.49534.d0.and.ajdep.le.50082.d0)       !  1-July-1994 to 1-January-1996
     *            idelt = 29
         if (ajdep.ge.50083.d0.and.ajdep.le.50629.d0)       !  96-97
     *            idelt = 30
         if (ajdep.ge.50630.d0.and.ajdep.le.51178.d0)       !  97-99
     *            idelt = 31
         if (ajdep.ge.51179.d0.and.ajdep.le.53735.d0)       !  1-January-1999- 1-January-2005
     *            idelt = 32
         if (ajdep.ge.53736.d0.and.ajdep.le.54831.d0)   ! 1-January-2005  Change on 21.01.2009
     *            idelt = 33                            ! By Oleg Titov
         if (ajdep.ge.54832.d0.and.ajdep.le.56108.d0)   ! Change on 07.05.2012
     *            idelt = 34                            ! By Oleg Titov
         if (ajdep.ge.56109.d0.and.ajdep.le.57203.d0)   ! Change on 16.10.2015
     *            idelt = 35                            ! By Oleg Titov
         if (ajdep.ge.57204.d0.and.ajdep.le.57753.d0)   ! Change on 22.01.2017
     *            idelt = 36                            ! By Oleg Titov
         if (ajdep.ge.57753.d0)                         ! Change on 22.01.2017
     *            idelt = 37                            ! By Oleg Titov

      WRITE (40,125) AJDEP,XP,YP,UT1UTC,idelt,DP,DE

c124   FORMAT (9x,f7.0,f9.6,f9.6,f10.7,14x,f9.5,f9.5)
c125   FORMAT ('EOP1  ',F7.1,1X,F8.6,F8.6,1X,F10.7,i4,f8.5,f8.5)

c
c  Change of reading/writing format of the EOP files  04.01.2008 (OT)
c

124   FORMAT (12x,f7.0,f11.6,f11.6,f12.7,12x,f11.6,f11.6)

125   FORMAT ('EOP1  ',F7.1,1X,F8.6,F8.6,1X,F10.7,i4,f9.6,f9.6)

      GO TO 2

      ENDIF

C   END OF READING EPHEMERIS
******************************************************************

99    CLOSE (40)
      CLOSE (41)
      CLOSE (42)

      RETURN

************************** ERROR MESSAGES ***********************

900   WRITE (*,*) ' ***** ERROR FOUND OPENING FILE ',EPHNAM
      GO TO 902
901   WRITE (*,*) ' ***** ERROR FOUND OPENING FILE ',EOPNAM

902   WRITE (*,103)
103   FORMAT (/' C-CONTINUE, Q-QUIT'/' -')

      READ (*,'(A1)') QASK
      IF (QASK.EQ.'C'.OR.QASK.EQ.'c') RETURN
      IF (QASK.EQ.'Q'.OR.QASK.EQ.'q') STOP
      GO TO 902

910   WRITE (*,*) ' ***** FILE ',EPHNAM,' DOES NOT CONTAIN THE'
      WRITE (*,*) '       EPHEMERIS FOR THE DATES ',ID1,' - ',ID2
      GO TO 902

920   WRITE (*,*) ' ***** FILE ',EOPNAM,' DOES NOT CONTAIN THE'
      WRITE (*,*) '       EARTH OR. PAR. FOR THE DATE ',ID1,' - ',ID2
      GO TO 902

      END

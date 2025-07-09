       PROGRAM DTHAND 

C      PROGRAM DTHAND  (DATA HANDLER)
C
C      THIS VERSION IS INTENDED TO BE COMPATIBLE WITH ALL COMPILERS OF
C   FORTRAN 77. 
C
C      THIS VERSION IS READY TO WORK ON A PC.
C
C      THIS PROGRAM IS DESIGNED FOR READING NGS-FORMAT VLBI DATA FILES,
C   SELECT THE DIFFERENT BASELINES AND CONVERT THE DATA FORMAT TO
C   OCCAM FORM
C
C      THE PROGRAM HAVE THE NEXT PARTS  (1990-NOV-16) :
C
C      1. ASK FOR AN IDENTIFICATION OF THE EXPERIMENT AND USER
C          1.A.- GET THE RUN CODE
C
C      2. ASK FOR THE CATALOG FILE TO BE USED
C
C      3. ASK FOR THE NAME OF THE DATA FILE CONTAINING THE DATA CARDS
C        FROM NGS-FORMAT (THE REST IS IGNORED)
C
C      4. CHECK IF ALL THE SOURCES AND STATIONS IN THE DATA FILE ARE IN
C        THE CATALOG.
C           4.A. A TABLE WITH THE COORDINATES IS CREATED FOR THE SOURCES
C               AND STATIONS OF THE EXPERIMENT.
C
C      6. SELECTION OF WHICH STATIONS SHOULD BE RUN
C
C      7. LOOK FOR THE DATA CARDS CORRESPONDING TO EACH BASELINE
C        SELECTED AND CONVERT THOSE CARDS INTO OCCAM FORMAT
C
C      8. CREATE DIRECT ACCESS FILES CONTAINING:
C           I.  SOURCES AND TIMES (PRECESION AND NUTATION)
C           II. STATIONS AND TIMES (STATION DEPENDENT CORRECTIONS)
C           III. BASELINES AND TIMES (BASELINE DEPENDENT CORRECTIONS)
C
C      THESE FILES ARE THE INPUT FOR INDEPENDENT PROGRAMS WHICH
C      CALCULATE THE CORRECTIONS TO BE READ BY KVLBI, KVLBI_IE, COLL, LSM
C
C
C     REQUIREMENTS TO RUN DTHANDLR:
C
C      . INPUT DATA FILE ON NGS-FORMAT.
C
C      . A CATALOG FILE, MUST INCLUDE INFORMATION ABOUT:
C
C        - HISTORY OF THE CATALOG FILE. CODE HH
C        - STATIONS COORDINATES. CODE TC  (TERRESTRIAL COORDINATES)
C        - SOURCES COORDINATES. CODE CC  (CELESTIAL COORDINATES)
C        - INFORMATION ABOUT THE FORMATS OF THE COORDINATES
C
C      OTHER POSSIBILITIES
C
C       . SELECT OBSERVATIONS BASED ON THE DATA QUALITY FLAG OR FORMAL
C        ERRORS OF THE OBSERVATIONS.
C        NOW THE PROGRAM ONLY CONSIDERS THE DATA WITH QUALITY FLAG 0
C
C   ISTAT.- PROGRAM'S STATUS FLAG:
C        ISTAT = 0     EVERYTHING IS WORKING PROPERLY
C        ISTAT = 990   PROGRAM FORCED TO QUIT
C        ISTAT = 1000  IN SOME ROUTINES IT SHOWS RUN TIME ERRORS
C
C
C REVISION  1990 NOVEMBER 19  (N.Z.L.)
C REVISION  1993 NOVEMBER 17  (NZL) EPHEMERIS BUILT AS PART 9
C REVISION  1997 MAY 15   (OT)  THE NUMBER OF STATIONS AND SOURCES
C                                   HAS BEEN INCREASED
C Revision 1999, September, 22 (OT)
C REVISION  2000 July  13 (OT) No baselines choice, list of stations
C                              was added for screen output
C REVISION 2001  March, 11 (OT)    COMMON /LOGIC/ was deleted
C                                  Subroutine CORROP was deleted
C
C REVISION 2001 April, 08 (OT)  Subroutine SHORTB was deleted
C                                    Subroutine SELBAS was deleted
C REVISION  2002  MAY 16  (OT) - Observations of sources out of
C                                    catalogue are deleted
C REVISION  2003  APRIL 16  (OT) - Position break
C REVISION  2003 JUNE 08 (OT) - Changes in OCCAM_N.FI and OPENING.FOR (dictio)
C REVISION  2005 MAY 04 (OT) - Reading and writing of the DBC code for
C                                   SINEX file submission
C REVISION  2006 AUGUST 30 (OT) - S1 declared as character
C
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
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INCLUDE 'OCCAM_N.FI'

      CHARACTER CATNAM*12,FILNAM*12,EXPID*8,BASFUL*18,s*1,s1*1
      CHARACTER FULSOR*8,FULSTA*8,name_no(400)*8      !  Tesmer - up to 300 sources
      CHARACTER mo*3, dbc_code*2

      COMMON /BASID0/ BASFUL(nbases),FULSOR(nsources),
     *  FULSTA(nstations),ISTA, ISOR, IOBS
      COMMON /BASID1/ NDIM, nbas
      COMMON /STATUS/ ISTAT
      COMMON /CATALG/ CATNAM
c+ jboehm, 2002, July 29
      common /lsmopt/ filnam
c- jboehm

C------------------------
C   GET PARAMETERS
C------------------------

      open (1,file='dbc_code',status='unknown')
      CALL PARAM

C
C     INITIALIZE THE STATUS VARIABLE ISTAT
C
      ISTAT = 0

C
C     PART  1. ASK FOR USER'S NAME, DATE AND TIME
C

      WRITE (*,100)
100   FORMAT (/' PLEASE, ENTER AN IDENTIFIER FOR THIS RUN  (UP TO '
     *,'8 CHARACTERS'/' --------',15X,'HIT "Q" TO QUIT')
      READ (*,200) EXPID
200   FORMAT (A8)
      IF (EXPID.EQ.'q       '.OR.EXPID.EQ.'Q       ') THEN
         ISTAT=990
         GO TO 990
      ENDIF

C
C     PART  2.  ASK FOR THE CATALOG FILE
C

22    CONTINUE
      WRITE (*,102)
102   FORMAT (/' PLEASE, ENTER THE NAME OF THE CATALOG FILE TO BE ',
     *'USED'/' --------.---',15X,' HIT "Q" TO QUIT')
      READ (*,201,ERR=22) CATNAM

201   FORMAT (A12)
      print 201, CATNAM


      IF (CATNAM.EQ.'            ') GO TO 22
      IF (CATNAM.EQ.'Q           '.OR.CATNAM.EQ.'q           ') THEN
          ISTAT=990
          GO TO 990
      ENDIF

c      OPEN (7,FILE=CATNAM,STATUS='OLD',ERR=1001)
c      CLOSE (7)

C
C     PART  3.  ASK FOR THE NAME OF THE DATA FILE.
C

23    CONTINUE
      WRITE (*,103)
103   FORMAT (/' PLEASE, ENTER THE DATA FILE NAME'/
     *' --------.---',15X,' HIT "Q" TO QUIT')

      READ (*,201,ERR=23) FILNAM

      IF (FILNAM.EQ.'            ') GO TO 23
      IF (FILNAM.EQ.'Q           '.OR.FILNAM.EQ.'q           ') THEN
          ISTAT=990
          GO TO 990
      ENDIF

      OPEN (4, FILE=FILNAM,STATUS='OLD',ERR=1000)

C  2004 June 14, jboehm s1 added  to read Zinovys ngs files
      read (4, 97) s1,s                                    !  28.04.2003
      rewind 4
      if (s.eq.'$') then
        read (4, 98) iye, mo, idd, dbc_code                !  28.04.2003
      else
        if (s1.eq.'B') then
          read (4, 96) iye, mo, idd, dbc_code
         else                                              !  16.04.2003
          read (4, 99) iye, mo, idd, dbc_code
        end if
      end if
      write (1,95) iye, mo, idd, dbc_code                  !  OT, 04.05.2005
      print *, iye, mo, idd

95    format (i2,a3,i2,a2)
96    format (33x,i2,a3,i2,a2)
97    format (28x,a1,5x,a1)
98    format (35x,i2,a3,i2,a2)
99    format (34x,i2,a3,i2,a2)

      if(mo.eq.'JAN') m1=1
      if(mo.eq.'FEB') m1=2
      if(mo.eq.'MAR') m1=3
      if(mo.eq.'APR') m1=4
      if(mo.eq.'MAY') m1=5
      if(mo.eq.'JUN') m1=6
      if(mo.eq.'JUL') m1=7
      if(mo.eq.'AUG') m1=8
      if(mo.eq.'SEP') m1=9
      if(mo.eq.'OCT') m1=10
      if(mo.eq.'NOV') m1=11
      if(mo.eq.'DEC') m1=12
      CALL DDOY (IYE,M1,IDD,0,0,IT,IDOY0)

      if (iye.gt.70) iye=1900+iye
      if (iye.le.70) iye=2000+iye

C
C     PART  4. CHECK SOURCES AND STATIONS WITH THE CATALOG
C

      CALL CHKCAT (NOBS,name_no,i_no,iye,idoy0,expid)   ! Change by Volker

C     CHECK THE STATUS RETURNED BY THE ROUTINE

      IF (ISTAT.EQ.990) GO TO 9999

C
C     PART  5.  CONVERSION TO OCCAM FORMAT
C

      CALL CONVER (EXPID,NOBS,IYEAR,name_no,i_no)

C
C     PART  6. BUILD THE EPHEMERIS FILE FOR THE EXPERIMENT
C

      CALL BUILDEPH (IYEAR)

C    CHECK THE STATUS CODE RETURNED BY THE ROUTINE

      IF (ISTAT.EQ.990) GO TO 990
      IF (ISTAT.NE.0) THEN
          WRITE (*,*)  ' ISTATUS = ',ISTAT
          STOP
      ENDIF

      GO TO 999

C
C     WARNING MESSAGES FROM THE PROGRAM TO THE USER
C

1000  WRITE (*,900) FILNAM
900   FORMAT(/' RUN TIME ERROR OPENING FILE... ',A12//
     *' PROBABLY YOUR DATA FILE HAS A DIFFERENT NAME.'/
     *' CHECK THAT NAME AND TRY AGAIN'/)
      GO TO 23

1001  WRITE (*,901) CATNAM
901   FORMAT(/' RUN TIME ERROR OPENING FILE... ',A12//
     *' PROBABLY YOUR CATALOG FILE HAS A DIFFERENT NAME.'/
     *' CHECK THAT NAME AND TRY AGAIN'/)
      GO TO 22

990   WRITE (*,991)
991   FORMAT (///' *** PROGRAM TERMINATED BY USER'///)
      GO TO 9999
999   WRITE (*,998)
998   FORMAT (///' PROGRAM TERMINATED SUCCESFULLY'///
     *' COME HERE AGAIN IF YOU HAVE MORE DATA TO CONVERT'///)
9999  CLOSE (4)
      CLOSE (5)
      STOP
      END

************************************************************************

      SUBROUTINE DDOY (IY,IM,ID,IH,IMI,ITIM,IDOY)
C
C     SUBROUTINE DDOY (IY,IM,ID,IH,IMI,ITIM,IDOY)
C
C     THIS SUBROUTINE COMPUTES THE DAY OF YEAR AND THE TIME IN A SPECIAL
C  FORMAT...  11 h  20 M. ---  113  (11.3 * 10)
C     USEFUL TO CREATE THE RUN IDENTIFICATION CODE AND FOR THE CONVERSION
C  TO OCCAM FORMAT.
C
C    LATEST REVISION  1990 NOVEMBER 19  (N.Z.L.)
C    Last revision 1999, September, 22 (OT)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      IFLG=0
      IF(MOD(IY,4).EQ.0) IFLG=1
      IF(IM.EQ.1) IDOY=ID
      IF(IM.EQ.2) IDOY=31+ID
      IF(IM.EQ.3) IDOY=59+ID+IFLG
      IF(IM.EQ.4) IDOY=90+ID+IFLG
      IF(IM.EQ.5) IDOY=120+ID+IFLG
      IF(IM.EQ.6) IDOY=151+ID+IFLG
      IF(IM.EQ.7) IDOY=181+ID+IFLG
      IF(IM.EQ.8) IDOY=212+ID+IFLG
      IF(IM.EQ.9) IDOY=243+ID+IFLG
      IF(IM.EQ.10) IDOY=273+ID+IFLG
      IF(IM.EQ.11) IDOY=304+ID+IFLG
      IF(IM.EQ.12) IDOY=334+ID+IFLG
      TIM=(IH+IMI/60.d0)*10.d0
      ITIM=TIM
      RETURN
      END

************************************************************************
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCC tesmer 17.04.04 CCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE CHKCAT (NOBS,name_no,i_no,iye,idoy0,expid)

C
C     SUBROUTINE CHKCAT
C
C     THIS SUBROUTINE CHECKS IF THE SOURCES AND STATIONS CONTAINED IN
C   THE DATA FILE ARE INCLUDED IN THE CATALOG, WITH THE HELP OF THE
C   SUBROUTINE CHCK
C
C   REVISION  NOVEMBER 15, 1993  (NZL) ADD VELOCITIES IN CATALOG
C   Revision 1999, September, 22 (OT)
C   REVISION 2001 April, 08 (OT)  Subroutine SELBAS was deleted
C                                 Array BASFIL is generating here
C   REVISION 2001 May, 15 (OT)
C   REVISION  2002  MAY 16  (OT) - Observations of sources out of
C                                       catalogue are deleted
C   REVISION  2003  APRIL 16  (OT) - Position break
C   LAST REVISION  2006  FEBRUARY 23  (OT) - Volker's change
C          The part of automatic creation of the KVLBI.OPT files revised
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INCLUDE 'OCCAM_N.FI'

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCC TESMER 22.09.04 allow up to 300 sources CCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      CHARACTER*8 STA1,STA2,SOUR,FULSTA(nstations),FULSOR(nsources),
     *     REF_ST, sta_list(nstations), BLANK, name_no(400)

      CHARACTER CATNAM*12, QQ*1, CC*1,AXTYPE*4, AXISTY(nstations)*4
      CHARACTER*18 BASFUL, bass, BLNK18
      CHARACTER*4 stcode, stc(nstations)
      DIMENSION nst(100)            !   The number of stations in catalogue
      LOGICAL NOVEL, AUTO, AGRAD
      DIMENSION CORSOR(NSOURCES,2), CORSTA(NSTATIONS,3),
     * VELSTA(NSTATIONS,3), OFFS(NSTATIONS)          !   Change for OFFS on 31-Oct-2006

c+ jboehm, 2002, July 29
      character filnam*12
c- jboehm
      CHARACTER EXPID*8   !  Volker Tesmer

      COMMON /CATALG/ CATNAM
      COMMON /STATUS/ ISTAT
      COMMON /COORD/ RASOR,DECSOR,XCOR,YCOR,ZCOR,XVEL,YVEL,ZVEL
      COMMON /BASID0/ BASFUL(NBASES),FULSOR,FULSTA, ISTA, ISOR, IOBS
      COMMON /BASID1/ NDIM, nbas
      COMMON /LOCCAT/ CORSOR,CORSTA,VELSTA,EPOCH

c+ jboehm, 2002, July 29
      common /lsmopt/ filnam
c- jboehm

C     INITIALIZE THE CHARACTER VARIABLES

      BLANK='        '
      BLNK18='                  '


c      open (433,file='dtau0.dat',status='unknown')

      DO  I=1,nstations
         FULSTA(I)=BLANK
      END DO

      DO  I=1,nsources
         FULSOR(I)=BLANK
      END DO

      DO I=1,nstations
         DO J=1,3
            CORSTA(I,J) = 0.d0
            VELSTA(I,J) = 0.d0
         ENDDO
      ENDDO

      DO I=1,NBASES
         BASFUL(I)=BLNK18
      ENDDO


      REWIND 4

C   GET THE EPOCH OF THE CATALOG

      NOVEL = .FALSE.

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCC tesmer 17.04.04 CCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      CALL CHCK ('        ','EP',EPOCH,AXTYPE,OFFST1,stcode,
     * name_no,i_no,iye,idoy0,br,expid,npp)

      IF (EPOCH.EQ.0.d0) THEN
         NOVEL=.TRUE.
         WRITE (*,341)
341      FORMAT (' *** NO EPOCH FOR THE CATALOG. NO VELOCITY',
     *   ' MODEL WILL BE APPLIED ***'/'  HIT ENTER TO CONTINUE')
         READ (*,'(I1)') ICONT
      ENDIF

C     THE PROCEDURE OF CHECKING WILL BE DONE ONLY ONCE FOR EACH STATION
C  OR SOURCE. THIS SUBROUTINE WILL CHECK THAT THE NUMBER OF STATIONS IS
C  NOT GREATER THAN 20 AND THAT THE SOURCES ARE LESS THAN 200, AND WILL
C  CALL THE SUBROUTINE CHCK TO DO THE ACTUAL CHECKING WITH THE CATALOG
C

      NOBS=0
      ISTA=0
    
      ISOR=0
      ndim=0
      i_no=0

      

1     CONTINUE

         READ (4,300,END=10) STA1,STA2,SOUR,ICODE
300      FORMAT (A8,2X,A8,2X,A8,50X,I2)

c         print 9000, STA1, STA2
9000     format (a8,2x,a8)

         IF (ICODE.NE.1) GO TO 1

         NOBS = NOBS+1
         IOBS = NOBS

C     CHECK STATIONS, CREATE A TABLE WITH NAMES AND A TABLE WITH COORDS.

         IFL1 = 0

         DO  I=1,nstations
            IF (STA1.EQ.FULSTA(I)) then
               IFL1=1
c               print  9000, Sta1, fulsta(i)
               

            end if

c            print *, ' ista = ', ista, nobs, sta1, sta2

         end do

         IF (IFL1.EQ.0) THEN

      
            ISTA=ISTA+1

c            print *, ifl1, ista

            IF (ISTA.GT.nstations) GO TO 1000
            FULSTA(ISTA)=STA1

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCC tesmer 17.04.04 CCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            CALL CHCK (STA1,'TC',epoch,AXTYPE,OFFST1,stcode,
     *                   name_no,i_no,iye,idoy0,br,expid,npp)

            nst(ista)=npp

            if (br.eq.1.d0) then

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCC tesmer 17.04.04 CCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
              CALL CHCK (STA1,'TB',epoch,AXTYPE,OFFST1,stcode,
     *                   name_no,i_no,iye,idoy0,br,expid,npp)

            end if

            IF (ISTAT.NE.0) GO TO 990

            CORSTA(ISTA,1)=XCOR
            CORSTA(ISTA,2)=YCOR
            CORSTA(ISTA,3)=ZCOR
            AXISTY(ISTA)=AXTYPE
            OFFS(ISTA)=OFFST1
            stc(ISTA)=stcode

            IF (.NOT.NOVEL) THEN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCC tesmer 17.04.04 CCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
               CALL CHCK (STA1,'TV',epoch,AXTYPE,OFFST1,stcode,
     *       name_no,i_no,iye,idoy0,br,expid,npp)

               IF (ISTAT.NE.0) GO TO 990
               VELSTA (ISTA,1) = XVEL
               VELSTA (ISTA,2) = YVEL
               VELSTA (ISTA,3) = ZVEL
            ENDIF

         END IF

         IFL2=0

         DO  I=1,nstations
            IF (STA2.EQ.FULSTA(I)) then
               IFL2=1
            end if
         END DO

         IF (IFL2.EQ.0) THEN

            ISTA=ISTA+1
            IF (ISTA.GT.nstations) GO TO 1000
            FULSTA(ISTA)=STA2

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCC tesmer 17.04.04 CCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            CALL CHCK (STA2,'TC',epoch,AXTYPE,OFFST1,stcode,
     *                      name_no,i_no,iye,idoy0,br,expid,npp)
            nst(ista)=npp

            if (br.eq.1.d0) then

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCC tesmer 17.04.04 CCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
               CALL CHCK (STA2,'TB',epoch,AXTYPE,OFFST1,stcode,
     *                      name_no,i_no,iye,idoy0,br,expid,npp)

            end if

            IF (ISTAT.NE.0) GO TO 990

            CORSTA(ISTA,1)=XCOR
            CORSTA(ISTA,2)=YCOR
            CORSTA(ISTA,3)=ZCOR
            AXISTY(ISTA)=AXTYPE
            OFFS(ISTA)=OFFST1
            stc(ISTA)=stcode

            IF (.NOT.NOVEL) THEN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCC tesmer 17.04.04 CCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
               CALL CHCK (STA2,'TV',epoch,AXTYPE,OFFST1,stcode,
     *                 name_no,i_no,iye,idoy0,br,expid,npp)

               IF (ISTAT.NE.0) GO TO 990
               VELSTA (ISTA,1) = XVEL
               VELSTA (ISTA,2) = YVEL
               VELSTA (ISTA,3) = ZVEL
            END IF

         ENDIF

         bass(1:8)=sta1
         bass(9:10)='  '
         bass(11:18)=sta2

         nbas = ista*(ista-1)/2
cc         print *, ' nbas  = ', nbas

         DO  I=1,nbas
            IF (bass.eq.basful(i)) go to 323
         END DO

         
         ndim=ndim+1
         basful(ndim)=bass

 323     continue

C     CHECK SOURCES, CREATE A TABLE WITH NAMES AND A TABLE WITH COORDS.

         islg=0
      
         DO  I=1,nsources
            IF (SOUR.EQ.FULSOR(I)) then
               islg=1
               exit
            end if
         END DO

         if (islg.eq.0) then

            ISOR=ISOR+1

c            print *, 'isor = ', isor, nsources

            IF (ISOR.GT.nsources) GO TO 1001
            FULSOR(ISOR)=SOUR

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCC tesmer 17.04.04 CCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

            CALL CHCK (SOUR,'CC',epoch,AXTYPE,OFFST1,stcode,
     *                   name_no,i_no,iye,idoy0,br,expid,npp)

            tt=dble(iye)+dble(idoy0)/365.2422

            IF (ISTAT.EQ.0) THEN
               CORSOR(ISOR,1) = RASOR
               CORSOR(ISOR,2) = DECSOR

c               write (433,444) dcos(rasor)*dcos(decsor),
c     *          dsin(rasor)*dcos(decsor), dsin(decsor)
c 444           format(3(2x,d22.16))

c
c   The model for radiosource
c   can be used if apriori positional variations are implemented
c
c            call unstable_source (tt, isor, sour, corsor)

c            print *, corsor(isor,1), corsor(isor,2)

c            call anom_acceleration (tt, isor, sour, corsor)

            call gal_acceleration (tt, isor, sour, corsor)

c            print *, corsor(isor,1), corsor(isor,2)

c            call seven_parameters (tt, isor, sour, corsor)
c            call second_polynom (tt, isor, sou	r, corsor)
c            call anom_acceleration_advanced (tt, isor, sour, corsor)
c            call grav_waves (tt, isor, sour, corsor)


            END IF

         end if

990      CONTINUE

         IF (ISTAT.EQ.990) RETURN

      GO TO 1

C     CREATION OF A DIRECT ACCESS FILE FOR THE STATIONS WITH THEIR
C   COORDINATES

10    CONTINUE

      CALL UNFOLD(20,'STACAT')

      WRITE (*,344)
344   FORMAT (/'  LIST OF VLBI STATIONS AT THE EXPERIMENT  '/)

      DO IJ=1,ISTA

          WRITE (*,345) IJ, FULSTA(IJ), NST(IJ)
  345     FORMAT (1x,I3,2x,a8,2X,I6)
          IJK=IJ
          IF (IJ.EQ.1) IJK=ISTA
          WRITE (20,REC=IJ) IJK,FULSTA(IJ),CORSTA(IJ,1),CORSTA(IJ,2),
     *    CORSTA(IJ,3),EPOCH,VELSTA(IJ,1),VELSTA(IJ,2),VELSTA(IJ,3),
     *    AXISTY(IJ),OFFS(IJ),stc(ij),NST(IJ)

      END DO

      CLOSE (20)

      GO TO 19

C     WARNING MESSAGES FROM THE PROGRAM TO THE USER

1000  WRITE (*,100)
100   FORMAT (///' ***** WARNING !!! *****'//' YOUR DATA FILE CONTAINS'
     *,' MORE THAN 20 STATIONS, LIMIT OF THIS PROGRAM'//' YOU HAVE TO'
     *,' REDUCE THE NUMBER OF STATIONS TO 20 BEFORE RUNNING ME AGAIN')
      ISTAT=990
      RETURN

1001  WRITE (*,101)  
101   FORMAT (///' ***** WARNING !!! *****'//' YOUR DATA FILE CONTAINS'
     *,' MORE THAN 200 SOURCES, LIMIT OF THIS PROGRAM'//' YOU HAVE TO'
     *,' REDUCE THE NUMBER OF SOURCES TO 100 BEFORE RUNNING ME AGAIN')
      write (*,*) isor, nsources  

      ISTAT=990
      RETURN

19    IF (ISTAT.EQ.1000) THEN
          WRITE (*,102)
102   FORMAT (////' *** ERROR MESSAGES FROM THE CHECKING ROUTINE ***'//
     *' IF THEY ARE NOT CORRECTED, SERIOUS PROBLEMS MAY HAPPEN ',
     *' WITH THE PERFORMANCE OF THIS PROGRAM'/
     *' HIT "Q" TO QUIT THE PROGRAM AND CORRECT THE MISSING ITEMS',
     *' (STRONGLY RECOMMENDED)'/' IF YOU JUST WANT TO ACCEPT THE',
     *' RISK, HIT "C"')
          READ (*,201) CC
201   FORMAT (A1)
          ISTAT=990
          IF (CC.EQ.'C'.OR.CC.EQ.'c') ISTAT=0
      ENDIF
      RETURN
      END


************************************************************************

      SUBROUTINE CHCK (NAME,CODE,EPOCH,AXTYPE,OFFST1,stcode,
     * name_no,i_no,iye,idoy0,br,expid, nst)                     

C
C     SUBROUTINE CHCK
C
C     THIS SUBROUTINE CHECKS STATIONS AND SOURCES WITH THE CATALOG TO
C  BE SURE THAT WE HAVE THEIR COORDINATES, AND SEND THESE TO THE
C  CALLING ROUTINE IN ORDER TO GENERATE A TABLE OF THEM.
C
C  CODES AVAILABLE
C    'EP' EPOCH OF COORDINATES
C    'TB' TERRESTRIAL COORDINATES OF STATIONS  before break
C    'TC' TERRESTRIAL COORDINATES OF STATIONS
C    'TV' TERRESTRIAL VELOCITIES OF STATIONS
C    'CC' CELESTIAL COORDINATES OF SOURCES
C
C   REVISION  1990 NOVEMBER 19  (N.Z.L.)
C   REVISION  1993 NOVEMBER 15  (N.Z.L.) ADD VELOCITIES AND EPOCH
C   REVISION  1999 JANUARY  16  (O.T.) FORMAT OF RADIOSOURCES READING
C   Revision 1999, September, 22 (OT)
C   REVISION  2002  MAY 16  (OT) - Observations of sources out of
C                                       catalogue are deleted
C   REVISION  2003  APRIL 16  (OT) - Position breaks (change in *.CAT also)
C   LAST REVISION  2006  FEBRUARY 23  (OT) - All changes for OCCAM 6.2 version
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      CHARACTER CODE*2,NAME*8,LCODE*2,LNAME*8,C0*1,CATNAM*12,LINE*96
      CHARACTER SIGN*1,AXTYPE*4, name_no(400)*8, stcode*4            ! Change by Volker
      CHARACTER EXPID*8                                              ! Change by Volker

      COMMON /CATALG/ CATNAM
      COMMON /COORD/ RA,DEC,XCOR,YCOR,ZCOR,XVEL,YVEL,ZVEL
      COMMON /STATUS/ ISTAT
      COMMON /MATH/ PI, TWOPI, ROG, ROS, ROZ, TOL

C     OPEN CATALOG FILE AND LOOK FOR THE ITEM

      OPEN (9,FILE=CATNAM,STATUS='OLD',ERR=1002)

      OPEN (2119,access='append',FILE='NEW_SOUR',STATUS='unknown')
      OPEN (2219,FILE='list_of_source',STATUS='unknown')

      INUM=0
      nstation = 0

      icode_s = 0      

      br=0.d0    !   16.04.2003

1     CONTINUE

      READ (9,200,END=1000,ERR=1001) LINE
      READ (LINE,210) LCODE,LNAME
200   FORMAT (A96)
210   FORMAT (A2,1X,A8)



      INUM=INUM+1
 
      if (LCODE.eq.'CC') icode_s = icode_s + 1

      IF (LCODE.NE.CODE) GO TO 1

      IF (LCODE.EQ.'EP') THEN
         READ (LINE,212) EPOCH
212      FORMAT(3X,f8.3)
         GO TO 999
      ENDIF

      IF (LCODE.EQ.'TC'.OR.LCODE.EQ.'TB') THEN
          nstation = nstation + 1
          IF (LNAME.EQ.NAME) then
		    nst = nstation
             nstation = 0
          end if
      end if

C
C     HERE WE READ THE COORDINATES OF THE RADIO SOURCES
C    AND PUT THEM IN A TABLE



      IF (LNAME.EQ.NAME) THEN
C
C     HERE WE READ THE COORDINATES OF THE RADIO SOURCES
C    AND PUT THEM IN A TABLE
C
C     SOMETHING SIMILAR IS DONE WITH THE STATIONS COORDINATES
C

       IF (LCODE.EQ.'CC') THEN

            READ (LINE,211) IRAH,IRAM,RAS,SIGN,IDECD,IDECM,DECS
            RA = (IRAH+IRAM/60.D0+RAS/3600.D0)*2.D0*PI/24.D0
            DEC = (IDECD+IDECM/60.D0+DECS/3600.D0)*2.D0*PI/360.D0
            IF (SIGN.EQ.'-') DEC=-DEC
211         FORMAT (12X,I2,1X,I2,1X,F12.9,1X,A1,I2,1X,I2,1X,F11.8)
          
            write (2219,2000) name, icode_s

2000        format (3x,a8,2x,i8)

         ENDIF

         IF (LCODE.EQ.'TC') THEN
            READ (LINE,213) XCOR,YCOR,ZCOR,AXTYPE,OFFST1,stcode
213         FORMAT (14X,3(F14.5,1X),A4,1x,F10.5,2x,a4)
         ENDIF

         IF (LCODE.EQ.'TB') THEN
            READ (LINE,215) XCOR,YCOR,ZCOR,AXTYPE,OFFST1,stcode,iy1,id1
215         FORMAT (14x,3(F14.5,1X),A4,1x,F10.5,2x,a4,8x,i4,1x,i3)
            if (iye.lt.iy1.or.(iye.le.iy1.and.idoy0.lt.id1)) br=1.d0
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccc  tesmer 16.04.04  cccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
              if (br.eq.1.) then
                print*, ' '
                print*, '(', name, ' before jump)'
              end if
              if (br.eq.0.) then
                print*, ' '
                print*, '(', name, ' after jump)'
              end if
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccc  tesmer 16.04.04  cccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
         ENDIF

         IF (LCODE.EQ.'TV') THEN
            READ (LINE,214) XVEL,YVEL,ZVEL

 214         format (11x,f8.4,1x,f8.4,1x,f8.4)   ! OT 08-Mar-2017

c214         FORMAT (13X,F7.5,2X,F7.5,2X,F7.5)  !  Change by Volker
c214         FORMAT (13X,F6.4,2X,F6.4,2X,F6.4)  !
         ENDIF

         GO TO 999

      ENDIF

      GO TO 1

C     WARNING MESSAGES FROM THE PROGRAM TO THE USER
C  IF CODE EQ 'EP'  JUST RETURN. WARNING MESSAGE WILL BE GIVEN
C  IN CHKCAT

1000  continue                    !   OT

      IF (CODE.EQ.'EP') THEN
         CLOSE (9)
         RETURN
      ENDIF

c
c    If radio source is not in catalogue
c


      IF (LCODE.EQ.'CC') THEN

         print 106, name, code, lcode
106      format (a8, 2x, a2, 2x, a2)

         if (CODE.EQ.'CC') THEN

            print *, 'i_no = ', i_no
            i_no=i_no+1
            name_no(i_no)=NAME
            WRITE (2119,107) iye, idoy0, NAME

         end if
107      format (' Year  ',i4,'  DOY  ',i4,2x,a8)
         go to 999

cccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccc  tesmer 16.04.04  cccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
         IF (CODE.EQ.'TB') GO TO 999
            print*, ' '
            print*, ' '
            print*, '***************************************'
            print*, '***************************************'
            print*, '***************************************'
            print*, ' '
            WRITE (*,108) NAME,i_no+1
            print*, '***************************************'
            print*, '***************************************'
            print*, '***************************************'
            print*, ' '
            print*, ' '
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccc  tesmer 16.04.04  cccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccc

         i_no=i_no+1
         name_no(i_no)=NAME

cccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccc  tesmer 16.04.04  cccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
108      format (' The source is out of catalogue',1x,a8,2x,i4/)
         WRITE (2119,*) NAME, ' ', expid
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccc  tesmer 16.04.04  cccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccc

         go to 999
      END IF

      WRITE (19,*)

      IF (LCODE.NE.'CC') THEN

         WRITE (*,100) NAME,CODE
100      FORMAT (///' ******** WARNING !!! ********'//
     *' ITEM ',A8,' IN CODE ',A2,' HAS NOT BEEN FOUND IN THE CATALOG'/
     */' WRITE DOWN THE NAME OF THAT ITEM AND HIT RETURN WHEN READY')

         ISTAT=1000

         READ (*,299) C0
299      FORMAT (A1)
         GO TO 999
      END IF
c         WRITE (*,100) NAME,CODE
c100      FORMAT (///' ******** WARNING !!! ********'//
c     *' ITEM ',A8,' IN CODE ',A2,' HAS NOT BEEN FOUND IN THE CATALOG'/
c     */' WRITE DOWN THE NAME OF THAT ITEM AND HIT RETURN WHEN READY')
c         ISTAT=1000
c         READ (*,299) C0
c299      FORMAT (A1)
c         GO TO 999

1001  WRITE (*,101) CATNAM
101   FORMAT (///' ******** RUN TIME ERROR READING THE CATALOG FILE'/
     */' IN FILE ... ',A12)
      ISTAT=990
      GO TO 999
1002  WRITE (*,102) CATNAM
102   FORMAT (///' ******** RUN TIME ERROR OPENING THE CATALOG FILE '
     *,'... ',A12//
     *' CHECK IF YOU ENTERED THE RIGHT NAME FOR IT, IT MUST EXIST !!')
      ISTAT=990
999   CLOSE (9)
      RETURN
      END

************************************************************************

      SUBROUTINE CONVER (EXPID,NOBS,IYEAR,name_no,i_no)
C
C     SUBROUTINE CONVER
C
C     THIS SUBROUTINE WILL CONVERT THE ORIGINAL NGS-FORMAT DATA FILE
C  INTO  N  DIFFERENT INPUT FILES IN OCCAM FORMAT
C
C     THE REAL CONVERSION IS DONE IN THE SUBROUTINE GETDAT
C
C REVISION  1990 NOVEMBER 19  (N.Z.L.)
C REVISION  1999 SEPTEMBER 22  (OLEG TITOV)
C REVISION 2000 JULY, 13 (OT) Output into BATCH file is changed
C REVISION  2002  MAY 16  (OT) - Observations of sources out of
C                                       catalogue are deleted
C LAST REVISION  2003 JUNE 08 (OT) - Changes in OCCAM_N.FI and OPENING.FOR (dictio)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INCLUDE 'OCCAM_N.FI'

      CHARACTER CATNAM*12,BASFUL*18
      character*8 FULSOR,FULSTA,EXPID,name_no(400)

      COMMON /BASID1/ NDIM, nbas
      COMMON /STATUS/ ISTAT
      COMMON /CATALG/ CATNAM
      COMMON /BASID0/ BASFUL(nbases),FULSOR(nsources),FULSTA(nstations),
     * ista, isor, iobs
      COMMON /LOCCAT/ CORSOR(nsources,2),CORSTA(nstations,3),
     * VELSTA(nstations,3),EPOCH

C     OPEN A BATCH FILE TO RUN OCCAM

      OPEN (11,FILE='BATCH')

      REWIND 11
      WRITE (11,103) EXPID,CATNAM
103   FORMAT (A8,2X,A12)

C
C     GENERATE THE DATA RECORDS
C     NGS data file is being read here
C
      print *, 'The full amount of observations is', nobs
      CALL GETDAT  (IYEAR,NOBS,name_no,i_no)

C    CHECK THE STATUS CODE

      IF (ISTAT.NE.0) RETURN

      WRITE (11,150) ISTA
150   FORMAT (1x,' THE NUMBER OF STATIONS IS ',i3/)

      DO I=1, nbas
         WRITE (11,151) I,BASFUL(I)
      END DO

151   FORMAT (I2,2X,A18)

      CLOSE (11)

      RETURN
      END



      subroutine gal_acceleration (tt, i, sour, corsor)


      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INCLUDE 'OCCAM_N.FI'

      CHARACTER*8  SOUR
      DIMENSION CORSOR(NSOURCES,2)
      COMMON /MATH/ PI, TWOPI, ROG, ROS, ROZ, TOL



c      OPEN (56, access='append',FILE='anon_acc.dat',STATUS='unknown')

      
      AA = 5.8d0   ! magnitude of the acceleration vector

      rag  = pi * 270.d0/180.d0  !  RA of the acceleration vector
      deg =  pi * (-21.d0)/180.d0  !  Dec of the acceleration vector

      a1 =  AA * dcos(rag) * dcos(deg)
      a2 =  AA * dsin(rag) * dcos(deg)
      a3 =  AA * dsin(deg)

c      write (*,*) '  Secular aberration drift is added '
c      write (*,50) a1, a2, a3

c 50    format (' a1 = ', f6.2, ' a2 = ', f6.2, ' a3 = ', f6.2)


      corr_a = 0.d0
      corr_d = 0.d0


      t0 = 2015.d0



      corr_a = -a1*dsin(corsor(i,1)) + a2*dcos(corsor(i,1))  ! microsec of arc/year

      corr_d = -a1*dsin(corsor(i,2))*dcos(corsor(i,1))       ! microsec of arc/year
     *         -a2*dsin(corsor(i,2))*dsin(corsor(i,1))
     *         +a3*dcos(corsor(i,2))


      corr_a = corr_a*(tt-t0)                    !   microsec of arc
      corr_d = corr_d*(tt-t0)                    !   microsec of arc

c      print *, tt-t0, corr_a, corr_d

c      write (56,7771) tt, corsor(i,1), corsor(i,2), corr_a, corr_d
c 7771 format (f10.4,4(2x,f10.5))

      corr_a = 1.0d-6*corr_a/ros                         !   rad
      corr_d = 1.0d-6*corr_d/ros                         !   rad

c      print *, tt-t0, corr_a, corr_d

c      print *, corsor(i,1), corsor(i,2), corr_a

ccc      CORSOR(i,1) = corsor(i,1) + corr_a
ccc      CORSOR(i,2) = corsor(i,2) + corr_d

c      print *, corsor(i,1), corsor(i,2)

      return
      end


      

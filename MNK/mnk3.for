                                       
C**********************************************************************          
C    
C   COLL IS A PROGRAM DESIGNED TO ANALYSE GEODETIC VLBI DATA
C   USING A LEAST SQUARES COLLOCATION TECHNIQUE. 
C 
C   IT IS INTEGRATED ON THE PACKAGE OCCAM AND USES THE OCCAM STANDARD
C   DATA FILES AS INPUT. IT ALSO NEEDS OTHER FILES WITH THE A PRIORI
C   INFORMATION AND THE OPTIONS TO BE USED 
C
C   Revision  07, April, 1998, (OT)
C   Revision  20, June, 1998, (OT)
C   Revision  14, March, 2001, (OT)
C   Revision  06, April, 2001, (OT)  COMMON /DATE/ deleted 
C   Revision  18, July, 2001, (OT)  Zenith distanse can be read from
C                                        COLLOCAT.OPT
C   Revision  18, July, 2001, (OT)  Maximum of observables is 3500
C                                        instead of 2000
C   Revision  29, November, 2001, (OT)  PSOU calculated as elements
C                  of partial derivatives matrixes for source coordinates
C                  adjustment
C   REVISION March 05, 2002, (OT)   812 radiosources in list,
C                                        3310 observations
C   REVISION May 20, 2002, (OT) Weighting of constraints,
C                  Limits for observations increased to 20000 (IDICTIO, IOBS)
C   REVISION May 21, 2002, (OT)  Clock offset is refered to middle of
C                      the session instead of begin;
C                      Also geodetic nutation is subtracted from result
C                      in case of IAU1980 nutation model implemented
C   REVISION June 11, 2002, (OT)   4000 observations instead of 3310
C   REVISION July 29, 2002, (OT) NNR, NNT approaches for ITRF;
C                                parameter chi for local parameters added
C   REVISION August 26, 2002, (OT)  3300 instead of 4000
C   LAST REVISION April 28, 2004, (OT)  5000 instead of 4000
C   Revision 24, February, 2005, (OT) Clock second derivatives on default;
C                                     Single clock break implementation
C   Revision 17, February, 2006, (OT) 1912 sources are available;
C                 Several lists of global/arc radiosources are optional
C                 Formal errors of the spherical component positions have
C                 been revised
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
C*********************************************************************

C
C   SPECIFICATIONS
C

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'OCCAM.FI'

      CHARACTER*3 CMODEL
      DIMENSION R2(nobs1), Z2(nobs1), H2(nobs1*npar)


 
C
C    COMMON LIST.  DESCRIPTION
C
C    * CINDEX .- LIST OF NUMBERS AND DIMENSIONS
C    * MATH   .- MATHEMATICAL CONSTANTS
C    * SGM    .- SIGMA STORED INFORMATION
C


      COMMON /CINDEX/ NBAS,NOBS,NSTM,NTIM,NSTA,NSEL,
     *               IBAS,ISTM,ITIM,ISTA,IRAND,NOUT
      COMMON /MATH/ PI, TWOPI, ROG, ROS, ROZ, TOL
      COMMON /OPT/ SGADD,ZEN_DIST,SGFAC(nmax),TSTFAC
      COMMON /PHYS/ C, FL, AH, AU, J20

C
C    GET PARAMETERS, OPEN FILES, INITIALIZE
C

c      open (55, file='test_col.dat', status='unknown')
      open (54, file='test_sou.dat', status='unknown')

      CALL PARAM
      CALL KOPEN
      CALL READOPT (CMODEL)
      CALL KINIT  (NR1)

C
C    READ STOCHASTIC PROCESSES INFORMATION
C    REORDER THE CONTENTS OF DICTIO
C

      CALL KTRANS

      CALL PREPARE  (CMODEL, NR1)

 
C    CLOSE FILES

      CALL KCLOSE

C
C    END THE PROGRAM
C
      STOP
      END


C**********************************************************************
C
C    SUBROUTINE KOPEN
C
C    THIS ROUTINE OPEN MOST OF THE NECCESSARY FILES TO RUN COLL
C
C   REVISION MARCH 17, 1998
C
C**********************************************************************

      SUBROUTINE KOPEN

      COMMON /CINDEX/NBAS,NOBS,NSTM,NTIM,NSTA,NSEL,
     *               IBAS,ISTM,ITIM,ISTA,IRAND,NOUT

C
C   OPEN FILES
C

      CALL UNFOLD (17,'BASTIM      ')
      CALL UNFOLD (18,'STATIM      ')
      CALL UNFOLD (19,'SORTIM      ')
      CALL UNFOLD (20,'STACAT      ')
      CALL UNFOLD (21,'DICTIO      ')
      OPEN (23,FILE='collocat.opt')

C
C    READ NUMBER OF RECORDS OF DIRECT ACCESS FILES
C

      read (17,rec=1) NBAS
      read (18,rec=1) NSTM
      READ (19,rec=1) NTIM
      read (20,rec=1) NSTA


      END


C**********************************************************************
C
C    SUBROUTINE KCLOSE
C
C    THIS ROUTINE CLOSES MOST OF THE FILES USED WHEN RUNNING COLL
C
C    REVISION MARCH 17, 1998
C
C**********************************************************************

      SUBROUTINE KCLOSE
      CLOSE (17) 
      CLOSE (18)
      CLOSE (19)
      CLOSE (20)
      CLOSE (21)
      CLOSE (25)
      END

C**********************************************************************
C
C    SUBROUTINE READOPT
C
C    THIS ROUTINE READS A PRIORI INFORMATION AND INITIALIZE THE
C    VARIABLES TO BE USED BY COLL
C
C    REVISION MARCH 17, 1998
C    LAST REVISION MARCH 14, 2001  (OT)
C
C**********************************************************************

      SUBROUTINE READOPT (CMODEL)

C   ESPECIFICATIONS

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'OCCAM.FI'

      CHARACTER CMODEL*3

      COMMON /OPT/ SGADD, ZEN_DIST, SGFAC(nmax), TSTFAC
      COMMON /PHYS/ C, FL, AH, AU, J20


C    READ FILE WITH SETUP DATA


      READ (*,307) CMODEL, ZEN_DIST, SGADD, TSTFAC
      write (*,307) cmodel, zen_dist, sgadd, tstfac
307   FORMAT (A3,4x,f3.0,1x,f6.3,2x,f5.2)

C     READ IN METER, PASS TO Seconds

      SGADD = SGADD/c

      RETURN

      END


C**********************************************************************
C
C    SUBROUTINE KINIT
C
C    THIS ROUTINE READS A PRIORI INFORMATION AND INITIALIZE THE
C    VARIABLES TO BE USED BY COLL
C
C    REVISION MARCH 17, 1998
C    REVISION January 21, 2001, (OT)   Troposphere gradients added
C    LAST REVISION MAY 05, 2002, (OT)  Check IPOL=0
C
C**********************************************************************

      SUBROUTINE KINIT (NR1)

C   ESPECIFICATIONS

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'OCCAM.FI'
 
c
c   Start of update by Oleg Titov (13-Nov-2009, Friday!!)
c   Extra blank field was added
c
      CHARACTER STAT(nmax)*8,STATCAT(nmax)*8,ANTENNA(nmax)*8,
     * ANTE1*8,ANTE2*8
c
c   End of update by Oleg Titov (13-Nov-2009, Friday!!)
c
      CHARACTER STC*4, AX*4
      DIMENSION IANTENNA(nmax),KANTENNA(nmax)
ccc      DIMENSION ECC(nmax,3), IEC(nmax), CONVER(nmax)

C   COMMON LIST. DESCRIPTION OF NEW COMMON AREAS
C
C   * CANTENNA.- INFORMATION ABOUT SELECTED ANTENNAS
C   * PHYS    .- PHYSICAL CONSTANTS

      COMMON /CINDEX/NBAS,NOBS,NSTM,NTIM,NSTA,NSEL,
     *               IBAS,ISTM,ITIM,ISTA,IRAND,NOUT
      COMMON /CANTENNA/ ANTENNA,IANTENNA,KANTENNA,NANTENNA
      COMMON /PHYS/ C, FL, AH, AU, J20
      COMMON /OPT/ SGADD, ZEN_DIST, SGFAC(nmax), TSTFAC
      COMMON /MATH/ PI, TWOPI, ROG, ROS, ROZ, TOL


C   SET VARIABLES TO ZERO

      NVER = 0

C   SET ARRAYS TO ZERO

      DO I=1,nmax
         ANTENNA(I)='        '
         IANTENNA(I)=0
         KANTENNA(I)=0
c         CONVER(I) = 0.D0
      END DO

C
C   GET WHICH ANTENNAS ARE TO BE SOLVED
C     READ INFORMATION FROM FILE 'STATIONS' WITH A LOOP
C

      NANTENNA=0

      READ (23,389,ERR=301,END=301) NANTENNA
389   FORMAT (13X,I2)


      DO IA=1,NANTENNA

c
c   Start of update by Oleg Titov (13-Nov-2009, Friday!!)
c   Extra blank field was added
c
ccc         READ (23,100,ERR=301,END=301) STAT(IA),SGFAC(IA),CVER,
ccc     *    CONVER(IA), ECC(IA,1),ECC(IA,2),ECC(IA,3)

c
c  18-Nov-2013
c 
 

         READ (23,100,ERR=301,END=301) STAT(IA),SGFAC(IA)



c
c   End of update by Oleg Titov (13-Nov-2009, Friday!!)
c

100      FORMAT (A8,2X,F6.3)

ccc100      FORMAT (A8,2X,F6.3,2X,A1,2X,F8.5,3F8.3)


c
c   End of update by Oleg Titov (18-Nov-2013, Friday!!)
c


C       IF ECCENTRICITIES
 
c
c Commented on 02-12-2011
c

c          IEC(IA)=0
c          DO II=1,3
c             IF(ECC(IA,II).NE.0.D0)THEN
c                IEC(IA)=1
c                GO TO 110
c             ENDIF
c          ENDDO
c
c110       CONTINUE

          SGFAC(IA) = SGFAC(IA)/C   !  [meter] -> [sec]


c
c   Start of update by Oleg Titov (13-Nov-2009, Friday!!)
c   Extra blank field was added
c
c          IF (STAT.NE.'        ') THEN

      END DO
c
c   End of update by Oleg Titov (13-Nov-2009, Friday!!)
c



C     CHECK IF THIS ANTENNA IS IN THE CATALOG

c
c   Start of update by Oleg Titov (13-Nov-2009, Friday!!)
c   Extra blank field was added
c
 
      DO I=1,NSTA

         READ (20,REC=I) JST,STATCAT(I),X,Y,Z,EP,VX,VY,VZ,AX,
     *       OFF,STC,npp1
             print 6767, JST, STATCAT(i)
6767         format (i4,2x,a8)

c             write (54,*) STATCAT(I), X, Y, Z, EP, OFF

c         print 34, ax, off, stc, npp1
c 34      format (a4,2x,f10.5,2x,a4,2x,i6)
      END DO

      print *, ' NANTENNA = ', NANTENNA , ' NSTA = ', NSTA

      DO IA=1,NANTENNA

         IFLAG=0

         DO I=1,NSTA

            IF (STATCAT(i).EQ.STAT(ia)) THEN

               ANTENNA(IA)=STAT(ia)
               IANTENNA(IA)=I
cc               IANT_NST(IA)=Npp1
               KANTENNA(I)=IA

c               print *, i, ia

               IFLAG=1

            ENDIF

         END DO

C    WARNING IF ANTENNA NOT IN CATALOG

         IF (IFLAG.EQ.0) then

             WRITE (*,*) ' ***** THE ANTENNA ',STAT,' NOT FOUND'

         END IF
      END DO

c
c   End of update by Oleg Titov (13-Nov-2009, Friday!!)
c

C   WRITE INFORMATION ON FILE KHISTORY

C  WARNING AND BREAK IF LESS THAN 2 ANTENNAS FOUND

      IF (NANTENNA.LT.2) THEN
         WRITE (*,*) NANTENNA, ' ANTENNAS FOUND. NOT ENOUGH TO PROCEED'
         STOP
      ENDIF

C   CHECK IF ANY DIRECTION IS TO BE CONSTRAINED


C
C    COMPUTE THE NUMBER OF PARAMETERS TO BE SOLVED FOR
C    WE PLAN TO SOLVE 9 PARAMETERS AT EACH STATION:
C       - THREE COORDINATES X,Y,Z
C       - CLOCK OFFSET CLOCK RATE and CLOCK SECOND RATE
C       - ATMOSPHERE
C       - TWO ATMOSPHERE GRADIENTS
C    FOR THE FIRST STATION WE WILL SOLVE ONLY THE ATMOSPHERE
C    and STATION COORDINATES (Six parameters)
C

      if (NANTENNA.eq.2) then 
              
          NR = (NANTENNA-1)*(NS-3) + (ns-6)   !  Update on 02-12-2011

c          NR = (NANTENNA-1)*(NS-3)    !  Update on 21-03-2023


      else

          NR = (NANTENNA-1)*NS + (ns-3)   !  Update on 02-12-2011

      end if
   
      print *, ' nantenna = ', nantenna, ' nr = ', nr, ' ns = ', ns

     

      if (NANTENNA.eq.2) then 

         NR1 = NR+2   ! for nutation

      ELSE

         NR1 = NR + 8   ! for nutation,EOP and rates

      END IF

      NR1 = NR1 + 3 !   for Sun's centre motion

      print *, ' nr1 = ', nr1

      RETURN

C   WARNING AND BREAK IS COLLOCAT.OPT IS WRONG OR NOT COMPLETE

301   WRITE (*,*)
     +   ' ***** The file COLLOCAT.OPT has errors. UNABLE TO PROCEED'
      STOP

      END
  
  
C**********************************************************************
C
C    SUBROUTINE prepare (cmodel, nr1)
C
C    THIS ROUTINE READS THE OBSERVABLES, THEIR WEIGHTING MATRIX AND
C    THE JACOBIAN MATRIX
C
C    REVISION MARCH 17, 1998
C    REVISION January 21, 2001, (OT)   Troposphere gradients added
C    REVISION September 14, 2001, (OT)   DELT instead of IDELT
C    REVISION May 05, 2002, (OT)   Output of local and stochastic parameters
C    REVISION July 29, 2002, (OT) NNR, NNT approaches for ITRF
C                                      parameter chi for local parameters added
C    REVISION October 22, 2002, (OT) 642 sources
C    REVISION December 08, 2003, (OT) Output
C    REVISION June 02, 2005, (OT) Statistic check,
C          change for Gilcreek model, n_unstable from the 'OCCAM.FI' file
C    REVISION June 07, 2005, (OT) Output for middle epoch calculation revised
C    REVISION July 21, 2005, (OT) Input for amount of unstable quasars
C                                      from the INPUT8.TXT file
C    REVISION March 22, 2006, (OT) Cable correction implemented;
C                                      One break per session
C    REVISION October 31, 2006, (OT) Baseline results output
C    REVISION November 07, 2008, (OT)  Session identificator
C                                outputed for source.dat
C   LAST REVISION December 02, 2009, (OT) Dry troposphere delay output 
C                              for individual sites
C
C**********************************************************************
    
      SUBROUTINE prepare (cmodel, nr1)
   
C   ESPECIFICATIONS
 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'OCCAM.FI'

      CHARACTER*8 ANTENNA(nmax), STAN, fulsor(nss), arcsour(nss)
      CHARACTER*3 CMODEL
      character*4 axtyp, stc
      CHARACTER*8 SOURCE(nobs1), sour, icrf107(n_unstable)
      CHARACTER*8 st_eli_1(nmax), st_eli_2(nmax), ss(nss)
	
      character*12 antfile, stfile, sphfile
      character*96 line
      character*21 basfile
      character dbc_code*2, mo*3                 !   07.11.2008

      integer*4 ibrkst(10), ibrktm(10), iis(nobs1)            !   24.11.2004
      integer*4 k_br(10), idat_br(10), irow_break(10)       !   22.11.2010
      integer*4 nsig_bas(nmax,nmax)                         !   17.09.2011
      integer*4 num_eli(nmax,2)                             !   17.09.2011
      integer*4 i_st1(nobs1), i_st2(nobs1)                  !   21.11.2013
      integer*2 isor, im, in 

      integer*4 inn(nobs1)                                  !   01.09.2016
      integer*4 indr(nobs1), ng(nss), num(nsou)             !   13.10.2022

   
      real*8 ut_br(10), year_br(10), ts(10)                 !   22.11.2010

 

      DIMENSION r2(nobs1),z2(nobs1),h2(nobs1*npar), numob(n_unstable)
      DIMENSION TIM(nobs1),IDICTIO(nobserv),IOBS(nobserv)
      DIMENSION PCOORD(2,8),PCOORR(2,8),PNUT(2),PPOL(3),PSOU(2),PPOL1(3)
      DIMENSION IANTENNA(nmax),KANTENNA(nmax),KPAR(2),cx(nmax,3)
      DIMENSION number(nss), ra(nss), de(nss)
   
           
      DIMENSION teta(nobs1) !  30 May 2014
      
      
      DIMENSION bb(nobs1)
 
      dimension is1(nobs1), is2(nobs1), year(nobs1),
     *       eps(nobs1), covm(npar*(npar+1)/2), covm1(npar*(npar+1)/2),
     *       corr(npar*(npar+1)/2), eps1(nobs1), eps2(nobs1)
      dimension y(npar), sy(npar)

      dimension p1(nobs1)
      dimension p2(nobs1), dd2(nmax), chik2(nmax)


      dimension wrms2(nmax), eps_res2(nmax), chik(nmax), wrms(nmax)

      dimension nsta_rew(nmax), eps_res(nmax), nk(nmax), dd(nmax)
      dimension dt(nobs1), ra2(n_unstable), de2(n_unstable)
      
     
      dimension cov1(npar)
 
      dimension z31_1(nobs1), z32_1(nobs1)
     
      dimension sig_bas(nantenna,nantenna), sigma_bas(nantenna,nantenna)       ! 17-Sep-2011
      
c      dimension corr_str(nobs1)                                             ! 1-Apr-2015 (OT)
      
    

      real*8 mjd, chisq, chi, s1, sd, e, disp_pol, disp_tr, disp_ut


C  COMMON LIST. DESCRIPTION OF NEW ITEMS
C   * CEXTRACT.- KEEPS INFORMATION ABOUT OBSERVABLES, WEIGHTS AND
C                PARTIAL DERIVATIVES READ FROM OCCAM DATA FILES

      COMMON /CANTENNA/ ANTENNA,IANTENNA,KANTENNA,NANTENNA
      COMMON /CINDEX/ NBAS,NOBS,NSTM,NTIM,NSTA,NSEL,
     *               IBAS,ISTM,ITIM,ISTA,IRAND,NOUT
      COMMON /DICTION/ TIM,IDICTIO,IOBS
      COMMON /CEXTRACT/ OC,SOBD,OCR,SOBR,PCOORD,PCOORR,
     * PNUT,PPOL,PSOU,PPOL1
      COMMON /OPT/ SGADD, ZEN_DIST, SGFAC(nmax), TSTFAC
      COMMON /PHYS/ C, FL, AH, AU, J20
      COMMON /MATH/ PI, TWOPI, ROG, ROS, ROZ, TOL
      COMMON  /G/   UTTJ, XP(30), YP(30), MJD(30), UT1KOR(30), IP

      COMMON /GRAV/ GMS 
                      
      open (1,file='dbc_code',status='unknown')           !  07.11.2008
      read (1,6) iy1, mo, idd, dbc_code                    !  07.11.2008
 6    format (i2.2,a3,i2.2,a2)                                !  07.11.2008
      write (*,6) iy1, mo, idd, dbc_code                    !  07.11.2008
 
      close (1)
c
c  INPUT FILE - THE CATALOGUE of ARC SOURCES
c
      read (*,*) num_unst
      print *, 'num_unst = ', num_unst

      if (num_unst.eq.102) OPEN(24,FILE='arc_102.CAT',status='unknown')   !  102 ICRF 'other' sources
      if (num_unst.eq.163) OPEN(24,FILE='arc_163.CAT',status='unknown')   !   163 unstable Feissel
      if (num_unst.eq.11)  OPEN(24,FILE='11.CAT', status='unknown')   !  for 1922-224
      if (num_unst.eq.13)  OPEN(24,FILE='13.CAT', status='unknown')   !  for estimation of gamma

      if (num_unst.eq.14)  OPEN(24,FILE='14.CAT', status='unknown')   !  14 high redshift
      if (num_unst.eq.1495)OPEN(24,FILE='arc_1495.CAT',status='unknown')  !  all distant and close not 'stable'
      if (num_unst.eq.1591)OPEN(24,FILE='arc_1591.CAT',status='unknown')  !


      OPEN (28,access='append',FILE='EOP.DAT',status='unknown')

 
      if (num_unst.eq.57)  OPEN(24,file='arc_57.cat',status='unknown')    !  for aus2021b
      if (num_unst.eq.87)  OPEN(24,file='arc_87.cat',status='unknown')    !  for aus2024a
      if (num_unst.eq.123)  OPEN(24,file='arc_123.cat',status='unknown')    !  for aus2024a
     
      if (num_unst.eq.8014)OPEN(24,FILE='arc_8014.CAT',status='unknown')  !  non 303  (for aus2021b)

      if (num_unst.eq.5485)                                               !  all sources for 2023
     *       open(24,file='arc_5485.cat',status='unknown')

      if (num_unst.eq.7667)                                               !  non 303  (for 2019)
     *       open(24,file='arc_7667.cat',status='unknown')

 
      if (num_unst.ne.11.and.num_unst.ne.163.and.num_unst.ne.123.
     *  and.num_unst.ne.57.and.num_unst.ne.87.and.num_unst.ne.8014.
     *  and.num_unst.ne.102.and.num_unst.ne.14.and.num_unst.ne.1495.
     *  and.num_unst.ne.1591.and.num_unst.ne.5485.and.num_unst.ne.7667)
     *   print *, 'ERROR in num_unst!! '

c      print *, ' num_unst = ', num_unst 
 
c
c  OUTPUT FILES
c

      open (12, access='append',file = 'WEIGHTS.DAT', status='unknown')
      open (14, access='append',file = 'STATIST.DAT', status='unknown')
 
      open (114, access='append',
     * file = 'STATIST_OF_BAD_SESSIONS.DAT', status='unknown')    !  1-April-2016

      OPEN (22,access='append',FILE='SOURCES.DAT',status='unknown')  !  arc sources output
      open (27, file = 'CHISQ', status='unknown')

      open (11, file = 'residuals.DAT', status='unknown')

      open (333, access='append',file='trop_clock.DAT',status='unknown')
      open (400, file = 'list_of_source', status = 'unknown')
    
      open (477, file='tempor.dat',status='unknown')


      open (88, access='append', file='baseline.dat',status='unknown')
c
      open (332,access='append',file='Sun_pos.dat',status='unknown')
 

      do i=1,nantenna

         antfile(1:8) = antenna(i)
         stfile(1:8) = antenna(i)
         sphfile(1:8) = antenna(i)
         antfile(9:12) = '.EPS'
         stfile(9:12) = '.PAR'
         sphfile(9:12) = '.SPH'

         ist = iantenna(i)
c         print *, i, ist, kantenna(ist)

         OPEN (30+ist,access='append',FILE=antfile,status='unknown')
         OPEN (60+ist,access='append',FILE=stfile,status='unknown')
         OPEN (90+ist,access='append',FILE=sphfile,status='unknown')

      end do

      do i=1, nmax

         nk(i) = 0

      end do

      do i=1,nobs1
 
         p1(i) = 0.d0
         p2(i) = 0.d0

      end do

      print *, ' nsta = ', nsta 


C   INITIALIZE TO ZERO ALL THE ARRAYS INVOLVED

      CALL KZERO (nfull)

c
c  No Clock breaks for reference station!!!
c  to be first in COLLOCAT.OPT
c
c
      CALL CLOCK_BREAKS (IDNM, IBRKST, IBRKTM)

      if (idnm.gt.0) print *, 'The number of breaks is ', IDNM


      if (idnm.gt.0) then

         print *, nr1, 3*idnm
         print *, 'Total number of pars after clock breaks = ',
     *   nr1+3*idnm

         nr1=nr1+3*idnm

         do i=1,idnm

            print *, ' ibrktm = ', ibrktm(i)
            ts(i) = 0.d0  

            k_br(i) = kantenna(ibrkst(i))

      READ (19,REC=IBRKTM(i)) j,source(k_br(i)),
     *                ra2000,de2000,Idat_br(i),UT_br(i)

            year_br(i) = dble(IDAT_BR(i)) + UT_BR(i)/TWOPI
            print 67, i, antenna(k_br(i)), year_br(i)
 67         format (i4,2x,a8,2x,f12.4)

         end do

      end if

      CALL BASELINE_ELIMINATED (IBI, ST_ELI_1, ST_ELI_2)

      do j=1,ibi

         do i=1,nsta

            if (antenna(i).eq.ST_ELI_1(j)) then

               num_eli(j,1) = i

            end if

            if (antenna(i).eq.ST_ELI_2(j)) then

               num_eli(j,2) = i

            end if

         end do

c         print *, ibi, ' num_eli = ', num_eli(j,1), num_eli(j,2)

      end do


      READ (19,REC=1) J,source(1),ra2000,de2000,Idat_1,UT1_1
c      print *, idat_1, ut1_1/twopi
      READ (19,REC=ntim) J,source(ntim),ra2000,de2000,Idat_ntim,UT1_ntim
 

      TIME_1 = dble(IDAT_1) + UT1_1/TWOPI
      TIME_NTIM = dble(IDAT_ntim) + UT1_ntim/TWOPI

      tt=(time_ntim + time_1)/2.d0


 
      isor=0                      

      nfull=nfull+6
      nf=nfull-6


      if (nantenna.eq.2) then

             
         nfull=nfull-6    !  NO NNT, NNR for ITRF

      end if


      print *, ' The number of observations is = ', nf


      ijj = 0

      do i=1,nss

       
         
         ng(i) = 0

         FULSOR(i) = '        '    !  02.09.2016


      end do

      num1 = 0  

   1  continue


         num1 = num1 + 1

         read (400,440,end=445) ss(num1)

 440     format(3x,a8)

         go to 1

 445  continue

      print *, ' number of radio sources is = ', num1 

      rewind 400 

      DO K = 1,num_unst

c         number(k) = 0

         arcsour(1:8) = '        '

         read (24,125) icrf107(k)
 
125      format (a8)
c225      format (3x,a8)
         
c         print 125, icrf107(k)

      END DO

      dda = 0.d0
      dde = 0.d0 

      close (400)


      DO IROW=1,nf 

         IRECORD = iobs(IROW)
        
         indr(irow) = 0

         READ (17,rec=irecord) JBAS,ISOUR,IST1,IST2
         READ (18,REC=IST1) JSTA,ISTA
         READ (18,REC=IST2) JSTA,ISTB
         READ (19,REC=isour) J,source(irow),ra2000,de2000,Idat0,UT0,DELT  ! Sources name

         if (irow.eq.1) TIME_0 = dble(IDAT0) + UT0/TWOPI
         if (irow.eq.nf) TIME_N = dble(IDAT0) + UT0/TWOPI

         TT_irow = IDAT0 + UT0/TWOPI

         t1= tt_irow - tt

         
         CALL ICRF2000 (source(irow), ng, num1,isor,ra2000, de2000,
     *                            fulsor, num, in, ra, de)

         inn(irow) = in  
        
c         print *, ' isour = ', isour, ibrktm(1), ibrktm(2)

         if (idnm.gt.0) then

            do i=1,idnm
               if (isour.eq.ibrktm(i)) then

			    irow_break(i) = irow
c                  print *, 'isor = ', isor, 'irow = ', irow

               end if
            end do

c            print *, idnm, i, irow, irow_break(i)

         end if 

         iis(irow)=isour

         is1(irow)=ista
         is2(irow)=istb

         i_st1(irow) = kantenna(is1(irow))      !        for collocation
         i_st2(irow) = kantenna(is2(irow))      !        for collocation

         IF (KANTENNA(ista).NE.0.AND.KANTENNA(istb).NE.0) THEN

            KPAR(1) = ista
            KPAR(2) = istb

            CALL EXTRACT (irecord, cable1, cable2, source(irow),
c+  OT  04.06.2010
     *       ra2000,de2000,
c-  OT  04.06.2010
     *       antenna(is1(irow)), antenna(is2(irow)),                                        !  pressure on 29-Jan-2013
     *   cmodel, part_gamma,
     *   teta(irow), bb(irow),      ! TETA added on 22 May, 2014; COS_A added on 28 December, 2014
     *   thdelay, freq0, z31, z32, parallax )     !                      ! 26-June-2015, thdelay added on 6-Nov-2015 for software cross-check, 19-Sep-2019

         par = parallax/dsin(teta(irow))
c
c  freq0 - first term for delay rate
c
 
            
 
              z31_1(irow) = z31
              z32_1(irow) = z32
                                            
c              p1(irow) = part_gamma
            

C   FILL THE OBSERVABLES VECTOR AND THE WEIGTHING MATRIX (DIAGONAL)
C     ADD EXTRA NOISE THAT MUST BE READ FROM AN OPTIONS FILE

            Z2(IROW) = OC          !  [sec]

            Z2(irow) = Z2(irow) - cable1*1.d-9      !  Cable correction 22.03.2006
            Z2(irow) = Z2(irow) + cable2*1.d-9      !  Cable correction 22.03.2006


            SGTOT = SOBD**2 + SGADD**2
                                                  

C    FILL THE NUTATION OFFSET PART OF THE JACOBIAN MATRIX

            if (nantenna.eq.2) then                     !  singlebase  !!!!

               H2 (IROW + (NR1-2)*nfull) = PNUT(1)
               H2 (IROW + (NR1-1)*nfull) = PNUT(2)

            else


               H2 (IROW +(NR1-11)*nfull) =-par*dcos(ra2000)*dcos(de2000)
               H2 (IROW +(NR1-10)*nfull) =-par*dsin(ra2000)*dcos(de2000)
               H2 (IROW +(NR1-9)*nfull) = -par*dsin(de2000)   !  Sun's centre shift

               H2 (IROW + (NR1-8)*nfull) = PNUT(1)
               H2 (IROW + (NR1-7)*nfull) = PNUT(2)

               H2 (IROW + (NR1-6)*nfull) = PPOL(1)
               H2 (IROW + (NR1-5)*nfull) = PPOL(2)
               H2 (IROW + (NR1-4)*nfull) = PPOL(3)

 

               H2 (IROW + (NR1-3)*nfull) = PPOL(1)*t1
               H2 (IROW + (NR1-2)*nfull) = PPOL(2)*t1
               H2 (IROW + (NR1-1)*nfull) = PPOL(3)*t1
             

            end if

C    DETERMINE THE COLUMNS FOR ANT1 AND ANT2 FOR THE JACOBIAN MATRIX

            TIMDT = dble(IDAT0) + UT0/TWOPI + (32.184+DELT)/86400.d0
            YEAR(irow) = timdt



            DO K=1,2

               IF (K.EQ.1) ISIG=-1
               IF (K.EQ.2) ISIG=+1

               NKANT = KANTENNA(KPAR(k))


               CALL JULDAT (IY,ID,IH,IM,S,IDAT0,UT0,1)

               IFAC = 365
               IF (MOD(IY,4).EQ.0) IFAC=366
               YEAR_1 = dble(IY)+(dble(ID)+dble(IH)/23.93447d0+
     *          dble(IM)/1440.d0+dble(S)/86400.d0)/(dble(IFAC)+0.2422d0)
               


C    FILL JACOBIAN MATRIX


               IF (NKANT.EQ.1) THEN                  !  Reference station


                 
c                   SGTOT = SGTOT + SGFAC(NKANT)**2

                  IJ=IROW

                  if (nantenna.eq.2) then   !  singlebase

c                     H2(IJ)        =ISIG*PCOORD(K,6)       !  atmosphere
c                     H2(IJ+nfull)  =ISIG*PCOORD(K,7)       !  N-S atmosphere grad.
c                     H2(IJ+2*nfull)=ISIG*PCOORD(K,8)       !  E-W atmosphere grad.

c                     tt1(irow,k) = 1.d9 * ISIG*PCOORD(K,6)       ! partials  [nsec/cm]
c                     tt2(irow,k) = 1.d9 * ISIG*PCOORD(K,7) 
c                     tt3(irow,k) = 1.d9 * ISIG*PCOORD(K,8) 

                  else

                     H2(IJ)        =ISIG*PCOORD(K,1)       !  X   [sec/m]
                     H2(IJ+nfull)  =ISIG*PCOORD(K,2)       !  Y   [sec/m]
                     H2(ij+2*nfull)=ISIG*PCOORD(K,3)       !  Z   [sec/m]


                  end if

               ELSE                                 !  No reference station

c                   SGTOT = SGTOT + SGFAC(NKANT)**2

                  if (nantenna.eq.2)  then

                     IJ = IROW + 3*nfull

                     ts0=(tim(ntim)+tim(1))/2.d0
                     dt(irow) = (tim(isour)-ts0)/43200.d0

                     H2(IJ + 3*nfull)=ISIG/(100.d0*c)               !  [sec/cm]   clock offset
                     H2(ij + 4*nfull)=isig*dt(irow)/(100.d0*c)      !  clock rate
                     H2(ij + 5*nfull)=isig*dt(irow)**2/(100.d0*c)   !  clock 2nd rate


                  else

                     IJ = IROW + (NKANT-2)*ns*nfull + 3 * nfull    !  6 with troposphere, 3 without troposphere

ccc                     IJ = IROW + (NKANT-2)*ns*nfull + 6 * nfull


C-------------------------------------------
C   PCOORD(K,1) = PARTIAL WITH RESPECT TO X
C   PCOORD(K,2) = PARTIAL WITH RESPECT TO Y
C   PCOORD(K,3) = PARTIAL WITH RESPECT TO Z
C-------------------------------------------
                    
                     H2(IJ)        =  ISIG*PCOORD(K,1)          !  X   [sec/m]
                     H2(IJ+nfull)  =  ISIG*PCOORD(K,2)          !  Y   [sec/m]
                     H2(ij+2*nfull)=  ISIG*PCOORD(K,3)          !  Z   [sec/m]

C--------------------------------------------------------------------
C   PARTIALS WITH RESPECT TO CLOCK OFFSET, CLOCK RATE AND ATMOSPHERE
C--------------------------------------------------------------------

                     if (idnm.gt.0) then                      !  Clock break found
                      
                        do i=1,idnm

                           if (nkant.eq.kantenna(ibrkst(i)) ) then


                              if (isour.le.ibrktm(i)) then

                                 ts(i) = (tim(ibrktm(i))+tim(1))/2.d0
 

                                 factor=43200.d0*(year_br(i)-year(1))
                                 dt(irow) = (tim(isour)-ts(i))/factor

      

                       H2(IJ+3*nfull)=ISIG/(100.d0*c)   !  clock offset   [sec/cm]
                       H2(ij+4*nfull)=isig*dt(irow)/(100.d0*c)           !  clock rate
                       H2(ij+5*nfull)=isig*dt(irow)**2/(100.d0*c)        !  clock 2nd rate

c
c     Partials for clock break - 3 (or more) columns before EOP (nr1-9, nr1-10, nr1-11, etc)
c


         	                 H2 (IROW + (NR1-(3*i+11))*nfull) = 0.d0      !  clock offset         
                                 H2 (IROW + (NR1-(3*i+10))*nfull) = 0.d0      !  clock rate
                                 H2 (IROW + (NR1-(3*i+9))*nfull) = 0.d0      !  clock 2nd rat

                             else if (isour.gt.ibrktm(i)) then

                               ts(i)=(tim(ibrktm(i))+tim(ntim))/2.d0
         
                               factor=43200.d0*(1.d0+year(1)-year_br(i))
                               dt(irow) = (tim(isour)-ts(i))/factor

                                 H2(IJ+3*nfull)=0.d0                   !  clock offset
                                 H2(ij+4*nfull)=0.d0                   !  clock rate
                                 H2(ij+5*nfull)=0.d0                   !  clock 2nd rate


 

c
c     Partials for clock break - 3 (or more) columns before EOP (nr1-9, nr1-10, nr1-11, etc)
c

           H2 (IROW + (NR1-(3*i+11))*nfull) = ISIG/(100.d0*c)             !  [sec/cm]clock offset
           H2 (IROW + (NR1-(3*i+10))*nfull) = isig*dt(irow)/(100.d0*c)    !  clock rate
           H2 (IROW + (NR1-(3*i+9))*nfull) = isig*dt(irow)**2/(100.d0*c) !  clock 2nd rate

                              end if

                           else                    !  dealing with problem station

                              ts0=(tim(ntim)+tim(1))/2.d0
                              dt(irow) = (tim(isour)-ts0)/43200.d0

                 H2(IJ+3*nfull)=ISIG/(100.d0*c)              !  [sec/cm] clock offset
                 H2(ij+4*nfull)=isig*dt(irow)/(100.d0*c)     !  clock rate
                 H2(ij+5*nfull)=isig*dt(irow)**2/(100.d0*c)  !  clock 2nd rate

                           end if

                        end do

                     else                                     ! No clock breaks

                        ts0=(tim(ntim)+tim(1))/2.d0
                        dt(irow) = (tim(isour)-ts0)/43200.d0

                        H2(IJ+3*nfull)=ISIG/(100.d0*c)               !  [sec/cm]   clock offset
                        H2(ij+4*nfull)=isig*dt(irow)/(100.d0*c)      !  clock rate
                        H2(ij+5*nfull)=isig*dt(irow)**2/(100.d0*c)   !  clock 2nd rate


                     end if


                  END IF  !   More than 2 stations

               END IF   !  No  reference stations

            END DO

   

            ns1 = irow + 2*(inn(irow)-1)*nfull                  !  01.09.2016
            ns2 = irow + (2*inn(irow)-1)*nfull                 !  01.09.2016
    
c            if (bb(irow).gt.3.d6) then           
           
               p1(irow) = psou(1)   !   [rad]
               p2(irow) = psou(2)   !   [rad]

c               number(inn(irow)) = number(inn(irow)) + 1

      
c            else

c              p1(irow) = 0.d0
c              p2(irow) = 0.d0

c            end if
         
            R2(IROW) = sgtot  !  [sec**2]

            
            if (temp1.eq.999.d0.or.temp2.eq.999.d0) then 

		     print *, ' temp1 =', temp1, ' temp2 = ',temp2
			 R2(IROW) = r2(irow)*1.d9

            end if

            

            if (mo.eq.'JUN') imo = 06
            if (mo.eq.'JUL') imo = 07
            iy2 = 2000 +iy1


            if (dabs(s-60.d0).le.0.01) then

                im = im+1
                s=0.d0

            end if
      
            
         

c      write (76,333) irow, iy2, imo, idd, ih, im, s, source(irow),
c     *  antenna(i_st1(irow)), antenna(i_st2(irow)), thdelay

c        if (irow.lt.10) print *, thdelay, 2.d0*part_gamma


c 333  format (i6,1x,i4,'/',i2.2,'/',i2.2,1x,i2.2,':',i2.2,':','0',
c     * f4.2,1x,a8, 1x,a8,1x,a8,1x,e21.14)



         ENDIF

         do k=1,num_unst

            if (fulsor(inn(irow)).eq.icrf107(k)) then
                
               indr(irow) = 1
                  
            end if

         end do


      ENDDO
                   	
 
      tt=(time_0+time_n)/2.d0

      

      if (idnm.ne.0) then
         do i=1,idnm
            print *, i, 'irow_break = ', irow_break(i), '  ts = ', ts(i)
         end do
      end if

      write (12,86) tt
 86   format (//'  NEW SESSION ', f12.4/)

c
c   Partial derivatives for arc sources as daily parameters
c


      do i = 1,isor

         number(i) = ng(num(i))   ! amount of observations sorted on source

         write (477,*) i, num(i), ng(num(i))

      end do

      imarc = 0


c      print *, ' ok 1', ' nf = ', nf, ' num_unst = ', num_unst


c      write (477,*)

      limit_on_observations = 4
      print *, ' limit of observations per source per experiment = ',
     *  limit_on_observations

      DO IROW = 1,nf

         if (indr(irow).eq.1) then
        
            if (number(inn(irow)).lt.limit_on_observations) then

                indr(irow) = 0

            end if

         end if

      END DO

                           
      narc = 0

      
      DO IROW = 1, nf

         if (indr(irow).eq.1) then

            islg = 0

            do kk=1,num_unst


               if (fulsor(inn(irow)).eq.arcsour(kk)) then

                 

                  islg=1
                  km=kk
                  exit

               end if
            end do

            if (islg.eq.0) then

               narc = narc + 1
               imarc = imarc + 1 
               arcsour(narc) = fulsor(inn(irow))

c               print *, irow, narc, arcsour(narc)

               
               km = imarc

               numob(narc) = number(inn(irow))

               ra2(narc)  = ra(inn(irow))
               de2(narc)  = de(inn(irow))

            end if

         end if

     
cc         if (narc.eq.0) print *, ' irow = ', irow
cc         if (narc.ne.0) print *, ' irow = ', irow, ' narc = ', narc, fulsor(inn(irow))

c
c   ARC sources as daily parameters
c 

         if (narc.gt.0) then

            ns_arc1 = irow + nr1*nfull + 2*(km-1) * nfull   !  RA for arc sources
            ns_arc2 = irow + nr1*nfull + (2*km-1) * nfull   !  DE for arc sources

            if (indr(irow).eq.1) then

c               if (bb(irow).le.6.d6) then

                  h2(ns_arc1) = p1(irow)  !  [radian]
                  h2(ns_arc2) = p2(irow)  !  [radian] 

c               else 

c                  h2(ns_arc1)=0.d0
c                  h2(ns_arc2)=0.d0              

c               end if      
                             
            else

               h2(ns_arc1)=0.d0
               h2(ns_arc2)=0.d0

            end if
         end if

      ENDDO




      do i=1,6

          r2(nf+i)=1.d-5      !  covariance for contrained obs

c          r2(nf+i) = 1.d0

      end do


c
c +OT    18-Sep-2011
c

c
c   Downweighting of eliminated baselines
c

      do i=1,nf

         
         do j=1,ibi


           if ((i_st1(i).eq.num_eli(j,1).and.i_st2(i).eq.num_eli(j,2))
     *     .or.(i_st1(i).eq.num_eli(j,2).and.i_st2(i).eq.num_eli(j,1)))
     *     then

              print *, j, num_eli(j,1), num_eli(j,2)

              r2(i) = r2(i)*1.d6

           end if

         end do


      end do




c
c -OT    18-Sep-2011
c

      if (nsta.ne.2) then   ! if not single-baseline
      
         do i=1,nsta

            READ (20,REC=i) JST,STAN,RX,RY,RZ,EPOCH,VX,VY,VZ
c            print 796, stan, vx, vy, vz
            

            ik=kantenna(i)
            if (ik.ne.0) then

               rx=rx+vx*(year_1-epoch)
               ry=ry+vy*(year_1-epoch)
               rz=rz+vz*(year_1-epoch)

               rr = rx*rx + ry*ry + rz*rz

               if (ik.eq.1) then


          if (stan.eq.'METSAHOV'.or.stan.eq.'DSS65A  '.or.stan.eq.
     *   'TIGOCONC'.or.stan.eq.'KUNMING '.or
     * .stan.eq.'PARKES  '.or.stan.eq.'ISHIOKA '.or
c     * or.stan.eq.'ONSA13SW'.or.stan.eq.'ONSA13NE'
c     * .or.stan.eq.'MACGO12M'.or.stan.eq.'RAEGSMAR'.or
     * .stan.eq.'GILCREEK'.or.stan.eq.'KOKEE12M'.or
     * .stan.eq.'KOGANEI '.or.stan.eq.'KASHIM11'.or
     * .stan.eq.'DSS36   '.or.stan.eq.'TSUKUB32'.or
     * .stan.eq.'KASHIM34'.or.stan.eq.'URUMQI  ')     then


c                  if (stan.eq.'GILCREEK'.and.year_1.ge.2002.84d0) then

                     H2(nf+1)=0.d0
                     H2(nf+2+nfull)=0.d0
                     H2(nf+3+2*nfull)=0.d0


                     H2(nf+4+nfull)=  0.d0
                     H2(nf+4+2*nfull)=  0.d0

                     H2(nf+5)=  0.d0
                     H2(nf+5+2*nfull)=  0.d0

                     H2(nf+6)=  0.d0
                     H2(nf+6+nfull)=  0.d0

                  else

                     H2(nf+1) = 1.d0
                     H2(nf+2+nfull) = 1.d0
                     H2(nf+3+2*nfull) = 1.d0


                     H2(nf+4+nfull)=  -ah*rz/rr
                     H2(nf+4+2*nfull)=  ah*ry/rr

                     H2(nf+5)=  ah*rz/rr
                     H2(nf+5+2*nfull)=  -ah*rx/rr

                     H2(nf+6)=  -ah*ry/rr
                     H2(nf+6+nfull)=  ah*rx/rr

                  end if
               else

c
c   No troposphere estimated
c
                  IJ1=nf+1 + (ik-2)*ns*nfull + 3*nfull
                  IJ2=nf+2 + (ik-2)*ns*nfull + 3*nfull
                  IJ3=nf+3 + (ik-2)*ns*nfull + 3*nfull
                  IJ4=nf+4 + (ik-2)*ns*nfull + 3*nfull
                  IJ5=nf+5 + (ik-2)*ns*nfull + 3*nfull
                  IJ6=nf+6 + (ik-2)*ns*nfull + 3*nfull

c
c  - Troposphere estimated
c


c                  IJ1=nf+1 + (ik-2)*ns*nfull + 6*nfull
c                  IJ2=nf+2 + (ik-2)*ns*nfull + 6*nfull
c                  IJ3=nf+3 + (ik-2)*ns*nfull + 6*nfull
c                  IJ4=nf+4 + (ik-2)*ns*nfull + 6*nfull
c                  IJ5=nf+5 + (ik-2)*ns*nfull + 6*nfull
c                  IJ6=nf+6 + (ik-2)*ns*nfull + 6*nfull



          if (stan.eq.'METSAHOV'.or.stan.eq.'DSS65A  '.or.stan.eq.
     *   'TIGOCONC'.or.stan.eq.'KUNMING '.or
     * .stan.eq.'PARKES  '.or.stan.eq.'ISHIOKA ' .or
c     * or.stan.eq.'ONSA13SW'.or.stan.eq.'ONSA13NE'
c     * .or.stan.eq.'MACGO12M'.or.stan.eq.'RAEGSMAR'.or
     * .stan.eq.'GILCREEK'.or.stan.eq.'KOKEE12M'.or
     * .stan.eq.'KOGANEI '.or.stan.eq.'KASHIM11'.or
     * .stan.eq.'DSS36   '.or.stan.eq.'TSUKUB32'.or
     * .stan.eq.'KASHIM34'.or.stan.eq.'URUMQI  ')     then


c                  if (stan.eq.'GILCREEK'.and.year_1.ge.2002.84d0) then

                     H2(ij1)=0.d0
                     H2(ij2+nfull)=0.d0
                     H2(ij3+2*nfull)=0.d0


                     H2(ij4+nfull)=  0.d0
                     H2(ij4+2*nfull)=  0.d0

                     H2(ij5)=  0.d0
                     H2(ij5+2*nfull)=  0.d0

                     H2(ij6)=  0.d0 
                     H2(ij6+nfull)=  0.d0

                  else

                     H2(IJ1)=1.d0
                     H2(IJ2+nfull)=1.d0
                     H2(IJ3+2*nfull)=1.d0


                     H2(IJ4+nfull)=  -ah*rz/rr
                     H2(IJ4+2*nfull)=  ah*ry/rr

                     H2(IJ5)=  ah*rz/rr
                     H2(IJ5+2*nfull)=  -ah*rx/rr

                     H2(IJ6)=  -ah*ry/rr
                     H2(IJ6+nfull)=  ah*rx/rr
                   end if
                end if
             end if

         end do
      end if

      nr2 = nr1 + 2*narc   !  Increasing amount of daily parameters due to arc sources

      print *, ' nr1 = ', nr1, ' nr2 = ', nr2, nfull


      do i=1,nantenna
         
         wrms(i) = 0.d0
         chik(i) = 0.d0
         eps_res(i) = 0.d0

      end do

      schi = 0.d0
      sw = 0.d0
      s00 = 0.d0
      s11 = 0.d0


      do i=1, nfull

c         write (477,777) (h2(nfull*(j-1)+i),j=1,nr2)

          write (477,778) c*z2(i), 
     * antenna(i_st1(i)), antenna(i_st2(i)), source(i)

      end do
c 777    format (200(1x,f18.12))
 778   format (f16.8,3(2x,a8))


      call mnk1 (nr1, nr2, nfull, year, r2, z2, h2, 
     *                         y, schi, sw, eps1, s00, s11, covm1)

c      print *, ' first mnk is finished '

       print *, ' s00 = ', s00, '  s11 = ', s11


   
      n_i=0
c
c     Downweighting of outliers after LSM
c
      do i=1, nf

         eps_res(i_st1(i)) = eps_res(i_st1(i)) + (eps1(i)**2)/r2(i)
         eps_res(i_st2(i)) = eps_res(i_st2(i)) + (eps1(i)**2)/r2(i)

         dd(i_st1(i)) = dd(i_st1(i))+1.d0/r2(i)
         dd(i_st2(i)) = dd(i_st2(i))+1.d0/r2(i)

         nk(i_st1(i)) = nk(i_st1(i))+1
         nk(i_st2(i)) = nk(i_st2(i))+1

      end do

      do i=1,nantenna

         nsta_rew(i) = 0
         wrms(i) = dsqrt(eps_res(i)/dd(i))                !  sec
         chik(i) = dsqrt(eps_res(i)/dble(nk(i)))

      end do



c      print *, ' schi = ', schi, ' sw = ', sw
c      print *, ' ibi  = ', ibi


      do i=1, nf


         if ((eps1(i)**2/r2(i)).gt.tstfac * schi) then

    
            r2(i) = r2(i)*100000.d0

            n_i=n_i+1

            nsta_rew(i_st1(i)) = nsta_rew(i_st1(i))+1
            nsta_rew(i_st2(i)) = nsta_rew(i_st2(i))+1

         end if

      end do

      print *, n_i, ' observations have been downweighted  after LSM '


c      do i=1,nantenna
c
c         write (12,88) antenna(i), nk(i), chik(i), 
c     *  100.d0 * c *wrms(i), nsta_rew(i)                        !  cm

c         print 88, antenna(i), nk(i), chik(i), 
c     *  100.d0*c*wrms(i),nsta_rew(i)
c
c      end do

c 88   format (1x,a8,4x,i4,2(3x,f8.4),3x,i4,' observations downweighted')


      write (*,444) isor, narc, isor-narc
444   format (' The number of sources in the session is ',i3,';'2x,i4,'
     * - Arc',2x,i4,' - Global')

      print *, ' schi = ', schi, ' sw = ', 1.d12 * sw, ' pks ',
     * ' sw = ', 100.d0 * c * sw, '  cm '
 
      schi_limit = 20
      sw_limit = 20

      print *, ' schi_limit = ', schi_limit,'  sw_limit =  ', sw_limit

      if (schi.gt.schi_limit.or.sw.gt.sw_limit) then

         print *, ' BAD SESSION '
         write (114,*) tt, 'schi = ', schi, 'sw = ', sw

          do i=1, nantenna

             write (114,87) antenna(i), nk(i), chik(i), wrms(i)
             print 87, antenna(i), nk(i), chik(i), wrms(i)

          end do
 
 87   format (1x,a8,4x,i4,2(3x,f8.4))

      else


         chi   = 0.d0
         s1 = 0.0d0
         sd = 0.0d0 
         e = 0.d0
         s0 = 0.0d0
        
         call mnk1 (nr1, nr2, nfull, year, r2, z2, h2, 
     *                         y, chi1, sd, eps2, s0, s1, covm) 


         s1 = s1 / ((100.d0 * c )**2) !  s1 converted from [1/s^2] to [1/cm^2]
       
         print  *, ' s0 = ', s0, '  s1 = ', s1

         print *, ' second LSM iteration has finished '


      do i=1,nf

          if ((i_st1(i).eq.2.and.i_st2(i).eq.3).or.
     *      (i_st1(i).eq.3.and.i_st2(i).eq.2))  then

           write (88,589) year(i), c * eps2(i), c * z2(i), iis(i),
     *    antenna(i_st1(i)), antenna(i_st2(i)), source(i)
 

          end if


  
       end do

  589  format (f12.6,2(2x,f17.4),2x,i6,3(2x,a8))


         do i=1,nantenna
         
            wrms2(i) = 0.d0
            chik2(i) = 0.d0
            eps_res2(i) = 0.d0

         end do                           

         do i=1, nf

            eps_res2(i_st1(i)) = eps_res2(i_st1(i)) + (eps2(i)**2)/r2(i)
            eps_res2(i_st2(i)) = eps_res2(i_st2(i)) + (eps2(i)**2)/r2(i)

            dd2(i_st1(i)) = dd2(i_st1(i)) + 1.d0/r2(i)
            dd2(i_st2(i)) = dd2(i_st2(i)) + 1.d0/r2(i)

         end do

         do i=1,nantenna

            wrms2(i) = dsqrt(eps_res2(i)/dd2(i))                !  sec
            chik2(i) = dsqrt(eps_res2(i)/dble(nk(i)))

         end do

         print *, ' second iteration output '       

         print *, ' residual written in baseline.dat file '

         do i=1,nantenna
         
            write (12,189) antenna(i), nk(i), chik2(i), 
     *         100.d0 * c *wrms2(i), nsta_rew(i)                    !  cm
            print 189, antenna(i), nk(i), chik2(i), 
     *        100.d0*c*wrms2(i), nsta_rew(i)

         end do

189   format (1x,a8,4x,i4,2(3x,f8.4),2x,i4,' observations downweighted')

c    eps2 in [sec]

         sd1 = 3.d10 * sd  !  [sec] -> [cm]
        
         print *, ' chi = ', chi1, ' weighted mean root square = ', sd1
c
c  17-Sep-2011
c
         do j=1, nantenna

	    do k=1, j-1

                sigma_bas(k,j) = 0.d0
	        sig_bas(k,j) = 0.d0
	        nsig_bas(k,j) = 0

	    end do

         end do
    
c
c   17-Sep-2011
c
   

         read (27,102) chisq_sum, s22, nf_sum, nr1_sum

         print *, chisq_sum, s22, nf_sum, nr1_sum
         print *, ' s0 = ', s0, '  s1 = ', s1
         print *, chi1, sd1, nf, nr2 

         chisq_sum = chisq_sum + s0
         s22 = s22 + s1 
     
         nf_sum = nf_sum + nf
         nr1_sum = nr1_sum + nr2

              do i=1, nr2

            ii=i*(i+1)/2

            do j=1,i
               ij=i*(i-1)/2+j
               jj=j*(j+1)/2
               corr(ij)=covm(ij)/dsqrt(covm(ii)*covm(jj))

            end do

         end do
                               
         rewind 27

         write (27,102) chisq_sum, s22, nf_sum, nr1_sum

 102      format (1x,f16.2,1x,f12.2,1x,i8,2x,i6)
 
         write (14,103) tt, nf, chi1, sd1, 1.d12 * sd, n_i, 
     * iy1, mo, idd, dbc_code

 103  format (f11.4,1x,i5,3(2x,f8.4),2x,i5,2x,i2.2,a3,i2.2,a2)
  
       
c
c     Calculation of the EOP
c
         call eop           !          Reading of EPHEM.DAT
         time_m = (year(1)+year(nf))/2.d0

         CALL LAGINT1 (  MJD, UT1KOR, IP,   time_m, UT10, utm1, utp1)
         CALL LAGINT1 (  MJD,     XP, IP,   time_m, X0, xm1, xp1 )
         CALL LAGINT1 (  MJD,     YP, IP,   time_m, Y0, ym1, yp1 )


c
c        Daily parameters
c
c         print *, ' nr1 before writing to 77', nr1

          do i=1,nr1

            sy(i) = dsqrt(covm(i*(i+1)/2))       
            
c            print *, i, y(i)
 
         end do

         
        
c
c        Daily parameters (EOP)       mas -> arcsec
c
         do i = nr1-7,nr1


            y(i)=y(i)/1.d3
c
c   Correction back for geodetic nutation in case of IAU1980 nutation model
c

            if (nutmod.eq.1) then

               call geodnut (tt, dpsig)
               y(nr1-7)=y(nr1-7)-dpsig

            end if

            sy(i) = chi1*sy(i)/1.d3   !   mas -> arcsec
 
         end do
c
c        Output for local parameters (station coordinates and EOP)
c

         corr_dp_de = corr((NR1-6)*(NR1-5)/2-1)

         corr_x_y   = corr((NR1-3)*(NR1-4)/2-1)

         corr_x_ut1 = corr((NR1-3)*(NR1-2)/2-2)

         corr_y_ut1 = corr((NR1-3)*(NR1-2)/2-1)

         xa1 = (xp1-xm1)                        !  Apriori for X-dot
         ya1 = (yp1-ym1)                        !  Apriori for Y-dot
         ua1 = -(utp1-utm1)                     !  Apriori for LOD


         
        write (332,2200) tt, 
     *   (1.496d11*y(nr1-11+i), 1.496d11*sy(nr1-11+i), i=1,3)


2200     format (f12.6,6(2x,f10.5))


         print *, ' Dx  ',  1.496d11 * y(nr1-10), 1.496d11*sy(nr1-10)
         print *, ' Dy  ',  1.496d11 * y(nr1-9), 1.496d11*sy(nr1-9)
         print *, ' Dz  ',  1.496d11 * y(nr1-8), 1.496d11*sy(nr1-8)
        


         nhour = 24

          write(28,105) tt,y(nr1-5), y(nr1-4),
     *  y(nr1-3), y(nr1-7), y(nr1-6), sy(nr1-5), sy(nr1-4),
     *  sy(nr1-3), sy(nr1-7), sy(nr1-6), sd1, nfull-6,
     *  nhour, xa1 + y(nr1-2), ya1+y(nr1-1), ua1+y(nr1),
     *  sy(nr1-2), sy(nr1-1), sy(nr1), iy1, mo, idd, dbc_code


c        write(277,1105) tt, c*y(nr1-2)/ros,c*y(nr1-1)/ros, c*y(nr1)/ros,
c     *  c*sy(nr1-2)/ros, c*sy(nr1-1)/ros, c*sy(nr1)/ros
     

 105  format (1x,f12.5,2(2(1x,f11.6),1x,f12.7,2(1x,f11.6)),2x,f10.5,
     * 2x,i6,2x,i4,6(1x,f11.6),2x,i2.2,a3,i2.2,a2)

c1105  format (f13.5,6(2x,f11.6))

c         write (77,*) tt


         do i=1,nantenna

            sr=0.d0
            sphi=0.d0
            slam=0.d0

            ist=iantenna(i)

            READ (20,REC=IST) JST,STAN,RX,RY,RZ,EPOCH,VX,VY,VZ

            CALL TRANSF (RX,RY,RZ,PHI,GLONG,ELHGT,R0,0)


            rx=(rx+vx*(year_1-epoch))      ! [m] 
            ry=(ry+vy*(year_1-epoch))      ! [m] 
            rz=(rz+vz*(year_1-epoch))      ! [m] 

            IY = EPOCH
            ID = (EPOCH-IY)*365.25d0
            IH = 0
            IM = 0
            S  = 0.d0

C      ADD VELOCITY MODEL

            CALL JULDAT (IY,ID,IH,IM,S,IDJEP,UTJEP,0)

            d_earthquake=2002.842d0    !   Earthquake epoch
            if (year_1.ge.d_earthquake) then

               TG = YEAR_1-D_EARTHQUAKE

c               if (STAN.EQ.'GILCREEK') write (77,998) stan, rx, ry, rz
c               if (STAN.EQ.'TSUKUB32') write (77,998) stan, rx, ry, rz

               if (STAN.EQ.'GILCREEK') then

                  tau=0.25d0

                  du = 0.8d0
                  dn = -4.78d0 - 2.37d0*(1.d0-dexp(-tg/tau))
                  dea = 2.26d0 + 0.86d0*(1.d0-dexp(-tg/tau))

c      if (STAN.EQ.'GILCREEK') write (77,998) stan,du,dn,dea,tg,tg/tau

                  offx= -dsin(phi)*dcos(glong)*dn - dsin(glong)*dea +
     *  dcos(phi)*dcos(glong)*du
                  offy= -dsin(phi)*dsin(glong)*dn  +dcos(glong)*dea +
     *  dcos(phi)*dsin(glong)*du
                  offz=  dcos(phi)*dn + dsin(phi)*du


c                 if (STAN.EQ.'GILCREEK') print *, rx, offx
c                  if (STAN.EQ.'GILCREEK') print *, ry, offy
c                  if (STAN.EQ.'GILCREEK') print *, rz, offz

                  RX=RX+offx
                  RY=RY+offy
                  RZ=RZ+offz


               end if

            end if
c-  OT

c            print *, rx,ry,rz


            if (i.eq.1) then

               rx_f = rx + y(1)
               ry_f = ry + y(2)
               rz_f = rz + y(3)


               s1x  = sy(1)
               s1y  = sy(2)
               s1z  = sy(3)


            else

               rx_f = rx + y(ns*(i-1)-2)
               ry_f = ry + y(ns*(i-1)-1)
               rz_f = rz + y(ns*(i-1))


cc               print *, y(ns*(i-1)-2)
cc               print *, y(ns*(i-1)-1)
cc               print *, y(ns*(i-1))


               s1x  = sy(ns*(i-1)-2)
               s1y  = sy(ns*(i-1)-1)
               s1z  = sy(ns*(i-1))



            end if


            cx(i,1) = rx_f               !   [m]
            cx(i,2) = ry_f               !   [m]
            cx(i,3) = rz_f               !   [m]

            CALL ESPHER (I,COVM,CX,R,PHI,DLAM,SR,SPHI,SLAM)

c
c  For output to SINEX file -- 18-Dec-2013
c
           
            s1x = chi1*s1x           !    m
            s1y = chi1*s1y           !    m
            s1z = chi1*s1z           !    m

      
            
            sr   = chi1*sr
            sphi = chi1*sphi
            slam = chi1*slam

c
c   Format changed from 198 to 197  (noted by N.Pavlovskaya on 19-Jan-2010)
c
c            write (60+ist,198) tt,year_1,rx_f,ry_f,rz_f,


            if (i.eq.1) then

                write (60+ist,197) tt,year_1,rx_f,ry_f,rz_f,
     *           s1x,s1y,s1z,nk(i), rx, ry, rz, y(1), y(2), y(3), ist

            else

                write (60+ist,197) tt,year_1,rx_f,ry_f,rz_f,
     *           s1x,s1y,s1z,nk(i), rx, ry, rz, 
     *           y(ns*(i-1)-2), y(ns*(i-1)-1), y(ns*(i-1)), ist

            end if


            if (i.eq.1) write (90+ist,198) tt,year_1,r,phi,dlam,
     *          sr,sphi,slam,nk(i)

            if (i.ne.1) write (90+ist,198) tt,year_1,r,phi,dlam,
     *          sr,sphi,slam, nk(i)




         end do

 98   format (1x,f15.2,1x,f11.2)
 97   format (1x,f11.6,1x,f11.6)
 96   format (1x,f15.7,1x,f11.7)
 197  format (1x,f12.5,1x,f10.5,2(6(1x,f17.5),2x,i6))   !  19-Jan-2010
 198  format (1x,f12.5,1x,f10.5,6(1x,f17.5),2x,i6)   !  16-Dec-2013
 
 299  format (a8,1x,3(1x,f15.4),3(f7.4))

 897  format (a8,3(2x,f14.2))
 898  format (8x,6(2x,f14.2))
 899  format (8x,3(2x,f14.2))
 998  format (a8,6(2x,f14.2))


c
c     WRITING TO FILES FOR BASELINES    31.10.2006 (OT)
c
         do i=2,nantenna
        
            do k=1,i-1

               dx = cx(i,1) - cx(k,1)                                  ! cm
               dy = cx(i,2) - cx(k,2)                                  ! cm
               dz = cx(i,3) - cx(k,3)                                  ! cm

               do jj=1, nr2
                  cov1(jj)=0.d0
               end do

               ipp = (NS-3)+(i-2)*ns       ! Change on 28.11.2005

               COV1(ipp+1) = dx
               COV1(ipp+2) = dy
               COV1(ipp+3) = dz

               IF (k.EQ.1) THEN
                  COV1(1) = -dx
                  COV1(2) = -dy
                  COV1(3) = -dz
               ELSE
                  kpp = (NS-3)+(k-2)*ns    ! Change on 28.11.2005
                  COV1(kpp+1) = -dx
                  COV1(kpp+2) = -dy
                  COV1(kpp+3) = -dz
               ENDIF

               CALL ABAT(COV1,1,nr2,covm,1,0,SBL,NRX,IERb)

c               base(i,k) = dsqrt( dx*dx + dy*dy + dz*dz)               ! m

c               print *, dsqrt(sbl)*chi1, base(i,k)                     ! mm
c               if (ierb.ne.0) print *, i, k, ierb

c               sbase(i,k) = dsqrt(sbl)*chi1/base(i,k)                  ! mcm

               
               ik  = i*(i-1)/2+k



            end do

            
         end do

       

         do i=1,nf
 
            

            write (30+is1(i),99) year(i), eps2(i)/dsqrt(r2(i)), 
     *  1.d2 * c * eps2(i),
     *  source(i), antenna(i_st1(i)), antenna(i_st2(i)), bb(i)
    

            write (30+is2(i),99) year(i), eps2(i)/dsqrt(r2(i)),
     *   1.d2 * c * eps2(i), 
     *   source(i), antenna(i_st1(i)), antenna(i_st2(i)), bb(i)
   
       

c
c       Statistic on individual baselines, added on 17.09.2011
c       Hotel Hilton, in Zurich Airport (we stuck with Johannes, when the
c       flight to Vienna was canceled)
c


c

            do j=1, nantenna

	         do k=1, j-1

                  if ((i_st1(i).eq.j.and.i_st2(i).eq.k).or.
     *                (i_st1(i).eq.k.and.i_st2(i).eq.j))  then

                      if (r2(i).ne.0) then

      sigma_bas(k,j) = sigma_bas(k,j) + eps2(i)**2/r2(i)
      nsig_bas(k,j) = nsig_bas(k,j) + 1

	sig_bas(k,j) = sig_bas(k,j) + 1.d0/r2(i)

	                end if


                  end if

	         end do

            end do

         end do


         

         do j=1, nantenna

            do k=1, j-1

c	         print *, sigma_bas(k,j), sig_bas(k,j)

               if (sig_bas(k,j).ne.0) then

			      sigma_bas(k,j) = dsqrt(sigma_bas(k,j)/sig_bas(k,j))

c	              print 445, antenna(k), antenna(j),
c     *              nsig_bas(k,j), sigma_bas(k,j), k,j
c 445                format (a8,2x,a8,2x,i5,2x,f10.5,4x,i4,2x,i4)

             
               end if

            end do



         end do


 99    format (f12.6,2x,f7.3,2x,f14.4,3(2x,a8),2x,f11.2)

399   format (4(2x,f12.6),2x,f13.7)



c
c  Output of the arc radio source daily estimates
c
         do i=nr1+1,nr2

            y(i) = y(i)*ros                          !  rad -> arcsec

            sy(i) = chi1 * dsqrt(covm(i*(i+1)/2))*ros  !  rad -> arcsec
     

         end do

         print *, 'nr1 = ', nr1, 'narc = ', narc, ' nr2 = ', nr2

         do i=1,narc

            i1 = (nr1+2*i-1)*(nr1+2*i)/2
            i2 = (nr1+2*i)*(nr1+2*i+1)/2

            corr1 = covm(i2-1)/dsqrt(covm(i1)*covm(i2))

            write (22,126) year_1, 
     * ra2(i), de2(i), arcsour(i), 
     * y(nr1+2*i-1)/15.d0,     sy(nr1+2*i-1)/15.d0,  !  sec of time
    
     *      y(nr1+2*i), sy(nr1+2*i), 
     *      dsqrt (y(nr1+2*i-1)**2 + y(nr1+2*i)**2),            !  deflection angle alpha               22 December 2014
     *      dsqrt (sy(nr1+2*i-1)**2 + sy(nr1+2*i)**2),          !  sigma of deflection angle alpha      22 December 2014
     *      corr1, 
     *      
     *      numob(i), ijj, time_m,
     *      iy1, mo, idd, dbc_code                   !  07.11.2008

         end do

      end if  


     

126   format (f12.5,2(2x,f14.8), 2x,a8,6(2x,f18.9),2x,f9.5,                !   TETA added on 22 May, 2014
     *2(2x,i4),2x,f11.5,1x,i2.2,a3,i2.2,a2)                       !  07.11.2008


c909   RETURN

999   WRITE (*,*)  '***** ERROR FOUND ON COLLOCAT.OPT *****'
      WRITE (*,*)  ' CONSTRAINED DIRECTION NOT OBSERVED'
      WRITE (*,*)  '  --- HIT ANY KEY TO CONTINUE --- '

      close (15)
      close (26)
      close (27)


      
      END


C**********************************************************************
C
C    SUBROUTINE KZERO
C
C    THIS ROUTINE INITIALIZES THE VARIABLES USED BY PREPARE
C
C   REVISION MARCH 17, 1998
C   REVISION May 20, 2002, (OT) Limits for observations increased to
C                                    20000 (IDICTIO, IOBS)
C   LAST REVISION July 07, 2003, (OT) Limits for observations increased to
C                                    32000 (IDICTIO, IOBS)
C
C**********************************************************************

      SUBROUTINE KZERO (nfull)

C   ESPECIFICATIONS

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'OCCAM.FI'

      CHARACTER ANTENNA(nmax)*8
      DIMENSION TIM(nobs1),IDICTIO(nobserv),IOBS(nobserv)
      DIMENSION IANTENNA(nmax),KANTENNA(nmax)

      COMMON /CANTENNA/ ANTENNA,IANTENNA,KANTENNA,NANTENNA

      COMMON /CINDEX/NBAS,NOBS,NSTM,NTIM,NSTA,NSEL,
     *               IBAS,ISTM,ITIM,ISTA,IRAND,NOUT
      COMMON /DICTION/ TIM,IDICTIO,IOBS

      NOBS=0
      DO ITIM=1,NTIM
         DO I=1,nobserv

            IF (IDICTIO(I).EQ.ITIM) THEN
               READ (17,rec=I) JBAS,ISOUR,IST1,IST2
               READ (18,REC=IST1) JSTA,ISTA
               READ (18,REC=IST2) JSTA,ISTB
               IF (KANTENNA(ISTA).NE.0.AND.KANTENNA(ISTB).NE.0) THEN
                  NOBS=NOBS+1
                  IOBS(NOBS) = i
               END IF
            ENDIF

         END DO
      END DO
      nfull=nobs

      END

C**********************************************************************
C
C    SUBROUTINE KTRANS
C
C    THIS ROUTINE READS OCCAM STANDARD DATA FILE DICTIO AND GENERATES
C      THE ARRAY IDICTIO.
C    THIS ARRAY CONTAIN FOR EACH TIME WHICH BASELINES HAVE BEEN OBSERVED
C      AT THAT TIME AND THE RECORD NUMBER OF SUCH OBSERVATIONS ON FILE
C      BASTIM
C    THE ROUTINE ALSO FILLS THE TIME ARRAY
C
C   REVISION MARCH 17, 1998
C   Revision January,21, 2001, (OT)   Filling of IDICTIO
C   LAST REVISION May 20, 2002, (OT) Limits for observations increased to
C                                    20000 (IDICTIO, IOBS)
C
C**********************************************************************

      SUBROUTINE KTRANS

C    ESPECIFICATIONS

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'OCCAM.FI'

      CHARACTER*8 SOURCE, ANTENNA(nmax)*8
      DIMENSION IANTENNA(nmax),KANTENNA(nmax)
      DIMENSION TIM(nobs1),IDICTIO(nobserv),NBASE(ntimmax),IOBS(nobserv)

      COMMON /DICTION/ TIM,IDICTIO,IOBS
      COMMON /CINDEX/NBAS,NOBS,NSTM,NTIM,NSTA,NSEL,
     *               IBAS,ISTM,ITIM,ISTA,IRAND,NOUT
      COMMON /MATH/ PI, TWOPI, ROG, ROS, ROZ, TOL

C   PUT ZEROS IN IDICTIO

      DO I=1,nobserv
         IDICTIO(I)=0
      END DO

C   LOOK FOR THE NUMBER OF BASELINES IN DICTIO

      READ (21,REC=1) NSEL

C   READ INITIAL TIME

      READ (19,REC=1) j,source,ra2000, de2000, Idat0,UT0


C   FILL THE TIME ARRAY WITH NUMBER OF SECONDS FROM FIRST OBSERVATION

      do i=1,ntim
         READ (19,REC=I) JP,source,ra2000, de2000, Idatj,UTJ
         TIM(I) = (IDATJ-IDAT0)*86400.d0 + (UTJ-UT0)/TWOPI*86400.D0
      end do


C   READ INFORMATION IN DICTIO AND FILL ARRAY IDICTIO

      DO I=1,NSEL
         READ (21,REC=I) K2,KST1,KST2,NBASE
         IF (KST1.NE.KST2) THEN      !       21.01.2001 (O.Titov)

C  FILL IDICTIO
C  Massive IDICTIO is filling up by observational indexes
C

            DO J=1,ntimmax    !  NTIMMAX in DTAU1.FOR   !!
               IF(NBASE(J).NE.0) THEN
                  READ (17,REC=NBASE(J)) JBAS,ISOUR
                  IDICTIO (NBASE(J)) = ISOUR
               END IF
            END DO
         END IF                      !       21.01.2001 (O.Titov)
      ENDDO

      RETURN

      END   


C**********************************************************************
C
C    SUBROUTINE EXTRACT  (IRECORD, CABL1, CABL2, CMODEL,
C     Z31, Z32, azimut1, azimut2, zdry1, zdry2, part_gamma)
C
C    THIS ROUTINE READS INFORMATION OF OBSERVABLES, SIGMAS AND
C      PARTIAL DERIVATIVES FROM OCCAM STANDARD DATA FILES
C
C   REVISION MARCH 17, 1998
C   REVISION MAY 04, 2002     3310 observations
C   REVISION MAY 13, 2002  (OT)   Calculation of EOP partial
C                  derivatives for delay rates has been added
C   REVISION MAY 04, 2002     4000 observations
C   REVISION March 22, 2006, (OT) Cable correction implemented
C   REVISION 2006 SEPTEMBER 18 (O.T.) Sign for thermal correction was
C   changed due to corresponding change in 'STATION.FOR' and part_gamma added
C   LAST REVISION December 02, 2009, (OT) Dry troposphere delay output
C                              for individual sites
C
C
C**********************************************************************
 
      SUBROUTINE EXTRACT (IRECORD, CABL1, CABL2, sour,
     * ra, de,
     * ant1, ant2, 
     *   CMODEL,  part_gamma,
     * teta, bb, thdelay, 
     * freq0, z31, z32, parallax)      !   19-Sep-2019
                                              

C   ESPECIFICATIONS

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'OCCAM.FI'

      character*3 cmodel
      character*8 sour, ant1, ant2
      DIMENSION PCOORD(2,8),PCOORR(2,8),AXKT1(2),AXKT2(2)
      DIMENSION CRX1(3),CRY1(3),CRZ1(3),CGLONG1(3),CPHI1(3),CELHGT1(3)
      DIMENSION CRX2(3),CRY2(3),CRZ2(3),CGLONG2(3),CPHI2(3),CELHGT2(3)
      DIMENSION cmniel1_w(2), cmniel1_d(2), a_grad1(2)
      DIMENSION cmniel2_w(2), cmniel2_d(2), a_grad2(2)
      DIMENSION PEOP1(3,3),PEOP2(3,3),PPOL(3),PPOL1(3)
      DIMENSION PNUT1(2),PNUT2(2),PNUT(2),PSOU(2),drda(3),drdd(3),v(3)
      dimension acc_der(3), ww1(3), ww2(3), p(3), b(3), bsw(3)

      dimension ppolar_rate1(3), ppolar_rate2(3)

      dimension azimut1(3), azimut2(3)
      
      COMMON /CINDEX/NBAS,NOBS,NSTM,NTIM,NSTA,NSEL,
     *               IBAS,ISTM,ITIM,ISTA,IRAND,NOUT
      COMMON /CEXTRACT/ OC,SOBD,OCR,SOBR,PCOORD,PCOORR,
     * PNUT,PPOL,PSOU,PPOL1
      COMMON /PHYS/ C, FL, AH, AU, J20
      COMMON /MATH/ PI, TWOPI, ROG, ROS, ROZ, TOL
      COMMON /OPT/ SGADD,ZEN_DIST,SGFAC(nmax),TSTFAC


c      open (323, file = 'occam.dat',status = 'unknown')


C   READ ON BASTIM THE RECORD "IRECORD", READ FROM IDICTIO
C   IF CONTAINS THE THEORETICAL AND OBSERVED DELAY AND THE
C   SIGMA OF THE LATTER

      READ (17,rec=IRECORD) JBAS,ISOUR,IST1,IST2,D1,SD1,R01,SR1,
     *           DI1,SDI1,RI1,SRI1,TP1,T0,TM1,
     *           TGravP1,TGrav0,TGravM1,part_gamma,teta,      !   TETA added on 22 May, 2014
     *           cos_A, bb, phi_A          !   cos_a  added on 28 Dec, 2014
     *           , AA, freq0, diff_a, RSun  !   25-Jan-2017
     *           , v(1), v(2), v(3), parallax


C   READ THE STATION INFORMATION AT THE TIME FOR BOTH ANTENNAS
C   THIS INCLUDES THE PARTIAL DERIVATIVES, THE ATMOSPHERIC MODELS
C   AND THE AXIS OFFSETS

      READ (18,REC=IST1) JSTA,ISTA,ISOUR,TEMP1,PRES1,HUMD,IFAC,CABL1,
     *                 CRX1,CRY1,CRZ1,RX1,RY1,RZ1,
     *                 CGLONG1,CPHI1,CELHGT1,Z1,HA1,
     *                 HSTA1,Z2,HA2,HSTA2,Z31,HA3,HSTA31,
     *                 ZDRY1,cmniel1_w,cmniel1_d,
     *                 w11, w12, w13,              !  geodetic velocity of first station
     *                 a11, a12, a13,              !  geodetic acceleration
     *                 AXKT1,
     *                 (PCOORD(1,I),I=1,5),(PCOORR(1,I),I=1,5),
     *                 PNUT1,PEOP1,
     *                 ppolar_rate1,            !  03-Jan-2020
     *            therm_d1,azimut1,a_grad1,
     *            dr1_ut1,                      !  01-Jan-2021
     *            z200_1,smfw3_1,azim1,elev1,
     *            dimfh1,dimfw1,   !  dry and wet IMF maping functions
     *            gimfh1,          !  dry gradient IMF mapping function
     *            grn1,gre1,       !  apriori north and east gradients
     *            iavail_1,
     *            drmfh1,drmfw1,    !  dry and wet ray traced mapping function
     *            iavailr_1,
     *            dqmfh1,dqmfw1,    !  dry and wet Vienna mapping function
     *            iavailq_1
     *            ,gmfh1,gmfw1              !  Global Mapping Function GMF
     *            ,grnh1,greh1,grnw1,grew1  ! apriori gradients
     *            ,iavailv1_1               ! VMF1
     *            ,vmf1h1,vmf1w1
     

      READ (18,REC=IST2) JSTA,ISTB,ISOUR,TEMP2,PRES2,HUMD,IFAC,CABL2,
     *                 CRX2,CRY2,CRZ2,RX2,RY2,RZ2,
     *                 CGLONG2,CPHI2,CELHGT2,Z1,HA1,
     *                 HSTA1,Z2,HA2,HSTA2,Z32,HA3,HSTA32,
     *                 ZDRY2,cmniel2_w,cmniel2_d,
     *                 w21, w22, w23,             !  geodetic velocity of second station
     *                 a21, a22, a23,             !  geodetic acceleration
     *                 AXKT2,
     *                 (PCOORD(2,I),I=1,5),(PCOORR(2,I),I=1,5),
     *                 PNUT2,PEOP2,
     *                 ppolar_rate2,              !  03-Jan-2020
     *                 therm_d2,azimut2,a_grad2,
     *            dr2_ut1,                        !  01-Jan-2021
     *            z200_2,smfw3_2,azim2,elev2,
     *            dimfh2,dimfw2,   !  dry and wet IMF maping functions
     *            gimfh2,          !  dry gradient IMF mapping function
     *            grn2,gre2,       !  apriori north and east gradients
     *            iavail_2,
     *            drmfh2,drmfw2,    !  dry and wet ray traced mapping function
     *            iavailr_2,
     *            dqmfh2,dqmfw2,    !  dry and wet Vienna mapping function
     *            iavailq_2
     *            ,gmfh2,gmfw2              !  Global Mapping Function GMF
     *            ,grnh2,greh2,grnw2,grew2  ! apriori gradients
     *            ,iavailv1_2               ! VMF1
     *            ,vmf1h2,vmf1w2
     

      if (temp1.eq.999.d0.or.temp2.eq.999.d0) then
         print *, 333, temp1, temp2
 333    format (f10.5,2x, f10.5)
      end if

c
c   Added on 04.06.2010
c

      DRDA(1) = -DSIN(ra)*DCOS(de)   !  ds/da
      DRDA(2) =  DCOS(ra)*DCOS(de)   !
      DRDA(3) =  0.D0                      !

      DRDD(1) = -DCOS(ra)*DSIN(de)   !  ds/dd
      DRDD(2) = -DSIN(ra)*DSIN(de)   !
      DRDD(3) =  DCOS(de)               !

      p(1) = dcos(ra)*dcos(de)
	p(2) = dsin(ra)*dcos(de)
	p(3) = dsin(de)



      dva=0.d0
	dvd=0.d0

      denom = 0.d0


      ww2(1) = w21
      ww2(2) = w22
      ww2(3) = w23

      ww1(1) = w11
      ww1(2) = w12
      ww1(3) = w13  
     

    

      freq0 = 0.d0

      b1 = 0.d0
      b2 = 0.d0

      bb2s = 0.d0

      b(1) = crx2(3) - crx1(3)
      b(2) = cry2(3) - cry1(3)
      b(3) = crz2(3) - crz1(3)

c      bsw(1) = (dotpr(b,p)*p(1) - b(1)) * 100.d0 ! cm
c      bsw(2) = (dotpr(b,p)*p(2) - b(2)) * 100.d0 ! cm
c      bsw(3) = (dotpr(b,p)*p(3) - b(3)) * 100.d0 ! cm

       bsw(1) = -100.d0 * dotpr(b,p)*v(1)/c   !  cm
       bsw(2) = -100.d0 * dotpr(b,p)*v(2)/c
       bsw(3) = -100.d0 * dotpr(b,p)*v(3)/c   


      DO I=1,3

c         print *, i, v(i), ww2(i)

         DVA = DVA + (V(i)+ww2(i))*Drda(I)
         DVD = DVD + (V(i)+ww2(i))*drdd(i)

         denom = denom + p(i)*(v(i)+ww2(i))/c

         b1 = b1 + b(i)*p(i)           

         b2 = b2 + ww2(i) * p(i)

         freq0 = freq0 +  (ww2(i)-ww1(i))*p(i)/c

      ENDDO

      denom = 1.d0 + denom

      bb2s = -b1*freq0/c


c
c     MULTIPLICATION ON SEC ( zenit distance )
c
c           BY O.TITOV
c
      z31_1=z31*180.d0/pi
      z32_1=z32*180.d0/pi

      if (z31_1.gt.zen_dist.or.z32_1.gt.zen_dist) then

         if (z32_1.le.z31_1) then
              sd1=sd1/dcos(z31)
         else
              sd1=sd1/dcos(z32)
         end if

      end if

C
C   CONVERT TO PROPER UNITS (SECONDS)
C

      OBD     = D1 *1.D-06
      SOBD    = SD1*1.D-09
      OBR     = R01*1.D-12   !  sec/sec
      SOBR    = SR1*1.D-12   !  sec/sec
      

C   GET TROPOSPHERIC MODEL. IT MUST BE SELECTED VIA AN OPTIONS FILE
C   GET PARTIAL DERIVATIVES FOR THE TROPOSPHERIC PARAMETER


c+ jboehm
c  if no ECMWF data is available

      if (iavail_1.eq.0) then
        dimfh1 = cmniel1_d(1)
        gimfh1 = cmniel1_d(1)
        dimfw1 = cmniel1_w(1)
      end if

      if (iavail_2.eq.0) then
        dimfh2 = cmniel2_d(1)
        gimfh2 = cmniel2_d(1)
        dimfw2 = cmniel2_w(1)
      end if

c  if no RMF data is available

      if (iavailr_1.eq.0) then
        drmfh1 = cmniel1_d(1)
        drmfw1 = cmniel1_w(1)
      end if

      if (iavailr_2.eq.0) then
        drmfh2 = cmniel2_d(1)
        drmfw2 = cmniel2_w(1)
      end if

c  if no VMF data is available

      if (iavailq_1.eq.0) then
        dqmfh1 = cmniel1_d(1)
        dqmfw1 = cmniel1_w(1)
      end if

      if (iavailq_2.eq.0) then
        dqmfh2 = cmniel2_d(1)
        dqmfw2 = cmniel2_w(1)
      end if

c  if no VMF1 data is available

      if (iavailv1_1.eq.0) then
        vmf1h1 = cmniel1_d(1)
        vmf1w1 = cmniel1_w(1)
      end if

      if (iavailv1_2.eq.0) then
        vmf1h2 = cmniel2_d(1)
        vmf1w2 = cmniel2_w(1)
      end if

c
c  ZDRY - in meters
c  TROPO - seconds
c 


      

      IF (CMODEL.EQ.'NIE') THEN
         TROPO = (ZDRY2*cmniel2_d(1) - ZDRY1*cmniel1_d(1))/C            ! sec
         TROPOR= (ZDRY2*cmniel2_d(2) - ZDRY1*cmniel1_d(2))/C      ! sec
         PCOORD(1,6) = cmniel1_w(1)
         PCOORD(2,6) = cmniel2_w(1)
      ENDIF

      IF (CMODEL.EQ.'VMF') THEN
         TROPO = (ZDRY2*dqmfh2-ZDRY1*dqmfh1)/C                        ! sec
         TROPOR= (ZDRY2*cmniel2_d(2)-ZDRY1*cmniel1_d(2))/C      ! sec
         PCOORD(1,6) = dqmfw1
         PCOORD(2,6) = dqmfw2
      ENDIF

      IF (CMODEL.EQ.'VM1') THEN
         TROPO = (ZDRY2*vmf1h2-ZDRY1*vmf1h1)/C                        ! sec
         TROPOR= (ZDRY2*cmniel2_d(2)-ZDRY1*cmniel1_d(2))/C      ! sec
         PCOORD(1,6) = vmf1w1
         PCOORD(2,6) = vmf1w2
      ENDIF


     
c
c  downweighting of low elevation observations, 15 Oct 2015
c

      EXS1 = 0.d0
      EXS2 = 0.d0

c
c     ADDITION of ERROR  ( zenit distance )
c
c           BY O.TITOV
c
      
cc      if (z31_1.gt.75.d0) then           !  >75 degrees downweighted 
cc 
cc         EXS1 = 10.d-12 * PCOORD(1,6)    !  mapping function * 10 pks formal error for 1st station; 01 Oct 2015
cc         SOBD = DSQRT (SOBD**2 + EXS1**2)
cc    
cc      end if
cc
cc      if (z32_1.gt.75.d0) then           ! >75 degrees downweighted 
cc
cc         EXS2 = 10.d-12 * PCOORD(2,6)    !  mapping function * 10 pks formal error for 2nd station; 01 Oct 2015
cc         SOBD = DSQRT (SOBD**2 + EXS2**2)
cc  
cc      end if


c 
c  For Station positions from [sec] to [meter]
c


      PCOORD(1,1) = PCOORD(1,1)/c         !  sec/m
      PCOORD(1,2) = PCOORD(1,2)/c
      PCOORD(1,3) = PCOORD(1,3)/c
      PCOORD(2,1) = PCOORD(2,1)/c
      PCOORD(2,2) = PCOORD(2,2)/c
      PCOORD(2,3) = PCOORD(2,3)/c

c
c  Conversion of units for troposphere derivatives
c  from '1' to sec/cm for troposphere
c

      PCOORD(1,6) = pcoord(1,6)/(100.d0*c)      !  sec/cm
      PCOORD(2,6) = pcoord(2,6)/(100.d0*c)      !  sec/cm

      PCOORD(1,7) = a_grad1(1)*pcoord(1,6)      !  North-South   1st station
      PCOORD(2,7) = a_grad2(1)*pcoord(2,6)      !  North-South   2nd station
      PCOORD(1,8) = a_grad1(2)*pcoord(1,6)      !  East-West     1st station
      PCOORD(2,8) = a_grad2(2)*pcoord(2,6)      !  East-West     2nd station

c      write (55,555) ista, sour, ant1, vmf1w1, pcoord(1,7),  pcoord(1,8)
c      write (55,555) istb, sour, ant2, vmf1w2, pcoord(2,7),  pcoord(2,8)
c 555  format (i4,2(2x,a8),3(2x,f8.4))


C   GET AXIS OFFSET CORRECTIONS

      AXIS  = AXKT2(1)-AXKT1(1)
      AXISR = AXKT2(2)-AXKT1(2)  !        nsec/sec

C   GET THERMIC DEFORMATION CORRECTIONS

      therm = therm_d1-therm_d2    !  Sign was changed on 18.09.2006 (OT)

C   ADD FALSE CLOCK MODEL TO THE DATA FOR TESTING

C   COMPUTE THEORETICAL DELAY ADDING THE GEOMETRICAL MODEL, THE
C  TROPOSPHERIC CONTRIBUTION AND THE AXIS OFFSET CORRECTION


      THDELAY = T0 + TROPO + AXIS + therm    !  to be checked for all types
ccc                                             !  of antennas
 

      THRATE = (TP1-TM1)/2.D0 + TROPOR + AXISR/1.d9  ! sec/sec  23 Sep 2015

      FRAD = 180.D0*3600.D0/PI                   !   206265 arcsec

C   GET OBSERVABLE MINUS CALCULUS

      OC = OBD - THDELAY                    !  sec

      OCR= OBR - THRATE
     
   


C   GET PARTIAL DERIVATIVES WITH RESPECT TO NUTATION

c      PNUT (1) =  (PNUT2(1)-PNUT1(1))/C
c      PNUT (2) =  (PNUT2(2)-PNUT1(2))/C

      PNUT (1) =  (PNUT2(1)-PNUT1(1))/(1.d3*c*FRAD)   ! [m] -> [sec/mas]
      PNUT (2) =  (PNUT2(2)-PNUT1(2))/(1.d3*c*FRAD)   ! [m] -> [sec/mas]


C   GET PARTIAL DERIVATIVES WITH RESPECT TO POLE MOTION AND UT1



      PPOL(1) =  (PEOP2(1,3)-PEOP1(1,3))/(1.d3 *c*FRAD)           !  m -> sec/mas
      PPOL(2) =  (PEOP2(2,3)-PEOP1(2,3))/(1.d3 *c*FRAD)           !  m -> sec/mas
      PPOL(3) =  15.d0 * (PEOP2(3,3)-PEOP1(3,3))/(1.d3*c*FRAD)   !  m -> sec/ms




C   GET PARTIAL DERIVATIVES on delay rates  WITH RESPECT TO EOP   13.05.2002


      PX1  = (PEOP1(1,2)-PEOP1(1,1))/2        !     X - comp  for 1st station
      PY1  = (PEOP1(2,2)-PEOP1(2,1))/2        !     Y - comp  for 1st station
      PUT1  = (PEOP1(3,2)-PEOP1(3,1))/2       !     UT1-UTC   for 1st station

      PX2  = (PEOP2(1,2)-PEOP2(1,1))/2        !     X - comp  for 2st station
      PY2  = (PEOP2(2,2)-PEOP2(2,1))/2        !     Y - comp  for 2st station
      PUT2  = (PEOP2(3,2)-PEOP2(3,1))/2       !     UT1-UTC   for 2st station

      PPOL1(1) =  (PX2-PX1)/(1.d3*c*FRAD)      !  sec/mas
      PPOL1(2) =  (PY2-PY1)/(1.d3*c*FRAD)      !  sec/mas
      PPOL1(3) =  15.d0 * (PUT2-PUT1)/(1.d3*c*FRAD)    !  sec/msec

C   GET PARTIAL DERIVATIVES WITH RESPECT TO SOURCE COORDINATES


c      PSOU (1) =  (PCOORD(2,4)-PCOORD(1,4))*100.d0     !  cm
c      PSOU (2) =  (PCOORD(2,5)-PCOORD(1,5))*100.d0     !  cm

      PSOU (1) =  (PCOORD(2,4)-PCOORD(1,4))/(c*denom)     !  m -> second
      PSOU (2) =  (PCOORD(2,5)-PCOORD(1,5))/(c*denom)     !  m -> second


      END

************************************************************************

      SUBROUTINE eop

C*************************************************************
C
C   SUBROUTINE eop
C
C     THIS SUBROUTINE READS THE TABULAR VALUES FOR THE EARTH
C   ORIENTATION PARAMETERS C   FROM THE FILE "EPHEM.DAT"
C
C     IF EOP ARE OF 1-DAY SERIES, IT REMOVES THE SHORT PERIOD
C   VARIATIONS CORRECTION FROM TABULAR DATA OF UT1-UTC
C
C     IF "EPHEM.DAT" DOES NOT EXIST IT FLAGS AN ERROR MESSAGE
C
C    REVISION 1991 SEPTEMBER 30 (N.Z.L.)
C    REVISION 1993 APRIL 22  READING NEW EPH. FORMAT (N.Z.L.)
C    REVISION 1993 NOV. 17  INCREASE DIMENSIONS (N.Z.L.)
C    REVISION 1997 APRIL 23; READING NEW EOP1 FORMAT
C    Last revision October, 4, 1999, (OT)
C*************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION MJD
      CHARACTER*250 LINE
      CHARACTER*5 BODY

      COMMON  /G/   UTTJ, XP(30), YP(30), MJD(30), UT1KOR(30), IP
      COMMON  /NUT/  de(30), dp(30)
      COMMON /MATH/ PI, TWOPI, ROG, ROS, ROZ, TOL

      OPEN (16,FILE='EPHEM.DAT',STATUS='OLD')
      REWIND 16
      IP  = 0

1     CONTINUE

      READ (16,201,END=2) LINE
      READ (LINE,'(A5)') BODY
201   FORMAT (A250)

C----------------------------------------------------------------
C  CHECK IF THE LINE CORRESPONDS TO EARTH ROTATION PARAMETERS
C  ONE-DAY SERIES WITH SHORT PERIOD UT1 VARIATIONS REMOVED (EOP1)
C  AND READ JULIAN DATE, X WOBBLE, Y WOBBLE AND UT1-UTC
C----------------------------------------------------------------

      IF (BODY.EQ.'EOP1 ') THEN
         IP = IP+1
         READ (LINE,211)
     *        MJD(IP),XP(IP),YP(IP),UT1KOR(IP),idelt,dp(ip),de(ip)
c211      FORMAT (6X,F7.1,1x,2(F8.6),1x,F10.7,i4,f8.5,f8.5)

c
c   Change of EPHEM.DAT reading format - 04.01.2008
c


211      FORMAT (6X,F7.1,1x,2(F8.6),1x,F10.7,i4,f9.6,f9.6)

C----------------------------------------------------------------
C MJD IS SUPPOSED TO BE IN TDT, WHICH IS NOT VERY CLEAR. PROBABLY
C IT IS REALLY IN UTC, BUT IT MUST BE CHECKED. CALC 7.4 ASSUMES IT
C TO BE CT AND THEN INTERPOLATE ALSO IN CT. THE REAL OCCAM VERSION
C SHOUL TAKE IT AS UTC AND INTERPOLATE AS UTC, AT LEAST UNTIL THIS
C POINT IS CROSS CHECKED PROPERLY.
C FOR THE TEST, OCCAM WILL FOLLOW CALC STRATEGY FOR A WHILE.
C----------------------------------------------------------------

      ENDIF

      GO TO 1

2     CONTINUE

      CLOSE (16)
      RETURN

      WRITE (*,*) ' ***** FILE EPHEM.DAT NOT FOUND. HIT CTRL+C TO ABORT'
      READ (*,'(I1)') KK

      RETURN
      END

C************************************************************************
C
C     SUBROUTINE geodnut (mjd, dpsig)
C
C     The subroutine calculate the geodetic nutation correction
C    in accordance with Fukushima (1991)
C
C    21 May 2002 (OT)
C    Last revision, 12 December 2005 (OT) - change the sign due to transition
C                                           to the mean equinox system
C
C*************************************************************
      subroutine geodnut (mjd, dpsig)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION LPRIME, MJD

      COMMON /PHYS/ C, FL, AH, AU, J20
      COMMON /MATH/ PI, TWOPI, ROG, ROS, ROZ, TOL

      T = (MJD - J20)/36525.0D0
      LPRIME = -0.00001149d0*T**4 - 0.000136d0*T**3 - 0.5532d0*T**2
     *+ 129596581.0481d0*T + 1287104.79305d0
      LPRIME = DMOD(LPRIME,1296000d0)
      LPRIME = DMOD(PI*LPRIME/648000.d0,TWOPI)
      DPSIG =  0.000153d0*DSIN(lprime) + 0.000002D0*DSIN(2.d0*lprime)

      return
      end

C************************************************************************
C
      SUBROUTINE LAGINT1 (X,Y,N,XINT,YOUT,ym1,yp1)
C
C*************************************************************
C
C ------- RECOMMENDED BY IERS !!!
C
C     THIS SUBROUTINE PERFORMS LAGRANGIAN INTERPOLATION
C     WITHIN A SET OF (X,Y) PAIRS TO GIVE THE Y
C     VALUE CORRESPONDING TO XINT.  THIS PROGRAM USES A
C     WINDOW OF 4 DATA POINTS TO PERFORM THE INTERPOLATION.
C     IF THE WINDOW SIZE NEEDS TO BE CHANGED, THIS CAN BE
C     DONE BY CHANGING THE INDICES IN THE DO LOOPS FOR
C     VARIABLES M AND J.
C
C     PARAMETERS ARE :
C     X     - ARRAY OF VALUES OF THE INDEPENDENT VARIABLE
C     Y     - ARRAY OF FUNCTION VALUES CORRESPONDING TO X
C     N     - NUMBER OF POINTS
C     XINT  - THE X-VALUE FOR WHICH ESTIMATE OF Y IS DESIRED
C     YOUT  - THE Y VALUE RETURNED TO CALLER
C
C     Distributed by Daniel Gambis (1997)
C
C
C*************************************************************
      REAL*8 X(N),Y(N),XINT,YOUT,TERM, YM1, YP1
      INTEGER N,I,J
C
      YOUT = 0.0D0
      DO 5 I = 1,N-1
        IF ( XINT .GE. X(I) .AND. XINT .LT. X(I+1) ) then
            K = I
            ym1=y(i)
            yp1=y(i+1)
        END IF
    5 CONTINUE
      IF ( K .LT. 2 ) K = 2
      IF ( K .GT. N-2 ) K = N-2
      DO 20 M = K-1,K+2
        TERM = Y(M)
        DO 10 J = K-1,K+2
          IF ( M .NE. J ) THEN
            TERM = TERM * (XINT - X(J))/(X(M) - X(J))
          END IF
   10   CONTINUE
        YOUT = YOUT + TERM
   20 CONTINUE
      RETURN
      END




C**********************************************************************
C      SUBROUTINE DMAX_ARR (DD,N,DMAX,MAX)
C
C      Maximum element of the vector
C
C**********************************************************************

      SUBROUTINE DMAX_ARR (DD,N,DMAX,MAX)
      real*8 dd(n), dmax
      integer n, max
      dmax=dd(1)
      max=1
      do i=2,n
         if (dd(i).gt.dmax) then
            max=i
            dmax=dd(i)
         end if
      end do
      return
      end

C**********************************************************************
c
c      SUBROUTINE CLOCK_BREAKS (IDNM, IBRKST, IBRKTM)
c
c
c    Clock breaks reading (OT)  - 24.11.2004
c    Update  (OT)  - 25.07.2005
c    Last update  (OT)  - 22.11.2010
c
C**********************************************************************

      SUBROUTINE CLOCK_BREAKS (IDNM, IBRKST, IBRKTM)
      integer*4 ibrkst(10), ibrktm(10)

      OPEN (178,FILE='breaks.res',STATUS='unknown',ERR=99)

      I = 0

322   CONTINUE




	   READ (178,310,END=330) idst, idtm
310      FORMAT (I3,1X,I4)
         print 310, idst, idtm

         if (idst.eq.0.or.idtm.eq.0) then
            print *, ' WARNING  !!!  Check the BREAKS.RES file !!'
         else
            I = I+1
            IBRKST(I) = IDST
            IBRKTM(I) = IDTM
            print *, i, ibrkst(i), ibrktm(i)
         end if


         if (idst.ne.0.and.idtm.ne.0) GO TO 322



330   CONTINUE

      idnm=i

      print *, ' idnm = ', idnm

      CLOSE (178)
      RETURN
99    WRITE (*,*) ' ***** WARNING: File "BREAKS.RES"   not present'

      RETURN
      END


c
c
c

      SUBROUTINE BASELINE_ELIMINATED (IBI, ST1, ST2)
      INCLUDE 'OCCAM.FI'

      CHARACTER*8 ST1(nmax), ST2(nmax)   !    Fixed on 31-Mar-2014

      OPEN (179,FILE='baseline.eli',STATUS='unknown',ERR=99)

      IbI = 0

322   CONTINUE

         IbI = IbI + 1

         READ (179,310,END=330) ST1(ibi), ST2(ibi)
310      FORMAT (A8,2X,A8)
         print 310, st1(ibi), st2(ibi)


         GO TO 322

330   CONTINUE

      ibi = ibi - 1

      print *, ' ibi = ', ibi

      RETURN

 99   PRINT *, ' ERROR in BASELINE.ELI '

      RETURN

      END


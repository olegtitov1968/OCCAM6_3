                                     
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
C
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
C   REVISION April 28, 2004, (OT)  5000 instead of 4000
C   Revision 24, February, 2005, (OT) Clock second derivatives on default;
C                                     Single clock break implementation
C   Revision 17, February, 2006, (OT) 1912 sources are available;
C                 Several lists of global/arc radiosources are optional
C                 Formal errors of the spherical component positions have
C                 been revised
C
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
307   FORMAT (A3,4x,f3.0,1x,f6.3,1x,f6.2)

C     READ IN picorad, PASS TO rad

      SGADD = SGADD/1.d12




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
      DIMENSION IANTENNA(nmax),KANTENNA(nmax),IANT_NST(nmax)
ccc      DIMENSION ECC(nmax,3), IEC(nmax), CONVER(nmax)

C   COMMON LIST. DESCRIPTION OF NEW COMMON AREAS
C
C   * CANTENNA.- INFORMATION ABOUT SELECTED ANTENNAS
C   * PHYS    .- PHYSICAL CONSTANTS

      COMMON /CINDEX/NBAS,NOBS,NSTM,NTIM,NSTA,NSEL,
     *               IBAS,ISTM,ITIM,ISTA,IRAND,NOUT
      COMMON /CANTENNA/ANTENNA,IANTENNA,KANTENNA,NANTENNA,IANT_NST
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




c
c   End of update by Oleg Titov (18-Nov-2013, Friday!!)
c


          SGFAC(IA) = SGFAC(IA)/1.d12


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
 
      print *, ' nmax = ', nmax 

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
               IANT_NST(IA)=Npp1
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
C    WE PLAN TO SOLVE 7 PARAMETERS AT EACH STATION:  (NS = 7)
C       - CLOCK RATE 
C       - ATMOSPHERE rate
C       - ATMOSPHERE delay
C       - TWO ATMOSPHERE GRADIENT rates
C       - TWO ATMOSPHERE GRADIENTS


C    FOR THE FIRST STATION WE WILL SOLVE ONLY THE ATMOSPHERE 
C     (6 parameters)
C

      NR = (NANTENNA-1)*NS + (ns-1)   !  Update on 04-01-2019

      NR1 = NR+6 ! for rotation and scale

    
c      NR1 = NR1 + 1  ! for UT1-UTC


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

      CHARACTER*8 ANTENNA(nmax),STAN,fulsor(nss),arcsour(n_unstable)
      CHARACTER*3 CMODEL
      character*4 axtyp, stc
      CHARACTER*8 SOURCE(nobserv), sour
      CHARACTER*8 st_eli_1(nmax), st_eli_2(nmax)  ! 18.09.2011
	
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

   
      real*8 ut_br(10), year_br(10), ts(10)                 !   22.11.2010
    
      dimension obs(nobs1)                          !   09-Dec-2013
     
      DIMENSION r2(nobs1),z2(nobs1),h2(nobs1*npar), numob(n_unstable)
     

      DIMENSION TIM(nobs1),IDICTIO(nobserv),IOBS(nobserv)
      DIMENSION PCOORD(2,8),PCOORR(2,8),PNUT(2),PPOL(8)
      DIMENSION IANTENNA(nmax),KANTENNA(nmax),KPAR(2)

      DIMENSION IANT_NST(nmax)
            

 
      dimension is1(nobs1), is2(nobs1), year(nobs1),
     *       eps(nobs1), z3(nobs1), covm(npar*(npar+1)/2),
     *   corr(npar*(npar+1)/2), eps1(nobs1), syst(nobs1), omm(nobs1)
      dimension y(npar), sy(npar)
     
      dimension chik(nmax),  wrms(nmax)

      dimension nsta_rew(nmax), eps_res(nmax), nk(nmax), dd(nmax)
      dimension dz2(nobs1), dz1(nobs1), bb(nobs1), z31(nobs1)
      dimension z32(nobs1), db_ut1(nobs1)

      dimension rqu(3), rqp1(3), rqm1(3), omega(nobserv)
      dimension rbar(3,3), vbar(3,3)

      dimension dww1(nobs1), dww2(nobs1), dww3(nobs1), dw(3), p(3)
      dimension p1(nobs1), p2(nobs1), p3(nobs1), dwp(nobs1)

      dimension z1(nobs1), x(npar), sx(npar),eps2(nobs1)
      dimension syst1(nobs1), corr1(npar*(npar+1)/2)
      

      
      real*8 mjd


C  COMMON LIST. DESCRIPTION OF NEW ITEMS
C   * CEXTRACT.- KEEPS INFORMATION ABOUT OBSERVABLES, WEIGHTS AND
C                PARTIAL DERIVATIVES READ FROM OCCAM DATA FILES

      COMMON /CANTENNA/ ANTENNA,IANTENNA,KANTENNA,NANTENNA,IANT_NST
      COMMON /CINDEX/ NBAS,NOBS,NSTM,NTIM,NSTA,NSEL,
     *               IBAS,ISTM,ITIM,ISTA,IRAND,NOUT
      COMMON /DICTION/ TIM,IDICTIO,IOBS
      COMMON /CEXTRACT/ OC,SOBD,OCR,SOBR,PCOORD,PCOORR,
     * PNUT,PPOL
      COMMON /OPT/ SGADD, ZEN_DIST, SGFAC(nmax), TSTFAC
      COMMON /PHYS/ C, FL, AH, AU, J20
      COMMON /MATH/ PI, TWOPI, ROG, ROS, ROZ, TOL


      COMMON /GRAV/ GMS
                      
c

c
c  OUTPUT FILES
c

      open (12, access='append',file = 'WEIGHTS.DAT', status='unknown')
      open (13, file='corr.dat', status = 'unknown')
      open (14, access='append',file = 'STATIST.DAT', status='unknown')
      open (99, file = 'gradients.dat', status = 'unknown')

 
      open (114, access='append',
     * file = 'STATIST_OF_BAD_SESSIONS.DAT', status='unknown')    !  1-April-2016
 
     
      open (77,	file='temp.dat',status='unknown')
      open (88, access='append', file='triangle.dat',status='unknown')
c      open (87, access='append', file='resid.dat',status='unknown')

      open (932,file='test.dat',status='unknown')
 

c      coeff=ros/(100.d0*c)      !  0.0000068755 arcsec*sec/cm

c      do i=1,nantenna

c         antfile(1:8) = antenna(i)
        
c         antfile(9:12) = '.EPS'
        

c         ist = iantenna(i)


c         OPEN (300+ist,access='append',FILE=antfile,status='unknown')
         

c         nk(i) = 0

c      end do

    

      open (1,file='dbc_code',status='unknown')           !  07.11.2008
      read (1,6) iy1, mo, idd, dbc_code                    !  07.11.2008
 6    format (i2.2,a3,i2.2,a2)                                !  07.11.2008
      write (*,6) iy1, mo, idd, dbc_code                    !  07.11.2008

      OPEN (28,access='append',FILE='EOP.DAT',status='unknown')
      OPEN (29,access='append',FILE='scale.DAT',status='unknown')


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

         print *, nr1, idnm
         print *, 'Total number of pars after clock breaks = ',

     *   nr1 + idnm   !  Only one parameter per station is added for a break

         nr1 = nr1 + idnm

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
      print *, ' idat_ntim, ut1_ntim/twopi '

      TIME_1 = dble(IDAT_1) + UT1_1/TWOPI
      TIME_NTIM = dble(IDAT_ntim) + UT1_ntim/TWOPI

      tt=(time_ntim + time_1)/2.d0
c      print *, ' time_1 = ', time_1
c      print *, ' time_ntim = ', time_ntim
c      print *, ' tt = ', tt


      isor=0                      


      nf = nfull           


      print *, ' The number of observations is = ', nf
      print *, ' nr1 = ', nr1

      ijj = 0

      do i=1,nss

        
         

         FULSOR(i) = '        '    !  02.09.2016


      end do 


      iomega = 0

      dda = 0.d0
      dde = 0.d0 


      omm_mean = 0.d0

      DO IROW=1,nf 

         

         IRECORD = iobs(IROW)
         READ (17,rec=irecord) JBAS,ISOUR,IST1,IST2
         READ (18,REC=IST1) JSTA,ISTA
         READ (18,REC=IST2) JSTA,ISTB
         READ (19,REC=isour) J,source(irow),ra2000,de2000,Idat0,UT0,DELT  ! Sources name
     *     ,nu, rab, deb, rqp1,
     *          raa, dea, rqm1,
     *          raap, deap, rqu, rbar, vbar,
     *          hsgp1, hp1, hsgm1, hm1, hsg, h, omega(irow)

         if (irow.eq.1) then 

                 TIME_0 = dble(IDAT0) + UT0/TWOPI
                 omega_0 = omega(irow)

         end if


         if (irow.eq.nf) then

             TIME_N = dble(IDAT0) + UT0/TWOPI
             omega_n = omega(irow)

         end if


         TT_irow = IDAT0 + UT0/TWOPI

         t1= tt_irow - tt


         if (idnm.gt.0) then

            do i=1,idnm
               if (isour.eq.ibrktm(i)) then

			    irow_break(i) = irow

               end if
            end do


         end if


         is1(irow)=ista
         is2(irow)=istb

         i_st1(irow) = kantenna(is1(irow))      !        for collocation
         i_st2(irow) = kantenna(is2(irow))      !        for collocation

         IF (KANTENNA(ista).NE.0.AND.KANTENNA(istb).NE.0) THEN

            KPAR(1) = ista
            KPAR(2) = istb

           CALL EXTRACT (tt_irow, irecord, source(irow),
     *       ra2000,de2000, antenna(i_st1(irow)), antenna(i_st2(irow)),  
     *       cmodel, z31(irow), z32(irow), bb(irow), 
     *       dz2(irow), dz1(irow), z1(irow), omm(irow) )    

c
c  freq0 - analytical model for delay rate
c

             
            z1(irow) = 1.d12 * z1(irow)    !   freq00
                                                                   

C   FILL THE OBSERVABLES VECTOR AND THE WEIGTHING MATRIX (DIAGONAL)
C     ADD EXTRA NOISE THAT MUST BE READ FROM AN OPTIONS FILE

            Z2(IROW) = 1.d12 * OCR          !  in prad   (delay rate)
            

            omm_mean = omm_mean + omm(irow)   !  4-June-2024
                      
            SGTOT = SOBR**2 + SGADD**2           !  SGADD to be checked

c            print *, diff_pot

c            print *, sobr, sgadd, dsqrt(sgtot)

C    FILL THE ROTATION

            
            H2 (IROW + (NR1-3)*nfull) = PPOL(1)
            H2 (IROW + (NR1-2)*nfull) = PPOL(2)
            H2 (IROW + (NR1-1)*nfull) = PPOL(3)

            H2 (IROW + (NR1-6)*nfull) = PPOL(4)
            H2 (IROW + (NR1-5)*nfull) = PPOL(5)
            H2 (IROW + (NR1-4)*nfull) = PPOL(6)

                                

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

               YEAR_1 = dble(IY)+(dble(ID-1)+dble(IH)/23.93447d0+
     *          dble(IM)/1440.d0+dble(S)/86400.d0)/(dble(IFAC)+0.2422d0)
              

               SGTOT = SGTOT + SGFAC(NKANT) **2 

C    FILL JACOBIAN MATRIX


               IF (NKANT.EQ.1) THEN                  !  Reference station

                                    
                  H2(irow) = ISIG*PCOORD(K,6)       !    tropo rate
                  H2(irow +  nfull) = ISIG*PCOORR(K,6)       !    tropo delay
                  H2(irow +  2 * nfull) = ISIG*PCOORD(K,7)       !  N-S atmosphere grad. rate
                  H2(irow +  3 * nfull) = ISIG*PCOORD(K,8)       !  E-W atmosphere grad. rate
                  H2(irow +  4 * nfull) = ISIG*PCOORR(K,7)       !  N-S atmosphere grad.
                  H2(irow +  5 * nfull) = ISIG*PCOORR(K,8)       !  E-W atmosphere grad.
 
               ELSE                                 !  No reference station

                  
                  IJ = IROW + (NKANT-2)*ns * nfull + (ns-1) *nfull
                  


C--------------------------------------------------------------------
C   PARTIALS WITH RESPECT TO CLOCK OFFSET, CLOCK RATE AND ATMOSPHERE
C--------------------------------------------------------------------

                  if (idnm.gt.0) then                      !  Clock break found

                     do i=1,idnm

                        if (nkant.eq.kantenna(ibrkst(i)) ) then


                           if (isour.le.ibrktm(i)) then
                         
    

                              H2(IJ + 6*nfull)=ISIG   !  clock offset


c
c     Partials for clock break - 3 (or more) columns before EOP (nr1-4, etc)
c


                              H2 (IROW + (NR1-(i+3))*nfull) = 0.d0      !  clock rate

                          else if (isour.gt.ibrktm(i)) then

                             

                              H2(IJ + 6*nfull)=0.d0                   !  clock rate


c
c     Partials for clock break - 3 (or more) columns before EOP (nr1-4, etc)
c

                              H2 (IROW + (NR1-(i+3))*nfull) = ISIG                !  clock rate


                           end if

                        else                    !  dealing with problem station

                           
                           H2(IJ + 6*nfull)=ISIG                   !  clock rate

                        end if

                     end do

                  else                                     ! No clock breaks
                    

                     H2(IJ + (ns-1) * nfull)=ISIG                   !  clock rate
   

                  end if

                  H2(IJ) = ISIG*PCOORD(K,6)            !  tropo rate   
                  H2(IJ + nfull) = ISIG*PCOORR(K,6)            !  tropo rate   
                  H2(IJ + 2 * nfull) = ISIG*PCOORD(K,7)  !  N-S atmosphere grad. rate
                  H2(IJ + 3 * nfull) = ISIG*PCOORD(K,8)  !  E-W atmosphere grad. rate
                  H2(IJ + 4 * nfull) = ISIG*PCOORR(K,7)  !  N-S atmosphere grad.
                  H2(IJ + 5 * nfull) = ISIG*PCOORR(K,8)  !  E-W atmosphere grad.

               END IF

            END DO
  
               
                         
c            R2(irow) = sgtot * 1.d24 +  (sgadd)**2 * 1.d24 ! for mnk

            r2(irow) = sgtot * 1.d24

            
            

            if (mo.eq.'JUN') imo = 06
            if (mo.eq.'JUL') imo = 07
            iy2 = 2000 +iy1


            if (dabs(s-60.d0).le.0.01) then

                im = im+1
                s=0.d0

            end if
      

         ENDIF


c         IF (irow.le.5) print *, z2(irow), h2(irow + (nr1-1)*nfull)


      ENDDO

      print *, ' ijj = ', ijj 
 
      tt = (time_0+time_n)/2.d0


      omega_s = 72921151.4671d0

      omm_mean = 1.d12 * omm_mean/dble(nf)

      print *, omega_s - omm_mean
      

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




      nr2 = nr1 

      print *, ' nr1 = ', nr1, ' nr2 = ', nr2, nfull


      do i=1,nantenna
         
         wrms(i) = 0.d0
         chik(i) = 0.d0
         eps_res(i) = 0.d0

      end do

      schi = 0.d0
      sw = 0.d0

      call mnk1 (nr1, nr2, nfull, year, r2, z2, h2, schi, sw, eps1,
     *  y, sy, syst, corr)

      print *, '  mnk has finished '

      n_i=0
c
c     Downweighting of outliers after LSM
c
      print *, ' nf = ', nf

      do i=1, nf

         eps_res(i_st1(i)) = eps_res(i_st1(i)) + (eps1(i)**2)/r2(i)
         eps_res(i_st2(i)) = eps_res(i_st2(i)) + (eps1(i)**2)/r2(i)

         dd(i_st1(i)) = dd(i_st1(i)) + 1.d0/r2(i)
         dd(i_st2(i)) = dd(i_st2(i)) + 1.d0/r2(i)

         nk(i_st1(i)) = nk(i_st1(i)) + 1
         nk(i_st2(i)) = nk(i_st2(i)) + 1

      end do

      do i=1,nantenna

         nsta_rew(i)=0
         wrms(i)=dsqrt(eps_res(i)/dd(i))
         chik(i)=dsqrt(eps_res(i)/dble(nk(i)))

      end do



      print *, ' schi = ', schi, ' sw = ', sw
  
      


      do i=1,nantenna

         write (12,88) antenna(i), nk(i), chik(i), wrms(i)
         print 88, antenna(i), nk(i), chik(i), wrms(i)

      end do

 88   format (1x,a8,4x,i4,2(3x,f8.4))


   


c      do i=1, nfull

c        write (77,777) (h2(nfull*(j-1)+i),j=1,nr1)
c         write (77,777) (h2(nfull*(j-1)+i),j=1,nr2)
c      end do

c777    format (200(1x,f18.6))


      if (schi.gt.100.d0.and.sw.gt.100.d0) then

         print *, ' BAD SESSION '
         write (114,*) tt, 'schi = ', schi, 'sw = ', sw

      else


         print *, ' nr1 = ', nr1, ' nf = ', nf, ' nr2 = ', nr2


c     Check for statistics - bad session do not corrupt full solution
c
                      
  
        write (14,103) tt, nf, schi, sw



 103  format (f11.4,1x,i5,2(2x,f8.4))

             
  
          s_omee = dsqrt(sy(nr1-2)**2 + sy(nr1-1)**2+sy(nr1)**2)
          omega1 = dsqrt(y(nr1-2)**2 + y(nr1-1)**2 + y(nr1)**2)

c          day_length = -86400.d0 * (omee - omega_s)/omega_s
c          day_length = -86400.d0 * omee/omega_s

           day_length = -1.d3 * 86400.d0 * (omega1-omega_s)/omega_s
         
          omm_mean = 365.d0 * omm_mean/366.d0

          write(28,105) tt, y(nr1), y(nr1-1), 
     *  y(nr1-2), sy(nr1), sy(nr1-1), sy(nr1-2),
     *  omega1, s_omee, omega1 - omega_s, day_length,
     *  iy1, mo, idd, dbc_code, year_1 

          write(29,106) tt, y(nr1-4), y(nr1-5), 
     *  y(nr1-6), sy(nr1-4), sy(nr1-5), sy(nr1-6),
     *  iy1, mo, idd, dbc_code, year_1 


 105  format (f12.5,10(1x,f14.4),2x,i2.2,a3,i2.2,a2,2x,f12.5)
 106  format (f12.5,6(1x,f14.4),2x,i2.2,a3,i2.2,a2,2x,f12.5)

         do i=1, nr1

             write (13,601) (corr(i*(i-1)/2+j),j=1,i)

601          format (100(1x,f5.2))

            


         end do





      end if  


  

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
      DIMENSION IANTENNA(nmax),KANTENNA(nmax),IANT_NST(nmax)

      COMMON /CANTENNA/ ANTENNA,IANTENNA,KANTENNA,NANTENNA,IANT_NST

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
C    SUBROUTINE EXTRACT  (tt, IRECORD, CMODEL,
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
 
      SUBROUTINE EXTRACT (tt, IRECORD, sour, ra, de, ant1, ant2,
     * CMODEL, z31, z32, bba, dz2, dz1, freq00, omm )

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
      DIMENSION PEOP1(3,3),PEOP2(3,3),PPOL(8)
     
      DIMENSION PNUT1(2),PNUT2(2),PNUT(2),drda(3),drdd(3),v(3)
      dimension acc_der(3), ww1(3), ww2(3), p(3), b(3), dw(3),
     * a1(3), a2(3), p1(3), da(3)

      DIMENSION ppolar_rate1(3), ppolar_rate2(3)
      DIMENSION azimut1(3), azimut2(3)

      
      
      COMMON /CINDEX/NBAS,NOBS,NSTM,NTIM,NSTA,NSEL,
     *               IBAS,ISTM,ITIM,ISTA,IRAND,NOUT
      COMMON /CEXTRACT/ OC,SOBD,OCR,SOBR,PCOORD,PCOORR,
     * PNUT,PPOL
      COMMON /PHYS/ C, FL, AH, AU, J20
      COMMON /MATH/ PI, TWOPI, ROG, ROS, ROZ, TOL
      COMMON /OPT/ SGADD,ZEN_DIST,SGFAC(nmax),TSTFAC


      open (323, file = 'occam.dat',status = 'unknown')


C   READ ON BASTIM THE RECORD "IRECORD", READ FROM IDICTIO
C   IF CONTAINS THE THEORETICAL AND OBSERVED DELAY AND THE
C   SIGMA OF THE LATTER

      READ (17,rec=IRECORD) JBAS,ISOUR,IST1,IST2,D1,SD1,R01,SR1,
     *           DI1,SDI1,RI1,SRI1,TP1,T0,TM1,
     *           TAUG_P,TAUG_0,TAUG_M,part_gamma,teta,      !   TETA added on 22 May, 2014
     *           cos_A, bb
     *           , phi_A, aa
     *           , freq0, diff_a, RSun
c     *          , acc_der1, acc_der2, acc_der3
     *           , v1, v2, v3
     

                  


C   READ THE STATION INFORMATION AT THE TIME FOR BOTH ANTENNAS
C   THIS INCLUDES THE PARTIAL DERIVATIVES, THE ATMOSPHERIC MODELS
C   AND THE AXIS OFFSETS

      READ (18,REC=IST1) JSTA,ISTA,ISOUR,TEMP1,PRES1,HUMD,IFAC,CABL1,
     *                 CRX1,CRY1,CRZ1,RX1,RY1,RZ1,
     *                 CGLONG1,CPHI1,CELHGT1,Z11,HA1,
     *                 HSTA1,Z21,HA2,HSTA2,Z31,HA3,HSTA31,
     *                 ZDRY1,
     *                 cmniel1_w,cmniel1_d,
     *                 w11, w12, w13,              !  CRF velocity of first station
     *                 a11, a12, a13,
     *                 AXKT1,
     *                 (PCOORD(1,I),I=1,5),(PCOORR(1,I),I=1,5),
     *                 PNUT1,PEOP1,
     *                 ppolar_rate1,
     *                 therm_d1,azimut1,a_grad1,
     *                 dr1_ut1,                      !   01.01.2021
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
     *            ,vmf1h1d,vmf1w1d          !  for delay rate
     

      READ (18,REC=IST2) JSTA,ISTB,ISOUR,TEMP2,PRES2,HUMD,IFAC,CABL2,
     *                 CRX2,CRY2,CRZ2,RX2,RY2,RZ2,
     *                 CGLONG2,CPHI2,CELHGT2,Z12,HA1,
     *                 HSTA1,Z22,HA2,HSTA2,Z32,HA3,HSTA32,
     *                 ZDRY2, 
     *                 cmniel2_w,cmniel2_d,
     *                 w21, w22, w23,             !  CRF velocity of second station
     *                 a21, a22, a23, 
     *                 AXKT2,
     *                 (PCOORD(2,I),I=1,5),(PCOORR(2,I),I=1,5),
     *                 PNUT2,PEOP2,
     *                 ppolar_rate2,
     *                 therm_d2,azimut2,a_grad2,
     *                 dr2_ut1,                   !  01.01.2021
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
     *            ,vmf1h2d,vmf1w2d          !  for delay rate
     

      
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



      a2(1) = a21
      a2(2) = a22
      a2(3) = a23

      a1(1) = a11
      a1(2) = a12
      a1(3) = a13  


      v(1) = v1
      v(2) = v2
      v(3) = v3
     

      vc = dotpr(v,v)/(c**2)            !  V**2/c**2
                                         
      

      tropor = 0.d0
      b1 = 0.d0
      b2 = 0.d0


      b(1) = crx2(3) - crx1(3)
      b(2) = cry2(3) - cry1(3)
      b(3) = crz2(3) - crz1(3)
 
      

      w0 = dsqrt(w11**2 + w12**2 + w13**2)
      r0 = dsqrt(crx1(3)**2 + cry1(3)**2 + crz1(3)**2)

      w0 = w0/dcos(cphi1(3))

      omm = w0/r0

    
      write (323,201) tt, 1.d12 * omm, w0, r0

201   format (f12.5,2x,f12.2,2x,f10.6,2x,f12.4,2x,f12.2)


      DO I=1,3

         
         
         dw(i) = ww2(i) - ww1(i)
         da(i) = a2(i) - a1(i)

 
         

      END DO

      bba = dsqrt(dotpr(b,b))
                       
      w11 = dsqrt(dotpr(ww1,ww1))
      w22 = dsqrt(dotpr(ww2,ww2))
  

      dag = 180.d0
      ddg =  0.d0
        
      dag = pi * dag/180.d0
      ddg = pi * ddg/180.d0

      dgd = 0.d0

c
c   DGD is cos of angular distance
c


      dgd = dcos(de)*dcos(ddg)*dcos(ra - dag) + dsin(de)*dsin(ddg)

      dgd = dacos(dgd)

      sdgd2 = dsin(dgd)**2/2.d0

      
      freq00 = -dotpr(dw,p) 
      


      
c
c     MULTIPLICATION ON SEC ( zenit distance )
c
c           BY O.TITOV
c
      z31_1=z31*180.d0/pi
      z32_1=z32*180.d0/pi

c      print *, z31_1, z32_1, zen_dist

      if (z31_1.gt.zen_dist.or.z32_1.gt.zen_dist) then

  
         if (z32_1.le.z31_1) then

c              print *, sr1

              sd1 = sd1/dcos(z31)
c              sr1 = 1.d0 * sr1/dcos(z31)      !   08-Jan-2019 (OT)
                          
c              print *, sr1

         else


                  
              sd1 = sd1/dcos(z32)
c              sr1 = 1.d0 * sr1/dcos(z32)     !    08-Jan-2019 (OT)



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

         TROPO = (ZDRY2 *cmniel2_d(1)-ZDRY1 *cmniel1_d(1))/C            ! sec
         TROPOR= (ZDRY2 *cmniel2_d(2)-ZDRY1 *cmniel1_d(2))/C      ! rad
         PCOORD(1,6) = cmniel1_w(1)
         PCOORD(2,6) = cmniel2_w(1)

      ENDIF

      IF (CMODEL.EQ.'VMF') THEN
         	
         TROPO = (ZDRY2 *dqmfh2-ZDRY1 *dqmfh1)/C                        ! sec


         TROPOR= (ZDRY2 *cmniel2_d(2)-ZDRY1 *cmniel1_d(2))/C      ! rad

         PCOORD(1,6) = dqmfw1
         PCOORD(2,6) = dqmfw2

      ENDIF



      IF (CMODEL.EQ.'VM1') THEN

         TROPO = (ZDRY2*vmf1h2 - ZDRY1*vmf1h1)/C                        ! sec
                      

         dz2 =  z22 -  z12
         dz1 =  z21 -  z11

         daz1 = azimut1(2) - azimut1(1)
         daz2 = azimut2(2) - azimut2(1)



   

         tropor2 = ZDRY2 * cmniel2_d(2) * dz2/2.d0     ! m/s
 
         tropor1 = ZDRY1 * cmniel1_d(2) * dz1/2.d0      ! m/s


         tropor = (tropor2 - tropor1)/c                       !  rad
  

         PCOORD(1,6) = vmf1w1      !   partial for tropo delay first station
         PCOORD(2,6) = vmf1w2      !   partial for tropo delay second station

         
         PCOORR(1,6) = cmniel1_w(2)   !   derivatives for delay rate station #1
         PCOORR(2,6) = cmniel2_w(2)   !   derivatives for delay rate station #2

         
      ENDIF


         
      PCOORD(1,7) = a_grad1(1)*pcoord(1,6)      !  North-South   1st station
      PCOORD(2,7) = a_grad2(1)*pcoord(2,6)      !  North-South   2nd station

   
      PCOORD(1,8) = a_grad1(2)*pcoord(1,6)      !  East-West     1st station
      PCOORD(2,8) = a_grad2(2)*pcoord(2,6)      !  East-West     2nd station


      PCOORR(1,7) = a_grad1(1)*pcoorr(1,6)*dz1      !  North-South   1st station
     *    - pcoord(1,6) * a_grad1(2) * daz1
     *    + pcoord(1,6) * a_grad1(1) * dz1/(dsin(z31)*dcos(z31)) 

      PCOORR(2,7) = a_grad2(1)*pcoorr(2,6)*dz2      !  North-South   2nd station
     *    - pcoord(2,6) * a_grad2(2) * daz2
     *    + pcoord(2,6) * a_grad2(1) * dz2/(dsin(z32)*dcos(z32)) 
   
      PCOORR(1,8) = a_grad1(2)*pcoorr(1,6)*dz1      !  East-West     1st station
     *    + pcoord(1,6) * a_grad1(1) * daz1
     *    + pcoord(1,6) * a_grad1(2) * dz1/(dsin(z31)*dcos(z31))

      PCOORR(2,8) = a_grad2(2)*pcoorr(2,6)*dz2      !  East-West     2nd station
     *    + pcoord(2,6) * a_grad2(1) * daz2
     *    + pcoord(2,6) * a_grad2(2) * dz2/(dsin(z32)*dcos(z32)) 


c      PCOORR(1,7) =       !  North-South   1st station
c     *            - pcoord(1,6) * a_grad1(2) 
c     *            + pcoord(1,6) * a_grad1(1)/dcos(z31)

c      PCOORR(2,7) =       !  North-South   2nd station
c     *            + pcoord(2,6) * a_grad2(2) 
c     *            + pcoord(2,6) * a_grad2(1)/dcos(z32)
   
c      PCOORR(1,8) =      !  East-West     1st station
c     *            - pcoord(1,6) * a_grad1(1) 
c     *            + pcoord(1,6) * a_grad1(2)/dcos(z31)

c      PCOORR(2,8) =       !  East-West     2nd station
c     *            + pcoord(2,6) * a_grad2(1) 
c     *            + pcoord(2,6) * a_grad2(2)/dcos(z32)


C   GET AXIS OFFSET CORRECTIONS

      AXIS  = AXKT2(1)-AXKT1(1)
      AXISR = (AXKT2(2)-AXKT1(2))/1.d9  !        nrad -> rad

      
C   GET THERMIC DEFORMATION CORRECTIONS

      therm = therm_d1-therm_d2    !  Sign was changed on 18.09.2006 (OT)

      

C   COMPUTE THEORETICAL DELAY ADDING THE GEOMETRICAL MODEL, THE
C  TROPOSPHERIC CONTRIBUTION AND THE AXIS OFFSET CORRECTION


      THDELAY = T0 + TROPO + AXIS + therm    !  to be checked for all types of antennas


      THRATE = (TP1-TM1)/2.D0    ! sec/sec = rad  23 Sep 2015

      TGEOM_1 = (TAUG_P - TAUG_M)/2.d0


      FRAD = 180.D0*3600.D0/PI                   !   206265 arcsec

C   GET OBSERVABLE MINUS CALCULUS

      OC = OBD - THDELAY                    !  sec

c      OCR = OBR - THRATE  - AXISR  - tropor                !  rad


      OCR = OBR - freq0 - axisr - tropor           ! rad

c      print *, obr, freq0, ocr


           
C   GET PARTIAL DERIVATIVES WITH RESPECT TO BARYCENTRE ACCELERATION
       

c      PPOL(1) =  - 1.d10 * (B(1) - P(1) * DOTPR(B,P))/(c*c)          !  sec**2/(1.d10*m)
c      PPOL(2) =  - 1.d10 * (B(2) - P(2) * DOTPR(B,P))/(c*c)          !  sec**2/(1.d10*m)
c      PPOL(3) =  - 1.d10 * (B(3) - P(3) * DOTPR(B,P))/(c*c)          !  sec**2/(1.d10*m)

c      print *, ppol(1), ppol(2), ppol(3)


C   GET PARTIAL DERIVATIVES WITH RESPECT TO ROTATION

     
  
      PV1 = P(1) * DOTPR(P,V) - V(1)
      PV2 = P(2) * DOTPR(P,V) - V(2)
      PV3 = P(3) * DOTPR(P,V) - V(3)



      PPOL(1) =  - (B(2) * P(3) - B(3) * P(2))/c    !  sec
      PPOL(2) =  - (B(3) * P(1) - B(1) * P(3))/c    !  sec
      PPOL(3) =  - (B(1) * P(2) - B(2) * P(1))/c    !  sec

      PPOL(4) =    B(3) * P(3)/c    !  sec
      PPOL(5) =    B(2) * P(2)/c    !  sec
      PPOL(6) =    B(1) * P(1)/c    !  sec


 

      db_ut1 = dr2_ut1 - dr1_ut1                    !  1/rad

      db_ut1 = twopi * 1.d12 * db_ut1 / 86400.d0    !  prad/sec

             
c      print *, dr1_ut1, dr2_ut1, db_ut1

        
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


      
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
C   Revision  18, July, 2001, (OT)  Maximum of observables is 3310
C                                        instead of 2000
C   Revision  29, November, 2001, (OT)  PSOU calculated as elements
C                  of partial derivatives matrixes for source coordinates
C                  adjustment
C   REVISION March 05, 2002, (OT)   812 radiosources in list,
C                                        3310 observations
C   REVISION May 20, 2002, (OT) Weighting of constraints,
C                  Limits for observations increased to 20000 (IDICTIO, IOBS)
C   Revision 11,June,2002, (OT)  Max. is 4000 observations instead of 3310
C                                     NFULL subroutine has been deleted
C   Revision 05, July, 2002, (OT) NNT and NNR approaches for ITRF treatment
C                                Direction and vertical constraints are
C                                not fixed more
C   Revision 10, July, 2002, (OT) COMMON-BLOCK /CKALMAN/ canceled
C   Revision 28, April, 2004, (OT) 5000 observations
C   Revision 22, July, 2004, (OT) Arc sources estimated as daily parameters;
C          Clock second term and troposphere gradients are estimated on default
C          250 daily parameters are possible (OCCAM.FI)
C   Revision 24, February, 2005, (OT) Clock second derivatives check;
C                                     Single clock break implementation
C   REVISION JUNE 25, 2007  (OT) - NNR parameter is read from file
C                                       partials for acceleration revised
C   REVISION JANUARY 19, 2009  (OT) - The error for partial derivatives
C                                          fixed
c   Revision: 05, September, 2009 (OT)  17 extra global parameters
c      (3 for Galaxy, 3 for rotation, 10 for Grav waves, 1 for gamma, 1 for parallax)
c
c   Revision: 25, June, 2010 (OT)  Large update on partials (Auckland)
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
C***************************************************************************
 
C 
C   SPECIFICATIONS
C

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'OCCAM.FI'

      CHARACTER CMODEL*3, NNR*6     !  Change of 25.07.2007

C
C    COMMON LIST.  DESCRIPTION
C
C    * CINDEX .- LIST OF NUMBERS AND DIMENSIONS
C    * MATH   .- MATHEMATICAL CONSTANTS
C    * SGM    .- SIGMA STORED INFORMATION
C

      COMMON /CINDEX/ NBAS,NOBS,NSTM,NTIM,NSTA,NSEL,
     *               IBAS,ISTM,ITIM,ISTA
      COMMON /MATH/ PI, TWOPI, ROG, ROS, ROZ, TOL
      COMMON /OPT/ SGADD,ZEN_DIST,SGFAC(nmax),TSTFAC
      COMMON /PHYS/ C, FL, AH, AU, J20

C
C    GET PARAMETERS, OPEN FILES, INITIALIZE
C

      CALL PARAM        !  Calls subroutine from ARC.FOR
      CALL KOPEN
      CALL READOPT (CMODEL, NNR)            !  Change of 25.07.2007
      CALL KINIT  (NR1)

C
C    READ STOCHASTIC PROCESSES INFORMATION
C    REORDER THE CONTENTS OF DICTIO
C

      CALL KTRANS
      CALL PREPARE (CMODEL, NNR, NR1)       !  Change of 25.07.2007


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
     *               IBAS,ISTM,ITIM,ISTA

C
C   OPEN FILES
C

      open (12, access='append',file = 'details.dat', status='unknown')
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
C    REVISION MARCH 14, 2001  (OT)
C    LAST REVISION JUNE 25, 2007  (OT) - NNR parameter is read from file
C
C**********************************************************************

      SUBROUTINE READOPT (CMODEL, NNR)

C   ESPECIFICATIONS

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'OCCAM.FI'

      CHARACTER CMODEL*3, NNR*6

      COMMON /OPT/ SGADD, ZEN_DIST, SGFAC(nmax),TSTFAC
      COMMON /PHYS/ C

C    READ FILE WITH SETUP DATA


      READ (*,307) CMODEL, ZEN_DIST, SGADD, TSTFAC, NNR
      write (*,307) cmodel, zen_dist, sgadd, tstfac, nnr
      write (12,307) cmodel, zen_dist, sgadd, tstfac, nnr
307   FORMAT (A3,4x,f3.0,1x,f6.3,2x,f5.2,2x,a6)

C     READ IN METERS, PASS TO SECONDS

      SGADD = SGADD/C

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

      CHARACTER STAT*8,STATCAT*8,ANTENNA(nmax)*8,ANTE1*8,ANTE2*8,CVER*1
      DIMENSION IANTENNA(nmax),KANTENNA(nmax),CONVER(nmax)

C   COMMON LIST. DESCRIPTION OF NEW COMMON AREAS
C
C   * CANTENNA.- INFORMATION ABOUT SELECTED ANTENNAS 
C   * PHYS    .- PHYSICAL CONSTANTS

      COMMON /CINDEX/NBAS,NOBS,NSTM,NTIM,NSTA,NSEL,
     *               IBAS,ISTM,ITIM,ISTA
      COMMON /CANTENNA/ANTENNA,IANTENNA,KANTENNA,NANTENNA
      COMMON /PHYS/ C, FL, AH, AU, J20
      COMMON /OPT/ SGADD,ZEN_DIST, SGFAC(nmax),TSTFAC
      COMMON /MATH/ PI, TWOPI, ROG, ROS, ROZ, TOL
c      COMMON /ECC/ ECC(nmax,3),IEC(nmax)


C   SET VARIABLES TO ZERO

      NVER = 0

C   SET ARRAYS TO ZERO

      DO I=1,nmax
         ANTENNA(I)='        '
         IANTENNA(I)=0
         KANTENNA(I)=0
         CONVER(I) = 0.D0
      END DO
 
C
C   GET WHICH ANTENNAS ARE TO BE SOLVED
C     READ INFORMATION FROM FILE 'STATIONS' WITH A LOOP
C

      NANTENNA=0

      READ (23,389,ERR=301,END=301) NANTENNA
389   FORMAT (13X,I2)

      DO IA=1,NANTENNA
         IFLAG=0


c   Start of update by Oleg Titov (25-June-2010, Auckland)

c         READ (23,100,ERR=302,END=302) STAT,SGFAC(IA),CVER,CONVER(IA),
c     *   ECC(IA,1),ECC(IA,2),ECC(IA,3)
c100      FORMAT (A8,2X,F6.3,2X,A1,2X,F8.5,3F8.3)

         READ (23,100,ERR=301,END=301) STAT,SGFAC(IA)



100      FORMAT (A8,2X,F6.4)

c
c   End of update by Oleg Titov (25-Jun-2010, Auckland)
c

C       IF ECCENTRICITIES

c          IEC(IA)=0
c          DO II=1,3
c             IF(ECC(IA,II).NE.0.D0)THEN
c                IEC(IA)=1
c                GO TO 110
c             ENDIF
c          ENDDO

c110       CONTINUE

          SGFAC(IA) = SGFAC(IA)/C
          CONVER(IA) = CONVER(IA)/C

          IF (STAT.NE.'       ') THEN

C     CHECK IF THIS ANTENNA IS IN THE CATALOG

          DO I=1,NSTA

             READ (20,REC=I) JST,STATCAT

             IF (STATCAT.EQ.STAT) THEN
                ANTENNA(IA)=STAT
                IANTENNA(IA)=I
                KANTENNA(I)=IA
                IFLAG=1

             ENDIF
          END DO

C    WARNING IF ANTENNA NOT IN CATALOG

          IF (IFLAG.EQ.0)

     *          WRITE (*,*) ' ***** THE ANTENNA ',STAT,' NOT FOUND'

          END IF

C   END OF LOOP

      ENDDO

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

      NR=(NANTENNA-1)*NS + (ns-3)

      NR1=NR
      NR1=NR1+2 ! for Nutation

      NR1=NR1+3 ! for EOP

      RETURN

C   WARNING AND BREAK IS COLLOCAT.OPT IS WRONG OR NOT COMPLETE

301   WRITE (*,*)
     +' File COLLOCAT.OPT has errors in ANTENNA line. UNABLE TO PROCEED'
302   WRITE (*,*)
     +' File COLLOCAT.OPT has errors in station line. UNABLE TO PROCEED'
303   WRITE (*,*)
     +'File COLLOCAT.OPT has errors in DIRECTION line.UNABLE TO PROCEED'
304   WRITE (*,*)
     +' File COLLOCAT.OPT has errors in ST1-ST2 line. UNABLE TO PROCEED'
      STOP

      END


C**********************************************************************
C
C   SUBROUTINE prepare (cmodel, nnr, nr1)
C
C   THIS ROUTINE READS THE OBSERVABLES, THEIR WEIGHTING MATRIX AND
C   THE JACOBIAN MATRIX
C
C   REVISION MARCH 17, 1998
C   REVISION January 21, 2001, (OT)   Troposphere gradients added
C   REVISION September 14, 2001, (OT)   DELT instead of IDELT
C   REVISION March 05, 2002, (OT)   812 radiosources in list
C   REVISION May 20, 2002, (OT) Weighting of constraints,
C                  Limits for observations increased to 20000 (IDICTIO, IOBS)
C   REVISION May 15, 2003, (OT)   810 radiosources in list
C   REVISION NOVEMBER 13, 2003, (OT)   783 radiosources in list
C   REVISION DECEMBER 01, 2003, (OT)   769 radiosources in list
C
C   REVISION JUNE 21, 2004, (OT)   NNR-constraints are for 207 ICRF sources
C   REVISION MAY 20, 2005, (OT)  Statistic test before writing
C   REVISION July 21, 2005, (OT) Input for amount of unstable and stable
C                                  quasars from the INPUT8.TXT file
C   REVISION March 01, 2006, (OT) More handy selection of unstable
C       radiosources;
C   REVISION March 22, 2006, (OT) Cable correction implemented;
C                                      One break per session
C   LAST REVISION October 04, 2007, (OT) Octople systematic added
C
C**********************************************************************

      SUBROUTINE prepare (CMODEL, NNR, NR1)       !  Change of 25.07.2007

C   ESPECIFICATIONS

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'OCCAM.FI'

      CHARACTER*4 AXTYP, stc
      CHARACTER CMODEL*3,ANTENNA(nmax)*8,STAN*8,FULSOR(nss)*8,NNR*6 !  Change of 25.07.2007
      CHARACTER*8 SOURCE, source_1
      CHARACTER*8 st_eli_1(nmax), st_eli_2(nmax)  ! 18.09.2011

      character*8 arcsour(n_unstable), ICRF107 (n_unstable)
      character*8 ICRF207(n_stable)
 
      integer*4 i_st1(nobs1), i_st2(nobs1)

      integer*4 ibrkst(10), ibrktm(10)                        !   24.11.2004

      integer*4 k_br(10), idat_br(10), irow_break(10)       !   22.11.2010
      real*8 ut_br(10), year_br(10), ts(10)                 !   22.11.2010
 
      integer*4 nsig_bas(nmax,nmax)                         !   17.09.2011
      integer*4 num_eli(nmax,2)                             !   17.09.2011
      dimension dt(nobs1)                                   !   10.11.2011
      integer*2 im, i_n, isor, num1
      integer*4 inn(nobs1), in, km(nobs1), iarc(nobs1), num2(nobs1)


      DIMENSION r2(nobs1),z2(nobs1),h2(nobs1*npar),part_gamma(nobs1)
     
      DIMENSION TIM(nobs1),IDICTIO(nobserv),IOBS(nobserv)
      DIMENSION PCOORD(2,8),PCOORR(2,8),PNUT(2),PPOL(3),PSOU(2),PPOL1(3)
      DIMENSION IANTENNA(nmax),KANTENNA(nmax),KPAR(2)
      DIMENSION num(nsou), p1(nobs1), p2(nobs1)
      dimension ng(nss), w(nss), number(nss)
      dimension ra(nss), de(nss), ca(nss), sa(nss), cd(nss), sd(nss)
      dimension nsta_rew(nmax), eps_res(nmax), nk(nmax), dd(nmax)
      dimension chik(nmax), wrms(nmax)

      dimension is1(nobs1),is2(nobs1),year(nobs1),eps(nobs1)
      dimension c8(nobs1*nsou*2), cm8(nsou*(2*nsou+1))

      dimension v8(2*nsou), part_parallax(nobs1)
      dimension pg(nobs1,n14-2)              !  for GAL(3)+rot(3)+gamma(1)+quadrupole(10)
c                                 
c     6.07.07
c
      dimension acc_der(3)    !  partials for absolute acceleration
      dimension b(3), p(3)
      
c
c
c
      

C  COMMON LIST. DESCRIPTION OF NEW ITEMS
C   * CEXTRACT.- KEEPS INFORMATION ABOUT OBSERVABLES, WEIGHTS AND
C                PARTIAL DERIVATIVES READ FROM OCCAM DATA FILES

      COMMON /CANTENNA/ ANTENNA,IANTENNA,KANTENNA,NANTENNA
      COMMON /CINDEX/ NBAS,NOBS,NSTM,NTIM,NSTA,NSEL,
     *               IBAS,ISTM,ITIM,ISTA
      COMMON /DICTION/ TIM,IDICTIO,IOBS
      COMMON /CEXTRACT/ OC,SOBD,OCR,SOBR,PCOORD,PCOORR,
     * PNUT,PPOL,PSOU,PPOL1
      COMMON /OPT/ SGADD,ZEN_DIST, SGFAC(nmax), TSTFAC
      COMMON /PHYS/ C, FL, AH, AU, J20
      COMMON /MATH/ PI, TWOPI, ROG, ROS, ROZ, TOL


      open (1,file='dbc_code',status='unknown')
      read (1,6) iy1, mo, idd, dbc_code
6     format (i2.2,a3,i2.2,a2)
      write (*,6) iy1, mo, idd, dbc_code

      close (1)

      open (15, access='append',file = 'INFORM.DAT', status='unknown')

      read (*,*) num_unst, num_sta

      print *, ' The number of unstable sources is ', num_unst
      print *, ' The number of   stable sources is ', num_sta
c      write (12,*) ' The number of unstable sources is ', num_unst
c      write (12,*) ' The number of   stable sources is ', num_sta

      if (num_sta.gt.n_stable) go to 9991
      if (num_unst.gt.n_unstable) go to 9992

      
      if (num_unst.eq.5389)OPEN(24,FILE='arc_5389.CAT',status='unknown')    !  02 Dec 2023
      if (num_unst.eq.57)  OPEN(24,FILE='arc_57.cat',status='unknown')      !  57 'specially handled' sources for ICRF3
      if (num_unst.eq.87)  open(24,file='arc_87.cat',status='unknown')      !  87 special handled
      if (num_unst.eq.2108)  OPEN(24,FILE='low_z.CAT',status='unknown')     !  2108 low z quasars (z<1)
 
      if (num_unst.eq.10.or.num_unst.eq.87.and.num_unst.eq.1495.
     *   or.num_unst.eq.5389.or.num_unst.eq.2108.or.num_unst.eq.57) then

          print *, '  num_unst = ', num_unst
      else

         print *, 'ERROR in num_unst!! '

      end if

      if (num_sta.eq.303) OPEN (25,FILE='303.CAT',status='unknown')   !   ICRF3 (303)


      do k=1,num_unst
         read (24,125) icrf107(k)
      end do

      do k=1,num_sta

         read (25,125) icrf207(k)

      end do
125   format (a8)

      if (num_sta.ne.303)
     *   print *, 'ERROR in num_sta!! ', ' num_sta = ', num_sta

      icc=0

      if (num_unst.eq.2108.and.num_sta.eq.303)  icc=1        ! ICRF2   lowz
      if (num_unst.eq.87. and.num_sta.eq.303)  icc=1        !  aus2024b
      if (num_unst.eq.57. and.num_sta.eq.303)  icc=1        !  ICRF3
      if (num_unst.eq.5389. and.num_sta.eq.303)  icc=1        !  post ICRF3

      if (icc.eq.0) print *, 'The combination of stable are unstable
     * radiosources is not consistent !!! '

      OPEN (26,FILE='col_matr.dat',access='direct',recl=24,
     *  form='formatted')
      OPEN (27,FILE='col_vec.dat',access='direct',recl=24,
     *  form='formatted')
      OPEN (28,access='append',FILE='PARAMS.DAT',status='unknown')
      OPEN (128,access='append',
     * FILE='BAD_DATA_PARAMS.DAT',status='unknown')

c      OPEN (29,FILE='TEST.DAT',status='unknown')

      open (30,access='append',file='statistic.dat',status='unknown')

C   INITIALIZE TO ZERO ALL THE ARRAYS INVOLVED

      CALL KZERO (nfull)

      nfull = nfull+6   !  NNT, NNR constraints for ITRF

      do i=1,nss  !            The number of radiosources

         FULSOR(i)='        '
       
         ng(i)=0

         num2(i) = 0

      end do

      print *, ' The number of radiosources in catalogue ', nss
c      write (12,*) ' The number of radiosources in catalogue ', nss

      do i=1,num_unst
         arcsour(i)='        '
      end do

      CALL CLOCK_BREAKS (IDNM, IBRKST, IBRKTM)

      if (idnm.gt.0) print *, 'The number of breaks is ', IDNM


      if (idnm.gt.0) then


         nr1=nr1+3*idnm

         do i=1,idnm

c            print *, ' ibrktm = ', ibrktm(i)

            k_br(i)=kantenna(ibrkst(i))

      READ (19,REC=IBRKTM(i)) j,source,
     *                ra2000,de2000,Idat_br(i),UT_br(i)

            year_br(i) = dble(IDAT_BR(i)) + UT_BR(i)/TWOPI
c            print 67, i, antenna(k_br(i)), year_br(i)
c 67         format (i4,2x,a8,2x,f12.4)

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

         print *, ibi, ' num_eli = ', num_eli(j,1), num_eli(j,2)

      end do



         narc=0


c
c
c   Change on 25.07.2007
c
      if (NNR.eq.'NNRYES') then

         nfull=nfull+3     !  NNR constraints for ICRF
         nf=nfull-6-3

         print *, ' NNR constraints impose the CRF '

      END IF

      if (NNR.eq.'NNR_NO') then

         nf=nfull-6  !   No NNR constraints for ICRF

         print *, ' NO NNR constraints impose the CRF '

      END IF

c
c    End of change on 25.07.2007
c


      if (nsta.eq.2) then
         nr1=nr1-6        !  No EOP and reference station for single-base network
         nfull=nfull-6    !  NO NNT, NNR for ITRF
      end if

      print *, ' nr1 = ',nr1, ' nf = ', nf, ' nfull = ', nfull

c      im = 0


      do i=1,nf

         inn(i) = 0
         km(i) = 0

      end do

      do i=1,nmax

        nk(i) = 0

      end do


      DO IROW=1,nf

         IRECORD = iobs(IROW)
         READ (17,rec=irecord) JBAS,ISOUR,IST1,IST2
         READ (18,REC=IST1) JSTA,ISTA
         READ (18,REC=IST2) JSTA,ISTB
         READ (19,REC=isour) J,source_1,ra2000,de2000,Idat0,UT0,DELT  ! Sources name

         if (irow.eq.1) TIME_0 = dble(IDAT0) + UT0/TWOPI
         if (irow.eq.nf) TIME_N = dble(IDAT0) + UT0/TWOPI


c         source_1 = source(irow)


         im = 0
         in = 0

                           
         CALL ICRF2000 (source_1, ng, isor, ra2000, de2000,
     *                 fulsor, num, num1, in, ra, de)


         inn(irow) = in

        

         num2(irow) = num1

c         print *, irow, source_1, num1

         do i = 1, num_unst

            if (source_1.eq.icrf107(i)) num2(irow) = 0 

         end do


         
         if (idnm.gt.0) then

            do i=1,idnm
               if (isour.eq.ibrktm(i)) then

			    irow_break(i) = irow

               end if
            end do

c            print *, idnm, i, irow, irow_break(i)

         end if

      
         do j=1,nsou  !            The number of radiosources in the file

            c8(irow+2*(j-1)*nfull) = 0.d0
            c8(irow+(2*j-1)*nfull) = 0.d0

         end do

         is1(irow)=ista
         is2(irow)=istb



         i_st1(irow)=kantenna(is1(irow))      !        for collocation
         i_st2(irow)=kantenna(is2(irow))      !        for collocation

         IF (KANTENNA(IS1(irow)).NE.0.
     *                AND.KANTENNA(IS2(irow)).NE.0) THEN

            KPAR(1) = IS1(irow)
            KPAR(2) = IS2(irow)

c
c  Partial derivatives for acceleration on 6.07.2007
c
            CALL EXTRACT (irecord, cable1, cable2,
     *       ra2000,de2000, cmodel, pgamma, b, p, parallax, teta)


C   FILL THE OBSERVABLES VECTOR AND THE WEIGTHING MATRIX (DIAGONAL)
C     ADD EXTRA NOISE THAT MUST BE READ FROM AN OPTIONS FILE

            part_gamma(irow) = pgamma



            part_parallax(irow) = parallax
            

            
            Z2(IROW)=OC    !   in seconds



            Z2(IROW)=Z2(irow)-cable1*1.d-9      !  Cable correction 22.03.2006
            Z2(IROW)=Z2(irow)+cable2*1.d-9      !  Cable correction 22.03.2006

            SGTOT = SOBD**2 + SGADD**2

C    FILL THE NUTATION OFFSET PART OF THE JACOBIAN MATRIX

            if (nsta.eq.2) then                     !  singlebase  !!!!
               H2 (IROW + (NR1-2)*nfull) = PNUT(1)
               H2 (IROW + (NR1-1)*nfull) = PNUT(2)
            else

	       H2 (IROW + (NR1-5)*nfull) = PNUT(1)
               H2 (IROW + (NR1-4)*nfull) = PNUT(2)
               H2 (IROW + (NR1-3)*nfull) = PPOL(1)
               H2 (IROW + (NR1-2)*nfull) = PPOL(2)
               H2 (IROW + (NR1-1)*nfull) = PPOL(3)

            end if

c    Observational epoch for collocation procedure
c
            TIMDT = dble(IDAT0) + UT0/TWOPI + (32.184+DELT)/86400.d0
            YEAR(irow) = timdt

C    DETERMINE THE COLUMNS FOR ANT1 AND ANT2 FOR THE JACOBIAN MATRIX



            DO K=1,2
               IF (K.EQ.1) ISIG=-1
               IF (K.EQ.2) ISIG=+1
               NKANT = KANTENNA(KPAR(k))


               CALL JULDAT (IY,ID,IH,IM,S,IDAT0,UT0,1)
               IFAC = 365
               IF (MOD(IY,4).EQ.0) IFAC=366

               YEAR_1 = dble(IY)+(dble(ID)+dble(IH)/23.93447d0+
     *        dble(IM)/1440.d0+dble(S)/86400.d0)/(dble(IFAC)+0.2422d0)


C    FILL JACOBIAN MATRIX


               IF (NKANT.EQ.1) THEN

                  SGTOT = SGTOT + SGFAC(NKANT)**2
      

                  IJ=IROW

                  if (nsta.eq.2) then

                     H2(IJ)        =ISIG*PCOORD(K,6)       !  atmosphere
                     H2(IJ+nfull)  =ISIG*PCOORD(K,7)       !  atmosphere grad.
                     H2(IJ+2*nfull)=ISIG*PCOORD(K,8)       !  atmosphere grad.

                  else

                     H2(IJ)        =ISIG*PCOORD(K,1)       !  X
                     H2(IJ+nfull)  =ISIG*PCOORD(K,2)       !  Y
                     H2(ij+2*nfull)=ISIG*PCOORD(K,3)       !  Z

                     H2(IJ+3*nfull)=ISIG*PCOORD(K,6)       !  atmosphere
                     H2(IJ+4*nfull)=ISIG*PCOORD(K,7)       !  atmosphere grad.
                     H2(IJ+5*nfull)=ISIG*PCOORD(K,8)       !  atmosphere grad.

                  end if

               ELSE

                  SGTOT = SGTOT + SGFAC(NKANT)**2



                  if (nsta.eq.2)  then
                     IJ=IROW + 3*nfull
                  else
                     IJ=IROW+(NKANT-2)*ns*nfull +6*nfull  
                  end if


C-------------------------------------------
C   PCOORD(K,1) = PARTIAL WITH RESPECT TO Z
C   PCOORD(K,2) = PARTIAL WITH RESPECT TO X
C   PCOORD(K,3) = PARTIAL WITH RESPECT TO Y
C-------------------------------------------

                  H2(IJ)        =ISIG*PCOORD(K,1)               !  X
                  H2(IJ+nfull)  =ISIG*PCOORD(K,2)               !  Y
                  H2(ij+2*nfull)=ISIG*PCOORD(K,3)               !  Z

C--------------------------------------------------------------------
C   PARTIALS WITH RESPECT TO CLOCK OFFSET, CLOCK RATE AND ATMOSPHERE
C--------------------------------------------------------------------

                  if (idnm.gt.0) then              !   if clock break presents

                     do i=1,idnm

                        if (nkant.eq.kantenna(ibrkst(i)) ) then


                           if (isour.le.ibrktm(i)) then

                              ts(i)=(tim(ibrktm(i))+tim(1))/2.d0
                              factor=43200.d0*(year_br(i)-year(1))
                              dt(irow) = (tim(isour)-ts(i))/factor

                              H2(IJ+3*nfull)=ISIG                   !  clock offset
                              H2(ij+4*nfull)=isig*dt(irow)          !  clock rate
                              H2(ij+5*nfull)=isig*dt(irow)**2       !  clock 2nd rate


		         H2 (IROW + (NR1-(3*i+5))*nfull) = 0.d0      !  clock offset
                         H2 (IROW + (NR1-(3*i+4))*nfull) = 0.d0      !  clock rate
                         H2 (IROW + (NR1-(3*i+3))*nfull) = 0.d0      !  clock 2nd rat



                           else if (isour.gt.ibrktm(i)) then

                            ts(i)=(tim(ibrktm(i))+tim(ntim))/2.d0
                            factor=43200.d0*(1.d0+year(1)-year_br(i))
                            dt(irow) = (tim(isour)-ts(i))/factor

                                  H2(IJ+3*nfull)=0.d0                   !  clock offset
                                  H2(ij+4*nfull)=0.d0                   !  clock rate
                                  H2(ij+5*nfull)=0.d0                   !  clock 2nd rate

                 H2 (IROW + (NR1-(3*i+5))*nfull) = ISIG                !  clock offset
                 H2 (IROW + (NR1-(3*i+4))*nfull) = isig*dt(irow)       !  clock rate
                 H2 (IROW + (NR1-(3*i+3))*nfull) = isig*dt(irow)**2    !  clock 2nd rat



                           end if

                        else                    !  dealing with problem station

                           ts0=(tim(ntim)+tim(1))/2.d0
                           dt(irow) = (tim(isour)-ts0)/43200.d0

                           H2(IJ+3*nfull)=ISIG                   !  clock offset
                           H2(ij+4*nfull)=isig*dt(irow)          !  clock rate
                           H2(ij+5*nfull)=isig*dt(irow)**2       !  clock 2nd rate

                        end if

                     end do

                  else                         !   No clock breaks

                     ts0=(tim(ntim)+tim(1))/2.d0
                     dt(irow) = (tim(isour)-ts0)/43200.d0

                     H2(IJ+3*nfull)=ISIG                   !  clock offset
                     H2(ij+4*nfull)=isig*dt(irow)          !  clock rate
                     H2(ij+5*nfull)=isig*dt(irow)**2       !  clock 2nd rate

                  end if

                  H2(IJ+6*nfull)=ISIG*PCOORD(K,6)  !  atmosphere           '1'
                  H2(IJ+7*nfull)=ISIG*PCOORD(K,7)  !  N-S atmosphere grad. '1'
                  H2(IJ+8*nfull)=ISIG*PCOORD(K,8)  !  E-W atmosphere grad. '1'

               END IF


            

            END DO


            if (num2(irow).ne.0) then

               ns1=irow + 2*(inn(irow)-1)*nfull    !  in -> inn(irow)   22.04.2016
               ns2=irow + (2*inn(irow)-1)*nfull    !  in -> inn(irow)   22.04.2016


                        
               c8(ns1) = psou(1)  ! partial derivatives
               c8(ns2) = psou(2)  ! partial derivatives

            end if

            p1(irow) = psou(1)
            
            p2(irow) = psou(2)

c
c
c     Revised on 11.07.2018 (OT) - epoch 2015.0 is used for reference
c

            dtq = (timdt - 57023.d0)/365.2422d0  !   2015.d0


            pa1 = -dsin(ra2000)                         !   A1
            pa2 =  dcos(ra2000)                         !   A2
            pa3 =  0.d0                                        !   A3

            pa4 = -dsin(de2000)*dcos(ra2000)            !   omega1
            pa5 = -dsin(de2000)*dsin(ra2000)            !   omega2
            pa6 =  dcos(de2000)                        !   omega3

c
c   Error for partial derivatives fixed on 19.01.2009
c
            pa7 =  dsin(de2000)*dcos(ra2000)            !   S21'
            pa8 = -dsin(de2000)*dsin(ra2000)            !   S21
            pa9 =  dcos(de2000)*dcos(2.d0*ra2000)       !   S22'
            pa10 = -dcos(de2000)*dsin(2.d0*ra2000)       !   S22   H12
            pa11 =  0.d0                                        !   S20



            pa12 = dsin(ra2000)*dcos(2.d0*de2000)        ! T21
            pa13 = dcos(ra2000)*dcos(2.d0*de2000)        ! T21'
c
c  Change to coeff 0.5 on 5-June-2009
c
            pa14 = - 0.5d0*dsin(2.d0*de2000)*dsin(2.d0*ra2000)            ! T22
            pa15 = - 0.5d0*dsin(2.d0*de2000)*dcos(2.d0*ra2000)            ! T22'
            pa16 =  dsin(2.d0*de2000)                    ! T20






            pd1 = -dsin(de2000)*dcos(ra2000)            !   A1
            pd2 = -dsin(de2000)*dsin(ra2000)            !   A2
            pd3 =  dcos(de2000)                        !   A3

            pd4 = -dsin(ra2000)                         !   omega1
            pd5 =  dcos(ra2000)                         !   omega2
            pd6 =  0.d0                                        !   omega3

            pd7 =  dsin(ra2000)*dcos(2.d0*de2000)       !   S21'
            pd8 =  dcos(ra2000)*dcos(2.d0*de2000)       !   S21
c
c  Change to coeff 0.5 on 5-June-2009
c
            pd9 = -0.5d0*dsin(2.d0*de2000)*dsin(2.d0*ra2000)  !   S22'
            pd10 = -0.5d0*dsin(2.d0*de2000)*dcos(2.d0*ra2000)  !   S22   H12
            pd11 =  dsin(2.d0*de2000)                   !   S20


c
c   Error for partial derivatives fixed on 19.01.2009
c
            pd12  = -dsin(de2000)*dcos(ra2000)           !   T21
            pd13 =  dsin(de2000)*dsin(ra2000)           !   T21'
            pd14 =  dcos(de2000)*dcos(2.d0*ra2000)      !   T22
            pd15 = -dcos(de2000)*dsin(2.d0*ra2000)      !   T22'
            pd16 =  0.d0                                       !   T20

c
c           Construction of vector spherical harmonic partials
c

ccc      For constant velocities

ccc            pg(irow,1)= psou(1)*pa1  + psou(2)*pd1     !  sin(d0)*cos(a0)
ccc            pg(irow,2)= psou(1)*pa2  + psou(2)*pd2     !  sin(d0)*sin(a0)
ccc            pg(irow,3)= psou(1)*pa3  + psou(2)*pd3     !  cos(d0)

ccc      For constant acceleration

            pg(irow,1)= dtq*(psou(1)*pa1  + psou(2)*pd1 )    !  sin(d0)*cos(a0)
            pg(irow,2)= dtq*(psou(1)*pa2  + psou(2)*pd2 )    !  sin(d0)*sin(a0)
            pg(irow,3)= dtq*(psou(1)*pa3  + psou(2)*pd3 )    !  cos(d0)

c
c     For vector direction to the solar centre
c

cc            bb = dsqrt(dotpr(b,b))   

cc            pg(irow,1)= (p(1) * dotpr(b,p) - b(1))/bb       !  
cc            pg(irow,2)= (p(2) * dotpr(b,p) - b(2))/bb      !  
cc            pg(irow,3)= (p(3) * dotpr(b,p) - b(3))/bb      !  




c
c       Rotation
c

            pg(irow,4)= dtq*(psou(1)*pa4  - psou(2)*pd4 )    !  T11
            pg(irow,5)= dtq*(psou(1)*pa5  - psou(2)*pd5 )    !  T11'
            pg(irow,6)= dtq*(psou(1)*pa6  - psou(2)*pd6 )    !  T10


            pg(irow,7) = dtq*(psou(1)*pa7  + psou(2)*pd7 )    !  S21
            pg(irow,8) = dtq*(psou(1)*pa8  + psou(2)*pd8 )    !  S21'
            pg(irow,9) = dtq*(psou(1)*pa9  + psou(2)*pd9 )    !  S22
            pg(irow,10)= dtq*(psou(1)*pa10 + psou(2)*pd10 )    !  S22'
            pg(irow,11)= dtq*(psou(1)*pa11 + psou(2)*pd11 )    !  S20

            pg(irow,12)= dtq*(psou(1)*pa12 - psou(2)*pd12 )    !   T21
            pg(irow,13)= dtq*(psou(1)*pa13 - psou(2)*pd13)    !   T21'
            pg(irow,14)= dtq*(psou(1)*pa14 - psou(2)*pd14)    !   T22
            pg(irow,15)= dtq*(psou(1)*pa15 - psou(2)*pd15)    !   T22'
            pg(irow,16)= dtq*(psou(1)*pa16 - psou(2)*pd16)    !   T20


            r2(IROW) = sgtot*(c*100.d0)**2   !   cm**2
 
         ENDIF

      ENDDO

      print *, ' isor = ' , isor


      do i=1,isor

         number(i) = ng(num(i))   ! amount of observations sorted on source
 

      end do


      if (idnm.ne.0) then
         do i=1,idnm
            print *, 'irow_break = ', irow_break(i), '  ts = ', ts(i)
         end do
      end if


      dt2 = year(nf)-year(1)
      tt = year(1)+dt2/2.d0

      nobs_unstab_sources = 0
      nobs_nounst_sources = 0
      nobs_stable_sources = 0

c      print *, 'num_unst = ', num_unst
      print *, 'dtq = ', dtq


  


c      print *, num_unst, num_sta, isor

      write (15,888) tt, isor, nf

      if (isor.gt.nsou) then
          print *, 'The number of sources exceeds NSOU'
          write (15,*) 'The number of sources exceeds NSOU', isor, nsou
          go to 999
      end if

      if (nr2.gt.npar) then
          print *, 'The number of parameters exceeds NPAR'
          write (15,*) 'The number of parameters exceeds NPAR',nr2,npar
          go to 999
      end if

 
      rewind 24


      

      imarc = 0
      iflag = 0
  
      DO IROW = 1,nf


         iarc (irow) =0

         do k=1,num_unst

            read (24,125) icrf107(k)

            if (fulsor(inn(irow)).eq.icrf107(k)) then    !  in -> inn(irow)

c
c   Limit of the number of observations of a single radio source to add the 
c   source to the list of reduced parameters.
c   Now the limit number is 3.
c   If the number is not sufficient, the second inversion may fail and
c   IER2 will be -1. Take care on that
c   (03-May-2016)
c                   
               if (number(inn(irow)).gt.3) then          !  in -> inn(irow)

                  iarc(irow) = 1

                  islg=0

                  do kk=1,num_unst
                      

                     if (fulsor(inn(irow)).eq.arcsour(kk)) then     !  in -> inn(irow)
                       

                        islg=1

                        km (irow) = kk



                        exit

                     end if

                  end do

                  if(islg.eq.0) then

                     narc=narc+1
                     imarc=imarc+1
                     arcsour(narc)=fulsor(inn(irow))            ! in -> inn(irow)

                     km (irow) =imarc
                  

                  end if

               end if

            end if


         end do

         rewind 24

      END DO  


    

      iglob = isor-narc   !  the number of global sources for the current session
      print *, 'num_unst = ', num_unst, '  num_sta  = ', num_sta
      print *, 'isor  = ', isor, ' narc = ', narc, ' iglob = ', iglob

      if (iflag.eq.1) write (15,443) iflag

443   format ('iflag = ', i4)


888   format (/'d',9x, 'NEW SESSION',1x,f12.5/
     *2x,' The number of radiosources is',i4/
     *2x,' The number of observations is',i10/)

890   format ('s',i4,2x,a8,2x,i4,2x,i5,2x,f5.2)



 
      DO IROW=1,nf

 
c
c    107 ICRF sources to be considered as 'arc' ones
c    (Fey, et al., Astonomical Journal, 2004)
c
c
c    163 ICRF sources to be considered as 'arc' ones
c    (Feissel, A&A, 2003)
c



         c8(irow+2*iglob*nfull)     = pg(irow,1)            ! for GAL
         c8(irow+(2*iglob+1)*nfull) = pg(irow,2)            ! for GAL
         c8(irow+(2*iglob+2)*nfull) = pg(irow,3)            ! for GAL
         c8(irow+(2*iglob+3)*nfull) = pg(irow,4)            ! for omega1
         c8(irow+(2*iglob+4)*nfull) = pg(irow,5)            ! for omega2
         c8(irow+(2*iglob+5)*nfull) = pg(irow,6)            ! for omega3
         c8(irow+(2*iglob+6)*nfull) = pg(irow,7)            ! for Grav  S21'
         c8(irow+(2*iglob+7)*nfull) = pg(irow,8)            ! for Grav  S21
         c8(irow+(2*iglob+8)*nfull) = pg(irow,9)            ! for Grav  S22'
         c8(irow+(2*iglob+9)*nfull) = pg(irow,10)           ! for Grav  S22
         c8(irow+(2*iglob+10)*nfull) = pg(irow,11)          ! for Grav S20
         c8(irow+(2*iglob+11)*nfull) = pg(irow,12)          ! for   T21
         c8(irow+(2*iglob+12)*nfull) = pg(irow,13)          ! for   T21'
         c8(irow+(2*iglob+13)*nfull) = pg(irow,14)          ! for   T22
         c8(irow+(2*iglob+14)*nfull) = pg(irow,15)          ! for   T22'
         c8(irow+(2*iglob+15)*nfull) = pg(irow,16)          ! for   T20

         c8(irow+(2*iglob+16)*nfull) = part_gamma(irow)     ! for gamma
         c8(irow+(2*iglob+17)*nfull) = part_parallax(irow)  ! for parallax 16.11.2023



         ns1 = irow + 2*(inn(irow)-1)*nfull         !  in -> inn(irow)
         ns2 = irow + (2*inn(irow)-1)*nfull         !  in -> inn(irow)

         do k=1,num_unst


            if (fulsor(inn(irow)).eq.icrf107(k)) then   !  in -> inn(irow)

               c8(irow+2*iglob*nfull)     = 0.d0          ! for GAL
               c8(irow+(2*iglob+1)*nfull) = 0.d0          ! for GAL
               c8(irow+(2*iglob+2)*nfull) = 0.d0          ! for GAL
               c8(irow+(2*iglob+3)*nfull) = 0.d0          ! omega1
               c8(irow+(2*iglob+4)*nfull) = 0.d0          ! omega2
               c8(irow+(2*iglob+5)*nfull) = 0.d0          ! omega3
               c8(irow+(2*iglob+6)*nfull) = 0.d0          ! for Grav
               c8(irow+(2*iglob+7)*nfull) = 0.d0          ! for Grav
               c8(irow+(2*iglob+8)*nfull) = 0.d0          ! for Grav
               c8(irow+(2*iglob+9)*nfull) = 0.d0          ! for Grav
               c8(irow+(2*iglob+10)*nfull) = 0.d0         ! for Grav
               c8(irow+(2*iglob+11)*nfull) = 0.d0         ! for   T21
               c8(irow+(2*iglob+12)*nfull) = 0.d0         ! for   T21'
               c8(irow+(2*iglob+13)*nfull) = 0.d0         ! for   T22
               c8(irow+(2*iglob+14)*nfull) = 0.d0         ! for   T22'
               c8(irow+(2*iglob+15)*nfull) = 0.d0         ! for   T20
               c8(irow+(2*iglob+16)*nfull) = 0.d0     ! for gamma
               c8(irow+(2*iglob+17)*nfull) = 0.d0     ! for parallax  16.11.2023

               c8 (ns1) = 0.d0  ! partial derivatives
               c8 (ns2) = 0.d0  ! partial derivatives


            end if

         end do

      ENDDO

      print *, ' OK after c8 '

      if (nsta.ne.2) then

         do i=1,3

            r2(nf+i) = 1.d0

            do j=1,iglob

               c8 (nf+i+2*(j-1)*nfull) = 0.d0
               c8 (nf+i+(2*j-1)*nfull) = 0.d0

            end do

         end do


         do i=4,6

            r2(nf+i) = 1.d0

            do j=1,iglob

              c8 (nf+i+2*(j-1)*nfull) = 0.d0
              c8 (nf+i+(2*j-1)*nfull) = 0.d0

            end do

         end do

         do i=1,nsta

            READ (20,REC=i) JST,STAN,RX,RY,RZ,EPOCH,VX,VY,VZ

            ik=kantenna(i)
            if (ik.ne.0) then

               rx=rx+vx*(year_1-epoch)
               ry=ry+vy*(year_1-epoch)
               rz=rz+vz*(year_1-epoch)

               rr = rx*rx + ry*ry + rz*rz

               if (ik.eq.1) then

          if (stan.eq.'METSAHOV'.or.stan.eq.'DSS65A  '.or.stan.eq.
     *   'TIGOCONC'.or.stan.eq.'KUNMING '.or
     * .stan.eq.'PARKES  '.or.stan.eq.'ISHIOKA '
c     * or.stan.eq.'ONSA13SW'.or.stan.eq.'ONSA13NE'
     * .or.stan.eq.'MACGO12M'.or.stan.eq.'RAEGSMAR'.or
     * .stan.eq.'GILCREEK'.or.stan.eq.'KOKEE12M'.or
     * .stan.eq.'KOGANEI '.or.stan.eq.'KASHIM11'.or
     * .stan.eq.'DSS36   '.or
     * .stan.eq.'TSUKUB32'.or.stan.eq.'URUMQI  '.or
     * .stan.eq.'KASHIM34')    then


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

                     H2(nf+1)=1.d0
                     H2(nf+2+nfull)=1.d0
                     H2(nf+3+2*nfull)=1.d0


                     H2(nf+4+nfull)=  -ah*rz/rr
                     H2(nf+4+2*nfull)=  ah*ry/rr

                     H2(nf+5)=  ah*rz/rr
                     H2(nf+5+2*nfull)=  -ah*rx/rr

                     H2(nf+6)=  -ah*ry/rr
                     H2(nf+6+nfull)=  ah*rx/rr

                  end if
               else

                  IJ1=nf+1 + (ik-2)*ns*nfull + 6*nfull
                  IJ2=nf+2 + (ik-2)*ns*nfull + 6*nfull
                  IJ3=nf+3 + (ik-2)*ns*nfull + 6*nfull
                  IJ4=nf+4 + (ik-2)*ns*nfull + 6*nfull
                  IJ5=nf+5 + (ik-2)*ns*nfull + 6*nfull
                  IJ6=nf+6 + (ik-2)*ns*nfull + 6*nfull


          if (stan.eq.'METSAHOV'.or.stan.eq.'DSS65A  '.or.stan.eq.
     *   'TIGOCONC'.or.stan.eq.'KUNMING '.or
     * .stan.eq.'PARKES  '.or.stan.eq.'ISHIOKA '
c     * or.stan.eq.'ONSA13SW'.or.stan.eq.'ONSA13NE'
     * .or.stan.eq.'MACGO12M'.or.stan.eq.'RAEGSMAR'.or
     * .stan.eq.'GILCREEK'.or.stan.eq.'KOKEE12M'.or
     * .stan.eq.'KOGANEI '.or.stan.eq.'KASHIM11'.or
     * .stan.eq.'DSS36   '.or
     * .stan.eq.'TSUKUB32'.or.stan.eq.'URUMQI  '.or
     * .stan.eq.'KASHIM34')                 then


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

      do i=1, nfull-3

         z2(i) = z2(i)*c*100.d0    !   sec -> cm

         

      end do

c
c     Moved here on 12-Jan-2016
c

      do i=1, nf

c
c +OT    18-Sep-2011
c

         do j=1,ibi


           if ((i_st1(i).eq.num_eli(j,1).and.i_st2(i).eq.num_eli(j,2))
     *     .or.(i_st1(i).eq.num_eli(j,2).and.i_st2(i).eq.num_eli(j,1)))
     *     then


              r2(i) = r2(i)*1.d6

           end if

         end do

c
c -OT    18-Sep-2011
c
                 

      end do

c
c    End of move on 12-Jan-2016
c

      print *, ' narc = ', narc

      nr2 = nr1 + 2*narc   !  Increasing amount of daily parameters due to arc sources

      print *, ' nr1 = ', nr1, ' nr2 = ', nr2, nfull


      DO irow = 1, nf

c
c   ARC sources as daily parameters
c
         if (narc.gt.0) then

            
            ns_arc1 = irow + nr1 * nfull + 2*(km(irow)-1)*nfull   !  RA for arc sources
            ns_arc2 = irow + nr1 * nfull + (2*km(irow)-1)*nfull   !  DE for arc sources

c            write (29,44) irow, km(irow), iarc(irow), ns_arc1, p1(irow)

c44          format (4(2x,i6),2x,f20.10)
   
            if(ns_arc1.ge.nobs1 * npar) iflag = 1         !  Something is wrong here

            if (iarc(irow).eq.1) then           !   enough observations to include the source

               h2(ns_arc1) = p1(irow)
               h2(ns_arc2) = p2(irow)

c            else

c               h2(ns_arc1)=0.d0
c               h2(ns_arc2)=0.d0


            end if

         end if

      ENDDO

c
c   call mnk1 is changed on 25.07.2007 (NNR added)
c                                       

c      do i=1, nfull


c         write (29,777) (c8(i+(j-1)*nfull),j=1,2*iglob + n14)

c         write (29,777) (h2(i+(j-1)*nfull),j=1,nr2), p1(i), p2(i)

c      end do
c777   format (220(1x,f10.6))

      schi = 0.d0
      sw = 0.d0

      if (nfull.ge.nobs1) go to 2999

      call mnk1 (nr1, nr2, nfull, year, r2, z2, h2, schi, sw, eps, nnr)
 

c
c     Downweighting of outliers after LSM
c
      do i=1, nf

         eps_res(i_st1(i))=eps_res(i_st1(i))+(eps(i)**2)/r2(i)
         eps_res(i_st2(i))=eps_res(i_st2(i))+(eps(i)**2)/r2(i)

         dd(i_st1(i)) = dd(i_st1(i)) + 1.d0/r2(i)
         dd(i_st2(i)) = dd(i_st2(i)) + 1.d0/r2(i)

         nk(i_st1(i)) = nk(i_st1(i))+1
         nk(i_st2(i)) = nk(i_st2(i))+1                               

      end do

      do i=1,nantenna
         nsta_rew(i)=0
         wrms(i) = dsqrt(eps_res(i)/dd(i))
         chik(i) = dsqrt(eps_res(i)/dble(nk(i)))
      end do

      print *, ' schi = ', schi, ' sw = ', sw
      n_i=0


c
c     Check for statistics - bad session do not corrupt full solution
c
c
c     More stringent rejection of bad data       !  1-April-2016
c
c     MOST OF SCHI and SW are less than 10!


      if ((schi.gt.10.and.sw.gt.10).or.iflag.eq.1) then

          write (128,*) timdt, schi, sw
                   
          do i=1, nantenna

             write (128,87) antenna(i), nk(i), chik(i), wrms(i)
             print 87, antenna(i), nk(i), chik(i), wrms(i)

          end do
         

          iflag = 2

          if (iflag.eq.2) go to 8777

           
          
          do i=1,nantenna

            write (12,88) antenna(i), nk(i), chik(i),wrms(i),nsta_rew(i)
c            print 88, antenna(i), nk(i), chik(i), wrms(i), nsta_rew(i)

          end do

      else

         print *, 'num_unst = ', num_unst, ' num_sta = ', num_sta, isor

c
c   Only good sessions are collected in 'INFORM.DAT'  - moved here on 15-Apr-2016
c

         do i=1,nantenna

            write (15,790) i, antenna(i), sgfac(i)

         end do
790      format (i3,1x,a8,2x,f16.6)

         write (15,*)

c         write (12,*)
c         write (12,*) year_1
 

         do i=1,isor

            w(i) = dble(ng(num(i)))/dble(nf)   !  Weights of constraints

c            if (i.le.5) print 890, i, fulsor(i)

            write (15,890) i, fulsor(i), ng(num(i)), num(i), w(i)

            sa(i)=dsin(ra(i))
            sd(i)=dsin(de(i))
            ca(i)=dcos(ra(i))
            cd(i)=dcos(de(i))

cc            number(i)=ng(num(i))   ! amount of observations sorted on source

            do k=1,num_unst

               if (fulsor(i).eq.icrf107(k)) then
                  nobs_unstab_sources = nobs_unstab_sources + number(i)
                  
               end if

            end do

            do k=1,num_sta

               
               if (fulsor(i).eq.icrf207(k)) then

                  
                  nobs_stable_sources = nobs_stable_sources + number(i)

c                  print *, k, number(i), nobs_stable_sources

               end if
            end do

         end do

         write (15,891) narc
891      format ('The number of arc sources',i4/)

   

      nobs_nounst_sources = nf-nobs_stable_sources-nobs_unstab_sources

      print *, 'number of obs of stable sources = ', nobs_stable_sources
      print *, 'number of obs of unstab sources = ', nobs_unstab_sources


         do i=1, nf

            if ((eps(i)**2/r2(i)).gt.tstfac*sw) then

               r2(i) = r2(i)*100000.d0
               n_i = n_i+1

               nsta_rew(i_st1(i)) = nsta_rew(i_st1(i))+1
               nsta_rew(i_st2(i)) = nsta_rew(i_st2(i))+1

            end if

         end do

         write (30,330) timdt, schi, sw, iy1, mo, idd

 330     format (f11.5,2x,f8.4,2x,f8.4,2x,i2.2,a3,i2.2)      

         do i=1,nantenna

            write (30,88) antenna(i), nk(i), chik(i),wrms(i),nsta_rew(i)
            print 88, antenna(i), nk(i), chik(i), wrms(i), nsta_rew(i)

         end do

 87   format (1x,a8,4x,i4,2(3xf8.4))
 88   format (1x,a8,4x,i4,2(3x,f8.4),3x,i4,' observations downweighted')

         print *, n_i, ' observations have been downweighted '


         write (12,444) isor, narc, iglob

444   format (' The number of sources in the session is ',i3,';'2x,i4,'
     *    - Arc',2x,i4,' - Global')

         

         if (iglob.lt.2) print *, 555
555   format (' Only one global radio source in this session ! ')



         if (iglob.lt.2) go to 8777

      
c
c     Parameter for NNR CRF constraint  (01.01.2015, Good start for 2015 year!)
c
c     LOOSE Constraint!! 3-June-2016
c     ON ALL Sources!!   3-June-2016
 

        constr_nnr_crf = 1.d8    !  26-July-2024   

c         constr_nnr_crf = 1.d6   !  01.01.2015   100 --> 1000000 on 03-Jun-2016

c         constr_nnr_crf = 1.d2   !   

c          constr_nnr_crf = 1.d-2   !   0.01 on 13-Jan-2025


         print *, ' LOOSE constraint for CRF NNR = ', constr_nnr_crf
         write (12,*) 'constraint for CRF NNR = ', constr_nnr_crf

c
c
c   Change on 25.07.2007
c
         if (NNR.eq.'NNRYES') then

            do irow=nfull-2, nfull

               z2(irow)=0.d0
               R2(IROW) = constr_nnr_crf     !    10000 on 27-Mar-2015

              do j=1,nr2

                  ij=irow+(j-1)*nfull
                  h2(ij)=0.d0

              end do

           end do

c
c    NNR-Constraints for CRF (on Gipson)
c
            do j=1,iglob  !  The number of global radiosources in session - 21-Jan-2013


c
c    Constraints are imposed for limited set of sources
c
c
c  Constraint on ALL Sources!!!   3-June-2016
c
c
c
ccccc               do k=1,num_sta  !  The number of sources fixed for constraints
   
ccccc                   if (fulsor(j).eq.icrf207(k)) then

c  Weghted NNR-constraints  (not used, but can be used)

c                  c8 (nfull-2+2*(j-1)*nfull) = -sd(j)*ca(j)*w(j)
c                  c8 (nfull-2+(2*j-1)*nfull) = sa(j)*w(j)
c                  c8 (nfull-1+2*(j-1)*nfull) = -sd(j)*sa(j)*w(j)
c                  c8 (nfull-1+(2*j-1)*nfull) = -ca(j)*w(j)
c                  c8 (nfull  +2*(j-1)*nfull) = cd(j)*w(j)

c  Unweghted NNR-constraints ! 20.08.2004

                     c8 (nfull-2+2*(j-1)*nfull) = -sd(j)*ca(j)
                     c8 (nfull-2+(2*j-1)*nfull) = sa(j)
                     c8 (nfull-1+2*(j-1)*nfull) = -sd(j)*sa(j)
                     c8 (nfull-1+(2*j-1)*nfull) = -ca(j)
                     c8 (nfull  +2*(j-1)*nfull) = cd(j)
                     c8 (nfull  +(2*j-1)*nfull) = 0.d0

c
c   CONSTRAINT on ALL Sources!  3-June-2016
c

ccccc                  end if

ccccc               end do

            end do
         END IF



c
c    End of change on 25.07.2007
c



c    No NNR for ICRF

c      do irow=nfull-2, nfull
c         z2(irow)=0.d0
c         R2(IROW) = 1000.d0
c         do j=1,nr2
c            ij=irow+(j-1)*nfull
c            h2(ij)=0.d0
c         end do
c      end do


         idj=int(timdt)
         CALL ddate (idj, iy, mon, iday)

         write (28,299) timdt, schi, sw, nfull, n_i, iy, mon, iday-1,
     *nobs_stable_sources, nobs_nounst_sources, nobs_unstab_sources,
     * iy1, mo, idd, dbc_code

299   format (f10.4,1x,'schi = ',f5.2,1x,'sw = ',f5.2,3(1x,i5),2(i3),
     *3(1x,i5),2x,i2.2,a3,i2.2,a2)



         CALL mnk2 (nr1, nr2, nfull, iglob,
     *         R2, Z2, H2, c8, nnr, cm8, v8, ier2)
      
         print *, 'Ier2 = ', ier2

c
c     Calculation for GLOBAL ADJUSTMENT
c
c

         print *, ' iglob = ', iglob, ' nss = ', nss, isor


         print *, ' total number of estimabales ', 2*nss + n14

         i_crit = 0

          do i = 2*iglob+4, 2*iglob+n14 -2

             if (dabs(v8(i)).ge.1.d6) then 

                 i_crit = -1

                 write (12,*) '  Warning  ', tt, timdt, i, v8(i)

             end if             
         
          end do   
       

         if (ier2.ne.-1.and.i_crit.ne.-1) then


         do i=1,iglob      !  iglob - the number of global radiosources in NGS file

C
C  ATTENTION !!   isor vs iglob
C

ccc         do i=1,isor         !  isor - full number of radio sources in NGS file

            i1=2*num(i)-1
            i2=2*num(i)

            ii1=2*i-1
            ii2=2*i

c            print *, ' i = ', i, i1, num(i)
c            print *, ' OK ', ii1, ii2
c            print *, ' num = ', num(i), 2*num(i)

cc            read (27,102,rec=i1) vv(i1)
cc            read (27,102,rec=i2) vv(i2)

            read (27,102,rec=i1) vvi1
            read (27,102,rec=i2) vvi2

c        write (12,*) i, i1, i2, ii1, ii2, num(i), 
c     *         fulsor(i), vvi1, vvi2, v8(ii1), v8(ii2)


            vvi1 = vvi1 + v8(ii1)
            vvi2 = vvi2 + v8(ii2)



c      end if

            write (27,102,rec=i1) vvi1
            write (27,102,rec=i2) vvi2



            do j=1,i

               j1=2*num(j)-1
               j2=2*num(j)
               jj1=2*j-1
               jj2=2*j

               if(i1.ge.j1) then
                  ij1=i1*(i1-1)/2+j1
                  read (26,101,rec=ij1) cc11
                  iijj1=ii1*(ii1-1)/2+jj1
                  cc11 =cc11 +cm8(iijj1)
                  write (26,101,rec=ij1) cc11
               end if

               if(i1.lt.j1) then
                  ij1=j1*(j1-1)/2+i1
                  read (26,101,rec=ij1) cc12
                  iijj1=ii1*(ii1-1)/2+jj1
                  cc12 = cc12+cm8(iijj1)
                  write (26,101,rec=ij1) cc12
               end if


               if(i2.ge.j2) then
                  ij2=i2*(i2-1)/2+j2
                  read (26,101,rec=ij2) cc22
                  iijj2=ii2*(ii2-1)/2+jj2
                  cc22 = cc22+cm8(iijj2)
                  write (26,101,rec=ij2) cc22
               end if
               if(i2.lt.j2) then
                  ij2=j2*(j2-1)/2+i2
                  read (26,101,rec=ij2) cc21
                  iijj2=ii2*(ii2-1)/2+jj2
                  cc21 = cc21 +cm8(iijj2)
                  write (26,101,rec=ij2) cc21
               end if


               if(i2.ge.j1)  then
                  ij3=i2*(i2-1)/2+j1
                  read (26,101,rec=ij3) cc33
                  iijj3=ii2*(ii2-1)/2+jj1
                  cc33=cc33+cm8(iijj3)
                  write (26,101,rec=ij3) cc33
               end if
               if(i2.lt.j1)  then
                  ij3=j1*(j1-1)/2+i2
                  read (26,101,rec=ij3) cc34
                  iijj3=ii2*(ii2-1)/2+jj1
                  cc34 = cc34+cm8(iijj3)
                  write (26,101,rec=ij3) cc34
               end if

               if (i.ne.j) then

                  if(i1.ge.j2) then
                     ij4=i1*(i1-1)/2+j2   !  non-diagonal elements of correlation matrix
                     read (26,101,rec=ij4) cc43
                     iijj4=ii1*(ii1-1)/2+jj2     !  non-diagonal elements of correlation matrix
                     cc43=cc43+cm8(iijj4)
                     write (26,101,rec=ij4) cc43
                  end if

                  if(i1.lt.j2) then
                     ij4=j2*(j2-1)/2+i1   !  non-diagonal elements of correlation matrix
                     read (26,101,rec=ij4) cc44
                     iijj4=ii1*(ii1-1)/2+jj2     !  non-diagonal elements of correlation matrix
                     cc44=cc44+cm8(iijj4)
                     write (26,101,rec=ij4) cc44
                  end if

               end if
            end do
         end do

c
c   14.05.2003    ! for the part for Galaxy rotation
c   19.07.2006      and for gamma
c   01.03.2007      and for grav waves
c   04.10.2007      and for octople
c   08.10.2007      and for T2
c   16.11.2023      and for parallax
c

         print *, ' first part of recording is fine '

         do i=1,n14   ! n14 = 18 = 3+3+5+5+1+1   (Gal + Rotation + GW + gamma + parallax)

            ir1 = 2 * nss + i  

            read (27,102,rec=ir1) vvir1

c
c    isor --> iglob  21-Jan-2013
c
c

            vvir1 = vvir1 + v8(2*iglob+i)   ! 2 * iglob + 18

          
            write (27,102, rec=ir1) vvir1
c
c    isor --> iglob  21-Jan-2013
c
c

            do j=1,iglob

               j1=2*num(j)-1
               j2=2*num(j)
               jj1=2*j-1
               jj2=2*j

               ij1=(2*nss+i)*(2*nss+i-1)/2+j1
               ij2=(2*nss+i)*(2*nss+i-1)/2+j2

c
c    isor --> iglob  21-Jan-2013
c
c


               iijj1=(2*iglob+i)*(2*iglob+i-1)/2+jj1
               iijj2=(2*iglob+i)*(2*iglob+i-1)/2+jj2




               read (26,101,rec=ij1) cc1
               read (26,101,rec=ij2) cc2

               cc1 = cc1 + cm8(iijj1)
               cc2 = cc2 + cm8(iijj2)

               write (26,101,rec=ij1) cc1
               write (26,101,rec=ij2) cc2
            end do

           do j=1,i                       !   06.06.2003

c
c    isor --> iglob  21-Jan-2013
c
c

               ij=(2*nss+i)*(2*nss+i-1)/2+2*nss+j

               iijj=(2*iglob+i)*(2*iglob+i-1)/2 + 2*iglob+j


               read (26,101,rec=ij) cc3
               cc3 = cc3 + cm8(iijj)
               write (26,101,rec=ij) cc3

            end do

         end do

c
c  END of BAD SESSION REJECTION
c
      end if

      end if

 101  format (d24.16)
 102  format (d24.16)
 103  format (300(1x,d18.10))
 104  format (d24.16)
 105  format (14(1x,i7))

      close (15)
      close (26)
      close (27)

      RETURN

9991  WRITE (*,*)  '***** ERROR ***** num_sta is more than n_stable '

      RETURN

9992  WRITE (*,*)  '***** ERROR ***** num_unst is more than n_unstable '

      RETURN

2999   WRITE (*,*)  '** the number of observations exceed the limit **'

999   WRITE (*,*) ' ERROR '

8777  continue

      RETURN

      END


C**********************************************************************
C
C    SUBROUTINE KZERO
C
C    THIS ROUTINE INITIALIZES THE VARIABLES USED BY PREPARE
C
C   REVISION MARCH 17, 1998
C   LAST REVISION May 20, 2002, (OT) Limits for observations increased to
C                                    20000 (IDICTIO, IOBS)
C   REVISION June 11, 2002, (OT) NFULL is calculating here
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
     *               IBAS,ISTM,ITIM,ISTA
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
      nfull = nobs

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
      DIMENSION IANTENNA(nmax),KANTENNA(nmax),TIM(nobs1)
      DIMENSION IDICTIO(nobserv),NBASE(ntimmax),IOBS(nobserv)

      COMMON /DICTION/ TIM,IDICTIO,IOBS
      COMMON /CINDEX/NBAS,NOBS,NSTM,NTIM,NSTA,NSEL,
     *               IBAS,ISTM,ITIM,ISTA
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
c         IBLN(I,1) = KST1
c         IBLN(I,2) = KST2
         IF (KST1.NE.KST2) THEN      !       21.01.2001 (O.Titov)

C  FILL IDICTIO
C  Massive IDICTIO is filling up by observational indexes
C

            DO J=1,ntimmax     !  NTIMMAX in DTAU1.FOR   !!
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
C
C    THIS ROUTINE READS INFORMATION OF OBSERVABLES, SIGMAS AND
C      PARTIAL DERIVATIVES FROM OCCAM STANDARD DATA FILES
C
C   REVISION MARCH 17, 1998
C   REVISION MAY 13, 2002  (OT)   Calculation of EOP partial
C                  derivatives for delay rates has been added
C   REVISION May 20, 2002, (OT) Limits for observations increased to
C                                    20000 (IDICTIO, IOBS)
C   REVISION June 11,2002, (OT) 4000 observations; 'AMB' has cut
C   REVISION March 22, 2006, (OT) Cable correction implemented
C   REVISION June 27, 2006, (OT) Partial derivatives for gamma implemented
C   REVISION 2006 SEPTEMBER 18 (O.T.) Sign for thermal correction was
C   changed due to corresponding change in r'STATION.FOR'
C   LAST REVISION 2010 JUNE 25 (O.T.) Large update in Auckland
C
C**********************************************************************

      SUBROUTINE EXTRACT (IRECORD, CABL1, CABL2,
     * ra, de, CMODEL, part_gamma, b, p, part_parallax, teta)

C   ESPECIFICATIONS

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'OCCAM.FI'

      character*3 cmodel
      DIMENSION PCOORD(2,8),PCOORR(2,8),AXKT1(2),AXKT2(2)
      DIMENSION CRX1(3),CRY1(3),CRZ1(3),CGLONG1(3),CPHI1(3),CELHGT1(3)
      DIMENSION CRX2(3),CRY2(3),CRZ2(3),CGLONG2(3),CPHI2(3),CELHGT2(3)
      DIMENSION cmniel1_w(2), cmniel1_d(2), a_grad1(2)
      DIMENSION cmniel2_w(2), cmniel2_d(2), a_grad2(2)
      DIMENSION PEOP1(3,3),PEOP2(3,3),PPOL(3),PPOL1(3)
      DIMENSION ppolar_rate1(3), ppolar_rate2(3)
      DIMENSION PNUT1(2),PNUT2(2),PNUT(2),PSOU(2),drda(3),drdd(3),v(3) ! 25.06.2010 Auckland
      dimension ve(3), ww1(3), ww2(3), p(3)       !  25.06.10 Auckland


      dimension b(3), azimut1(3), azimut2(3)

      COMMON /CINDEX/NBAS,NOBS,NSTM,NTIM,NSTA,NSEL,
     *               IBAS,ISTM,ITIM,ISTA
      COMMON /CEXTRACT/ OC,SOBD,OCR,SOBR,PCOORD,PCOORR,
     * PNUT,PPOL,PSOU,PPOL1
      COMMON /PHYS/ C, FL, AH, AU, J20
      COMMON /MATH/ PI, TWOPI, ROG, ROS, ROZ, TOL
      COMMON /OPT/ SGADD, ZEN_DIST, SGFAC(nmax), tstfac

C   READ ON BASTIM THE RECORD "IRECORD", READ FROM IDICTIO
C   IF CONTAINS THE THEORETICAL AND OBSERVED DELAY AND THE
C   SIGMA OF THE LATTER

c      open (1, file = 'test1.dat', status = 'unknown')
c      open (2, file = 'test2.dat', status = 'unknown')

c      READ (17,rec=IRECORD) JBAS,ISOUR,IST1,IST2,D1,SD1,R01,SR1,
c     *        DI1,SDI1,RI1,SRI1,TP1,T0,TM1,TGP1,TG0,TGM1,part_gamma,     ! 27.06.2006

      READ (17,rec=IRECORD) JBAS,ISOUR,IST1,IST2,D1,SD1,R01,SR1,
     *        DI1,SDI1,RI1,SRI1,TP1,T0,TM1,TGravP1,TGrav0,TGravM1,
     *        part_gamma, teta, cos_A, bb, phia, AA, freq0,
     *        DIFFA, RSUN
     *        ,ve(1), ve(2), ve(3), part_parallax         ! 16 November 2023




C   READ THE STATION INFORMATION AT THE TIME FOR BOTH ANTENNAS
C   THIS INCLUDES THE PARTIAL DERIVATIVES, THE ATMOSPHERIC MODELS
C   AND THE AXIS OFFSETS

      READ (18,REC=IST1) JSTA,ISTA,ISOUR,TEMP,PRES,HUMD,IFAC,CABL1,
     *                 CRX1,CRY1,CRZ1,RX1,RY1,RZ1,
     *                 CGLONG1,CPHI1,CELHGT1,Z1,HA1,
     *                 HSTA1,Z2,HA2,HSTA2,Z31,HA3,HSTA3,
     *                 ZDRY1,cmniel1_w,cmniel1_d,
     *                 w11, w12, w13,              !  geodetic velocity of first station
     *                 a11, a12, a13,              !  geodetic acceleration
     *                 AXKT1,
     *                 (PCOORD(1,I),I=1,5),(PCOORR(1,I),I=1,5),
     *                 PNUT1,PEOP1,
     *                 ppolar_rate1,
     *                 therm_d1,azimut1,a_grad1,
     *            dr1_ut1,                      ! 01-01-2021
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
     

      READ (18,REC=IST2) JSTA,ISTB,ISOUR,TEMP,PRES,HUMD,IFAC,CABL2,
     *                 CRX2,CRY2,CRZ2,RX2,RY2,RZ2,
     *                 CGLONG2,CPHI2,CELHGT2,Z1,HA1,
     *                 HSTA1,Z2,HA2,HSTA2,Z32,HA3,HSTA3,
     *                 ZDRY2,cmniel2_w,cmniel2_d,
     *                 w21, w22, w23,             !  geodetic velocity of second station
     *                 a21, a22, a23,             !  geodetic acceleration
     *                 AXKT2,
     *                 (PCOORD(2,I),I=1,5),(PCOORR(2,I),I=1,5),
     *                 PNUT2,PEOP2,
     *                 ppolar_rate2,
     *                 therm_d2,azimut2,a_grad2,
     *            dr2_ut1,                          !  01-01-2021  
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

      b(1) = 0.d0
      b(2) = 0.d0
      b(3) = 0.d0



      dva=0.d0
      dvd=0.d0
      denom = 0.d0


      ww2(1) = w21
      ww2(2) = w22
      ww2(3) = w23

      DO I=1,3

c         print *, ww2(i)

         DVA = DVA + (V(i) + ww2(i)) *Drda(I)
         DVD = DVD + (V(i) + ww2(i)) *drdd(i)

         denom = denom + p(i)*(v(i)+ww2(i))/c

      
      ENDDO

      B(1) = CRX2(3) - CRX1(3)
      B(2) = CRY2(3) - CRY1(3)
      B(3) = CRZ2(3) - CRZ1(3)


      denom = 1.d0 + denom

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
      OBR     = R01*1.D-12
      SOBR    = SR1*1.D-12

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

      

      IF (CMODEL.EQ.'NIE') THEN
         TROPO = (ZDRY2*cmniel2_d(1)-ZDRY1*cmniel1_d(1))/C            ! sec
         TROPOR= (ZDRY2*cmniel2_d(2)-ZDRY1*cmniel1_d(2))*1.D09/C      ! sec
         PCOORD(1,6) = cmniel1_w(1)
         PCOORD(2,6) = cmniel2_w(1)
      ENDIF

      IF (CMODEL.EQ.'VMF') THEN
         TROPO = (ZDRY2*dqmfh2-ZDRY1*dqmfh1)/C                        ! sec
         TROPOR= (ZDRY2*cmniel2_d(2)-ZDRY1*cmniel1_d(2))*1.D09/C      ! sec
         PCOORD(1,6) = dqmfw1
         PCOORD(2,6) = dqmfw2
      ENDIF

      IF (CMODEL.EQ.'VM1') THEN
         TROPO = (ZDRY2*vmf1h2-ZDRY1*vmf1h1)/C                        ! sec
         TROPOR= (ZDRY2*cmniel2_d(2)-ZDRY1*cmniel1_d(2))*1.D09/C      ! sec
         PCOORD(1,6) = vmf1w1
         PCOORD(2,6) = vmf1w2
      ENDIF


      

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
c  For troposphere gradient partial derivatives
c

      PCOORD(1,7) = a_grad1(1)*pcoord(1,6)      !  Troposphere
      PCOORD(2,7) = a_grad2(1)*pcoord(2,6)      !        Gradients
      PCOORD(1,8) = a_grad1(2)*pcoord(1,6)      !            partial
      PCOORD(2,8) = a_grad2(2)*pcoord(2,6)      !                derivatives


C   GET AXIS OFFSET CORRECTIONS

      AXIS  = AXKT2(1)-AXKT1(1)
      AXISR = AXKT2(2)-AXKT1(2)

C   GET THERMIC DEFORMATION CORRECTIONS

      therm = therm_d1-therm_d2   !  Sign was changed on 18.09.2006 (OT)

C   ADD FALSE CLOCK MODEL TO THE DATA FOR TESTING

C   COMPUTE THEORETICAL DELAY ADDING THE GEOMETRICAL MODEL, THE
C  TROPOSPHERIC CONTRIBUTION AND THE AXIS OFFSET CORRECTION


      THDELAY = T0 + TROPO + AXIS + therm    !  to be checked for all types
ccc                                             !  of antennas

      THRATE = (TP1-TM1)/2.D0 +TROPOR+AXISR

C   GET OBSERVABLE MINUS CALCULUS

      OC = OBD - THDELAY                     !  seconds

      OCR= OBR - THRATE

C   GET PARTIAL DERIVATIVES WITH RESPECT TO NUTATION

c      PNUT (1) =  (PNUT2(1)-PNUT1(1))/C
c      PNUT (2) =  (PNUT2(2)-PNUT1(2))/C

      PNUT (1) =  (PNUT2(1)-PNUT1(1))*100.d0   ! 'cm'
      PNUT (2) =  (PNUT2(2)-PNUT1(2))*100.d0   ! 'cm'

C   GET PARTIAL DERIVATIVES WITH RESPECT TO POLE MOTION AND UT1

c      PPOL (1) = (PEOP2(1,3)-PEOP1(1,3))/C        !  sec
c      PPOL (2) = (PEOP2(2,3)-PEOP1(2,3))/C        !  sec
c      PPOL (3) = (PEOP2(3,3)-PEOP1(3,3))/C        !  sec


      PPOL (1) =  (PEOP2(1,3)-PEOP1(1,3))*100.d0           !  cm
      PPOL (2) =  (PEOP2(2,3)-PEOP1(2,3))*100.d0           !  cm
      PPOL (3) =  (PEOP2(3,3)-PEOP1(3,3))*100.d0           !  cm

C   GET PARTIAL DERIVATIVES on delay rates  WITH RESPECT TO EOP   13.05.2002

      PX1  = (PEOP1(1,2)-PEOP1(1,1))/2        !     X - comp  for 1st station
      PY1  = (PEOP1(2,2)-PEOP1(2,1))/2        !     Y - comp  for 1st station
      PUT1  = (PEOP1(3,2)-PEOP1(3,1))/2       !     UT1-UTC   for 1st station

      PX2  = (PEOP2(1,2)-PEOP2(1,1))/2        !     X - comp  for 2st station
      PY2  = (PEOP2(2,2)-PEOP2(2,1))/2        !     Y - comp  for 2st station
      PUT2  = (PEOP2(3,2)-PEOP2(3,1))/2       !     UT1-UTC   for 2st station

      PPOL1 (1) =  (PX2-PX1)*100.d0      !  cm
      PPOL1 (2) =  (PY2-PY1)*100.d0      !  cm
      PPOL1 (3) =  (PUT2-PUT1)*100.d0    !  cm

c
c   Update for partials on 25-June-2010  in Auckland
c
c
C   GET PARTIAL DERIVATIVES WITH RESPECT TO SOURCE COORDINATES

      PSOU (1) =  (PCOORD(2,4)-PCOORD(1,4))/(c*denom)     ! [m]-> [seconds]
      PSOU (2) =  (PCOORD(2,5)-PCOORD(1,5))/(c*denom)     ! [m]-> [seconds]

      second_term_da =  sc_t*dva/(denom)**2
      second_term_dd =  sc_t*dvd/(denom)**2

c      print *, psou(1), second_term_da

      PSOU (1) =  PSOU(1) + second_term_da     
      PSOU (2) =  PSOU(2) + second_term_dd     

c      print *, psou(1), psou(2)

c
c   End update of partials on 25-June-2010 in Auckland
c

      Return
      END



C**********************************************************************
c
c      SUBROUTINE CLOCK_BREAKS (IDNM, IBRKST, IBRKTM)
c
c
c    Clock breaks reading (OT)  - 24.11.2004
c    Last update  (OT)  - 25.07.2005
c
C**********************************************************************

      SUBROUTINE CLOCK_BREAKS (IDNM, IBRKST, IBRKTM)
      integer*4 ibrkst(10), ibrktm(10)

      OPEN (178,FILE='breaks.res',STATUS='OLD',ERR=99)

      IDNM = 0
322   CONTINUE
         READ (178,310,END=330) IDST,IDTM
310      FORMAT (I3,1X,I4)
         if (idst.eq.0.or.idtm.eq.0) then
            print *, ' WARNING  !!!  Check the BREAKS.RES file !!'
         else
            IDNM = IDNM+1
            IBRKST(IDNM) = IDST
            IBRKTM(IDNM) = IDTM
         end if
         if (idst.ne.0.and.idtm.ne.0) GO TO 322
330   CONTINUE

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

      CHARACTER*8 ST1(nmax), ST2(nmax)           ! nbase -> nmax on 03.04.2014

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

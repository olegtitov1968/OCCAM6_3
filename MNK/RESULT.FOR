C**********************************************************************
C
C    SUBROUTINE RESULT
C
C    THIS ROUTINE COMPUTES THE FINAL PARAMETER VALUES AND THEIR
C      STANDARD DEVIATIONS
C    IT WRITES ON OUTPUT FILE RESULT.KAL THIS INFORMATION, AS WELL
C      AS THE NORMALIZED VARIANCE-COVARIANCE MATRIX (IF OPTION
C      SELECTED)
C    A FILE WITH INFORMATION FOR A MULTIEXPERIMENT ANALYSIS IS ALSO
C      CREATED HERE
C
C    
C   LAST REVISION: 01 April 2021 (OT)  Clock break to SINEX file
C
C   Copyright (C) 2021 Oleg Titov, Geoscience Australia, 2001-2021
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
C**********************************************************************

      SUBROUTINE RESULT (RX, RY, RZ, x00, xm1, xp1, y00, ym1, yp1,
     * ut10, utm1, utp1, norm, vect, est, ngr2,
     * ntot, varf, sw, sq_sum, ut1_1, idat_1, ut1_nf,
     * idat_nf, arcsour, narc, ra_narc, de_narc)

C   ESPECIFICATIONS

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'OCCAM.FI'

      CHARACTER ANTENNA(nmax)*8,CMODEL*3,SOURCE*8,stan(nmax)*8,
     * dom*9,st_des*22, lin1*28
      CHARACTER mo*3, filen1*3, filen2*24, dbc_code*2, filename*31
      CHARACTER*4 ant(nmax),axtyp, stc,code(nmax)

      CHARACTER*8 nut_mod, prec_mod
      CHARACTER*8 arcsour(narc)
      character*8 s_ivs, s_iers, ss_iers(7970)
      character*8 s_icrf, ss_icrf(7970)

      DIMENSION IANTENNA(nmax),KANTENNA(nmax), dolg(nmax), shir(nmax)

      dimension sinex(npar*npar), svec(npar)        !    26.11.2010
      dimension est(npar)                           !    02.04.2012

      real*8 mjd
      real*8 norm(npar*(npar+1)/2), vect(npar)
      real*8 rx(nmax), ry(nmax), rz(nmax)
	real*8 rx0(nantenna), ry0(nantenna), rz0(nantenna)
	real*8 vx0(nantenna), vy0(nantenna), vz0(nantenna)
      dimension ra_narc(n_unstable), de_narc(n_unstable)

      integer*4 icode(narc)

      integer*4 nums(npar)                      !    25.11.2010

      COMMON /CANTENNA/ANTENNA,IANTENNA,KANTENNA,NANTENNA,IANT_NST
      COMMON /PHYS/ C, FL, AH, AU, J20
      COMMON /MATH/ PI, TWOPI, ROG, ROS, ROZ, TOL

      istan = 1


C   SEND REPORT ON OUTLIERS FOUND AND IF THEY WERE CORRECTED OR NOT

c      print *, ' ut1_1 = ', ut1_1, ' idat_1 = ', idat_1
c      print *, ' ut1_nf = ', ut1_nf, ' idat_nf = ', idat_nf

      open (11,file='dbc_code',status='unknown')
c

      read (11,96) iy, mo, idd, dbc_code                   !  04.05.2005
96    format (i2,a3,i2,a2)
      write (*,96) iy, mo, idd, dbc_code

      if (dbc_code(2:2).eq.' ')  dbc_code(2:2)='_'

      filen1(1:3)=mo
      filen2(1:2)=dbc_code
      filen2(3:24)='_aus2020a.snx'

      WRITE (filename,97) IY, filen1, IDD, filen2
 97   FORMAT (I2.2,a3,I2.2,a24)

      open (554,file=filename, status='unknown')
      open (55,file='icrf2',status='unknown') 
      open (65,file='eop_apriori.txt',status='unknown')


C   READ TIME AND WRITE IN OUTPUT FILE FOR KMEXP


      CALL JULDAT (IY,ID,IH,IM,S,IDAT_1,UT1_1,1)
      IFAC = 365
      IF (MOD(IY,4).EQ.0) IFAC=366
      YEAR = IY + (ID+IH/24.+IM/1440.+S/86400.)/IFAC

      sp= IH/24.+ IM/1440. + s/86400.
      iss=int(86400.0*sp)

      CALL JULDAT (IY_1,ID_1,IH_1,IM_1,S_1,IDAT_nf,UT1_nf,1)

      iss_1=int(86400.0*(IH_1/24.+IM_1/1440.+s_1/86400.))

      if(iy.lt.2000) iyy=iy-1900
      if(iy.ge.2000) iyy=iy-2000
      if(iy_1.lt.2000) iyy_1=iy_1-1900
      if(iy_1.ge.2000) iyy_1=iy_1-2000

      ih_1_24=ih_1+24

      if (ih+ih_1_24.ge.48) then
         idat_m = idat_1+1
         ut_m=(ut1_1+ut1_nf+twopi)/2.d0-twopi
      else
         idat_m = idat_1                                        	
         ut_m=(ut1_1+ut1_nf+twopi)/2.d0
      end if

      CALL JULDAT (IY_m,ID_m,IH_m,IM_m,S_m,IDAT_m,UT_m,1)
      iss_m=int(86400.0*(IH_m/24.+IM_m/1440.+s_m/86400.))
      if(iy_m.lt.2000) iyy_m=iy_m-1900
      if(iy_m.ge.2000) iyy_m=iy_m-2000


      n_pars_to_be = nantenna * 3 + 8 + 2 * narc 

c
c     The Number of pre-reduced parameters must be added as Sabine Bachmann adviced
c     10-Jan-2017
c     So, there are 6 parameters for each antenna but the reference one
c     and three  parameters for the reference antenna!

      n_pars_to_be_plus = n_pars_to_be + (nantenna - 1) * 6 + 3


c
c     SINEX OUTPUT
c
      in=2

      write (554,233) IYY,ID,iss,IYY_1,ID_1,iss_1, n_pars_to_be, in


233   format ('%=SNX 2.00 AUS ','00:000:00000',1x,'AUS',1x,i2.2,':',
     *i3.3,':',i5.5,1x,i2.2,':',i3.3,':',i5.5,' R ',i5.5,1x,i1,' S',
     *' E',' C')

      write (554,234)
234   format ('*======================================================')
      write (554,235)
235   format ('+FILE/REFERENCE')
      write (554,236)
236   format (1x,'DESCRIPTION       GEOSCIENCE AUSTRALIA  - VLBI')
      write (554,237)
237   format (1x,'OUTPUT            DAILY STATION COORDINATES   ')
      write (554,238)
238   format (1x,'CONTACT           oleg.titov@ga.gov.au ')
      write (554,239)
239   format (1x,'SOFTWARE          OCCAM v. 6.3  ')
      write (554,240)
240   format (1x,'HARDWARE          PENTIUM-III  ')
      write (554,2350)
2350  format ('-FILE/REFERENCE')
      write (554,234)

      write (554,241)
241   format ('+SITE/ID')
      write (554,242)
242   format ('*CODE PC   DOMES   T  Station Description       LONG
     *       LAT       HEIGHT')

      DO I=1,NANTENNA

         IST = IANTENNA(I)
         ANT(I)=ANTENNA(I)(1:4)

         READ (20,REC=IST) JST,STAN(i),RX0(i),RY0(i),RZ0(i),EPOCH,
     *    VX0(i),VY0(i),VZ0(i),axtyp,of,stc,phi, glong

         call snx (stc,code(i),dom,st_des,lod,lom,slo,lad,lam,sla, hei)

         write (554,243) code(i),dom,st_des,lod,lom,slo,lad,lam,sla,hei
243      format (1x,a4,2x,'A',1x,a9,' R',1x,a22,1x,i3,1x,i2,1x,f4.1,
     *    1x,i3,1x,i2,1x,f4.1,1x,f7.1)

         dolg(i)=glong
         shir(i)=phi

      END DO
      write (554,2431)
2431  format ('-SITE/ID')

      write (554,234)
      write (554,2432)
2432  format ('+SOLUTION/EPOCHS')

      isol=1

      DO I=1,NANTENNA
         write (554,2433) code(i),isol,IYY,ID,iss,IYY_1,ID_1,
     *    iss_1, iyy_m, id_m, iss_m
      END DO

2433  format (1x,a4,2x,'A',4x,i1,1x,'R',3(1x,i2.2,':',i3.3,':',i5.5))
      write (554,2434)
2434   format ('-SOLUTION/EPOCHS')
      write (554,234)

      nut_mod = 'iau2000a'
      prec_mod = 'iers1996'


      write (554,2435) nut_mod
2435  format ('*',1x,a8,1x,'MODEL IAU2000A',56x)

      write (554,2436) prec_mod
2436  format ('*',1x,a8,1x,'MODEL IERS1996',56x)

      write (554,24361)
24361 format ('+SOURCE/ID ')

      write (554,24362)
24362 format (' Code IERS des ICRF designation Comments ')


      do i=1,7970

cc         read (55,998) s_ivs, s_icrf, s_iers

         read (55,998) s_ivs

         do j=1,narc

            if (arcsour(j).eq.s_ivs) then

               ss_icrf(j) = s_ivs

               ss_iers(j) = s_ivs

cc               if(s_iers(1:1).ne.'-') then
cc
cc                 ss_iers(j) = s_iers
cc
cc               else
cc
cc                  ss_iers(j) = arcsour(j)
cc
cc               end if
            
            end if


         end do

      end do

ccc  998   format (a8,20x,a10,2x,a8)
998   format (3x,a8)    



      do i=1,narc

         icode(i)=i
         write (554,2437) icode(i), ss_iers(i), ss_icrf(i)

      end do
2437  format (1x,i4,1x,a8,1x,a8,68x)

      write (554,24363)
24363 format ('-SOURCE/ID ')


      write (554,244)
244   format ('+SOLUTION/APRIORI ')
      write (554,245)
245   format ('*INDEX TYPE',19x,'EPOCH',4x,'UNITS',9x,'ESTIMATES',10x,
     *'STD_DEV')

      s1=0.d0
      s0=0.d0

      ic=1
      in=2

      DO I=1,NANTENNA

C    GET A PRIORI COORDINATES OF STATIONS


c
c        Writing of a priori parameters
c

         CALL NO_BREAK (i, RX, RY, RZ, code, isol, iyy_m,
     *   id_m, iss_m, ic, in, s1, ii)


      END DO

      xx=x00*1000.d0
      yy=y00*1000.d0
      ut=ut10*1000.d0

      xa1 = (xp1-xm1)*1000.d0                !  Apriori for X-dot
      ya1 = (yp1-ym1)*1000.d0                !  Apriori for Y-dot
      ua1 = -(utp1-utm1)*1000.d0             !  Apriori for LOD


      write (554,264) ii+3, isol, iyy_m, id_m, iss_m, in, s0, s1
      write (554,265) ii+4, isol, iyy_m, id_m, iss_m, in, s0, s1
      write (554,266) ii+5, isol, iyy_m, id_m, iss_m, in, xx, s1
      write (554,267) ii+6, isol, iyy_m, id_m, iss_m, in, yy, s1
      write (554,268) ii+7, isol, iyy_m, id_m, iss_m, in, ut, s1
      write (554,269) ii+8, isol, iyy_m, id_m, iss_m, in, xa1, s1
      write (554,270) ii+9, isol, iyy_m, id_m, iss_m, in, ya1, s1
      write (554,271) ii+10, isol, iyy_m, id_m, iss_m, in, ua1, s1

266     format (1x,i5,1x,'XPO   ',' ____  A',1x,i4,1x,i2.2,':',i3.3,':',
     *   i5.5,1x,'mas',2x,i1,1x,e21.15,1x,e11.6)
267     format (1x,i5,1x,'YPO   ',' ____  A',1x,i4,1x,i2.2,':',i3.3,':',
     *   i5.5,1x,'mas',2x,i1,1x,e21.15,1x,e11.6)
268     format (1x,i5,1x,'UT    ',' ____  A',1x,i4,1x,i2.2,':',i3.3,':',
     *   i5.5,1x,'ms ',2x,i1,1x,e21.15,1x,e11.6)
269     format (1x,i5,1x,'XPOR  ',' ____  A',1x,i4,1x,i2.2,':',i3.3,':',
     *   i5.5,1x,'masD',1x,i1,1x,e21.15,1x,e11.6)
270     format (1x,i5,1x,'YPOR  ',' ____  A',1x,i4,1x,i2.2,':',i3.3,':',
     *   i5.5,1x,'masD',1x,i1,1x,e21.15,1x,e11.6)
271     format (1x,i5,1x,'LOD   ',' ____  A',1x,i4,1x,i2.2,':',i3.3,':',
     *   i5.5,1x,'ms ',2x,i1,1x,e21.15,1x,e11.6)


      do i=1,narc

         write (554,2438) ii+11+2*(i-1), icode(i),
     *    isol, iyy_m, id_m, iss_m, in, ra_narc(i), s1
         write (554,2439) ii+11+2*i-1, icode(i),
     *    isol, iyy_m, id_m, iss_m, in, de_narc(i), s1

      end do

2438  format (1x,i5,1x,'RS_RA ',1x,i4,2x,'A',1x,i4,1x,i2.2,':',i3.3,':',
     *   i5.5,1x,'rad',2x,i1,1x,e21.15,1x,e11.6)
2439  format (1x,i5,1x,'RS_DE ',1x,i4,2x,'A',1x,i4,1x,i2.2,':',i3.3,':',
     *   i5.5,1x,'rad',2x,i1,1x,e21.15,1x,e11.6)


c
c     Number of parameters to be estimated
c


      write (554,2440)
2440  format ('-SOLUTION/APRIORI ')

      write (554,234)
      write (554,2441)
2441  format ('+SOLUTION/STATISTICS ')
      write (554,2442) float(ntot)
2442  format (1x,'NUMBER OF OBSERVATIONS       ',1x,e22.15)
      write (554,2443) float(n_pars_to_be_plus)
2443  format (1x,'NUMBER OF UNKNOWNS           ',1x,e22.15)
      write (554,2444) float(ntot - n_pars_to_be_plus)
2444  format (1x,'NUMBER OF DEGREES OF FREEDOM ',1x,e22.15)
      write (554,2445) varf**2
2445  format (1x,'VARIANCE FACTOR              ',1x,e22.15)
      write (554,24451) sw * c
24451 format (1x,'SQUARE SUM OF RESIDUALS (VTVP)',1x,e22.15)
      write (554,24452) sq_sum * c 
24452 format (1x,'WEIGHTED SQUARE SUM OF O-C   ',1x,e22.15)

      write (554,2446)
2446  format ('-SOLUTION/STATISTICS ')
      write (554,234)

      FRAD = 180.D0*3600.D0/PI


264   format (1x,i5,1x,'NUT_X ',' ____  A',1x,i4,1x,i2.2,':',i3.3,':',
     *i5.5,1x,'mas',2x,i1,1x,e21.15,1x,e11.6)
265   format (1x,i5,1x,'NUT_Y ',' ____  A',1x,i4,1x,i2.2,':',i3.3,':',
     *i5.5,1x,'mas',2x,i1,1x,e21.15,1x,e11.6)

c      print *, ' Apriori values have been written to SINEX file '

c      do i=1,ngr2
c         print *, vect(i)
c      end do


      if (istan.eq.1) then

         write (554,234)

C   WRITE COVARIANCE IF SELECTED

         ik=0
         ic=1
         in=2

         write (554,2741)
2741     format ('+SOLUTION/NORMAL_EQUATION_VECTOR')
         write (554,2739)
2739     format ('*INDEX TYPE',19x,'EPOCH',4x,'UNITS',9x,'ESTIMATES')



         DO I=1,NANTENNA

c
c        Writing of a priori parameters
c

            CALL NO_BREAK_VEC (i, vect, code, isol, iyy_m,

     *          id_m, iss_m, ic, in, ii)


         END DO

      write (554,2640) ii+3, isol, iyy_m, id_m, iss_m, in, vect(ii+3)
      write (554,2650) ii+4,isol, iyy_m, id_m, iss_m, in, vect(ii+4)
      write (554,2660) ii+5,isol, iyy_m, id_m, iss_m, in, vect(ii+5)
      write (554,2670) ii+6,isol, iyy_m, id_m, iss_m, in, vect(ii+6)
      write (554,2680) ii+7,isol, iyy_m, id_m, iss_m, in, vect(ii+7)
      write (554,2690) ii+8,isol, iyy_m, id_m, iss_m, in, vect(ii+8)
      write (554,2700) ii+9,isol, iyy_m, id_m, iss_m, in, vect(ii+9)
      write (554,2710) ii+10,isol, iyy_m, id_m, iss_m, in, vect(ii+10)

2640   format (1x,i5,1x,'NUT_X ',' ____  A',1x,i4,1x,i2.2,':',i3.3,':',
     * i5.5,1x,'mas',2x,i1,1x,e21.15)
2650   format (1x,i5,1x,'NUT_Y ',' ____  A',1x,i4,1x,i2.2,':',i3.3,':',
     * i5.5,1x,'mas',2x,i1,1x,e21.15)
2660   format (1x,i5,1x,'XPO   ',' ____  A',1x,i4,1x,i2.2,':',i3.3,':',
     * i5.5,1x,'mas',2x,i1,1x,e21.15)
2670   format (1x,i5,1x,'YPO   ',' ____  A',1x,i4,1x,i2.2,':',i3.3,':',
     * i5.5,1x,'mas',2x,i1,1x,e21.15)
2680   format (1x,i5,1x,'UT    ',' ____  A',1x,i4,1x,i2.2,':',i3.3,':',
     * i5.5,1x,'ms ',2x,i1,1x,e21.15)
2690   format (1x,i5,1x,'XPOR  ',' ____  A',1x,i4,1x,i2.2,':',i3.3,':',
     * i5.5,1x,'masD',1x,i1,1x,e21.15)
2700   format (1x,i5,1x,'YPOR  ',' ____  A',1x,i4,1x,i2.2,':',i3.3,':',
     * i5.5,1x,'masD',1x,i1,1x,e21.15)
2710   format (1x,i5,1x,'LOD   ',' ____  A',1x,i4,1x,i2.2,':',i3.3,':',
     * i5.5,1x,'ms ',2x,i1,1x,e21.15)


      do i=1,narc

         write (554,2711) ii+11+2*(i-1), icode(i),
     *    isol, iyy_m, id_m, iss_m, in, vect(ii+11+2*(i-1))
         write (554,2712) ii+11+2*i-1, icode(i),
     *    isol, iyy_m, id_m, iss_m, in, vect(ii+11+2*i-1)
         
c         print *, i, ii, vect(ii+11+2*(i-1))

      end do

2711  format (1x,i5,1x,'RS_RA ',1x,i4,2x,'A',1x,i4,1x,i2.2,':',i3.3,':',
     *   i5.5,1x,'rad',2x,i1,1x,e21.15,1x,e11.6)
2712  format (1x,i5,1x,'RS_DE ',1x,i4,2x,'A',1x,i4,1x,i2.2,':',i3.3,':',
     *   i5.5,1x,'rad',2x,i1,1x,e21.15,1x,e11.6)


         write (554,2742)
2742     format ('-SOLUTION/NORMAL_EQUATION_VECTOR ')

         write (554,234)

         write (554,275)
275      format ('+SOLUTION/NORMAL_EQUATION_MATRIX U')

c         print *, ' ngr2 = ', ngr2, ' narc = ', narc


c         write (77,*)

         do i=1, ngr2

cc         do i=1,ngr2 + 2*narc

c            print *, i, norm(i*(i+1)/2)

c             print *, i, vect(ii+11+2*i-1)


c              write (77,*) (norm(i*(i-1)/2+j),j=1,i)

cc            do j=i,ngr2 + 2*narc,3

            do j=i,ngr2,3

               l=2

cc               if ((ngr2+2*narc-j).lt.3) l = ngr2+2*narc-j

               if ((ngr2-j).lt.3) l = ngr2 - j

               write (554,276) i, j, (norm((j+k)*((j+k)-1)/2+i),k=0,l)

            end do

         end do

276      format (1x,i5,1x,i5,300(1x,e21.14))

         write (554,277)
277      format ('-SOLUTION/NORMAL_EQUATION_MATRIX U')

         write (554,234)


         ik=0
         ic=1
         in=2

c         write (554,2751)
c2751     format ('+SOLUTION/ESTIMATE ')
c
c         DO I=1,NANTENNA
c
c
c        Writing of a priori parameters
c
c
c            CALL NO_BREAK_VEC (i, est, code, isol, iyy_m,
c     *          id_m, iss_m, ic, in, ii)
c
c
c         END DO
c
c      write (554,2640) ii+3, isol, iyy_m, id_m, iss_m, in, est(ii+3)
c      write (554,2650) ii+4,isol, iyy_m, id_m, iss_m, in, est(ii+4)
c      write (554,2660) ii+5,isol, iyy_m, id_m, iss_m, in, est(ii+5)
c      write (554,2670) ii+6,isol, iyy_m, id_m, iss_m, in, est(ii+6)
c      write (554,2680) ii+7,isol, iyy_m, id_m, iss_m, in, est(ii+7)
c      write (554,2690) ii+8,isol, iyy_m, id_m, iss_m, in, est(ii+8)
c      write (554,2700) ii+9,isol, iyy_m, id_m, iss_m, in, est(ii+9)
c      write (554,2710) ii+10,isol, iyy_m, id_m, iss_m, in, est(ii+10)
c
c
c      do i=1,narc
c         write (554,2711) ii+11+2*(i-1), icode(i),
c     *    isol, iyy_m, id_m, iss_m, in, est(ii+11+2*(i-1))
c         write (554,2712) ii+11+2*i-1, icode(i),
c     *    isol, iyy_m, id_m, iss_m, in, est(ii+11+2*i-1)
c      end do
c
c
c         write (554,2752)
c2752     format ('-SOLUTION/ESTIMATE ')

         write (554,234)
         write (554,278)
278      format ('%ENDSNX')

      else !  (istan.eq.0)

          print *, ' ERROR - STATION POSITIONS ARE NOT ESTIMATED !! '

      end if

      END




      SUBROUTINE NO_BREAK (i, RX, RY, RZ, code, isol, iyy_m, id_m,
     * iss_m, ic, in, s1, ii)


      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'OCCAM.FI'

      real*8 rx(nmax), ry(nmax), rz(nmax)
      CHARACTER*4 ant(nmax),axtyp, stc,code(nmax)


      ik = 0

      if (i.eq.1) then
        write (554,251) ik+1,code(i),isol,iyy_m,id_m,iss_m,ic,RX(i),s1
        write (554,252) ik+2,code(i),isol,iyy_m,id_m,iss_m,ic,RY(i),s1
        write (554,253) ik+3,code(i),isol,iyy_m,id_m,iss_m,ic,RZ(i),s1

      else

        ii=(ns-6)*(i-2)+4

        write (554,251) ii  ,code(i),isol,iyy_m,id_m,iss_m,ic,RX(i),s1
        write (554,252) ii+1,code(i),isol,iyy_m,id_m,iss_m,ic,RY(i),s1
        write (554,253) ii+2,code(i),isol,iyy_m,id_m,iss_m,ic,RZ(i),s1
      end if

251   format (1x,i5,1x,'STAX  ',1x,a4,1x,' A',1x,i4,1x,i2.2,':',i3.3,
     *':',i5.5,1x,'m',4x,i1,1x,e21.15,1x,e11.6)
252   format (1x,i5,1x,'STAY  ',1x,a4,1x,' A',1x,i4,1x,i2.2,':',i3.3,
     *':',i5.5,1x,'m',4x,i1,1x,e21.15,1x,e11.6)
253   format (1x,i5,1x,'STAZ  ',1x,a4,1x,' A',1x,i4,1x,i2.2,':',i3.3,
     *':',i5.5,1x,'m',4x,i1,1x,e21.15,1x,e11.6)

      RETURN
      END








      SUBROUTINE NO_BREAK_VEC (i, vect, code, isol, iyy_m, id_m,
     * iss_m, ic, in, ii)


      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'OCCAM.FI'

      real*8 vect(npar)
      CHARACTER*4 ant(nmax),axtyp, stc,code(nmax)


      ik = 0

      if (i.eq.1) then
         write (554,2610) ik+1,code(i),isol,iyy_m,id_m,iss_m,ic,vect(1)
         write (554,2620) ik+2,code(i),isol,iyy_m,id_m,iss_m,ic,vect(2)
         write (554,2630) ik+3,code(i),isol,iyy_m,id_m,iss_m,ic,vect(3)
      else

         ii = 3*(i-2)+4

        write (554,2610) ii  ,code(i),isol,iyy_m,id_m,iss_m,ic,vect(ii)
        write (554,2620)ii+1,code(i),isol,iyy_m,id_m,iss_m,ic,vect(ii+1)
        write (554,2630)ii+2,code(i),isol,iyy_m,id_m,iss_m,ic,vect(ii+2)
       end if

2610  format (1x,i5,1x,'STAX  ',1x,a4,1x,' A',1x,i4,1x,i2.2,':',i3.3,
     *':',i5.5,1x,'m',4x,i1,1x,e21.15)
2620  format (1x,i5,1x,'STAY  ',1x,a4,1x,' A',1x,i4,1x,i2.2,':',i3.3,
     *':',i5.5,1x,'m',4x,i1,1x,e21.15)
2630  format (1x,i5,1x,'STAZ  ',1x,a4,1x,' A',1x,i4,1x,i2.2,':',i3.3,
     *':',i5.5,1x,'m',4x,i1,1x,e21.15)


      RETURN
      END




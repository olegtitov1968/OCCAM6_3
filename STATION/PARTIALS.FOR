      SUBROUTINE PARTIA  (RQU,HSTA,RBAR,VBAR,RX,RY,RZ,SN,
     *      PREC,SS,ALPHA,DELTA,EROT,STAT,SIDGRE_dot,dsdxw_t,
     * dsdyw_t,dsdut_t, ppolar_rate)

C    SUBROUTINE PARTIA
C
C    THIS SUBROUTINE COMPUTE THE PARTIAL DERIVATIVES OF DELAYS AND
C  RATES WITH RESPECT TO THE PARAMETERS
C
C  REVISION  1990 NOVEMBER 28  (N.Z.L.)
C  REVISION 2000  July 14 (OT)  option for tides has been changed
C  Revision 2001 April 17 (O.T.)  Remark by Johannes Boehm
C  Revision 2002 MAY 13 (O.T.), Partial derivatives for EOP rates
C                         Change in common-block /PART/, /PDATA/
C  Revision 2006 MARCH 24 (O.T.), Changes for OCCAM 6.2 - All new partials
C  Revision 2007 FEBRUARY 02 (O.T.), HSTA array added, partials for sources
C                                         updated.
C
c   Last Revision: 12 April 2021 (OT)
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
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      LOGICAL OPTTD,OPTOC,OPTAT,OPTPT,OPTLOV
      CHARACTER*8 STAT
      DIMENSION PCOORD(5),PCOORR(5)
      DIMENSION CGLONG(3),CPHI(3),CELHGT(3),CRX(3),CRY(3),CRZ(3)
      DIMENSION HA(3),Z(3),RQU(3),VBAR(3,3), DE(3), HSTA(3)
      DIMENSION RS(3), v(3), SS(3,3), EROT(3,3), V2(3)
      DIMENSION RBAR(3,3), dist(3)
      DIMENSION dsdut_t(3),dsdxw_t(3),dsdyw_t(3), ppolar_rate(3)

      DIMENSION SN(3,3), PREC(3,3)                                !   24-Mar-2014

      
      COMMON /PART / PCOORD,PCOORR,PNUTAT(2),PPOLAR(3,3)
      COMMON /PDATA/ SD,CD,SIDGRE,DKDPSI(3),DKDEPS(3),
     * DSDXW(3,3),DSDYW(3,3), DSDUT(3,3), DKDPSI_v(3),DKDEPS_v(3)
      COMMON /COCDEL/ CGLONG,CPHI,CELHGT,CRX,CRY,CRZ,HA,Z,DE
      COMMON /OPT/ OPTTD,OPTOC,OPTAT,OPTPT,OPTLOV
      COMMON /PHYS/ C, FL, AH, AU, J20
c      COMMON /BAR1/ RBAR0(3)
      COMMON /PARM/ OMEGA          !   added on 02.02.2007
      COMMON /GRAV/ GMS            !   added on 02.02.2007

C   BEWARE! THE MODULE TO OBTAIN THE PARTIAL DERIVATIVES WITH RESPECT
C  TO THE SOURCES MUST BE CHANGED TO A CLEARER FORM

C     COMPUTE THE PARTIALS WITH RESPECT TO THE STATION COORDINATES
C  (GEOMETRICAL MODEL)

      CINV = 1.D0/C
      vrqu=0.d0

c      open (1, file = 'partials.dat', status = 'unknown')

      do i=1,3

         v(i)=vbar(i,3)
         dist(i)=rbar(i,3)

         v2(i)=0.d0
      end do

      V2(1) = - OMEGA * RY
      V2(2) =   OMEGA * RX
      V2(3) =   0.D0

      VEV2=  DOTPR(v,v2)                     !  02.02.20007


      do i=1,3

         vrqu = vrqu + RQU(I)*(v(i) + v2(i))

      end do

      RS(1)=rx                       !
      RS(2)=ry                       !   Coordinates before rotation!
      RS(3)=rz                       !

c      print *, rx, crx(3)


      vmod=dsqrt(dotpr(v,v))
      RSUN= DSQRT( DOTPR(dist,dist) )
      VB= DOTPR(v,rs)/(c*c)

      drp  = 1.d0 + cinv*vrqu
      drp2 = 1.d0 + 0.5d0*cinv*vrqu



c      cof=1.d0-2.d0*GMS/(RSUN*c*c)-vmod**2/(2.d0*c*c)             !
      cof=1.d0-2.d0*GMS/(RSUN*c*c)-vmod**2/(2.d0*c*c)-vev2/(c*c)   !  02.02.2007


      CALL PSTATI (SS,EROT,RQU,VBAR,cof,drp,drp2,STAT,SIDGRE,sidgre_dot)

C     COMPUTE THE PARTIALS WITH RESPECT TO SOURCE COORDINATES
C  (GEOMETRICAL MODEL)

c+ OT  04.06.2010

      CALL PSOUR2 (RQU,SS,CRX,CRY,CRZ,V,VB,v2,ALPHA,DELTA,
     * SIDGRE,SIDGRE_dot,cof,drp)

      pcoord(4) = pcoord(4)/dcos(delta)
      pcoorr(4) = pcoorr(4)/dcos(delta)

c- OT  04.06.2010

C     COMPUTE THE PARTIALS WITH RESPECT TO NUTATION

      CALL PNUT (CRX,CRY,CRZ,cof,drp)

C     COMPUTE THE PARTIALS WITH RESPECT TO EOP

      CALL PEOP (RQU,SN,PREC,cof,drp,vb, sidgre_dot, dsdxw_t,
     * dsdyw_t,dsdut_t,  ppolar_rate)

C     COMPUTE THE PARTIALS WITH RESPECT TO LOVE NUMBERS

ccc   IF (OPTLOV)     CALL PLOVE

      RETURN
      END

************************************************************************

      SUBROUTINE PSTATI (SS,EROT,RQU,VBAR,cof,drp,drp2,STAT,
     * SIDGRE,sidgre_dot)

C
C     SUBROUTINE PSTATI
C
C    THIS ROUTINE COMPUTES THE PARTIAL DERIVATIVES OF DELAY AND RATE
C  WITH RESPECT TO STATION COORDINATES
C
C  REVISION  1993 JUNE 30 (NZL)  CORRECTED PARTIALS RESPECT TO
C                                      STATIONS (POLAR MOTION INC.)
C  Revision 2002 MAY 13 (O.T.), Partial derivatives for EOP rates
C                         Change in common-block /PART/, /PDATA/
C  Revision 2006 MARCH 24 (O.T.), Changes for OCCAM 6.2
C  Last revision 2006 August 28 (O.T.) /COCDEL/   is off
C
C-----------------------------------------------------------------------

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      CHARACTER*8 STAT
      DIMENSION PCOORD(5),PCOORR(5)
      DIMENSION RQU(3),VBAR(3,3),RV(3),EROT(3,3),ss(3,3)
      dimension p1(3), p2(3), p3(3), v1(3), v2(3), v3(3), p2_dot(3)
      dimension p3_dot(3)
      COMMON /PART / PCOORD,PCOORR,PNUTAT(2),PPOLAR(3,3)
      COMMON /PARM/ OMEGA
      COMMON /PHYS/ C, FL, AH, AU, J20


      STH = DSIN(SIDGRE)
      CTH = DCOS(SIDGRE)



C
C   COMPUTE THE TERMS for VECTOR of PARTIALS
C


      do i=1,3

         p1(i)=0.d0
         v1(i)=0.d0

         do j=1,3

            p1(i)=p1(i)+rqu(j)*ss(j,i)            !    s*P*N
            v1(i)=v1(i)+vbar(j,3)*ss(j,i)/c       !    V_e*P*N

         end do

      end do


      p2(1)= p1(1)*CTH+p1(2)*STH                !    s*P*N*R3(-H)
      p2(2)=-p1(1)*sth+p1(2)*cth                !
      p2(3)= p1(3)                              !


c
c   for delay rate
c
      p2_dot(1)= (-p1(1)*STH + p1(2)*CTH) * sidgre_dot              !    s*P*N*dR3(-H)
      p2_dot(2)= (-p1(1)*CTH - p1(2)*STH) * sidgre_dot              !
      p2_dot(3)= 0.d0                           !


      v2(1)= v1(1)*CTH+v1(2)*STH                !    v*P*N*R3(-H)
      v2(2)=-v1(1)*STH+v1(2)*CTH               !
      v2(3)= v1(3)                              !

      do i=1,3

         p3(i)=0.d0
         p3_dot(i)=0.d0
         v3(i)=0.d0

         do j=1,3

            p3(i) = p3(i) + p2(j)*erot(j,i)            !    s*P*N*R3(-H)*R1(y)*R2(x)
            p3_dot(i) = p3_dot(i) + p2_dot(j)*erot(j,i)    !    s*P*N*dR3(-H)*R1(y)*R2(x)
            v3(i) = v3(i) + v2(j)*erot(j,i)            !    v*P*N*R3(-H)*R1(y)*R2(x)

        end do

      end do


C   COMPUTE THE PARTIALS

     
      PCOORD(1) = -(p3(1)*cof + v3(1) * (1.d0 + drp2)) / drp     ! ddx
      PCOORD(2) = -(p3(2)*cof + v3(2) * (1.d0 + drp2)) / drp     ! ddy
      PCOORD(3) = -(p3(3)*cof + v3(3) * (1.d0 + drp2)) / drp     ! ddz
    

      
      PCOORR(1) = -p3_dot(1)*cof / drp     ! ddx for delay rate
      PCOORR(2) = -p3_dot(2)*cof / drp     ! ddy for delay rate
      PCOORR(3) = -p3_dot(3)*cof / drp     ! ddz for delay rate
     

    
c      print *, ' cof = ',cof, ' drp = ', drp
c      write (15,*) p3(1), p3(2), p3(3)
c      write (15,*) v3(1), v3(2), v3(3)
c      print *, '  new ', pcoord(1), pcoord(2), pcoord(3)
c      print *, pcoorr(1), pcoorr(2), pcoorr(3)

      RETURN
      END

***********************************************************************

      SUBROUTINE PSOUR2 (RQU,SS,CRX,CRY,CRZ,V,VB,v2,ALPHA,DELTA,
     * SIDGRE,SIDGRE_dot,cof,drp)

C
C    SUBROUTINE PSOUR2
C
C    THIS SUBROUTINE COMPUTES THE PARTIAL DERIVATIVES WITH RESPECT
C  TO SOURCE COORDINATES PCOORD(4), (5)
C
C   NO PARTIALS WILL BE CALCULATED FOR THE DELAY RATE AT THIS MOMENT
C
C    IN ORDER TO ESTIMATE THE COORDINATES AT J2000.0, THE PARTIAL
C  DERIVATIVES MUST TAKE ADVANTAGE OF THE PRECESION AND NUTATION
C  MATRICES COMPUTED IN PN
C
C    THIS IS A FIRST CODE TO INCORPORATE THE PARTIALS WITH RESPECT TO
C  J2000 COORDINATES
C
C  REVISION 1990 NOVEMBER 19  (N.Z.L.)
C  Revision 2002 MAY 13 (O.T.), Partial derivatives for EOP rates
C                         Change in common-block /PART/, /PDATA/
C  Revision 2002 June 19 (O.T.), Rotation matrix R3 moved to STATION.FOR
C  Revision 2006 MARCH 24 (O.T.), Changes for OCCAM 6.2
C  Revision 2006 August 28 (O.T.) /COCDEL/   is off
C                        CRX,CRY,CRZ - from SUBROUTINE line
C  Revision 2006 August 28 (O.T.) VB array added, partials updated
C  Last revision 2010 June 04 (O.T.) Many changes
C
C-----------------------------------------------------------------------
C

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION V(3),E(3,3),ERK(3),RQU(3)
      DIMENSION DRDA(3),DRDD(3),DTDA(3),DTDD(3),P_TRF(3), P_TRF1(3)
      DIMENSION CRX(3), CRY(3), CRZ(3), SS(3,3), V2(3)

      DIMENSION ERK_dot(3), p_trf_dot(3)  !   03-Sep-2019


      COMMON /PART / PCOORD(5),PCOORR(5)
      COMMON /PARM/ OMEGA
      COMMON /PHYS/ C, FL, AH, AU, J20
      COMMON /GRAV/ GMS, GMM, GME, GM_Mer, GM_Ven, GM_Mar, GM_Jup,
     * GM_Sat, GM_Ura, GM_Nep

C
C    PARTIAL DERIVATIVES FOR THE SOURCE POSITIONS
C    THESE PARTIALS ARE WITH RESPECT TO ARC INSTEAD OF ANGLE
C
C   THE VALUES PCOOR?(4) ARE THE DERIVATIVES WITH RESPECT TO RA
C  DIVIDED BY C

      STH = DSIN(SIDGRE)
      CTH = DCOS(SIDGRE)

c
c     ERK = R3(-H)*CRX = R3(-H)*R1(y)*R2(x)*r_st
c


      ERK(1) = CRX(3)*CTH - CRY(3)*STH
      ERK(2) = CRX(3)*STH + CRY(3)*CTH
      ERK(3) = CRZ(3)

c
c     Partial on ERK to calculate the partial on source for delay rate
c
c     03-Sep-2019
c

      ERK_dot(1) = (-CRX(3)*STH - CRY(3) * CTH) * sidgre_dot
      ERK_dot(2) = ( CRX(3)*CTH - CRY(3) * STH) * sidgre_dot
      ERK_dot(3) =  0.d0



      do i=1,3

         p_trf(i)=0.d0
         p_trf_dot (i) = 0.d0

         do j=1,3

            p_trf(i) = p_trf(i) + ss(i,j)*erk(j)    !  P*N*R3(-H)*R1(y)*R2(x)*r
            p_trf_dot(i) = p_trf_dot(i) + ss (i,j) * erk_dot(j)

         end do

      end do

      DRDA(1) = -DSIN(ALPHA)*DCOS(DELTA)   !  ds/da
      DRDA(2) =  DCOS(ALPHA)*DCOS(DELTA)   !
      DRDA(3) =  0.D0                      !

      DRDD(1) = -DCOS(ALPHA)*DSIN(DELTA)   !  ds/dd
      DRDD(2) = -DSIN(ALPHA)*DSIN(DELTA)   !
      DRDD(3) =  DCOS(DELTA)               !


c
c   Added on 02.02.2007
c

      DO I=1,3

         DVA = DVA + V(i)*Drda(I)
         DVD = DVD + V(i)*Drdd(I)

      ENDDO

      DVA=DVA/c
      DVD=DVD/c


      BS = cof*dotpr(p_trf,rqu)

      pcoord(4) = 0.d0
      pcoord(5) = 0.d0

      pcoorr(4) = 0.d0
      pcoorr(5) = 0.d0

      do i=1,3

         pcoord(4) = pcoord(4) + drda(i)*p_trf(i)   !   [m]
         pcoord(5) = pcoord(5) + drdd(i)*p_trf(i)   !   [m]

         pcoorr(4) = pcoorr(4) + drda(i) * p_trf_dot(i)  ! [m/s]
         pcoorr(5) = pcoorr(5) + drdd(i) * p_trf_dot(i)  ! [m/s]    

	end do

      pcoord(4) = - pcoord(4)*cof
      pcoord(5) = - pcoord(5)*cof

c      pcoorr(4) = - pcoorr(4)
c      pcoorr(5) = - pcoorr(5)

c      pcoord(4) = ( -DOTPR(DRDA,P_TRF1) - DVA )/drp
c      pcoord(5) = ( -DOTPR(DRDD,P_TRF1) - DVD )/drp

c      write (1,*) ' cof = ',cof, ' drp = ', drp
c      write (1,100) pcoord(4), pcoord(5)
c      write (1,101) alpha, delta
c      write (1,104) drda(1), drda(2), drda(3)
c      write (1,105) drdd(1), drdd(2), drdd(3)
c	write (1,*) v(1),v(2),v(3), dva, dvd
c	write (1,102) pcoord(4) - pcoord(5)
c      write (1,103) p_trf(1), p_trf(2), p_trf(3)
c      print 100,  pcoord(4), pcoord4
c      print 101,  pcoord(5), pcoord5
c100   format (3(2x,d22.16))
c101   format (' alpha =  ', d22.16,2x, ' delta =', d22.16)
c104   format (' drda =  ', 3(2x,d22.16))
c105   format (' drdd =  ', 3(2x,d22.16))
c 102   format (' p0 =  ', d22.16,2x, ' p1 =', d22.16)
c 103   format (' p_trf =  ', d22.16,2x, d22.16, 2x, d22.16)
c101   format (' p5 =  ', d22.16, 2x, d22.16)


      RETURN
      END

************************************************************************

      SUBROUTINE PNUT (CRX,CRY,CRZ,cof,drp)

C
C    SUBROUTINE PNUT
C
C    THIS SUBROUTINE COMPUTES THE PARTIAL DERIVATIVES WITH RESPECT
C  TO NUTATION OFFSETS IN EPSILON AND PSI
C
C  REVISION 1991 NOVEMBER 12  (N.Z.L.)
C
C  Revision 2002 MAY 13 (O.T.), Partial derivatives for EOP rates
C                         Change in common-block /PART/, /PDATA/
C  Revision 2006 MARCH 24 (O.T.), Changes for OCCAM 6.2
C  Last revision 2006 August 28 (O.T.) /COCDEL/   is off
C                        CRX,CRY,CRZ - from SUBROUTINE line
C
C-----------------------------------------------------------------------

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION CRX(3),CRY(3),CRZ(3)

      COMMON /PART / PCOORD(5),PCOORR(5),PNUTAT(2),PPOLAR(3,3)
      COMMON /PARM/ OMEGA
      COMMON /PDATA/ SD,CD,SIDGRE,DKDPSI(3),DKDEPS(3),
     * DSDXW(3,3),DSDYW(3,3), DSDUT(3,3),DKDPSI_v(3),DKDEPS_v(3)

C   AS PARTIAL DERIVATIVE WE WILL CONSIDER THE - ST * DKD

c
c    Multiplication on R1(y)*R2(x)*TRS
c
c

      P1 = -(CRX(3)*DKDEPS(1) + CRY(3)*DKDEPS(2) +
     *              CRZ(3)*DKDEPS(3))*cof/drp

      P1_v = -(CRX(3)*DKDEPS_v(1) + CRY(3)*DKDEPS_v(2) +
     *              CRZ(3)*DKDEPS_v(3))/drp


      P2 = -(CRX(3)*DKDPSI(1) + CRY(3)*DKDPSI(2) +
     *              CRZ(3)*DKDPSI(3))*cof/drp

      P2_v = -(CRX(3)*DKDPSI_v(1) + CRY(3)*DKDPSI_v(2) +
     *              CRZ(3)*DKDPSI_v(3))/drp

      pnutat(1) = p1 + p1_v
      pnutat(2) = p2 + p2_v


c      print *, pnutat(1), p1, p1_v
c      print *, pnutat(2), p2, p2_v

      RETURN
      END

************************************************************************

      SUBROUTINE PEOP (RQU,SN,PREC,cof,drp,vb, sidgre_dot,
     * dsdxw_t, dsdyw_t,dsdut_t, ppolar_rate)

C
C    SUBROUTINE PEOP
C
C    THIS SUBROUTINE COMPUTES THE PARTIAL DERIVATIVES WITH RESPECT
C  TO POLAR WOBBLE X, Y AND UT1
C
C  REVISION 1993 MAY 7  (N.Z.L.)
C
C  Revision 2002 MAY 13 (O.T.), Partial derivatives for EOP rates
C                         Change in common-block /PART/, /PDATA/
C  Last revision 2006 MARCH 24 (O.T.), Changes for OCCAM 6.2
C
C-----------------------------------------------------------------------

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION RQU(3),PCOORD(5),PCOORR(5)
      dimension px1(3), px2(3), px3(3), py1(3), py2(3), py3(3)
      dimension pu1(3), pu2(3), pu3(3), SN(3,3), PREC(3,3)
      dimension pnx1(3), pnx2(3), pnx3(3), pny1(3), pny2(3), pny3(3)
      dimension ppx1(3), ppx2(3), ppx3(3), ppy1(3), ppy2(3), ppy3(3)
      dimension pnu1(3), pnu2(3), pnu3(3), ppu1(3), ppu2(3), ppu3(3)


      dimension dsdxw_t(3), dsdyw_t(3), dsdut_t(3), ppolar_rate(3)
      dimension peopx(3), peopy(3), peoput(3), px(3), py(3), pu(3)

     

      COMMON /PART / PCOORD,PCOORR,PNUTAT(2),PPOLAR(3,3)
      COMMON /PDATA/ SD,CD,SIDGRE,DKDPSI(3),DKDEPS(3),
     * DSDXW(3,3),DSDYW(3,3), DSDUT(3,3),DKDPSI_v(3),DKDEPS_v(3)

C   AS PARTIAL DERIVATIVE WE WILL CONSIDER THE - DSD * K


      do i=1,3

         px1(i) = dsdxw(1,i)       !  R3(-H)*R1(y)*dR2(x)/dx*b
         px2(i) = dsdxw(2,i)
         px3(i) = dsdxw(3,i)
         py1(i) = dsdyw(1,i)       !  R3(-H)*dR1(y)/dy*R2(x)*b
         py2(i) = dsdyw(2,i)
         py3(i) = dsdyw(3,i)
         pu1(i) = dsdut(1,i)       !  dR3(-H)/dH*R1(y)*R2(x)*b
         pu2(i) = dsdut(2,i)
         pu3(i) = dsdut(3,i)

      end do

c      print *, 'prec'
c      print *, prec(1,1), prec(2,1), prec(3,1)
c      print *, prec(1,2), prec(2,2), prec(3,2)
c      print *, prec(1,3), prec(2,3), prec(3,3)
c      print *

c      print *, 'nut'
c      print *, sn(1,1), sn(2,1), sn(3,1)
c      print *, sn(1,2), sn(2,2), sn(3,2)
c      print *, sn(1,3), sn(2,3), sn(3,3)
c      print *

c      do i=1,3
c         print *, rqu(i)
c      end do

      do i=1,3

         pnx1(i)=0.d0
         pnx2(i)=0.d0
         pnx3(i)=0.d0
         pny1(i)=0.d0
         pny2(i)=0.d0
         pny3(i)=0.d0
         pnu1(i)=0.d0
         pnu2(i)=0.d0
         pnu3(i)=0.d0

         peopx(i) = 0.d0
         peopy(i) = 0.d0
         peoput(i) = 0.d0


         do j=1,3

            pnx1(i) = pnx1(i) + sn(i,j)*px1(j)    !  N*R3(-H)*R1(y)*dR2(x)/dx*b
            pnx2(i) = pnx2(i) + sn(i,j)*px2(j)
            pnx3(i) = pnx3(i) + sn(i,j)*px3(j)
            pny1(i) = pny1(i) + sn(i,j)*py1(j)    !  N*R3(-H)*dR1(y)/dy*R2(x)*b
            pny2(i) = pny2(i) + sn(i,j)*py2(j)
            pny3(i) = pny3(i) + sn(i,j)*py3(j)
            pnu1(i) = pnu1(i) + sn(i,j)*pu1(j)    !  N*dR3(-H)/dH*R1(y)*R2(x)*b
            pnu2(i) = pnu2(i) + sn(i,j)*pu2(j)
            pnu3(i) = pnu3(i) + sn(i,j)*pu3(j)
        
            peopx(i) =   peopx(i) + sn(i,j) * dsdxw_t(j) !  delay rate
            peopy(i) =   peopy(i) + sn(i,j) * dsdyw_t(j) !  delay rate
            peoput(i) = peoput(i) + sn(i,j) * dsdut_t(j) !  delay rate

         end do

         

      end do

       


      do i=1,3

         ppx1(i)=0.d0
         ppx2(i)=0.d0
         ppx3(i)=0.d0
         ppy1(i)=0.d0
         ppy2(i)=0.d0
         ppy3(i)=0.d0
         ppu1(i)=0.d0
         ppu2(i)=0.d0
         ppu3(i)=0.d0

         px(i)=0.d0
         py(i)=0.d0
         pu(i)=0.d0


         do j=1,3

            ppx1(i) = ppx1(i) + prec(i,j)*pnx1(j) !  P*N*R3(-H)*R1(y)*dR2(x)/dx*b
            ppx2(i) = ppx2(i) + prec(i,j)*pnx2(j)
            ppx3(i) = ppx3(i) + prec(i,j)*pnx3(j)
            ppy1(i) = ppy1(i) + prec(i,j)*pny1(j) !  P*N*R3(-H)*dR1(y)/dy*R2(x)*b
            ppy2(i) = ppy2(i) + prec(i,j)*pny2(j)
            ppy3(i) = ppy3(i) + prec(i,j)*pny3(j)
            ppu1(i) = ppu1(i) + prec(i,j)*pnu1(j) !  P*N*dR3(-H)/dH*R1(y)*R2(x)*b
            ppu2(i) = ppu2(i) + prec(i,j)*pnu2(j)
            ppu3(i) = ppu3(i) + prec(i,j)*pnu3(j)

            px(i) =   px(i) + prec(i,j) * peopx(j) !  delay rate
            py(i) =   py(i) + prec(i,j) * peopy(j) !  delay rate
            pu(i) =   pu(i) + prec(i,j) * peoput(j) !  delay rate

        end do
      end do


      PPOLAR(1,1) = (-DOTPR(RQU,PPX1)*cof-vb)/drp   !  partial derivative x one second before
      PPOLAR(2,1) = (-DOTPR(RQU,PPY1)*cof-vb)/drp   !  partial derivative y one second before
      PPOLAR(3,1) = (-DOTPR(RQU,PPU1)*cof-vb)/drp   !  partial derivative ut one second before

      PPOLAR(1,2) = (-DOTPR(RQU,PPX2)*cof-vb)/drp   !  partial derivative x one second after
      PPOLAR(2,2) = (-DOTPR(RQU,PPY2)*cof-vb)/drp   !  partial derivative y one second after
      PPOLAR(3,2) = (-DOTPR(RQU,PPU2)*cof-vb)/drp   !  partial derivative ut one second after

      PPOLAR(1,3) = (-DOTPR(RQU,PPX3)*cof-vb)/drp   !  partial derivative x in time
      PPOLAR(2,3) = (-DOTPR(RQU,PPY3)*cof-vb)/drp   !  partial derivative y in time
      PPOLAR(3,3) = (-DOTPR(RQU,PPU3)*cof-vb)/drp   !  partial derivative ut in time

      PPOLAR_rate(1) = -DOTPR(RQU,PX) * sidgre_dot   !  partial derivative x-dot in time 
      PPOLAR_rate(2) = -DOTPR(RQU,PY) * sidgre_dot   !  partial derivative y-dot in time
      PPOLAR_rate(3) = -DOTPR(RQU,PU) * sidgre_dot   !  partial derivative LOD in time


      RETURN
      END

************************************************************************

c      SUBROUTINE PLOVE
C
C     SUBROUTINE PLOVE
C
C  ESTA RUTINA CALCULA LAS DERIVADAS PARCIALES DE LAS COORDENADAS
C CILINDRICAS DE LA ESTACION, Z, EQ, LONG,  EN UN SISTEMA DE
C REFERENCIA TERRESTRE, CON RESPECTO DE LOS NUMEROS DE LOVE Y SHIDA
C H Y K, Y TAMBIEN DEL NUMERO DE LOVE PARA MAREA K1 HK1
C
C  REVISION 1990 NOVEMBER 19   (N.Z.L.)
C  REVISION 1997 APRIL, 26    (OT) IERS 96 Conventions
C  Last revision 2002 MAY 13 (O.T.), Partial derivatives for EOP rates
C                         Change in common-block /PART/, /PDATA/
C
C-----------------------------------------------------------------------
c
c      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c      DOUBLE PRECISION L02, L22, L3, LI_1, LI_2, L1_1, L1_2
c
c      DIMENSION PCOORD(5),PCOORR(5),PLOVED(3),PLOVER(3)
c      DIMENSION CGLONG(3),CPHI(3),CELHGT(3),CRX(3),CRY(3),CRZ(3)
c      DIMENSION HA(3),Z(3)
c
c      COMMON /PART / PCOORD,PCOORR,PNUTAT(2),PPOLAR(3,3)
c      COMMON /PDATA/ SD,CD,SIDGRE,DKDPSI(3),DKDEPS(3),
c     * DSDXW(3,3),DSDYW(3,3), DSDUT(3,3),DKDPSI_v(3),DKDEPS_v(3)
c      COMMON /COCDEL/ CGLONG,CPHI,CELHGT,CRX,CRY,CRZ,HA,Z,DE
c      COMMON /DATEb/ IDATJ,UTJ,XJ,YJ,IDELT
ccc      COMMON /LOVE/ H2,H2K1,L2
c      COMMON /LOVE/ H02, L02, H22, L22, H3, L3, HI_1, LI_1, HI_2, LI_2,
c     *              L1_1, L1_2
c      COMMON /MATH/ PI, TWOPI, ROG, ROS, ROZ, TOL
c      H02K1 = 0.5203d0
c
C   F FACTOR DE CONVERSION DE METROS A CENTIMETROS
c
c      F = 1.D2
c      RX=CRX(3)
c      RY=CRY(3)
c      RZ=CRZ(3)
c
C   PRERARACION DE LA ENTRADA DE LA RUTINA ETMB
C  P ES LA LATITUD GEOCENTRICA, R0 LA LONGITUD DEL RADIO VECTOR Y
C  G LA LONGITUD.
c
c      P = DATAN (RZ/(DSQRT(RX**2+RY**2)))
c      R0 = DSQRT (RX**2+RY**2+RZ**2)
c      G = CGLONG(3)
c      ELHGT=CELHGT(3)
c      PHI=CPHI(3)
c
C    PARAMETROS DE TIEMPO
C  T00 ES EL NUMERO DE SIGLOS JULIANOS DESDE 1900.0
c
c      DATJ = IDATJ + UTJ/TWOPI
c      T00 = (DATJ - 15019.5D0) / 36525.D0
c
C   COMIENZO DEL DO QUE CALCULA EL VALOR DE LAS COORDENADAS PARA
C  X+0.0002 Y X-0.0002 DE FORMA QUE SE OBTENGA LA DERIVADA
C  PARCIAL NUMERICA (X REPRESENTA A LOS NUMEROS DE LOVE)
c
c      DO 1 ILOVE = 1,3
c
C  SELECCIONAR EL NUMERO DE LOVE A VARIAR

c         IF (ILOVE.EQ.1) H02 = H02 + 0.2D-3
c         IF (ILOVE.EQ.2) H02K1 = H02K1 + 0.2D-3
c         IF (ILOVE.EQ.3) L02 = L02 + 0.2D-3

C   LLAMADA A LA RUTINA ETMB
C  WHEN POSSIBLE IT MUST BE CHANGED TO THE NEW TIDES ROUTINE

c         CALL ETMB (T00,UTJ,P,G,R0,PHI,ELHGT,DR,DPHI,DLAM)

C    CALCULAR LA VARIACION

c       RN = R0 + DR/F
c       PN = P + DPHI / (RN*F)

c       CF = DCOS (PN)
c       SF = DSIN (PN)

C     OBTENCION DE LAS COORDENADAS EN EL SISTEMA CILINDRICO

c       GNP = G + DLAM / (RN*CF*F)

c       ZP = RN * SF

c       EQP = RN * CF

C   CALCULAR AHORA LO MISMO PARA X-0.0002

c         IF (ILOVE.EQ.1) H02 = H02 - 0.4D-3
c         IF (ILOVE.EQ.2) H02K1 = H02K1 - 0.4D-3
c         IF (ILOVE.EQ.3) L02 = L02 - 0.4D-3

C   LLAMADA A LA RUTINA ETMB

c         CALL ETMB (T00,UTJ,P,G,R0,PHI,ELHGT,DR,DPHI,DLAM)

C    CALCULAR LA VARIACION

c       R0N = R0 + DR/F
c       PN = P + DPHI / (RN*F)

c       CF = DCOS (PN)
c       SF = DSIN (PN)

C     OBTENCION DE LAS COORDENADAS EN EL SISTEMA CILINDRICO

c       GNM = G + DLAM / (RN*CF*F)

c       ZM = RN * SF

c       EQM = RN * CF

C     CALCULO DE  (F(X+0.2D-3) - F(X-0.2D-3)) / 0.4D-3 , PARCIAL
C  NUMERICA

c       PARZ = (ZP - ZM) / 0.4D-3

c       PAREQ = (EQP - EQM) / 0.4D-3

c       PARGN = (GNP - GNM) / 0.4D-3

C     COMPROBAR A QUE NUMERO DE LOVE PERTENECEN ESTAS DERIVADAS
C     CALCULAR LA PARCIAL REAL (DTAU/DH = SUM DTAU/DR * DR/DH)

c      PLOVED(ILOVE) = PCOORD(1)*PARZ+PCOORD(2)*PAREQ*RN+PCOORD(3)*
c     *   PARGN*RN*DCOS(PN)

c      PLOVER(ILOVE) = PCOORR(1)*PARZ+PCOORR(2)*PAREQ*RN+PCOORR(3)*
c     *   PARGN*RN*DCOS(PN)

c1     CONTINUE

c      RETURN
c      END



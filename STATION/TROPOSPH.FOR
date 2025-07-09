      SUBROUTINE CDELAY (ZENDRY, istcat, isour, 
     * iav_imf,iav_rmf,iav_vmf,
     * iav_vm1, iav_grad, YESMAR, YESAZD, YESAGR)
      
C**********************************************************************
C
C    SUBROUTINE CDELAY
C
C    THIS ROUTINE COMPUTES THE MAPPING FUNCTIONS OF THE TROPOSPHERICAL
C  MODELS TO CORRECT THE OBSERVED DELAYS
C    NIELL MAPPING FUNCTIONS ARE USED
C    DRY DELAY (ZENDRY) IS CALCULATED USING MARINI's CODE
C
C
C  REVISION 1990 NOVEMBER 27 (N.Z.L.)
C  Revision 1997 April 09 (O.T.); Mapping function of A.Niell
C  Revision 2000 November 20 (O.T.)
C  Revision 2001 January 31 (O.T.)  Only mapping functions of A.Niell are usable
C  Revision 2001 April 17 (O.T.)  Remark by Johannes Boehm
c  revision 2002 July 29 by jboehm
c  revision 2003 May 19 by jboehm: iavail, iavailr and iavailq
c                                       were changed to integer arrays
C  revision 2005 Sep 5 by jboehm: Global Mapping Function GMF was added
C  Revision 2006 February 21 (O.T.)  Chages by Volker, Johannes, Oleg
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
C**********************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION HA(3),Z(3),DE(3)
      DIMENSION CGLONG(3),CPHI(3),CELHGT(3),CX(3),CY(3),CZ(3)
      DIMENSION CGLONG0(3),CPHI0(3),CELHGT0(3),HA0(3)

      

c+ jboehm
      dimension tmjdimf(8),adz200(20,8),adsmfw3(20,8),
     .          adazim(20,8),adelev(20,8),dz200(8),dsmfw3(8),
     .          dazim(8),delev(8),iavail(20)
      dimension tmjdrmf(8),
     .          adah(20,8),adbh(20,8),adch(20,8),
     .          adaw(20,8),adbw(20,8),adcw(20,8),
     .          dah(8),dbh(8),dch(8),daw(8),dbw(8),dcw(8),
     .          iavailr(20)
      dimension tmjdqmf(8),
     .          adahq(20,8),adawq(20,8),
     .          dahq(8),dawq(8),iavailq(20)
c+ jboehm, 2005-Oct-28
      dimension tmjdvmf(8),
     .          adahv(20,8),adawv(20,8),
     .          dahv(8),dawv(8),iavailv(20)
     *          ,pressure(20,8),temperat(20,8), pres1(8),temp1(8)              !  24-Mar-2014
c- jboehm, 2005-Oct-28

c+ jboehm 2005-Oct-24
      dimension tmjdgra(8),
     .          adgrnh(20,8),adgreh(20,8),
     .          adgrnw(20,8),adgrew(20,8),
     .          dgrnh(8),dgreh(8),dgrnw(8),dgrew(8),iavailg(20)
c- jboehm 2005-Oct-24


      dimension hd1(20,8), wd1(20,8), hd(8), wd(8)

      character*8 stat, lin*32, li2*18
      character*1 YESMAR, YESAZD, YESAGR
c- jboehm

      COMMON /CORR/ CTIDE(3),cmniel_w(2),cmniel_d(2),AXKT(2),therm_d
      COMMON /COCDEL/ CGLONG0,CPHI0,CELHGT0,CX,CY,CZ,HA0,Z,DE
      COMMON /COCDEL2/ CGLONG,CPHI,CELHGT,HA
      COMMON /METDAT/ TDRY,PRES,TWET,KW
c+ jboehm
      common /corrimf/ tmjdimf,adz200,adsmfw3,adazim,adelev,iavail,
     .                 azimuth
      common /corrrmf/ tmjdrmf,adah,adbh,adch,adaw,adbw,adcw,iavailr
      common /corrqmf/ tmjdqmf,adahq,adawq,iavailq
      common /corrout/ z200,smfw3,azim,elev,dimfh,dimfw,gimfh,grn,gre,
     .                 drmfh,drmfw,dqmfh,dqmfw
c- jboehm
c+ jboehm
      COMMON /PHYS/ C, FL, AH, AU, J20
      COMMON /MATH/ PI, TWOPI, ROG, ROS, ROZ, TOL
c- jboehm
C+ jboehm, 2005-Sep-5
      common /corrgmf/ gmfh,gmfw
C- jboehm, 2005-Sep-5
C+ jboehm, 2005-Oct-24
      common /corrgra/ tmjdgra,adgrnh,adgreh,adgrnw,adgrew,iavailg,
     .                 grnh,greh,grnw,grew
C- jboehm, 2005-Oct-24
c+ jboehm, 2005-Oct-28, OT, 5-Nov-2019

      common /corrvmf/ tmjdvmf,adahv,adawv,hd1,wd1,   !  7-Nov-2019 (OT)
     * iavailv,vmf1h,vmf1w,vmf1hd
     * ,vmf1wd,pressure,temperat                                                         !  24-Mar-2014

c- jboehm, 2005-Oct-28, 5-Nov-2019

C
C  CALL CORRECTIONS DUE TO TROPOSPHERE MODELS. THESE CORRECTIONS
C  ARE CHANGES IN THE OBSERVED DELAY AND DELAY RATE, BUT AT THIS STAGE
C  THEY WILL JUST BE THE MAPPING FUNCTION VALUES WHICH WILL BE
C  MULTIPLIED BY THE ZENITH TROPOSPHERIC DELAY IN THE MAIN PROGRAM
C

      CLT = DCOS (CPHI(3))
      SHA = DSIN (HA(3))
      CDEC = DCOS (DE(3))

c+ isour instead of 1
      READ (19,REC=isour) J,source,ra2000, de2000, Idat0, UT0
      READ (20,REC=ISTCAT) JP,STAT
c- jboehm
 
      CALL JULDAT (IY,ID,IH,IM,S,IDAT0,UT0,1)
      doy = ID+IH/24.d0+IM/1440.d0+S/86400.d0

C
C   MARINI TROPOSPHERIC MODEL (ZENDRY is calculated)
C

c
c   Height here is in meters
c

      F1 = 1.D0 - 0.00266D0 * DCOS(2.D0*CPHI(3)) - 0.28D-6*CELHGT(3)

  

c
c     Saastamoinen's correction
c

ccc      C3 = CELHGT(3)/1.d3    ! meter -> km 

ccc      BB = 1.1555d0 - 0.1565d0 * C3 + 0.0078d0 * C3**2  !   mbar
    
ccc      print *, pres, bb, z(3)

ccc      A1 = 0.0022768D0 * (PRES - BB * Dtan(z(3))**2) 
         
      A1 = 0.0022768D0 * PRES 



      ZENDRY = A1/F1             !  Dry zenith delay using Marini's approach
      

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc TESMER 23.02.05 NOW USE const. hydr. ZD directly           ccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         if (YESAZD.eq.'N'.or.YESAZD.eq.'n') then

c+ otitov 22.03.2006
c           if (isour.eq.1) print*, ' APRIORI CONST. ZD USED !!!'
c- otitov 22.03.2006

           rewind(44)
687        continue
           READ (44,688,END=686) li2
688        FORMAT(A18)
           IF (li2(1:8).EQ.STAT) then

             read(li2(10:18),'(d9.2)') ZENDRY
             

           else
             goto 687
           end if
686        continue

         else

c+ otitov 22.03.2006
c           if (isour.eq.1) print*, 'APRIORI ZD as usual!'
c- otitov 22.03.2006

         end if
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc TESMER 23.02.05 NOW USE const. hydr. ZD directly           ccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c+ jboehm

c  get mjd of observation
      tmjd = idat0 + ut0/TWOPI

c  create vectors for Lagrange interpolation
      do i = 1,8
        if (iav_imf.ne.0) then               !  OT
           dz200(i)  = adz200(istcat,i)
           dsmfw3(i) = adsmfw3(istcat,i)
           dazim(i)  = adazim(istcat,i)
           delev(i)  = adelev(istcat,i)
        end if
c  for rmf
        if (iav_rmf.ne.0) then               !  OT
           dah(i)    = adah(istcat,i)
           dbh(i)    = adbh(istcat,i)
           dch(i)    = adch(istcat,i)
           daw(i)    = adaw(istcat,i)
           dbw(i)    = adbw(istcat,i)
           dcw(i)    = adcw(istcat,i)
        end if
c  for qmf

        if (iav_vmf.ne.0) then               !  OT
           dahq(i)    = adahq(istcat,i)
           dawq(i)    = adawq(istcat,i)
        end if

c+ jboehm, 2005-Oct-28

        if (iav_vm1.ne.0) then               !  OT    21-Feb-2006

           dahv(i)   = adahv(istcat,i)
           dawv(i)   = adawv(istcat,i)

           hd(i) = hd1(istcat,i)             !  OT 7-Nov-2019
           wd(i) = wd1(istcat,i)             !     7-Nov-2019

cc           print *, i, hd(i), wd(i)
           

           temp1(i)   = temperat(istcat,i)   !  OT    24-Mar-2014
           pres1(i)   = pressure(istcat,i)   !  OT    24-Mar-2014

        end if                               !  OT    21-Feb-2006
c- jboehm, 2005-Oct-28
c+ jboehm, 2005-Oct-24
        if (iav_grad.ne.0) then              !  OT    21-Feb-2006
           dgrnh(i)  = adgrnh(istcat,i)
           dgreh(i)  = adgreh(istcat,i)
           dgrnw(i)  = adgrnw(istcat,i)
           dgrew(i)  = adgrew(istcat,i)
        end if                               !  OT    21-Feb-2006
c- jboehm, 2005-Oct-24
      end do


c  Lagrange interpolation
      if (iav_imf.ne.0) then                              !  OT

         call lagint4 (tmjdimf,dz200, 8,tmjd,z200)
         call lagint4 (tmjdimf,dsmfw3,8,tmjd,smfw3)
         call lagint4 (tmjdimf,dazim, 8,tmjd,azim)
         call lagint4 (tmjdimf,delev, 8,tmjd,elev)

c  calculate imf

         call imf (cphi(3),celhgt(3),z(3),z200,smfw3,dimfh,dimfw)
c
c  calculate gradient mapping function
c      print *, azimuth, azim, azimuth-azim

         zd = z(3) - (PI/2.d0 - elev)*dcos(azimuth - azim)
         call imf (cphi(3),celhgt(3),zd,z200,smfw3,gimfh,gimfw)

         grn = -(PI/2.d0 - elev)*dcos(0.0d0   - azim)*zendry
         gre = -(PI/2.d0 - elev)*dcos(PI/2.d0 - azim)*zendry

      end if
c  for rmf
      if (iav_rmf.ne.0) then                              !  OT
         call lagint4 (tmjdrmf,dah,8,tmjd,ad)
         call lagint4 (tmjdrmf,dbh,8,tmjd,bd)
         call lagint4 (tmjdrmf,dch,8,tmjd,cd)
         call lagint4 (tmjdrmf,daw,8,tmjd,aw)
         call lagint4 (tmjdrmf,dbw,8,tmjd,bw)
         call lagint4 (tmjdrmf,dcw,8,tmjd,cw)
c  calculate rmf
         call rmf (z(3),ad,bd,cd,aw,bw,cw,celhgt(3),drmfh,drmfw)
      end if
c  for qmf
      if (iav_vmf.ne.0) then                              !  OT
         call lagint4 (tmjdqmf,dahq,8,tmjd,adq)
         call lagint4 (tmjdqmf,dawq,8,tmjd,awq)
c  calculate qmf
         call qmf (adq,awq,cphi(3),celhgt(3),z(3),dqmfh,dqmfw)
      end if
c+ jboehm, 2005-Oct-28

      if (iav_vm1.ne.0) then                              !  OT    21-Feb-2006

         call lagint4 (tmjdvmf,dahv,8,tmjd,adv)
         call lagint4 (tmjdvmf,dawv,8,tmjd,awv)
         call lagint4 (tmjdvmf,temp1,8,tmjd,temp2)        !  OT    24-Mar-2014
         call lagint4 (tmjdvmf,pres1,8,tmjd,pres2)        !  OT    24-Mar-2014

         call lagint4 (tmjdvmf,hd,8,tmjd,hdelay)        !  OT    24-Mar-2014
         call lagint4 (tmjdvmf,wd,8,tmjd,wdelay)        !  OT    24-Mar-2014



ccc         print *, 'zendry old = ', zendry

         A2 = 0.0022768D0 * PRES2                         !  OT    24-Mar-2014





      end if                                              !  OT    21-Feb-2006
c- jboehm, 2005-Oct-28
c+ jboehm 2005-Oct-24
      if (iav_vm1.ne.0) then                              !  OT    21-Feb-2006

         call lagint4 (tmjdgra,dgrnh,8,tmjd,grnh)
         call lagint4 (tmjdgra,dgreh,8,tmjd,greh)
         call lagint4 (tmjdgra,dgrnw,8,tmjd,grnw)
         call lagint4 (tmjdgra,dgrew,8,tmjd,grew)

c+ empirically determined factors are applied
         grnh = grnh/1.9d0
         greh = greh/1.9d0
         grnw = grnw/1.4d0
         grew = grew/1.4d0
      end if                                              !  OT    21-Feb-2006
c- jboehm 2005-Oct-24

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc TESMER 26.07.05 read apriori gradients from GSFC from file ccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C NOW INITIALISE THE APRIORI GRADIENTS FROM gsfc_dao_9095.mgr
C UNITS ARE [mm]
c
c IT IS NOT VERY CLEAR IF EVERYTHING IS CALCULATED VERY RIGHT,
c BUT IT SEEMS TO BE OK... HAS TO BE CHECKED IN DETAIL...

      if (YESAGR.eq.'N'.or.YESAGR.eq.'n') then

        if (isour.eq.1) print*, ' APRIORI GRADIENTS APPLIED !!!'

        rewind(43)
697     continue
        READ (43,698,END=696) lin
698     FORMAT(A32)
        IF (lin(1:8).EQ.STAT) then
          read(lin(12:20),'(E9.2)') grnh   ! comment hannes
          read(lin(24:32),'(E9.2)') greh   ! comment hannes
        else
          grnh=0.d0     ! comment hannes
          greh=0.d0     ! comment hannes
          goto 697
        end if
696     continue

        grnw = 0.d0   ! comment hannes
        grew = 0.d0   ! comment hannes

c        grn=grn*0.001d0    ! [mm] -> [m]                ! comment hannes
c        gre=gre*0.001d0    ! [mm] -> [m]                ! comment hannes
c        azim=atan2(gre,grn)                             ! comment hannes
c        elev=PI/2.d0 + grn/(dcos(azim)*zendry)          ! comment hannes
c                                                        ! comment hannes
c        zd = z(3) - (PI/2.d0-elev)*dcos(azimuth-azim)   ! comment hannes
c        z(3) = zd                                       ! comment hannes

      end if

c+ jboehm, 2005-Oct-28
      if (YESAGR.eq.'Y'.or.YESAGR.eq.'y') then

        grnh = 0.d0
        greh = 0.d0
        grnw = 0.d0
        grew = 0.d0

      end if
c- jboehm, 2005-Oct-28

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc TESMER 26.07.05 read apriori gradients from GSFC from file ccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c  calculate imf
      call imf (cphi(3),celhgt(3),z(3),z200,smfw3,dimfh,dimfw)

c  calculate rmf
      call rmf (z(3),ad,bd,cd,aw,bw,cw,celhgt(3),drmfh,drmfw)

c  calculate vmf
      call vmf (adq,awq,cphi(3),celhgt(3),z(3),dqmfh,dqmfw)

c  calculate vmf1

      call vmf1 (adv,awv,doy,cphi(3),z(3),vmf1h,vmf1w,vmf1hd,vmf1wd) ! comment hannes
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc TESMER 23.02.05 NOW USE apriori gradients from GSFC        ccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
Cc  calculate gradient mapping function
C      zd = z(3) - (PI/2.d0 - elev)*dcos(azimuth - azim)
C      call imf (cphi(3),celhgt(3),zd,z200,smfw3,gimfh,gimfw)
C
Cc- jboehm
C
Cc+ jboehm


c      print *, ' dry = ', hdelay, zendry, iav_vm1
c      zendry = hdelay

c  calculate gradient mapping function

        zd = z(3) - (PI/2.d0 - elev)*dcos(azimuth - azim)
        call imf (cphi(3),celhgt(3),zd,z200,smfw3,gimfh,gimfw)



c      end if                                      ! comment hannes
c- jboehm, 2005-Oct-28

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc TESMER 23.02.05 NOW USE apriori gradients from GSFC        ccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



c      write (12,'(7f16.6)')
c     .       cphi(3),celhgt(3),z(3),z200,smfw3,dimfh,dimfw

c  set not available parameters to -9.d0
      if (iavail(istcat).eq.0) then
        z200  = -9.d0
        smfw3 = -9.d0
        azim  = -9.d0
        elev  = -9.d0
        dimfh = -9.d0
        dimfw = -9.d0
        gimfh = -9.d0
        grn   = -9.d0
        gre   = -9.d0
      end if

      if (iavailr(istcat).eq.0) then
        drmfh = -9.d0
        drmfw = -9.d0
      end if

      if (iavailq(istcat).eq.0) then
        dqmfh = -9.d0
        dqmfw = -9.d0
      end if

c- jboehm

c
c    WET NIELL MAPPING FUNCTION
c
      CALL cmniell_w ( Z(3), cphi(3), cmniel_w)

c
c    DRY NIELL MAPPING FUNCTION
c
      CALL cmniell_d ( Z(3), doy, cphi(3), celhgt(3), cmniel_d)



C+ jboehm, 2005-Sep-5
C    Global Mapping Function (Boehm, Niell, Schuh, 2005)

      CALL gmf (tmjd,cphi(3),cglong(3),celhgt(3),Z(3),gmfh,gmfw)

C- jboehm, 2005-Sep-5


      RETURN
      END

************************************************************************


      subroutine cmniell_w (z, latitude, wmf)

* new aen 930517 Routine to compute the new wmf2.0 mapping function which
*                depends only on latitude.


      integer*4 i

*   a,b,c       - the a,b,and c coefficients in the continued fraction
*                 form of Marini
*   beta        - intermediate term in calculation
*   gamma       - intermediate term in calculation
*   sine        - Sine of elevation angle
*   cose        - Cos of elevation angle
*   wmf(1)      - wet delay mapping function
*   wmf(2)      - d_wet_mapping_function/d_elevation
*   topcon      - Constant of top of mapping fuinction to ensure
*                 that value is 1.0000 at zenith

      real*8 a,b,c, beta, cose, wmf(2), gamma, sine, topcon

*   latitude   - latitude (degrees)
*   l          - absolute latitude
*   dl         - incremental latitude from last lat_wmf
*   elev       - elevation (degrees)
*   dl,da,db,dc  - used for interpolation

      real*8 lat_wmf(5), abc_w2p0(5,3)
      real*8 dl, da, db, dc
      real*8 latitude, l, elev, z, deg2rad, latitude_deg


*   define parameters used for calculating coefficients.

      data lat_wmf / 15, 30, 45, 60, 75/

*   coefficients are from fits to raytraces of the standard atmospheres
*   for July for latitudes 15, 45, 60, and 75 degrees latitude and for
*   January for 30 degrees latitude (930517).

      data abc_w2p0 /

     . 5.8021897e-4,5.6794847e-4,5.8118019e-4,5.9727542e-4,6.1641693e-4,
     . 1.4275268e-3,1.5138625e-3,1.4572752e-3,1.5007428e-3,1.7599082e-3,
     . 4.3472961e-2,4.6729510e-2,4.3908931e-2,4.4626982e-2,5.4736038e-2/


      deg2rad = 3.1415926535897932D0/180.d0

C
C     Transformation from radians to degrees
C
      latitude_deg=latitude/deg2rad
      elev=3.1415926535897932D0/2.d0 - z
      elev=elev/deg2rad

      l = dabs(latitude_deg)


*   Coefficients for the continued fraction expansion for each latitude.

*   for latitudes less than 15 degrees:

      if (l .le. lat_wmf(1)) then
         a = abc_w2p0(1,1)
         b = abc_w2p0(1,2)
         c = abc_w2p0(1,3)
      endif

*   for latitudes between 15 and 75  degrees:

      do i = 1,4
          if (l .gt. lat_wmf(i) .and. l .le. lat_wmf(i+1)) then
             dl = (l-lat_wmf(i))/(lat_wmf(i+1)-lat_wmf(i))
             da  =   abc_w2p0(i+1,1)-abc_w2p0(i,1)
             a   =   abc_w2p0(i,1) + dl*da
c     write(*,'(" dl,da ,a  ",6e15.6)')
c    .            dl,da ,a

             db  =   abc_w2p0(i+1,2)-abc_w2p0(i,2)
             b   =   abc_w2p0(i,2) + dl*db
c     write(*,'(" dl,db ,b ",6e15.6)')
c    .            dl,db ,b

             dc  =   abc_w2p0(i+1,3)-abc_w2p0(i,3)
             c   =   abc_w2p0(i,3) + dl*dc
c     write(*,'(" dl,dc ,c ",6e15.6)')
c    .            dl,dc ,c

          endif
      end do

*   for latitudes greater than 75 degrees:

      if (l .ge. lat_wmf(5)) then
         a = abc_w2p0(5,1)
         b = abc_w2p0(5,2)
         c = abc_w2p0(5,3)
      endif

*   Now the coefficients exist; calculate the mapping function, wmf(1),
*       and the change of mapping function with elevation,
*       dwmf/d_el =wmf(2).
*   To calculate the delay-rate correction, d_tau/dt:
*       d_tau/dt = d_tau_zen/dt * wmf(1) + tau_zen * dwmf/d_el * d_el/dt


c      a = a + 0.0005d0

      sine  = sin( elev * deg2rad)
      cose  = cos( elev * deg2rad)
      beta  = b/( sine + c )
      gamma = a/( sine + beta)


   
      topcon = (1.d0 + a/(1.d0 + b/(1.d0 + c)))


      wmf(1) = topcon / ( sine + gamma )

      wmf(2) = -topcon / ( sine + gamma )**2 *
     .         ( cose - a/( sine + beta)**2 * cose *
     .         ( 1.d0 - b/( sine + c )**2 ) )


      
      return
      end


      subroutine cmniell_d (z,doy,latitude,height,hmf)


*     Routine to compute the hydrostatic mapping function nhmf2 which
*     depends on DOY (day of year) and station position (latitude
*     and height above geoid; use ellipsoid height for now).

* 931007 aen NEW nhmf2
* 951129 aen MOD nmfh2p1 Add derivative of height correction wrt elevation to hmf(2).
*                NOTE change in spelling of subroutine from nhmf to nmfh.
* 960530 aen MOD correction: add doy term for |elevations| .ge. 75 degrees.
* 960814 aen MOD change correction for southern latitude so day number (doy) isn't changed.

      integer*4 i

*   a,b,c       - the a,b,and c coeffiecents in the continued fraction
*                 form of Marini
*   beta        - intermediate term in calculation
*   gamma       - intermediate term in calculation
*   sine        - sine of elevation angle
*   cose        - cos of elevation angle
*   hmf(1)      - delay mapping function
*   hmf(2)      - d_mapping_function/d_elevation (dhmf2/d_el)
*   topcon      - constant of top of mapping function to ensure
*                 that value is 1.0000 at zenith

      real*8 a,b,c, beta, cose, hmf(2), gamma, sine, topcon

*   height     - height of site above geoid (meters)
*   hs_km      - Height of site in kms.
*   latitude   - latitude (degrees)
*   l          - absolute latitude
*   dl         - incremental latitude from last lat_hmf
*   elev       - elevation (degrees)
*   epoch      - if Julian date of observation is known for the observation,
*              - then epoch can be used to get day of year.
*              - (if epoch  is passed as argument, then un-comment the
*                 line converting epoch to doy.)
*   doy        - days since Dec 31
*   doy_atm    - doy for atmosphere relative to Jan 28.
*   doyr_atm   - doy_atm in radians;
*   cost       - cosine(day of year)
*   doy2rad    - convert doy to radians

      real*8 doy, latitude, height, elev, latitude_deg
      real*8 hs_km, l, dl, doy_atm, doyr_atm, cost
      real*8 doy2rad, deg2rad, z

*   lat_hmf     - latitudes at which coefficients are defined (5).
*   abc_avg     - continued fraction coefficients at latitudes lat_hmf
*   abc_amp     - amplitude of annual variation of abc_avg
*   daavg, daamp, etc - incremental values for interpolation
*   aavg,  aamp,  etc - average and amplitude at latitude

      real*8 lat_hmf(5)
      real*8 abc_avg(5,3), abc_amp(5,3)
      real*8 daavg, daamp, dbavg, dbamp, dcavg, dcamp
      real*8 aavg,  aamp,  bavg,  bamp,  cavg,  camp

*   a_ht, b_ht, c_ht - parameters for continued fraction for height corr'n.
*   dhcc_del    - derivative of height correction coefficient with elevation
*   dht_corr_del - derivative of height correction with elevation

      real*8 a_ht, b_ht, c_ht, ht_corr_coef, ht_corr
      real*8 dhcc_del, dht_corr_del

*   define parameters used for calculating coefficients.

      data lat_hmf / 15, 30, 45, 60, 75/

      data abc_avg /
     .1.2769934e-3,1.2683230e-3,1.2465397e-3,1.2196049e-3,1.2045996e-3,
     .2.9153695e-3,2.9152299e-3,2.9288445e-3,2.9022565e-3,2.9024912e-3,
     .62.610505e-3,62.837393e-3,63.721774e-3,63.824265e-3,64.258455e-3/

      data abc_amp /
     .  0.0,   1.2709626e-5, 2.6523662e-5, 3.4000452e-5, 4.1202191e-5,
     .  0.0,   2.1414979e-5, 3.0160779e-5, 7.2562722e-5, 11.723375e-5,
     .  0.0,   9.0128400e-5, 4.3497037e-5, 84.795348e-5, 170.37206e-5/

      data a_ht / 2.53e-5/
     .     b_ht / 5.49e-3/
     .     c_ht / 1.14e-3/

*   conversions:

      doy2rad = 2.d0*3.1415926535897932d0/365.25d0
      deg2rad = 3.1415926535897932D0/180.d0

*   convert height in meters to kilometers

      hs_km  = height/1000.d0

*   If Julian date is used for epoch, then calculate day of year;
*      use 1980 Jan 0 as reference epoch.

*     doy = epoch - 2444238.5

* mod aen 930517 Use phase of 28 days (winter extremum corresponds to Jan 28)
*                based on least-square fit to
*                raytrace of radiosonde data for DRT, ELP, ALB, CHH, FAI,
*                MUN, and LIH.
*
      doy_atm  = doy - 28.d0

C
C     Transformation from radians to degrees
C
      latitude_deg=latitude/deg2rad
      elev=3.1415926535897932D0/2.d0-z
      elev=elev/deg2rad

*   to account for the six month difference in seasons between hemispheres,
*   add 365.25/2 days to doy_atm if station is in the southern hemisphere.

      l = dabs(latitude_deg)
      if (latitude .lt. 0) doy_atm = doy_atm + 365.25d0/2.d0

      doyr_atm = doy_atm * doy2rad

      cost = dcos(doyr_atm)

*   Coefficients for the continued fraction expansion for each latitude.

*   for latitudes less than 15 degrees:

      if (l .le. lat_hmf(1)) then
         a = abc_avg(1,1)
         b = abc_avg(1,2)
         c = abc_avg(1,3)
      endif

*   for latitudes between 15 and 75  degrees:

      do i = 1,4
          if (l .gt. lat_hmf(i) .and. l .le. lat_hmf(i+1)) then
             dl = (l-lat_hmf(i))/(lat_hmf(i+1)-lat_hmf(i))
             daavg =   abc_avg(i+1,1)-abc_avg(i,1)
             daamp =   abc_amp(i+1,1)-abc_amp(i,1)
             aavg  =   abc_avg(i,1) + dl*daavg
             aamp  =   abc_amp(i,1) + dl*daamp
             a     = aavg - aamp*cost
c     write(*,'(" dl,daavg,daamp,aavg,aamp,a ",6e15.6)')
c    .            dl,daavg,daamp,aavg,aamp,a

             dbavg =   abc_avg(i+1,2)-abc_avg(i,2)
             dbamp =   abc_amp(i+1,2)-abc_amp(i,2)
             bavg  =   abc_avg(i,2) + dl*dbavg
             bamp  =   abc_amp(i,2) + dl*dbamp
             b     = bavg - bamp*cost
c     write(*,'(" dl,dbavg,dbamp,bavg,bamp,b ",6e15.6)')
c    .            dl,dbavg,dbamp,bavg,bamp,b

             dcavg =   abc_avg(i+1,3)-abc_avg(i,3)
             dcamp =   abc_amp(i+1,3)-abc_amp(i,3)
             cavg  =   abc_avg(i,3) + dl*dcavg
             camp  =   abc_amp(i,3) + dl*dcamp
             c     = cavg - camp*cost
c     write(*,'(" dl,dcavg,dcamp,cavg,camp,c ",6e15.6)')
c    .            dl,dcavg,dcamp,cavg,camp,c

          endif
      end do

*   for latitudes greater than 75 degrees:

      if (l .ge. lat_hmf(5)) then
caen 960530 correct the next three lines:
c        a = abc_avg(5,1)
c        b = abc_avg(5,2)
c        c = abc_avg(5,3)
         a = abc_avg(5,1) - abc_amp(5,1)*cost
         b = abc_avg(5,2) - abc_amp(5,2)*cost
         c = abc_avg(5,3) - abc_amp(5,3)*cost
      endif

*   Now the coefficients exist; calculate for the sea level part
*   the mapping function, hmf(1), and the derivative wrt elevation
*   dhmf/d_el = hmf(2).

*   To get delay-rate correction d_tau/dt:
*      d_tau/dt = d_tau-zen/dt*hmf(1) + tau-zen*hmf(2)*d_el/dt
*      where  hmf(2)=dhmf/d_el

c      a = a + 0.0005d0
c      b = b - 0.0010d0
c      c = c - 0.0100d0

      sine   = dsin(elev * deg2rad)
      cose   = dcos(elev * deg2rad)
      beta   = b/( sine + c )
      gamma  = a/( sine + beta)

      topcon = (1.d0 + a/(1.d0 + b/(1.d0 + c)))

      hmf(1) =     topcon / ( sine + gamma )

      hmf(2) =     -topcon*cose / ( sine + gamma )**2 *
     .            ( 1.d0 - a/ ( sine + beta)**2 *
     .            ( 1.d0 - b/ ( sine + c   )**2 ) )

c     write(*,'("sine, cose, beta, gamma, topcon = ", 5f10.5)')
c    .           sine, cose, beta, gamma, topcon
c     write(*,'("hmf(1), hmf(2) = ", 2f10.4)') hmf(1), hmf(2)
c     write(*,'("hmf(1), hmf(2) = ", 2f10.4)') hmf

*   Apply height correction to mapping function and derivative wrt elevation:
*
*      1) height correction coefficient is
*         1/sine(elev) - continued fraction(a_ht,b_ht,c_ht).
*      2) height correction is ht_corr_coef times height in km.
*      3) height correction to derivative wrt elevation is (derivative of
*         height correction coefficient wrt elevation)*height in km.

      beta   = b_ht/( sine + c_ht )
      gamma  = a_ht/( sine + beta)
      topcon = (1.d0 + a_ht/(1.d0 + b_ht/(1.d0 + c_ht)))

      ht_corr_coef = 1/sine - topcon/(sine + gamma)
      ht_corr      = ht_corr_coef * hs_km
      hmf(1)       = hmf(1)     + ht_corr

*    951129 The derivative of the height correction wrt elevation is added
*    to hmf(2) after Chris Jacobs pointed out the magnitude of the term.

      dhcc_del   = -cose/sine**2
     .             +topcon*cose / ( sine + gamma)**2 *
     .            ( 1.d0 - a_ht/ ( sine + beta)**2 *
     .            ( 1.d0 - b_ht/ ( sine + c_ht)**2) )
      dht_corr_del = dhcc_del * hs_km

      hmf(2)       = hmf(2)     + dht_corr_del

c     write(*,'("ht_corr_coef, ht_corr, hs_km, hmf(1) = ", 4f15.6)')
c    .           ht_corr_coef, ht_corr, hs_km, hmf(1)

      return
      end



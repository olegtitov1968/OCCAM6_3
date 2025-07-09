cc

      subroutine getimfpar(flat,flon,tmjd,hyd,wet,az,el,an,ae)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  This subroutine calculates the hydrostatic (z200) and wet (smfw3) parameters
c  for the isobaric mapping functions by A. Niell. It also provides the
c  angles (azimuth and elevation) for the gradients (tilting of the 200 mbar
c  surface).
c
c  input parameters
c  ----------------
c  flat:  latitude radians
c  flon:  longitude in radians
c  tmjd:  modified Julian date
c
c  output parameters
c  -----------------
c  hyd: hydrostatic parameter (z200), height of the 200 mbar surface
c  wet: wet parameter (smfw3)
c  el:  elevation corresponding to the normal to the 200 mbar surface
c       in radians (always lower than PI/2)
c  az:  azimuth corresponding to the elevation in radians
c  an:  tilting to the north direction in radians
c  ae:  tilting to the east direction in radians
c
c  jboehm, 2002 May 8
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit double precision (a-h,o-z)

      dimension z200(4),smfw(4),tmjdr(4),
     .          z200n(4),z200s(4),z200e(4),z200w(4)
      character*10 fnameh(4),fnamew(4)

      COMMON /PHYS/ C, FL1, AH, AU, J20
      COMMON /MATH/ PI, TWOPI, ROG, ROS, ROZ, TOL

      idj = int(tmjd)
      utj = mod(tmjd,1.d0)*TWOPI

c  get year, day of year and hour of the observation

      s0 = 0.d0
      
      call juldat (iy,id,ih,imin,s,idj,utj,1)

      if (ih.eq.24) then
        id = id + 1
        ih = 0
      end if

c  Lagrange interpolation is used with two epochs before and two after the obs.
c  get the second epoch for the interpolation
      iy2 = mod(iy,100)
      id2 = id
      ih2 = ih - mod(ih,6)
      write (fnameh(2),'(a1,i2.2,i3.3,a2,i2.2)') 'h',iy2,id2,'.h',ih2
      write (fnamew(2),'(a1,i2.2,i3.3,a2,i2.2)') 'w',iy2,id2,'.h',ih2
      call juldat (iy,id2,ih2,0,s0,idj,utj,0)
      tmjdr(2) = idj + utj/TWOPI

c  get the first epoch for the interpolation
      iy1 =iy2
      if (ih2.eq.0) then
        id1 = id2 - 1
        ih1 = 18
       else
        id1 = id2
        ih1 = ih2 - 6
      end if
      write (fnameh(1),'(a1,i2.2,i3.3,a2,i2.2)') 'h',iy1,id1,'.h',ih1
      write (fnamew(1),'(a1,i2.2,i3.3,a2,i2.2)') 'w',iy1,id1,'.h',ih1
      call juldat (iy,id1,ih1,0,s0,idj,utj,0)
      tmjdr(1) = idj + utj/TWOPI
         
c  get the third epoch for the interpolation
      iy3 =iy2
      if (ih2.eq.18) then
        id3 = id2 + 1
        ih3 = 00
       else
        id3 = id2
        ih3 = ih2 + 6
      end if
      write (fnameh(3),'(a1,i2.2,i3.3,a2,i2.2)') 'h',iy3,id3,'.h',ih3
      write (fnamew(3),'(a1,i2.2,i3.3,a2,i2.2)') 'w',iy3,id3,'.h',ih3
      call juldat (iy,id3,ih3,0,s0,idj,utj,0)
      tmjdr(3) = idj + utj/TWOPI
         
c  get the fourth epoch for the interpolation
      iy4 =iy3
      if (ih3.eq.18) then
        id4 = id3 + 1
        ih4 = 00
       else
        id4 = id3
        ih4 = ih3 + 6
      end if
      write (fnameh(4),'(a1,i2.2,i3.3,a2,i2.2)') 'h',iy4,id4,'.h',ih4
      write (fnamew(4),'(a1,i2.2,i3.3,a2,i2.2)') 'w',iy4,id4,'.h',ih4
      call juldat (iy,id4,ih4,0,s0,idj,utj,0)
      tmjdr(4) = idj + utj/TWOPI
         
c  get the interpolated values for the coordinates given
      do i = 1,4
        call intpt1(fnameh(i),flat,flon,z200(i))
        call intpt1(fnamew(i),flat,flon,smfw(i))
      end do

c  get the interpolated values 0.5 deg north, south, east and west of the point
      dl = 0.5d0/180.0d0*PI
      do i = 1,4
        call intpt1(fnameh(i),flat+dl,flon,z200n(i))
        call intpt1(fnameh(i),flat-dl,flon,z200s(i))
        call intpt1(fnameh(i),flat,flon-dl,z200w(i))
        call intpt1(fnameh(i),flat,flon+dl,z200e(i))
      end do

c  get the interpolated values
c      call LAGINT1 (tmjdr,z200, 4,tmjd,hyd)
c      call LAGINT1 (tmjdr,smfw, 4,tmjd,wet)
c      call LAGINT1 (tmjdr,z200n,4,tmjd,hydn)
c      call LAGINT1 (tmjdr,z200s,4,tmjd,hyds)
c      call LAGINT1 (tmjdr,z200e,4,tmjd,hyde)
c      call LAGINT1 (tmjdr,z200w,4,tmjd,hydw)

c  get some parameters of the ellipsoid
c  flattening (a-b)/a
      fl = 1/fl1
c  semiminor axis
      bh = ah*(1 - fl)
c  excentricity squared
      e2 = 2*fl - fl**2
c  2nd excentricity squared
      ep2 = e2/(1 - e2)
c  a*a/b
      ch = ah**2/bh
c  V
      V = dsqrt(1 + ep2*(dcos(flat))**2)
c  Meridiankruemmungsradius
      dM = ch/V**3
c  Querkruemmungsradius
      dN = ch/V

c  lengths of one degree in the height of z200 along the meridian and the parallel
      dlatlen = (dM + hyd)*PI/180.d0
      dlonlen = (dN + hyd)*dcos(flat)*PI/180.d0

c  get angle in sn and we direction
      an = datan2((hyds - hydn),dlatlen)
      ae = datan2((hydw - hyde),dlonlen)

c  get gradient angles (azimuth and elevation)
      az = datan2(dsin(ae),dsin(an))
      el = datan(dsqrt((dsin(an))**2 + (dsin(ae))**2))
      if (az.lt.(-1.0d0*PI)) az = az + 2.0d0*PI
      el = PI/2.0d0 - el

      return
      end subroutine

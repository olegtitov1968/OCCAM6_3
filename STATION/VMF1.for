      subroutine vmf1 (ah,aw, doy, dlat,zd,vmf1h,vmf1w, 
     * vmf1h_dot, vmf1w_dot)

C     This subroutine determines the VMF1 (Vienna Mapping Functions 1)
C     Reference: Boehm, J., B. Werl, H. Schuh: Troposphere Mapping Functions
C                for GPS and VLBI from ECMWF operational analysis data
C                accepted by JGR, 2005
C
C     input data
C     ----------
C     ah:   hydrostatic coefficient a (www.hg.tuwien.ac.at/~ecmwf1)
C     aw:   wet coefficient a         (www.hg.tuwien.ac.at/~ecmwf1)
C     dmjd: modified julian date
C     dlat: latitude in radians
C     zd:   zenith distance in radians
C
C     output data
C     -----------
C     vmf1h: hydrostatic mapping function
C     vmf1w: wet mapping function
C
C     Johannes Boehm, 2005 October 2
C
c     vmf1h: hydrostatic mapping function for delay rate  (OT)  5-Nov-2019
C     vmf1w: wet mapping function for delay rate          (OT)  5-Nov-2019

      implicit double precision (a-h,o-z)

      pi = 3.141592653589d0

C     reference day is 28 January 1980
C     this is taken from Niell (1996) to be consistent

c      doy = dmjd  - 44239.d0 + 1 - 28

c      print *, dmjd -44239.d0 -27.d0 , doy


c
c     44239 - 1 Jan 1980
c

      bh = 0.0029d0
      c0h = 0.062d0

      if (dlat.lt.0.d0) then   ! southern hemisphere
          phh  = pi
          c11h = 0.007d0
          c10h = 0.002d0
      else                     ! northern hemisphere
          phh  = 0.d0
          c11h = 0.005d0
          c10h = 0.001d0
      end if

      ch = c0h + ((dcos(doy/365.d0*2.d0*pi + phh)+1.d0)*c11h/2.d0
     .     + c10h)*(1.d0-dcos(dlat))


      sine   = dsin(pi/2.d0 - zd)
      cose =   dcos(pi/2.d0 - zd)

      beta   = bh/( sine + ch  )
      gamma  = ah/( sine + beta)
      topcon = (1.d0 + ah/(1.d0 + bh/(1.d0 + ch)))
      vmf1h   = topcon/(sine+gamma)
c
c  Oleg Titov 5-Nov-2019
c
      vmf1h_dot = -topcon*cose/ (sine + gamma)**2*
     * (1.d0 - ah/(sine + beta)**2 *
     * (1.d0 - bh/(sine + ch )**2))


      bw = 0.00146d0
      cw = 0.04391d0
      betaw   = bw/( sine + cw )
      gammaw  = aw/( sine + betaw)
      topconw = (1.d0 + aw/(1.d0 + bw/(1.d0 + cw)))
      vmf1w   = topcon/(sine+gammaw)
c
c   Oleg Titov 5-Nov-2019
c
      vmf1w_dot = -topconw * cose/ (sine + gammaw)**2*
     * (1.d0 - aw/(sine + betaw)**2*
     * (1.d0 - bw/(sine + cw )**2))
     

      end subroutine

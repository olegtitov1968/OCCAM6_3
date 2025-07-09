cc

      subroutine qmf (ah,aw,dlat,height,zd,dqrmfh,dqrmfw)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  This subroutine uses ECMWF data to generate the
c  'quasi-raytraced' mapping funtions qmf
c
c  input data
c  ----------
c  dlat:   latitude in radians
c  height: ellipsoidal height in meters
c  zd:     zenith distance in radians
c  ah:     hydrostatic coefficient
c  aw:     wet coefficient
c
c  output data
c  -----------
c  dqrmfh:  hydrostatic mf
c  dqrmfw:  wet mf
c
c  jboehm, 2003 Feb 27
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit double precision (a-h,o-z)

      COMMON /MATH/ PI, TWOPI, ROG, ROS, ROZ, TOL

c  coefficients for the determination of the mapping functions

      zd1  = zd  *180.0d0/PI

      bh  = 0.002905d0

      cm0 = 0.0634d0
      cm1 = 0.0014d0
      ch  = cm0 + cm1*dcos(2.d0*dlat)

c  imfh

      elev = 90.d0 - zd1

      sine   = dsin(elev * pi/180.d0)
      cose   = dcos(elev * pi/180.d0)
      beta   = bh/( sine + ch )
      gamma  = ah/( sine + beta)
      topcon = (1.d0 + ah/(1.d0 + bh/(1.d0 + ch)))

      dqrmfh = topcon/(sine+gamma)

c  height correction

      a_ht = 2.53d-5
      b_ht = 5.49d-3
      c_ht = 1.14d-3

      hs_km  = height/1000.d0

      beta         = b_ht/( sine + c_ht )
      gamma        = a_ht/( sine + beta)
      topcon       = (1.d0 + a_ht/(1.d0 + b_ht/(1.d0 + c_ht)))
      ht_corr_coef = 1/sine - topcon/(sine + gamma)
      ht_corr      = ht_corr_coef * hs_km
      dqrmfh       = dqrmfh + ht_corr

c  wet mapping function

c coefficients by Niell 96

      bw = 0.00146
      cw = 0.04391

c      aw = aw + 2.0e-4

      sine  = dsin( elev * pi/180.d0)
      cose  = dcos( elev * pi/180.d0)
      beta  = bw/( sine + cw )
      gamma = aw/( sine + beta)
      topcon = (1.d0 + aw/(1.d0 + bw/(1.d0 + cw)));

      dqrmfw = topcon/(sine+gamma)

      return
      end

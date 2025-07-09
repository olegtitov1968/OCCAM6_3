cc

      subroutine rmf (zd,ah,bh,ch,aw,bw,cw,height,drmfh,drmfw)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  This subroutine uses raytracing data to generate the
c  mapping funtions fmf (hydrostatic drmfh and wet drmfw).
c
c  input data
c  ----------
c  zd:     zenith distance in radians
c  ah,bh,ch,aw,bw,cw: hyd and wet coefficients
c
c  output data
c  -----------
c  drmfh:  hydrostatic mf
c  drmfw:  wet mf
c
c  jboehm, 2003 Feb 18
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit double precision (a-h,o-z)

      COMMON /MATH/ PI, TWOPI, ROG, ROS, ROZ, TOL

c  coefficients for the determination of the mapping functions

      zd1   = zd  *180.0d0/PI

c  rmfh

      elev = 90.d0 - zd1

      sine   = sin(elev * pi/180.d0);
      cose   = cos(elev * pi/180.d0);
      beta   = bh/( sine + ch );
      gamma  = ah/( sine + beta);
      topcon = (1.d0 + ah/(1.d0 + bh/(1.d0 + ch)));

      drmfh = topcon/(sine+gamma);

c  height correction

      a_ht = 2.53e-5;
      b_ht = 5.49e-3;
      c_ht = 1.14e-3;

      hs_km  = height/1000.d0;

      beta         = b_ht/( sine + c_ht );
      gamma        = a_ht/( sine + beta);
      topcon       = (1.d0 + a_ht/(1.d0 + b_ht/(1.d0 + c_ht)));
      ht_corr_coef = 1/sine - topcon/(sine + gamma);
      ht_corr      = ht_corr_coef * hs_km;

      drmfh = drmfh + ht_corr

c  wet mapping function

      sine  = sin( elev * pi/180.d0);
      cose  = cos( elev * pi/180.d0);
      beta  = bw/( sine + cw );
      gamma = aw/( sine + beta);
      topcon = (1.d0 + aw/(1.d0 + bw/(1.d0 + cw)));

      drmfw = topcon/(sine+gamma);

      return
      end

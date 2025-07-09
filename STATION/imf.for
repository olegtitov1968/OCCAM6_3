cc

      subroutine imf (dlat,height,zd,z200,smfw3,dimfh,dimfw)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  This subroutine uses ECMWF data to generate the isobaric
c  Niell mapping funtions imf (hydrostatic dimfh and wet dimfw).
c
c  input data
c  ----------
c  dlat:   latitude in radians
c  height: ellipsoidal height in meters
c  zd:     zenith distance in radians
c  z200:   parameter for imfh
c  smfw3:  parameter for imfw
c
c  output data
c  -----------
c  dimfh:  hydrostatic path delay
c  dimfw:  wet path delay
c
c  jboehm, 2002 March 28
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit double precision (a-h,o-z)

      COMMON /MATH/ PI, TWOPI, ROG, ROS, ROZ, TOL

c  coefficients for the determination of the mapping functions

      dlat1 = dlat*180.0d0/PI
      zd1   = zd  *180.0d0/PI

      a00    = 0.00124;
      a01    = 4.0e-5;
      dlata0 = 2.0;
      dadz0  = 7.4e-8;
      dadz1  = -1.6e-8;
      dlatd0 = 0.0;
      bm0    = 0.002905;
      cm0    = 0.0634;
      cm1    = 0.0014;
      dlatc0 = 0.0;
      z0     = 11836.0;
      z1     = 619.0;
      dlatz0 = 3.0;

      alat =      a00 +   a01*cos(2*(dlat1-dlata0)*pi/180.d0);
      dadzlat = dadz0 + dadz1*cos(2*(dlat1-dlatd0)*pi/180.d0);
      zm =         z0 +    z1*cos(2*(dlat1-dlatz0)*pi/180.d0);

      a = alat + dadzlat*(z200-zm);
      b = bm0;
      c = cm0 + cm1*cos(2*(dlat1-dlatc0)*pi/180);

c  imfh

      elev = 90.d0 - zd1

      sine   = sin(elev * pi/180.d0);
      cose   = cos(elev * pi/180.d0);
      beta   = b/( sine + c );
      gamma  = a/( sine + beta);
      topcon = (1.d0 + a/(1.d0 + b/(1.d0 + c)));

      dimfh = topcon/(sine+gamma);

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
      dimfh        = dimfh + ht_corr;

c  wet mapping function as of 030307

      a0   =  6.8827e-004;
      a_ht = -1.6580e-007;
      dads = -2.0795e-004;
      b0 =    1.3503e-003;
      dbds =  1.8882e-004;
      c0 =    3.9647e-002;
      dcds =  4.8581e-003;
      smfw0 = 15.5;

      aw  = (smfw3-smfw0)*dads + a0 + a_ht*height;
      bw  = (smfw3-smfw0)*dbds + b0;
      cw  = (smfw3-smfw0)*dcds + c0;

      sine  = sin( elev * pi/180.d0);
      cose  = cos( elev * pi/180.d0);
      beta  = bw/( sine + cw );
      gamma = aw/( sine + beta);
      topcon = (1.d0 + aw/(1.d0 + bw/(1.d0 + cw)));

      dimfw = topcon/(sine+gamma);

      return
      end

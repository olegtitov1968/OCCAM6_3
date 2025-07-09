cxxxxxx
c   subroutines to compute the diurnal and semidiurnal variations
c   in the earth orientation from the version of Richard Ray's ocean
c   tide model that was listed in IERS Technical Note 21, July 1996.
c   This code includes the variations from 71 diurnal and semidiurnal
c   terms instead of the 8 that are listed in the report.
cxxxxxx

      subroutine eanes (time, corx, cory, corut1)
c
c...Purpose: to compute the diurnal and semidiurnal variations
c            in EOP (x,y,UT1) from ocean tides
c
c...Coded by: Richard Eanes, UT/CSR, Feb 1997
c
c...Input: time = Modified Julian Date
c
c...Output: eop = (delta_x, delta_y, delta_UT1)
c                 microarcsec for x and y, microsec for UT1
c
      real*8 eop(3), orthow(12,3), h(12), alpha, time
      real*8 corx, cory, corut1
c
c...diurnal and semidiurnal orthoweights fit to the 8 constituents
c   listed in IERS Technical Note 21, July 1996 which are from the
c   paper "Diurnal and Semidiurnal Variations in the Earth's Rotation
c   Rate Induced by Ocean Tides"  by Ray,R.D., Steinberg,D.J.,
c   Chao,B.F., and Cartwright, D.E., Science, 264, pp. 830-832.
c
      data orthow /
     . -6.77832,-14.86323,  0.47884, -1.45303,  0.16406,  0.42030,
     .  0.09398, 25.73054, -4.77974,  0.28080,  1.94539, -0.73089,
     . 14.86283, -6.77846,  1.45234,  0.47888, -0.42056,  0.16469,
     . 15.30276, -4.30615,  0.07564,  2.28321, -0.45717, -1.62010,
     . -1.76335,  1.03364, -0.27553,  0.34569, -0.12343, -0.10146,
     . -0.47119,  1.28997, -0.19336,  0.02724,  0.08955,  0.04726/
c
c...compute the partials of the tidal variations to the orthoweights
      call cnmtx (time, h, alpha)
c
c...compute eop changes
      do k=1,3
         eop(k) = 0.
         do j=1,12
            eop(k) = eop(k) + h(j)*orthow(j,k)
         enddo
         corx = eop(1) * 1.0d-6         !  to arcseconds
         cory = eop(2) * 1.0d-6         !  to arcseconds
         corut1 = eop(3) * 1.0d-6       !  to seconds of time
      end do
c
      return
      end



      subroutine cnmtx (dmjd, h, alpha)
c
c...Purpose: To compute the time dependent part of the second degree
c            diurnal and semidiurnal tidal potential from the dominant
c            spectral lines in the Cartwright-Tayler-Edden harmonic
c            decomposition
c
c...Coded by: Richard Eanes, UT/CSR, Feb 1997
c
c...Input: dmjd = modified julian date
c
c...Output: h = vector of length 12 with partials of the tidal
c           variation with respect to the orthoweights
c
      parameter (nlines=71)
      character*7 numarg(nlines)
      real*8 h(12), alpha, dt60, pinm , dmjd, twopi
      dimension nj(nlines),mj(nlines)
      real*8 hs(nlines),phase(nlines),freq(nlines)
      real*8 anm(2:3,0:3,-1:1), bnm(2:3,0:3,-1:1)
      real*8 p(0:2,2),q(0:2,2),sp(6,2)
c
c...the orthotide weight factors
      data ((sp(i,m),i=1,6),m=1,2) /
     . 0.0298d0, 0.1408d0, +0.0805d0, 0.6002d0, +0.3025d0, 0.1517d0,
     . 0.0200d0, 0.0905d0, +0.0638d0, 0.3476d0, +0.1645d0, 0.0923d0/
c
      data twopi /6.2831853071796d0/
      data dt / 2.d0 /
      data nmax /2/
c
c...tidal potential model for 71 diurnal and semidiurnal lines
c
      data d1960/37076.5d0/
      data (nj(j),mj(j),hs(j),phase(j),freq(j),numarg(j),j=1,15)
     ./2, 1,  -1.94, 9.0899831, 5.18688050, '117.655',            ! 29.073
     . 2, 1,  -1.25, 8.8234208, 5.38346657, '125.745',            ! 28.011
     . 2, 1,  -6.64,12.1189598, 5.38439079, '125.755',            ! 28.006
     . 2, 1,  -1.51, 1.4425700, 5.41398343, '127.545',            ! 27.853
     . 2, 1,  -8.02, 4.7381090, 5.41490765, '127.555',            ! 27.848
     . 2, 1,  -9.47, 4.4715466, 5.61149372, '135.645',            ! 26.872
     . 2, 1, -50.20, 7.7670857, 5.61241794, '135.655',            ! Q1
     . 2, 1,  -1.80,-2.9093042, 5.64201057, '137.445',            ! 26.727
     . 2, 1,  -9.54, 0.3862349, 5.64293479, '137.455',            ! 26.723
     . 2, 1,   1.52,-3.1758666, 5.83859664, '145.535',            ! 25.827
     . 2, 1, -49.45, 0.1196725, 5.83952086, '145.545',            ! 25.823
     . 2, 1,-262.21, 3.4152116, 5.84044508, '145.555',            ! O1  25.819
     . 2, 1,   1.70,12.8946194, 5.84433381, '145.755',            ! 25.802
     . 2, 1,   3.43, 5.5137686, 5.87485066, '147.555',            ! 25.668
     . 2, 1,   1.94, 6.4441883, 6.03795537, '153.655'/            ! 24.975
      data (nj(j),mj(j),hs(j),phase(j),freq(j),numarg(j),j=16,30)
     ./2, 1,   1.37,-4.2322016, 6.06754801, '155.445',            ! 24.853
     . 2, 1,   7.41,-0.9366625, 6.06847223, '155.455',            ! 24.849
     . 2, 1,  20.62, 8.5427453, 6.07236095, '155.655',            ! 24.833
     . 2, 1,   4.14,11.8382843, 6.07328517, '155.665',            ! 24.829
     . 2, 1,   3.94, 1.1618945, 6.10287781, '157.455',            ! 24.709
     . 2, 1,  -7.14, 5.9693878, 6.24878055, '162.556',            ! 24.132
     . 2, 1,   1.37,-1.2032249, 6.26505830, '163.545',            ! 24.069
     . 2, 1,-122.03, 2.0923141, 6.26598252, '163.555',            ! P1 24.066
     . 2, 1,   1.02,-1.7847596, 6.28318449, '164.554',            ! 24.000
     . 2, 1,   2.89, 8.0679449, 6.28318613, '164.556',            ! 23.999
     . 2, 1,  -7.30, 0.8953321, 6.29946388, '165.545',            ! 23.939
     . 2, 1, 368.78, 4.1908712, 6.30038810, '165.555',            ! K1 23.935
     . 2, 1,  50.01, 7.4864102, 6.30131232, '165.565',            ! 23.931
     . 2, 1,  -1.08,10.7819493, 6.30223654, '165.575',            ! 23.927
     . 2, 1,   2.93, 0.3137975, 6.31759007, '166.554'/            ! 23.869
      data (nj(j),mj(j),hs(j),phase(j),freq(j),numarg(j),j=31,45)
     ./2, 1,   5.25, 6.2894282, 6.33479368, '167.555',            ! 23.804
     . 2, 1,   3.95, 7.2198478, 6.49789839, '173.655',            ! 23.207
     . 2, 1,  20.62,-0.1610030, 6.52841524, '175.455',            ! 23.098
     . 2, 1,   4.09, 3.1345361, 6.52933946, '175.465',            ! 23.095
     . 2, 1,   3.42, 2.8679737, 6.72592553, '183.555',            ! 22.420
     . 2, 1,   1.69,-4.5128771, 6.75644239, '185.355',            ! 22.319
     . 2, 1,  11.29, 4.9665307, 6.76033111, '185.555',            ! 22.306
     . 2, 1,   7.23, 8.2620698, 6.76125533, '185.565',            ! 22.303
     . 2, 1,   1.51,11.5576089, 6.76217955, '185.575',            ! 22.300
     . 2, 1,   2.16, 0.6146566, 6.98835826, '195.455',            ! 21.578
     . 2, 1,   1.38, 3.9101957, 6.98928248, '195.465',            ! 21.575
     . 2, 2,   1.80,20.6617051,11.45675174, '225.855',            ! 13.162
     . 2, 2,   4.67,13.2808543,11.48726860, '227.655',            ! 13.127
     . 2, 2,  16.01,16.3098310,11.68477889, '235.755',            ! 12.905
     . 2, 2,  19.32, 8.9289802,11.71529575, '237.555'/            ! 12.872
      data (nj(j),mj(j),hs(j),phase(j),freq(j),numarg(j),j=46,60)
     ./2, 2,   1.30, 5.0519065,11.73249771, '238.554',            ! 12.853
     . 2, 2,  -1.02,15.8350306,11.89560406, '244.656',            ! 12.677
     . 2, 2,  -4.51, 8.6624178,11.91188181, '245.645',            ! 12.659
     . 2, 2, 120.99,11.9579569,11.91280603, '245.655',            ! 12.658
     . 2, 2,   1.13, 8.0808832,11.93000800, '246.654',            ! 12.640
     . 2, 2,  22.98, 4.5771061,11.94332289, '247.455',            ! 12.626
     . 2, 2,   1.06, 0.7000324,11.96052486, '248.454',            ! 12.608
     . 2, 2,  -1.90,14.9869335,12.11031632, '253.755',            ! 12.452
     . 2, 2,  -2.18,11.4831564,12.12363121, '254.556',            ! 12.438
     . 2, 2, -23.58, 4.3105437,12.13990896, '255.545',            ! 12.4215
     . 2, 2, 631.92, 7.6060827,12.14083318, '255.555',            ! M2  12.4206
     . 2, 2,   1.92, 3.7290090,12.15803515, '256.554',            ! 12.403
     . 2, 2,  -4.66,10.6350594,12.33834347, '263.655',            ! 12.222
     . 2, 2, -17.86, 3.2542086,12.36886033, '265.455',            ! 12.192
     . 2, 2,   4.47,12.7336164,12.37274905, '265.655'/            ! 12.188
      data (nj(j),mj(j),hs(j),phase(j),freq(j),numarg(j),j=61,71)
     ./2, 2,   1.97,16.0291555,12.37367327, '265.665',            ! 12.187
     . 2, 2,  17.20,10.1602590,12.54916865, '272.556',            ! 12.016
     . 2, 2, 294.00, 6.2831853,12.56637061, '273.555',            ! S2  12.000
     . 2, 2,  -2.46, 2.4061116,12.58357258, '274.554',            ! 11.984
     . 2, 2,  -1.02, 5.0862033,12.59985198, '275.545',            ! 11.968
     . 2, 2,  79.96, 8.3817423,12.60077620, '275.555',            ! K2  11.967
     . 2, 2,  23.83,11.6772814,12.60170041, '275.565',            ! 11.966
     . 2, 2,   2.59,14.9728205,12.60262463, '275.575',            ! 11.965
     . 2, 2,   4.47, 4.0298682,12.82880334, '285.455',            ! 11.754
     . 2, 2,   1.95, 7.3254073,12.82972756, '285.465',            ! 11.753
     . 2, 2,   1.17, 9.1574019,13.06071921, '295.555'/            ! 11.546
c
c...compute the time dependent potential matrix
c
      do k=-1,1
         dt60 = (dmjd - k*dt) - d1960
         anm(2,2,k) = 0.d0
         bnm(2,2,k) = 0.d0
         anm(2,1,k) = 0.d0
         bnm(2,1,k) = 0.d0
         do j=1,nlines
            n = nj(j)
            m = mj(j)
            pinm = mod(n+m,2)*twopi/4.d0
            alpha = phase(j) + freq(j)*dt60 - pinm
            alpha = mod(alpha,twopi)
            anm(n,m,k) = anm(n,m,k) + hs(j)*cos(alpha)
            bnm(n,m,k) = bnm(n,m,k) - hs(j)*sin(alpha)
         enddo
      enddo
c
c...orthogonalize the response terms
c
      do m = 1,2
        ap = anm(2,m,1) + anm(2,m,-1)
        am = anm(2,m,1) - anm(2,m,-1)
        bp = bnm(2,m,1) + bnm(2,m,-1)
        bm = bnm(2,m,1) - bnm(2,m,-1)
        p(0,m) = sp(1,m)*anm(2,m,0)
        p(1,m) = sp(2,m)*anm(2,m,0) - sp(3,m)*ap
        p(2,m) = sp(4,m)*anm(2,m,0) - sp(5,m)*ap + sp(6,m)*bm
        q(0,m) = sp(1,m)*bnm(2,m,0)
        q(1,m) = sp(2,m)*bnm(2,m,0) - sp(3,m)*bp
        q(2,m) = sp(4,m)*bnm(2,m,0) - sp(5,m)*bp - sp(6,m)*am
        anm(2,m,-1) = p(0,m)
        anm(2,m, 0) = p(1,m)
        anm(2,m, 1) = p(2,m)
        bnm(2,m,-1) = q(0,m)
        bnm(2,m, 0) = q(1,m)
        bnm(2,m, 1) = q(2,m)
      enddo
c
c...fill partials vector
      j = 1
      do n=2,nmax
         do m = 1,n
            do k = -1,1
               h(j)   = anm(n,m,k)
               h(j+1) = bnm(n,m,k)
               j = j + 2
             enddo
         enddo
      enddo
c
      return
      end


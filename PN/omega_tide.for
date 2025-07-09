      SUBROUTINE OMEGA_TIDE (RJD, COM)

c   This subroutine implements the Ray model for diurnal/subdirunal
c   tides.  It uses the Simon et al. Fundamental Arguments.  The
c   corrections in Omega are in rad.  This correction should be added to "average"
c   Omega value to get estimates of the instantaneous value.


      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DOUBLE PRECISION  L,        LPRIME

      HALFPI = 1.5707963267948966d0
      PI0 = 2.d0 * HALFPI 

      TWOPI = 4.d0 * HALFPI

      T = (RJD - 51544.5D0)/36525.0D0
      L = -0.00024470d0*T**4 + 0.051635d0*T**3 + 31.8792d0*T**2
     .  + 1717915923.2178d0*T + 485868.249036d0
      L = DMOD(L,1296000d0)
      LPRIME = -0.00001149d0*T**4 - 0.000136d0*T**3 - 0.5532d0*T**2
     .  + 129596581.0481d0*T + 1287104.79305d0
      LPRIME = DMOD(LPRIME,1296000d0)
      CAPF = 0.00000417d0*T**4 - 0.001037d0*T**3 - 12.7512d0*T**2
     .  + 1739527262.8478d0*T + 335779.526232d0
      CAPF = DMOD(CAPF,1296000d0)
      CAPD = -0.00003169d0*T**4 + 0.006593d0*T**3 - 6.3706d0*T**2
     .  + 1602961601.2090d0*T + 1072260.70369d0
      CAPD = DMOD(CAPD,1296000d0)
      OMEGA = -0.00005939d0*T**4 + 0.007702d0*T**3 + 7.4722d0*T**2
     .  - 6962890.2665d0*T + 450160.398036d0
      OMEGA = DMOD(OMEGA,1296000d0)
      THETA = (67310.54841d0 +
     .        (876600d0*3600d0 + 8640184.812866d0)*T +
     .         0.093104d0*T**2 -
     .         6.2d-6*T**3)*15.0d0 + 648000.0d0



      ARG7 = DMOD((-L - 2.0D0*CAPF - 2.0D0*OMEGA + THETA)             ! Q1
     .     * PI0/648000.0D0,TWOPI) - HALFPI
      ARG1 = DMOD((-2.0d0*CAPF - 2.0d0*OMEGA + THETA)                 ! O1
     .     * PI0/648000.0D0,TWOPI) - HALFPI
      ARG2 = DMOD((-2.0d0*CAPF + 2.0d0*CAPD - 2.0d0*OMEGA + THETA)    ! P1
     .     * PI0/648000.0D0,TWOPI) - HALFPI
      ARG3 = DMOD(THETA * PI0/648000.0D0,TWOPI)   + HALFPI            ! K1
     

      ARG4 = DMOD((-L - 2.0d0*CAPF - 2.0D0*OMEGA + 2.0d0*THETA)       ! N2
     .     * PI0/648000.0D0,TWOPI)
      ARG5 = DMOD((-2.0D0*CAPF - 2.0D0*OMEGA + 2.0d0*THETA)           ! M2
     .     * PI0/648000.0D0,TWOPI)
      ARG6 = DMOD((-2.0d0*CAPF + 2.0d0*CAPD - 2.0d0*OMEGA + 2.0d0*THETA)  ! S2
     .     * PI0/648000.0D0,TWOPI)
      ARG8 = DMOD((2.0d0*THETA)     * PI0/648000.0D0,TWOPI)               ! K2 

      COM =     12.d0*DSIN(ARG7) - 24.D0*DCOS(ARG7)
     .       +  60.D0*DSIN(ARG1) - 79.D0*DCOS(ARG1)
     .       +  15.D0*DSIN(ARG2) - 27.D0*DCOS(ARG2)
     .       +  46.D0*DSIN(ARG3) - 94.D0*DCOS(ARG3)
     .       -  38.D0*DSIN(ARG4) + 16.D0*DCOS(ARG4)
     .       - 166.D0*DSIN(ARG5) + 74.D0*DCOS(ARG5)
     .       -  81.D0*DSIN(ARG6) +  4.D0*DCOS(ARG6)
     .       -  21.D0*DSIN(ARG8) +  4.D0*DCOS(ARG8)
 
      COM = COM * 1.0d-15  !  rad/s
      

      RETURN
      END




      SUBROUTINE RAY (RJD,CORX,CORY,CORT)
c   This subroutine implements the Ray model for diurnal/subdirunal
c   tides.  It uses the Simon et al. Fundamental Arguments.  The
c   corrections in x and y are in units of sec. of arc and UT1-UTC
c   in sec. of time.  These corrections should be added to "average"
c   EOP values to get estimates of the instantaneous values.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DOUBLE PRECISION  L,        LPRIME
      HALFPI = 1.5707963267948966d0

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
     .     * 3.14159265D0/648000.0D0,6.28318530718D0) - HALFPI
      ARG1 = DMOD((-2.0d0*CAPF - 2.0d0*OMEGA + THETA)                 ! O1
     .     * 3.14159265D0/648000.0D0,6.28318530718D0) - HALFPI
      ARG2 = DMOD((-2.0d0*CAPF + 2.0d0*CAPD - 2.0d0*OMEGA + THETA)    ! P1
     .     * 3.14159265D0/648000.0D0,6.28318530718D0) - HALFPI
      ARG3 = DMOD(THETA * 3.14159265D0/648000.0D0,6.28318530718D0)    ! K1
     .     + HALFPI
      ARG4 = DMOD((-L - 2.0d0*CAPF - 2.0D0*OMEGA + 2.0d0*THETA)       ! N2
     .     * 3.14159265D0/648000.0D0,6.28318530718D0)
      ARG5 = DMOD((-2.0D0*CAPF - 2.0D0*OMEGA + 2.0d0*THETA)           ! M2
     .     * 3.14159265D0/648000.0D0,6.28318530718D0)
      ARG6 = DMOD((-2.0d0*CAPF + 2.0d0*CAPD - 2.0d0*OMEGA + 2.0d0*THETA)  ! S2
     .     * 3.14159265D0/648000.0D0,6.28318530718D0)
      ARG8 = DMOD((2.0d0*THETA)                                           ! K2
     .     * 3.14159265D0/648000.0D0,6.28318530718D0)
      CORX = - 0.026D0*DSIN(ARG7) + 0.006D0*DCOS(ARG7)
     .       - 0.133D0*DSIN(ARG1) + 0.049D0*DCOS(ARG1)
     .       - 0.050D0*DSIN(ARG2) + 0.025D0*DCOS(ARG2)
     .       - 0.152D0*DSIN(ARG3) + 0.078D0*DCOS(ARG3)
     .       - 0.057D0*DSIN(ARG4) - 0.013D0*DCOS(ARG4)
     .       - 0.330D0*DSIN(ARG5) - 0.028D0*DCOS(ARG5)
     .       - 0.145D0*DSIN(ARG6) + 0.064D0*DCOS(ARG6)
     .       - 0.036D0*DSIN(ARG8) + 0.017D0*DCOS(ARG8)
      CORY = - 0.006D0*DSIN(ARG7) - 0.026D0*DCOS(ARG7)
     .       - 0.049D0*DSIN(ARG1) - 0.133D0*DCOS(ARG1)
     .       - 0.025D0*DSIN(ARG2) - 0.050D0*DCOS(ARG2)
     .       - 0.078D0*DSIN(ARG3) - 0.152D0*DCOS(ARG3)
     .       + 0.011D0*DSIN(ARG4) + 0.033D0*DCOS(ARG4)
     .       + 0.037D0*DSIN(ARG5) + 0.196D0*DCOS(ARG5)
     .       + 0.059D0*DSIN(ARG6) + 0.087D0*DCOS(ARG6)
     .       + 0.018D0*DSIN(ARG8) + 0.022D0*DCOS(ARG8)
      CORT = + 0.0245D0*DSIN(ARG7) + 0.0503D0*DCOS(ARG7)
     .       + 0.1210D0*DSIN(ARG1) + 0.1605D0*DCOS(ARG1)
     .       + 0.0286D0*DSIN(ARG2) + 0.0516D0*DCOS(ARG2)
     .       + 0.0864D0*DSIN(ARG3) + 0.1771D0*DCOS(ARG3)
     .       - 0.0380D0*DSIN(ARG4) - 0.0154D0*DCOS(ARG4)
     .       - 0.1617D0*DSIN(ARG5) - 0.0720D0*DCOS(ARG5)
     .       - 0.0759D0*DSIN(ARG6) - 0.0004D0*DCOS(ARG6)
     .       - 0.0196D0*DSIN(ARG8) - 0.0038D0*DCOS(ARG8)
      CORX = CORX * 1.0d-3
      CORY = CORY * 1.0d-3
      CORT = CORT * 0.1d-3
      RETURN
      END




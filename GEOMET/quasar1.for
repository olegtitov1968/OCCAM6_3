C
C  Subroutine QUASAR1
C
C--------------------------------------------------------------------
C
C  Purpose: calculate theoretical values for VLBI observables
C           when observing infinitely distant objects.
C
C           This subroutine is based on the VLBI model
C             by Sergei A. Klioner published in Proc. of the
C             AGU Chapman Conference on Geodetic VLBI:
C             Monitoring Global Change, NOAA Technical Report
C             NOS 137 NGS 49, 188-202.
C           The accuracy of observables calculated by the subroutine
C             is 1 ps in time delay and 1 fs/s in delay rate.
C           We suppose herewith that
C
C             1) the observers are situated on the Earth's surface
C                (not on the Earth satellite)
C             2) geocentric reference system is kinematically
C                non-rotating (according to the recommendations
C                of IAU, 1991)
C             3) the source is infinitely distant (neither parallax
C                nor proper motion is accounted for)
C             4) the angular distances between the source and
C                gravitating bodies are not less than the values
C                specified below:
C
C                      Sun         1 degree
C                      Mercury     2'
C                      Venus      55'
C                      Mars        5'
C                      Jupiter    23 degrees  3'
C                      Saturn      3 degrees 27'
C                      Uran       15'
C                      Neptune    11'
C                      Moon       42'
C
C
C  INPUT PARAMETERS:
C
C
C    EARTH  =.true.  Gravitational term due to the Earth is applied.
C           =.false. The term is neglected.
C
C    PPN    =.true.  The post-post-Newtonian term due to the Sun
C                      is applied.
C           =.false. The term is neglected.
C
C    RATE   =.true.  Delay rate will be calculated analytically.
C           =.false. Delay rate will not be calculated (input
C                      parameters dtau and dtaugr will not be
C                      affected).
C
C    VARK   =.true.  When calculating delay rate the derivative of
C                    the unit vector directed to the source is
C                    considered. Some VLBI packages (e.g. OCCAM,
C                    Kiev-GR1) account for precession and nutation
C                    by computing the position of the source at
C                    the moment of observation. In this case
C                    the derivative of the unit vector to the source
C                    due to precession and nutation must be considered.
C           =.false. The position of the source is supposed to be
C                    constant.
C
C    SCALEDX =.true.  The spatial coordinates are scaled x=x_GCRS*(1-L_G),
C                     where x_GCRS are thje spatial coordinates of the IAU.
C                     This choice corresponds to the model from IERS Standards
C                     (1992), but contradicts the IAU resolutions
C
C            =.false. The spatial coodinates are those of the GCRS of the IAU
C                     without any additional scaling
C
C             NOTE that this subroutine implements a simplified version
C              of the scaling which neglects terms of order LG^2=10^(-18),
C              which is quite reasonable for both delay and delay rate
C              for Earth-born baselines and the accuracies specified above.
C
C
C    K      unit vector from the barycenter to the source
C (*)DK     derivative of the vector K (see above)
C
C    W1     geocentric coordinates of the first station
C (*)V1     geocentric velocity of the first station
C    W2     geocentric coordinates of the second station
C    V2     geocentric velocity of the second station
C (*)A2     geocentric acceleration of the second station
C
C    XE     barycentric coordinates of the Earth
C    VE     barycentric velocity of the Earth
C (*)AE     barycentric acceleration of the Earth
C    XS     barycentric coordinates of the Sun
C
C           All these quantities must be calculated at the moment of
C             reception of the signal at the first station.
C           The star (*) means that the parameter is used only when
C             delay rate is computed analytically (RATE=.true.)
C
C    GMS    G * the mass of the Sun
C    GME    G * the barycentric mass of the Earth
C             (geocentric value can also be used at the adopted
C              level of accuracy)
C
C    C      the light velocity
C
C
C  OUTPUT PARAMETERS:
C
C    tau    theoretical time delay (including gravitational part)
C    taugr  theoretical gravitational time delay
C    dtau   theoretical delay rate
C    dtaugr theoretical gravitational delay rate
C    pgamm  partial derivative for gamma
C
C
C  EXTERNAL SUBROUTINES:
C
C    SLOG(X)   double precision function which compute DLOG(1+X)
C
C         Revision: 30 March 2001
C  Programmer:      Sergei A. Klioner
C  Latest revision: 05 October 2006 (OT) - partial derivative for gamma
C                    Variable EARTH deleted
C--------------------------------------------------------------------
C



      SUBROUTINE quasar1 (tau,taugr,dtau,dtaugr,PPN,RATE,VARK,
     ,                  SCALEDX,
     ,         K,DK,W1,V1,W2,V2,A2,XE,VE,AE,XS,GMS,GME,C,pgamm,teta)
c
c    Added on 17.12.2007 (OT)
c
    
      IMPLICIT none

C
C Input parameters:
C

      LOGICAL PPN,RATE,VARK,SCALEDX
      DOUBLE PRECISION K(3),DK(3),W1(3),V1(3),W2(3),V2(3),A2(3),
     ,                 XE(3),VE(3),AE(3),XS(3)
      DOUBLE PRECISION GMS,GME,C

C
C Output parameters:
C

      DOUBLE PRECISION tau,taugr,dtau,dtaugr
      DOUBLE PRECISION pgamm, pgamm1, du1, teta, cos_teta      ! 05-Oct-2006    (OT)
  

C
C Local variables:
C
C Index
C

      INTEGER i

C
C B       is geocentric baseline vector
C Bbar    is barycentric baseline vector
C xi      are barycentric coordinates of i-th station
C dxi     are barycentric velocities of i-th station minus
C             the Earth's velocity at the moment t_1
C VE2     is velocity of the Earth at the moment t_2
C V22     is velocity of the second station at the moment u_2
C W22     is coordinates of the second station at the moment u_2
C

      DOUBLE PRECISION B(3),Bbar(3),dx1(3),dx2(3),x1(3),x2(3),
     ,                 VE2(3),V22(3),W22(3)

C
C additional constant: L_G=6.969290134e-10
C

      DOUBLE PRECISION LG

C
C auxiliary variables
C

      DOUBLE PRECISION Cinv,Cinv2,Cinv3,
     ,                 aux1,aux2,aux3,aux4

C
C dt     is the difference of barycentric coordinate moments
C           corresponding to the receptions of the signal by the
C           second and first stations
C du     is the analogous difference of geocentric coordinate moments
C tgrav  is gravitational time delay without aberrational term
C dtdui  are partial derivatives dt_i/du_i
C dtgri  are partial derivatives dt_{gr}/dt_i
C dt21   is partial derivative dt_2/dt_1
C dt21gr is the part of dt21 induced by gravitational time delay
C
 
      DOUBLE PRECISION dt,du,tgrav,
     ,                 dtdu1,dtdu2,dtgr1,dtgr2,dt21,dt21gr

C
C The letter "a" at the beginning of the identifier means that the
C   corresponding local variable is absolute value of the vector
C   whose identifier coincides with the rest of the identifier
C   in question. The variables whose identifies begin with "a2" are
C   square of the relevant absolute values.
C The symbol "_" in the middle of the identifier means that the
C   variable is a scalar product of the vectors whose names are
C   situated on both sides from "_".
C
C ES     is the distance from the Sun to the geocenter
C U      is Newtonian potential of the Sun evaluated at the
C        geocenter
C rSi    are the vectors from the Sun to the i-th station
C K_dx12=K_dx1-K_dx2
C

      DOUBLE PRECISION arS1,arS2,aW1,aW2,aW22,a2VE,a2VE2,
     ,                 K_B,K_Bbar,K_dx12,K_dx2,K_rS1,K_rS2,K_V,
     ,                 K_W1,K_W2,
     ,                 B_VE,VE_V2,VE_V1,VE_W1,VE_W2,DK_B,
     ,                 ES,U

C
C External function
C

      DOUBLE PRECISION SLOG

C
C initialize tgrav (workaround for the bug with Earth .eq. .false.)
C

      tgrav=0.0d0

C
C value of the scaling constant L_G
C

      LG=6.969290134d-10
C
C Calculate auxiliary constants
C

      Cinv=1.0d0/C
      Cinv2=Cinv*Cinv
      Cinv3=Cinv2*Cinv

C
C Calculate geocentric baseline vector
C

      DO 5 i=1,3
         B(i)=W2(i)-W1(i)
   5  CONTINUE

C
C Calculate barycentric baseline vector: Eq.(4.15)
C
C           Compute several scalar products:
C

      B_VE=0.d0
      ES=0.d0

      DO 10 i=1,3
         B_VE=B_VE+B(i)*VE(i)
         ES=ES+(XE(i)-XS(i))**2
  10  CONTINUE

      ES=dsqrt(ES)

C
C           U is the Newtonian potential of the Sun evaluated
C             at the geocenter
C

      U=GMS/ES

C
C           Evaluate Eq.(4.15)
C

c      print *, ' ok ' 

      DO 15 i=1,3
         Bbar(i)=B(i)*(1.d0-Cinv2*U)-
     -           Cinv2*B_VE*(0.5d0*VE(i)+V2(i))
  15  CONTINUE

C
C
C Calculate gravitational time delay: Eqs.(4.8)-(4.12)
C
C           The post-Newtonian term due to the Earth: Eq.(4.9)
C
C                         Geocentric distances of the stations
C

      aW1=0.0d0
      aW2=0.0d0

      DO 20 i=1,3
         aW1=aW1+W1(i)**2
         aW2=aW2+W2(i)**2
  20  CONTINUE

      aW1=dsqrt(aW1)
      aW2=dsqrt(aW2)

      

C
C                       Scalar products K*Wi
C
      K_W1=0.0d0
      K_W2=0.0d0

      DO 25 i=1,3
         K_W1=K_W1+K(i)*W1(i)
         K_W2=K_W2+K(i)*W2(i)
  25  CONTINUE

C
C                         Compute the term...
C

      tgrav=tgrav+2.d0*GME*Cinv3*
     *      SLOG((K_W1/aW1-K_W2/aW2)/(1.0d0+K_W2/aW2))


        
C
C           The post-Newtonian term due to the Sun
C
C                         Barycentric coordinates of the stations:
C                          Eq.(2.2)
C

      VE_W1=0.0d0
      VE_W2=0.0d0

      DO 30 i=1,3
         VE_W1=VE_W1+VE(i)*W1(i)
         VE_W2=VE_W2+VE(i)*W2(i)
  30  CONTINUE

      DO 35 i=1,3
         x1(i)=XE(i)+W1(i)*(1.d0-Cinv2*U)-0.5d0*Cinv2*VE_W1*VE(i)
         x2(i)=XE(i)+W2(i)*(1.d0-Cinv2*U)-0.5d0*Cinv2*VE_W2*VE(i)
  35  CONTINUE

C
C                         Compute the terms under logarithm
C                          in Eq.(4.9)
C

      arS1=0.d0
      arS2=0.d0
      K_rS1=0.d0
      K_rS2=0.d0

      DO 40 i=1,3
         aux1=x1(i)-XS(i)
         aux2=x2(i)-XS(i)
         arS1 = arS1 + aux1**2
         arS2 = arS2 + aux2**2
         K_rS1= K_rS1+K(i)*aux1
         K_rS2= K_rS2+K(i)*aux2        ! for teta
  40  CONTINUE

      arS1=dsqrt(arS1)
      arS2=dsqrt(arS2)
   
      cos_teta = -K_rS2/arS2
      teta = dacos(cos_teta)

C
C                         Compute the term using Eqs.(4.9) and (5.3)
C                            to avoid the loss of accuracy
C
      aux1=0.d0
      aux2=arS1+arS2
      DO 45 i=1,3
         aux1=aux1-Bbar(i)*((x1(i)+x2(i)-2.0d0*XS(i))/aux2+K(i))
  45  CONTINUE

      tgrav=tgrav+ 2.d0*GMS*Cinv3 * SLOG(aux1/(arS2+K_rS2))

C
C           The post-post-Newtonian term due to the Sun
C             (4.12)
C

      IF (PPN) THEN

         aux2=-4.d0*aux1/(arS1+K_rS1)/(arS2+K_rS2)+
     +        0.25d0*(K_rS2/arS2**2-K_rS1/arS1**2)+
     +        3.75d0*(dacos(K_rS2/arS2)/dsqrt(arS2**2-K_rS2**2)-
     -                dacos(K_rS1/arS1)/dsqrt(arS1**2-K_rS1**2))
         tgrav=tgrav+ Cinv*(GMS*Cinv2)**2*aux2

      END IF

C
C Barycentric coordinate time delay: Eq.(4.2)
C

      K_Bbar=0.d0
      K_V  =0.d0

      DO 50 i=1,3
         K_Bbar=K_Bbar+K(i)*Bbar(i)
         K_V=K_V+K(i)*(VE(i)+V2(i))
  50  CONTINUE

      aux1=K_V*Cinv

      pgamm = 0.5d0*tgrav/(1.d0+aux1)    !   partial derivative for gamma
      
      aux2=aux1*aux1

      taugr=tgrav*(1.0d0-aux1)

      dt=-Cinv*K_Bbar*(1.0d0-aux1+aux2)+taugr

c      dt=-Cinv*K_Bbar*(1.0d0-aux1+aux2)


C
C Geocentric coordinate time delay: Eq.(4.18)
C

      K_B=0.0d0
      aux1=0.0d0

      DO 55 i=1,3
         K_B=K_B+K(i)*B(i)
         aux1=aux1+VE(i)*(0.5d0*VE(i)+V2(i))
  55  CONTINUE



c      du1=dt-Cinv2*B_VE+Cinv3*K_B*aux1

      pgamm1 = pgamm + Cinv3*K_B*U/(1.d0 + K_V*Cinv)      !   partial derivative for gamma    18-Jun-2015

c      print *, pgamm1, pgamm, Cinv3*K_B*U/(1.d0 +K_V*Cinv)

      du=dt-Cinv2*B_VE+Cinv3*K_B*(aux1+U)

c      print *, du, du1, du - du1



C
C Proper time delay expressed in TT.
C

      IF (SCALEDX) THEN
C
C The geocentric time scale is TT (not TCG), the baseline is measured
C   in TT meters. Therefore, scaling term is not needed. And observed
C   time delay is simply
C
      tau=du

      ELSE

C
C The delay is expressed in TT, but the spatial coordinates in
C   the non-scaled coodinates of the Geocentric Celestial Reference
C   System (GCRS) of the IAU.
C
      tau=du*(1-LG)
      taugr=taugr*(1-LG)

      ENDIF

C
C Compute observed delay rate analytically if requested
C

      IF (RATE) THEN

C
C Compute barycentric velocities of the i-th station at the
C   moment t_i: Eq.(4.16)
C
C           The first station. To avoid loss of accuracy when
C             subtracting the velocities of the stations we
C             will not add the velocity of the Earth at the
C             moment t_1 to both velocities.
C
      a2VE=0.0d0
      VE_V1=0.0d0

      DO 60 i=1,3
         a2VE=a2VE+VE(i)**2
         VE_V1=VE_V1+VE(i)*V1(i)
  60  CONTINUE

      DO 65 i=1,3
         dx1(i)=-VE(i)*Cinv2*0.5d0*VE_V1+
     +          V1(i)*(1.0d0-Cinv2*(0.5d0*a2VE+2.0d0*U+VE_V1))
  65  CONTINUE

C
C           The second station.
C           We need the velocity at the moment t_2.
C           Therefore, we need to compute the position and velocity
C           of the second station at the moment t_2, and the velocity
C           of the Earth at the same moment. We can neglect
C           herewith the variability of the accelerations of the
C           station and the Earth.
C

      DO 70 i=1,3
         V22(i)=V2(i)+du*A2(i)
         W22(i)=W2(i)+du*(V2(i)+0.5d0*du*A2(i))
         VE2(i)=VE(i)+dt*AE(i)
  70  CONTINUE

C
C           The variation of the potential of the Sun
C             can be neglected.
C

      a2VE2=0.0d0
      VE_V2=0.0d0

      DO 75 i=1,3
         a2VE2=a2VE2+VE2(i)**2
         VE_V2=VE_V2+VE2(i)*V22(i)
  75  CONTINUE

      DO 80 i=1,3
         dx2(i)=-VE2(i)*Cinv2*0.5d0*VE_V2+
     +          dt*AE(i)+
     +          V22(i)*(1.0d0-Cinv2*(0.5d0*a2VE2+2.0d0*U+VE_V2))
  80  CONTINUE

C
C Compute partial derivatives of post-Newtonian gravitational
C   time delay: Eq.(4.13)
C   The velocity of the Sun is neglected.
C

      aW22=0.0d0

      DO 85 i=1,3
         aW22=aW22+W22(i)**2
  85  CONTINUE

      aW22=dsqrt(aW22)

      aux1=0.0d0
      aux2=0.0d0
      aux3=0.0d0
      aux4=0.0d0

      DO 90 i=1,3
         aux1=aux1+((x1(i)-XS(i))/arS1+K(i))*(dx1(i)+VE(i))
         aux2=aux2+(W1(i)/aW1+K(i))*V1(i)
         aux3=aux3+((x2(i)-XS(i))/arS2+K(i))*(dx2(i)+VE(i))
         aux4=aux4+(W22(i)/aW22+K(i))*V22(i)
  90  CONTINUE

      dtgr1= 2.0d0*Cinv3*(GMS*aux1/(arS1+K_rS1)+
     +                    GME*aux2/(aW1+K_W1))
      dtgr2=-2.0d0*Cinv3*(GMS*aux3/(arS2+K_rS2)+
     +                    GME*aux4/(aW2+K_W2))

C
C Compute derivative of the moment t_2 w.r.t. t_1 minus 1:
C   Eq.(4.7)
C

      K_dx12=0.0d0
      K_dx2=0.0d0

      DO 95 i=1,3
         K_dx12=K_dx12+K(i)*(dx1(i)-dx2(i))
         K_dx2=K_dx2+K(i)*(dx2(i)+VE(i))
  95  CONTINUE

      dt21gr=(dtgr1+dtgr2)*(1.0d0-Cinv*K_dx2)+
     +        dtgr2*Cinv*K_dx12

      dt21=Cinv *K_dx12*
     *     (1.0d0-Cinv*K_dx2*(1.0d0-Cinv*K_dx2))+dt21gr

C
C Compute delay rate due to possible variability of
C   the direction to the source (see description of
C   the parameter input VARK).
C

      IF (VARK) THEN

      DK_B=0.0d0

      DO 100 i=1,3
         DK_B=DK_B+DK(i)*(Bbar(i)+dt*(dx2(i)+VE(i)))
 100  CONTINUE

      dt21=dt21-Cinv*DK_B

      END IF

C
C Compute derivatives of the moments t_i w.r.t. u_i minus 1:
C   Eq.(4.22)
C

      aux1=0.0d0
      aux2=0.0d0
      aux3=0.0d0
      aux4=0.0d0

      DO 105 i=1,3

         aux1=aux1+AE(i)*W1(i)
         aux2=aux2+VE(i)*V1(i)
         aux3=aux3+AE(i)*W22(i)
         aux4=aux4+VE2(i)*V22(i)

 105  CONTINUE

      dtdu1=Cinv2*(0.5d0*a2VE +U+aux1+aux2)
      dtdu2=Cinv2*(0.5d0*a2VE2+U+aux3+aux4)

C
C Compute observable delay rate.
C   The term due to the difference in the heights of the
C     stations =Cinv2 * 9.8 m/c * (h2-h1) is neglected.
C
C                      Compute the terms in Eq.(4.21)
C                        due to dt_i/du_i
C

      aux1=(dtdu1-dtdu2)/(1.0d0+dtdu2)

      IF (SCALEDX) THEN
C
C The geocentric time scale is TT (not TCG), the baseline is measured
C   in TT meters. Therefore, scaling term is not needed. And observed
C   time delay is simply
C
C
C no scaling is necessary here
C                      Eq.(4.21)
C
      dtaugr=dt21gr+dt21gr*aux1
      dtau  =dt21+aux1+dt21*aux1

c      print *, dtau, dtaugr

      ELSE

C
C The delay is expressed in TT, but the spatial coordinates in
C   the non-scaled coodinates of the Geocentric Celestial Reference
C   System (GCRS) of the IAU.
C
      dtau=(dt21+aux1+dt21*aux1)*(1-LG)
      dtaugr=(dt21gr+dt21gr*aux1)*(1-LG)

      END IF

      END IF
      RETURN
      END

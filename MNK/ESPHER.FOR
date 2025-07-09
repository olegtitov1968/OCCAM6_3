***********************************************************************
C
C    SUBROUTINE ESPHER
C
C    THIS ROUTINE COMPUTES THE FINAL PARAMETER VALUES IN ESPHERICAL
C      COORDINATES AND THEIR STANDARD DEVIATIONS
C
C    REVISION APRIL 18, 1991
C    Revision July, 22, 1998, (OT) The number of parameters has been
C                                      increased
C    REVISION OCTOBER 04, 1999 (OT)
C    REVISION AUGUST 13, 2001 (OT)  Error with DATAN fixed
C    REVISION AUGUST 22, 2001 (OT)  IGRAD parameter aded
C    LAST REVISION January 25, 2006 (OT)  NS, NMAX and NPAR called from OCCAM.FI
C
C**********************************************************************

      SUBROUTINE ESPHER (I,COVM,CX,R,DPHI,DDLAM,SR,SPHI,SLAM)

C   ESPECIFICATIONS

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'OCCAM.FI'

      DIMENSION CX(nmax,3),A(3),CV(6),COVM(npar*(npar+1)/2)


      COMMON /PHYS/ C, FL, AH, AU, J20
      COMMON /MATH/ PI, TWOPI, ROG, ROS, ROZ, TOL

C  GET ESPHERICAL COORDINATES

      R = DSQRT(CX(I,1)**2+CX(I,2)**2+CX(I,3)**2)     !  [m]

      EQR = DSQRT (CX(I,1)**2+CX(I,2)**2)

      PHI = DATAN2(CX(I,3),EQR)
      DPHI = PHI*AH               !  [m]

      DLAM = DATAN (CX(i,2)/CX(i,1))

      phi1=phi*2.d4
      iphi1=int(phi1)
      phi2=dble(iphi1)/2.d4


      DDLAM = DLAM*AH*DCOS(PHI2)  !  [m]

C  GET SIGMAS OF SPHERICAL COORDINATES

      SR = 0.D0
      SPHI = 0.D0
      SLAM = 0.D0

      IF (I.EQ.1) then
         do j=1,3
            do k=1,j
               IJK = J*(J-1)/2+K
               CV(IJK) = COVM(IJK)
            end do
         end do
      ELSE
         IPAR = (I-2)*ns+(ns-3)
         DO J=1,3
            DO K=1,J
               IPS = (IPAR+J-1)*(IPAR+J)/2 + IPAR+K
               IJK = J*(J-1)/2+K
               CV(IJK) = COVM(IPS)
            ENDDO
         ENDDO
      END IF

C    SIGMA OF R

      DO J=1,3
         A(J) = CX(I,J)/R
      ENDDO

      CALL ABAT (A,3,1,CV,1,1,SR,NRX,IER)
      SR = DSQRT(SR)

C    SIGMA OF PHI

      DO K=1,2
         A(K) = -CX(I,3)*CX(I,K)/(R**2*EQR)
      ENDDO
      A(3) = EQR/R**2

      CALL ABAT (A,3,1,CV,1,1,SPHI,NRX,IER)
      SPHI = DSQRT(SPHI)*AH               ! [m]

C    SIGMA OF LAM

      A(1) = -CX(I,2)/EQR**2
      A(2) =  CX(I,1)/EQR**2
      A(3) = 0.D0

      CALL ABAT (A,3,1,CV,1,1,SLAM,NRX,IER)
      SLAM = DSQRT(SLAM)*AH*DCOS(PHI)       ! [m]

      END




      SUBROUTINE FEED (STAT,AXTYP,CORZ,HSTA,CPHI,Z,DE,FEED_R)

C************************************************************************
C
C   THIS ROUTINE COMPUTES THE CORRECTIONS TO THE DELAY RATE DUE TO THE FEED ROTATION
C
C       15-SEP-2019 (OT)
C
C************************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      CHARACTER AXTYP*4,STAT*8,STT*1
      DIMENSION Z(3),DE(3), CPHI(3), CORZ(3)
      DIMENSION HSTA(3)
      COMMON /PHYS/ C, FL, AH, AU, J20
      COMMON /MATH/ PI, TWOPI, ROG, ROS, ROZ, TOL

c
c     For delay calculation only IK=3 is necessary
c     For delay rate calculation only IK=1,2 is necessary
c



      DO IK=1,3

         SZ = DSIN (Z(IK)-CORZ(IK))
         CZ = DCOS (Z(IK)-CORZ(IK))
         CD = DCOS (DE(IK))
         SD = DSIN (DE(IK))

C  TREATING THE CASE OF AN AZIMUTH-ELEVATION MOUNTING
C----------------------------------------------------------

         IF (AXTYP.EQ.'AZEL')  THEN

             A1 = DCOS(cphi(ik))*DSIN(hsta(ik))
             A2 = DSIN(cphi(ik))*dcos(de(ik)) - 
     *             DCOS(cphi(ik))*dsin(de(ik)) *DCOS(hsta(ik))

             TAN_FEED = A1/A2
             FEED0 = DATAN(TAN_FEED)

C  TREATING THE CASE OF AN EQUATORIAL ANTENNA
C----------------------------------------------------------

         ELSE IF (AXTYP.EQ.'EQUA') THEN

             FEED0 = 0.d0


C  TREATING THE CASE OF AN X-Y ANTENNA WITH FIXED AXIS NORTH-SOUTH
C-----------------------------------------------------------------

         ELSE IF (AXTYP.EQ.'X-Y1'.OR.AXTYP.EQ.'X-YN') THEN

             A1 = -DSIN(cphi(ik))*DSIN(hsta(ik))
             A2 = DCOS(cphi(ik))*dcos(de(ik)) +
     *             DSIN(cphi(ik))*dsin(de(ik)) *DCOS(hsta(ik))

             TAN_FEED = A1/A2
             FEED0 = DATAN(TAN_FEED)

            

C  TREATING THE CASE OF AN X-Y ANTENNA WITH FIXED AXIS EAST-WEST
C---------------------------------------------------------------

         ELSE IF (AXTYP.EQ.'X-Y2'.OR.AXTYP.EQ.'X-YE') THEN

             A1 = DCOS(hsta(ik))
             A2 = dsin(de(ik)) *Dsin(hsta(ik))

             TAN_FEED = A1/A2
             FEED0 = DATAN(TAN_FEED)
     

         ELSE

            WRITE (*,321) AXTYP,STAT



321   FORMAT (//' ***** WARNING *****  AXIS MODEL ',A4,' FOR STATION',
     *1X,A8,' NOT AVAILABLE'//' HIT "Q" TO QUIT, "C" TO CONTINUE'//)
            READ (*,'(A1)') STT
            IF (STT.NE.'C'.AND.STT.NE.'c') STOP
         ENDIF

C AXIS MODEL FOR RICHMOND
C-----------------------------

         
         IF (IK.EQ.1) FEED1 = FEED0/8.4d9       !   for phase delay [sec]
         IF (IK.EQ.2) FEED2 = FEED0/8.4d9       !   for phase delay [sec]

      end do
         
      FEED_R = (FEED2 - FEED1) / 2.D0    !   correction to the delay rate [rad]

      RETURN
      END
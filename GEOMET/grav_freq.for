      SUBROUTINE GRAV_FREQ (R1,R2,S,VE,B,Jmap,GM_Body, 
     *   RB, RBody, DV2V1, Freq_Grav, Freq_1, Freq_2)                                                        


CCORRECTION DUE TO RELATIVISTIC LIGHT BENDING (TGRAV) for Sun, Moon and planets
C-----------------------------------------------------------------------
C THIS IS THE ALGORITHM developed in 2013
C
c   Last Revision: 12 April 2021 (OT)
C
C   Copyright (C) Oleg Titov, Geoscience Australia, 2001-2021
C
C   This program is free software: you can redistribute it and/or modify it
C   under the terms of the GNU General Public License as published by the
C   Free Software Foundation, either version 3 of the License, or (at your
C   option) any later version.
C
C   This program is distributed in the hope that it will be useful, but
C   WITHOUT ANY WARRANTY; without even the implied WARRANTY of 
C   MERCHANTABILITY of FITNESS FOR A PARTICULAR PURPOSE.
C   See the GNU General Public License for more details.
C
C   You should have received a copy of the GNU General Public License
C   along with the program. If not, see <http://www.gnu.org/license/>.
C
C
C
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      IMPLICIT DOUBLE PRECISION  (A-H,O-Z)

      COMMON /PHYS/ C, FL, AH, AU, J20
      COMMON /MATH/ PI, TWOPI, ROG, ROS, ROZ, TOL

      DIMENSION S(3), R1(3), R2(3), RB(3), B(3), B2(3), VE(3)
      DIMENSION R2MRBody(3), R1MRBody(3), R2MRBodyV(3), DV2V1(3)
      DIMENSION  RN2(3), SSN(3)

      CHARACTER SOURCE*8
    

      CINV  = 1.D0 / C
      CINV2 = CINV * CINV
      CINV3 = CINV * CINV2

c       R1, R2 - geocentic positions of the first and second sites
c       RB - geocentric position of the gravity body
c       

           
      R2_Dist = 0.d0
      BETA = 0.d0
      phi = 0.d0
      RBod2 = 0.d0  

      phi_V = 0.d0
     
      SIN_TETA = 0.d0  
      COS_TETA = 0.d0
      SIN_PHI = 0.d0
      TETA = 0.d0
      

      COS_A =0.d0
      COS_PHI = 0.d0
      COS_PSI = 0.d0
     

      COS_A_V =0.d0
      COS_PHI_V = 0.d0
      COS_PSI_V = 0.d0
     

      Freq_1 = 0.d0
      Freq_2 = 0.d0
      Freq_Grav = 0.d0

   

      
      do i=1,3
    
   

        R2MRBody(I) = R2(I) - RB(I)    ! 6-Oct-06
        R1MRBody(I) = R1(I) - RB(I)    ! 6-Oct-06
      
        
      end do

      VSUN = DSQRT(DOTPR(VE,VE))  !  Barycentric velocity of the Earth

 

      BB = DSQRT (DOTPR(B,B))    !  Baseline length

      VV2 = DSQRT(DOTPR(DV2V1, DV2V1))  ! lenght of the w2-w1 vector

      freq_b = - DOTPR(DV2V1,S)/C 


      RBod2 = DSQRT (DOTPR(R2MRBody,R2MRBody))   ! Corrected distance from the body to site #2
      RBod1 = DSQRT (DOTPR(R1MRBody,R1MRBody))   ! Corrected distance from the body to site #1



      COS_TETA = -DOTPR(R2MRBody,S)/RBod2         ! Angle between the body (wrt of the second station) and source (teta)
      COS_TETA_1 = -DOTPR(R1MRBody,S)/RBod1       ! Angle between the body (wrt of the first station) and source (teta)


      DO I = 1, 3


         RN2(I) = R2MRBODY(I)/RBOD2  !   Unit vector R2/|R2|

      END DO 

      SN = DOTPR(RN2,S)

      DO I=1,3

         SSN(I) = RN2(I) - S(I)*SN

      END DO

      SUM = DOTPR(DV2V1,SSN)

c      if (jmap.eq.3) print *, SUM, SUM/(1.d0 - cos_teta)


      COS_PHI  = DOTPR(B,S)/BB                 ! Angle between the baseline vector and source (phi)

      COS_PHI_V  = DOTPR(DV2V1,S)/VV2                 ! Angle between the w2-w1 vector and source (phi')



      COS_PSI  = DOTPR(R2MRBody,B)/(BB*RBod2)        ! Angle between the baseline vector and body (psi) as from station #2


      COS_PSI_V  = DOTPR(R2MRBody,DV2V1)/(VV2*RBod2)        ! Angle between the w2-w1 vector and body (psi) as from station #2

      TETA = DACOS(COS_TETA)
       
      PHI = DACOS(COS_PHI)

      PHI_V = DACOS(COS_PHI_V)

      R2_Dist = DSQRT(DOTPR(R2,R2))         ! geocentric distance to station #2

           
c
c     SINUS of the Angles PHI and TETA
c

      SIN_TETA = DSIN(TETA)
      SIN_PHI = DSIN(PHI)     
      SIN_PHI_V = DSIN(PHI_V)     


      SIN_PSI = DSQRT (1.d0 - COS_PSI**2)

      COS_A = - (COS_PSI + COS_PHI*COS_TETA)/(SIN_TETA * SIN_PHI)

c
c    Angle A for w2-w1 vector
c
      COS_A_V = -(COS_PSI_V + COS_PHI_V*COS_TETA)/(SIN_TETA * SIN_PHI_V)

     
c
c   First term with barycentric velocity vector
c


c      print *, VSun


cc      Freq_1 = -2.d0*GM_Body * cinv3 * bb * VSun*sin_phi*COS_A/
cc     *  ( RBody**2 *(1.d0 - COS_TETA))

      Freq_1 = -2.d0*GM_Body * cinv3 * bb * sin_phi*COS_A/
     *  ( RBody *(1.d0 - COS_TETA))
 


c
c   Second term with vector (w2-w1)
c

      
      Freq_2 = 2.d0*GM_Body * cinv3*VV2 * sin_phi_V * sin_teta*COS_A_V/
     *  ( RBody *(1.d0 - COS_TETA))
    
     
               	
      Freq_Grav = Freq_1 * VSun * sin_teta/RBody + Freq_2 !  radians


      RETURN

      END
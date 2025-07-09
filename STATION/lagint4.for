C     Last change:  JB   28 Nov 2001    7:41 am
      SUBROUTINE LAGINT4 (X,Y,N,XINT,YOUT)

C************************************************************************
C
C     SUBROUTINE LAGINT (X,Y,N,XINT,YOUT)
C
C*************************************************************
C
C ------- RECOMMENDED BY IERS !!!
C
C     THIS SUBROUTINE PERFORMS LAGRANGIAN INTERPOLATION
C     WITHIN A SET OF (X,Y) PAIRS TO GIVE THE Y
C     VALUE CORRESPONDING TO XINT.  THIS PROGRAM USES A
C     WINDOW OF 4 DATA POINTS TO PERFORM THE INTERPOLATION.
C     IF THE WINDOW SIZE NEEDS TO BE CHANGED, THIS CAN BE
C     DONE BY CHANGING THE INDICES IN THE DO LOOPS FOR
C     VARIABLES M AND J.
C
C     PARAMETERS ARE :
C     X     - ARRAY OF VALUES OF THE INDEPENDENT VARIABLE
C     Y     - ARRAY OF FUNCTION VALUES CORRESPONDING TO X
C     N     - NUMBER OF POINTS
C     XINT  - THE X-VALUE FOR WHICH ESTIMATE OF Y IS DESIRED
C     YOUT  - THE Y VALUE RETURNED TO CALLER
C
C     Distributed by Daniel Gambis (1997)
C
C
C*************************************************************
      REAL*8 X(N),Y(N),XINT,YOUT,TERM
      INTEGER N,I,J
C
      YOUT = 0.0D0
      DO 5 I = 1,N-1
        IF ( XINT .GE. X(I) .AND. XINT .LT. X(I+1) ) K = I
    5 CONTINUE
      IF ( K .LT. 2 ) K = 2
      IF ( K .GT. N-2 ) K = N-2
      DO 20 M = K-1,K+2
        TERM = Y(M)
        DO 10 J = K-1,K+2
          IF ( M .NE. J ) THEN
            TERM = TERM * (XINT - X(J))/(X(M) - X(J))
          END IF
   10   CONTINUE
        YOUT = YOUT + TERM
   20 CONTINUE
      RETURN
      END


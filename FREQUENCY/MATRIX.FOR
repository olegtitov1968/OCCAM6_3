C ********************************************************************
C SUBROUTINE ABAT FOR MATRIX PRODUCT  (AT*B*A) OR (A*B*AT)
C ********************************************************************

C USAGE        CALL ABAT (A,NRA,NCA,B,MB,MPROD,X,NRX,IER)
C
C ARGUMENTS    A  - GENERAL MATRIX OF DIMENSIONS NRA BY NCA (INPUT)
C              NRA  - NUMBER OF ROWS IN A (INPUT)
C              NCA  - NUMBER OF COLUMNS IN A (INPUT)
C
C              B  - MATRIX OF DIMENSIONS NCA BY NCA IF MPROD=0 OR
C                                        NRA BY NRA IF MPROD=1
C              MB - TYPE OF MATRIX B (INPUT)
C                   1 IF B SYMMETRIC
C                   2 IF B DIAGONAL
C
C              MPROD - IF EQUAL TO 0, COMPUTE A * B * AT
C                    - IF EQUAL TO 1, COMPUTE AT * B * A
C
C              X  - SYMMETRIC MATRIX OF DIMENSIONS NRX BY NRX (OUTPUT)
C              NRX  - NUMBER OF ROWS IN X (OUTPUT)
C
C              IER- ERROR PARAMETER (OUTPUT)
C                    IER=300 INDICATES THAT INPUT MATRIX HAS WRONGS DI-
C                            MENSIONS.
C                   OTHER SOURCES FOR AN IER NON EQUAL TO ZERO ARE:
C                    IERS IN LOC
C
C OTHERS ROUTINES REQUIRED - LOC
C
C LATEST REVISION     - JANUARY 2, 1991
C

      SUBROUTINE ABAT (A,NRA,NCA,B,MB,MPROD,X,NRX,IER)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(1),B(1),X(1)
      IER=0

C     COMPUTE NRX AND NCX
      IF (MPROD.EQ.0) NRX=NRA
      IF (MPROD.EQ.1) NRX=NCA

C     PUT ZEROS INTO THE X
      DO 1 I=1,NRX
      DO 2 J=I,NRX
         IJ=(J*(J-1))/2 + I
         X(IJ)=0.D0
2     CONTINUE
1     CONTINUE

C      ACCUMULATE RESULTS

      DO 10 I=1,NRX
         DO 20 J=I,NRX
            IJ=(J*(J-1))/2 + I
            IF (MPROD.EQ.0) THEN
               IF (MB.EQ.1) THEN
                  DO K=1,NCA
                     IK = (K-1)*NRA+I
                     IF (A(IK).NE.0.D0) THEN
                        DO L=1,NCA
                           JL = (L-1)*NRA + J
                           IF (A(JL).NE.0.D0) THEN
                              IF (L.GE.K) KL = (L*(L-1))/2 + K
                              IF (L.LT.K) KL = (K*(K-1))/2 + L
                              X(IJ) = X(IJ)+A(IK)*B(KL)*A(JL)
                           ENDIF
                        ENDDO
                     ENDIF
                  ENDDO

               ELSE IF (MB.EQ.2) THEN
                  DO K=1,NCA
                     IK = (K-1)*NRA+I
                     IF (A(IK).NE.0.D0) THEN
                        JK = (K-1)*NRA + J
                        IF (A(JK).NE.0.D0) X(IJ)=X(IJ)+A(IK)*B(K)*A(JK)
                     ENDIF
                  ENDDO

               ENDIF

            ELSE IF (MPROD.EQ.1) THEN
               IF (MB.EQ.1) THEN
                  DO K=1,NRA
                     KI = (I-1)*NRA+K
                     IF (A(KI).NE.0.D0) THEN
                        DO L=1,NRA
                           LJ = (J-1)*NRA + L
                           IF (A(LJ).NE.0.D0) THEN
                              IF (L.GE.K) KL = (L*(L-1))/2 + K
                              IF (L.LT.K) KL = (K*(K-1))/2 + L
                              X(IJ) = X(IJ)+A(KI)*B(KL)*A(LJ)
                           ENDIF
                        ENDDO
                     ENDIF
                  ENDDO

               ELSE IF (MB.EQ.2) THEN
                  DO K=1,NRA
                     KI = (I-1)*NRA+K
                     IF (A(KI).NE.0.D0) THEN
                        KJ = (J-1)*NRA + K
                        IF (A(KJ).NE.0.D0) X(IJ)=X(IJ)+A(KI)*B(K)*A(KJ)
                     ENDIF
                  ENDDO

               ENDIF

            ENDIF

20       CONTINUE
10    CONTINUE

      END

C ********************************************************************

C USAGE        CALL ATBC (A,NRA,NCA,B,MB,C,NCC,X,NRX,NCX,IER)
C                                                          T
C PURPOSE      FAST COMPUTATION OF THREE MATRICES PRODUCT A B C
C              OUTPUT IN GENERAL MODE
C
C ARGUMENTS    A  - GENERAL MATRIX OF DIMENSIONS NRA BY NCA (INPUT)
C              NRA  - NUMBER OF ROWS IN A (INPUT)
C              NCA  - NUMBER OF COLUMNS IN A (INPUT)
C
C              B  - MATRIX OF DIMENSIONS NRA BY NRA (INPUT)
C              MB - 1 IF SYMMETRIC AND 2 IF DIAGONAL
C
C              C  - GENERAL MATRIX OF DIMENSIONS NRA BY NCC (INPUT)
C              NCC  - NUMBER OF COLUMNS IN C (INPUT)
C
C              X  - GENERAL MATRIX OF DIMENSIONS NRX BY NCX (OUTPUT)
C              NRX  - NUMBER OF ROWS IN X (OUTPUT)
C              NCX  - NUMBER OF COLUMNS IN X (OUTPUT)
C
C              IER- ERROR PARAMETER (OUTPUT)
C                    IER=300 INDICATES THAT INPUT MATRIX HAS WRONGS DI-
C                            MENSIONS.
C                    IER=310 MEANS THAT MODA, MODB,MODC OR MODX ARE
C                            .GT.3 .OR. .LT.0
C
C                   OTHER SOURCES FOR AN IER NON EQUAL TO ZERO ARE:
C                    IERS IN LOC
C
C OTHERS ROUTINES REQUIRED - LOC
C
C LATEST REVISION     - JANUARY 3, 1991
C

      SUBROUTINE ATBC (A,NRA,NCA,B,MB,C,NCC,X,NRX,NCX,IER)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(1),B(1),C(1),X(1)
      IER=0
C     COMPUTE NRX AND NCX
      NRX=NCA
      NCX=NCC

C     PUT ZEROS INTO THE X
      DO 1 I=1,NRX
         DO 2 J=1,NCX
            ij=(j-1)*nrx+i
            X(IJ)=0.D0
2        CONTINUE
1     CONTINUE

C      ACCUMULATE RESULTS

      DO 10 I=1,NRX
         DO 20 J=1,NCX
            ij=(j-1)*nrx+i
            DO 30 K=1,nra
               ki=(i-1)*nra+k
               IF (A(KI).NE.0.D0) THEN
                  IF(MB.EQ.1) THEN
                     DO L=1,NRA
                        lj=(j-1)*nra+l
                        IF (C(LJ).NE.0.D0) THEN
                           if (k.ge.l) kl=(k*(k-1))/2+l
                           if (k.lt.l) kl=(l*(l-1))/2+k
                           X(IJ)=X(IJ)+A(ki)*B(kl)*c(lj)
                        ENDIF
                     ENDDO
                  ELSE IF (MB.EQ.2) THEN
                     KJ=(J-1)*NRA+K
                     IF (C(KJ).NE.0.D0) X(IJ) = X(IJ)+A(KI)*B(K)*C(KJ)
                  ENDIF
               ENDIF
30          CONTINUE
20       CONTINUE
10    CONTINUE

      RETURN

      END

C ********************************************************************
C SUBROUTINE DSDT FOR MATRIX PRODUCT  (D*S*DT )
C ********************************************************************

C USAGE        CALL DSDT (A,NRA,B,MB,X,NRX,IER)
C
C ARGUMENTS    A  - DIAGONAL MATRIX OF DIMENSIONS NRA (INPUT)
C              NRA  - NUMBER OF ROWS AND COLUMNS IN A (INPUT)
C
C              B  - SYMMETRIC OR DIAGONAL MATRIX OF DIMENSIONS NRA BY NRA
C              MB - TYPE OF MATRIX B (INPUT)
C                   1 IF B SYMMETRIC
C                   2 IF B DIAGONAL
C
C              X  - SYMMETRIC MATRIX OF DIMENSIONS NRX BY NRX (OUTPUT)
C              NRX  - NUMBER OF ROWS IN X (OUTPUT)
C
C              IER- ERROR PARAMETER (OUTPUT)
C                    IER=300 INDICATES THAT INPUT MATRIX HAS WRONGS DI-
C                            MENSIONS.
C                   OTHER SOURCES FOR AN IER NON EQUAL TO ZERO ARE:
C                    IERS IN LOC
C
C OTHERS ROUTINES REQUIRED - LOC
C
C CREATED - JULY 13, 1993
C

      SUBROUTINE DSDT (A,NRA,B,MB,X,NRX,IER)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(1),B(1),X(1)
      IER=0

C     COMPUTE NRX
      NRX=NRA

C     PUT ZEROS INTO THE X
      DO 1 I=1,NRX
         DO 2 J=I,NRX
            IJ=(J*(J-1))/2 + I
            X(IJ)=0.D0
2        CONTINUE
1     CONTINUE

C      ACCUMULATE RESULTS

      IF (MB.EQ.1) THEN

	 DO I=1,NRX
            DO J=I,NRX
               IJ=(J*(J-1))/2 + I

               X(IJ) = B(IJ) * A(I) * A(J)
            ENDDO
         ENDDO

      ELSE IF (MB.EQ.2) THEN

         DO I=1,NRX
	    X(I) = B(I) * A(I) * A(I)
         ENDDO
      ENDIF

      END

C ********************************************************************
C SUBROUTINE FIPR FOR MATRIX PRODUCT  (FI*P*FIT)
C ********************************************************************

C USAGE        CALL FIPR (A,NRA,B,MB,X,NRX,IER)
C
C ARGUMENTS    A  - TRANSITION MATRIX OF DIMENSIONS NRA BY NCA (INPUT)
C              NRA  - NUMBER OF ROWS AND COLUMNS IN A (INPUT)
C
C              B  - MATRIX OF DIMENSIONS NRA BY NRA
C              MB - TYPE OF MATRIX B (INPUT)
C                   1 IF B SYMMETRIC
C                   2 IF B DIAGONAL
C
C              X  - SYMMETRIC MATRIX OF DIMENSIONS NRX BY NRX (OUTPUT)
C              NRX  - NUMBER OF ROWS IN X (OUTPUT)
C
C              IER- ERROR PARAMETER (OUTPUT)
C                    IER=300 INDICATES THAT INPUT MATRIX HAS WRONGS DI-
C                            MENSIONS.
C                   OTHER SOURCES FOR AN IER NON EQUAL TO ZERO ARE:
C                    IERS IN LOC
C
C OTHERS ROUTINES REQUIRED - LOC
C
C LATEST REVISION     - JANUARY 2, 1991
C

      SUBROUTINE FIPR (A,NRA,B,MB,X,NRX,IER)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(1),B(1),X(1)
      IER=0

C     COMPUTE NRX
      NRX=NRA

C     PUT ZEROS INTO THE X
      DO 1 I=1,NRX
         DO 2 J=I,NRX
            IJ=(J*(J-1))/2 + I
            X(IJ)=0.D0
2        CONTINUE
1     CONTINUE

C      ACCUMULATE RESULTS

      DO 10 I=1,NRX
         DO 20 J=I,NRX
            IJ=(J*(J-1))/2 + I

            IF (MB.EQ.1) THEN
               DO K=I,MIN(I+1,NRX)
                  IK = (K-1)*NRA+I
                  IF (A(IK).NE.0.D0) THEN
                     DO L=J,MIN(J+1,NRA)
                        JL = (L-1)*NRA + J
                        IF (A(JL).NE.0.D0) THEN
                           IF (L.GE.K) KL = (L*(L-1))/2 + K
                           IF (L.LT.K) KL = (K*(K-1))/2 + L
                           X(IJ) = X(IJ)+A(IK)*B(KL)*A(JL)
                        ENDIF
                     ENDDO
                  ENDIF
               ENDDO

            ELSE IF (MB.EQ.2) THEN
               DO K=I,MIN(I+1,NRX)
                  IF (J.EQ.K.OR.J+1.EQ.K) THEN
                     IK = (K-1)*NRA+I
                     IF (A(IK).NE.0.D0) THEN
                        JK = (K-1)*NRA + J
                        IF (A(JK).NE.0.D0) X(IJ)=X(IJ)+A(IK)*B(K)*A(JK)
                     ENDIF
                  ENDIF
               ENDDO

            ENDIF

20       CONTINUE
10    CONTINUE

      END


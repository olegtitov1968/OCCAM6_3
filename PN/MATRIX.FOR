C  CREATED  JULY 13, 1993   General matrix operation routines
C    NZ               Previously contained in leastt.for and kalman.for

C *********************************************************************
C   SUBROUTINE ABX FOR GENERAL TWO-MATRIX MULTIPLICATION (A*B)
C *********************************************************************

C USAGE        CALL ABX (A,NRA,NCA,MODA,TA,B,NRB,NCB,MODB,TB,X,NRX,NCX,
C             +MODX,IER)
C
C ARGUMENTS    A  - MATRIX OF DIMENSIONS NRA BY NCA (INPUT)
C              NRA  - NUMBER OF ROWS IN A (INPUT)
C              NCA  - NUMBER OF COLUMNS IN A (INPUT)
C              MODA - MODE OF MATRIX A (INPUT)
C                     IF MODA = 0 THEN THE MATRIX A IS A GENERAL MATRIX
C                     IF MODA = 1 THEN THE MATRIX A IS SYMMETRIC
C                     IF MODA = 2 THEN THE MATRIX A IS DIAGONAL
C                     IF MODA = 3 THEN A IS AN ESCALAR MATRIX
C              TA - TYPE OF MATRIX A (INPUT)
C                   IF TA = 'N' THE SUBROUTINE USES THE MATRIX A
C                   IF TA = 'T' THE SUBROUTINE USES THE TRANSPONSE OF A
C
C              B  - MATRIX OF DIMENSIONS NRB BY NCB (INPUT)
C              NRB  - NUMBER OF ROWS IN B (INPUT)
C              NCB  - NUMBER OF COLUMNS IN B (INPUT)
C              MODB - MODE OF MATRIX B (SAME CODE AS MODA) (INPUT)
C              TB - TYPE OF MATRIX B (INPUT)
C                    IF TB = 'N' THE SUBROUTINE USES THE MATRIX B
C                    IF TB = 'T' THE SUBROUTINE USES THE TRANSPONSE OF B
C
C              X  - MATRIX OF DIMENSIONS NRX BY NCX (OUTPUT)
C              NRX  - NUMBER OF ROWS IN X (OUTPUT)
C              NCX  - NUMBER OF COLUMNS IN X (OUTPUT)
C              MODX - MODE OF MATRIX X (SAME CODE AS MODA) (INPUT)
C
C              IER- ERROR PARAMETER (OUTPUT)
C                    IER=200 INDICATES THAT THE DIMENSIONS OF INPUT
C                            MATRIX ARE WRONG.
C                    IER=210 MEANS THAT MODA,MODB OR MODX.GT.3.OR..LT.0
C                   OTHER SOURCES FOR AN IER NON EQUAL TO ZERO ARE:
C                    IERS IN LOCB
C
C OTHERS ROUTINES REQUIRED - LOCB
C
C LATEST REVISION     - JUNE 27, 1989
C

      SUBROUTINE ABX (A,NRA,NCA,MODA,TA,B,NRB,NCB,MODB,TB,X,NRX,NCX,
     +MODX,IER)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER TA*1,TB*1
      DIMENSION A(1),B(1),X(1)
      IER=0
      IF (MODA.GT.3.OR.MODA.LT.0) IER=210
      IF (MODB.GT.3.OR.MODB.LT.0) IER=210
      IF (MODX.GT.3.OR.MODX.LT.0) IER=210
      IF (IER.GE.100) RETURN

      IF (TA.EQ.'N'.AND.TB.EQ.'N') THEN
           IF(NCA.NE.NRB) THEN
           IER = 200
           RETURN
           ENDIF
      ELSE IF (TA.EQ.'N'.AND.TB.EQ.'T') THEN
           IF(NCA.NE.NCB) THEN
           IER = 200
           RETURN
           ENDIF
      ELSE IF (TA.EQ.'T'.AND.TB.EQ.'N') THEN
           IF(NRA.NE.NRB) THEN
           IER = 200
           RETURN
           ENDIF
      ELSE IF (TA.EQ.'T'.AND.TB.EQ.'T') THEN
           IF(NRA.NE.NCB) THEN
           IER = 200
           RETURN
           ENDIF
      ENDIF
C     COMPUTE NRX AND NCX

      NRX=NRA
      IF (TA.EQ.'T') NRX=NCA
      NCX=NCB
      IF (TB.EQ.'T') NCX=NRB


C     PUT ZEROS INTO THE X
      DO 1 I=1,NRX
      DO 2 J=1,NCX
      IF ((MODX.EQ.1).AND.(I.LT.J)) GO TO 2
      IF ((MODX.EQ.2).AND.(I.NE.J)) GO TO 2

           call locb (I,J,IJ,NRX,NCX,MODX,IER)
           X(IJ)=0.D0
2     CONTINUE
1     CONTINUE
C      ACCUMALTE RESULTS
      DO 10 I=1,NRX
      DO 20 J=1,NCX
      IF ((MODX.EQ.1).AND.(I.LT.J)) GO TO 20
      IF ((MODX.EQ.2).AND.(I.NE.J)) GO TO 20
      call locb (I,J,IJ,NRX,NCX,MODX,IER)
      NDUMM1=NCA
      IF (TA.EQ.'T') NDUMM1=NRA
      DO 30 K=1,NDUMM1
      IF (TA.EQ.'N') call locb (I,K,IPOSA,NRA,NCA,MODA,IER)
      IF (TA.EQ.'T') call locb (K,I,IPOSA,NRA,NCA,MODA,IER)
      IF (TB.EQ.'N') call locb (K,J,IPOSB,NRB,NCB,MODB,IER)
      IF (TB.EQ.'T') call locb (J,K,IPOSB,NRB,NCB,MODB,IER)
      IF((IJ.NE.0).AND.(IPOSA.NE.0).AND.(IPOSB.NE.0))THEN
            X(IJ)=X(IJ)+A(IPOSA)*B(IPOSB)

      ENDIF
30    CONTINUE
20    CONTINUE
10    CONTINUE
      RETURN

      END


C ********************************************************************
C SUBRUTINE ABCX FOR GENERAL THREE-MATRIX PRODUCT  (A*B*C)
C ********************************************************************

C USAGE        CALL ABCX (A,NRA,NCA,MODA,TA,B,NRB,NCB,MODB,TB,
C                         C,NRC,NCC,MODC,TC,X,NRX,NCX,MODX,IER)
C
C ARGUMENTS    A  - MATRIX OF DIMENSIONS NRA BY NCA (INPUT)
C              NRA  - NUMBER OF ROWS IN A (INPUT)
C              NCA  - NUMBER OF COLUMNS IN A (INPUT)
C              MODA - MODE OF MATRIX A (INPUT)
C                     IF MODA = 0 THEN THE MATRIX A IS A GENERAL MATRIX
C                     IF MODA = 1 THEN THE MATRIX A IS SYMMETRIC
C                     IF MODA = 2 THEN THE MATRIX A IS DIAGONAL
C                     IF MODA = 3 THEN A IS AN ESCALAR MATRIX
C              TA - TYPE OF MATRIX A (INPUT)
C                   IF TA = 'N' THE SUBROUTINE USES THE MATRIX A
C                   IF TA = 'T' THE SUBROUTINE USES THE TRANSPONSE OF A
C
C              B  - MATRIX OF DIMENSIONS NRB BY NCB (INPUT)
C              NRB  - NUMBER OF ROWS IN B (INPUT)
C              NCB  - NUMBER OF COLUMNS IN B (INPUT)
C              MODB - MODE OF MATRIX B (SAME CODE AS MODA) (INPUT)
C              TB - TYPE OF MATRIX B (INPUT)
C                   IF TB = 'N' THE SUBROUTINE USES THE MATRIX B
C                   IF TB = 'T' THE SUBROUTINE USES THE TRANSPONSE OF B
C
C              C  - MATRIX OF DIMENSIONS NRC BY NCC (INPUT)
C              NRC  - NUMBER OF ROWS IN C (INPUT)
C              NCC  - NUMBER OF COLUMNS IN C (INPUT)
C              MODC - MODE OF MATRIX C (SAME CODE AS MODA) (INPUT)
C              TC - TYPE OF MATRIX C (INPUT)
C                   IF TC = 'N' THE SUBROUTINE USES THE MATRIX C
C                   IF TC = 'T' THE SUBROUTINE USES THE TRANSPONSE OF C
C
C              X  - MATRIX OF DIMENSIONS NRX BY NCX (OUTPUT)
C              NRX  - NUMBER OF ROWS IN X (OUTPUT)
C              NCX  - NUMBER OF COLUMNS IN X (OUTPUT)
C              MODX - MODE OF MATRIX X (SAME CODE AS MODA) (INPUT)
C
C              IER- ERROR PARAMETER (OUTPUT)
C                    IER=300 INDICATES THAT INPUT MATRIX HAS WRONGS DI-
C                            MENSIONS.
C                    IER=310 MEANS THAT MODA, MODB,MODC OR MODX ARE
C                            .GT.3 .OR. .LT.0
C                   OTHER SOURCES FOR AN IER NON EQUAL TO ZERO ARE:
C                    IERS IN LOCB
C
C OTHERS ROUTINES REQUIRED - LOCB
C
C LATEST REVISION     - JUNE 27, 1989
C


      SUBROUTINE ABCX (A,NRA,NCA,MODA,TA,B,NRB,NCB,MODB,TB,
     +C,NRC,NCC,MODC,TC,X,NRX,NCX,MODX,IER)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER TA*1,TB*1,TC*1
      DIMENSION A(1),B(1),C(1),X(1)
      IER=0
      IF (MODA.GT.3.OR.MODA.LT.0) IER=310
      IF (MODB.GT.3.OR.MODB.LT.0) IER=310
      IF (MODC.GT.3.OR.MODC.LT.0) IER=310
      IF (MODX.GT.3.OR.MODX.LT.0) IER=310

      IF (TA.EQ.'N'.AND.TB.EQ.'N') THEN
           IF(NCA.NE.NRB) THEN
           IER = 300
           RETURN
           ENDIF
      ELSE IF (TA.EQ.'N'.AND.TB.EQ.'T') THEN
           IF(NCA.NE.NCB) THEN
           IER = 300
           RETURN
           ENDIF
      ELSE IF (TA.EQ.'T'.AND.TB.EQ.'N') THEN
           IF(NRA.NE.NRB) THEN
           IER = 300
           RETURN
           ENDIF
      ELSE IF (TA.EQ.'T'.AND.TB.EQ.'T') THEN
           IF(NRA.NE.NCB) THEN
           IER = 300
           RETURN
           ENDIF
      ENDIF

      IF (TB.EQ.'N'.AND.TC.EQ.'N') THEN
           IF(NCB.NE.NRC) THEN
           IER = 300
           RETURN
           ENDIF
      ELSE IF (TB.EQ.'N'.AND.TC.EQ.'T') THEN
           IF(NCB.NE.NCC) THEN
           IER = 300
           RETURN
           ENDIF
      ELSE IF (TB.EQ.'T'.AND.TC.EQ.'N') THEN
           IF(NRB.NE.NRC) THEN
           IER = 300
           RETURN
           ENDIF
      ELSE IF (TB.EQ.'T'.AND.TC.EQ.'T') THEN
           IF(NRB.NE.NCC) THEN
           IER = 300
           RETURN
           ENDIF
      ENDIF
C     COMPUTE NRX AND NCX
      NRX=NRA
      IF (TA.EQ.'T') NRX=NCA
      NCX=NCC
      IF (TC.EQ.'T') NCX=NRC

C     PUT ZEROS INTO THE X
      DO 1 I=1,NRX
      DO 2 J=1,NCX
      IF ((MODX.EQ.1).AND.(I.LT.J)) GO TO 2
      IF ((MODX.EQ.2).AND.(I.NE.J)) GO TO 2
           call locb (I,J,IJ,NRX,NCX,MODX,IER)
           X(IJ)=0.D0
2     CONTINUE
1     CONTINUE
C      ACCUMULATE RESULTS
      DO 10 I=1,NRX
      DO 20 J=1,NCX
      IF ((MODX.EQ.1).AND.(I.LT.J)) GO TO 20
      IF ((MODX.EQ.2).AND.(I.NE.J)) GO TO 20
      call locb (I,J,IJ,NRX,NCX,MODX,IER)
      NDUMM1=NCA
      NDUMM2=NCB
      IF (TA.EQ.'T') NDUMM1=NRA
      IF (TB.EQ.'T') NDUMM2=NRB
      DO 30 K=1,NDUMM1
      IF (TA.EQ.'N') call locb (I,K,IPOSA,NRA,NCA,MODA,IER)
      IF (TA.EQ.'T') call locb (K,I,IPOSA,NRA,NCA,MODA,IER)
      DO  40 L=1,NDUMM2
      IF (TB.EQ.'N') call locb (K,L,IPOSB,NRB,NCB,MODB,IER)
      IF (TB.EQ.'T') call locb (L,K,IPOSB,NRB,NCB,MODB,IER)
      IF (TC.EQ.'N') call locb (L,J,IPOSC,NRC,NCC,MODC,IER)
      IF (TC.EQ.'T') call locb (J,L,IPOSC,NRC,NCC,MODC,IER)
      IF((IJ.NE.0).AND.(IPOSA.NE.0).AND.(IPOSB.NE.0).AND.
     +(IPOSC.NE.0))THEN
            X(IJ)=X(IJ)+A(IPOSA)*B(IPOSB)*C(IPOSC)
      ENDIF
40    CONTINUE
30    CONTINUE
20    CONTINUE
10    CONTINUE

      RETURN

      END

C ********************************************************************

C USAGE        CALL ABCXp (A,NRA,NCA,B,NRB,NCB,C,NRC,NCC,
C                            X,NRX,NCX,igen,IER)
C                                                          T
C PURPOSE      FAST COMPUTATION OF THREE MATRICES PRODUCT A B C
C              OUTPUT IN GENERAL MODE
C
C ARGUMENTS    A  - GENERAL MATRIX OF DIMENSIONS NRA BY NCA (INPUT)
C              NRA  - NUMBER OF ROWS IN A (INPUT)
C              NCA  - NUMBER OF COLUMNS IN A (INPUT)
C
C              B  - MATRIX OF DIMENSIONS NRB BY NCB (INPUT)
C              NRB  - NUMBER OF ROWS IN B (INPUT)
C              NCB  - NUMBER OF COLUMNS IN B (INPUT)
C
C              C  - GENERAL MATRIX OF DIMENSIONS NRB BY NCB (INPUT)
C              NRC  - NUMBER OF ROWS IN B (INPUT)
C              NCC  - NUMBER OF COLUMNS IN B (INPUT)
C
C              X  - MATRIX OF DIMENSIONS NRX BY NCX (OUTPUT)
C              NRX  - NUMBER OF ROWS IN X (OUTPUT)
C              NCX  - NUMBER OF COLUMNS IN X (OUTPUT)
C              IGEN - INPUT PARAMETER. SET THE CHARACTER OF THE OUTPUT
C                     MATRIX TO INCREASE COMPUTATION SPEED
C                     1 IF X IS SYMMETRIC
C                     0 IF X IS GENERAL
C
C              IER- ERROR PARAMETER (OUTPUT)
C                    IER=300 INDICATES THAT INPUT MATRIX HAS WRONGS DI-
C                            MENSIONS.
C                    IER=310 MEANS THAT MODA, MODB,MODC OR MODX ARE
C                            .GT.3 .OR. .LT.0
C
C                   OTHER SOURCES FOR AN IER NON EQUAL TO ZERO ARE:
C                    IERS IN LOCB
C
C OTHERS ROUTINES REQUIRED - LOCB
C
C LATEST REVISION     - JUNE 27, 1989
C


      SUBROUTINE ABCXp (A,NRA,NCA,B,NRB,NCB,C,NRC,NCC,
     *                   X,NRX,NCX,igen,IER)

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
2     CONTINUE
1     CONTINUE
C      ACCUMULATE RESULTS
      DO 10 I=1,NRX
      i1=i
      if (igen.eq.0) i1=1
      DO 20 J=i1,NCX
      ij=(j-1)*nrx+i
      DO 30 K=1,nra
      ki=(i-1)*nra+k
      DO  40 L=1,nrc
         if (k.ge.l) kl=(k*(k-1))/2+l
         if (k.lt.l) kl=(l*(l-1))/2+k
         lj=(j-1)*nra+l
            X(IJ)=X(IJ)+A(ki)*B(kl)*c(lj)
40    CONTINUE
30    CONTINUE
      if (igen.eq.1) THEN
        ji=(i-1)*nrx+j
        x(ji)=x(ij)
      ENDIF
20    CONTINUE
10    CONTINUE

      RETURN

      END


C ********************************************************************
C SUBROUTINE LOCB
C ********************************************************************

C USAGE        call locb (I,J,IJ,N,M,MS,IER)
C
C PURPOSE      COMPUTE THE VECTOR SUBSCRIPT POSITION FOR AN ELEMENT
C              IN A MATRIX OF SPECIFIED STORAGE MODE
C
C ARGUMENTS    I - ROW SUBSCRIPT OF THE MATRIX ELEMENT (INPUT)
C              J - COLUMN SUBSCRIPT OF THE MATRIX ELEMENT (INPUT)
C              IJ - SUBSCRIPT OF THE VECTOR (OUTPUT)
C              N  - NUMBER OF ROWS IN MATRIX (INPUT)
C              M  - NUMBER OF COLUMNS IN MATRIX (INPUT)
C              MS - MODE OF THE MATRIX  (INPUT)
C                     MS=0   GENERAL
C                     MS=1   SYMMETRIC
C                     MS=2   DIAGONAL
C                     MS=3   ESCALAR
C              IER- ERROR PARAMETER (OUTPUT)
C                    IER=104 INDICATES THAT I*J.GT.N*M
C                    IER=105 INDICATES THAT N.NE.M (WHEN THE MATRIX IS
C                            SYMMETRIC OR DIAGONAL)
C                    IER=106 INDICATES THAT N.NE.M (WHEN THE MATRIX IS
C                            ESCALAR)
C REQUIRED SUBROUTINES - NONE REQUIRED.
C
C LATEST REVISION     - JUNE 27, 1989
C

      SUBROUTINE LOCb (I,J,IJ,N,M,MS,IER)

       IER=0

       IF(I*J.GT.N*M) THEN
       IER=104
       RETURN
       ENDIF


      IF(MS.EQ.0) THEN
        IJ=N*(J-1)+I
      ENDIF

      IF(MS.EQ.1) THEN
        IF(I.LT.J) THEN
          IJ=I+(J*J-J)/2

        IF(N.NE.M) THEN
        IER=105
       RETURN
       ENDIF
        ENDIF
        IF(I.GE.J) THEN
          IJ=J+(I*I-I)/2
        ENDIF
       ENDIF

      IF(MS.EQ.2) THEN
        IF(N.NE.M) THEN
        IER=105
       RETURN
       ENDIF
        IF(I.EQ.J) THEN
          IJ=I
        ENDIF
        IF(I.NE.J) THEN
          IJ=0
        ENDIF
      ENDIF
      IF(MS.EQ.3) THEN
        IF(N.NE.M) THEN
        IER=106
       RETURN
       ENDIF
        IF(I.EQ.J) THEN
          IJ=1
        ENDIF
        IF(I.NE.J) THEN
          IJ=0
        ENDIF
      ENDIF

      RETURN
      END

C *******************************************************************
C     SUBROUTINE ADD
C *******************************************************************

      subroutine add (a,nra,nca,moda,plus,b,nrb,ncb,modb,x,nrx,ncx,
     +modx,ier)

      implicit double precision (a-h,o-z)
      character plus*1
      dimension A(1),B(1),X(1)
      IER=0
      if (nra.ne.nrb) ier=401
      if (nca.ne.ncb) ier=402
      if (ier.ne.0) return

      if (nra.ne.nrx) ier=403
      if (nca.ne.ncx) ier=404
      if (ier.ne.0) return

C      ACCUMALTE RESULTS
      DO 10 I=1,NRX
      DO 20 J=1,NCX
      IF ((MODX.EQ.1).AND.(I.LT.J)) GO TO 20
      IF ((MODX.GT.1).AND.(I.NE.J)) GO TO 20
         call locb (I,J,IJ,NRX,NCX,MODX,IER)
         IF (IJ.EQ.0) GO TO 20
         call locb (I,J,IPOSA,NRA,NCA,MODA,IER)
         call locb (I,J,IPOSB,NRB,NCB,MODB,IER)
            TEMP=0.D0
            if (iposa.ne.0) TEMP=a(iposa)
            if (iposb.ne.0) then
                 if (plus.eq.'+') TEMP=TEMP+b(iposb)
                 if (plus.eq.'-') TEMP=TEMP-b(iposb)
            endif
            X(IJ)=TEMP
20    CONTINUE
10    CONTINUE
      RETURN
      END

C ********************************************************************
C SUBRUTINE ABCXG FOR THREE GENERAL MATRIX PRODUCT  (A*B*C)
C ********************************************************************

C USAGE        CALL ABCXG (A,NRA,NCA,TA,B,NRB,NCB,TB,
C                         C,NRC,NCC,TC,X,NRX,NCX,MODX,IER)
C
C ARGUMENTS    A  - GENERAL MATRIX OF DIMENSIONS NRA BY NCA (INPUT)
C              NRA  - NUMBER OF ROWS IN A (INPUT)
C              NCA  - NUMBER OF COLUMNS IN A (INPUT)
C              TA - TYPE OF MATRIX A (INPUT)
C                   IF TA = 'N' THE SUBROUTINE USES THE MATRIX A
C                   IF TA = 'T' THE SUBROUTINE USES THE TRANSPONSE OF A
C
C              B  - GENERAL MATRIX OF DIMENSIONS NRB BY NCB (INPUT)
C              NRB  - NUMBER OF ROWS IN B (INPUT)
C              NCB  - NUMBER OF COLUMNS IN B (INPUT)
C              TB - TYPE OF MATRIX B (INPUT)
C                   IF TB = 'N' THE SUBROUTINE USES THE MATRIX B
C                   IF TB = 'T' THE SUBROUTINE USES THE TRANSPONSE OF B
C
C              C  - GENERAL MATRIX OF DIMENSIONS NRC BY NCC (INPUT)
C              NRC  - NUMBER OF ROWS IN C (INPUT)
C              NCC  - NUMBER OF COLUMNS IN C (INPUT)
C              TC - TYPE OF MATRIX C (INPUT)
C                   IF TC = 'N' THE SUBROUTINE USES THE MATRIX C
C                   IF TC = 'T' THE SUBROUTINE USES THE TRANSPONSE OF C
C
C              X  - MATRIX OF DIMENSIONS NRX BY NCX (OUTPUT)
C              NRX  - NUMBER OF ROWS IN X (OUTPUT)
C              NCX  - NUMBER OF COLUMNS IN X (OUTPUT)
C              MODX - MODE OF MATRIX X (SAME CODE AS MODA) (INPUT)
C
C              IER- ERROR PARAMETER (OUTPUT)
C                    IER=300 INDICATES THAT INPUT MATRIX HAS WRONGS DI-
C                            MENSIONS.
C                    IER=310 MEANS THAT MODA, MODB,MODC OR MODX ARE
C                            .GT.3 .OR. .LT.0
C                   OTHER SOURCES FOR AN IER NON EQUAL TO ZERO ARE:
C                    IERS IN LOCB
C
C OTHERS ROUTINES REQUIRED - LOCB
C
C LATEST REVISION     - JUNE 27, 1989
C


      SUBROUTINE ABCXG (A,NRA,NCA,TA,B,NRB,NCB,TB,
     +C,NRC,NCC,TC,X,NRX,NCX,MODX,IER)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER TA*1,TB*1,TC*1
      DIMENSION A(1),B(1),C(1),X(1)
      IER=0
      IF (MODX.GT.3.OR.MODX.LT.0) IER=310

C     COMPUTE NRX AND NCX
      NRX=NRA
      IF (TA.EQ.'T') NRX=NCA
      NCX=NCC
      IF (TC.EQ.'T') NCX=NRC

C     PUT ZEROS INTO THE X
      DO 1 I=1,NRX
         DO 2 J=1,NCX
            IF ((MODX.EQ.1).AND.(I.LT.J)) GO TO 2
            IF ((MODX.EQ.2).AND.(I.NE.J)) GO TO 2
            call locb (I,J,IJ,NRX,NCX,MODX,IER)
            X(IJ)=0.D0
2        CONTINUE
1     CONTINUE

C      ACCUMULATE RESULTS

      DO 10 I=1,NRX
         DO 20 J=1,NCX
            IF ((MODX.EQ.1).AND.(I.LT.J)) GO TO 20
            IF ((MODX.EQ.2).AND.(I.NE.J)) GO TO 20
            call locb (I,J,IJ,NRX,NCX,MODX,IER)
            NDUMM1=NCA
            NDUMM2=NCB
            IF (TA.EQ.'T') NDUMM1=NRA
            IF (TB.EQ.'T') NDUMM2=NRB
            DO 30 K=1,NDUMM1
               IF (TA.EQ.'N') IPOSA=I+(K-1)*NRA
               IF (TA.EQ.'T') IPOSA=K+(I-1)*NRA
               DO  40 L=1,NDUMM2
                  IF (TB.EQ.'N') IPOSB=K+(L-1)*NRB
                  IF (TB.EQ.'T') IPOSB=L+(K-1)*NRB
                  IF (TC.EQ.'N') IPOSC=L+(J-1)*NRC
                  IF (TC.EQ.'T') IPOSC=J+(L-1)*NRC
                  IF((IJ.NE.0).AND.(IPOSA.NE.0).AND.(IPOSB.NE.0).AND.
     +               (IPOSC.NE.0))THEN
                     X(IJ)=X(IJ)+A(IPOSA)*B(IPOSB)*C(IPOSC)
                  ENDIF
40             CONTINUE
30          CONTINUE
20       CONTINUE
10    CONTINUE

      RETURN

      END

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
C                    IERS IN LOCB
C
C OTHERS ROUTINES REQUIRED - LOCB
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
C                    IERS IN LOCB
C
C OTHERS ROUTINES REQUIRED - LOCB
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
C                    IERS IN LOCB
C
C OTHERS ROUTINES REQUIRED - LOCB
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
C                    IERS IN LOCB
C
C OTHERS ROUTINES REQUIRED - LOCB
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


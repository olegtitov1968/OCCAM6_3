      SUBROUTINE UNFOLD (LU,CNAME)
C
C     THIS ROUTINE TAKES THE PLACE OF AN OPEN STATEMENT FOR THE FILES
C   'SORTIM','STATIM','BASTIM','STACAT' AND 'DICTIO'
C
c  REVISION:  22, FEBRUARY, 2001 (jboehm)
c  REVISION:  27, July, 2001 (OT)
c  REVISION:  13, May, 2002 (OT) - STATIM changed due to increasing of
c                  PPOLAR array - 9 elememts instead of 3  (+48=6*8)
c  REVISION:  21, May, 2002 (OT) - SORTIM changed due to adding of
c                  NUT -  (+4)
c  REVISION:  08, June, 2003 (OT) - DICTIO changed due to increasing
c                          of ntimmax
c  REVISION:  24, September, 2003 (JB) - STATIM changed to 800
c  REVISION:  03, March, 2006 (OT)
c  REVISION:  19, July, 2006 (OT) In BASTIM 144 changed to 152 due to gamma
c  LAST REVISION:  30, August, 2006 (OT) In SORTIM 1184 changed to 1052
c                                        In DICTIO 12028 changed to 8028
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      CHARACTER*6 CNAME

      IF (LU.EQ.17) OPEN (17,FILE='BASTIM',ACCESS='DIRECT',RECL=296,  ! acc_der
     *FORM='UNFORMATTED')
c
c   Last change on 17.11.2005 due to more coordinates in STATION.FOR
c

      IF (LU.EQ.18) OPEN (18,FILE='STATIM',ACCESS='DIRECT',RECL=1052,
     *FORM='UNFORMATTED')

c
c   Last change on 14.11.2005 due to aberration term at the nutation partials
c
      IF (LU.EQ.19) OPEN (19,FILE='SORTIM',ACCESS='DIRECT',RECL=1286,
     *FORM='UNFORMATTED')

      IF (LU.EQ.20) OPEN (20,FILE='STACAT',ACCESS='DIRECT',RECL=128,
     *FORM='UNFORMATTED')
c
c   Change by Johannes

      IF (LU.EQ.21) OPEN (21,FILE='DICTIO',ACCESS='DIRECT',RECL=20028,
     *FORM='UNFORMATTED')

      RETURN
      END

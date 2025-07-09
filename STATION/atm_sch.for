c+ jboehm
      subroutine atm_sch(model,stan,i,tmjddec,dxatm,dyatm,dzatm)
c- jboehm

c This subroutine opens stations dependent files with corrections
c due to atmospheric loading computed by HG Scherneck
c ftp: gere.oso.chalmers.se/pub/hgs/apload (129.16.208.42)
c Means of corrections have to be corrected out of file
c atmean97 (original apmean_upto971031_olfg.aplc)
c
c about 10 values around the session day are collected, and have
c to interpolated later in station.for for each observation.
c
c model    string3 to determine which model should be used
c          'atm','oce','soi' or 'sno'
c stan     station's name IN
c i        record in file 20 for station IN
c tmjddec  vector with times of corrections OUT
c dxatm    ... OUT
c dyatm    ... OUT
c dzatm    ... OUT
c
c units are meters.
c
c      Tesmer 9.11.2000
c   LAST REVISION 24.01.2001 (OT)   STANOP  has 19 characters instead of 12
C                                   to read the data from separate directory

c 020516 jboehm made the program usable for other models as well
c               all changes after
c               c+ jboehm and before c- jboehm
c               distinct names for atm are kept for simplicity
c               mind: the mean files are called ATMMEAN or OCEMEAN or ..
c                     (not atmean97 !!)
C
C*******************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      CHARACTER*8 STAN,QUASAR,STAT,stmean
      character*8 sta_list
      CHARACTER*19 stanop
      character*12 name
c+ jboehm ATM, SNO, SOI or OCE
      character*3 model
c- jboehm

      dimension tmjddec(30), dratm(30), deatm(30), dnatm(30)
      dimension dxatm(30), dyatm(30), dzatm(30)

      COMMON /MATH/ PI, TWOPI, ROG, ROS, ROZ, TOL
      COMMON /PHYS/ C, FL1, AH, AU, J20

c ...
      do j=1,30

        tmjddec(j)=0.
        dratm(j)=0.
        deatm(j)=0.
        dnatm(j)=0.
        dxatm(j)=0.
        dyatm(j)=0.
        dzatm(j)=0.

      end do


      dr = 0.d0
      de = 0.d0
      dn = 0.d0

      ameanr = 0.d0
      ameane = 0.d0
      ameann = 0.d0

c construct stationspecific correctionfile stationname.txt
c+ jboehm
      write(name,60) stan,'.',model
c- jboehm
60    format(8a)


c find real filename from stationame.txt
      ilfnam=7
c+ jboehm
      stanop(1:7)='../'//model//'/'
c- jboehm
      stanop(8:19)='            '

         sta_list=name(1:8)

         do k=1,8
         if (sta_list(k:k).eq.' ') then
           ij=0
           do kk=k+1,8
             if (sta_list(kk:kk).ne.' ') then
               ij=ij+1
             end if
           end do
           if (ij.eq.0) then
             goto 199
           else
             do kk=k,8-ij
               sta_list(kk:kk)='_'
             end do
           end if
         end if
         end do
199     continue
        name(1:8)=sta_list

      do infnam=1,12
         if(name(infnam:infnam).ne.' ') then
           ilfnam=ilfnam+1
           stanop(ilfnam:ilfnam)=name(infnam:infnam)
         endif
      enddo

c open real stationspecific correctionfile station.txt
      open(45,file=stanop,status='old',err=999)

c get time of first obs in mjd.dec -> tmjd1st
      read(18,rec=1) jsta, istcat, isour
      read(19,rec=isour) jsor,quasar,ra2000,de2000,
     *                   idatj,utj0
      read(20,rec=i) jp,stat,rx,ry,rz,epoch,vx,vy,vz

C******** OT+   04.01.2006
c
c      IY = EPOCH
c      ID = (EPOCH-IY)*365.25
c      IH = 0
c      IM = 0
c      S  = 0.d0
c
c      CALL JULDAT (IY,ID,IH,IM,S,IDJEP,UTJEP,0)
c      TS = (IDATJ-IDJEP)/365.25d0
c
c      rx=rx+vx*ts
c      ry=ry+vy*ts
c      rz=rz+vz*ts
C
C******** OT-

c+ jboehm
c      write(*,*)' -> reading ',model,' file ',stanop
c- jboehm
c      write(*,*) ' '

      tmjd1st=dble(idatj)+utj0/twopi

c now find the right points of time in the file

c     read first line only to check if any corrections are there
      read(45,*,err=999,end=999) tmjdatm
      if(tmjd1st.lt.tmjdatm) go to 999


c loop over the whole file
      ikorr=1
      rewind(45)
      do icoun=1,20000  ! this reads atmloadingdata ~15 years if neces

        read(45,*,err=999,end=999) tmjdatm,dr,de,dn
c        print *, 'atm_sch; tmjdatm = ', tmjdatm
        inttatm=tmjdatm ! -> integer

        if(inttatm.eq.idatj-1) then
          if(tmjd1st-tmjdatm.le.1) then
            tmjddec(ikorr)=tmjdatm
            dratm(ikorr)=dr
            deatm(ikorr)=de
            dnatm(ikorr)=dn
            ikorr=ikorr+1
          end if
        end if

        if(inttatm.eq.idatj) then
          tmjddec(ikorr)=tmjdatm
          dratm(ikorr)=dr
          deatm(ikorr)=de
          dnatm(ikorr)=dn
          ikorr=ikorr+1
        end if

        if(inttatm.eq.idatj+1) then
          tmjddec(ikorr)=tmjdatm
          dratm(ikorr)=dr
          deatm(ikorr)=de
          dnatm(ikorr)=dn
          ikorr=ikorr+1
        end if

        if(inttatm.eq.idatj+2) then
          if(tmjdatm-tmjd1st.ge.2) then
            tmjddec(ikorr)=tmjdatm
            dratm(ikorr)=dr
            deatm(ikorr)=de
            dnatm(ikorr)=dn
            ikorr=ikorr+1
          end if
        end if

        if(inttatm.gt.idatj+2) go to 800

c close loop over the whole file
      end do

800   continue
      ikorr=ikorr-1

c+ jboehm
      or = 0.d0
      oe = 0.d0
      on = 0.d0
c- jboehm

c compute corrections as dxatm, dyatm, dzatm

      If (ikorr.gt.3) then
c open file with mean values
c+ jboehm
        open(46,file='../'//model//'/'//model//'MEAN',err=999)
c- jboehm
        do iiii=1,200
          read(46,90,end=600) stmean,ameanr,ameane,ameann
90        format(a8,4x,1x,f9.5,1x,f9.5,1x,f9.5)
          if(stmean.eq.stan) then
            or=ameanr   ! [m]
            oe=ameane   ! [m]
            on=ameann   ! [m]
            go to 600
          end if
        end do
        close(46)
        goto 999

600     continue
        close(46)

c tranform a-priori rx,ry,rz to spherical coords rr,re,rn [m]
        aprr=dsqrt(rx*rx+ry*ry+rz*rz)       ! r [m]
        aprer=datan2(rz,dsqrt(rx*rx+ry*ry)) ! phi [rad]
        aprnr=datan2(ry,rx)            ! lambda [rad]

        fac_lam=dcos(aprer)

        apre=aprer*ah          ! [rad] -> [m]
        aprn=aprnr*ah*fac_lam  ! [rad] -> [m]

        do iko=1,ikorr   ! open loop over corrections

c add the mean correction
          coffr=dratm(iko)-or ! should be the right sign
          coffe=deatm(iko)-oe
          coffn=dnatm(iko)-on


c add the final correction to the a-priori coordinates
          corrr=aprr+coffr   ! [m]
          corre=apre+coffe   ! [m]
          corrn=aprn+coffn   ! [m]

c transform the correction to cartesians dxatm, dyatm, dzatm [m]
          correr=corre/ah           ! [rad]
          corrnr=corrn/ah/fac_lam   ! [rad]

          corrx=corrr*dcos(correr)*dcos(corrnr)  ! [m]
          corry=corrr*dcos(correr)*dsin(corrnr)  ! [m]
          corrz=corrr*dsin(correr)               ! [m]

C create final values

          dxatm(iko)=corrx-rx	              ! [m]
          dyatm(iko)=corry-ry	              ! [m]
          dzatm(iko)=corrz-rz	              ! [m]

c          print *, dxatm(iko), dyatm(iko), dzatm(iko)

        end do ! close loop over corrections

        close(45)

      end if

      if (ikorr.le.3) then
c+ jboehm
c        write(*,*) 'just ',ikorr,' ',model,' points for ',stan
c- jboehm
c        write(*,*) 'corrections will be zero'
c        write(*,*) ' '
        go to 999
      end if

      go to 500

c+ jboehm
999   write(*,*) 'No ',model,' data for ',stan
c- jboehm
c      write(*,*) '=> corrections will be zero'
c      write(*,*) ' '
c      go to 499

c998   write(*,*) 'data in file for ',stan
c      write(*,*) 'exists up from ',tmjdatm
c      write(*,*) '1st obs is at  ',tmjd1st
c      write(*,*) '=> corrections will be zero !!!'
c      write(*,*) ' '
c      go to 499

c+ jboehm
c997   write(*,*)'Sorry, no mean ',model,' for ',stan
c- jboehm
c      write(*,*) '=> corrections will be zero'
c      write(*,*) ' '
c      go to 499

c+ jboehm
c995   write(*,*) 'file ',model,'MEAN not found'
c- jboehm
c      write(*,*) '=> corrections will be zero'
c      write(*,*) ' '
c      go to 499

c+ jboehm
c992   write(*,*) 'unknown troubles with ',model,' loading for ',stan
c- jboehm
c      write(*,*) '=> corrections will be zero'
c      write(*,*) ' '
c      go to 499

c+ jboehm
c990   write(*,*) model,' of range at ',stan
c- jboehm
c      write(*,*) 'last correction is at   ',tmjdatm
c      write(*,*) 'first observation is at ',tmjd1st
c      write(*,*) '=> corrections will be zero'
c      write(*,*) ' '


499   continue



      do iii=1,5

        tmjddec(iii)=0.
        dxatm(iii)=0.
        dyatm(iii)=0.
        dzatm(iii)=0.

      end do

500   continue


      

c output for checking
c      do j=1,ikorr
c        write(*,*) stan,tmjddec(j),dxatm(j),dyatm(j),dzatm(j)
c      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine atmo_sch (datj,istatm,atmsta,atmmjd,
     *                     atmx,atmy,atmz,
     *                     sta,dxat_i,dyat_i,dzat_i)

c subroutine to interpolate some corrections due to atmo loading
c computed by scherneck for all timepointas of observations
c
c called by cposit (=correc)
c
c IN:   sta                 name of station to be handled
c       datj                time of obs in mjd.dec
c       istatm 	            nstations with corrections
c       atmsta(1:8*istatm)  names of stations with corr in row
c       atmmjd(istatm,j)    j times of corrections of stations
c       atmx(istatm,j)      j dx corrections of stations
c       atmy(istatm,j)      j dy corrections of stations
c       atmz(istatm,j)      j dz corrections of stations
c
c OUT:  dxat_i
c       dyat_i
c       dzat_i
c
c tesmer 10.11.00

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'OCCAM_N.FI'   ! To call 'nstations'

      character*(8*nstations) atmsta
      character*8 sta,stacor
      dimension atmmjd(nstations,30)
      dimension atmx(nstations,30),atmy(nstations,30),atmz(nstations,30)
      dimension timeatm(30), dxatm(30), dyatm(30), dzatm(30)

c find the stations name

      ns=0
      do i=1,istatm
        if(sta(1:8).eq.atmsta(((i*8)-7):(i*8))) then
          ns=i
        end if
      end do

c no correction for station ?
      if (ns.eq.0) go to 100

c put info in nice vectors
      ncor=0
      do j=1,30
        if(atmmjd(ns,j).ne.0.) ncor=j
      end do

      do j=1,ncor
        timeatm(j)=atmmjd(ns,j)
        dxatm(j)=atmx(ns,j)
        dyatm(j)=atmy(ns,j)
        dzatm(j)=atmz(ns,j)
      end do

c do the interpolation
c in:  datj
c      timeatm(ncor)
c      dxatm(ncor)
c      dyatm(ncor)
c      dzatm(ncor)
c
c out: dxat_i=0.
c      dyat_i=0.
c      dzat_i=0.

      call lagint_1 (timeatm,dxatm,ncor,datj,dxat_i)
      call lagint_1 (timeatm,dyatm,ncor,datj,dyat_i)
      call lagint_1 (timeatm,dzatm,ncor,datj,dzat_i)

      go to 200

100   continue
      dxat_i=0.
      dyat_i=0.
      dzat_i=0.

200   continue


c      print *, datj, dxat_i, dyat_i

c output for checking
c      write(*,300) sta,datj,dxat_i
c300   format(a8,1x,f16.8,1x,f16.13)
      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C************************************************************
C
      SUBROUTINE LAGINT_1 (X,Y,N,XINT,YOUT)
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
c      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
c      INTERPOLATION IS CHANGED FROM 4 DATAPONITS TO 6
c      FOR THE MOMENT. IT IS NOT CLEAR WHAT IS THE BEST
c      APPROACH (ALSO FOR EOP INTERPOLATION). MAYBE IT
c      DOES NOT MATTER AT ALL, BECAUSE THE INTERPOLATION
c      IS ALWAYS DONE 'AROUND' THE DATA POINT ITSELF.
c
c      TESMER 10.11.00
c
c      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
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

      DO 20 M = K-1,K+2    ! original
c      DO 20 M = K-2,K+3

        TERM = Y(M)

        DO 10 J = K-1,K+2  ! original
c        DO 10 J = K-2,K+3

          IF ( M .NE. J ) THEN
            TERM = TERM * (XINT - X(J))/(X(M) - X(J))
          END IF
   10   CONTINUE
        YOUT = YOUT + TERM
   20 CONTINUE
      RETURN
      END





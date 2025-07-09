C     Last change:  JB    5 Dec 2001    8:27 am

      subroutine read_gra (stan,tmjddec,grnh,greh,grnw,grew,iavail)

C***********************************************************************

      implicit double precision (a-h,o-z)

      CHARACTER*8 stan, quasar
      character*8 sta_list

      dimension tmjddec(8),grnh(8),greh(8),grnw(8),grew(8)

      COMMON /MATH/ PI, TWOPI, ROG, ROS, ROZ, TOL

c  initialize
      iavail = 0
      do i = 1,8
        tmjddec(i) = 0.
        grnh(i)   = 0.
        greh(i)   = 0.
        grnw(i)   = 0.
        grew(i)   = 0.
      end do

c  filename of rmf data
      do ista = 1,8
        if (stan(ista:ista).eq.' ') goto 99
      end do
99    ista = ista - 1

         sta_list=stan

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
199     k=k-1
        ista=k

c  open file with imf data for the station
      open(45,file=('../GRA/'//sta_list(1:ista)//'.GRA'),
     .        status='old',err=999)

c  get time of first observation in mjd
      read(18,rec=1)     jsta,istcat,isour
      read(19,rec=isour) jsor,quasar,ra2000,de2000,
     *                   idatj,utj0
      tmjd1st = dble(idatj) + utj0/twopi

c  read first line to get information

      read(45,*,err=999,end=999) tmjd0,tmjd1,tint

      if ( ( (tmjd1st     -2.d0*tint).lt.tmjd0 ) .or.
     .     ( (tmjd1st+1.d0+2.d0*tint).gt.tmjd1 ) ) go to 999

c loop over the whole file
      do
        read (45,*,err=999,end=999) tmjd,rnh,reh,rnw,rew
        if (tmjd.gt.(tmjd1st-2.0d0*tint)) then
          tmjddec(1)= tmjd
          grnh(1) = rnh
          greh(1) = reh
          grnw(1) = rnw
          grew(1) = rew
          do i = 2,8
            read (45,*,err=999,end=999)
     .            tmjddec(i),grnh(i),greh(i),grnw(i),grew(i)
          end do
          iavail = 1
          goto 111
        end if
      end do
111   continue
      goto 990

999   write(*,*) 'No GRA data for ',stan,tmjd

c      goto 990
c998   write(*,*) 'GRA FILE OUT OF RANGE FOR ',stan

990   continue

      end subroutine
C     Last change:  JB    5 Dec 2001    8:27 am

      subroutine read_imf (stan,tmjddec,dz200,dsmfw3,dazim,delev,iavail)

C***********************************************************************

      implicit double precision (a-h,o-z)

      CHARACTER*8 stan, quasar
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCC                changes begin            CCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      character*8 sta_list
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCC                changes end              CCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      dimension tmjddec(8), dz200(8), dsmfw3(8), dazim(8), delev(8)

      COMMON /MATH/ PI, TWOPI, ROG, ROS, ROZ, TOL

c  initialize
      iavail = 0
      do i = 1,8
        tmjddec(i) = 0.
        dz200(i)   = 0.
        dsmfw3(i)  = 0.
        dazim(i)   = 0.
        delev(i)   = 0.
      end do

c  filename of imf data
      do ista = 1,8
        if (stan(ista:ista).eq.' ') goto 99
      end do
99    ista = ista - 1

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCC                changes begin            CCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         sta_list=stan

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
199     k=k-1
        ista=k

c  open file with imf data for the station
      open(45,file=('../IMF/'//sta_list(1:ista)//'.IMF'),
     .        status='old',err=999)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCC                changes end              CCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c  get time of first observation in mjd

      read(18,rec=1)     jsta,istcat,isour
      read(19,rec=isour) jsor,quasar,ra2000,de2000,
     *                   idatj,utj0
      tmjd1st = dble(idatj) + utj0/twopi
c  read first line to get information
      read(45,*,err=999,end=999) tmjd0,tmjd1,tint
      if ( ( (tmjd1st     -2.d0*tint).lt.tmjd0 ) .or.
     .     ( (tmjd1st+1.d0+2.d0*tint).gt.tmjd1 ) ) go to 999

c loop over the whole file
      do
        read (45,*,err=999,end=999) tmjd,dz2,dsm,daz,del
        if (tmjd.gt.(tmjd1st-2.0d0*tint)) then
          tmjddec(1)= tmjd
          dz200(1)  = dz2
          dsmfw3(1) = dsm
          dazim(1)  = daz
          delev(1)  = del
          do i = 2,8
            read (45,*,err=999,end=999)
     .            tmjddec(i),dz200(i),dsmfw3(i),dazim(i),delev(i)
          end do
          iavail = 1
          goto 111
        end if
      end do
111   continue
      goto 990

999   continue

c      goto 990
c998   write(*,*) 'IMF FILE OUT OF RANGE FOR ',stan

990   continue

      end subroutine



C     Last change:  JB    5 Dec 2001    8:27 am

      subroutine read_rmf (stan,tmjddec,dah,dbh,dch,daw,dbw,dcw,iavail)

c  28,January, 2004 (OT) - no screen output of missing data

C***********************************************************************

      implicit double precision (a-h,o-z)

      CHARACTER*8 stan, quasar

      dimension tmjddec(8),dah(8),dbh(8),dch(8),daw(8),dbw(8),dcw(8)

      COMMON /MATH/ PI, TWOPI, ROG, ROS, ROZ, TOL

c  initialize
      iavail = 0
      do i = 1,8
        tmjddec(i) = 0.
        dah(i)   = 0.
        dbh(i)   = 0.
        dch(i)   = 0.
        daw(i)   = 0.
        dbw(i)   = 0.
        dcw(i)   = 0.
      end do

c  filename of rmf data
      do ista = 1,8
        if (stan(ista:ista).eq.' ') goto 99
      end do
99    ista = ista - 1

c  open file with rmf data for the station
      open(45,file=('../RMF/'//stan(1:ista)//'.RMF'),
     .        status='old',err=999)

c  get time of first observation in mjd
      read(18,rec=1)     jsta,istcat,isour
      read(19,rec=isour) jsor,quasar,ra2000,de2000,
     *                   idatj,utj0
      tmjd1st = dble(idatj) + utj0/twopi

c  read first line to get information
      read(45,*,err=999,end=999) tmjd0,tmjd1,tint
      if ( ( (tmjd1st     -2.d0*tint).lt.tmjd0 ) .or.
     .     ( (tmjd1st+1.d0+2.d0*tint).gt.tmjd1 ) ) go to 999

c loop over the whole file
      do
        read (45,*,err=999,end=999) tmjd,ah,bh,ch,aw,bw,cw
        if (tmjd.gt.(tmjd1st-2.0d0*tint)) then
          tmjddec(1)= tmjd
          dah(1) = ah
          dbh(1) = bh
          dch(1) = ch
          daw(1) = aw
          dbw(1) = bw
          dcw(1) = cw
          do i = 2,8
            read (45,*,err=999,end=999)
     .            tmjddec(i),dah(i),dbh(i),dch(i),daw(i),dbw(i),dcw(i)
          end do
          iavail = 1
          goto 111
        end if
      end do
111   continue
      goto 990

999   continue !  (OT) 28.01.2004

ccc   write(*,*) 'No RMF data for ',stan   ! (OT) 28.01.2004

c      goto 990
c998   write(*,*) 'RMF FILE OUT OF RANGE FOR ',stan

990   continue

      end subroutine



C     Last change:  JB    5 Dec 2001    8:27 am

      subroutine read_vm1 (stan,tmjddec,dah,daw,hd,wd,
     *        press1,tempe1,iavail)

C***********************************************************************

      implicit double precision (a-h,o-z)

      CHARACTER*8 stan, quasar
      character*8 sta_list

      dimension tmjddec(8),dah(8),daw(8),press1(8),tempe1(8)
      dimension hd(8), wd(8)      !  7-Nov-2019

      COMMON /MATH/ PI, TWOPI, ROG, ROS, ROZ, TOL

c  initialize
      iavail = 0


      do i = 1,8
        tmjddec(i) = 0.d0
        dah(i)   = 0.d0
        daw(i)   = 0.d0
        press1(i) = 0.d0   !   24-Mar-2014
        tempe1(i) = 0.d0   !   24-Mar-2014
      end do


      
c  filename of rmf data
      do ista = 1,8
        if (stan(ista:ista).eq.' ') goto 99
      end do
99    ista = ista - 1

         sta_list=stan

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
199     k=k-1
        ista=k

c  open file with imf data for the station

      open(45,file=('../VM1/'//sta_list(1:ista)//'.VM1'),
     .        status='old',err=999)

c  get time of first observation in mjd

      read(18,rec=1)     jsta,istcat,isour
      read(19,rec=isour) jsor,quasar,ra2000,de2000,
     *                   idatj,utj0


      tmjd1st = dble(idatj) + utj0/twopi

      
c  read first line to get information

      read (45,*,err=999,end=999) tmjd0,tmjd1,tint

c      print *, ' tint = ', tint
c      print *, tmjd0, tmjd1

      if (tmjd1.eq.0.d0) print 1999

      if ( ( (tmjd1st     -2.d0*tint).lt.tmjd0 ) .or.
     .     ( (tmjd1st+1.d0+2.d0*tint).gt.tmjd1 ) ) go to 999

c loop over the whole file
      do
cc        read (45,*,err=999,end=999) tmjd,ah,aw



        read (45,600,err=999,end=999) tmjd,ah,aw,hd1,wd1,pres1,temp1       !  24-Mar-2014  pressure and temperature are read

                
        if (tmjd.gt.(tmjd1st-2.0d0*tint)) then


          tmjddec(1)= tmjd
          dah(1) = ah
          daw(1) = aw

          hd(1) = hd1
          wd(1) = wd1

          do i = 2,8

c
c  7-Nov-2019 Reading of hd and wd added
c


            read (45,600,err=999,end=999)
     .       tmjddec(i),dah(i),daw(i),hd(i),wd(i),press1(i),tempe1(i)
           
          
          end do

          iavail = 1

          goto 111

        end if
      end do
111   continue
      goto 990

600   format (f8.2,2x,f10.8,2x,f10.8,2(1x,f7.4),9x,f7.2,2x,f7.2)

999   write(*,*) 'No VM1 data for ',stan,tmjd

c      goto 990
c998   write(*,*) 'VMF FILE OUT OF RANGE FOR ',stan

990   continue

1999  format (' Problem with VM1 file reading!! Stop the code, idiot! ')

      end subroutine
C     Last change:  JB    5 Dec 2001    8:27 am

      subroutine read_vmf (stan,tmjddec,dah,daw,iavail)

C***********************************************************************

      implicit double precision (a-h,o-z)

      CHARACTER*8 stan, quasar
      character*8 sta_list
      dimension tmjddec(8),dah(8),daw(8)

      COMMON /MATH/ PI, TWOPI, ROG, ROS, ROZ, TOL

c  initialize
      iavail = 0
      do i = 1,8
        tmjddec(i) = 0.
        dah(i)   = 0.
        daw(i)   = 0.
      end do

c  filename of rmf data
      do ista = 1,8
        if (stan(ista:ista).eq.' ') goto 99
      end do
99    ista = ista - 1

         sta_list=stan

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
199     k=k-1
        ista=k

c  open file with imf data for the station
      open(45,file=('../VMF/'//sta_list(1:ista)//'.VMF'),
     .        status='old',err=999)

c  get time of first observation in mjd
      read(18,rec=1)     jsta,istcat,isour
      read(19,rec=isour) jsor,quasar,ra2000,de2000,
     *                   idatj,utj0
      tmjd1st = dble(idatj) + utj0/twopi

c  read first line to get information
      read(45,*,err=999,end=999) tmjd0,tmjd1,tint

      if ( ( (tmjd1st     -2.d0*tint).lt.tmjd0 ) .or.
     .     ( (tmjd1st+1.d0+2.d0*tint).gt.tmjd1 ) ) go to 999

c loop over the whole file
      do
        read (45,*,err=999,end=999) tmjd,ah,aw
        if (tmjd.gt.(tmjd1st-2.0d0*tint)) then
          tmjddec(1)= tmjd
          dah(1) = ah
          daw(1) = aw
          do i = 2,8
            read (45,*,err=999,end=999)
     .            tmjddec(i),dah(i),daw(i)
          end do
          iavail = 1
          goto 111
        end if
      end do
111   continue
      goto 990

999   write(*,*) 'No VMF data for ',stan,tmjd

c      goto 990
c998   write(*,*) 'VMF FILE OUT OF RANGE FOR ',stan

990   continue

      end subroutine





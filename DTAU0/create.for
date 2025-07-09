

      subroutine create_vm1 (ista, sta_list_orig, iye, idoy0)

C  ista      number of stations
C  sta_list  list with station names (a8)
C  iye       year of experiment (yyyy)
C
C  This yearly files for VMF1 have to be available in ../VM1_ALL/ .
C  This routine extracts station-wise yearly files and
C  puts them into ../VM1/ .
C  2005-Oct-28 , jboehm

Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit double precision (a-h,o-z)

      parameter (n0 = 30000)  

      COMMON /PHYS/ C, FL1, AH, AU, J20
      COMMON /MATH/ PI, TWOPI, ROG, ROS, ROZ, TOL

      character*8 sta_list(20),sta_list_orig(20),str8
      character fname*30,str90*90, str_do(n0)*90

      call JULDAT (iye,idoy0,0,0,0.d0,IDJ,UTJ,0)

      tmjd1 = idj - 5.d0
      tmjd2 = idj + 5.d0

C  get the filename and open it

      write (fname,100) iye

100   format ('../VM1_ALL/vlbi',i4.4,'.vm1')            !  Change by OT on 21.02.2006

      do k = 1,30
        if (fname(k:k).eq.' ') goto 99
      end do

99    k = k - 1

      open (11,file = fname(1:k),status='old',err=898)


      read (11,'(a90)') str90

      do

        read (11,'(a90)') str90

        if (str90(1:1).ne.'#') goto 399

      end do

399   print *,'creating auxiliary file'
      ifound = 0

cc      open (12,file = 'auxiliar.dat')

      km = 0

      do

        read (11,'(A90)',end=799) str90
        read (str90,'(9x,f8.2)') tmjd

        if (ifound.eq.0) tmjd01 = tmjd

        ifound = 1

        

        if ((tmjd.ge.tmjd1).and.(tmjd.le.tmjd2)) then


          km = km + 1

          str_do(km) = str90 
c          print *, km, str_do(km)

cc          write (12,'(A90)') str90

          if (km.ge.n0) print *, km, n0

        end if

        tmjd02 = tmjd

      end do



799   continue

      if (tmjd01.gt.tmjd1) tmjd1 = tmjd01
      if (tmjd02.lt.tmjd2) tmjd2 = tmjd02

c      print *, tmjd1, tmjd2

      close (11)
cc      close (12)





cc      open (11,file = 'auxiliar.dat')
     
C  loop over stations


      do i = 1,ista

c initialise sta_list


        sta_list(i)=sta_list_orig(i)

cc        rewind(11)

        ifound = 0

        do k = 1,8

          if (sta_list(i)(k:k).eq.' ') then
            ij=0
            do kk=k+1,8
              if (sta_list(i)(kk:kk).ne.' ') then
                ij=ij+1
              end if
            end do


            if (ij.eq.0) then
              goto 199
            else
              do kk=k,8-ij
                sta_list(i)(kk:kk)='_'
              end do
            end if


          end if

        end do

199     k = k - 1

        print *,'for ',sta_list(i)

        open (12,file=('../VM1/'//sta_list(i)(1:k)//'.VM1'))

        write (12,'(F8.2,1X,F8.2,1X,F5.2)') tmjd1,tmjd2,0.25d0

C  write the file

cc        rewind(11)

c        print *, km, km1

        do j=1,km

cc          read (11,'(A90)',end=899) str90

 
 
            if (str_do(j)(1:8).eq.sta_list(i)(1:8)) then

c              print *, km, j, str_do(j)(1:8), sta_list(i)(1:8)
 

                write (12,'(A80)') str_do(j)(11:90)
   
            end if

        end do


899     continue

        close (12)

      end do

cc      close (11)

898   continue

      end subroutine
 
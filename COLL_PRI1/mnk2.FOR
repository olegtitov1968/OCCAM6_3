      subroutine mnk2 (npar0, npar1, nobs, iglob,
     * cov, l, a, c, nnr, covm8, v8, ier2 ) 

c***************************************************************************
c
c   Least squares collocation method
c
c
C   Revision  18, July, 2001, (OT)  Maximum of observables is 3310
C                                        instead of 2000
c   Revision: 26.06.2002 (OT)  4000 obs., 101 parameters for NNT network
c   Revision: 10, July, 2002 (OT)  Maximum 11 stations (101 parameters) !!
c   Revision: 15, May, 2003 (OT)
c   Revision: 31, May, 2004 (OT)  Reading OCCAM.FI
c                                      20 stations, 5000 observations
c   Revision: 31, May, 2004 (OT)  Optimization
c   Revision: 01, March, 2007 (OT)  Nine extra global parameters
c                       (3 for Galaxy, 5 for Grav and 1 for gamma) instead of four
c   Revision: 25, June, 2007 (OT)  NNR added
c   Last revision: 18, March, 2008 (OT)  14 extra global parameters
c      (3 for Galaxy, 10 for Grav and 1 for gamma)
c   Last revision: 05, September, 2009 (OT)  17 extra global parameters
c      (3 for Galaxy, 3 for rotation, 10 for Grav and 1 for gamma)
c
c   IER2 added to the list of output paarmeters
c
c
c  (C) Oleg Titov, Saint-Petersburg University, 1998-2001
c  (C) Oleg Titov, Geoscience Australia, 2001-2016
c
c
c**************************************************************************

      implicit double precision (a-h,o-z)
      include 'OCCAM.FI'

      real*8 l(nobs1), a(nobs1*npar), a1(nobs1*npar), c(2*nobs1*nsou)
      real*8 a2(nobs1*npar), a3(2*nsou*nobs1)
    
      real*8 covm(npar*(npar+1)/2), cov(nobs1)
      real*8 q(nobs1*(nobs1+1)/2)
  
      real*8 covm8(nsou*(2*nsou+1)), v8(2*nsou)
     
      character nnr*6         !                25.07.2007



c
c  Change on 25.07.2007
c
      if (nnr.eq.'NNRYES') then
         nobs_9=nobs-9
         nobs_3=nobs-3
      end if

      if (nnr.eq.'NNR_NO') then
         nobs_9=nobs-6   !  No NNR-constraints
         nobs_3=nobs     !
      end if
c
c  End of change on 25.07.2007
c

c      open (661, file='test_matr.dat', status = 'unknown')


      if(nobs.gt.nobs1) print *, ' error!!! The value n is too large ! '
      if(npar1.gt.npar) print *, 'error!!! The value npar is too large!'

    

      print *, 'Least squares method in progress ... ', npar1, nobs

c
c
c      T   -1
c     B * Q  * B
c
      do i=1,nobs      

         do j=1,npar1

            nji=nobs*(j-1)+i
            
            a1(nji) = a(nji)/cov(i)

         end do

      end do


      do i=1,npar1

         do j=1,i

            ij=i*(i-1)/2+j

            covm(ij)=0.d0

            do k=1,nobs          !

               nik=nobs*(i-1)+k
               njk=nobs*(j-1)+k

               covm(ij) = covm(ij) + a(nik)*a(njk)/cov(k)

            end do

            

         end do


      end do

      eps2=0.0000001d0
c
c                                   T   -1
c     invertion of the matrix Q1 = B * Q  * B
c
      call dsinv ( covm, npar1, eps2, ier2 )

      print *, ' IER2 = ',ier2

      if (ier2.ne.0) write (15,*) 'ier2 = ', ier2


c
c           -1    T   -1
c         Q1  *  B * Q
c

      do i=1,npar1

         do k=1,nobs

            nki=nobs*(i-1)+k

            a2(nki)=0.d0

            do j=1,npar1
               if(j.le.i) then
                  ji=i*(i-1)/2+j
               else
                  ji=j*(j-1)/2+i
               end if
               a2(nki) = a2(nki) + covm(ji)*a1(nobs*(j-1)+k)
            end do

         end do
      end do
c
c           -1   -1       -1    T   -1
c     Z =  Q  - Q * B * Q1  *  B * Q      !     Projectional matrix
c
      do i=1,nobs
         do j=1,i


            ij=i*(i-1)/2+j

            if (i.eq.j) then 

                q(ij) = 1.d0 / cov(i)


            else

                q(ij) = 0.d0

            end if 
                 

            do k=1,npar1

               kj=nobs*(k-1)+j
               ki=nobs*(k-1)+i

               q(ij) = q(ij) - a1(kj)*a2(ki)

            end do

         end do
      end do

c      T             T
c     C * Z * C,    C * Z * l
c
c
      print *, ' iglob  = ', iglob, 2*iglob + n14, 'nobs = ', nobs

      do i=1,nobs

         do j=1,2*iglob+n14      !  18.03.2008 iglob - The number of radiosources in one session

            nji=nobs*(j-1)+i

            a3(nji)=0.d0

            do k=1,nobs
               if(k.le.i) then
                  ki=i*(i-1)/2+k
               else
                  ki=k*(k-1)/2+i
               end if

               jk=nobs*(j-1)+k                !       T

               a3(nji)=a3(nji)+c(jk)*q(ki)    !      C * Z

            end do

         end do
      end do


      print *, 2*iglob + n14

      do i=1,2*iglob + n14        !  18.03.2008  iglob - the number of radio sources in one session

         v8(i)=0.d0

         do j=1,i

            ij=i*(i-1)/2+j
            covm8(ij)=0.d0

            do k=1,nobs
               ik=nobs*(i-1)+k
               jk=nobs*(j-1)+k
               covm8(ij)=covm8(ij)+a3(ik)*c(jk)
            end do

         end do

         do k=1,nobs

            ik=nobs*(i-1)+k
            v8(i)=v8(i)+a3(ik)*l(k)

         end do

c         print *, i, v8(i)

      end do

      return

      end

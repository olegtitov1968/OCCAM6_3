      subroutine mnk1 (npar0, np, nobs, time,
     *                cov, l, a, schi, sw, eps, nnr)
c***********************************************************************
c
c     Conventional least squares method (MNK) without separation
c     on segments
c
c     INPUT PARAMETERS
c
c     NOBS - the number of observations
c     NPAR - the number of parameters under adjustment
c     A - Matrix of partial derivatives
c     L - O-C vector
c     COV - vector of covariances (from quadratic formal errors)
c
c     OUTPUT
c
c     X - vector of parameters
c     CORR - matrix of correlations (diagonal)
c     SCHI - chi-squared value
c     SW - W.R.M.S value
c     EPS - vector of post-fit residuals
c
c   Version: 14.10.98
c   Revision: 30.07.2001
c   Revision: 03.01.2002
C   REVISION June 11, 2002, (OT)   4000 observations instead of 3310
c   Revision: 26.06.2002 (OT)  101 parameters for NNT network
c   Revision: 10.07.2002 (OT)  SUBS1, SUBS2 deleted
c   Revision: 27.05.2003 (OT)
c   Revision: 31.05.2004 (OT)  Reading OCCAM.FI
c   Revision: 25.01.2007 (OT)  Matrix CORR is commented off
c   LAST Revision: 25.07.2007 (OT)  NNR is added
c
c  (C) Oleg Titov, Saint-Petersburg University, 1998-2000
c  (C) Oleg Titov, Geoscience Australia, 2001-2007
c
c
c***********************************************************************

      implicit double precision (a-h,o-z)

      INCLUDE 'OCCAM.FI'

      real*8 l(nobs1), time(nobs1), a(nobs1*npar), vv(nobs1)
      real*8 x(npar),  eps(nobs1), s1, s0

      real*8 covm(npar*(npar+1)/2), cov(nobs1), sw, schi

      character nnr*6         !                25.07.2007
c      real*8 corr(npar*(npar+1)/2)

      if(nobs.gt.nobs1) then
        print *, 'The number of observations exceeds the limit ! '
        print *, ' nobs1 = ', nobs1
      end if
      if(np.gt.npar) then
         print *, 'The number of parameters exceeds the limit ! '
      end if

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

      print *, npar0, np, nobs

      do i=1,nobs_3
         if (cov(i).eq.0) write (29,*) i, cov(i)
      end do

c
c      T   -1         T   -1
c     A * Q  * A   ; A * Q  * l
c
      do i=1,np
         vv(i)=0.d0
         do j=1,i

            covm(i*(i-1)/2+j)=0.d0

            do k=1,nobs_3

c               if (cov(k).eq.0) print *, 'step ', k, ' cov = 0 '
               covm(i*(i-1)/2+j)=covm(i*(i-1)/2+j)+
     *         a(nobs*(i-1)+k)*a(nobs*(j-1)+k)/cov(k)

            end do
         end do
         do k=1,nobs_3
            vv(i)=vv(i)+a(nobs*(i-1)+k)*l(k)/cov(k)
         end do

         if (i.gt.npar0) then

            crit1=covm(i*(i+1)/2)*covm((i-1)*i/2)
            crit2=covm(i*(i+1)/2-1)**2

            if (crit2.ge.crit1) then
c               print *, 'mnk', i, crit1, crit2
               covm(i*(i+1)/2-1)=0.d0
            end if

         end if

c	print *, i, covm(i*(i+1)/2)
      end do
      eps2=0.0000001d0
c                                   T   -1
c     invertion of the matrix Q1 = A * Q  * A
c
      call dsinv ( covm, np, eps2, ier )
      print *, ' IER = ',ier
c
c           -1    T   -1
c     X = Q1  *  A * Q  * l
c
      do i=1,np
         x(i)=0.d0
         do j=1,np
            if(j.le.i) then
               ji=i*(i-1)/2+j
            else
               ji=j*(j-1)/2+i
            end if
            x(i)=x(i)+covm(ji)*vv(j)
         end do
      end do

c
c     eps = l - A * x
c

      do i=1,nobs_3
         eps(i)=l(i)
         do j=1,np
            eps(i)=eps(i)-a(nobs*(j-1)+i)*x(j)
         end do
      end do
c
c     OUTPUT:  Vector of parameters
c
c      do i=1,np
c         ii=i*(i+1)/2
c         do j=1,i
c            jj=j*(j+1)/2
c            corr(i*(i-1)/2+j)=covm(i*(i-1)/2+j)/dsqrt(covm(ii)*covm(jj))
c         end do
c      end do

      s0=0.d0
      s1=0.d0

      do i=1,nobs_9
         s0=s0+eps(i)**2/cov(i)
         s1=s1+1.d0/cov(i)
      end do

      schi=dsqrt(s0/dble(nobs_9-np))
      sw=dsqrt(s0/s1)

      return
      end

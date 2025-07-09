      subroutine mnk1 (npar0, npar1, nobs, time, cov, l, a, 
     * schi, sw, eps, x, sx, syst, corr)

c***********************************************************************
c
c     Conventional least squares method (MNK) without separation
c     on segments
c
c     INPUT PARAMETERS
c
c     NOBS1 - the number of observations
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
c   Revision: 31.05.2004 (OT)  Preparation of the OCCAM 6.0 version
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
c
c***********************************************************************

      implicit double precision (a-h,o-z)

      INCLUDE 'OCCAM.FI'
      real*8 l(nobs1), time(nobs1), a(nobs1*npar), vv(nobs1)
      real*8 x(npar), corr(npar*(npar+1)/2), eps(nobs1), s1, s0
      real*8 covm(npar*(npar+1)/2), cov(nobs1), sw, schi, syst(nobs1)
      real*8 sx(npar)

      if(nobs.gt.nobs1) then
        print *, 'The number of observations exceeds the limit ! '
      end if

      np=npar0

      if(np.gt.npar) then
         print *, 'The number of parameters exceeds the limit ! '
      end if

      
      

c
c      T   -1         T   -1
c     A * Q  * A   ; A * Q  * l
c
      do i=1,npar1

         vv(i)=0.d0

         do j=1,i
            covm(i*(i-1)/2+j)=0.d0
            do k=1,nobs
               if (cov(k).eq.0) print *, 'step ', k, ' cov = 0 '
               covm(i*(i-1)/2+j)=covm(i*(i-1)/2+j)+
     *         a(nobs*(i-1)+k)*a(nobs*(j-1)+k)/cov(k)
            end do
         end do

         do k=1,nobs
            vv(i)=vv(i)+a(nobs*(i-1)+k)*l(k)/cov(k)
         end do

      end do

      eps2=0.0000001d0

c      print *, ' npar1 = ', npar1
c      do i=1,npar1
c         write (77,*) i, covm(i*(i+1)/2)
c      end do

c                                   T   -1
c     invertion of the matrix Q1 = A * Q  * A
c
      call dsinv ( covm, npar1, eps2, ier )
      print *, ' IER = ',ier
c
c           -1    T   -1
c     X = Q1  *  A * Q  * l
c
      do i=1,npar1

         x(i)=0.d0

         do j=1,npar1

            if(j.le.i) then
               ji=i*(i-1)/2+j
            else
               ji=j*(j-1)/2+i
            end if

            x(i)=x(i)+covm(ji)*vv(j)

         end do


      end do

      print *, ' OK 1 ', nobs, npar1, x(npar1)

c
c     eps = l - A * x
c

         
      do i=1,nobs

         eps(i)=0.d0
         syst(i) = 0.d0
         
         
         do j=1,npar1

            syst(i) = syst(i) + a(nobs*(j-1)+i)*x(j)

                      
         end do

        
         eps(i) = l(i) - syst(i)


      end do
c
c     OUTPUT:  Vector of parameters
c
      do i=1,npar1

         ii=i*(i+1)/2

         do j=1,i
            jj=j*(j+1)/2
            corr(i*(i-1)/2+j)=covm(i*(i-1)/2+j)/dsqrt(covm(ii)*covm(jj))
         end do

      end do

      s0=0.d0
      s1=0.d0
      schi =0.d0
      sw = 0.d0

      

      do i=1,nobs

         s0=s0+eps(i)**2/cov(i)
         s1=s1+1.d0/cov(i)

      end do

      schi = dsqrt(s0/dble(nobs-npar1))

      do i=1,npar1
         
         sx(i) = schi * dsqrt(covm(i*(i+1)/2))

       
      end do
    


      sw=dsqrt(s0/s1)

      return
	end

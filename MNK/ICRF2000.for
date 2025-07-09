      subroutine ICRF2000 (source, ng, num1, isor, ra2000, de2000, 
     *                                fulsor, num, in, ra, de)
C**********************************************************************
C
C    SUBROUTINE ICRF2000
C
C    THIS ROUTINE READS the CATALOGUE ICRF2000 to create the list of
C    sources for every session
C
C    REVISION 27, MAY, 2003 (OT)
C    REVISION 13, NOVEMBER, 2003 (OT)
C    REVISION 27, NOVEMBER, 2003 (OT)
C
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
C**********************************************************************

      INCLUDE 'OCCAM.FI'

      integer*2 number0, im, in, isor, islg
      integer*4 ng(nss), num(nsou)

      character*8 source, ss, fulsor(nss)
      real*8 ra(nss), de(nss), ra2000, de2000

      open (401, file='list_of_source', status='unknown')  

      number0 = 0
      islg = 0

      DO I=1,isor

         IF (SOURCE.EQ.FULSOR(I)) then     !  this source is not new in this session
           
           islg=1
           in=i

           exit           
         end if
      END DO

      do i=1, num1 

         read (401,100) ss

         if (source.eq.ss) then

            number0 = i        ! the order number of this source in ICRF2
            ng(i)=ng(i)+1      ! The amount of observations for one radiosource

c            print *, i, number0, ng(i)
            exit

         end if

      end do




      if (islg.eq.0) then       !  new source in this session

         ISOR=ISOR+1            !  for NUM(isor)
         im=im+1                !  for in

         FULSOR(ISOR) = SOURCE
         ra(isor)=ra2000
         de(isor)=de2000

c         print *, isor, ra(isor)
        
         num(isor)=number0      ! the order number of this source in ICRF2
         in=im                  ! the order number of this source in this session

      end if

      if (number0.eq.0) write (*,101) i, source
 100  format(3x, a8)

 101  format ('icrf2000   ', i4,2x,a8)

      close (401)

      return
      end


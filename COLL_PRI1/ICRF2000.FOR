      subroutine ICRF2000 (source,  ng, isor, ra2000, de2000,
     *             fulsor, num, number, in, ra, de)
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
C    REVISION 20, MAY, 2005 (OT)
C    REVISION 01, MARCH, 2006 (OT)  maximum 250  sources per session is available
C    LAST REVISION 13, FEBRUARY, 2007 (OT) Mean epoch for each sources is read
C                         from ICRF2000
C
C**********************************************************************


      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'OCCAM.FI'

      integer*2 number, im, islg, isor
      integer*4 ng(nss), num(nsou), in
      character*8 source, ss, fulsor(nss)
      real*8 ra(nss), de(nss)

      open (10, file='ICRF2', status='unknown')   ! Check the number of sources (nss) in OCCAM.FI!!!

      number = 0 
      im = 0 
 
c      print *, ' nss in ICRF2000.FOR =  ', nss


      do i=1,nss   !   The number of lines in icrf2000

         read (10,100) ss

         if (source.eq.ss) then	

            number = i

            ng(i)=ng(i)+1      ! The amount of observations for one radiosource
            exit

         end if
      end do

c      print *, nss

      islg=0


      DO I=1,nss

         IF (SOURCE.EQ.FULSOR(I)) then

           islg=1
           in=i

           exit

         end if

      END DO


      if (islg.eq.0) then

         ISOR=ISOR+1               !  An order number of the source in the current session

         im = im + 1 
      
  
         FULSOR(ISOR) = SOURCE



         ra(isor) = ra2000
         de(isor) = de2000
         
         num(isor)=number          ! An order number of the source in the master catalogue

c         print *, ' im  = ', im, source, isor

 

      

         if (im.eq.1) then

            in = isor

         else

            in = im                     ! A number of sources in list

         end if

c         print *, ' in = ', in, isor, im


       end if

      

      if (number.eq.0) write (*,101) i, source, number

 100  format(3x,a8)
 101  format (i4,2x,a8,2x,i6)

      close (10)

      return
      end


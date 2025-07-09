************************************************************************

      SUBROUTINE kvlbi_opt (ista, fulsta, sta_list, ref_st)
C
C     SUBROUTINE kvlbi_opt
C
C     THIS SUBROUTINE WILL CREATE 'kvlbi.opt' FILE AUTOMATICALLY
C
C     1997 JULY 14  (OT)
c
c     2000-5-17, jboehm: sta_list becomes an output parameter also
C     Revision 2000, July, 13 (OT)  Zenith distance in .OPT
C     Revision 2001, May, 21 (OT)  Change in KVLBI.opt file making
C     Last revision 2006, February, 23 (OT)  Change in KVLBI.opt file making
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER FULSTA(20)*8, STA_LIST(20)*8, ref_st*8
      character down*1, model*6, v(20)*8, yes*3, no*2
      real*8 ad_sig(20), v_c(20), en(20), ee(20), eh(20)


C     OPEN A 'kvlbi.opt' FILE TO RUN Kalman filter code

      open (33,file='kvlbi.opt')

c
c   The flags for choice of the reference station
c
      k_flag=0

      do i=1,ista
         if(fulsta(i).eq.ref_st) then
            j=1
            sta_list(j)=fulsta(i)
            k_flag=1
         else
           if (k_flag.eq.0) then
              j=i+1
           else
              j=i
           end if
           sta_list(j)=fulsta(i)
         end if
         v_c(j)=0.d0
         en(j)=0.d0
         ee(j)=0.d0
         eh(j)=0.d0
      end do
c
c
c     A priori math.expectations and dispersions
c
      sig=0.01
      tr_disp=1.
      tr_grad_disp=0.01
      cl0_disp=10000.
      cl1_disp=0.01
      disp_nut=1.
      disp_EOP=1.
c
c     Random walk dispersions
c
      tr_dis=0.1
      cl0_dis=1.
      cl1_dis=0.

      write (33,101)  yes
      write (33,102)  no
      write (33,*)
      write (33,103)  cl0_disp
      write (33,104)  cl1_disp
      write (33,105)  tr_disp
      write (33,106)  tr_grad_disp
      write (33,107)  disp_nut
      write (33,108)  disp_eop
      write (33,*)
      write (33,109)  ista

      do i=1, ista
         if (i.eq.1) then
            write (33,110) sta_list(i), sig, v_c(i), en(i), ee(i),
     *       eh(i), tr_dis
         else
            write (33,111) sta_list(i), sig, v_c(i), en(i), ee(i),
     *       eh(i), tr_dis, cl0_dis, cl1_dis
         end if
      end do

 101  format ('EOP:',26x,a3)
 102  format ('STATION COORDINATES:',10x,a2)
 103  format ('Apriori disp for clock offset:',3x,f8.2,'   (m**2)')
 104  format ('Apriori disp for clock rate  :',3x,f8.2,'   (m**2)')
 105  format ('Apriori disp for wet trop. delay:',f8.2,'   (m**2)')
 106  format ('Apriori disp for trop. gradients:',f8.2,'   (m**2)')
 107  format ('Apriori disp for nutation       :',f8.2,' (arcsec**2)')
 108  format ('Apriori disp for EOP            :',f8.2,' (arcsec**2)')
 109  format ('NO. STATIONS',1x,i2,'   CONST N_EC E_EC H_EC   WET   CL_O
     *  CL_1' )
 110  format (a8,2x,f6.3,3x,5(1x,f4.2),2x,f5.2)
 111  format (a8,2x,f6.3,3x,5(1x,f4.2),2x,f5.2,2x,f4.1,1x,f8.6)


      CLOSE (33)

      RETURN
      END



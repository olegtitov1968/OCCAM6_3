      subroutine snx (stc, code, domes, st_des, lod, lom, slo,
     *                                       lad, lam, sla, hei)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  Program reads specific information from SNX.CAT file to provide an output
c  to the SINEX file
c
c  LAST CORRECTION - 20.05.2004 (OT)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      real*8 slo, sla, hei
      character stc*4, code*4, domes*9, st_des*22
      character*75 line
      open (1, file='snx.cat',  status='unknown')

      i=0

2     continue

      read (1,100,end=999) line
100   format (a75)
      read (line,101) code
101   format (1x,a4)

      i=i+1

      if (code.eq.stc) then

         read (line,102) code,domes,st_des,lod,lom,slo,lad,lam,sla,hei
102      format (1x,a4,4x,a9,3x,a22,1x,i3,1x,i2,1x,f4.1,1x,i3,1x,i2,1x,
     *          f4.1,1x,f7.1)

         go to 999
      end if

      go to 2

999   continue
      rewind 1

      return
      end

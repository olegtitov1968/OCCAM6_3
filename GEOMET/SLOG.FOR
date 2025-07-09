C
C    Function SLOG
C
C--------------------------------------------------------------------
C
C  Purpose: Compute double precision natural logarithm of 1.0d0+X
C           without the loss of accuracy
C
C           The use of built-in function DLOG may cause the loss
C             of accuracy when X is small. In this case we make
C             use of a series expansion of LN(1+X).
C
C  INPUT PARAMETER:  X
C
C  OUTPUT PARAMETER: SLOG (the value of the function)
C
C  Latest revision: 27 Dec 1991
C  Programmer:      Sergei A. Klioner
C--------------------------------------------------------------------
C
      double precision function slog(X)
      implicit none
C
C  Parameter
C
      double precision X
C
C  Local variables
C
      double precision s1,f,y
C
C  Parameter X must be greater then -1
C
      if (X.le.-1.d0) then
      print 414
414   format(1x, 'SLOG: domain error')
      stop
C
C If abs(X)<0.1, use a series expansion rather than the built-in
C   function DLOG
C
      else if (dabs(X).le.0.1d0) then
      slog=X
      s1=0.d0
      f=slog
      y=1.d0
10    f=-f*X
      s1=slog
      y=y+1.d0
      slog=slog+f/y
      if (slog.NE.s1) goto 10
C
C If abs(X)>0.1, the use of DLOG is safe
C
      else
      slog=dlog(1.d0+X)
      endif
      return
      end

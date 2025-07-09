C *******************************************************************
C
C USAGE : call DDate (idj,iy,mon,iday)
C
C PROGRAMMER: Oleg Titov
C
C ARGUMENTS: IY    - Year
C            ID    - Day of year
C            IH    - Hour
C            IM    - Minutes
C
C REQUIRED ROUTINES : None
C
C LAST REVISION: 24, February 2002 (Oleg Titov)
C
C***************************************************************************

      SUBROUTINE DDATE (idj,iy,mon,iday)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      COMMON /PHYS/ C, FL1, AH, AU, J20
      COMMON /MATH/ PI, TWOPI, ROG, ROS, ROZ, TOL


C   PASS FROM JULIAN DATE TO YEAR + DOY


        IC0 = IDJ-J20
        IC=IC0/1461

C   CHANGE THE NUMBER OF 4 YEARS CYCLES IF DATE BEFORE 2000

        IF (IC0.LT.0) IC=IC-1

        IR=IC0-IC*1461
        IY = 2000+IC*4
        IF (IR.LT.366) ID = IR+1

        IF (IR.GE.366) THEN
            IR = IR-366
            IY1 = IR/365
            IY = IY + 1 + IY1
            ID = IR - IY1*365 + 1
        ENDIF

        if (id.le.31) then
           mon=1
           iday=id
        end if
        if (mod(iy,4).eq.0) then
           if (id.gt.31.and.id.le.60) then
              mon=2
              iday=id-31
           end if
           if (id.gt.60.and.id.le.91) then
              mon=3
              iday=id-60
           end if
           if (id.gt.91.and.id.le.121) then
              mon=4
              iday=id-91
           end if
           if (id.gt.121.and.id.le.152) then
              mon=5
              iday=id-121
           end if
           if (id.gt.152.and.id.le.182) then
              mon=6
              iday=id-152
           end if
           if (id.gt.182.and.id.le.213) then
              mon=7
              iday=id-182
           end if
           if (id.gt.213.and.id.le.244) then
              mon=8
              iday=id-213
           end if
           if (id.gt.244.and.id.le.274) then
              mon=9
              iday=id-244
           end if
           if (id.gt.274.and.id.le.305) then
              mon=10
              iday=id-274
           end if
           if (id.gt.305.and.id.le.335) then
              mon=11
              iday=id-305
           end if
           if (id.gt.335.and.id.le.366) then
              mon=12
              iday=id-335
           end if
        else
           if (id.gt.31.and.id.le.59) then
              mon=2
              iday=id-31
           end if
           if (id.gt.59.and.id.le.90) then
              mon=3
              iday=id-59
           end if
           if (id.gt.90.and.id.le.120) then
              mon=4
              iday=id-90
           end if
           if (id.gt.120.and.id.le.151) then
              mon=5
              iday=id-120
           end if
           if (id.gt.151.and.id.le.181) then
              mon=6
              iday=id-151
           end if
           if (id.gt.181.and.id.le.212) then
              mon=7
              iday=id-181
           end if
           if (id.gt.212.and.id.le.243) then
              mon=8
              iday=id-212
           end if
           if (id.gt.243.and.id.le.273) then
              mon=9
              iday=id-243
           end if
           if (id.gt.273.and.id.le.304) then
              mon=10
              iday=id-273
           end if
           if (id.gt.304.and.id.le.334) then
              mon=11
              iday=id-304
           end if
           if (id.gt.334.and.id.le.365) then
              mon=12
              iday=id-334
           end if
        end if

      RETURN
      END


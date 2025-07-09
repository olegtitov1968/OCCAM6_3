      subroutine rot1 (rqu, s)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      dimension rqu(3), rq(3), tr(3,3)
      tr(1,1)=dcos(s)
      tr(2,1)=-dsin(s)
      tr(3,1)=0.d0
      tr(1,2)=dsin(s)
      tr(2,2)=dcos(s)
      tr(3,2)=0.d0
      tr(1,3)=0.d0
      tr(2,3)=0.d0
      tr(3,3)=1.d0
      do i=1,3
         rq(i)=0.d0
         do j=1,3
            rq(i)=rq(i)+tr(i,j)*rqu(j)
         end do
      end do
      do i=1,3
         rqu(i)=rq(i)
      end do
      return
      end

      subroutine rot2 (rqu, s)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      dimension rqu(3), rq(3), tr(3,3)
      tr(1,1)=dcos(s)
      tr(2,1)=0.d0
      tr(3,1)=dsin(s)
      tr(1,2)=0.d0
      tr(2,2)=1.d0
      tr(3,2)=0.d0
      tr(1,3)=-dsin(s)
      tr(2,3)=0.d0
      tr(3,3)=dcos(s)
      do i=1,3
         rq(i)=0.d0
         do j=1,3
            rq(i)=rq(i)+tr(i,j)*rqu(j)
         end do
      end do
      do i=1,3
         rqu(i)=rq(i)
      end do
      return
      end

      subroutine rot3 (rqu, s)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      dimension rqu(3), rq(3), tr(3,3)
      tr(1,1)=1.d0
      tr(2,1)=0.d0
      tr(3,1)=0.d0
      tr(1,2)=0.d0
      tr(2,2)=dcos(s)
      tr(3,2)=-dsin(s)
      tr(1,3)=0.d0
      tr(2,3)=dsin(s)
      tr(3,3)=dcos(s)
      do i=1,3
         rq(i)=0.d0
         do j=1,3
            rq(i)=rq(i)+tr(i,j)*rqu(j)
         end do
      end do
      do i=1,3
         rqu(i)=rq(i)
      end do
      return
      end


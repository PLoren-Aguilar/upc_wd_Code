      subroutine indexx(n,list,arr,indx)
c************************************************************
c                                                           *
c  Indexes particles using a quicksort algorithm.           *
c                                                           *
c  1) Description of variables in calling sequence:         *
c  ================================================         *
c  a) INPUT                                                 *
c  --------                                                 *
c  n        : number of active nodes or particles           *
c  list     : array containing particles addresses          *
c  arr      : array to be sorted                            *
c  indx     : index                                         *
c                                                           *
c  2) Quantities computed in common blocks:                 *
c  ========================================                 *
c  none                                                     *
c                                                           * 
c  3) Subroutine called                                     *
c  ====================                                     *
c  none                                                     *
c                                                           *
c************************************************************
c
      implicit none
      real*8 arr(*),fm,fa,fc,fmi,fx,a
      integer n,iq,l,ir,m,nstack
      integer i,j,indxt,jstack,p
      parameter (m=12,nstack=50,fm=7875.0d0,fa=211.0d0,fc=1663.0d0,
     &          fmi=1.0d0/fm)
c
      integer list(n), istack(nstack), indx(n)
 
      jstack=0
      l=1
      ir=n
      fx=0.0d0
      do 34 j=1,n
         indx(j)=j
34    continue
10    if(ir-l.lt.m)then
         do 36 j=l+1,ir
            indxt=indx(j)
            a=arr(list(indxt))
            do 35 i=j-1,1,-1
               if(arr(list(indx(i))).le.a)go to 11
               indx(i+1)=indx(i)
35          continue
            i=0
11          indx(i+1)=indxt
36       continue
!         if(jstack.eq.0)return
! CAMBIO1!!
         if(jstack.eq.0)then
           return
         end if
         ir=istack(jstack)
         l=istack(jstack-1)
         jstack=jstack-2
      else
         i=l
         j=ir
         fx=mod(fx*fa+fc,fm)
         iq=l+(ir-l+1)*(fx*fmi)
         indxt=indx(iq)
         a=arr(list(indxt))
         indx(iq)=indx(l)
20       continue
21       if(j.gt.0)then
            if(a.lt.arr(list(indx(j))))then
               j=j-1
               go to 21
            endif
         endif
         if(j.le.i)then
            indx(i)=indxt
            go to 30
         endif
         indx(i)=indx(j)
         i=i+1
22       if(i.le.n)then
            if(a.gt.arr(list(indx(i))))then
               i=i+1
               go to 22
            endif
         endif
         if(j.le.i)then
            indx(j)=indxt
            i=j
            go to 30
         endif
         indx(j)=indx(i)
         j=j-1
         go to 20
30       jstack=jstack+2
         if(jstack.gt.nstack)then
            write(*,*)'jstack greater than nstack!'
            stop
         end if
         if(ir-i.ge.i-l)then
            istack(jstack)=ir
            istack(jstack-1)=i+1
            ir=i-1
         else
            istack(jstack)=i-1
            istack(jstack-1)=l
            l=i+1
         endif
      endif
      go to 10

      end

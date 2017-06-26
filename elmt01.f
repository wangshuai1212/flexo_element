c$Id:$
      subroutine elmt01(d,ul,xl,ix,tl,s,p,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c        Shuai Wang, Darmstadt, Germany, 22.06.2017                                   
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose:
c     simple solid 2d 
c     IGA method,2D

c     Inputs:
c            Stiffness matrix

c     Outputs: stress strain component
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none
      
      include   'bdata.h'
      include   'cdata.h'
      include   'eldata.h'
      include   'eltran.h'
      include   'errchk.h'
      include   'hdata.h'
      include   'iofile.h'
      include   'prstrs.h'
      include   'comblk.h'
      include   'tdata.h'
      include   'pointer.h'
      include   'cnurb.h'
      include   'qudshp.h'
      include   'ptdat6.h'


      integer ix(*), ndf , ndm , nst , isw
      real*8  d(*), ul(ndf,nen,*), xl(ndm,*), tl(*), s(nst,nst), p(nst)
      
      real*8  sig(3),eps(3),gradeps(3,2),E,mu,alpha,cc(3,3)
      real*8  Biu(3,2),Bju(3,2),Btc(2,3)
      
      real*8  xsj,cxsj
      
      integer i,j,k,gpInd,iInd,jInd,nthreads,i1,j1
      
      real*8 valp(9),body(2)
      
      go to(1,2,3,3,2,3,2,3,2,2,2,2,3,3,3), isw
      return
      

1     continue
      call dinput(d(1),2)
      E = d(1)
      mu = d(2)
      call dinput(d(3),1)
      alpha = d(3)
      call dinput(d(4),1)
      nthreads = int(d(4))
      call dinput(d(5),2)
      d(190) = int(d(5))
      d(191) = int(d(6))
      call dinput(d(7),2)
      body(1)=d(7)
      body(2)=d(8)
      
      d(189) = 1.d0
      
      write(iow,2000) E,mu,alpha

      write(iow,2001) (nint(d(i)), i=189,191)
      
      return
      
2     continue
      return
      
3     continue

      E = d(1)
      mu = d(2)
      alpha = d(3)
      nthreads = int(d(4))
      d(190) = int(d(5))
      d(191) = int(d(6))
      d(189) = 1.d0
      body(1)=d(7)
      body(2)=d(8)
      
      cc=0.d0
      cc(1,1)=E/(1+mu)/(1-2*mu)*(1-mu)
      cc(1,2)=E/(1+mu)/(1-2*mu)*mu
      cc(3,3)=E/(1+mu)/2.d0
      cc(2,2)=cc(1,1)
      cc(2,1)=cc(1,2)    !plane strain
        
      
      s =0.d0
      p =0.d0
      
      call quadr2d(d)     !Gauss int.
      
      
      do gpInd=1,lint
        
      
        call interps2d(gpInd,xl,ix,ndm,nel,.false.)
        xsj = jac(gpInd)
        cxsj = xsj * ctan(1)
        
        sig=0.d0
        eps=0.d0
        gradeps=0.d0
        
        do i=1,nel
          
          eps(1) = eps(1) + shp2(1,i,gpInd) * ul(1,i,1)
          eps(2) = eps(2) + shp2(2,i,gpInd) * ul(2,i,1)
          eps(3) = eps(3) + shp2(1,i,gpInd) * ul(2,i,1)
     x                    + shp2(2,i,gpInd) * ul(1,i,1)
          
          gradeps(1,1) = gradeps(1,1) + shps2(1,i,gpInd) * ul(1,i,1)
          gradeps(1,2) = gradeps(1,2) + shps2(3,i,gpInd) * ul(1,i,1)
          gradeps(2,1) = gradeps(2,1) + shps2(3,i,gpInd) * ul(2,i,1)
          gradeps(2,2) = gradeps(2,2) + shps2(2,i,gpInd) * ul(2,i,1)
          gradeps(3,1) = gradeps(3,1) + shps2(3,i,gpInd) * ul(1,i,1)
     x                                + shps2(1,i,gpInd) * ul(2,i,1)
          gradeps(3,2) = gradeps(3,2) + shps2(2,i,gpInd) * ul(1,i,1)
     x                                + shps2(3,i,gpInd) * ul(2,i,1)
          
        enddo
                
        do i=1,3
          do j=1,3
           sig(i)=sig(i)+cc(i,j)*eps(j)
          enddo
        enddo
        
        if (mod(isw,3).eq.0) then
          do iInd=1,nel
          
          Biu=0.d0
          Biu(1,1)=shp2(1,iInd,gpInd)
          Biu(2,2)=shp2(2,iInd,gpInd)
          Biu(3,1)=shp2(2,iInd,gpInd)
          Biu(3,2)=shp2(1,iInd,gpInd)
          
          Btc=0.d0
          
          do i=1,2
          do j=1,3
            do k=1,3
          Btc(i,j)=Btc(i,j)+Biu(k,i)*cc(k,j)
          enddo
          enddo
          enddo
          
          do i1=1,2
          do k=1,3
          p((iInd-1)*ndf+i1) =   p((iInd-1)*ndf+i1)
     x            - sig(k)* Biu(k,i1)*xsj 
     x            + shp2(3,iInd,gpInd)*body(i1)*xsj 
          end do !k
          enddo !i1
          
          
          if (isw.eq.3) then   ! compute tangent
          
          do jInd=1,nel
          Bju=0.d0
          Bju(1,1)=shp2(1,jInd,gpInd)
          Bju(2,2)=shp2(2,jInd,gpInd)
          Bju(3,1)=shp2(2,jInd,gpInd)
          Bju(3,2)=shp2(1,jInd,gpInd)
          
          do i1=1,2
          do j1=1,2
          do k=1,3
            s((iInd-1)*ndf+i1,(jInd-1)*ndf+j1)
     x   =  s((iInd-1)*ndf+i1,(jInd-1)*ndf+j1)
     x    + Btc(i1,k)*Bju(k,j1)*cxsj
          enddo
          enddo
          enddo
          
          
          enddo !jInd
          
          endif !isw=3
          
          end do !iInd
        
        endif !if (mod(isw,3).eq.0)
        
                   if (isw.eq.8) then
                   valp(1) = eps(1)
                   valp(2) = eps(2)
                   valp(3) = eps(3)
                   valp(4) = gradeps(1,1)
                   valp(5) = gradeps(1,2)
                   valp(6) = gradeps(2,1)
                   valp(7) = gradeps(2,2)
                   valp(8) = gradeps(3,1)
                   valp(9) = gradeps(3,2)
                   
             call  lumpToNodes(ix,nel,valp,9,shp2,
     x              xsj,hr(nph),
     x              hr(nph+numnp),numnp)
                   endif
                   

         
      
      end do !gpInd

      
 
2000  format(' Isogeometric element for high order elastic',/,
     x       ' E    = ',e12.4,'        mu =  ', e12.4,/,
     x       ' alpha= ',e12.4)

2001  format(' NURBS flag: ',1i8,' Integration pattern: ',2i8) 
      end SUBROUTINE elmt01

c$Id:$
      subroutine elmt03(d,ul,xl,ix,tl,s,p,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2014: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c     Modification log                                Date (dd/mm/year)
c       Original version                                    01/11/2006
c        Shuai Wang, Darmstadt, Germany, 22.06.2017                                   
c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose:
c     ferro + elastic
c     DOF: p1 p2 phi u1 u2
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
      real*8  Biu(3,2),Bju(3,2)
      real*8  dsigde(3),dsigdu(3,2),ddelde(2),ddeldu(2,2)
      
      real*8  xsj,cxsj
      
      integer i,j,k,gpInd,iInd,jInd,nthreads,i1,j1
      
      real*8 valp(9),body(2)
      
c     define parameter for ferro and electric variables      
      real*8  Efield(2),Dfield(2),kappa,phi
      real*8  d31,d33,d51,b31,b33,b51,bbb(2,3)
      
      real*8  eps0pol(2)
      
      real*8 pol(2),gradpol(4),pol0,a(8),G11,G12,G44
      
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
      call dinput(d(9),1)
      kappa=d(9)
      call dinput(d(10),3)
      d31=d(10)
      d33=d(11)
      d51=d(12)
      call dinput(d(13),3)
      G11=d(13)
      G12=d(14)
      G44=d(15)
      call dinput(d(16),8)
      a(1)=d(16)
      a(2)=d(17)
      a(3)=d(18)
      a(4)=d(19)
      a(5)=d(20)
      a(6)=d(21)
      a(7)=d(22)
      a(8)=d(23)
      
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
      kappa=d(9)
      d31=d(10)
      d33=d(11)
      d51=d(12)
      G11=d(13)
      G12=d(14)
      G44=d(15)
      a(1)=d(16)
      a(2)=d(17)
      a(3)=d(18)
      a(4)=d(19)
      a(5)=d(20)
      a(6)=d(21)
      a(7)=d(22)
      a(8)=d(23)
            
      cc=0.d0
      cc(1,1)=E/(1+mu)/(1-2*mu)*(1-mu)
      cc(1,2)=E/(1+mu)/(1-2*mu)*mu
      cc(3,3)=E/(1+mu)/2.d0
      cc(2,2)=cc(1,1)
      cc(2,1)=cc(1,2)    !plane strain
      
      call d33_to_b33_2d03(d31,d33,d51,b31,b33,b51,cc)
      
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
        phi=0.d0
        Efield=0.d0
        Dfield=0.d0
        pol=0.d0
        gradpol=0.d0
        
        do i=1,nel
        
          pol(1) = pol(1) + shp2(3,i,gpInd)*ul(1,j,1)
          pol(2) = pol(2) + shp2(3,i,gpInd)*ul(2,j,1)
          
          Efield(1)=Efield(1)-shp2(1,i,gpInd) * ul(3,i,1)
          Efield(2)=Efield(2)-shp2(2,i,gpInd) * ul(3,i,1)
          
          eps(1) = eps(1) + shp2(1,i,gpInd) * ul(4,i,1)
          eps(2) = eps(2) + shp2(2,i,gpInd) * ul(5,i,1)
          eps(3) = eps(3) + shp2(1,i,gpInd) * ul(5,i,1)
     x                    + shp2(2,i,gpInd) * ul(4,i,1)
          
          gradeps(1,1) = gradeps(1,1) + shps2(1,i,gpInd) * ul(4,i,1)
          gradeps(1,2) = gradeps(1,2) + shps2(3,i,gpInd) * ul(4,i,1)
          gradeps(2,1) = gradeps(2,1) + shps2(3,i,gpInd) * ul(5,i,1)
          gradeps(2,2) = gradeps(2,2) + shps2(2,i,gpInd) * ul(5,i,1)
          gradeps(3,1) = gradeps(3,1) + shps2(3,i,gpInd) * ul(4,i,1)
     x                                + shps2(1,i,gpInd) * ul(5,i,1)
          gradeps(3,2) = gradeps(3,2) + shps2(2,i,gpInd) * ul(4,i,1)
     x                                + shps2(3,i,gpInd) * ul(5,i,1)
          
          gradpol(1) = gradpol(1) + shp2(1,i,gpInd)*ul(1,j,1)
          gradpol(2) = gradpol(2) + shp2(2,i,gpInd)*ul(2,j,1)
          gradpol(3) = gradpol(3) + shp2(2,i,gpInd)*ul(1,j,1)
          gradpol(4) = gradpol(4) + shp2(1,i,gpInd)*ul(2,j,1)

        enddo
        

        pol0=1.d0
        
        call b33_to_bbb03(b31,b33,b51,pol,pol0,bbb)
        print*,b31,bbb
        

                
        do i=1,3
          do j=1,3
           sig(i)=sig(i)+cc(i,j)*eps(j)
          enddo
          do j=1,2
           sig(i)=sig(i)-bbb(j,i)*Efield(j)
          enddo
        enddo
        
        do i=1,2
          Dfield(i)=Dfield(i)+kappa*Efield(i)
          do j=1,3
           Dfield(i)=Dfield(i)+bbb(i,j)*eps(j)
          enddo
        enddo
        
        
        if (mod(isw,3).eq.0) then
          do iInd=1,nel
          
          Biu=0.d0
          Biu(1,1)=shp2(1,iInd,gpInd)
          Biu(2,2)=shp2(2,iInd,gpInd)
          Biu(3,1)=shp2(2,iInd,gpInd)
          Biu(3,2)=shp2(1,iInd,gpInd)
          
!           Btc=0.d0
!           
!           do i=1,2
!           do j=1,3
!             do k=1,3
!           Btc(i,j)=Btc(i,j)+Biu(k,i)*cc(k,j)
!           enddo
!           enddo
!           enddo
          
          do k=1,2
           p((iInd-1)*ndf+3) =   p((iInd-1)*ndf+3)
     x            - shp2(k,iInd,gpInd)*Dfield(k)*xsj
          enddo
          
          
          do i1=1,2
          do k=1,3
          p((iInd-1)*ndf+3+i1) =   p((iInd-1)*ndf+3+i1)
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
          
          dsigde=0.d0
          dsigdu=0.d0
          do i1=1,3
            do j1=1,2
             dsigde(i1)= dsigde(i1)+bbb(j1,i1)*shp2(j1,jInd,gpInd)
               do k=1,3
                   dsigdu(i1,j1)=dsigdu(i1,j1)+cc(i1,k)*Bju(k,j1)
                 end do
              end do
           end do
          
         ddelde=0.d0 
         ddeldu=0.d0
         do i1=1,2
            do j1=1,2
            ddelde(i1)=ddelde(i1)-kappa*shp2(i1,jInd,gpInd)
               do k=1,3
                 ddeldu(i1,j1)=ddeldu(i1,j1)+bbb(i1,k)*Bju(k,j1)
                  end do
               end do
            end do  
           
          
c           Kuu
          do i1=1,2
          do j1=1,2
          do k=1,3
            s((iInd-1)*ndf+3+i1,(jInd-1)*ndf+3+j1)
     x   =  s((iInd-1)*ndf+3+i1,(jInd-1)*ndf+3+j1)
     x     +Biu(k,i1)*dsigdu(k,j1)*cxsj
          enddo
          enddo
          enddo
          
!           Kue  - sig(k)* Biu(k,i1)*xsj
         
          do i1=1,2
          do k=1,3
            s((iInd-1)*ndf+3+i1,(jInd-1)*ndf+3)
     x   =  s((iInd-1)*ndf+3+i1,(jInd-1)*ndf+3)
     x    + Biu(k,i1)*dsigde(k)*cxsj
          enddo
          enddo
!           
!           Keu
          do i1=1,2
          do k=1,2
            s((iInd-1)*ndf+3,(jInd-1)*ndf+3+i1)
     x   =  s((iInd-1)*ndf+3,(jInd-1)*ndf+3+i1)
     x    + shp2(k,iInd,gpInd)*ddeldu(k,i1)*cxsj
          enddo
          enddo

          
          
!           Kphiphi
         do k=1,2
           s((iInd-1)*ndf+3,(jInd-1)*ndf+3)
     x   =  s((iInd-1)*ndf+3,(jInd-1)*ndf+3)
     x    - shp2(k,iInd,gpInd)*ddelde(k)*cxsj
         enddo
           
           
          
          enddo !jInd
          
          endif !isw=3
          
          end do !iInd
        
        endif !if (mod(isw,3).eq.0)
        
                   if (isw.eq.8) then
                   valp(1) = Dfield(1)
                   valp(2) = Dfield(2)
                   valp(3) = Efield(1)
                   valp(4) = Efield(2)
                   valp(5) = sig(1)
                   valp(6) = sig(2)
                   valp(7) = eps(1)
                   valp(8) = eps(2)
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
      end SUBROUTINE elmt03
      
      
      
          subroutine d33_to_b33_2d03(d31,d33,d51,b31,b33,b51,cc)
          real*8,intent(in):: d31,d33,d51,cc(3,3)
          real*8,intent(out):: b31,b33,b51
          

         b31=d31*cc(1,1)+d33*cc(1,2)
         b33=d31*cc(1,2)+d33*cc(1,1)
         b51=d51*cc(3,3)

          end subroutine d33_to_b33_2d03
          
         subroutine b33_to_bbb03(b31,b33,b51,pol,pol0,bbb)
         real*8,intent(in):: b31,b33,b51,pol(2),pol0
         real*8,intent(out):: bbb(2,3)
         
         real*8 abspol,abspol2,abspoli,apoldp0,pol1q,pol2q
         real*8 nvec(2),nn(2,2),nnn(2,2,2),krnn(2,2),brot(2,2,2)
         
        abspol=dsqrt(pol(1)*pol(1)+pol(2)*pol(2))
        abspol2=abspol**2.d0
        abspoli=1.d0/abspol
        apoldp0=abspol/pol0
        pol1q=pol(1)**2.d0
        pol2q=pol(2)**2.d0
        
       
        
        if(abspol.ne.0.d0) then
          nvec(1)=pol(1)/abspol
          nvec(2)=pol(2)/abspol
          nn(1,1)=nvec(1)*nvec(1)
          nn(1,2)=nvec(1)*nvec(2)
          nn(2,1)=nn(1,2)
          nn(2,2)=nvec(2)*nvec(2)
          nnn(1,1,1)=nn(1,1)*nvec(1)
          nnn(1,1,2)=nn(1,1)*nvec(2)
          nnn(1,2,1)=nnn(1,1,2)
          nnn(1,2,2)=nn(1,2)*nvec(2)
          nnn(2,1,1)=nnn(1,1,2)
          nnn(2,1,2)=nnn(1,2,2)
          nnn(2,2,1)=nnn(1,2,2)
          nnn(2,2,2)=nn(2,2)*nvec(2)
          krnn(1,1)=1.d0-nn(1,1)
          krnn(1,2)=-nn(1,2)
          krnn(2,1)=-nn(2,1)
          krnn(2,2)=1.d0-nn(2,2)
          
         do k=1,2
          do i=1,2
            do j=1,2
         brot(k,i,j)=(b33*nnn(i,j,k)+b31*krnn(i,j)*nvec(k)
     x                      +b51*.5d0*(krnn(k,i)*nvec(j)
     x                      +krnn(k,j)*nvec(i)))*apoldp0
          enddo
          enddo
          enddo
          
         bbb(1,1)=brot(1,1,1)
         bbb(1,2)=brot(1,2,2)
         bbb(1,3)=brot(1,1,2)
         bbb(2,1)=brot(2,1,1)
         bbb(2,2)=brot(2,2,2)
         bbb(2,3)=brot(2,1,2)
          
       else
       
       bbb=0.d0
           
       endif   

         end subroutine b33_to_bbb03
      

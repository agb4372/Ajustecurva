      program mincua
      implicit none
      real x,y,e,indy,indyx,indyx2,a0,a0x,a0x2
      real a1,a1x,a1x2,a2,a2x,a2x2,p
      real deta0,deta1,deta2,det
      real res0,res1,res2
      integer nmax,n,u,i,io
      parameter (nmax=100,u=20)
      dimension x(nmax),y(nmax),e(nmax)
      character*50 f
      write(*,*)'Ingrese el nombre del archivo'
      read(*,*)f
      n=0
      
c     --------------------------------
c     Modulo de lectura
c     --------------------------------
      open(u,file=f,status='old')
      do i=1,nmax
         read(u,*,iostat=io)p
         
         if(io .eq. 0) then
            n=n+1   !CANTIDAD DE ELEMENTOS
            else
            continue
            end if
         end do
         
      write(*,*)n
      close (u)
      
      open(u,file=f,status='old')
      do i=1,n
         read(u,*)x(i),y(i),e(i)
         write(*,*)x(i),y(i),e(i)
      end do
      
      
      close(u)

c     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     AJUSTE POR MINIMOS CUADRADOS

      
c     --------------------------------
c     Subrutinas de calculo
c     -------------------------------- 
      call termy(x,y,n,indy,indyx,indyx2)
      call terma0(x,n,a0,a0x,a0x2)
      call terma1(x,n,a1,a1x,a1x2)
      call terma2(x,n,a2,a2x,a2x2)
c     ---------------------------------
c     Imprime el sistema de ecuaciones
c     ---------------------------------
      write(*,*)'------------------------'
      write(*,*)a2,'X**2',a1,'X',a0,'=',indy
      write(*,*)a2x,'X**2',a1x,'X',a0x,'=',indyx
      write(*,*)a2x2,'X**2',a1x2,'X',a0x2,'=',indyx2
c     ---------------------------------
c     Determinacion de cada coeficiente
c     ---------------------------------
c     --------------------------------

      call determ(a2,a2x,a2x2,a1,a1x,a1x2,a0,a0x,a0x2,det)
      call detra2(indy,indyx,indyx2,a1,a1x,a1x2,a0,a0x,a0x2,deta2)
      call detra1(a2,a2x,a2x2,indy,indyx,indyx2,a0,a0x,a0x2,deta1)
      call detra0(a2,a2x,a2x2,a1,a1x,a1x2,indy,indyx,indyx2,deta0)
      res0=deta0/det
      res1=deta1/det
      res2=deta2/det
      write(*,*)'---------------'
      write(*,*)'Ecuacion por minimos cudrados:'
      write(*,*)res2,'X**2',res1,'X',res0

      write(*,*)'---------------'
c     ---------------------------------
c     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     ---------------------------------

c     --------------------------------
c     AJUSTE MEDIANTE { u= Icm/(xcm - x) +m(xcm - x)
c     --------------------------------
      call ajuste(x,y,n,e)

      end

c     -------------------------------
c     Definicion de la subrutina del termino =y
c     ------------------------------- 

      subroutine termy(x,y,n,indy,indyx,indyx2)
      real x,y,indy,indyx,indyx2
      integer i,n
      dimension x(n),y(n)
      indy=0.
      indyx=0.
      indyx2=0.
      do i=1,n
         indy=indy+y(i)
         indyx=indyx+(y(i)*x(i))
         indyx2=indyx2+(y(i)*(x(i)**2))
      end do
      end

c     ------------------------------
c     Definicion del termino independiente
c     Factores del termino independiente a0
c     -----------------------------
      subroutine terma0(x,n,a0,a0x,a0x2)
      real x,a0,a0x,a0x2
      integer i,n
      dimension x(n)
      a0=0.
      a0x=0.
      a0x2=0.
      a0=n                      !Hay tantos a0 como puntos en la base de datos

      do i=1,n
         a0x=a0x+x(i)
         a0x2=a0x2+(x(i)**2)
      
      end do
      end

c     ------------------------------
c     Subrutina del termino a1
c     ------------------------------
      subroutine terma1(x,n,a1,a1x,a1x2)
      real x,a1,a1x,a1x2
      integer i,n
      dimension x(n)
      a1=0.
      a1x=0.
      a1x2=0.
      do i=1,n
         a1=a1+x(i)
         a1x=a1x+(x(i)**2)
         a1x2=a1x2+(x(i)**3)
      end do
      end
c     ------------------------------
c     Subrutina del termino a2
c     ------------------------------ 
      subroutine terma2(x,n,a2,a2x,a2x2)
      real x,a2,a2x,a2x2
      integer i,n
      dimension x(n)
      a2=0.
      a2x=0.
      a2x2=0.

      do i=1,n
         a2=a2+(x(i)**2)
         a2x=a2x+(x(i)**3)
         a2x2=a2x2+(x(i)**4)
      end do
      end

      subroutine determ(a2,a2x,a2x2,a1,a1x,a1x2,a0,a0x,a0x2,det)
      real a2,a2x,a2x2,a1,a1x,a1x2,a0,a0x,a0x2,det
      det=a2*((a1x*a0x2)-(a1x2*a0x))-a1*((a2x*a0x2)-(a2x2*a0x))+a0*((a2x&
     &*a1x2)-(a1x*a2x2))
      end

      subroutine detra2(indy,indyx,indyx2,a1,a1x,a1x2,a0,a0x,a0x2,deta2)
      real indy,indyx,indyx2,a1,a1x,a1x2,a0,a0x,a0x2,deta0
      deta2=indy*((a1x*a0x2)-(a0x*a1x2))-a1*((indyx*a0x2)-(indyx2*a0x)) &
     &+a0*((indyx*a1x2)-(indyx2*a1x))
      end

      subroutine detra1(a2,a2x,a2x2,indy,indyx,indyx2,a0,a0x,a0x2,deta1)
      real a2,a2x,a2x2,indy,indyx,indyx2,a0,a0x,a0x2,deta1
      deta1=a2*((indyx*a0x2)-(a0x*indyx2))-indy*((a2x*a0x2)-(a0x*a2x2)) &
     &+a0*((a2x*indyx2)-(a2x2*indyx))
      end

      subroutine detra0(a2,a2x,a2x2,a1,a1x,a1x2,indy,indyx,indyx2,deta0)
      real a2,a2x,a2x2,a1,a1x,a1x2,indy,indyx,indyx2,deta0
      deta0=a2*((a1x*indyx2)-(a1x2*indyx))-a1*((a2x*indyx2)-(a2x2*indyx)&
     &)+indy*((a2x*a1x2)-(a2x2*a1x))
      end
c     ------------------------------------------------------------------
c     COMPARACION ENTRE u y Icm/(x-xcm) + m(x-xcm) 
c     ------------------------------------------------------------------
      
      
      subroutine ajuste(x,y,n,e)
      real x,y,r,Icm,m,dif,e
      integer n,i,j,g,k,d,f
      dimension x(n),y(n),e(n)
      parameter (g=30)
      character*20 I1,Xcm1,k1
      character*200 fun,fun2
      Xcm=0.
      Icm=0.
      k=0
      open(g,file='fun.dat')
     
       do i=1,n
         do j=1,1000

         Xcm=0.
         Xcm=Xcm+0.001*j

              do f=1,1000
                 icm=0.
                 Icm=0.001*f
                 
         
                 r=Icm/(-x(i)+Xcm) + 1.834*(-x(i)+Xcm)
                 m=y(i)
                 dif=m-r
         
         
         if((abs(dif).le.e(i)).and.(icm.le.0.16).and.(xcm.le.0.6)) then
            k=k+1
            
            write(I1,'(a20)')Icm
            write(Xcm1,'(a20)')Xcm
            write(K1,'(a20)')k
            

            write(fun,*)'f',k,'(x)=',Icm,'/(-x+', Xcm,')','+','(-x+',Xcm&
     &  ,          ')','*1.834'
            fun2=trim(fun)
         else
            continue
         end if
     
      
      end do
      end do
      end do
      
      
      
      
      end
      

      
         
         

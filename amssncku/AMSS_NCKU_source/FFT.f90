

#if 0
program checkFFT
use dfport
implicit none
double precision::x
integer,parameter::N=256
double precision,dimension(N*2)::p
double precision,dimension(N/2)::s
integer::ncount,j,idum
character(len=8)::tt
tt=clock()
idum=iachar(tt(8:8))-48
p=0.0
open(77,file='prime.dat',status='unknown')
loop1:do ncount=1,N
   x=ran(idum)
   p(2*ncount-1)=x
   write(77,'(f15.3)')x
enddo loop1
close(77)
call four1(p,N,1)
do j=1,N/2
  s(j)=p(2*j)*p(2*j)+p(2*j-1)*p(2*j-1)
enddo
x=0.0
do j=1,N/2
  x=x+s(j)
enddo
s=s/x
open(77,file='power.dat',status='unknown')
do j=1,N/2
  write(77,'(2(1x,f15.3))')dble(j-1)/dble(N),s(j)
enddo
close(77)
end program checkFFT
#endif

!-------------
SUBROUTINE four1(dataa,nn,isign)
implicit none
INTEGER::isign,nn
double precision,dimension(2*nn)::dataa
INTEGER::i,istep,j,m,mmax,n
double precision::tempi,tempr
DOUBLE PRECISION::theta,wi,wpi,wpr,wr,wtemp
n=2*nn
j=1
do i=1,n,2
  if(j.gt.i)then
     tempr=dataa(j)
     tempi=dataa(j+1)
     dataa(j)=dataa(i)
     dataa(j+1)=dataa(i+1)
     dataa(i)=tempr
     dataa(i+1)=tempi
  endif
  m=nn
1 if ((m.ge.2).and.(j.gt.m)) then
  j=j-m
  m=m/2
goto 1
  endif
j=j+m
enddo
mmax=2
2  if (n.gt.mmax) then
     istep=2*mmax
     theta=6.28318530717959d0/(isign*mmax)
     wpr=-2.d0*sin(0.5d0*theta)**2
     wpi=sin(theta)
     wr=1.d0
     wi=0.d0
     do m=1,mmax,2
       do i=m,n,istep
         j=i+mmax
         tempr=sngl(wr)*dataa(j)-sngl(wi)*dataa(j+1)
         tempi=sngl(wr)*dataa(j+1)+sngl(wi)*dataa(j)
         dataa(j)=dataa(i)-tempr
         dataa(j+1)=dataa(i+1)-tempi
         dataa(i)=dataa(i)+tempr
         dataa(i+1)=dataa(i+1)+tempi
       enddo
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
     enddo
mmax=istep
goto 2
endif
return
END SUBROUTINE four1

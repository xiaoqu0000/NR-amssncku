

#include "macrodef.fh"

#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif    
!---------------------------------------------------------------------------------------------------
! copy a point of data into data target for vertext center code
!---------------------------------------------------------------------------------------------------
  subroutine pointcopy(wei,llbout,uubout,ext_out,data_out,xx,yy,zz,dv)
  implicit none
  integer,intent(in) :: wei
  integer,dimension(3),intent(in) ::ext_out
  real*8,dimension(3) :: llbout,uubout
  real*8,dimension(ext_out(1),ext_out(2),ext_out(3)),intent(inout)::data_out
  real*8,intent(in) :: xx,yy,zz,dv

  real*8,dimension(3) :: ho
  integer :: i,j,k

!sanity check
  if(wei.ne.3)then
     write(*,*)"fmisc.f90::pointcopy: this routine only surport 3 dimension"
     write(*,*)"dim = ",wei
     stop
  endif

!!!   
  if(any(ext_out == 1))then
       write(*,*)"fmisc.f90::pointcopy: meets iolated points for out data"
       write(*,*) llbout,uubout
       stop
  else
    ho = (uubout-llbout)/(ext_out-1)
  endif    
  i = idint((xx-llbout(1))/ho(1)+0.4)+1
  j = idint((yy-llbout(2))/ho(2)+0.4)+1
  k = idint((zz-llbout(3))/ho(3)+0.4)+1

  if(i<1 .or. i>ext_out(1) .or. &
     j<1 .or. j>ext_out(2) .or. &
     k<1 .or. k>ext_out(3) )then
     write(*,*)"i,j,k = ",i,j,k
     write(*,*)"ext = ",ext_out
     stop
  endif
  if(dabs(llbout(1)+(i-1)*ho(1)-xx)>ho(1)/2 .or. &
     dabs(llbout(2)+(j-1)*ho(2)-yy)>ho(2)/2 .or. &
     dabs(llbout(3)+(k-1)*ho(3)-zz)>ho(3)/2 )then    
     write(*,*)"fmisc.f90::pointcopy: llbout = ",llbout
     write(*,*)"fmisc.f90::pointcopy: ho = ",ho
     write(*,*)"fmisc.f90::pointcopy: x,y,z = ",llbout(1)+(i-1)*ho(1),llbout(2)+(j-1)*ho(2),llbout(3)+(k-1)*ho(3)
     write(*,*)"fmisc.f90::pointcopy: point = ",xx,yy,zz
     stop
   endif

  data_out(i,j,k)=dv

  return

  end subroutine pointcopy
!---------------------------------------------------------------------------------------------------
! copy a part of data from data source, for vertex center code
!---------------------------------------------------------------------------------------------------
  subroutine copy(wei,llbout,uubout,ext_out,data_out,llbin,uubin,ext_in,data_in,lcopy,ucopy)
  implicit none
  integer,intent(in) :: wei
  integer,dimension(3),intent(in) ::ext_out,ext_in
  real*8,dimension(3),intent(in) :: lcopy,ucopy
  real*8,dimension(3) :: llbout,uubout,llbin,uubin
  real*8,dimension(ext_out(1),ext_out(2),ext_out(3)),intent(inout)::data_out
  real*8,dimension(ext_in(1),ext_in(2),ext_in(3)),intent(in)::data_in

  real*8,dimension(3) :: ho,hi
  integer,dimension(3) :: illo,iuuo,illi,iuui

!sanity check
  if(wei.ne.3)then
     write(*,*)"fmisc.f90::copy: this routine only surport 3 dimension"
     write(*,*)"dim = ",wei
     stop
  endif

!!!   
  if(any(ext_out == 1))then
    if(any(ext_in == 1))then
       write(*,*)"fmisc.f90::copy: meets iolated points for both in and out data"
       write(*,*) llbin,uubin
       write(*,*) llbout,uubout
       stop
     else
       hi = (uubin-llbin)/(ext_in-1)
       ho = hi
    endif
  else
    ho = (uubout-llbout)/(ext_out-1)
    if(any(ext_in == 1))then
      hi = ho
    else
      hi = (uubin-llbin)/(ext_in-1)
      if(any(abs(hi-ho) > min(hi,ho)/2))then
         write(*,*)"fmisc.f90::copy: meets copy reqest for different numerical grid"
         write(*,*)hi,ho
         stop
      endif
    endif
  endif    
  illo = idint((lcopy-llbout)/ho+0.4)+1
  iuuo = ext_out - idint((uubout-ucopy)/ho+0.4)
  illi = idint((lcopy-llbin)/hi+0.4)+1
  iuui = ext_in - idint((uubin-ucopy)/hi+0.4)

  if(any(llbout-lcopy>ho/2) .or. any(ucopy-uubout>ho/2))then
     write(*,*)"fmisc.f90::copy: llbout = ",llbout
     write(*,*)"fmisc.f90::copy: uubout = ",uubout
     write(*,*)"fmisc.f90::copy: llbcopy = ",lcopy
     write(*,*)"fmisc.f90::copy: uubcopy = ",ucopy
     write(*,*)"fmisc.f90::copy: ho = ",ho
     write(*,*)llbout-lcopy,ucopy-uubout
     stop
  elseif(any(llbin -lcopy>hi/2) .or. any(ucopy-uubin >hi/2))then
     write(*,*)"fmisc.f90::copy: llbin = ",llbin
     write(*,*)"fmisc.f90::copy: uubin = ",uubin
     write(*,*)"fmisc.f90::copy: llbcopy = ",lcopy
     write(*,*)"fmisc.f90::copy: uubcopy = ",ucopy
     stop
  elseif(any(illo<1) .or. any(illi<1) .or. any(illo-iuuo>0) .or. any(illi-iuui>0) .or. &
         any(iuui-ext_in>0) .or. any(iuuo-ext_out>0))then
     write(*,*)"fmisc.f90::copy: illi = ",illi
     write(*,*)"fmisc.f90::copy: iuui = ",iuui
     write(*,*)"fmisc.f90::copy: illo = ",illo
     write(*,*)"fmisc.f90::copy: iuuo = ",iuuo     
     write(*,*)"fmisc.f90::copy: llbout = ",llbout
     write(*,*)"fmisc.f90::copy: uubout = ",uubout
     write(*,*)"fmisc.f90::copy: llbin = ",llbin
     write(*,*)"fmisc.f90::copy: uubin = ",uubin
     write(*,*)"fmisc.f90::copy: llbcopy = ",lcopy
     write(*,*)"fmisc.f90::copy: uubcopy = ",ucopy
     stop
   endif

  data_out(illo(1):iuuo(1),illo(2):iuuo(2),illo(3):iuuo(3))=data_in(illi(1):iuui(1),illi(2):iuui(2),illi(3):iuui(3))

  return

  end subroutine copy

!---------------------------------------------------------------------------------------
subroutine symmetry_bd(ord,extc,func,funcc,SoA)
  implicit none

!~~~~~~> input arguments
  integer,intent(in) :: ord
  integer,dimension(3),   intent(in) :: extc
  real*8, dimension(extc(1),extc(2),extc(3)),intent(in ):: func
  real*8, dimension(-ord+1:extc(1),-ord+1:extc(2),-ord+1:extc(3)),intent(out):: funcc
  real*8, dimension(1:3), intent(in) :: SoA

  integer::i

  funcc = 0.d0
  funcc(1:extc(1),1:extc(2),1:extc(3)) = func
   do i=0,ord-1
      funcc(-i,1:extc(2),1:extc(3)) = funcc(i+2,1:extc(2),1:extc(3))*SoA(1)
   enddo
   do i=0,ord-1
      funcc(:,-i,1:extc(3)) = funcc(:,i+2,1:extc(3))*SoA(2)
   enddo
   do i=0,ord-1
      funcc(:,:,-i) = funcc(:,:,i+2)*SoA(3)
   enddo

end subroutine symmetry_bd

subroutine symmetry_tbd(ord,extc,func,funcc,SoA)
  implicit none

!~~~~~~> input arguments
  integer,intent(in) :: ord
  integer,dimension(3),   intent(in) :: extc
  real*8, dimension(extc(1),extc(2),extc(3)),intent(in ):: func
  real*8, dimension(-ord+1:extc(1)+ord,-ord+1:extc(2)+ord,-ord+1:extc(3)+ord),intent(out):: funcc
  real*8, dimension(1:3), intent(in) :: SoA

  integer::i

  funcc = 0.d0
  funcc(1:extc(1),1:extc(2),1:extc(3)) = func
   do i=0,ord-1
      funcc(-i,1:extc(2),1:extc(3)) = funcc(i+2,1:extc(2),1:extc(3))*SoA(1)
      funcc(extc(1)+1+i,1:extc(2),1:extc(3)) = funcc(extc(1)-1-i,1:extc(2),1:extc(3))*SoA(1)
   enddo
   do i=0,ord-1
      funcc(:,-i,1:extc(3)) = funcc(:,i+2,1:extc(3))*SoA(2)
      funcc(:,extc(2)+1+i,1:extc(3)) = funcc(:,extc(2)-1-i,1:extc(3))*SoA(2)
   enddo
   do i=0,ord-1
      funcc(:,:,-i) = funcc(:,:,i+2)*SoA(3)
      funcc(:,:,extc(3)+1+i) = funcc(:,:,extc(3)-1-i)*SoA(3)
   enddo

end subroutine symmetry_tbd

subroutine symmetry_stbd(ord,extc,func,funcc,SoA)
  implicit none

!~~~~~~> input arguments
  integer,intent(in) :: ord
  integer,dimension(3),   intent(in) :: extc
  real*8, dimension(extc(1),extc(2),extc(3)),intent(in ):: func
  real*8, dimension(-ord+1:extc(1)+ord,-ord+1:extc(2)+ord,extc(3)),intent(out):: funcc
  real*8, dimension(2), intent(in) :: SoA

  integer::i

  funcc = 0.d0
  funcc(1:extc(1),1:extc(2),1:extc(3)) = func
   do i=0,ord-1
      funcc(-i,1:extc(2),1:extc(3)) = funcc(i+2,1:extc(2),1:extc(3))*SoA(1)
      funcc(extc(1)+1+i,1:extc(2),1:extc(3)) = funcc(extc(1)-1-i,1:extc(2),1:extc(3))*SoA(1)
   enddo
   do i=0,ord-1
      funcc(:,-i,1:extc(3)) = funcc(:,i+2,1:extc(3))*SoA(2)
      funcc(:,extc(2)+1+i,1:extc(3)) = funcc(:,extc(2)-1-i,1:extc(3))*SoA(2)
   enddo

end subroutine symmetry_stbd

subroutine symmetry_sntbd(ord,extc,func,funcc,SoA,actd)
  implicit none

!~~~~~~> input arguments
  integer,intent(in) :: ord,actd
  integer,dimension(3),   intent(in) :: extc
  real*8, dimension(extc(1),extc(2),extc(3)),intent(in ):: func
  real*8, dimension(-ord+1:extc(1)+ord,-ord+1:extc(2)+ord,extc(3)),intent(out):: funcc
  real*8, intent(in) :: SoA

  integer::i

  funcc = 0.d0
  funcc(1:extc(1),1:extc(2),1:extc(3)) = func
  if(actd==0)then
   do i=0,ord-1
      funcc(-i,1:extc(2),1:extc(3)) = funcc(i+2,1:extc(2),1:extc(3))*SoA
      funcc(extc(1)+1+i,1:extc(2),1:extc(3)) = funcc(extc(1)-1-i,1:extc(2),1:extc(3))*SoA
   enddo
  elseif(actd==1)then
   do i=0,ord-1
      funcc(1:extc(1),-i,1:extc(3)) = funcc(1:extc(1),i+2,1:extc(3))*SoA
      funcc(1:extc(1),extc(2)+1+i,1:extc(3)) = funcc(1:extc(1),extc(2)-1-i,1:extc(3))*SoA
   enddo
  else
    write(*,*)"symmetry_sntbd: not recognized actd = ",actd
  endif

end subroutine symmetry_sntbd


subroutine d2dump(wei,llb,uub,ext,data_in,data_out,gord,SoA)
  implicit none
  integer,             intent(in) :: wei,gord
  integer,dimension(3),intent(in) :: ext
  real*8, dimension(3),intent(in) :: SoA
  real*8, dimension(3) :: llb,uub
  real*8, dimension(ext(1),ext(2),ext(3)),intent(in)   ::data_in
  real*8, dimension(ext(1),ext(2)),       intent(inout)::data_out

  real*8  :: dZ
  integer :: i,j,k

!sanity check
  if(wei.ne.3)then
     write(*,*)"fmisc.f90::copy: this routine only surport 3 dimension"
     write(*,*)"dim = ",wei
     stop
  endif

  dZ = (uub(3)-llb(3))/(ext(3)-1)
  k = idint((0-llb(3))/dZ+0.4)+1
     
  if(k < 1)then
      write(*,*) "d2dump: something must be wrong"
      return
  endif

  data_out(i,j) = data_in(i,j,k)

end subroutine d2dump

#else
#ifdef Cell
!subroutine interp_2 support cell center only
!-----------------------------------------------------------------------------
!
! Interpolate function f using weights Delx, Dely and Delz
!
!-----------------------------------------------------------------------------

  subroutine interp_2(ex,f,f_int,il,iu,jl,ju,kl,ku,Dx,Dy,Dz,&
                                     ordn,SoA,symmetry)
  implicit none

!~~~~~~> Input parameters:

  integer,                             intent(in) :: ex(1:3), symmetry
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: f
  real*8,                              intent(out):: f_int
  integer,                             intent(in) :: il,iu,jl,ju,kl,ku,ordn
  real*8,                              intent(in) :: Dx,Dy,Dz,SoA(3)

!~~~~~~> Other parameters:

  integer :: j,imin,jmin,kmin
  real*8, dimension(1:ordn) :: x1a
  real*8, dimension(1:ordn,1:ordn,1:ordn) :: ya
  real*8, parameter :: ONE=1.d0
  integer, parameter :: NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2
  real*8 :: ddy,symX,symY,symZ

  symX = SoA(1)
  symY = SoA(2)
  symZ = SoA(3)
 
  imin = lbound(f,1)
  jmin = lbound(f,2)
  kmin = lbound(f,3)

  forall( j = 1:ordn ) x1a(j) = ( j - 1 )* ONE

  ya(2:ordn,2:ordn,2:ordn) = f(il+1:iu,jl+1:ju,kl+1:ku)
   
  if( il < imin .and. symmetry < OCTANT ) then
      write(*,*) 'Error in interp_2!!!'
      stop
  endif
  if( il < imin ) then
    ya(1,2:ordn,2:ordn) = f(imin,jl+1:ju,kl+1:ku)* symX
  else
    ya(1,2:ordn,2:ordn) = f(il  ,jl+1:ju,kl+1:ku)
  endif
  
  if( jl < jmin .and. symmetry < OCTANT ) then
     write(*,*) 'Error in interp_2!!!'
     stop
  endif
   
  if( jl < jmin ) then
    ya(2:ordn,1,2:ordn) = f(il+1:iu,jmin,kl+1:ku)* symY
  else
    ya(2:ordn,1,2:ordn) = f(il+1:iu,jl,kl+1:ku)
  endif

  if( kl < kmin .and. symmetry < EQUATORIAL ) then
    write(*,*)  'Error in interp_2!!!'
    stop
  endif

  if( kl < kmin ) then
   ya(2:ordn,2:ordn,1) = f(il+1:iu,jl+1:ju,kmin)* symZ
  else
   ya(2:ordn,2:ordn,1) = f(il+1:iu,jl+1:ju,kl  )
  endif

  if( il < imin .and. jl < jmin ) then
   ya(1,1,2:ordn) = f(imin,jmin,kl+1:ku)* symX * symY
  else if( il >= imin .and. jl < jmin ) then
   ya(1,1,2:ordn) = f(il,jmin,kl+1:ku)* symY
  else if( il < imin .and. jl >= jmin ) then
   ya(1,1,2:ordn) = f(imin,jl,kl+1:ku)* symX
  else 
   ya(1,1,2:ordn) = f(il,jl,kl+1:ku)
  endif
   
  if( il < imin .and. kl < kmin ) then
   ya(1,2:ordn,1) = f(imin,jl+1:ju,kmin)* symX * symZ
  else if( il >= imin .and. kl < kmin ) then
   ya(1,2:ordn,1) = f(il,jl+1:ju,kmin)* symZ
  else if( il < imin .and. kl >= kmin ) then
   ya(1,2:ordn,1) = f(imin,jl+1:ju,kl)* symX
  else
   ya(1,2:ordn,1) = f(il,jl+1:ju,kl)
  endif

  if( jl < jmin .and. kl < kmin ) then
    ya(2:ordn,1,1) = f(il+1:iu,jmin,kmin)* symY * symZ
  else if( jl >= jmin .and. kl < kmin ) then
    ya(2:ordn,1,1) =   f(il+1:iu,jl,kmin)* symZ
  else if( jl < jmin .and. kl >= kmin ) then
    ya(2:ordn,1,1) =   f(il+1:iu,jmin,kl)* symY
  else
   ya(2:ordn,1,1) = f(il+1:iu,jl,kl)
  endif

  if( il < imin ) then
   if( jl < jmin .and. kl < kmin) then
    ya(1,1,1) = f(imin,jmin,kmin)* symX * symY * symZ
   else if( jl >= jmin .and. kl < kmin ) then
    ya(1,1,1) = f(imin,jl,kmin)* symX * symZ
   else if( jl < jmin .and. kl >= kmin ) then
    ya(1,1,1) = f(imin,jmin,kl)* symX * symY
   else
    ya(1,1,1) = f(imin,jl,kl)* symX
   endif
  else
   if( jl < jmin .and. kl < kmin) then
    ya(1,1,1) = f(il,jmin,kmin)* symY * symZ
   else if( jl >= jmin .and. kl < kmin ) then
    ya(1,1,1) = f(il,jl,kmin)* symZ
   else if( jl < jmin .and. kl >= kmin ) then
    ya(1,1,1) = f(il,jmin,kl)* symY
   else
    ya(1,1,1) = f(il,jl,kl)
   endif
  endif
  
  call polin3(x1a,x1a,x1a,ya,Dx,Dy,Dz,f_int,ddy,ordn)

  if(.not.(dabs(f_int).ge.0))then
    write(*,*)"find nan in interp_2:",f_int,"inputs are:"
!    write(*,*)ya
!    write(*,*)"-----------------------------------------"
!    write(*,*)f(il:iu,jl:ju,kl:ku)
    write(*,*)Dx,Dy,Dz,symx,symy,symz,ordn
    write(*,*)il,iu,jl,ju,kl,ku,ex,symmetry
  endif

  return

  end subroutine interp_2
!---------------------------------------------------------------------------------------------------
! copy a point of data into data target for vertext center code
!---------------------------------------------------------------------------------------------------
  subroutine pointcopy(wei,llbout,uubout,ext_out,data_out,xx,yy,zz,dv)
  implicit none
  integer,intent(in) :: wei
  integer,dimension(3),intent(in) ::ext_out
  real*8,dimension(3) :: llbout,uubout
  real*8,dimension(ext_out(1),ext_out(2),ext_out(3)),intent(inout)::data_out
  real*8,intent(in) :: xx,yy,zz,dv

  real*8,dimension(3) :: ho
  integer :: i,j,k

!sanity check
  if(wei.ne.3)then
     write(*,*)"fmisc.f90::pointcopy: this routine only surport 3 dimension"
     write(*,*)"dim = ",wei
     stop
  endif

!!!   
  ho = (uubout-llbout)/ext_out
  i = idint((xx-llbout(1))/ho(1)+0.4)+1
  j = idint((yy-llbout(2))/ho(2)+0.4)+1
  k = idint((zz-llbout(3))/ho(3)+0.4)+1

  if(i<1 .or. i>ext_out(1) .or. &
     j<1 .or. j>ext_out(2) .or. &
     k<1 .or. k>ext_out(3) )then
     write(*,*)"i,j,k = ",i,j,k
     write(*,*)"ext = ",ext_out
     stop
  endif
  if(dabs(llbout(1)+(i-0.5)*ho(1)-xx)>ho(1)/2 .or. &
     dabs(llbout(2)+(j-0.5)*ho(2)-yy)>ho(2)/2 .or. &
     dabs(llbout(3)+(k-0.5)*ho(3)-zz)>ho(3)/2 )then    
     write(*,*)"fmisc.f90::pointcopy: llbout = ",llbout
     write(*,*)"fmisc.f90::pointcopy: ho = ",ho
     write(*,*)"fmisc.f90::pointcopy: x,y,z = ",llbout(1)+(i-0.5)*ho(1),llbout(2)+(j-0.5)*ho(2),llbout(3)+(k-0.5)*ho(3)
     write(*,*)"fmisc.f90::pointcopy: point = ",xx,yy,zz
     stop
   endif

  data_out(i,j,k)=dv

  return

  end subroutine pointcopy
!---------------------------------------------------------------------------------------------------
! copy a part of data from data source, for cell center code
!---------------------------------------------------------------------------------------------------
  subroutine copy(wei,llbout,uubout,ext_out,data_out,llbin,uubin,ext_in,data_in,lcopy,ucopy)
  implicit none
  integer,intent(in) :: wei
  integer,dimension(3),intent(in) ::ext_out,ext_in
  real*8,dimension(3),intent(in) :: lcopy,ucopy
  real*8,dimension(3) :: llbout,uubout,llbin,uubin
  real*8,dimension(ext_out(1),ext_out(2),ext_out(3)),intent(inout)::data_out
  real*8,dimension(ext_in(1),ext_in(2),ext_in(3)),intent(in)::data_in

  real*8,dimension(3) :: ho,hi
  integer,dimension(3) :: illo,iuuo,illi,iuui

!sanity check
  if(wei.ne.3)then
     write(*,*)"fmisc.f90::copy: this routine only surport 3 dimension"
     write(*,*)"dim = ",wei
     stop
  endif

!!!   
  ho = (uubout-llbout)/ext_out
  hi = (uubin-llbin)/ext_in
  illo = idint((lcopy-llbout)/ho+0.4)+1
  iuuo = ext_out - idint((uubout-ucopy)/ho+0.4)
  illi = idint((lcopy-llbin)/hi+0.4)+1
  iuui = ext_in - idint((uubin-ucopy)/hi+0.4)

  if(any(llbout-lcopy>ho/2) .or. any(ucopy-uubout>ho/2))then
     write(*,*)"fmisc.f90::copy: llbout = ",llbout
     write(*,*)"fmisc.f90::copy: uubout = ",uubout
     write(*,*)"fmisc.f90::copy: llbcopy = ",lcopy
     write(*,*)"fmisc.f90::copy: uubcopy = ",ucopy
     write(*,*)"fmisc.f90::copy: ho = ",ho
     write(*,*)llbout-lcopy,ucopy-uubout
     stop
  elseif(any(llbin -lcopy>hi/2) .or. any(ucopy-uubin >hi/2))then
     write(*,*)"fmisc.f90::copy: llbin = ",llbin
     write(*,*)"fmisc.f90::copy: uubin = ",uubin
     write(*,*)"fmisc.f90::copy: llbcopy = ",lcopy
     write(*,*)"fmisc.f90::copy: uubcopy = ",ucopy
     stop
  elseif(any(illo<1) .or. any(illi<1) .or. any(illo-iuuo>0) .or. any(illi-iuui>0) .or. &
         any(iuui-ext_in>0) .or. any(iuuo-ext_out>0))then
     write(*,*)"fmisc.f90::copy: illi = ",illi
     write(*,*)"fmisc.f90::copy: iuui = ",iuui
     write(*,*)"fmisc.f90::copy: illo = ",illo
     write(*,*)"fmisc.f90::copy: iuuo = ",iuuo     
     write(*,*)"fmisc.f90::copy: llbout = ",llbout
     write(*,*)"fmisc.f90::copy: uubout = ",uubout
     write(*,*)"fmisc.f90::copy: llbin = ",llbin
     write(*,*)"fmisc.f90::copy: uubin = ",uubin
     write(*,*)"fmisc.f90::copy: llbcopy = ",lcopy
     write(*,*)"fmisc.f90::copy: uubcopy = ",ucopy
     stop
   endif

  data_out(illo(1):iuuo(1),illo(2):iuuo(2),illo(3):iuuo(3))=data_in(illi(1):iuui(1),illi(2):iuui(2),illi(3):iuui(3))

  return

  end subroutine copy

!---------------------------------------------------------------------------------------
subroutine symmetry_bd(ord,extc,func,funcc,SoA)
  use sft_trace_mod
  use omp_lib
  implicit none
  integer(8) :: ts0_thr
  integer :: tid

!~~~~~~> input arguments
  integer,intent(in) :: ord
  integer,dimension(3),   intent(in) :: extc
  real*8, dimension(extc(1),extc(2),extc(3)),intent(in ):: func
  real*8, dimension(-ord+1:extc(1),-ord+1:extc(2),-ord+1:extc(3)),intent(out):: funcc
  real*8, dimension(1:3), intent(in) :: SoA

  integer :: i, j, k


! initialize ghost zone array to zero
  do k=-ord+1,extc(3)
   do j=-ord+1,extc(2)
    do i=-ord+1,extc(1)
     funcc(i,j,k) = 0.d0
    enddo
   enddo
  enddo

! copy interior values
  do k=1,extc(3)
   do j=1,extc(2)
    do i=1,extc(1)
     funcc(i,j,k) = func(i,j,k)
    enddo
   enddo
  enddo

! fill x ghost zones
  do i=0,ord-1
   do k=1,extc(3)
    do j=1,extc(2)
     funcc(-i,j,k) = funcc(i+1,j,k)*SoA(1)
    enddo
   enddo
  enddo

! fill y ghost zones
  do i=0,ord-1
   do k=1,extc(3)
    do j=-ord+1,extc(1)
     funcc(j,-i,k) = funcc(j,i+1,k)*SoA(2)
    enddo
   enddo
  enddo

! fill z ghost zones
  do i=0,ord-1
   do k=-ord+1,extc(2)
    do j=-ord+1,extc(1)
     funcc(j,k,-i) = funcc(j,k,i+1)*SoA(3)
    enddo
   enddo
  enddo


end subroutine symmetry_bd

! Inner version: orphaned OMP DO directives only, no OMP PARALLEL.
! Must be called from within an existing OMP PARALLEL region.
subroutine symmetry_bd_inner(ord,extc,func,funcc,SoA)
  use sft_trace_mod
  use omp_lib
  implicit none
  integer(8) :: ts0_thr
  integer :: tid

!~~~~~~> input arguments
  integer,intent(in) :: ord
  integer,dimension(3),   intent(in) :: extc
  real*8, dimension(extc(1),extc(2),extc(3)),intent(in ):: func
  real*8, dimension(-ord+1:extc(1),-ord+1:extc(2),-ord+1:extc(3)),intent(out):: funcc
  real*8, dimension(1:3), intent(in) :: SoA

  integer :: i, j, k, kk
  integer(8) :: ts1
  real*8 :: sx, sy, sz, sxy, sxz, syz, sxyz

  tid = omp_get_thread_num()
  ts0_thr = sft_get_ts()

  sx = SoA(1)
  sy = SoA(2)
  sz = SoA(3)
  sxy = sx * sy
  sxz = sx * sz
  syz = sy * sz
  sxyz = sx * sy * sz

! ── Loop 1: interior-z slabs (k=1:nz) ──────────────────────────
! Each k-slab fills the full x,y range of funcc:
!   interior copy  (i=1:nx,  j=1:ny )  ── pure memcpy
!   x-face ghost   (i=-ord+1:0, j=1:ny )  ── *sx
!   y-face ghost   (i=1:nx,  j=0:-ord+1)  ── *sy
!   xy-edge ghost  (i=-ord+1:0, j=0:-ord+1)  ── *sxy
!$OMP DO COLLAPSE(1)
  do k=1,extc(3)
    do j=1,extc(2)
      do i=-ord+1,0
        funcc(i,j,k) = func(1-i,j,k) * sx
      enddo
      do i=1,extc(1)
        funcc(i,j,k) = func(i,j,k)
      enddo
    enddo
    do j=0,ord-1
      do i=-ord+1,0
        funcc(i,-j,k) = func(1-i,j+1,k) * sxy
      enddo
      do i=1,extc(1)
        funcc(i,-j,k) = func(i,j+1,k) * sy
      enddo
    enddo
  enddo
!$OMP END DO NOWAIT
  ts1 = sft_get_ts()
  call sft_trace("sym_bd_inner:main", ts0_thr, ts1, 0, tid)
  ts0_thr = ts1

! ── Loop 2: ghost-z slabs (k=0,-1,...,-ord+1) ──────────────────
! Same x,y structure but every value carries an extra *sz factor:
!   z-face ghost   (i=1:nx,  j=1:ny )  ── *sz
!   xz-edge ghost  (i=-ord+1:0, j=1:ny )  ── *sxz
!   yz-edge ghost  (i=1:nx,  j=0:-ord+1)  ── *syz
!   xyz-corner     (i=-ord+1:0, j=0:-ord+1)  ── *sxyz
!$OMP DO COLLAPSE(1)
  do kk=0,ord-1
    do j=1,extc(2)
      do i=-ord+1,0
        funcc(i,j,-kk) = func(1-i,j,kk+1) * sxz
      enddo
      do i=1,extc(1)
        funcc(i,j,-kk) = func(i,j,kk+1) * sz
      enddo
    enddo
    do j=0,ord-1
      do i=-ord+1,0
        funcc(i,-j,-kk) = func(1-i,j+1,kk+1) * sxyz
      enddo
      do i=1,extc(1)
        funcc(i,-j,-kk) = func(i,j+1,kk+1) * syz
      enddo
    enddo
  enddo
!$OMP END DO NOWAIT
  ts1 = sft_get_ts()
  call sft_trace("sym_bd_inner:ghost_z", ts0_thr, ts1, 0, tid)
  ts0_thr = ts1

end subroutine symmetry_bd_inner

subroutine symmetry_tbd(ord,extc,func,funcc,SoA)
  implicit none

!~~~~~~> input arguments
  integer,intent(in) :: ord
  integer,dimension(3),   intent(in) :: extc
  real*8, dimension(extc(1),extc(2),extc(3)),intent(in ):: func
  real*8, dimension(-ord+1:extc(1)+ord,-ord+1:extc(2)+ord,-ord+1:extc(3)+ord),intent(out):: funcc
  real*8, dimension(1:3), intent(in) :: SoA

  integer::i

  funcc = 0.d0
  funcc(1:extc(1),1:extc(2),1:extc(3)) = func
   do i=0,ord-1
      funcc(-i,1:extc(2),1:extc(3)) = funcc(i+1,1:extc(2),1:extc(3))*SoA(1)
      funcc(extc(1)+1+i,1:extc(2),1:extc(3)) = funcc(extc(1)-i,1:extc(2),1:extc(3))*SoA(1)
   enddo
   do i=0,ord-1
      funcc(:,-i,1:extc(3)) = funcc(:,i+1,1:extc(3))*SoA(2)
      funcc(:,extc(2)+1+i,1:extc(3)) = funcc(:,extc(2)-i,1:extc(3))*SoA(2)
   enddo
   do i=0,ord-1
      funcc(:,:,-i) = funcc(:,:,i+1)*SoA(3)
      funcc(:,:,extc(3)+1+i) = funcc(:,:,extc(3)-i)*SoA(3)
   enddo

end subroutine symmetry_tbd

subroutine symmetry_stbd(ord,extc,func,funcc,SoA)
  implicit none

!~~~~~~> input arguments
  integer,intent(in) :: ord
  integer,dimension(3),   intent(in) :: extc
  real*8, dimension(extc(1),extc(2),extc(3)),intent(in ):: func
  real*8, dimension(-ord+1:extc(1)+ord,-ord+1:extc(2)+ord,extc(3)),intent(out):: funcc
  real*8, dimension(2), intent(in) :: SoA

  integer::i

  funcc = 0.d0
  funcc(1:extc(1),1:extc(2),1:extc(3)) = func
   do i=0,ord-1
      funcc(-i,1:extc(2),1:extc(3)) = funcc(i+1,1:extc(2),1:extc(3))*SoA(1)
      funcc(extc(1)+1+i,1:extc(2),1:extc(3)) = funcc(extc(1)-i,1:extc(2),1:extc(3))*SoA(1)
   enddo
   do i=0,ord-1
      funcc(:,-i,1:extc(3)) = funcc(:,i+1,1:extc(3))*SoA(2)
      funcc(:,extc(2)+1+i,1:extc(3)) = funcc(:,extc(2)-i,1:extc(3))*SoA(2)
   enddo

end subroutine symmetry_stbd

subroutine symmetry_sntbd(ord,extc,func,funcc,SoA,actd)
  implicit none

!~~~~~~> input arguments
  integer,intent(in) :: ord,actd
  integer,dimension(3),   intent(in) :: extc
  real*8, dimension(extc(1),extc(2),extc(3)),intent(in ):: func
  real*8, dimension(-ord+1:extc(1)+ord,-ord+1:extc(2)+ord,extc(3)),intent(out):: funcc
  real*8, intent(in) :: SoA

  integer::i

  funcc = 0.d0
  funcc(1:extc(1),1:extc(2),1:extc(3)) = func
  if(actd==0)then
   do i=0,ord-1
      funcc(-i,1:extc(2),1:extc(3)) = funcc(i+1,1:extc(2),1:extc(3))*SoA
      funcc(extc(1)+1+i,1:extc(2),1:extc(3)) = funcc(extc(1)-i,1:extc(2),1:extc(3))*SoA
   enddo
  elseif(actd==1)then
   do i=0,ord-1
      funcc(1:extc(1),-i,1:extc(3)) = funcc(1:extc(1),i+1,1:extc(3))*SoA
      funcc(1:extc(1),extc(2)+1+i,1:extc(3)) = funcc(1:extc(1),extc(2)-i,1:extc(3))*SoA
   enddo
  else
    write(*,*)"symmetry_sntbd: not recognized actd = ",actd
  endif

end subroutine symmetry_sntbd

subroutine d2dump(wei,llb,uub,ext,data_in,data_out,gord,SoA)
  implicit none
  integer,intent(in) :: wei,gord
  integer,dimension(3),intent(in) ::ext
  real*8,dimension(3),intent(in) :: SoA
  real*8,dimension(3) :: llb,uub
  real*8,dimension(ext(1),ext(2),ext(3)),intent(in)::data_in
  real*8,dimension(ext(1),ext(2)),intent(inout)::data_out

  real*8 :: dZ
  integer :: i,j,k

!sanity check
  if(wei.ne.3)then
     write(*,*)"fmisc.f90::copy: this routine only surport 3 dimension"
     write(*,*)"dim = ",wei
     stop
  endif

  dZ = (uub(3)-llb(3))/ext(3)
  k = idint((0-llb(3))/dZ+0.4)+1

  select case (gord)
  case (2)
     if(k > 2)then
         do i=1,ext(1)
         do j=1,ext(2)
             data_out(i,j) = 0.5625d0*(data_in(i,j,k)+data_in(i,j,k-1))-0.0625d0*(data_in(i,j,k+1)+data_in(i,j,k-2))
         enddo
         enddo
     else if(k == 1)then
         do i=1,ext(1)
         do j=1,ext(2)
             data_out(i,j) = 0.5625d0*(data_in(i,j,k)+SoA(3)*data_in(i,j,k))-0.0625d0*(data_in(i,j,k+1)+SoA(3)*data_in(i,j,k+1))
         enddo
         enddo
     else
      write(*,*) "d2dump: something must be wrong, k = ",k
      return
     endif
  case (3)
     if(k > 3)then
         do i=1,ext(1)
         do j=1,ext(2)
             data_out(i,j) = 0.5859375d0*(data_in(i,j,k)+data_in(i,j,k-1))   &
                           -0.9765625d-1*(data_in(i,j,k+1)+data_in(i,j,k-2)) &
                           +0.1171875d-1*(data_in(i,j,k+2)+data_in(i,j,k-3))
         enddo
         enddo
     else if(k == 1)then
         do i=1,ext(1)
         do j=1,ext(2)
             data_out(i,j) = 0.5859375d0*(data_in(i,j,k)+SoA(3)*data_in(i,j,k))   &
                           -0.9765625d-1*(data_in(i,j,k+1)+SoA(3)*data_in(i,j,k+1)) &
                           +0.1171875d-1*(data_in(i,j,k+2)+SoA(3)*data_in(i,j,k+2))
         enddo
         enddo
     else
      write(*,*) "d2dump: something must be wrong, k = ",k
      return
     endif
  case (4)
     if(k > 4)then
         do i=1,ext(1)
         do j=1,ext(2)
             data_out(i,j) = 0.5981445312d0*(data_in(i,j,k)+data_in(i,j,k-1))   &
                            -0.1196289063d0*(data_in(i,j,k+1)+data_in(i,j,k-2)) &
                           +0.2392578125d-1*(data_in(i,j,k+2)+data_in(i,j,k-3)) &
                           -0.2441406250d-2*(data_in(i,j,k+3)+data_in(i,j,k-4))
         enddo
         enddo
     else if(k == 1)then
         do i=1,ext(1)
         do j=1,ext(2)
             data_out(i,j) = 0.5981445312d0*(data_in(i,j,k)+SoA(3)*data_in(i,j,k))   &
                            -0.1196289063d0*(data_in(i,j,k+1)+SoA(3)*data_in(i,j,k+1)) &
                           +0.2392578125d-1*(data_in(i,j,k+2)+SoA(3)*data_in(i,j,k+2)) &
                           -0.2441406250d-2*(data_in(i,j,k+3)+SoA(3)*data_in(i,j,k+3))
         enddo
         enddo
     else
      write(*,*) "d2dump: something must be wrong, k = ",k
      return
     endif
  case (5)
     if(k > 5)then
         do i=1,ext(1)
         do j=1,ext(2)
             data_out(i,j) = 0.6056213378d0*(data_in(i,j,k)+data_in(i,j,k-1))   &
                            -0.1345825196d0*(data_in(i,j,k+1)+data_in(i,j,k-2)) &
                           +0.3460693359d-1*(data_in(i,j,k+2)+data_in(i,j,k-3)) &
                           -0.6179809571d-2*(data_in(i,j,k+3)+data_in(i,j,k-4)) &
                           +0.5340576171d-3*(data_in(i,j,k+4)+data_in(i,j,k-5))
         enddo
         enddo
     else if(k == 1)then
         do i=1,ext(1)
         do j=1,ext(2)
             data_out(i,j) = 0.6056213378d0*(data_in(i,j,k)+SoA(3)*data_in(i,j,k))   &
                            -0.1345825196d0*(data_in(i,j,k+1)+SoA(3)*data_in(i,j,k+1)) &
                           +0.3460693359d-1*(data_in(i,j,k+2)+SoA(3)*data_in(i,j,k+2)) &
                           -0.6179809571d-2*(data_in(i,j,k+3)+SoA(3)*data_in(i,j,k+3)) &
                           +0.5340576171d-3*(data_in(i,j,k+4)+SoA(3)*data_in(i,j,k+4))
         enddo
         enddo
     else
      write(*,*) "d2dump: something must be wrong, k = ",k
      return
     endif
  case default
    write(*,*) "d2dump: not recognized ord = ",gord
    return
  end select

end subroutine d2dump

#else
#error Not define Vertex nor Cell
#endif  
#endif
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!--------------------------------------------------------------------------------------
! calculate L2norm  
  subroutine l2normhelper(ex, X, Y, Z,xmin,ymin,zmin,xmax,ymax,zmax,&
                          f,f_out,gw)

  implicit none
!~~~~~~> Input parameters:
  integer,intent(in ):: ex(1:3)
  real*8, intent(in ):: X(1:ex(1)),Y(1:ex(2)),Z(1:ex(3)),xmin,ymin,zmin,xmax,ymax,zmax
  integer,intent(in)::gw
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: f
  real*8, intent(out) :: f_out
!~~~~~~> Other variables:

  real*8, parameter :: ZEO = 0.D0
  real*8            :: dX, dY, dZ
  integer::imin,jmin,kmin
  integer::imax,jmax,kmax
  integer::i,j,k

  dX = X(2) - X(1)
  dY = Y(2) - Y(1)
  dZ = Z(2) - Z(1)

! for ghost zone
   imin = gw+1
   jmin = gw+1
   kmin = gw+1

   imax = ex(1) - gw
   jmax = ex(2) - gw
   kmax = ex(3) - gw

!for patch boundary (i.e., not ghost boundary)

if(dabs(X(ex(1))-xmax) < dX) imax = ex(1)
if(dabs(Y(ex(2))-ymax) < dY) jmax = ex(2)
if(dabs(Z(ex(3))-zmax) < dZ) kmax = ex(3)
if(dabs(X(1)-xmin) < dX) imin = 1
if(dabs(Y(1)-ymin) < dY) jmin = 1
if(dabs(Z(1)-zmin) < dZ) kmin = 1

f_out = sum(f(imin:imax,jmin:jmax,kmin:kmax)*f(imin:imax,jmin:jmax,kmin:kmax))

f_out = f_out*dX*dY*dZ

  return

  end subroutine l2normhelper
!--------------------------------------------------------------------------------------
! calculate L2norm especially for shell Blocks
  subroutine l2normhelper_sh(ex, X, Y, Z,xmin,ymin,zmin,xmax,ymax,zmax,&
                          f,f_out,gw,ogw,Symmetry)

  implicit none
!~~~~~~> Input parameters:
  integer,intent(in ):: ex(1:3),Symmetry
  real*8, intent(in ):: X(1:ex(1)),Y(1:ex(2)),Z(1:ex(3)),xmin,ymin,zmin,xmax,ymax,zmax
  integer,intent(in)::gw,ogw
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: f
  real*8, intent(out) :: f_out
!~~~~~~> Other variables:

  real*8, parameter :: ZEO = 0.D0
  real*8            :: dX, dY, dZ
  integer::imin,jmin,kmin
  integer::imax,jmax,kmax
  integer::i,j,k

  real*8 :: PIo4

  PIo4 = dacos(-1.d0)/4.d0

  dX = X(2) - X(1)
  dY = Y(2) - Y(1)
  dZ = Z(2) - Z(1)

! for ghost zone
   imin = gw+1
   jmin = gw+1
   kmin = gw+1

   imax = ex(1) - gw
   jmax = ex(2) - gw
   kmax = ex(3) - gw

!for patch boundary (i.e., not ghost boundary)

if(dabs(X(ex(1))-xmax) < dX)then
  if(X(ex(1))-PIo4 > dX)then
    imax = ex(1)-ogw             ! for overlap zone
  else
    imax = ex(1)
  endif
endif
if(dabs(Y(ex(2))-ymax) < dY)then
  if(Y(ex(2))-PIo4 > dY)then
    jmax = ex(2)-ogw             ! for overlap zone
  else
    jmax = ex(2)
  endif
endif
if(dabs(Z(ex(3))-zmax) < dZ) kmax = ex(3)

if(dabs(X(1)-xmin) < dX)then
  if(X(1)+PIo4 < dX)then
    imin = 1+ogw             ! for overlap zone
  else
    imin = 1
  endif
endif
if(dabs(Y(1)-ymin) < dY)then
  if(Y(1)+PIo4 < dY)then
    jmin = 1+ogw             ! for overlap zone
  else
    jmin = 1
  endif
endif
if(dabs(Z(1)-zmin) < dZ) kmin = 1

!for Symmetry ghost points
if(Symmetry==1)then
  if(dabs(ymin+gw*dY)<dY.and.Y(1)<0.d0) jmin = gw+1
  if(dabs(ymax-gw*dY)<dY.and.Y(ex(2))>0.d0) jmax = ex(2)-gw
endif
if(Symmetry==2)then
  if(dabs(xmin+gw*dX)<dX.and.X(1)<0.d0) imin = gw+1
  if(dabs(ymin+gw*dY)<dY.and.Y(1)<0.d0) jmin = gw+1
endif

f_out = sum(f(imin:imax,jmin:jmax,kmin:kmax)*f(imin:imax,jmin:jmax,kmin:kmax))

f_out = f_out*dX*dY*dZ

  return

  end subroutine l2normhelper_sh
!--------------------------------------------------------------------------------------
! calculate L2norm especially for shell Blocks
! use root mean sqrt method
  subroutine l2normhelper_sh_rms(ex, X, Y, Z,xmin,ymin,zmin,xmax,ymax,zmax,&
                          f,f_out,gw,ogw,Symmetry,Nout)

  implicit none
!~~~~~~> Input parameters:
  integer,intent(in ):: ex(1:3),Symmetry
  real*8, intent(in ):: X(1:ex(1)),Y(1:ex(2)),Z(1:ex(3)),xmin,ymin,zmin,xmax,ymax,zmax
  integer,intent(in)::gw,ogw
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: f
  real*8, intent(out) :: f_out
  integer,intent(out) :: Nout
!~~~~~~> Other variables:

  real*8, parameter :: ZEO = 0.D0
  real*8            :: dX, dY, dZ
  integer::imin,jmin,kmin
  integer::imax,jmax,kmax
  integer::i,j,k

  real*8 :: PIo4

  PIo4 = dacos(-1.d0)/4.d0

  dX = X(2) - X(1)
  dY = Y(2) - Y(1)
  dZ = Z(2) - Z(1)

! for ghost zone
   imin = gw+1
   jmin = gw+1
   kmin = gw+1

   imax = ex(1) - gw
   jmax = ex(2) - gw
   kmax = ex(3) - gw

!for patch boundary (i.e., not ghost boundary)

if(dabs(X(ex(1))-xmax) < dX)then
  if(X(ex(1))-PIo4 > dX)then
    imax = ex(1)-ogw             ! for overlap zone
  else
    imax = ex(1)
  endif
endif
if(dabs(Y(ex(2))-ymax) < dY)then
  if(Y(ex(2))-PIo4 > dY)then
    jmax = ex(2)-ogw             ! for overlap zone
  else
    jmax = ex(2)
  endif
endif
if(dabs(Z(ex(3))-zmax) < dZ) kmax = ex(3)

if(dabs(X(1)-xmin) < dX)then
  if(X(1)+PIo4 < dX)then
    imin = 1+ogw             ! for overlap zone
  else
    imin = 1
  endif
endif
if(dabs(Y(1)-ymin) < dY)then
  if(Y(1)+PIo4 < dY)then
    jmin = 1+ogw             ! for overlap zone
  else
    jmin = 1
  endif
endif
if(dabs(Z(1)-zmin) < dZ) kmin = 1

!for Symmetry ghost points
if(Symmetry==1)then
  if(dabs(ymin+gw*dY)<dY.and.Y(1)<0.d0) jmin = gw+1
  if(dabs(ymax-gw*dY)<dY.and.Y(ex(2))>0.d0) jmax = ex(2)-gw
endif
if(Symmetry==2)then
  if(dabs(xmin+gw*dX)<dX.and.X(1)<0.d0) imin = gw+1
  if(dabs(ymin+gw*dY)<dY.and.Y(1)<0.d0) jmin = gw+1
endif

f_out = sum(f(imin:imax,jmin:jmax,kmin:kmax)*f(imin:imax,jmin:jmax,kmin:kmax))

f_out = f_out

Nout = (imax-imin+1)*(jmax-jmin+1)*(kmax-kmin+1)

  return

  end subroutine l2normhelper_sh_rms
! locating the position of NaN
  subroutine ScaneNaN(ext,X,Y,Z,f)
  implicit none
  integer,dimension(3),intent(in) :: ext
  real*8,dimension(ext(1)) :: X
  real*8,dimension(ext(2)) :: Y
  real*8,dimension(ext(3)) :: Z
  real*8,dimension(ext(1),ext(2),ext(3)) :: f

  integer :: i,j,k

  do k=1,ext(3)
  do j=1,ext(2)
  do i=1,ext(1)
     if(abs(f(i,j,k)) .ne. abs(f(i,j,k))) write(*,*)X(i),Y(j),Z(k),f(i,j,k)
  enddo
  enddo
  enddo

  end subroutine ScaneNaN
! fortran version writefile
  subroutine fwritefile(time,nx,ny,nz,xmin,xmax,ymin,ymax,zmin,zmax,NN,filename,data_out)
  implicit none

  real*8,intent(in) :: time,xmin,xmax,ymin,ymax,zmin,zmax
  integer,intent(in) :: nx,ny,nz,NN
  real*8,dimension(nx,ny,nz),intent(in) :: data_out
  Character(Len=NN) :: filename

 
!  Open( 12 , File = filename,form='BINARY', access="SEQUENTIAL",status="replace",action='Write')
  integer :: recl_val
  INQUIRE(IOLENGTH=recl_val) time,nx,ny,nz,xmin,xmax,ymin,ymax,zmin,zmax,data_out
  Open( 12 , File = filename,form='UNFORMATTED', access="DIRECT",status="replace",action='Write',RECL=recl_val)
  Write( 12 ) time,nx,ny,nz,xmin,xmax,ymin,ymax,zmin,zmax,data_out

  Close( 12 )

  end subroutine fwritefile
!--------------------------------------------------------------------------
!
! average for interface 
!
!--------------------------------------------------------------------------
  subroutine average(ext,f1,f2,fout)
  use sft_trace_mod
  use omp_lib
  implicit none
  integer,dimension(3),   intent(in) :: ext
  real*8, dimension(ext(1),ext(2),ext(3)),intent(in):: f1,f2
  real*8, dimension(ext(1),ext(2),ext(3)),intent(out):: fout

  real*8,parameter::HLF=0.5d0
  integer :: i,j,k
  integer(8) :: ts0_thr
  integer :: tid

!$OMP PARALLEL PRIVATE(i,j,k,ts0_thr,tid) SHARED(f1,f2,fout,ext)
  tid = omp_get_thread_num()
  ts0_thr = sft_get_ts()
!$OMP DO COLLAPSE(1)
  do k=1,ext(3)
    do j=1,ext(2)
      do i=1,ext(1)
        fout(i,j,k) = HLF*(f1(i,j,k)+f2(i,j,k))
      enddo
    enddo
  enddo
!$OMP END DO
  call sft_trace("average", ts0_thr, sft_get_ts(), 0, tid)
!$OMP END PARALLEL

  return

  end subroutine average
!-----------------------------------------------------------------------------  
  subroutine average3(ext,f1,f2,fout)
  use sft_trace_mod
  use omp_lib
  implicit none
  integer,dimension(3),   intent(in) :: ext
  real*8, dimension(ext(1),ext(2),ext(3)),intent(in):: f1,f2
  real*8, dimension(ext(1),ext(2),ext(3)),intent(out):: fout
! f1 ----------                ^
!                fout ------p  | t
! f2 ----------                |
! 2 points, 1st order interpolation
! 1   2
! f2  f1
! *---*--> t
!    ^
! f=3/4*f_1 + 1/4*f_2

  real*8,parameter::C1=0.75d0,C2=0.25d0
  integer :: i,j,k
  integer(8) :: ts0_thr
  integer :: tid

!$OMP PARALLEL PRIVATE(i,j,k,ts0_thr,tid) SHARED(f1,f2,fout,ext)
  tid = omp_get_thread_num()
  ts0_thr = sft_get_ts()
!$OMP DO COLLAPSE(1)
  do k=1,ext(3)
    do j=1,ext(2)
      do i=1,ext(1)
        fout(i,j,k) = C1*f1(i,j,k)+C2*f2(i,j,k)
      enddo
    enddo
  enddo
!$OMP END DO
  call sft_trace("average3", ts0_thr, sft_get_ts(), 0, tid)
!$OMP END PARALLEL

  return

  end subroutine average3
!-----------------------------------------------------------------------------  
  subroutine average2(ext,f1,f2,f3,fout)
  use sft_trace_mod
  use omp_lib
  implicit none
  integer,dimension(3),   intent(in) :: ext
  real*8, dimension(ext(1),ext(2),ext(3)),intent(in):: f1,f2,f3
  real*8, dimension(ext(1),ext(2),ext(3)),intent(out):: fout
! f1 ----------                ^
!                fout ------   |
! f2 ----------                | t
!                              |
! f3 ----------                |
! 3 points, 2nd order interpolation
! 1   2   3
! f3  f2  f1
! *---*---*--> t
!       ^
! f=3/8*f_1 + 3/4*f_2 - 1/8*f_3

  real*8,parameter::C1=3.d0/8.d0,C2=3.d0/4.d0,C3=-1.d0/8.d0
  integer :: i,j,k
  integer(8) :: ts0_thr
  integer :: tid

!$OMP PARALLEL PRIVATE(i,j,k,ts0_thr,tid) SHARED(f1,f2,f3,fout,ext)
  tid = omp_get_thread_num()
  ts0_thr = sft_get_ts()
!$OMP DO COLLAPSE(1)
  do k=1,ext(3)
    do j=1,ext(2)
      do i=1,ext(1)
        fout(i,j,k) = C1*f1(i,j,k)+C2*f2(i,j,k)+C3*f3(i,j,k)
      enddo
    enddo
  enddo
!$OMP END DO
  call sft_trace("average2", ts0_thr, sft_get_ts(), 0, tid)
!$OMP END PARALLEL

  return

  end subroutine average2
!-----------------------------------------------------------------------------  
  subroutine average2p(ext,f1,f2,f3,fout)
  use sft_trace_mod
  use omp_lib
  implicit none
  integer,dimension(3),   intent(in) :: ext
  real*8, dimension(ext(1),ext(2),ext(3)),intent(in):: f1,f2,f3
  real*8, dimension(ext(1),ext(2),ext(3)),intent(out):: fout
! f1 ----------                ^
!                fout ------p  |
! f2 ----------                | t
!                              |
! f3 ----------                |
! 3 points, 2nd order interpolation
! 1   2   3
! f3  f2  f1
! *---*---*--> t
!        ^
! f=21/32*f_1 + 7/16*f_2 - 3/32*f_3

  real*8,parameter::C1=5.d0/3.2d1,C2=1.5d1/1.6d1,C3=-3.d0/3.2d1
  integer :: i,j,k
  integer(8) :: ts0_thr
  integer :: tid

!$OMP PARALLEL PRIVATE(i,j,k,ts0_thr,tid) SHARED(f1,f2,f3,fout,ext)
  tid = omp_get_thread_num()
  ts0_thr = sft_get_ts()
!$OMP DO COLLAPSE(1)
  do k=1,ext(3)
    do j=1,ext(2)
      do i=1,ext(1)
        fout(i,j,k) = C1*f1(i,j,k)+C2*f2(i,j,k)+C3*f3(i,j,k)
      enddo
    enddo
  enddo
!$OMP END DO
  call sft_trace("average2p", ts0_thr, sft_get_ts(), 0, tid)
!$OMP END PARALLEL

  return

  end subroutine average2p
!-----------------------------------------------------------------------------  
  subroutine average2m(ext,f1,f2,f3,fout)
  use sft_trace_mod
  use omp_lib
  implicit none
  integer,dimension(3),   intent(in) :: ext
  real*8, dimension(ext(1),ext(2),ext(3)),intent(in):: f1,f2,f3
  real*8, dimension(ext(1),ext(2),ext(3)),intent(out):: fout
! f1 ----------                ^
!                fout ------m  |
! f2 ----------                | t
!                              |
! f3 ----------                |
! 3 points, 2nd order interpolation
! 1   2   3
! f3  f2  f1
! *---*---*--> t
!      ^
! f=5/32*f_1 + 15/16*f_2 - 3/32*f_3

  real*8,parameter::C1=5.d0/3.2d1,C2=1.5d1/1.6d1,C3=-3.d0/3.2d1
  integer :: i,j,k
  integer(8) :: ts0_thr
  integer :: tid

!$OMP PARALLEL PRIVATE(i,j,k,ts0_thr,tid) SHARED(f1,f2,f3,fout,ext)
  tid = omp_get_thread_num()
  ts0_thr = sft_get_ts()
!$OMP DO COLLAPSE(1)
  do k=1,ext(3)
    do j=1,ext(2)
      do i=1,ext(1)
        fout(i,j,k) = C1*f1(i,j,k)+C2*f2(i,j,k)+C3*f3(i,j,k)
      enddo
    enddo
  enddo
!$OMP END DO
  call sft_trace("average2m", ts0_thr, sft_get_ts(), 0, tid)
!$OMP END PARALLEL

  return

  end subroutine average2m
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine lowerboundset(ex,chi0,TINNY)
  implicit none

!~~~~~~% Input parameters:

  integer ,intent(in):: ex(1:3)
  real*8  ,intent(in):: TINNY
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) ::chi0

  where(chi0 < TINNY) chi0 = TINNY

  return   

  end subroutine lowerboundset
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!global interpolation with given index and coeffients
!------------------------------------------
!fortran version of Wigner d function
!Eq.(42) of PRD 77, 024027 (2008)
!we consider only theta in [0,pi]
!------------------------------------------
  function fWigner_d_function(l,m,s,costheta) result(gont)
  implicit none
  integer,intent(in) :: l,m,s
  real*8,intent(in) :: costheta

  real*8 :: gont

  integer :: t,C1,C2
  real*8 :: ffact,vv,sinht,cosht

  C1=max(0,m-s)
  C2=min(l+m,l-s)
  vv=0
  sinht=dsqrt((1.d0-costheta)/2.d0)
  cosht=dsqrt((1.d0+costheta)/2.d0);
  if(C1/2*2==C1)then
    do t=C1,C2,2
       vv=vv+cosht**(2*l+m-s-2*t)*sinht**(2*t+s-m)/(ffact(l+m-t)*ffact(l-s-t)*ffact(t)*ffact(t+s-m))
    enddo
    do t=C1+1,C2,2
       vv=vv-cosht**(2*l+m-s-2*t)*sinht**(2*t+s-m)/(ffact(l+m-t)*ffact(l-s-t)*ffact(t)*ffact(t+s-m))
    enddo
  else
    do t=C1,C2,2
       vv=vv-cosht**(2*l+m-s-2*t)*sinht**(2*t+s-m)/(ffact(l+m-t)*ffact(l-s-t)*ffact(t)*ffact(t+s-m))
    enddo
    do t=C1+1,C2,2
       vv=vv+cosht**(2*l+m-s-2*t)*sinht**(2*t+s-m)/(ffact(l+m-t)*ffact(l-s-t)*ffact(t)*ffact(t+s-m))
    enddo
  endif
  
  gont = vv*dsqrt(ffact(l+m)*ffact(l-m)*ffact(l+s)*ffact(l-s))

  return

  end function fWigner_d_function
!----------------------------------
  function ffact(N) result(gont)
  implicit none
  integer,intent(in) :: N

  real*8 :: gont

  integer :: i

! sanity check
  if(N < 0)then
     write(*,*) "ffact: error input for factorial"
     return
  endif

  gont = 1.d0
  do i=1,N
     gont = gont*i
  enddo

  return

  end function ffact
!---------------------------  
!Eq.(41) of PRD 77, 024027 (2008)
!----------------------------------
  function Yslm(s,l,m,the,phi) result(gont)
  implicit none
  integer,intent(in) :: s,l,m
  real*8,intent(in) :: the,phi

  double complex :: gont

  real*8 :: fWigner_d_function,PI,rp

  PI = dacos(-1.d0)

  rp = fWigner_d_function(l,m,s,dcos(the))
  rp = rp*dsqrt((2*l+1.d0)/4.d0/PI)
  if(s/2*2.ne.s) rp = -rp

  gont = dcmplx(dcos(m*phi),dsin(m*phi))

  gont = rp*gont

  return

  end function Yslm
!------------------------------------------------------------------------------------  
subroutine set_value(ext,data_out,rr)

  IMPLICIT NONE

  integer, intent(in) :: ext(3)
  REAL*8, DIMENSION(ext(1),ext(2),ext(3)), intent(out) :: data_out
  REAL*8, intent(in) :: rr

  data_out = rr

  return
end subroutine set_value
subroutine add_value(ext,data_out,rr)

  IMPLICIT NONE

  integer, intent(in) :: ext(3)
  REAL*8, DIMENSION(ext(1),ext(2),ext(3)), intent(inout) :: data_out
  REAL*8, intent(in) :: rr

  data_out = data_out + rr

  return
end subroutine add_value
! copy array2 to array1  
subroutine array_copy(ext,data1,data2)

  IMPLICIT NONE

  integer, intent(in) :: ext(3)
  REAL*8, DIMENSION(ext(1),ext(2),ext(3)), intent(out) :: data1
  REAL*8, DIMENSION(ext(1),ext(2),ext(3)), intent(in) :: data2

  data1 = data2

  return
  end subroutine array_copy
! add array2 to array1  
subroutine array_add(ext,data1,data2)

  IMPLICIT NONE

  integer, intent(in) :: ext(3)
  REAL*8, DIMENSION(ext(1),ext(2),ext(3)), intent(inout) :: data1
  REAL*8, DIMENSION(ext(1),ext(2),ext(3)), intent(in) :: data2

  data1 = data1 + data2

  return
  end subroutine array_add
! subtract array2 from array1  
subroutine array_subtract(ext,data1,data2)

  IMPLICIT NONE

  integer, intent(in) :: ext(3)
  REAL*8, DIMENSION(ext(1),ext(2),ext(3)), intent(inout) :: data1
  REAL*8, DIMENSION(ext(1),ext(2),ext(3)), intent(in) :: data2

  data1 = data1 - data2

  return
  end subroutine array_subtract
! find out the maximum
subroutine find_maximum(ext,X,Y,Z,fun,val,pos,llb,uub)

  implicit none

  integer,intent(in) :: ext(3),llb(3),uub(3)
  real*8 :: X(ext(1)),Y(ext(2)),Z(ext(3))
  REAL*8, DIMENSION(ext(1),ext(2),ext(3)), intent(in) :: fun
  real*8,intent(out) :: val,pos(3)

  integer :: i,j,k,ii,jj,kk
  real*8 :: tmp

  tmp = 0.d0

  ii=1
  jj=1
  kk=1

  do k=llb(3)+1,ext(3)-uub(3)
  do j=llb(2)+1,ext(2)-uub(2)
  do i=llb(1)+1,ext(1)-uub(1)
      if(dabs(fun(i,j,k)) > tmp)then
         tmp = dabs(fun(i,j,k))
         ii = i
         jj = j
         kk = k
      endif
  enddo
  enddo
  enddo

  pos(1) = X(ii)
  pos(2) = Y(jj)
  pos(3) = Z(kk)
  val = tmp

  return

end subroutine

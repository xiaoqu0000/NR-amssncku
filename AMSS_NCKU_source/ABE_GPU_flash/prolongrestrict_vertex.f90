

! Because of overlap determination, source region is always larger than target
! region

#include "macrodef.fh"

#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif   

!--------------------------------------------------------------------------
!
! Prepare the data on coarse level for prolong
! valid for all finite difference order
!--------------------------------------------------------------------------

  subroutine prolongcopy3(wei,llbc,uubc,extc,func,&
                      llbf,uubf,exto,funo,&
                      llbp,uubp,SoA,Symmetry)
  implicit none

!~~~~~~> input arguments
  integer,intent(in) :: wei
!                                       coarse    fine       coarse
  real*8,dimension(3),   intent(in) :: llbc,uubc,llbf,uubf,llbp,uubp
  integer,dimension(3),  intent(in) :: extc,exto
  real*8, dimension(extc(1),extc(2),extc(3)),intent(in)   :: func
! both bounds ghost_width
  real*8, dimension(exto(1)+2*ghost_width,exto(2)+2*ghost_width,exto(3)+2*ghost_width),intent(out):: funo
  real*8, dimension(1:3), intent(in) :: SoA
  integer,intent(in)::Symmetry

!~~~~~~> local variables

  real*8,dimension(1-ghost_width:extc(1),1-ghost_width:extc(2),1-ghost_width:extc(3)) :: fh
  real*8, dimension(1:3) :: base
  integer,dimension(3) :: lbc,ubc,lbf,ubf,lbp,ubp,lbpc,ubpc,cxI
  integer :: i,j,k

  integer::imini,imaxi,jmini,jmaxi,kmini,kmaxi
  integer::imino,imaxo,jmino,jmaxo,kmino,kmaxo

  real*8,dimension(3) :: CD,FD

  if(wei.ne.3)then
     write(*,*)"prolongrestrict.f90::prolongcopy3: this routine only surport 3 dimension"
     write(*,*)"dim = ",wei
     stop
  endif

! it's possible a iolated point for target but not for source
  CD = (uubc-llbc)/(extc-1)
  FD = CD/2

!take care the mismatch of the two segments of grid
  do i=1,3
     if(llbc(i) <= llbf(i))then
        base(i) = llbc(i)
     else
        j=idint((llbc(i)-llbf(i))/FD(i)+0.4)
        if(j/2*2 == j)then
           base(i) = llbf(i)
        else
           base(i) = llbf(i) - CD(i)/2
        endif
     endif
  enddo

!!! function idint:
!If A is of type REAL and |A| < 1, INT(A) equals 0. If |A| \geq 1, 
!then INT(A) equals the largest integer that does not exceed the range of A 
!and whose sign is the same as the sign of A.

    lbf = idint((llbf-base)/FD+0.4)+1
    ubf = idint((uubf-base)/FD+0.4)+1
    lbc = idint((llbc-base)/CD+0.4)+1
    ubc = idint((uubc-base)/CD+0.4)+1
    lbp = idint((llbp-base)/FD+0.4)+1
    lbpc = idint((llbp-base)/CD+0.4)+1
    ubp = idint((uubp-base)/FD+0.4)+1
    ubpc = idint((uubp-base)/CD+0.4)+1
!sanity check
!|*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*| 
!|x===============x===============x===============x========|
!                 ^                               ^
  imini=lbpc(1)-lbc(1) + 1 - ghost_width
  imaxi=ubpc(1)-lbc(1) + 1 + ghost_width
  jmini=lbpc(2)-lbc(2) + 1 - ghost_width
  jmaxi=ubpc(2)-lbc(2) + 1 + ghost_width
  kmini=lbpc(3)-lbc(3) + 1 - ghost_width
  kmaxi=ubpc(3)-lbc(3) + 1 + ghost_width

  cxI(1) = imaxi-imini+1
  cxI(2) = jmaxi-jmini+1
  cxI(3) = kmaxi-kmini+1
  if(any(cxI.ne.exto+2*ghost_width).or. &
     imaxi.gt.extc(1)+1.or.jmaxi.gt.extc(2)+1.or.kmaxi.gt.extc(3)+1)then
          write(*,*)"error in prolongationcopy3 for"
          if(any(cxI.ne.exto+2*ghost_width))then
             write(*,*) cxI,exto+2*ghost_width
             return
          endif
          write(*,*)"from"
          write(*,*)llbc,uubc
          write(*,*)lbc,ubc
          write(*,*)"to"
          write(*,*)llbf,uubf
          write(*,*)lbf,ubf
          write(*,*)"want"
          write(*,*)llbp,uubp
          write(*,*)lbp,ubp,lbpc,ubpc
          if(imini.lt.1) write(*,*)"imini = ",imini
          if(jmini.lt.1) write(*,*)"jmini = ",jmini       
          if(kmini.lt.1) write(*,*)"kmini = ",kmini     
          if(imaxi.gt.extc(1)) write(*,*)"imaxi = ",imaxi,"extc(1) = ",extc(1) 
          if(jmaxi.gt.extc(2)) write(*,*)"jmaxi = ",jmaxi,"extc(2) = ",extc(2) 
          if(kmaxi.gt.extc(3)) write(*,*)"kmaxi = ",kmaxi,"extc(3) = ",extc(3) 
          return
  endif

! because some point needs 2*ghost_width
! while   some point needs 2*ghost_width-1
! so we use 0 to fill empty points
  if(imini < 1.or.jmini < 1.or.kmini < 1)then
    if(imini<1.and.dabs(llbp(1))>CD(1)) write(*,*)"prolongcopy3 warning: ",llbp(1)
    if(jmini<1.and.dabs(llbp(2))>CD(2)) write(*,*)"prolongcopy3 warning: ",llbp(2)
    if(kmini<1.and.dabs(llbp(3))>CD(3)) write(*,*)"prolongcopy3 warning: ",llbp(3)
    call symmetry_bd(ghost_width,extc,func,fh,SoA)
    if(imaxi<=extc(1).and.jmaxi<=extc(2).and.kmaxi<=extc(3))then
      funo = fh(imini:imaxi,jmini:jmaxi,kmini:kmaxi)
    else
      funo = 0.d0
      cxI = 0
      if(imaxi>extc(1))then
        cxI(1) = 1
        imaxi = extc(1)
      endif
      if(jmaxi>extc(2))then
        cxI(2) = 1
        jmaxi = extc(2)
      endif
      if(kmaxi>extc(3))then
        cxI(3) = 1
        kmaxi = extc(3)
      endif
      funo(1:exto(1)+2*ghost_width-cxI(1), &
           1:exto(2)+2*ghost_width-cxI(2), &
           1:exto(3)+2*ghost_width-cxI(3)) = fh(imini:imaxi,jmini:jmaxi,kmini:kmaxi)
    endif
  else
    if(imaxi<=extc(1).and.jmaxi<=extc(2).and.kmaxi<=extc(3))then
      funo = func(imini:imaxi,jmini:jmaxi,kmini:kmaxi)
    else
      funo = 0.d0
      cxI = 0
      if(imaxi>extc(1))then
        cxI(1) = 1
        imaxi = extc(1)
      endif
      if(jmaxi>extc(2))then
        cxI(2) = 1
        jmaxi = extc(2)
      endif
      if(kmaxi>extc(3))then
        cxI(3) = 1
        kmaxi = extc(3)
      endif
      funo(1:exto(1)+2*ghost_width-cxI(1), &
           1:exto(2)+2*ghost_width-cxI(2), &
           1:exto(3)+2*ghost_width-cxI(3)) = func(imini:imaxi,jmini:jmaxi,kmini:kmaxi)
    endif
  endif

  return

  end subroutine prolongcopy3
!=================================================================================================
!--------------------------------------------------------------------------
!
! Prolong data throug mix data of fine and coarse levels
!--------------------------------------------------------------------------

  subroutine prolongmix3(wei,llbf,uubf,extf,funf,&
                      llbc,uubc,exti,funi,&
                      llbp,uubp,SoA,Symmetry, &
                      illb,iuub)
  implicit none

!~~~~~~> input arguments
  integer,intent(in) :: wei
!                                       coarse      fine     coarse   fine (real inner points)
  real*8,dimension(3),   intent(in) :: llbc,uubc,llbf,uubf,llbp,uubp,illb,iuub
  integer,dimension(3),  intent(in) :: exti,extf
  real*8, dimension(extf(1),extf(2),extf(3)),intent(inout)   :: funf
! lower bound ghost_width; upper bound ghost_width-1  
  real*8, dimension(exti(1)+2*ghost_width,exti(2)+2*ghost_width,exti(3)+2*ghost_width),intent(in):: funi
  real*8, dimension(1:3), intent(in) :: SoA
  integer,intent(in)::Symmetry

!~~~~~~> local variables

  real*8,dimension(1-ghost_width:extf(1),1-ghost_width:extf(2),1-ghost_width:extf(3)) :: fh
  real*8, dimension(1:3) :: base
  integer,dimension(3) :: lbc,ubc,lbf,ubf,lbp,ubp,lbpc,ubpc,ilb,iub
  integer :: i,j,k,n,ii,jj,kk

  integer::imino,imaxo,jmino,jmaxo,kmino,kmaxo

  real*8,dimension(3) :: CD,FD
  integer,dimension(3) :: cxI,cxB,cxT,fg

  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2

  real*8,dimension(2*ghost_width,2*ghost_width,2*ghost_width) :: ya
  real*8,dimension(2*ghost_width) :: X,Y,Z
  real*8, dimension(2*ghost_width,2*ghost_width) :: tmp2
  real*8, dimension(2*ghost_width) :: tmp1
  real*8 :: ddy

#if (ghost_width == 2)
  real*8, parameter :: C1=-1.d0/16,C2=9.d0/16
#elif (ghost_width == 3)
  real*8, parameter :: C1=3.d0/2.56d2,C2=-2.5d1/2.56d2,C3=7.5d1/1.28d2
#elif (ghost_width == 4)
  real*8, parameter :: C1=-5.d0/2.048d3,C2=4.9d1/2.048d3,C3=-2.45d2/2.048d3,C4=1.225d3/2.048d3
#elif (ghost_width == 5)
  real*8, parameter :: C1=3.5d1/6.5536d4,C2=-4.05d2/6.5536d4,C3=5.67d2/1.6384d4
  real*8, parameter :: C4=-2.205d3/1.6384d4,C5=1.9845d4/3.2768d4
#endif  

  if(wei.ne.3)then
     write(*,*)"prolongrestrict.f90::prolongmix3: this routine only surport 3 dimension"
     write(*,*)"dim = ",wei
     stop
  endif

! it's possible a iolated point for target but not for source
  FD = (uubf-llbf)/(extf-1)
  CD = FD*2.d0

!take care the mismatch of the two segments of grid
  do i=1,3
     if(llbc(i) <= llbf(i))then
        base(i) = llbc(i)
     else
        j=idint((llbc(i)-llbf(i))/FD(i)+0.4)
        if(j/2*2 == j)then
           base(i) = llbf(i)
        else
           base(i) = llbf(i) - CD(i)/2
        endif
     endif
  enddo

!!! function idint:
!If A is of type REAL and |A| < 1, INT(A) equals 0. If |A| \geq 1, 
!then INT(A) equals the largest integer that does not exceed the range of A 
!and whose sign is the same as the sign of A.

    lbf = idint((llbf-base)/FD+0.4)+1
    ubf = idint((uubf-base)/FD+0.4)+1
    lbc = idint((llbc-base)/CD+0.4)+1
    ubc = idint((uubc-base)/CD+0.4)+1
    lbp = idint((llbp-base)/FD+0.4)+1
    lbpc = idint((llbp-base)/CD+0.4)+1
    ubp = idint((uubp-base)/FD+0.4)+1
    ubpc = idint((uubp-base)/CD+0.4)+1
    ilb = idint((illb-base)/FD+0.4)+1
    iub = idint((iuub-base)/FD+0.4)+1
!sanity check
  imino=lbp(1)-lbf(1) + 1
  imaxo=ubp(1)-lbf(1) + 1
  jmino=lbp(2)-lbf(2) + 1
  jmaxo=ubp(2)-lbf(2) + 1
  kmino=lbp(3)-lbf(3) + 1
  kmaxo=ubp(3)-lbf(3) + 1

!sanity check
!|*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*| 
!|x===============x===============x===============x===============x===============x|
!                 ^                                               ^
! ghost_width for both sides
  lbpc = lbpc - ghost_width
  ubpc = ubpc + ghost_width
! index for real inner points  
  ilb = ilb - lbf+1
  iub = iub - lbf+1

! because of domain division by parallelization
  ilb = max(ilb,1)
  iub = min(iub,extf)

  if(imino.lt.1.or.jmino.lt.1.or.kmino.lt.1.or.&
     imaxo.gt.extf(1).or.jmaxo.gt.extf(2).or.kmaxo.gt.extf(3))then
          write(*,*)"error in prolongmix3 for"
          write(*,*)"from"
          write(*,*)llbc,uubc
          write(*,*)lbc,ubc
          write(*,*)"to"
          write(*,*)llbf,uubf
          write(*,*)lbf,ubf
          write(*,*)base,FD
          write(*,*)"want"
          write(*,*)llbp,uubp
          write(*,*)lbp,ubp
          if(imino.lt.1) write(*,*)"imino = ",imino
          if(jmino.lt.1) write(*,*)"jmino = ",jmino       
          if(kmino.lt.1) write(*,*)"kmino = ",kmino
          if(imaxo.gt.extf(1)) write(*,*)"imaxo = ",imaxo,"extf(1) = ",extf(1) 
          if(jmaxo.gt.extf(2)) write(*,*)"jmaxo = ",jmaxo,"extf(2) = ",extf(2) 
          if(kmaxo.gt.extf(3)) write(*,*)"kmaxo = ",kmaxo,"extf(3) = ",extf(3) 
          return
  endif

  if(Symmetry > NO_SYMM .and. dabs(illb(3)) < FD(3)) ilb(3) = 1-ghost_width
  if(Symmetry > EQ_SYMM .and. dabs(illb(1)) < FD(1)) ilb(1) = 1-ghost_width
  if(Symmetry > EQ_SYMM .and. dabs(illb(2)) < FD(2)) ilb(2) = 1-ghost_width

  if(any(ilb<1))then
    call symmetry_bd(ghost_width,extf,funf,fh,SoA)
  else
    fh(1:extf(1),1:extf(2),1:extf(3)) = funf
  endif

  do k=kmino,kmaxo
  do j=jmino,jmaxo
  do i=imino,imaxo
       cxI(1) = i
       cxI(2) = j
       cxI(3) = k

! for fine level we use cxI-ghost_width,....cxI,....cxI+ghost_width-1       
       cxB = max(cxI-ghost_width  ,ilb)
       cxT = min(cxI+ghost_width-1,iub)
! change to coarse level reference
!|*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*| 
!|x===============x===============x===============x===============x===============x|
       cxI = (cxI+lbf)/2
! change to array index      
       cxI = cxI - lbpc + 1

       ya = funi(cxI(1)-ghost_width+1:cxI(1)+ghost_width,cxI(2)-ghost_width+1:cxI(2)+ghost_width,cxI(3)-ghost_width+1:cxI(3)+ghost_width)

       fg = 0
       if(cxT(1)>=i.and.cxB(1)<=i) fg(1) = 1
       if(cxT(2)>=j.and.cxB(2)<=j) fg(2) = 1
       if(cxT(3)>=k.and.cxB(3)<=k) fg(3) = 1

       if(cxT(1)>=cxB(1) .and. cxT(2)>=cxB(2) .and. cxT(3)>=cxB(3).and. sum(fg).eq.2)then
          if(any(cxB<1-ghost_width).or.any(cxT>extf))then
             write(*,*) "error in prolongmix3: "
             if(any(cxB<1-ghost_width)) write(*,*) cxB,1-ghost_width
             if(any(cxT>extf)         ) write(*,*) cxT,extf,iuub,uubf
             stop
          endif

! fix the wanted point at (0,0,0), set FD = 1       
       ii=i+lbf(1)-1
       jj=j+lbf(2)-1
       kk=k+lbf(3)-1

       if(ii/2*2==ii)then
         do n=1,ghost_width
           X(ghost_width-n+1) = -1.d0-(n-1)*2
           X(ghost_width+n  ) =  1.d0+(n-1)*2
         enddo
       else
         do n=1,ghost_width
           X(ghost_width-n+1) = -(n-1)*2.d0
           X(ghost_width+n  ) =   n   *2.d0
         enddo
       endif

       if(jj/2*2==jj)then
         do n=1,ghost_width
           Y(ghost_width-n+1) = -1.d0-(n-1)*2
           Y(ghost_width+n  ) =  1.d0+(n-1)*2
         enddo
       else
         do n=1,ghost_width
           Y(ghost_width-n+1) = -(n-1)*2.d0
           Y(ghost_width+n  ) =   n   *2.d0
         enddo
       endif

       if(kk/2*2==kk)then
         do n=1,ghost_width
           Z(ghost_width-n+1) = -1.d0-(n-1)*2
           Z(ghost_width+n  ) =  1.d0+(n-1)*2
         enddo
       else
         do n=1,ghost_width
           Z(ghost_width-n+1) = -(n-1)*2.d0
           Z(ghost_width+n  ) =   n   *2.d0
         enddo
       endif

! i=>(ghost_width,0), i-ghost_width=>(1,1-ghost_width)
       do n=cxB(1)+ghost_width-i+1,cxT(1)+ghost_width-i+1
          X(n) = n-ghost_width 
       enddo

       do n=cxB(2)+ghost_width-j+1,cxT(2)+ghost_width-j+1
          Y(n) = n-ghost_width 
       enddo

       do n=cxB(3)+ghost_width-k+1,cxT(3)+ghost_width-k+1
          Z(n) = n-ghost_width 
       enddo

! because of the mismatch of points for fine level and coarse level 
! we have to deal in this way

! for x direction
       if(fg(1) .eq. 0)then

#if (ghost_width == 2)
         if(kk/2*2==kk)then
             tmp2= C1*(ya(:,:,1)+ya(:,:,4))+C2*(ya(:,:,2)+ya(:,:,3))
         else
             tmp2= ya(:,:,2)
         endif
         if(jj/2*2==jj)then
             tmp1= C1*(tmp2(:,1)+tmp2(:,4))+C2*(tmp2(:,2)+tmp2(:,3))
         else
             tmp1= tmp2(:,2)
         endif
#elif (ghost_width == 3)
         if(kk/2*2==kk)then
             tmp2= C1*(ya(:,:,1)+ya(:,:,6))+C2*(ya(:,:,2)+ya(:,:,5))+C3*(ya(:,:,3)+ya(:,:,4))
         else
             tmp2= ya(:,:,3)
         endif
         if(jj/2*2==jj)then
             tmp1= C1*(tmp2(:,1)+tmp2(:,6))+C2*(tmp2(:,2)+tmp2(:,5))+C3*(tmp2(:,3)+tmp2(:,4))
         else
             tmp1= tmp2(:,3)
         endif
#elif (ghost_width == 4)
         if(kk/2*2==kk)then
             tmp2= C1*(ya(:,:,1)+ya(:,:,8))+C2*(ya(:,:,2)+ya(:,:,7)) &
                  +C3*(ya(:,:,3)+ya(:,:,6))+C4*(ya(:,:,4)+ya(:,:,5))
         else
             tmp2= ya(:,:,4)
         endif
         if(jj/2*2==jj)then
             tmp1= C1*(tmp2(:,1)+tmp2(:,8))+C2*(tmp2(:,2)+tmp2(:,7)) &
                  +C3*(tmp2(:,3)+tmp2(:,6))+C4*(tmp2(:,4)+tmp2(:,5))
         else
             tmp1= tmp2(:,4)
         endif
#elif (ghost_width == 5)
         if(kk/2*2==kk)then
             tmp2= C1*(ya(:,:,1)+ya(:,:,10))+C2*(ya(:,:,2)+ya(:,:,9)) &
                  +C3*(ya(:,:,3)+ya(:,:,8 ))+C4*(ya(:,:,4)+ya(:,:,7)) &
                  +C5*(ya(:,:,5)+ya(:,:,6 ))
         else
             tmp2= ya(:,:,5)
         endif
         if(jj/2*2==jj)then
             tmp1= C1*(tmp2(:,1)+tmp2(:,10))+C2*(tmp2(:,2)+tmp2(:,9)) &
                  +C3*(tmp2(:,3)+tmp2(:,8 ))+C4*(tmp2(:,4)+tmp2(:,7)) &
                  +C5*(tmp2(:,5)+tmp2(:,6 ))
         else
             tmp1= tmp2(:,5)
         endif
#endif

         tmp1(cxB(1)+ghost_width-i+1:cxT(1)+ghost_width-i+1) = fh(cxB(1):cxT(1),j,k)

         call polint(X,tmp1,0.d0,funf(i,j,k),ddy,2*ghost_width)

! for y direction
       elseif (fg(2) .eq. 0)then

#if (ghost_width == 2)
         if(kk/2*2==kk)then
             tmp2= C1*(ya(:,:,1)+ya(:,:,4))+C2*(ya(:,:,2)+ya(:,:,3))
         else
             tmp2= ya(:,:,2)
         endif
         if(ii/2*2==ii)then
             tmp1= C1*(tmp2(1,:)+tmp2(4,:))+C2*(tmp2(2,:)+tmp2(3,:))
         else
             tmp1= tmp2(2,:)
         endif
#elif (ghost_width == 3)
         if(kk/2*2==kk)then
             tmp2= C1*(ya(:,:,1)+ya(:,:,6))+C2*(ya(:,:,2)+ya(:,:,5))+C3*(ya(:,:,3)+ya(:,:,4))
         else
             tmp2= ya(:,:,3)
         endif
         if(ii/2*2==ii)then
             tmp1= C1*(tmp2(1,:)+tmp2(6,:))+C2*(tmp2(2,:)+tmp2(5,:))+C3*(tmp2(3,:)+tmp2(4,:))
         else
             tmp1= tmp2(3,:)
         endif
#elif (ghost_width == 4)
         if(kk/2*2==kk)then
             tmp2= C1*(ya(:,:,1)+ya(:,:,8))+C2*(ya(:,:,2)+ya(:,:,7)) &
                  +C3*(ya(:,:,3)+ya(:,:,6))+C4*(ya(:,:,4)+ya(:,:,5))
         else
             tmp2= ya(:,:,4)
         endif
         if(ii/2*2==ii)then
             tmp1= C1*(tmp2(1,:)+tmp2(8,:))+C2*(tmp2(2,:)+tmp2(7,:)) &
                  +C3*(tmp2(3,:)+tmp2(6,:))+C4*(tmp2(4,:)+tmp2(5,:))
         else
             tmp1= tmp2(4,:)
         endif
#elif (ghost_width == 5)
         if(kk/2*2==kk)then
             tmp2= C1*(ya(:,:,1)+ya(:,:,10))+C2*(ya(:,:,2)+ya(:,:,9)) &
                  +C3*(ya(:,:,3)+ya(:,:,8 ))+C4*(ya(:,:,4)+ya(:,:,7)) &
                  +C5*(ya(:,:,5)+ya(:,:,6 ))
         else
             tmp2= ya(:,:,5)
         endif
         if(ii/2*2==ii)then
             tmp1= C1*(tmp2(1,:)+tmp2(10,:))+C2*(tmp2(2,:)+tmp2(9,:)) &
                  +C3*(tmp2(3,:)+tmp2(8 ,:))+C4*(tmp2(4,:)+tmp2(7,:)) &
                  +C5*(tmp2(5,:)+tmp2(6 ,:))
         else
             tmp1= tmp2(5,:)
         endif
#endif

         tmp1(cxB(2)+ghost_width-j+1:cxT(2)+ghost_width-j+1) = fh(i,cxB(2):cxT(2),k)

         call polint(Y,tmp1,0.d0,funf(i,j,k),ddy,2*ghost_width)

! for z direction
       else

#if (ghost_width == 2)
         if(ii/2*2==ii)then
             tmp2= C1*(ya(1,:,:)+ya(4,:,:))+C2*(ya(2,:,:)+ya(3,:,:))
         else
             tmp2= ya(2,:,:)
         endif
         if(jj/2*2==jj)then
             tmp1= C1*(tmp2(1,:)+tmp2(4,:))+C2*(tmp2(2,:)+tmp2(3,:))
         else
             tmp1= tmp2(2,:)
         endif
#elif (ghost_width == 3)
         if(ii/2*2==ii)then
             tmp2= C1*(ya(1,:,:)+ya(6,:,:))+C2*(ya(6,:,:)+ya(5,:,:))+C3*(ya(3,:,:)+ya(4,:,:))
         else
             tmp2= ya(3,:,:)
         endif
         if(jj/2*2==jj)then
             tmp1= C1*(tmp2(1,:)+tmp2(6,:))+C2*(tmp2(2,:)+tmp2(5,:))+C3*(tmp2(3,:)+tmp2(4,:))
         else
             tmp1= tmp2(3,:)
         endif
#elif (ghost_width == 4)
         if(ii/2*2==ii)then
             tmp2= C1*(ya(1,:,:)+ya(8,:,:))+C2*(ya(2,:,:)+ya(7,:,:)) &
                  +C3*(ya(3,:,:)+ya(6,:,:))+C4*(ya(4,:,:)+ya(5,:,:))
         else
             tmp2= ya(4,:,:)
         endif
         if(jj/2*2==jj)then
             tmp1= C1*(tmp2(1,:)+tmp2(8,:))+C2*(tmp2(2,:)+tmp2(7,:)) &
                  +C3*(tmp2(3,:)+tmp2(6,:))+C4*(tmp2(4,:)+tmp2(5,:))
         else
             tmp1= tmp2(4,:)
         endif
#elif (ghost_width == 5)
         if(ii/2*2==ii)then
             tmp2= C1*(ya(1,:,:)+ya(10,:,:))+C2*(ya(2,:,:)+ya(9,:,:)) &
                  +C3*(ya(3,:,:)+ya(8 ,:,:))+C4*(ya(4,:,:)+ya(7,:,:)) &
                  +C5*(ya(5,:,:)+ya(6 ,:,:))
         else
             tmp2= ya(5,:,:)
         endif
         if(jj/2*2==jj)then
             tmp1= C1*(tmp2(1,:)+tmp2(10,:))+C2*(tmp2(2,:)+tmp2(9,:)) &
                  +C3*(tmp2(3,:)+tmp2(8 ,:))+C4*(tmp2(4,:)+tmp2(7,:)) &
                  +C5*(tmp2(5,:)+tmp2(6 ,:))
         else
             tmp1= tmp2(5,:)
         endif
#endif

         tmp1(cxB(3)+ghost_width-k+1:cxT(3)+ghost_width-k+1) = fh(i,j,cxB(3):cxT(3))

         call polint(Z,tmp1,0.d0,funf(i,j,k),ddy,2*ghost_width)

       endif
      
       else

       ii=i+lbf(1)-1
       jj=j+lbf(2)-1
       kk=k+lbf(3)-1

#if (ghost_width == 2)
         if(kk/2*2==kk)then
             tmp2= C1*(ya(:,:,1)+ya(:,:,4))+C2*(ya(:,:,2)+ya(:,:,3))
         else
             tmp2= ya(:,:,2)
         endif
         if(jj/2*2==jj)then
             tmp1= C1*(tmp2(:,1)+tmp2(:,4))+C2*(tmp2(:,2)+tmp2(:,3))
         else
             tmp1= tmp2(:,2)
         endif
         if(ii/2*2==ii)then
             funf(i,j,k)= C1*(tmp1(1)+tmp1(4))+C2*(tmp1(2)+tmp1(3))
         else
             funf(i,j,k)= tmp1(2)
         endif
#elif (ghost_width == 3)
         if(kk/2*2==kk)then
             tmp2= C1*(ya(:,:,1)+ya(:,:,6))+C2*(ya(:,:,2)+ya(:,:,5))+C3*(ya(:,:,3)+ya(:,:,4))
         else
             tmp2= ya(:,:,3)
         endif
         if(jj/2*2==jj)then
             tmp1= C1*(tmp2(:,1)+tmp2(:,6))+C2*(tmp2(:,2)+tmp2(:,5))+C3*(tmp2(:,3)+tmp2(:,4))
         else
             tmp1= tmp2(:,3)
         endif
         if(ii/2*2==ii)then
             funf(i,j,k)= C1*(tmp1(1)+tmp1(6))+C2*(tmp1(2)+tmp1(5))+C3*(tmp1(3)+tmp1(4))
         else
             funf(i,j,k)= tmp1(3)
         endif
#elif (ghost_width == 4)
         if(kk/2*2==kk)then
             tmp2= C1*(ya(:,:,1)+ya(:,:,8))+C2*(ya(:,:,2)+ya(:,:,7)) &
                  +C3*(ya(:,:,3)+ya(:,:,6))+C4*(ya(:,:,4)+ya(:,:,5))
         else
             tmp2= ya(:,:,4)
         endif
         if(jj/2*2==jj)then
             tmp1= C1*(tmp2(:,1)+tmp2(:,8))+C2*(tmp2(:,2)+tmp2(:,7)) &
                  +C3*(tmp2(:,3)+tmp2(:,6))+C4*(tmp2(:,4)+tmp2(:,5))
         else
             tmp1= tmp2(:,4)
         endif
         if(ii/2*2==ii)then
             funf(i,j,k)= C1*(tmp1(1)+tmp1(8))+C2*(tmp1(2)+tmp1(7)) &
                         +C3*(tmp1(3)+tmp1(6))+C4*(tmp1(4)+tmp1(5))
         else
             funf(i,j,k)= tmp1(4)
         endif
#elif (ghost_width == 5)
         if(kk/2*2==kk)then
             tmp2= C1*(ya(:,:,1)+ya(:,:,10))+C2*(ya(:,:,2)+ya(:,:,9)) &
                  +C3*(ya(:,:,3)+ya(:,:,8 ))+C4*(ya(:,:,4)+ya(:,:,7)) &
                  +C5*(ya(:,:,5)+ya(:,:,6 ))
         else
             tmp2= ya(:,:,5)
         endif
         if(jj/2*2==jj)then
             tmp1= C1*(tmp2(:,1)+tmp2(:,10))+C2*(tmp2(:,2)+tmp2(:,9)) &
                  +C3*(tmp2(:,3)+tmp2(:,8 ))+C4*(tmp2(:,4)+tmp2(:,7)) &
                  +C5*(tmp2(:,5)+tmp2(:,6 ))
         else
             tmp1= tmp2(:,5)
         endif
         if(ii/2*2==ii)then
             funf(i,j,k)= C1*(tmp1(1)+tmp1(10))+C2*(tmp1(2)+tmp1(9)) &
                         +C3*(tmp1(3)+tmp1(8 ))+C4*(tmp1(4)+tmp1(7)) &
                         +C5*(tmp1(5)+tmp1(6 ))
         else
             funf(i,j,k)= tmp1(5)
         endif
#endif
       endif

  enddo
  enddo
  enddo

  return

  end subroutine prolongmix3
!///////////////////////////////////////////////////////////////////////////////////////////////
!--------------------------------------------------------------------------------------
!
! Restrict from finner grids to coarser grids ignore the boundary point
! this routine is valid for all orders finite difference
!
! 1   2   3   4
! *---*---*---*
!     ^
! COPY directly!
!--------------------------------------------------------------------------------------

  subroutine restrict3(wei,llbc,uubc,extc,func,&
                       llbf,uubf,extf,funf,&
                       llbr,uubr,SoA,Symmetry)
  implicit none

!~~~~~~> input arguments
  integer,intent(in) :: wei
!                                       coarse    fine       coarse
  real*8,dimension(3),   intent(in) :: llbc,uubc,llbf,uubf,llbr,uubr
  integer,dimension(3),  intent(in) :: extc,extf
  real*8, dimension(extc(1),extc(2),extc(3)),intent(inout):: func
  real*8, dimension(extf(1),extf(2),extf(3)),intent(in):: funf
  real*8, dimension(1:3), intent(in) :: SoA
  integer,intent(in)::Symmetry

!~~~~~~> local variables

  real*8, dimension(1:3) :: base
  integer,dimension(3) :: lbc,ubc,lbf,ubf,lbr,ubr,lbrf,ubrf
  integer,dimension(3) :: cxB,cxT,cxI
  integer :: i,j,k
  real*8,  dimension(4,4,4) :: ya
  real*8, dimension(4,4) :: tmp2
  real*8, dimension(4) :: tmp1
  real*8, parameter :: C1=-1.d0/1.6d1,C2=9.d0/1.6d1

  integer::imini,imaxi,jmini,jmaxi,kmini,kmaxi
  integer::imino,imaxo,jmino,jmaxo,kmino,kmaxo
  real*8,dimension(3) :: CD,FD

  if(wei.ne.3)then
     write(*,*)"prolongrestrict.f90::restrict3: this routine only surport 3 dimension"
     write(*,*)"dim = ",wei
     stop
  endif
! it's possible a iolated point for target but not for source
  FD = (uubf-llbf)/(extf-1)
  CD = 2*FD

!take care the mismatch of the two segments of grid
  do i=1,3
     if(llbc(i) <= llbf(i))then
        base(i) = llbc(i)
     else
        j=idint((llbc(i)-llbf(i))/FD(i)+0.4)
        if(j/2*2 == j)then
           base(i) = llbf(i)
        else
           base(i) = llbf(i) - CD(i)/2
        endif
     endif
  enddo
!!! function idint:
!If A is of type REAL and |A| < 1, INT(A) equals 0. If |A| \geq 1, 
!then INT(A) equals the largest integer that does not exceed the range of A 
!and whose sign is the same as the sign of A.

    lbf = idint((llbf-base)/FD+0.4)+1
    ubf = idint((uubf-base)/FD+0.4)+1
    lbc = idint((llbc-base)/CD+0.4)+1
    ubc = idint((uubc-base)/CD+0.4)+1
    lbr = idint((llbr-base)/CD+0.4)+1
    lbrf = idint((llbr-base)/FD+0.4)+1
    ubr = idint((uubr-base)/CD+0.4)+1
    ubrf = idint((uubr-base)/FD+0.4)+1

!sanity check
  imino=lbr(1)-lbc(1) + 1
  imaxo=ubr(1)-lbc(1) + 1
  jmino=lbr(2)-lbc(2) + 1
  jmaxo=ubr(2)-lbc(2) + 1
  kmino=lbr(3)-lbc(3) + 1
  kmaxo=ubr(3)-lbc(3) + 1

  imini=lbrf(1)-lbf(1) + 1
  imaxi=ubrf(1)-lbf(1) + 1
  jmini=lbrf(2)-lbf(2) + 1
  jmaxi=ubrf(2)-lbf(2) + 1
  kmini=lbrf(3)-lbf(3) + 1
  kmaxi=ubrf(3)-lbf(3) + 1

  if(imino.lt.1.or.jmino.lt.1.or.kmino.lt.1.or.&
     imini.lt.1.or.jmini.lt.1.or.kmini.lt.1.or.&
     imaxo.gt.extc(1).or.jmaxo.gt.extc(2).or.kmaxo.gt.extc(3).or.&
     imaxi.gt.extf(1).or.jmaxi.gt.extf(2).or.kmaxi.gt.extf(3))then
          write(*,*)"error in restrict for"
          write(*,*)"mino = ",imino,jmino,kmino
          write(*,*)"maxo = ",imaxo,jmaxo,kmaxo
          write(*,*)"extc = ",extc
          write(*,*)"CD = ",CD
          write(*,*)"mini = ",imini,jmini,kmini
          write(*,*)"maxi = ",imaxi,jmaxi,kmaxi
          write(*,*)"extf = ",extf
          write(*,*)"FD = ",FD
          write(*,*)"from"
          write(*,*)lbf,ubf,extf
          write(*,*)"to"
          write(*,*)lbc,ubc,extc
          write(*,*)"want"
          write(*,*)lbr,ubr,lbrf,ubrf
          write(*,*)"llbf = ",llbf
          write(*,*)"uubf = ",uubf
          write(*,*)"llbc = ",llbc
          write(*,*)"uubc = ",uubc
          write(*,*)"llbr = ",llbr
          write(*,*)"uubr = ",uubr
          stop
  endif

!~~~~~~> restriction start...
  do k = kmino,kmaxo
   do j = jmino,jmaxo
    do i = imino,imaxo

       cxI(1) = i
       cxI(2) = j
       cxI(3) = k
! change to fine level reference
!|*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*| 
!|x===============x===============x===============x========|
       cxI = 2*(cxI+lbc-1) - 1
! change to array index      
       cxI = cxI - lbf + 1

       func(i,j,k)= funf(cxI(1),cxI(2),cxI(3))
    enddo
   enddo
  enddo
  
  return

  end subroutine restrict3
!===========================================================================================

! for different finite differnce order
#if (ghost_width == 2)
! 2nd order
!--------------------------------------------------------------------------
!
! Prolongation from coarser grids to finer grids
! 4 points, 3rd order interpolation
! 1   2   3   4
! *---*---*---*
!       ^
! f=-1/16*f_1 + 9/16*f_2
!   -1/16*f_4 + 9/16*f_3
!--------------------------------------------------------------------------

  subroutine prolong3(wei,llbc,uubc,extc,func,&
                      llbf,uubf,extf,funf,&
                      llbp,uubp,SoA,Symmetry)
  implicit none

!~~~~~~> input arguments
  integer,intent(in) :: wei
!                                       coarse    fine       coarse
  real*8,dimension(3),   intent(in) :: llbc,uubc,llbf,uubf,llbp,uubp
  integer,dimension(3),  intent(in) :: extc,extf
  real*8, dimension(extc(1),extc(2),extc(3)),intent(in)   :: func
  real*8, dimension(extf(1),extf(2),extf(3)),intent(inout):: funf
  real*8, dimension(1:3), intent(in) :: SoA
  integer,intent(in)::Symmetry

!~~~~~~> local variables

  real*8, dimension(1:3) :: base
  integer,dimension(3) :: lbc,ubc,lbf,ubf,lbp,ubp,lbpc,ubpc
  integer,dimension(3) :: cxB,cxT,cxI
  integer :: i,j,k,ii,jj,kk
  real*8,  dimension(4,4,4) :: ya
  real*8, dimension(4,4) :: tmp2
  real*8, dimension(4) :: tmp1

  real*8, parameter :: C1=-1.d0/16,C2=9.d0/16
  real*8, parameter :: C4=C1,C3=C2

  integer::imini,imaxi,jmini,jmaxi,kmini,kmaxi
  integer::imino,imaxo,jmino,jmaxo,kmino,kmaxo
  logical::decide3d

  real*8,dimension(3) :: CD,FD

  if(wei.ne.3)then
     write(*,*)"prolongrestrict.f90::prolong3: this routine only surport 3 dimension"
     write(*,*)"dim = ",wei
     stop
  endif

! it's possible a iolated point for target but not for source
  CD = (uubc-llbc)/(extc-1)
  FD = CD/2

!take care the mismatch of the two segments of grid
  do i=1,3
     if(llbc(i) <= llbf(i))then
        base(i) = llbc(i)
     else
        j=idint((llbc(i)-llbf(i))/FD(i)+0.4)
        if(j/2*2 == j)then
           base(i) = llbf(i)
        else
           base(i) = llbf(i) - CD(i)/2
        endif
     endif
  enddo

!!! function idint:
!If A is of type REAL and |A| < 1, INT(A) equals 0. If |A| \geq 1, 
!then INT(A) equals the largest integer that does not exceed the range of A 
!and whose sign is the same as the sign of A.

    lbf = idint((llbf-base)/FD+0.4)+1
    ubf = idint((uubf-base)/FD+0.4)+1
    lbc = idint((llbc-base)/CD+0.4)+1
    ubc = idint((uubc-base)/CD+0.4)+1
    lbp = idint((llbp-base)/FD+0.4)+1
    lbpc = idint((llbp-base)/CD+0.4)+1
    ubp = idint((uubp-base)/FD+0.4)+1
    ubpc = idint((uubp-base)/CD+0.4)+1
!sanity check
  imino=lbp(1)-lbf(1) + 1
  imaxo=ubp(1)-lbf(1) + 1
  jmino=lbp(2)-lbf(2) + 1
  jmaxo=ubp(2)-lbf(2) + 1
  kmino=lbp(3)-lbf(3) + 1
  kmaxo=ubp(3)-lbf(3) + 1

  imini=lbpc(1)-lbc(1) + 1
  imaxi=ubpc(1)-lbc(1) + 1
  jmini=lbpc(2)-lbc(2) + 1
  jmaxi=ubpc(2)-lbc(2) + 1
  kmini=lbpc(3)-lbc(3) + 1
  kmaxi=ubpc(3)-lbc(3) + 1

  if(imino.lt.1.or.jmino.lt.1.or.kmino.lt.1.or.&
     imini.lt.1.or.jmini.lt.1.or.kmini.lt.1.or.&
     imaxo.gt.extf(1).or.jmaxo.gt.extf(2).or.kmaxo.gt.extf(3).or.&
     imaxi.gt.extc(1)-1.or.jmaxi.gt.extc(2)-1.or.kmaxi.gt.extc(3)-1)then
          write(*,*)"error in prolongation for"
          write(*,*)"from"
          write(*,*)llbc,uubc
          write(*,*)lbc,ubc
          write(*,*)"to"
          write(*,*)llbf,uubf
          write(*,*)lbf,ubf
          write(*,*)"want"
          write(*,*)llbp,uubp
          write(*,*)lbp,ubp,lbpc,ubpc
          if(imini.lt.1) write(*,*)"imini = ",imini
          if(jmini.lt.1) write(*,*)"jmini = ",jmini       
          if(kmini.lt.1) write(*,*)"kmini = ",kmini      
          if(imino.lt.1) write(*,*)"imino = ",imino       
          if(jmino.lt.1) write(*,*)"jmino = ",jmino       
          if(kmino.lt.1) write(*,*)"kmino = ",kmino   
          if(imaxi.gt.extc(1)) write(*,*)"imaxi = ",imaxi,"extc(1) = ",extc(1) 
          if(jmaxi.gt.extc(2)) write(*,*)"jmaxi = ",jmaxi,"extc(2) = ",extc(2) 
          if(kmaxi.gt.extc(3)) write(*,*)"kmaxi = ",kmaxi,"extc(3) = ",extc(3) 
          if(imaxo.gt.extf(1)) write(*,*)"imaxo = ",imaxo,"extf(1) = ",extf(1) 
          if(jmaxo.gt.extf(2)) write(*,*)"jmaxo = ",jmaxo,"extf(2) = ",extf(2) 
          if(kmaxo.gt.extf(3)) write(*,*)"kmaxo = ",kmaxo,"extf(3) = ",extf(3) 
          return
  endif
!~~~~~~> prolongation start...
  do k = kmino,kmaxo
   do j = jmino,jmaxo
    do i = imino,imaxo
! change to coarse level reference
!|*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*| 
!|x===============x===============x===============x========|
!       if(i/2*2 == i)then
!         cxI(1) = (i+lbf(1)-1)/2
!       else
!         cxI(1) = (i+lbf(1)-1)/2+1
!       endif      
!       if(j/2*2 == j)then
!         cxI(2) = (j+lbf(2)-1)/2
!       else
!         cxI(2) = (j+lbf(2)-1)/2+1
!       endif
!       if(k/2*2 == k)then
!         cxI(3) = (k+lbf(3)-1)/2
!       else
!         cxI(3) = (k+lbf(3)-1)/2+1
!       endif
! above code segment is equivalent to 
       cxI(1) = i
       cxI(2) = j
       cxI(3) = k
       cxI = (cxI+lbf)/2
! change to array index      
       cxI = cxI - lbc + 1

       ii=i+lbf(1)-1
       jj=j+lbf(2)-1
       kk=k+lbf(3)-1

       if(any(cxI+2 > extc)) write(*,*)"error in prolong"
       if(ii/2*2==ii)then
         if(jj/2*2==jj)then
           if(kk/2*2==kk)then
!  due to ghost zone, we can deal with symmetry boundary like this                   
             if(cxI(1)>1.and.cxI(2)>1.and.cxI(3)>1)then
                tmp2= C1*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)-1)+&
                      C2*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)  )+&
                      C3*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)+1)+&
                      C4*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)+2)
             else
                cxB=cxI-1
                cxT=cxI+2
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,4,Symmetry))then
                        write(*,*)"prolong3 position index: ",i+lbf(1)-1,j+lbf(2)-1,k+lbf(3)-1
                        return
                endif
                tmp2= C1*ya(:,:,1)+C2*ya(:,:,2)+C3*ya(:,:,3)+C4*ya(:,:,4)
             endif
             tmp1= C1*tmp2(:,1)+C2*tmp2(:,2)+C3*tmp2(:,3)+C4*tmp2(:,4)
             funf(i,j,k)= C1*tmp1(1)+C2*tmp1(2)+C3*tmp1(3)+C4*tmp1(4)
           else
             if(cxI(1)>1.and.cxI(2)>1.and.cxI(3)>1)then
                tmp2= func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3))
             else
                cxB=cxI-1
                cxT=cxI+2
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,4,Symmetry))then
                        write(*,*)"prolong3 position index: ",i+lbf(1)-1,j+lbf(2)-1,k+lbf(3)-1
                        return
                endif
                tmp2= ya(:,:,2)
             endif
             tmp1= C1*tmp2(:,1)+C2*tmp2(:,2)+C3*tmp2(:,3)+C4*tmp2(:,4)
             funf(i,j,k)=  C1*tmp1(1)+C2*tmp1(2)+C3*tmp1(3)+C4*tmp1(4)
           endif
         else
           if(kk/2*2==kk)then
             if(cxI(1)>1.and.cxI(2)>1.and.cxI(3)>1)then
                tmp2= C1*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)-1)+&
                      C2*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)  )+&
                      C3*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)+1)+&
                      C4*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)+2)
             else
                cxB=cxI-1
                cxT=cxI+2
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,4,Symmetry))then
                        write(*,*)"prolong3 position index: ",i+lbf(1)-1,j+lbf(2)-1,k+lbf(3)-1
                        return
                endif
                tmp2= C1*ya(:,:,1)+C2*ya(:,:,2)+C3*ya(:,:,3)+C4*ya(:,:,4)
             endif
             tmp1= tmp2(:,2)
             funf(i,j,k)= C1*tmp1(1)+C2*tmp1(2)+C3*tmp1(3)+C4*tmp1(4)
           else
             if(cxI(1)>1.and.cxI(2)>1.and.cxI(3)>1)then
                tmp2= func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3))
             else
                cxB=cxI-1
                cxT=cxI+2
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,4,Symmetry))then
                        write(*,*)"prolong3 position index: ",i+lbf(1)-1,j+lbf(2)-1,k+lbf(3)-1
                        return
                endif
                tmp2= ya(:,:,2)
             endif
             tmp1= tmp2(:,2)
             funf(i,j,k)=  C1*tmp1(1)+C2*tmp1(2)+C3*tmp1(3)+C4*tmp1(4)
           endif
         endif
       else
         if(jj/2*2==jj)then
           if(kk/2*2==kk)then               
             if(cxI(1)>1.and.cxI(2)>1.and.cxI(3)>1)then
                tmp2= C1*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)-1)+&
                      C2*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)  )+&
                      C3*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)+1)+&
                      C4*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)+2)
             else
                cxB=cxI-1
                cxT=cxI+2
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,4,Symmetry))then
                        write(*,*)"prolong3 position index: ",i+lbf(1)-1,j+lbf(2)-1,k+lbf(3)-1
                        return
                endif
                tmp2= C1*ya(:,:,1)+C2*ya(:,:,2)+C3*ya(:,:,3)+C4*ya(:,:,4)
             endif
             tmp1= C1*tmp2(:,1)+C2*tmp2(:,2)+C3*tmp2(:,3)+C4*tmp2(:,4)
             funf(i,j,k)= tmp1(2)
           else
             if(cxI(1)>1.and.cxI(2)>1.and.cxI(3)>1)then
                tmp2= func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3))
             else
                cxB=cxI-1
                cxT=cxI+2
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,4,Symmetry))then
                        write(*,*)"prolong3 position index: ",i+lbf(1)-1,j+lbf(2)-1,k+lbf(3)-1
                        return
                endif
                tmp2= ya(:,:,2)
             endif
             tmp1= C1*tmp2(:,1)+C2*tmp2(:,2)+C3*tmp2(:,3)+C4*tmp2(:,4)
             funf(i,j,k)= tmp1(2)
           endif
         else
           if(kk/2*2==kk)then
             if(cxI(1)>1.and.cxI(2)>1.and.cxI(3)>1)then
                tmp2= C1*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)-1)+&
                      C2*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)  )+&
                      C3*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)+1)+&
                      C4*func(cxI(1)-1:cxI(1)+2,cxI(2)-1:cxI(2)+2,cxI(3)+2)
             else
                cxB=cxI-1
                cxT=cxI+2
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,4,Symmetry))then
                        write(*,*)"prolong3 position index: ",i+lbf(1)-1,j+lbf(2)-1,k+lbf(3)-1
                        return
                endif
                tmp2= C1*ya(:,:,1)+C2*ya(:,:,2)+C3*ya(:,:,3)+C4*ya(:,:,4)
             endif
             tmp1= tmp2(:,2)
             funf(i,j,k)= tmp1(2)
           else
             funf(i,j,k)= func(cxI(1),cxI(2),cxI(3))
           endif
         endif
       endif
    enddo
   enddo
  enddo

  return

  end subroutine prolong3

#elif (ghost_width == 3)
! fourth order code
!--------------------------------------------------------------------------
!
! Prolongation from coarser grids to finer grids
! 6 points, 5th order interpolation
! 1   2   3   4   5   6
! *---*---*---*---*---*
!           ^
! f=3/256*(f_1+f_6) - 25/256*(f_2+f_5) + 75/128*(f_3+f_4)
!--------------------------------------------------------------------------

  subroutine prolong3(wei,llbc,uubc,extc,func,&
                      llbf,uubf,extf,funf,&
                      llbp,uubp,SoA,Symmetry)
  implicit none

!~~~~~~> input arguments
  integer,intent(in) :: wei
!                                       coarse    fine       coarse
  real*8,dimension(3),   intent(in) :: llbc,uubc,llbf,uubf,llbp,uubp
  integer,dimension(3),  intent(in) :: extc,extf
  real*8, dimension(extc(1),extc(2),extc(3)),intent(in)   :: func
  real*8, dimension(extf(1),extf(2),extf(3)),intent(inout):: funf
  real*8, dimension(1:3), intent(in) :: SoA
  integer,intent(in)::Symmetry

!~~~~~~> local variables

  real*8, dimension(1:3) :: base
  integer,dimension(3) :: lbc,ubc,lbf,ubf,lbp,ubp,lbpc,ubpc
  integer,dimension(3) :: cxB,cxT,cxI
  integer :: i,j,k,ii,jj,kk
  real*8, dimension(6,6) :: tmp2
  real*8, dimension(6) :: tmp1

  real*8, parameter :: C1=3.d0/2.56d2,C2=-2.5d1/2.56d2,C3=7.5d1/1.28d2

  integer::imini,imaxi,jmini,jmaxi,kmini,kmaxi
  integer::imino,imaxo,jmino,jmaxo,kmino,kmaxo
  real*8, dimension(-1:extc(1),-1:extc(2),-1:extc(3))   :: funcc

  real*8,dimension(3) :: CD,FD

  if(wei.ne.3)then
     write(*,*)"prolongrestrict.f90::prolong3: this routine only surport 3 dimension"
     write(*,*)"dim = ",wei
     stop
  endif

! it's possible a iolated point for target but not for source
  CD = (uubc-llbc)/(extc-1)
  FD = CD/2

!take care the mismatch of the two segments of grid
  do i=1,3
     if(llbc(i) <= llbf(i))then
        base(i) = llbc(i)
     else
        j=idint((llbc(i)-llbf(i))/FD(i)+0.4)
        if(j/2*2 == j)then
           base(i) = llbf(i)
        else
           base(i) = llbf(i) - CD(i)/2
        endif
     endif
  enddo

!!! function idint:
!If A is of type REAL and |A| < 1, INT(A) equals 0. If |A| \geq 1, 
!then INT(A) equals the largest integer that does not exceed the range of A 
!and whose sign is the same as the sign of A.

    lbf = idint((llbf-base)/FD+0.4)+1
    ubf = idint((uubf-base)/FD+0.4)+1
    lbc = idint((llbc-base)/CD+0.4)+1
    ubc = idint((uubc-base)/CD+0.4)+1
    lbp = idint((llbp-base)/FD+0.4)+1
    lbpc = idint((llbp-base)/CD+0.4)+1
    ubp = idint((uubp-base)/FD+0.4)+1
    ubpc = idint((uubp-base)/CD+0.4)+1
!sanity check
  imino=lbp(1)-lbf(1) + 1
  imaxo=ubp(1)-lbf(1) + 1
  jmino=lbp(2)-lbf(2) + 1
  jmaxo=ubp(2)-lbf(2) + 1
  kmino=lbp(3)-lbf(3) + 1
  kmaxo=ubp(3)-lbf(3) + 1

  imini=lbpc(1)-lbc(1) + 1
  imaxi=ubpc(1)-lbc(1) + 1
  jmini=lbpc(2)-lbc(2) + 1
  jmaxi=ubpc(2)-lbc(2) + 1
  kmini=lbpc(3)-lbc(3) + 1
  kmaxi=ubpc(3)-lbc(3) + 1

  if(imino.lt.1.or.jmino.lt.1.or.kmino.lt.1.or.&
     imini.lt.1.or.jmini.lt.1.or.kmini.lt.1.or.&
     imaxo.gt.extf(1).or.jmaxo.gt.extf(2).or.kmaxo.gt.extf(3).or.&
     imaxi.gt.extc(1)-2.or.jmaxi.gt.extc(2)-2.or.kmaxi.gt.extc(3)-2)then
          write(*,*)"error in prolongation for"
          write(*,*)"from"
          write(*,*)llbc,uubc
          write(*,*)lbc,ubc
          write(*,*)"to"
          write(*,*)llbf,uubf
          write(*,*)lbf,ubf
          write(*,*)"want"
          write(*,*)llbp,uubp
          write(*,*)lbp,ubp,lbpc,ubpc
          if(imini.lt.1) write(*,*)"imini = ",imini
          if(jmini.lt.1) write(*,*)"jmini = ",jmini       
          if(kmini.lt.1) write(*,*)"kmini = ",kmini      
          if(imino.lt.1) write(*,*)"imino = ",imino       
          if(jmino.lt.1) write(*,*)"jmino = ",jmino       
          if(kmino.lt.1) write(*,*)"kmino = ",kmino   
          if(imaxi.gt.extc(1)) write(*,*)"imaxi = ",imaxi,"extc(1) = ",extc(1) 
          if(jmaxi.gt.extc(2)) write(*,*)"jmaxi = ",jmaxi,"extc(2) = ",extc(2) 
          if(kmaxi.gt.extc(3)) write(*,*)"kmaxi = ",kmaxi,"extc(3) = ",extc(3) 
          if(imaxo.gt.extf(1)) write(*,*)"imaxo = ",imaxo,"extf(1) = ",extf(1) 
          if(jmaxo.gt.extf(2)) write(*,*)"jmaxo = ",jmaxo,"extf(2) = ",extf(2) 
          if(kmaxo.gt.extf(3)) write(*,*)"kmaxo = ",kmaxo,"extf(3) = ",extf(3) 
          return
  endif
  
  call symmetry_bd(2,extc,func,funcc,SoA)

!~~~~~~> prolongation start...
  do k = kmino,kmaxo
   do j = jmino,jmaxo
    do i = imino,imaxo
! change to coarse level reference        v
!|*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*| 
!|x===============x===============x===============x===============x===============x|
       cxI(1) = i
       cxI(2) = j
       cxI(3) = k
       cxI = (cxI+lbf)/2
! change to array index      
       cxI = cxI - lbc + 1

       ii=i+lbf(1)-1
       jj=j+lbf(2)-1
       kk=k+lbf(3)-1

       if(any(cxI+3 > extc)) write(*,*)"error in prolong"
       if(ii/2*2==ii)then
         if(jj/2*2==jj)then
           if(kk/2*2==kk)then
!  due to ghost zone, we can deal with symmetry boundary like this       
             tmp2= C1*(funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)-2)+funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)+3))&
                  +C2*(funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)-1)+funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)+2))&
                  +C3*(funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)  )+funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)+1))
             tmp1= C1*(tmp2(:,1)+tmp2(:,6))+C2*(tmp2(:,2)+tmp2(:,5))+C3*(tmp2(:,3)+tmp2(:,4))
             funf(i,j,k)= C1*(tmp1(1)+tmp1(6))+C2*(tmp1(2)+tmp1(5))+C3*(tmp1(3)+tmp1(4))
           else
             tmp2= funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3))            
             tmp1= C1*(tmp2(:,1)+tmp2(:,6))+C2*(tmp2(:,2)+tmp2(:,5))+C3*(tmp2(:,3)+tmp2(:,4))
             funf(i,j,k)= C1*(tmp1(1)+tmp1(6))+C2*(tmp1(2)+tmp1(5))+C3*(tmp1(3)+tmp1(4))
           endif
         else
           if(kk/2*2==kk)then
             tmp2= C1*(funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)-2)+funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)+3))&
                  +C2*(funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)-1)+funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)+2))&
                  +C3*(funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)  )+funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)+1))
             tmp1= tmp2(:,3)
             funf(i,j,k)= C1*(tmp1(1)+tmp1(6))+C2*(tmp1(2)+tmp1(5))+C3*(tmp1(3)+tmp1(4))
           else
             tmp2= funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3))  
             tmp1= tmp2(:,3)
             funf(i,j,k)= C1*(tmp1(1)+tmp1(6))+C2*(tmp1(2)+tmp1(5))+C3*(tmp1(3)+tmp1(4))
           endif
         endif
       else
         if(jj/2*2==jj)then
           if(kk/2*2==kk)then    
             tmp2= C1*(funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)-2)+funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)+3))&
                  +C2*(funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)-1)+funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)+2))&
                  +C3*(funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)  )+funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)+1))
             tmp1= C1*(tmp2(:,1)+tmp2(:,6))+C2*(tmp2(:,2)+tmp2(:,5))+C3*(tmp2(:,3)+tmp2(:,4))
             funf(i,j,k)= tmp1(3)
           else
             tmp2= funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3))  
             tmp1= C1*(tmp2(:,1)+tmp2(:,6))+C2*(tmp2(:,2)+tmp2(:,5))+C3*(tmp2(:,3)+tmp2(:,4))
             funf(i,j,k)= tmp1(3)
           endif
         else
           if(kk/2*2==kk)then
             tmp2= C1*(funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)-2)+funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)+3))&
                  +C2*(funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)-1)+funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)+2))&
                  +C3*(funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)  )+funcc(cxI(1)-2:cxI(1)+3,cxI(2)-2:cxI(2)+3,cxI(3)+1))
             tmp1= tmp2(:,3)
             funf(i,j,k)= tmp1(3)
           else
             funf(i,j,k)= funcc(cxI(1),cxI(2),cxI(3))
           endif
         endif
       endif
    enddo
   enddo
  enddo

  return

  end subroutine prolong3

#elif (ghost_width == 4)
! sixth order code
!--------------------------------------------------------------------------
!
! Prolongation from coarser grids to finer grids
! 8 points, 7th order interpolation
! 1   2   3   4   5   6   7   8
! *---*---*---*---*---*---*---*
!               ^
! f=-5/2048*(f_1+f_8) + 49/2048*(f_2+f_7) - 245/2048*(f_3+f_6) + 1225/2048*(f_4+f_5)
!--------------------------------------------------------------------------

  subroutine prolong3(wei,llbc,uubc,extc,func,&
                      llbf,uubf,extf,funf,&
                      llbp,uubp,SoA,Symmetry)
  implicit none

!~~~~~~> input arguments
  integer,intent(in) :: wei
!                                       coarse    fine       coarse
  real*8,dimension(3),   intent(in) :: llbc,uubc,llbf,uubf,llbp,uubp
  integer,dimension(3),  intent(in) :: extc,extf
  real*8, dimension(extc(1),extc(2),extc(3)),intent(in)   :: func
  real*8, dimension(extf(1),extf(2),extf(3)),intent(inout):: funf
  real*8, dimension(1:3), intent(in) :: SoA
  integer,intent(in)::Symmetry

!~~~~~~> local variables

  real*8, dimension(1:3) :: base
  integer,dimension(3) :: lbc,ubc,lbf,ubf,lbp,ubp,lbpc,ubpc
  integer,dimension(3) :: cxB,cxT,cxI
  integer :: i,j,k,ii,jj,kk
  real*8,  dimension(8,8,8) :: ya
  real*8, dimension(8,8) :: tmp2
  real*8, dimension(8) :: tmp1

  real*8, parameter :: C1=-5.d0/2.048d3,C2=4.9d1/2.048d3,C3=-2.45d2/2.048d3,C4=1.225d3/2.048d3

  integer::imini,imaxi,jmini,jmaxi,kmini,kmaxi
  integer::imino,imaxo,jmino,jmaxo,kmino,kmaxo
  logical::decide3d

  real*8,dimension(3) :: CD,FD

  if(wei.ne.3)then
     write(*,*)"prolongrestrict.f90::prolong3: this routine only surport 3 dimension"
     write(*,*)"dim = ",wei
     stop
  endif

! it's possible a iolated point for target but not for source
  CD = (uubc-llbc)/(extc-1)
  FD = CD/2

!take care the mismatch of the two segments of grid
  do i=1,3
     if(llbc(i) <= llbf(i))then
        base(i) = llbc(i)
     else
        j=idint((llbc(i)-llbf(i))/FD(i)+0.4)
        if(j/2*2 == j)then
           base(i) = llbf(i)
        else
           base(i) = llbf(i) - CD(i)/2
        endif
     endif
  enddo

!!! function idint:
!If A is of type REAL and |A| < 1, INT(A) equals 0. If |A| \geq 1, 
!then INT(A) equals the largest integer that does not exceed the range of A 
!and whose sign is the same as the sign of A.

    lbf = idint((llbf-base)/FD+0.4)+1
    ubf = idint((uubf-base)/FD+0.4)+1
    lbc = idint((llbc-base)/CD+0.4)+1
    ubc = idint((uubc-base)/CD+0.4)+1
    lbp = idint((llbp-base)/FD+0.4)+1
    lbpc = idint((llbp-base)/CD+0.4)+1
    ubp = idint((uubp-base)/FD+0.4)+1
    ubpc = idint((uubp-base)/CD+0.4)+1
!sanity check
  imino=lbp(1)-lbf(1) + 1
  imaxo=ubp(1)-lbf(1) + 1
  jmino=lbp(2)-lbf(2) + 1
  jmaxo=ubp(2)-lbf(2) + 1
  kmino=lbp(3)-lbf(3) + 1
  kmaxo=ubp(3)-lbf(3) + 1

  imini=lbpc(1)-lbc(1) + 1
  imaxi=ubpc(1)-lbc(1) + 1
  jmini=lbpc(2)-lbc(2) + 1
  jmaxi=ubpc(2)-lbc(2) + 1
  kmini=lbpc(3)-lbc(3) + 1
  kmaxi=ubpc(3)-lbc(3) + 1

  if(imino.lt.1.or.jmino.lt.1.or.kmino.lt.1.or.&
     imini.lt.1.or.jmini.lt.1.or.kmini.lt.1.or.&
     imaxo.gt.extf(1).or.jmaxo.gt.extf(2).or.kmaxo.gt.extf(3).or.&
     imaxi.gt.extc(1)-3.or.jmaxi.gt.extc(2)-3.or.kmaxi.gt.extc(3)-3)then
          write(*,*)"error in prolongation for"
          write(*,*)"from"
          write(*,*)llbc,uubc
          write(*,*)lbc,ubc
          write(*,*)"to"
          write(*,*)llbf,uubf
          write(*,*)lbf,ubf
          write(*,*)"want"
          write(*,*)llbp,uubp
          write(*,*)lbp,ubp,lbpc,ubpc
          if(imini.lt.1) write(*,*)"imini = ",imini
          if(jmini.lt.1) write(*,*)"jmini = ",jmini       
          if(kmini.lt.1) write(*,*)"kmini = ",kmini      
          if(imino.lt.1) write(*,*)"imino = ",imino       
          if(jmino.lt.1) write(*,*)"jmino = ",jmino       
          if(kmino.lt.1) write(*,*)"kmino = ",kmino   
          if(imaxi.gt.extc(1)) write(*,*)"imaxi = ",imaxi,"extc(1) = ",extc(1) 
          if(jmaxi.gt.extc(2)) write(*,*)"jmaxi = ",jmaxi,"extc(2) = ",extc(2) 
          if(kmaxi.gt.extc(3)) write(*,*)"kmaxi = ",kmaxi,"extc(3) = ",extc(3) 
          if(imaxo.gt.extf(1)) write(*,*)"imaxo = ",imaxo,"extf(1) = ",extf(1) 
          if(jmaxo.gt.extf(2)) write(*,*)"jmaxo = ",jmaxo,"extf(2) = ",extf(2) 
          if(kmaxo.gt.extf(3)) write(*,*)"kmaxo = ",kmaxo,"extf(3) = ",extf(3) 
          return
  endif
!~~~~~~> prolongation start...
  do k = kmino,kmaxo
   do j = jmino,jmaxo
    do i = imino,imaxo
! change to coarse level reference                        v
!|*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*| 
!|x===============x===============x===============x===============x===============x===============x===============x|
       cxI(1) = i
       cxI(2) = j
       cxI(3) = k
       cxI = (cxI+lbf)/2
! change to array index      
       cxI = cxI - lbc + 1

       ii=i+lbf(1)-1
       jj=j+lbf(2)-1
       kk=k+lbf(3)-1

       if(any(cxI+4 > extc)) write(*,*)"error in prolong"
       if(ii/2*2==ii)then
         if(jj/2*2==jj)then
           if(kk/2*2==kk)then
!  due to ghost zone, we can deal with symmetry boundary like this                  
             if(cxI(1)>3.and.cxI(2)>3.and.cxI(3)>3)then
               tmp2= C1*(func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)-3)+func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+4))&
                    +C2*(func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)-2)+func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+3))&
                    +C3*(func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)-1)+func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+2))&
                    +C4*(func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)  )+func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+1))
             else
               cxB=cxI-3
               cxT=cxI+4
               if(decide3d(extc,func,func,cxB,cxT,SoA,ya,8,Symmetry))then
                 write(*,*)"prolong3 position: "
                 write(*,*)llbf(1)+(i-1)*FD(1),llbf(2)+(j-1)*FD(2),llbf(3)+(k-1)*FD(3)
                 write(*,*)"llbf = ",llbf
                 stop
               endif
               tmp2=C1*(ya(:,:,1)+ya(:,:,8))+C2*(ya(:,:,2)+ya(:,:,7))+C3*(ya(:,:,3)+ya(:,:,6))+C4*(ya(:,:,4)+ya(:,:,5))
             endif
             tmp1= C1*(tmp2(:,1)+tmp2(:,8))+C2*(tmp2(:,2)+tmp2(:,7))+C3*(tmp2(:,3)+tmp2(:,6))+C4*(tmp2(:,4)+tmp2(:,5))
             funf(i,j,k)= C1*(tmp1(1)+tmp1(8))+C2*(tmp1(2)+tmp1(7))+C3*(tmp1(3)+tmp1(6))+C4*(tmp1(4)+tmp1(5))
           else
             if(cxI(1)>3.and.cxI(2)>3.and.cxI(3)>3)then
                tmp2= func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3))
             else
                cxB=cxI-3
                cxT=cxI+4
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,8,Symmetry))then
                 write(*,*)"prolong3 position: "
                 write(*,*)llbf(1)+(i-1)*FD(1),llbf(2)+(j-1)*FD(2),llbf(3)+(k-1)*FD(3)
                 write(*,*)"llbf = ",llbf
                 stop
                endif
                tmp2= ya(:,:,4)
             endif
             tmp1= C1*(tmp2(:,1)+tmp2(:,8))+C2*(tmp2(:,2)+tmp2(:,7))+C3*(tmp2(:,3)+tmp2(:,6))+C4*(tmp2(:,4)+tmp2(:,5))
             funf(i,j,k)= C1*(tmp1(1)+tmp1(8))+C2*(tmp1(2)+tmp1(7))+C3*(tmp1(3)+tmp1(6))+C4*(tmp1(4)+tmp1(5))
           endif
         else
           if(kk/2*2==kk)then
             if(cxI(1)>3.and.cxI(2)>3.and.cxI(3)>3)then
               tmp2= C1*(func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)-3)+func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+4))&
                    +C2*(func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)-2)+func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+3))&
                    +C3*(func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)-1)+func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+2))&
                    +C4*(func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)  )+func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+1))
             else
                cxB=cxI-3
                cxT=cxI+4
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,8,Symmetry))then
                 write(*,*)"prolong3 position: "
                 write(*,*)llbf(1)+(i-1)*FD(1),llbf(2)+(j-1)*FD(2),llbf(3)+(k-1)*FD(3)
                 write(*,*)"llbf = ",llbf
                 stop
                endif
                tmp2=C1*(ya(:,:,1)+ya(:,:,8))+C2*(ya(:,:,2)+ya(:,:,7))+C3*(ya(:,:,3)+ya(:,:,6))+C4*(ya(:,:,4)+ya(:,:,5))
             endif
             tmp1= tmp2(:,4)
             funf(i,j,k)= C1*(tmp1(1)+tmp1(8))+C2*(tmp1(2)+tmp1(7))+C3*(tmp1(3)+tmp1(6))+C4*(tmp1(4)+tmp1(5))
           else
             if(cxI(1)>3.and.cxI(2)>3.and.cxI(3)>3)then
                tmp2= func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3))
             else
                cxB=cxI-3
                cxT=cxI+4
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,8,Symmetry))then
                 write(*,*)"prolong3 position: "
                 write(*,*)llbf(1)+(i-1)*FD(1),llbf(2)+(j-1)*FD(2),llbf(3)+(k-1)*FD(3)
                 write(*,*)"llbf = ",llbf
                 stop
                endif
                tmp2= ya(:,:,4)
             endif
             tmp1= tmp2(:,4)
             funf(i,j,k)= C1*(tmp1(1)+tmp1(8))+C2*(tmp1(2)+tmp1(7))+C3*(tmp1(3)+tmp1(6))+C4*(tmp1(4)+tmp1(5))
           endif
         endif
       else
         if(jj/2*2==jj)then
           if(kk/2*2==kk)then               
             if(cxI(1)>3.and.cxI(2)>3.and.cxI(3)>3)then
               tmp2= C1*(func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)-3)+func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+4))&
                    +C2*(func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)-2)+func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+3))&
                    +C3*(func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)-1)+func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+2))&
                    +C4*(func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)  )+func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+1))
             else
                cxB=cxI-3
                cxT=cxI+4
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,8,Symmetry))then
                 write(*,*)"prolong3 position: "
                 write(*,*)llbf(1)+(i-1)*FD(1),llbf(2)+(j-1)*FD(2),llbf(3)+(k-1)*FD(3)
                 write(*,*)"llbf = ",llbf
                 stop
                endif
                tmp2=C1*(ya(:,:,1)+ya(:,:,8))+C2*(ya(:,:,2)+ya(:,:,7))+C3*(ya(:,:,3)+ya(:,:,6))+C4*(ya(:,:,4)+ya(:,:,5))
             endif
             tmp1= C1*(tmp2(:,1)+tmp2(:,8))+C2*(tmp2(:,2)+tmp2(:,7))+C3*(tmp2(:,3)+tmp2(:,6))+C4*(tmp2(:,4)+tmp2(:,5))
             funf(i,j,k)= tmp1(4)
           else
             if(cxI(1)>3.and.cxI(2)>3.and.cxI(3)>3)then
                tmp2= func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3))
             else
                cxB=cxI-3
                cxT=cxI+4
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,8,Symmetry))then
                 write(*,*)"prolong3 position: "
                 write(*,*)llbf(1)+(i-1)*FD(1),llbf(2)+(j-1)*FD(2),llbf(3)+(k-1)*FD(3)
                 write(*,*)"llbf = ",llbf
                 stop
                endif
                tmp2= ya(:,:,4)
             endif
             tmp1= C1*(tmp2(:,1)+tmp2(:,8))+C2*(tmp2(:,2)+tmp2(:,7))+C3*(tmp2(:,3)+tmp2(:,6))+C4*(tmp2(:,4)+tmp2(:,5))
             funf(i,j,k)= tmp1(4)
           endif
         else
           if(kk/2*2==kk)then
             if(cxI(1)>3.and.cxI(2)>3.and.cxI(3)>3)then
               tmp2= C1*(func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)-3)+func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+4))&
                    +C2*(func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)-2)+func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+3))&
                    +C3*(func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)-1)+func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+2))&
                    +C4*(func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)  )+func(cxI(1)-3:cxI(1)+4,cxI(2)-3:cxI(2)+4,cxI(3)+1))
             else
                cxB=cxI-3
                cxT=cxI+4
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,8,Symmetry))then
                 write(*,*)"prolong3 position: "
                 write(*,*)llbf(1)+(i-1)*FD(1),llbf(2)+(j-1)*FD(2),llbf(3)+(k-1)*FD(3)
                 write(*,*)"llbf = ",llbf
                 stop
                endif
                tmp2=C1*(ya(:,:,1)+ya(:,:,8))+C2*(ya(:,:,2)+ya(:,:,7))+C3*(ya(:,:,3)+ya(:,:,6))+C4*(ya(:,:,4)+ya(:,:,5))
             endif
             tmp1= tmp2(:,4)
             funf(i,j,k)= tmp1(4)
           else
             funf(i,j,k)= func(cxI(1),cxI(2),cxI(3))
           endif
         endif
       endif
    enddo
   enddo
  enddo

  return

  end subroutine prolong3

#elif (ghost_width == 5)
! eighth order code
!--------------------------------------------------------------------------
!
! Prolongation from coarser grids to finer grids
! 10 points, 9th order interpolation
! 1   2   3   4   5   6   7   8   9   10
! *---*---*---*---*---*---*---*---*---*
!                   ^
! f=35/65536(f_1+f_10)-405/65536*(f_2+f_9) + 567/16384*(f_3+f_8) - 2205/16384*(f_4+f_7) + 19845/32768*(f_5+f_6)
!--------------------------------------------------------------------------

  subroutine prolong3(wei,llbc,uubc,extc,func,&
                      llbf,uubf,extf,funf,&
                      llbp,uubp,SoA,Symmetry)
  implicit none

!~~~~~~> input arguments
  integer,intent(in) :: wei
!                                       coarse    fine       coarse
  real*8,dimension(3),   intent(in) :: llbc,uubc,llbf,uubf,llbp,uubp
  integer,dimension(3),  intent(in) :: extc,extf
  real*8, dimension(extc(1),extc(2),extc(3)),intent(in)   :: func
  real*8, dimension(extf(1),extf(2),extf(3)),intent(inout):: funf
  real*8, dimension(1:3), intent(in) :: SoA
  integer,intent(in)::Symmetry

!~~~~~~> local variables

  real*8, dimension(1:3) :: base
  integer,dimension(3) :: lbc,ubc,lbf,ubf,lbp,ubp,lbpc,ubpc
  integer,dimension(3) :: cxB,cxT,cxI
  integer :: i,j,k,ii,jj,kk
  real*8,  dimension(10,10,10) :: ya
  real*8, dimension(10,10) :: tmp2
  real*8, dimension(10) :: tmp1

  real*8, parameter :: C1=3.5d1/6.5536d4,C2=-4.05d2/6.5536d4,C3=5.67d2/1.6384d4
  real*8, parameter :: C4=-2.205d3/1.6384d4,C5=1.9845d4/3.2768d4

  integer::imini,imaxi,jmini,jmaxi,kmini,kmaxi
  integer::imino,imaxo,jmino,jmaxo,kmino,kmaxo
  logical::decide3d

  real*8,dimension(3) :: CD,FD

  if(wei.ne.3)then
     write(*,*)"prolongrestrict.f90::prolong3: this routine only surport 3 dimension"
     write(*,*)"dim = ",wei
     stop
  endif

! it's possible a iolated point for target but not for source
  CD = (uubc-llbc)/(extc-1)
  FD = CD/2

!take care the mismatch of the two segments of grid
  do i=1,3
     if(llbc(i) <= llbf(i))then
        base(i) = llbc(i)
     else
        j=idint((llbc(i)-llbf(i))/FD(i)+0.4)
        if(j/2*2 == j)then
           base(i) = llbf(i)
        else
           base(i) = llbf(i) - CD(i)/2
        endif
     endif
  enddo

!!! function idint:
!If A is of type REAL and |A| < 1, INT(A) equals 0. If |A| \geq 1, 
!then INT(A) equals the largest integer that does not exceed the range of A 
!and whose sign is the same as the sign of A.

    lbf = idint((llbf-base)/FD+0.4)+1
    ubf = idint((uubf-base)/FD+0.4)+1
    lbc = idint((llbc-base)/CD+0.4)+1
    ubc = idint((uubc-base)/CD+0.4)+1
    lbp = idint((llbp-base)/FD+0.4)+1
    lbpc = idint((llbp-base)/CD+0.4)+1
    ubp = idint((uubp-base)/FD+0.4)+1
    ubpc = idint((uubp-base)/CD+0.4)+1
!sanity check
  imino=lbp(1)-lbf(1) + 1
  imaxo=ubp(1)-lbf(1) + 1
  jmino=lbp(2)-lbf(2) + 1
  jmaxo=ubp(2)-lbf(2) + 1
  kmino=lbp(3)-lbf(3) + 1
  kmaxo=ubp(3)-lbf(3) + 1

  imini=lbpc(1)-lbc(1) + 1
  imaxi=ubpc(1)-lbc(1) + 1
  jmini=lbpc(2)-lbc(2) + 1
  jmaxi=ubpc(2)-lbc(2) + 1
  kmini=lbpc(3)-lbc(3) + 1
  kmaxi=ubpc(3)-lbc(3) + 1

  if(imino.lt.1.or.jmino.lt.1.or.kmino.lt.1.or.&
     imini.lt.1.or.jmini.lt.1.or.kmini.lt.1.or.&
     imaxo.gt.extf(1).or.jmaxo.gt.extf(2).or.kmaxo.gt.extf(3).or.&
     imaxi.gt.extc(1)-4.or.jmaxi.gt.extc(2)-4.or.kmaxi.gt.extc(3)-4)then
          write(*,*)"error in prolongation for"
          write(*,*)"from"
          write(*,*)llbc,uubc
          write(*,*)lbc,ubc
          write(*,*)"to"
          write(*,*)llbf,uubf
          write(*,*)lbf,ubf
          write(*,*)"want"
          write(*,*)llbp,uubp
          write(*,*)lbp,ubp,lbpc,ubpc
          if(imini.lt.1) write(*,*)"imini = ",imini
          if(jmini.lt.1) write(*,*)"jmini = ",jmini       
          if(kmini.lt.1) write(*,*)"kmini = ",kmini      
          if(imino.lt.1) write(*,*)"imino = ",imino       
          if(jmino.lt.1) write(*,*)"jmino = ",jmino       
          if(kmino.lt.1) write(*,*)"kmino = ",kmino   
          if(imaxi.gt.extc(1)) write(*,*)"imaxi = ",imaxi,"extc(1) = ",extc(1) 
          if(jmaxi.gt.extc(2)) write(*,*)"jmaxi = ",jmaxi,"extc(2) = ",extc(2) 
          if(kmaxi.gt.extc(3)) write(*,*)"kmaxi = ",kmaxi,"extc(3) = ",extc(3) 
          if(imaxo.gt.extf(1)) write(*,*)"imaxo = ",imaxo,"extf(1) = ",extf(1) 
          if(jmaxo.gt.extf(2)) write(*,*)"jmaxo = ",jmaxo,"extf(2) = ",extf(2) 
          if(kmaxo.gt.extf(3)) write(*,*)"kmaxo = ",kmaxo,"extf(3) = ",extf(3) 
          return
  endif
!~~~~~~> prolongation start...
  do k = kmino,kmaxo
   do j = jmino,jmaxo
    do i = imino,imaxo
! change to coarse level reference
!|*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*--- ---*| 
!|x===============x===============x===============x========|
       cxI(1) = i
       cxI(2) = j
       cxI(3) = k
       cxI = (cxI+lbf)/2
! change to array index      
       cxI = cxI - lbc + 1

       ii=i+lbf(1)-1
       jj=j+lbf(2)-1
       kk=k+lbf(3)-1

       if(any(cxI+5 > extc)) write(*,*)"error in prolong"
       if(ii/2*2==ii)then
         if(jj/2*2==jj)then
           if(kk/2*2==kk)then
!  due to ghost zone, we can deal with symmetry boundary like this                                   
             if(cxI(1)>4.and.cxI(2)>4.and.cxI(3)>4)then
               tmp2= C1*(func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-4)+func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+5))&
                    +C2*(func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-3)+func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+4))&
                    +C3*(func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-2)+func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+3))&
                    +C4*(func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-1)+func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+2))&
                    +C5*(func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)  )+func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+1))
             else
               cxB=cxI-4
               cxT=cxI+5
               if(decide3d(extc,func,func,cxB,cxT,SoA,ya,10,Symmetry))then
                write(*,*)"prolong3 position index: ",i+lbf(1)-1,j+lbf(2)-1,k+lbf(3)-1
                stop
               endif
               tmp2=C1*(ya(:,:,1)+ya(:,:,10))+C2*(ya(:,:,2)+ya(:,:,9))+C3*(ya(:,:,3)+ya(:,:,8)) &
                   +C4*(ya(:,:,4)+ya(:,:, 7))+C5*(ya(:,:,5)+ya(:,:,6))
             endif
             tmp1= C1*(tmp2(:,1)+tmp2(:,10))+C2*(tmp2(:,2)+tmp2(:,9))+C3*(tmp2(:,3)+tmp2(:,8)) &
                  +C4*(tmp2(:,4)+tmp2(:, 7))+C5*(tmp2(:,5)+tmp2(:,6))
             funf(i,j,k)= C1*(tmp1(1)+tmp1(10))+C2*(tmp1(2)+tmp1(9))+C3*(tmp1(3)+tmp1(8)) &
                         +C4*(tmp1(4)+tmp1( 7))+C5*(tmp1(5)+tmp1(6))
           else
             if(cxI(1)>4.and.cxI(2)>4.and.cxI(3)>4)then
                tmp2= func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3))
             else
                cxB=cxI-4
                cxT=cxI+5
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,10,Symmetry))then
                        write(*,*)"prolong3 position index: ",i+lbf(1)-1,j+lbf(2)-1,k+lbf(3)-1
                        return
                endif
                tmp2= ya(:,:,5)
             endif
             tmp1= C1*(tmp2(:,1)+tmp2(:,10))+C2*(tmp2(:,2)+tmp2(:,9))+C3*(tmp2(:,3)+tmp2(:,8)) &
                  +C4*(tmp2(:,4)+tmp2(:, 7))+C5*(tmp2(:,5)+tmp2(:,6))
             funf(i,j,k)= C1*(tmp1(1)+tmp1(10))+C2*(tmp1(2)+tmp1(9))+C3*(tmp1(3)+tmp1(8)) &
                         +C4*(tmp1(4)+tmp1( 7))+C5*(tmp1(5)+tmp1(6))
           endif
         else
           if(kk/2*2==kk)then
             if(cxI(1)>4.and.cxI(2)>4.and.cxI(3)>4)then
               tmp2= C1*(func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-4)+func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+5))&
                    +C2*(func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-3)+func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+4))&
                    +C3*(func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-2)+func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+3))&
                    +C4*(func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-1)+func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+2))&
                    +C5*(func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)  )+func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+1))
             else
                cxB=cxI-4
                cxT=cxI+5
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,10,Symmetry))then
                        write(*,*)"prolong3 position index: ",i+lbf(1)-1,j+lbf(2)-1,k+lbf(3)-1
                        return
                endif
               tmp2=C1*(ya(:,:,1)+ya(:,:,10))+C2*(ya(:,:,2)+ya(:,:,9))+C3*(ya(:,:,3)+ya(:,:,8)) &
                   +C4*(ya(:,:,4)+ya(:,:, 7))+C5*(ya(:,:,5)+ya(:,:,6))
             endif
             tmp1= tmp2(:,5)
             funf(i,j,k)= C1*(tmp1(1)+tmp1(10))+C2*(tmp1(2)+tmp1(9))+C3*(tmp1(3)+tmp1(8)) &
                         +C4*(tmp1(4)+tmp1( 7))+C5*(tmp1(5)+tmp1(6))
           else
             if(cxI(1)>4.and.cxI(2)>4.and.cxI(3)>4)then
                tmp2= func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3))
             else
                cxB=cxI-4
                cxT=cxI+5
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,10,Symmetry))then
                        write(*,*)"prolong3 position index: ",i+lbf(1)-1,j+lbf(2)-1,k+lbf(3)-1
                        return
                endif
                tmp2= ya(:,:,5)
             endif
             tmp1= tmp2(:,5)
             funf(i,j,k)= C1*(tmp1(1)+tmp1(10))+C2*(tmp1(2)+tmp1(9))+C3*(tmp1(3)+tmp1(8)) &
                         +C4*(tmp1(4)+tmp1( 7))+C5*(tmp1(5)+tmp1(6))
           endif
         endif
       else
         if(jj/2*2==jj)then
           if(kk/2*2==kk)then               
             if(cxI(1)>4.and.cxI(2)>4.and.cxI(3)>4)then
               tmp2= C1*(func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-4)+func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+5))&
                    +C2*(func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-3)+func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+4))&
                    +C3*(func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-2)+func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+3))&
                    +C4*(func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-1)+func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+2))&
                    +C5*(func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)  )+func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+1))
             else
                cxB=cxI-4
                cxT=cxI+5
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,10,Symmetry))then
                        write(*,*)"prolong3 position index: ",i+lbf(1)-1,j+lbf(2)-1,k+lbf(3)-1
                        return
                endif
               tmp2=C1*(ya(:,:,1)+ya(:,:,10))+C2*(ya(:,:,2)+ya(:,:,9))+C3*(ya(:,:,3)+ya(:,:,8)) &
                   +C4*(ya(:,:,4)+ya(:,:, 7))+C5*(ya(:,:,5)+ya(:,:,6))
             endif
             tmp1= C1*(tmp2(:,1)+tmp2(:,10))+C2*(tmp2(:,2)+tmp2(:,9))+C3*(tmp2(:,3)+tmp2(:,8)) &
                  +C4*(tmp2(:,4)+tmp2(:, 7))+C5*(tmp2(:,5)+tmp2(:,6))
             funf(i,j,k)= tmp1(5)
           else
             if(cxI(1)>4.and.cxI(2)>4.and.cxI(3)>4)then
                tmp2= func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3))
             else
                cxB=cxI-4
                cxT=cxI+5
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,10,Symmetry))then
                        write(*,*)"prolong3 position index: ",i+lbf(1)-1,j+lbf(2)-1,k+lbf(3)-1
                        return
                endif
                tmp2= ya(:,:,5)
             endif
             tmp1= C1*(tmp2(:,1)+tmp2(:,10))+C2*(tmp2(:,2)+tmp2(:,9))+C3*(tmp2(:,3)+tmp2(:,8)) &
                  +C4*(tmp2(:,4)+tmp2(:, 7))+C5*(tmp2(:,5)+tmp2(:,6))
             funf(i,j,k)= tmp1(5)
           endif
         else
           if(kk/2*2==kk)then
             if(cxI(1)>4.and.cxI(2)>4.and.cxI(3)>4)then
               tmp2= C1*(func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-4)+func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+5))&
                    +C2*(func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-3)+func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+4))&
                    +C3*(func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-2)+func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+3))&
                    +C4*(func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)-1)+func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+2))&
                    +C5*(func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)  )+func(cxI(1)-4:cxI(1)+5,cxI(2)-4:cxI(2)+5,cxI(3)+1))
             else
                cxB=cxI-4
                cxT=cxI+5
                if(decide3d(extc,func,func,cxB,cxT,SoA,ya,10,Symmetry))then
                        write(*,*)"prolong3 position index: ",i+lbf(1)-1,j+lbf(2)-1,k+lbf(3)-1
                        return
                endif
               tmp2=C1*(ya(:,:,1)+ya(:,:,10))+C2*(ya(:,:,2)+ya(:,:,9))+C3*(ya(:,:,3)+ya(:,:,8)) &
                   +C4*(ya(:,:,4)+ya(:,:, 7))+C5*(ya(:,:,5)+ya(:,:,6))
             endif
             tmp1= tmp2(:,5)
             funf(i,j,k)= tmp1(5)
           else
             funf(i,j,k)= func(cxI(1),cxI(2),cxI(3))
           endif
         endif
       endif
    enddo
   enddo
  enddo

  return

  end subroutine prolong3
#endif

#else
#ifndef Cell
#error Not define Vertex nor Cell
#endif
#endif

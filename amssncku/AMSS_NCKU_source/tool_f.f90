subroutine fout3D(ftag,matrix,msize)
  
  implicit none
  integer,intent(in ):: ftag,msize(1:3)
  double precision,intent(in),dimension(msize(1),msize(2),msize(3))::matrix
  integer mat_size
  mat_size = msize(1)*msize(2)*msize(3)
  
  call writefile_f(ftag,matrix,mat_size)
  return

end subroutine fout3D

subroutine fout1D(ftag,matrix,msize)
  
  implicit none
  integer,intent(in ):: ftag,msize
  double precision,intent(in),dimension(1:msize)::matrix
 
  call writefile_f(ftag,matrix,msize)
  return

end subroutine fout1D

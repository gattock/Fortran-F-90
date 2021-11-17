program heatgaussf90
implicit none

integer :: i,j,k
integer, parameter :: a=40, b=60
real, parameter :: x=4, y=6, dx=x/a, dy=y/b

real*4 :: c(a*b, a*b+1 ) !linear system with external condition column(+1)

! i,j are the indexes for the elements of the metal sheet 
! the indexes for the matrix containing all the coefficient are p,k
! is possible to obtain p=(i-1)*b+j 
! k is the index for all the equation of the matrix c

!generate matrix of zeroes of size a*b x a*b+1
do k=1,a*b
  do p=1,a*b+1
c(k,p)=0
  end do
end do

!boundary 300°K conditions:
  do i=2,(a-1)
	c( (i-1)*b+1 , (i-1)*b+1 ) =1 !j=1
	c( (i-1)*b+1 , a*b+1 ) =300 !j=1
	c( i*b , i*b )=1 	      !j=b
	c( i*b , a*b+1)=300 	      !j=b
  end do
!boundary 500°K conditions:
  do j=2, (b-1)
	c( j, j )=1 			!i=1
	c( j, a*b+1)=500		!i=1
	c( (a-1)*b+j , (a-1)*b+j )= 1	!i=a
	c( (a-1)*b+j , a*b+1) =500	!i=a
  end do

!imposition of T steady state of not boundary elements:
! T(i,j) - 1/4*( T(i+1,j)+T(i-1,j)+T(i,j+1)+T(i,j-1))=0

  do i=2, (a-1)
	do j=2, (b-1)
	c( (i-1)*b+j , (i-1)*b+j )=1
	c( (i-1)*b+j , (i-1)*b+j-1 )=-0.25
	c( (i-1)*b+j , (i-1)*b+j+1 )=-0.25
	c( (i-1)*b+j , (i-2)*b+j )=-0.25
	c( (i-1)*b+j , i*b+j )=-0.25
	end do
  end do

!Gauss-Jordan elimination 
!working on this,programming the subroutine is taking more time than expected

!Reading values of the T(i,j) from the last column of the matrix g=c(gauss-reduced)
  do i=2,(a-1)
	do j=2, (b-1)
	T(i,j)= c( (i-1)*b+j , a*b+1)
	end do
  end do

!plot with python
open(10, file = 'out.dat')

do i=2,(a-1)
	do j=2, (b-1)
	write(10, *) T(i,j)
	end do
  end do
  close(10)

end program heatgaussf90

program heatgaussf90
implicit none

integer :: i,j,k
integer, parameter :: a=20, b=10
real, parameter :: x=4, y=6, dx=x/a, dy=y/b

real*4 :: c(a*b, a*b+1 ) !linear system with external condition column(+1)
real*4 :: G(a*b, a*b+1 )
real*4 :: T(a,b)

! i,j are the indexes for the elements of the metal sheet 
! the indexes for the matrix containing all the coefficient are p,k
! is possible to obtain p=(i-1)*b+j 
! k is the index for all the equation of the matrix c

!generate matrix of zeroes of size a*b x a*b+1
do i=1,a*b
  do j=1,a*b+1
c(i,j)=0
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

	c( (i-1)*b+j , (i-1)*b+j-1 )=-0.25 !check formula
	c( (i-1)*b+j , (i-1)*b+j+1 )=-0.25

	c( (i-1)*b+j , (i-2)*b+j )=-0.25
	c( (i-1)*b+j , i*b+j )=-0.25
	end do
  end do

!Gauss-Jordan elimination 
call gauss_sub (a*b,a*b+1,c,G)
call jordan_sub (a*b,a*b+1,G,c)

!Reading values of the T(i,j) from the last column of the matrix g=c(gauss-reduced)
  do i=2,(a-1)
	do j=2, (b-1)
	T(i,j)= c( (i-1)*b+j , a*b+1)
	end do
  end do

call print_matrix (a,b,T)

end program heatgaussf90




!_______________________________________
SUBROUTINE gauss_sub (a,b,LS,G)
implicit none
integer,intent(in) :: a,b !a rows, b columns
real, intent(in) :: LS(a,b) !LinSist,
real, intent(out) :: G(a,b) !gaussj.ed
real :: Gin(a,b) !gauss-input, 
integer :: i
		call matrix_substitution(a,b,LS,Gin)
	do i=1,a
	call pivot_the_rows (a,b,Gin,i,G)
	call matrix_substitution(a,b,G,Gin)
	call one_pivot_survival (a,b,Gin,i,G)
	call matrix_substitution(a,b,G,Gin)
	end do
END SUBROUTINE gauss_sub
!________________________________________________
SUBROUTINE jordan_sub (nr,nc,G,GJout)
implicit none
integer,intent(in) :: nr,nc
real, intent(in) :: G(nr,nc) !in this sub G will not be G after 1 i-round
real, intent(out) :: GJout(nr,nc)
integer :: i,j
		do i=1,nr-1
call row_jordaning (nr,nc,i,G,GJout)
call matrix_substitution(nr,nc,GJout,G)
		end do
END SUBROUTINE jordan_sub
!________________________________________________
SUBROUTINE row_jordaning(nr,nc,i,G,GJout)
implicit none
integer,intent(in) :: nr,nc,i
real, intent(in) :: G(nr,nc)
real, intent(out) :: GJout(nr,nc)
integer :: j
		call matrix_substitution(nr,nc,G,GJout)
GJout(nr-i,nc)=G(nr-i,nc)-G(nr,nc)*G(nr-i,nc-1)
	do j=1,i
GJout(nr-i,nc-j)=G(nr-i,nc-j)-G(nr-j+1,nc-j)*G(nr-i,nc-j)
	call matrix_substitution(nr,nc,GJout,G)
	end do
END SUBROUTINE row_jordaning
!________________________________________________
SUBROUTINE print_matrix (nr,nc,mat)
implicit none
integer,intent(in) :: nr,nc
real, intent(in) :: mat(nr,nc)
integer :: i,j
do i=1,nr
	!print *, (mat(i,j), j=1,nc )
	write(*,80) (mat(i,j),j=1,nc )
	80 FORMAT ('',10F6.1)
end do
print *, '-------------------------------------------'
END SUBROUTINE print_matrix
!_______________________________________________
SUBROUTINE matrix_substitution (nr,nc,A,Aupdate)
implicit none
integer,intent(in) :: nr,nc
real, intent(in) :: A(nr,nc) 
real :: Aupdate(nr,nc) 
integer :: i,j
	do i=1,nr
	do j=1,nc
Aupdate(i,j)=A(i,j)
	end do
	end do
END SUBROUTINE matrix_substitution
!_______________________________________________
SUBROUTINE pivot_the_rows (nr,nc,matrix,index,newMat)
implicit none
integer,intent(in) :: nr,nc,index
real, intent(in) :: matrix(nr,nc)
real, intent(out) :: newMat(nr,nc)
integer :: j,k
		do j=index,nr
	do k=index,nc
newMat(j,k)=matrix(j,k)/matrix(j,index)
	end do
		end do
END SUBROUTINE pivot_the_rows
!_______________________________________________
SUBROUTINE one_pivot_survival (nr,nc,matrix,index,newMat)
implicit none
integer,intent(in) :: nr,nc,index
real, intent(in) :: matrix(nr,nc)
real,intent(out) :: newMat(nr,nc)
integer :: j,k
		do j=index+1,nr
	do k=index,nc
newMat(j,k)=matrix(j,k)-matrix(index,k)
	end do
		end do
END SUBROUTINE one_pivot_survival
!_______________________________________________
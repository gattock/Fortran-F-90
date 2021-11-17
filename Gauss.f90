program gauss 
implicit none

integer, parameter :: a=5,b=6 !a rows, b columns
real :: LS(a,b),Gin(a,b),G(a,b) !LinSist, gauss-input, gaussed
integer :: i,h
		call rand_matrix(a,b,LS)
		call matrix_substitution(a,b,LS,Gin)
	do i=1,a
call print_matrix(a,b,Gin)
	call pivot_the_rows (a,b,Gin,i,G)
	call matrix_substitution(a,b,G,Gin)
call print_matrix(a,b,Gin)
	call one_pivot_survival (a,b,Gin,i,G)
	call matrix_substitution(a,b,G,Gin)
	end do
call print_matrix(a,b,Gin)
end program
!________________________________________________
SUBROUTINE rand_matrix (nr,nc,mat)
implicit none
integer, intent(in) :: nr,nc
real ,intent(out) :: mat(nr,nc)
integer :: i,j
real :: ciccio
	do i=1,nr	!filling matrix with random integers
		do j=1,nc
call random_number(ciccio)
mat(i,j)=int(ciccio*10)
		end do
	end do
END SUBROUTINE rand_matrix
!_______________________________________________
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
real, intent(out) :: Aupdate(nr,nc)
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
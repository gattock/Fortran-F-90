program prova_gauss
implicit none
integer, parameter :: imax=4,jmax=5
real :: LinSis(imax,jmax), gaussed(imax,jmax)

call rand_matrix(imax,jmax,LinSis)
call gauss_sub (imax,jmax,LinSis,gaussed)
call print_matrix (imax,jmax,gaussed)

end program
!_______________________________________
SUBROUTINE gauss_sub (a,b,LS,G)
implicit none
integer,intent(in) :: a,b !a rows, b columns
real, intent(in) :: LS(a,b) !LinSist,
real, intent(out) :: G(a,b) !gaussed
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

!this one must be ereased if LinSis is coming from real problem

SUBROUTINE rand_matrix (nr,nc,mat)
implicit none
integer, intent(in) :: nr,nc
real ,intent(out) :: mat(nr,nc)
integer :: i,j
real :: ciccio
	do i=1,nr	!filling matrix with random integers
		do j=1,nc
call random_number(ciccio)
mat(i,j)=int(ciccio*10)+1
		end do
	end do
END SUBROUTINE rand_matrix
!_______________________________________________
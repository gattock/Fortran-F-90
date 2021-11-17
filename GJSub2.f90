program prova_gauss_jordan
implicit none
integer, parameter :: imax=3,jmax=4
real :: LinSis(imax,jmax), G(imax,jmax), GJ(imax,jmax)
call rand_matrix(imax,jmax,LinSis)
call print_matrix(imax,jmax,LinSis)
call gauss_sub (imax,jmax,LinSis,G)
call print_matrix(imax,jmax,G)
call jordan_sub (imax,jmax,G,GJ)
call print_matrix (imax,jmax,GJ)
end program
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
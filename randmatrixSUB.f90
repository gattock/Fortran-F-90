!call rand_matrix(a,b,C)
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
mat(i,j)=int(ciccio*10)+1
		end do
	end do
END SUBROUTINE rand_matrix

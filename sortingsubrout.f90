program prova
implicit none
integer, parameter :: a=20
integer :: i
real*4 :: v(a),numero,w(a)

do i=1,a
	call random_number(numero)
	v(i)=int(numero*10)
end do

call sorting_vector (a,v,w)

end program

!_______________________________________________
SUBROUTINE sorting_vector (size,v,w)
implicit none

integer,intent(in) :: size
integer :: temp,i,j,check
real*4 :: v(size),numero
real*4, intent(out) :: w(size)
	
		do i=1,size
check=0
	do j=i,size
if (v(i)<v(j)) then
	temp=v(i)
	v(i)=v(j)
	v(j)=temp
end if

if (w(j)/=v(j)) then
check=1
end if
	end do
do j=1,size
w(j)=v(j)
end do

if (check==0) then
stop
end if	


	write(*,90) v
	90 FORMAT ('v=',20F4.0) !update: aF8.0 !!
	write(*,100) w
	100 FORMAT ('w=',20F4.0) !update: aF8.0 !!
	print *, '________________________'

		end do
END SUBROUTINE sorting_vector
!_______________________________________________
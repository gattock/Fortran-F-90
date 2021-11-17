program beam_prog
  implicit none

  integer :: i, j
  real*4    :: d
  real*4, parameter :: dx = 0.01, L = 10.
  integer, parameter   :: n = int(L/dx)
  real*4    :: t = 0.01
  real*4, parameter :: E = 200.e9
  real*4    :: beam(3,n) ! Beam parameters 1: positions, 2: height, 3: inertia
  real*4    :: q(n) ! Load 
  real*4    :: s(n), m(n) ! Shear, momentum
  real*4    :: chi(n) ! Curvature
  real*4    :: theta(n) ! Slope the beam
  real*4    :: y(n) ! Vertical displacement

!  n = int(L/dx) 
  beam(1,:) = 0.
  beam(2,:) = 0.1
  q(:) = 1.e2
  d = 0.

  do i = 1, n
    beam(1,i) = d
    beam(3,i) = t*beam(2,i)**3/12.
    d = d + dx
  end do

  print *, "Beam 1"
  print *, beam(1,:)
  
  do i = 1, n
    do j = i-1, n
      s(i) = s(i) + q(j) * (beam(1,j) - beam(1,j-1))
    end do
  end do
  do i = 1, n
    do j = i-1, n
      m(i) = m(i) + s(j) * (beam(1,j) - beam(1,j-1))
    end do
    chi(i) = m(i)/(E*beam(3,i))
  end do

  do i = 2, n
    do j = 2, i
      theta(i) = theta(i) + chi(j) * (beam(1,j) - beam(1,j-1))
    end do
  end do
!  theta(1)= theta(1) + chi(j) * (beam(1,j) - beam(1,j-1))

  print *, "Theta"
  print *, theta(:)

  do i = 1, n
    do j = 2, i
      y(i) = y(i) + theta(j) * (beam(1,j) - beam(1,j-1))
    end do
  end do

  open(10, file = 'out.dat')
  do i = 1, n
    write(10, *) beam(1,i), theta(i), y(i)
  end do
  close(10)

end program beam_prog




program heat
implicit none

integer :: i,j,k
integer, parameter :: a=40,b=60,cy=1000
real*4, parameter :: dt=0.005,tmax=cy*dt,x=4,y=6,dx=x/a,dy=y/b,q=10
!dt=time discretization -- tmax=simulation time -- cy=num cycles
!q=heat transmission coefficient a(T)
!a,b= element numbers on x,y -- x,y=side sizes [m]  --  dx,dy=element size [m]

real*4 :: T(a,b,cy) !temperature of i,j element in time t 
real*4 :: dTa(a,b,cy) ! dT/dx of i,j....
real*4 :: dTb(a,b,cy) ! dT/dy of i,j....
real*4 :: lap(a,b,cy) !laplacian of i,j...
real*4 :: dTt(a,b,cy) ! dT/dt of i,j...

!boundary temperature begin (assumed all plate 300°K except long borders at 500°K)
do j=1,b
	do i=2,a-1
	T(i,j,1)=300
        end do
T(1,j,1)=500
T(a,j,1)=500
end do !boundary temperature end


do k= 1,cy-1 !along timeeee
   
   !partial derivative matrix, Laplacian, calcs...
   do i=2,a-1 !along short side x-->a
		do j=2,b-1 !along long side y-->b
                dTa(i,j,k)= ( T(i+1,j,k)-T(i-1,j,k) )/(2*dx)
		dTb(i,j,k)= ( T(i,j+1,k)-T(i,j-1,k) )/(2*dy)
     lap(i,j,k)= (T(i+1,j,k)+T(i-1,j,k)+T(i,j+1,k)+T(i,j-1,k)-4*T(i,j,k)) / (dx*dy)
                end do
   end do !end ...calcs...

   !dT/dt
   do i=2,a-1
      do j=2,b-1
         dTt(i,j,k)=q*lap(i,j,k)
         T(i,j,k+1)=T(i,j,k)+dTt(i,j,k)*dt
      end do
   end do

   print *, "T°sample=",T(a/4,b/4,k)
   
end do !end along timeeeee

end program heat

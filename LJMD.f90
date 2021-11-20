program LennardF90

implicit none

real*4, parameter :: dt=0.01, eps=1, sig=0.1 ,side=1, tmax=0.05 ! dtime, 24*epsilon, sigma, box side, timemax
integer, parameter :: n=3, m=3		 ! n particles, m=3Dimensions
integer, parameter :: smax=tmax/dt +1
integer i,j,k,s
real*4 :: rsq,l,p,t

	!declaration of matrixes, parameters and variables
real*4 :: X(n,m,smax) !position tensor (n=which particle,m= 3 components, s= time ascissa)
real*4 :: A(n,m,smax) !are memorized only position and acceleration of each particle over the time, speed and force local variables->forgotten
real*4 :: VER(n,n,m) !versor direction of force
real*4 :: V(n,m)
real*4 :: F(n,n)
real*4 :: MASS(n) !masses vector

			print*, "smax=cyclesnumber=",smax
t=0
s=1	

	!fill X and V of random numbers
DO i= 1,n
DO j= 1,m
	call random_number (l) 
	call random_number (p)
	X(i,j,1)=l*side !random initial position
	V(i,j)=l*p/5	!small random inital velocity
END DO 
MASS(i)=1
END DO

		DO s=1,smax-1 !TIME
		
	!interactions and cinematics (rsq=r**2)
DO i=1,n-1	
	DO k=i+1,n
	rsq=0
		DO j=1,m
		rsq=rsq+(X(i,j,s)-X(k,j,s))**2		!r^2 between i(th) and k(th) particle
		END DO
		print*, "r^2=", rsq, "i,k=",i,k
		DO j=1,m
		VER(i,k,j)=(X(i,j,s)-X(k,j,s))/(rsq**0.5) !versor of force (direction)
		END DO
		print*, "VER(",i,k,")=", (VER(i,k,j),j=1,m)
		print*, "-----------------"
	F(i,k)=eps*(2*(sig**12)/(rsq**6.5)-(sig**6)/(rsq**3.5))		!modulus of force between i(th) and k(th) particle
	
	END DO
END DO
	DO i=1,n
			print*, "F(i,k)=" , (F(i,k),k=1,n)
	END DO
	
DO i=1,n
A(i,j,s)=0
	DO j=1,m
		DO k=i+1,n		
		A(i,j,s)=A(i,j,s)+F(i,k)*VER(i,k,j)/MASS(i)
		END DO
	V(i,j)=V(i,j)+A(i,j,s)*dt
	X(i,j,s+1)=X(i,j,s)+V(i,j)*dt
	
	!if OUT OF THE BOX, THEN elastic collision with boundaries
		IF (X(i,j,s+1)<=0 .OR. X(i,j,s+1)>=side) THEN
			V(i,j)=-V(i,j)
		END IF	
	END DO
		print*, "A(",i,")=(", (A(i,j,s),j=1,m)
		print*, "X(",i,")=(", (X(i,j,s),j=1,m)
		print*, "V(",i,")=(", (V(i,j),j=1,m)
END DO
	t=t+dt
print*, "----------------------------------"
print*, "t=",t
		END DO !time cycle


			!go to 101
!export
open(10,file='out.dat')

DO s=1,smax
	write(10,*) "s=",s
   	DO i=1,n
		
  
      write(10,*) "X(", i,")=", (X(i,j,s),j=1,m)
      write(10,*) "A(", i,")=", (A(i,j,s),j=1,m)
		
   END DO
write(10,*) "-----------------------------------"
END DO
close(10)

			!101 continue

END PROGRAM LennardF90

function readArray(n) result(c)
	implicit none 
	integer, intent(in) :: n
	real(8) c(n,n)
	integer :: i,j
	open(15, file = 'cond_100000.txt',status='old')
	do i = 1,n
			read(15,*) (c(i,j), j = 1,n)
	end do
	close(15)
end function readArray

function readB(n) result(c)
	implicit none 
	integer, intent(in) :: n
	real(8) c(n)
	integer :: i,j
	open(15, file = 'cond_100000.txt',status='old')
	do i = 1,(n+1)
		read(15,*) (c(j), j = 1,n )
	end do
	close(15)
end function readB



real(8) function normINF(B,n)
	implicit none
	integer, intent(in) :: n
	integer i
	real(8) B(n), y
	y = dabs(B(1))
	do i = 2 ,n
		if (dabs(B(i)) > y) then
		y = dabs(B(i))
		end if
	end do
	normINF = y
end function normINF

real(8) function MatrixNorm(C, n) 
	implicit none
	integer, intent(in) :: n 
	real(8) C(n, n), maximum, summa  
	integer i, j

	maximum = 0  
	do I = 1, n 
		summa = 0 
		do j = 1, n 
			summa = summa + dabs(C(i,j)) 
		end do 
		if (summa > maximum) then 
			maximum = summa 
		end if 
	end do 
	MatrixNorm = maximum 
end function

function findX(U,y,n)  result(x)
	implicit none
	integer,intent(in) :: n
	real(8) :: U(n,n), x(n), y(n),summa
	integer :: i, j
	summa = 0
	x(n) = y(n)/U(n,n)
	do i = 1, n-1
		summa = 0
		do j = 1 , i
			summa = summa + U(n-i,n-j+1)*x(n-j+1)
		end do
		x(n-i) = (y(n-i) - summa)/U(n-i,n-i)
	end do
end function findX


subroutine rotate(A,B,C,B_C,n)
	implicit none
	integer :: n,j,i,k
	real(8) :: A(n,n), B(n), C(n,n), B_C(n), c_1, s_1, Y(n)
	do i = 1,n
		do j = 1,n
			C(i,j) = A(i,j)
		end do
	end do
	do j = 1,n
		B_C(j) = B(j)
	end do
	do k = 1, n-1
		do j = k+1, n
			c_1 = C(k,k)/(sqrt(C(k,k)**2+C(j,k)**2))
			s_1 = C(j,k)/(sqrt(C(k,k)**2+C(j,k)**2))
			do i = k , n
			Y(i) = C(k,i)
			C(k,i) = Y(i)*c_1 + s_1*C(j,i)
			C(j,i) = -s_1*Y(i) + c_1*C(j,i)
			end do
			Y(k) = B_C(k)
			B_C(k) = c_1*Y(k) + s_1*B_C(j)
			B_C(j) = -s_1*Y(k)+c_1*B_C(j)
		end do
	end do
end subroutine


subroutine normNEV(A,B, X, n, cond)
implicit none
integer, intent(in) :: n
real(8) :: X(n),X2(n),dX(n),A(n,n),A2(n,n),dA(n,n),C2(n,n), B(n), B_C2(n),normINF, cond, norm_dx, norm_dA, MatrixNorm, eq1,eq2
integer :: i,j
interface
	function findX(U,y,n)  result(x)
	implicit none
	integer,intent(in) :: n
	real(8) :: U(n,n), x(n), y(n)
	end function findX
end interface
do i = 1, N
	do j = 1,n 
	 	 call random_number(dA(i,j)) 
	 	 dA(i,j) = A(i,j)*(2 - dA(i,j))/100
	 	 A2(i,j) = A(i,j) + dA(i,j)
	end do
end do
call rotate(A2, B, C2, B_C2, n)
X2 = findX(C2, B_C2, n)
do i = 1, n
	dX(i) = dabs(X2(i) - X(i))
end do
norm_dx = normINF(dX,n)
norm_dA = MatrixNorm(dA, n)
eq1 = cond*norm_dA/(MatrixNorm(A, n))
eq2 = 1-cond*norm_dA/(MatrixNorm(A, n))
print*, eq2
if (norm_dx/(normINF(X,n)) .LE. eq1/eq2) then
	write(*, '(a   e32.16  a   e32.16)') ' TRUE: ', norm_dx/(normINF(X,n)), '  <=  ', eq1/eq2
else
	write(*, '(a   e32.16  a   e32.16)') ' FALSE: ', norm_dx/(normINF(X,n)), ' not <=  ', eq1/eq2
end if
end subroutine normNEV
	
	
	
program test
	implicit none
	integer :: n,i,j,k
	real(8), allocatable :: A(:,:), B(:), C(:,:), B_C(:), X(:)
	real(8) :: cond, summa,proma
	interface
		function readArray(n) result(c)
		implicit none
		integer, intent(in) :: n
		real(8) c(n,n)
		end function readArray
		
		function readB(n) result(c)
		implicit none 
		integer, intent(in) :: n
		real(8) c(n)
		end function readB
		
		function findX(U,y,n)  result(x)
		implicit none
		integer,intent(in) :: n
		real(8) :: U(n,n), x(n), y(n)
		end function findX
	end interface
	write(*,*) ' Enter N '
	read(*,*) n
	allocate(A(n,n))
	allocate(B(n))
	allocate(X(n))
	allocate(B_C(n))
	allocate(C(n,n))
	cond =   193854.37984180983
	A = readArray(n)
	B = readB(n)
	call rotate(A,B,C,B_C,n)
	X = findX(C, B_C, n)
	write(*,*) (X(j), j = 1,n)
	!call normNEV(A,B, X, n, cond)
	open(20, file = '5_M_solve_3.txt',status='replace')
	write(20,*) (x(j), j = 1,n)
	close(20)
	deallocate(A)
	deallocate(B)
	deallocate(C)
	deallocate(B_C)
	deallocate(X)
end program
	

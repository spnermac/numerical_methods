function dot(H,G,n) result(M)
	implicit none
	integer, intent(in) :: n
	real(8) :: H(n,n), G(n,n), M(n,n), summ
	integer :: i,j,k
	summ = 0
	do i = 1, n
		do j = 1,n
			do k = 1, n
			summ = summ + H(i,k)*G(k,j)
			end do
			M(i,j) = summ
			summ = 0
		end do
	end do
end function dot

function transpos(U,n) result(Ut)
	implicit none
	integer, intent(in) :: n
	real(8) :: U(n,n),Ut(n,n)
	integer :: i,j
	do i = 1, n
		do j = 1,n
		Ut(i,j) = U(j,i)
		end do
	end do
end function transpos
	
function readArray(n) result(c)
	implicit none 
	integer, intent(in) :: n
	real(8) c(n,n)
	integer :: i,j
	open(15, file = '5matrix_cond424812.txt',status='old')
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
	open(15, file = '5matrix_cond424812.txt',status='old')
	do i = 1,(n+1)
		read(15,*) (c(j), j = 1,n )
	end do
	close(15)
end function readB

function findU(A,n) result(U)
	implicit none
	integer, intent(in) :: n
	real(8) :: U(n,n), A(n,n), summa_sq,summa_2,prop, prop2, R
	integer :: i,j,t,t2
	summa_sq = 0
	summa_2 = 0
	prop = 0
	do i = 1,n
		do j = i, n 
			if (j .EQ. i) then
				if (i .NE. 1) then
					do t = 1, i-1
						R = U(t,i)**2
						summa_sq = summa_sq + R
						prop = summa_sq
					end do
				else
						prop = 0
				end if
				U(i,j) = sqrt(A(i,j)-prop)
				summa_sq = 0
			else
			if (j .NE. 1) then
				do t2 = 1, j-1
					summa_2 = summa_2 + U(t2,i)*U(t2,j)
					prop2 = summa_2
				end do
			else
				prop2 = 0
			end if
				U(i,j) = (A(i,j) - prop2)/U(i,i)
			end if
			summa_2 = 0
		end do
	end do
end function findU

function findY(U,b,n)  result(y)
	implicit none
	integer,intent(in) :: n
	real(8) :: U(n,n), y(n), b(n),summa, proma
	integer :: i, j
	summa = 0
	do i = 1,n
		if (i .EQ. 1) then
			y(i) = b(i)/U(i,i)
		else
			do j = 1, i-1
				summa = summa + U(i,j)*y(j)
				proma = summa
			end do
			y(i) = (b(i) - proma)/U(i,i)	
			summa = 0
		end if
	end do
end function findY

function findX(U,y,n)  result(x)
	implicit none
	integer,intent(in) :: n
	real(8) :: U(n,n), x(n), y(n),summa, proma
	integer :: i, j
	summa = 0
	x(n) = y(n)/U(n,n)
	do i = 1, n-1
		do j = 1 , i
			summa = summa + U(n-i,n-j+1)*x(n-j+1)
			proma = summa
		end do
		x(n-i) = (y(n-i) - proma)/U(n-i,n-i)
		summa = 0
	end do
end function findX

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


subroutine normNEV(U,Ut,B,X, n, cond)
implicit none
integer, intent(in) :: n
real(8) :: B(n), X(n), dB(n), dX(n),B2(n), U(n,n), Ut(n,n), X2(n), Y2(n), normINF,cond, norm_dx, norm_db
integer :: i,j
interface
	function findX(U,y,n)  result(x)
	implicit none
	integer,intent(in) :: n
	real(8) :: U(n,n), x(n), y(n)
	end function findX

	function findY(U,b,n)  result(y)
	implicit none
	integer,intent(in) :: n
	real(8) :: U(n,n), y(n), b(n)
	end function findY
end interface
	do I = 1, N 
	  call random_number(dB(I)) 
	  dB(I) = 2*(0.5 - dB(I))*B(I) 
	  B2(I) = B(I) + dB(I) 
end do
do i = 1,n
	!write(*,*) (U(i,j), j = 1, n)
end do
Y2 = findY(Ut,B2,n)
X2 = findX(U,Y2,n)
do i = 1, n
	dX(i) = dabs(X2(i) - X(i))
	! print*, Y2(i)
end do
norm_dx = normINF(dX,n)
norm_db = normINF(dB, n)
if (norm_dx/(normINF(X,n)) .LE. cond*(norm_db/(normINF(B,n)))) then
	write(*, '(a   e32.16  a   e32.16)') ' TRUE: ', norm_dx/(normINF(X,n)), '  <=  ', cond*(norm_db/(normINF(B,n)))
else
	write(*, '(a   e32.16  a   e32.16)') ' FALSE: ', norm_dx/(normINF(X,n)), '  <=  ', cond*(norm_db/(normINF(B,n)))
end if
end subroutine normNEV
	
program test
	implicit none
	integer :: n,i,j,k
	real(8), allocatable :: A(:,:), B(:), U(:,:), Ut(:,:), Y(:), X(:), cond
	interface
		function findX(U,y,n)  result(x)
		implicit none
		integer,intent(in) :: n
		real(8) :: U(n,n), x(n), y(n)
		end function findX
	
		function findY(U,b,n)  result(y)
		implicit none
		integer,intent(in) :: n
		real(8) :: U(n,n), y(n), b(n)
		end function findY
			
		function dot(H,G,n) result(M)
		implicit none
		integer, intent(in) :: n
		real(8) :: H(n,n), G(n,n), M(n,n)
		end function dot
	
		function transpos(U,n) result(Ut)
		implicit none
		integer, intent(in) :: n
		real(8) :: U(n,n),Ut(n,n)
		end function transpos
	
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
		
		function findU(A,n) result(U)
		implicit none
		integer, intent(in) :: n
		real(8) :: U(n,n), A(n,n), summa_sq,summa_2
		end function findU
	end interface
	write(*,*) ' Enter N '
	read(*,*) n
	allocate(A(n,n))
	allocate(U(n,n))
	allocate(Ut(n,n))
	allocate(Y(n))
	allocate(X(n))
	allocate(B(n))
	cond = 424812.1786362939
	A = readArray(n)
	B = readB(n)
	U = findU(A,n)
	Ut = transpos(U,n)
	!A = dot(U,Ut,n)
	Y = findY(Ut,B,n)
	X = findX(U,y,n)
	call normNEV(U,Ut,B,X, n,cond)
	open(20, file = '5_M_solve.txt',status='new')
	write(20,*) (x(j), j = 1,n)
	close(20)
	deallocate(A)
	deallocate(U)
	deallocate(Ut)
	deallocate(B)
	deallocate(Y)
	deallocate(X)
end program
	

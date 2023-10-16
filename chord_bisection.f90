real(8) function f(x, func_num)
	implicit none
	real(8) :: x
	integer :: func_num
	if (func_num .EQ. 0) then
		f = 2*x**4 + 8*x**3 + 8*x**2 -1
	else 
	 	f = 2*atan(x) - x + 3
	end if
	end function f


real(8) function chord_meth(a, b, func_num, eps, error_ch)
	implicit none
	real(8) ::  f, a, b, x, x0, eps, error_ch, max_1, min_1
	integer :: func_num, i,n
	x = a
	x0 = 0
	if (func_num .EQ. 0) then
		max_1 = 8*b**3+24*b**2+16*b
		min_1 = 8*a**3+24*a**2+16*a
	else
		max_1 = 2/(1+a**2)-1
		min_1 = 2/(1+b**2)-1
	end if
	do while (abs(x - x0)*dabs((max_1-min_1)/(max_1)) > eps)
		x0 = x
		x = x - f(x, func_num)*(x - b)/(f(x, func_num)-f(b, func_num))
	end do
	chord_meth = x
	if (func_num .EQ. 0) then
		error_ch = dabs(x-0.3065629648763764)
	else
		error_ch = dabs(x-5.80012972)
	end if
		write(*,'(a e32.16)') ' The error to the equation with CHORDS METHOD is: ' , error_ch
	end function chord_meth
	
real(8) function bisec_meth(a, b, eps, func_num, error_bis)
	implicit none
	real(8) :: f, x, a, b, c, eps, error_bis
	integer :: func_num
	do while ((b-a) > 2*eps)
		c = (b+a)*0.5
		if (f(a, func_num)*f(c, func_num) < 0) then
			b = c
		else 
			a = c
		end if
	end do
	bisec_meth = c
	if (func_num .EQ. 0) then
		error_bis = dabs(c-0.3065629648763764)
	else
		error_bis = dabs(c-5.80012972)
	end if
		write(*,'(a e32.16)') ' The error to the equation with BISECTION METHOD is: ' , error_bis
	end function bisec_meth

program solv_eq
	implicit none
	real(8) :: bisec_meth, a, b, eps, solve_bis, chord_meth, solve_ch, error_bis, error_ch,f
	integer :: func_num
		a = 5
		b = 7
		write(*,*) ' Enter the number of function(0 is for algebral function; 1 for transcendental function): num =  '
		read(*,*) func_num
		write(*,*) ' Enter the first search boundary: a = '
		read(*,*) a
		write(*,*) ' Enter the second search boundary: b = '
		read(*,*) b
		write(*,*) ' Enter the precision: eps =  '
		read(*,*) eps
		solve_ch = chord_meth(a,b, func_num, eps, error_ch)
		solve_bis = bisec_meth(a, b, eps, func_num, error_bis)
		write(*,'(a e12.7)') ' The solution to the equation with CHORD METHOD is: x = ', solve_ch
		write(*,'(a e12.7)') ' The solution to the equation with BISECTION METHOD is: x = ', solve_bis
	end program

module helper
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none
    contains

    function f(x) result(x3)
        real(dp), intent(in)  :: x
        real(dp) :: x3
        x3 = x*x*x
    end function

    subroutine trap(a,b,n,sum,x_values,f_values)
        real(dp), intent(inout) :: sum
        real(dp), intent(in)    :: a,b
        integer , intent(in)    :: n
        real(dp), intent(inout) :: x_values(:),f_values(:)
        real(dp) :: dx
        integer :: i
        dx = (b-a)/n
        sum = 0.5*(f(a)+f(b)) ! end point contribution
        do i=1,n-1
            x_values(i) = a+(i)*dx
            f_values(i) = f(x_values(i))
            sum = sum + f_values(i)
        enddo
        sum = sum * dx
    end subroutine

end module

program main
    use helper
    implicit none
    integer :: n
    real(dp) :: a = 0., b = 1.0, result, result_an
    real(dp), allocatable :: x_values(:), f_values(:)
    integer :: i, num_args, stat
    character(len=64) :: arg
    num_args = command_argument_count()
    call get_command_argument(1,arg)
    read(arg,*,iostat=stat) n
    allocate(x_values(n), f_values(n))

    write(*,*)"numerical integration of f(x)=x^3 from a=",a," to b=",b,' with ',n,' points'
    call trap(a,b,n,result,x_values,f_values)
    result_an = 0.25*(b**4 - a**4)

    do i=1,n-1
        write(*,'(a,i4,a,f8.4,a,i4,a,f8.4)')'x[',i,']=',x_values(i),', f(x[',i,'])=',f_values(i)
    enddo

    write(*,'(a,f8.4)')"         result=",result
    write(*,'(a,f8.4)')"analytic result=",result_an
    write(*,'(a,e8.2)')"      abs error=",abs(result - result_an)
end program

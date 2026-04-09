module helper
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none
    contains

    function f(x) result(x3)
        implicit none
        real(dp), intent(in)  :: x
        real(dp) :: x3
        x3 = x*x*x
    end function

    subroutine trap(a,b,n,sum)
        implicit none
        real(dp), intent(inout) :: sum
        real(dp), intent(in)    :: a,b
        integer , intent(in)    :: n
        real(dp) :: dx
        integer :: i
        dx = (b-a)/n
        sum = 0.5*(f(a)+f(b)) ! end point contribution
        do i=1,n-1
            sum = sum + f(a+i*dx)
        enddo
        sum = sum * dx
    end subroutine

end module

program main
    use helper
    implicit none
    real(dp) :: a = 0., b = 1.0, result, result_an
    integer  :: n = 1000

    write(*,*)"numerical integration of f(x)=x^3 from a=",a," to b=",b
    call trap(a,b,n,result)
    result_an = 0.25*(b**4 - a**4);
    write(*,'(a,f8.4)')"         result=",result
    write(*,'(a,f8.4)')"analytic result=",result_an
    write(*,'(a,e8.2)')"      abs error=",abs(result - result_an)
end program
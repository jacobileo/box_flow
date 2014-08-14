module boxflow_solvers
    implicit none
    integer, parameter :: DBL=selected_real_kind(20, 100), MAXIMUM_TIME_OF_RECURSION=5000

contains

    pure subroutine Interior_omega(w, p, dt, dx, dy, m, n, nu)
        implicit none
        integer, intent(in) :: m, n
        real(kind=DBL), dimension(m,n), intent(inout) :: w !omega
        real(kind=DBL), dimension(m,n), intent(in) :: p !psi
        real(kind=DBL), intent(in) :: dt, dx, dy, nu
        real(kind=DBL) :: Dxp, Dxw, Dyp, Dyw, Dhw
        integer :: i, j

        do i=2,m-1
            do j=2,n-1
                Dxp = 0.5D0*(p(i+1,j) - p(i-1,j))/dx
                Dxw = 0.5D0*(w(i+1,j) - w(i-1,j))/dx
                Dyp = 0.5D0*(p(i,j+1) - p(i,j-1))/dy
                Dyw = 0.5D0*(w(i,j+1) - w(i,j-1))/dy
                Dhw = (w(i+1,j)-2.D0*w(i,j)+w(i-1,j))/dx**2 + (w(i,j+1)-2.D0*w(i,j)+w(i,j-1))/dy**2
                w(i,j) = w(i,j) + dt*(Dxp*Dyw - Dyp*Dxw + nu*Dhw)
            end do
        end do
    end subroutine Interior_omega

    subroutine SOR (T, z, dx, dy, w, eol, m, n)
        !===============================
        ! T     The temperature field
        ! b     \frac{\Delta x}{\Delta y}
        ! w     factor of SOR
        ! eol   error of limit
        ! r     time of recursion
        ! m,n   dimension of matrix T
        !==============================
        implicit none
        integer, intent(in) :: n, m    
        real(kind=DBL), intent(in) :: w, eol, dx, dy
        real(kind=DBL), intent(inout), dimension(m, n) :: T
        real(kind=DBL), intent(in), dimension(m, n) :: z
        !integer, intent(out) :: r
        integer :: r
        real(kind=DBL), dimension(m, n) :: T_tmp
        real(kind=DBL) :: b2, tmp
        integer :: i, j, k

        b2 = dx**2 / dy**2
        !write(*,*) 'Solving for w =', w, ' with n = ', n
        T_tmp = T
        iteration: do k=1, MAXIMUM_TIME_OF_RECURSION
            do i=2,m-1
                do j=2,n-1
                    tmp = T(i+1,j) + T(i-1,j) + b2*(T(i,j+1)+T(i,j-1)) +z(i,j)*dx*dx
                    T(i,j) = (1-w) * T(i,j) + 0.5D0 * w / (1+b2) * tmp
                end do
            end do
            if(sum(abs(T-T_tmp)) < eol) then
                r = k
                write(*,*) '        recursion for ', r, 'times.'
                exit iteration
            end if
            T_tmp = T
        end do iteration
        if(k > MAXIMUM_TIME_OF_RECURSION - 2) then
            write(*,*) 'Failed to converge.'
        end if
    end subroutine SOR

    pure subroutine Vortex_Boundary_Condition(omega, psi, dx, dy, m, n)
        implicit none
        integer, intent(in) :: m, n
        real(kind=DBL), dimension(m,n), intent(inout) :: omega
        real(kind=DBL), dimension(m,n), intent(in) :: psi
        real(kind=DBL), intent(in) :: dx, dy
        real(kind=DBL) :: u, x
        integer :: i, j

        do i=2,m-1
            !u = Velocity(real(i,DBL)/real(n,DBL))
            x = real(i,DBL) / real(n,DBL)
            u = 1.6D1 * x**2 * (1-x)**2
            omega(i,n) = -1.D0*(2D0*psi(i,n-1) + u*dy)/dy**2
            omega(i,1) = -1.D0*(2D0*psi(i,2))/dy**2
        end do
        do j=2,n-1
            omega(1,j) = -1.D0*(2D0*psi(2,j))/dx**2
            omega(n,j) = -1.D0*(2D0*psi(n-1,j))/dx**2
        end do

        !contains

        !    pure function Velocity(x)
        !        implicit none
        !        real(kind=DBL), intent(in) :: x
        !        real(kind=DBL) :: Velocity
        !        Velocity = 1.6D1 * x**2 * (1-x)**2
        !    end function Velocity

        end subroutine Vortex_Boundary_Condition

    end module boxflow_solvers


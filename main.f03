program main
    use boxflow_solvers
    implicit none

    integer, parameter :: GRID_MAX=100, ITERATION_MAX=1000000
    integer :: i, j, p, q
    real(kind=DBL), dimension(GRID_MAX, GRID_MAX) :: omega, psi, omega_tmp, psi_tmp
    
    real(kind=DBL), parameter :: dt=1.D-3, dx=1.D0/real(GRID_MAX, DBL), &
        dy=1.D0/real(GRID_MAX, DBL), w=1.94D0, eol=1.0D-5, nu=1.D-5, eol2=1.0D-3

    character(len=30) :: filename, filename_head, filename_tail
    character(len=20) :: formatstr="(???(F10.6, ' '))"
    filename_head = 'output_psi_'
    write(formatstr(2:4), '(i3)') GRID_MAX

    forall(i=1:GRID_MAX, j=1:GRID_MAX)
        psi(i,j) = 0.D0
        omega(i,j) = 0.D0
    end forall

    iteration: do i=1,ITERATION_MAX
        write(*,*) 'The ', i, 'th step of time forward.'
        psi_tmp = psi
        omega_tmp = omega
        call Interior_omega (omega, psi, dt, dx, dy, GRID_MAX, GRID_MAX, nu)
        write(*,*) '        STEP 1'
        call SOR (psi, omega, dx, dy, w, eol, GRID_MAX, GRID_MAX)
        write(*,*) '        STEP 2'
        call Vortex_Boundary_Condition (omega, psi, dx, dy, GRID_MAX, GRID_MAX)
        write(*,*) '        STEP 3'

        if(mod(i,1000) == 0) then
            j = i/1000
            filename_tail = IntToStr(j)
            filename = trim(filename_head) // trim(filename_tail) // '.txt'
            open(unit=j, file=filename, status='new', action='write')
            write(j,fmt=formatstr) ((psi(p,q), q=1,GRID_MAX), p=1,GRID_MAX)
            close(unit=j)
        end if
        
        if(i > ITERATION_MAX - 2) then
            write(*,*) 'Fail to converge.'
        else if(sum(abs(omega-omega_tmp)) < eol2 .and. sum(abs(psi-psi_tmp)) < eol2) then
            exit iteration
        end if
    end do iteration

    contains

        pure function IntToStr(j)
            integer, intent(in) :: j
            character(len=5) :: IntToStr
            integer :: i, n

            IntToStr = ''
            n = j
            do i=1,5
                if(n==0) then
                    exit
                end if
                IntToStr = achar(mod(n,10)+iachar('0')) // IntToStr
                n = n / 10
            end do
        end function IntToStr

end program main


! =============================================================================
! Program: brownian_settling
! Phase 3 – Ermak–McCammon integration loop for a settling Brownian sphere.
!
! Algorithm (overdamped Langevin / Euler–Maruyama):
!
!   x(t+dt) = x(t) +            sqrt(2 D dt) * Rx
!   y(t+dt) = y(t) +            sqrt(2 D dt) * Ry
!   z(t+dt) = z(t) - U_s*dt  + sqrt(2 D dt) * Rz
!
! where Rx, Ry, Rz ~ N(0,1) independently.
! The deterministic term -U_s*dt is the sedimentation drift (z is upward).
! This is equivalent to (D/kT) * F_ext * dt with F_ext = -F_gravity_net.
!
! Output: trajectory.dat  (columns: step, x, y, z  [all in SI metres])
! =============================================================================
program brownian_settling
    use physical_params
    use rng_module
    implicit none

    ! ---- Simulation parameters (tune these) ---------------------------------
    integer,  parameter :: N_steps   = 500000    ! total integration steps
    integer,  parameter :: log_every = 500        ! write output every this many steps
    real(8),  parameter :: dt        = 1.0d-2     ! time step [s]
    integer,  parameter :: rng_seed  = 42         ! reproducibility

    ! ---- Local variables ----------------------------------------------------
    real(8) :: x, y, z            ! position [m]
    real(8) :: rx, ry, rz         ! Gaussian random displacements
    real(8) :: noise              ! sqrt(2 D dt) – precomputed
    real(8) :: t_now              ! current time [s]
    integer :: i, funit
    integer, allocatable :: seed(:)
    integer :: seed_size

    ! ---- Initialise ---------------------------------------------------------
    call init_params()

    ! Seed the RNG for reproducibility
    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    seed = rng_seed
    call random_seed(put=seed)
    deallocate(seed)

    ! Precompute noise prefactor
    noise = sqrt(2.0d0 * D * dt)

    ! Initial position at the origin
    x = 0.0d0;  y = 0.0d0;  z = 0.0d0

    ! ---- Open output file ---------------------------------------------------
    open(newunit=funit, file='trajectory.dat', status='replace', action='write')
    write(funit, '(A)') '# Ermak-McCammon BD simulation: settling sphere'
    write(funit, '(A,ES12.4,A,I0,A,I0)') &
        '# dt = ', dt, ' s  |  N_steps = ', N_steps, '  |  log_every = ', log_every
    write(funit, '(A)') '#      step          time[s]            x[m]            y[m]            z[m]'
    write(funit, '(I10, 4ES16.8)') 0, 0.0d0, x, y, z

    ! ---- Main integration loop ----------------------------------------------
    write(*,'(A)') '  Running settling simulation ...'

    do i = 1, N_steps
        call gaussian_3(rx, ry, rz)

        x = x             + noise * rx
        y = y             + noise * ry
        z = z - U_s * dt  + noise * rz

        if (mod(i, log_every) == 0) then
            t_now = dble(i) * dt
            write(funit, '(I10, 4ES16.8)') i, t_now, x, y, z
        end if
    end do

    close(funit)

    write(*,'(A,F8.1,A)') '  Done.  Total simulated time: ', dble(N_steps)*dt, ' s'
    write(*,'(A,ES10.3,A)') '  Expected z drift:            ', -U_s * dble(N_steps) * dt, ' m'
    write(*,'(A)') '  Output written to trajectory.dat'
    write(*,*)

end program brownian_settling

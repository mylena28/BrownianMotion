! =============================================================================
! Program: msd_validation
! Phase 4 – Mean Squared Displacement (MSD) validation.
!
! Runs a purely diffusive BD simulation (gravity OFF) and computes the
! time-averaged MSD:
!
!   MSD(tau) = < |r(t + tau) - r(t)|^2 >_t
!
! For a free Brownian particle theory predicts:
!   MSD_1D(tau) = 2 D tau       (per axis, e.g. <Dx^2>)
!   MSD_3D(tau) = 6 D tau       (full 3-D)
!
! Outputs:
!   no_gravity_traj.dat   – trajectory (for visual check)
!   msd_output.dat        – tau, MSD_x, MSD_y, MSD_z, MSD_3D, MSD_theory_3D
! =============================================================================
program msd_validation
    use physical_params
    use rng_module
    implicit none

    ! ---- Simulation parameters ----------------------------------------------
    integer,  parameter :: N_steps   = 100000   ! steps to integrate
    real(8),  parameter :: dt        = 1.0d-2    ! time step [s]
    integer,  parameter :: N_lags    = 2000      ! number of lag values to compute
    integer,  parameter :: rng_seed  = 137

    ! ---- Storage ------------------------------------------------------------
    real(8), allocatable :: px(:), py(:), pz(:)   ! stored trajectory
    real(8) :: noise
    real(8) :: rx, ry, rz
    real(8) :: msd_x, msd_y, msd_z, msd_3d
    real(8) :: dx, dy, dz
    integer :: i, k, n_avg
    integer :: funit_traj, funit_msd
    integer, allocatable :: seed(:)
    integer :: seed_size
    integer :: lag_step      ! lag in number of steps

    ! ---- Initialise ---------------------------------------------------------
    call init_params()

    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    seed = rng_seed
    call random_seed(put=seed)
    deallocate(seed)

    noise = sqrt(2.0d0 * D * dt)

    allocate(px(0:N_steps), py(0:N_steps), pz(0:N_steps))
    px(0) = 0.0d0;  py(0) = 0.0d0;  pz(0) = 0.0d0

    ! ---- Integrate (no gravity) --------------------------------------------
    write(*,'(A)') '  Running no-gravity simulation for MSD validation ...'

    do i = 1, N_steps
        call gaussian_3(rx, ry, rz)
        px(i) = px(i-1) + noise * rx
        py(i) = py(i-1) + noise * ry
        pz(i) = pz(i-1) + noise * rz
    end do

    ! ---- Write no-gravity trajectory (log every 1000 steps) ----------------
    open(newunit=funit_traj, file='no_gravity_traj.dat', status='replace', action='write')
    write(funit_traj,'(A)') '# No-gravity BD trajectory'
    write(funit_traj,'(A)') '#      step          time[s]            x[m]            y[m]            z[m]'
    do i = 0, N_steps, 1000
        write(funit_traj,'(I10,4ES16.8)') i, dble(i)*dt, px(i), py(i), pz(i)
    end do
    close(funit_traj)

    ! ---- Compute time-averaged MSD -----------------------------------------
    ! Lag values: 1, 2, ..., N_lags  (with N_lags << N_steps for good averaging)
    ! The maximum lag is N_lags steps = N_lags * dt seconds.
    ! Using N_lags <= N_steps/10 ensures each lag has >= 90% of the trajectory
    ! contributing to the average.

    write(*,'(A)') '  Computing time-averaged MSD ...'

    open(newunit=funit_msd, file='msd_output.dat', status='replace', action='write')
    write(funit_msd,'(A)') '# MSD validation: free Brownian particle'
    write(funit_msd,'(A,ES12.4,A)') '# D = ', D, ' m^2/s'
    write(funit_msd,'(A)') &
        '#    tau[s]         MSD_x[m^2]      MSD_y[m^2]      MSD_z[m^2]' // &
        '      MSD_3D[m^2]   Theory_3D[m^2]'

    do k = 1, N_lags
        lag_step = k * (N_steps / (10 * N_lags))   ! subsample lags evenly
        if (lag_step < 1) lag_step = k              ! fall back for small N_steps
        if (lag_step > N_steps / 4) exit            ! hard cap: max lag = N/4

        n_avg = N_steps - lag_step
        msd_x = 0.0d0;  msd_y = 0.0d0;  msd_z = 0.0d0

        do i = 0, n_avg - 1
            dx = px(i + lag_step) - px(i)
            dy = py(i + lag_step) - py(i)
            dz = pz(i + lag_step) - pz(i)
            msd_x = msd_x + dx*dx
            msd_y = msd_y + dy*dy
            msd_z = msd_z + dz*dz
        end do

        msd_x = msd_x / dble(n_avg)
        msd_y = msd_y / dble(n_avg)
        msd_z = msd_z / dble(n_avg)
        msd_3d = msd_x + msd_y + msd_z

        write(funit_msd,'(6ES16.8)') &
            dble(lag_step) * dt, &   ! tau [s]
            msd_x, msd_y, msd_z, &
            msd_3d, &
            6.0d0 * D * dble(lag_step) * dt    ! theory: 6 D tau
    end do

    close(funit_msd)
    deallocate(px, py, pz)

    write(*,'(A)') '  MSD written to msd_output.dat'
    write(*,'(A)') '  No-gravity trajectory written to no_gravity_traj.dat'
    write(*,*)

end program msd_validation

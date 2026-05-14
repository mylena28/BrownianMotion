! =============================================================================
! Module: rng_module
! Phase 2 – Gaussian random number generation via the Box–Muller transform.
!
! Fortran's intrinsic random_number() returns uniform U[0,1] deviates.
! Brownian motion requires N(0,1) deviates; Box–Muller achieves this exactly:
!
!   Given u1, u2 ~ Uniform(0,1):
!     z1 = sqrt(-2 ln u1) * cos(2 pi u2)    <- N(0,1)
!     z2 = sqrt(-2 ln u1) * sin(2 pi u2)    <- N(0,1), independent of z1
! =============================================================================
module rng_module
    implicit none
    real(8), parameter :: PI_RNG = 4.0d0 * atan(1.0d0)

contains

    ! -------------------------------------------------------------------------
    ! box_muller: generate one pair of independent N(0,1) deviates.
    ! u1 is drawn in a loop to guard against the log(0) singularity.
    ! -------------------------------------------------------------------------
    subroutine box_muller(z1, z2)
        real(8), intent(out) :: z1, z2
        real(8) :: u1, u2

        do
            call random_number(u1)
            if (u1 > 0.0d0) exit   ! avoid log(0)
        end do
        call random_number(u2)

        z1 = sqrt(-2.0d0 * log(u1)) * cos(2.0d0 * PI_RNG * u2)
        z2 = sqrt(-2.0d0 * log(u1)) * sin(2.0d0 * PI_RNG * u2)
    end subroutine box_muller

    ! -------------------------------------------------------------------------
    ! gaussian_3: produce three independent N(0,1) numbers (Rx, Ry, Rz).
    ! Two Box–Muller calls yield 4 values; we use 3 and discard 1.
    ! -------------------------------------------------------------------------
    subroutine gaussian_3(rx, ry, rz)
        real(8), intent(out) :: rx, ry, rz
        real(8) :: z1, z2

        call box_muller(z1, z2)
        rx = z1
        ry = z2
        call box_muller(z1, z2)
        rz = z1          ! z2 discarded; all three are independent
    end subroutine gaussian_3

end module rng_module

! =============================================================================
! Module: physical_params
! Phase 1 – Physical constants and derived hydrodynamic quantities for a
! settling Brownian sphere in a Newtonian fluid.
! =============================================================================
module physical_params
    implicit none

    ! -------------------------------------------------------------------------
    ! Fundamental constants
    ! -------------------------------------------------------------------------
    real(8), parameter :: PI    = 4.0d0 * atan(1.0d0)
    real(8), parameter :: kB    = 1.380649d-23   ! Boltzmann constant   [J/K]
    real(8), parameter :: g_acc = 9.81d0          ! Gravitational accel. [m/s^2]

    ! -------------------------------------------------------------------------
    ! Fluid / particle properties (edit these to change the physical problem)
    ! -------------------------------------------------------------------------
    real(8), parameter :: T     = 298.15d0   ! Temperature               [K]
    real(8), parameter :: mu    = 1.0d-3     ! Dynamic viscosity (water) [Pa·s]
    real(8), parameter :: a     = 1.0d-6     ! Sphere radius             [m]
    real(8), parameter :: rho_p = 2000.0d0   ! Particle density (SiO2)   [kg/m^3]
    real(8), parameter :: rho_f = 1000.0d0   ! Fluid density (water)     [kg/m^3]

    ! -------------------------------------------------------------------------
    ! Derived quantities – computed at runtime by init_params()
    ! -------------------------------------------------------------------------
    real(8) :: D       ! Stokes–Einstein diffusion coefficient  [m^2/s]
    real(8) :: U_s     ! Sedimentation (settling) velocity      [m/s]
    real(8) :: gamma_d ! Stokes drag coefficient                [N·s/m]

contains

    ! -------------------------------------------------------------------------
    ! Compute derived quantities and print a summary to stdout.
    ! Must be called once before the integration loop.
    ! -------------------------------------------------------------------------
    subroutine init_params()
        real(8) :: F_ext, V_sphere

        ! Stokes drag:  gamma = 6 pi mu a
        gamma_d  = 6.0d0 * PI * mu * a

        ! Stokes–Einstein:  D = k_B T / gamma
        D        = kB * T / gamma_d

        ! Net gravitational force (buoyancy-corrected):  F = V (rho_p - rho_f) g
        V_sphere = (4.0d0 / 3.0d0) * PI * a**3
        F_ext    = V_sphere * (rho_p - rho_f) * g_acc

        ! Sedimentation velocity:  U = F / gamma
        U_s      = F_ext / gamma_d

        write(*,'(/,A)') repeat('=', 55)
        write(*,'(A)')   '  Ermak–McCammon Brownian Dynamics: Settling Sphere'
        write(*,'(A)')   repeat('=', 55)
        write(*,'(A,ES12.4,A)') '  Diffusion coefficient  D   = ', D,       ' m^2/s'
        write(*,'(A,ES12.4,A)') '  Sedimentation velocity U   = ', U_s,     ' m/s'
        write(*,'(A,ES12.4,A)') '  Drag coefficient       gam = ', gamma_d, ' N·s/m'
        write(*,'(A,F10.4)')    '  Peclet number       Pe=Ua/D = ', U_s * a / D
        write(*,'(A)')   repeat('=', 55)
        write(*,*)
    end subroutine init_params

end module physical_params

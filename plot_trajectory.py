"""
plot_trajectory.py – Phase 4 visualisation
Plots:
  1. 3-D settling trajectory  (trajectory.dat)
  2. MSD vs lag time, linear scale, with theory line  (msd_output.dat)
  3. MSD log-log plot to confirm the diffusive slope = 1
Run after:
    make all_run
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D   # noqa: F401  (needed for 3-D projection)
import sys
import os

# ---------------------------------------------------------------------------
# Helper: load a whitespace-/comment-delimited data file safely
# ---------------------------------------------------------------------------
def load(fname):
    if not os.path.isfile(fname):
        print(f"[ERROR] File not found: {fname}")
        print("  Run  'make all_run'  first.")
        sys.exit(1)
    return np.loadtxt(fname, comments='#')


# ---------------------------------------------------------------------------
# 1. 3-D settling trajectory
# ---------------------------------------------------------------------------
traj = load('trajectory.dat')
step_col, t_col, x_col, y_col, z_col = 0, 1, 2, 3, 4

x  = traj[:, x_col]
y  = traj[:, y_col]
z  = traj[:, z_col]
t  = traj[:, t_col]

# Convert to micrometres for readability
x_um = x * 1e6
y_um = y * 1e6
z_um = z * 1e6

fig = plt.figure(figsize=(14, 6))

# --- 3-D view ---------------------------------------------------------------
ax3d = fig.add_subplot(121, projection='3d')
sc = ax3d.scatter(x_um, y_um, z_um, c=t, cmap='viridis', s=3, alpha=0.8)
ax3d.set_xlabel('x  [µm]')
ax3d.set_ylabel('y  [µm]')
ax3d.set_zlabel('z  [µm]')
ax3d.set_title('3-D Settling Trajectory\n(colour = time)')
cbar = fig.colorbar(sc, ax=ax3d, shrink=0.6, pad=0.1)
cbar.set_label('time  [s]')

# --- z(t) panel – shows steady drift with superimposed noise ----------------
ax_zt = fig.add_subplot(122)
ax_zt.plot(t, z_um, lw=0.5, color='royalblue', label='z(t)  simulation')

# Theoretical mean position
U_s_label = (z_um[-1] - z_um[0]) / (t[-1] - t[0])   # slope from data
t_line = np.array([t[0], t[-1]])
ax_zt.plot(t_line, z_um[0] + U_s_label * t_line,
           'r--', lw=1.5, label=r'$-U_s\,t$  (mean drift)')

ax_zt.set_xlabel('time  [s]')
ax_zt.set_ylabel('z  [µm]')
ax_zt.set_title('Vertical position  z(t)\nSteady settling + Brownian fluctuations')
ax_zt.legend(fontsize=9)
ax_zt.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('trajectory_3d.png', dpi=150)
print('Saved: trajectory_3d.png')
plt.show()


# ---------------------------------------------------------------------------
# 2 & 3. MSD validation
# ---------------------------------------------------------------------------
if not os.path.isfile('msd_output.dat'):
    print("[INFO] msd_output.dat not found – skipping MSD plots.")
    print("  Run  './msd_validate'  to generate it.")
    sys.exit(0)

msd = load('msd_output.dat')
# columns: tau, MSD_x, MSD_y, MSD_z, MSD_3D, Theory_3D
tau        = msd[:, 0]
msd_x      = msd[:, 1]
msd_y      = msd[:, 2]
msd_z      = msd[:, 3]
msd_3d     = msd[:, 4]
theory_3d  = msd[:, 5]

fig2, axes = plt.subplots(1, 2, figsize=(13, 5))

# ---- Linear MSD plot -------------------------------------------------------
ax = axes[0]
ax.plot(tau, msd_3d * 1e12,  lw=1.5, color='steelblue', label='MSD 3D (simulation)')
ax.plot(tau, msd_x  * 1e12,  lw=1.0, color='tomato',    ls='--', label=r'$\langle\Delta x^2\rangle$')
ax.plot(tau, msd_y  * 1e12,  lw=1.0, color='seagreen',  ls='--', label=r'$\langle\Delta y^2\rangle$')
ax.plot(tau, msd_z  * 1e12,  lw=1.0, color='orchid',    ls='--', label=r'$\langle\Delta z^2\rangle$')
ax.plot(tau, theory_3d * 1e12, lw=2.0, color='black', ls=':', label=r'Theory: $6D\tau$')
ax.set_xlabel(r'$\tau$  [s]')
ax.set_ylabel(r'MSD  [µm²]')
ax.set_title('Mean Squared Displacement\n(free diffusion – no gravity)')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# ---- Log-log MSD plot – slope = 1 confirms normal diffusion ---------------
ax2 = axes[1]
mask = (tau > 0) & (msd_3d > 0) & (theory_3d > 0)
ax2.loglog(tau[mask], msd_3d[mask]    * 1e12, lw=1.5, color='steelblue', label='MSD 3D (simulation)')
ax2.loglog(tau[mask], theory_3d[mask] * 1e12, lw=2.0, color='black', ls=':', label=r'Theory: $6D\tau$  (slope = 1)')

# Reference slope-1 line
tau_ref = np.array([tau[mask][0], tau[mask][-1]])
y_ref   = theory_3d[mask][0] * 1e12 * (tau_ref / tau[mask][0])
ax2.loglog(tau_ref, y_ref, 'r--', lw=1.2, alpha=0.6, label='Slope = 1 guide')

ax2.set_xlabel(r'$\tau$  [s]')
ax2.set_ylabel(r'MSD  [µm²]')
ax2.set_title('MSD log–log plot\nSlope = 1 ⟹ normal (Fickian) diffusion')
ax2.legend(fontsize=9)
ax2.grid(True, which='both', alpha=0.3)

plt.tight_layout()
plt.savefig('msd_validation.png', dpi=150)
print('Saved: msd_validation.png')
plt.show()

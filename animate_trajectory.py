"""
animate_trajectory.py
Animates the settling Brownian sphere trajectory from trajectory.dat.

Layout:
  Left  – 3-D view: trail grows over time, sphere marker moves, view rotates.
  Right – z(t) panel: red dot and vertical line track the current time.

Outputs:
  trajectory_settling.gif   (Pillow writer, ~200 frames)
  trajectory_settling.mp4   (FFmpeg writer, all frames, 25 fps)

Run after:  make run
"""

import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D   # noqa: F401

# ---------------------------------------------------------------------------
# Load trajectory
# ---------------------------------------------------------------------------
fname = 'trajectory.dat'
if not os.path.isfile(fname):
    print(f"[ERROR] {fname} not found – run  'make run'  first.")
    sys.exit(1)

data = np.loadtxt(fname, comments='#')
x_m = data[:, 2];  y_m = data[:, 3];  z_m = data[:, 4]
t   = data[:, 1]

# Convert to micrometres for readability
x = x_m * 1e6;  y = y_m * 1e6;  z = z_m * 1e6
N = len(t)

print(f"Loaded {N} trajectory points  (t = {t[0]:.0f} … {t[-1]:.0f} s)")

# ---------------------------------------------------------------------------
# Subsample: use N_frames evenly-spaced indices for the GIF
# The MP4 uses all N points.
# ---------------------------------------------------------------------------
N_frames_gif = 200
gif_indices  = np.round(np.linspace(0, N - 1, N_frames_gif)).astype(int)

# ---------------------------------------------------------------------------
# Axis limits (fixed throughout animation)
# ---------------------------------------------------------------------------
pad_xy = max(abs(x).max(), abs(y).max()) * 0.15 + 2
x_lim  = [x.min() - pad_xy, x.max() + pad_xy]
y_lim  = [y.min() - pad_xy, y.max() + pad_xy]
z_lim  = [z.min() - abs(z.max() - z.min()) * 0.05,
          z.max() + abs(z.max() - z.min()) * 0.05]

# Colormap for trail (viridis, mapped to time)
norm   = plt.Normalize(t.min(), t.max())
cmap   = cm.viridis


# ===========================================================================
# Figure factory: creates a fresh figure + axes + static artists
# (called once per output format so each writer gets its own figure)
# ===========================================================================
def make_figure():
    fig = plt.figure(figsize=(13, 6), dpi=100)
    fig.patch.set_facecolor('#0e1117')

    ax3 = fig.add_subplot(121, projection='3d')
    ax2 = fig.add_subplot(122)

    for ax in (ax2,):
        ax.set_facecolor('#161b22')
        for sp in ax.spines.values():
            sp.set_color('#30363d')
        ax.tick_params(colors='#8b949e')
        ax.xaxis.label.set_color('#8b949e')
        ax.yaxis.label.set_color('#8b949e')
        ax.title.set_color('#c9d1d9')

    ax3.set_facecolor('#161b22')
    ax3.tick_params(colors='#8b949e')
    ax3.xaxis.pane.fill = False
    ax3.yaxis.pane.fill = False
    ax3.zaxis.pane.fill = False
    ax3.xaxis.pane.set_edgecolor('#30363d')
    ax3.yaxis.pane.set_edgecolor('#30363d')
    ax3.zaxis.pane.set_edgecolor('#30363d')
    ax3.tick_params(colors='#8b949e')

    # --- 3-D static: ghost of the complete future path ---------------------
    ax3.plot(x, y, z, lw=0.4, color='#30363d', alpha=0.5, zorder=1)

    ax3.set_xlim(x_lim);  ax3.set_ylim(y_lim);  ax3.set_zlim(z_lim)
    ax3.set_xlabel('x  [µm]', color='#8b949e', labelpad=6)
    ax3.set_ylabel('y  [µm]', color='#8b949e', labelpad=6)
    ax3.set_zlabel('z  [µm]', color='#8b949e', labelpad=6)
    ax3.view_init(elev=20, azim=30)

    # Dynamic artists on the 3-D axes
    trail_line, = ax3.plot([], [], [], lw=1.2, color='#58a6ff',
                           alpha=0.85, zorder=2)
    sphere_dot  = ax3.scatter([], [], [], s=80, c='#f78166',
                              zorder=5, depthshade=False)
    title3d = ax3.set_title('', color='#c9d1d9', fontsize=10, pad=8)

    # --- z(t) panel --------------------------------------------------------
    ax2.plot(t, z, lw=0.6, color='#30363d', label='full path')
    ax2.set_xlim(t[0], t[-1])
    ax2.set_ylim(z_lim)
    ax2.set_xlabel('time  [s]')
    ax2.set_ylabel('z  [µm]')
    ax2.set_title('Vertical position  z(t)')
    ax2.grid(True, alpha=0.15, color='#30363d')

    vline, = ax2.plot([t[0], t[0]], z_lim, color='#f78166',
                      lw=1.2, ls='--', alpha=0.7)
    pt2d,  = ax2.plot([], [], 'o', color='#f78166', ms=6, zorder=5)

    # Colourbar (attached to z(t) panel to show time progression)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax2, shrink=0.7, pad=0.03)
    cbar.set_label('time  [s]', color='#8b949e')
    cbar.ax.yaxis.set_tick_params(color='#8b949e')
    plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='#8b949e')

    fig.tight_layout(pad=2.0)

    return fig, ax3, ax2, trail_line, sphere_dot, title3d, vline, pt2d


# ===========================================================================
# Update function (shared between GIF and MP4 renders)
# ===========================================================================
def make_update(indices, trail_line, sphere_dot, title3d, vline, pt2d,
                ax3, rotate=True):
    """Return a closure over the chosen index array."""
    def update(frame):
        k = indices[frame]

        # Growing trail coloured by time
        trail_line.set_data(x[:k+1], y[:k+1])
        trail_line.set_3d_properties(z[:k+1])

        # Sphere position
        sphere_dot._offsets3d = (np.array([x[k]]),
                                 np.array([y[k]]),
                                 np.array([z[k]]))

        # Rotate 3-D view: one full rotation over all frames
        if rotate:
            azim = 30 + 360 * frame / len(indices)
            ax3.view_init(elev=20, azim=azim)

        # Title
        title3d.set_text(
            f'3-D Settling Trajectory\n'
            f't = {t[k]:.0f} s   z = {z[k]:.1f} µm')

        # z(t) panel tracking
        vline.set_xdata([t[k], t[k]])
        pt2d.set_data([t[k]], [z[k]])

        return trail_line, sphere_dot, vline, pt2d

    return update


# ===========================================================================
# Render GIF
# ===========================================================================
print("\nRendering GIF  (200 frames) ...")
fig_g, ax3_g, ax2_g, tl_g, sd_g, tt_g, vl_g, p2_g = make_figure()

update_gif = make_update(gif_indices, tl_g, sd_g, tt_g, vl_g, p2_g,
                         ax3_g, rotate=True)
ani_gif = animation.FuncAnimation(fig_g, update_gif,
                                  frames=N_frames_gif,
                                  interval=80,    # ms between frames
                                  blit=False)

gif_writer = animation.PillowWriter(fps=12)
ani_gif.save('trajectory_settling.gif', writer=gif_writer, dpi=100)
plt.close(fig_g)
print("  Saved: trajectory_settling.gif")


# ===========================================================================
# Render MP4 (all trajectory points, 25 fps)
# ===========================================================================
print(f"\nRendering MP4  ({N} frames at 25 fps) ...")
fig_m, ax3_m, ax2_m, tl_m, sd_m, tt_m, vl_m, p2_m = make_figure()

all_indices  = np.arange(N)
update_mp4   = make_update(all_indices, tl_m, sd_m, tt_m, vl_m, p2_m,
                           ax3_m, rotate=True)
ani_mp4 = animation.FuncAnimation(fig_m, update_mp4,
                                  frames=N,
                                  interval=40,
                                  blit=False)

mp4_writer = animation.FFMpegWriter(fps=25, bitrate=1800,
                                    extra_args=['-vcodec', 'libx264',
                                                '-pix_fmt', 'yuv420p'])
ani_mp4.save('trajectory_settling.mp4', writer=mp4_writer, dpi=100)
plt.close(fig_m)
print("  Saved: trajectory_settling.mp4")

print("\nDone.")

"""Figures for the PB3D tutorial (page_tutorial), from the runs of the tutorial.

Inputs, relative to the folder it runs in:
  wout_cbm18a_cut_SI.nc     the tutorial's VMEC file (Examples/wout.nc)
  r/mem4000/                PB3D with Examples/PB3D.input, then POST with Examples/POST.input
  r/f4_fixed_keep/          PB3D with BC_style = 1 1, then POST with only plot_E_rec
  r/post2d/                 POST on the output of r/mem4000 for figure 3 (settings in the tutorial)
  r/post3d/                 POST on the output of r/mem4000 for figure 4 (settings in the tutorial)
Output: fig/tutorial_*.png, copied to Documentation/images/
Usage: python make_tutorial_figures.py   (python with numpy, matplotlib, h5py, netCDF4)
"""
import os
import numpy as np
import h5py
import netCDF4
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

os.makedirs("fig", exist_ok=True)

# colours: reference palette of the dataviz method (light surface)
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#9a9893"
GRID = "#e4e2dc"
BLUE, ORANGE, RED = "#2a78d6", "#eb6834", "#e34948"
BLUE_RAMP = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]
DIVERGING = LinearSegmentedColormap.from_list(
    "bgr", ["#0d366b", "#256abf", "#86b6ef", "#f0efec", "#f2a3a2", "#e34948", "#8f2423"])

plt.rcParams.update({
    "font.size": 10, "axes.titlesize": 10.5, "axes.labelsize": 10,
    "text.color": INK, "axes.labelcolor": INK, "axes.edgecolor": MUTED,
    "xtick.color": INK2, "ytick.color": INK2, "axes.titlecolor": INK,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "savefig.dpi": 160, "savefig.bbox": "tight", "figure.facecolor": "white",
})


def wout(fn):
    d = netCDF4.Dataset(fn)
    return {k: np.array(d.variables[k][:]) for k in d.variables}


def surface_RZ(w, js, th):
    xm = w["xm"]
    R = (w["rmnc"][js][:, None] * np.cos(np.outer(xm, th))).sum(0)
    Z = (w["zmns"][js][:, None] * np.sin(np.outer(xm, th))).sum(0)
    return R, Z


cut = wout("wout_cbm18a_cut_SI.nc")
th = np.linspace(0, 2 * np.pi, 721)
psi_n = cut["chi"] / cut["chi"][-1]                 # normalized poloidal flux

# ---------------------------------------------------------------- figure 1: equilibrium
fig = plt.figure(figsize=(9.2, 4.3))
gs = fig.add_gridspec(2, 2, width_ratios=[0.9, 1], hspace=0.12, wspace=0.32)
ax = fig.add_subplot(gs[:, 0])
ax.grid(False)
Rb, Zb = surface_RZ(cut, -1, th)
for p in np.arange(0.1, 1.0, 0.1):
    js = np.argmin(abs(psi_n - p))
    R, Z = surface_RZ(cut, js, th)
    ax.plot(R, Z, color=BLUE_RAMP[3], lw=0.8)
ax.plot(Rb, Zb, color=BLUE_RAMP[8], lw=2)
ax.plot(cut["raxis_cc"][0], 0, marker="+", color=INK, ms=9, mew=1.5)
ax.set_aspect("equal")
ax.set_xlim(2.0, 5.0); ax.set_ylim(-1.55, 1.55)
ax.set_xlabel("$R$ (m)"); ax.set_ylabel("$Z$ (m)")
ax.set_title("flux surfaces, every 0.1 in $\\psi$", loc="left")
axp = fig.add_subplot(gs[0, 1])
axq = fig.add_subplot(gs[1, 1], sharex=axp)
axp.plot(psi_n, cut["presf"] / 1e3, color=BLUE, lw=2)
axp.set_ylabel("$p$ (kPa)")
axp.set_ylim(0, 31)
axp.tick_params(labelbottom=False)
axp.set_title("pressure and safety factor", loc="left")
axq.plot(psi_n, abs(1 / cut["iotaf"]), color=BLUE, lw=2)
axq.set_ylim(0, 3.3)
axq.set_ylabel("$|q|$")
axq.set_xlabel("normalized poloidal flux $\\psi$")
axq.set_xlim(0, 1)
fig.savefig("fig/tutorial_equilibrium.png")
plt.close(fig)

# ---------------------------------------------------------------- figure 2: Fourier harmonics
f = h5py.File("r/mem4000/PB3D_out.h5", "r")
vec = f["sol_R_1/RE_sol_vec/var"][:] + 1j * f["sol_R_1/IM_sol_vec/var"][:]
n_mod, n_r = 20, 400
vec = vec.reshape((n_mod, n_r), order="F")
mtab = np.loadtxt("r/mem4000/Data/modes_m_sol.dat")
r = mtab[:, 0]
mm = mtab[:, n_mod:].astype(int)
amp = {}
for i in range(n_r):
    for k in range(n_mod):
        amp.setdefault(mm[i, k], np.full(n_r, np.nan))[i] = abs(vec[k, i])
mmax = max(np.nanmax(a) for a in amp.values())
ms = sorted(amp)
cmap = LinearSegmentedColormap.from_list("blues", BLUE_RAMP[2:])
fig, ax = plt.subplots(figsize=(8.4, 3.9))
nq = 10 * abs(np.interp(r, psi_n, 1 / cut["iotaf"]))
for m in ms:
    a = amp[m] / mmax
    c = cmap((m - ms[0]) / (ms[-1] - ms[0]))
    ax.plot(r, a, color=c, lw=1.3)
    if m % 5 == 0:
        i = np.nanargmax(a)
        if a[i] > 0.12:
            ax.text(r[i], a[i] + 0.025, f"$m={m}$", ha="center", va="bottom", fontsize=8.5, color=INK2)
ax.set_xlim(0.3, 1.0); ax.set_ylim(0, 1.12)
ax.set_xlabel("normalized poloidal flux $\\psi$")
ax.set_ylabel("$|X_m|$ / largest value")
ax.set_title("Poloidal Fourier harmonics $X_m$ of the normal displacement, $n = 10$", loc="left")
sec = ax.secondary_xaxis("top", functions=(lambda x: 10 * abs(np.interp(x, psi_n, 1 / cut["iotaf"])),
                                            lambda y: np.interp(y, 10 * abs(1 / cut["iotaf"]), psi_n)))
sec.set_xlabel("$n|q|$", color=INK2)
sec.set_ticks([15, 20, 25, 30])
fig.savefig("fig/tutorial_harmonics.png")
plt.close(fig)

# ---------------------------------------------------------------- figure 3: mode and vacuum potential
xs = h5py.File("r/post2d/Data/1_sol_X_RE.h5", "r")
RX, ZX, VX = (xs[k][:, 0, :] for k in ("R_1", "Z_1", "var_1"))
vp = h5py.File("r/post2d/Data/1_vac_pot_RE.h5", "r")
RV, ZV, VV = (vp[k][0] for k in ("R_1", "Z_1", "var_1"))
inside = h5py.File("r/post2d/Data/vac_in_plasma.h5", "r")["var_1"][0] > 0.5
VV = np.ma.masked_where(inside, VV)
fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.2, 4.6), gridspec_kw=dict(wspace=0.08))
lim = abs(VX).max()
p1 = a1.pcolormesh(RX, ZX, VX, cmap=DIVERGING, norm=TwoSlopeNorm(0, -lim, lim), shading="gouraud", rasterized=True)
aV = abs(VV.compressed()); lim2 = aV[aV < 2 * np.percentile(aV, 99.9)].max()   # two points next to the boundary lie 3 to 5 times above the rest
p2 = a2.pcolormesh(RV, ZV, VV, cmap=DIVERGING, norm=TwoSlopeNorm(0, -lim2, lim2), shading="gouraud", rasterized=True)
Ri, Zi = RX[0], ZX[0]
for a in (a1, a2):
    a.plot(Rb, Zb, color=INK, lw=1.2)
    a.set_aspect("equal"); a.grid(False)
    a.set_xlim(1.9, 5.9); a.set_ylim(-1.9, 1.9)
    a.set_xlabel("$R$ (m)")
    for s in a.spines.values():
        s.set_visible(False)
a1.plot(Ri, Zi, color=INK2, lw=0.8, ls=(0, (3, 2)))
a1.text(3.55, 0, "$\\psi < 0.3$:\nnot solved", ha="center", va="center", fontsize=8.5, color=INK2)
a1.set_ylabel("$Z$ (m)")
a2.tick_params(labelleft=False)
a1.set_title("$X = \\vec{\\xi}\\cdot\\nabla\\psi$ in the plasma", loc="left")
a2.set_title("vacuum potential outside the plasma", loc="left")
for a, p in ((a1, p1), (a2, p2)):
    cb = fig.colorbar(p, ax=a, orientation="horizontal", fraction=0.05, pad=0.13, ticks=[])
    cb.outline.set_visible(False)
    cb.ax.text(0, -0.9, "negative", transform=cb.ax.transAxes, ha="left", va="top", fontsize=8.5, color=INK2)
    cb.ax.text(1, -0.9, "positive", transform=cb.ax.transAxes, ha="right", va="top", fontsize=8.5, color=INK2)
    cb.ax.text(0.5, -0.9, "real part, arbitrary units", transform=cb.ax.transAxes, ha="center", va="top", fontsize=8.5, color=INK2)
fig.savefig("fig/tutorial_mode.png")
plt.close(fig)

# ---------------------------------------------------------------- figure 4: energy terms
def energy(fn):
    lines = [l for l in open(fn) if not l.startswith("#")]
    row = lambda i: np.array(lines[i].split(), float)
    ev, ratio, E_pot, E_kin = row(0)
    return ev, ratio, row(4) / E_kin


ev, ratio, vals = energy("r/mem4000/POST_out_EN_R1.txt")
ev_f, ratio_f, vals_f = energy("r/f4_fixed_keep/POST_out_EN_R1.txt")
names = ["line bending, normal", "line bending, geodesic", "ballooning, normal", "ballooning, geodesic",
         "kink, normal", "kink, geodesic", "vacuum", "sum"]
free = np.append(vals, vals.sum())
fixed = np.append(vals_f, vals_f.sum())
fig, ax = plt.subplots(figsize=(7.8, 4.7))
y = np.arange(len(names))[::-1].astype(float)
h = 0.36
ax.barh(y + h / 2 + 0.02, free, height=h, color=BLUE, label="free boundary (this tutorial)")
ax.barh(y - h / 2 - 0.02, fixed, height=h, color=ORANGE, label="fixed boundary")
for yi, v in zip(y, free):
    ax.text(v + (0.004 if v >= 0 else -0.004), yi + h / 2 + 0.02, f"{v:+.3f}", va="center",
            ha="left" if v >= 0 else "right", fontsize=8, color=INK)
ax.text(fixed[-1] - 0.004, y[-1] - h / 2 - 0.02, f"{fixed[-1]:+.3f}", va="center", ha="right", fontsize=8, color=INK)
ax.set_yticks(y, names)
ax.get_yticklabels()[-1].set_fontweight("bold")
ax.axvline(0, color=INK2, lw=0.8)
ax.axhline(0.5, color=GRID, lw=0.8)
ax.grid(axis="y", visible=False)
ax.set_xlim(-0.235, 0.16)
ax.set_xlabel("contribution to $\\lambda = E_\\mathrm{pot}/E_\\mathrm{kin}$")
ax.set_title("Terms of the potential energy of the eigenmode, divided by its kinetic energy", loc="left")
ax.legend(frameon=False, loc="lower right", fontsize=8.5)
ax.annotate("", xy=(-0.1, -0.175), xytext=(-0.003, -0.175), xycoords=("data", "axes fraction"),
            textcoords=("data", "axes fraction"), arrowprops=dict(arrowstyle="->", color=INK2, lw=0.8))
ax.annotate("", xy=(0.1, -0.175), xytext=(0.003, -0.175), xycoords=("data", "axes fraction"),
            textcoords=("data", "axes fraction"), arrowprops=dict(arrowstyle="->", color=INK2, lw=0.8))
ax.text(-0.105, -0.175, "destabilizing", transform=ax.get_xaxis_transform(), ha="right", va="center", fontsize=8.5, color=INK2)
ax.text(0.105, -0.175, "stabilizing", transform=ax.get_xaxis_transform(), ha="left", va="center", fontsize=8.5, color=INK2)
fig.savefig("fig/tutorial_energy.png")
plt.close(fig)

# ---------------------------------------------------------------- figure 5: mode on the torus
s3 = h5py.File("r/post3d/Data/1_sol_X_RE.h5", "r")
X3, Y3, Z3, V3 = (s3[k][-1] for k in ("X_1", "Y_1", "Z_1", "var_1"))   # outer surface, (zeta, theta)
lim3 = abs(V3).max()
norm3 = TwoSlopeNorm(0, -lim3, lim3)
fig = plt.figure(figsize=(7.6, 4.4))
ax = fig.add_axes([0, -0.06, 1, 1.0], projection="3d")
ax.plot_surface(X3, Y3, Z3, facecolors=DIVERGING(norm3(V3)), rstride=1, cstride=1,
                linewidth=0, antialiased=False, shade=False)
L3 = max(abs(X3).max(), abs(Y3).max()); Lz = abs(Z3).max()
ax.set_xlim(-L3, L3); ax.set_ylim(-L3, L3); ax.set_zlim(-Lz, Lz)
ax.set_box_aspect((L3, L3, Lz), zoom=1.35)
ax.view_init(elev=38, azim=-62)
ax.set_axis_off()
fig.text(0.02, 0.95, "$X = \\vec{\\xi}\\cdot\\nabla\\psi$ on the plasma boundary, $n = 10$ (real part)", fontsize=10.5,
         color=INK, va="top")
with plt.rc_context({"savefig.bbox": "standard"}):   # the tight box of a 3-D axes cuts the torus
    fig.savefig("fig/tutorial_mode_3d.png", dpi=200)
plt.close(fig)
print("done", sorted(os.listdir("fig")))
print("free: lambda", ev, "E_pot/E_kin", ratio, "sum", vals.sum(), "| fixed:", ev_f, ratio_f, vals_f.sum())
print("free terms", np.round(vals, 4)); print("fixed terms", np.round(vals_f, 4))

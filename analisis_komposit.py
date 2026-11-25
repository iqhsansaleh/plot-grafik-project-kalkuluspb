"""Visualisasi 3D dan validasi komputasi analisis struktur sayap komposit.

Mencakup:
1. Analisis turunan parsial, gradien, turunan berarah, dan bidang singgung t(x, y).
2. Uji titik kritis dengan determinan Hessian D dan optimasi Lagrange P(x, y, z)
   pada permukaan elipsoida.
3. Integral lipat dalam koordinat polar untuk volume V, massa M, dan pusat massa.
4. Plot 3D Matplotlib (permukaan, bidang singgung, peta kontur gradien,
   elipsoida tekanan aerodinamika, dan kerapatan massa) serta verifikasi numerik.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
from scipy import integrate, optimize
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401


def thickness(x: np.ndarray | float, y: np.ndarray | float) -> np.ndarray | float:
    """Fungsi ketebalan sayap komposit t(x, y) dalam cm."""
    return 15.0 - 0.8 * x**2 - 0.8 * y**2 + 0.3 * x * y


def grad_thickness(x: np.ndarray | float, y: np.ndarray | float) -> tuple[np.ndarray | float, np.ndarray | float]:
    """Turunan parsial pertama (t_x, t_y) dari fungsi ketebalan."""
    tx = -1.6 * x + 0.3 * y
    ty = -1.6 * y + 0.3 * x
    return tx, ty


def density(x: np.ndarray | float, y: np.ndarray | float) -> np.ndarray | float:
    """Kerapatan variabel material komposit rho(x, y) dalam kg/m^2."""
    return 12.0 + 2.0 * np.sqrt(x**2 + y**2)


def pressure(x: np.ndarray | float, y: np.ndarray | float, z: np.ndarray | float) -> np.ndarray | float:
    """Fungsi tekanan aerodinamika P(x, y, z)."""
    return 3.0 * x - 2.0 * y + 4.0 * z


def run_numerical_verification() -> dict[str, float]:
    """Jalankan perhitungan analitik dan verifikasi numerik seluruh soal."""
    print("=" * 72)
    print("1. ANALISIS TURUNAN PADA t(x, y) = 15 - 0.8x^2 - 0.8y^2 + 0.3xy")
    print("=" * 72)
    x0, y0 = 1.0, 2.0
    z0 = float(thickness(x0, y0))
    tx0, ty0 = grad_thickness(x0, y0)
    u = np.array([-3.0, 4.0])
    u_hat = u / np.linalg.norm(u)
    dir_deriv = float(np.dot([tx0, ty0], u_hat))

    # Verifikasi turunan berarah dengan beda hingga terpusat
    eps = 1e-7
    dir_deriv_fd = (
        float(thickness(x0 + eps * u_hat[0], y0 + eps * u_hat[1]))
        - float(thickness(x0 - eps * u_hat[0], y0 - eps * u_hat[1]))
    ) / (2.0 * eps)

    print(f"t(1, 2)                 = {z0:.6f} cm")
    print(f"nabla t(1, 2)           = <{tx0:.6f}, {ty0:.6f}>")
    print(f"Vektor satuan u_hat     = <{u_hat[0]:.6f}, {u_hat[1]:.6f}>")
    print(f"D_u t(1, 2) [analitik]  = {dir_deriv:.6f}")
    print(f"D_u t(1, 2) [numerik]   = {dir_deriv_fd:.6f}")
    d_const = z0 - tx0 * x0 - ty0 * y0
    print(f"Bidang singgung         : z = {d_const:.2f} + ({tx0:.2f})x + ({ty0:.2f})y")

    print("\n" + "=" * 72)
    print("2. OPTIMASI HESSIAN & PENGALI LAGRANGE")
    print("=" * 72)
    txx, tyy, txy = -1.6, -1.6, 0.3
    hess_d = txx * tyy - txy**2
    print(f"Titik kritis interior   = (0, 0) -> di luar anulus R (r^2 = 0 < 1)")
    print(f"Hessian D               = {hess_d:.4f} > 0, t_xx = {txx:.4f} < 0 (Maksimum lokal R^2)")
    print(f"t(0, 0)                 = {thickness(0.0, 0.0):.4f} cm")

    t_max_R = 15.0 - 1.0**2 * (0.8 - 0.15)
    t_min_R = 15.0 - 3.0**2 * (0.8 + 0.15)
    print(f"Pada batas R: t_max     = {t_max_R:.4f} cm (r=1, theta=pi/4, 5pi/4)")
    print(f"Pada batas R: t_min     = {t_min_R:.4f} cm (r=3, theta=3pi/4, 7pi/4)")

    p_exact_max = 3.0 * np.sqrt(35.0)
    pt_max_exact = np.array([3.0, -1.0, 1.0]) * (np.sqrt(35.0) / 5.0)
    pt_min_exact = -pt_max_exact

    cons = {"type": "eq", "fun": lambda v: v[0] ** 2 + 2 * v[1] ** 2 + 4 * v[2] ** 2 - 21.0}
    res_max = optimize.minimize(lambda v: -pressure(v[0], v[1], v[2]), [1.0, -1.0, 1.0], constraints=cons, tol=1e-12)
    res_min = optimize.minimize(lambda v: pressure(v[0], v[1], v[2]), [-1.0, 1.0, -1.0], constraints=cons, tol=1e-12)

    print(f"P_max [analitik]        = {p_exact_max:.8f} di {np.round(pt_max_exact, 6)}")
    print(f"P_max [scipy]           = {-res_max.fun:.8f} di {np.round(res_max.x, 6)}")
    print(f"P_min [analitik]        = {-p_exact_max:.8f} di {np.round(pt_min_exact, 6)}")
    print(f"P_min [scipy]           = {res_min.fun:.8f} di {np.round(res_min.x, 6)}")

    print("\n" + "=" * 72)
    print("3. INTEGRAL LIPAT & VERIFIKASI NUMERIK (ANULUS 1 <= r <= 3)")
    print("=" * 72)
    v_exact = 88.0 * np.pi
    m_exact = (392.0 / 3.0) * np.pi

    v_num, v_err = integrate.dblquad(
        lambda r, theta: thickness(r * np.cos(theta), r * np.sin(theta)) * r,
        0.0,
        2.0 * np.pi,
        1.0,
        3.0,
    )
    m_num, m_err = integrate.dblquad(
        lambda r, theta: density(r * np.cos(theta), r * np.sin(theta)) * r,
        0.0,
        2.0 * np.pi,
        1.0,
        3.0,
    )
    mx_num, _ = integrate.dblquad(
        lambda r, theta: (r * np.cos(theta)) * density(r * np.cos(theta), r * np.sin(theta)) * r,
        0.0,
        2.0 * np.pi,
        1.0,
        3.0,
    )
    my_num, _ = integrate.dblquad(
        lambda r, theta: (r * np.sin(theta)) * density(r * np.cos(theta), r * np.sin(theta)) * r,
        0.0,
        2.0 * np.pi,
        1.0,
        3.0,
    )
    x_bar_num = mx_num / m_num
    y_bar_num = my_num / m_num

    print(f"Volume V [analitik]     = 88*pi = {v_exact:.8f} (satuan luas * cm)")
    print(f"Volume V [dblquad]      = {v_num:.8f} (est. error = {v_err:.2e})")
    print(f"Massa  M [analitik]     = (392/3)*pi = {m_exact:.8f} kg")
    print(f"Massa  M [dblquad]      = {m_num:.8f} kg (est. error = {m_err:.2e})")
    print("Pusat massa [analitik]  = (0, 0)")
    print(f"Pusat massa [dblquad]   = ({x_bar_num:.2e}, {y_bar_num:.2e})")

    return {
        "z0": z0,
        "tx0": float(tx0),
        "ty0": float(ty0),
        "dir_deriv": dir_deriv,
        "p_max": float(p_exact_max),
        "v_exact": float(v_exact),
        "v_num": float(v_num),
        "m_exact": float(m_exact),
        "m_num": float(m_num),
    }


def create_thickness_figure(fig: Figure | None = None, dpi: int = 120) -> tuple[Figure, Axes3D]:
    """Bangun figure 3D permukaan t(x, y), bidang singgung, & peta kontur gradien."""
    r = np.linspace(1.0, 3.0, 80)
    theta = np.linspace(0.0, 2.0 * np.pi, 160)
    R_grid, T_grid = np.meshgrid(r, theta)
    X = R_grid * np.cos(T_grid)
    Y = R_grid * np.sin(T_grid)
    Z = thickness(X, Y)

    x0, y0 = 1.0, 2.0
    z0 = float(thickness(x0, y0))
    tx0, ty0 = grad_thickness(x0, y0)
    u_hat = np.array([-3.0, 4.0]) / 5.0

    if fig is None:
        fig = plt.figure(figsize=(14.5, 6.6), dpi=dpi)
    fig.suptitle(
        r"Analisis Ketebalan Sayap Komposit $t(x,y) = 15 - 0.8x^2 - 0.8y^2 + 0.3xy$ pada Anulus $1 \leq x^2+y^2 \leq 9$",
        fontsize=12.5,
        fontweight="bold",
        y=0.98,
    )

    ax1 = fig.add_subplot(1, 2, 1, projection="3d")
    surf = ax1.plot_surface(
        X, Y, Z, cmap="viridis", edgecolor="none", alpha=0.88, antialiased=True
    )

    px = np.linspace(x0 - 1.1, x0 + 1.1, 20)
    py = np.linspace(y0 - 1.1, y0 + 1.1, 20)
    PX, PY = np.meshgrid(px, py)
    PZ = z0 + tx0 * (PX - x0) + ty0 * (PY - y0)
    ax1.plot_surface(PX, PY, PZ, color="tomato", alpha=0.45, edgecolor="darkred", linewidth=0.25)

    ax1.scatter([x0], [y0], [z0], color="red", s=65, zorder=10, label=r"Titik singgung $(1, 2, 11.6)$")
    ax1.scatter([0.0], [0.0], [15.0], color="gold", edgecolor="black", s=75, marker="^", zorder=10, label=r"Titik kritis $(0, 0, 15) \notin R$")

    t_b = np.linspace(0, 2 * np.pi, 200)
    for rad, lbl in [(1.0, r"Batas $r=1, r=3$"), (3.0, None)]:
        xb, yb = rad * np.cos(t_b), rad * np.sin(t_b)
        ax1.plot(xb, yb, thickness(xb, yb), color="black", linewidth=1.3, label=lbl)

    ax1.set_title("Grafik Permukaan 3D & Bidang Singgung", fontsize=11, pad=10)
    ax1.set_xlabel("$x$")
    ax1.set_ylabel("$y$")
    ax1.set_zlabel("$t(x, y)$ (cm)", labelpad=8)
    ax1.view_init(elev=28, azim=42)
    ax1.legend(loc="upper left", fontsize=8.5, framealpha=0.9)
    fig.colorbar(surf, ax=ax1, shrink=0.62, pad=0.08, label="$t(x, y)$ (cm)")

    ax2 = fig.add_subplot(1, 2, 2)
    levels = np.linspace(6.4, 14.5, 18)
    cf = ax2.contourf(X, Y, Z, levels=levels, cmap="viridis", alpha=0.85)
    cs = ax2.contour(X, Y, Z, levels=levels, colors="white", linewidths=0.6, alpha=0.7)
    ax2.clabel(cs, inline=True, fontsize=7.5, fmt="%.1f")

    rq = np.linspace(1.2, 2.8, 4)
    tq = np.linspace(0.0, 2.0 * np.pi, 16, endpoint=False)
    RQ, TQ = np.meshgrid(rq, tq)
    XQ, YQ = RQ * np.cos(TQ), RQ * np.sin(TQ)
    TXQ, TYQ = grad_thickness(XQ, YQ)
    ax2.quiver(XQ, YQ, TXQ, TYQ, color="black", alpha=0.55, scale=45, width=0.004, label=r"Medan $\nabla t(x,y)$")

    ax2.scatter([x0], [y0], color="red", s=55, zorder=10)
    ax2.quiver(
        [x0], [y0], [tx0], [ty0],
        angles="xy", scale_units="xy", scale=3.0, color="red", width=0.008, zorder=11,
        label=r"$\nabla t(1,2) = \langle -1, -2.9 \rangle$ (skala 1/3)",
    )
    ax2.quiver(
        [x0], [y0], [u_hat[0]], [u_hat[1]],
        angles="xy", scale_units="xy", scale=1.0, color="cyan", edgecolor="black", linewidth=0.6, width=0.008, zorder=11,
        label=r"$\hat{\mathbf{u}} = \langle -0.6, 0.8 \rangle$",
    )

    ax2.plot(np.cos(t_b), np.sin(t_b), "k-", linewidth=1.6)
    ax2.plot(3 * np.cos(t_b), 3 * np.sin(t_b), "k-", linewidth=1.6)
    ax2.scatter([0], [0], color="gold", edgecolor="black", marker="^", s=70, zorder=10, label=r"Titik kritis $(0,0)$")

    ax2.set_title(r"Peta Kontur & Vektor Gradien pada $R$", fontsize=11, pad=10)
    ax2.set_xlabel("$x$")
    ax2.set_ylabel("$y$")
    ax2.set_aspect("equal")
    ax2.set_xlim(-3.4, 3.4)
    ax2.set_ylim(-3.6, 3.6)
    ax2.grid(True, linestyle="--", alpha=0.35)
    ax2.legend(loc="upper right", fontsize=8, framealpha=0.94)
    fig.colorbar(cf, ax=ax2, shrink=0.82, label="$t(x, y)$ (cm)")

    fig.tight_layout(rect=[0.02, 0, 1, 0.94])
    return fig, ax1


def create_ellipsoid_figure(fig: Figure | None = None, dpi: int = 120) -> tuple[Figure, Axes3D]:
    """Bangun figure 3D elipsoida x^2 + 2y^2 + 4z^2 = 21 dengan warna tekanan P(x, y, z)."""
    u = np.linspace(0.0, 2.0 * np.pi, 100)
    v = np.linspace(0.0, np.pi, 60)
    U, V = np.meshgrid(u, v)

    a = np.sqrt(21.0)
    b = np.sqrt(21.0 / 2.0)
    c = np.sqrt(21.0 / 4.0)

    X = a * np.cos(U) * np.sin(V)
    Y = b * np.sin(U) * np.sin(V)
    Z = c * np.cos(V)
    P = pressure(X, Y, Z)

    pt_max = np.array([3.0, -1.0, 1.0]) * (np.sqrt(35.0) / 5.0)
    pt_min = -pt_max
    p_max = 3.0 * np.sqrt(35.0)

    if fig is None:
        fig = plt.figure(figsize=(10.5, 7.2), dpi=dpi)
    ax = fig.add_subplot(1, 1, 1, projection="3d")

    norm = colors.Normalize(vmin=-p_max, vmax=p_max)
    face_colors = cm.coolwarm(norm(P))
    ax.plot_surface(
        X, Y, Z, facecolors=face_colors, rstride=1, cstride=1, alpha=0.68, antialiased=True, linewidth=0
    )

    mappable = cm.ScalarMappable(norm=norm, cmap="coolwarm")
    mappable.set_array(P)
    fig.colorbar(mappable, ax=ax, shrink=0.72, pad=0.08, label=r"Tekanan $P(x,y,z) = 3x - 2y + 4z$")

    ax.scatter(
        [pt_max[0]], [pt_max[1]], [pt_max[2]],
        color="gold", edgecolor="black", linewidth=1.2, s=130, zorder=15, depthshade=False,
        label=rf"$P_{{\max}} = 3\sqrt{{35}} \approx {p_max:.3f}$ di $({pt_max[0]:.2f}, {pt_max[1]:.2f}, {pt_max[2]:.2f})$",
    )
    ax.scatter(
        [pt_min[0]], [pt_min[1]], [pt_min[2]],
        color="cyan", edgecolor="black", linewidth=1.2, s=130, zorder=15, depthshade=False,
        label=rf"$P_{{\min}} = -3\sqrt{{35}} \approx {-p_max:.3f}$ di $({pt_min[0]:.2f}, {pt_min[1]:.2f}, {pt_min[2]:.2f})$",
    )

    grad_p = np.array([3.0, -2.0, 4.0])
    grad_p_unit = grad_p / np.linalg.norm(grad_p) * 1.8
    ax.quiver(
        pt_max[0], pt_max[1], pt_max[2],
        grad_p_unit[0], grad_p_unit[1], grad_p_unit[2],
        color="black", linewidth=2.0, label=r"Arah $\nabla P = \langle 3, -2, 4 \rangle$",
    )

    ax.set_title(
        r"Optimasi Lagrange Tekanan Aerodinamika $P(x,y,z)$ pada Elipsoida $x^2 + 2y^2 + 4z^2 = 21$",
        fontsize=12,
        fontweight="bold",
        pad=14,
    )
    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")
    ax.set_zlabel("$z$")
    ax.set_box_aspect([a, b, c])
    ax.view_init(elev=22, azim=35)
    ax.legend(loc="upper left", fontsize=9, framealpha=0.92)

    fig.tight_layout()
    return fig, ax


def create_density_figure(fig: Figure | None = None, dpi: int = 120) -> tuple[Figure, Axes3D]:
    """Bangun figure 3D kerapatan rho(x, y) dan peta kerapatan 2D + pusat massa."""
    r = np.linspace(1.0, 3.0, 80)
    theta = np.linspace(0.0, 2.0 * np.pi, 160)
    R_grid, T_grid = np.meshgrid(r, theta)
    X = R_grid * np.cos(T_grid)
    Y = R_grid * np.sin(T_grid)
    RHO = density(X, Y)

    if fig is None:
        fig = plt.figure(figsize=(14.5, 6.4), dpi=dpi)
    fig.suptitle(
        r"Kerapatan Material Komposit $\rho(x,y) = 12 + 2\sqrt{x^2+y^2}$ ($\mathrm{kg/m^2}$) & Pusat Massa pada Anulus $R$",
        fontsize=12.5,
        fontweight="bold",
        y=0.98,
    )

    ax1 = fig.add_subplot(1, 2, 1, projection="3d")
    surf = ax1.plot_surface(X, Y, RHO, cmap="plasma", edgecolor="none", alpha=0.9, antialiased=True)
    ax1.set_title(r"Permukaan 3D Kerapatan $\rho(x, y)$", fontsize=11, pad=10)
    ax1.set_xlabel("$x$")
    ax1.set_ylabel("$y$")
    ax1.set_zlabel(r"$\rho(x, y)$ ($\mathrm{kg/m^2}$)")
    ax1.view_init(elev=27, azim=40)
    fig.colorbar(surf, ax=ax1, shrink=0.62, pad=0.08, label=r"$\rho$ ($\mathrm{kg/m^2}$)")

    ax2 = fig.add_subplot(1, 2, 2)
    cf = ax2.contourf(X, Y, RHO, levels=20, cmap="plasma")
    t_b = np.linspace(0, 2 * np.pi, 200)
    ax2.plot(np.cos(t_b), np.sin(t_b), "k-", linewidth=1.6)
    ax2.plot(3 * np.cos(t_b), 3 * np.sin(t_b), "k-", linewidth=1.6)
    ax2.scatter([0], [0], color="cyan", edgecolor="black", s=110, marker="*", zorder=10, label=r"Pusat massa $(\bar{x}, \bar{y}) = (0, 0)$")
    ax2.axhline(0, color="gray", linestyle="--", linewidth=0.8, alpha=0.6)
    ax2.axvline(0, color="gray", linestyle="--", linewidth=0.8, alpha=0.6)

    ax2.set_title(r"Distribusi Simetri Radial & Pusat Massa ($M = \frac{392\pi}{3}$ kg)", fontsize=11, pad=10)
    ax2.set_xlabel("$x$")
    ax2.set_ylabel("$y$")
    ax2.set_aspect("equal")
    ax2.set_xlim(-3.4, 3.4)
    ax2.set_ylim(-3.4, 3.4)
    ax2.legend(loc="upper right", fontsize=9, framealpha=0.92)
    fig.colorbar(cf, ax=ax2, shrink=0.82, label=r"$\rho$ ($\mathrm{kg/m^2}$)")

    fig.tight_layout(rect=[0.02, 0, 1, 0.94])
    return fig, ax1


def save_and_optionally_show(show: bool) -> None:
    """Simpan seluruh gambar PNG dan tampilkan jendela interaktif jika diminta."""
    out_dir = Path(__file__).resolve().parent
    builders = [
        (out_dir / "plot_ketebalan_3d_kontur.png", create_thickness_figure),
        (out_dir / "plot_elipsoida_lagrange_3d.png", create_ellipsoid_figure),
        (out_dir / "plot_kerapatan_pusat_massa_3d.png", create_density_figure),
    ]

    print("\nFile visualisasi berhasil disimpan:")
    for path, builder in builders:
        fig, _ = builder(dpi=200)
        fig.savefig(path, dpi=220, bbox_inches="tight", pad_inches=0.2)
        print(f" - {path}")
        if not show:
            plt.close(fig)

    if show:
        plt.show()


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot grafik 3D dan validasi komputasi kalkulus sayap komposit.")
    parser.add_argument("--show", action="store_true", help="Tampilkan jendela plot interaktif Matplotlib.")
    args = parser.parse_args()

    if not args.show:
        matplotlib.use("Agg")

    run_numerical_verification()
    save_and_optionally_show(show=args.show)


if __name__ == "__main__":
    main()

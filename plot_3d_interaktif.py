"""Penampil Grafik 3D Interaktif Matplotlib (TkAgg + Tkinter).

Memungkinkan pengguna memutar ketiga grafik 3D ke berbagai arah menggunakan
mouse (seret kiri untuk putar 360 derajat, seret kanan untuk zoom),
maupun memilih sudut pandang preset dan slider elevasi/azimut secara langsung.
"""

from __future__ import annotations

import argparse
import tkinter as tk
from tkinter import ttk
import matplotlib
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

from analisis_komposit import (
    create_density_figure,
    create_ellipsoid_figure,
    create_thickness_figure,
)


def build_interactive_viewer(auto_close_ms: int | None = None) -> None:
    """Bangun dan jalankan aplikasi penampil 3D interaktif berbasis Tkinter."""
    matplotlib.use("TkAgg")

    root = tk.Tk()
    root.title("Visualisasi 3D Interaktif - Analisis Struktur Sayap Komposit")
    root.geometry("1280x780")

    header = ttk.Frame(root, padding=(10, 6))
    header.pack(fill=tk.X)
    ttk.Label(
        header,
        text=(
            "Kontrol 3D: Seret tombol KIRI mouse pada grafik 3D untuk memutar ke segala arah  |  "
            "Seret tombol KANAN untuk zoom  |  Gunakan tombol sudut pandang / slider di bawah."
        ),
        font=("Segoe UI", 10, "bold"),
    ).pack(side=tk.LEFT)

    notebook = ttk.Notebook(root)
    notebook.pack(fill=tk.BOTH, expand=True, padx=8, pady=(0, 6))

    tab_specs = [
        ("1. Ketebalan t(x, y) & Gradien", create_thickness_figure, (28.0, 42.0)),
        ("2. Elipsoida Lagrange P(x, y, z)", create_ellipsoid_figure, (22.0, 35.0)),
        ("3. Kerapatan rho(x, y) & Pusat Massa", create_density_figure, (27.0, 40.0)),
    ]

    for title, builder, (default_elev, default_azim) in tab_specs:
        tab = ttk.Frame(notebook)
        notebook.add(tab, text=title)

        fig = Figure(figsize=(12.2, 6.2), dpi=105)
        _, ax3d = builder(fig=fig)

        canvas = FigureCanvasTkAgg(fig, master=tab)
        canvas_widget = canvas.get_tk_widget()
        canvas_widget.pack(fill=tk.BOTH, expand=True)

        # Pastikan rotasi mouse 3D aktif penuh pada TkAgg
        ax3d.mouse_init(rotate_btn=1, zoom_btn=3, pan_btn=2)
        canvas.draw()

        controls = ttk.Frame(tab, padding=(8, 4))
        controls.pack(fill=tk.X)

        toolbar = NavigationToolbar2Tk(canvas, controls, pack_toolbar=False)
        toolbar.update()
        toolbar.pack(side=tk.LEFT)

        elev_var = tk.DoubleVar(value=default_elev)
        azim_var = tk.DoubleVar(value=default_azim)

        def make_sync_from_mouse(ax=ax3d, ev=elev_var, av=azim_var):
            def sync_view(_event: object = None) -> None:
                ev.set(round(float(ax.elev), 1))
                av.set(round(float(ax.azim), 1))

            return sync_view

        canvas.mpl_connect("button_release_event", make_sync_from_mouse())

        def make_updater(ax=ax3d, c=canvas, ev=elev_var, av=azim_var):
            def update_view(_event: object = None) -> None:
                ax.view_init(elev=ev.get(), azim=av.get())
                c.draw_idle()

            return update_view

        updater = make_updater()

        def make_preset(elev: float, azim: float, ax=ax3d, c=canvas, ev=elev_var, av=azim_var):
            def apply_preset() -> None:
                ev.set(elev)
                av.set(azim)
                ax.view_init(elev=elev, azim=azim)
                c.draw_idle()

            return apply_preset

        btn_frame = ttk.Frame(controls)
        btn_frame.pack(side=tk.RIGHT)

        ttk.Label(btn_frame, text="Elevasi:").pack(side=tk.LEFT, padx=(6, 2))
        ttk.Scale(
            btn_frame, from_=-90, to=90, variable=elev_var, orient=tk.HORIZONTAL, length=110, command=updater
        ).pack(side=tk.LEFT)

        ttk.Label(btn_frame, text="Azimut:").pack(side=tk.LEFT, padx=(8, 2))
        ttk.Scale(
            btn_frame, from_=-180, to=180, variable=azim_var, orient=tk.HORIZONTAL, length=130, command=updater
        ).pack(side=tk.LEFT)

        presets = [
            ("Reset", default_elev, default_azim),
            ("Atas (XY)", 90.0, -90.0),
            ("Depan (XZ)", 0.0, -90.0),
            ("Samping (YZ)", 0.0, 0.0),
        ]
        for label, p_elev, p_azim in presets:
            ttk.Button(
                btn_frame,
                text=label,
                command=make_preset(p_elev, p_azim),
            ).pack(side=tk.LEFT, padx=3)

    if auto_close_ms is not None:
        root.after(auto_close_ms, root.destroy)

    root.mainloop()


def main() -> None:
    parser = argparse.ArgumentParser(description="Buka penampil grafik 3D interaktif Matplotlib.")
    parser.add_argument(
        "--test-ms",
        type=int,
        default=None,
        help="Tutup jendela otomatis setelah N milidetik (untuk pengujian).",
    )
    args = parser.parse_args()
    build_interactive_viewer(auto_close_ms=args.test_ms)


if __name__ == "__main__":
    main()

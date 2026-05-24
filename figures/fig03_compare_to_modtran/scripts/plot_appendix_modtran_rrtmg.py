"""
Generate the radiation-comparison figure (MODTRAN vs original/modified RRTMG).

NOTE on the filename: "appendix" is historical -- this figure originated in the
supplementary material and was later moved to the main text. The name is kept
for traceability with the source repository (transport-paper fork,
layer_toy_model/scripts/plot_appendix_modtran_rrtmg.py).

Reads a cached rf_results_day{DDD}.npz (one per sampled day, produced upstream by
comparison_to_modtran.py, which ran climlab + parsed MODTRAN TAPE6 files) and
writes a single-panel PNG of daily-integrated SW radiative forcing per unit
burden, as a function of latitude, comparing:

    - MODTRAN, hourly-sampled and daily-integrated   (black solid)   -- "3D" / ground truth
    - original RRTMG, hourly-sampled, daily-integrated (green solid)  -- "3D"
    - original RRTMG, single call with <mu>_daily      (green dashed) -- "2D"
    - modified RRTMG, single call with <mu>_daily      (purple dashed)-- "2D"

("3D" = the way a radiation code is driven inside a 3D GCM, i.e. one call per
hour with longitudinal/diurnal sampling; "2D" = the way a zonally-averaged 2D
transport model must drive it, i.e. one call per latitude with the daily-averaged
cosine of the solar zenith angle as the only solar-geometry input. The figure
also has access to a "modified RRTMG, 3D" curve in the cached data; it is omitted
from the plot because it is visually almost coincident with MODTRAN.)

Scenario (encoded in the cached .npz): silica layer, r = 250 nm (d = 500 nm),
10 Tg uniformly distributed in p in [23, 83] hPa (~18-22 km), annual-mean
background atmosphere, clear sky.

Dependencies: numpy, matplotlib only (see requirements.txt).

Run from anywhere:
    python scripts/plot_appendix_modtran_rrtmg.py [--day 91] [--material silica]
                                                  [--ymin -0.85] [--ymax 0.05]
                                                  [--outdir DIR]

Default output: <bundle>/figures/appendix_modtran_rrtmg_day{DDD}.png

Available cached days (silica): 21 52 81 91 111 141 172 202 233 264 294 326 356
The figure used in the paper is day 91.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from _style import apply
apply()

M_TG = 10.0  # total layer mass [Tg] used when producing the cached .npz


def daily_integrate(mat_extended: np.ndarray, hours_mat_extended: np.ndarray) -> np.ndarray:
    """Daily-average a (lat, hour) quantity by trapezoidal integration over valid hours,
    normalized by the full 24-hour (2*pi rad) cycle so that night-time zeros are counted.
    Matches the convention used to produce the cached .npz."""
    nlat = mat_extended.shape[0]
    out = np.zeros(nlat)
    two_pi = 2 * np.pi
    for ilat in range(nlat):
        ind = np.where(~np.isnan(hours_mat_extended[ilat, :]))[0]
        if len(ind) >= 2:
            out[ilat] = np.trapezoid(mat_extended[ilat, ind], hours_mat_extended[ilat, ind]) / two_pi
    return out


def load_day(material: str, day: int) -> dict:
    here = Path(__file__).resolve().parent
    npz_path = here.parent / "data" / material / f"rf_results_day{day:03d}.npz"
    if not npz_path.exists():
        raise FileNotFoundError(
            f"cached data not found: {npz_path}\n"
            f"(available materials/days are under {here.parent / 'data'})"
        )
    return dict(np.load(npz_path))


def build_figure(material: str, day: int, out_dir: Path, lat_lim=(-85, 85), ylim=None) -> Path:
    d = load_day(material, day)

    lat = d["lat_vect"]
    hours = d["hours_mat_extended"]

    rf_modtran = daily_integrate(d["rf_mat_modtran_extended"], hours) / M_TG
    rf_old_inst = daily_integrate(d["model_rf_mat_old_rrtmg_extended"], hours) / M_TG
    rf_new_inst = daily_integrate(d["model_rf_mat_extended"], hours) / M_TG  # modified RRTMG, 3D (not plotted)

    rf_old_avg = d["model_rf_mat_old_rrtmg_daily_avg"] / M_TG
    rf_new_avg = d["model_rf_mat_new_rrtmg_daily_avg"] / M_TG

    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.tick_params(labelsize=11)
    ax.plot(lat, rf_modtran, "k-", lw=2.0, label="MODTRAN")
    ax.plot(lat, rf_old_inst, color="#1F77B4", ls="-", lw=1.5,
            label="original RRTMG 3D")
    ax.plot(lat, rf_old_avg, color="#1F77B4", ls="--", lw=1.5,
            label="original RRTMG 2D")
    # ax.plot(lat, rf_new_inst, color="#C44E52", ls="-", lw=1.5,
    #         label="modified RRTMG 3D")   # omitted: visually almost coincident with MODTRAN
    ax.plot(lat, rf_new_avg, color="#C44E52", ls="--", lw=1.5,
            label="modified RRTMG 2D")

    ax.axhline(0.0, color="0.6", lw=0.8)
    ax.set_xlabel("Latitude [deg]")
    ax.set_ylabel(r"$RF_{\mathrm{SW}}$  [$W\,m^{-2}\,Tg^{-1}$]")
    ax.set_xlim(lat_lim)
    if ylim:
        ax.set_ylim(ylim)
    ax.grid(True, alpha=0.4)
    ax.legend(loc="upper center", fontsize=11, framealpha=0.9)

    out_dir.mkdir(parents=True, exist_ok=True)
    out_path_png = out_dir / f"appendix_modtran_rrtmg_day{day:03d}.png"
    out_path_pdf = out_dir / f"appendix_modtran_rrtmg_day{day:03d}.pdf"
    fig.tight_layout()
    fig.savefig(out_path_png, dpi=300)
    fig.savefig(out_path_pdf)
    plt.close(fig)
    return out_path_png


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--day", type=int, default=91, help="calendar day of the cached run (default 91)")
    p.add_argument("--material", default="silica", help="aerosol material subfolder under data/ (default silica)")
    p.add_argument("--ymin", type=float, default=-0.85, help="y-axis lower limit (default -0.85)")
    p.add_argument("--ymax", type=float, default=0.05, help="y-axis upper limit (default 0.05)")
    p.add_argument("--outdir", default=None, help="output directory (default: <bundle>/figures)")
    args = p.parse_args()

    here = Path(__file__).resolve().parent
    out_dir = Path(args.outdir).expanduser().resolve() if args.outdir else (here.parent / "figures")
    out = build_figure(args.material, args.day, out_dir, ylim=(args.ymin, args.ymax))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()

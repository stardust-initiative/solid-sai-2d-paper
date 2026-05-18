#!/usr/bin/env python
r"""Preprocess transport-run aerosol profiles into SARF-ready NetCDFs.

The 2-D transport runner writes a per-injection-point ``.npz`` carrying a
``total_tracers`` field of shape ``(n_bins, n_snapshots, lat, lev)``. The
SARF case driver, however, wants a flat per-bin mass-mixing-ratio field on
the ``(lat, lev)`` grid. This script performs that reduction.

Two modes
---------
``--mode annual`` averages the last 12 monthly snapshots (one year) of
each run into a single annual-mean profile, written as
``<lat>deg_<plev>hpa.nc``.

``--mode seasonal`` splits year 5 of a 5-year run into the four
climatological seasons -- DJF/MAM/JJA/SON -- by averaging the three
matching monthly snapshots, written as
``<lat>deg_<plev>hpa_<season>.nc``. The four seasonal means average back
to the annual mean, so seasonal SARF stays directly comparable with the
annual sweep.

Both modes write, per transport-run directory:

  * ``radius_mapping.npz``  -- maps each ``bin_<i>`` variable to its
    effective particle radius (= equivalent-volume diameter / 2)
  * one ``.nc`` per injection point (per season in seasonal mode)

and a global ``manifest.csv`` consumed by ``run_sarf_sweep``.

Rate normalisation: every output profile is normalised to the canonical
20 Tg/yr total injection. The original single-injection nocoag runs were
performed at 10 Tg/yr and are scaled by 2.0; coag runs and the symmetric
pair-study runs were already at 20 Tg/yr (factor 1.0).

Two ways to use it
-------------------
Command line::

    conda run -n climlab_stardust_ext_env python scripts/preprocess_profiles.py \
        --mode annual --material silica

As a Python module::

    from preprocess_profiles import PreprocessConfig, preprocess_profiles
    cfg = PreprocessConfig(mode='seasonal', dirs='symmetric')
    rows = preprocess_profiles(cfg)
"""

import argparse
import csv
import os
import re
from dataclasses import dataclass
from glob import glob

import numpy as np
import xarray as xr


# Repo layout: this file lives in transport-paper-umbrella/scripts/.
_SCRIPTS_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.dirname(_SCRIPTS_DIR)

# Material -> bulk density [kg/m^3], for the equivalent-sphere mass
# conversion (mass per particle = 4/3 pi r^3 rho).
RHO = {
    'silica':  2200.0,
    'sulfate': 1770.0,
    'calcite': 2710.0,
}

# Seasonal snapshot indices from year 5 of a 5-year monthly-output run
# starting Jan 1 (61 snapshots: index 0 is the initial state, 1..60 are
# end-of-month for months 1..60). The 12 indices below cover every
# snapshot of year 5 exactly once.
SEASON_SNAPSHOTS = {
    'DJF': [49, 50, 60],   # Jan, Feb, Dec of year 5
    'MAM': [51, 52, 53],
    'JJA': [54, 55, 56],
    'SON': [57, 58, 59],
}

# Transport-run directory -> (output label, material, monomer diameter
# [nm], coag flag, injection schedule, rate factor).
#
# rate factor normalises each run to the canonical 20 Tg/yr: the original
# single-injection nocoag runs were 10 Tg/yr (factor 2.0); coag runs and
# the symmetric pair-study runs were already 20 Tg/yr (factor 1.0).
DIR_SPECS = {
    # ---- symmetric pair study (20 Tg/yr, both schedules) ----
    'for_paper_2026_04_diameter05_symmetric':
        ('silica_d500_nocoag_sym', 'silica', 500, False, 'sym', 1.0),
    'for_paper_2026_04_diameter05_symmetric_seasons':
        ('silica_d500_nocoag_symseas', 'silica', 500, False, 'symseas', 1.0),
    'for_paper_2026_04_diameter05_symmetric_20Tg_coag':
        ('silica_d500_coag_sym', 'silica', 500, True, 'sym', 1.0),
    'for_paper_2026_04_diameter05_symmetric_seasons_20Tg_coag':
        ('silica_d500_coag_symseas', 'silica', 500, True, 'symseas', 1.0),
    'for_paper_2026_04_diameter03_rho2710_symmetric':
        ('calcite_d300_nocoag_sym', 'calcite', 300, False, 'sym', 1.0),
    'for_paper_2026_04_diameter03_rho2710_symmetric_seasons':
        ('calcite_d300_nocoag_symseas', 'calcite', 300, False, 'symseas', 1.0),
    'for_paper_2026_04_diameter03_rho2710_symmetric_20Tg_coag':
        ('calcite_d300_coag_sym', 'calcite', 300, True, 'sym', 1.0),
    'for_paper_2026_04_diameter03_rho2710_symmetric_seasons_20Tg_coag':
        ('calcite_d300_coag_symseas', 'calcite', 300, True, 'symseas', 1.0),
    # ---- original single-injection sweep ----
    # silica (no _rho suffix)
    'for_paper_2026_04_diameter03':
        ('silica_d300_nocoag_single', 'silica', 300, False, 'single', 2.0),
    'for_paper_2026_04_diameter05':
        ('silica_d500_nocoag_single', 'silica', 500, False, 'single', 2.0),
    'for_paper_2026_04_diameter07':
        ('silica_d700_nocoag_single', 'silica', 700, False, 'single', 2.0),
    'for_paper_2026_04_diameter03_20Tg_coag':
        ('silica_d300_coag_single', 'silica', 300, True, 'single', 1.0),
    'for_paper_2026_04_diameter05_20Tg_coag':
        ('silica_d500_coag_single', 'silica', 500, True, 'single', 1.0),
    'for_paper_2026_04_diameter07_20Tg_coag':
        ('silica_d700_coag_single', 'silica', 700, True, 'single', 1.0),
    # sulfate (rho1770; nocoag only)
    'for_paper_2026_04_diameter03_rho1770':
        ('sulfate_d300_nocoag_single', 'sulfate', 300, False, 'single', 2.0),
    'for_paper_2026_04_diameter05_rho1770':
        ('sulfate_d500_nocoag_single', 'sulfate', 500, False, 'single', 2.0),
    'for_paper_2026_04_diameter07_rho1770':
        ('sulfate_d700_nocoag_single', 'sulfate', 700, False, 'single', 2.0),
    # calcite (rho2710)
    'for_paper_2026_04_diameter03_rho2710':
        ('calcite_d300_nocoag_single', 'calcite', 300, False, 'single', 2.0),
    'for_paper_2026_04_diameter05_rho2710':
        ('calcite_d500_nocoag_single', 'calcite', 500, False, 'single', 2.0),
    'for_paper_2026_04_diameter07_rho2710':
        ('calcite_d700_nocoag_single', 'calcite', 700, False, 'single', 2.0),
    'for_paper_2026_04_diameter03_rho2710_20Tg_coag':
        ('calcite_d300_coag_single', 'calcite', 300, True, 'single', 1.0),
    'for_paper_2026_04_diameter05_rho2710_20Tg_coag':
        ('calcite_d500_coag_single', 'calcite', 500, True, 'single', 1.0),
    'for_paper_2026_04_diameter07_rho2710_20Tg_coag':
        ('calcite_d700_coag_single', 'calcite', 700, True, 'single', 1.0),
}


# ===========================================================================
# Configuration
# ===========================================================================

@dataclass
class PreprocessConfig:
    """Configuration for one preprocessing pass.

    ``mode`` selects annual-mean vs four-season reduction. ``dirs`` /
    ``material`` select which transport-run directories are processed.
    """
    mode: str = 'annual'              # 'annual' | 'seasonal'
    input_root: str = os.path.join(_REPO_ROOT, 'data', 'transport_runs')
    out_dir: str = None               # default depends on mode (see resolve)
    # annual mode: comma-separated material list, or 'all'
    material: str = 'silica'
    # seasonal mode: 'all' | 'symmetric' | 'single' | comma-list of dir names
    dirs: str = 'single'
    dry_run: bool = False

    def resolve_out_dir(self):
        """Return the output directory, defaulting on ``mode``."""
        if self.out_dir is not None:
            return self.out_dir
        sub = 'preprocessed' if self.mode == 'annual' else 'preprocessed_seasonal'
        return os.path.join(_REPO_ROOT, 'data', sub)


def parse_config():
    """Parse command-line arguments into a :class:`PreprocessConfig`."""
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    d = PreprocessConfig()
    p.add_argument('--mode', choices=['annual', 'seasonal'], default=d.mode,
                   help="'annual': one year-mean profile per run; "
                        "'seasonal': four climatological-season profiles")
    p.add_argument('--input-root', default=d.input_root,
                   help='directory holding the transport-run subdirectories')
    p.add_argument('--out-dir', default=None,
                   help='output directory (default: data/preprocessed[_seasonal])')
    p.add_argument('--material', default=d.material,
                   help='annual mode: comma-separated material list '
                        '(silica,sulfate,calcite) or "all"')
    p.add_argument('--dirs', default=d.dirs,
                   help='seasonal mode: "all" | "symmetric" | "single" | '
                        'comma-separated explicit directory names')
    p.add_argument('--dry-run', action='store_true',
                   help='enumerate outputs without reading npz or writing nc')
    a = p.parse_args()
    return PreprocessConfig(
        mode=a.mode, input_root=a.input_root, out_dir=a.out_dir,
        material=a.material, dirs=a.dirs, dry_run=a.dry_run,
    )


# ===========================================================================
# Helpers
# ===========================================================================

def parse_npz_filename(fname):
    """Map an npz filename's trailing tag to ``(signed_lat, plev_hPa)``.

    Matches ``..._<lat>{N|S}_<plev>hpa.npz``; works for both the original
    single-injection sweep and the symmetric pair study (the token before
    the lat -- sym / symseas / nocoag / coag -- is irrelevant here).
    """
    m = re.search(r'_(\d+)([NS])_(\d+)hpa\.npz$', fname)
    if m is None:
        return None
    sign = +1 if m.group(2) == 'N' else -1
    return sign * int(m.group(1)), int(m.group(3))


def _column_burden_Tg(mmr_sum, lat, lev):
    """Lat-weighted total column burden [Tg] of a summed-over-bins mmr field.

    ``mmr_sum`` has shape ``(lat, lev)``.
    """
    mid = 0.5 * (lev[:-1] + lev[1:])
    lev_bounds = np.concatenate([[0.0], mid, [1013.25]])
    dp_pa = 1e2 * np.diff(lev_bounds)
    col = np.sum(mmr_sum * dp_pa[None, :] / 9.80665, axis=-1)   # kg/m^2
    R = 6.371e6
    dlat = np.deg2rad(2.0)
    area_lat = 2 * np.pi * R ** 2 * dlat * np.cos(np.deg2rad(lat))
    return float((col * area_lat).sum() / 1e9)


def _bin_metadata(sample):
    """Return ``(bin_names, radius_mapping)`` from a loaded sample npz."""
    bin_diameters = np.asarray(sample['aerosol_sizes'])   # equiv-vol diameter [m]
    bin_radii = bin_diameters / 2.0
    bin_names = [f'bin_{i:02d}_{int(round(d * 1e9))}nm'
                 for i, d in enumerate(bin_diameters)]
    radius_mapping = {name: float(r) for name, r in zip(bin_names, bin_radii)}
    return bin_names, radius_mapping


# ===========================================================================
# Per-directory processing
# ===========================================================================

def _process_dir(input_dir, output_dir, material, d_nm, coag, schedule,
                 rate_factor, mode, input_root, dry_run=False):
    """Process one transport-run directory in the requested mode."""
    npz_files = sorted(glob(os.path.join(input_dir, '*.npz')))
    if not npz_files:
        print(f'  [skip] no npz files in {input_dir}')
        return []

    rho = RHO[material]
    os.makedirs(output_dir, exist_ok=True)

    bin_names, radius_mapping = _bin_metadata(
        np.load(npz_files[0], allow_pickle=True))
    if not dry_run:
        np.savez(os.path.join(output_dir, 'radius_mapping.npz'),
                 **radius_mapping)

    # In annual mode, the whole run is reduced to one snapshot group; in
    # seasonal mode, to four. A unified iteration over named groups
    # covers both.
    if mode == 'annual':
        groups = {'': None}            # '' = no season suffix; idx set below
    else:
        groups = dict(SEASON_SNAPSHOTS)

    rows = []
    for npz_path in npz_files:
        meta = parse_npz_filename(os.path.basename(npz_path))
        if meta is None:
            print(f'  [skip] cannot parse {os.path.basename(npz_path)}')
            continue
        inj_lat, inj_plev = meta
        sign = 'p' if inj_lat >= 0 else 'm'

        d = None if dry_run else np.load(npz_path, allow_pickle=True)

        for season, snap_idx in groups.items():
            tag = f'_{season}' if season else ''
            out_name = f'{sign}{abs(inj_lat)}deg_{inj_plev}hpa{tag}.nc'
            out_path = os.path.join(output_dir, out_name)
            if dry_run:
                rows.append({'material': material, 'diameter_nm': d_nm,
                             'coag': coag, 'schedule': schedule,
                             'season': season or 'Annual', 'inj_lat': inj_lat,
                             'inj_plev_hPa': inj_plev,
                             'output_path': out_path, 'status': 'dry_run'})
                continue

            tt = d['total_tracers']    # (n_bins, n_snapshots, lat, lev)
            if tt.shape[0] != len(bin_names):
                raise ValueError(
                    f'bin count mismatch in {npz_path}: '
                    f'{tt.shape[0]} vs expected {len(bin_names)}')

            if mode == 'annual':
                mmr = tt[..., -12:, :, :].mean(axis=-3)        # (n_bins, lat, lev)
            else:
                mmr = tt[:, snap_idx, :, :].mean(axis=1)       # (n_bins, lat, lev)
            mmr = rate_factor * mmr
            # clip negatives (cubic-extrapolation artifacts in the solver)
            mmr = np.where(mmr > 0.0, mmr, 0.0)

            lat = np.asarray(d['lat'])
            lev = np.asarray(d['lev'])
            data_vars = {name: (('lat', 'lev'), mmr[i])
                         for i, name in enumerate(bin_names)}
            attrs = {
                'material': material,
                'monomer_diameter_nm': d_nm,
                'coag': int(coag),
                'schedule': schedule,
                'inj_lat_deg': inj_lat,
                'inj_plev_hPa': inj_plev,
                'rate_Tg_yr_total': 20.0,
                'rate_factor_applied': rate_factor,
                'particle_density_kg_m3': rho,
                'source_npz': os.path.relpath(npz_path, input_root),
            }
            if season:
                attrs['season'] = season
                attrs['season_snapshots'] = str(snap_idx)
            else:
                attrs['avg_period_months'] = 12
            xr.Dataset(data_vars, coords={'lat': lat, 'lev': lev},
                       attrs=attrs).to_netcdf(out_path)

            burden_Tg = _column_burden_Tg(mmr.sum(axis=0), lat, lev)
            rows.append({'material': material, 'diameter_nm': d_nm,
                         'coag': coag, 'schedule': schedule,
                         'season': season or 'Annual', 'inj_lat': inj_lat,
                         'inj_plev_hPa': inj_plev, 'burden_Tg': burden_Tg,
                         'n_bins': len(bin_names),
                         'output_path': out_path, 'status': 'ok'})

    n_ok = sum(1 for r in rows if r['status'] == 'ok')
    print(f'  wrote {n_ok} files to {output_dir}')
    return rows


def _annual_dir_label(material, d_nm, coag):
    coag_lbl = 'coag' if coag else 'nocoag'
    return f'{material}_d{d_nm:03d}_{coag_lbl}'


# ===========================================================================
# Entry
# ===========================================================================

def preprocess_profiles(cfg):
    """Run one preprocessing pass and return the manifest rows.

    Parameters
    ----------
    cfg : PreprocessConfig
        The preprocessing configuration.

    Returns
    -------
    list of dict
        One row per written (or, with ``dry_run``, planned) output file.
    """
    out_dir = cfg.resolve_out_dir()
    os.makedirs(out_dir, exist_ok=True)
    manifest_rows = []

    if cfg.mode == 'annual':
        # Annual mode walks every directory matching DIR_SPECS' single-
        # injection entries, filtered by material.
        if cfg.material == 'all':
            wanted = set(RHO)
        else:
            wanted = set(cfg.material.split(','))
        for dirname, spec in DIR_SPECS.items():
            label, material, d_nm, coag, schedule, rate_factor = spec
            if schedule != 'single' or material not in wanted:
                continue
            full = os.path.join(cfg.input_root, dirname)
            if not os.path.isdir(full):
                continue
            out_sub = os.path.join(out_dir,
                                   _annual_dir_label(material, d_nm, coag))
            print(f'[{material} d={d_nm}nm '
                  f'{"coag" if coag else "nocoag"}] {dirname} -> {out_sub}')
            manifest_rows += _process_dir(
                full, out_sub, material, d_nm, coag, schedule, rate_factor,
                'annual', cfg.input_root, dry_run=cfg.dry_run)
    else:
        # Seasonal mode selects directories by --dirs.
        sym = [k for k, v in DIR_SPECS.items() if v[4] in ('sym', 'symseas')]
        single = [k for k, v in DIR_SPECS.items() if v[4] == 'single']
        if cfg.dirs == 'all':
            wanted_dirs = sym + single
        elif cfg.dirs == 'symmetric':
            wanted_dirs = sym
        elif cfg.dirs == 'single':
            wanted_dirs = single
        else:
            wanted_dirs = cfg.dirs.split(',')
        for dirname in wanted_dirs:
            if dirname not in DIR_SPECS:
                print(f'[skip] unknown directory {dirname}')
                continue
            full = os.path.join(cfg.input_root, dirname)
            if not os.path.isdir(full):
                print(f'[skip] directory does not exist: {full}')
                continue
            label, material, d_nm, coag, schedule, rate_factor = \
                DIR_SPECS[dirname]
            out_sub = os.path.join(out_dir, label)
            print(f'[{label}] {dirname} -> {out_sub}')
            manifest_rows += _process_dir(
                full, out_sub, material, d_nm, coag, schedule, rate_factor,
                'seasonal', cfg.input_root, dry_run=cfg.dry_run)

    # --- write manifest --------------------------------------------------
    if manifest_rows:
        manifest_path = os.path.join(out_dir, 'manifest.csv')
        keys = list(manifest_rows[0].keys())
        with open(manifest_path, 'w') as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            w.writerows(manifest_rows)
        print(f'\nwrote manifest with {len(manifest_rows)} rows '
              f'-> {manifest_path}')
        ok = [r for r in manifest_rows if r.get('status') == 'ok']
        if ok:
            burdens = [r['burden_Tg'] for r in ok]
            print(f'burden across {len(ok)} ok cases: '
                  f'min={min(burdens):.2f} median={np.median(burdens):.2f} '
                  f'max={max(burdens):.2f} Tg')
    return manifest_rows


def main():
    cfg = parse_config()
    preprocess_profiles(cfg)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

#!/usr/bin/env python
r"""Build a converged SARF reference model on the transport-paper data grid.

The SARF (stratospheric-aerosol radiative-forcing) pipeline computes a
forcing as the difference between a perturbed and a reference radiative
state. This script builds that reference: an ERA5-pinned, transport-free
radiative-convective column model integrated to a true radiative fixed
point, saved as a NetCDF that ``run_sarf_case`` later loads.

The reference is built on the lat/lev grid extracted from one of the
transport-run ``.npz`` files (90 latitudes, 50 log-spaced pressure
levels) so the aerosol profiles need no regridding when SARF is computed.

The model is the ``era5ref`` protocol: ``minimal_for_sarf=True`` (only the
radiative core), ``disable_transport=True``, troposphere pinned to the
ERA5 seasonal-mean T/q/Ts. It is integrated for ``--n-cycle`` spin-up
cycles of ``--t-cycle`` days each, then one long ``--t-avg`` averaging
cycle. The default ``--n-cycle`` (240, ~20 model years) is large enough
to reach a genuine fixed point even at the polar stratopause where
radiative timescales are long.

cold-n20: the published SARF results use ``n_rrtmg_repeat = 20`` (the
cloud-overlap MCICA ensemble size; the rev6 model default is 5). The
``--n-rrtmg`` value is passed straight into ``get_ref`` at construction
and recorded in the output NetCDF attrs so ``run_sarf_case`` can assert
the reference and perturbed models agree.

Two ways to use it
-------------------
Command line::

    conda run -n climlab_stardust_ext_env python scripts/build_sarf_ref.py \
        --season Annual --n-cycle 240

As a Python module::

    from build_sarf_ref import RefConfig, build_sarf_ref
    cfg = RefConfig(season='DJF', n_cycle=240)
    result = build_sarf_ref(cfg)
"""

import argparse
import json
import os
import subprocess
import time
from dataclasses import dataclass, asdict

import numpy as np
import xarray as xr

from climate_runs_ext import load_project_config
from climate_runs_ext.reference_model import get_ref
from climate_runs_ext.utils.era5_data import (
    SeasonTypes, era5_annual_initial_state,
)
from climate_runs_ext.utils.state_io import iteratively_update_internal
from climate_runs_ext.utils.model_control import fix_Ts, fix_q, fix_Tatm_trop


# Repo layout: this file lives in transport-paper-umbrella/scripts/.
_SCRIPTS_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.dirname(_SCRIPTS_DIR)


# ===========================================================================
# Configuration
# ===========================================================================

@dataclass
class RefConfig:
    """Configuration for one SARF-reference build.

    The reference is fully specified by the CO2 concentration, the
    insolation season (which also selects the ERA5 month list used for
    the troposphere pinning), the spin-up schedule, the tropopause
    pressure below which T is held to ERA5, and the RRTMG cloud-overlap
    ensemble size ``n_rrtmg``.
    """
    # grid source: a transport-run npz from which lat/lev are read
    grid_npz: str = os.path.join(
        _REPO_ROOT, 'data', 'transport_runs',
        'for_paper_2026_04_mean_winds_diameter05_nocoag_1N_68hpa.npz')
    out_dir: str = os.path.join(_REPO_ROOT, 'data', 'sarf_ref')
    # physics
    co2_ppm: float = 420.0
    season: str = 'Annual'            # Annual | DJF | MAM | JJA | SON
    p_trop_hPa: float = 175.0
    n_rrtmg: int = 20                 # RRTMG cloud-overlap ensemble size
    # spin-up schedule
    n_cycle: int = 240                # spin-up cycles (~20 model years)
    t_cycle: float = 30.0             # days per spin-up cycle
    t_avg: float = 365.0              # days in the final averaging cycle
    # behaviour
    force: bool = False               # rebuild even if the output exists


def parse_config():
    """Parse command-line arguments into a :class:`RefConfig`."""
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    d = RefConfig()
    p.add_argument('--grid-npz', default=d.grid_npz,
                   help='transport-run npz to read the lat/lev grid from')
    p.add_argument('--out-dir', default=d.out_dir,
                   help='directory for the reference NetCDF + provenance')
    p.add_argument('--co2', dest='co2_ppm', type=float, default=d.co2_ppm,
                   help='CO2 concentration [ppm]')
    p.add_argument('--season', default=d.season,
                   choices=['Annual', 'DJF', 'MAM', 'JJA', 'SON'],
                   help='insolation season + ERA5 month list; Annual writes '
                        'the legacy filename, others append a _<season> tag')
    p.add_argument('--p-trop-hPa', type=float, default=d.p_trop_hPa,
                   help='pressure below which Tatm is pinned to ERA5 [hPa]')
    p.add_argument('--n-rrtmg', type=int, default=d.n_rrtmg,
                   help='RRTMG cloud-overlap ensemble size (published SARF '
                        'uses 20; the rev6 model default is 5)')
    p.add_argument('--n-cycle', type=int, default=d.n_cycle,
                   help='number of spin-up cycles (the original converged '
                        'seasonal refs used ~240)')
    p.add_argument('--t-cycle', type=float, default=d.t_cycle,
                   help='days per spin-up cycle')
    p.add_argument('--t-avg', type=float, default=d.t_avg,
                   help='days in the final averaging cycle')
    p.add_argument('--force', action='store_true',
                   help='rebuild even if the output NetCDF already exists')
    a = p.parse_args()
    return RefConfig(
        grid_npz=a.grid_npz, out_dir=a.out_dir, co2_ppm=a.co2_ppm,
        season=a.season, p_trop_hPa=a.p_trop_hPa, n_rrtmg=a.n_rrtmg,
        n_cycle=a.n_cycle, t_cycle=a.t_cycle, t_avg=a.t_avg, force=a.force,
    )


# ===========================================================================
# Helpers
# ===========================================================================

def _git_info(path):
    """Return {sha, branch, dirty} for a git repo, or 'unknown' on failure."""
    try:
        sha = subprocess.check_output(
            ['git', '-C', path, 'rev-parse', 'HEAD'],
            stderr=subprocess.DEVNULL).decode().strip()
        branch = subprocess.check_output(
            ['git', '-C', path, 'rev-parse', '--abbrev-ref', 'HEAD'],
            stderr=subprocess.DEVNULL).decode().strip()
        dirty = bool(subprocess.check_output(
            ['git', '-C', path, 'status', '--porcelain'],
            stderr=subprocess.DEVNULL).decode().strip())
        return {'sha': sha, 'branch': branch, 'dirty': dirty}
    except Exception:
        return {'sha': 'unknown', 'branch': 'unknown', 'dirty': None}


def _lat_weighted(values, lat):
    """Cosine-latitude-weighted global mean of a 1-D lat array."""
    w = np.cos(np.deg2rad(np.asarray(lat)))
    return float((np.asarray(values).squeeze() * w).sum() / w.sum())


# ===========================================================================
# Build
# ===========================================================================

def build_sarf_ref(cfg, verbose=True):
    """Build one SARF reference model and write it to ``cfg.out_dir``.

    Parameters
    ----------
    cfg : RefConfig
        The reference configuration.
    verbose : bool
        Print a header and per-cycle convergence metrics.

    Returns
    -------
    dict
        ``ref_nc``      -- path of the written reference NetCDF
        ``log_file``    -- path of the per-cycle convergence log
        ``prov_file``   -- path of the provenance JSON
        ``n_rrtmg``     -- the RRTMG ensemble size that was used
        ``skipped``     -- True if the output already existed (no --force)
    """
    log = print if verbose else (lambda *a, **k: None)

    suffix = '' if cfg.season == 'Annual' else f'_{cfg.season}'
    os.makedirs(cfg.out_dir, exist_ok=True)
    ref_nc = os.path.join(cfg.out_dir,
                          f'model_ref_era5ref_minimal_notransp{suffix}.nc')
    log_file = os.path.join(cfg.out_dir, f'build_ref{suffix}.log')
    prov_file = os.path.join(cfg.out_dir, f'provenance{suffix}.json')

    if os.path.exists(ref_nc) and not cfg.force:
        log(f'[build_sarf_ref] {ref_nc} exists. Use --force to rebuild.')
        return {'ref_nc': ref_nc, 'log_file': log_file,
                'prov_file': prov_file, 'n_rrtmg': cfg.n_rrtmg,
                'skipped': True}

    # --- pull the data grid from a transport-run npz ---------------------
    sample = np.load(cfg.grid_npz, allow_pickle=True)
    data_lat = np.asarray(sample['lat'])
    data_lev = np.asarray(sample['lev'])
    nlat, nlev = len(data_lat), len(data_lev)
    log('=' * 70)
    log('SARF reference builder (era5ref, minimal, transport-free)')
    log('=' * 70)
    log(f'  grid          : {nlat} lat ({data_lat[0]:+.0f}..{data_lat[-1]:+.0f}'
        f') x {nlev} lev ({data_lev[0]:.3f}..{data_lev[-1]:.1f} hPa)')
    log(f'  CO2           : {cfg.co2_ppm:.0f} ppm')
    log(f'  season        : {cfg.season}')
    log(f'  n_rrtmg       : {cfg.n_rrtmg}  (cloud-overlap ensemble size)')
    log(f'  spin-up       : {cfg.n_cycle} x {cfg.t_cycle:.0f} d '
        f'(~{cfg.n_cycle * cfg.t_cycle / 365.0:.1f} yr) + '
        f'{cfg.t_avg:.0f} d averaging')
    log(f'  p_trop        : {cfg.p_trop_hPa:.0f} hPa')
    log()

    # --- build the model -------------------------------------------------
    project_cfg = load_project_config()
    model = get_ref(
        project_cfg, season_str=cfg.season,
        lat=data_lat, lev=data_lev, nlev=nlev,
        CO2=cfg.co2_ppm * 1e-6,
        n_rrtmg_repeat=cfg.n_rrtmg,
        minimal_for_sarf=True, disable_transport=True,
    )

    # --- pin troposphere to ERA5 (era5ref protocol) ----------------------
    log(f'  injecting ERA5 {cfg.season}-mean T/q/Ts ...')
    months = SeasonTypes.months_dict[cfg.season]
    era5_init = era5_annual_initial_state(model.state['Tatm'].domain,
                                          months, project_cfg)
    model.state['Tatm'][:] = era5_init['Tatm']
    model.state['Ts'][:] = era5_init['Ts']
    model.state['q'][:] = era5_init['q']
    iteratively_update_internal(model, model.state)
    model.compute_diagnostics()

    log(f'  pinning Ts, q, Tatm below p_trop={cfg.p_trop_hPa:.0f} hPa ...')
    fix_Ts(model)
    fix_q(model)
    fix_Tatm_trop(model, p_trop_hPa=cfg.p_trop_hPa)

    # --- integrate to convergence ----------------------------------------
    log()
    log('  integrating to a radiative fixed point ...')
    lat = np.asarray(model.Tatm.domain.lat.points)
    lines = []
    t0 = time.time()
    for k in range(cfg.n_cycle + 1):
        t_cyc = cfg.t_cycle if k < cfg.n_cycle else cfg.t_avg
        model.integrate_days(t_cyc + 1e-9)
        model.compute_diagnostics()
        t = model.timeave

        def ta(key):
            v = t.get(key, 0.0)
            return (np.asarray(v).squeeze()
                    if hasattr(v, 'shape') and v.shape != () else 0.0)

        asr, olr, ohu = ta('ASR'), ta('OLR'), ta('ohu')
        eb = (_lat_weighted(asr - olr - ohu, lat)
              if hasattr(asr, 'shape') else 0.0)
        msg = (f'k={k:4d}  ASR={_lat_weighted(asr, lat):7.3f}  '
               f'OLR={_lat_weighted(olr, lat):7.3f}  EB={eb:+7.4f}')
        lines.append(msg)
        if verbose and (k % 20 == 0 or k == cfg.n_cycle):
            log(f'  {msg}  ({time.time() - t0:.0f}s)')

    # --- save ------------------------------------------------------------
    out_xr = model.to_xarray(diagnostics=True, timeave=True)
    # Persist the build configuration so the SARF case driver can verify
    # its perturbed model is built consistently with this reference.
    out_xr.attrs['season_str'] = cfg.season
    out_xr.attrs['p_trop_hPa'] = cfg.p_trop_hPa
    out_xr.attrs['n_rrtmg_repeat'] = int(cfg.n_rrtmg)
    out_xr.attrs['co2_ppm'] = cfg.co2_ppm
    out_xr.to_netcdf(ref_nc)
    with open(log_file, 'w') as f:
        f.write('\n'.join(lines) + '\n')

    # --- provenance ------------------------------------------------------
    provenance = {
        'timestamp': time.strftime('%Y-%m-%dT%H:%M:%S%z'),
        'script': os.path.abspath(__file__),
        'config': asdict(cfg),
        'grid': {'nlat': int(nlat), 'nlev': int(nlev),
                 'lat_range': [float(data_lat[0]), float(data_lat[-1])],
                 'lev_range': [float(data_lev[0]), float(data_lev[-1])]},
        'repos': {'climate_runs_ext': _git_info(
            os.path.dirname(os.path.dirname(__import__(
                'climate_runs_ext').__file__)))},
    }
    with open(prov_file, 'w') as f:
        json.dump(provenance, f, indent=2)

    log()
    log(f'  saved reference  -> {ref_nc}')
    log(f'  saved log        -> {log_file}')
    log(f'  saved provenance -> {prov_file}')
    log('=' * 70)
    return {'ref_nc': ref_nc, 'log_file': log_file, 'prov_file': prov_file,
            'n_rrtmg': cfg.n_rrtmg, 'skipped': False}


def main():
    cfg = parse_config()
    build_sarf_ref(cfg, verbose=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

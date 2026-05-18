#!/usr/bin/env python
r"""Parallel SARF sweep over a preprocessed-profile manifest.

Enumerates every case in a ``manifest.csv`` written by
``preprocess_profiles``, spawns a pool of worker processes, and computes
one SARF per case via ``run_sarf_case.run_one``. Cases whose output NetCDF
already exists are skipped unless ``--force`` is given. A per-case results
CSV, a sweep log, and a provenance JSON are written to the output tree.

Two modes
---------
``--mode annual`` sweeps the annual-mean manifest against a single
annual SARF reference. ``--mode seasonal`` sweeps the four-season
manifest, selecting the matching seasonal reference
(``model_ref_era5ref_minimal_notransp_<SEASON>.nc``) per case.

cold-n20: ``--n-rrtmg`` (default 20, the published-results value) is
plumbed through to every case. ``run_sarf_case`` asserts each reference
NetCDF was built at the same ensemble size.

Two ways to use it
-------------------
Command line::

    conda run -n climlab_stardust_ext_env python scripts/run_sarf_sweep.py \
        --mode annual --material silica --nproc 30

As a Python module::

    from run_sarf_sweep import SweepConfig, run_sweep
    cfg = SweepConfig(mode='seasonal', material='silica')
    results = run_sweep(cfg)
"""

import argparse
import csv
import json
import multiprocessing
import os
import subprocess
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, asdict


# Repo layout: this file lives in transport-paper-umbrella/scripts/.
_SCRIPTS_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.dirname(_SCRIPTS_DIR)


# ===========================================================================
# Configuration
# ===========================================================================

@dataclass
class SweepConfig:
    """Configuration for one SARF sweep.

    ``mode`` selects the annual-mean or four-season manifest. The
    ``material`` / ``diameter`` / ``coag`` / ``schedule`` / ``season``
    fields filter which manifest rows are run. The remaining fields are
    the per-case integration schedule and physics, forwarded verbatim to
    ``run_sarf_case``.
    """
    mode: str = 'annual'              # 'annual' | 'seasonal'
    # input / output roots
    pre_root: str = None              # default depends on mode
    ref_dir: str = os.path.join(_REPO_ROOT, 'data', 'sarf_ref')
    out_root: str = None              # default depends on mode
    # case filters
    material: str = 'silica'          # material name, or 'all' (annual only)
    diameter: int = None              # monomer diameter [nm]; None = all
    coag: str = 'both'                # 'coag' | 'nocoag' | 'both'
    schedule: str = 'all'             # seasonal only: all|symmetric|single|...
    season: str = 'all'               # seasonal only: DJF|MAM|JJA|SON|all
    # per-case integration schedule
    n_cycle: int = 24
    t_cycle_days: float = 30.0
    t_avg_days: float = 365.0
    # per-case physics
    p_trop_hPa: float = 175.0
    co2_ppm: float = 420.0
    n_rrtmg: int = 20                 # RRTMG cloud-overlap ensemble size
    # parallelism / behaviour
    nproc: int = max(1, multiprocessing.cpu_count() - 2)
    force: bool = False
    limit: int = None                 # run only the first N cases

    def resolve_pre_root(self):
        if self.pre_root is not None:
            return self.pre_root
        sub = 'preprocessed' if self.mode == 'annual' else 'preprocessed_seasonal'
        return os.path.join(_REPO_ROOT, 'data', sub)

    def resolve_out_root(self):
        if self.out_root is not None:
            return self.out_root
        sub = 'sarf_annual' if self.mode == 'annual' else 'sarf_seasonal'
        return os.path.join(_REPO_ROOT, 'output', sub)


def parse_config():
    """Parse command-line arguments into a :class:`SweepConfig`."""
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    d = SweepConfig()
    p.add_argument('--mode', choices=['annual', 'seasonal'], default=d.mode,
                   help="'annual' sweeps the annual-mean manifest; "
                        "'seasonal' sweeps the four-season manifest")
    p.add_argument('--pre-root', default=None,
                   help='preprocessed-profile root (default: '
                        'data/preprocessed[_seasonal])')
    p.add_argument('--ref-dir', default=d.ref_dir,
                   help='directory holding the SARF reference NetCDFs')
    p.add_argument('--out-root', default=None,
                   help='output root (default: output/sarf_[annual|seasonal])')
    p.add_argument('--material', default=d.material,
                   help='material name, or "all" (annual mode only)')
    p.add_argument('--diameter', type=int, default=d.diameter,
                   help='filter by monomer diameter [nm]')
    p.add_argument('--coag', choices=['coag', 'nocoag', 'both'],
                   default=d.coag)
    p.add_argument('--schedule', default=d.schedule,
                   help='seasonal mode: all|symmetric|single|comma-list')
    p.add_argument('--season', choices=['DJF', 'MAM', 'JJA', 'SON', 'all'],
                   default=d.season, help='seasonal mode: which season(s)')
    p.add_argument('--n-cycle', type=int, default=d.n_cycle)
    p.add_argument('--t-cycle-days', type=float, default=d.t_cycle_days)
    p.add_argument('--t-avg-days', type=float, default=d.t_avg_days)
    p.add_argument('--p-trop-hPa', type=float, default=d.p_trop_hPa)
    p.add_argument('--co2', dest='co2_ppm', type=float, default=d.co2_ppm)
    p.add_argument('--n-rrtmg', type=int, default=d.n_rrtmg,
                   help='RRTMG cloud-overlap ensemble size (published SARF '
                        'uses 20; must match the reference build)')
    p.add_argument('--nproc', type=int, default=d.nproc)
    p.add_argument('--force', action='store_true',
                   help='recompute cases whose output already exists')
    p.add_argument('--limit', type=int, default=d.limit,
                   help='run only the first N cases (for testing)')
    a = p.parse_args()
    return SweepConfig(
        mode=a.mode, pre_root=a.pre_root, ref_dir=a.ref_dir,
        out_root=a.out_root, material=a.material, diameter=a.diameter,
        coag=a.coag, schedule=a.schedule, season=a.season,
        n_cycle=a.n_cycle, t_cycle_days=a.t_cycle_days,
        t_avg_days=a.t_avg_days, p_trop_hPa=a.p_trop_hPa, co2_ppm=a.co2_ppm,
        n_rrtmg=a.n_rrtmg, nproc=a.nproc, force=a.force, limit=a.limit,
    )


# ===========================================================================
# Helpers
# ===========================================================================

def _case_filename(inj_lat, inj_plev_hPa, season=None):
    """Per-case output filename, with an optional season suffix."""
    sign = 'p' if int(inj_lat) >= 0 else 'm'
    tag = f'_{season}' if season else ''
    return f'{sign}{abs(int(inj_lat))}deg_{int(inj_plev_hPa)}hpa{tag}.nc'


def _case_subdir(material, diameter_nm, coag, schedule=None):
    """Per-case output subdirectory name."""
    coag_lbl = 'coag' if (str(coag).lower() in ('true', '1')) else 'nocoag'
    base = f'{material}_d{int(diameter_nm):03d}_{coag_lbl}'
    return f'{base}_{schedule}' if schedule else base


def _git_info(path):
    """Return {sha, branch} for a git repo, or 'unknown' on failure."""
    try:
        sha = subprocess.check_output(
            ['git', '-C', path, 'rev-parse', 'HEAD'],
            stderr=subprocess.DEVNULL).decode().strip()
        branch = subprocess.check_output(
            ['git', '-C', path, 'rev-parse', '--abbrev-ref', 'HEAD'],
            stderr=subprocess.DEVNULL).decode().strip()
        return {'sha': sha, 'branch': branch}
    except Exception:
        return {'sha': 'unknown', 'branch': 'unknown'}


def _worker(task):
    """Compute one SARF case in a fresh worker process.

    Each worker pins every BLAS / OpenMP / RRTMG-Fortran library to a
    single thread: without this, N workers each spawn one OpenMP thread
    per physical core, oversubscribing the box into thrashing.
    """
    for var in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
                'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS',
                'VECLIB_MAXIMUM_THREADS'):
        os.environ[var] = '1'
    # Imported here, inside the worker, so the single-thread pinning above
    # is in effect before any numerical library initialises.
    sys.path.insert(0, _SCRIPTS_DIR)
    from run_sarf_case import SarfCaseConfig, run_one
    try:
        case_cfg = SarfCaseConfig(
            preprocessed_nc=task['preprocessed_nc'],
            ref_nc=task['ref_nc'],
            out_nc=task['out_nc'],
            radius_mapping=task['radius_mapping'],
            n_cycle=task['n_cycle'],
            t_cycle_days=task['t_cycle_days'],
            t_avg_days=task['t_avg_days'],
            p_trop_hPa=task['p_trop_hPa'],
            co2_ppm=task['co2_ppm'],
            n_rrtmg=task['n_rrtmg'],
        )
        return {**task, **run_one(case_cfg, verbose=False)}
    except Exception as e:
        return {**task, 'status': 'failed',
                'error': f'{type(e).__name__}: {e}',
                'traceback': traceback.format_exc()}


# ===========================================================================
# Case enumeration
# ===========================================================================

def _enumerate_cases(cfg):
    """Build the per-case task list from the manifest, applying filters."""
    pre_root = cfg.resolve_pre_root()
    out_root = cfg.resolve_out_root()
    manifest_path = os.path.join(pre_root, 'manifest.csv')
    if not os.path.exists(manifest_path):
        raise FileNotFoundError(
            f'manifest not found: {manifest_path}; '
            f'run preprocess_profiles.py --mode {cfg.mode} first')

    annual_ref = os.path.join(
        cfg.ref_dir, 'model_ref_era5ref_minimal_notransp.nc')

    if cfg.schedule == 'all':
        schedule_filter = {'sym', 'symseas', 'single'}
    elif cfg.schedule == 'symmetric':
        schedule_filter = {'sym', 'symseas'}
    elif cfg.schedule == 'single':
        schedule_filter = {'single'}
    else:
        schedule_filter = {s.strip() for s in cfg.schedule.split(',')}

    cases = []
    with open(manifest_path) as f:
        for row in csv.DictReader(f):
            if row.get('status') != 'ok':
                continue
            if cfg.material != 'all' and row['material'] != cfg.material:
                continue
            if cfg.diameter is not None and \
               int(row['diameter_nm']) != cfg.diameter:
                continue
            row_coag = (row['coag'] == 'True')
            if cfg.coag == 'coag' and not row_coag:
                continue
            if cfg.coag == 'nocoag' and row_coag:
                continue

            season = row.get('season', 'Annual')
            schedule = row.get('schedule', 'single')
            if cfg.mode == 'seasonal':
                if schedule not in schedule_filter:
                    continue
                if cfg.season != 'all' and season != cfg.season:
                    continue
                ref_nc = os.path.join(
                    cfg.ref_dir,
                    f'model_ref_era5ref_minimal_notransp_{season}.nc')
                subdir = _case_subdir(row['material'], row['diameter_nm'],
                                      row['coag'], schedule)
                out_name = _case_filename(row['inj_lat'],
                                          row['inj_plev_hPa'], season)
            else:
                ref_nc = annual_ref
                subdir = _case_subdir(row['material'], row['diameter_nm'],
                                      row['coag'])
                out_name = _case_filename(row['inj_lat'],
                                          row['inj_plev_hPa'])

            out_dir = os.path.join(out_root, subdir)
            out_nc = os.path.join(out_dir, out_name)
            if os.path.exists(out_nc) and not cfg.force:
                continue
            os.makedirs(out_dir, exist_ok=True)

            preprocessed_nc = row['output_path']
            cases.append({
                'material': row['material'],
                'diameter_nm': int(row['diameter_nm']),
                'coag': row_coag,
                'schedule': schedule,
                'season': season,
                'inj_lat': int(row['inj_lat']),
                'inj_plev_hPa': int(row['inj_plev_hPa']),
                'burden_Tg': float(row.get('burden_Tg', 0.0)),
                'preprocessed_nc': preprocessed_nc,
                'radius_mapping': os.path.join(
                    os.path.dirname(preprocessed_nc), 'radius_mapping.npz'),
                'ref_nc': ref_nc,
                'out_nc': out_nc,
                'n_cycle': cfg.n_cycle,
                't_cycle_days': cfg.t_cycle_days,
                't_avg_days': cfg.t_avg_days,
                'p_trop_hPa': cfg.p_trop_hPa,
                'co2_ppm': cfg.co2_ppm,
                'n_rrtmg': cfg.n_rrtmg,
            })

    # Verify each referenced reference NetCDF exists before launching.
    missing = sorted({c['ref_nc'] for c in cases
                      if not os.path.exists(c['ref_nc'])})
    if missing:
        raise FileNotFoundError(
            'missing SARF reference NetCDF(s): ' + ', '.join(missing) +
            '; build them with build_sarf_ref.py')

    if cfg.limit is not None:
        cases = cases[:cfg.limit]
    return cases


# ===========================================================================
# Run
# ===========================================================================

def run_sweep(cfg, verbose=True):
    """Run the full SARF sweep and return the per-case result records.

    Parameters
    ----------
    cfg : SweepConfig
        The sweep configuration.
    verbose : bool
        Print per-case progress.

    Returns
    -------
    list of dict
        One record per completed case (each carries either ``status``
        ``'computed'`` with ``forcing_W_m2`` or ``'failed'`` with an
        ``error``).
    """
    log = print if verbose else (lambda *a, **k: None)
    out_root = cfg.resolve_out_root()
    os.makedirs(out_root, exist_ok=True)
    log_path = os.path.join(out_root, 'sweep.log')
    results_csv = os.path.join(out_root, 'results.csv')
    prov_path = os.path.join(out_root, 'provenance.json')

    cases = _enumerate_cases(cfg)
    log(f'[sweep] mode={cfg.mode}: {len(cases)} cases to compute, '
        f'nproc={cfg.nproc}')
    if not cases:
        log('[sweep] nothing to do; outputs already exist '
            '(--force to recompute)')
        return []

    # --- provenance ------------------------------------------------------
    provenance = {
        'timestamp': time.strftime('%Y-%m-%dT%H:%M:%S%z'),
        'script': os.path.abspath(__file__),
        'config': asdict(cfg),
        'n_cases': len(cases),
        'repos': {'climate_runs_ext': _git_info(
            os.path.dirname(os.path.dirname(__import__(
                'climate_runs_ext').__file__)))},
    }
    with open(prov_path, 'w') as f:
        json.dump(provenance, f, indent=2)

    # --- run -------------------------------------------------------------
    completed = []
    t_start = time.time()
    with open(log_path, 'w') as log_f, \
            ProcessPoolExecutor(max_workers=cfg.nproc) as pool:
        futures = {pool.submit(_worker, c): c for c in cases}
        for n_done, fut in enumerate(as_completed(futures), start=1):
            res = fut.result()
            elapsed = time.time() - t_start
            tag = (f'[{n_done}/{len(cases)}] {res["material"]} '
                   f'd{res["diameter_nm"]} '
                   f'{"coag" if res["coag"] else "nocoag"} '
                   f'{res["schedule"]} {res["season"]} '
                   f'{res["inj_lat"]:+d}/{res["inj_plev_hPa"]}hpa')
            if res.get('status') == 'computed':
                msg = (f'{tag}  forcing={res["forcing_W_m2"]:+.4f} W/m^2 '
                       f'mass={res["total_mass_Tg"]:.2f}Tg '
                       f'(elapsed {elapsed / 60:.1f}min)')
            else:
                msg = (f'{tag}  FAILED: {res.get("error", "?")} '
                       f'(elapsed {elapsed / 60:.1f}min)')
            log(msg)
            log_f.write(msg + '\n')
            log_f.flush()
            completed.append(res)

    # --- write results CSV ----------------------------------------------
    keys = ['material', 'diameter_nm', 'coag', 'schedule', 'season',
            'inj_lat', 'inj_plev_hPa', 'burden_Tg', 'forcing_W_m2',
            'total_mass_Tg', 'avg_D_m', 'status', 'out_nc', 'error']
    with open(results_csv, 'w') as f:
        w = csv.DictWriter(f, fieldnames=keys, extrasaction='ignore')
        w.writeheader()
        for r in completed:
            w.writerow({k: r.get(k, '') for k in keys})

    n_ok = sum(1 for r in completed if r.get('status') == 'computed')
    n_fail = len(completed) - n_ok
    log(f'\n[sweep] done: {n_ok} ok, {n_fail} failed; '
        f'wall={(time.time() - t_start) / 60:.1f} min')
    log(f'[sweep] results CSV: {results_csv}')
    log(f'[sweep] log: {log_path}')
    return completed


def main():
    cfg = parse_config()
    completed = run_sweep(cfg, verbose=True)
    n_fail = sum(1 for r in completed if r.get('status') != 'computed')
    return 1 if n_fail else 0


if __name__ == '__main__':
    raise SystemExit(main())

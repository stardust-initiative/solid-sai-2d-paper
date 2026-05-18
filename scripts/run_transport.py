#!/usr/bin/env python
r"""Zonal-mean 2-D aerosol-transport runner (latitude x pressure), with coagulation.

The clean, config-driven driver for the stratospheric-aerosol transport
calculation of the transport paper. Per particle-size bin it assembles a
coupled transport model::

    TwoDimensionalAdvectionDiffusion  (residual winds + eddy diffusion + sedimentation)
    + ParticleSource                  (aerosol injection)
    + ParticleSink                    (removal below the tropopause)

The bins are coupled by a single ``Coagulation`` process and driven by
time-interpolated winds / diffusivities supplied by ``AtmosphericData``.

It is built on the ``climlab_stardust_extension`` and ``climate_runs_ext``
packages and stock climlab.

Two ways to use it
------------------
Command line::

    conda run -n climlab_stardust_ext_env python scripts/run_transport.py \
        --injection-geometry symmetric --inject-lat 30 --months 12 \
        --output run.npz

As a Python module -- construct a ``TransportConfig`` and call
``run_transport``; this is the friendly path for building a scenario sweep::

    from run_transport import TransportConfig, run_transport

    for lat in (10, 20, 30):
        for geom in ('single', 'symmetric'):
            cfg = TransportConfig(inject_lat=lat, injection_geometry=geom,
                                  rho_p=1800.0, months=12)
            results = run_transport(cfg, verbose=False)
            ...  # results['total_mass'], results['total_tracers'], ...

The defaults are a SHORT proof run (25x25 grid, 4 bins, ~2 months); a
paper-scale run is a larger grid with ``--months 60``.
"""

import argparse
import copy
import json
import os
import time as _time
from dataclasses import dataclass, asdict
from datetime import datetime

import numpy as np
import climlab
from climlab import couple
from climlab import constants as const
from scipy import interpolate

# --- component packages ----------------------------------------------------
from climlab_stardust_extension.dynamics import (
    TwoDimensionalAdvectionDiffusion,
    ParticleSink,
    ParticleSource,
    AtmosphericData,
)
from climlab_stardust_extension.microphysics import (
    Coagulation,
    sedimentation_velocity,
)
from climate_runs_ext import load_project_config
from climate_runs_ext.utils.transport_data import get_atm_data
from climate_runs_ext.utils.data_loading import load_xr_from_repo


# months covered by each named season (matches the extension's _SeasonalTime)
_SEASON_MONTHS = {
    'winter': [12, 1, 2], 'spring': [3, 4, 5],
    'summer': [6, 7, 8], 'fall': [9, 10, 11],
}


# ===========================================================================
# Configuration
# ===========================================================================

@dataclass
class TransportConfig:
    """Configuration for one 2-D aerosol-transport run.

    All of the paper's transport scenarios are expressed through these
    fields: injection geometry (single point vs symmetric N/S pair),
    injection timing (continuous vs seasonal), the injection point, the
    material (density ``rho_p``, monomer diameter ``d0``, fractal
    parameters), and the coagulation toggle.

    ``inject_tg`` is the *annual delivered mass* [Tg/yr]: under seasonal
    timing the instantaneous rate is scaled up so that the gated pulse
    still delivers this annual total.
    """
    # grid
    nlat: int = 25
    nlev: int = 25
    lev_top: float = 1.0          # top-of-domain pressure [hPa]
    lev_bot: float = 1000.0       # bottom-of-domain pressure [hPa]
    # size bins (fractal aggregates of monomer diameter d0)
    nbins: int = 4
    d0: float = 0.3e-6            # monomer diameter [m]
    rho_p: float = 2200.0         # material density [kg/m^3]
    Df: float = 1.6               # fractal dimension
    kf: float = 1.0               # fractal prefactor
    # injection
    inject_tg: float = 20.0       # annual delivered mass [Tg/yr]
    inject_lat: float = 15.0      # injection latitude [deg]
    inject_lev: float = 65.0      # injection pressure [hPa]
    injection_geometry: str = 'single'     # 'single' | 'symmetric' | 'antisymmetric'
    injection_timing: str = 'continuous'   # 'continuous' | 'seasonal' (single/symmetric)
    inject_months: str = 'winter'          # season name or month list, e.g. '6,7,8'
    # time
    months: int = 2
    timestep: float = 3600.0               # transport timestep [s]
    atm_timestep: float = 3600.0 * 12.0    # driving-field refresh cadence [s]
    t0: str = '2024-01-01'
    # atmospheric driver fields
    time_type: int = 1            # 0: annual mean, 1: monthly climatology
    snr_limit: float = 1.0        # Kzz signal-to-noise floor threshold
    data_year: int = -1           # -1: average all years, else a year
    # microphysics
    coagulation: bool = True
    # output
    output: str = None            # optional .npz path


def parse_config():
    """Parse command-line arguments into a :class:`TransportConfig`."""
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    # --- grid ---
    p.add_argument('--nlat', type=int, default=25)
    p.add_argument('--nlev', type=int, default=25)
    p.add_argument('--lev-top', type=float, default=1.0,
                   help='top-of-domain pressure [hPa]')
    p.add_argument('--lev-bot', type=float, default=1000.0,
                   help='bottom-of-domain pressure [hPa]')
    # --- size bins ---
    p.add_argument('--nbins', type=int, default=4)
    p.add_argument('--d0', type=float, default=0.3e-6,
                   help='monomer diameter [m]')
    p.add_argument('--rho-p', type=float, default=2200.0,
                   help='material density [kg/m^3]')
    p.add_argument('--Df', type=float, default=1.6,
                   help='fractal dimension')
    p.add_argument('--kf', type=float, default=1.0,
                   help='fractal prefactor')
    # --- injection ---
    p.add_argument('--inject-tg', type=float, default=20.0,
                   help='annual delivered mass [Tg/yr]')
    p.add_argument('--inject-lat', type=float, default=15.0,
                   help='injection latitude [deg]')
    p.add_argument('--inject-lev', type=float, default=65.0,
                   help='injection pressure [hPa]')
    p.add_argument('--injection-geometry',
                   choices=['single', 'symmetric', 'antisymmetric'],
                   default='single',
                   help="'single': one point source; 'symmetric': a mirrored "
                        "N/S pair splitting the rate; 'antisymmetric': a north "
                        "point active in --inject-months and its southern "
                        "mirror active in the complementary months, same rate")
    p.add_argument('--injection-timing', choices=['continuous', 'seasonal'],
                   default='continuous',
                   help="single/symmetric only -- 'seasonal' gates injection "
                        "to --inject-months, scaling the rate so the annual "
                        "total is preserved")
    p.add_argument('--inject-months', type=str, default='winter',
                   help="months of injection: a season name "
                        "(winter|spring|summer|fall) or a comma-separated "
                        "month-number list, e.g. '6,7,8'")
    # --- time ---
    p.add_argument('--months', type=int, default=2)
    p.add_argument('--timestep', type=float, default=3600.0,
                   help='transport timestep [s]')
    p.add_argument('--atm-timestep', type=float, default=3600.0 * 12.0,
                   help='driving-field refresh cadence [s]')
    p.add_argument('--t0', type=str, default='2024-01-01',
                   help='simulation start date (YYYY-MM-DD)')
    # --- atmospheric driver fields ---
    p.add_argument('--time-type', type=int, default=1,
                   help='0: annual mean, 1: monthly climatology')
    p.add_argument('--snr-limit', type=float, default=1.0)
    p.add_argument('--data-year', type=int, default=-1,
                   help='-1: average all years, else a calendar year')
    # --- physics toggles / output ---
    p.add_argument('--no-coagulation', action='store_true',
                   help='disable inter-bin coagulation')
    p.add_argument('--output', type=str, default=None,
                   help='optional .npz path to save the full result fields')
    a = p.parse_args()
    return TransportConfig(
        nlat=a.nlat, nlev=a.nlev, lev_top=a.lev_top, lev_bot=a.lev_bot,
        nbins=a.nbins, d0=a.d0, rho_p=a.rho_p, Df=a.Df, kf=a.kf,
        inject_tg=a.inject_tg, inject_lat=a.inject_lat, inject_lev=a.inject_lev,
        injection_geometry=a.injection_geometry,
        injection_timing=a.injection_timing,
        months=a.months, timestep=a.timestep, atm_timestep=a.atm_timestep,
        t0=a.t0, time_type=a.time_type, snr_limit=a.snr_limit,
        data_year=a.data_year, coagulation=not a.no_coagulation,
        output=a.output, inject_months=a.inject_months,
    )


# ===========================================================================
# Grid + state
# ===========================================================================

def build_grid_and_state(cfg):
    """Build the latitude x pressure grid and zero per-bin tracer states."""
    lev_bound = np.logspace(np.log10(cfg.lev_top * 100.0),
                            np.log10(cfg.lev_bot * 100.0), cfg.nlev + 1) / 100.0
    lat_bound = np.linspace(-90.0, 90.0, cfg.nlat + 1)
    lev = 0.5 * (lev_bound[:-1] + lev_bound[1:])
    lat = 0.5 * (lat_bound[:-1] + lat_bound[1:])

    # A climlab column-state only used to borrow a correctly-shaped (lat, lev)
    # array template; the transport tracers are independent state variables.
    template = climlab.column_state(lev=lev, lat=lat, water_depth=10.0)
    zero2d = 0.0 * template['Tatm']

    bin_names = [f"Si_{i}" for i in range(1, cfg.nbins + 1)]
    state = {name: 0.0 + zero2d for name in bin_names}

    return {
        'lev_bound': lev_bound, 'lat_bound': lat_bound,
        'lev': lev, 'lat': lat,
        'bin_names': bin_names, 'state': state,
    }


def cell_centered_temperature(project_cfg, lat, lev):
    """Interpolate the monthly-mean ERA5 zonal temperature onto (lat, lev)."""
    fields = load_xr_from_repo('Monthly_Zonal_Variables_2008_2017', project_cfg)
    fields = fields.sortby('latitude')
    temperature = fields['T'].mean(dim='month')
    interp = interpolate.RectBivariateSpline(
        fields.latitude, fields.level, temperature.T, kx=1, ky=1,
    )
    return interp(lat, lev)


# ===========================================================================
# Size-bin distribution (fractal aggregates: cores = monomer count)
# ===========================================================================

def build_size_bins(cfg):
    """Geometric ladder of monomer counts -> per-bin aggregate diameters.

    All injected mass is placed in the smallest bin (single monomer); the
    coagulation process moves mass up the ladder over the run.
    """
    cores = 2 ** np.arange(cfg.nbins)
    sizes = cfg.d0 * cores ** (1.0 / 3.0)          # volume-equivalent diameter
    fractions = np.zeros(cfg.nbins)
    fractions[0] = 1.0                              # inject into smallest bin
    return cores, sizes, fractions


# ===========================================================================
# Sedimentation velocity per bin
# ===========================================================================

def build_sedimentation(cfg, lat, lev_bound, project_cfg, cores):
    """Sedimentation pressure-velocity [Pa/s] on cell-vertical-bounds, per bin.

    ``sedimentation_velocity`` takes the monomer diameter ``d0`` and the
    aggregate's monomer count, NOT an aggregate diameter -- so each bin is
    just a different ``cores_number``.
    """
    fields = load_xr_from_repo('Monthly_Zonal_Variables_2008_2017', project_cfg)
    fields = fields.sortby('latitude')
    temperature = fields['T'].mean(dim='month')
    interp = interpolate.RectBivariateSpline(
        fields.latitude, fields.level, temperature.T, kx=1, ky=1,
    )
    T_on_bounds = interp(lat, lev_bound)
    rho_on_bounds = lev_bound[None, :] * 1e2 / T_on_bounds / const.Rd

    v_sed = []
    for n_cores in cores:
        omega = sedimentation_velocity(
            cfg.d0, n_cores, cfg.rho_p, T_on_bounds, rho_on_bounds,
            Df=cfg.Df, kf=cfg.kf,
        )
        # sedimentation_velocity returns a negative pressure-tendency;
        # W_sedimentation must be positive for downward settling (toward
        # higher pressure), so negate -- matching the legacy runner.
        v_sed.append(-omega)
    return v_sed


# ===========================================================================
# Injection sources
# ===========================================================================

def _parse_months(spec):
    """Parse a months specification into a sorted list of month numbers.

    Accepts a season name (``winter``/``spring``/``summer``/``fall``) or a
    comma-separated list of month numbers, e.g. ``'6,7,8'``.
    """
    spec = str(spec).strip().lower()
    if spec in _SEASON_MONTHS:
        return sorted(_SEASON_MONTHS[spec])
    try:
        months = sorted(int(x) for x in spec.split(','))
    except ValueError:
        raise ValueError(f"cannot parse months spec {spec!r}")
    if not months or any(m < 1 or m > 12 for m in months):
        raise ValueError(f"months spec {spec!r} must be 1-12")
    return months


def build_sources_config(cfg):
    """Assemble the ParticleSource config for the configured injection scenario.

    Three injection geometries:

    * ``single``        -- one point source at ``inject_lat``
    * ``symmetric``     -- a mirrored N/S pair at +/-|inject_lat|, each
      carrying half the rate, both on the same schedule
    * ``antisymmetric`` -- a north point active during ``inject_months`` and
      its southern mirror active during the complementary months, each at
      the full rate

    For ``single``/``symmetric`` the timing axis applies: ``continuous`` is
    constant injection; ``seasonal`` gates injection to ``inject_months`` and
    scales the rate by ``12 / (number of active months)``.

    ``inject_tg`` is the annual delivered mass [Tg/yr]; every scenario is
    constructed so the realised annual total equals it.
    """
    rate_kg_s = cfg.inject_tg * 1e9 / (365.0 * 86400.0)   # Tg/yr -> kg/s
    lat = abs(cfg.inject_lat)

    if cfg.injection_geometry == 'antisymmetric':
        north = _parse_months(cfg.inject_months)
        south = [m for m in range(1, 13) if m not in north]
        # the two month-sets tile the year, so each source at the full rate
        # delivers a combined annual total of exactly inject_tg
        return [
            {'name': 'injection north', 'space_type': 'single_grid_point',
             'time_type': 'by_month', 'month_list': north,
             'rate': rate_kg_s, 'point_source': [lat, cfg.inject_lev]},
            {'name': 'injection south', 'space_type': 'single_grid_point',
             'time_type': 'by_month', 'month_list': south,
             'rate': rate_kg_s, 'point_source': [-lat, cfg.inject_lev]},
        ]

    # --- single / symmetric ---------------------------------------------
    if cfg.injection_timing == 'seasonal':
        months = _parse_months(cfg.inject_months)
        rate_kg_s *= 12.0 / len(months)            # preserve the annual total
        timing = {'time_type': 'by_month', 'month_list': months}
    else:
        timing = {'time_type': 'const'}

    if cfg.injection_geometry == 'symmetric':
        lats = [lat, -lat]
        per_source_rate = rate_kg_s / 2.0          # split across hemispheres
    else:
        lats = [cfg.inject_lat]
        per_source_rate = rate_kg_s

    sources = []
    for k, lat_inj in enumerate(lats):
        src = {
            'name': f'injection {k + 1}',
            'space_type': 'single_grid_point',
            'rate': per_source_rate,
            'point_source': [lat_inj, cfg.inject_lev],
        }
        src.update(timing)
        sources.append(src)
    return sources


# ===========================================================================
# Per-bin transport model
# ===========================================================================

def build_bin_transport(name, bin_state, cfg, fields_2d, v_sed,
                         bin_fraction, get_time, temperature_cells,
                         sources_config):
    """Couple advection-diffusion + source + sink for a single size bin."""
    advdiff = TwoDimensionalAdvectionDiffusion(
        name=f'aerosol Transport 2D {name}', state=bin_state,
        timestep=cfg.timestep,
        U=fields_2d['vlat'], W=fields_2d['vlev'], W_sedimentation=v_sed,
        Kyy=fields_2d['kyy'], Kzz=fields_2d['kzz'], Kyz=fields_2d['kyz'],
        diagnostic_name_suffix=f'_{name}',
    )
    sink = ParticleSink(
        name=f'transport sink {name}', state=bin_state, timestep=cfg.timestep,
        tropopause_p=fields_2d['tropopause'],
    )
    # ParticleSource scales each source dict's 'rate' in place by bin_fraction,
    # so every bin must receive its own deep copy of the source list.
    source = ParticleSource(
        name=f'transport source {name}', state=bin_state, timestep=cfg.timestep,
        bin_fraction=bin_fraction, current_time=get_time,
        diagnostic_name_suffix=f'_{name}',
        sources_config=copy.deepcopy(sources_config),
        temperature=temperature_cells,
    )
    return couple([advdiff, source, sink], name=f'total transport {name}')


# ===========================================================================
# Run
# ===========================================================================

def run_transport(cfg, verbose=True):
    """Run one 2-D aerosol-transport calculation.

    Parameters
    ----------
    cfg : TransportConfig
        The run configuration.
    verbose : bool
        Print a header, per-month progress, and a diagnostics summary.

    Returns
    -------
    dict
        Result fields, with array axis 0 = size bin and axis 1 = time
        (time index 0 is the initial zero state):

        ``lat``, ``lev``        -- grid cell centers
        ``cores``, ``sizes``    -- per-bin monomer count and diameter
        ``diag_days``           -- elapsed-time axis [days]
        ``total_mass``          -- per-bin burden [kg]            (nbins, ntime)
        ``total_source``        -- per-bin cumulative injected [kg]
        ``total_column_density``-- per-bin column density   (nbins, ntime, nlat)
        ``total_tracers``       -- per-bin mixing-ratio field
                                   (nbins, ntime, nlat, nlev)
        ``coag_tendencies``     -- coagulation tendencies (None if disabled)
        ``config_json``         -- the run configuration, for provenance
        ``physical``            -- bool: output finite and mass-conserving
    """
    log = print if verbose else (lambda *a, **k: None)
    t0 = datetime.strptime(cfg.t0, '%Y-%m-%d')

    log('=' * 70)
    log('2-D aerosol-transport runner (zonal mean, latitude x pressure)')
    log('=' * 70)
    log(f'  grid          : {cfg.nlat} lat x {cfg.nlev} lev')
    log(f'  size bins     : {cfg.nbins}  (d0 = {cfg.d0 * 1e6:.3f} um, '
        f'rho_p = {cfg.rho_p:.0f} kg/m^3)')
    log(f'  coagulation   : {"ON" if cfg.coagulation else "OFF"}')
    if cfg.injection_geometry == 'antisymmetric':
        log(f'  injection     : {cfg.inject_tg:.1f} Tg/yr | antisymmetric | '
            f'north in [{cfg.inject_months}], south in the complement')
    else:
        _t = (f' ({cfg.inject_months})'
              if cfg.injection_timing == 'seasonal' else '')
        log(f'  injection     : {cfg.inject_tg:.1f} Tg/yr | '
            f'{cfg.injection_geometry} geometry | '
            f'{cfg.injection_timing} timing{_t}')
    log(f'  inject point  : {cfg.inject_lat:.0f} deg, {cfg.inject_lev:.0f} hPa'
        + ('  (+ N/S mirror)'
           if cfg.injection_geometry in ('symmetric', 'antisymmetric')
           else ''))
    log(f'  integration   : {cfg.months} month(s), dt = {cfg.timestep:.0f} s')
    log()

    # --- project config + grid -------------------------------------------
    project_cfg = load_project_config()
    grid = build_grid_and_state(cfg)
    lat, lev = grid['lat'], grid['lev']
    state = grid['state']
    bin_names = grid['bin_names']

    temperature_cells = cell_centered_temperature(project_cfg, lat, lev)
    log(f'  temperature   : {temperature_cells.min():.1f} - '
        f'{temperature_cells.max():.1f} K (cell centers)')

    # --- atmospheric driver fields ---------------------------------------
    log('  loading transport driver fields (offline cache) ...')
    vlev, vlat, kzz, kyy, kyz, tropopause = get_atm_data(
        project_cfg, time_type=cfg.time_type,
        snr_limit=cfg.snr_limit, data_year=cfg.data_year,
    )
    atm_data = AtmosphericData(
        param_configs={
            'kyy': {'data': kyy, 'method': 'linear', 'grid_type': 'bounds*centers'},
            'kzz': {'data': kzz, 'method': 'linear', 'grid_type': 'centers*bounds'},
            'kyz': {'data': kyz, 'method': 'linear', 'grid_type': 'bounds*bounds'},
            'vlat': {'data': vlat, 'method': 'linear', 'grid_type': 'bounds*centers'},
            'vlev': {'data': vlev, 'method': 'linear', 'grid_type': 'centers*bounds'},
            'tropopause': {'data': tropopause, 'method': 'linear', 'grid_type': 'centers'},
        },
        t_0=t0, time_type=cfg.time_type, state=state,
        timestep=cfg.atm_timestep, name='atm data for transport',
    )
    fields_2d = {k: atm_data._current_data[k].data
                 for k in ('kyy', 'kzz', 'kyz', 'vlat', 'vlev', 'tropopause')}
    get_time = atm_data.get_current_time
    log(f'  tropopause    : {fields_2d["tropopause"].min():.0f} - '
        f'{fields_2d["tropopause"].max():.0f} hPa')

    # --- size bins + sedimentation ---------------------------------------
    cores, sizes, fractions = build_size_bins(cfg)
    log(f'  aggregate dia : {", ".join(f"{s*1e6:.3f}" for s in sizes)} um')
    v_sed_all = build_sedimentation(cfg, lat, grid['lev_bound'],
                                    project_cfg, cores)

    # --- injection source(s) + per-bin transport models ------------------
    sources_config = build_sources_config(cfg)
    bin_models = [
        build_bin_transport(name, {name: state[name]}, cfg, fields_2d,
                            v_sed_all[i], fractions[i], get_time,
                            temperature_cells, sources_config)
        for i, name in enumerate(bin_names)
    ]

    # --- coagulation couples the bins together ---------------------------
    processes = bin_models + [atm_data]
    coag = None
    if cfg.coagulation:
        coag = Coagulation(
            name='coagulation process', state=state, d0=cfg.d0,
            cores={name: int(cores[i]) for i, name in enumerate(bin_names)},
            timestep=cfg.timestep, temperature=temperature_cells,
            rho_p=cfg.rho_p, Df=cfg.Df, kf=cfg.kf,
        )
        processes.append(coag)
    transport_model = couple(processes, name='all bins transport model')

    # --- integration loop (records the full diagnostic field set) --------
    log()
    log('  integrating ...')
    nbins, nj, ni = cfg.nbins, cfg.nlat, cfg.nlev
    total_mass = [[np.array([0.0])] for _ in range(nbins)]
    total_source = [[0.0] for _ in range(nbins)]
    total_column_density = [[np.zeros(nj)] for _ in range(nbins)]
    total_tracers = [[np.zeros([nj, ni])] for _ in range(nbins)]
    diag_days = [0.0]

    wall0 = _time.time()
    for m in range(cfg.months):
        transport_model.integrate_days(365.0 / 12.0 + 1e-6)
        for i, name in enumerate(bin_names):
            total_mass[i].append(
                transport_model.timeave[f'total_tracer_mass_{name}'])
            total_source[i].append(
                transport_model.diagnostics[f'total_tracer_source_{name}'][0])
            total_column_density[i].append(
                transport_model.timeave[f'column_density_{name}'])
            total_tracers[i].append(transport_model.timeave[name])
        diag_days.append(transport_model.time['days_elapsed'])

        burden = sum(float(np.atleast_1d(total_mass[i][-1])[0])
                     for i in range(nbins)) / 1e9
        injected = sum(total_source[i][-1] for i in range(nbins)) / 1e9
        log(f'    month {m + 1:2d}: burden = {burden:9.4f} Tg | '
            f'injected = {injected:9.4f} Tg  [{_time.time() - wall0:7.1f} s]')

    total_mass = np.array(total_mass)[:, :, 0]            # (nbins, ntime)
    total_source = np.array(total_source)                 # (nbins, ntime)
    total_column_density = np.array(total_column_density)  # (nbins, ntime, nj)
    total_tracers = np.array(total_tracers)               # (nbins, ntime, nj, ni)
    diag_days = np.array(diag_days)

    # --- summary ---------------------------------------------------------
    final_burden = total_mass[:, -1].sum() / 1e9
    final_injected = total_source[:, -1].sum() / 1e9
    finite = bool(np.all(np.isfinite(total_mass)))
    mass_closed = bool(0.0 <= final_burden <= final_injected + 1e-9)
    physical = finite and mass_closed
    log()
    log('  summary')
    log('  ' + '-' * 60)
    log(f'    total injected       : {final_injected:.4f} Tg')
    log(f'    final burden         : {final_burden:.4f} Tg')
    log(f'    removed (trop. sink) : {final_injected - final_burden:.4f} Tg')
    log(f'    per-bin burden [Tg]  : '
        + ', '.join(f'{x:.4f}' for x in total_mass[:, -1] / 1e9))
    if coag is not None:
        log(f'    mass in bins > 1     : '
            f'{total_mass[1:, -1].sum() / 1e9:.4f} Tg (coagulation)')
    log(f'    all output finite    : {finite}')
    log(f'    mass-conserving      : {mass_closed}')
    log(f'    wall time            : {_time.time() - wall0:.1f} s')
    log(f'  RESULT: {"PASS" if physical else "FAIL - see diagnostics above"}')
    log('=' * 70)

    return {
        'lat': lat, 'lev': lev, 'cores': cores, 'sizes': sizes,
        'diag_days': diag_days,
        'total_mass': total_mass, 'total_source': total_source,
        'total_column_density': total_column_density,
        'total_tracers': total_tracers,
        'coag_tendencies': coag.tendencies if coag is not None else None,
        'config_json': json.dumps(asdict(cfg)),
        'physical': physical,
    }


def main():
    cfg = parse_config()
    results = run_transport(cfg, verbose=True)
    if cfg.output:
        np.savez(cfg.output, **{k: v for k, v in results.items()
                                if v is not None})
        print(f'  results saved to {os.path.abspath(cfg.output)}')
    return 0 if results['physical'] else 1


if __name__ == '__main__':
    raise SystemExit(main())

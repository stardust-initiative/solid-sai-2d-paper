#!/usr/bin/env python
r"""Compute the SARF of one preprocessed transport-run aerosol profile.

SARF -- stratospheric-aerosol radiative forcing -- is the top-of-atmosphere
radiation imbalance induced by an aerosol layer when the surface, the
water-vapour field, and the tropospheric temperature are all held fixed.
Only the stratosphere is free to respond, so the result isolates the
direct radiative effect of the aerosol.

This driver computes one case:

  1. load the SARF reference NetCDF (built by ``build_sarf_ref``) and
     rebuild the matching ``model_ref``;
  2. load the preprocessed per-bin aerosol mass-mixing-ratio profile,
     convert it to volume mixing ratio per size bin, and build the
     aerosol optics (``AerosolsOptDepTables`` /
     ``get_radiation_with_aerosols_params``);
  3. build the perturbed model ``model_pert`` with that aerosol optics;
  4. pin both models -- ``fix_Ts``, ``fix_q``, ``fix_Tatm_trop`` -- so
     only the stratosphere responds (the SARF constraint set);
  5. integrate ``model_pert`` to a radiative fixed point;
  6. save the perturbed-minus-reference anomaly NetCDF, with the
     lat-weighted TOA forcing in ``forcing_W_m2``.

cold-n20: the published SARF results use ``n_rrtmg_repeat = 20`` (the
RRTMG cloud-overlap MCICA ensemble size; the rev6 model default is 5).
``--n-rrtmg`` is passed straight into ``get_ref`` for both the reference
and the perturbed model, and the loaded reference NetCDF's
``n_rrtmg_repeat`` attr is asserted to match -- so a SARF case can never silently mix a
reference built at one ensemble size with a perturbed model at another.

Two ways to use it
-------------------
Command line::

    conda run -n climlab_stardust_ext_env python scripts/run_sarf_case.py \
        --preprocessed-nc data/preprocessed/silica_d500_coag/p10deg_68hpa.nc \
        --ref-nc data/sarf_ref/model_ref_era5ref_minimal_notransp.nc \
        --out-nc out/p10deg_68hpa.nc

As a Python module (this is the path ``run_sarf_sweep`` uses)::

    from run_sarf_case import SarfCaseConfig, run_one
    cfg = SarfCaseConfig(preprocessed_nc=..., ref_nc=..., out_nc=...)
    result = run_one(cfg)
"""

import argparse
import json
import os
import sys
import time
import traceback
from dataclasses import dataclass

import numpy as np
import xarray as xr

from climate_runs_ext import load_project_config
from climate_runs_ext.reference_model import get_ref
from climate_runs_ext.utils.state_io import (
    update_model_from_xr, iteratively_update_internal,
)
from climate_runs_ext.utils.model_control import fix_Ts, fix_q, fix_Tatm_trop
from climate_runs_ext.utils.era5_data import lat_avg
from climlab_stardust_extension.radiation.optical_depth_tables_aerosols import (
    aerosol_instance, AerosolsOptDepTables,
    get_radiation_with_aerosols_params,
)
from climlab_stardust_extension.utils.constants import n_avogadro
from climlab.utils import constants as const


# Preprocessor material name -> key in config.json:aerosols_table_dict.
# AerosolsOptDepTables does a case-insensitive membership check, so these
# values must match the config keys modulo case.
MATERIAL_TABLE_KEY = {
    'silica':  'silica',
    'sulfate': 'sulfate',
    'calcite': 'calcite average ray',
}


# ===========================================================================
# Configuration
# ===========================================================================

@dataclass
class SarfCaseConfig:
    """Configuration for one SARF case.

    The aerosol profile (``preprocessed_nc`` + its companion
    ``radius_mapping``) and the SARF reference (``ref_nc``) are the
    inputs; ``out_nc`` is the anomaly output. The remaining fields are
    the integration schedule and the physical knobs -- which must be
    consistent with how ``ref_nc`` was built (the ``n_rrtmg_repeat``
    consistency is checked at run time).
    """
    preprocessed_nc: str = None       # per-bin aerosol mmr profile
    ref_nc: str = None                # SARF reference NetCDF
    out_nc: str = None                # anomaly output NetCDF
    radius_mapping: str = None        # default: companion of preprocessed_nc
    # integration schedule
    n_cycle: int = 24                 # spin-up cycles
    t_cycle_days: float = 30.0        # days per spin-up cycle
    t_avg_days: float = 365.0         # days in the final averaging cycle
    # physics (must match the reference build)
    p_trop_hPa: float = 175.0         # tropopause pin pressure
    co2_ppm: float = 420.0
    n_rrtmg: int = 20                 # RRTMG cloud-overlap ensemble size

    def resolve_radius_mapping(self):
        """Return the radius-mapping path, defaulting to the companion file."""
        if self.radius_mapping is not None:
            return self.radius_mapping
        return os.path.join(os.path.dirname(self.preprocessed_nc),
                            'radius_mapping.npz')


def parse_config():
    """Parse command-line arguments into a :class:`SarfCaseConfig`."""
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    d = SarfCaseConfig()
    p.add_argument('--preprocessed-nc', required=True,
                   help='preprocessed per-bin aerosol mmr NetCDF')
    p.add_argument('--ref-nc', required=True,
                   help='SARF reference NetCDF (from build_sarf_ref)')
    p.add_argument('--out-nc', required=True,
                   help='output anomaly NetCDF (perturbed minus reference)')
    p.add_argument('--radius-mapping', default=None,
                   help='bin-name -> radius npz (default: radius_mapping.npz '
                        'next to --preprocessed-nc)')
    p.add_argument('--n-cycle', type=int, default=d.n_cycle,
                   help='number of spin-up cycles')
    p.add_argument('--t-cycle-days', type=float, default=d.t_cycle_days,
                   help='days per spin-up cycle')
    p.add_argument('--t-avg-days', type=float, default=d.t_avg_days,
                   help='days in the final averaging cycle')
    p.add_argument('--p-trop-hPa', type=float, default=d.p_trop_hPa,
                   help='tropopause pin pressure [hPa]')
    p.add_argument('--co2', dest='co2_ppm', type=float, default=d.co2_ppm,
                   help='CO2 concentration [ppm]')
    p.add_argument('--n-rrtmg', type=int, default=d.n_rrtmg,
                   help='RRTMG cloud-overlap ensemble size (published SARF '
                        'uses 20; must match the reference build)')
    a = p.parse_args()
    return SarfCaseConfig(
        preprocessed_nc=a.preprocessed_nc, ref_nc=a.ref_nc, out_nc=a.out_nc,
        radius_mapping=a.radius_mapping, n_cycle=a.n_cycle,
        t_cycle_days=a.t_cycle_days, t_avg_days=a.t_avg_days,
        p_trop_hPa=a.p_trop_hPa, co2_ppm=a.co2_ppm, n_rrtmg=a.n_rrtmg,
    )


# ===========================================================================
# Helpers
# ===========================================================================

def _ref_season(ref_xr, ref_nc):
    """Determine the insolation season the reference was built with.

    Looked up from the ``season_str`` attr (current builders set this),
    falling back to a filename-suffix parse, then to ``Annual``. Building
    the perturbed model with a season different from the reference's
    causes a polar-night stratopause artifact, so this must be exact.
    """
    season = ref_xr.attrs.get('season_str', None)
    if season is not None:
        return season
    base = os.path.basename(ref_nc).rstrip('.nc')
    for s in ('DJF', 'MAM', 'JJA', 'SON'):
        if base.endswith(f'_{s}'):
            return s
    return 'Annual'


def _column_burden(mmr, data_lat, data_lev):
    """Lat-weighted total column burden [Tg] of a ``(lat, lev)`` mmr field."""
    mid = 0.5 * (data_lev[:-1] + data_lev[1:])
    lev_b = np.concatenate([[0.0], mid, [1013.25]])
    dp_pa = 1e2 * np.diff(lev_b)
    col = np.sum(mmr * dp_pa[None, :] / 9.80665, axis=-1)   # kg/m^2
    R = 6.371e6
    dlat = np.deg2rad(2.0)
    area = 2 * np.pi * R ** 2 * dlat * np.cos(np.deg2rad(data_lat))
    return float((col * area).sum() / 1e9)


# ===========================================================================
# Run
# ===========================================================================

def run_one(cfg, verbose=True):
    """Compute the SARF for one preprocessed case.

    Parameters
    ----------
    cfg : SarfCaseConfig
        The case configuration.
    verbose : bool
        Print progress and per-cycle energy-balance diagnostics.

    Returns
    -------
    dict
        ``forcing_W_m2``  -- lat-weighted TOA SARF
        ``total_mass_Tg`` -- total aerosol burden
        ``avg_D_m``       -- mass-weighted mean aerosol diameter
        ``out_nc``        -- path of the written anomaly NetCDF
        ``status``        -- ``'computed'``
    """
    project_cfg = load_project_config()
    radius_mapping_npz = cfg.resolve_radius_mapping()

    # --- load reference --------------------------------------------------
    ref_xr = xr.open_dataset(cfg.ref_nc).load()
    data_lat = ref_xr.lat.values
    data_lev = ref_xr.lev.values
    nlev = len(data_lev)

    # cold-n20: the reference must have been built at the same RRTMG
    # ensemble size we are about to build the perturbed model with.
    # ``n_rrtmg_repeat`` is the get_ref kwarg name and the attr written
    # by build_sarf_ref.py and carried by every legacy reference NetCDF.
    ref_n_rrtmg = ref_xr.attrs.get('n_rrtmg_repeat', None)
    if ref_n_rrtmg is not None:
        assert int(ref_n_rrtmg) == int(cfg.n_rrtmg), (
            f'n_rrtmg mismatch: reference NetCDF was built with '
            f'n_rrtmg={int(ref_n_rrtmg)} but this case requests '
            f'n_rrtmg={int(cfg.n_rrtmg)}')
    elif verbose:
        print('[run_sarf_case] WARNING: reference NetCDF has no '
              'n_rrtmg_repeat attr; cannot verify cloud-ensemble '
              'consistency', flush=True)

    season = _ref_season(ref_xr, cfg.ref_nc)
    if verbose:
        print(f'[run_sarf_case] season={season}, n_rrtmg={cfg.n_rrtmg}, '
              f'ref={os.path.basename(cfg.ref_nc)}', flush=True)

    kwargs_ref = dict(
        season_str=season, lat=data_lat, lev=data_lev, nlev=nlev,
        CO2=cfg.co2_ppm * 1e-6, n_rrtmg_repeat=cfg.n_rrtmg,
        minimal_for_sarf=True, disable_transport=True,
    )
    model_ref = get_ref(project_cfg, **kwargs_ref)
    update_model_from_xr(model_ref, ref_xr, do_compute=True)
    rad = model_ref.subprocess['Atmosphere'].subprocess['Radiation']
    coszen = rad.coszen
    domain = rad.Tatm.domain
    ref_state = {k: model_ref.timeave[k] for k in model_ref.state.keys()}

    # --- load the preprocessed aerosol profile ---------------------------
    state_xr = xr.open_dataset(cfg.preprocessed_nc)
    rho_particle = float(state_xr.attrs['particle_density_kg_m3'])
    material_raw = state_xr.attrs['material']
    material = MATERIAL_TABLE_KEY.get(material_raw, material_raw)
    rm_data = np.load(radius_mapping_npz)
    radius_mapping = {k: float(rm_data[k]) for k in rm_data.files}

    m_air = 1e-3 * const.molecular_weight['dry air'] / n_avogadro

    aerosol_instance_list = []
    mass_tot_Tg = 0.0
    avg_D = 0.0
    for name, r_m in radius_mapping.items():
        m_p = (4.0 / 3.0) * np.pi * r_m ** 3 * rho_particle
        mmr = state_xr[name].values                  # (lat, lev)
        if mmr.shape != (len(data_lat), len(data_lev)):
            mmr = mmr.T
        assert mmr.shape == (len(data_lat), len(data_lev)), (
            f'{name}: shape {mmr.shape} vs expected '
            f'({len(data_lat)}, {len(data_lev)})')
        mmr = np.where(mmr > 0.0, mmr, 0.0)
        vmr = m_air / m_p * mmr
        aerosol_instance_list.append(aerosol_instance(material, r_m, vmr))
        dm = _column_burden(mmr, data_lat, data_lev)
        mass_tot_Tg += dm
        avg_D += 2 * r_m * dm
    if mass_tot_Tg > 0:
        avg_D /= mass_tot_Tg
    inj_lat = int(state_xr.attrs.get('inj_lat_deg', 0))
    inj_plev = int(state_xr.attrs.get('inj_plev_hPa', 0))
    state_xr.close()

    if verbose:
        print(f'[run_sarf_case] material={material} rho={rho_particle:g} '
              f'n_bins={len(aerosol_instance_list)} '
              f'mass={mass_tot_Tg:.2f}Tg avg_D={avg_D * 1e9:.0f}nm', flush=True)

    # --- build the perturbed model --------------------------------------
    aer = AerosolsOptDepTables(
        aerosol_instance_list=aerosol_instance_list,
        domain=domain, coszen=coszen,
        **project_cfg['aerosols_input_dict'],
    )
    rad_with_aero = get_radiation_with_aerosols_params(
        rad.state, aer, rad.coszen)

    model_pert = get_ref(project_cfg, **kwargs_ref,
                         rad_with_aero_param_dict=rad_with_aero)
    update_model_from_xr(model_pert, ref_xr, do_compute=True)
    iteratively_update_internal(model_pert, ref_state)

    # --- assert the bypass-RRTMG SW path is engaged ---------------------
    sw_proc = (model_pert.subprocess['Atmosphere'].subprocess['Radiation']
               .subprocess['SW'])
    assert int(np.asarray(sw_proc.add_aero_layer).ravel()[0]) == 1, (
        f'bypass-RRTMG SW not engaged: add_aero_layer={sw_proc.add_aero_layer}')
    r_mu_arr = np.asarray(sw_proc.r_mu)
    t_mu_arr = np.asarray(sw_proc.t_mu)
    assert r_mu_arr.ndim == 3, (
        f'bypass r_mu wrong shape: {r_mu_arr.shape} (expected (spec, lat, lev))')
    assert r_mu_arr.max() > 0.0, (
        f'bypass r_mu trivially zero everywhere: max={r_mu_arr.max():.3e}')
    assert t_mu_arr.min() < 1.0 - 1e-6, (
        f'bypass t_mu trivially unity everywhere: min={t_mu_arr.min():.6f}')
    if verbose:
        print(f'[run_sarf_case] bypass OK: add_aero_layer=1, '
              f'r_mu max={r_mu_arr.max():.3e}, '
              f't_mu min={t_mu_arr.min():.6f}', flush=True)

    # --- apply the SARF constraint set ----------------------------------
    # Surface, water vapour, and tropospheric temperature are pinned, so
    # only the stratosphere responds to the aerosol.
    fix_Ts(model_pert)
    fix_q(model_pert)
    fix_Tatm_trop(model_pert, p_trop_hPa=cfg.p_trop_hPa)

    # --- integrate to a radiative fixed point ---------------------------
    t_start = time.time()
    for k in range(cfg.n_cycle + 1):
        t_cyc = cfg.t_cycle_days if k < cfg.n_cycle else cfg.t_avg_days
        model_pert.integrate_days(t_cyc + 1e-9)
        model_pert.compute_diagnostics()
        if verbose and (k % 6 == 0 or k == cfg.n_cycle):
            asr = np.asarray(model_pert.timeave['ASR'])
            olr = np.asarray(model_pert.timeave['OLR'])
            ohu = np.asarray(model_pert.timeave.get('ohu', 0.0))
            eb = lat_avg((asr - olr - ohu)[:, 0], data_lat)
            print(f'[run_sarf_case] k={k:3d} EB={eb:+7.4f} '
                  f'(elapsed {time.time() - t_start:.1f}s)', flush=True)

    # --- save the perturbed-minus-reference anomaly ---------------------
    pert_xr = model_pert.to_xarray(diagnostics=True, timeave=True)
    diag = pert_xr - ref_xr

    dASR = np.asarray(diag.ASR.values).squeeze()
    dOLR = np.asarray(diag.OLR.values).squeeze()
    forcing = float(lat_avg(dASR - dOLR, data_lat))

    diag.attrs.update({
        'forcing_W_m2': forcing,
        'total_mass_Tg': mass_tot_Tg,
        'avg_D_m': avg_D,
        'mode': 'sarf',
        'n_cycle': cfg.n_cycle,
        't_cycle_days': cfg.t_cycle_days,
        't_avg_days': cfg.t_avg_days,
        'n_rrtmg': int(cfg.n_rrtmg),
        'co2_ppm': cfg.co2_ppm,
        'p_trop_hPa': cfg.p_trop_hPa,
        'preprocessed_nc': cfg.preprocessed_nc,
        'radius_mapping_npz': radius_mapping_npz,
        'ref_nc': cfg.ref_nc,
        'inj_lat_deg': inj_lat,
        'inj_plev_hPa': inj_plev,
    })

    out_dir = os.path.dirname(cfg.out_nc)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    diag.to_netcdf(cfg.out_nc)
    if verbose:
        print(f'[run_sarf_case] saved {cfg.out_nc}  '
              f'forcing={forcing:.4f} W/m^2', flush=True)
    return {
        'forcing_W_m2': forcing,
        'total_mass_Tg': mass_tot_Tg,
        'avg_D_m': avg_D,
        'out_nc': cfg.out_nc,
        'status': 'computed',
    }


def main():
    cfg = parse_config()
    try:
        result = run_one(cfg, verbose=True)
        print(json.dumps(result, indent=2))
    except Exception as e:
        print(f'FAILED: {type(e).__name__}: {e}')
        traceback.print_exc()
        return 1
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

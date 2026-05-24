# Material-comparison figure — panels (a) & (b)

Self-contained bundle for the 2×2 material-comparison figure
(silica / calcite / sulfate, monomer diameter **d = 500 nm**). The **top row**
(panels **a** and **b**) is produced here; the **bottom row** (panels **c**, **d**)
is left as empty placeholder axes for the co-author to fill in.

## Files

| file | what it is | needs |
|---|---|---|
| `mie_d500_data.npz` | **cached** Mie spectral slices at d = 500 nm for the three materials (β̄, q_sca, q_ext, q_abs on a solar grid and a thermal grid), plus the spectral-weighting parameters | — (just data) |
| `make_panels_ab.py` | builds the 2×2 figure; populates (a)/(b); leaves (c)/(d) as placeholders | `numpy`, `matplotlib` |
| `compute_mie_d500.py` | (re)computes the Mie curves from the refractive-index tables — exposes `compute_all(ri_dir)` and a CLI to refresh the `.npz` | `numpy`, `pandas`, `scipy`, `miepython` + the `Materials_optical_tables` checkout |
| `figure_2x2_material_panels.{pdf,png}` | current output, 300 dpi — panel (b) **log y** | — |
| `figure_2x2_material_panels_linB.{pdf,png}` | same figure but panel (b) has a **linear y** axis | — |

Both versions are written on every run of `make_panels_ab.py` (panels (a),
(c), (d) are identical between them).

## Run

```bash
python make_panels_ab.py            # -> figure_2x2_material_panels.{pdf,png}
```

By default this just reads the bundled `mie_d500_data.npz` — only `numpy` and
`matplotlib` are needed. To recompute the curves from the refractive-index
tables instead:

```bash
python make_panels_ab.py --recompute --ri-dir /path/to/Materials_optical_tables
# or, to refresh the cache file itself:
python compute_mie_d500.py --ri-dir /path/to/Materials_optical_tables
```

If `--recompute` is requested (or the cache file is absent) but the calculation
packages are missing, the script **falls back to `mie_d500_data.npz` and prints
a warning** rather than erroring — so it always produces the figure as long as
the cache is present.

## What panels (a) and (b) show

- **(a)** volume-normalised effective back-scatter
  `σ_b/V = (3/4)·β̄·q_sca / a`  (a = particle radius = 250 nm), units **µm⁻¹**,
  versus wavelength over the solar range (0–3 µm), **linear y**.
- **(b)** volume-normalised effective absorption
  `σ_abs/V = (3/4)·q_abs / a`, units **µm⁻¹**, versus wavelength over the
  thermal range (1–50 µm), **log y**.

`β̄` is the exact expression used throughout the transport-paper Mie comparison:
`β̄ = (1/2π) ∫_{-1}^{+1} P(µ) arccos(µ) dµ = (1/2) ∫_0^π P(cosΘ) sinΘ (Θ/π) dΘ`,
with `P` the full Mie unpolarised phase function (4π-normalised) — no
asymmetry-parameter / two-stream approximation. (`σ/V` is density-independent;
the bulk densities are recorded in the `.npz` for reference only.)

Each panel's legend entry is `<material>  (<spectral-weighted average> µm⁻¹)`,
where the average is over the panel's weighting function: panel (a) — solar
Planck, T = 6000 K, zeroed below λ = 0.25 µm; panel (b) — piecewise Planck,
brightness temperatures `T = [255, 220, 290, 255] K` for the four thermal
windows split at `[600, 800, 1250] cm⁻¹`. Those weighting functions are also
drawn on each panel as a grey semi-transparent fill + black line **behind the
material curves, peaking at ~90 % of the panel height** (the overlay is a
*shape*, not a value — its absolute height is arbitrary; on the log-y panel (b)
it sits on a hidden linear twin axis).

Colours: silica `#0072B2`, calcite `#E69F00`, sulfate `#D55E00` — three of the
Okabe–Ito colour-blind-safe palette (blue / orange / vermillion).

Current spectral-weighted averages (µm⁻¹):

| | panel (a) ⟨σ_b/V⟩ (solar) | panel (b) ⟨σ_abs/V⟩ (thermal) |
|---|---|---|
| silica  | 1.30 | 0.384 |
| calcite | 1.67 | 0.113 |
| sulfate | 1.26 | 0.201 |

## Completing (c) and (d)

In `make_panels_ab.py`, `axes[1, 0]` (= `axC`) and `axes[1, 1]` (= `axD`)
are the placeholders. Replace the "reserved for co-author" block with your
plotting code; the figure size, fonts, dpi, panel labels, and `constrained_layout`
are already set so the rows will match.

# Public API (Grok Build)

NASA CEA rocket + collisionless plume, callable from a Grok Build web app on `https://*.grok.me` (and local preview). This is the contract Grok Build should call. The desktop Vite UI is unchanged.

Base URL: the Render / Railway host, e.g. `https://YOUR-SERVICE.onrender.com`. All routes are under `/api`.

## Endpoints

| Method | Path | Auth | Body |
|---|---|---|---|
| `GET` | `/api/health` | never (load-balancer friendly) | — |
| `GET` | `/api/catalog` | `X-API-Key` if `API_KEY` set | — |
| `POST` | `/api/mixture-preview` | same | JSON mixture |
| `POST` | `/api/solve` | same | JSON (below) |
| `POST` | `/api/characteristics` | same | JSON (below) |
| `OPTIONS` | `/api/*` | never (CORS preflight) | — |

**Solve URL path:** `POST /api/solve`

## Example: thesis point IPG6-S (2)

IPG6-S, 100 Pa, 23 MJ/kg assigned enthalpy, pure O₂, 37/20/40 mm nozzle.

```http
POST /api/solve
Content-Type: application/json
X-API-Key: <only if the host set API_KEY>
```

```json
{
  "pinj_Pa": 100,
  "hinj_MJ_kg": 23,
  "mixture": { "O2": 1.0 },
  "basis": "mole",
  "d_c_mm": 37,
  "d_t_mm": 20,
  "d_e_mm": 40,
  "nozzle_name": "IPG6-S",
  "xmax_m": 2.0,
  "ymax_m": 1.0,
  "nx": 65,
  "ny": 65
}
```

`nx` and `ny` are clamped to **17–80** (phone-safe grid). Sending 97×81 is accepted and capped at 80×80.

Optional `p_tank_Pa` (default **10.0**, `gt=0.05`, `lt=2e5`) is the ambient pressure used only for the continuum shock overlay. Old clients that omit it keep working.

Thesis PWA clients send `"plume_mode": "collisionless"` and do not rely on `p_tank_Pa` for the field. Advanced clients send `"auto"` or `"sudden_freeze"` and may send `p_tank_Pa`.

Response top keys: `cea` (stations, frozen `exit`, geometry, mixture, `mdot_mg_s`, `hinj_MJ_kg`, `delta_h_MJ_kg`, `converged`) and `plume` (`H`, `S0`, `T0`, `U0`, `n0`, `n_ratio`, `u`, `v`, `t_ratio`, `h_tot_MJ_kg`, `h_tot_ratio` arrays, plus `contours` and `probe`). Pointwise interpolation of the plume grid is still client-side; the optional calorimeter disk is `plume.probe` (JSON `null` if `probe_x_m` is omitted). CEA already reports `exit.p_Pa` (`p_e`); nozzle pressure ratio is `NPR = p_e / p_tank`.

### Input modes

`mode` on `/api/solve` (default `"enthalpy"`):

- `"enthalpy"` (alias `"point"`): set `pinj_Pa` and `hinj_MJ_kg` (assigned-enthalpy rocket problem). CEA returns `mdot_mg_s`. Clicking a point on the characteristics plot is this mode — no extra invert.
- `"generator"`: set `pinj_Pa` and `mdot_mg_s` (what the operator actually sets: injection pressure + MFC mass flow). The API **inverts CEA** for `hinj_MJ_kg` (mdot falls as enthalpy rises at fixed pinj). Do **not** send `power_W` for this; `hinj = P/ṁ` is the wrong inversion. Response includes inverted `cea.hinj_MJ_kg` / `cea.delta_h_MJ_kg` and `requested_mdot_mg_s`.

Generator example (IPG6-S O2, ~thesis point 2):

```json
{
  "mode": "generator",
  "pinj_Pa": 100,
  "mdot_mg_s": 13,
  "mixture": { "O2": 1.0 },
  "basis": "mole",
  "d_c_mm": 37,
  "d_t_mm": 20,
  "d_e_mm": 40,
  "nozzle_name": "IPG6-S",
  "nx": 49,
  "ny": 49
}
```

### Total enthalpy field

Frozen suggestion on the collisionless grid:

`h_static ≈ href + (h_exit − href) · (T/T0)` then `h_tot = h_static + ½(u²+v²)`.

Arrays: `plume.h_tot_MJ_kg`, `plume.h_tot_ratio` (`h_tot / hinj`). Contours under `plume.contours.h_tot`.

### Mixture preview

```http
POST /api/mixture-preview
Content-Type: application/json
```

```json
{ "mixture": { "N2": 0.79, "O2": 0.21 }, "basis": "mole" }
```

Returns `MW`, `R`, `h_ref_MJ_kg`, mole/mass fractions. Air as `{"Air": 1}` expands to N₂/O₂ 79/21; CEA is never sent the reactant name `"Air"`.

`GET /api/catalog` lists facilities, gas chips, thesis presets (IPG6-S (1)(2)(3), IPG3 O#01, IPG4 Burghaus).


## Characteristics field (`POST /api/characteristics`)

One hinj sweep at a single `pinj_ref` (~20–30 NASA CEA rocket calls, not a 2-D grid). Returns isolines for the (pinj, hinj) operating map and chamber/exit mole fractions vs hinj.

Approximation: at fixed hinj, composition and T are only weakly p-dependent, so `mdot ≈ k(h) · pinj` with `k(h) = mdot(h, pinj_ref) / pinj_ref`.

- Solid isolines: constant `mdot` (default 2, 5, 8, 13, 20, 30, 50 mg/s). `pinj(h) = mdot_target / k(h)`.
- Dashed isolines: coupled power `P = mdot · (hinj − href)` in W (default 50, 150, 300, 450, 600). `mdot = P / ((hinj − href)·1e6)`, then `pinj = mdot / k(h)`.
- Kinks on the mdot isolines are detected from the sweep (finite differences on chamber mole fractions): start/end of parent-molecule dissociation and start of ionization. For pure O2 they typically sit near ~2, ~21, ~28 MJ/kg. For N2/CO2/He mixes the parent is the dominant molecule (atomic gases skip dissociation kinks).
- Failed CEA hinj points are skipped.

Plot axes for the (pinj, hinj) map: pinj 0–250 Pa, hinj 0–40 MJ/kg. Companion composition plot: chamber mole fractions vs hinj at `pinj_ref` (composition is nearly independent of pinj).

Click a point `(pinj_Pa, hinj_MJ_kg)` then `POST /api/solve` with those two fields (`mode` omitted, `"enthalpy"`, or `"point"`). Same geometry and mixture as the sweep.

```http
POST /api/characteristics
Content-Type: application/json
```

```json
{
  "pinj_ref_Pa": 100,
  "mixture": { "O2": 1.0 },
  "basis": "mole",
  "d_c_mm": 37,
  "d_t_mm": 20,
  "d_e_mm": 40,
  "nozzle_name": "IPG6-S",
  "hinj_max": 40,
  "n_h": 29,
  "ions": true
}
```

Optional body fields: `gas`, `he_mole_frac` (legacy mix), `hinj_min` (default `href+0.2`), `mdot_mg_s_lines`, `power_W_lines`. `n_h` defaults to 29 and is capped at 41.

Response shape:

```json
{
  "pinj_ref_Pa": 100.0,
  "href_MJ_kg": 0.0,
  "geometry": { "name": "IPG6-S", "d_c_mm": 37.0, "d_t_mm": 20.0, "d_e_mm": 40.0 },
  "hinj_MJ_kg": [0.2, "..."],
  "chamber": {
    "T": ["..."],
    "MW": ["..."],
    "mdot_mg_s": ["..."],
    "x": { "O2": ["..."], "O": ["..."], "O+": ["..."], "e-": ["..."] }
  },
  "exit": {
    "T0": ["..."],
    "MW": ["..."],
    "x": { "O2": ["..."], "O": ["..."], "O+": ["..."], "e-": ["..."] }
  },
  "kinks": [
    { "hinj_MJ_kg": 2.0, "kind": "dissociation_start", "label": "O2 dissociation start" }
  ],
  "mdot_isolines": [
    { "mdot_mg_s": 13.0, "pinj_Pa": ["..."], "hinj_MJ_kg": ["..."] }
  ],
  "power_isolines": [
    { "power_W": 300.0, "pinj_Pa": ["..."], "hinj_MJ_kg": ["..."] }
  ],
  "axes": { "pinj_Pa": [0, 250], "hinj_MJ_kg": [0, 40] },
  "notes": [
    "mdot isolines use mdot ≈ k(h)·pinj at fixed hinj (composition weakly p-dependent).",
    "Kinks mark composition change: energy into dissociation/ionization, T and pinj rise more slowly."
  ]
}
```

`chamber.x` / `exit.x` always include O2, O, O+, e- and any other species whose mole fraction exceeds 1e-3 somewhere on the sweep. Extra keys `k_kg_s_Pa`, `ions`, `parent_molecule` are also returned.

## CORS

Allowed `Origin` values (reflected on the response):

- `https://*.grok.me` and `https://grok.me`
- `https://grok.com` and `https://*.grok.com`
- `http://localhost`, `http://localhost:<port>`
- `http://127.0.0.1`, `http://127.0.0.1:<port>`
- `https://*.onrender.com`

Preflight: `OPTIONS /api/*` with `Access-Control-Request-Method` / `Access-Control-Request-Headers`. Allowed request headers include `Content-Type` and `X-API-Key`. Other origins get no `Access-Control-Allow-Origin`.

Grok Build fetch sketch:

```js
const res = await fetch(`${API_BASE}/api/solve`, {
  method: "POST",
  headers: {
    "Content-Type": "application/json",
    ...(API_KEY ? { "X-API-Key": API_KEY } : {}),
  },
  body: JSON.stringify({
    pinj_Pa: 100,
    hinj_MJ_kg: 23,
    mixture: { O2: 1.0 },
    basis: "mole",
    d_c_mm: 37,
    d_t_mm: 20,
    d_e_mm: 40,
    nozzle_name: "IPG6-S",
    nx: 65,
    ny: 65,
  }),
});
```

## Environment

| Variable | Required | Meaning |
|---|---|---|
| `PORT` | set by Render/Railway | Listen port. Default **8765**. Binds `0.0.0.0`. |
| `API_KEY` | optional | If **set and non-empty**, every `/api/*` request except `GET /api/health` and `OPTIONS` must send header `X-API-Key` with that value. 401 otherwise. Leave unset for an open demo. |
| `PYTHONPATH` | Docker/Procfile set this to `backend` | So `uvicorn app:app` finds `app.py`. |

## Deploy

Docker (NASA CEA is installed and import-checked at build):

```bash
docker build -t irs-plume-api .
docker run --rm -p 8765:8765 -e PORT=8765 irs-plume-api
# optional: -e API_KEY=...
```

The image / Render start command:

```text
PYTHONPATH=backend python3 -m uvicorn app:app --app-dir backend --host 0.0.0.0 --port $PORT
```

Set `API_KEY` in the Render dashboard if you want to require `X-API-Key` — do not commit secrets. Confirm CEA with `GET /api/health` → `cea_version` like `3.3.3`.

## Tests

From the repo root:

```text
pip install -r requirements-test.txt
PYTHONPATH=backend python -m pytest tests -q
```

Most cases fixture a CEA exit state (no live NASA CEA). One smoke test runs `solve_operating_point` if `import cea` works. CI is `.github/workflows/test.yml`. Bug reports use `.github/ISSUE_TEMPLATE/bug.yml`.


### Plume mode

`plume_mode` on `/api/solve` (default `"auto"` when omitted, so old clients do not break). **Knudsen number is the only Auto trigger.**

`Kn_exit = λ/H` at the lip. `KN_CRIT = 0.05` (Boyd `Kn_GLL` / Bird P-order).

Phone / PWA layers (this API does not enforce them; the client chooses `plume_mode`):

- **Thesis** (PWA default): always send `"collisionless"` explicitly. Hard override: `p_tank_Pa` is ignored for the field (no barrel, no disk). `p_e_Pa`, `p_tank_Pa`, and `npr` may still be echoed as diagnostics if a tank pressure was sent.
- **Advanced**: send `"auto"` or `"sudden_freeze"`, plus optional `p_tank_Pa`. Auto uses the Kn trigger and may apply the shock overlay. `sudden_freeze` forces the continuum core and overlay (subject to the freeze/Kn veto).

Mode chips:

- `"collisionless"`: Khasawneh–Cai 2-D free-molecular jet from the exit slit. Hard override of Auto and of `p_tank`. No shock overlay, even if NPR is huge or `Kn_exit` is low.
- `"sudden_freeze"`: planar isentropic source flow (collisions on) until Boyd `Kn_GLL = λ/R` reaches 0.05, then translational T freezes and density continues as 1/R. Outside the Prandtl–Meyer vacuum cone the collisionless jet is used. Chip overrides Auto. Continuum core plus shock overlay if `p_tank_Pa` allows.
- `"auto"` (API default when the field is omitted): `Kn_exit >= 0.05` → collisionless (unchanged Khasawneh–Cai). `Kn_exit < 0.05` → continuum core (sudden-freeze source flow) **plus** the shock overlay below if `p_tank` allows. Response also has `plume.plume_mode_requested`.

Response already includes `plume.mode`, `plume.kn_gll_exit`, `plume.r_freeze_m`, `plume.kn_crit`.

### Shock overlay (continuum only)

Engineering model, not Euler / Navier–Stokes / DSMC. Applied only on the continuum path.

`NPR = p_e / p_tank` with `p_e = cea.exit.p_Pa`.

- **Underexpanded** (`p_e > p_tank`): Mach-disk station (axisymmetric) `x_m / D_e = 0.67 * sqrt(p_e / p_tank)` (Crist 1966 / Addy 1981 / Ashkenas–Sherman sonic-orifice family; `D_e = 2H`). Barrel is a smooth curve from the lip `(x=0, y=H)` to the triple point `(x_m, ~0.3–0.4 x_m)`, then a normal disk across `|y| < r_tp`. Rankine–Hugoniot jump on the disk (normal shock at the pre-disk Mach from the isentropic core). Downstream: subsonic, higher `n` and `T`, `U` drops.
- **Overexpanded** (`p_tank > p_e`): lip oblique shock from θ–β–M (pressure ratio → `β`). Core is deflected inward. If the wave would Mach-reflect, a small disk is placed; otherwise regular reflection.
- **Matched** (`|log NPR| < ~0.05`): no shock, sudden-freeze only.
- **Freeze veto:** if `r_freeze` is closer than `x_m`, or Kn at the would-be disk is `>= 0.05`, the shock is **not** applied (`shock_applied=false`, `shock_reason=freeze_before_disk`). Vacuum-like ambient pressures then look like today's sudden-freeze.

Optional-safe extras on `plume`:

| Field | Meaning |
|---|---|
| `p_tank_Pa` | Ambient pressure used for NPR |
| `p_e_Pa` | CEA exit pressure |
| `npr` | `p_e / p_tank` |
| `regime` | `underexpanded` \| `overexpanded` \| `matched` \| `vacuum` |
| `x_mach_disk_m` | Disk axial station, or `null` if no disk |
| `r_triple_m` | Triple-point radius, or `null` |
| `shock_applied` | `true` if the overlay mutated the field |
| `shock_reason` | Why it was applied or skipped (`freeze_before_disk`, `collisionless`, `matched`, …) |
| `barrel_xy` | `[[x, y], …]` lip → triple, upper half `y >= 0` (empty if no barrel) |
| `disk_y0`, `disk_y1` | Disk span at `x_mach_disk_m`, or `null` |
| `r_freeze_m`, `kn_gll_exit`, `kn_crit` | Already present |

Grid extras on `plume`:

- `mach`: local Mach \(U / \sqrt{\gamma R T}\) with \(T = (T/T_0) T_0\).
- `e_kin_eV`: directed particle kinetic energy \( \tfrac12 \bar m U^2 \) in eV (mixture-mean mass).
- `e_O_eV`: same for atomic oxygen, \( \tfrac12 m_O U^2 \).

### Optional sample disk (`plume.probe`)

All three body fields are optional. If `probe_x_m` is omitted, `plume.probe` is `null` (old clients unchanged).

| Field | Default | Limits | Meaning |
|---|---|---|---|
| `probe_x_m` | omitted | `> 0`, `< xmax_m` | Station of a circular disk on the axis (`y = 0`). Face in the *yz* plane (normal along *x*), intercepting the jet on the upstream face. |
| `probe_r_mm` | `20` | `2`–`80` | Disk radius. |
| `probe_Tw_K` | `300` | `200`–`2000` | Wall temperature (water-cooled calorimeter default). |

The incident state is sampled from the **existing** plume field at `(probe_x_m, 0)` *before* inserting a body: `n_inf`, `T_inf`, `U_inf`. The object Knudsen number is `kn_obj = λ / probe_r` at that state, with `λ = 1 / (√2 π d_hs² n)` and `KN_CRIT = 0.05`.

**Regime**

- `kinetic` (`model`: `khasawneh_diffuse`) if `plume.mode` is `collisionless` **or** `kn_obj >= 0.05`.
- `continuum` (`model`: `newtonian_billig`) if `kn_obj < 0.05` **and** the field is sudden-freeze / auto-continuum (`plume.mode == sudden_freeze`).

**Kinetic / free-molecular.** Khasawneh–Cai incident flux on the forward face of a drifting Maxwellian equal to the sampled state; fully diffuse re-emission at `Tw` with accommodation `α = 1`. Wall pressure `p_w` is incident plus re-emitted normal momentum. Heat `q_w` is incident translational energy minus wall re-emission. No chemistry catalysis in v1. No bow shock is drawn on this path. Collisionless geometric shadow (zero / incident-only downstream of the disk in `|y| < R`, `x > x_p`) is **omitted in v1**; the returned `n_ratio`/`u`/`v`/`t_ratio` arrays are the undisturbed field.

**Continuum.** Engineering blunt-face closure, not Euler/NS/DSMC:

- Billig (1967) sphere-like standoff on a blunt face: `Δ/R ≈ 0.143 exp(3.24/M²)`, plus the hyperbolic shock with vertex radius `R_c/R = 1.143 exp[0.54/(M−1)^{1.2}]`. Optional draw-data `bow_xy` is the upper-half polyline in plume `(x, y)`.
- Modified Newtonian on the forward face: `p = p_inf + (p_t2 − p_inf) cos²θ`. Stagnation `p_stag` is the Rayleigh–Pitot `p_t2`; for a flat face `θ = 0` so the area-average `p_w` equals `p_stag`.
- Heat: Sutton–Graves stagnation `q_stag` for a sphere of radius `R` and cold wall `Tw`, `q = K (p_s/R)^{1/2} (h_s − h_w)` with `p_s` in atm (NASA TR R-376 eq. 33). Face average `q_w = (2/3) q_stag` (frontal-area mean of a `q/q_s = cos θ` distribution on an equivalent spherical nose). Not a 2-D NS body.

Response `plume.probe` (or `null`):

```json
{
  "x_m": 0.4,
  "r_mm": 20.0,
  "Tw_K": 300.0,
  "kn_obj": 1.2,
  "regime": "kinetic",
  "n_inf": 1.0e18,
  "T_inf": 800.0,
  "U_inf": 2500.0,
  "p_w_Pa": 12.3,
  "q_w_W_m2": 4.5e4,
  "p_stag_Pa": null,
  "q_stag_W_m2": null,
  "bow_xy": [],
  "model": "khasawneh_diffuse",
  "notes": ["Translational heat only; no chemistry catalysis in v1."]
}
```

`p_stag_Pa` and `q_stag_W_m2` are continuum-only (`null` in kinetic). `bow_xy` is `[[x, y], ...]` or `[]`. `T0` on the plume is nozzle-exit translational temperature.

Example with a disk:

```json
{
  "pinj_Pa": 100,
  "hinj_MJ_kg": 23,
  "mixture": { "O2": 1.0 },
  "d_c_mm": 37,
  "d_t_mm": 20,
  "d_e_mm": 40,
  "plume_mode": "collisionless",
  "probe_x_m": 0.4,
  "probe_r_mm": 20,
  "probe_Tw_K": 300,
  "nx": 33,
  "ny": 33
}
```

### Citations

- Khasawneh, K., Liu, H., and Cai, C., “Highly rarefied two-dimensional jet impingement on a flat plate,” *Phys. Fluids* **22**, 117101 (2010). doi:[10.1063/1.3490409](https://doi.org/10.1063/1.3490409)
- Cai, C. and Boyd, I. D., “Theoretical and Numerical Study of Free Molecular-Flow Problems,” *J. Spacecraft and Rockets* **44**(3) 619–624 (2007). doi:[10.2514/1.25893](https://doi.org/10.2514/1.25893)
- Billig, F. S., “Shock-wave shapes around spherical- and cylindrical-nosed bodies,” *J. Spacecraft and Rockets* **4**(6) 822–823 (1967). doi:[10.2514/3.28969](https://doi.org/10.2514/3.28969)
- Sutton, K. and Graves, R. A., Jr., “A general stagnation-point convective-heating equation for arbitrary gas mixtures,” NASA TR R-376 (1971). [NTRS 19720003329](https://ntrs.nasa.gov/citations/19720003329)

## References

- NASA CEA (Gordon & McBride / CEA2)
- Bird, G. A. (1970). Breakdown of translational and rotational equilibrium in gaseous expansions. *AIAA J.* 8(11), 1998–2003. doi:10.2514/3.6037
- Bird, G. A. (1994). *Molecular Gas Dynamics and the Direct Simulation of Gas Flows*. Oxford.
- Boyd, I. D., Chen, G., & Candler, G. V. (1995). Predicting failure of the continuum fluid equations in transitional hypersonic flows. *Phys. Fluids* 7, 210–219. doi:10.1063/1.868720
- Cai, C., & Boyd, I. D. (2007). Theoretical and numerical study of free molecular-flow problems. *J. Spacecraft Rockets* 44(3), 619–624. doi:10.2514/1.25893
- Khasawneh, K. R., Liu, H., & Cai, C. (2010). Highly rarefied two-dimensional jet impingement on a flat plate. *Phys. Fluids* 22. doi:10.1063/1.3490409
- Billig, F. S. (1967). Shock-wave shapes around spherical- and cylindrical-nosed bodies. *J. Spacecraft Rockets* 4(6), 822–823. doi:10.2514/3.28969
- Sutton, K., & Graves, R. A., Jr. (1971). A general stagnation-point convective-heating equation for arbitrary gas mixtures. NASA TR R-376. https://ntrs.nasa.gov/citations/19720003329
- Crist, S., Sherman, P. M., & Glass, D. R. (1966). Study of the highly underexpanded sonic jet. *AIAA J.* 4(1), 68–71. doi:10.2514/3.3386
- Addy, A. L. (1981). Effects of axisymmetric sonic nozzle geometry on Mach disk characteristics. *AIAA J.* 19(1), 121–122. doi:10.2514/3.7751
- Ashkenas, H., & Sherman, F. S. (1966). The structure and utilization of supersonic free jets in low density wind tunnels. *Rarefied Gas Dynamics*.
- Albini, F. A. (1965). Approximate computation of underexpanded jet structure. *AIAA J.* 3(8), 1535–1537. doi:10.2514/3.3194
- Boynton, F. P. (1967). Highly underexpanded jet structure — exact and approximate calculations. *AIAA J.* 5(9), 1703–1704. doi:10.2514/3.4283
- Dettleff, G. (1991). Plume flow and plume impingement in space technology. *Prog. Aerosp. Sci.* 28, 1–71. doi:10.1016/0376-0421(91)90008-R
- Rankine–Hugoniot relations; Prandtl–Meyer turning (standard)


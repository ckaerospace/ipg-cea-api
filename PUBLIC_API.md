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

Response top keys: `cea` (stations, frozen `exit`, geometry, mixture, `mdot_mg_s`, `hinj_MJ_kg`, `delta_h_MJ_kg`, `converged`) and `plume` (`H`, `S0`, `T0`, `U0`, `n0`, `n_ratio`, `u`, `v`, `t_ratio`, `h_tot_MJ_kg`, `h_tot_ratio` arrays, plus `contours`). There is no top-level `probe`; interpolate the plume grid client-side.

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


### Plume mode

`plume_mode` on `/api/solve` (default `"auto"`):

- `"collisionless"`: Khasawneh–Cai 2-D free-molecular jet from the exit slit (original).
- `"sudden_freeze"`: planar isentropic source flow (collisions on) until Boyd `Kn_GLL = λ/R` reaches 0.05, then translational T freezes and density continues as 1/R. Outside the Prandtl–Meyer vacuum cone the collisionless jet is used. Response `plume.mode`, `plume.kn_gll_exit`, `plume.r_freeze_m`.

- `"auto"` (default): if `Kn_exit = λ/H` at the lip is below 0.05, use sudden-freeze; otherwise collisionless. Response also has `plume.plume_mode_requested`.

Grid extras on `plume`:

- `mach`: local Mach \(U / \sqrt{\gamma R T}\) with \(T = (T/T_0) T_0\).
- `e_kin_eV`: directed particle kinetic energy \( \tfrac12 \bar m U^2 \) in eV (mixture-mean mass).
- `e_O_eV`: same for atomic oxygen, \( \tfrac12 m_O U^2 \). Probe also returns `e_th_eV` = \( \tfrac32 kT \).


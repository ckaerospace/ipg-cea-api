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

- `"enthalpy"`: set `pinj_Pa` and `hinj_MJ_kg` (thesis assigned-enthalpy rocket problem). CEA returns `mdot_mg_s`.
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

## CORS

Allowed `Origin` values (reflected on the response):

- `https://*.grok.me` and `https://grok.me`
- `https://grok.com` and `https://*.grok.com`
- `http://localhost`, `http://localhost:<port>`
- `http://127.0.0.1`, `http://127.0.0.1:<port>`

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

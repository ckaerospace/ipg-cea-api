"""
2-D collisionless free-molecular jet into vacuum.

Primary sources (trusted over thesis typesetting if they disagree):
  Cai & Boyd, J. Spacecraft and Rockets 44(3) 2007 (thesis ref. 42)
  Khasawneh, Cai & Wei, AIAA 2010-0986 / related gaskinetic papers
  Cai, Aerospace 4(1) 5 (2017) — compact 2-D jet formulae
  Cai & Cai, Fluids 6(7) 250 (2021) — A(t), B(t), C(t) integrand form

The operating-point notes eqs. (4.1)–(4.9) match the Fluids 2021
integrand form after expanding A, B, C.  HTML conversions of some
open-access papers drop √π vs π; the stable A/B/C form is used here.

Geometry: planar slit of half-height H.  Applied to the IPG6-S round
exit as in the thesis (H = D_E / 2).  Angles:

    θ1 = atan2(Y - H, X),   θ2 = atan2(Y + H, X)

Speed ratio:  S0 = U0 / sqrt(2 R T0)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import erf, erfcx

SQRT_PI = np.sqrt(np.pi)
SQRT_2 = np.sqrt(2.0)

# Default Gauss–Legendre order.  64 is well past contour-plot convergence
# for S0 up to ~6; tests use 96.
_DEFAULT_QUAD = 64


def fig48_speed_ratio() -> float:
    """S0 for thesis Fig. 4.8: U0 = sqrt(5/3 R T0) ⇒ S0 = sqrt(5/6)."""
    return float(np.sqrt(5.0 / 6.0))


def _w(S0: float, theta: np.ndarray) -> np.ndarray:
    """exp(-S0² sin²θ) (1 + erf(S0 cos θ)), overflow-safe."""
    theta = np.asarray(theta, dtype=np.float64)
    s = np.sin(theta)
    t = S0 * np.cos(theta)
    S0sq = S0 * S0
    out = np.empty_like(theta)
    pos = t >= 0.0
    if np.any(pos):
        out[pos] = np.exp(-S0sq * s[pos] ** 2) * (1.0 + erf(t[pos]))
    neg = ~pos
    if np.any(neg):
        # exp(-S0² sin²θ) (1+erf(t)) = exp(-S0²) erfcx(-t)  for t < 0
        out[neg] = np.exp(-S0sq) * erfcx(-t[neg])
    return out


def _integrands(S0: float, theta: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (I_A, I_Bcos, I_C) = exp(-S0²) × {A, B cosθ, C}.

    A(t) = 1 + √π t e^{t²} [1+erf(t)]
    B(t) = t + √π (1/2 + t²) e^{t²} [1+erf(t)]
    C(t) = 3/4 + t²/2 + √π (t + t³/2) e^{t²} [1+erf(t)]
    with t = S0 cos θ.

    Multiplying by exp(-S0²) up front keeps every term bounded.
    """
    c = np.cos(theta)
    t = S0 * c
    S0sq = S0 * S0
    eS = np.exp(-S0sq)
    w = _w(S0, theta)  # exp(-S0² sin²θ) (1+erf(t))

    I_A = eS + SQRT_PI * t * w
    I_B = t * eS + SQRT_PI * (0.5 + t * t) * w
    I_C = eS * (0.75 + 0.5 * t * t) + SQRT_PI * (t + 0.5 * t**3) * w
    I_Bcos = I_B * c
    return I_A, I_Bcos, I_C


def _gamma3(S0: float, theta: np.ndarray) -> np.ndarray:
    """Thesis (4.9): γ3(θ) = exp(-S0² sin²θ) cosθ [1+erf(S0 cosθ)]."""
    return _w(S0, theta) * np.cos(theta)


@dataclass(frozen=True)
class PlumePoint:
    n_ratio: float
    u: float  # m/s, streamwise (X)
    v: float  # m/s, transverse (Y)
    t_ratio: float
    speed: float


class CollisionlessPlume:
    """Khasawneh–Cai 2-D collisionless jet from a slit of half-height H."""

    def __init__(
        self,
        S0: float,
        T0: float,
        R_specific: float,
        U0: float,
        n0: float,
        H: float,
        quad: int = _DEFAULT_QUAD,
    ) -> None:
        if S0 < 0:
            raise ValueError("S0 must be non-negative")
        if H <= 0:
            raise ValueError("H must be positive")
        self.S0 = float(S0)
        self.T0 = float(T0)
        self.R = float(R_specific)
        self.U0 = float(U0)
        self.n0 = float(n0)
        self.H = float(H)
        self.quad = int(quad)
        self._xi, self._wi = leggauss(self.quad)
        self.thermal = SQRT_2 * np.sqrt(self.R * self.T0)  # sqrt(2 R T0)
        # n/n0 floor avoids 1/n blow-ups in empty vacuum
        self._n_floor = 1e-16

    @classmethod
    def from_exit(
        cls,
        T0: float,
        R_specific: float,
        U0: float,
        n0: float,
        H: float,
        quad: int = _DEFAULT_QUAD,
    ) -> "CollisionlessPlume":
        thermal = SQRT_2 * np.sqrt(R_specific * T0)
        S0 = U0 / thermal if thermal > 0 else 0.0
        return cls(S0, T0, R_specific, U0, n0, H, quad=quad)

    def angles(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """θ1, θ2 via atan2 so X → 0+ is well defined."""
        x = np.asarray(x, dtype=np.float64)
        y = np.asarray(y, dtype=np.float64)
        # Keep a tiny positive X so the exit plane is approached from downstream.
        x_safe = np.maximum(x, 1e-12 * self.H)
        th1 = np.arctan2(y - self.H, x_safe)
        th2 = np.arctan2(y + self.H, x_safe)
        return th1, th2

    def _integrate(self, th1: np.ndarray, th2: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Gauss–Legendre of (I_A, I_B cosθ, I_C) on [θ1, θ2]."""
        th1 = np.asarray(th1, dtype=np.float64)
        th2 = np.asarray(th2, dtype=np.float64)
        half = 0.5 * (th2 - th1)
        mid = 0.5 * (th2 + th1)
        theta = half[..., None] * self._xi + mid[..., None]
        wgt = half[..., None] * self._wi
        I_A, I_Bcos, I_C = _integrands(self.S0, theta)
        return (
            np.sum(wgt * I_A, axis=-1),
            np.sum(wgt * I_Bcos, axis=-1),
            np.sum(wgt * I_C, axis=-1),
        )

    def evaluate(self, x: np.ndarray, y: np.ndarray) -> dict[str, np.ndarray]:
        """Vectorized field at arrays x, y (broadcastable).

        Returns dict of n_ratio, u, v, t_ratio, speed (same shape as broadcast).
        """
        x_arr, y_arr = np.broadcast_arrays(
            np.asarray(x, dtype=np.float64),
            np.asarray(y, dtype=np.float64),
        )
        shape = x_arr.shape
        th1, th2 = self.angles(x_arr, y_arr)
        int_A, int_Bcos, int_C = self._integrate(th1.ravel(), th2.ravel())

        n_ratio = (int_A / (2.0 * np.pi)).reshape(shape)
        n_ratio = np.maximum(n_ratio, 0.0)
        n_safe = np.maximum(n_ratio, self._n_floor)
        inv_n = 1.0 / n_safe

        # U / sqrt(2RT0) = [exp(-S0²)/(2π)] (n0/n) ∫ B cosθ dθ
        # I_B already includes exp(-S0²), so:
        u_hat = (inv_n / (2.0 * np.pi)) * int_Bcos.reshape(shape)

        # Closed form (4.4) / Fluids (3): numerically stabler than ∫ B sinθ
        g1 = _gamma3(self.S0, th1)
        g2 = _gamma3(self.S0, th2)
        v_hat = (inv_n / (4.0 * SQRT_PI)) * (g1 - g2)

        u = u_hat * self.thermal
        v = v_hat * self.thermal

        # T/T0 = (2/3) (n0/n) / π  ∫ I_C dθ  − (U²+V²)/(3 R T0)
        kinetic = (u * u + v * v) / (3.0 * self.R * self.T0)
        t_ratio = (2.0 / 3.0) * inv_n / np.pi * int_C.reshape(shape) - kinetic
        # Rare empty-vacuum noise; translational T cannot be negative.
        t_ratio = np.clip(t_ratio, 0.0, 5.0)

        # Downstream only; mask points inside the nozzle (X<0) as empty.
        inside = x_arr < 0.0
        if np.any(inside):
            n_ratio = n_ratio.copy()
            u = u.copy()
            v = v.copy()
            t_ratio = t_ratio.copy()
            n_ratio[inside] = 0.0
            u[inside] = 0.0
            v[inside] = 0.0
            t_ratio[inside] = 0.0

        speed = np.hypot(u, v)
        return {
            "n_ratio": n_ratio,
            "u": u,
            "v": v,
            "t_ratio": t_ratio,
            "speed": speed,
        }

    def point(self, x: float, y: float) -> PlumePoint:
        f = self.evaluate(np.array(x), np.array(y))
        return PlumePoint(
            n_ratio=float(f["n_ratio"]),
            u=float(f["u"]),
            v=float(f["v"]),
            t_ratio=float(f["t_ratio"]),
            speed=float(f["speed"]),
        )

    def grid(
        self,
        x: Iterable[float],
        y: Iterable[float],
    ) -> dict[str, np.ndarray]:
        """Evaluate on a tensor-product mesh.  x, y are 1-D axes.

        Returns 2-D arrays indexed [iy, ix] (row = Y, col = X) plus the axes.
        """
        x_ax = np.asarray(list(x), dtype=np.float64)
        y_ax = np.asarray(list(y), dtype=np.float64)
        xx, yy = np.meshgrid(x_ax, y_ax, indexing="xy")
        field = self.evaluate(xx, yy)
        field["x"] = x_ax
        field["y"] = y_ax
        field["X"] = xx
        field["Y"] = yy
        return field


def marching_squares(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    levels: Iterable[float],
) -> dict[str, list[list[list[float]]]]:
    """Lightweight contour extractor.  Returns {level: [[x,y], ...] polylines}."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    cs = plt.contour(x, y, z, levels=sorted(set(float(v) for v in levels)))
    out: dict[str, list[list[list[float]]]] = {}
    try:
        all_paths = cs.get_paths() if hasattr(cs, "get_paths") else None
        all_levels = cs.levels
        # matplotlib >= 3.8: cs.allsegs
        segs = getattr(cs, "allsegs", None)
        if segs is not None:
            for lev, paths in zip(all_levels, segs):
                polylines = []
                for p in paths:
                    if len(p) >= 2:
                        polylines.append(p.astype(float).tolist())
                out[f"{lev:g}"] = polylines
        elif all_paths is not None:
            # single collection fallback
            out["all"] = [p.vertices.astype(float).tolist() for p in all_paths]
    finally:
        plt.close("all")
    return out

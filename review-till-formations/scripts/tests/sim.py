"""Simulated counties with known answers, in the same shape io.prepare returns."""
import numpy as np
import pandas as pd

from tillreview import variables as V
from tillreview.io import Data

# Broad-group proportions in the order of V.LEAVES:
# light, dark, red, pc_other, carbonate, pz_other, shale, k_other
SHALY = [0.36, 0.05, 0.01, 0.00, 0.22, 0.01, 0.30, 0.05]
CRYSTAL = [0.62, 0.20, 0.12, 0.01, 0.05, 0.00, 0.00, 0.00]
CARB = [0.50, 0.08, 0.02, 0.00, 0.38, 0.01, 0.005, 0.005]
UNITS = {"A": (SHALY, [42.0, 20.0]), "B": (CRYSTAL, [60.0, 14.0]), "C": (CARB, [33.0, 21.0])}


def make(holes, units=UNITS, n_grains=250, AL=40.0, AW=200.0, sd_between=2.0, sd_within=2.0, seed=1):
    """holes: list of boreholes, each a list of (true unit, coded unit or None, n samples), top down.

    Each borehole draws its own composition around the unit's (spread AL, sd_between); each
    sample draws around its borehole's (spread AW, sd_within). n_grains may be a number or a
    function of the sample index.
    """
    rng = np.random.default_rng(seed)
    rows, leaf, tex = [], [], []
    for h, runs in enumerate(holes):
        depth = 5.0
        local = {}
        for true, coded, n in runs:
            pi, mu = units[true]
            if true not in local:
                local[true] = (rng.dirichlet(AL * np.maximum(np.array(pi), 1e-4)),
                               np.array(mu) + rng.normal(0, sd_between, 2))
            theta, m = local[true]
            for _ in range(n):
                ng = n_grains(len(rows)) if callable(n_grains) else n_grains
                leaf.append(rng.multinomial(ng, rng.dirichlet(AW * np.maximum(theta, 1e-5))))
                tex.append(np.round(m + rng.normal(0, sd_within, 2)))
                rows.append({"hole": f"H{h:02d}", "true": true, "code": coded, "depth": depth})
                depth += 5.0
    S = pd.DataFrame(rows)
    leaf = np.array(leaf, float)
    raw = pd.DataFrame(0.0, index=S.index, columns=V.GRAIN_COLS)
    for k, g in enumerate(V.LEAVES):
        raw[V.GROUPS[g][0]] = leaf[:, k]
    hid = S.hole.str[1:].astype(int)
    S = S.assign(row=np.arange(len(S)) + 2, sample_num=np.arange(len(S)), stratunit=S.code, in_county=True,
                 x=1000.0 * hid, y=0.0, surface_elev=1000.0, deposittype="T", elev=1000.0 - S.depth, till=True,
                 n_grains=leaf.sum(axis=1), nonprov_share=0.0, leached="N", wet_color="2.5Y 4/2", dry_color="",
                 counted_by="X", depth_top=np.nan, depth_bottom=np.nan, secondary=0.0, incomplete=0, type="RS",
                 rotasonic=True, gravel_pct=3.0, silt_pct=0.0, notes="", bad_depth=False, missing_xy=False,
                 zero_grains=False, no_texture=False, low_count=leaf.sum(axis=1) < 100, unusable=False,
                 formation=S.code, level="formation")
    return Data(S, raw, leaf, np.array(tex, float), [])


def units_of(data):
    return np.array([c if isinstance(c, str) else None for c in data.S.code], dtype=object)

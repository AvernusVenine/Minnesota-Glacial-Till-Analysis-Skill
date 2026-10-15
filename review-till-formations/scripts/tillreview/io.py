"""Loading, the area under review (county + buffer) and step 1 (bad samples)."""
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from . import variables as V

TILL = ("T", "T*")
FT_PER_M = 3.28084
NOTE_COLS = ["notes1", "notes2", "notes3", "notes4", "notes5", "grain_note1", "grain_note2",
             "grain_note3", "grain_note4", "previnterpt", "strat_work"]


@dataclass
class Data:
    """Samples of the area under review. S has one row per sample, index 0..N-1."""
    S: pd.DataFrame
    raw: pd.DataFrame      # grain columns, numeric, NaN -> 0
    leaf: np.ndarray       # (N, 8) counts of the broad groups
    tex: np.ndarray        # (N, 2) sand %, clay %; NaN where there is no texture
    codes_not_in_lookup: list

    def subset(self, mask):
        idx = np.flatnonzero(mask)
        return Data(self.S.iloc[idx].reset_index(drop=True), self.raw.iloc[idx].reset_index(drop=True),
                    self.leaf[idx], self.tex[idx], self.codes_not_in_lookup)


def load_csv(path):
    df = pd.read_csv(path, low_memory=False, dtype={"relateid": str, "stratunit": str})
    return df


def load_codes(path):
    codes = pd.read_csv(path, dtype=str).fillna("")
    return codes.set_index("code")


def check_rollups(df, tol=0.5):
    """Share of rows where the groups reproduce each QDI rollup column (M1)."""
    out = {}
    for col, groups in V.ROLLUPS.items():
        if col not in df.columns:
            continue
        mine = df[[c for g in groups for c in V.GROUPS[g]]].apply(pd.to_numeric, errors="coerce").fillna(0).sum(axis=1)
        theirs = pd.to_numeric(df[col], errors="coerce")
        ok = theirs.notna()
        out[col] = float(((mine[ok] - theirs[ok]).abs() <= tol).mean()) if ok.any() else None
    return out


def prepare(df, county, codes, buffer_km=25.0, low_count=100):
    """Rows of every borehole in `county` or within buffer_km of one of its boreholes.

    Adds the columns used downstream and the step-1 flags. Nothing is dropped here except
    boreholes outside the area.
    """
    df = df.copy()
    for c in ("POINT_X", "POINT_Y", "sampledepth", "elevation"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    holes = df.groupby("relateid").agg(x=("POINT_X", "first"), y=("POINT_Y", "first"),
                                       cty=("cty_name", "first"))
    inside = holes[holes.cty == county]
    if inside.empty:
        raise ValueError(f"no boreholes with cty_name = {county!r}")
    located = inside.dropna(subset=["x", "y"])
    other = holes[(holes.cty != county)].dropna(subset=["x", "y"])
    near = []
    if buffer_km > 0 and len(located) and len(other):
        d, _ = cKDTree(located[["x", "y"]].to_numpy()).query(other[["x", "y"]].to_numpy())
        near = list(other.index[d <= buffer_km * 1000.0])
    S = df[df.relateid.isin(set(inside.index) | set(near))]
    csv_row = S.index.to_numpy() + 2          # line in the CSV file (header = line 1)
    S = S.reset_index(drop=True)

    raw = S[V.GRAIN_COLS].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    leaf = V.leaf_counts(raw)
    out = pd.DataFrame({
        "row": csv_row, "sample_num": S.sample_num if "sample_num" in S.columns else "",
        "hole": S.relateid.astype(str), "stratunit": S.stratunit,
        "in_county": (S.cty_name == county).to_numpy(), "x": S.POINT_X, "y": S.POINT_Y,
        "surface_elev": S.elevation, "deposittype": S.deposittype.fillna("").astype(str),
    })
    depth = S.sampledepth.copy()
    metres = S.depth_unit.astype(str).str.lower().str.startswith("m")
    depth[metres & (depth >= 0)] *= FT_PER_M
    out["depth"] = depth
    out["elev"] = out.surface_elev - out.depth
    out["till"] = out.deposittype.isin(TILL)
    out["n_grains"] = leaf.sum(axis=1)
    nonprov = S[V.NON_PROVENANCE].apply(pd.to_numeric, errors="coerce").fillna(0).sum(axis=1)
    out["nonprov_share"] = (nonprov / (nonprov + out.n_grains).replace(0, np.nan)).fillna(0.0)
    for c in ["leached", "wet_color", "dry_color", "counted_by", "depth_top", "depth_bottom",
              "secondary", "incomplete", "type", "rotasonic", "gravel_pct", "silt_pct"]:
        out[c] = S[c] if c in S.columns else np.nan
    out["notes"] = S[[c for c in NOTE_COLS if c in S.columns]].fillna("").astype(str).agg(" ".join, axis=1).str.lower()

    sand = pd.to_numeric(S.sand_pct, errors="coerce")
    silt = pd.to_numeric(S.silt_pct, errors="coerce")
    clay = pd.to_numeric(S.clay_pct, errors="coerce")
    total = sand + silt + clay
    has_tex = total.notna() & ((total - 100).abs() <= 2)   # all-zero rows were never measured
    tex = np.column_stack([sand.where(has_tex), clay.where(has_tex)]).astype(float)

    # Step 1 flags
    out["bad_depth"] = ~(out.depth >= 0)
    out["missing_xy"] = out.x.isna() | out.y.isna() | (out.x == 0) | (out.y == 0)
    out["zero_grains"] = out.n_grains == 0
    out["no_texture"] = ~has_tex.to_numpy()
    out["low_count"] = (out.n_grains > 0) & (out.n_grains < low_count)
    out["unusable"] = out.bad_depth | out.missing_xy | (out.zero_grains & out.no_texture)

    code = out.stratunit.fillna("").astype(str).str.strip().str.upper()
    known = code.isin(codes.index)
    out["code"] = code.where(code != "", None)
    out["formation"] = [codes.formation_code[c] if k else (c or None) for c, k in zip(code, known)]
    out["level"] = [codes.level[c] if k else ("formation" if c else None) for c, k in zip(code, known)]
    missing = sorted(set(code[(code != "") & ~known & out.till & out.in_county]))
    return Data(out, raw, leaf, tex, missing)


def bad_sample_table(data):
    S = data.S[data.S.in_county & data.S.till]
    rows = []
    for flag, effect in [("bad_depth", "excluded: depth negative or missing"),
                         ("missing_xy", "excluded: no coordinates"),
                         ("zero_grains", "texture only: no grains counted"),
                         ("no_texture", "grains only: texture missing or not summing to 100 +/- 2")]:
        for _, r in S[S[flag]].iterrows():
            rows.append({"row": r.row, "hole": r.hole, "depth_ft": r.depth, "stratunit": r.code,
                         "problem": flag, "effect": effect})
    return pd.DataFrame(rows, columns=["row", "hole", "depth_ft", "stratunit", "problem", "effect"])

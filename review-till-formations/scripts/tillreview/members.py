"""M13: which individual grain types join the member model of one formation."""
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import chi2

from . import variables as V
from .model import LOG_A_MAX, LOG_A_MIN, dm_logpmf


def _bb_loglik(x, n, p):
    """Beta-binomial log likelihood at the best shared concentration, for per-sample means p."""
    xx = np.column_stack([x, n - x])
    pp = np.clip(np.column_stack([p, 1 - p]), 1e-6, 1)
    res = minimize_scalar(lambda t: -dm_logpmf(xx, np.exp(t) * pp).sum(), bounds=(LOG_A_MIN, LOG_A_MAX),
                          method="bounded")
    return -res.fun


def bh(p):
    """Benjamini-Hochberg q-values."""
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = np.minimum.accumulate((p[o] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    out = np.empty(len(p))
    out[o] = np.minimum(q, 1.0)
    return out


def select_grains(raw, unit, idxs, q_max=0.05):
    """Test each grain type (as a share of its broad group) for a difference between members.

    Returns the extra branches for the member model and a table with every grain tested.
    """
    idxs = np.asarray(idxs, int)
    members = np.asarray(unit, dtype=object)[idxs]
    rows = []
    for group, grains in V.GROUPS.items():
        if len(grains) < 2:
            continue
        n = raw[grains].to_numpy(float)[idxs].sum(axis=1)
        for g in grains:
            x = raw[g].to_numpy(float)[idxs]
            keep = n > 0
            if x[keep].sum() == 0 or x[keep].sum() == n[keep].sum() or len(set(members[keep])) < 2:
                continue
            xk, nk, mk = x[keep], n[keep], members[keep]
            p0 = np.full(len(xk), (xk.sum() + 0.5) / (nk.sum() + 1))
            p1 = np.array([(xk[mk == m].sum() + 0.5) / (nk[mk == m].sum() + 1) for m in mk])
            G = max(0.0, 2 * (_bb_loglik(xk, nk, p1) - _bb_loglik(xk, nk, p0)))
            row = {"grain": g, "group": group, "G": G, "p": float(chi2.sf(G, len(set(mk)) - 1)),
                   "confusable": V.is_confusable(g)}
            for m in sorted(set(mk)):      # the grain as a percentage of its group, per member
                row[f"pct_in_{m}"] = round(100 * xk[mk == m].sum() / max(nk[mk == m].sum(), 1), 1)
            rows.append(row)
    table = pd.DataFrame(rows)
    if table.empty:
        return [], table
    table["q"] = bh(table.p)
    table["selected"] = table.q < q_max
    nodes = []
    for group in V.GROUPS:
        sel = list(table.grain[(table.group == group) & table.selected])
        if len(sel) == len(V.GROUPS[group]):
            sel = sel[:-1]          # the last grain is determined by the others
        if sel:
            nodes.append(V.member_node(group, sel))
    return nodes, table.sort_values("q").reset_index(drop=True)

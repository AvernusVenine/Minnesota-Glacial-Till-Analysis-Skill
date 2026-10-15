"""M3, M4: the one order every borehole must follow, and the samples that break it."""
from itertools import combinations

import numpy as np
import pandas as pd


def precedence(hole, unit, depth):
    """W[a][b] = boreholes where unit a's median depth is shallower than unit b's."""
    df = pd.DataFrame({"hole": hole, "unit": unit, "depth": depth}).dropna()
    med = df.groupby(["hole", "unit"]).depth.median()
    units = sorted(df.unit.unique())
    W = pd.DataFrame(0, index=units, columns=units)
    for _, m in med.groupby(level=0):
        m = m.droplevel(0)
        for a, b in combinations(m.index, 2):
            if m[a] < m[b]:
                W.loc[a, b] += 1
            elif m[b] < m[a]:
                W.loc[b, a] += 1
    return W, df.groupby("unit").depth.median()


def best_order(W, tiebreak=None):
    """Order (top first) that the fewest boreholes violate. Exact for <= 16 units.

    tiebreak: a value per unit (smaller = higher in the order) used only between orders
    with the same number of violations, e.g. median depth.
    """
    units = list(W.index)
    k = len(units)
    w = W.to_numpy(float)
    if tiebreak is not None:
        t = np.array([tiebreak.get(u, 0.0) for u in units], float)
        t = (t - t.min()) / (np.ptp(t) + 1e-9)
        w = w + 1e-6 * (t[:, None] > t[None, :])     # tiny cost for a deeper unit placed above
    if k <= 16:
        best = {0: (0.0, [])}
        for mask in range(1, 1 << k):
            cand = None
            for v in range(k):
                if mask >> v & 1:
                    prev = mask ^ (1 << v)
                    # v goes below every unit in prev: a violation wherever v was seen above one
                    cost = best[prev][0] + sum(w[v, a] for a in range(k) if prev >> a & 1)
                    if cand is None or cost < cand[0] - 1e-12:
                        cand = (cost, best[prev][1] + [v])
            best[mask] = cand
        order = best[(1 << k) - 1][1]
    else:
        order = list(np.argsort(-(w.sum(axis=1) - w.sum(axis=0)), kind="stable"))
        improved = True
        while improved:
            improved = False
            for i in range(k - 1):
                a, b = order[i], order[i + 1]
                if w[b, a] > w[a, b]:
                    order[i], order[i + 1] = b, a
                    improved = True
    order = [units[i] for i in order]
    viol = int(sum(W.loc[b, a] for i, a in enumerate(order) for b in order[i + 1:]))
    return order, viol


def pair_table(W, order):
    rows = []
    for i, a in enumerate(order):
        for b in order[i + 1:]:
            agree, against = int(W.loc[a, b]), int(W.loc[b, a])
            rows.append({"upper": a, "lower": b, "adjacent": order.index(b) == i + 1,
                         "boreholes_agree": agree, "boreholes_disagree": against,
                         "order_observed": agree + against > 0})
    return pd.DataFrame(rows)


def violations(ranks):
    """M4 for one borehole. ranks: timeline rank of each coded sample, top to bottom.

    Returns a status per sample: 'ok' (in every largest order-consistent subset),
    'violation' (in none) or 'ambiguous' (in some: either it or another sample is wrong).
    """
    r = list(ranks)
    n = len(r)
    f, cf = [1] * n, [1] * n          # longest non-decreasing run ending at i, and how many
    for i in range(n):
        for j in range(i):
            if r[j] <= r[i]:
                if f[j] + 1 > f[i]:
                    f[i], cf[i] = f[j] + 1, cf[j]
                elif f[j] + 1 == f[i]:
                    cf[i] += cf[j]
    b, cb = [1] * n, [1] * n          # the same, starting at i
    for i in range(n - 1, -1, -1):
        for j in range(i + 1, n):
            if r[i] <= r[j]:
                if b[j] + 1 > b[i]:
                    b[i], cb[i] = b[j] + 1, cb[j]
                elif b[j] + 1 == b[i]:
                    cb[i] += cb[j]
    L = max(f) if n else 0
    T = sum(cf[i] for i in range(n) if f[i] == L)
    out = []
    for i in range(n):
        if f[i] + b[i] - 1 < L:
            out.append("violation")
        elif cf[i] * cb[i] == T:
            out.append("ok")
        else:
            out.append("ambiguous")
    return out


def find_violations(S, unit, order):
    """Status of every coded row of the county against `order`. Returns a Series on S.index."""
    rank = {u: i for i, u in enumerate(order)}
    status = pd.Series("", index=S.index, dtype=object)
    use = S.in_county & ~S.bad_depth & pd.Series(unit, index=S.index).isin(order)
    for _, g in S[use].groupby("hole"):
        g = g.assign(rank=[rank[unit[i]] for i in g.index]).sort_values(["depth", "rank"], kind="stable")
        status[g.index] = violations(g["rank"])
    return status

"""M11: are two timeline-adjacent units one unit or two?"""
import numpy as np
import pandas as pd

from .within_hole import VIEWS


def prob(log_bf):
    return float(1.0 / (1.0 + np.exp(-np.clip(log_bf, -700, 700))))


def shared_params(model, ua, ub, na, nb):
    """Spreads common to both models, so the comparison is about the means only."""
    wa = na / (na + nb)
    P = {"AW": {}, "AL": {}, "pi": {}}
    for nd in model.nodes:
        (aw, al), (bw, bl) = model.spread(nd.name, ua), model.spread(nd.name, ub)
        P["AW"][nd.name] = float(np.exp(wa * np.log(aw) + (1 - wa) * np.log(bw)))
        P["AL"][nd.name] = float(np.exp(wa * np.log(al) + (1 - wa) * np.log(bl)))
    P["SW"] = wa * model.SW.get(ua, model.SW0) + (1 - wa) * model.SW.get(ub, model.SW0)
    P["SL"] = wa * model.SL.get(ua, model.SL0) + (1 - wa) * model.SL.get(ub, model.SL0)
    return P


def with_means(model, P, idxs, weights=None):
    Q = dict(P, pi={j: model.pool_pi(j, idxs, P["AW"][j], P["AL"][j], weights) for j in P["AW"]})
    Q["mu"] = model.pool_mu(idxs, P["SW"], P["SL"], weights)
    return Q


def group_ll(model, depth, idxs, weights, P):
    """Sum over boreholes of the joint log probability of that borehole's samples, per view."""
    idxs = np.asarray(idxs, int)
    weights = np.ones(len(idxs)) if weights is None else np.asarray(weights, float)
    tot = {}
    for h in dict.fromkeys(model.hole[idxs]):
        m = model.hole[idxs] == h
        o = np.argsort(depth[idxs[m]], kind="stable")
        terms = model.chain(idxs[m][o], P, weights[m][o])
        for j, t in terms.items():
            tot[j] = tot.get(j, 0.0) + t.sum()
    return tot


def dof(model):
    d = {nd.name: nd.k - 1 for nd in model.nodes}
    d["texture"] = 2
    d["overall"] = sum(model.weight[j] * d[j] for j in model.overall) + 2
    return d


def penalty(d, ha, hb):
    """Log of the unit-information factor (1 + n)^(d/2) for comparing two groups of ha and hb
    boreholes: n = ha*hb/(ha+hb), the number of boreholes the difference of means rests on."""
    return 0.5 * d * np.log1p(ha * hb / max(ha + hb, 1))


def one_vs_two(model, depth, ua, ub, ia, ib, wa=None, wb=None, n_holes=None):
    """Log Bayes factor (one unit : two units), overall and for each view alone."""
    ia, ib = np.asarray(ia, int), np.asarray(ib, int)
    wa = np.ones(len(ia)) if wa is None else wa
    wb = np.ones(len(ib)) if wb is None else wb
    P = shared_params(model, ua, ub, len(ia), len(ib))
    both = np.concatenate([ia, ib])
    wboth = np.concatenate([wa, wb])
    # a sample present on both sides (an uncoded one) is one sample of the single unit
    uniq, inv = np.unique(both, return_inverse=True)
    wuniq = np.bincount(inv, wboth)
    one = group_ll(model, depth, uniq, wuniq, with_means(model, P, uniq, wuniq))
    la = group_ll(model, depth, ia, wa, with_means(model, P, ia, wa))
    lb = group_ll(model, depth, ib, wb, with_means(model, P, ib, wb))
    two = {j: la.get(j, 0.0) + lb.get(j, 0.0) for j in one}
    ha, hb = n_holes or (len(set(model.hole[ia])), len(set(model.hole[ib])))
    d = dof(model)
    out = {}
    delta = model.combine(two) - model.combine(one)
    out["overall"] = -delta + penalty(d["overall"], ha, hb)
    for v in VIEWS:
        out[v] = -(two[v] - one[v]) + penalty(d[v], ha, hb)
    return out


def hellinger(p, q):
    return float(np.sqrt(0.5 * np.sum((np.sqrt(p) - np.sqrt(q)) ** 2)))


def separation_reference(model, leaf, min_grains=100):
    """Pooled broad-group composition of every (borehole, unit) group, and the distances
    between groups of the same unit in different boreholes (how far apart one unit gets)."""
    comps = {}
    for u, idx in model.members.items():
        for h in set(model.hole[idx]):
            x = leaf[idx[model.hole[idx] == h]].sum(axis=0)
            if x.sum() >= min_grains:
                comps[(u, h)] = x / x.sum()
    keys = list(comps)
    same = [hellinger(comps[a], comps[b]) for i, a in enumerate(keys) for b in keys[i + 1:] if a[0] == b[0]]
    return comps, np.array(same)


def merge_table(model, S, order, leaf, W, uncoded=None):
    """One row per pair of units adjacent in the timeline."""
    depth = S.depth.to_numpy()
    comps, same = separation_reference(model, leaf)
    rows = []
    for a, b in zip(order[:-1], order[1:]):
        if a not in model.units or b not in model.units:
            continue
        ia, ib = model.members_county[a], model.members_county[b]
        holes_a, holes_b = set(model.hole[ia]), set(model.hole[ib])
        lb = one_vs_two(model, depth, a, b, ia, ib)
        row = {"upper": a, "lower": b, "n_upper": len(ia), "n_lower": len(ib),
               "boreholes_upper": len(holes_a), "boreholes_lower": len(holes_b),
               "boreholes_with_both_stacked": len(holes_a & holes_b),
               "order_observed_in_boreholes": int(W.loc[a, b] + W.loc[b, a]) if a in W.index and b in W.index else 0,
               "pct_one_unit": 100 * prob(lb["overall"]),
               "one_vs_two_ratio": float(np.exp(np.clip(lb["overall"], -700, 700)))}
        for v in VIEWS:
            row[f"pct_one_unit_{v}"] = 100 * prob(lb[v])
        hs = [hellinger(comps[(a, x)], comps[(b, y)]) for x in holes_a for y in holes_b
              if (a, x) in comps and (b, y) in comps]
        row["separation_vs_same_code_pct"] = (100 * float(np.mean(same <= np.median(hs)))
                                             if hs and len(same) else np.nan)
        if uncoded is not None and len(uncoded):
            ua = uncoded[uncoded.unit.isin([a, b])].pivot_table(index="idx", columns="unit", values="posterior",
                                                                aggfunc="first").reindex(columns=[a, b]).fillna(0)
            cap = model.cfg["uncoded_cap"]
            tot = ua[a] + ua[b]
            ua = ua[tot > 0.01]
            tot = tot[tot > 0.01]
            wab = np.minimum(cap, tot)
            xa, xb = (wab * ua[a] / tot), (wab * ua[b] / tot)
            ea, eb = xa[xa > 0], xb[xb > 0]
            lu = one_vs_two(model, depth, a, b, np.concatenate([ia, ea.index.to_numpy(int)]),
                            np.concatenate([ib, eb.index.to_numpy(int)]),
                            np.concatenate([np.ones(len(ia)), ea.to_numpy()]),
                            np.concatenate([np.ones(len(ib)), eb.to_numpy()]),
                            n_holes=(len(holes_a), len(holes_b)))
            row["pct_one_unit_with_uncoded"] = 100 * prob(lu["overall"])
            row["uncoded_samples_used"] = int(len(ua))
            row["uncoded_weight_used"] = float(wab.sum())
        rows.append(row)
    return pd.DataFrame(rows)

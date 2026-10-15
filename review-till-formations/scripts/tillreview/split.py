"""M12: is one coded unit really two, one above the other in every borehole?"""
import numpy as np
import pandas as pd
from scipy.special import logsumexp

from .merge import prob
from .within_hole import VIEWS, split_ll


def _params(model, P, idxs, weights):
    Q = dict(P, pi={j: model.pool_pi(j, idxs, P["AW"][j], P["AL"][j], weights) for j in P["AW"]})
    Q["mu"] = model.pool_mu(idxs, P["SW"], P["SL"], weights)
    return Q


def split_unit(model, depth, u, idxs, weights=None, min_part=3, max_iter=20):
    """Best two-part vertical split of unit u and its log Bayes factor (two : one).

    idxs: samples of the unit (coded, plus any weighted uncoded ones). Each borehole gets a
    boundary position; a borehole may hold only the upper or only the lower part.
    "One unit" is the drift baseline: the same boundary positions, averaged the same way, with
    the borehole's local composition restarting at the boundary, but one typical composition for
    both parts. Only a difference in typical composition can favour two units.
    Returns None when no borehole has four or more coded samples to start from.
    """
    idxs = np.asarray(idxs, int)
    weights = np.ones(len(idxs)) if weights is None else np.asarray(weights, float)
    P = model.unit_params(u)
    holes = {}
    for h in dict.fromkeys(model.hole[idxs]):
        m = model.hole[idxs] == h
        o = np.argsort(depth[idxs[m]], kind="stable")
        holes[h] = (idxs[m][o], weights[m][o])
    one = {h: split_ll(model, ix, P, P, w) for h, (ix, w) in holes.items()}
    ll_one = {h: {j: L[j][len(holes[h][0])] for j in L} for h, L in one.items()}
    ll_drift = {h: {j: logsumexp(L[j]) - np.log(len(holes[h][0]) + 1) for j in L} for h, L in one.items()}

    starts = []
    for h, (ix, w) in holes.items():
        coded = np.flatnonzero(w >= 1.0)
        if len(coded) >= 4:
            gain = one[h]["overall"] - ll_one[h]["overall"]
            ks = [k for k in range(1, len(ix)) if (w[:k] >= 1).sum() >= 2 and (w[k:] >= 1).sum() >= 2]
            if ks:
                k = max(ks, key=lambda k: gain[k])
                starts.append((gain[k], h, k))
    if not starts:
        return None
    best = None
    for _, h0, k0 in sorted(starts, reverse=True)[:3]:
        cut = {h0: k0}
        up, low = holes[h0][0][:k0], holes[h0][0][k0:]
        wu, wl = holes[h0][1][:k0], holes[h0][1][k0:]
        L = {}
        for _ in range(max_iter):
            P1, P2 = _params(model, P, up, wu), _params(model, P, low, wl)
            L = {h: split_ll(model, ix, P1, P2, w) for h, (ix, w) in holes.items()}
            new = {h: int(L[h]["overall"].argmax()) for h in holes}
            up = np.concatenate([holes[h][0][:k] for h, k in new.items()])
            low = np.concatenate([holes[h][0][k:] for h, k in new.items()])
            wu = np.concatenate([holes[h][1][:k] for h, k in new.items()])
            wl = np.concatenate([holes[h][1][k:] for h, k in new.items()])
            if (wu >= 1).sum() < min_part or (wl >= 1).sum() < min_part:
                L = None
                break
            if new == cut:
                break
            cut = new
        if L is None:
            continue
        ll_two = sum(logsumexp(L[h]["overall"]) - np.log(len(holes[h][0]) + 1) for h in holes)
        if best is None or ll_two > best["ll_two"]:
            best = {"ll_two": ll_two, "cut": cut, "L": L, "up": up, "low": low, "wu": wu, "wl": wl,
                    "P1": P1, "P2": P2}
    if best is None:
        return {"found": False, "reason": f"no split leaves {min_part}+ coded samples on each side"}
    # Score: every borehole is predicted from typical compositions estimated from the OTHER boreholes
    # only (leave one borehole out), under "two units" (upper and lower compositions from the other
    # boreholes' parts) and under the drift baseline (one composition from the other boreholes, local
    # composition free to restart at the boundary). Boundary positions are averaged over, not picked.
    # Because no borehole helps estimate the compositions it is judged against, chance differences
    # found by the search for a split cannot inflate the result, and no penalty term is needed.
    H = len(holes)
    up_set, wu_all = best["up"], best["wu"]
    low_set, wl_all = best["low"], best["wl"]
    ind = [nd for nd in model.nodes if nd.view == "individual"]
    groups = {"overall": None, **{v: [v] for v in VIEWS}}
    if ind:
        groups["individual"] = [nd.name for nd in ind]
        nc = [nd.name for nd in ind if not nd.confusable]
        if nc:
            groups["individual_no_confusable"] = nc
    two = {g: 0.0 for g in groups}
    base = {g: 0.0 for g in groups}
    lat = fixed = 0.0
    for h, (ix, w) in holes.items():
        o_up, o_low = ~np.isin(up_set, ix), ~np.isin(low_set, ix)
        o_all = ~np.isin(idxs, ix)
        Ph = _params(model, P, idxs[o_all], weights[o_all])
        P1h = _params(model, P, up_set[o_up], wu_all[o_up]) if (wu_all[o_up] >= 1).any() else Ph
        P2h = _params(model, P, low_set[o_low], wl_all[o_low]) if (wl_all[o_low] >= 1).any() else Ph
        L2 = split_ll(model, ix, P1h, P2h, w)
        L1 = split_ll(model, ix, Ph, Ph, w)
        n1 = np.log(len(ix) + 1)
        for g, names in groups.items():
            a = L2["overall"] if names is None else sum(L2[j] for j in names)
            b = L1["overall"] if names is None else sum(L1[j] for j in names)
            two[g] += logsumexp(a) - n1
            base[g] += logsumexp(b) - n1
        lat += logsumexp(L2["overall"][[0, -1]]) - np.log(2)
        fixed += L1["overall"][len(ix)]
    out = {"found": True, "log_bf_two": {g: two[g] - base[g] for g in groups}}
    # lateral control: the same two groups with every borehole kept whole, against one unit kept whole;
    # vertical against lateral compares the evidence each gives for two typical compositions.
    out["log_bf_lateral"] = lat - fixed
    out["log_bf_vertical_vs_lateral"] = (two["overall"] - base["overall"]) - (lat - fixed)
    out.update(cut=best["cut"], holes=holes, upper=best["up"][best["wu"] >= 1], lower=best["low"][best["wl"] >= 1],
               P1=best["P1"], P2=best["P2"], n_boreholes=H)
    return out


def describe(model, P):
    """A unit's proportions in the reported categories, as percentages."""
    out = {}
    for j in ("age", "lithology", "colour"):
        for name, v in zip(model.node[j].parts, P["pi"][j]):
            out[f"{name}_pct"] = round(100 * float(v), 1)
    out["sand_pct"], out["clay_pct"] = round(float(P["mu"][0]), 1), round(float(P["mu"][1]), 1)
    out["silt_pct"] = round(100 - out["sand_pct"] - out["clay_pct"], 1)
    return out


def same_direction(model, res):
    """Boreholes holding both parts in which the upper part differs from the lower part in the same
    direction as the split overall, per view: (same, compared). Random drift scatters both ways."""
    out = {}
    for v in VIEWS:
        same = n = 0
        for h, k in res["cut"].items():
            ix, w = res["holes"][h]
            if not 0 < k < len(ix):
                continue
            up, low = ix[:k][w[:k] >= 1], ix[k:][w[k:] >= 1]
            if v == "texture":
                up, low = up[model.has_tex[up]], low[model.has_tex[low]]
                if not len(up) or not len(low):
                    continue
                d = model.tex[up].mean(axis=0) - model.tex[low].mean(axis=0)
                D = res["P1"]["mu"] - res["P2"]["mu"]
            else:
                cu, cl = model.cnt[v][up].sum(axis=0), model.cnt[v][low].sum(axis=0)
                if cu.sum() <= 0 or cl.sum() <= 0:
                    continue
                d = cu / cu.sum() - cl / cl.sum()
                D = res["P1"]["pi"][v] - res["P2"]["pi"][v]
            n += 1
            same += int(float(d @ D) > 0)
        out[v] = (same, n)
    return out


def split_table(model, S, uncoded=None, min_part=3):
    depth = S.depth.to_numpy()
    rows, bounds = [], []
    for u in model.units:
        idxs = model.members_county[u]
        row = {"unit": u, "n_samples": len(idxs), "n_boreholes": len(set(model.hole[idxs]))}
        res = split_unit(model, depth, u, idxs, min_part=min_part)
        if res is None or not res["found"]:
            row["result"] = "not testable: " + ("no borehole with 4+ samples" if res is None else res["reason"])
            rows.append(row)
            continue
        both = [h for h, k in res["cut"].items() if 0 < k < len(res["holes"][h][0])]
        row.update(result="tested", pct_two_units=100 * prob(res["log_bf_two"]["overall"]),
                   two_vs_one_ratio=float(np.exp(np.clip(res["log_bf_two"]["overall"], -700, 700))),
                   lateral_vs_one_ratio=float(np.exp(np.clip(res["log_bf_lateral"], -700, 700))),
                   vertical_vs_lateral_ratio=float(np.exp(np.clip(res["log_bf_vertical_vs_lateral"], -700, 700))),
                   n_upper=len(res["upper"]), n_lower=len(res["lower"]), boreholes_showing_both=len(both))
        for v in VIEWS + ["individual", "individual_no_confusable"]:
            if v in res["log_bf_two"]:
                row[f"pct_two_units_{v}"] = 100 * prob(res["log_bf_two"][v])
        for v, (same, n) in same_direction(model, res).items():
            row[f"same_direction_{v}"] = f"{same} of {n}" if n else ""
        for part in ("upper", "lower"):
            codes = S.code[res[part]].value_counts()
            row[f"codes_{part}"] = ", ".join(f"{c} {n}" for c, n in codes.items())
        for part, P in (("upper", res["P1"]), ("lower", res["P2"])):
            row.update({f"{part}_{k}": v for k, v in describe(model, P).items()})
        if uncoded is not None and len(uncoded):
            e = uncoded[(uncoded.unit == u) & (uncoded.weight > 0)]
            if len(e):
                ru = split_unit(model, depth, u, np.concatenate([idxs, e.idx.to_numpy(int)]),
                                np.concatenate([np.ones(len(idxs)), e.weight.to_numpy()]), min_part=min_part)
                if ru and ru["found"]:
                    row["pct_two_units_with_uncoded"] = 100 * prob(ru["log_bf_two"]["overall"])
                row["uncoded_samples_used"] = len(e)
        rows.append(row)
        for h, k in res["cut"].items():
            ix = res["holes"][h][0]
            bounds.append({"unit": u, "hole": h, "n_samples": len(ix), "n_upper": k, "n_lower": len(ix) - k,
                           "boundary_ft": (f"{depth[ix[k - 1]]:g}-{depth[ix[k]]:g}" if 0 < k < len(ix)
                                           else "all lower" if k == 0 else "all upper")})
    return pd.DataFrame(rows), pd.DataFrame(bounds)

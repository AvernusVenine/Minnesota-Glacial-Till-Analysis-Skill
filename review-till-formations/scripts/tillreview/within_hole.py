"""M8-M10: each sample against the units its borehole allows, contact positions, mixing."""
import numpy as np
import pandas as pd

from .model import dm_logpmf, mvn_logpdf

VIEWS = ["age", "lithology", "colour", "texture"]


def ratio(d):
    return float(np.exp(np.clip(d, -700, 700)))


def hole_neighbours(S, unit, order, usable):
    """Unit coded directly above and below each row, among usable coded rows of its borehole."""
    rank = {u: r for r, u in enumerate(order)}
    above, below = {}, {}
    coded = np.array([u in rank for u in unit])
    for _, g in S[S.in_county & ~S.bad_depth].groupby("hole"):
        g = g.sort_values("depth", kind="stable")
        idx = list(g.index)
        ref = [i for i in idx if usable[i] and coded[i]]
        for i in idx:
            d = S.depth[i]
            up = [unit[r] for r in ref if S.depth[r] < d and unit[r] != unit[i]]
            dn = [unit[r] for r in ref if S.depth[r] > d and unit[r] != unit[i]]
            above[i], below[i] = (up[-1] if up else None), (dn[0] if dn else None)
    return above, below


def candidates(u, up, dn, order):
    """Units a sample may be compared with: above, below and any the timeline puts between."""
    rank = {x: r for r, x in enumerate(order)}
    lo = rank[up] if up is not None else rank.get(u)
    hi = rank[dn] if dn is not None else rank.get(u)
    if lo is None or hi is None:
        lo = hi = lo if lo is not None else hi
    if lo is None:
        return []
    lo, hi = min(lo, hi), max(lo, hi)
    out = []
    for c in order[lo:hi + 1]:
        if c != u:
            out.append((c, "above" if c == up else "below" if c == dn else "between"))
    return out


def view_ratios(model, lp_alt, lp_own):
    """M8 breakdown: overall (weighted views, M1) and each view alone."""
    out = {"ratio_overall": ratio(model.total(lp_alt) - model.total(lp_own))}
    for v in VIEWS:
        out["ratio_" + v] = ratio(lp_alt[v] - lp_own[v])
    ind = [nd for nd in model.nodes if nd.view == "individual"]
    if ind:
        out["ratio_individual"] = ratio(sum(lp_alt[nd.name] - lp_own[nd.name] for nd in ind))
        out["ratio_individual_no_confusable"] = ratio(
            sum(lp_alt[nd.name] - lp_own[nd.name] for nd in ind if not nd.confusable))
    return out


def sample_table(model, S, unit, order, usable, typ, rng, R=999):
    """M8: one row per coded county till sample and allowed alternative unit."""
    above, below = hole_neighbours(S, unit, order, usable)
    rows = []
    for i in S.index[S.in_county & S.till & ~S.unusable]:
        u = unit[i]
        if u not in order:
            continue
        base = {"idx": i, "row": S.row[i], "hole": S.hole[i], "depth_ft": S.depth[i],
                "stratunit": S.code[i], "unit": u, "n_grains": int(S.n_grains[i]), "unit_above": above.get(i), "unit_below": below.get(i),
                "p_typical_of_coded": typ.p.get(i, np.nan)}
        if u not in model.units:
            rows.append({**base, "candidate": "", "relation": "coded unit has no clean samples to compare with"})
            continue
        pred_own = model.sample_params(i, u)
        lp_own = model.logp(i, pred_own)
        cands = candidates(u, above.get(i), below.get(i), order)
        if not cands:
            rows.append({**base, "candidate": "", "relation": "no other unit allowed"})
        for c, rel in cands:
            if c not in model.units:
                continue
            pred = model.sample_params(i, c)
            lp = model.logp(i, pred)
            r = view_ratios(model, lp, lp_own)
            rows.append({**base, "candidate": c, "relation": rel, **r,
                         "prob_candidate": r["ratio_overall"] / (1 + r["ratio_overall"]),
                         "p_typical_of_candidate": model.typicality(i, pred, rng, R)["p"]})
    return pd.DataFrame(rows)


def _view_pi(model, leaf_pi, j):
    """Proportions of view j implied by eight broad-group proportions."""
    from .variables import leaf_map
    v = leaf_pi @ leaf_map(model.node[j])
    return v / v.sum()


def mixing(model, i, predA, predB, step=0.01):
    """M10: best blend of units A and B for sample i, overall (weighted views) and per view.

    The blend is made of the broad groups, as when two tills physically mix, and each view's
    proportions are read from it. An individual-grain view blends within its group, weighted by
    how much of that group each unit brings. Concentrations blend on the log scale.
    """
    from .variables import LEAVES
    la = predA["alpha"]["leaf"] / predA["alpha"]["leaf"].sum()
    lb = predB["alpha"]["leaf"] / predB["alpha"]["leaf"].sum()
    grid = np.arange(0.0, 1.0 + step / 2, step)
    L = {j: np.zeros(len(grid)) for j in model.overall + ["texture"]}
    for t, lam in enumerate(grid):
        leaf = lam * la + (1 - lam) * lb
        for j in model.overall:
            if model.n[j][i] <= 0:
                continue
            nd = model.node[j]
            pa, pb = predA["alpha"][j], predB["alpha"][j]
            if nd.view == "individual":
                g = LEAVES.index(nd.within)
                wa, wb = lam * la[g], (1 - lam) * lb[g]
                pi = (wa * pa / pa.sum() + wb * pb / pb.sum()) / max(wa + wb, 1e-12)
            else:
                pi = _view_pi(model, leaf, j)
            A = np.exp(lam * np.log(pa.sum()) + (1 - lam) * np.log(pb.sum()))
            L[j][t] = dm_logpmf(model.cnt[j][i], A * np.maximum(pi, 1e-12))
        if model.has_tex[i]:
            L["texture"][t] = mvn_logpdf(model.tex[i], lam * predA["mean"] + (1 - lam) * predB["mean"],
                                         lam * predA["cov"] + (1 - lam) * predB["cov"])
    tot = model.combine(L)
    k = int(tot.argmax())
    out = {"mix_fraction_upper": float(grid[k]), "mix_vs_best_pure": ratio(tot[k] - max(tot[0], tot[-1])),
           "upper_vs_lower": ratio(tot[-1] - tot[0])}
    for v in VIEWS:
        out["mix_fraction_upper_" + v] = float(grid[int(L[v].argmax())])
        out["mix_vs_best_pure_" + v] = ratio(L[v].max() - max(L[v][0], L[v][-1]))
    return out


def split_ll(model, idxs, P_up, P_low, weights=None):
    """Log probability, per view, of every position k: idxs[:k] upper unit, idxs[k:] lower."""
    idxs = list(idxs)
    w = None if weights is None else list(weights)
    up = model.chain(idxs, P_up, w)
    low = model.chain(idxs[::-1], P_low, None if w is None else w[::-1])
    out = {}
    for j in up:
        a = np.concatenate([[0.0], np.cumsum(up[j])])
        b = np.concatenate([[0.0], np.cumsum(low[j])])[::-1]
        out[j] = a + b
    out["overall"] = model.combine(out)
    return out


def contact_table(model, S, unit, order, clean, rng):
    """M9 (+ M10 for the samples involved): every contact between two units stacked in a borehole.

    Returns the contacts, the mixing table, and one row per layer (sample) that the best-supported
    contact would move to the other unit, with that layer's own evidence for each unit.
    """
    rows, mixes, layers = [], [], []
    R = model.cfg["simulations"]
    for h, g in S[clean & S.in_county].groupby("hole"):
        g = g.sort_values("depth", kind="stable")
        idx = [i for i in g.index if unit[i] in order]
        runs = []
        for i in idx:
            if runs and unit[runs[-1][-1]] == unit[i]:
                runs[-1].append(i)
            else:
                runs.append([i])
        for A, B in zip(runs[:-1], runs[1:]):
            a, b = unit[A[0]], unit[B[0]]
            seq, k0, N = A + B, len(A), len(A) + len(B)
            local = {}
            params = []
            for u, own in ((a, A), (b, B)):
                elsewhere = len(model.members_county[u]) > len(own)
                local[u] = not elsewhere
                params.append(model.unit_params(u, exclude=tuple(own) if elsewhere else ()))
            L = split_ll(model, seq, params[0], params[1])
            ks = np.arange(1, N)
            ov = L["overall"][1:N]
            best = int(ks[ov.argmax()])
            post = np.exp(ov - ov.max())
            post /= post.sum()
            d = S.depth
            moved = seq[best:k0] if best < k0 else seq[k0:best]
            row = {"hole": h, "upper_unit": a, "lower_unit": b, "n_upper": len(A), "n_lower": len(B),
                   "geologist_contact_ft": f"{d[seq[k0 - 1]]:g}-{d[seq[k0]]:g}",
                   "best_contact_ft": f"{d[seq[best - 1]]:g}-{d[seq[best]]:g}",
                   "samples_that_change": len(moved),
                   "change_to": "" if best == k0 else (b if best < k0 else a),
                   "ratio_best_vs_geologist": ratio(L["overall"][best] - L["overall"][k0]),
                   "prob_geologist_position": float(post[k0 - 1]), "prob_best_position": float(post[best - 1])}
            for v in VIEWS:
                row["ratio_" + v] = ratio(L[v][best] - L[v][k0])
            row["mean_from_this_borehole_only"] = ", ".join(u for u, f in local.items() if f)
            row["changed_depths_ft"] = " ".join(f"{d[i]:g}" for i in moved)
            row["changed_rows"] = " ".join(str(S.row[i]) for i in moved)
            rows.append(row)
            for i in dict.fromkeys([seq[k0 - 1], seq[k0]] + moved):
                # each unit as logged in this borehole, without the sample itself
                pa = model.sample_params(i, a, others=np.array([o for o in A if o != i], int))
                pb = model.sample_params(i, b, others=np.array([o for o in B if o != i], int))
                mx = mixing(model, i, pa, pb)
                mixes.append({"row": S.row[i], "hole": h, "depth_ft": d[i],
                              "contact_interval_ft": row["geologist_contact_ft"], "stratunit": S.code[i],
                              "upper_unit": a, "lower_unit": b, **mx})
                if i not in moved:
                    continue
                coded, proposed = (a, b) if i in A else (b, a)
                p_c, p_p = (pa, pb) if i in A else (pb, pa)
                r = view_ratios(model, model.logp(i, p_p), model.logp(i, p_c))
                layers.append({"hole": h, "row": S.row[i], "depth_ft": d[i], "stratunit": S.code[i],
                               "coded_unit": coded, "proposed_unit": proposed, "n_grains": int(S.n_grains[i]),
                               **{("ratio_proposed_vs_coded" if k == "ratio_overall" else k): v
                                  for k, v in r.items()},
                               "prob_proposed": r["ratio_overall"] / (1 + r["ratio_overall"]),
                               "p_typical_of_coded": model.typicality(i, p_c, rng, R)["p"],
                               "p_typical_of_proposed": model.typicality(i, p_p, rng, R)["p"],
                               "mix_fraction_upper": mx["mix_fraction_upper"], "mix_vs_best_pure": mx["mix_vs_best_pure"],
                               "geologist_contact_ft": row["geologist_contact_ft"],
                               "best_contact_ft": row["best_contact_ft"],
                               "contact_ratio_best_vs_geologist": row["ratio_best_vs_geologist"],
                               "layers_moved_at_contact": len(moved)})
    return pd.DataFrame(rows), pd.DataFrame(mixes), pd.DataFrame(layers)

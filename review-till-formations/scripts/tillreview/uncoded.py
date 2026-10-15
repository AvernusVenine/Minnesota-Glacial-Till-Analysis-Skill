"""M14: how much an uncoded till sample may count toward each unit."""
import numpy as np
import pandas as pd


def uncoded_weights(model, S, unit, order, usable, target, rng, cap=0.5, pad_ft=15.0, p_min=0.01, R=499,
                    within_parent=False):
    """Weights for the samples in `target` (boolean on S: uncoded or member-unknown till samples).

    A unit is allowed for a sample only if it keeps the borehole in timeline order. Where the
    borehole gives no position at all (no coded row above or below, and the sample is not
    already coded to the units' formation: within_parent=False), the sample's elevation must
    also lie within the unit's county elevation range (+/- pad_ft).
    weight = min(cap, probability of the unit among the allowed units, equal priors), and 0
    if the sample is not typical (p < p_min) of its most probable unit.
    """
    rank = {u: r for r, u in enumerate(order)}
    units = [u for u in order if u in model.units]
    elev = S.elev.to_numpy()
    span = {u: (np.nanmin(elev[model.members_county[u]]) - pad_ft, np.nanmax(elev[model.members_county[u]]) + pad_ft)
            for u in units}
    coded = np.array([u in rank for u in unit]) & usable
    rows = []
    for h, g in S[S.in_county & ~S.unusable].groupby("hole"):
        ref = [i for i in g.index if coded[i]]
        for i in g.index[target[g.index]]:
            if S.n_grains[i] <= 0 or str(S.leached[i]).upper() == "Y":
                continue
            up = [rank[unit[r]] for r in ref if S.depth[r] < S.depth[i]]
            dn = [rank[unit[r]] for r in ref if S.depth[r] > S.depth[i]]
            lo, hi = (max(up) if up else 0), (min(dn) if dn else len(order) - 1)
            anchored = within_parent or bool(up or dn)
            allowed = [u for u in units if lo <= rank[u] <= hi
                       and (anchored or not np.isfinite(elev[i]) or span[u][0] <= elev[i] <= span[u][1])]
            if not allowed:
                continue
            preds = {u: model.sample_params(i, u) for u in allowed}
            lp = np.array([model.total(model.logp(i, preds[u])) for u in allowed])
            post = np.exp(lp - lp.max())
            post /= post.sum()
            best = allowed[int(post.argmax())]
            p = model.typicality(i, preds[best], rng, R)["p"]
            for u, q in zip(allowed, post):
                rows.append({"idx": i, "row": S.row[i], "hole": h, "depth_ft": S.depth[i],
                             "stratunit": S.code[i], "unit": u, "posterior": float(q),
                             "weight": float(min(cap, q)) if p >= p_min and q >= 0.01 else 0.0,
                             "most_likely_unit": best, "p_typical_of_most_likely": p,
                             "units_allowed": len(allowed)})
    return pd.DataFrame(rows, columns=["idx", "row", "hole", "depth_ft", "stratunit", "unit", "posterior",
                                       "weight", "most_likely_unit", "p_typical_of_most_likely", "units_allowed"])

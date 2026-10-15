"""The review steps, run once per loop (formations, or the members of one formation)."""
import numpy as np
import pandas as pd

from geological_priors import miscount, processes

from . import merge, split, timeline, uncoded, within_hole
from .model import Model


def typicality_table(model, unit, idxs, rng, R):
    rows = {}
    for i in idxs:
        t = model.typicality(i, model.sample_params(i, unit[i]), rng, R)
        rows[i] = {k: t[k] for k in ("p", "p_grains", "p_texture", "worst_branch", "branch_p")}
    return pd.DataFrame.from_dict(rows, orient="index",
                                  columns=["p", "p_grains", "p_texture", "worst_branch", "branch_p"])


def withheld_views(model, typ, idxs, p_cut):
    """M7: which measurements of each outlier are withheld.

    A view is atypical if its own p is below p_cut. If some measured views are atypical and at least
    one is not, only the atypical ones are withheld and the rest of the sample is still used. If no
    single view is atypical (the sample is odd only taken as a whole) or every measured view is,
    the whole sample is withheld. The hidden eight-group composition goes with the grain views.
    Returns {sample: list of withheld views, or None for the whole sample}.
    """
    out = {}
    for i in idxs:
        bp = typ.branch_p[i] or {}
        bad = [j for j, p in bp.items() if p < p_cut]
        out[i] = bad if bad and len(bad) < len(bp) else None
    return out


def explain_outliers(model, data, unit, idxs, rng, R):
    """Step 4: why each withheld sample may be atypical (priors sections 1 and 2)."""
    S = data.S
    ctx = processes.hole_context(S)
    rows = []
    for i in idxs:
        pred = model.sample_params(i, unit[i])
        typ = model.typicality(i, pred, rng, R)
        found = processes.check(model, S, data.leaf, i, pred, typ, rng, ctx)
        found += miscount.check(model, data.raw, data.leaf, i, pred, typ, rng)
        base = {"idx": i, "row": S.row[i], "hole": S.hole[i], "depth_ft": S.depth[i], "stratunit": S.code[i],
                "unit": unit[i], "n_grains": int(S.n_grains[i]), "p_typical": typ["p"],
                "p_grains": typ["p_grains"], "p_texture": typ["p_texture"], "least_typical_view": typ["worst_branch"]}
        if not found:
            found = [{"explanation": "unexplained", "status": "", "cite": "", "detail": ""}]
        rows += [{**base, **f} for f in found]
    return pd.DataFrame(rows)


def distribution_tables(model, rng, n=250, R=4000):
    """Step 5: each unit in the reported categories, with the range expected of a sample."""
    dist, spread = [], []
    for u in model.units:
        P = model.unit_params(u)
        idx = model.members_county[u]
        size = {"n_samples": len(idx), "n_boreholes": len(set(model.hole[idx]))}
        for j in ("age", "lithology", "colour"):
            AW, AL = P["AW"][j], P["AL"][j]
            rw, rl = 1 / (AW + 1), 1 / (AL + 1)
            AT = 1 / (rw + rl - rw * rl) - 1
            sims = {}
            for name, A in (("new_borehole", AT), ("same_borehole", AW)):
                x = rng.multinomial(n, rng.dirichlet(np.maximum(A * P["pi"][j], 1e-12), size=R))
                sims[name] = np.percentile(100 * x / n, [5, 95], axis=0)
            for k, part in enumerate(model.node[j].parts):
                dist.append({"unit": u, **size, "view": j, "category": part, "mean_pct": 100 * P["pi"][j][k],
                             "new_borehole_5": sims["new_borehole"][0, k], "new_borehole_95": sims["new_borehole"][1, k],
                             "same_borehole_5": sims["same_borehole"][0, k], "same_borehole_95": sims["same_borehole"][1, k]})
        tot, wit = P["SL"] + P["SW"], P["SW"]
        tex = {"sand": (P["mu"][0], tot[0, 0], wit[0, 0]), "clay": (P["mu"][1], tot[1, 1], wit[1, 1]),
               "silt": (100 - P["mu"].sum(), tot.sum(), wit.sum())}
        for part, (m, vt, vw) in tex.items():
            dist.append({"unit": u, **size, "view": "texture", "category": part, "mean_pct": m,
                         "new_borehole_5": m - 1.645 * vt ** 0.5, "new_borehole_95": m + 1.645 * vt ** 0.5,
                         "same_borehole_5": m - 1.645 * vw ** 0.5, "same_borehole_95": m + 1.645 * vw ** 0.5})
        for nd in model.nodes:
            if not nd.overall:
                continue
            AW, AL = P["AW"][nd.name], P["AL"][nd.name]
            spread.append({"unit": u, "view": nd.name, "weight": model.weight[nd.name],
                           "A_within_borehole": AW, "A_between_boreholes": AL,
                           "A_within_all_units": model.pooled[nd.name][0],
                           "A_between_all_units": model.pooled[nd.name][1]})
        spread.append({"unit": u, "view": "texture", "weight": model.weight["texture"],
                       "sd_sand_within": P["SW"][0, 0] ** 0.5, "sd_clay_within": P["SW"][1, 1] ** 0.5,
                       "sd_sand_between": P["SL"][0, 0] ** 0.5, "sd_clay_between": P["SL"][1, 1] ** 0.5})
    return pd.DataFrame(dist), pd.DataFrame(spread)


def run_loop(data, unit, nodes, cfg, rng, order=None, exclude=None, unknown=None, stop_after=None):
    """Steps 2-8 for one set of units.

    unit: label per sample (None = no unit at this level). order: confirmed timeline, or
    None to use the derived one. exclude: samples kept out of this loop. unknown: samples
    with no unit that may be used as evidence (M14); default every uncoded till sample.
    """
    S = data.S
    unit = np.asarray(unit, dtype=object)
    coded = np.array([u is not None for u in unit])
    exclude = np.zeros(len(S), bool) if exclude is None else np.asarray(exclude, bool)
    till_ok = (S.till & ~S.unusable).to_numpy() & ~exclude
    county = S.in_county.to_numpy()
    R = cfg["simulations"]
    out = {}

    # Step 2: timeline (M3)
    units = sorted(set(unit[county & till_ok & coded]))
    in_units = np.isin(unit, units) & coded
    rows = county & ~S.bad_depth.to_numpy() & in_units & ~exclude
    W, med = timeline.precedence(S.hole[rows], unit[rows], S.depth[rows])
    derived, n_viol = timeline.best_order(W, med.to_dict())
    confirmed = order is not None
    if confirmed:
        missing = [u for u in units if u not in order]
        if missing:
            raise ValueError(f"confirmed timeline lacks units present in the county: {missing}")
        order = [u for u in order if u in units]
        n_viol = int(sum(W.loc[b, a] for i, a in enumerate(order) for b in order[i + 1:]))
    else:
        order = derived
    out.update(order=order, derived_order=derived, timeline_confirmed=confirmed, order_violating_boreholes=n_viol,
               timeline=pd.DataFrame({"position": range(1, len(order) + 1), "unit": order,
                                      "median_depth_ft": [med.get(u, np.nan) for u in order],
                                      "boreholes": [S.hole[rows & (unit == u)].nunique() for u in order]}),
               timeline_pairs=timeline.pair_table(W, order), W=W)
    if stop_after == "timeline" or len(order) == 0:
        return out

    # Step 3: violations (M4)
    status = timeline.find_violations(S, np.where(exclude, None, unit), order)
    bad_order = status.isin(["violation", "ambiguous"]).to_numpy()
    v = S[bad_order]
    out["violations"] = pd.DataFrame({"idx": v.index, "row": v.row, "hole": v.hole, "depth_ft": v.depth, "stratunit": v.code,
                                      "unit": unit[bad_order], "deposittype": v.deposittype,
                                      "status": status[bad_order]})

    # Steps 4-5: robust first fit, withhold outliers (M7), refit (M5)
    eligible = till_ok & in_units & ~bad_order
    first = Model(data, unit, eligible, nodes, cfg)
    typ0 = typicality_table(first, unit, np.flatnonzero(eligible), rng, R)
    outlier = np.zeros(len(S), bool)
    outlier[typ0.index[typ0.p < cfg["outlier_p"]]] = True
    partial = withheld_views(first, typ0, np.flatnonzero(outlier), cfg["outlier_p"])
    whole = np.zeros(len(S), bool)
    whole[[i for i, v in partial.items() if v is None]] = True
    views = {j: np.zeros(len(S), bool) for j in [nd.name for nd in nodes] + ["texture"]}
    for i, v in partial.items():
        for j in v or []:
            views[j][i] = True
            if j in first.grain_views and "leaf" in views:
                views["leaf"][i] = True
    fit = eligible & ~whole
    model = Model(data, unit, fit, nodes, cfg, withheld_views=views)
    out["model"] = model
    out["weights"] = pd.DataFrame(model.weight_table())
    reviewed = np.flatnonzero(county & till_ok & in_units)
    typ = typicality_table(model, unit, reviewed, rng, R)
    out["outliers"] = explain_outliers(model, data, unit, np.flatnonzero(outlier & county), rng, R)
    if len(out["outliers"]):
        label = {i: ("all" if v is None else ", ".join(v)) for i, v in partial.items()}
        kept = {i: ("" if v is None else ", ".join(j for j in list(typ0.branch_p[i]) if j not in v))
                for i, v in partial.items()}
        o = out["outliers"]
        o.insert(o.columns.get_loc("least_typical_view") + 1, "withheld_views", o.idx.map(label))
        o.insert(o.columns.get_loc("withheld_views") + 1, "views_kept", o.idx.map(kept))
    out["distributions"], out["spread"] = distribution_tables(model, rng)
    part_label = {i: f"outlier ({', '.join(v)} withheld)" for i, v in partial.items() if v is not None}
    withheld = np.where(bad_order, "timeline violation", np.where(whole, "outlier", ""))
    withheld = np.array([part_label.get(i, w) for i, w in enumerate(withheld)], dtype=object)

    # Step 6: within boreholes (M8-M10)
    usable = (status == "ok").to_numpy()
    samples = within_hole.sample_table(model, S, unit, order, usable, typ, rng, R)
    if len(samples):
        samples.insert(7, "withheld", withheld[samples.idx.to_numpy()])
    out["samples"] = samples
    out["contacts"], out["mixing"], out["contact_layers"] = within_hole.contact_table(model, S, unit, order, fit, rng)

    # M14 evidence, then steps 7-8 (M11, M12)
    target = (S.till & ~S.unusable & S.in_county).to_numpy() & ~exclude & (~coded if unknown is None else unknown)
    unc = uncoded.uncoded_weights(model, S, unit, order, usable, target, rng, cfg["uncoded_cap"],
                                  cfg["uncoded_pad_ft"], cfg["outlier_p"], within_parent=unknown is not None)
    out["uncoded"] = unc
    out["merges"] = merge.merge_table(model, S, order, data.leaf, W, unc)
    out["splits"], out["split_boundaries"] = split.split_table(model, S, unc, cfg["split_min_samples"])
    out["counts"] = {"reviewed_samples": int(len(reviewed)), "clean_samples": int((fit & county).sum()),
                     "outliers": int((outlier & county).sum()),
                     "outliers_partly_used": int(sum(1 for i, v in partial.items() if v is not None and county[i])),
                     "timeline_violations": int(bad_order.sum()),
                     "buffer_samples_used_for_spread": int((fit & ~county).sum())}
    return out

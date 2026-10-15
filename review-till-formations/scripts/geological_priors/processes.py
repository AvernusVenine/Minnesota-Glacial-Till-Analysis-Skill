"""Checks from references/geological_priors.md section 1 for a sample that is atypical of its unit.

Each check follows the rule common to all processes:
  1. the altered parts moved in the expected direction,
  2. the composition with those parts removed is typical of the unit in that borehole,
  3. an independent indicator is present.
1 + 2 + 3 = supported, 1 + 2 = possible, 3 alone = indicator only.
"""
import re

import numpy as np
from scipy.stats import betabinom

from tillreview.model import dm_logpmf
from tillreview.variables import LEAVES

P_TYPICAL = 0.01
P_DIRECTION = 0.05     # one-sided: the altered groups are lower (or higher) than the unit predicts


def chroma(colour):
    """Munsell chroma from e.g. '2.5Y 5/6' (None if unreadable)."""
    m = re.search(r"/\s*(\d+(?:\.\d+)?)", str(colour))
    return float(m.group(1)) if m else None


def leaf_pi(model, pred):
    """Expected broad-group proportions of a sample under predictive parameters `pred`."""
    a = pred["alpha"]["leaf"]
    return a / a.sum()


def reduced_p(model, leaf, i, pred, rng, drop=(), merge=(), R=999):
    """Typicality of sample i's broad groups after dropping groups or pooling pairs of groups."""
    pi, x = leaf_pi(model, pred), leaf[i].astype(float)
    keep = [k for k, g in enumerate(LEAVES) if g not in drop]
    names = [LEAVES[k] for k in keep]
    pi, x = pi[keep], x[keep]
    for a, b in merge:
        if a in names and b in names:
            ka, kb = names.index(a), names.index(b)
            pi[ka] += pi[kb]; x[ka] += x[kb]
            pi, x = np.delete(pi, kb), np.delete(x, kb)
            names.pop(kb)
    if x.sum() <= 0:
        return np.nan
    alpha = pred["alpha"]["leaf"].sum() * pi / pi.sum()
    sims = rng.multinomial(int(x.sum()), rng.dirichlet(np.maximum(alpha, 1e-12), size=R))
    return float((1 + np.sum(dm_logpmf(sims, alpha) <= dm_logpmf(x, alpha) + 1e-9)) / (R + 1))


def moved(model, leaf, i, pred, groups, side):
    """Rule 1: one-sided probability that the groups' share is this low ('low') or this high
    ('high') under the unit's prediction for this borehole. True if below P_DIRECTION."""
    k = [LEAVES.index(g) for g in groups]
    n, x = int(leaf[i].sum()), int(leaf[i][k].sum())
    if n == 0:
        return False
    A, pi = pred["alpha"]["leaf"].sum(), float(leaf_pi(model, pred)[k].sum())
    a, b = A * pi, A * (1 - pi)
    p = betabinom.cdf(x, n, a, b) if side == "low" else betabinom.sf(x - 1, n, a, b)
    return bool(p < P_DIRECTION)


def _share(leaf, i, groups):
    n = leaf[i].sum()
    return leaf[i][[LEAVES.index(g) for g in groups]].sum() / n if n else np.nan


def _status(direction, typical_without, indicator):
    if direction and typical_without and indicator:
        return "supported"
    if direction and typical_without:
        return "possible"
    if indicator:
        return "indicator only"
    return None


def hole_context(S):
    """Non-till rows and counters of each borehole, looked up by the checks."""
    ctx = {}
    for h, g in S.groupby("hole"):
        ctx[h] = {"nontill": g.loc[~g.till, ["depth", "deposittype", "code"]],
                  "counters": g.loc[g.till, "counted_by"].dropna().astype(str)}
    return ctx


def check(model, S, leaf, i, pred, typ, rng, ctx):
    """Explanations for sample i being atypical of its unit. Returns a list of dicts."""
    out = []
    notes, depth = str(S.notes[i]), S.depth[i]
    expected = leaf_pi(model, pred)
    exp = lambda groups: expected[[LEAVES.index(g) for g in groups]].sum()
    near = ctx[S.hole[i]]["nontill"]
    grains_odd = typ["p_grains"] < P_TYPICAL if np.isfinite(typ["p_grains"]) else False
    texture_odd = typ["p_texture"] < P_TYPICAL if np.isfinite(typ["p_texture"]) else False

    def add(name, cite, status, detail):
        if status:
            out.append({"explanation": name, "status": status, "cite": cite, "detail": detail})

    def without(groups):
        p = reduced_p(model, leaf, i, pred, rng, drop=groups)
        return bool(np.isfinite(p) and p >= P_TYPICAL), p

    # section 1.C: not geological, rule these out first
    if S.low_count[i]:
        add("low grain count", "priors 1.C", "supported", f"{int(S.n_grains[i])} provenance grains")
    if S.nonprov_share[i] > 0.25:
        add("non-provenance grains", "priors 1.C", "supported",
            f"{100 * S.nonprov_share[i]:.0f}% unknown + secondary + other")
    counters = ctx[S.hole[i]]["counters"]
    if len(counters) > 1 and str(S.counted_by[i]) not in ("nan", counters.mode().iloc[0]):
        add("different counter from the rest of the borehole", "priors 1.C, 2", "indicator only",
            f"counted by {S.counted_by[i]}")
    top, bottom = S.depth_top[i], S.depth_bottom[i]
    if np.isfinite(top) and np.isfinite(bottom) and bottom > top:
        inside = near[(near.depth > top) & (near.depth < bottom)]
        if len(inside):
            add("sample interval spans a change of deposit", "priors 1.C", "indicator only",
                f"interval {top:g}-{bottom:g} ft")

    if grains_odd:
        # 1.1 leaching: carbonate lost
        carb = ["carbonate", "k_other"]
        ok, p = without(carb)
        ind = str(S.leached[i]).upper() == "Y" or bool(re.search(r"leach|noncalc|no carbonate", notes))
        add("leaching", "priors 1.1", _status(moved(model, leaf, i, pred, ["carbonate"], "low"), ok, ind),
            f"carbonate {100 * _share(leaf, i, ['carbonate']):.0f}% vs {100 * exp(['carbonate']):.0f}% expected; "
            f"p without carbonate = {p:.3f}")
        # 1.2 oxidation: shale, pyrite, lignite lost
        ok, p = without(["shale", "k_other"])
        c = chroma(S.wet_color[i])
        ind = (c is not None and c >= 3) or bool(re.search(r"oxid|limonite|goethite|iron oxide|stain", notes))
        add("oxidation", "priors 1.2", _status(moved(model, leaf, i, pred, ["shale"], "low"), ok, ind),
            f"shale {100 * _share(leaf, i, ['shale']):.0f}% vs {100 * exp(['shale']):.0f}% expected; "
            f"wet chroma {c}; p without shale = {p:.3f}")
        # 1.8 bedrock incorporation: Cretaceous or Paleozoic gained at the base
        below = near[(near.depth > depth) & (near.depth <= depth + 10)]
        ind = bool(((below.deposittype == "B") | below.code.fillna("").str.startswith(("K", "P", "O", "C", "D"))).any())
        for name, groups in (("Cretaceous", ["shale", "k_other"]), ("Paleozoic", ["carbonate", "pz_other"])):
            ok, p = without(groups)
            add(f"bedrock incorporation ({name})", "priors 1.8",
                # position defines this process: without bedrock just below it is not a candidate
                _status(True, ok, ind) if ind and moved(model, leaf, i, pred, groups, "high") else None,
                f"{name} {100 * _share(leaf, i, groups):.0f}% vs {100 * exp(groups):.0f}% expected; "
                f"p without it = {p:.3f}")
    # 1.4 secondary minerals
    sec = S.secondary[i] if np.isfinite(S.secondary[i]) else 0.0
    share = sec / (sec + S.n_grains[i]) if sec + S.n_grains[i] > 0 else 0.0
    noted = re.search(r"siderite|concretion|caliche|cement|gypsum|selenite|mangan|till ball", notes)
    if share > 0.10 or noted:
        add("secondary minerals", "priors 1.4", "indicator only",
            f"secondary = {100 * share:.0f}% of grains" + (f"; notes mention '{noted.group(0)}'" if noted else ""))
    # 1.5 soil formation
    if depth <= 5 or re.search(r"\bsoil\b|horizon|\bbt\b|root|organic|paleosol", notes):
        add("soil formation", "priors 1.5", "indicator only", f"depth {depth:g} ft")
    if texture_odd and not grains_odd:
        # 1.10 sorted sediment / 1.11 washing: texture shifted, grains typical
        sorted_near = near[(near.depth >= depth - 5) & (near.depth <= depth + 5) & near.deposittype.isin(["O", "L"])]
        kinds = ", ".join(sorted(set(sorted_near.deposittype)))
        ind = len(sorted_near) > 0 or bool(re.search(r"sand lens|stringer|laminat|interbed|sandy till", notes))
        add("incorporation of sorted sediment", "priors 1.10", _status(True, True, ind),
            f"texture atypical (p = {typ['p_texture']:.3f}), grains typical; sorted sediment within 5 ft: {kinds or 'none'}")
        if re.search(r"washed|\blag\b|colluv|slopewash|flow till|frost|cryo", notes):
            add("washing or surface reworking", "priors 1.11", "supported", "texture atypical, grains typical, noted")
        if re.search(r"loess|silt cap", notes):
            add("loess mixed into the top", "priors 1.6", "supported", "texture atypical, grains typical, noted")
    return out

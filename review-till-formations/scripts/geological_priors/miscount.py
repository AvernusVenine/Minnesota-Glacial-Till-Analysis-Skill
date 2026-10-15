"""Checks from references/geological_priors.md section 2: grains counters confuse.

A miscount is a possible explanation when pooling the two confused grain types (so that
it no longer matters which of the two a grain was called) makes the sample typical.
"""
import numpy as np

from tillreview.variables import CONFUSABLE, GROUPS

from .processes import P_TYPICAL, reduced_p


def group_of(grain):
    return next(g for g, cols in GROUPS.items() if grain in cols)


def check(model, raw, leaf, i, pred, typ, rng):
    """Returns a list of explanations for sample i (empty if its grains are typical)."""
    if not (np.isfinite(typ["p_grains"]) and typ["p_grains"] < P_TYPICAL):
        return []
    out = []
    for a, b in CONFUSABLE:
        if b is None:
            # confused with any grain: the sample must be typical once that grain's group is ignored
            if raw[a].iloc[i] <= 0:
                continue
            p = reduced_p(model, leaf, i, pred, rng, drop=[group_of(a)])
            detail = f"{int(raw[a].iloc[i])} {a} counted; p ignoring {group_of(a)} = {p:.3f}"
        else:
            ga, gb = group_of(a), group_of(b)
            if ga == gb:
                continue      # same broad group: the pair is already pooled at this level
            if raw[a].iloc[i] + raw[b].iloc[i] <= 0:
                continue
            p = reduced_p(model, leaf, i, pred, rng, merge=[(ga, gb)])
            detail = f"p with {ga} and {gb} pooled = {p:.3f}"
        if np.isfinite(p) and p >= P_TYPICAL:
            out.append({"explanation": f"possible miscount: {a} vs {b or 'another grain'}",
                        "status": "possible", "cite": "priors 2", "detail": detail})
    return out


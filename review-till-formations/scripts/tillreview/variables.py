"""M1, M2: the variables every statistic is computed on.

Grain counts are used as raw counts in four views, each its own composition (M1): age,
lithology and crystalline colour, plus (member loops) the individual grains of a group.
Texture is sand % and clay % as measured (M2). Nothing is log-transformed, so true zeros
need no replacement. A hidden eight-group composition (`leaf`) is kept only to simulate
samples (M1 weights, M7) and for the priors checks; it never enters a reported number.
"""
from dataclasses import dataclass, field

import numpy as np

# The eight non-overlapping groups the nine reported categories are summed from.
GROUPS = {
    "light": ["felsic", "quartzite", "clear_quartz"],
    "dark": ["mafic_igneous", "metasedvol"],
    "red": ["iron", "red_volcanic", "pc_arkosic", "pc_qtz_arenite"],
    "pc_other": ["pc_other"],
    "carbonate": ["carbonate"],
    "pz_other": ["paleo_ss", "paleo_shale", "paleo_other"],
    "shale": ["gray_shale", "speckled_shale"],
    "k_other": ["limestone", "inoceramus", "pyrite", "lignite", "ostrandersand", "cret_other"],
}
LEAVES = list(GROUPS)
GRAIN_COLS = [c for cols in GROUPS.values() for c in cols]
AGE = {"precambrian": ["light", "dark", "red", "pc_other"],
       "paleozoic": ["carbonate", "pz_other"],
       "cretaceous": ["shale", "k_other"]}
# Not provenance indicators: never part of a grain statistic.
NON_PROVENANCE = ["unknown", "secondary", "other_other"]
# QDI rollup columns the groups must reproduce (checked by io.check_rollups).
ROLLUPS = {"precambrian": AGE["precambrian"], "paleozoic": AGE["paleozoic"],
           "cret_total": AGE["cretaceous"], "shale": ["shale"]}
# Grain pairs counters confuse (geological_priors.md section 2). None = any other grain.
CONFUSABLE = [("felsic", "clear_quartz"), ("felsic", "carbonate"),
              ("mafic_igneous", "metasedvol"), ("inoceramus", None), ("red_volcanic", None)]
TEXTURE = ["sand_pct", "clay_pct"]


def _cols(groups):
    return [c for g in groups for c in GROUPS[g]]


@dataclass
class Node:
    """One view: a set of parts that split a total. `parts` maps part name -> grain columns.

    overall = True if the view enters the overall number (with its weight, M1); the hidden
    `leaf` composition does not. within = for an individual-grain view, the broad group whose
    grains it splits.
    """
    name: str
    parts: dict
    overall: bool = True
    view: str = ""
    within: str = ""
    confusable: bool = False
    note: str = field(default="")

    @property
    def k(self):
        return len(self.parts)


def formation_nodes():
    """The three grain views (M1) and the hidden eight-group composition used for simulation."""
    return [
        Node("age", {a: _cols(g) for a, g in AGE.items()}, view="age"),
        Node("lithology", {"carbonate": GROUPS["carbonate"],
                           "crystalline": _cols(AGE["precambrian"]),
                           "shale": GROUPS["shale"]}, view="lithology"),
        Node("colour", {g: GROUPS[g] for g in ("light", "dark", "red")}, view="colour"),
        Node("leaf", dict(GROUPS), overall=False, view="leaf"),
    ]


def leaf_map(node):
    """(8, k) matrix: which broad groups make up each part of a view (None if not built from groups)."""
    M = np.zeros((len(LEAVES), node.k))
    for gi, g in enumerate(LEAVES):
        for pk, cols in enumerate(node.parts.values()):
            if set(GROUPS[g]) <= set(cols):
                M[gi, pk] = 1.0
    return M if M.any() and node.view != "individual" else None


def is_confusable(grain):
    return any(grain in pair for pair in CONFUSABLE)


def member_node(group, grains):
    """M13: the selected individual grains of one group as a view of that group's grains."""
    parts = {g: [g] for g in grains}
    rest = [c for c in GROUPS[group] if c not in grains]
    if rest:
        parts["rest_of_" + group] = rest
    return Node("grains_" + group, parts, view="individual", within=group,
                confusable=any(is_confusable(g) for g in grains))


def node_counts(raw, node):
    """Counts of each part of a node. raw: DataFrame of grain columns. Returns (N, k)."""
    return np.column_stack([raw[cols].to_numpy(float).sum(axis=1) for cols in node.parts.values()])


def leaf_counts(raw):
    return np.column_stack([raw[cols].to_numpy(float).sum(axis=1) for cols in GROUPS.values()])

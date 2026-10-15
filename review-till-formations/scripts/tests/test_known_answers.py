"""Known-answer tests: each plants something in simulated data and checks it is recovered."""
import numpy as np
import pandas as pd
import pytest

import sim
from tillreview import DEFAULTS, glossary, report, timeline, variables
from tillreview.model import Model, dm_logpmf, dm_pred_alpha
from tillreview.pipeline import run_loop

CFG = {**DEFAULTS, "simulations": 199}
NODES = variables.formation_nodes


def model_of(data):
    unit = sim.units_of(data)
    return Model(data, unit, np.array([u is not None for u in unit]), NODES(), CFG), unit


def ratio(model, i, alt, own):
    lp_a = model.logp(i, model.sample_params(i, alt))
    lp_o = model.logp(i, model.sample_params(i, own))
    return {j: float(np.exp(lp_a[j] - lp_o[j])) for j in lp_a} | {
        "overall": float(np.exp(model.total(lp_a) - model.total(lp_o)))}


# ---------------------------------------------------------------- M1
def test_dm_probabilities_sum_to_one():
    alpha = np.array([2.0, 0.5, 0.01])
    xs = np.array([[a, b, 6 - a - b] for a in range(7) for b in range(7 - a)])
    assert np.exp(dm_logpmf(xs, alpha)).sum() == pytest.approx(1.0, abs=1e-9)


# ---------------------------------------------------------------- M3, M4
def test_timeline_recovers_order_and_finds_violation():
    hole = ["h1"] * 3 + ["h2"] * 3 + ["h3"] * 2 + ["h4"] * 3
    unit = ["A", "B", "C", "A", "B", "C", "B", "C", "B", "A", "C"]       # h4 has A below B
    depth = [1, 2, 3, 1, 2, 3, 1, 2, 1, 2, 3]
    W, med = timeline.precedence(hole, unit, depth)
    order, viol = timeline.best_order(W, med.to_dict())
    assert order == ["A", "B", "C"] and viol == 1
    assert timeline.violations([0, 1, 0, 1, 2]) == ["ok", "ambiguous", "ambiguous", "ok", "ok"]
    assert timeline.violations([0, 0, 2, 0, 0, 1]) == ["ok", "ok", "violation", "ok", "ok", "ok"]


# ---------------------------------------------------------------- M8
def test_zero_count_gives_finite_ratio_without_replacement():
    """A sample with no shale: strong evidence against the shaly unit, none between two shale-free units."""
    data = sim.make([[("A", "A", 8), ("B", "B", 8), ("C", "C", 8)] for _ in range(5)])
    model, unit = model_of(data)
    i = next(k for k in np.flatnonzero(unit == "B") if data.leaf[k][6] == 0 and data.leaf[k][7] == 0)
    r_shaly = ratio(model, i, "A", "B")
    assert np.isfinite(r_shaly["overall"]) and r_shaly["age"] < 1e-6
    # between two units that both lack shale, a zero shale count says nothing
    x = np.array([0.0, 120.0, 30.0])                     # shale, crystalline, carbonate
    a1, a2 = np.array([0.01, 60.0, 20.0]), np.array([0.01, 60.0, 20.0]) * 1.0
    a2[0] = 0.02
    assert 0.5 < np.exp(dm_logpmf(x, a2) - dm_logpmf(x, a1)) <= 1.0 + 1e-9


# ---------------------------------------------------------------- M1 weights
def test_weights_count_shared_grains_once_and_overall_is_the_weighted_product():
    data = sim.make([[("A", "A", 8), ("B", "B", 8), ("C", "C", 8)] for _ in range(6)], seed=29)
    model, unit = model_of(data)
    w = model.weight
    assert 0.35 < w["age"] < 0.7 and 0.35 < w["lithology"] < 0.7     # mostly the same grains
    assert w["colour"] > w["age"] and w["texture"] == 1.0
    i = int(np.flatnonzero(unit == "B")[0])
    r = ratio(model, i, "C", "B")
    by_hand = sum(w[v] * np.log(r[v]) for v in ("age", "lithology", "colour", "texture"))
    assert np.log(r["overall"]) == pytest.approx(by_hand, rel=1e-6, abs=1e-6)


def test_low_count_sample_is_less_decisive():
    data = sim.make([[("A", "A", 8), ("C", "C", 8)] for _ in range(5)], n_grains=lambda k: 60 if k == 3 else 300)
    model, unit = model_of(data)
    big = next(k for k in np.flatnonzero(unit == "A") if k != 3)
    assert ratio(model, 3, "C", "A")["age"] > ratio(model, big, "C", "A")["age"] * 100


def test_wide_view_contributes_little():
    """Unit W varies widely in crystalline colour between boreholes; unit T is tight. The same
    colour difference moves the ratio far less when the coded unit is the wide one."""
    wide = {"T": (sim.CRYSTAL, [60.0, 14.0]), "W": (sim.CRYSTAL, [60.0, 14.0]), "O": (sim.CARB, [33.0, 21.0])}
    rng = np.random.default_rng(3)
    data = sim.make([[("T", "T", 6), ("W", "W", 6), ("O", "O", 6)] for _ in range(8)], units=wide, AL=400.0)
    unit = sim.units_of(data)
    for h in set(data.S.hole):                      # give each borehole of W its own colour mix
        m = ((data.S.hole == h) & (unit == "W")).to_numpy()
        pc = data.leaf[m][:, :3].sum(axis=1)
        mix = rng.dirichlet([2.0, 2.0, 2.0])
        new = np.array([rng.multinomial(int(n), mix) for n in pc])
        data.leaf[m, :3] = new
        for k, g in enumerate(variables.LEAVES[:3]):
            data.raw.loc[m, variables.GROUPS[g][0]] = new[:, k]
    model = Model(data, unit, np.ones(len(unit), bool), NODES(), CFG)
    assert model.AL["colour"]["W"] < model.AL["colour"]["T"] / 3
    odd = np.array([120.0, 60.0, 40.0])             # a colour mix unlike the tight unit's
    typical = np.round(220 * model.unit_params("T")["pi"]["colour"])

    def drop(u):                                    # log probability lost by being odd, in a new borehole
        P = model.unit_params(u)
        a = dm_pred_alpha(P["pi"]["colour"], P["AW"]["colour"], P["AL"]["colour"], 0.0, 0.0)
        return dm_logpmf(typical, a) - dm_logpmf(odd, a)
    assert drop("W") < drop("T") / 2


# ---------------------------------------------------------------- M7
def test_typicality_is_calibrated_and_flags_a_planted_sample():
    data = sim.make([[("A", "A", 10), ("B", "B", 10)] for _ in range(8)], seed=5)
    data.leaf[4] = data.leaf[15]                    # a crystalline composition coded as shaly A
    for k, g in enumerate(variables.LEAVES):
        data.raw.loc[4, variables.GROUPS[g][0]] = data.leaf[4, k]
    model, unit = model_of(data)
    rng = np.random.default_rng(0)
    p = np.array([model.typicality(i, model.sample_params(i, unit[i]), rng, 199)["p"] for i in range(len(unit))])
    assert p[4] <= 0.01
    assert np.mean(np.delete(p, 4) < 0.05) < 0.15


def test_an_outlier_keeps_the_views_that_are_ordinary():
    """M7: a sample odd only in texture keeps its grains; a sample odd only in its grains keeps its texture."""
    data = sim.make([[("A", "A", 10), ("B", "B", 10)] for _ in range(8)], seed=31)
    tex_odd, grain_odd = 3, 25                      # both coded A
    data.tex[tex_odd] = [90.0, 2.0]                 # far from A's 42 % sand, 20 % clay
    data.leaf[grain_odd] = data.leaf[15]            # a crystalline B grain count coded as shaly A
    for k, g in enumerate(variables.LEAVES):
        data.raw.loc[grain_odd, variables.GROUPS[g][0]] = data.leaf[grain_odd, k]
    res = run_loop(data, sim.units_of(data), NODES(), CFG, np.random.default_rng(0), order=["A", "B"])
    o = res["outliers"].drop_duplicates("row").set_index("row")
    m = res["model"]
    r_tex, r_grain = data.S.row[tex_odd], data.S.row[grain_odd]
    assert o.loc[r_tex, "withheld_views"] == "texture"
    assert m.fit[tex_odd] and m.excl["texture"][tex_odd] and not m.excl["age"][tex_odd]
    assert "texture" not in o.loc[r_grain, "withheld_views"] and "texture" in o.loc[r_grain, "views_kept"]
    assert m.fit[grain_odd] and m.excl["age"][grain_odd] and m.excl["leaf"][grain_odd]
    assert not m.excl["texture"][grain_odd]


# ---------------------------------------------------------------- M9, pipeline
def test_miscoded_contact_samples_are_found():
    """Two samples at the top of B are coded A in one borehole: the contact must move up by two."""
    holes = [[("A", "A", 8), ("B", "B", 8)] for _ in range(5)]
    holes[2] = [("A", "A", 8), ("B", "A", 2), ("B", "B", 6)]
    data = sim.make(holes, seed=7)
    res = run_loop(data, sim.units_of(data), NODES(), CFG, np.random.default_rng(0), order=["A", "B"])
    planted = set(data.S.row[(data.S.hole == "H02") & (data.S.true == "B") & (data.S.code == "A")])
    c = res["contacts"].set_index("hole").loc["H02"]
    sm = res["samples"]
    assert set(sm.row[(sm.ratio_overall > 100)]) == planted
    # outliers are withheld before contacts are tested, so the contact is tested without them
    assert c.samples_that_change == 0 or set(map(int, c.changed_rows.split())) <= planted
    assert set(res["outliers"].row) >= planted


def test_contact_moves_when_miscoded_samples_are_not_extreme():
    """A and a similar unit A2: the miscoded samples are not outliers, so the contact test finds them."""
    near = {"A": (sim.SHALY, [42.0, 20.0]),
            "A2": ([0.44, 0.06, 0.01, 0.00, 0.22, 0.01, 0.22, 0.04], [46.0, 18.0])}
    holes = [[("A", "A", 10), ("A2", "A2", 10)] for _ in range(5)]
    holes[1] = [("A", "A", 10), ("A2", "A", 3), ("A2", "A2", 7)]
    data = sim.make(holes, units=near, seed=11)
    res = run_loop(data, sim.units_of(data), NODES(), {**CFG, "outlier_p": 0.0}, np.random.default_rng(0),
                   order=["A", "A2"])
    c = res["contacts"].set_index("hole")
    assert c.loc["H01", "samples_that_change"] == 3 and c.loc["H01", "change_to"] == "A2"
    assert c.loc["H01", "ratio_best_vs_geologist"] > 1e6
    # one row per moved layer, each with its own evidence for the proposed unit
    ly = res["contact_layers"]
    ly = ly[ly.hole == "H01"]
    planted = set(data.S.row[(data.S.hole == "H01") & (data.S.true == "A2") & (data.S.code == "A")])
    assert set(ly.row) == planted and (ly.coded_unit == "A").all() and (ly.proposed_unit == "A2").all()
    assert (ly.ratio_proposed_vs_coded > 1).all() and ly.p_typical_of_proposed.notna().all()


def test_correct_contacts_rarely_move_across_simulated_counties():
    """False alarms: in 10 counties with every contact correctly coded, a strong suggestion (> 20)
    to move a contact appears in at most one county."""
    near = {"A": (sim.SHALY, [42.0, 20.0]),
            "A2": ([0.44, 0.06, 0.01, 0.00, 0.22, 0.01, 0.22, 0.04], [46.0, 18.0])}
    worst = []
    for seed in range(200, 210):
        data = sim.make([[("A", "A", 10), ("A2", "A2", 10)] for _ in range(5)], units=near, seed=seed)
        res = run_loop(data, sim.units_of(data), NODES(), {**CFG, "outlier_p": 0.0}, np.random.default_rng(0),
                       order=["A", "A2"])
        worst.append(res["contacts"].ratio_best_vs_geologist.max())
    assert np.median(worst) < 5 and np.sum(np.array(worst) > 20) <= 1


# ---------------------------------------------------------------- M11
def test_one_unit_split_by_borehole_merges_and_distinct_units_do_not():
    holes = [[("A", "A1" if h % 2 else "A2", 6), ("B", "B", 6)] for h in range(10)]
    data = sim.make(holes, seed=13)
    res = run_loop(data, sim.units_of(data), NODES(), CFG, np.random.default_rng(0), order=["A1", "A2", "B"])
    m = res["merges"].set_index(["upper", "lower"])
    assert m.loc[("A1", "A2"), "pct_one_unit"] > 90
    assert m.loc[("A2", "B"), "pct_one_unit"] < 1


def test_one_unit_with_two_codes_mostly_comes_out_as_one_unit_across_counties():
    """False alarms: a single unit coded A1 or A2 by borehole, in 10 simulated counties. The overall is the
    weighted product of the views, so chance leanings in several views can combine: most counties, not all,
    come out as one unit (validation record in statistical_methods.md)."""
    pct = []
    for seed in range(300, 310):
        holes = [[("A", "A1" if h % 2 else "A2", 6), ("B", "B", 6)] for h in range(10)]
        data = sim.make(holes, seed=seed)
        res = run_loop(data, sim.units_of(data), NODES(), CFG, np.random.default_rng(0), order=["A1", "A2", "B"])
        pct.append(res["merges"].set_index(["upper", "lower"]).loc[("A1", "A2"), "pct_one_unit"])
    assert np.median(pct) > 60 and np.sum(np.array(pct) < 10) <= 2


# ---------------------------------------------------------------- M12
def test_two_units_under_one_code_split_and_a_single_unit_does_not():
    holes = [[("A", "X", 6), ("C", "X", 6), ("B", "B", 6)] for _ in range(6)]
    data = sim.make(holes, seed=17)
    res = run_loop(data, sim.units_of(data), NODES(), {**CFG, "outlier_p": 0.0}, np.random.default_rng(0),
                   order=["X", "B"])
    s = res["splits"].set_index("unit")
    assert s.loc["X", "pct_two_units"] > 99 and s.loc["X", "boreholes_showing_both"] == 6
    b = res["split_boundaries"]
    assert (b[b.unit == "X"].n_upper == 6).all()
    assert s.loc["X", "vertical_vs_lateral_ratio"] > 100


def test_a_single_unit_rarely_splits_across_simulated_counties():
    """False alarms: a homogeneous unit in 10 counties. Judged out of sample (each borehole against the
    other boreholes), a unit with nothing to split reads 'cannot decide, leans one' rather than 'certainly
    one'; any high value rests on at most two boreholes that hold both parts, which the report shows."""
    pct, both = [], []
    for seed in range(100, 110):
        holes = [[("A", "X", 6), ("C", "X", 6), ("B", "B", 6)] for _ in range(6)]
        data = sim.make(holes, seed=seed)
        res = run_loop(data, sim.units_of(data), NODES(), {**CFG, "outlier_p": 0.0}, np.random.default_rng(0),
                       order=["X", "B"])
        r = res["splits"].set_index("unit").loc["B"]
        pct.append(r.pct_two_units)
        both.append(r.boreholes_showing_both)
    pct, both = np.array(pct), np.array(both)
    assert np.median(pct) < 30
    assert (both[pct > 50] <= 2).all()


def _drift_county(seed):
    """One unit X whose composition drifts down every core in a random direction: each borehole's top and
    bottom are the same base composition nudged independently. No consistent second composition exists."""
    rng = np.random.default_rng(seed)
    units, holes = {"B": sim.UNITS["B"]}, []
    for h in range(8):
        for part in ("t", "b"):
            units[f"{part}{h}"] = (list(rng.dirichlet(60 * np.array(sim.SHALY) + 0.5)),
                                   list(np.array([42.0, 20.0]) + rng.normal(0, 3, 2)))
        holes.append([(f"t{h}", "X", 5), (f"b{h}", "X", 5), ("B", "B", 6)])
    return sim.make(holes, units=units, seed=seed)


def test_drift_down_cores_is_not_a_split():
    """Composition that changes down each core in random directions must not come out as two units."""
    pct = []
    for seed in range(700, 710):
        data = _drift_county(seed)
        res = run_loop(data, sim.units_of(data), NODES(), {**CFG, "outlier_p": 0.0}, np.random.default_rng(0),
                       order=["X", "B"])
        pct.append(res["splits"].set_index("unit").loc["X", "pct_two_units"])
    assert np.median(pct) < 5 and np.sum(np.array(pct) > 50) <= 1


def test_a_small_but_consistent_split_is_still_found():
    """Power: two similar units (A over A2 in every borehole, same code) still come out as two."""
    near = {"A": (sim.SHALY, [42.0, 20.0]), "B": (sim.CRYSTAL, [60.0, 14.0]),
            "A2": ([0.44, 0.06, 0.01, 0.00, 0.22, 0.01, 0.22, 0.04], [46.0, 18.0])}
    pct = []
    for seed in range(500, 506):
        data = sim.make([[("A", "X", 5), ("A2", "X", 5), ("B", "B", 6)] for _ in range(8)], units=near, seed=seed)
        res = run_loop(data, sim.units_of(data), NODES(), {**CFG, "outlier_p": 0.0}, np.random.default_rng(0),
                       order=["X", "B"])
        pct.append(res["splits"].set_index("unit").loc["X", "pct_two_units"])
    assert np.median(pct) > 90 and min(pct) > 50


def test_a_difference_between_boreholes_is_flagged_as_lateral_not_vertical():
    holes = [[("A" if h % 2 else "C", "X", 1 if h else 4), ("B", "B", 6)] for h in range(16)]
    data = sim.make(holes, seed=23)
    res = run_loop(data, sim.units_of(data), NODES(), {**CFG, "outlier_p": 0.0}, np.random.default_rng(0),
                   order=["X", "B"])
    s = res["splits"].set_index("unit")
    assert s.loc["X", "vertical_vs_lateral_ratio"] < 1


# ---------------------------------------------------------------- column guide
def test_every_output_column_has_a_plain_language_explanation():
    data = sim.make([[("A", "A", 6), ("C", "C", 6), ("B", "B", 6)] for _ in range(5)], seed=19)
    res = run_loop(data, sim.units_of(data), NODES(), CFG, np.random.default_rng(0))
    missing = [(key, c) for key, _ in report.TABLES if res.get(key) is not None
               for c in res[key].columns if c != "idx" and not glossary.explain(key, c)]
    assert not missing
    member_only = ["ratio_individual", "ratio_individual_no_confusable", "G", "p", "q", "selected", "pct_in_QNUH"]
    assert all(glossary.explain("samples", c) or glossary.explain("grain_selection", c) for c in member_only)
    assert all(glossary.explain("splits", c) for c in
               ["pct_two_units_individual", "pct_two_units_individual_no_confusable"])

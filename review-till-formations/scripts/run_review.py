"""Review the till codes of one county.

    python run_review.py config.json --stop-after timeline   # derive the timelines, confirm them
    python run_review.py config.json                         # full review

config.json keys:
    data      path to the QDI export (CSV)
    county    cty_name to review
    out       output folder
    codes     path to stratunit_codes.csv (default: ../references/stratunit_codes.csv)
    timeline  {"formation": [top ... bottom], "members": {"QNUU": [top ... bottom], ...}}
              confirmed orders; a loop without one uses the derived order and the report says so
    plus any key of tillreview.DEFAULTS to override it
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from tillreview import DEFAULTS, io, members, report, variables  # noqa: E402
from tillreview.pipeline import run_loop  # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("config")
    ap.add_argument("--stop-after", choices=["timeline"])
    args = ap.parse_args(argv)
    with open(args.config, encoding="utf-8") as f:
        user = json.load(f)
    cfg = {**DEFAULTS, **{k: v for k, v in user.items() if k in DEFAULTS}}
    out = user["out"]
    os.makedirs(out, exist_ok=True)
    tl = user.get("timeline", {})
    codes = io.load_codes(user.get("codes", os.path.join(HERE, "..", "references", "stratunit_codes.csv")))
    df = io.load_csv(user["data"])
    data = io.prepare(df, user["county"], codes, cfg["buffer_km"], cfg["low_count"])
    S = data.S
    rng = np.random.default_rng(cfg["seed"])

    bad = io.bad_sample_table(data)
    bad.to_csv(os.path.join(out, "01_bad_samples.csv"), index=False)
    sections, ledger_loops, results = [], {}, {}

    # Formation loop
    formation = np.array([f if isinstance(f, str) and f else None for f in S.formation], dtype=object)
    res = run_loop(data, formation, variables.formation_nodes(), cfg, rng, order=tl.get("formation"),
                   stop_after=args.stop_after)
    report.write_loop(res, os.path.join(out, "formation"))
    results["formation"] = res
    sections.append(report.loop_section("Formation loop", res, "formation", cfg))
    ledger_loops["formation"] = {k: res.get(k) for k in ("order", "derived_order", "timeline_confirmed", "counts")}
    if "model" in res:
        ledger_loops["formation"]["view_weights"] = res["model"].weight

    # Member loops: original codes, minus samples flagged in the formation loop
    flagged = np.zeros(len(S), bool)
    if "samples" in res and len(res["samples"]) and "ratio_overall" in res["samples"]:
        sm = res["samples"]
        flagged[sm.idx[sm.ratio_overall >= cfg["member_exclusion_ratio"]].to_numpy(int)] = True
        flagged[res["violations"].idx.to_numpy(int)] = True
    is_member = (S.level == "member").to_numpy()
    for F in res.get("order", []):
        in_F = (formation == F)
        unit = np.where(in_F & is_member, S.code.to_numpy(dtype=object), None)
        county_till = (S.in_county & S.till & ~S.unusable).to_numpy()
        if len(set(unit[county_till & in_F & is_member & ~flagged])) < 2:
            continue
        nodes, table = variables.formation_nodes(), None
        if args.stop_after != "timeline":
            fit_idx = np.flatnonzero((S.till & ~S.unusable).to_numpy() & in_F & is_member & ~flagged)
            extra, table = members.select_grains(data.raw, unit, fit_idx, cfg["member_grain_q"])
            nodes += extra
        rm = run_loop(data, unit, nodes, cfg, rng, order=tl.get("members", {}).get(F), exclude=flagged,
                      unknown=in_F & ~is_member, stop_after=args.stop_after)
        rm["grain_selection"] = table
        report.write_loop(rm, os.path.join(out, f"members_{F}"))
        results[f"members_{F}"] = rm
        sections.append(report.loop_section(f"Member loop: {F}", rm, f"members_{F}", cfg))
        ledger_loops[f"members_{F}"] = {k: rm.get(k) for k in ("order", "derived_order", "timeline_confirmed", "counts")}
        if "model" in rm:
            ledger_loops[f"members_{F}"]["view_weights"] = rm["model"].weight
        ledger_loops[f"members_{F}"]["excluded_from_formation_loop"] = int((flagged & in_F & S.in_county.to_numpy()).sum())

    cty = S[S.in_county & S.till]
    meta = {"county": user["county"], "data": user["data"], "seed": cfg["seed"], "buffer_km": cfg["buffer_km"],
            "config_path": args.config, "cfg": cfg, "timeline_config": tl,
            "county_till_samples": int(len(cty)), "county_coded_till": int(cty.code.notna().sum()),
            "area_holes": int(S.hole.nunique()), "rollups": io.check_rollups(df[df.cty_name == user["county"]]),
            "codes_not_in_lookup": data.codes_not_in_lookup,
            "bad_summary": bad.groupby(["problem", "effect"]).size().reset_index(name="samples"),
            "versions": {"python": sys.version.split()[0], "numpy": np.__version__, "pandas": pd.__version__},
            "stopped_after": args.stop_after}
    report.write_report(os.path.join(out, "report.md"), meta, sections)
    report.write_guide(os.path.join(out, "column_guide.md"), bad, results)
    report.write_ledger(os.path.join(out, "ledger.json"), {**{k: v for k, v in meta.items() if k != "bad_summary"},
                                                           "loops": ledger_loops})
    print(f"wrote {out}/report.md")
    for name, lp in ledger_loops.items():
        state = "confirmed" if lp["timeline_confirmed"] else "DERIVED, not confirmed"
        print(f"{name}: {' > '.join(lp['order'] or [])}  [{state}]")


if __name__ == "__main__":
    main()

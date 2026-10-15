"""CSV tables, the Markdown report and the ledger."""
import json
import os

import numpy as np
import pandas as pd

from . import glossary

TABLES = [("timeline", "02_timeline"), ("timeline_pairs", "02_timeline_pairs"), ("violations", "03_violations"),
          ("outliers", "04_outliers"), ("distributions", "05_distributions"), ("spread", "05_spread"),
          ("samples", "06_sample_likelihoods"), ("contacts", "06_contacts"),
          ("contact_layers", "06_contact_layers"), ("mixing", "06_mixing"),
          ("merges", "07_merge_pairs"), ("splits", "08_split_candidates"),
          ("split_boundaries", "08_split_boundaries"), ("uncoded", "09_uncoded_evidence"),
          ("grain_selection", "10_member_grain_selection"), ("weights", "11_view_weights")]
VIEW_COLS = ["ratio_overall", "ratio_age", "ratio_lithology", "ratio_colour", "ratio_texture"]
# column order of the unit distributions table: texture, lithology, age, colour
DIST_ORDER = ["sand", "silt", "clay", "crystalline", "carbonate", "shale",
              "precambrian", "paleozoic", "cretaceous", "light", "dark", "red"]


def write_loop(res, folder):
    os.makedirs(folder, exist_ok=True)
    for key, name in TABLES:
        if key in res and res[key] is not None:
            res[key].drop(columns=["idx"], errors="ignore").to_csv(os.path.join(folder, name + ".csv"), index=False)


def _fmt(v):
    if v is None or (isinstance(v, float) and not np.isfinite(v)):
        return ""
    if isinstance(v, (bool, np.bool_)):
        return "yes" if v else "no"
    if isinstance(v, (float, np.floating)):
        return f"{v:.3g}"
    return str(v)


def _pct(v):
    if v is None or not np.isfinite(v):
        return ""
    return "<0.01" if v < 0.01 else ">99.99" if v > 99.99 else f"{v:.3g}"


def md_table(df, cols=None, n=20, rename=None):
    if df is None or len(df) == 0:
        return "_None._\n"
    cols = [c for c in (cols or df.columns) if c in df.columns]
    head = [(rename or {}).get(c, c) for c in cols]
    lines = ["| " + " | ".join(head) + " |", "|" + "---|" * len(cols)]
    for _, r in df[cols].head(n).iterrows():
        lines.append("| " + " | ".join(_pct(v) if c.startswith("pct_") else _fmt(v) for c, v in zip(cols, r)) + " |")
    if len(df) > n:
        lines.append(f"\n_{len(df) - n} more rows in the CSV._")
    return "\n".join(lines) + "\n"


def table(key, df, cols=None, n=20):
    """A table followed by the plain-language meaning of the columns it shows."""
    if df is None or len(df) == 0:
        return md_table(df)
    shown = [c for c in (cols or df.columns) if c in df.columns and c != "idx"]
    return md_table(df, shown, n) + "\n" + glossary.column_notes(key, shown)


def loop_section(title, res, folder, cfg):
    s = [f"## {title}\n", f"Tables: `{folder}/`\n"]
    if not res.get("order"):
        return "\n".join(s + ["No coded units to review.\n"])
    conf = "confirmed by the user" if res["timeline_confirmed"] else "**derived, NOT yet confirmed by the user**"
    s.append(f"### Timeline (M3) — {conf}\n")
    s.append("Top to bottom: " + " → ".join(res["order"]) + f". Borehole-pair disagreements with this order: "
             f"{res['order_violating_boreholes']}.\n")
    if not res["timeline_confirmed"] and res["derived_order"] != res["order"]:
        s.append("Derived order: " + " → ".join(res["derived_order"]) + "\n")
    pairs = res["timeline_pairs"]
    s.append("Pairs with a borehole that disagrees, and adjacent pairs never seen in one borehole "
             "(their order is not observed):\n")
    s.append(table("timeline_pairs",
                   pairs[(pairs.boreholes_disagree > 0) | (pairs.adjacent & ~pairs.order_observed)], n=30))
    if "samples" not in res:
        return "\n".join(s)
    c = res["counts"]
    part = c.get("outliers_partly_used", 0)
    s.append(f"Reviewed {c['reviewed_samples']} coded till samples; {c['clean_samples']} used in the unit "
             f"descriptions, {c['outliers']} outliers ({part} of them still used for their ordinary views, "
             f"{c['outliers'] - part} withheld entirely) and {c['timeline_violations']} timeline violations "
             f"withheld; {c['buffer_samples_used_for_spread']} buffer samples lent spread (M5).\n")
    s.append("### Timeline violations (M4)\n")
    s.append(table("violations", res["violations"]))
    s.append(f"### Outliers (M7, p < {cfg['outlier_p']}; explanations from priors sections 1-2)\n")
    s.append("Only the views in which an outlier is itself atypical are withheld; its other views still describe "
             "its unit. A sample odd only as a whole, or odd in every view, is withheld entirely.\n")
    s.append(table("outliers", res["outliers"], ["row", "hole", "depth_ft", "stratunit", "p_typical",
                                        "least_typical_view", "withheld_views", "views_kept", "explanation",
                                        "status", "cite", "detail"], n=30))
    s.append("### Unit distributions (M5), mean % by category\n")
    d = res["distributions"]
    if len(d):
        piv = d.pivot_table(index="unit", columns="category", values="mean_pct", sort=False).round(1)
        piv = piv[[c for c in DIST_ORDER if c in piv.columns] + [c for c in piv.columns if c not in DIST_ORDER]]
        s.append(md_table(piv.reindex(res["order"]).dropna(how="all").reset_index(), n=50))
        s.append("Each value is the unit's typical percentage of that category. Sand + silt + clay add to 100, "
                 "as do crystalline + carbonate + shale, Precambrian + Paleozoic + Cretaceous, and "
                 "light + dark + red.\n")
        s.append("Ranges for a sample in a new borehole and within one borehole are in `05_distributions.csv`; "
                 "spreads per view in `05_spread.csv`.\n")
    if res.get("weights") is not None:
        s.append("### Weights of the views in the overall number (M1)\n")
        s.append("The overall number multiplies the view numbers, each raised to its weight, so that grains "
                 "two views both count are counted once. Overlaps are measured on samples simulated from the "
                 "fitted units, where only shared counting can make two views agree.\n")
        s.append(table("weights", res["weights"], n=20))
    s.append("### Samples most likely to belong to another allowed unit (M8)\n")
    s.append("`ratio_overall` = how many times more likely the candidate unit is than the coded one. It is the "
             "view ratios multiplied together after each is raised to its weight (table above), so age and "
             "lithology, which count mostly the same grains, together count about once.\n")
    sm = res["samples"]
    if len(sm) and "ratio_overall" in sm:
        top = sm[sm.ratio_overall.notna()].sort_values("ratio_overall", ascending=False)
        extra = [x for x in ("ratio_individual",) if x in top]
        s.append(table("samples", top, ["row", "hole", "depth_ft", "stratunit", "candidate", "relation", "n_grains"]
                          + VIEW_COLS + extra + ["p_typical_of_coded", "p_typical_of_candidate", "withheld"], n=25))
    else:
        s.append("_No sample had another allowed unit._\n")
    s.append("### Contacts (M9) and mixing (M10)\n")
    ct = res["contacts"]
    if len(ct):
        s.append(table("contacts", ct.sort_values("ratio_best_vs_geologist", ascending=False),
                          ["hole", "upper_unit", "lower_unit", "n_upper", "n_lower", "geologist_contact_ft",
                           "best_contact_ft", "samples_that_change", "change_to", "ratio_best_vs_geologist",
                           "prob_geologist_position"] + VIEW_COLS[1:], n=25))
        mx = res["mixing"]
        s.append("Samples at those contacts as mixtures (fraction of the upper unit; fit relative to the better "
                 "pure unit):\n")
        s.append(table("mixing", mx.sort_values("mix_vs_best_pure", ascending=False), n=15))
    else:
        s.append("_No borehole has two units stacked among clean samples._\n")
    s.append("### Adjacent units: one unit or two? (M11)\n")
    s.append("Percentages assume equal prior odds. `separation_vs_same_code_pct` = share of same-code borehole "
             "pairs that are closer together than this pair of units.\n")
    s.append(table("merges", res["merges"], n=40))
    s.append("### Units: one unit or an upper and a lower one? (M12)\n")
    sp = res["splits"]
    keep = [c for c in sp.columns if not (c.startswith("upper_") or c.startswith("lower_"))]
    s.append(table("splits", sp, keep, n=40))
    if res.get("grain_selection") is not None:
        s.append("### Individual grains in the member model (M13)\n")
        s.append(table("grain_selection", res["grain_selection"], n=30))
    return "\n".join(s)


def write_report(path, meta, sections):
    head = [f"# Till code review: {meta['county']} County\n",
            f"Data: `{meta['data']}`. Seed {meta['seed']}. Buffer {meta['buffer_km']} km. "
            "Methods M1-M14: `references/statistical_methods.md`. Priors: `references/geological_priors.md`. "
            "The original `stratunit` codes are never changed; everything below is evidence for a geologist.\n",
            "## How to read the numbers\n", glossary.HOW_TO_READ,
            "Every table below is followed by the meaning of its columns. `column_guide.md` in this folder "
            "explains every column of every CSV.\n",
            "## Data and bad samples (step 1)\n",
            f"{meta['county_till_samples']} till samples in the county, {meta['county_coded_till']} coded; "
            f"{meta['area_holes']} boreholes in the area including the buffer.\n",
            md_table(meta["bad_summary"]),
            f"Rollup check (M1), share of rows where the groups reproduce the QDI column: {meta['rollups']}.\n"]
    if meta["codes_not_in_lookup"]:
        head.append(f"**Codes not in `stratunit_codes.csv` (treated as their own formation; confirm):** "
                    f"{', '.join(meta['codes_not_in_lookup'])}\n")
    tail = ["## Settings\n", md_table(pd.DataFrame({"setting": list(meta["cfg"]), "value": [str(v) for v in meta["cfg"].values()]}), n=99),
            "## Reproduce\n", f"`python scripts/run_review.py {meta['config_path']}`\n"]
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(head + sections + tail))


def write_guide(path, bad, loops):
    """column_guide.md: every column of every table written in this run."""
    frames = {"bad_samples": [bad]}
    for res in loops.values():
        for key, _ in TABLES:
            if res.get(key) is not None:
                frames.setdefault(key, []).append(res[key])
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(glossary.guide([("bad_samples", "01_bad_samples")] + TABLES, frames))


def write_ledger(path, meta):
    with open(path, "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2, default=str)

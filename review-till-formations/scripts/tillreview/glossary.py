"""Plain-language meaning of every column in every output table.

Written for a geologist with no statistics background: what the number is, then how to
read it. The report prints these under each table, `column_guide.md` lists them all, and a
test fails if an output column has no entry.
"""
import re

HOW_TO_READ = """\
Three kinds of number appear throughout. Nothing here is a verdict: each number is evidence for a geologist
to weigh against what was seen in the core.

- **"Times more likely" (columns starting `ratio_`).** The grains and texture of a sample are compared with
  what two units normally look like. A value of 1 means the sample fits both equally well. A value of 50 means
  it fits the second unit 50 times better than the first. A value of 0.02 means the opposite: it fits the first
  unit 50 times better (1 ÷ 0.02). Values between about 0.3 and 3 say very little either way. Very large values
  are written like 2.4e+07, which means 24,000,000.
- **"How ordinary is this sample" (columns starting `p_typical`).** A score from 0 to 1. It is the share of
  genuine samples of that unit, in that borehole, that would look at least as unusual as this one. 0.5 is a
  perfectly ordinary sample. 0.01 means only 1 genuine sample in 100 would look this odd. The smallest value
  the method can give is 0.001. A ratio only means "belongs to the other unit" when the sample is ordinary for
  that other unit; if it is odd for both, it fits neither.
- **"Chance out of 100" (columns starting `pct_one_unit` or `pct_two_units`).** How likely an explanation is,
  starting from the position that both explanations were equally likely before looking at the data. Near 100 or
  near 0 is a strong statement; 30 to 70 means the data cannot decide. These chances are only as good as the
  number of samples and boreholes behind them, which are always listed beside them.

**The four views.** Every comparison is given overall and for four views on their own: **age** (Precambrian /
Paleozoic / Cretaceous grains), **lithology** (carbonate / crystalline / shale), **crystalline colour** (light /
dark / red) and **texture** (sand / silt / clay); member loops add the individual grains. Each view is judged
on its own. The overall multiplies the view numbers together after raising each to a weight: age and lithology
count mostly the same grains, so each gets about half a vote, while texture (a separate measurement) gets a full
one. The weights are listed in each loop's "Weights" table. Use the overall to spot the clear cases and the view
numbers to see what drives them.

**Why a low grain count matters less.** A sample with 80 grains counted gives a less certain picture than one
with 300, so its numbers are automatically pulled toward "cannot tell".

**Withheld samples.** Samples that break the timeline or look unlike their own unit are kept out of the
description of each unit, so they cannot distort it. They are still compared with the units around them.
"""

VIEW = {"age": "age (Precambrian / Paleozoic / Cretaceous)", "lithology": "lithology (carbonate / crystalline / shale)",
        "colour": "crystalline colour (light / dark / red)", "texture": "texture (sand / silt / clay)"}

COMMON = {
    "row": "Line number of the sample in the data file, so it can be looked up.",
    "hole": "Borehole identifier (`relateid`).",
    "depth_ft": "Sample depth in feet.",
    "stratunit": "The code the geologist assigned. It is never changed by this review.",
    "unit": "The unit being reviewed in this loop: the formation in the formation loop, the member in a member loop.",
    "n_grains": "Number of provenance grains counted in the sample. Below about 100 the sample says little; "
                "its numbers are automatically made less decisive.",
    "n_samples": "Number of clean coded samples of the unit in the county. Fewer than about 10 means every "
                 "number for this unit is weak.",
    "n_boreholes": "Number of boreholes those samples come from. A unit seen in one or two boreholes may only "
                   "reflect those boreholes.",
    "upper": "The unit higher in the timeline (younger).",
    "lower": "The unit lower in the timeline (older).",
    "upper_unit": "The unit above the contact in this borehole.",
    "lower_unit": "The unit below the contact in this borehole.",
    "n_upper": "Number of samples in the upper unit or part.",
    "n_lower": "Number of samples in the lower unit or part.",
}

TABLES = {
    "bad_samples": {
        "problem": "Why the sample cannot be used normally.",
        "effect": "What the review does about it: leaves it out, or uses only its texture or only its grains.",
    },
    "timeline": {
        "position": "Place in the timeline, 1 = uppermost (youngest).",
        "median_depth_ft": "Typical depth at which the unit is found in the county. Only used to break ties.",
        "boreholes": "Number of county boreholes in which the unit is coded.",
    },
    "timeline_pairs": {
        "adjacent": "Yes if the two units are next to each other in the timeline.",
        "boreholes_agree": "Boreholes where the upper unit really is found above the lower one. More = the order "
                           "is better established.",
        "boreholes_disagree": "Boreholes where the two are the other way round. Any value above 0 means the "
                              "codes in those boreholes conflict with the timeline.",
        "order_observed": "No = the two units are never coded in the same borehole, so their order is assumed, "
                          "not seen. A geologist must confirm it.",
    },
    "violations": {
        "deposittype": "QDI deposit type of the row (T = till, O = outwash, L = lake sediment).",
        "status": "`violation` = this sample is out of order however the borehole is read. `ambiguous` = either "
                  "this sample or another one in the borehole is out of order, and the data cannot say which.",
    },
    "outliers": {
        "p_typical": "How ordinary the sample is for the unit it is coded as, using grains and texture together "
                     "(0 to 1). Every sample in this table is below 0.01: fewer than 1 genuine sample in 100 "
                     "would look this odd.",
        "p_grains": "The same score using only the grain count. Small = the grains are what is unusual.",
        "p_texture": "The same score using only sand, silt and clay. Small = the texture is what is unusual.",
        "least_typical_view": "The view in which the sample is most unusual: `age` (Precambrian / Paleozoic / "
                                "Cretaceous), `lithology` (carbonate / crystalline / shale), `colour` (light / dark "
                                "/ red), an individual-grain view (`grains_...`, member loops) or `texture`.",
        "withheld_views": "The views kept out of every unit description and comparison: those in which the "
                          "sample is itself atypical (p below the outlier cut-off). `all` = the whole sample, "
                          "because it is odd only taken as a whole or odd in every view.",
        "views_kept": "The views of this sample that are still used, because they are ordinary for its unit. "
                      "For example, a leached sample can keep its texture while its grain views are withheld.",
        "explanation": "A possible reason the sample looks odd without being a different unit, or `unexplained`. "
                       "Unexplained samples are the ones most worth checking for a wrong code.",
        "status": "`supported` = the composition changed the way that process predicts, the rest of the sample "
                  "looks normal, and the record independently points to it. `possible` = the composition fits "
                  "but nothing in the record confirms it. `indicator only` = the record points to the process "
                  "but the composition does not show it.",
        "cite": "Section of `geological_priors.md` the explanation comes from.",
        "detail": "The values behind the explanation.",
    },
    "distributions": {
        "view": "Which view the category belongs to (age, lithology, colour or texture).",
        "category": "The grain or texture category.",
        "mean_pct": "The unit's typical percentage of that category in the county.",
        "new_borehole_5": "Low end of what a 250-grain sample of this unit would normally show in a borehole "
                          "that has not been drilled yet (9 samples in 10 fall between the low and high end).",
        "new_borehole_95": "High end of the same range. A wide range means the unit varies a lot from place to place.",
        "same_borehole_5": "Low end of what samples normally show within a single borehole, relative to that "
                           "borehole's own average.",
        "same_borehole_95": "High end of the same range. It is usually much narrower than the new-borehole range: "
                            "a unit is more uniform down one core than across the county.",
    },
    "spread": {
        "view": "The view: `age`, `lithology`, `colour`, an individual-grain view (`grains_...`, member "
                  "loops) or `texture`.",
        "weight": "How much this view counts in the overall number (see the Weights table).",
        "A_within_borehole": "How uniform the unit is down a single borehole in this view. Larger = more "
                             "uniform. Around 20 is loose, 200 is tight.",
        "A_between_boreholes": "How uniform the unit is from one borehole to the next. Larger = more uniform. A "
                               "small value means the unit naturally varies here, so this view counts for "
                               "little when judging its samples.",
        "A_within_all_units": "The same within-borehole figure averaged over all units, for comparison.",
        "A_between_all_units": "The same between-borehole figure averaged over all units, for comparison.",
        "sd_sand_within": "Typical difference in sand % between samples of the unit in one borehole.",
        "sd_clay_within": "Typical difference in clay % between samples of the unit in one borehole.",
        "sd_sand_between": "Typical difference in average sand % between boreholes of the unit.",
        "sd_clay_between": "Typical difference in average clay % between boreholes of the unit.",
    },
    "samples": {
        "withheld": "Blank = a normal sample. `outlier` or `timeline violation` = kept out of the unit "
                    "descriptions, but still compared here. `outlier (texture withheld)` and the like = only "
                    "those views are kept out; the sample's other views still describe its unit.",
        "unit_above": "The different unit coded directly above this sample in its borehole, if any.",
        "unit_below": "The different unit coded directly below this sample in its borehole, if any.",
        "candidate": "The other unit the sample is being compared with. Only units that keep the borehole in "
                     "timeline order are considered.",
        "relation": "`above` / `below` = the candidate is the unit coded above or below in this borehole. "
                    "`between` = a unit the timeline places between those two, which the borehole may have skipped.",
        "p_typical_of_coded": "How ordinary the sample is for the unit it is coded as (0 to 1). 0.5 = ordinary; "
                              "0.01 = 1 genuine sample in 100 would look this odd.",
        "p_typical_of_candidate": "How ordinary the sample would be as a sample of the candidate unit. If this "
                                  "is also small, the sample fits neither unit and a large ratio only means "
                                  "'less bad'.",
        "ratio_overall": "How many times better the sample's grains and texture fit the candidate unit than the "
                         "unit it is coded as, in this borehole. 1 = fits both equally. 50 = fits the candidate "
                         "50 times better, so the code deserves a second look. 0.02 = fits its own code 50 times "
                         "better.",
        "prob_candidate": "The same comparison as a chance from 0 to 1 that the sample is the candidate, if those "
                          "two units were the only possibilities and equally likely beforehand. 0.5 = cannot tell.",
        "ratio_individual": "Member loops only: the same comparison using only the individual grain types that "
                            "separate the members of this formation.",
        "ratio_individual_no_confusable": "As `ratio_individual`, leaving out grain types that counters commonly "
                                          "confuse. If the two differ a lot, the result may come from how the "
                                          "sample was counted.",
    },
    "contacts": {
        "geologist_contact_ft": "Depth interval in which the contact lies according to the codes (between the "
                                "last sample of the upper unit and the first of the lower).",
        "best_contact_ft": "Depth interval where the grains and texture best support the contact. The same as "
                           "the geologist's when the data agree.",
        "samples_that_change": "Number of samples that would switch unit if the contact moved there. 0 = the "
                               "data agree with the logged contact.",
        "change_to": "The unit those samples would switch to.",
        "ratio_best_vs_geologist": "How many times better the best contact position fits than the logged one. "
                                   "1 = they are the same position. Below about 20 is a weak suggestion; "
                                   "thousands or more is strong.",
        "prob_geologist_position": "Of all the places the contact could sit between these two units, the share "
                                   "of support (0 to 1) for where it was logged. 0.9 = well supported. "
                                   "0.05 = another depth is clearly preferred.",
        "prob_best_position": "The same share for the best position. Low values for both mean several positions "
                              "are about equally good.",
        "mean_from_this_borehole_only": "Lists a unit that occurs in no other county borehole. Its normal "
                                        "composition is then known only from this borehole, so the comparison is weak.",
        "changed_depths_ft": "Depths of the samples that would switch unit.",
        "changed_rows": "Data-file line numbers of those samples.",
    },
    "contact_layers": {
        "coded_unit": "The unit the layer is coded as (in the formation loop, its formation).",
        "proposed_unit": "The unit the layer would belong to if the contact moved to the best-supported depth.",
        "ratio_proposed_vs_coded": "How many times better this layer's grains and texture fit the proposed unit "
                                   "than the unit it is coded as, each judged as it was logged in this borehole "
                                   "(the layer itself left out). 1 = fits both equally. 50 = fits the proposed "
                                   "unit 50 times better. Below 1 = this layer on its own favours its code, even "
                                   "if the contact as a whole moves.",
        "prob_proposed": "The same comparison as a chance from 0 to 1 that the layer belongs to the proposed "
                         "unit, if only these two units were possible and equally likely beforehand.",
        "p_typical_of_coded": "How ordinary the layer is for the unit it is coded as, in this borehole (0 to 1). "
                              "0.5 = ordinary; 0.01 = 1 genuine sample in 100 would look this odd.",
        "p_typical_of_proposed": "How ordinary the layer would be as a sample of the proposed unit. If this is "
                                 "also small, the layer fits neither unit and a large ratio only means 'less bad'.",
        "mix_fraction_upper": "The blend of the upper and lower unit that best matches the layer, as the fraction "
                              "of the upper unit (M10). Well away from 0 and 1 = the layer may be a mixture.",
        "mix_vs_best_pure": "How many times better that blend fits than either unit alone. Never below 1; large "
                            "values with a blend well away from 0 and 1 point to a mixed layer, which should be "
                            "reported as mixed rather than moved.",
        "geologist_contact_ft": "Depth interval of the contact as logged.",
        "best_contact_ft": "Depth interval where the grains and texture best support the contact.",
        "contact_ratio_best_vs_geologist": "For the whole contact: how many times better the best position fits "
                                           "than the logged one, using all samples of both units in the borehole.",
        "layers_moved_at_contact": "How many layers this contact move would switch in total.",
        "ratio_individual": "Member loops only: the same comparison using only the individual grain types that "
                            "separate the members of this formation.",
        "ratio_individual_no_confusable": "As `ratio_individual`, leaving out grain types that counters commonly "
                                          "confuse. If the two differ a lot, the result may come from how the "
                                          "layer was counted.",
    },
    "mixing": {
        "contact_interval_ft": "The logged contact the sample was tested against, as a depth range: from the last "
                               "clean till sample of the upper unit to the first clean till sample of the lower unit. "
                               "Non-till layers and withheld samples are skipped, so a wide range means the two "
                               "tills do not touch here and the sample cannot have picked up the other till "
                               "directly.",
        "mix_fraction_upper": "The blend of the two units that best matches the sample, as the fraction of the "
                              "upper unit. 0.5 = looks like half of each, as when the base of a till has picked "
                              "up the till below. Near 1 or 0 = looks like the upper or lower unit alone, not a mixture.",
        "mix_vs_best_pure": "How many times better that blend fits than the better of the two units alone. It is "
                            "never below 1. Up to about 5 is expected by chance; large values with a fraction "
                            "well away from 0 and 1 point to a genuinely mixed sample.",
        "upper_vs_lower": "How many times better the sample fits the upper unit than the lower one. Above 1 "
                          "favours the upper unit, below 1 the lower.",
    },
    "weights": {
        "view": "A view of the sample: `age`, `lithology`, `colour`, an individual-grain view (`grains_...`, "
                "member loops) or `texture`.",
        "weight": "How much the view counts in the overall number. 1 = a full, independent vote. About 0.5 = "
                  "half a vote, because another view counts mostly the same grains. Texture and individual-grain "
                  "views share no counts with any other view and always count 1.",
    },
    "merges": {
        "boreholes_upper": "Boreholes in which the upper unit is coded.",
        "boreholes_lower": "Boreholes in which the lower unit is coded.",
        "boreholes_with_both_stacked": "Boreholes where the two are coded one above the other. A geologist saw a "
                                       "contact there, which is evidence against them being one unit.",
        "order_observed_in_boreholes": "Boreholes in which the order of the two can be seen. 0 = their order is assumed.",
        "pct_one_unit": "The chance out of 100 that the two codes describe one unit, judged against how much a "
                        "single unit normally varies between boreholes. Near 100 = the data cannot tell them "
                        "apart. Near 0 = they differ more than one unit would. Weak if either unit has only a "
                        "few samples or boreholes.",
        "one_vs_two_ratio": "The same result as 'times more likely': how many times better 'one unit' explains "
                            "the data than 'two units'. Above 1 favours one unit, below 1 favours two.",
        "separation_vs_same_code_pct": "How far apart the two units are compared with boreholes that share one "
                                       "code. 50 = about as far apart as two typical boreholes of the same "
                                       "unit. Near 100 = further apart than almost any two boreholes of one unit.",
        "pct_one_unit_with_uncoded": "The chance out of 100 again, this time also using nearby uncoded till "
                                     "samples as evidence. The number without them is the main one.",
        "uncoded_samples_used": "How many uncoded samples were used.",
        "uncoded_weight_used": "Their combined weight. Each counts for at most half a coded sample.",
    },
    "splits": {
        "result": "`tested`, or why the unit could not be tested (too few samples in any one borehole).",
        "pct_two_units": "The chance out of 100 that the code covers two units, one above the other, in the same "
                         "order wherever both are present. Near 100 = the upper and lower parts differ "
                         "consistently. Near 0 = one unit. Weak if `boreholes_showing_both` is 1. It is judged "
                         "against one unit whose composition may drift down each core, so a core that simply "
                         "changes down its length does not count as two units; only upper and lower parts with "
                         "different typical compositions do. Read it with `vertical_vs_lateral_ratio` and the "
                         "`same_direction_` columns.",
        "two_vs_one_ratio": "The same result as 'times more likely': above 1 favours two units, below 1 one unit.",
        "lateral_vs_one_ratio": "Lateral check: how many times better the same two groups fit when every borehole "
                                "is kept whole (each borehole entirely upper or entirely lower), compared with one "
                                "unit. This is a difference between areas or projects, not one unit above another.",
        "vertical_vs_lateral_ratio": "How much more evidence there is for two typical compositions when the parts "
                                     "are upper and lower than when the same two groups are whole boreholes. Well "
                                     "above 1 = the evidence comes from changes down the boreholes. Near or below 1 = "
                                     "the split is really a difference between boreholes (area, project or "
                                     "counter), not two stacked units, however high `pct_two_units` is.",
        "pct_two_units_individual": "Member loops only: the chance out of 100 that the code covers two units, "
                                    "judged on the individual grain types that separate the members alone.",
        "pct_two_units_individual_no_confusable": "As `pct_two_units_individual`, leaving out grain types that "
                                                  "counters commonly confuse. If the two differ a lot, the split may "
                                                  "come from how samples were counted.",
        "boreholes_showing_both": "Boreholes that contain both the upper and the lower part. If this is 1, the "
                                  "split rests on a single core and may be a local change.",
        "codes_upper": "The original codes of the samples in the upper part, with counts. If the two parts carry "
                       "different member codes, the split is confirming the members.",
        "codes_lower": "The original codes of the samples in the lower part, with counts.",
        "pct_two_units_with_uncoded": "The chance out of 100 again, also using nearby uncoded till samples. The "
                                      "number without them is the main one.",
        "uncoded_samples_used": "How many uncoded samples were used.",
    },
    "split_boundaries": {
        "boundary_ft": "Depth interval where the split falls in this borehole, or `all upper` / `all lower` if "
                       "the borehole holds only one part.",
    },
    "uncoded": {
        "unit": "A unit the sample could belong to without breaking the timeline order in its borehole.",
        "posterior": "The chance (0 to 1) that the sample is this unit, among the units its position allows.",
        "weight": "How much the sample counts toward this unit in the one-unit / two-unit tests. At most 0.5 (half "
                  "a coded sample); 0 if the sample does not look like any allowed unit.",
        "most_likely_unit": "The allowed unit the sample most resembles. For a sample coded only to formation, "
                            "this is the member it most resembles.",
        "p_typical_of_most_likely": "How ordinary the sample is for that unit (0 to 1). Below 0.01 it is not used.",
        "units_allowed": "How many units the sample's position allows. With 1, the position alone decided.",
    },
    "grain_selection": {
        "grain": "An individual grain type.",
        "group": "The broad group it belongs to.",
        "G": "How strongly the grain's share of its group differs between the members. 0 = no difference; larger "
             "= a clearer difference.",
        "p": "How easily a difference that large could arise by chance if the members were really the same "
             "(0 to 1). Small = unlikely to be chance.",
        "q": "The same, adjusted for the fact that many grain types were tested at once. This is the value the "
             "selection uses.",
        "selected": "Yes if the grain is used to tell this formation's members apart (q below 0.05).",
        "confusable": "Yes if counters commonly confuse this grain with another (geological_priors.md section 2). "
                      "A difference in such a grain may reflect who counted the sample.",
    },
}

PATTERNS = [
    (r"^ratio_(age|lithology|colour|texture)$",
     lambda m: f"The same 'times more likely' comparison using only {VIEW[m.group(1)]}. Shows whether this "
               f"view supports or opposes the overall number. The overall is the view numbers multiplied "
               f"together, each raised to its weight."),
    (r"^same_direction_(age|lithology|colour|texture)$",
     lambda m: f"Of the boreholes that hold both parts, how many show the upper part differing from the lower part "
               f"in {VIEW[m.group(1)]} the same way as the split overall. A real split points the same way in "
               f"nearly every borehole; drift down the cores scatters both ways (about half)."),
    (r"^overlap_with_(\w+)$",
     lambda m: f"How much this view's evidence moves together with the {m.group(1)} view's when the true unit is "
               f"known (0 to 1): 0 = they count different grains, near 1 = they count the same grains."),
    (r"^mix_fraction_upper_(age|lithology|colour|texture)$",
     lambda m: f"The best blend of the two units judged on {VIEW[m.group(1)]} alone, as the fraction of the "
               f"upper unit."),
    (r"^mix_vs_best_pure_(age|lithology|colour|texture)$",
     lambda m: f"How many times better the best blend fits than either unit alone, judged on {VIEW[m.group(1)]} "
               f"alone. Never below 1. A large value in one view only, with the others near 1, is weak evidence "
               f"of mixing."),
    (r"^pct_one_unit_(age|lithology|colour|texture)$",
     lambda m: f"The chance out of 100 that the two are one unit, judged on {VIEW[m.group(1)]} alone."),
    (r"^pct_two_units_(age|lithology|colour|texture)$",
     lambda m: f"The chance out of 100 that the code covers two units, judged on {VIEW[m.group(1)]} alone. The "
               f"views can each favour a split at different depths, so they can be high while the overall is low."),
    (r"^(upper|lower)_(\w+)_pct$",
     lambda m: f"Typical percentage of {m.group(2)} in the {m.group(1)} part of the proposed split."),
    (r"^pct_in_(\w+)$",
     lambda m: f"The grain as a percentage of its broad group in member {m.group(1)}."),
]


def explain(table, column):
    """Plain-language meaning of a column of an output table, or None if there is none."""
    text = TABLES.get(table, {}).get(column) or COMMON.get(column)
    if text:
        return text
    for pattern, make in PATTERNS:
        m = re.match(pattern, column)
        if m:
            return make(m)
    return None


def column_notes(table, columns):
    """Markdown list explaining the given columns of a table."""
    lines = ["**Columns**\n"]
    for c in columns:
        text = explain(table, c)
        if text:
            lines.append(f"- `{c}`: {text}")
    return "\n".join(lines) + "\n"


def guide(table_names, frames):
    """The full column guide: every column of every table that was written."""
    out = ["# Column guide\n", "What every column in the output tables means, in plain language.\n",
           "## How to read the numbers\n", HOW_TO_READ]
    for key, name in table_names:
        cols = list(dict.fromkeys(c for f in frames.get(key, []) for c in f.columns if c != "idx"))
        if cols:
            out += [f"## `{name}.csv`\n", column_notes(key, cols)]
    return "\n".join(out)

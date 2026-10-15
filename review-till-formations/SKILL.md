---
name: "review-till-formations"
description: "Review till formations and members previously identified by geologists based on texture and 1-2mm grain count data, anchored on rotasonic cores, with one chronological order every borehole must follow."
---

# Till formation review

Statistically review geologist labeled till units based on raw borehole data, one county at a time.  Every conclusion 
must be **repeatable** (a geologist should be able to follow a clean trail of exact statistical methods and logic based 
decisions, WHEN they were used, WHY they were used, and HOW they were used. especially when a geological prior made a 
choice, rather than pure statistics) and **cited** (statistical methods must be cited by number found in 
`references/statistical_methods.md` and decisions based on geological priors must be cited by number found in 
`references/geological_priors.md`).

## The data situation

- **Rotasonic cores** (`rotasonic = True`, type RS): few, deep, often to bedrock. Most formation codes (`stratunit`) 
  were defined using these cores, and they typically see **every** till in the general vicinity
- **Shallow holes** (giddings GD, hand samples HS): many, typically < 25 ft, mostly uncoded.  They give lateral coverage 
  of the uppermost tills only.
- **Grain counts** are hand counts of an unknown composition: multinomial within a sample,
  overdispersed between samples. **Texture** is weight %, a composition.

Grain counts come in both individual counts (24 separate grains) and 9 broad categories.  When identifying formations 
geologists almost **exclusively** use some or all of the 9 broad categories + the 3 texture categories, almost never 
the individual counts.  This helps reduce human-error with similar looking grains, but as a consequence can omit useful 
information.

Additionally, a units within-core spread is typically much smaller than a unit's lateral spread.

True zeros are common (no shale, no red grains, no carbonate), so grain counts are analysed as counts and texture
as measured percentages. Nothing is log-transformed and no zero is ever replaced (M1, M2).

## Two loops

The review runs twice, and the two loops are **not** the same:

| | Formation loop | Member loop (once per formation) |
|---|---|---|
| Units | member codes rolled up to their formation (`references/stratunit_codes.csv`) | the member codes of one formation |
| Grain variables | the broad categories only (M1) | the broad categories **plus** the individual grains that measurably separate that formation's members (M13) |
| Samples | every coded till sample | the formation's samples, minus those the formation loop found 20× or more likely to be another formation, or out of order |
| Formation-level codes (QNUU, QLHU, …) | a unit | "member unknown": never wrong at member level, used as evidence only (M14) |
| Timeline | order of formations | order of members inside the formation, equally strict |

Both loops use texture. The formation loop never uses individual grain counts.

## How numbers are reported

- **No verdicts or tiers.** Report the exact number and let the geologist decide: "this sample is 31× more
  likely QHWU than QNUM", "QLHS and QLHM are less than 0.01 % likely to be one unit".
- **Always the overall number and the four views**: age, lithology, crystalline colour, texture. Each view is
  modelled on its own. The overall is the view numbers multiplied together after each is raised to its weight
  (M1): age and lithology count mostly the same grains, so they share about one vote; texture and the individual
  grains always count fully. The weights are in each loop's "Weights" table (`11_view_weights.csv`). Use the
  overall to spot clear cases and the views for the deeper look. In the member loop add the fifth number,
  `ratio_individual`.
- A category that is innately wide in a formation automatically counts for little in that formation, because
  each unit has its own spread per view (M5). There is no table of "what each formation cares about": it is
  read from `05_spread.csv`.
- Percentages ("X % likely one unit") assume equal prior odds. State that, and give the ratio beside it.
- Never report a ratio without the sample's typicality (`p_typical_of_coded`, `p_typical_of_candidate`) and
  grain total. A sample atypical of both units fits neither.
- The original `stratunit` is never changed. Everything is evidence for a geologist.

## Explaining the numbers to a geologist

The reader is a geologist who may have no statistics background. Every number you present must come with what
it is and what it means for the code under review, in plain words.

- Never leave a statistical term unexplained. Do not write "likelihood ratio", "posterior", "p-value",
  "Bayes factor", "concentration" or "dispersion" without the plain restatement below, and prefer the plain
  restatement alone.
- Say what the number means for this sample or unit, not only its value: "the sample at 92 ft fits Lake Henry
  855 times better than New Ulm, and looks like an ordinary Lake Henry sample" rather than "ratio 855, p 0.03".
- Always say what a number rests on: grains counted, samples and boreholes.
- Every table you show must be followed by the meaning of its columns. `report.md` already does this under each
  table, and `column_guide.md` in the output folder explains every column of every CSV. Use that wording
  (source: `scripts/tillreview/glossary.py`); do not invent a different explanation.

The three kinds of number, and how to explain them:

| Columns | What it is, in plain words | How to read it |
|---|---|---|
| `ratio_overall`, `ratio_best_vs_geologist`, `one_vs_two_ratio`, `two_vs_one_ratio`, `mix_vs_best_pure`, `upper_vs_lower` | "Times more likely": how many times better the grains and texture fit one explanation than the other | 1 = fits both equally. 50 = fits the second 50 times better. 0.02 = fits the first 50 times better. Between about 0.3 and 3 says very little. 2.4e+07 means 24,000,000 |
| `ratio_age`, `ratio_lithology`, `ratio_colour`, `ratio_texture`, `ratio_individual`, `mix_*_<view>` | The same comparison using one view only | Shows which evidence drives the overall number. The overall = each view number raised to its weight, multiplied; age and lithology share about one vote because they count mostly the same grains |
| `weight`, `overlap_with_*` | How much a view counts in the overall, and how much its evidence repeats another view's | 1 = a full independent vote; about 0.5 = half a vote |
| `p_typical`, `p_typical_of_coded`, `p_typical_of_candidate`, `p_grains`, `p_texture` | "How ordinary is this sample" for a unit, from 0 to 1: the share of genuine samples of that unit in that borehole that would look at least as unusual | 0.5 = perfectly ordinary. 0.01 = only 1 genuine sample in 100 would look this odd. Small for both units = the sample fits neither |
| `pct_one_unit*`, `pct_two_units*` | "Chance out of 100" that two codes are one unit, or that one code covers two, if both answers were equally likely beforehand | Near 100 or near 0 is a strong statement; 30 to 70 means the data cannot decide. Only as good as the samples and boreholes behind it |
| `prob_geologist_position`, `prob_best_position`, `prob_candidate`, `posterior` | A share of support from 0 to 1 among the possibilities considered | 0.9 = well supported. 0.05 = another possibility is clearly preferred |
| `mix_fraction_upper` | The blend of the two units that best matches the sample | 0.5 = looks like half of each. Near 0 or 1 = not a mixture |
| `A_within_borehole`, `A_between_boreholes`, `sd_*` | How tightly a unit holds its composition, down one borehole and from borehole to borehole | Larger A = more uniform. A small between-borehole value means the unit naturally varies there, so that part counts for little |
| `n_grains`, `n_samples`, `n_boreholes`, `boreholes_*` | What the result rests on | Under about 100 grains, 10 samples or 3 boreholes, treat the numbers beside them as weak |

## Running a review

Scripts are in `scripts/` (numpy, scipy, pandas; pytest for the tests). Before the first run on a machine:

    python -m pytest -q scripts/tests

Do not analyse real data until the tests pass.

Copy `assets/config.example.json`, set `data`, `county` and `out`, then:

    python scripts/run_review.py config.json --stop-after timeline     # steps 1-2
    python scripts/run_review.py config.json                           # everything

A county of about 650 till samples takes about two minutes. The same config and seed give the same numbers.

## Review process

Each loop follows the same steps. Cite the method (M#) or prior (§#) wherever a number or decision appears.

0) **Scope** Confirm the county and data file with the user. If the report lists codes that are not in
   `references/stratunit_codes.csv`, ask what formation they belong to and add confirmed ones to the table.
1) **Identify bad samples within the dataset** Find and filter out samples that can cause problems with statistical methods.
   - Zero grain count samples: Only use these when looking at **Texture** data, never include in grain count statistics
   - Sample depth negative or NaN: depth is indeterminate and as such unusable
   - Missing X/Y coordinates
   - Texture all zero or not summing to 100 ± 2: never measured, use grains only (M2)

   Depths recorded in metres are converted to feet. Output: `01_bad_samples.csv`.
2) **Create geological timeline** (M3) Create a timeline of formations based on existing geological codes. All formations MUST
   follow this ordering. A borehole can SKIP over formations, but never re-order them.

   **Stop here.** Run with `--stop-after timeline` and show the user, for the formation loop and each member loop:
   the derived order, every pair with a borehole that disagrees, and every adjacent pair never seen in one
   borehole (its order is not observed, only assumed). Record the confirmed or corrected orders in the config
   under `timeline` and run again. A report whose timeline says "derived, NOT yet confirmed" is provisional and
   must be presented as such.
3) **Find timeline violations** (M4) In each borehole, the smallest set of samples whose removal restores the
   order. `violation` = wrong in every reading; `ambiguous` = either it or another sample is wrong. Both are
   reported and **withheld** from the distributions. Output: `03_violations.csv`.
4) **Identify sample outliers** (M7) A sample is an outlier if fewer than 1 in 100 genuine samples of its unit, in
   its borehole, would be as improbable. Each outlier is then explained:
   - Check if the sample has a low grain count 
   - Check if the sample has had any geological processes applied to it (`references/geological_priors.md`, §1)
   - Check if the sample was possibly miscounted (`references/geological_priors.md`, §2)

   or recorded as **unexplained**. A process is `supported` (direction + typical without the altered parts +
   an independent indicator), `possible` (no indicator) or `indicator only`, as §1 defines.
   Only what is atypical is **withheld**: each view of an outlier (age, lithology, colour, individual grains,
   texture) has its own p, and views below the cut-off are kept out of every statistical method while the
   sample's ordinary views are still used. A leached sample loses its grain views but keeps its texture; a
   sample with sand mixed in loses its texture but keeps its grains. The whole sample is withheld only when it is
   odd as a whole with no single odd view, or odd in every view (`withheld_views`, `views_kept`). Unexplained
   outliers are the most likely miscodes: look them up in step 6.
   Output: `04_outliers.csv`.
5) **Compute current formation distributions** (M5) With violations and withheld outlier measurements omitted, each unit's proportions
   and its two spreads: between boreholes and within a borehole. Spread borrows from boreholes within the buffer
   (25 km); means use the county only. Output: `05_distributions.csv`, `05_spread.csv`.
6) **Test samples and contacts within boreholes** (M8, M9, M10)
   - *Samples* (M8): every coded sample against the unit coded above it, the unit below it, and any unit the
     timeline puts between them. `ratio_overall` = how many times more likely the candidate is than the code.
     Low confidence is a ratio near or above 1; there is no cut-off.
   - *Contacts* (M9): reassigning an edge sample is the same as moving the contact, so every contact between two
     stacked units is also tested as a whole: the best-supported position against the geologist's. For a closer
     look at one contact, `06_contact_layers.csv` gives each sample it would move with that sample's own ratio.
   - *Mixing* (M10): before presenting a contact move, check each sample that would change for mixing (§1.9)
     and for spanning the contact (§1.C). A mixed sample is reported as mixed, not reassigned.

   A lone odd sample (M8) and a misplaced contact (M9) are different findings; report both.
   Output: `06_sample_likelihoods.csv`, `06_contacts.csv`, `06_contact_layers.csv`, `06_mixing.csv`.
7) **Find adjacent units that may be one** (M11) Every pair adjacent in the timeline: how likely the two are one
   unit rather than two, judged against how much one unit varies between boreholes. Boreholes where the two are
   stacked with a contact are evidence against and must be stated. Output: `07_merge_pairs.csv`.
8) **Find units that may be two** (M12) Every unit: how likely it is an upper and a lower unit, in the same order
   in every borehole that shows both. Each borehole is judged against compositions from the other boreholes,
   and against one unit whose composition may drift down the core, so drift and chance do not make a split.
   Vertical splits only; a difference between areas is lateral variation and is not proposed.
   Output: `08_split_candidates.csv`, `08_split_boundaries.csv`.

Steps 7 and 8 are given with and without uncoded till samples as evidence (M14). The number without is primary.

## Reading the results

- `report.md` in the output folder has the top rows of every table, each followed by the plain-language meaning
  of its columns; the CSVs in `formation/` and `members_<code>/` have all rows. `column_guide.md` explains every
  column of every CSV. `ledger.json` records the config, seed, versions and timelines. Samples are identified by
  `row` (line in the data file), `hole` and `depth_ft`.
- A formation with distinct members is expected to show as a split in the formation loop. Check `codes_upper`
  and `codes_lower`: if the split follows the member codes, it confirms the members, not a new formation.
- A unit that truly holds two sub-units is wide, so its sub-unit samples can be withheld as outliers in step 4
  before step 8 sees them. Several unexplained outliers in one unit that lean the same way are a possible
  split: say so.
- Always read a split's `pct_two_units` with `vertical_vs_lateral_ratio`. Near or below 1, the "split" is a
  grouping of whole boreholes (area, project or counter), which is lateral variation: report it as such, never
  as two stacked units. Check `08_split_boundaries.csv` and the boreholes' `project` / `counted_by`.
- Read a split with its `same_direction_<view>` counts: a real split points the same way in nearly every
  borehole showing both parts; about half means the parts do not differ consistently. A homogeneous unit
  typically reads 5-30 % "two units" (cannot decide, leans one), not 0 %.
- `boreholes_showing_both = 1` on a split, or a unit with one or two samples in a merge, means the number rests
  on almost nothing. Say what it rests on.
- Read every mixing row with its `contact_interval_ft`. A wide interval means non-till layers or withheld samples
  separate the two tills: the sample cannot have picked up the other till directly (§1.9), so a high
  `mix_vs_best_pure` there means "intermediate in composition", not "mixed at the contact". A sample whose depth
  lies outside the interval is one the best-supported contact would move.
- A contact row with `mean_from_this_borehole_only` compares the unit with itself in that borehole; treat its
  ratio as weak.
- When a member finding depends on `felsic` vs `clear_quartz`, `mafic_igneous` vs `metasedvol`, `inoceramus` or
  `red_volcanic`, compare `ratio_individual` with `ratio_individual_no_confusable` and cite priors §2.

## Report to the user

Write it so a geologist with no statistics background can follow every number (see "Explaining the numbers to a
geologist"). Open with the "How to read the numbers" text from `report.md`, and follow every table with the
meaning of its columns.

1. **Summary**: the 5-8 strongest findings, one sentence each, in plain words, with the number, what it means
   and what it rests on.
2. **Data**: samples, exclusions (step 1), codes not in the lookup table.
3. **Timeline** per loop: order, whether confirmed, disagreements, pairs whose order is not observed.
4. **Violations and outliers**, with each outlier's explanation and citation.
5. **Samples and contacts**: highest ratios first, each with the four views, typicality and mixing.
6. **Merge and split candidates**, with borehole counts and the with/without-uncoded numbers.
7. **Priors invoked**: every place a geological prior, not a statistic, decided something.
8. **Settings and how to reproduce** (from `report.md`).

## Things not to do

- Do not compute statistics on grain percentages or log-ratios; counts only (M1). Do not replace zeros.
- Do not add the nine categories as nine variables, or multiply the four view numbers together without their
  weights (M1).
- Do not use the county-wide spread to judge a sample against its own borehole, or the reverse (M6 handles both).
- Do not let a withheld measurement into a distribution: a whole withheld sample, or the withheld views of an
  outlier that is used in part.
- Do not turn a ratio into a verdict ("miscoded", "should be merged"). Give the number and the evidence.
- Do not use individual grain counts in the formation loop.
- Do not propose a relabel that breaks the timeline.
- Do not present a derived timeline as confirmed, or change a default without recording it in the config.
- Do not present a geological judgement as a statistical result, or the reverse.
- Do not show a table or a number without its plain-language meaning, or add an output column without adding it
  to `scripts/tillreview/glossary.py` (a test fails if one is missing).

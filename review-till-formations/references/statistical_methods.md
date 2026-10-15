# Statistical methods reference

Every number in a report comes from one of the methods below. Cite them by number (M1-M14). Each is implemented
in `scripts/tillreview/` (module named in each entry) and run by `scripts/run_review.py`. Defaults are in the
last table; every one can be overridden in the run config and is printed in the report.

Contents: M1 Grain variables · M2 Texture variables · M3 Timeline · M4 Timeline violations · M5 Unit model ·
M6 Probability of a sample under a unit · M7 Typicality and outliers · M8 Sample likelihood ratio ·
M9 Contact position · M10 Mixing · M11 One unit or two · M12 Vertical split · M13 Member grain selection ·
M14 Uncoded evidence · Defaults · References · Validation record

This file is the technical description. The plain-language meaning of every output column, written for a reader
with no statistics background, is in `scripts/tillreview/glossary.py`; it is printed under each table of
`report.md` and in full in `column_guide.md` in every output folder.

**No log-ratios anywhere.** Zeros are real in this data: among till samples with 100+ grains, 40 % have no
shale, 39 % no red grains, 11 % no carbonate, and 99.6 % have a zero in at least one broad group. A log-ratio
needs every zero replaced by an invented value. Grain counts are therefore modeled as counts (M1) and texture
as measured percentages (M2).

---

## M1 Grain variables: four views and their weights (`variables.py`, `model.py: _fit_weights`)

**When:** every grain statistic. **Why:** counts handle zeros and counting uncertainty with no adjustment. Each
view geologists read is modelled on its own, with its own spread per unit, so a category that is innately wide
in a unit carries little weight for that unit, and every view number can be read and checked alone.

**How.** Eight non-overlapping groups, which reproduce the QDI rollup columns (`precambrian`, `paleozoic`,
`cret_total`, `shale`; checked on every run and printed in the report):

| Group | Grain columns |
|---|---|
| light | felsic, quartzite, clear_quartz |
| dark | mafic_igneous, metasedvol |
| red | iron, red_volcanic, pc_arkosic, pc_qtz_arenite |
| pc_other | pc_other |
| carbonate | carbonate |
| pz_other | paleo_ss, paleo_shale, paleo_other |
| shale | gray_shale, speckled_shale |
| k_other | limestone, inoceramus, pyrite, lignite, ostrandersand, cret_other |

`unknown`, `secondary`, `other_other` are not provenance grains and `chert` is in no group: none enters a statistic.

**The views.** Each view is its own Dirichlet-multinomial (Mosimann 1962) on its own counts:

| View | Parts | Counts it uses |
|---|---|---|
| Age | Precambrian / Paleozoic / Cretaceous | all provenance grains |
| Lithology | carbonate / crystalline (= all Precambrian) / shale | all except pz_other and k_other |
| Crystalline colour | light / dark / red | light + dark + red |
| Individual grains (member loops, M13) | selected grains of one group / rest of that group | that group |
| Texture | sand / silt / clay | M2 |

For counts x of total n and parameters α = A·π,

    log P(x) = lnΓ(n+1) − Σ lnΓ(x_k+1) + lnΓ(A) − lnΓ(n+A) + Σ [lnΓ(x_k+α_k) − lnΓ(α_k)]

π = the unit's proportions in that view, A = concentration (larger = tighter). A view whose total is zero
contributes a factor of 1. Fewer grains give a flatter probability, so a 79-grain sample automatically weighs
less than a 300-grain one. The eight groups are also kept as a hidden composition (`leaf`), fitted the same way.
It is used only to simulate samples (weights below, M7), for the priors checks and to build mixtures (M10). It
never enters a reported number.

**The overall number: a weighted product of the views.** Age and lithology count mostly the same grains (carbonate
is most of Paleozoic, crystalline is all of Precambrian, shale most of Cretaceous), and colour splits the
Precambrian count that age also uses. Multiplying the four view numbers would count that evidence more than once.
The overall is therefore a composite likelihood with overlap-adjusted weights (Varin, Reid & Firth 2011; Pauli,
Racugno & Ventura 2011):

    log overall = Σ_views w_v · log(view number)            e.g. overall ratio = Π ratio_v ^ w_v

    w_v = 1 / Σ_u max(ρ_vu, 0)                              (ρ_vv = 1)

ρ_vu is the overlap of views v and u: the rank correlation of their log ratios **when the true unit is known**.
Each loop simulates 2000 samples from every unit (eight groups, a new borehole, grain totals drawn from the unit's
real samples), scores each against every other unit and the area average in every view, centres each comparison,
and pools. With the truth fixed, only counting shared between views can make two views move together; units that
really differ in several ways do not inflate ρ (agreement between independent views is corroboration and keeps
full weight). Texture and the individual-grain views share no counts with any other view: ρ = 0, w = 1. Views that
are independent get w = 1 each and the overall is their plain product; two views that are the same evidence get
about ½ each and together count once. The weights and overlaps are printed per loop (`11_view_weights.csv`) and
recorded in the ledger.

The same weighted sum is used wherever several views are combined: M7, M8, M9, M10 and, with each view's own
penalty, M11 and M12. Wherever the views are scored at the same positions (M8, M9, M10, M11) the overall number
equals the reported view numbers raised to their weights and multiplied, and can be checked by hand.

*Limitations:* one set of weights per loop: the real overlap varies somewhat with the units compared and the grain
total. Because the overall mirrors the views, chance leanings in several views combine: a single unit given two
codes by borehole comes out as one unit in most, not all, simulated counties (Validation record).

## M2 Texture variables (`variables.py`, `model.py`)

**When:** every texture statistic. **How:** a bivariate normal on sand % and clay % as recorded (silt is 100
minus the two). A true 0 % is an ordinary value. Gravel is reported, never modelled.

A row where sand, silt and clay are all zero, or do not sum to 100 ± 2, has no texture (905 till rows statewide
are all-zero: never measured). Such a sample is judged on grains only.

*Limitation:* the normal ignores the 0-100 bounds. This matters only for units whose sand or clay sits within
about two standard deviations of 0 or 100.

## M3 Timeline (`timeline.py: precedence, best_order, pair_table`)

**When:** step 2, once for formations and once per formation for its members. **Why:** one order that every
borehole must follow; a borehole may skip units, never reorder them.

**How.** Units are the codes that occur on till samples in the county. Coded non-till rows (outwash, lake
sediment) of those units also count for position. For each borehole and each pair of units in it, the unit with
the shallower median depth is "above". W[a,b] = number of boreholes with a above b. The timeline is the order
that minimises the number of borehole-pair disagreements, Σ W[b,a] over pairs placed a-above-b (linear ordering
problem; Grötschel et al. 1984), solved exactly by dynamic programming over subsets for up to 16 units and by
adjacent swaps beyond that. Ties are broken by median depth.

Reported per pair: boreholes agreeing, boreholes disagreeing, and whether the pair is ever seen in one
borehole. **A pair never seen together has no observed order**: its position rests on the tie-break and must be
confirmed by the user. The derived order is used only until the user confirms or corrects it in the config.

## M4 Timeline violations (`timeline.py: violations, find_violations`)

**When:** step 3. **How.** In each borehole the coded rows are read top to bottom as a sequence of timeline
ranks. The largest subset that is in order (longest non-decreasing subsequence) is kept. A sample is

- `ok` if it is in every largest in-order subset,
- `violation` if it is in none,
- `ambiguous` if it is in some: either it or another sample is wrong, and the data cannot say which.

Violations and ambiguous samples are reported and withheld from every distribution. They are still given
likelihood ratios in M8.

## M5 Unit model: two spreads per unit (`model.py: Model`)

**When:** step 5 and everything after. **Why:** a unit varies much more between boreholes than within one, so
one spread cannot serve both "is this sample like its neighbours" and "is this borehole like the unit".

**How, grains.** For each unit u and view j:

    borehole composition   θ_uh ~ Dirichlet(A_L · π_u)        A_L = between-borehole (lateral) concentration
    sample counts          x    ~ DM(n, A_W · θ_uh)           A_W = within-borehole concentration

- **Spreads** are estimated from every clean coded till sample of the unit in the county **and in boreholes
  within the buffer distance** (default 25 km) of a county borehole. (A_W, A_L) maximise the sum of leave-one-out
  predictive log probabilities (M6): each sample is predicted from the other samples of its unit in its
  borehole and from the unit's other boreholes. They are therefore the spreads that make M6 best calibrated.
- **Shrinkage.** The same estimate over all units together gives (A_W0, A_L0). A unit's own estimate is the
  maximum of its log likelihood plus a normal penalty on log A centred on the all-unit value with standard
  deviation `shrink_tau` (1.0). A unit with many samples keeps its own spread; one with few is pulled toward the
  all-unit value. A unit seen in one borehole has no lateral information and takes A_L0.
- **Proportions** π_u come from **county samples only**: a borehole-weighted mean of the borehole compositions
  (weight 1/(ρ_L + 1/E_h), ρ_L = 1/(A_L+1), E_h = effective grains in the borehole), so one deep core does not
  dominate, plus `prior_grains` (1) grain at the area composition so no proportion is exactly zero.

Effective grains of a sample for locating its borehole's composition: n_eff = n(A_W+1)/(A_W+n).

**How, texture.** y = μ_u + b_uh + e, b ~ N(0, S_L), e ~ N(0, S_W). S_W is the pooled within-(borehole, unit)
covariance; S_L is the covariance of borehole means minus S_W/m (method of moments), both shrunk toward the
all-unit value as if it were `shrink_df` (10) samples, with a floor of 0.25 %² (texture is recorded to the
nearest percent). μ_u is the precision-weighted mean of county borehole means.

Reported per unit (`05_distributions.csv`): mean % of every category in the four views, the 5-95 % range of a
250-grain sample in a new borehole (both spreads) and within one borehole (A_W only). `05_spread.csv` lists
A_W and A_L per view with the all-unit values and the view's weight (M1).

*Limitation:* the two-level model is fitted by a moment-matching approximation (below), not by full
integration over borehole compositions.

## M6 Probability of a sample under a unit (`model.py: sample_params, logp, chain`)

**When:** every likelihood in M7-M14. **Why:** a sample must be judged against the unit *as it appears in that
borehole* where the borehole has other samples of the unit, and against the unit's wider lateral spread where
it does not.

**How.** Let S = Σ n_eff·(x/n) over the other samples of unit u in the borehole and E = Σ n_eff. Then

    a      = A_L·π_u + S                       expected proportions = a / Σa
    ρ_post = 1/(A_L + E + 1)                   uncertainty left about this borehole's composition
    ρ      = ρ_post + ρ_W − ρ_post·ρ_W          ρ_W = 1/(A_W + 1)
    sample ~ DM(n, (1/ρ − 1) · a/Σa)

With no other samples this is the unit's county proportions at the full lateral + within spread. With many it
is the borehole's own composition at the within-borehole spread. Texture is the normal equivalent: the borehole
effect has posterior precision S_L⁻¹ + m·S_W⁻¹, and the predictive covariance is its inverse plus S_W.

**The tested sample is always left out**, of its borehole's evidence and of its own unit's proportions.

A set of samples (M9, M11, M12) is scored sequentially: each sample given those before it in the same borehole
(`chain`), which is the joint probability of the set.

## M7 Typicality and outliers (`model.py: typicality`, `pipeline.py`)

**When:** step 4. **How.** Under M6 for the sample's coded unit, R = 999 samples with the same grain total are
simulated (seeded) as the eight groups, so the views keep the counts they share. Each simulated sample is scored
in every view and combined with the M1 weights, exactly as the observed sample is. p = (1 + number of simulated
samples no more probable than the observed one)/(R + 1). The same p is computed for grains alone, for texture
alone (exact: Mahalanobis distance against χ²(2)), and per view (with that view's observed total); the view with
the smallest p is reported as "least typical".

**Outlier rule (the only binary rule in the review besides the member-loop exclusion): p < 0.01 → outlier.**
The model is fitted once on all eligible samples, outliers are identified, and it is refitted once. **Only the
measurements that are atypical are withheld.** Each outlier's views (age, lithology, colour, individual-grain
views, texture) have their own p (above); a view with p < 0.01 is atypical. If some of the sample's measured views
are atypical and at least one is not, only the atypical views are withheld: they are left out of every spread
(M5), proportion or mean (M5) and borehole evidence (M6, M9, M11, M12), and the sample's ordinary views are still
used. A process that alters one part of a sample (leaching removes carbonate grains, sorted sediment shifts the
texture) then costs only that part. When any of age, lithology or colour is withheld, the hidden eight-group
composition is withheld with it. The whole sample is withheld when no single view is atypical (it is odd only
taken as a whole) or when every measured view is. Reported per outlier: `withheld_views` and `views_kept`.
A sample is always judged with all its views when it is itself the one being tested (M7, M8). Each outlier is
then explained by `scripts/geological_priors/` (low count, non-provenance grains, a process from priors §1, a
miscount from priors §2), or recorded as **unexplained**. Unexplained outliers are the ones most likely to be
miscoded.

Process checks use two supporting tests on the eight groups, with the sample's M6 prediction:
- *direction* (priors §1, rule 1): one-sided beta-binomial p < 0.05 that the altered groups are as low (or as
  high) as observed;
- *typical with the altered parts removed* (rule 2), and *typical with a confused pair pooled* (priors §2): the
  same simulated p on the remaining groups, p ≥ 0.01.

*Limitations:* p cannot be smaller than 1/(R+1) = 0.001. A unit that truly contains two members is wide, so
real sub-unit samples can be withheld as outliers before M12 sees them: read a cluster of unexplained outliers
that share a direction as a possible split.

## M8 Sample likelihood ratio (`within_hole.py: sample_table, view_ratios`)

**When:** step 6, for every coded till sample in the county, including withheld ones.

**Allowed alternatives:** the unit coded directly above it in the borehole, the unit directly below, and every
unit the timeline places between those two (a thin unit the borehole may have skipped). Nothing that would
break the order. Above and below are read from rows that are `ok` under M4.

**How.** ratio = P(sample | candidate unit, M6) / P(sample | coded unit, M6). Reported:

| Column | Meaning |
|---|---|
| `ratio_overall` | the views' ratios raised to their weights and multiplied (M1): "the candidate is this many times more likely" |
| `ratio_age`, `ratio_lithology`, `ratio_colour`, `ratio_texture` | each view alone (M1) |
| `prob_candidate` | ratio/(1+ratio): probability of the candidate if only these two were possible, equal priors |
| `p_typical_of_coded`, `p_typical_of_candidate` | M7 under each |

A ratio is only as meaningful as the typicality beside it: a sample atypical of both units (both p small) fits
neither, and a large ratio then says "less bad", not "belongs".

There are **no verdict tiers**. For orientation only, Kass & Raftery (1995) call 3-20 "positive", 20-150
"strong" and above 150 "very strong".

*Limitation:* neighbouring samples in a core are correlated beyond what the borehole effect captures, so
several samples leaning the same way are somewhat less than independent pieces of evidence.

## M9 Contact position (`within_hole.py: contact_table, split_ll`)

**When:** step 6, for every pair of units stacked in a borehole (clean samples only). **Why:** relabelling
samples at a contact is the same as moving the contact; testing the contact pools the evidence and can never
produce an out-of-order result.

**How.** The N samples of the two units, top to bottom. For each position k = 1…N−1, the upper k samples are
scored as one borehole of the upper unit and the rest as one borehole of the lower unit (M6, sequential). Unit
proportions exclude this borehole when the unit occurs elsewhere in the county; otherwise the row is marked
`mean_from_this_borehole_only`.

Reported per contact (`06_contacts.csv`): geologist's and best contact as depth intervals, the samples that
would change and to which unit, `ratio_best_vs_geologist`, the same ratio per view, and the probability of the
geologist's position and of the best position among all positions (uniform prior).

**Also reported per layer (`06_contact_layers.csv`, for a closer look; the report shows the per-contact table).** Every sample the best contact
would move gets its own row: its code and the proposed unit, and its own evidence. The sample is judged
against each unit as logged in that borehole, without the sample itself (M6, the other samples of that unit in
the borehole as evidence): `ratio_proposed_vs_coded` (overall, M1 weights) and each view alone,
`prob_proposed` at equal prior odds, `p_typical_of_coded` and `p_typical_of_proposed` (M7), and its best blend of
the two units (M10). The contact-level ratio is repeated beside it. A contact can move as a whole while one of
its layers, alone, favours its own code; a layer ordinary for neither unit, or best explained as a blend, is
not evidence that it belongs to the proposed unit.

## M10 Mixing (`within_hole.py: mixing`)

**When:** step 6, for the two samples either side of each contact and every sample M9 would move. **Why:**
priors §1.9: the base of a till can incorporate the till below, and a sample can span a contact. A mixture
must be reported as a mixture, not moved.

**How.** Mixing is linear in proportions of the eight groups, as when two tills physically mix. For λ from 0 to 1
in steps of 0.01 the expected composition is λ·(upper unit) + (1−λ)·(lower unit) in the eight groups, and each
view's proportions are read from that blend (an individual-grain view blends within its group, weighted by how
much of the group each unit brings); concentrations blend on the log scale and the texture mean and covariance
mix linearly. Each view is scored on its own and the overall is the M1 weighted sum. Reported: the logged
contact interval the sample was tested against (`contact_interval_ft`, last clean upper to first clean lower
till sample; non-till layers and withheld samples between them are skipped, so a wide interval means the tills
do not touch and §1.9 cannot apply), the best λ
(`mix_fraction_upper`), how many times more likely the best mixture is than the better pure unit
(`mix_vs_best_pure`), upper versus lower, and for each view alone its own best λ and mixture ratio
(`mix_fraction_upper_<view>`, `mix_vs_best_pure_<view>`). A mixture seen in one view only is weak evidence.

*Limitation:* the mixture has one more free parameter than a pure unit, so `mix_vs_best_pure` ≥ 1 always.
Values of a few are expected by chance; read it together with λ (a λ near 0 or 1 is not a mixture).

## M11 One unit or two (`merge.py`)

**When:** step 7, for every pair of units adjacent in the timeline.

**How.** Two models of the pair's clean county samples, with the same spreads in both (the two units' spreads
averaged, weighted by sample count) so that only the means differ:

- *one unit:* one set of proportions and one texture mean; a borehole holding both codes is one borehole;
- *two units:* each code has its own.

Where the two codes are stacked, "one unit" keeps one local composition through the logged contact, so a sharp
change at the contact is evidence for two units. This is deliberate and differs from M12: the contact was placed
by a geologist from the core, not searched for. Letting "one unit" restart there (tested 2026-10-14) discarded
that evidence and judged formations only against their between-borehole spread, which for a formation with
very different members (New Ulm) is wide enough to absorb real differences.

Each borehole's samples are scored jointly (M6), so the difference is judged against lateral variation and a
unit sampled in one borehole cannot look distinct merely because that borehole is. With Δ = log L(two) −
log L(one) = Σ_views w_v·Δ_v (M1 weights), d = Σ_views w_v·d_v, d_v the extra means of view v (2 per
3-part view, 2 for texture, parts − 1 for an individual-grain view), and H_a, H_b the boreholes of each unit,

    log BF(one : two) = −Δ + (d/2)·ln(1 + H_a·H_b/(H_a + H_b))

the unit-information approximation to the Bayes factor (Kass & Wasserman 1995; Kass & Raftery 1995) with the
number of boreholes the difference of means rests on. `pct_one_unit` = BF/(1+BF) **at equal prior odds**;
`one_vs_two_ratio` = BF. Each view alone gets the same calculation with its own Δ_v and d_v, so the overall
ratio is the view ratios raised to their weights and multiplied.

Also reported: boreholes of each unit, boreholes where they are stacked (a geologist saw a contact: evidence
against merging), and `separation_vs_same_code_pct`: the percentage of same-code borehole pairs in the area
whose pooled compositions are closer (Hellinger distance on the eight groups) than the median borehole pair of
these two units. This percentile is descriptive only.

*Limitations:* means are plug-in estimates; a unit with one or two samples can fit its own mean too well, and
the penalty only partly offsets that. Always read the percentage with the sample and borehole counts.

## M12 Vertical split (`split.py`)

**When:** step 8, for every unit. Only vertical splits: an upper and a lower part in the same order in every
borehole. A borehole may hold only one part.

**How.** Every borehole gets a boundary position k (0 = all lower, all = all upper). Starting from the borehole
whose samples gain most from being treated as two groups (four or more coded samples, two each side; the best
three starts are tried), alternate: (1) upper and lower typical compositions from the current assignment, (2)
each borehole's best boundary given those, until nothing changes. This finds the candidate split.

**Score: out of sample, against a drift baseline.** Each borehole h is then predicted from typical compositions
estimated from the **other boreholes only** (leave one borehole out):

    two units:  upper and lower compositions from the other boreholes' upper and lower parts
    one unit:   one composition from the other boreholes; the borehole's local composition may still restart
                at the boundary (the drift baseline)

    log L_h(model) = ln[ mean over positions k of exp(log L_h(k)) ]          (boundary averaged, not picked)
    log BF(two : one) = Σ_h [ log L_h(two) − log L_h(one) ]

with every log L the M1 weighted sum over views. This is a cross-validated predictive comparison (a pseudo-Bayes
factor; Geisser & Eddy 1979). Two problems of a split search are removed:
- **Drift.** Real cores change down their length. Because the baseline may restart at the boundary too, a core
  that only drifts gains nothing from "two units"; only parts with different typical compositions do.
- **Selection.** The split is found by searching for the boundaries and compositions that fit best. Judged on
  the same boreholes, chance differences would look real. No borehole helps estimate the compositions it is
  judged against, so they cannot; no penalty term is needed.
A part found in no other borehole is predicted from the one-unit composition, so a split that rests on one
borehole cannot score on that borehole.

`pct_two_units` at equal prior odds, and the same per view (in member loops also the M13 individual-grain views
together, with and without confusable grains). Reported: samples in each part, boreholes showing both parts,
the original codes in each part (in the formation loop this shows whether the split follows the member codes),
each part's composition, the boundary depth per borehole (`08_split_boundaries.csv`), and per view
`same_direction_<view>`: of the boreholes that hold both parts, how many differ upper-to-lower in the same
direction as the split overall (sign of the dot product of the borehole's upper-minus-lower proportions with the
split's upper-minus-lower typical proportions). A real split points the same way in nearly every borehole; drift
scatters about half each way. This count uses no model.

**Lateral check (reported only, does not change `pct_two_units`).** With the same out-of-sample compositions,
each borehole is also scored kept whole: ln[½ exp(log L_h(all lower)) + ½ exp(log L_h(all upper))], against one
unit kept whole. Reported as `lateral_vs_one_ratio` and `vertical_vs_lateral_ratio` (the evidence for two
compositions as upper and lower parts, divided by the evidence for them as groups of whole boreholes). Because a
borehole may hold only one part, the "split" can be a grouping of boreholes by area, project or counter; a
`vertical_vs_lateral_ratio` near or below 1 means the evidence does not come from order down the boreholes.

*Limitations:* `boreholes_showing_both` = 1 or 2 means the split rests on almost no borehole that shows the
order. An out-of-sample comparison cannot give strong support to "one unit" when nothing differs: a homogeneous
unit reads "cannot decide, leans one" (median 16 % in simulation) rather than near 0 %. The views can each
favour a split at different depths while the overall does not, and each borehole's position is averaged once for
all views together, so the overall is not exactly the weighted product of the view numbers. Each part needs
`split_min_samples` (3) coded samples.

## M13 Member grain selection (`members.py`)

**When:** member loop only, before step 2. **Why:** members may be separated by individual grain types or by
smaller shifts inside a broad group, but most individual types are rare or noisy, and some are confused by
counters.

**How.** For each grain type in a group of two or more types, its count as a share of its group is compared
between the coded members: beta-binomial with one mean versus one mean per member (shared concentration),
G = 2·Δ log L against χ²(members − 1), then Benjamini-Hochberg q across all grains tested. A grain with
q < `member_grain_q` (0.05) is **selected**. The selected grains of a group become one extra view of that
group's grains (selected grains + "rest of group"). It splits counts no other view looks at, so its M1 weight
is 1.

Reported (`10_member_grain_selection.csv`): G, p, q, the grain's % of its group in each member, and
`confusable` (the grain is in a priors §2 pair). Every member result carries `ratio_individual` (the extra
views together) and `ratio_individual_no_confusable` (without views containing a confusable grain).
If a finding disappears without the confusable grains, say so and cite priors §2.

## M14 Uncoded evidence (`uncoded.py`)

**When:** steps 7-8. Uncoded till samples in the formation loop; samples coded only to formation level
("member unknown") in the member loop.

**How.** A unit is allowed for a sample only if it keeps the borehole in timeline order (between the coded rows
above and below). Elevation is used only when the borehole gives no position at all: an uncoded sample in a
borehole with no coded row above or below it must also lie within the unit's county elevation range ±
`uncoded_pad_ft` (15 ft). A member-unknown sample is already placed by its formation code, so every member the
order allows is a candidate and elevation is never used. Among allowed units, with equal priors, posterior ∝ M6 probability. weight = min(`uncoded_cap` (0.5),
posterior), and 0 if the sample is not typical (M7 p < 0.01) of its most probable unit or is flagged leached.
M11 and M12 are reported **with and without** this evidence. `09_uncoded_evidence.csv` lists every sample's
allowed units, posteriors and weights; for member-unknown samples `most_likely_unit` is the member it most
resembles.

*Limitation:* the weights come from the same unit model the tests then examine, so uncoded evidence tends to
make units look more distinct, and many shallow samples can outweigh a few coded ones. The without-uncoded
number is the primary one.

---

## Defaults

| Setting | Default | Method |
|---|---|---|
| `buffer_km` | 25 | M5 |
| `low_count` | 100 provenance grains | M7 explanation |
| `outlier_p` | 0.01 | M7 |
| `simulations` | 999 | M7 |
| `prior_grains` | 1 | M5 |
| `shrink_tau` | 1.0 (log scale) | M5 |
| `shrink_df` | 10 | M5 |
| `member_exclusion_ratio` | 20 | member loop: a sample with `ratio_overall` ≥ this for another formation, or an M4 violation, is left out |
| `uncoded_cap`, `uncoded_pad_ft` | 0.5, 15 ft | M14 |
| `split_min_samples` | 3 | M12 |
| `member_grain_q` | 0.05 | M13 |
| `seed` | 20261012 | all simulation |

## References

- Benjamini Y, Hochberg Y (1995) Controlling the false discovery rate. *JRSS B* 57:289-300.
- Geisser S, Eddy WF (1979) A predictive approach to model selection. *JASA* 74:153-160.
- Grötschel M, Jünger M, Reinelt G (1984) A cutting plane algorithm for the linear ordering problem. *Oper. Res.* 32:1195-1220.
- Kass RE, Raftery AE (1995) Bayes factors. *JASA* 90:773-795.
- Kass RE, Wasserman L (1995) A reference Bayesian test for nested hypotheses and its relationship to the Schwarz criterion. *JASA* 90:928-934.
- Pauli F, Racugno W, Ventura L (2011) Bayesian composite marginal likelihoods. *Statistica Sinica* 21:149-164.
- Varin C, Reid N, Firth D (2011) An overview of composite likelihood methods. *Statistica Sinica* 21:5-42.
- Mosimann JE (1962) On the compound multinomial distribution, the multivariate β-distribution, and correlations among proportions. *Biometrika* 49:65-82.

## Validation record

Known-answer tests on simulated counties: `python -m pytest -q scripts/tests` (19 tests; the last checks that
every output column has a plain-language explanation). They check that
Dirichlet-multinomial probabilities sum to 1; the timeline recovers a planted order and M4 separates violations
from ambiguous samples; a zero count gives a finite ratio with no replacement value, decisive against a unit
rich in that group and silent between two units that lack it; the M1 weights give age and lithology about half
a vote each and texture a full one, and the overall equals the weighted product of the views; a low-count sample
is less decisive; a unit that is wide in one view loses far less probability for an odd value in that view;
typicality flags a planted sample and stays near its nominal rate otherwise; an outlier odd only in its
texture keeps its grain views and one odd only in its grains keeps its texture; planted miscoded contact samples are
found, with one layer row per moved sample; one unit given two codes by borehole comes out as one unit and two
distinct units do not; two units under one code are split at the right depth in every borehole, and two similar
units still are; composition that drifts down cores in random directions is not a split; a difference between
boreholes is flagged as lateral.

False alarms over 10 simulated counties each (weighted views; the earlier tree model in brackets):
- a correctly coded contact: largest `ratio_best_vs_geologist` median 1.0 (1.0), above 20 in 0 of 10 (1 of 10);
- a single homogeneous unit tested for a split (out-of-sample score, M12): `pct_two_units` median 16 (tree:
  0.02); the two values above 50 (66 and 99) each rest on 1-2 boreholes showing both parts;
- one unit drifting down its cores in random directions: median 0.3, above 50 in 1 of 10 (in-sample score
  without the drift baseline: 100 in 10 of 10);
- two similar units stacked in every borehole (power): 100 in 4 of 6 counties, 71 and 94 in the others;
- one unit given two codes by borehole: `pct_one_unit` median 79 (96), below 10 in 2 of 10 (0 of 10). The
  overall mirrors the views, so mild chance leanings toward "two" in several views combine. Read a merge result
  with its view numbers.
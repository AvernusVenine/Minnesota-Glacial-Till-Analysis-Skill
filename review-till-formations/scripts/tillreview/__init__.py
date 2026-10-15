"""Statistics for the review-till-formations skill. Methods are numbered as in
references/statistical_methods.md (M1-M14)."""

DEFAULTS = {
    "buffer_km": 25.0,            # M5: boreholes this close to the county lend their spread
    "low_count": 100,             # provenance grains below which a count is "low"
    "outlier_p": 0.01,            # M7: withheld if less typical than this
    "simulations": 999,           # M7
    "prior_grains": 1.0,          # M5: weight of the area composition in a unit's proportions
    "shrink_tau": 1.0,            # M5: how far a unit's grain spread may sit from the all-unit value (log scale)
    "shrink_df": 10.0,            # M5: all-unit texture spread counts as this many samples
    "member_exclusion_ratio": 20.0,   # samples this much more likely another formation skip the member loop
    "uncoded_cap": 0.5,           # M14
    "uncoded_pad_ft": 15.0,       # M14
    "split_min_samples": 3,       # M12
    "member_grain_q": 0.05,       # M13
    "seed": 20261012,
}

# Reproduction materials

Six self-contained notebooks reproducing the paper's key figures and
robustness checks:

- **`reproduce_figure2.ipynb`** - Figure 2: minimum land per adult needed
  to reach a given level of mean nutrient adequacy (MAR) across 18
  countries, and how that compares to land currently cultivated by
  surveyed households.
- **`reproduce_figure3_b12.ipynb`** - Figure 3: how much of each country's
  land requirement is spent on vitamin B12 alone, and a leave-one-out
  sweep across every nutrient showing B12 is a genuine outlier.
- **`reproduce_figure4_and_sensitivity.ipynb`** - Figure 4: the adequacy-
  vs-land curve with and without vitamin B12 supplied off-farm, and how
  much animal-product land collapses once it is; plus a threshold-
  robustness check (does the "less land than currently cultivated" finding
  survive moving the adequacy target, and how does B12's land share vary
  with it).
- **`reproduce_eggcap_sensitivity.ipynb`** - supplementary: whether B12's
  dominance is an artifact of the specific egg-production land cap, tested
  by sweeping that cap from 0.25x to 4x baseline.
- **`reproduce_regional_yield_robustness.ipynb`** - robustness check: each
  country's average-yield-optimal product basket, re-evaluated under real
  subnational/yearly yield variation (HarvestStat-Africa data) instead of
  the single country-level average used everywhere else.
- **`reproduce_nar_distribution.ipynb`** - the household nutrient-
  adequacy-ratio (NAR) distribution figure: which nutrients real surveyed
  households tend to fall short on, and by how much.

Each notebook only needs the `.py` files and `data/` in this repository -
no other setup, no notebook depends on another. Outputs are included as
committed, so figures and verification tables are visible without running
anything; running any one of them yourself takes at most a few minutes.

## Interactive dashboard

For interactive exploration of the same underlying model - adjustable
nutrient constraints, country selection, household composition by
age/gender, and a per-country B12 drill-down - see:

**[Dashboard](https://foodsec-dashboard.politerock-163c694a.switzerlandnorth.azurecontainerapps.io)**

## Method, briefly

A linear program allocates land across crop and animal products to
maximize mean adequacy ratio (MAR) - the average, across 17 nutrients, of
how much of each nutrient's dietary requirement is met, capped at 100% per
nutrient - subject to: a caloric intake range, a macronutrient (protein/
carbohydrate) energy-share range, per-nutrient upper intake limits, and
per-country caps on land assignable to each animal product (derived from
survey data on livestock holdings per person). Crop yields are taken from
FAOSTAT; nutrient composition per crop/animal product and dietary
reference intakes are from the sources cited in the notebooks and the
paper's supplementary materials. The B12-counterfactual notebooks
additionally use `functions_supply.py`, which extends the core model with
an "assume this nutrient is supplied off-farm" variant (the nutrient still
counts toward mean adequacy, but its land-linked production requirement is
dropped). The regional-robustness notebook additionally uses subnational
yield data (HarvestStat-Africa) and three helper functions already in
`functions.py` (`preprocess_region_yield`, `create_full_dataframe`,
`fill_yield`) to build a per-region-per-year yield grid, falling back to
national-level FAOSTAT yield wherever no regional match exists.

## Data

- `data/original/FAOSTAT_data_en_3-25-2026.csv` - crop yields by country
  (FAOSTAT).
- `data/original/FINAL_results_10-12-2024.csv`,
  `data/original/CLEAN_nut_requirements_HH_10.12.24.csv`,
  `data/original/hh_full_production_02.08.csv`,
  `data/indicator_data.csv` - household-level survey data (production,
  nutrient requirements, land cultivated). Two columns were blanked before
  inclusion here - a GPS-coordinate column and internal survey
  dataset-linkage IDs - neither is used by any analysis here; the raw
  files simply shouldn't carry real respondents' location data.
- `data/processed/hvstat_crops_matched_long.csv` - subnational crop yields
  by region and year (HarvestStat-Africa), used only by the regional-
  robustness notebook.
- `data/crop_composition_table.csv`,
  `data/animal_composition_table_with_carrying_capacity.xlsx`,
  `data/food_comp_DRI.xlsx` - nutrient composition and dietary reference
  intake reference tables.
- `data/processed/country_mar_wlc.csv` - the paper's own precomputed
  MAR-vs-land curves, used in `reproduce_figure2.ipynb` to verify that
  reproduction matches the published results exactly, not just
  approximately.
- `data/processed/all_heads_per_person_95percentiles.csv`,
  `data/processed/crop_yield_country.csv` - supporting reference data
  (country-level aggregates, not household-level).

## Requirements

```
pip install -r requirements.txt
```

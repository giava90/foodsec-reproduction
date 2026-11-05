# Task: climate-change robustness check of the land-sparing / nutrient-adequacy model

You are working autonomously overnight. Nobody will answer questions. Make
reasonable decisions, document every one in `climate_robustness/DECISIONS.md`,
and keep going.

## Context

This repository reproduces a paper ("A single micronutrient sets the limit of
nutrition-sensitive land sparing across sub-Saharan Africa"). A linear program
(`functions.py`, scipy `linprog`) allocates land across crop and animal
products to maximise mean adequacy ratio (MAR, 17 nutrients) for an adult,
and finds the **minimum land per adult needed to reach MAR = 1** (and other
MAR levels) per country, using FAOSTAT average crop yields (`data/original/
FAOSTAT_data_en_3-25-2026.csv` -> `df_crop_yield`). Read `README.md` and the
notebooks first. `reproduce_regional_yield_robustness.ipynb` is the closest
template: it fixes the average-yield optimal basket and re-evaluates it under
other yields. `reproduce_figure2.ipynb` shows the min-land-vs-MAR curves.

## Research question (from the PI)

1. **Robustness of today's solution under future yields.** Take each
   country's optimal land allocation computed with current (FAOSTAT) yields
   (land per product, hectares). Keep that allocation fixed and evaluate it
   under projected future yields: what MAR does it reach, which nutrients
   fall short, how much worse does it get? Equivalent view: how much extra
   land would the current basket composition need to stay at MAR = 1.
2. **Re-optimised solution under future yields (safety check).** Re-solve the
   LP with future yields: minimum land to reach MAR = 1 (and the MAR-vs-land
   curve as in Figure 2). Compare with today: change in minimum land, change
   in the basket composition (which products enter or leave, land shares),
   whether the paper's headline findings still hold (vitamin B12 being the
   binding nutrient that sets the land limit; minimum land below currently
   cultivated land - see Figure 2/4 notebooks).
3. **Where improvement is possible.** Contrast (1) and (2): how much of the
   climate damage can be recovered just by changing the crop mix? Which
   countries/products drive that? Is there a "climate-robust" basket (e.g.
   the solution optimised under the worst scenario, evaluated under current
   and all future scenarios), and what does robustness cost in land?

## Data provided (in `data/climate/`)

**Timing note:** `isimip3b_country_annual.csv.gz` is still being produced on
a local machine and is pushed automatically to this branch when ready
(expected within ~30 min of your start). Start with everything that does not
need it (read the repo, reproduce/validate the baseline, crop mapping, the
Joshi cross-check data). Run `git pull` every ~10 minutes until the file
appears. If it is still missing after 60 minutes, produce it yourself with
`python data/climate/process_isimip.py data/climate` (needs `xarray netCDF4
requests`; it reads `data/climate/yield_datasets.json` and needs
`countrymasks-fractional_30arcmin.nc` and
`landuse-15crops_2015soc_annual_2015_2100.nc` downloaded from
https://files.isimip.org/ISIMIP3b/InputData/geo_conditions/countrymasks/ and
https://files.isimip.org/ISIMIP3b/InputData/socioeconomic/landuse/2015soc/
into `data/climate/`; output lands in `data/climate/out/`; do NOT commit the
.nc files). If ISIMIP is unreachable from your sandbox, do the analysis with
the Joshi literature data only and say so prominently.

`isimip3b_country_annual.csv.gz` - prepared from ISIMIP3b / GGCMI phase 3
global gridded crop model outputs (`yield-<crop>-noirr`, rainfed, 2015soc,
2015-2100, 0.5 deg), aggregated to the 19 countries of this study as a
harvested-area-weighted mean (weights: ISIMIP 2015soc landuse-15crops,
rainfed+irrigated area, times the ISIMIP fractional country mask).
Columns:

- `model` - crop model (up to 13 GGCMs: acea, crover, cygma1p74, dssat-pythia,
  epic-iiasa, isam, ldndc, lpj-guess, lpjml, pdssat, pepic, promet,
  simplace-lintul5; coverage differs by crop)
- `gcm` - climate model (gfdl-esm4, ipsl-cm6a-lr, mpi-esm1-2-hr, mri-esm2-0,
  ukesm1-0-ll)
- `scen` - ssp126 or ssp585
- `co2` - `default` (transient CO2, includes CO2 fertilisation) or `2015co2`
  (CO2 held at 2015 level, i.e. no fertilisation) - not all models ran both
- `crop` - mai (maize), ri1 (rice), soy (soybean), swh (spring wheat), mil
  (millet), sor (sorghum), bea (beans/pulses), pot (potato), cas (cassava)
- `iso3`, `year`
- `yield_areaw` - area-weighted country mean yield, t/ha dry matter
- `yield_unw` - unweighted mean over all country cells (fallback when the
  crop area weight is 0 / missing)
- `crop_area_ha_weight` - sum of weights (0 means the crop is essentially not
  grown there in the land-use data; treat those changes as low-confidence)

`PROVENANCE.md` has the exact processing. The processing script is
`data/climate/process_isimip.py`.

`joshi2026/yield_prediction_data.csv` (+ `metadata.xlsx` column dictionary,
`references.xlsx`) - independent literature meta-analysis "Projected Climate
Impacts on Crop Yields in Africa" (Harvard Dataverse doi:10.7910/DVN/SELOFD,
CC BY 4.0): 2378 published projections, column `Climate impacts (%)`, with
crop, country (some regional labels like "East Africa", misspellings like
"Cote D'lvoire", "Eygpt"), time slice / future mid-point, scenario (mostly
CMIP3/5 RCPs), CO2 and adaptation flags. Use it as (a) a **cross-check** of
the ISIMIP change factors for the overlapping crops/countries (compare the
no-adaptation, mid-century entries; report agreement/disagreement), and (b)
an **alternative factor source** for crops ISIMIP lacks (yam, sweet potato,
tomato, sunflower, sugarcane, coffee, barley...) - e.g. a sensitivity run where
these crops take the literature median change (by crop, pooled over
countries if a country has too few points) instead of the proxy. Filter to
non-adapted entries (`Adaptation` == "No") for the main comparison.

**Use relative changes, not absolute yields.** Crop models are not calibrated
to national yields. Build a change factor per (model, gcm, scen, co2, crop,
country) as mean(future window) / mean(baseline window), with baseline =
2015-2034 (closest to the FAOSTAT period used by the model; note in the
report that this baseline already contains some warming) and future windows
2041-2060 (headline, "2050") and 2081-2100 ("2090"). Multiply FAOSTAT yields
by the factor. Summarise across the ensemble (median, and an
inter-quartile/likely range, e.g. 25th/75th percentile across model x GCM),
and keep the worst plausible case (e.g. 10th percentile). Guard against
outliers (tiny baselines -> huge ratios): clip or drop factors where the
baseline mean is near zero, and document it.

You decide which scenario is the headline after looking at the data (the PI
explicitly left that to you). A sensible default: 2050 SSP5-8.5 and SSP1-2.6,
with and without CO2 fertilisation, ensemble median as headline, 25th/75th
and 10th percentile as uncertainty / worst case.

## Crops without a projection

The model has ~70 crops; ISIMIP covers 9. Map each model crop to the closest
ISIMIP crop and flag the mapping, e.g.: cereals (maize->mai, rice->ri1,
wheat/barley/oats/teff->swh, millet/pearl_millet/fingermillet/fonio->mil,
sorghum->sor), pulses/legumes (beans, cowpea, pigeon_pea, groundnut,
bambara_nut, chick_peas, lentils, peas...)->bea (groundnut/soya_bean->soy is
also defensible), roots/tubers (cassava->cas, irish_potato->pot, sweet_potato,
yam, taro, arrow_root->cas), everything else (vegetables, fruits, oil/cash
crops) -> the ensemble mean change across all ISIMIP crops for that country
("proxy"). Write the full mapping table to `climate_robustness/
crop_mapping.csv` with a `mapping_type` column (direct / proxy-similar /
proxy-mean). Run a **sensitivity** where all unmapped (proxy-mean) crops keep
0% change, and report whether conclusions change.

Animal products: the model's animal land productivity is fixed. Keep it
fixed in the main analysis (state this clearly as a limitation); optionally a
sensitivity where animal productivity scales with the mean crop/forage change.

## Deliverables (all under `climate_robustness/`)

1. `climate_yield_factors.csv` - tidy table of change factors (country,
   model crop, isimip crop, mapping_type, scenario, co2, window, median, p10,
   p25, p75, p90, n_members).
2. `reproduce_climate_robustness.ipynb` - a self-contained, executed
   notebook in the same style as the other `reproduce_*.ipynb` notebooks
   (markdown explaining each step, outputs committed, runs from the repo root
   with only repo files). Any heavy helper code can go in
   `climate_robustness/climate_helpers.py`.
3. Figures (PNG, also embedded in the notebook), at least:
   - map/heatmap of yield change factors by country x crop group (headline
     scenario);
   - per-country MAR achieved by today's fixed optimal allocation under
     future yields (distribution across ensemble members / scenarios, like
     the regional robustness boxplot), plus which nutrients become limiting;
   - minimum land to reach MAR = 1: today vs re-optimised future (with
     uncertainty), and the fixed-basket land needed to restore MAR = 1;
   - basket composition change (today vs future-optimal), and the
     recoverable fraction of climate damage by re-optimising.
4. `results_summary.csv` - one row per country x scenario with the key
   numbers.
5. `REPORT.md` - a concise scientific write-up (aim for something the PI can
   drop into a supplementary note): question, data and methods (including
   every assumption and limitation: rainfed only, 2015soc fixed harvest
   areas, uncalibrated models -> relative changes only, crop mapping proxies,
   animals fixed, no adaptation besides crop mix, baseline 2015-2034),
   results with numbers, and a clear answer to: *are the paper's headline
   conclusions robust to climate change?* Also include a short section "What
   I did" listing your steps, and a "Caveats / things the PI should check".
6. `REPORT.html` - the same report as a single self-contained HTML page
   (figures embedded as base64 PNG, inline CSS, no external resources,
   readable on a phone). This will be published as a web page.
7. Update the top-level `README.md` with a short entry for the new notebook
   (one bullet in the notebook list, consistent with the others).
8. Last step: create `climate_robustness/DONE` containing a 5-10 line plain
   text summary of the main findings, commit, and push. A watcher is waiting
   for this file.

## Rules

- **Anonymity is mandatory.** This repository is shown to double-blind
  reviewers. Before ANY commit run
  `git config user.name "anonymous" && git config user.email "anonymous@example.com"`.
  Never write person names, institutions, cities, e-mail addresses, or user
  names into files. Before every push, grep the diff for
  names/institutions/emails/local paths and fix any hit.
- Work on branch `climate-robustness` (already checked out / exists on the
  remote). Never push to `main`. If pushing to `climate-robustness` is
  rejected, push to `claude/climate-robustness` instead and write that in
  `DONE`.
- Commit and push **incrementally** (after the factor table, after the
  robustness evaluation, after re-optimisation, after the report) with
  descriptive messages, so partial work survives if you run out of time.
- Do not modify existing notebooks, `functions.py`, `functions_supply.py`, or
  data files outside `data/climate/` and `climate_robustness/` (except the one
  README bullet). If you need a variant of a model function, put it in
  `climate_robustness/climate_helpers.py`.
- Validate first: before the climate runs, reproduce today's baseline
  minimum-land-to-MAR=1 per country and check it against
  `data/processed/country_mar_wlc.csv` / the Figure 2 notebook output. If it
  doesn't match, fix that first (do not proceed on a broken baseline).
- Prefer cheap, vectorised runs. If an LP fails, use the same
  `failure="remove_constraints"` fallback and count failures; do not silently
  drop.
- Be honest: if something didn't work or a result is fragile, say so in the
  report prominently.
- **Time box: about 3.5 hours total.** At ~2h45m stop starting new analyses
  and finish the report, HTML, README bullet, DONE file, and push. A
  complete, honest, smaller analysis beats an unfinished big one.

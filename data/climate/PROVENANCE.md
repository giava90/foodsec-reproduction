# Provenance of `data/climate/`

`joshi2026/` - unmodified copy (files renamed only) of the Harvard Dataverse
dataset "Projected Climate Impacts on Crop Yields in Africa",
doi:10.7910/DVN/SELOFD, CC BY 4.0: a compilation of published crop-yield
projections for Africa.

# Provenance: `isimip3b_country_annual.csv.gz`

**Source.** ISIMIP3b agriculture sector output (GGCMI phase 3), accessed via
https://data.isimip.org and https://files.isimip.org. Variables
`yield-<crop>-noirr` (rainfed yield, t ha-1 yr-1 dry matter), files
`<model>_<gcm>_w5e5_<ssp>_2015soc_<sens>_yield-<crop>-noirr_global_annual-gs_2015_2100.nc`.
Licence: CC BY 4.0 (ISIMIP terms of use). Key references: Jägermeyr et al.
2021, *Nature Food* 2:873-885 (GGCMI phase 3 / CMIP6 crop projections);
Frieler et al. 2024, *GMD* (ISIMIP3 protocol).

**Selection.** All ISIMIP3b agriculture datasets with
`climate_scenario` in {ssp126, ssp585}, `irrigation` = noirr,
`crop` in {mai, ri1, soy, swh, mil, sor, bea, pot, cas}, both
`sens_scenario` = default (transient CO2) and 2015co2 (CO2 fixed at 2015
level), all crop models and all five bias-adjusted CMIP6 GCMs available
(909 datasets). Excluded: irrigated runs (irrigation is a small share of
cropland in the study countries), second-season rice (ri2), winter wheat
(wwh), historical runs (baseline is taken from 2015-2034 of the same run).

**Aggregation.** For each file and each of the 19 study countries:

```
yield_areaw[year] = sum_cells( y * A_crop * f_country ) / sum_cells( A_crop * f_country )
```

over cells with a valid yield, where `A_crop` = 2015 harvested area
(rainfed + irrigated) of the matching crop type in ISIMIP
`landuse-15crops_2015soc_annual_2015_2100.nc` (maize; rice; oil_crops_soybean;
temperate_cereals for wheat; tropical_cereals for millet and sorghum; pulses
for beans; temperate_roots for potato; tropical_roots for cassava) and
`f_country` = ISIMIP `countrymasks-fractional_30arcmin.nc`. `yield_unw` uses
`f_country` alone as the weight. Fill values (>1e19) and negatives were set
to missing. Script: `process_isimip.py` (this folder).

**Countries (ISO3).** BFA Burkina Faso, BDI Burundi, COM Comoros, CIV Côte
d'Ivoire, COD Democratic Republic of the Congo, ETH Ethiopia, GHA Ghana, KEN
Kenya, MWI Malawi, MLI Mali, NER Niger, NGA Nigeria, RWA Rwanda, SEN Senegal,
SLE Sierra Leone, ZAF South Africa, UGA Uganda, TZA United Republic of
Tanzania, ZMB Zambia.

**Caveats.** Absolute yields from global gridded crop models are not
calibrated to national statistics; use relative changes only. Comoros is a
handful of 0.5-degree cells. Some crops have near-zero area weight in some
countries (see `crop_area_ha_weight`).

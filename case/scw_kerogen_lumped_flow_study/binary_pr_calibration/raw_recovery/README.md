# Raw experimental table recovery

This directory separates **table recovery evidence** from rows that are admissible for EOS regression.

## Rule

A paper abstract, plot, fitted BIP, dataset count, or min/max range is not converted into a synthetic experimental row. Only pointwise values traceable to an original table (or an explicitly labelled digitization when the project intentionally permits that lower evidence level) may enter the regression-observation file.

## H2O + n-hexane

Primary source: Tian, Michelberger and Franck (1991), DOI `10.1016/S0021-9614(05)80063-1`.

The original paper reports five n-hexane mixture compositions and phase-boundary/PVT data over roughly 550–700 K and 20–240 MPa. A GERG-2004 data audit independently confirms that the Yiling et al. rows belong to the `n-C6H14 + H2O` system and gives the per-composition dataset counts/ranges preserved in `nhexane_tian1991_dataset_index.csv`.

The accessible audit identifies approximately `x_H2O = 0.39, 0.50, 0.60, 0.70, 0.80`; the article abstract rounds the first value when describing the five compositions as 0.40–0.80. The pointwise `(T,p,V_m)` rows are **not** reconstructed from these ranges.

Status: `RAW_POINT_ROWS_NOT_YET_RECOVERED`.

## H2O + dodecane

Primary source: Stevenson, LaBracio, Beaton and Thies (1994), DOI `10.1016/0378-3812(94)87016-0`.

The original article measured VLE/LLE compositions and critical phenomena from about 600–660 K and to about 31 MPa. A later correlation paper explicitly reports the Stevenson dataset counts used at the two temperatures:

- 603.6 K VLE: 18 points;
- 603.6 K LLE: 14 points;
- 633.0 K VLE: 14 points.

Those counts are preserved in `dodecane_stevenson1994_dataset_index.csv`. They are useful for checking whether a future transcription is complete, but they are not themselves experimental rows.

Status: `RAW_POINT_ROWS_NOT_YET_RECOVERED`.

## Existing exact rows

The repository already contains exact Stevenson water–squalane coexistence rows in `../h2o_squalane_reference.csv`; these remain the only target-window pointwise coexistence data currently admitted to the generic regression harness.

## Promotion rule

When an original n-hexane or dodecane table is recovered, every row must record the source table/page and original units, and the number of imported rows must match the recovered table count before the system can move from `DATA_BLOCKED` to `FIT_READY`.

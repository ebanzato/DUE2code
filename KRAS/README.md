# KRAS single-cell analysis with DUE2


The analysis uses the R package [DUE2](https://github.com/ebanzato/DUE2) to test for
differences in mean expression and associations among the genes of the KEGG
*Ras signaling pathway*, across the time points of KRAS G12C inhibition.


## Content

| File | Description |
|---|---|
| `KRAS.R` | Loads the raw counts, filters cells and genes, computes size factors, maps genes to ENTREZ IDs, selects the KEGG pathways containing KRAS, and saves the resulting dataset. |
| `H0_kras.R` | Simulation under the null hypothesis: group labels are randomly permuted at t0, and the empirical FDR of the negative binomial and Poisson versions of DUE2 is compared. |
| `KRAS_TS.R` | Analysis of the time sequence (t0 vs t4, t4 vs t24, t24 vs t72) with DUE2 and tables of the significant differences. |
| `names_converter_norm_finale.csv` | Ensembl to ENTREZ ID / HGNC symbol conversion table, obtained from biomaRt. |


## Data

The data are publicly available on GEO under accession
[GSE137912](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137912)
(cell lines H358, H2122 and SW1573 at 0, 4, 24 and 72 hours).

Download the following files into a single folder:

- `GSE137912_cell.annotation.csv.gz`
- for each sample listed below, the three files `*.matrix.mtx.gz`, `*.barcodes.tsv.gz` and `*.genes.tsv.gz`:


Copy `names_converter_norm_finale.csv` from this repository into the same folder.
Alternatively, the table can be regenerated with biomaRt by uncommenting the
corresponding block in `KRAS.R` (results may differ slightly depending on the
Ensembl release).



## How to run

Each script has a `Paths` section at the top: set the paths there before running it.

1. **`KRAS.R`**: set `path` (folder with the data) and `out_file` (where to save
   the processed dataset, e.g. `KRAS_dataset.RData`).
2. **`H0_kras.R`**: set `data_file` (the `out_file` of step 1) and `res_file`.
3. **`KRAS_TS.R`**: set `data_file` (the `out_file` of step 1) and `res_file`.

Steps 2 and 3 are independent and both require only the output of step 1.
Once the results have been saved, the analysis in steps 2 and 3 can be skipped by
uncommenting the `load(res_file)` line.





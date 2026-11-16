# skNIR Haploid Kernel Sorting in Field and Sweet Corn

Data and analysis code for:

> Sharma S, Gustin JL, Frei UK, Settles AM, Lübberstedt T, Resende MFR, Hershberger J (2026). Single-kernel near-infrared spectroscopy enables haploid kernel sorting in field and sweet corn using high-oil haploid inducers across diverse donor-inducer combinations. *Theoretical and Applied Genetics* 139:332. https://doi.org/10.1007/s00122-026-05424-1

- **Paper:** https://doi.org/10.1007/s00122-026-05424-1
- **Preprint:** https://doi.org/10.64898/2026.07.23.740370
- **Data (Dryad):** https://doi.org/10.5061/dryad.cfxpnvxp1

If you use this code or data, please cite the paper above.

## Overview

Rapid and accurate haploid kernel identification remains a major bottleneck in maize doubled haploid production, particularly in sweet corn germplasm where the widely used R1-nj kernel color marker is frequently suppressed, resulting in unreliable visual haploid sorting. To address this challenge, we evaluated an automated single-kernel near-infrared (skNIR) spectroscopy-based sorting system across 12 induction populations derived from two field corn and two sweet corn donor backgrounds and four high-oil haploid inducer lines. We compared two skNIR-based haploid sorting approaches: predicted oil content-based methods and spectra-based multivariate classification methods.

## Repository contents

| Path | Contents |
|---|---|
| `Data_skNIR_Manuscript_June27.2026.csv` | Kernel-level data (see variables below) |
| `Final R Codes_June28.2026/skNIR_RScript1.R` | Oil-based sorting methods and multivariate classification (PLS-DA, RF, SVM); generates Figure 3 and Tables 3–5, Supplementary Tables 3–5 |
| `Final R Codes_June28.2026/skNIR_RScript2.R` | All remaining analyses; generates Figures 2, 4, 5, Supplementary Figures 1–2, Tables 1–2, and Supplementary Tables 1–2 |
| `Final_Results_csvfiles_skNIR_June28.2026/Results_PART1_*` | CSV outputs of `skNIR_RScript1.R` |
| `Final_Results_csvfiles_skNIR_June28.2026/Results_PART2_*` | CSV outputs of `skNIR_RScript2.R` |
| `High-res_figures_June22.2026/` | Manuscript figures (Figures 1–5, Supplementary Figures 1–3) |

## Data dictionary

`Data_skNIR_Manuscript_June27.2026.csv`. Each row is one kernel. Missing values are coded as `NA`.

| Variable | Description |
|---|---|
| `Index` | Unique identifier for each kernel sample |
| `source` | Data source; all samples were collected at Iowa State University (ISU) |
| `population` | Population ID (1–12): induction crosses of four donor (female) parents with four inducer (male) parents |
| `plate` | ID of the 24-well microtiter plate used to track kernel samples |
| `vis.score` | Visual classification by experienced technical staff: 1 = haploid, 4 = hybrid, 5 = defective |
| `dna.score.final` | True kernel class from genomic DNA markers: 1 = haploid, 4 = hybrid, 5 = visually defective (same as `vis.score`). Used as the reference standard for evaluating all classification methods |
| `Donor` | Female parent |
| `Inducer` | Male (haploid inducer) parent |
| `pred.oil` | Kernel oil content (%) predicted by the skNIR sorter; used for the two oil-based classification methods |
| `nmr.oil` | Kernel oil content (%) measured by nuclear magnetic resonance (NMR); used to validate `pred.oil` |
| `X940`–`X1640` | Absorbance at wavelengths 940–1640 nm; used for the three multivariate classification methods |

## Running the code

Analyses were run in R v4.4.0. Install the required packages:

```r
install.packages(c("dplyr", "tidyr", "waves", "caret", "ggplot2", "ggtext",
                   "emmeans", "car", "FSA", "multcompView", "rstatix",
                   "cowplot", "ggrepel"))
```

**Before running**, update the `read.csv()` path near the top of each script (line 35 in `skNIR_RScript1.R`, line 31 in `skNIR_RScript2.R`) to point to `Data_skNIR_Manuscript_June27.2026.csv` on your machine. Output CSVs and figures are written to the working directory.

## Contact

- Project lead: Suman Sharma
- PI: Jenna Hershberger ([Hershberger Lab](https://hershlab.org))

Project start: August 2024

Project PI: Dr. Jenna Hershberger

Project Lead: Suman Sharma

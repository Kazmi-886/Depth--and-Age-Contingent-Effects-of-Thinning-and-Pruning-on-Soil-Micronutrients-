# Soil micronutrient depth profiles under thinning and pruning in *Larix principis-rupprechtii* plantations (Saihanba, China)

> **Repository purpose:** reproducible analyses and supporting materials for a manuscript examining how thinning and pruning interact with stand development and soil depth to structure Fe–Cu–Mn–Zn profiles in cold-temperate larch plantations.

## Overview
Micronutrient distributions in plantation soils are governed by coupled controls of organic-matter turnover, redox-sensitive transformations, and vertical transport. However, the extent to which routine silvicultural operations (notably thinning and pruning) modify *depth-resolved* micronutrient profiles remains insufficiently resolved for cold-temperate conifer plantations. This study evaluates soil **Fe, Cu, Mn, and Zn** across a factorial thinning × pruning experiment in *Larix principis-rupprechtii* plantations at Saihanba Forest Farm (Hebei Province, northern China), explicitly accounting for **stand age** and **soil depth**.

## Study system and design
- **Site:** Mandian area, Saihanba Mechanical Forest Farm, Hebei, China (42°02′–42°36′ N, 116°51′–117°39′ E; 1650–1830 m a.s.l.)
- **Stand ages:** 22, 26, and 42 years
- **Factorial treatments:** 5 thinning levels (T0–T4) × 3 pruning levels (P0–P2) per age class
- **Soil depths:** 0–10 cm, 10–20 cm, and 20–40 cm
- **Analytical unit:** plot × depth (after laboratory replicate aggregation)
- **Total observations:** 135 (3 ages × 15 TP combinations × 3 depths)

## Main findings (manuscript summary)
- **Stand age and soil depth were the dominant controls** on micronutrient concentrations, with management effects generally smaller and context dependent.
- **Consistent surface enrichment** was observed for **Fe, Mn, and Zn**, with topsoil typically higher than deeper layers.
- **Cu showed an age-dependent vertical reorganisation**, including an inversion in older stands where deeper layers exceeded surface concentrations.
- **Management signals were conditional**, emerging in specific age–depth domains (e.g., pruning-associated shifts in Cu and Mn profiles; age-dependent sensitivity of Fe and Cu to canopy interventions).
- **Multivariate profiling (Fe–Cu–Mn–Zn) indicated weak treatment separation** once stand age and depth were accounted for, suggesting limited evidence for a strong “treatment fingerprint” in bulk micronutrient pools.

## Methods summary
### Laboratory determination
- Air-dried, sieved (≤2 mm) mineral soil; concentrations reported on a dry-mass basis.
- Acid digestion followed by **flame atomic absorption spectrophotometry (FAAS)** for Fe, Cu, Mn, and Zn (with standard QA/QC elements: blanks, calibration checks, and batch controls).

### Statistical analysis (reproducible workflow in this repository)
- Element-wise modelling on the **log scale** to address right-skewness and enable multiplicative interpretation.
- **Mixed-effects models (LMMs)** with a unit-level random intercept for repeated depth observations where supported; **linear-model fallback** where the random effect was not supported.
- Planned **depth contrasts** (0–10 vs 10–20; 0–10 vs 20–40) with **BH–FDR** multiplicity control.
- Moderation analyses using intensity-based thinning/pruning parameterisation (+10% units) to test **Age × management** and **Depth × management** structure.
- Multivariate analysis of joint Fe–Cu–Mn–Zn profiles using **PERMANOVA** and conditional ordination (partial dbRDA) to evaluate whether thinning–pruning combinations produce coordinated multielement shifts beyond age–depth structuring.

## Repository structure (recommended)















## Reproducibility
### Option A (R; recommended if your analysis is in R)
1. Install R (≥ 4.2) and RStudio.
2. Restore packages (recommended):
   - Use `renv` for reproducibility:
     ```r
     install.packages("renv")
     renv::restore()
     ```
3. Run scripts in order:
   - `scripts/01_data_qc.R`
   - `scripts/02_models_univariate.R`
   - `scripts/03_models_multivariate.R`
   - `scripts/04_figures.R`
   - `scripts/05_tables.R`

### Key R packages (indicative)
`lme4`, `lmerTest`, `emmeans`, `performance` (R²), `vegan` (PERMANOVA/dbRDA), `ggplot2`, `dplyr`, `readr`

> Update this list to match your actual session info and lockfile.

## Data availability
Data sharing may be constrained by site agreements or institutional policies. If raw data cannot be public, this repository can still provide:
- a de-identified/aggregated dataset sufficient to reproduce models, **or**
- a synthetic example dataset plus complete scripts, **or**
- scripts that run on data available upon reasonable request.

(Replace this section with your journal-compliant statement.)

## How to cite
Add the final citation after acceptance. Suggested placeholder:
> Author(s). Year. *Title*. Journal (in review / accepted). DOI: xxx

Include a `CITATION.cff` file to enable GitHub’s citation button.

## License
Choose a license appropriate for your target journal and data policy (e.g., MIT for code; CC BY 4.0 for documentation). If data are restricted, clearly separate licensing for code vs. data.

## Contact
For questions, issues, or requests related to data access, please open an issue or contact the corresponding author.















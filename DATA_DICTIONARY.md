# Data dictionary and coding conventions

## Common abbreviations

| Code | Definition |
|---|---|
| DTF | Days to flowering |
| DTM | Days to maturity |
| HSW | Hundred-seed weight |
| PH | Plant height |
| PPP | Pods per plant |
| YPP / YPPlnt | Yield per plant; `YPPlnt` is the internal analysis label and `YPP` is the manuscript label |
| PA | Predictive ability |
| LD | Linkage disequilibrium |
| GWAS | Genome-wide association study |
| BRR | Bayesian ridge regression |
| GxE | Genotype-by-environment interaction |

## Cross-validation schemes

| Code | Target scenario |
|---|---|
| CV1 | Untested genotypes in observed environments |
| CV2 | Sparse testing: tested genotypes missing records in observed environments |
| CV0 | Tested genotypes in untested environments |
| CV00 | Untested genotypes in untested environments |

## Models

| Model | Linear predictor |
|---|---|
| M1 | E + L |
| M2 | E + L + G |
| M3 | E + L + G + GxE |
| M4 | E + L + G + GxE + S |

`E` is the environment incidence term, `L` is the line incidence term, `G` is
the genomic main-effect kernel, `GxE` is the genomic interaction kernel, and
`S` is a BRR-shrunk regression term constructed from fold-specific meta-GWAS
ranked markers. The primary M4 analysis uses 500 markers; 100 and 1,000 are
sensitivity analyses.

## Marker panels

| Internal code | Manuscript label | Definition |
|---|---|---|
| `chr` / `genome-wide` | Genome-wide | 409,128 SNPs after final QC |
| `tags` / `haplotype-tagged` | Haplotype-tagged | 9,845 LD-based tag SNPs |

For tagged-panel M4, a selected genome-wide SNP is retained if it is a tag.
Otherwise, the strongest tag proxy within 500 kb is used only if pairwise
`r2 >= 0.50`; duplicate proxies are collapsed.

## Fold-level predictive-ability files

Files named `PA_fold_level_M1_M4_top{K}.csv` contain:

| Column | Definition |
|---|---|
| Model | M1, M2, M3, or M4 |
| Panel | `genome-wide` or `haplotype-tagged` |
| Trait | Trait code |
| CV | Cross-validation scheme |
| Repetition | Replicate partition, 1-10 |
| Fold | Fold within repetition, 1-5 |
| PA | Mean of eligible within-environment Pearson correlations in the fold |

An environment contributes to PA only when it has at least five non-missing
validation observations. Eligible environment correlations are averaged
equally to obtain one PA value per fold.

## Fold-specific meta-GWAS ranking

| Column | Definition |
|---|---|
| SNP | Marker identifier |
| Z_meta | Signed, sample-size-weighted Stouffer statistic |
| P_meta | Two-sided probability derived from `Z_meta`, used for ranking |
| N_environments | Number of training environments contributing evidence |
| min_P | Smallest environment-specific GWAS P-value |
| consistent_sign | Whether all available non-zero effect estimates have the same sign |
| FDR_meta | Benjamini-Hochberg-adjusted `P_meta`, reported descriptively |
| Rank | Fold-specific meta-GWAS rank |

Because correlations among environment-specific association statistics were
not explicitly modelled, combined probabilities are used to rank markers and
not as formal inferential P-values.

## Prediction files

| Column | Definition |
|---|---|
| testing | Row index in the analysis phenotype table |
| Individual | Genotype identifier |
| environment | Site-year environment identifier |
| y | Observed phenotype/BLUE |
| yHat | Predicted value |

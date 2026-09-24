# Genotype data supplied with the revised manuscript

This archive contains the two processed marker panels used for genomic prediction. It replaces the previously supplied intermediate 596,202-marker matrix, which was not the analysis dataset described in the manuscript.

## Files

| File | Description |
|---|---|
| `genome_wide_genotypes_190x409128.rds` | Mean-imputed additive-dosage matrix for 190 chickpea accessions and 409,128 SNPs. The first column, `Taxa`, contains accession identifiers; remaining columns are SNPs coded as allele dosages (0, 1 or 2 before imputation, with marker means used for missing calls). |
| `genome_wide_marker_map.tsv` | Complete map for the genome-wide panel: SNP identifier, chromosome, physical position (bp), PLINK A1 allele and PLINK A2 allele. |
| `genome_wide_marker_list.txt` | Ordered list of the 409,128 genome-wide SNP identifiers. |
| `haplotype_tagged_genotypes_190x9845.rds` | Mean-imputed additive-dosage matrix for the same 190 accessions and 9,845 LD-based tag SNPs, in the same `Taxa`-first format. |
| `haplotype_tagged_marker_map.tsv` | Complete map for the tagged panel: SNP identifier, chromosome, physical position (bp), PLINK A1 allele and PLINK A2 allele. |
| `haplotype_tagged_marker_list.txt` | Ordered list of the 9,845 tag-SNP identifiers. |
| `panel_validation_summary.tsv` | Panel dimensions, chromosome labels and number of missing values after imputation. |
| `checksums_md5.tsv` | MD5 checksums for archive-integrity verification. |

RDS files can be read in R using `readRDS()`. The map and summary files are tab-delimited UTF-8 text files. Marker-map rows and genotype-matrix marker columns use the same SNP identifiers.

## Genome-wide panel processing

The source VCF contained 3,941,492 variants. Variants were restricted to the eight assembled chickpea chromosomes (Ca1-Ca8), leaving 2,470,880 SNPs before sample- and marker-level filtering. PLINK v1.90b6.21 was used with `--mind 0.2`, `--geno 0.1`, and `--maf 0.05`. Nineteen accessions exceeding the sample-missingness threshold were removed. The resulting panel contained 190 accessions and 409,128 SNPs. Missing genotype dosages in the retained panel were imputed to the corresponding marker mean before genomic prediction and GWAS input construction.

Genome-wide SNP counts by chromosome are:

| Chromosome | SNPs |
|---|---:|
| Ca1 | 60,702 |
| Ca2 | 36,667 |
| Ca3 | 27,699 |
| Ca4 | 127,525 |
| Ca5 | 30,504 |
| Ca6 | 58,582 |
| Ca7 | 48,216 |
| Ca8 | 19,233 |
| **Total** | **409,128** |

No unplaced-contig or scaffold markers are included in the supplied processed panel.

## Haplotype-tagged panel processing

The 9,845-marker panel was derived from the same quality-controlled Ca1-Ca8 source panel and the same 190 accessions. Tag SNPs were selected using the LD-based two-pass procedure described in the Methods, with the stated pairwise LD criterion. The final tag list was extracted from the processed genome-wide panel. Missing tag-SNP dosages were imputed to the corresponding marker mean before construction of the genomic relationship matrix.

## Genotype coding and imputation

PLINK additive dosages record the number of copies of the allele represented in each exported marker column. Because mean imputation can produce non-integer values, the supplied RDS matrices may contain continuous dosage values at originally missing calls. Imputation used genotype data only and did not use phenotypes or cross-validation outcomes.

The `A1` and `A2` columns in the marker maps reproduce the allele fields from the corresponding PLINK BIM files. They should not be interpreted as ancestral and derived alleles.

## Reproducibility checks

After extraction and imputation:

- both matrices must contain the same 190 accession identifiers;
- the genome-wide matrix must contain 409,128 marker columns;
- the tagged matrix must contain 9,845 marker columns;
- all markers must map to Ca1-Ca8;
- neither supplied matrix should contain missing values; and
- marker identifiers must be unique within each panel.

The archive can be regenerated from the verified analysis files by running `Rscript prepare_submission_genotype_archive.R` from the `haplotype.analysis` directory.

## Obsolete file warning

The historical project object `G_haploview_tagged.rds` contains 237,483 markers and does **not** represent the final 9,845-marker panel. It must not be included in the manuscript data archive or described as the haplotype-tagged prediction matrix.

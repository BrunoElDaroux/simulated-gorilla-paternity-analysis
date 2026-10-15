# Simulated Mountain Gorilla Paternity, Kinship & Population Genetics

**Computational genetics simulation project | Python · Jupyter · Synthetic genetic data**

![Python](https://img.shields.io/badge/Python-3.10%2B-blue) ![Jupyter](https://img.shields.io/badge/Jupyter-Notebooks-orange) ![Data](https://img.shields.io/badge/Data-100%25%20synthetic-purple)

## ⚠️ Data Transparency

The datasets used in this project is **synthetic** and was created strictly for **educational, analytical, and portfolio purposes**. It does not represent real gorilla genetic samples, genotypes, pedigrees, or verified paternity records.

The genetic structure, biological context, and assumptions were informed by published research on mountain gorilla population genetics and social behavior. However, **all paternity assignments, kinship estimates, reproductive-skew measurements, and other findings presented in this repository are generated from synthetic data and should not be interpreted as real-world findings**.

This project demonstrates how genetic data could be analyzed to investigate paternity, kinship, and reproductive patterns in mountain gorillas. **If a comparable study were conducted using appropriate real-world genetic data from an authorized institution, the analysis could follow a similar approach and produce results structured like those demonstrated here.**

Any future research involving real institutional or field genetic data would require the appropriate **data-access permission, research authorization, and applicable ethical or regulatory approvals**.

## Project overview

Can a genotype-based workflow recover simulated parentage and characterize kinship and genetic diversity in a structured population? This project demonstrates genotype QC, LOD-based paternity assignment, reproductive-skew analysis, relatedness estimation, and population-genetic summaries.

**Marker limitation:** The 15 STR-like locus names resemble a human forensic marker panel; this repository does **not** establish that the markers are validated for mountain gorillas.

## Workflow and formulas

### 1. Parentage assignment — LOD score

The custom routine uses a **CERVUS-style** likelihood comparison; it is not a validated implementation of CERVUS.

$$
\mathrm{LOD}(c)=\sum_{\ell=1}^{L}\log_{10}\left[\frac{P(G_o^{\ell}\mid G_m^{\ell},G_c^{\ell})}{P(G_o^{\ell}\mid G_m^{\ell})}\right]
$$

Here, \(G_o\), \(G_m\), and \(G_c\) denote offspring, mother, and candidate-father genotypes at locus \(\ell\). Scores compare the support for a candidate-father hypothesis across loci.

### 2. Reproductive skew — B

The current code calculates an **equal-sharing deviation score**:

$$
B_{\mathrm{simple}}=\sum_{i=1}^{k}\left(p_i-\frac{1}{k}\right)^2
$$

where \(p_i\) is male \(i\)'s share of assigned paternities and \(k\) is the number of sires included. This is a simple skew score, **not the full exposure-adjusted Nonacs B-index**. A commonly used Nonacs formulation accounts for each male's exposure time:

$$
B_{\mathrm{Nonacs}}=\sum_{i=1}^{N}\left(p_i-\frac{t_i}{T}\right)^2-\frac{1}{R}\left(1-\frac{1}{\bar N}\right)
$$

where \(p_i\) is individual \(i\)'s share of reproductive success, \(t_i\) is exposure time, \(T=\sum_i t_i\), \(R\) is total reproductive success, and \(\bar N=T/t_{\max}\) is the exposure-adjusted group-size term. The reported skew values below are retained from the current notebook output, but should be treated as **diagnostic values pending correction and validation of the implementation**.

### 3. Genetic-diversity metrics

Expected heterozygosity:

$$H_e=1-\sum_{a=1}^{A}p_a^2$$

Observed heterozygosity:

$$H_o=\frac{\text{heterozygous genotypes}}{\text{genotypes scored}}$$

Polymorphic information content (multiallelic form):

$$\mathrm{PIC}=1-\sum_{i=1}^{A}p_i^2-\sum_{i=1}^{A-1}\sum_{j=i+1}^{A}2p_i^2p_j^2$$

A simplified inbreeding-coefficient summary and heterozygosity-based differentiation measure are:

$$F_{IS}=1-\frac{H_o}{H_e}\qquad\qquad F_{ST}=\frac{H_T-H_S}{H_T}$$

In this repository, FIS and FST are descriptive outputs from the implemented estimators; the FST calculation is a **Nei-style approximation**, not a fully validated population-genetic inference workflow.

### 4. Relatedness and group comparisons

The notebooks calculate Queller–Goodnight-style and Lynch–Ritland relatedness estimates, allele-sharing summaries, within-/between-group comparisons, and marker-level diversity. Pairwise records share individuals, so conventional tests on all pairwise values may violate independence assumptions.

## Results — all simulated

### Result 1 — Genotype quality control

| Metric | Reported result |
|---|---:|
| Individuals before QC | 63 |
| Individuals retained | 63/63 |
| Individuals with >20% missing loci | 0 |
| Loci evaluated / retained | 15 / 15 |
| Mean missing genotype rate | 5.40% |
| Highest locus missingness | D8S1179: 15.87% (10/63 records) |
| Loci exceeding 20% missingness | 0 |

**Interpretation:** All 63 synthetic individuals and 15 loci passed the stated missingness filter. This describes completeness of the generated matrix, not laboratory quality or marker validity.

### Result 2 — Ground-truth paternity versus recovered assignments

The simulation's constructed paternity labels were:

| Simulated ground-truth category | Offspring | Share |
|---|---:|---:|
| Dominant / rank-1 sire | 12 | 60% |
| Subordinate sire | 2 | 10% |
| `UNSAMPLED` father | 6 | 30% |
| **Total** | **20** | **100%** |

The analysis pipeline recovered:

| Pipeline outcome | Offspring | Share of all offspring |
|---|---:|---:|
| Assigned a candidate father | 4 | 20% |
| Unresolved | 16 | 80% |
| Assigned father ranked dominant | 4 | 20% |
| Assigned father ranked subordinate | 0 | 0% |

| Offspring | Assigned father | LOD | Reported confidence | Matches simulated truth? |
|---|---|---:|---|---|
| JUV019 | SB011 | 6.3399 | High | Yes |
| JUV020 | SB011 | 0.5851 | Low | Yes |
| JUV029 | SB021 | 0.9689 | Low | Yes |
| JUV041 | SB033 | 1.1094 | Medium | Yes |

Confidence summary: **1 high (5%), 1 medium (5%), 2 low (10%), and 16 unresolved (80%)**. Of the 14 offspring whose true father was represented among candidates, **4/14 (28.6%)** were assigned and **10/14 (71.4%)** remained unresolved. All four recovered calls matched the stored simulated labels, but two were low-confidence and no subordinate sire was recovered despite two subordinate-sired offspring in the ground truth.

**Interpretation:** The main result is incomplete recovery, not evidence that subordinate paternity is absent. The pipeline assigned only 20% of all offspring.

### Result 3 — Dominance rank among recovered sires

| Statistic | Result |
|---|---:|
| Assigned offspring included | 4 |
| Mean rank of assigned sires | 1.00 |
| Mean rank among candidate males | 1.67 |
| Mann–Whitney U | 12 |
| Reported p-value | 0.0526 |
| Significant at \(\alpha=0.05\)? | No |

All four recovered fathers were rank 1 in the simulation. The reported test did not meet the 0.05 threshold, and the small number of assignments limits interpretation.

### Result 4 — Reproductive-skew diagnostic

| Metric from current notebook | Value |
|---|---:|
| Observed B score | 0.041667 |
| Reference value labelled \(E(B)\) | 0.222222 |
| Reported B/reference ratio | 0.188× |
| Permutation replicates | 2,000 |
| Reported permutation p-value | 0.7425 |

The reported p-value is not significant at 0.05. These values are preserved as printed by the current notebook; the implementation must be corrected and rerun before interpreting them as a validated Nonacs B-index.

### Result 5 — Overall genetic diversity

| Metric | Simulated result |
|---|---:|
| Mean expected heterozygosity, \(H_e\) | 0.8124 ± 0.0501 |
| Mean observed heterozygosity, \(H_o\) | 0.7014 ± 0.0570 |
| Mean PIC | 0.7896 |
| Mean \(F_{IS}\) | 0.1364 |
| Loci summarized | 15 |
| Individuals genotyped | 63 |

Mean observed heterozygosity was lower than expected heterozygosity in the overall summary. These are properties of the synthetic genotypes, not estimates of real gorilla diversity or inbreeding.

### Result 6 — Marker-level diversity and differentiation

`p(HWE)` below is the notebook's simplified chi-square approximation, not an exact multi-allelic Hardy–Weinberg test.

| Locus | N genotypes | Alleles | Ho | He | FIS | PIC | Approx. p(HWE) | FST |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| D2S1338 | 59 | 13 | 0.7119 | 0.8862 | 0.1968 | 0.8760 | 0.0000 | 0.0548 |
| D3S1358 | 60 | 7 | 0.6500 | 0.8249 | 0.2120 | 0.8012 | 0.0004 | 0.0933 |
| D5S818 | 57 | 9 | 0.7895 | 0.8270 | 0.0454 | 0.8059 | 0.4535 | 0.0970 |
| D7S820 | 58 | 9 | 0.7759 | 0.8291 | 0.0642 | 0.8087 | 0.2817 | 0.1116 |
| D8S1179 | 53 | 9 | 0.7547 | 0.8491 | 0.1111 | 0.8312 | 0.0551 | 0.0783 |
| D13S317 | 62 | 7 | 0.6774 | 0.8037 | 0.1571 | 0.7784 | 0.0123 | 0.0725 |
| D16S539 | 62 | 7 | 0.7581 | 0.8258 | 0.0821 | 0.8023 | 0.1594 | 0.0757 |
| D18S51 | 59 | 11 | 0.7119 | 0.8065 | 0.1174 | 0.7906 | 0.0657 | 0.0404 |
| D19S433 | 61 | 6 | 0.5738 | 0.7044 | 0.1854 | 0.6727 | 0.0254 | 0.0980 |
| D21S11 | 59 | 11 | 0.7288 | 0.8532 | 0.1458 | 0.8367 | 0.0069 | 0.0839 |
| vWA | 61 | 8 | 0.6557 | 0.7386 | 0.1122 | 0.7100 | 0.1405 | 0.0866 |
| FGA | 60 | 12 | 0.7000 | 0.8750 | 0.2000 | 0.8631 | 0.0000 | 0.0902 |
| CSF1PO | 61 | 7 | 0.7049 | 0.7734 | 0.0886 | 0.7367 | 0.2010 | 0.0679 |
| TPOX | 62 | 6 | 0.6452 | 0.7591 | 0.1501 | 0.7213 | 0.0359 | 0.0741 |
| TH01 | 60 | 7 | 0.6833 | 0.8306 | 0.1773 | 0.8090 | 0.0024 | 0.1006 |

Mean FST-style estimate: **0.0817**; range: **0.0404 (D18S51) to 0.1116 (D7S820)**. These are descriptive values from the synthetic group structure.

### Result 7 — Diversity by simulated group

| Group label | Individuals | Mean Ho | Mean He |
|---|---:|---:|---:|
| Susa | 10 | 0.6585 | 0.7641 |
| Amahoro | 10 | 0.7572 | 0.7261 |
| Umubano | 12 | 0.7037 | 0.7664 |
| Sabyinyo | 9 | 0.7057 | 0.7298 |
| Hirwa | 13 | 0.6997 | 0.7653 |
| Kwitonda | 9 | 0.6667 | 0.7253 |
| **Total** | **63** | — | — |

Amahoro's reported mean Ho (0.7572) exceeded its mean He (0.7261); the other five group summaries had mean Ho below mean He. All are synthetic summaries.

### Result 8 — Within-group versus between-group relatedness

| Pair category | Pair count | Mean rQG | Median rQG | SD |
|---|---:|---:|---:|---:|
| Within group | 306 | 0.0203 | -0.0012 | 0.1074 |
| Between groups | 1,647 | -0.0153 | -0.0184 | 0.0701 |

Reported Mann–Whitney result: **U = 294,336, p = 0.000001**. The mean difference was **0.0356** relatedness units, with within-group values higher on average. Because pairs reuse individuals, the nominal p-value may violate the test's independence assumption and should not be treated as definitive inference.

### Result 9 — Relatedness by simulated group

| Group label | Within-group pairs | Mean rQG | Median rQG | SD | Maximum rQG | Pairs with rQG > 0.25 |
|---|---:|---:|---:|---:|---:|---:|
| Susa | 45 | 0.0108 | -0.0157 | 0.1213 | 0.3337 | 6.7% |
| Amahoro | 45 | 0.0307 | 0.0072 | 0.1029 | 0.3046 | 2.2% |
| Umubano | 66 | 0.0172 | -0.0046 | 0.1030 | 0.3681 | 6.1% |
| Sabyinyo | 36 | 0.0320 | 0.0073 | 0.1068 | 0.3179 | 8.3% |
| Hirwa | 78 | 0.0123 | -0.0149 | 0.1098 | 0.3715 | 3.8% |
| Kwitonda | 36 | 0.0307 | 0.0087 | 0.1017 | 0.3582 | 2.8% |

Sabyinyo had the highest reported mean within-group relatedness (0.0320); Susa had the lowest (0.0108). Hirwa had the highest reported maximum (0.3715). These are not rankings of real gorilla groups.

## Repository structure

```text
simulated-gorilla-paternity-analysis/
├── data/
├── notebooks/
│   ├── 01_data_import_and_qc.ipynb
│   ├── 02_population_genetics.ipynb
│   ├── 03_parentage_analysis.ipynb
│   ├── 04_kinship_analysis.ipynb
│   └── 05_results_and_conservation.ipynb
├── scripts/
├── outputs/
│   ├── figures/
│   └── tables/
├── requirements.txt
├── environment.yml
└── README.md
```

The results notebook reports **33 generated output files**: **17 figures**, **13 analytical CSV tables**, and **3 processed CSV files** (`individuals_qc.csv`, `genotypes_qc.csv`, and `pairwise_relatedness.csv`). The auxiliary `qc_loci.txt` is listed separately and is not included in the 33-file count.

## Limitations

- All data and findings are synthetic.
- Marker suitability for mountain gorillas is not established.
- Parentage recovery is low: **4/20 (20%)**.
- The parentage procedure is custom and CERVUS-style, not validated CERVUS software.
- The skew calculation, HWE approximation, FST-style estimator, and relatedness significance test require validation before inferential use.
- No conservation recommendations for real populations can be drawn from these outputs.

## References

1. Marshall TC, et al. (1998). Statistical confidence for likelihood-based paternity inference. *Molecular Ecology*, 7, 639–655.
2. Nonacs P. (2000). Measuring and using skew in the study of social behaviour. *The American Naturalist*, 156, 577–589.
3. Queller DC, Goodnight KF. (1989). Estimating relatedness using genetic markers. *Evolution*, 43, 258–275.
4. Lynch M, Ritland K. (1999). Estimation of pairwise relatedness with molecular markers. *Genetics*, 152, 1753–1766.
5. Bradley BJ, et al. (2004). Dispersal patterns, phenotype, and genotype in free-ranging gorillas. *American Journal of Primatology*, 62, 1–14.
6. Nsubuga AM, et al. (2010). Gorilla microsatellite loci. *Molecular Ecology Resources*, 10(2), 397–399.
7. [NIST STR Fact Sheets](https://strbase-archive.nist.gov/str_fact.htm).
8. [15-locus AmpFlSTR Identifiler marker panel](https://www.sciencedirect.com/science/article/pii/S1875176808001650).

---

**Author:** Bruno Rudahongerwa Rwarahoze  
*Former Data Entry Intern at Dian Fossey Gorilla Fund (DFGF)*

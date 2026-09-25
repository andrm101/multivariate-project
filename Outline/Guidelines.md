# DSK805 - Course Project Assignment

**Instructor:** Henry Kirveslahti
**Date:** 2026-04-09

---

## Gene Set Enrichment Analysis

Gene set enrichment analysis (GSEA) [1] is a statistical method specifically designed to analyze gene expression data.

The goal of this type of analysis is to use prior knowledge of genetic pathways and to systematically study their relation to phenotypic variation. In other words, to relate meaningful genetic signal to characteristics of individuals that we can directly observe.

In this project we use the gene set enrichment analysis to study 5 genetic pathways and their relation to breast cancer metastasis. We also compare these results against the more traditional analysis of group differences based on Hotelling's T²-statistic.

---

## The Data Set

For this project we use data from two sources: a) 50 Hallmark genetic pathways from the Molecular Signatures Database [2] and b) the GSE2034 Breast cancer data set [3] from the Gene Expression Omnibus.

Links to the raw data are available, but since the data are quite messy and involve combining information from two sources with slightly different conventions, they have been pre-cleaned. The cleaned data comprises 3 `.RDS` files:

- **`gene_expressions.RDS`** — Contains the gene expression data frame with abundancy reads for *g* = 12,643 genes for *n* = 286 patients. **NB:** The data frame is *g × n* (rows index features, columns index samples — contrary to the usual course setup).
- **`group.RDS`** — A binary vector indicating which of the *n* patients belong to the two groups: relapse (*n₁* = 69) and no-relapse (*n₂* = 217), indicating breast cancer bone metastasis vs. no metastasis, respectively.
- **`HML.RDS`** — Contains 50 Hallmark genetic pathways. These are expert-chosen specific genetic pathways (sets of genes) that share a common biological theme or are associated with some phenotype. From a technical perspective, this is a named list of lists. Each element contains a list of genes corresponding to a subset of rows of `gene_expressions.RDS`. Your job is to analyze 5 of these pathways, assigned to you via the `genetic_lottery` randomization procedure.

---

## The Tools

To help with the analysis, you are given pre-written scripts. You do not need to write the analysis code from scratch.

- **`hotelling.R`** — Used to compute one- and two-sample Hotelling's T² tests.
- **`GSEA_functions.R`** — Contains 3 functions corresponding to the 3 steps of GSEA:
  1. `rank_genes` — Ranks genes based on their Signal-to-Noise Ratio and returns them in descending order.
  2. `compute_enrichment_score` — Computes the enrichment score for a gene set S, given ranked genes.
  3. `approximate_null_distribution` — Generates samples from the null distribution for the null hypothesis that gene set S is not differentially expressed between the two groups.
  - If you want to use other methods for steps 1–3, that is fine, but you will need to write or cite the code and describe what it does in sufficient detail.
- **`dks805utils.R`** — Utility functions including the Yeo-Johnson power transform to try to force the data to be normally distributed. Requires the `car` package.
- **`lottery.R`** — Used to determine the gene sets assigned to your group.

You can call the functions from these scripts directly by placing them in your working folder and running:

```r
source('hotelling.R')
source('GSEA_functions.R')
source('lottery.R')
source('dks805utils.R')
```

---

## Hints

- Some pathways may contain genes not present in the expression matrix. The statistically principled approach is to ignore them. For example, if your pathway contains 200 genes but only 184 are in the expression matrix, use those 184 genes (and use *p* = 184 for a Hotelling's T² test).
- Take a look at the `.RMD` file for the lectures on April 9 — it contains useful R commands for subsetting, and the `.pdf` files describe the GSEA steps.
- Fix your group early. The exact assignment (pathways given by `genetic_lottery`) depends on group members.
- Avoid dividing work pathway-by-pathway; coordinate roles instead (e.g., one person sets up the report structure, another handles Hotelling's T², a third sets up GSEA).
- Computing the null distribution can be computationally costly. Aim for at least **1,000 iterations**. Consider running on **ucloud** if you need more computational power.
- The GSEA makes no specific distributional requirements. However, for Hotelling's T² p-values, the data should be normally distributed with similar variances across groups. If assumptions are violated:
  - a) Try to transform the data (e.g., Yeo-Johnson transform).
  - b) Use a permutation test.
  - c) Document the violated assumptions and account for them in your results.

  Option c) is much better than ignoring the problem entirely. Note: testing equality of variance is not part of this project.

---

## Your Job

Following the development in the lecture slides, perform a gene set enrichment analysis for 5 genetic pathways determined by your group composition. Compare the results against Hotelling's T² test. Summarize your findings in a brief report, to be submitted on **itslearning by April 30**.

### Report Structure

Your report should contain the following:

1. **Project header** (see below)
2. **Generative AI statement** (separate form — see attached template for instructions; not counted toward page limit)
3. **Brief introduction** to what you are doing
4. **Description of methods** (naming what you did; no need for full mathematical details)
5. **Assumptions** behind the methods and whether they are met (include plots where relevant)
6. **Findings**, contrasting results from GSEA and Hotelling's T² approach, and discussing any relevant background from the literature on the genetic pathway

Submit your report as a **`.pdf` file**. Put any supporting code in a separate file. The Generative AI form should also be a separate file (e.g., `.docx`).

---

## Assignment Rules

- Open book, open everything assignment.
- Group size is **at most 5**.
- Collaboration across teams is welcome, but note that other groups may work on entirely different pathways.
- A discussion board on itslearning is available.
- **Maximum report length: 6 pages** (shorter is better, as long as all components are present).
- Additional details can go in an appendix, but the report should be self-contained.
- **Deadline: April 30**, submitted on itslearning.
- **Warning:** Tampering with the `genetic_lottery` function is considered **exam fraud**.

---

## The Project Header

On top of the first page of your project, include a header as follows:

| Names | Name1     | Name2 | Name3 |
|-------|-----------|-------|-------|
| Ids   | id1       | id2   | id3   |
| Hash  | 5921      |       |       |
| Set 1 | HALLMARK_PANCREAS_BETA_CELLS         | | |
| Set 2 | HALLMARK_NOTCH_SIGNALING             | | |
| Set 3 | HALLMARK_INFLAMMATORY_RESPONSE       | | |
| Set 4 | HALLMARK_TNFA_SIGNALING_VIA_NFKB     | | |
| Set 5 | HALLMARK_GLYCOLYSIS                  | | |

Replace `Name1`, `Name2`, `Name3` with your actual names, and `id1`, `id2`, `id3` with your actual SDU email prefixes. The rest is generated by `genetic_lottery()` from `lottery.R`.

### Specific Rules About the Header

- Use the provided `.RMD` to generate the header (see `assignment.RMD`). If you cannot use LaTeX with R, you may construct the header manually, but you **must** use the `genetic_lottery()` function from `lottery.R` to compute the hash and obtain the gene sets.
- In that case, provide the full input and output, like:

```r
source('lottery.R')
x <- c('id1', 'id2', 'id3')  # list of your IDs
print(x)
## [1] "id1" "id2" "id3"
print(genetic_lottery(x))
## [[1]]
## [1] 5921
## [[2]]
## [1] 19 10 22 7 28
## [[3]]
## [1] "HALLMARK_GLYCOLYSIS"        "HALLMARK_COAGULATION"
## [3] "HALLMARK_HYPOXIA"           "HALLMARK_APOPTOSIS"
## [5] "HALLMARK_KRAS_SIGNALING_DN"
```

- Include both the input and the full output of this script.
- Your **ID** is your SDU email address prefix (everything before the `@` symbol).
- Ids may be provided in any order.

---

## Grading Rubric

The assignment is graded on 3 dimensions (weights are approximate and act multiplicatively — avoid poor performance in any single area):

| Dimension | Weight | Description |
|-----------|--------|-------------|
| **Content** | 40% | Relevance, completeness, and depth of statistical analysis (did you do everything asked?) |
| **Accuracy** | 40% | Correctness of methods, calculations, and interpretations (do you understand what you did?) |
| **Clarity** | 20% | Organization, explanation, and presentation (is the report easy to read, on time, and well formatted?) |

The **target audience** is another student on this course.

### Expected Learning Outcomes Evaluated

- Understand and identify problems that can be solved using multivariate techniques.
- Perform a practical data analysis with techniques from the course.
- Perform programming relevant to course content in R.
- Summarize the results of an analysis in a statistical report.

---

## Generative AI Statement

Filling out the Generative AI statement (`GAI_form.docx`) is **optional but highly recommended**. See detailed instructions in the form. You must comply with generative AI guidelines whether or not you fill in the form.

Good uses of AI in this project include: helping make prettier graphics, and having AI explain what the Hallmark pathways mean intuitively. This is a project in multivariate statistics, not molecular biology.

---

## References

[1] A. Subramanian, P. Tamayo, V. K. Mootha, S. Mukherjee, B. L. Ebert, M. A. Gillette, A. Paulovich, S. L. Pomeroy, T. R. Golub, E. S. Lander, J. P. Mesirov, *Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles*, Proceedings of the National Academy of Sciences, 102(43):15545–15550, 2005.

[2] A. Liberzon, C. Birger, H. Thorvaldsdottir, M. Ghandi, J. P. Mesirov, P. Tamayo, *The Molecular Signatures Database (MSigDB) hallmark gene set collection*, Cell Systems, 1(6):417–425, 2015.

[3] M. J. van de Vijver, Y. He, L. J. van 't Veer, H. Dai, A. A. M. Hart, D. W. Voskuil, G. J. Schreiber, H. L. Peterse, C. Roberts, M. J. Marton, M. Parrish, D. Atsma, A. Witteveen, A. M. Glas, L. Delahaye, T. van der Velde, H. Bartelink, S. Rodenhuis, E. J. T. Rutgers, S. H. Friend, *A gene-expression signature as a predictor of survival in breast cancer*, The New England Journal of Medicine, 347:1999–2009, 2002.

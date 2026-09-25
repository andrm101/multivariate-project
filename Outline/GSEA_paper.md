# Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles

**Authors:** Aravind Subramanian, Pablo Tamayo, Vamsi K. Mootha, Sayan Mukherjee, Benjamin L. Ebert, Michael A. Gillette, Amanda Paulovich, Scott L. Pomeroy, Todd R. Golub, Eric S. Lander, and Jill P. Mesirov

**Contributed by:** Eric S. Lander, August 2, 2005

**Source:** PNAS, October 25, 2005, vol. 102, no. 43, pp. 15545–15550

Although genomewide RNA expression analysis has become a routine tool in biomedical research, extracting biological insight from such information remains a major challenge. Here, we describe a powerful analytical method called Gene Set Enrichment Analysis (GSEA) for interpreting gene expression data. The method derives its power by focusing on gene sets, that is, groups of genes that share common biological function, chromosomal location, or regulation. We demonstrate how GSEA yields insights into several cancer-related data sets, including leukemia and lung cancer. Notably, where single-gene analysis finds little similarity between two independent studies of patient survival in lung cancer, GSEA reveals many biological pathways in common. The GSEA method is embodied in a freely available software package, together with an initial database of 1,325 biologically defined gene sets.

**Keywords:** microarray

Genomewide expression analysis with DNA microarrays has become a mainstay of genomics research (1, 2). The challenge no longer lies in obtaining gene expression profiles, but rather in interpreting the results to gain insights into biological mechanisms.

In a typical experiment, mRNA expression profiles are generated for thousands of genes from a collection of samples belonging to one of two classes, for example, tumors that are sensitive vs. resistant to a drug. The genes can be ordered in a ranked list L, according to their differential expression between the classes. The challenge is to extract meaning from this list.

A common approach involves focusing on a handful of genes at the top and bottom of L (i.e., those showing the largest difference) to discern telltale biological clues. This approach has a few major limitations.

1. After correcting for multiple hypotheses testing, no individual gene may meet the threshold for statistical significance, because the relevant biological differences are modest relative to the noise inherent to the microarray technology.
2. Alternatively, one may be left with a long list of statistically significant genes without any unifying biological theme. Interpretation can be daunting and ad hoc, being dependent on a biologist’s area of expertise.
3. Single-gene analysis may miss important effects on pathways. Cellular processes often affect sets of genes acting in concert. An increase of 20% in all genes encoding members of a metabolic pathway may dramatically alter the flux through the pathway and may be more important than a 20-fold increase in a single gene.
4. When different groups study the same biological system, the list of statistically significant genes from the two studies may show distressingly little overlap (3).

To overcome these analytical challenges, a method called Gene Set Enrichment Analysis (GSEA) was developed that evaluates microarray data at the level of gene sets. The gene sets are defined based on prior biological knowledge, e.g., published information about biochemical pathways or coexpression in previous experiments. The goal of GSEA is to determine whether members of a gene set S tend to occur toward the top (or bottom) of the list L, in which case the gene set is correlated with the phenotypic class distinction.

A preliminary version of GSEA was used to analyze data from muscle biopsies from diabetics vs. healthy controls (4). The method revealed that genes involved in oxidative phosphorylation show reduced expression in diabetics, although the average decrease per gene is only 20%. The results from this study have been independently validated by other microarray studies (5) and by in vivo functional studies (6).

Given this success, GSEA was developed into a robust technique for analyzing molecular profiling data. Its characteristics and performance were studied and the original method was substantially revised and generalized for broader applicability.

This paper provides a full mathematical description of the GSEA methodology and illustrates its utility by applying it to several diverse biological problems. A software package, called GSEA-P, and an initial inventory of gene sets (Molecular Signature Database, MSigDB) were also created.

## Methods

### Overview of GSEA

GSEA considers experiments with genomewide expression profiles from samples belonging to two classes, labeled 1 or 2. Genes are ranked based on the correlation between their expression and the class distinction by using any suitable metric.

Given an a priori defined set of genes S (e.g., genes encoding products in a metabolic pathway, located in the same cytogenetic band, or sharing the same GO category), the goal of GSEA is to determine whether the members of S are randomly distributed throughout L or primarily found at the top or bottom. Sets related to the phenotypic distinction are expected to show the latter distribution.

There are three key elements of the GSEA method:

#### Step 1: Calculation of an Enrichment Score

An enrichment score (ES) reflects the degree to which a set S is overrepresented at the extremes (top or bottom) of the entire ranked list L. The score is calculated by walking down the list L, increasing a running-sum statistic when encountering a gene in S and decreasing it when encountering genes not in S. The magnitude of the increment depends on the correlation of the gene with the phenotype. The enrichment score is the maximum deviation from zero encountered in the random walk; it corresponds to a weighted Kolmogorov–Smirnov-like statistic.

#### Step 2: Estimation of Significance Level of ES

The statistical significance (nominal P value) of the ES is estimated by using an empirical phenotype-based permutation test procedure that preserves the complex correlation structure of the gene expression data. Specifically, the phenotype labels are permuted and the ES of the gene set is recomputed for the permuted data, generating a null distribution for the ES. The empirical nominal P value of the observed ES is then calculated relative to this null distribution.

Importantly, permutation of class labels preserves gene-gene correlations and thus provides a more biologically reasonable assessment of significance than would be obtained by permuting genes.

#### Step 3: Adjustment for Multiple Hypothesis Testing

When an entire database of gene sets is evaluated, the estimated significance level is adjusted to account for multiple hypothesis testing. First, the ES for each gene set is normalized to account for the size of the set, yielding a normalized enrichment score (NES). Then the proportion of false positives is controlled by calculating the false discovery rate (FDR) corresponding to each NES. The FDR is the estimated probability that a set with a given NES represents a false positive finding.

The implementation differs in several important ways from the preliminary version. In the original implementation, the running-sum statistic used equal weights at every step, which yielded high scores for sets clustered near the middle of the ranked list. These sets do not represent biologically relevant correlation with the phenotype. This issue was addressed by weighting the steps according to each gene’s correlation with a phenotype.

The original implementation used familywise-error rate (FWER) for multiple-testing correction, but this criterion was so conservative that many applications yielded no statistically significant results. Because the primary goal is hypothesis generation, the method instead uses FDR to focus on controlling the probability that each reported result is a false positive.

Based on statistical analysis and empirical evaluation, GSEA shows broad applicability. It can detect subtle enrichment signals and preserves the original oxidative phosphorylation result in normal samples (P = 0.008, FDR = 0.04). This methodology was implemented in a software tool called GSEA-P.

### The Leading-Edge Subset

Gene sets can be defined in many ways, but not all members of a gene set typically participate in a biological process. It is often useful to extract the core members of high-scoring gene sets that contribute to the ES.

The leading-edge subset is defined as those genes in the gene set S that appear in the ranked list L at, or before, the point where the running sum reaches its maximum deviation from zero. The leading-edge subset can be interpreted as the core of a gene set that accounts for the enrichment signal.

Examination of the leading-edge subset can reveal a biologically important subset within a gene set, as shown in the analysis of p53 status in cancer cell lines. High-scoring gene sets can also be grouped on the basis of shared leading-edge subsets, revealing which gene sets correspond to the same biological processes and which represent distinct processes.

### Variations of the GSEA Method

The paper focuses mainly on using GSEA to analyze a ranked gene list reflecting differential expression between two classes, each represented by a large number of samples. However, the method can be applied to ranked gene lists arising in other settings.

For small data sets with too few samples to allow rigorous evaluation of significance by permuting class labels, a P value can be estimated by permuting the genes while maintaining gene set size. This is not strictly accurate because it ignores gene-gene correlations and may overestimate significance, but it can still be useful for hypothesis generation.

Genes may also be ranked based on how well their expression correlates with a given target pattern, such as the expression pattern of a particular gene. Approximate P values can again be estimated by permutation of genes.

### An Initial Catalog of Human Gene Sets

GSEA evaluates a query microarray data set by using a collection of gene sets. An initial catalog of 1,325 gene sets, called MSigDB 1.0, was created and consisted of four types of sets:

- **Cytogenetic sets (C1, 319 gene sets):** 24 sets for the 24 human chromosomes and 295 sets corresponding to cytogenetic bands.
- **Functional sets (C2, 522 gene sets):** 472 sets containing genes involved in specific metabolic and signaling pathways from eight manually curated databases, and 50 sets containing genes coregulated in response to perturbations.
- **Regulatory-motif sets (C3, 57 gene sets):** based on 57 commonly conserved regulatory motifs in promoter regions of human genes.
- **Neighborhood sets (C4, 427 gene sets):** defined by expression neighborhoods centered on cancer-related genes.

### GSEA-P Software and MSigDB Gene Sets

To facilitate the use of GSEA, resources were developed and made available from the Broad Institute upon request. These included the GSEA-P software, MSigDB 1.0, and documentation.

The software was made available as:

- A platform-independent desktop application with a graphical user interface.
- Programs in R and Java for advanced users.
- An analytic module in the GenePattern microarray analysis package.
- A future web-based GSEA server for running analyses online.

## Results

The ability of GSEA to provide biologically meaningful insights was explored in six examples for which considerable background information was available. In each case, significantly associated gene sets were searched from one or both of the subcatalogs C1 and C2.

### Male vs. Female Lymphoblastoid Cells

mRNA expression profiles were generated from lymphoblastoid cell lines derived from 15 males and 17 females, and gene sets correlated with the distinctions “male/female” and “female/male” were identified.

When cytogenetic gene sets (C1) were tested, the male/female comparison yielded chromosome Y and the two Y bands with at least 15 genes (Yp11 and Yq11), as expected. For the female/male comparison, enrichment for chromosome X bands was not expected because most X-linked genes are subject to dosage compensation.

When functional gene sets (C2) were considered, three biologically informative sets were identified. One consisted of genes escaping X inactivation, showing the expected enrichment in female cells. Two additional sets consisted of genes enriched in reproductive tissues (testis and uterus), which remained significant even when restricted to autosomal genes.

### p53 Status in Cancer Cell Lines

Gene expression patterns from the NCI-60 collection of cancer cell lines were examined to identify targets of the transcription factor p53. Of 50 cell lines with reported p53 mutational status, 17 were classified as normal and 33 as carrying mutations.

GSEA identified five functional gene sets correlated with normal p53 function. These included genes in the p53 signaling pathway, downstream targets of p53, radiation-induced genes known to involve p53, hypoxia-induced genes acting through a p53-mediated pathway, and genes encoding heat shock-protein signaling pathways.

The complementary analysis identified one significant gene set in p53 mutant cells: the Ras signaling pathway. Two additional near-threshold sets involved the Ngf and Igf1 signaling pathways. Examination of leading-edge subsets showed a common subgroup containing MAP2K1, RAF1, ELK1, and PIK3CA, pointing to up-regulation of this MAPK subpathway as a key distinction between p53 mutant and wild-type tumors.

### Acute Leukemias

Acute lymphoid leukemia (ALL) and acute myeloid leukemia (AML) were studied by comparing gene expression profiles from 24 ALL patients and 24 AML patients.

Application of GSEA to cytogenetic gene sets (C1) yielded five gene sets enriched in ALL: chr6q21, chr5q31, chr13q14, chr14q32, and chr17q23. These regions were interpreted in terms of known leukemia biology. For example, 5q31 is consistent with known AML cytogenetics, 13q14 contains the RB locus, and 14q32 contains the Ig heavy chain locus and likely reflects lineage-specific expression.

The reciprocal AML/ALL analysis yielded no significantly enriched bands, likely reflecting the relative infrequency of deletions in ALL.

### Comparing Two Studies of Lung Cancer

To test robustness, two independent lung adenocarcinoma studies were reanalyzed: one from Boston (n = 62) and one from Michigan (n = 86), both with clinical outcomes classified as good or poor.

No genes in either study were strongly associated with outcome at a 5% significance level after correcting for multiple hypothesis testing. Comparing the top 100 genes associated with poor outcome in each study revealed only 12 genes in common, with weak statistical significance.

GSEA revealed much greater similarity. The poor-outcome gene set from Boston showed strong enrichment in the Michigan data (NES = 1.90, P < 0.001), and the poor-outcome gene set from Michigan was enriched in the Boston data (NES = 2.13, P < 0.001).

GSEA also identified significant functional gene sets associated with poor outcome: 8 in the Boston data and 11 in the Michigan data (FDR ≤ 0.25). About half of the significant gene sets were shared between the two studies, and several others were clearly related to the same biological process.

Across the Boston, Michigan, and Stanford studies, two broad themes emerged:

1. **Rapid cellular proliferation:** including gene sets related to Ras activation, cell cycle, hypoxia, angiogenesis, glycolysis, and carbohydrate metabolism.
2. **Amino acid biosynthesis and mTOR signaling:** including increased amino acid biosynthesis, mTor signaling, and up-regulation of genes down-regulated by amino acid deprivation and rapamycin treatment.

The analysis suggested much greater consistency across lung cancer data sets by using GSEA than by single-gene analysis, and it generated more compelling biological hypotheses for further exploration.

## Discussion

Traditional strategies for gene expression analysis have focused on identifying individual genes that differ between two states. Although useful, such approaches fail to detect biological processes distributed across entire gene networks and subtle at the level of individual genes.

GSEA addresses this by analyzing gene sets rather than single genes. It was initially used to discover metabolic pathways altered in human diabetes and later applied to diffuse large B cell lymphoma, nutrient-sensing pathways in prostate cancer, and cross-species expression comparisons.

GSEA has several advantages over single-gene methods:

- It eases interpretation of large-scale experiments by identifying pathways and processes.
- It can boost the signal-to-noise ratio when members of a gene set exhibit strong cross-correlation.
- Leading-edge analysis helps define gene subsets that elucidate the results.

Other pathway-analysis tools typically test whether a group of differentially expressed genes is enriched for a pathway or ontology term by using overlap statistics such as the cumulative hypergeometric distribution. GSEA differs in two important ways:

1. It considers all genes in an experiment, not only those above an arbitrary cutoff.
2. It assesses significance by permuting class labels, preserving gene-gene correlations and providing a more accurate null model.

The paper argues that the real power of GSEA lies in its flexibility. The molecular signature database can include sets based on biological pathways, chromosomal location, upstream cis motifs, responses to drug treatment, or expression profiles from previous microarray studies. Additional sets can be created through perturbation studies, computational analysis, and biological annotation.

## Appendix: Mathematical Description of Methods

### Inputs to GSEA

1. Expression data set D with N genes and k samples.
2. Ranking procedure to produce gene list L, including a correlation (or other ranking metric) and a phenotype or profile of interest C.
3. An exponent p to control the weight of the step.
4. Independently derived gene set S of NH genes (e.g., a pathway, a cytogenetic band, or a GO category).

### Enrichment Score ES(S)

1. Rank order the N genes in D to form L = {g1, ..., gN} according to the correlation, r(gj) = rj, of their expression profiles with C.
2. Evaluate the fraction of genes in S (“hits”) weighted by their correlation and the fraction of genes not in S (“misses”) present up to a given position i in L.

The enrichment score ES is the maximum deviation from zero of \(P_{hit} - P_{miss}\). For a randomly distributed S, ES(S) will be relatively small, but if it is concentrated at the top or bottom of the list, or otherwise nonrandomly distributed, then ES(S) will be correspondingly high.

When p = 0, ES(S) reduces to the standard Kolmogorov–Smirnov statistic; when p = 1, genes in S are weighted by their correlation with C normalized by the sum of the correlations over all genes in S. The paper uses p = 1.

### Estimating Significance

1. Randomly assign the original phenotype labels to samples, reorder genes, and recompute ES(S).
2. Repeat for 1,000 permutations and create a histogram of the corresponding enrichment scores ESNULL.
3. Estimate the nominal P value for S from ESNULL using the positive or negative portion of the distribution corresponding to the sign of the observed ES(S).

### Multiple Hypothesis Testing

1. Determine ES(S) for each gene set in the collection or database.
2. For each S and 1,000 fixed permutations of the phenotype labels, reorder the genes in L and determine ES(S, π).
3. Normalize ES(S, π) and the observed ES(S), separately rescaling positive and negative scores to yield normalized scores NES(S, π) and NES(S).
4. Compute FDR by controlling the ratio of false positives to the total number of gene sets attaining a fixed level of significance.

## References

1. Schena, M., Shalon, D., Davis, R. W. & Brown, P. O. (1995) *Science* 270, 467–470.
2. Lockhart, D. J., Dong, H., Byrne, M. C., Follettie, M. T., Gallo, M. V., Chee, M. S., Mittmann, M., Wang, C., Kobayashi, M., Horton, H., et al. (1996) *Nat. Biotechnol.* 14, 1675–1680.
3. Fortunel, N. O., Otu, H. H., Ng, H. H., Chen, J., Mu, X., Chevassut, T., Li, X., Joseph, M., Bailey, C., Hatzfeld, J. A., et al. (2003) *Science* 302, 393, author reply 393.
4. Mootha, V. K., Lindgren, C. M., Eriksson, K. F., Subramanian, A., Sihag, S., Lehar, J., Puigserver, P., Carlsson, E., Ridderstrale, M., Laurila, E., et al. (2003) *Nat. Genet.* 34, 267–273.
5. Patti, M. E., Butte, A. J., Crunkhorn, S., Cusi, K., Berria, R., Kashyap, S., Miyazaki, Y., Kohane, I., Costello, M., Saccone, R., et al. (2003) *Proc. Natl. Acad. Sci. USA* 100, 8466–8471.
6. Petersen, K. F., Dufour, S., Befroy, D., Garcia, R. & Shulman, G. I. (2004) *N. Engl. J. Med.* 350, 664–671.
7. Hollander, M. & Wolfe, D. A. (1999) *Nonparametric Statistical Methods* (Wiley, New York).
8. Benjamini, Y., Drai, D., Elmer, G., Kafkafi, N. & Golani, I. (2001) *Behav. Brain Res.* 125, 279–284.
9. Reiner, A., Yekutieli, D. & Benjamini, Y. (2003) *Bioinformatics* 19, 368–375.
10. Lamb, J., Ramaswamy, S., Ford, H. L., Contreras, B., Martinez, R. V., Kittrell, F. S., Zahnow, C. A., Patterson, N., Golub, T. R. & Ewen, M. E. (2003) *Cell* 114, 323–334.
11. Xie, X., Lu, J., Kulbokas, E. J., Golub, T. R., Mootha, V., Lindblad-Toh, K., Lander, E. S. & Kellis, M. (2005) *Nature* 434, 338–345.
12. Plath, K., Mlynarczyk-Evans, S., Nusinow, D. A. & Panning, B. (2002) *Annu. Rev. Genet.* 36, 233–278.
13. Carrel, L., Cottle, A. A., Goglin, K. C. & Willard, H. F. (1999) *Proc. Natl. Acad. Sci. USA* 96, 14440–14444.
14. Disteche, C. M., Filippova, G. N. & Tsuchiya, K. D. (2002) *Cytogenet. Genome Res.* 99, 36–43.
15. Olivier, M., Eeles, R., Hollstein, M., Khan, M. A., Harris, C. C. & Hainaut, P. (2002) *Hum. Mutat.* 19, 607–614.
16. Armstrong, S. A., Staunton, J. E., Silverman, L. B., Pieters, R., den Boer, M. L., Minden, M. D., Sallan, S. E., Lander, E. S., Golub, T. R. & Korsmeyer, S. J. (2002) *Nat. Genet.* 30, 41–47.
17. Zhao, N., Stoffel, A., Wang, P. W., Eisenbart, J. D., Espinosa, R., 3rd, Larson, R. A. & Le Beau, M. M. (1997) *Proc. Natl. Acad. Sci. USA* 94, 6948–6953.
18. Barbouti, A., Hoglund, M., Johansson, B., Lassen, C., Nilsson, P. G., Hagemeijer, A., Mitelman, F. & Fioretos, T. (2003) *Cancer Res.* 63, 1202–1206.
19. Tanaka, K., Arif, M., Eguchi, M., Guo, S. X., Hayashi, Y., Asaoku, H., Kyo, T., Dohy, H. & Kamada, N. (1999) *Leukemia* 13, 1367–1373.
20. Morelli, C., Karayianni, E., Magnanini, C., Mungall, A. J., Thorland, E., Negrini, M., Smith, D. I. & Barbanti-Brodano, G. (2002) *Oncogene* 21, 7266–7276.
21. Mrozek, K., Heerema, N. A. & Bloomfield, C. D. (2004) *Blood Rev.* 18, 115–136.
22. Bhattacharjee, A., Richards, W. G., Staunton, J., Li, C., Monti, S., Vasa, P., Ladd, C., Beheshti, J., Bueno, R., Gillette, M., et al. (2001) *Proc. Natl. Acad. Sci. USA* 98, 13790–13795.
23. Beer, D. G., Kardia, S. L., Huang, C. C., Giordano, T. J., Levin, A. M., Misek, D. E., Lin, L., Chen, G., Gharib, T. G., Thomas, D. G., et al. (2002) *Nat. Med.* 8, 816–824.
24. Garber, M. E., Troyanskaya, O. G., Schluens, K., Petersen, S., Thaesler, Z., Pacyna-Gengelbach, M., van de Rijn, M., Rosen, G. D., Perou, C. M., Whyte, R. I., et al. (2001) *Proc. Natl. Acad. Sci. USA* 98, 13784–13789.
25. Smith, L. L., Coller, H. A. & Roberts, J. M. (2003) *Nat. Cell Biol.* 5, 474–479.
26. Acker, T. & Plate, K. H. (2002) *J. Mol. Med.* 80, 562–575.
27. Peng, T., Golub, T. R. & Sabatini, D. M. (2002) *Mol. Cell. Biol.* 22, 5575–5584.
28. Boffa, D. J., Luan, F., Thomas, D., Yang, H., Sharma, V. K., Lagman, M. & Suthanthiran, M. (2004) *Clin. Cancer Res.* 10, 293–300.
29. Monti, S., Savage, K. J., Kutok, J. L., Feuerhake, F., Kurtin, P., Mihm, M., Wu, B., Pasqualucci, L., Neuberg, D., Aguiar, R. C., et al. (2004) *Blood* 105, 1851–1861.
30. Majumder, P. K., Febbo, P. G., Bikoff, R., Berger, R., Xue, Q., McMahon, L. M., Manola, J., Brugarolas, J., McDonnell, T. J., Golub, T. R., et al. (2004) *Nat. Med.* 10, 594–601.
31. Sweet-Cordero, A., Mukherjee, S., Subramanian, A., You, H., Roix, J. J., Ladd-Acosta, C., Mesirov, J., Golub, T. R. & Jacks, T. (2005) *Nat. Genet.* 37, 48–55.
32. Doniger, S. W., Salomonis, N., Dahlquist, K. D., Vranizan, K., Lawlor, S. C. & Conklin, B. R. (2003) *Genome Biol.* 4, R7.
33. Zhong, S., Storch, K. F., Lipan, O., Kao, M. C., Weitz, C. J. & Wong, W. H. (2004) *Appl. Bioinformatics* 3, 261–264.
34. Berriz, G. F., King, O. D., Bryant, B., Sander, C. & Roth, F. P. (2003) *Bioinformatics* 19, 2502–2504.

# Meeting 1 — Step-by-Step Setup Guide
## DSK805 | Gene Set Enrichment Analysis | 4-Person Team

**Duration:** ~1.5 hours  
**Goal:** Everyone unblocked, data loaded, lottery run, shared helper written, computation started.  
**Companion file:** `Meeting1_Setup.R` — run this script together, step by step.

---

## Before the Meeting (5 min each, done individually)

Each team member should do this **before** showing up:

1. Make sure R (≥ 4.0) and RStudio are installed and working.
2. Install required packages by running in R:
   ```r
   install.packages(c("car", "digest", "knitr"))
   ```
3. Place all project files in the **same working directory** and confirm they are there:
   - `gene_expressions.RDS`
   - `group.RDS`
   - `HML.RDS`
   - `hotelling.R`
   - `GSEA_functions.R`
   - `lottery.R`
   - `dsk805utils.R`
4. Set that folder as your working directory in RStudio:  
   *Session → Set Working Directory → To Source File Location*  
   or run `setwd("path/to/your/folder")`.

---

## Step 1 — Source All Scripts and Check Packages (10 min)

Run the following at the top of `Meeting1_Setup.R`. All four team members should see **no errors**.

```r
# ── Step 1: Package check ─────────────────────────────────────────────────────
required_pkgs <- c("car", "digest", "knitr")
missing_pkgs  <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]

if (length(missing_pkgs) > 0) {
  cat("Installing missing packages:", paste(missing_pkgs, collapse = ", "), "\n")
  install.packages(missing_pkgs)
} else {
  cat("✓ All required packages are installed.\n")
}

# ── Step 2: Source the provided scripts ──────────────────────────────────────
source("hotelling.R")      # loads: one_sample_test(), two_sample_test()
source("GSEA_functions.R") # loads: rank_genes(), compute_enrichment_score(),
                           #        approximate_null_distribution(), plot_gsea_heatmap()
source("lottery.R")        # loads: genetic_lottery(), HML vector
source("dsk805utils.R")    # loads: yj_transform()   [requires car]

cat("✓ All scripts sourced successfully.\n")
```

**What to look for:** The console should show no red error messages. If `dsk805utils.R` fails, run `install.packages("car")` and try again.

---

## Step 2 — Load the Data and Verify Dimensions (10 min)

```r
# ── Step 3: Load the three RDS files ─────────────────────────────────────────
gene_expr <- readRDS("gene_expressions.RDS")  # g × n data frame
group     <- readRDS("group.RDS")             # binary vector, length n
HML_data  <- readRDS("HML.RDS")              # named list of 50 pathways

# ── Step 4: Verify dimensions ────────────────────────────────────────────────
cat("── gene_expressions.RDS ──\n")
cat("  Dimensions (genes × patients):", nrow(gene_expr), "×", ncol(gene_expr), "\n")
cat("  Expected:                      12,643 × 286\n\n")

cat("── group.RDS ──\n")
cat("  Length:", length(group), "  (expected: 286)\n")
cat("  Class: ", class(group), "\n")
cat("  Table:\n"); print(table(group))
cat("  Expected: relapse = 69, no-relapse = 217\n\n")

cat("── HML.RDS ──\n")
cat("  Number of pathways:", length(HML_data), "  (expected: 50)\n")
cat("  First few pathway names:\n")
print(head(names(HML_data), 5))
```

**Expected output:**
```
── gene_expressions.RDS ──
  Dimensions (genes × patients): 12643 × 286
  Expected:                      12,643 × 286

── group.RDS ──
  Length: 286  (expected: 286)
  Table:
  no-relapse    relapse
         217         69

── HML.RDS ──
  Number of pathways: 50  (expected: 50)
```

> **⚠ Important orientation note:** `gene_expr` is **genes × patients** (rows = genes, columns = patients).  
> This is the opposite of the usual course setup. The GSEA functions expect this format directly.  
> For Hotelling's T², you will need to transpose subsets to get **patients × genes**.

---

## Step 3 — Run the Genetic Lottery (10 min)

> **Do this once, together, as a group.** The lottery result depends on your member IDs and must be recorded verbatim in the report header.

```r
# ── Step 5: Run the genetic lottery ──────────────────────────────────────────
# Replace with your actual SDU email prefixes (everything before @student.sdu.dk)
# Order does NOT matter — the function sorts them internally.
ids <- c("id1", "id2", "id3", "id4")   # ← FILL IN YOUR ACTUAL IDs

cat("Input IDs:\n"); print(ids)

lottery_result <- genetic_lottery(ids)

cat("\nFull output of genetic_lottery():\n")
print(lottery_result)

cat("\n── Summary ──\n")
cat("Hash value:    ", lottery_result[[1]], "\n")
cat("Pathway indices:", lottery_result[[2]], "\n")
cat("Assigned pathways:\n")
for (i in seq_along(lottery_result[[3]])) {
  cat(sprintf("  Set %d: %s\n", i, lottery_result[[3]][i]))
}
```

**Copy the full printed output into the report header section.** Both the input (`print(ids)`) and the full output (`print(lottery_result)`) must appear in the report.

**Your assigned pathways (confirmed from the assignment):**
```
Set 1: HALLMARK_PANCREAS_BETA_CELLS
Set 2: HALLMARK_NOTCH_SIGNALING
Set 3: HALLMARK_INFLAMMATORY_RESPONSE
Set 4: HALLMARK_TNFA_SIGNALING_VIA_NFKB
Set 5: HALLMARK_GLYCOLYSIS
```

```r
# Store the pathway names in a named vector for easy reference throughout the project
pathway_names <- lottery_result[[3]]
names(pathway_names) <- paste0("Set", 1:5)
print(pathway_names)
```

---

## Step 4 — Write the Shared Gene-Subsetting Helper (15 min)

This is the **most important output of Meeting 1**. Both the Hotelling Analyst and GSEA Analyst must use the exact same gene subsets. Writing this helper once and sharing it prevents any inconsistency.

```r
# ── Step 6: Shared gene-subsetting helper ────────────────────────────────────
#
# get_pathway_genes()
#   Given a pathway name, returns the gene names that exist in BOTH the pathway
#   definition AND the gene expression matrix. Genes in the pathway but absent
#   from the matrix are silently dropped (this is statistically correct — see
#   assignment hints).
#
# Arguments:
#   pathway_name  : character, e.g. "HALLMARK_GLYCOLYSIS"
#   expr_matrix   : the gene_expressions data frame (g × n)
#   hallmark_data : the HML list loaded from HML.RDS
#
# Returns a character vector of valid gene names.

get_pathway_genes <- function(pathway_name, expr_matrix, hallmark_data) {

  # Retrieve the full gene list for this pathway
  pathway_genes_full <- unlist(hallmark_data[[pathway_name]])

  if (is.null(pathway_genes_full)) {
    stop(paste("Pathway not found in HML data:", pathway_name))
  }

  # Find which of those genes are present as row names in the expression matrix
  available_genes <- rownames(expr_matrix)
  valid_genes <- intersect(pathway_genes_full, available_genes)

  # Report the overlap
  n_full  <- length(pathway_genes_full)
  n_valid <- length(valid_genes)
  n_lost  <- n_full - n_valid

  cat(sprintf(
    "  %-45s | %d genes in pathway | %d in matrix | %d dropped\n",
    pathway_name, n_full, n_valid, n_lost
  ))

  return(valid_genes)
}

# ── Step 7: Apply to all 5 pathways and inspect overlaps ─────────────────────
cat("Gene overlap summary for all 5 assigned pathways:\n")
cat(strrep("-", 80), "\n")

pathway_genes_list <- lapply(pathway_names, function(pname) {
  get_pathway_genes(pname, gene_expr, HML_data)
})
names(pathway_genes_list) <- names(pathway_names)

cat(strrep("-", 80), "\n")
cat("\nGene counts per pathway (p used in Hotelling's T²):\n")
print(sapply(pathway_genes_list, length))
```

**Expected output (approximate — exact numbers depend on your gene names):**
```
Gene overlap summary for all 5 assigned pathways:
────────────────────────────────────────────────────────────────────────────────
  HALLMARK_PANCREAS_BETA_CELLS              | ... genes in pathway | ... in matrix | ... dropped
  HALLMARK_NOTCH_SIGNALING                  | ... genes in pathway | ... in matrix | ... dropped
  HALLMARK_INFLAMMATORY_RESPONSE            | ... genes in pathway | ... in matrix | ... dropped
  HALLMARK_TNFA_SIGNALING_VIA_NFKB          | ... genes in pathway | ... in matrix | ... dropped
  HALLMARK_GLYCOLYSIS                       | ... genes in pathway | ... in matrix | ... dropped
```

> **Note for report:** These exact gene counts (`n_valid`) are what you call **p** in your Hotelling's T² test. Record them now.

---

## Step 5 — Prepare the Group Index (5 min)

The `rank_genes()` function needs `inds` as a **logical (TRUE/FALSE) vector**. Prepare and verify it here so both analysts use exactly the same object.

```r
# ── Step 8: Prepare the group boolean index ───────────────────────────────────
# group.RDS may be a factor or character — convert to logical
# TRUE  = relapse    (group 1, n1 = 69)
# FALSE = no-relapse (group 2, n2 = 217)

relapse_inds <- group == "relapse"   # adjust the string if your group uses different labels

cat("Group index summary:\n")
cat("  TRUE  (relapse):    ", sum(relapse_inds),  "\n")
cat("  FALSE (no-relapse): ", sum(!relapse_inds), "\n")
cat("  Total patients:     ", length(relapse_inds), "\n")
cat("  Matches ncol(gene_expr):", length(relapse_inds) == ncol(gene_expr), "\n")

# ── If group labels differ, inspect them first ────────────────────────────────
# Uncomment these lines if the TRUE/FALSE counts above look wrong:
# cat("Unique values in group.RDS:\n"); print(unique(group))
# cat("Table:\n"); print(table(group))
```

---

## Step 6 — Sanity-Check Each Function on One Pathway (15 min)

Run this together so everyone can see a working GSEA pipeline before splitting off. Use `HALLMARK_GLYCOLYSIS` as the test case — it's a well-known pathway and likely to show signal.

```r
# ── Step 9: Test GSEA pipeline on HALLMARK_GLYCOLYSIS ─────────────────────────
test_pathway  <- "Set5"   # HALLMARK_GLYCOLYSIS
test_genes    <- pathway_genes_list[[test_pathway]]
test_name     <- pathway_names[test_pathway]

cat("Testing GSEA pipeline on:", test_name, "\n")
cat("Using", length(test_genes), "genes.\n\n")

# Step 9a: Rank all genes by Signal-to-Noise Ratio
ranked <- rank_genes(gene_expr, relapse_inds)
# ranked$diff  = named numeric vector of SNR scores (length = nrow(gene_expr))
# ranked$genes = full gene expression matrix reordered by SNR descending

cat("rank_genes() output:\n")
cat("  Length of ranked diff vector:", length(ranked$diff), "\n")
cat("  Top 5 genes by SNR:\n")
print(head(ranked$diff, 5))
cat("  Bottom 5 genes by SNR:\n")
print(tail(ranked$diff, 5))

# Step 9b: Compute the enrichment score
es_result <- compute_enrichment_score(ranked$diff, test_genes, p = 1)
# es_result$ES     = running enrichment score vector (length = nrow(gene_expr))
# es_result$ES_max = the enrichment score (the max absolute deviation)

cat("\ncompute_enrichment_score() output:\n")
cat("  Enrichment Score (ES_max):", round(es_result$ES_max, 4), "\n")
cat("  Direction:",
    ifelse(es_result$ES_max > 0,
           "Positive → pathway enriched in RELAPSE group",
           "Negative → pathway enriched in NO-RELAPSE group"), "\n")

# Step 9c: Plot the running ES curve for this pathway
plot(es_result$ES,
     type = "l",
     col  = ifelse(es_result$ES_max > 0, "firebrick", "steelblue"),
     lwd  = 2,
     main = paste("Enrichment Score —", test_name),
     xlab = "Gene rank (high SNR → low SNR)",
     ylab = "Running Enrichment Score")
abline(h = 0, lty = 2, col = "grey50")
abline(v = which.max(abs(es_result$ES)), lty = 3, col = "orange", lwd = 1.5)
legend("topright",
       legend = paste("ES =", round(es_result$ES_max, 4)),
       bty = "n")
```

```r
# Step 9d: Sanity-check Hotelling's T² on the same pathway
# IMPORTANT: gene_expr is g × n → we need to transpose to n × p for Hotelling

subset_mat <- t(gene_expr[test_genes, ])  # now: 286 patients × p genes

data_relapse    <- subset_mat[ relapse_inds, ]   # n1 × p
data_norelapse  <- subset_mat[!relapse_inds, ]   # n2 × p

cat("\ntwo_sample_test() dimension check:\n")
cat("  data_relapse:   ", nrow(data_relapse),   "×", ncol(data_relapse),   "\n")
cat("  data_norelapse: ", nrow(data_norelapse), "×", ncol(data_norelapse), "\n")
cat("  p =", ncol(data_relapse), "(dimension used in T² test)\n\n")

hotelling_result <- two_sample_test(data_relapse, data_norelapse)

cat("Hotelling's T² result for", test_name, ":\n")
cat("  T² statistic:", round(hotelling_result$T2,    4), "\n")
cat("  F statistic: ", round(hotelling_result$Fstat, 4), "\n")
cat("  df1 =", hotelling_result$df1, " df2 =", hotelling_result$df2, "\n")
cat("  p-value:     ", round(hotelling_result$pvalue, 6), "\n")
cat("  Decision at α=0.05:",
    ifelse(hotelling_result$pvalue < 0.05,
           "REJECT H0 — significant group difference",
           "Fail to reject H0 — no significant difference"), "\n")
```

> **⚠ Hotelling's T² note:** With large `p` (many genes), the test may fail because `n1 + n2 - 2 < p` (more dimensions than observations). If that happens, the covariance matrix `S` is singular and `solve(S)` will error. See the troubleshooting section below.

---

## Step 7 — Kick Off the Null Distribution Computation (5 min)

> **Start this immediately after the sanity check.** It is your main computational bottleneck. With 1,000 iterations × 5 pathways, expect 10–60 minutes depending on pathway size and your machine.

```r
# ── Step 10: Launch null distribution for all 5 pathways ─────────────────────
# Run this at the END of Meeting 1, then let it run overnight or on ucloud.
# Results are saved to .RDS so you never need to recompute them.

N_ITERS <- 1000   # minimum required; use 5000 if time allows

cat("Starting null distribution computation...\n")
cat("Iterations per pathway:", N_ITERS, "\n")
cat("Pathways to process:   ", length(pathway_names), "\n\n")

null_distributions <- list()

for (i in seq_along(pathway_names)) {
  pname  <- pathway_names[i]
  pgenes <- pathway_genes_list[[names(pathway_names)[i]]]

  cat(sprintf("[%d/5] %s (%d genes)... ", i, pname, length(pgenes)))
  t_start <- proc.time()

  null_distributions[[pname]] <- approximate_null_distribution(
    expr_gene = gene_expr,
    inds      = relapse_inds,
    geneset   = pgenes,
    iters     = N_ITERS,
    p         = 1
  )

  elapsed <- round((proc.time() - t_start)[["elapsed"]], 1)
  cat("done in", elapsed, "seconds.\n")
}

# Save to disk — never run this again
saveRDS(null_distributions, "null_distributions.RDS")
cat("\n✓ Null distributions saved to null_distributions.RDS\n")
```

To load the saved results in a future session:
```r
null_distributions <- readRDS("null_distributions.RDS")
```

---

## Step 8 — Assign Deliverables and Close (5 min)

Go around the table and confirm each person's deliverable for Meeting 2:

| Role | Person | Deliverable by Meeting 2 (Apr 26–27) |
|------|--------|---------------------------------------|
| Report Lead | A | Report skeleton `.Rmd` with header, all 6 sections stubbed out |
| Hotelling Analyst | B | T² test + Q-Q plots for all 5 pathways; decision on assumption handling |
| GSEA Analyst | C | ES curves + empirical p-values for all 5 pathways using saved null distributions |
| Interpretation Lead | D | 2–3 sentence biology background for each of the 5 pathways |

**Set the Meeting 2 date before leaving.**

---

## Troubleshooting

### `solve(S)` fails — singular covariance matrix in Hotelling's T²

This happens when `p > n1 + n2 - 2` (more genes than effective observations). It is common for large pathways.

**Remedies:**
```r
# Option A: Check if p > n1 + n2 - 2
p_dim  <- length(test_genes)
n_total <- nrow(data_relapse) + nrow(data_norelapse)
cat("p =", p_dim, " | n1+n2-2 =", n_total - 2, "\n")
cat("Singular?", p_dim >= n_total - 2, "\n")

# Option B: Use a permutation-based Hotelling's T²
# Instead of the F-distribution p-value, permute group labels and compute empirical p-value
hotelling_permutation_pvalue <- function(data_full, group_inds, n_perm = 1000) {
  obs_result <- two_sample_test(data_full[ group_inds, ],
                                data_full[!group_inds, ])
  obs_T2     <- obs_result$T2

  perm_T2 <- replicate(n_perm, {
    perm_inds <- sample(group_inds)
    r <- two_sample_test(data_full[ perm_inds, ],
                         data_full[!perm_inds, ])
    r$T2
  })
  mean(perm_T2 >= obs_T2)  # empirical p-value
}

# Option C: Use only the top-k genes by variance to reduce p
top_k     <- 50  # keep the 50 most variable genes in the pathway
gene_vars <- apply(subset_mat, 2, var)
top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:min(top_k, ncol(subset_mat))]
subset_reduced <- subset_mat[, top_genes]
# then re-run two_sample_test on the reduced matrix
```

### `approximate_null_distribution` is too slow

```r
# Option A: Use the parallel package to run all 5 pathways simultaneously
library(parallel)
n_cores <- max(1, detectCores() - 1)
cat("Using", n_cores, "cores\n")

null_distributions <- mclapply(
  pathway_names,
  function(pname) {
    pgenes <- get_pathway_genes(pname, gene_expr, HML_data)
    approximate_null_distribution(gene_expr, relapse_inds, pgenes, iters = N_ITERS)
  },
  mc.cores = n_cores
)
names(null_distributions) <- pathway_names
saveRDS(null_distributions, "null_distributions.RDS")

# Option B: Run on ucloud — upload gene_expressions.RDS, group.RDS, HML.RDS,
# GSEA_functions.R, and Meeting1_Setup.R, then run just the null distribution loop.
```

### Group label not matching "relapse"

```r
# If relapse_inds is all FALSE or all TRUE, your group labels may use different strings
cat("Unique values in group:\n"); print(unique(group))
cat("Full table:\n"); print(table(group))

# Adjust accordingly, e.g.:
# relapse_inds <- group == 1          # if group is 0/1 integer
# relapse_inds <- group == TRUE       # if group is already logical
# relapse_inds <- as.logical(group)   # if group is "TRUE"/"FALSE" character
```

---

## End of Meeting 1 Checklist

Before everyone leaves, confirm:

- [ ] All four members can source all scripts with no errors
- [ ] Data dimensions verified: 12,643 × 286 gene matrix, 286-length group vector
- [ ] Genetic lottery run together; full input/output saved somewhere shared
- [ ] `get_pathway_genes()` helper written and tested on all 5 pathways
- [ ] Gene counts per pathway (p values for Hotelling) recorded
- [ ] Sanity check passed: ES curve plotted for HALLMARK_GLYCOLYSIS
- [ ] Hotelling's T² ran for HALLMARK_GLYCOLYSIS (or singular issue identified)
- [ ] Null distribution computation started and running in background
- [ ] Meeting 2 date set

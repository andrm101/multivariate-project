# DSK805 — GSEA Analysis
# anman25 | mhaus25 | cesko25 | dabay25  (hash: 45503)
# Subteam: GSEA + Report Lead

#load repo
library(car)
source("C:/Users/andre/Desktop/Sandbox/Multivariate Project/hotelling.R")
source("C:/Users/andre/Desktop/Sandbox/Multivariate Project/GSEA_functions.R")
source("C:/Users/andre/Desktop/Sandbox/Multivariate Project/lottery.R")
source("C:/Users/andre/Desktop/Sandbox/Multivariate Project/dsk805utils.R")

#load Data
gene_expr <- readRDS("C:/Users/andre/Desktop/Sandbox/Multivariate Project/gene_expressions.RDS") # g × n: 12,643 genes × 286 patients
group <- readRDS("C:/Users/andre/Desktop/Sandbox/Multivariate Project/group.RDS") # "relapse" / "no-relapse", length 286
HML_data <- readRDS("C:/Users/andre/Desktop/Sandbox/Multivariate Project/HML.RDS") # 50 hallmark pathways

relapse <- group == "relapse"  # TRUE = relapse (n1=69), FALSE = no-relapse (n2=217)

# Lottery assignments (hash 45503):
pathways <- c(
  Set1 = "HALLMARK_APICAL_JUNCTION",
  Set2 = "HALLMARK_INFLAMMATORY_RESPONSE",
  Set3 = "HALLMARK_E2F_TARGETS",
  Set4 = "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
  Set5 = "HALLMARK_CHOLESTEROL_HOMEOSTASIS"
)

#Gene subsetting 
# Keep only genes present in both the pathway definition and the expression matrix.
# Dropped genes are ignored
pathway_genes <- lapply(pathways, function(pw) {
  intersect(unlist(HML_data[[pw]]), rownames(gene_expr))
})

p_dims <- sapply(pathway_genes, length)
cat("Gene counts per pathway (= p for Hotelling's T²):\n"); print(p_dims)

#GSEA Step 1: Rank all genes by Signal-to-Noise Ratio 
# Computed once and reused across all 5 pathways.
# SNR = (μ_relapse − μ_no-relapse) / (σ_relapse + σ_no-relapse)
ranked <- rank_genes(gene_expr, relapse)

#GSEA Step 2: Enrichment scores 
es_all <- lapply(pathway_genes, function(g) {
  compute_enrichment_score(ranked$diff, g, p = 1)
})
obs_ES <- sapply(es_all, function(r) r$ES_max)
cat("\nObserved Enrichment Scores:\n"); print(round(obs_ES, 4))

#Figure 1: ES curves
pdf("fig_ES_curves.pdf", width = 10, height = 7)
par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
for (i in seq_along(pathways)) {
  col_i <- ifelse(obs_ES[i] > 0, "firebrick", "steelblue")
  plot(es_all[[i]]$ES, type = "l", lwd = 2, col = col_i,
       main  = gsub("HALLMARK_", "", pathways[i]),
       xlab  = "Gene rank (high SNR → low SNR)",
       ylab  = "Running Enrichment Score",
       cex.main = 0.9)
  abline(h = 0, lty = 2, col = "grey60")
  abline(v = which.max(abs(es_all[[i]]$ES)), lty = 3, col = "orange", lwd = 1.2)
  rug(which(names(ranked$diff) %in% pathway_genes[[i]]),
      col = adjustcolor(col_i, 0.3), ticksize = 0.04)
  legend("topright", bty = "n", cex = 0.8,
         legend = sprintf("ES = %+.4f", obs_ES[i]))
}
dev.off()
cat("Saved: fig_ES_curves.pdf\n")

#GSEA Step 3: Null distributions (permutation test)
# H0: pathway S is not differentially expressed between the two groups.
# Permute group labels 1000 times and recompute ES each time.
# Empirical p-value = fraction of null |ES| >= |observed ES| (two-sided).

N_PERM <- 1000

if (file.exists("C:/Users/andre/Desktop/Sandbox/Multivariate Project/null_dists.RDS")) {
  null_dists <- readRDS("C:/Users/andre/Desktop/Sandbox/Multivariate Project/null_dists.RDS")
  cat("Loaded saved null distributions.\n")
} else {
  cat("Computing null distributions (", N_PERM, "permutations × 5 pathways)...\n")
  null_dists <- setNames(vector("list", length(pathways)), names(pathways))
  for (i in seq_along(pathways)) {
    cat(sprintf("  [%d/5] %s... ", i, gsub("HALLMARK_", "", pathways[i])))
    t0 <- proc.time()
    null_dists[[i]] <- approximate_null_distribution(
      gene_expr, relapse, pathway_genes[[i]], iters = N_PERM
    )
    cat(sprintf("done (%.0f sec)\n", (proc.time() - t0)[["elapsed"]]))
  }
  saveRDS(null_dists, "null_dists.RDS")
}

# Empirical p-values (two-sided)
emp_pvals <- mapply(function(obs, null) mean(abs(null) >= abs(obs)),
                    obs_ES, null_dists)
cat("\nEmpirical p-values (two-sided, N =", N_PERM, "):\n")
print(data.frame(ES = round(obs_ES, 4), p_value = round(emp_pvals, 4)))

#Figure 2: Null distribution histograms (report figure)
pdf("fig_null_dists.pdf", width = 10, height = 7)
par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
for (i in seq_along(pathways)) {
  hist(null_dists[[i]], breaks = 35, col = "lightgrey", border = "white",
       main = gsub("HALLMARK_", "", pathways[i]),
       xlab = "ES under H0", cex.main = 0.9)
  abline(v =  obs_ES[i], col = "firebrick", lwd = 2)
  abline(v = -obs_ES[i], col = "firebrick", lwd = 1.2, lty = 2)
  legend("topright", bty = "n", cex = 0.78,
         legend = c(sprintf("ES = %+.4f", obs_ES[i]),
                    sprintf("p  =  %.4f",  emp_pvals[i])))
}
dev.off()
cat("Saved: fig_null_dists.pdf\n")

# Normality assessment (Lecture 6 Approach 2: PCA + Shapiro-Wilk) 
# For Hotelling's T² to be valid, data within each group should be MVN.
# Approach: project each pathway subset onto its principal components,
# then apply Shapiro-Wilk per PC — PCs are uncorrelated (independent under MVN),
# which avoids the overlapping-tests problem of testing raw genes one by one.

assess_normality <- function(genes, expr, group_idx, pathway_label) {
  X <- t(expr[genes, ]) # n × p (patients × genes)
  pca <- prcomp(X, center = TRUE, scale. = FALSE)
  varexp <- pca$sdev^2 / sum(pca$sdev^2)
  k80 <- which(cumsum(varexp) >= 0.80)[1] # PCs needed for ≥80% variance
  scores <- pca$x[, 1:k80, drop = FALSE]

  # Shapiro-Wilk per PC, per group
  sw_r <- apply(scores[ group_idx, , drop=FALSE], 2, function(x) shapiro.test(x)$p.value)
  sw_n <- apply(scores[!group_idx, , drop=FALSE], 2, function(x) shapiro.test(x)$p.value)

  any_fail <- any(c(sw_r, sw_n) < 0.05)
  cat(sprintf("  %-30s  k=%2d PCs  normality: %s\n",
              pathway_label, k80, ifelse(any_fail, "VIOLATED", "plausible")))

  list(pca = pca, k80 = k80, scores = scores, varexp = varexp,
       sw_relapse = sw_r, sw_norelapse = sw_n, X = X,
       violated = any_fail)
}

cat("\nNormality assessment (PCA + Shapiro-Wilk on first k PCs, α = 0.05):\n")
norm_res <- setNames(
  mapply(assess_normality,
         genes         = pathway_genes,
         pathway_label = gsub("HALLMARK_", "", pathways),
         MoreArgs      = list(expr = gene_expr, group_idx = relapse),
         SIMPLIFY      = FALSE),
  names(pathways)
)

#Figure 3: QQ plots for first PC, per group 
pdf("fig_QQ_plots.pdf", width = 12, height = 8)
par(mfrow = c(2, 5), mar = c(4, 4, 3, 1))
for (i in seq_along(pathways)) {
  nr <- norm_res[[i]]
  label <- gsub("HALLMARK_", "", pathways[i])
  for (grp in list(list(idx = relapse, col = "firebrick", name = "Relapse"),
                   list(idx = !relapse, col = "steelblue", name = "No-relapse"))) {
    qqnorm(nr$scores[grp$idx, 1], pch = 16, cex = 0.6, col = grp$col,
           main = sprintf("%s\nPC1 — %s", label, grp$name), cex.main = 0.75)
    qqline(nr$scores[grp$idx, 1], col = "black", lwd = 1.5)
    sw_p <- shapiro.test(nr$scores[grp$idx, 1])$p.value
    legend("topleft", bty = "n", cex = 0.75,
           legend = sprintf("SW p = %.4f", sw_p))
  }
}
dev.off()
cat("Saved: fig_QQ_plots.pdf\n")

# Yeo-Johnson transform (apply if normality is violated)
# yj_transform() from dsk805utils.R finds the optimal λ per gene.
# Re-run normality check after transforming — if still violated, document and
# proceed with Hotelling results interpreted with appropriate caution.

violated_sets <- names(which(sapply(norm_res, function(r) r$violated)))
if (length(violated_sets) > 0) {
  cat("\nApplying Yeo-Johnson to violated pathways:", paste(violated_sets, collapse=", "), "\n")
  yj_norm <- setNames(vector("list", length(violated_sets)), violated_sets)

  for (key in violated_sets) {
    genes <- pathway_genes[[key]]
    X_raw <- t(gene_expr[genes, ])
    X_yj <- tryCatch(apply(X_raw, 2, yj_transform),
                      error = function(e) { cat("  YJ failed for", key, "\n"); NULL })
    if (!is.null(X_yj)) {
      yj_norm[[key]] <- assess_normality(
        genes = genes,
        expr = gene_expr,
        group_idx = relapse,
        pathway_label = paste0(gsub("HALLMARK_", "", pathways[key]), " [YJ]")
      )
    }
  }
} else {
  cat("\nNo pathways require Yeo-Johnson transformation.\n")
  yj_norm <- list()
}

# ── Results summary table ─────────────────────────────────────────────────────
results <- data.frame(
  Pathway = gsub("HALLMARK_", "", pathways),
  p_genes = p_dims,
  ES = round(obs_ES, 4),
  Direction = ifelse(obs_ES > 0, "Relapse", "No-relapse"),
  GSEA_p = round(emp_pvals, 4),
  Sig = ifelse(emp_pvals < 0.05, "*", ""),
  Normality = sapply(norm_res, function(r) ifelse(r$violated, "Violated", "Plausible")),
  stringsAsFactors = FALSE
)
rownames(results) <- paste0("Set", 1:5)

cat("\n"); print(results)

# Save report
saveRDS(list(
  pathway_names = pathways,
  pathway_genes = pathway_genes,
  p_dims = p_dims,
  ranked = ranked,
  obs_ES = obs_ES,
  emp_pvals = emp_pvals,
  norm_res = norm_res,
  results = results
), "gsea_results.RDS")
cat("Saved: gsea_results.RDS\n")

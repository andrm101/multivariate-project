# DSK805 — Team Meeting Outline
## 4-Person Team | 3 Meetings | Deadline: April 30

---

## Team Roles

Divide responsibilities by **method, not by pathway** (as the guidelines recommend — pathway-by-pathway splits mean everyone rewrites the same setup code).

| # | Role | Owner |
|---|------|-------|
| 1 | **Report Lead** | Person A |
| 2 | **Hotelling Analyst** | Person B |
| 3 | **GSEA Analyst** | Person C |
| 4 | **Interpretation & Viz Lead** | Person D |

---

## Meeting 1 — Kickoff & Setup
**When:** April 23–24 (Day 1–2) | **Duration:** ~1.5 hours

### Goal
Get everyone unblocked and running. No results needed yet — just a clean, shared foundation.

### Agenda

**1. Run `genetic_lottery` together (15 min)**
- One person runs the lottery live so everyone sees the output.
- Record the full input/output — this goes directly into the report header.
- Confirm the 5 assigned pathways:
  - `HALLMARK_PANCREAS_BETA_CELLS`
  - `HALLMARK_NOTCH_SIGNALING`
  - `HALLMARK_INFLAMMATORY_RESPONSE`
  - `HALLMARK_TNFA_SIGNALING_VIA_NFKB`
  - `HALLMARK_GLYCOLYSIS`

**2. Verify everyone's R environment (20 min)**
- All four members source and run the provided scripts:
  ```r
  source('hotelling.R')
  source('GSEA_functions.R')
  source('lottery.R')
  source('dks805utils.R')  # requires the `car` package
  ```
- Load all three `.RDS` files and confirm dimensions:
  - `gene_expressions.RDS` → 12,643 rows × 286 columns
  - `group.RDS` → 286-length binary vector (69 relapse, 217 no-relapse)
  - `HML.RDS` → named list of 50 hallmark pathways

**3. Handle gene overlap issue (10 min)**
- Decide on subsetting logic: for each pathway, intersect the gene names with rows in `gene_expressions.RDS` and use only those.
- Person B (Hotelling) needs this for `p` (dimension in T²).
- Person C (GSEA) needs this as the gene set `S` passed to enrichment functions.
- Write a shared 3-line helper function for this so both analysts use identical subsets.

**4. Assign deliverables per role (15 min)**

| Role | By Meeting 2 |
|------|-------------|
| **Report Lead (A)** | Draft the report skeleton in `.Rmd` (header, intro placeholder, sections 1–6) |
| **Hotelling (B)** | T² test + normality diagnostics for all 5 pathways |
| **GSEA (C)** | GSEA pipeline running; null distribution started (≥1,000 iterations minimum) |
| **Interpretation (D)** | 1–2 paragraph biological background for each of the 5 pathways |

**5. Logistics (10 min)**
- Agree on a shared folder or Git repo for code and the `.Rmd` report file.
- Decide who handles ucloud if the null distribution computation is too slow locally.
- Set Meeting 2 date (target: April 26–27).

---

## Meeting 2 — Progress Check & Integration
**When:** April 26–27 (Day 4–5) | **Duration:** ~2 hours

### Goal
Review all partial results, catch assumption violations early, and start writing the findings.

### Agenda

**1. Hotelling T² — results review (30 min) — Person B presents**
- Show the T² test statistic and p-value for each of the 5 pathways.
- Walk through normality diagnostics: Q-Q plots and/or Mardia's skewness/kurtosis test.
- Decide together which of the 3 remedies to apply per pathway if assumptions are violated:
  - a) Apply Yeo-Johnson transform via `dsk805utils.R`
  - b) Use a permutation version of Hotelling's T²
  - c) Document the violation clearly and account for it in the interpretation
- **Note:** All three options are valid. Option c) is better than ignoring it. Be explicit.

**2. GSEA — progress check (30 min) — Person C presents**
- Confirm `rank_genes` output looks right: genes ranked by signal-to-noise ratio descending.
- Show enrichment score curves (running-sum walk) for at least 1–2 pathways.
- Check if null distribution sampling has reached 1,000 iterations. If not, discuss ucloud or parallelization.
- Compute observed ES vs. null distribution → empirical p-value for each pathway.

**3. Pathway backgrounds (15 min) — Person D presents**
- Brief preview of what each pathway does biologically — even 2–3 sentences per pathway is enough.
- Flag any known literature linking that pathway to breast cancer metastasis (this enriches the findings section).
- Good pathways to look up specifically:
  - **TNFA_SIGNALING_VIA_NFKB** — well-documented pro-inflammatory, pro-metastatic role in breast cancer.
  - **INFLAMMATORY_RESPONSE** — broad but highly relevant to tumor microenvironment.
  - **GLYCOLYSIS** — Warburg effect; metabolic reprogramming is a known cancer hallmark.
  - **NOTCH_SIGNALING** — cell fate, known role in cancer stem cells and metastasis.
  - **PANCREAS_BETA_CELLS** — less obvious connection to breast cancer; interesting to discuss if enriched or not.

**4. First draft of findings (30 min) — whole team**
- Fill in the findings section of the report together. Contrast GSEA vs. Hotelling's T² per pathway:
  - Do both methods agree on which pathways are differentially expressed?
  - If they disagree, why might that be? (GSEA uses the full ranked list; Hotelling's T² is multivariate but dimensionality-limited.)
- Start writing while the results are fresh — don't leave this for Meeting 3.

**5. Report review (15 min) — Person A leads**
- Check page count: you have max 6 pages.
- Verify all required sections are present:
  1. Project header (with lottery output)
  2. Generative AI statement (separate file)
  3. Introduction
  4. Methods description
  5. Assumption checks (with plots)
  6. Findings
- Assign remaining writing tasks for each person before Meeting 3.

---

## Meeting 3 — Final Assembly & Submission
**When:** April 29 (Day 7) | **Duration:** ~1.5 hours

### Goal
Polish, verify, and submit. No new analysis — only writing, checking, and formatting.

### Agenda

**1. Compile the full report (30 min) — Person A drives**
- Knit the `.Rmd` to PDF. Fix any LaTeX errors.
- Verify the project header renders correctly (genetic lottery output must be shown verbatim).
- Make sure figures are properly captioned and labeled.
- Confirm the page limit: cut or move to appendix anything that overruns 6 pages.

**2. Cross-check all results (20 min) — Persons B and C**
- Person B re-runs the Hotelling T² on the final, confirmed subsets to make sure nothing changed.
- Person C confirms enrichment scores and p-values are computed on the same gene subsets as Hotelling's T².
- Verify both analysts used the same pathway intersections (the shared helper function from Meeting 1).

**3. Final read-through (20 min) — Person D leads**
- Read the entire report aloud or silently, checking for:
  - Clarity: would another DSK805 student understand this without your help?
  - Consistency between methods and findings sections.
  - Any missing citations (the three required references are in the assignment).
- Check the Generative AI form is filled (optional but recommended) and is a separate `.docx`.

**4. Submission (10 min) — Person A submits**
- Upload to itslearning: report as `.pdf`, code as a separate file, GenAI form as `.docx`.
- Everyone confirms they can see the submission before leaving the meeting.

---

## Tips for an Excellent Case Study

### Statistical rigor
- **Don't skip assumption checking.** The grading weights Content and Accuracy equally at 40% each. Stating "we checked normality with Q-Q plots and found moderate deviation in pathway X, so we applied Yeo-Johnson transform" is far better than a silent pass.
- **Report exact p-values, not just ≥/≤ thresholds.** For both Hotelling and GSEA, show the numerical p-value alongside your decision.
- **Use at least 1,000 null permutations for GSEA.** Fewer iterations produce unstable p-values. If you have time, 5,000+ is even better — and this is worth the compute cost.

### Structure and clarity (20% of grade)
- **Write for a fellow DSK805 student.** Don't assume the reader knows GSEA; give a one-paragraph intuition of what a high enrichment score means.
- **Contrast the methods explicitly.** The assignment's centerpiece is comparing GSEA vs. Hotelling's T². Don't just list results side-by-side — explain *why* they might agree or disagree. GSEA considers the full ranked list; Hotelling's T² is sensitive to dimensionality and distributional assumptions.
- **Use figures strategically.** A running-sum (enrichment) plot per pathway and Q-Q plots for normality checks cover the assumption and findings requirements without ballooning the page count.

### Biology (the easy 10% that most teams skip)
- You don't need to be a molecular biologist, but one paragraph connecting your statistical findings to known biology impresses. For example: if TNFA_SIGNALING_VIA_NFKB is significantly enriched in the relapse group, note that NF-κB is a well-known driver of breast cancer invasiveness — this makes the finding interpretable.
- Use AI (per the assignment's own suggestion) to help explain the biology of the five pathways in plain language. Cite any sources you use.

### Pitfalls to avoid
- **Don't divide work pathway-by-pathway.** The setup code is identical across all 5 pathways — duplicating it wastes time and introduces inconsistencies.
- **Don't skip the project header.** It is generated by `genetic_lottery()` and is part of the grading. Tampering with the function is flagged as exam fraud.
- **Don't ignore violated assumptions.** A sentence saying "assumptions are violated, we ignore this" is penalized. Acknowledging and accounting for them (even via option c) is rewarded.
- **Don't run too few null samples.** Fewer than 1,000 makes your empirical p-values unreliable. Start the computation early (at Meeting 1 or right after) — it is the main bottleneck.
- **Start the null distribution computation immediately after Meeting 1.** With 5 pathways × 1,000+ permutations each, this is your longest compute job. Don't leave it for the night before the deadline.

---

## Quick Timeline Summary

| Date | Milestone |
|------|-----------|
| Apr 23–24 | **Meeting 1:** Setup, lottery run, roles assigned, computation starts |
| Apr 24–26 | Person B: Hotelling + diagnostics / Person C: GSEA running / Person D: literature |
| Apr 26–27 | **Meeting 2:** Review results, draft findings together |
| Apr 27–29 | Finish writing, polish figures, assemble report |
| Apr 29 | **Meeting 3:** Final assembly, cross-check, submit |
| **Apr 30** | **Deadline on itslearning** |

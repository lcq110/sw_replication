# SW replication files

Shang and Woźniak (2026), *Identification Verification for Structural Vector Autoregressions with Sparse Heterogeneous Markov Switching Heteroskedasticity*, [doi:10.48550/arXiv.2603.16035](https://doi.org/10.48550/arXiv.2603.16035).

Replication scripts and inputs for this paper. `SOURCES.md` records the source version and local adaptations; `VERSIONS.txt` lists the validated package versions. The GPL license is retained.

Install `bsvars`, `Rcpp`, `RcppArmadillo`, `RcppProgress`, `RcppTN`, `permute`, `mgcv`, `xtable`, `mvtnorm`, and `zoo`. Run the commands below from the repository root. The package versions used for validation are recorded in `VERSIONS.txt`.

## Empirical application and figures

1. `data/` contains the unchanged author `bpss21.rda` and exogenous regime labels. `estimations/model.R` applies the author's variable transformations and ten-lag specifications.
2. Estimate the four baseline models:

   ```r
   for (model in c("hmsh4", "msh10", "exh", "sv"))
     system2(file.path(R.home("bin"), "Rscript"), c("estimations/estimate.R", model))
   ```

   Each fit saves `estimations/results/e_<model>.rda`. Defaults keep the author's first-stage burn-in of 1,000,000 draws and 50,000 posterior draws. HMSH/MSH/EXH also run one 200,000-draw Step 2 pass; its historical repetition count is not recorded.
3. Choose the monetary-policy reference row explicitly for each fitted posterior. The paper selects the shock that contributes most to the impact forecast error variance of the federal funds rate. The author's original rows were HMSH=8, MSH=3, EXH=8, SV=10; they are not automatically valid after a new estimation. Check the row before calling, for example:

   ```r
   system2(file.path(R.home("bin"), "Rscript"),
           c("estimations/normalise.R", "hmsh4", "8"))
   ```

   Repeat for the other fits with their rows. In every draw, the shock closest to the reference is placed in row 1 and signed to match it (the reference has a positive interest-rate coefficient). This produces `_norm.rda`, `_post.rda` (SDDR and conditional standard deviations), and `_irfs.rda` (60-month responses). The local first-row alignment adapter fixes the row-index mismatch in the author's C++ routine and permutes the shock-specific state fields consistently; see `SOURCES.md`.
4. Before assembling `tab:sddr`, run prior cases 1 and 3 from the author
   sensitivity scripts (s=1 and transition-prior diagonal 10):

   ```r
   for (model in c("hmsh4", "msh10")) for (prior_case in c(1, 3))
     system2(file.path(R.home("bin"), "Rscript"),
             c("estimations/estimate.R", model, prior_case))
   ```

   Check the reference row separately for each posterior, then call, for
   example, `Rscript estimations/normalise.R hmsh4_1 8`. The old author's
   row numbers are not guaranteed for a new fit.
5. Generate the figure layouts and the complete baseline/sensitivity table:

   ```r
   source("estimations/out_csds.R")
   source("estimations/out_irfs.R")
   source("estimations/out_sddr.R")
   ```

   Outputs: `estimations/out_csds.pdf`, `out_irfs1.pdf`, `out_irfs2.pdf`, and
   `out_sddr.tex`, matching `fig:csds`, `fig:irfs1`, `fig:irfs2`, `tab:sddr`.

## Monte Carlo experiments (manuscript Section 3)

1. Set the working directory to `simulations/`, create `dgps/` and `results/`, and run `dgps260.R` and `dgps780.R`. They generate 100 datasets for each of 11 volatility scenarios at each sample size, plus the true B matrix and standard deviations.
2. Run the eight `bsvar_{hmsh2,msh2,exh,sv}_{260,780}.R` files for dataset IDs 1:1100. Each script receives one ID as its command-line argument; for example:

   ```r
   system2(file.path(R.home("bin"), "Rscript"), c("bsvar_hmsh2_260.R", "1"))
   ```

   These are the author's ordinary R scripts with their own C++ kernels, 10,000 burn-in and 10,000 posterior draws. A complete sequential run is a substantial research computation; the files supplied here do not require a cluster scheduler.
3. After all datasets and models finish, run `bsvar_ms_wrap.R` to collect the results, then `Rscript tables.R > results/tables.txt`. The collection and table code are split from the author's single original file; the statistics are unchanged.
4. `results/tables.txt` contains the original LaTeX tables and printed RMSE/rejection-rate summaries (`tab:mceffB`, `tab:mcass`). It requires the full study of 100 datasets per scenario; a short API test cannot populate the paper's tables.

## Forecasting study (manuscript Section 5)

1. From the repository root, enter `forecasting/` and run the author's full-sample burn-in scripts:

   ```r
   setwd("forecasting")
   source("burnin/burnin_hmsh.R")
   source("burnin/burnin_msh.R")
   ```

   They retain one million burn-in draws per model and save the last draws
   under `burnin/starting_values/`. Each forecast script truncates these
   states to the current estimation window, as in the author implementation.
2. Run the 28 ordinary author R files sequentially, preserving their seeds,
   117 origins and 5,000 burn-in/10,000 posterior draws:

   ```r
   for (id in 1:28)
     system2(file.path(R.home("bin"), "Rscript"), sprintf("bsvars%02d.R", id))
   setwd("..")
   for (model in c("sv", "exh"))
     system2(file.path(R.home("bin"), "Rscript"),
             c("forecasting/baselines.R", model))
   source("forecasting/score.R")
   ```

   `baselines.R` uses the same plain warm-start/expanding-window procedure
   for SV and EXH, p=10, Gaussian disturbances, and prediction-project EXH
   prior sigma_s=1/sigma_nu=3. Full-sample initialization uses seed 20241110;
   the existing window seed formula is extended with model IDs 29/30.
   Data values/order and regime labels match the project's prediction input.
   The new 30-model scores go to `forecasting/outputs/origin_scores-rerun.csv`;
   the frozen `inputs/origin_scores.csv` is preserved. Numerical equality to
   the published baseline scores has not been established.
3. Regenerate the complete paper figure from the frozen 30-model scores:

   ```r
   setwd("forecasting")
   source("make_forecasting_results.R")
   ```

   Outputs include `outputs/figures/lps_by_regimes.pdf`, `.png`,
   `outputs/tables/tab_forecast_metrics.tex` and
   `outputs/results/model_summary.csv`. This reproduces `fig:lps_regimes`
   from existing scores. After a complete rerun, generate its figure and
   table separately:

   ```r
   system2(file.path(R.home("bin"), "Rscript"),
           c("make_forecasting_results.R", "outputs/origin_scores-rerun.csv",
             "outputs/rerun"))
   ```

   The optional arguments select the score input/output directory; figure
   calculations are unchanged. No ECS, SLURM, deployment or parallel
   scheduling code is included.

The source data and scores end in October 2025, whereas the manuscript text
says September; this discrepancy needs author confirmation. Short execution
checks do not establish convergence or reproduce the paper's numbers.
The machine-dependent Section 4 speed comparison is not included.

# Type-I error of diff_methylsig() under the null, for the degrees of freedom
# of its t approximation.
#
# diff_methylsig() compares the likelihood ratio statistic D to a squared t
# distribution. The paper (Park et al. 2014, doi:10.1093/bioinformatics/btu339,
# section 2.3) uses p = sum_k H(k - i)(J_k1 + J_k2 - 2) degrees of freedom.
# Since 2015 (commit 5896c23) the code uses p + 2, because the t approximation
# was too conservative with small samples. This script estimates the false
# positive rate of both, and of the chi-square approximation, when there is no
# difference between groups.
#
# The null data are resampled from real data, so that methylation levels,
# dispersion, coverage, missingness, and the spacing of CpGs are realistic:
#
# * Loci and positions: CpGs of bsseqData's BS.cancer.ex with coverage 5-500
#   in at least 2 samples per group.
# * True methylation per locus: the methylation of all 6 samples pooled, so
#   the groups have no difference.
# * True dispersion per locus: diff_methylsig()'s estimate on the real data
#   (cancer vs normal, no local information), which is the dispersion around
#   the group means.
# * Coverage of each simulated sample: the coverage of a random real sample.
# * Methylated reads: beta-binomial with the true methylation and dispersion.
#
# Each setting is tested with diff_methylsig(), which returns D and p + 2.
# Rates are reported for all tested loci, for loci with true methylation
# between 0.1 and 0.9, and for overdispersed loci (true theta below the
# 1e6 upper bound of diff_methylsig()'s estimate, i.e. not binomial).
#
# Usage, from the repository root, in the devel image with methylSig installed:
#   Rscript simulations/df_calibration.R [output directory]
# The QQ-plot data go in <output directory>/null_pvalues.

suppressPackageStartupMessages({
    library(methylSig)
    library(bsseq)
})

args = commandArgs(trailingOnly = TRUE)
out_dir = if (length(args) > 0) args[1] else 'simulations'

# Loci of BS.cancer.ex to start from (before filtering). N_LOCI overrides it,
# e.g. for a quick test
n_loci = as.integer(Sys.getenv('N_LOCI', 250000))
samples_per_group = c(2, 3, 4, 5, 6, 8, 10)
local_settings = list(
    'none' = list(local_window_size = 0, local_disp = TRUE, local_meth = TRUE),
    'disp' = list(local_window_size = 200, local_disp = TRUE, local_meth = FALSE),
    'disp_meth' = list(local_window_size = 200, local_disp = TRUE, local_meth = TRUE))
alphas = c(0.05, 0.01, 0.001)
# diff_methylsig() forks a worker per core with mclapply(), and each worker
# can end up copying the parent's memory, so keep this small in Docker (8 GB)
n_cores = as.integer(Sys.getenv('N_CORES', 4))
dir.create(file.path(out_dir, 'null_pvalues'), recursive = TRUE, showWarnings = FALSE)

set.seed(20260926)

#####################################
# Truth from the real data

data(BS.cancer.ex, package = 'bsseqData')
real = BS.cancer.ex[seq(n_loci)]
rm(BS.cancer.ex)
real = filter_loci_by_coverage(real, min_count = 5, max_count = 500)
real = suppressMessages(filter_loci_by_group_coverage(
    real, 'Type', c('cancer' = 2, 'normal' = 2)))

real_fit = suppressMessages(diff_methylsig(
    bs = real,
    group_column = 'Type',
    comparison_groups = c('case' = 'cancer', 'control' = 'normal'),
    disp_groups = c('case' = TRUE, 'control' = TRUE),
    local_window_size = 0,
    n_cores = n_cores))

# Keep the loci that could be tested on the real data
real = real[match(granges(real_fit), granges(real))]
real_cov = as.matrix(getCoverage(real, type = 'Cov'))
real_meth = as.matrix(getCoverage(real, type = 'M'))
true_meth = rowSums(real_meth) / rowSums(real_cov)
true_theta = real_fit$disp_est
rm(real_fit)
invisible(gc())

cat(sprintf('%s loci. True methylation: %.0f%% at 0 or 1, %.0f%% between 0.1 and 0.9. True dispersion (theta): median %.1f, %.0f%% at the maximum (binomial).\n',
    length(real), 100 * mean(true_meth %in% c(0, 1)),
    100 * mean(true_meth > 0.1 & true_meth < 0.9), median(true_theta),
    100 * mean(true_theta >= 1e6)))

#####################################
# Null data and tests

simulate_null = function(n) {
    n_samples = 2 * n
    source_samples = sample(ncol(real_cov), n_samples, replace = TRUE)
    cov = unname(real_cov[, source_samples, drop = FALSE])

    # Beta-binomial with mean true_meth and theta = alpha + beta
    alpha = true_meth * true_theta
    beta = (1 - true_meth) * true_theta
    p = matrix(stats::rbeta(length(cov), alpha, beta), nrow = nrow(cov))
    meth = matrix(stats::rbinom(length(cov), cov, p), nrow = nrow(cov))

    sample_names = paste0('s', seq(n_samples))
    bs = BSseq(
        Cov = cov,
        M = meth,
        gr = granges(real),
        pData = data.frame(
            group = rep(c('a', 'b'), each = n),
            row.names = sample_names),
        sampleNames = sample_names)

    # As a user would, require at least 2 covered samples per group
    suppressMessages(filter_loci_by_group_coverage(
        bs, 'group', c('a' = 2, 'b' = 2)))
}

results = list()
for (n in samples_per_group) {
    bs = simulate_null(n)
    bs_idx = match(granges(bs), granges(real))

    for (setting in names(local_settings)) {
        started = Sys.time()
        fit = suppressMessages(do.call(diff_methylsig, c(list(
            bs = bs,
            group_column = 'group',
            comparison_groups = c('case' = 'a', 'control' = 'b'),
            disp_groups = c('case' = TRUE, 'control' = TRUE),
            t_approx = TRUE,
            n_cores = n_cores), local_settings[[setting]])))
        tested = bs_idx[match(granges(fit), granges(bs))]
        tested_mid = true_meth[tested] > 0.1 & true_meth[tested] < 0.9
        tested_overdispersed = true_theta[tested] < 1e6

        D = pmax(fit$log_lik_ratio, 0)
        pvalues = list(
            't_df_paper' = 2 * stats::pt(-sqrt(D), fit$df - 2),
            't_df_code' = 2 * stats::pt(-sqrt(D), fit$df),
            'chisq' = stats::pchisq(D, 1, lower.tail = FALSE))

        for (test in names(pvalues)) {
            for (subset in c('all', 'meth_10_90', 'overdispersed')) {
                p = pvalues[[test]]
                if (subset == 'meth_10_90') p = p[tested_mid]
                if (subset == 'overdispersed') p = p[tested_overdispersed]
                for (alpha in alphas) {
                    results[[length(results) + 1]] = data.frame(
                        samples_per_group = n,
                        local = setting,
                        test = test,
                        loci = subset,
                        n_loci = length(p),
                        median_df_paper = stats::median(fit$df - 2),
                        alpha = alpha,
                        false_positives = sum(p < alpha),
                        rate = mean(p < alpha),
                        rate_over_alpha = mean(p < alpha) / alpha)
                }
            }
        }

        # Keep the p-values for QQ-plots
        saveRDS(
            data.frame(D = D, df_code = fit$df, meth_10_90 = tested_mid,
                true_meth = true_meth[tested], true_theta = true_theta[tested]),
            file.path(out_dir, 'null_pvalues', sprintf('null_n%s_%s.rds', n, setting)))

        n_tested = length(tested)
        rm(fit, pvalues, D, tested, tested_mid, tested_overdispersed)
        invisible(gc())

        cat(sprintf('n = %2s per group, local = %-9s: %s loci, %.0f s\n',
            n, setting, n_tested,
            as.numeric(difftime(Sys.time(), started, units = 'secs'))))
    }
}

results = do.call(rbind, results)
utils::write.csv(results, file.path(out_dir, 'df_calibration_results.csv'), row.names = FALSE)

writeLines(
    c(sprintf('methylSig %s', utils::packageVersion('methylSig')), utils::capture.output(utils::sessionInfo())),
    file.path(out_dir, 'df_calibration_session.txt'))

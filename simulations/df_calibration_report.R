# Builds df_calibration.html from the results of df_calibration.R and
# df_calibration_template.html, by embedding the rates and thinned QQ-plot
# points as JSON.
#
# Usage, from the repository root, after df_calibration.R:
#   Rscript simulations/df_calibration_report.R [simulation directory]

args = commandArgs(trailingOnly = TRUE)
sim_dir = if (length(args) > 0) args[1] else 'simulations'

results = utils::read.csv(file.path(sim_dir, 'df_calibration_results.csv'))

# QQ-plot points: about 150 ranks evenly spaced in -log10(expected p), plus the
# 10 smallest p-values, for each sample size, setting, and test
qq = list()
for (f in list.files(file.path(sim_dir, 'null_pvalues'), pattern = '\\.rds$', full.names = TRUE)) {
    m = regmatches(basename(f), regexec('null_n(\\d+)_(.+)\\.rds', basename(f)))[[1]]
    d = readRDS(f)
    D = d$D
    pvalues = list(
        't_df_paper' = 2 * stats::pt(-sqrt(D), d$df_code - 2),
        't_df_code' = 2 * stats::pt(-sqrt(D), d$df_code),
        'chisq' = stats::pchisq(D, 1, lower.tail = FALSE))
    N = length(D)
    ranks = unique(sort(c(round(10^seq(0, log10(N), length.out = 150)), seq(10))))
    for (test in names(pvalues)) {
        p = sort(pvalues[[test]])
        qq[[length(qq) + 1]] = data.frame(
            n = as.integer(m[2]),
            local = m[3],
            test = test,
            exp = round(-log10((ranks - 0.5) / N), 3),
            obs = round(-log10(pmax(p[ranks], 1e-300)), 3))
    }
}
qq = do.call(rbind, qq)

data_json = jsonlite::toJSON(
    list(results = results, qq = qq),
    dataframe = 'rows', digits = NA, auto_unbox = TRUE)

template = readLines(file.path(sim_dir, 'df_calibration_template.html'), warn = FALSE)
template = paste(template, collapse = '\n')
stopifnot(grepl('/*__DATA__*/null', template, fixed = TRUE))
page = sub('/*__DATA__*/null', data_json, template, fixed = TRUE)
writeLines(page, file.path(sim_dir, 'df_calibration.html'))
cat(sprintf('Wrote %s (%.0f KB)\n', file.path(sim_dir, 'df_calibration.html'), nchar(page, type = 'bytes') / 1024))

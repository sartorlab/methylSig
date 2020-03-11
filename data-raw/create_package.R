# docker run --interactive --tty --rm --volume /Users/rcavalca/Projects:/Projects rcavalcante/bioconductor_docker:RELEASE_3_10

library(devtools)

# Description fields
description = list(
    Title = 'MethylSig: Differential Methylation Testing for WGBS and RRBS Data',
    Version = '0.99.0',
    Date = '2020-02-28',
    `Authors@R` = 'c(
        person(given = "Yongseok",
               family = "Park",
               role = c("aut"),
               email = "yongpark@pitt.edu"),
        person(given = "Raymond G.",
               family = "Cavalcante",
               role = c("aut", "cre"),
               email = "rcavalca@umich.edu"))',
    Description = 'MethylSig is a package for testing for differentially methylated
        cytosines (DMCs) or regions (DMRs) in whole-genome bisulfite sequencing
        (WGBS) or reduced representation bisulfite sequencing (RRBS) experiments.
        MethylSig uses a beta binomial model to test for significant differences
        between groups of samples. Several options exist for either site-specific
        or sliding window tests, and variance estimation.',
    BugReports = 'https://github.com/sartorlab/methylSig/issues',
    biocViews = 'DNAMethylation, DifferentialMethylation, Epigenetics, Regression, MethylSeq',
    License = 'GPL-3',
    Depends = 'R (>= 3.6)'
)

# Create package
path = '/Projects/methylSig'
create_package(path, fields = description)
activate_project(path)
# use_description(fields = description) # For updating

# Build ignore
build_ignore_files = c('README.md', '.travis.yml', '.git', '.gitignore')
use_build_ignore(files = build_ignore_files)

# Data
use_data_raw(name = '01-create_cov_files')
use_data_raw(name = '02-create_bsseq_rda')
use_data_raw(name = '03-create_internal_rda')

# Documentation
use_readme_md()
use_news_md()
use_package_doc()
use_vignette(name = 'using-methylSig', title = 'Using methylSig')

# Travis
use_travis()
use_travis_badge(ext = 'org')

# Coverage
use_coverage(type = 'coveralls')

# Testing
use_testthat()

# Package dependencies
use_package('bsseq', type = 'Imports')

# R files and test files
use_r('filter_loci_by_coverage')
use_r('filter_loci_by_snps')
use_r('tile_by_windows')
use_r('tile_by_regions')
use_r('filter_loci_by_group_coverage')
use_r('diff_binomial')
use_r('diff_methylsig')
use_r('diff_dss_fit')
use_r('diff_dss_test')
use_r('annotate_diff')
use_r('visualize_diff')
use_r('region_enrichment_diff')

use_test('filter_loci_by_coverage')
use_test('filter_loci_by_snps')
use_test('tile_by_windows')
use_test('tile_by_regions')
use_test('filter_loci_by_group_coverage')
use_test('diff_binomial')
use_test('diff_methylsig')
use_test('diff_dss_fit')
use_test('diff_dss_test')
use_test('annotate_diff')
use_test('visualize_diff')
use_test('region_enrichment_diff')

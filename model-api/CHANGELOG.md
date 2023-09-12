# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Adding...
- Multicollinearity tests

### Fixing
- Modularization of various utilities (presently bundled together in serializers.r)

## [0.7.0] - 2020-05-20
### Added
- Finite population correction
- Hierarchical decentering option
- Generalized Poisson for underdispersed counts
- More granular control opts for the two likelihoods in the hurdle model
- Contrasts code / derived quantities
- Options for deflections (effects) coding

### Fixed
- Censoring and truncation objects (`is.censored` and `censor.limit.vec`) and
  updated PPCs for models involving censored data

## [0.6.0] - 2019-05-29
### Added
- Support, where appropriate, for hierarchical variances specification
- Summaries that distinguish between sampled sites and all possible sites
- 'Superscripts' in folder dev/ for running jobs on a cluster

### Fixed
- Stratum-level summaries (no longer a mean of means)
- Default priors, esp. for models involving the inverse logit 'link'

## [0.5.0] - 2019-02-26
### Added
- Complete ZIP archive of code, data, config in the GCS bucket
- Trend by stratum, in addition to park-level trend
- Plots of coefficient estimates
- Somewhat better default figure dimensions
- New outputs for ordinal models (latent normal and latent beta)
- `g0` and `g0-g1` group level effects for the 'zeros' component of the beta
  hurdle model

### Fixed
- Inits caching -- inits are now cached _only_ when our `gelman_diag` < 2 or it beats the last benchmark
- Trimmed the file tree some to avoid too-long file paths
- Somewhat more appropriate variances structure tags in output directory names

## [0.4.2] - 2019-02-14 ("Happy Valentine's")
### Added
- Ordinal model for latent normal outcomes (e.g., soil stability)

### Fixed
- Proper DIC calculations for ordinal models with additional covariates

## [0.4.1] - 2018-11-15
### Added
- Improved inits handling (see [Issue 47](https://gitlab.com/nps-ima/uplands-ci/issues/47))

### Removed
- Bayesian $`R^2`$ (see [Issue 48](https://gitlab.com/nps-ima/uplands-ci/issues/48))

## [0.4.0] - 2018-09-12
### Added
- Bayesian R squared
- Status estimates conditioned on the time-varying values of covariates
- Zero-inflated binomial (inits for z need to be set at 0!)
- Proper posterior predictive checks for ordinal regressions
- Organization (subdirectories) for outputs
- A QA/QC review pipeline using an Azure VM and GCS:
  <https://gitlab.com/nps-ima/azure-vm>

### Fixed
- Default priors for coefficients associated with additional covariates in
  models involving the inverse-logit link
- Scaling of our `X.pred` JAGS object (using the mean and standard deviation
  applied to `X`)

## [0.3.1] - 2018-05-31
### Added
- Ability to declare the reference level for categorical
  covariates in the config files
- Somewhat more robust control flow for fits (via `while`, `tryCatch` and
  preemptive `return`'s in our rowwise fit routine)

### Fixed
- NPS-compatible plotting theme and palettes
- Separate priors for inverse logit (pg. 95 of Hobbs and Hooten 2015)
- 'Bad' inits, and control flow for default vs. JAGS-generated inits

## [0.3.0] - 2018-05-10
### Added
- Writing calling script and bare-bones metadata to outputs directory
- Extensive tabular output corresponding to the existing graphical summaries
- Trace plot functionality for 'fixed' (hard-coded) parameters (e.g., deviance)
- Lookup table in the README to outline the ways in which we are linking models
  to data
- Credible intervals not just for the mean, but for new observations
- Derived quantities for trend that are 'agnostic' WRT the deterministic
  function used in a given model
- More consistent labeling for axes in graphical output
- Better handling of sampling event date info

### Fixed
- Missing link in the beta models for ocular cover estimates

## [0.2.0] - 2018-04-20
### Added
- Support for modeling continuous proportion or cover data with a beta
  likelihood
- Support for specifying inits

### Removed
- The `write_mod_summaries_to_disk` function previously invoked in main calling
  script to save related model summary information to a single file

## [0.1.1] - 2018-04-18
### Fixed
- References/pointers to data were previously missing

### Added
- Instructions (in the README) for accessing data via Google Cloud Storage

### Security
- Using an ignored shell script to set object ACLs, giving access only to NPS
  colleagues

## 0.1.0 - 2018-04-18
### Added
- Cross OS support for the shell script used to compile JAGS model descriptions

[Unreleased]: https://gitlab.com/nps-ima/uplands-ci/compare/v0.7.0...HEAD
[0.7.0]: https://gitlab.com/nps-ima/uplands-ci/compare/v0.6.0...v0.7.0
[0.6.0]: https://gitlab.com/nps-ima/uplands-ci/compare/v0.5.0...v0.6.0
[0.5.0]: https://gitlab.com/nps-ima/uplands-ci/compare/v0.4.2...v0.5.0
[0.4.2]: https://gitlab.com/nps-ima/uplands-ci/compare/v0.4.1...v0.4.2
[0.4.1]: https://gitlab.com/nps-ima/uplands-ci/compare/v0.4.0...v0.4.1
[0.4.0]: https://gitlab.com/nps-ima/uplands-ci/compare/v0.3.1...v0.4.0
[0.3.1]: https://gitlab.com/nps-ima/uplands-ci/compare/v0.3.0...v0.3.1
[0.3.0]: https://gitlab.com/nps-ima/uplands-ci/compare/v0.2.0...v0.3.0
[0.2.0]: https://gitlab.com/nps-ima/uplands-ci/compare/v0.1.1...v0.2.0
[0.1.1]: https://gitlab.com/nps-ima/uplands-ci/compare/v0.1.0...v0.1.1

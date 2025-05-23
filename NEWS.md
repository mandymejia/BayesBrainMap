# 9.0

* Clean up dependencies
* Patch for `dual_reg_parc`

# 8.0

* Updates to prior calculations
    * Use effective sample size 
* Add prewhitening option to FC prior ICA
* Implement rescaling to BrainMap
* Default engagement type is `>` rather than `abs >`

# 7.1

* use `u` and `z` for activtions, rather than `gamma`

# 7.0

* no PARDISO requirement
* correct package documentation
* adjustments to FC, VB
* add scrubbing

# 6.0

* engagements: gray medial wall
* prep for adding FC prior

# 5.0

?

# 4.0

`estimate_prior`
* add custom smoothing
* add parallel over subjects

# 3.0

Fixes and cleaning up

# 2.0

Data pre-processing, for both prior estimation and TemplateICA
* `norm_BOLD` replaces `scale_BOLD`, since now it detrends the data too.
* Do not scale across space by default. Instead, scale across space only where necessary (one of the regressions in dual regression, and PCA-based operations).
* Option to detrend.
* Compute scale after centering and detrending.

Dual regression
* Option to normalize A matrix

Template estimation
* Merge `estimate_prior.cifti` and `estimate_prior.nifti` into `estimate_prior`
* Option to denoise each scan before computing DR
* If using pseudo retest data, split each subject's scan after denoising and detrending, not before.
* Option for different variance prior (non-negative)
* Option to obtain DR results
* Return parameters, including data normalization choices.
* `print` and `summary` for priors

TemplateICA
* Merge `BrainMap.cifti` and `BrainMap.nifti` into `BrainMap`
* Option to return DR results in `estimate_prior`
* Option to skip dimension reduction
* Option to use parallel computing for iterating over voxels.
* Option to provide `"xifti"` or `"nifti"` objects directly
* Option to provide multiple scans
* For denoising: replace `maxQ` with `Q2_max` to simplify logic
* Fix cases where `time_inds` and `resamp_res` are provided.
* Return parameters, including data normalization choices.
* Replace subject IC variance with the standard error 
* `summary` and `plot` for engagements

Misc
* Implement sign matching
* Refine `group_ICA`
* Refine `make_mesh`
* Remove INLA required version (bring back?)

For developers
* Laid some groundwork for estimating priors using >2 scans per subject: `var_decomp` can handle >2 scans per subject (but need to test it!)
* New code prior calculation

# 1.2

Update functions that work with `"xifti"` objects to conserve time and memory.

# 1.0

* Package maintenance:
    * format DESCRIPTION
    * change import to importFrom where possible
    * add README.Rmd
    * add travis and appveyor
    * delete Rhistory
    * clean up gitignore & buildignore
    * add internal keyword to functions not exported
    * add titles to functions
    * some formatting to function descriptions
    * replace equal signs with arrows
    * add CITATION

* Only require INLA for spatial modeling.
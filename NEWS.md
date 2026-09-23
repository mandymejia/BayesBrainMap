# 3.0

* Update scaling arguments. Use new `fMRItools` functions.
* Add scrubbing, nuisance regression, temporal filtering, and dropping volumes in a consistent manner between dual regression functions, `estimate_prior`, and `fit_BBM`. 
* Bug fixes.

# 2.0

### Rename:

The three main functions:

* (`estimate_prior` is the same)
* `BrainMap` -> `fit_BBM`
* `engagements` -> `id_engagements`

Also,

* `tvar_method` -> `var_method`

### Other changes

* Implement `FC_updateA`.
* Various improvements to the functions and documentation for calculating engagements and plotting
* Set default engagement level to `z=2` instead of `z=3`.
* FC BrainMap results structured differently in the object made by `fit_BBM`.

# 1.0 

Initialize package
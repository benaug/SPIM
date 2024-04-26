# SPIM
 Models fit:  2-Flank spatial partial identity model (SPIM), categorical SPIM, conventional and generalized categorical SMR. 
 
 No longer maintained, here for reproducibility reasons. Use nimble samplers also on github.

Note, the spatial mark-resight samplers in this package do not correctly handle "marked but not identified" and "unknown mark status" samples. Samplers on github are correct.

#4/26/24 Disclaimer: Do not use bernoulli observation model for catSPIM or SMR, the y/ID update isn't correct. (The bernoulli obsmod is correct in 2-flank SPIM). A correct version for catSPIM is available in nimble here: https://github.com/benaug/categorical-SPIM.
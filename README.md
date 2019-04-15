# crecs
## Comfy Multivariate Regression and Cross-validation for Spatial Objects

This R package is for multivariate regression and classification with spatial objects. It was developed for species-distribution modelling, but there's nothing biology-specific to it. Use it on spatially distributed continious quantities (as well as count and presence/absence data), featuring 
[random forest](http://cran.r-project.org/web/packages/randomForest/index.html),
[gbm.step](http://cran.r-project.org/web/packages/dismo/index.html), and
[gbm](http://cran.r-project.org/web/packages/gbm/index.html), 
[gam](http://cran.r-project.org/web/packages/mgcv/index.html), and
[maxent](http://www.cs.princeton.edu/~schapire/maxent/). It provides means for crossvalidation, various model performance metrics and helper functions for null models, spatial autocorrelation, etc., and a number of convenience functions.

----


### Installation
just issue
```
devtools::install_github("janhoo/crecs")
```
you may need to (`install.packages("devtools")`) first.

You are good to go.






### License

This package is licensed to you under the terms of the [GNU AFFERO GENERAL PUBLIC LICENSE](http://choosealicense.com/licenses/agpl-3.0/) version 3.0 or higher.

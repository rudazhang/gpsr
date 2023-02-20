# Gaussian process subspace regression, an R package [DEFUNCT: Development has moved to [https://github.com/UQUH/gpsr](https://github.com/UQUH/gpsr)]

This is a prototypical implementation of *GPS* in the R programming language.
For the original research article documenting the method, see [the Citation section](#citation).

## Installation

Option 1: Install the developing verion via `devtools`.

``` R
if (!("devtools" %in% installed.packages()[,"Package"])) {
    install.packages("devtools")
}
devtools::install_github("rudazhang/gpsr")
```

Option 2: Install from a bundled package.

First download a bundled package from [releases](https://github.com/rudazhang/gpsr/releases), then

``` R
## Run a line similar to the following.
install.packages("~/Downloads/gpsr_0.0.0.9000.tar.gz", type="source")
```

## Example Use

After installing the package, you can load it via:

``` R
library(gpsr)
```

Example scripts are included in `script/` under the installed directory. Find it via:

``` R
path <- system.file("script", package = "gpsr")
dir(path, full.names=TRUE)
```


## Citation

- Ruda Zhang, Simon Mak, and David Dunson.
  Gaussian Process Subspace Prediction for Model Reduction.
  SIAM Journal on Scientific Computing, 2022. https://epubs.siam.org/doi/10.1137/21M1432739

BibTeX citation:
``` bibtex
@Article{ZhangRD2022gps,
  author        = {Zhang, Ruda and Mak, Simon and Dunson, David},
  title         = {Gaussian Process Subspace Prediction for Model Reduction},
  journal       = {SIAM Journal on Scientific Computing},
  year          = {2022},
  volume        = {44},
  number        = {3},
  pages         = {A1428-A1449},
  doi           = {10.1137/21M1432739},
}
```

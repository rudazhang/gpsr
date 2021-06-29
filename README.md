# Gaussian process subspace regression, an R package

This is a prototypical implementation of *GP-subspace* in the R programming language.
For the original research article documenting the method, see [the Citation section](#citation).

## Installation

Option 1: Install the developing verion via `devtools`.

``` R
## Skip the first line if you have already installed `devtools`.
install.packages(devtools)
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

(arXiv paper to be added.)

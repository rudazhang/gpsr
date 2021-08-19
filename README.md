# Gaussian process subspace regression, an R package

This is a prototypical implementation of *GPS* in the R programming language.
For the original research article documenting the method, see [the Citation section](#citation).

## Installation

Option 1: Install the developing verion via `devtools`.

``` R
if (!("devtools" %in% installed.packages()[,"Package"])) {
    install.packages(devtools)
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
  Gaussian Process Subspace Regression for Model Reduction.
  arXiv, 2021.

BibTeX citation:
``` bibtex
@misc{ZhangRD2021gps,
  title={Gaussian Process Subspace Regression for Model Reduction},
  author={Ruda Zhang and Simon Mak and David Dunson},
  year={2021},
  eprint={2107.04668},
  archivePrefix={arXiv},
  primaryClass={math.ST},
  url={https://arxiv.org/abs/2107.04668},
}
```

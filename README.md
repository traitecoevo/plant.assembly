# {plant.assembly}: A package for modelling community assembly using the {plant} model trait ecology and evolution

<!-- badges: start -->
<!-- [![R-CMD-check](https://github.com/traitecoevo/plant.assembly/workflows/R-CMD-check/badge.svg)](https://github.com/traitecoevo/plant.assembly/master) -->
<!-- [![Codecov test coverage](https://codecov.io/gh/traitecoevo/plant.assembly/branch/master/graph/badge.svg)](https://codecov.io/gh/traitecoevo/plant?branch=master) -->
<!-- badges: end -->

The [{plant}](https://github.com/traitecoevo/plant) package for R is an extensible framework for modelling size-, patch- and trait-structured demography. The [{plant.assembly}](https://github.com/traitecoevo/plant.assembly) package extends the capabilities of {plant} to model community assembly using various algorithms. 

Current capabilities include

- stochastic assembly. 
- fitmax assembly

Envisioned future capabilities:

- More tools for 1D analysis, e.g. PIPs
- Solve 1-species attractor 
- use emulators to approximate fitness landscapes

## Installation

**Requirements**

- You must be using R 4.1.0 or newer. At this stage the package is not on CRAN. Your options for installing are described below.

**Option 1, using `remotes::install_github`**

The `{plant.assembly]` package can be installed direct from GitHub using the [`remotes`](https://cran.r-project.org/web/packages/remotes/index.html) package:

```r
remotes::install_github("traitecoevo/plant.assembly", dependencies=TRUE)
```

To install a specific (older) release, decide for the version number that you want to install in https://github.com/traitecoevo/plant.assembly/releases  e.g.

```r
remotes::install_github("traitecoevo/plant.assembly@v1.0.0", dependencies=TRUE)
```

with `"v1.0.0"` replaced by the appropriate version number.

**Option 2, building from source**

If familiar with [git](https://git-scm.com/) you might find it easiest to build `plant` directly from the source code. This is most useful if developing new models or strategies, or to contribute new features.

First, clone the `plant` repository

```
git clone https://github.com/traitecoevo/plant
```

Open an R session in the folder, then to install dependencies run

```
devtools::install_deps()
```

Then to compile the project

```
devtools::install()
```
or 

```
devtools::load_all()
```

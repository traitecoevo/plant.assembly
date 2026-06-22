# {plant.assembly}: Modelling community assembly with the {plant} model

<!-- badges: start -->
<!-- [![R-CMD-check](https://github.com/traitecoevo/plant.assembly/workflows/R-CMD-check/badge.svg)](https://github.com/traitecoevo/plant.assembly/master) -->
<!-- [![Codecov test coverage](https://codecov.io/gh/traitecoevo/plant.assembly/branch/master/graph/badge.svg)](https://codecov.io/gh/traitecoevo/plant?branch=master) -->
<!-- badges: end -->

The [{plant}](https://github.com/traitecoevo/plant) package for R is an extensible framework for modelling size-, patch- and trait-structured demography. The [{plant.assembly}](https://github.com/traitecoevo/plant.assembly) package extends {plant} to model community assembly and trait evolution, by repeatedly invading trait space, solving each community to demographic equilibrium, and following selection gradients to evolutionary attractors.

> **Status:** under active development. The fitness and equilibrium machinery was
> recently moved out of {plant} into this package, so a compatible recent {plant}
> is required (see `.plant-interface-version`). Developer notes, the working/not-working
> map, and build/test commands live in `CLAUDE.md`.

Current capabilities include:

- solving a community to demographic equilibrium (`community_demography`)
- 1D selection gradients and viable trait bounds
- stochastic assembly and maximum-fitness ("fitmax") assembly for one trait
- solving the 1-species, 1-trait evolutionary attractor

Envisioned future capabilities:

- emulators (Gaussian-process surrogates) to approximate fitness landscapes
- more tools for 1D analysis, e.g. pairwise invasibility plots (PIPs)
- attractor solving in 2D / 3D trait space


## Installation

**Requirements**

- You must be using R 4.5.0 or newer. At this stage the package is not on CRAN. Your options for installing are described below.

**Option 1, using `remotes::install_github`**

The `{plant.assembly}` package can be installed direct from GitHub using the [`remotes`](https://cran.r-project.org/web/packages/remotes/index.html) package:

```r
remotes::install_github("traitecoevo/plant.assembly", dependencies=TRUE)
```

To install a specific (older) release, decide for the version number that you want to install in https://github.com/traitecoevo/plant.assembly/releases  e.g.

```r
remotes::install_github("traitecoevo/plant.assembly@v1.0.0", dependencies=TRUE)
```

with `"v1.0.0"` replaced by the appropriate version number.

**Option 2, building from source**

If familiar with [git](https://git-scm.com/) you might find it easiest to build `plant.assembly` directly from the source code. This is most useful if developing new models or strategies, or to contribute new features.

First, clone the `plant.assembly` repository

```
git clone https://github.com/traitecoevo/plant.assembly
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

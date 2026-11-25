
<!-- README.md is generated from README.Rmd. Please edit that file and run devtools::build_readme() -->

# SRJPEmodel

Functions and documentation for the Spring Run Juvenile Production
Estimate (SRJPE) forecast model. Input data comes from the
[SRJPEdata](https://github.com/SRJPE/SRJPEdata) package, and the
forecast is run in [runSRJPE](https://github.com/SRJPE/runSRJPE).

The package currently contains the BT-SPAS-X juvenile abundance submodel
and the functions for storing model fits. The stock-recruit,
survival/travel-time, in-season and PLAD submodels are being developed
in their own repositories and will be added back as they are finalized.
Earlier versions of their code are in `archive/` and at the tag
`legacy-pre-migration`.

## Installation

``` r
# install.packages("pak")
pak::pak("SRJPE/SRJPEmodel")
```

## Submodels

- [BT-SPAS-X juvenile
  abundance](https://srjpe.github.io/SRJPEmodel/articles/btspas_submodel_documentation.html)
- Stock-recruit, survival, in-season and PLAD: articles are superseded
  until their code is added back.

## Development

``` r
devtools::document()     # regenerate NAMESPACE and man/
devtools::test()         # run tests
devtools::check()        # R CMD check
devtools::build_readme() # render this README
```

Style standards for this package are in
`planning/style_standards.md`.

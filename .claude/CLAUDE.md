# CLAUDE.md

Follow the FlowWest R standards (`FlowWest/infra`, `standards/r-coding-standards.md` and `standards/r-data-analysis-data-science-standards.md`) and this package's standards in `planning/style_standards_DRAFT.md`.

Keep functions and documentation simple. Match the code next to what you're changing, and change only what was asked.

Before proposing changes, run:

``` r
devtools::document()
devtools::test()
devtools::check()
```

`archive/` holds old code kept for reference during migration. It isn't part of the package; don't edit it.

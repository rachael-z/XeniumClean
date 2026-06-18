# Contributing to XeniumClean

## Setting up the package for development

You'll need `devtools`, `roxygen2`, and `testthat`:

```r
install.packages(c("devtools", "roxygen2", "testthat", "pkgdown"))
```

## Generating documentation

The `.R` files use roxygen2 comments. To regenerate the `man/` directory and
`NAMESPACE` from them:

```r
devtools::document()
```

Do this every time you change a function signature, add a function, or edit a
roxygen block. The version currently in this repo has the `NAMESPACE` written
by hand. After your first `devtools::document()` run, the `man/` directory
will be populated and `NAMESPACE` will be regenerated automatically.

## Building the package

```r
devtools::install()         # install locally for testing
devtools::check()           # run R CMD check (catches most issues)
devtools::test()            # run testthat suite
```

## Building the pkgdown site

```r
pkgdown::build_site()
```

The site lands in `docs/`. Push that subfolder to a `gh-pages` branch, or use
the `usethis::use_pkgdown_github_pages()` helper to set up GitHub Pages
automatically.

## Code style

- Function naming: `CamelCase` for exported user-facing functions (matches the
  Seurat convention readers will be used to). Internal helpers start with `.`.
- Argument naming: lowercase with dots (e.g. `group.by`, `gene.sets`,
  `coords.source`) — also matching Seurat.
- Roxygen blocks: title on first line, blank line, description, then
  `@param`, `@return`, `@examples`, `@export`. Mark internal helpers with
  `@keywords internal` and no `@export`.

## Adding a new feature

1. Branch off `main`: `git checkout -b feature/short-name`
2. Add or modify functions in `R/`
3. Add roxygen comments
4. Add at least one test in `tests/testthat/`
5. Update `NEWS.md`
6. Run `devtools::document()` and `devtools::check()`
7. Open a pull request

## Reporting bugs

Open an issue at https://github.com/rzemek/XeniumClean/issues with:

- Output of `sessionInfo()`
- The minimum reproducible example (small dataset preferred, dput-able if
  possible)
- The full error message and traceback

## Suggested improvements people may want to add

- Support for non-Seurat object types (SpatialExperiment, AnnData via
  zellkonverter)
- GPU-accelerated neighbour search for very large sections
- An interactive Shiny diagnostic dashboard for QC
- Built-in references for common tissues (download on first use)
- Per-cell removal rate threshold to flag potential outlier cells

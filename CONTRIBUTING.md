# Contributing to causaliflower

Open `causaliflower.Rproj` and use `devtools::load_all()` while editing.
This loads the package's functions together.


## Source style

- Keep the existing `.r` and `.R` filenames.
- Use the filename banner at the top of each source file. Describe the
  file's role in one or two short lines, without a final full stop.
- Leave two blank lines between sections and between the banner's filename
  and description. Put code or documentation directly below subsection headings.
- Keep constants near their main use, with comments on the assignment line.
  Do not wrap a short definition just to fit its comment. Keep related
  constants together when one depends on another.
- Put a comment about a whole function on its declaration line. Put comments
  about individual steps beside those steps.
- Leave the final full stop off single-sentence comments. Use normal
  punctuation for multiple sentences and roxygen documentation.


## Checking a change

Run the relevant tests while editing, then the full suite before release.

```r
devtools::load_all()
devtools::test(filter = "identification") # Use the relevant test file pattern
devtools::test()
```

Edit function help in the roxygen source. Use `@noRd` for internal functions
that do not need a help page. For 0.0.2, regenerate with roxygen2 7.3.3.
The RStudio project also runs roxygen when building, so check its version.

```r
stopifnot(as.character(packageVersion("roxygen2")) == "7.3.3")
roxygen2::roxygenise()
devtools::check(document = FALSE, args = "--as-cran")
```

Read the generated help and built vignette, including output and plots.
The HTML vignette needs Pandoc. Review the diff before submitting a change.

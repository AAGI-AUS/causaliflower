# Contributing to causaliflower

Bug reports, clearer documentation, examples and code changes are welcome.

To report a problem, [open an issue](https://github.com/AAGI-AUS/causaliflower/issues)
with a small reproducible example, what you expected, what happened, and your R
and package versions. For a larger change to the public interface, scientific
assumptions or dependencies, discuss the approach in an issue before starting.

## Working on code

See the [installation notes](README.md#installation) for package requirements.
From the package directory, load the code and run the tests relevant to your change:

```r
devtools::load_all()
devtools::test(filter = "identification")
```

Replace `identification` with the relevant test file pattern. RStudio users can
open `causaliflower.Rproj` to work in the package directory.

Keep changes focused and follow the surrounding style. Add regression coverage
for bug fixes and update documentation when behaviour changes. Edit function
help in its Roxygen source under `R/`. The [style guide](dev/STYLE.md) explains
documentation tools and the conventions used by maintainers and coding assistants.
Formatting details can be resolved during review.

## Sending a change

Review your diff, then open a pull request explaining the problem, the change
and how you checked it. Link any related issue. Mention checks you could not run
and any limitations you know about. Maintainers handle the broader package and
release checks described in the style guide.

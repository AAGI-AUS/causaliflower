# Style and maintenance guide

This guide records the conventions used by causaliflower's maintainers and
coding assistants. For the contribution workflow, see [CONTRIBUTING.md](../CONTRIBUTING.md).

Treat formatting and placement conventions as defaults for the code you change.
Follow nearby code when it differs. Use judgement when a convention makes a
change harder to understand, and explain substantial departures in review.
Keep formatting changes local to the work.

Changes should preserve existing behaviour and scientific meaning unless the
change is intentional and explained. Check the implementation before documenting
it. Coding assistants should read the relevant code and tests, keep edits within
the requested scope, and report checks actually run and any remaining uncertainty.

## Source organisation

Keep the existing `.r` and `.R` filenames. Within a file, place exported functions
before their dedicated internal helpers where the existing structure allows it.
Keep domain algorithms with their implementation and nested functions with the
code whose local state they use.

Shared helpers live in `R/helper_funcs.r`, under their main caller's filename.
Its sections use 80-hash banners and alphabetical filename order. Within a
section, follow caller order and put a helper's own helpers directly below it.
The public `string_to_dag()` function remains in the `helper_funcs.r` section,
followed by its parser helpers. Check callers, callbacks and exports before
moving a function. Moving it is a separate change from documenting it.

## Layout and names

- Follow the existing filename banner. Its role description is usually one or
  two short lines, without a final full stop.
- Use two blank lines between sections and between the banner's filename and
  description. Put code or documentation directly below subsection headings.
- Keep constants near their main use and dependent constants together. Put a
  useful comment beside the assignment, without wrapping a short definition
  solely to accommodate the comment.
- Use braces for `if` and `else` branches and retain explicit returns. Name local
  results for what they contain. A direct return needs no extra variable.
- Prefer descriptive domain names. Retain conventional method arguments,
  mathematical coordinates and generic utility names when they remain clear.

## Comments and prose

Add comments to explain reasons or non-obvious constraints. A comment about the
whole function belongs on its declaration line, and a comment about a step beside
that step. Remove comments that only repeat the code or Roxygen documentation.
Small nested callbacks can use local comments when an explanation is useful.

Omit the final full stop from single-sentence comments. Use normal punctuation
for multiple sentences and Roxygen documentation. Prefer plain verbs, concrete
nouns and British English. Avoid semicolon joins in prose. Preserve semicolons
that are part of executable syntax or data.

## Function documentation

Edit function help in its Roxygen source under `R/`. Document top-level internal
functions with their purpose, parameter meanings and return form. Use `@noRd`
when a help page is not needed. Reuse parameter documentation only when the
meanings match. Small nested callbacks rarely need a separate Roxygen block.

State logical returns directly as `TRUE` when the condition holds and `FALSE`
otherwise. Explain what happens when an option is unset. Describe empty results,
errors, mutation and ordering where callers depend on them. Read the function
and its callers before claiming that it guarantees a property.

When regenerating help, use the roxygen2 version recorded in DESCRIPTION:

```r
roxygen_version <- read.dcf("DESCRIPTION", fields = "Config/roxygen2/version")[[1L]]
stopifnot(as.character(packageVersion("roxygen2")) == roxygen_version)
roxygen2::roxygenise()
```

The RStudio project also runs roxygen when building, so check its version before
using that build action. Review changes to `man/`, NAMESPACE and DESCRIPTION.

## Checks for maintainers

Choose checks to match the change. Code changes need relevant tests and coverage
of the affected behaviour. Documentation changes need accurate examples, working
links and a reading of the affected help or vignette output. Changes to generated
help require inspection after regeneration. The HTML vignette needs Pandoc.

Run the full suite and package checks before release, and earlier when the scope
of a change calls for them:

```r
devtools::test()
devtools::check(document = FALSE, args = "--as-cran")
```

With `document = FALSE`, help must already be current. Read the check output,
including warnings and skips. For changed vignettes, inspect the built text,
results and plots. Record missing dependencies or checks that could not run.
Passing tests is one part of reviewing correctness and scientific claims.

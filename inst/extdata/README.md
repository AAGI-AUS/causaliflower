# Output availability register

`apsim-output-availability.csv` retains the seven column schema read by the
query engine. It contains no rows because the bundled data do not distribute
individual output locators. An empty register makes no claims of verified
output availability.

The columns are `crop`, `node_id`, `source_expression`, `status`,
`report_tested`, `evidence` and `note`. The reader selects the requested crop
and retains only records with status `source_verified_public` and a nonempty
source expression. Tests use a separate synthetic register to check those filters.

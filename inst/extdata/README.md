# Output availability register

`apsim-output-availability.csv` has the seven columns read by the query engine.
It contains column names but no records, so it does not certify that any
individual APSIM output is available.

The columns are `crop`, `node_id`, `source_expression`, `status`,
`report_tested`, `evidence` and `note`. The reader selects the requested crop
and retains only records with status `source_verified_public` and a nonempty
source expression. Tests use a separate synthetic register to check those filters.

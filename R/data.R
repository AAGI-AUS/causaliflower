# ==============================================================================
# data.R


# Roxygen documentation for the ESCDAGs lookup data
# ==============================================================================


#' Lookup table for ESCDAGs
#'
#' Causal criteria for assessing graph edges as presented in Ferguson et al., 2020, 'Evidence Synthesis for Constructing Directed Acyclic Graphs (ESC-DAGs): a novel and systematic method for building directed acyclic graphs'.
#'
#' About ESC-DAGs causal criteria ('Evidence synthesis for constructing directed acyclic graphs')
#' Each directed edge is assessed for three causal criteria: temporality, face validity and recourse to theory.
#' These are informed by the classic Bradford Hill viewpoints, compatible with the 'inference to the best explanation' approach advocated by Krieger and Davey Smith.
#' If a posited causal relationship satisfies all causal criteria, a counterfactual thought experiment derived from the potential outcomes framework can be used to further explicate and assess assumptions.
#' For more details, see Ferguson et al. (2020), DOI: https://doi.org/10.1093/ije/dyz150.
#' @format A data frame with 3 rows and five columns:
#' \describe{
#'   \item{name}{Causal criterion name or title.}
#'   \item{question}{Causal criterion phrased as a question.}
#'   \item{description}{Short (sentence/paragraph) explanation provided to users in case the question is not self-explanatory.}
#'   \item{source}{Journal article or reference.}
#'   \item{required}{Whether failing the criterion results in dropping an edge, e.g. if "yes" and a user inputs "n", the hypothesized causal relationship is rejected and the edge automatically removed. If "no", a user's answer is recorded, but it does not automatically remove the edge.}
#' }
#' @source Ferguson KD, McCann M, Katikireddi SV, Thomson H, Green MJ, Smith DJ,
#'   Lewsey JD (2020). 'Evidence synthesis for constructing directed acyclic
#'   graphs (ESC-DAGs): a novel and systematic method for building directed
#'   acyclic graphs.' International Journal of Epidemiology, 49(1), 322-329.
#'   DOI: https://doi.org/10.1093/ije/dyz150. The 'description' and 'source'
#'   columns quote p. 326 of this open-access article, distributed under the
#'   Creative Commons Attribution License 4.0
#'   (http://creativecommons.org/licenses/by/4.0/).
"ESCDAGs"


#' Readable structural causal model slice for DAG-by-query
#'
#' `apsim_dag_nodes` contains 576 quantities for barley, wheat and maize,
#' including shared quantities and the SoilWater family. Readable keys identify
#' quantities. Numbered concept keys connect structural references. Recovery
#' interpretations are included. Model path identifiers and individual output
#' locators are not distributed.
#'
#' @format A data frame with 576 rows and 22 character columns:
#' - `node_id`, `concept_id`, `node_name_readable`: quantity key, concept key and label.
#' - `tier`, `level`, `structural_category`: retained source classifications.
#' - `composition_of`, `function_of`, `init_of`: parent references, including crop
#'   qualifiers and recorded lags.
#' - `observable`, `measurable`, `proxy_measurable`: retained observation and measurement flags.
#' - `apsim_obtainable`: reserved output locator field, empty in this slice.
#' - `rs_obtainable`, `eo_obtainable`: retained source availability metadata.
#' - `crop_type`: applicable crops separated by semicolons, or `agnostic`.
#' - `modifiable_factor`, `provenance`: modification flag and provenance category.
#' - `hidden_node`, `path_dependent`, `temporal_class`: hidden detail and temporal metadata.
#' - `recovered_node`: `yes`, `yes:user` or `no:` followed by a reason.
#'
#' @details Boolean flags retain their source character values. The `cf_source_order`
#' attribute records ordering intervals and the quantity to concept mapping.
#' Queries retain this ordering metadata only while the identity mapping matches.
#' @docType data
"apsim_dag_nodes"


#' Aliases for the bundled query quantities
#'
#' `node_aliases` is a named list for the 576 quantities, containing 1,529 keys,
#' labels and plain language aliases. Conflicting aliases and model path aliases
#' are excluded.
#' @format A named list for the 576 quantities, containing 1,529 keys,
#'   labels and plain language aliases.
#' @docType data
"node_aliases"


#' Default management quantities for the bundled slice
#'
#' `management_defaults` contains three core quantity keys and an empty
#' `per_crop` list reserved for crop specific defaults.
#' @format Three core quantity keys and an empty `per_crop` list
#'   reserved for crop specific defaults.
#' @docType data
"management_defaults"


#' Developmental phase partitions for the bundled crops
#'
#' `crop_phases` describes developmental stages for barley, wheat and maize.
#' Its 32 rows contain `crop`, `phase_index` starting at zero,
#' `apsim_child_index`, `phase_name`, `start_stage`, `end_stage`, `phase_type`,
#' `terminal` and `canonical_stage`. The source type field is empty and model
#' class names are omitted.
#' @format A data frame with 32 rows and nine columns.
#' @docType data
"crop_phases"


#' Developmental phase membership of bundled quantities
#'
#' `phase_of` records active intervals for developmental stage queries. Its
#' 150 rows have seven character columns: `node_id`, `concept_id`, `crop_type`,
#' `active_from`, `active_to`, `archetype` and `basis`. These identify the quantity,
#' applicable crops, recorded interval, temporal category and supporting description.
#' @format A data frame with 150 rows and seven character columns.
#' @docType data
"phase_of"

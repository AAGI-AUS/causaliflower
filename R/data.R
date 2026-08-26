
#' Lookup table for ESCDAGs
#'
#' Causal criteria for assessing graph edges as presented in Ferguson et al., 2020, 'Evidence Synthesis for Constructing Directed Acyclic Graphs (ESC-DAGs): a novel and systematic method for building directed acyclic graphs'.
#'
#' About ESC-DAGs causal criteria ('Evidence synthesis for constructing directed acyclic graphs')
#' Each directed edge is assessed for three causal criteria: temporality; face-validity; and recourse to theory.
#' These are informed by the classic Bradford Hill viewpoints, compatible with the 'inference to the best explanation' approach advocated by Krieger and Davey Smith.
#' If a posited causal relationship satisfies all causal criteria, a counterfactual thought experiment derived from the potential outcomes framework can be used to further explicate and assess assumptions.
#' For more details, see Ferguson et al. (2020), DOI: https://doi.org/10.1093/ije/dyz150).
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

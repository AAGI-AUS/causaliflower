
# R/lookups.R

#' Lookup table for ESCDAGs
#'
#' Causal criteria for assessing graph edges as presented in Ferguson et al., 2020, 'Evidence Synthesis for Constructing Directed Acyclic Graphs (ESC-DAGs): a novel and systematic method for building directed acyclic graphs'.
#'
#' About ESC-DAGs causal criteria ('Evidence synthesis for constructing directed acyclic graphs')
#' Each directed edge is assessed for three causal criteria: temporality; face-validity; and recourse to theory.
#' These are informed by the classic Bradford Hill viewpoints, compatible with the ‘inference to the best explanation’ approach advocated by Krieger and Davey Smith.
#' If a posited causal relationship satisfies all causal criteria, a counterfactual thought experiment derived from the potential outcomes framework can be used to further explicate and assess assumptions.
#' For more details, see Ferguson et al. (2020), DOI: https://doi.org/10.1093/ije/dyz150).
#' @format A named list of 50 entries, each with three elements:
#' \describe{
#'   \item{name}{Causal criterion name or title.}
#'   \item{question}{Causal criterion phrased as a question.}
#'   \item{description}{Short (sentence/paragraph) explanation provided to users in case the question is not self-explanatory.}
#'   \item{source}{Journal article or reference.}
#'   \item{required}{Whether failing the criterion results in dropping an edge, e.g. if "yes" and a user inputs "n", the hypothesized causal relationship is rejected and the edge automatically removed. If "no", a user's answer is recorded, but it does not automatically remove the edge.}
#' }
"ESCDAGs"

#' Sentinel-2 spectral indices list
#'
#' A named list of spectral indices compatible with the
#' \href{https://github.com/awesome-spectral-indices/awesome-spectral-indices}{Awesome Spectral Indices} catalogue,
#' for use with \code{\link{get_rs_data}} and \code{\link{get_rs_raster}}.
#'
#' @format A named list where each element is a list with the following elements:
#' \describe{
#'   \item{assets}{Character vector of Sentinel-2 band names.}
#'   \item{fun}{A function that takes a SpatRaster and returns a single-layer SpatRaster.}
#'   \item{res}{Recommended resolution (highest integer value used by assets).}
#' }
"s2_index_list"


#' Landsat spectral indices list
#'
#' A named list of spectral indices compatible with Landsat 5 & 8 sources, for use with \code{\link{get_rs_data}} and \code{\link{get_rs_raster}}.
#'
#' @format A named list where each element is a list with the following elements:
#' \describe{
#'   \item{assets}{Character vector of Landast band names.}
#'   \item{fun}{A function that takes a SpatRaster and returns a single-layer SpatRaster.}
#'   \item{res}{Recommended resolution (highest integer value used by assets).}
#' }
"landsat_index_list"

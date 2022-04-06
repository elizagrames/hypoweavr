#' Synonyms grouped by term
#'
#' A dataset of terms to replace synonyms with, used in the example.
#'
#' @format A data frame 230 entries of 2 variables:
#' \describe{
#' \item{Term}{The replacement/standardized term}
#' \item{Synonyms}{Synonyms that should be standardized}
#' }
"ontology"



#' Study characteristics and hypotheses
#'
#' A dataset of the study characteristics/metadata and hypotheses
#' written as plain text that can be converted to a graphical model,
#' based on the case study systematic review.
#'
#' @format A data frame with 145 rows of 8 variables:
#' \describe{
#'   \item{Reference}{the citation for the study}
#'   \item{Article.type}{the type of publication, e.g. thesis or journal article}
#'   \item{Years.of.data.collection}{years during which data was collected}
#'   \item{Coordinates}{approximate study geographic coordinates}
#'   \item{Surrounding.landscape}{type of habitat surrounding the forest site}
#'   \item{Forest.type}{forest classification where study took place}
#'   \item{Focal.taxa}{four-letter alpha codes for bird species}
#'   \item{Pathways}{implied causal relationships written as plain text}
#' }
"studies"

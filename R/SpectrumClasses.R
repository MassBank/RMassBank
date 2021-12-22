#' @import MSnbase
#' @importFrom Biobase classVersion
#' @import S4Vectors
NULL

#' @title RMassBank Representation of an MSMS Spectrum
#'
#' @description This extends the \code{Spectrum2} class of the \code{MSnbase}
#' package and introduces further slots that are used to store information
#' during the \code{RMassBank} workflow.
#'
#' @slot satellite logical
#' If \code{TRUE}, the corresponding peak was removed as satellite.
#' @slot low logical
#' If \code{TRUE}, the corresponding peak was removed
#' because it failed the intensity cutoff.
#' @slot rawOk logical
#' If \code{TRUE}, the peak passed satellite and low-intensity cutoff removal.
#' @slot good logical
#' If \code{TRUE}, a formula could be found for the peak
#' and the peak passed all filter criteria. (see the
#' \code{RMassBank} vignette or the documentation of \code{\link{analyzeMsMs}}#' for details on filter settings)
#' @slot mzCalc numeric
#' The mz value calculated from the found formula for each peak (if any)
#' @slot formula character
#' The formula found for each peak.
#' \code{\link[rcdk]{generate.formula}} is used
#' for formula-fitting
#' @slot dbe numeric
#' The number of double bond equivalents.
#' This is calculated from the found formula for each peak (if any)
#' @slot formulaCount integer
#' The number of different formulae found for each peak.
#' Note: A peak for which multiple formulas were found will appear
#' multiple times. Hence there may be multiple entries in the \code{formula}
#' , \code{dppm} and \code{mzCalc} slot for the same mz value.
#' @slot formulaSource character "analyze" or "reanalysis"
#' Shows whether the current formula for the peak was determined by normal
#' analysis ("analyze") or by reanalysis of a failpeak ("reanalysis")
#' @slot dppm numeric
#' The ppm deviation of the mz value from the found formula (if any).
#' @slot dppmBest numeric
#' The ppm deviation of the mz value from the best formula found.
#' @slot ok logical one-element vector
#' If this is \code{TRUE}, the spectrum was successfully processed
#' with at least one resulting peak.
#' Otherwise, one of the following cases applies:
#' \itemize{
#' \item All peaks failed the intensity cutoff
#' i.e. the whole spectrum contains low intensity peaks, only.
#' \item All peaks were marked as satellites.
#' \item All peaks in the spectrum have a lower intensity than the value
#' given in the \code{specOkLimit} filter setting. (see the \code{RMassBank}
#' vignette or the documentation of \code{\link{analyzeMsMs}})
#' \item The precursor ion formula is invalid (see \code{\link{is.valid.formula}})
#' \item The spectrum is empty.
#' \item No molecular formula could be found for any of the peaks.
#' \item All peaks failed the \code{dbeMinLimit} criterion. (see the
#' \code{RMassBank} vignette or the documentation of \code{\link{analyzeMsMs}})
#' }
#' @slot info list
#' Spectrum identifying information
#' (collision energy, resolution, collision mode) from the \code{spectraList}
#' @slot properties data.frame
#' This is used as a flexible placeholder to store additional properties
#' for each peak throughout the workflow. After the last step of the
#' \code{mbWorkflow}, this will typically contain columns \code{mzRaw},
#' \code{noise}, \code{formulaMultiplicity}, \code{bestMultiplicity}
#' and \code{filterOK}. However, new columns may be added on demand
#' (see \code{\link{property<-}})
#' @seealso \code{\link[rcdk]{generate.formula}}, \code{\link{property<-}}
#' \code{\link{analyzeMsMs}}, \code{\link{generate.formula}},
#' \code{\link{is.valid.formula}}
#' @exportClass RmbSpectrum2
.RmbSpectrum2 <- setClass("RmbSpectrum2",
		representation = representation(
				satellite="logical",
				low="logical",
				rawOK ="logical",
				good = "logical",
				mzCalc = "numeric",
				formula = "character",
				dbe = "numeric",
				formulaCount = "integer",
				formulaSource = "character",
				dppm = "numeric",
				dppmBest = "numeric",
				ok = "logical",
				info = "list",
				properties = "data.frame"
		),
		contains=c("Spectrum2"),
		prototype = prototype(
				satellite = logical(),
				low = logical(),
				rawOK = logical(),
				good = logical(),
				mzCalc = numeric(),
				formula = character(),
				dbe = numeric(),
				formulaCount = integer(),
				formulaSource = character(),
				dppm = numeric(),
				dppmBest = numeric(),
				ok = logical(),
				info = list(),
				properties = data.frame(),
				new("Versioned", versions=c(classVersion("Spectrum2"), RmbSpectrum2 = "0.1.2"))
		),
)

#' @exportClass RmbSpectrum2List
.RmbSpectrum2List <- setClass("RmbSpectrum2List", contains="SimpleList",
		prototype=prototype(elementType="RmbSpectrum2"))
#
#setAs("ANY", "RmbSpectrum2List", function(from) {
#			coerceToSimpleList(from)
#		})

#' @exportClass RmbSpectraSet
.RmbSpectraSet <- setClass("RmbSpectraSet",
		representation = representation(
				parent = "Spectrum1",
				children = "RmbSpectrum2List",
				# These are done as slots and not as S4 functions, because they are set during the workflow
				# in "checking" steps. It's easier.
				found = "logical",
				complete = "logical",
				empty = "logical",
				formula = "character",
				id = "character",
				mz = "numeric",
				name = "character",
				mode = "character",
        smiles = "character"
				#annotations = "list"
				),
    contains=c("Versioned"),
		prototype = prototype(
				parent = new("Spectrum1"),
				children = new("RmbSpectrum2List"),
				found = FALSE,
				complete = NA,
				empty = NA,
				formula = character(),
				id = character(),
				mz = numeric(),
				name = character(),
				mode = character(),
        smiles = character(),
        new("Versioned", versions=c(RmbSpectraSet = "0.1.2"))
    # version 0.1.1: introduced versioning and SMILES slot
    # version 0.1.2: feed polarity to parent and children
		)
);

#' @exportClass RmbSpectraSetList
.RmbSpectraSetList <- setClass("RmbSpectraSetList", contains="SimpleList",
		prototype=prototype(elementType="RmbSpectraSet"))




setGeneric("getData",	function(s) standardGeneric("getData"))
setGeneric("setData",	function(s, df, ...) standardGeneric("setData"))


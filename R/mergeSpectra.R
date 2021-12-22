# TODO: Add comment
# 
# Author: stravsmi
###############################################################################
NULL


#' Merge peaks for spectra merging, FT shoulder elimination etc.
#' 
#' This procedure first sorts peaks by intensity (descending sort)
#' and then starts iterating over the peaks, removing all entries
#' that deviate "sufficiently far" from the currently selected peak.
#' See the Details section for a full explanation and information on
#' how to fine-tune peak removal.
#' 
#' Three parameters must be passed to \code{mergePeaks} for
#' peak-removal control in this order:
#' - cutoff_dppm_limit
#' - cutoff_absolute_limit
#' - cutoff_intensity_limit
#' The method iterates through the peaks, beginning with the
#' highest-intensity peak and in each step removes all other
#' peaks that fulfill conditions 1 AND 2 relative to the selected peak
#' 1. Their m/z value does not deviate too far from the one of the selected peak.
#' i.e. if the selected peak is p and the checked peak is c, it holds that
#' EITHER
#' |p$mz - c$mz| <= cutoff_absolute_limit
#' OR
#' |p$mz - c$mz| <= ppm(p$mz, cutoff_dppm_limit, p=TRUE)
#' (see \code{\link{ppm}})
#' 2. Their intensity is much smaller than the one of the selected peak, i.e.
#' c$mz < cutoff_intensity_limit * p$mz
#' for a suitable cutoff_intensity_limit between 0 and 1.
#'
#' @param peaks data.frame, matrix or RmbSpectrum2
#' The peak-table to be merged. In case of an \code{RmbSpectrum2}-object,
#' peaks are retrieved and updated via \code{\link{getData}}
#' and \code{\link{setData}}, respectively
#' @param ... 3 numeric values
#' These define cutoff limits (see details)
#' @return object of the same class as peaks
#' The result contains a reduced peak-table ordered by m/z
#' @examples \dontrun{mergePeaks(spectrum, 10, 0.5, 0.05)}
#' @seealso \code{\link{getData}}, \code{\link{setData}}, \code{\link{ppm}}
#' @rdname mergePeaks
#' @export
setGeneric("mergePeaks",	function(peaks, ...) standardGeneric("mergePeaks"))

#' Merge multiple spectra into one
#'
#' This method takes a collection of \code{RmbSpectrum2} objects
#' and merges them into a single \code{RmbSpectrum2} object
#'
#' Information from all spectra is retrieved via \code{\link{getData}}
#' combined with \code{rbind} and placed into the new spectrum with
#' \code{\link{setData}}
#'
#' @usage mergeSpectra(spectra, ...)
#' @param spectra \code{RmbSpectrum2List}
#' A list of \code{RmbSpectrum2} objects to be merged
#' @param ... NOTHING
#' (This parameter is reserved for future implementations of the generic)
#' @return A single \code{RmbSpectrum2} object
#' containing the merged information
#' @seealso \code{\link{getData}}, \code{\link{setData}}
#' @rdname mergeSpectra
#' @export
setGeneric("mergeSpectra", function(spectra, ...) standardGeneric("mergeSpectra"))

mergePeaks.df <- function(peaks, dppm, dabs, int)
{
	cutoff_int_limit <- int
	cutoff_mz_limit <- dabs
	cutoff_ppm_limit <- dppm
	# Order by intensity (descending)
	peaks_o <- peaks[order(peaks$intensity, decreasing=TRUE),,drop=FALSE]
	n <- 1
	# As long as there are peaks left AND the last peak is small enough (relative
	# to selected), move to the next peak
	while(n < nrow(peaks_o))
	{
		if(peaks_o[nrow(peaks_o),"intensity"] >= cutoff_int_limit *peaks_o[n,"intensity"])
			break
		# remove all peaks within cutoff_mz_limit (std. m/z = 0.5) which have intensity
		# of less than 5% relative to their "parent" peak
		#
		peaks_o <- peaks_o[ !(
							(
								((peaks_o$mz > peaks_o[n,"mz"] - cutoff_mz_limit)
									& (peaks_o$mz < peaks_o[n,"mz"] + cutoff_mz_limit))
								|  ((peaks_o$mz > peaks_o[n,"mz"] - ppm(peaks_o[n, "mz"], cutoff_ppm_limit, p=TRUE))
									& (peaks_o$mz < peaks_o[n,"mz"] + ppm(peaks_o[n, "mz"], cutoff_ppm_limit, p=TRUE)))
								)
							& (peaks_o$intensity < cutoff_int_limit * peaks_o[n,"intensity"])
							
							),,drop=FALSE]		 
		n <- n+1
	}
	return(peaks_o[order(peaks_o$mz),,drop=FALSE])
}

#' @rdname mergePeaks
setMethod("mergePeaks", "data.frame", function(peaks, ...)
		{
			mergePeaks.df(peaks, ...)
		})

#' @rdname mergePeaks
setMethod("mergePeaks", "matrix", function(peaks, ...)
		{
			mergePeaks.df(peaks, ...)
		})

#' @rdname mergePeaks
setMethod("mergePeaks", "RmbSpectrum2", function(peaks, ...)
		{
			df <- getData(peaks)
			df <- mergePeaks.df(df, ...)
			peaks <- setData(peaks, df)
			return(peaks)
		})

#' @rdname mergePeaks
setMethod("mergePeaks", "Spectrum", function(peaks, ...)
		{
			df <- as.data.frame(peaks)
			df <- mergePeaks.df(df, ...)
			peaks@mz <- df[,1]
			peaks@intensity <- df[,2]
			peaks@peaksCount <- nrow(df)
			return(peaks)
		})


mergeSpectra.RmbSpectrum2List <- function(spectra)
{
	if(length(spectra) == 0)
		stop("Merging empty spectra lists doesn't work because of missing metadata!")
	df.all <- lapply(spectra, getData)
	df <- do.call(rbind, df.all)
	spectrum <- spectra[[1]]
	spectrum <- setData(spectrum, df)
	spectrum@collisionEnergy <- -1
	spectrum@merged <- length(spectra)
	return(spectrum)
}

#' @rdname mergeSpectra
setMethod("mergeSpectra", "RmbSpectrum2List", function(spectra, ...)
			mergeSpectra.RmbSpectrum2List(spectra, ...))





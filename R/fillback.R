# TODO: Add comment
# 
# Author: stravsmi
###############################################################################

#' Fill back reanalyzed / refiltered peak info into spectra
#' 
#' This method takes the info which is added to the aggregated table in the reanalysis and 
#' multiplicity filtering steps of the workflow, and adds it back into the spectra.
#' 
#' @param o msmsWorkspace, RmbSpectraSet or RmbSpectrum2
#' The object information is filled back into. If applied to an RmbSpectraSet, information is added to all its RmbSpectrum2 children. If applied to the whole msmsWorkspace, information is added to all SpectraSets.
#' @param id character or missing
#' The id of the parent RmbSpectraSet if applied to RmbSpectrum2
#' @param aggregated data.frame or missing
#' The aggregated table of the parent msmsWorkspace if applied to RmbSpectraSet or RmbSpectrum2
#' @return o msmsWorkspace, RmbSpectraSet or Rmbspectrum2
#' The same object that was given as input with new information filled into it
#' @rdname fillback
#' @export
setGeneric("fillback", function(o, id, aggregated) standardGeneric("fillback"))

#' @rdname fillback
setMethod("fillback", c("msmsWorkspace", "missing", "missing"), function(o)
    {
      for(i in seq_len(length(o@spectra)))
        o@spectra[[i]] <- fillback(o@spectra[[i]], o@aggregated)
      o
    })

#' @rdname fillback
setMethod("fillback", c("RmbSpectraSet", "missing", "data.frame"), function(o, aggregated)
    {
      for(i in seq_len(length(o@children)))
        o@children[[i]] <- fillback(o@children[[i]], o@id, aggregated)
      o
    })

#' @rdname fillback
setMethod("fillback", c("RmbSpectrum2", "character", "data.frame"), function(o, id, aggregated)
    {
      .fillback(o, id, aggregated)  
    })

.fillback <- function(o, id, aggregated)
{
  peaks <- selectPeaks(aggregated,
      (cpdID == id) & (scan == o@acquisitionNum))
  curPeaks <- getData(o)
  # Check that our data processing assumptions are correct: the peaks that are rawOK are the peaks
  # that are in the aggregate table. (This is a very rough check: nrow identical)
  stopifnot(nrow(peaks)==nrow(curPeaks[curPeaks$rawOK,,drop=FALSE]))
  # If dppmBest is in the aggregate table, drop the one in the current spectrum, because the refiltering
  # can remove peaks with better dppm!
  if(("dppmBest" %in% colnames(peaks)) & ("dppmBest" %in% colnames(curPeaks)))
    curPeaks <- curPeaks[,which(colnames(curPeaks) != "dppmBest"),drop=FALSE]
  # Find all columns that we don't have yet in curPeaks
  colnames(peaks)[colnames(peaks) == "mzFound"] <- "mz"
  colsRemove <- c("scan", "cpdID", "parentScan", "index")
  colsRemoveIndex <- which(colnames(peaks) %in% colsRemove)
  peaks <- peaks[,-colsRemoveIndex,drop=FALSE]
  colsNew <- setdiff(colnames(peaks), colnames(curPeaks))
  for(col in colsNew)
  {
    if(col != "dppmBest")
      o <- addProperty(o, col, class(peaks[,col]))
  }
  peaksNew <- merge(curPeaks, peaks, all=TRUE)
  # Check that no stray "new peaks" were added by incorrect merging. If this happens, we have to write cleaner code
  stopifnot(nrow(peaksNew) == nrow(curPeaks))
  o <- setData(o, peaksNew)
  #browser()
  return(o)
}

# TODO: Add comment
# 
# Author: stravsmi
###############################################################################
#' @importFrom Biobase isVersioned
#' @importFrom Biobase isCurrent

.updateObject.RmbSpectrum2 <- setMethod("updateObject", signature(object="RmbSpectrum2"), function(object, ..., verbose = FALSE) 
		{
      callNextMethod()
			w <- object
			if(isVersioned(w))
				if(all(isCurrent(w)))
					return(w)
			# get RmbSpectrum2 version
			if(!isVersioned(w))
				v <- "0.0.0"
			else
				v <- classVersion(w)["RmbSpectrum2"]
			if(v < "0.1.1")
			{
				slot(w, "properties", check=FALSE) <- data.frame()
				classVersion(w)["RmbSpectrum2"] <- "0.1.1"
			}
			if(v < "0.1.2")
			{
				w <- .updateObject.RmbSpectrum2.formulaSource(w)
				classVersion(w)["RmbSpectrum2"] <- "0.1.2"
			}
			return(w)
		})



#' Add formulaSource column to spectrum.
#' 
#' TODO: consider whether to add functionality to move reanalysis stuff from legacy data back in.
#' 
#' @param w RmbSpectrum2
#' The object to be updated
#' @return The updated RmbSpectrum2 object
#' 
#' @author stravsmi
#' @export
.updateObject.RmbSpectrum2.formulaSource <- function(w)
{
	slot(w, "formulaSource", check=FALSE) <- character()
	if(length(w@formula) > 0)
	{
		w@formulaSource <- rep(NA, length(w@formula))
		w@formulaSource[!is.na(w@formula)] <- "analyze" 
	}
	w
}

.updateObject.RmbSpectraSet <- setMethod("updateObject", signature(object="RmbSpectraSet"), function(object, ..., verbose = FALSE) 
    {
      # first update the RmbSpectrum2s via the automated method:
      object <- updateObjectFromSlots(object, ..., verbose = verbose)
      
      w <- object
      if(isVersioned(w))
        if(all(isCurrent(w)))
          return(w)
      # get RmbSpectraSet version
      if(!isVersioned(w))
        v <- "0.0.0"
      else
        v <- classVersion(w)["RmbSpectraSet"]
      if(v < "0.1.1")
      {
        slot(w, "smiles", check=FALSE) <- character()
      }
      if(v < "0.1.2")
      {
        w <- .updateObject.RmbSpectraSet.updatePolarity(w)
        classVersion(w)["RmbSpectraSet"] <- "0.1.2"
      }
      
      return(w)
    })


.updateObject.RmbSpectraSet.updatePolarity <- function(w)
{
  w@parent@polarity <- .polarity[[w@mode]]
  for(n in seq_len(length(w@children)))
  {
    w@children[[n]]@polarity <- .polarity[[w@mode]]
  }
  w
}

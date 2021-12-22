# TODO: Add comment
# 
# Author: stravsmi
###############################################################################


#' Select peaks from aggregate table
#' 
#' Selects peaks from aggregate table according to different criteria.
#' 
#' @param o \code{msmsWorkspace} or aggregate \code{data.frame} from a workspace.
#' @param good if \code{TRUE}, include good (matched within filter criteria) peaks.
#' @param bad if \code{TRUE}, include bad (not matched within filter criteria) peaks. Note: \code{good} and \code{bad} can be
#' 			combined, both are returned in that case.
#' @param cleaned if \code{TRUE}, return only peaks which passed the noise filter. Note: If the noise filter was not applied, the
#' 			parameter has no effect. Also, a \code{noise} column is in any case added to the output, even if not present before.
#' @param best if \code{TRUE}, only select the best match for each peak (i.e. the formula with smallest delta ppm). Otherwise multiple
#' 			matches can be returned.
#' @param ... no additional parameters
#' @return Peak dataframe according to the specified criteria.
#' 
#' @author stravsmi
#' @export
setGeneric("selectPeaks", function(o, ...) standardGeneric("selectPeaks"))

#' Add and initialize dataframe column
#' 
#' Adds a new column of a defined type to a \code{data.frame} and initializes it to a value.
#' The advantage of doing this over adding it with \code{$} or \code{[,""]} is that the case
#' \code{nrow(o) == 0} is adequately handled and doesn't raise an error. 
#' 
#' @param o \code{data.frame} to add the column to
#' @param name Name of the new column
#' @param type Data type of the new column
#' @param value Initial value of the new column (\code{NA} if not given)
#' @return Expanded data frame.
#' 
#' @author stravsmi
#' @export
setGeneric("addProperty", function(o, name, type, value=NA) standardGeneric("addProperty"))

#' Get a property of an RmbSpectrum2 object
#'
#' This searches the 'properties' slot of the object
#' and returns a column with matching name (if found)
#' or NULL otherwise.
#'
#' @param o \code{RmbSpectrum2}
#' @param property character
#' The name of a property
#' @return The corresponding column of \code{o@properties}
#' @rdname property
#' @export
setGeneric("property", function(o, property) standardGeneric("property"))

#' Replacement function to set properties of an RmbSpectrum2 object
#'
#' Update the 'properties' slot of the given object.
#' If the column you want to update does not exist yet and
#' \code{addNew = FALSE} (default), this will cause a warning
#' and the object will not be changed
#'
#' Please note that this is a replacement method, meaning that
#' \code{property(o, property) <- value}
#' can be used as a short-hand for the equivalent
#' \code{o <- 'property<-'(o, property, value)}
#'
#' @usage property(o, property, addNew=FALSE, class="") <- value
#' @param o \code{RmbSpectrum2}
#' The object whos 'properties' slot should be updated
#' @param property character
#' The name of the column in the 'properties' data frame to be updated
#' @param addNew logical, Default: FALSE
#' Whether or not a new column should be added in case a column of the
#' given name does not exist yet.
#' @param class character or missing
#' The class of the entries for the column to be added/updated
#' @param value ANY
#' The value(s) to be written into the column
#' @return The \code{RmbSpectrum2} object with an updated 'properties' slot
#' @rdname property-set
#' @export
setGeneric("property<-", function(o, property, addNew = FALSE, class="", value) standardGeneric("property<-"))

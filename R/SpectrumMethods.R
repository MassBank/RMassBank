# TODO: Add comment
# 
# Author: stravsmi
###############################################################################

#' Get data frame with all present peak data 
#' 
#' Returns a data frame with columns for all non-empty slots in a \code{RmbSpectrum2} object. Note that \code{MSnbase::Spectrum} has
#' a method \code{as.data.frame}, however that one will return only mz, intensity. This function is kept separate to ensure downwards
#' compatibility since it returns more columns than MSnbase \code{as.data.frame}.
#' @name getData
#' @aliases getData,RmbSpectrum2-method
#' 
#' @param s The \code{RmbSpectrum2} object to extract data from.
#' @return A data frame with columns for every set slot.
#' 
#' @author stravsmi
#' @docType methods
#' @export
setMethod("getData", c("RmbSpectrum2"), function(s)
		{
			peaks <- s@peaksCount
			cols <- c("mz", "intensity", "satellite", "low", "rawOK", "good", "mzCalc", "formula", "dbe", "formulaCount", "formulaSource", "dppm", "dppmBest")
			slotLength <- unlist(lapply(cols, function(col) length(slot(s, col))))
			cols.isFilled <- slotLength == peaks
			cols.filled <- cols[cols.isFilled]
			data <- lapply(cols.filled, function(col) slot(s, col))
			data$stringsAsFactors <- FALSE
			df <- do.call(data.frame,data)
			colnames(df) <- cols.filled		
			if(nrow(s@properties) == peaks)
				df <- cbind(df, s@properties)
			else if(nrow(s@properties) > 0)
				stop("Incorrect number of rows in properties frame.")
			return(df)
		})




#' Set \code{RmbSpectrum2} data from data.frame
#' 
#' Sets all slots which are present as columns in the given dataframe. Optionally cleans the object, i.e. empties slots not defined in the data frame.
#' 
#' @name setData
#' @aliases setData,RmbSpectrum2,data.frame-method
#' 
#' @param s The \code{RmbSpectrum2} object to modify
#' @param df The data frame with new data
#' @param clean \code{TRUE} if slots which aren't present as columns in the data frame should be cleared.
#' @return The modified \code{RmbSpectrum2}.
#' 
#' @author stravsmi
#' @docType methods
#' @export
setMethod("setData", c("RmbSpectrum2", "data.frame"), function(s, df, clean = TRUE)
		{
			s <- .setData.main(s, df, clean)
			s <- .setData.properties(s, df, clean)
			s
		})

.setData.main <- function(s, df, clean = TRUE)
{
	cols <- c("mz", "intensity", "satellite", "low", "rawOK", "good", "mzCalc", "formula", "dbe", "formulaCount", "formulaSource", "dppm", "dppmBest")
	types <- c("mz" = "numeric", "intensity" = "numeric", "satellite" = "logical", "low" = "logical",
			"rawOK" = "logical", "good" = "logical", "mzCalc" = "numeric", "formula" = "character", 
			"dbe" = "numeric", "formulaCount" = "integer", "formulaSource" = "character", "dppm" = "numeric", "dppmBest" = "numeric"
	)
	s@peaksCount <- as.integer(nrow(df))
	cols.inDf <- (cols %in% colnames(df))
	cols.df <- cols[cols.inDf]
	for(col in cols.df)
	{
		slot(s, col) <- as(df[,col], types[[col]])
	}
	cols.notinDf <- !(cols.inDf)
	cols.no <- cols[cols.notinDf]
	if(clean)
	{
		for(col in cols.no)
		{
			slot(s, col) <- new(class(slot(s, col)))
		}
	}
	s
}

.setData.properties <- function(s, df, clean = TRUE)
{
	# first set everything that RmbSpectrum2 setData does already
	# then find which columns can be set for properties and always clean (do no remove properties!)
	
  # Find all properties which have a column in the new data (df)
	cols <- colnames(s@properties)
	cols.inDf <- (cols %in% colnames(df))
	
	#newDf <- s@properties[rep(NA, s@peaksCount),,drop=FALSE]
  # cols.df contains all columns in the new df that are a property
	cols.df <- cols[cols.inDf]
	newDf <- df[,cols.df,drop=FALSE]
	
  if(!clean)
  {
    # Properties which do not have a column in the data frame:
    # (Note: we do not care about columns in the dataframe that don't correspond to a property. These just go lost)
    # (Maybe we should warn about this?)
    cols.notinDf <- !(cols.inDf)
    cols.no <- cols[cols.notinDf]
    
    for(col in cols.no)
      newDf[,col] <- rep(as(NA, class(newDf[,col])), nrow(newDf))
    # reorder columns
    newDf <- newDf[,cols,drop=FALSE]
  }
  else
  {
    newDf <- newDf[,cols[cols.inDf],drop=FALSE]
  }
  
	s@properties <- newDf
	s
}

setMethod("initialize", "RmbSpectrum2", function(.Object, ...,
				satellite = .Object@satellite,
				low = .Object@low,
				rawOK = .Object@rawOK,
				good = .Object@good,
				mzCalc = .Object@mzCalc,
				formula = .Object@formula,
				dbe = .Object@dbe,
				formulaCount = .Object@formulaCount,
				formulaSource = .Object@formulaSource,
				dppm = .Object@dppm,
				dppmBest = .Object@dppmBest,
				ok = .Object@ok,
				info = .Object@info,
				properties = .Object@properties
				) {
			## do work of initialization
			callNextMethod(.Object, ..., 
					satellite = satellite,
					low = low,
					rawOK = rawOK,
					good = good,
					mzCalc = mzCalc,
					formula = formula,
					dbe = dbe,
					formulaCount = formulaCount,
					formulaSource = formulaSource,
					dppm = dppm,
					dppmBest = dppmBest,
					ok = ok,
					info = info,
					properties = properties)
		})



#' @export
#' @describeIn selectPeaks A method to filter spectra to the specified peaks
setMethod("selectPeaks", c("RmbSpectrum2"), function(o, filter, ..., enclos=parent.frame(2))
		{
			if(missing(filter))
				return(o)
			df <- getData(o)
			f <- substitute(filter)
			df <- df[eval(f, df, enclos) & !is.na(eval(f, df, enclos)),,drop=FALSE]
			o <- setData(o, df)
			o
		})

#' @export
setMethod("selectPeaks", c("Spectrum"), function(o, filter, ..., enclos=parent.frame(2))
		{
			if(missing(filter))
				return(o)
			
			df <- as.data.frame(o)
			f <- substitute(filter)
			df <- df[eval(f, df, enclos),,drop=FALSE]
			o@mz <- df[,1]
			o@intensity <- df[,2]
			o
		})


#' @export
#' @describeIn selectPeaks A method to filter spectra to the specified peaks
setMethod("selectPeaks", c("RmbSpectrum2List"), function(o, ..., enclos=parent.frame(2))
		{
			for(n in seq_len(length(o)))
				o[[n]] <- selectPeaks(o[[n]], ..., enclos=enclos)
			return(o)
		})

#' @export 
setMethod("normalize", c(object="RmbSpectrum2"), function(object, ..., scale=999, precision=0, slot="intensity")
		{
			intensity <- .normalize.int(object, ...) * scale
			if(precision >= 0) intensity <- round(intensity, precision)
			if(slot != "intensity")
			{
				property(object, slot, TRUE, "numeric") <- intensity
			}
			else
				object@intensity <- intensity
			object
		})


.normalize.int <- function(object, ...)
{
	spec <- selectPeaks(object, ...)
	maxint <- max(spec@intensity)
	intensity <- object@intensity / maxint
	return(intensity)
}

#' @export
setMethod("normalize", c("RmbSpectrum2List"), function(object, ...)
		{
			s <- lapply(object, function(s) normalize(s, ...))
			for(n in seq_len(length(object)))
				object[[n]] <- s[[n]]
			return(object)
		})



setMethod("+", c("Spectrum", "numeric"), function(e1, e2) 
		{
			e1@mz <- e1@mz + e2
			return(e1)
		}) 



setMethod("-", c("Spectrum", "numeric"), function(e1, e2) 
		{
			e1@mz <- e1@mz - e2
			return(e1)
		}) 

setMethod("+", c("RmbSpectraSet", "ANY"), function(e1, e2)
		{
			e1@parent <- e1@parent + e2
			for(n in seq_len(length(e1@children)))
				e1@children[[n]] <- e1@children[[n]] + e2
			e1
		})

setMethod("-", c("RmbSpectraSet", "ANY"), function(e1, e2)
		{
			e1@parent <- e1@parent - e2
			for(n in seq_len(length(e1@children)))
				e1@children[[n]] <- e1@children[[n]] - e2
			e1
		})


setMethod("+", c("RmbSpectrum2List", "ANY"), function(e1, e2)
		{
			for(n in seq_len(length(e1)))
				e1[[n]] <- e1[[n]] + e2
			e1
		})


setMethod("-", c("RmbSpectrum2List", "ANY"), function(e1, e2)
		{
			for(n in seq_len(length(e1)))
				e1[[n]] <- e1[[n]] - e2
			e1
		})

#' @describeIn addProperty Add a new column to the RmbSpectrum2 properties
#'
#' @export 
setMethod("addProperty", c("RmbSpectrum2", "character", "character", "ANY"), function(o, name, type, value=NA)
		{
			if(ncol(o@properties) == 0)
				o@properties <- data.frame(row.names = seq_len(o@peaksCount)) 
			o@properties[,name] <- as(rep(value, o@peaksCount), type)
			o
		})

#setGeneric("setData",	function(s, df, ...) standardGeneric("setData"))


#' @rdname property
setMethod("property", c("RmbSpectrum2", "character"), function(o, property)
		{
			if(property %in% colnames(o@properties))
				return(o@properties[,property])
			else
				# We can't use FALSE or NA, since it could be confused with a 1-length logical FALSE or 1-length ANY NA  
				return(NULL)
		})


.propertySet <- function(o, property, addNew = FALSE, class="", value)
{
	if(class == "") class <- class(value)
	if(!(property %in% colnames(o@properties)) & !addNew)
	{
		warning("Trying to set inexistent property.")
		return(o)
	}
	else if(!(property %in% colnames(o@properties)))
		o <- addProperty(o, property, class)
	o@properties[,property] <- value
	return(o)
}


#' @rdname property-set
setMethod("property<-", signature(o="RmbSpectrum2", property="character", addNew="logical", class="character", value="ANY"), .propertySet )

#' @rdname property-set
setMethod("property<-", signature(o="RmbSpectrum2", property="character", addNew="missing", class="character", value="ANY"), .propertySet )

#' @rdname property-set
setMethod("property<-", signature(o="RmbSpectrum2", property="character", addNew="logical", class="missing", value="ANY"), .propertySet)

#' @rdname property-set
setMethod("property<-", signature(o="RmbSpectrum2", property="character", addNew="missing", class="missing", value="ANY"), .propertySet )


.fillSlots <- function(o, slotNames)
{
  for(entry in slotNames)
  {
    if(length(slot(o, entry)) != length(o@mz))
      slot(o, entry) <- rep(new(class(slot(o, entry)),NA), length(o@mz))
  }
  return(o)
}


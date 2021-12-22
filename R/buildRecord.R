# TODO: Add comment
# 
# Author: stravsmi
###############################################################################

#' @import assertthat


#' @title Build MassBank records
#' 
#' @description Takes a spectra block for a compound, as returned from
#' \code{\link{analyzeMsMs}}, and an aggregated cleaned peak table, together
#' with a MassBank information block, as stored in the infolists and loaded via
#' \code{\link{loadInfolist}}/\code{\link{readMbdata}} and processes them to a
#' MassBank record
#' 
#' @usage buildRecord(o, ..., cpd, mbdata, analyticalInfo, additionalPeaks)
#' @param o \code{RmbSpectraSet} or \code{RmbSpectrum2}
#' The spectra (or single spectrum) should be taken from a compound after analysis (\code{\link{analyzeMsMs}}).
#' Note that \bold{peaks are not read from this
#' object anymore}: Peaks come from the \code{aggregated} dataframe (and from
#' the global \code{additionalPeaks} dataframe; cf. \code{\link{addPeaks}} for
#' usage information.)
#' @param ...
#' keyword arguments for intensity normalization and peak selection (see \code{\link{normalize}} and \code{\link{selectPeaks}})
#' @param cpd \code{RmbSpectraSet} or missing
#' In case o is an \code{RmbSpectrum2}, this represents the \code{RmbSpectraSet} it belongs to
#' @param mbdata list
#' The information data block for the record header, as stored in
#' \code{mbdata_relisted} after loading an infolist.
#' @param additionalPeaks data.frame
#' If present, a table with additional peaks to add into the spectra.
#' 		As loaded with \code{\link{addPeaks}}.
#' @return An object of the same type as was used for the input with new information added to it
#' @author Michael Stravs
#' @seealso \code{\link{mbWorkflow}}, \code{\link{addPeaks}},
#' \code{\link{gatherCompound}}, \code{\link{toMassbank}}
#' @references MassBank record format:
#' \url{http://www.massbank.jp/manuals/MassBankRecord_en.pdf}
#' @rdname buildRecord
#' @export
setGeneric("buildRecord", function(o, ..., cpd, mbdata, analyticalInfo, additionalPeaks) standardGeneric("buildRecord"))

.buildRecord.RmbSpectraSet <- function(cpd, ..., mbdata = list(), additionalPeaks = NULL)
{
	# gather the individual spectra data
	analyticalInfo <- getAnalyticalInfo(cpd)

	# Go through all child spectra, and add metadata to all the info slots
	# Pass them the AC_LC and AC_MS data, which are added at the right place
	# directly in there.
	allSpectra <- lapply(cpd@children, function(s)
				buildRecord(s, ..., cpd=cpd, mbdata=mbdata, analyticalInfo=analyticalInfo, 
            additionalPeaks=additionalPeaks))
	allSpectra <- allSpectra[which(!is.na(allSpectra))]
	if(length(allSpectra) > 0)
		cpd@children <- as(allSpectra, "SimpleList")
	else
		cpd@children <- new("RmbSpectrum2List")
  cpd
}

#' @rdname buildRecord
setMethod("buildRecord", "RmbSpectraSet", function(o, ..., mbdata = list(), additionalPeaks = NULL)
      .buildRecord.RmbSpectraSet(cpd=o, ..., mbdata = mbdata, additionalPeaks = additionalPeaks)
    )


.addGenericInfo <- function(ac, annotations, search_string=c("^AC\\$MASS_SPECTROMETRY_", "^AC\\$CHROMATOGRAPHY_")) {
	# Note: For whatever reason, recursivity is inverted for the unlist
	# function, meaning that recursive=FALSE actually leads to the
	# behaviour expected when setting recursive=TRUE, which is desired
	# here, because nested lists exist. See help(unlist)

	properties <- names(unlist(annotations, recursive=FALSE))
	presentProperties <- names(ac)
	
	theseProperties <- grepl(x = properties, pattern = search_string)
	properties2     <- gsub(x = properties, pattern = search_string,
	  replacement = "")
	theseProperties <- theseProperties &
	  !(properties2 %in% presentProperties)
	theseProperties <- theseProperties &
	  (unlist(annotations, recursive=FALSE) != "NA")
	ac[properties2[theseProperties]] <-
	  unlist(annotations, recursive=FALSE)[theseProperties]
	return(ac)
}

# For each compound, this function creates the "lower part" of the MassBank record, i.e.
# everything that comes after AC$INSTRUMENT_TYPE.

#' Compose data block of MassBank record
#' 
#' \code{gatherCompound} composes the data blocks (the "lower half") of all
#' MassBank records for a compound, using the annotation data in the RMassBank
#' options, spectrum info data from the \code{analyzedSpec}-type record and the
#' peaks from the reanalyzed, multiplicity-filtered peak table. It calls
#' \code{gatherSpectrum} for each child spectrum.
#' 
#' The returned data blocks are in format \code{list( "AC\$MASS_SPECTROMETRY" =
#' list('FRAGMENTATION_MODE' = 'CID', ...), ...)} etc.
#' 
#' @aliases gatherCompound gatherSpectrum
#' @usage gatherCompound(spec, aggregated, additionalPeaks = NULL, retrieval="standard")
#' 
#' 		gatherSpectrum(spec, msmsdata, ac_ms, ac_lc, aggregated, 
#'	 		additionalPeaks = NULL, retrieval="standard")
#' @param spec A \code{RmbSpectraSet} object, representing a compound with multiple spectra.
#' @param aggregated An aggregate peak table where the peaks are extracted from.
#' @param msmsdata A \code{RmbSpectrum2} object from the \code{spec} spectra set, representing a single spectrum to give a record.
#' @param ac_ms,ac_lc Information for the AC\$MASS_SPECTROMETRY and
#' AC\$CHROMATOGRAPHY fields in the MassBank record, created by
#' \code{gatherCompound} and then fed into \code{gatherSpectrum}.
#' @param additionalPeaks If present, a table with additional peaks to add into the spectra.
#' 		As loaded with \code{\link{addPeaks}}.
#' @param retrieval A value that determines whether the files should be handled either as "standard",
#' if the compoundlist is complete, "tentative", if at least a formula is present or "unknown"
#' if the only know thing is the m/z
#' @return \code{gatherCompound} returns a list of tree-like MassBank data
#' blocks. \code{gatherSpectrum} returns one single MassBank data block or
#' \code{NA} if no useful peak is in the spectrum. 
#' @note Note that the global table \code{additionalPeaks} is also used as an
#' additional source of peaks.
#' @author Michael Stravs
#' @seealso \code{\link{mbWorkflow}}, \code{\link{buildRecord}}
#' @references MassBank record format:
#' \url{http://www.massbank.jp/manuals/MassBankRecord_en.pdf}
#' @examples \dontrun{
#'      myspectrum <- w@@spectra[[1]]
#' 		massbankdata <- gatherCompound(myspectrum, w@@aggregated)
#' 		# Note: ac_lc and ac_ms are data blocks usually generated in gatherCompound and
#' 		# passed on from there. The call below gives a relatively useless result :)
#' 		ac_lc_dummy <- list()
#' 		ac_ms_dummy <- list() 
#' 		justOneSpectrum <- gatherSpectrum(myspectrum, myspectrum@@child[[2]],
#' 			ac_ms_dummy, ac_lc_dummy, w@@aggregated)
#' }
#' 
#' 
#' @export
getAnalyticalInfo <- function(cpd = NULL)
{
  ai <- list()
	# define positive or negative, based on processing mode.
  if(!is.null(cpd))
	  mode <- .ionModes[[cpd@mode]]
	
  # again, these constants are read from the options:
  ai[['AC$INSTRUMENT']] <- getOption("RMassBank")$annotations$instrument
  ai[['AC$INSTRUMENT_TYPE']] <- getOption("RMassBank")$annotations$instrument_type

	# for format 2.01
	ac_ms <- list();
	ac_ms[['MS_TYPE']] <- getOption("RMassBank")$annotations$ms_type
	ac_ms[['ION_MODE']] <- mode
	ac_ms[['IONIZATION']] <- getOption("RMassBank")$annotations$ionization
	
	# This list could be made customizable.
	ac_lc <- list();
  if(!is.null(cpd))
	  rt  <- cpd@parent@rt / 60
	ac_lc[['COLUMN_NAME']] <- getOption("RMassBank")$annotations$lc_column
	ac_lc[['FLOW_GRADIENT']] <- getOption("RMassBank")$annotations$lc_gradient
	ac_lc[['FLOW_RATE']] <- getOption("RMassBank")$annotations$lc_flow
	ac_lc[['RETENTION_TIME']] <- sprintf("%.3f min", rt)  
	lc_solvents <- getOption("RMassBank")$annotations$lc_solvents
	ac_lc[['SOLVENT A']] <- lc_solvents$lc_solvent_a
	ac_lc[['SOLVENT B']] <- lc_solvents$lc_solvent_b
	if(length(lc_solvents) > 2)
		ac_lc[['SOLVENT C']] <- lc_solvents$lc_solvent_c
	
	ac_ms <- .addGenericInfo(ac_ms, getOption('RMassBank')$annotations,
	  search_string="^AC\\$MASS_SPECTROMETRY_")
	ac_lc <- .addGenericInfo(ac_lc, getOption('RMassBank')$annotations,
	  search_string="^AC\\$CHROMATOGRAPHY_")
	return(list( ai=ai, ac_lc=ac_lc, ac_ms=ac_ms))
}


#' @rdname buildRecord
setMethod("buildRecord", "RmbSpectrum2", function(o, ..., cpd = NULL, mbdata = list(), analyticalInfo = list(), additionalPeaks = NULL)
      .buildRecord.RmbSpectrum2(spectrum = o, cpd=cpd, mbdata=mbdata, analyticalInfo=analyticalInfo, additionalPeaks=additionalPeaks, ...)
)

.buildRecord.RmbSpectrum2 <- function(spectrum, ..., cpd = NULL, mbdata = list(), analyticalInfo = list(), additionalPeaks = NULL)
{

	if(length(analyticalInfo$ac_ms) > 0)
		ac_ms=analyticalInfo$ac_ms
	else
		ac_ms=list()

	if(length(analyticalInfo$ac_lc) > 0)
		ac_lc=analyticalInfo$ac_lc
	else
		ac_lc=list()

	if(length(mbdata) == 0)
	{
		if(is.null(cpd))
			mbdata <- gatherDataMinimal.spectrum(spectrum)
		else
			mbdata <- gatherDataMinimal.cpd(cpd)
	}

	if(length(analyticalInfo$ai) > 0)
		mbdata <- c(mbdata, analyticalInfo$ai)

	# If the spectrum is not filled, return right now. All "NA" spectra will
	# not be treated further.
	# If step 2 was not performed, instead, spectrum@ok is empty and we want to export it, so proceed.
	if(length(spectrum@ok) > 0)
	{
		if(spectrum@ok == FALSE)
			return(NA)
	}
	# get data
	scan <- spectrum@acquisitionNum


	# Further fill the ac_ms datasets, and add the ms$focused_ion with spectrum-specific data:
	ac_ms[['FRAGMENTATION_MODE']] <- spectrum@info$mode
	#ac_ms['PRECURSOR_TYPE'] <- precursor_types[spec$mode]
	if(length(spectrum@info$ce) > 0)
		ac_ms[['COLLISION_ENERGY']] <- spectrum@info$ce
	else
		ac_ms[['COLLISION_ENERGY']] <- spectrum@collisionEnergy
	ac_ms[['RESOLUTION']] <- spectrum@info$res

	# Calculate exact precursor mass with Rcdk, and find the base peak from the parent
	# spectrum. (Yes, that's what belongs here, I think.)

	ms_fi <- list()
	if(!is.null(cpd))
	{
		ms_fi[['BASE_PEAK']] <- round(mz(cpd@parent)[which.max(intensity(cpd@parent))],4)
		ms_fi[['PRECURSOR_M/Z']] <- round(cpd@mz,4)
		ms_fi[['PRECURSOR_TYPE']] <- .precursorTypes[cpd@mode]

		if(all(!is.na(spectrum@precursorIntensity), 
		   spectrum@precursorIntensity != 0, 
		   spectrum@precursorIntensity != 100, na.rm = TRUE))
			ms_fi[['PRECURSOR_INTENSITY']] <- round(spectrum@precursorIntensity, 2)
	}

	# Add scan range to AC$MS, if present
	if (all(c("scanWindowUpperLimit", "scanWindowLowerLimit") %in%
	  names(spectrum@info))) {
		ac_ms[['MASS_RANGE_M/Z']] <- paste(
		  floor(spectrum@info$scanWindowLowerLimit),
		  ceiling(spectrum@info$scanWindowUpperLimit),
		  sep='-')
	}

	# Create the "lower part" of the record.  

	# Add the AC$MS, AC$LC info.
	if(getOption("RMassBank")$use_version == 2)
	{
		if(length(ac_ms) >0)
			mbdata[["AC$MASS_SPECTROMETRY"]] <- ac_ms
		if(length(ac_lc) >0)
			mbdata[["AC$CHROMATOGRAPHY"]] <- ac_lc
	}
	else
	{
		# Fix for MassBank data format 1, where ION_MODE must be renamed to MODE
		ac <- c(ac_ms, ac_lc)
		if(length(ac) > 0)
		{
			mbdata[["AC$ANALYTICAL_CONDITION"]] <- ac
			names(mbdata[["AC$ANALYTICAL_CONDITION"]])[[
			  which(names(mbdata[["AC$ANALYTICAL_CONDITION"]]) == "ION_MODE")
			]] <- "MODE"
		}
	}
	# Add the MS$FOCUSED_ION info.
	if(length(ms_fi) > 0)
		mbdata[["MS$FOCUSED_ION"]] <- ms_fi

	## The SPLASH is a hash value calculated across all peaks
	## http://splash.fiehnlab.ucdavis.edu/
	## Has to be temporarily added as "PK$SPLASH" in the "lower" part
	## of the record, but will later be moved "up" when merging parts in compileRecord()  

	# the data processing tag :)
	# Change by Tobias:
	# I suggest to add here the current version number of the clone due to better distinction between different makes of MB records
	# Could be automatised from DESCRIPTION file?
	if(getOption("RMassBank")$use_rean_peaks)
		processingComment <- list("REANALYZE" = "Peaks with additional N2/O included")
	else
		processingComment <- list()
	mbdata[["MS$DATA_PROCESSING"]] <- c(
	  getOption("RMassBank")$annotations$ms_dataprocessing,
	  processingComment,
	  list("WHOLE" = paste("RMassBank", packageVersion("RMassBank")))
	)

	if(length(spectrum@info$ces) > 0)
		mbdata[['RECORD_TITLE_CE']] <- spectrum@info$ces
	else
		mbdata[['RECORD_TITLE_CE']] <- spectrum@collisionEnergy

	# Mode of relative scan calculation: by default it is calculated relative to the
	# parent scan. If a corresponding option is set, it will be calculated from the first
	# present child scan in the list.

	if(!is.null(cpd))
	{
		relativeScan <- "fromParent"
		if(!is.null(getOption("RMassBank")$recomputeRelativeScan))
			if(getOption("RMassBank")$recomputeRelativeScan == "fromFirstChild")
				relativeScan <- "fromFirstChild"
		if(relativeScan == "fromParent")
			subscan <- spectrum@acquisitionNum - cpd@parent@acquisitionNum #relative scan
		else if(relativeScan == "fromFirstChild")
		{
			firstChild <- min(unlist(lapply(cpd@children,function(d) d@acquisitionNum)))
			subscan <- spectrum@acquisitionNum - firstChild + 1
		}
	}


	# Here is the right place to fix the name of the INTERNAL ID field.
	if(!is.null(getOption("RMassBank")$annotations$internal_id_fieldname))
	{
		id.col <- which(names(mbdata[["COMMENT"]]) == "ID")
		if(length(id.col) > 0)
		{
			names(mbdata[["COMMENT"]])[[id.col]] <-
			getOption("RMassBank")$annotations$internal_id_fieldname
		}
	}
	# get mode parameter (for accession number generation) depending on version 
	# of record definition
	# Generate the title and then delete the temprary RECORD_TITLE_CE field used before
	mbdata[["RECORD_TITLE"]] <- .parseTitleString(mbdata)
	mbdata[["RECORD_TITLE_CE"]] <- NULL
	userSettings = getOption("RMassBank")
	# Include project tag, if present
	if("project" %in% names(userSettings))
	{
		mbdata[["PROJECT"]] <- userSettings$project
	}
	# Use 'simple', 'standard' or 'selfDefined' accessionBuilder
	# depending on user input
	if("accessionBuilderType" %in% names(userSettings))
	{
		assert_that(userSettings$accessionBuilderType %in% c(
		  "standard", "simple", "selfDefined"),
		  msg=paste("accessionNumberType must be one of",
		  "'standard', 'simple' or 'selfDefined'"))
		mbdata[['ACCESSION']] <- switch(
		  userSettings$accessionBuilderType,
		  simple = .simpleAccessionBuilder(subscan),
		  standard = .standardAccessionBuilder(cpd, subscan),
		  selfDefined = .selfDefinedAccessionBuilder(cpd, spectrum,
		  subscan)
		)
	}
	else
	{
		mbdata[['ACCESSION']] <- .standardAccessionBuilder(cpd, subscan)
	}

	spectrum@info <- mbdata

	spectrum <- renderPeaks(spectrum, cpd=cpd, additionalPeaks=additionalPeaks, ...)

	return(spectrum)
}

.simpleAccessionBuilder <- function(subscan)
{
	userSettings = getOption("RMassBank")
	assert_that('accessionNumberStart' %in% names(userSettings),
	  msg=paste("accessionBuilderType is 'simple', but",
	  "accessionNumberStart is not provided.",
	  "You may set the value of accessionBuilderType to",
	  "'standard' or 'selfDefined' to use a different accessionBuilder.", 
	  "For detailed explanations of accessionBuilders, check out the",
	  "'Settings' section of the RMassBank vignette."))
	sprintf('%s%06d', userSettings$annotations$entry_prefix,
	  subscan + userSettings$accessionNumberStart)
}

.standardAccessionBuilder <- function(cpd, subscan)
{
	userSettings = getOption("RMassBank")
	shift <- userSettings$accessionNumberShifts[[cpd@mode]]
	sprintf("%s%04d%02d", userSettings$annotations$entry_prefix,
	  as.numeric(cpd@id), subscan+shift)
}
	
.selfDefinedAccessionBuilder <- function(cpd, spectrum, subscan)
{
	#This is a wrapper for the user-defined accessionBuilder
	userSettings = getOption("RMassBank")
	assert_that("accessionBuilderFile" %in% names(userSettings),
	  msg=paste("accessionBuilderType is 'selfDefined', but",
	  "accessionBuilderFile is not provided.",
	  "You may set the value of accessionBuilderType",
	  "to 'standard' or 'simple' to use a different accessionBuilder.", 
	  "For detailed explanations of accessionBuilders, check out the",
	  "'Settings' section of the RMassBank vignette."))
	source(userSettings$accessionBuilderFile)
	#The file must contain a function called 'accessionBuilder'
	#with arguments cpd, spectrum, subscan
	assert_that(exists("accessionBuilder"),
	  msg=paste('No accessionBuilder defined in',
	  userSettings$accessionBuilderFile))
	assert_that(class(accessionBuilder)=='function',
	  msg='accessionBuilder must be a function')
	assert_that(has_args(accessionBuilder,
	  c('cpd', 'spectrum', 'subscan'), exact=TRUE),
	  msg=paste('accessionBuilder must have function arguments',
	  'cpd, spectrum, subscan in this order'))
	accessionBuilder(cpd, spectrum, subscan)
}

renderPeaks <- function(spectrum, ..., cpd = NULL, additionalPeaks = NULL)
{
	# Select all peaks which belong to this spectrum (correct cpdID and scan no.)
	# from peaksOK
	# Note: Here and below it would be easy to customize the source of the peaks.
	# Originally the peaks came from msmsdata$childFilt, and the subset
	# was used where dppm == dppmBest (because childFilt still contains multiple formulas)
	# per peak.
  spectrum <- .fillSlots(spectrum, c("good", "dppm", "dppmBest", "mzCalc", "formula", "formulaCount"))
  peaks <- getData(spectrum)
  property(spectrum, "best", addNew=TRUE, "logical") <- (peaks$good %in% TRUE) & (peaks$dppm == peaks$dppmBest)
  spectrum <- normalize(spectrum, 999, slot="intrel", ...)
  peaks <- getData(selectPeaks(spectrum, ...))
  # filterOK is the final criterion for selection, it includes both reanalyzed and original matches.
  # If there was no peak filtering performed, use best | matchedReanalysis (which gets both regular and reanalyzed matches)
  # To get peaks without the reanalyzed matches, use best
  # rawOK gives the unfiltered, not denoised spectrum.
  # Any other condition can be used (also for example intrel > 50)
  
  if(!getOption("RMassBank")$use_rean_peaks)
    peaks <- peaks[peaks$formulaSource == "analyze",,drop=FALSE]
  
  
	# No peaks? Aha, bye
	if(nrow(peaks) == 0)
		return(NA)
	
	# Calculate relative intensity and make a formatted m/z to use in the output
	# (mzSpec, for "spectrum")
	#peaks$intrel <- floor(peaks$intensity / max(peaks$intensity) * 999)

	# reorder peaks after addition of the reanalyzed ones


	# copy the peak table to the annotation table. (The peak table will then be extended
	# with peaks from the global "additional_peaks" table, which can be used to add peaks
	# to the spectra by hand.

	
	# Here add the additional peaks if there are any for this compound!
	# They are added without any annotation.
  if(is.null(additionalPeaks))
    additionalPeaks <- data.frame()
  
  if(!is.null(cpd))
  {
    if(ncol(additionalPeaks) > 0)
    {
      # select the peaks from the corresponding spectrum which were marked with "OK=1" in the table.
      spec_add_peaks <- additionalPeaks[
          (!is.na(additionalPeaks$OK)) &
              (additionalPeaks$OK == 1) & 
              (as.character(additionalPeaks$cpdID) == cpd@id) &
              (additionalPeaks$scan == spectrum@acquisitionNum),
          c("mzFound", "intensity")]
      # If there are peaks to add:
      if(nrow(spec_add_peaks)>0)
      {
        colnames(spec_add_peaks) <- c("mz", "intensity")
        # bind tables together. First add in NA fillers for all columns not in spec_add_peaks
        for(column in setdiff(colnames(peaks), colnames(spec_add_peaks)))
          spec_add_peaks[,column] <- new(class(peaks[,column]), NA)
        #print(spec_add_peaks)
        peaks <- rbind(peaks, spec_add_peaks[,colnames(peaks),drop=FALSE])
        # recalculate rel.int.  and reorder list
      }
    }
  }
  
  # recalculate relative intensity with the newly added peaks. Note that possibly this leads to spectra
  # with intrel < what was specified in the original filter!
  peaks$intrel <- floor(peaks$intensity / max(peaks$intensity) * 999)
  peaks <- peaks[order(peaks$mz),]
  
	# build annotation
  spectrum <- setData(spectrum, peaks) 

	return(spectrum)
}

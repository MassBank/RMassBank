## For generating the NAMESPACE
#' @import mzR

#' @import Rcpp 
#' @import XML 
#' @import methods 
#' @import mzR 
#' @import rcdk 
#' @import rjson 
#' @import yaml 
#' @import digest
NULL # This is required so that roxygen knows where the first manpage starts


# # importClassesFrom mzR ## Causes error 
# # importMethodsFrom mzR 

#' Extract MS/MS spectra for specified precursor
#' 
#' Extracts MS/MS spectra from LC-MS raw data for a specified precursor, specified
#' either via the RMassBank compound list (see \code{\link{loadList}}) or via a mass.
#' 
#' Different versions of the function get the data from different sources. Note that 
#' 		findMsMsHR and findMsMsHR.direct differ mainly in that findMsMsHR opens a file
#' 		whereas findMsMs.direct uses an open file handle - both are intended to be used
#' 		in a full process which involves compound lists etc. In contrast, findMsMsHR.mass
#' 		is a low-level function which uses the mass directly for lookup and is intended for
#' 		use as a standalone function in unrelated applications.
#' 
#' @note \code{findMsMs.direct} is deactivated
#' 
## # @usage findMsMsHR(fileName, cpdID, mode="pH",confirmMode =0, useRtLimit = TRUE,
## # 		ppmFine = getOption("RMassBank")$findMsMsRawSettings$ppmFine,
## # 		mzCoarse = getOption("RMassBank")$findMsMsRawSettings$mzCoarse,
## # 		fillPrecursorScan = getOption("RMassBank")$findMsMsRawSettings$fillPrecursorScan,
## # 		rtMargin = getOption("RMassBank")$rtMargin,
## # 		deprofile = getOption("RMassBank")$deprofile,
## # 		headerCache = NULL,
## # 		peaksCache = NULL)
## # 		
## # findMsMsHR.mass(msRaw, mz, limit.coarse, limit.fine, rtLimits = NA, maxCount = NA,
## # 		headerCache = NULL, fillPrecursorScan = FALSE,
## # 		deprofile = getOption("RMassBank")$deprofile, peaksCache = NULL, cpdID = NA)
#' 
#' 
#' @aliases findMsMsHR.mass findMsMsHR
#' @param fileName The file to open and search the MS2 spectrum in.
#' @param msRaw The opened raw file (mzR file handle) to search the MS2 spectrum in. Specify either this
#' 			or \code{fileName}.
#' @param cpdID The compound ID in the compound list (see \code{\link{loadList}})
#' 			to use for formula lookup. Note: In \\code{findMsMsHR.mass}, this is entirely optional and
#' 			used only in case a warning must be displayed; compound lookup is done via mass only. 
#' @param mz The mass to use for spectrum search.
#' @param ppmFine The limit in ppm to use for fine limit (see below) calculation.
#' @param mzCoarse The coarse limit to use for locating potential MS2 scans:
#'			this tolerance is used when finding scans with a suitable precursor
#' 			ion value.  
#' @param limit.fine The fine limit to use for locating MS2 scans: this tolerance
#' 			is used when locating an appropriate analyte peak in the MS1 precursor
#' 			spectrum. 
#' @param limit.coarse Parameter in \code{findMsMsHR.mass} corresponding to \code{mzCoarse}.
#' 			(The parameters are distinct to clearly conceptually distinguish findMsMsHR.mass
#' 			(a standalone useful function) from the cpdID based functions (workflow functions).)
#' @param mode The processing mode (determines which ion/adduct is searched):
#' 			\code{"pH", "pNa", "pM", "pNH4", "mH", "mM", "mFA"} for different ions 
#' 			([M+H]+, [M+Na]+, [M]+, [M+NH4]+, [M-H]-, [M]-, [M+FA]-). 
#' @param confirmMode Whether to use the highest-intensity precursor (=0), second-
#' 			highest (=1), third-highest (=2)...
#' @param useRtLimit Whether to respect retention time limits from the compound list.
#' @param rtLimits \code{c(min, max)}: Minimum and maximum retention time to use
#' 			when locating the MS2 scans. 
#' @param headerCache If present, the complete \code{mzR::header(msRaw)}. Passing
#' 			this value is useful if spectra for multiple compounds should be 
#' 			extracted from the same mzML file, since it avoids getting the data
#' 			freshly from \code{msRaw} for every compound.
#' @param peaksCache If present, the complete output of \code{mzR::peaks(msRaw)}. This speeds up the lookup
#' 			if multiple compounds should be searched in the same file.
#' @param maxCount The maximal number of spectra groups to return. One spectra group
#' 			consists of all data-dependent scans from the same precursor whose precursor
#' 			mass matches the specified search mass.
#' @param fillPrecursorScan If \code{TRUE}, the precursor scan will be filled from MS1 data.
#' 			To be used for data where the precursor scan is not stored in the raw data.
#' @param rtMargin	The retention time tolerance to use.
#' @param deprofile	Whether deprofiling should take place, and what method should be
#' 			used (cf. \code{\link{deprofile}}) 
#' @param diaWindows A data frame with columns \code{precursorMz}, \code{mzMin}, \code{mzMax} which specifies the precursor and 
#'      window size of each window for DIA acquisition.
#' @return	An \code{RmbSpectraSet} (for \code{findMsMsHR}). Contains parent MS1 spectrum (\code{@@parent}), a block of dependent MS2 spectra ((\code{@@children})
#' 			and some metadata (\code{id},\code{mz},\code{name},\code{mode} in which the spectrum was acquired.
#' 
#' 			For \code{findMsMsHR.mass}: a list of \code{RmbSpectraSet}s as defined above, sorted
#' 			by decreasing precursor intensity.
#' 
#' @examples \dontrun{
#' 			loadList("mycompoundlist.csv")
#' 			# if Atrazine has compound ID 1:
#' 			msms_atrazine <- findMsMsHR(fileName = "Atrazine_0001_pos.mzML", cpdID = 1, mode = "pH")
#' 			# Or alternatively:
#' 			msRaw <- openMSfile("Atrazine_0001_pos.mzML")
#' 			msms_atrazine <- findMsMsHR(msRaw=msRaw, cpdID = 1, mode = "pH")
#' 			# Or directly by mass (this will return a list of spectra sets):
#' 			mz <- findMz(1)$mzCenter
#' 			msms_atrazine_all <- findMsMsHR.mass(msRaw, mz, 1, ppm(msRaw, 10, p=TRUE))
#' 			msms_atrazine <- msms_atrazine_all[[1]]
#' }
#' @author Michael A. Stravs, Eawag <michael.stravs@@eawag.ch>
#' @seealso findEIC
#' @export
findMsMsHR <- function(fileName = NULL, msRaw = NULL, cpdID, mode="pH",confirmMode =0, useRtLimit = TRUE,
		ppmFine = getOption("RMassBank")$findMsMsRawSettings$ppmFine,
		mzCoarse = getOption("RMassBank")$findMsMsRawSettings$mzCoarse,
		fillPrecursorScan = getOption("RMassBank")$findMsMsRawSettings$fillPrecursorScan,
		rtMargin = getOption("RMassBank")$rtMargin,
		deprofile = getOption("RMassBank")$deprofile,
		headerCache = NULL,
		peaksCache = NULL,
		enforcePolarity = getOption("RMassBank")$enforcePolarity,
		diaWindows = getOption("RMassBank")$findMsMsRawSettings$diaWindows)
{
  retrieval <- findLevel(cpdID,TRUE)
  # old behaviour: do not enforce polarity
  if(is.null(enforcePolarity))
    enforcePolarity <- FALSE
  
  if(enforcePolarity)
    polarity <- .polarity[[mode]]
  else
    polarity <- NA
	# access data directly for finding the MS/MS data. This is done using
	# mzR.
	if(!is.null(fileName) & !is.null(msRaw))
		stop("Both MS raw data and MS filename given. Only one can be handled at the same time.")
	if(!is.null(fileName))
		msRaw <- openMSfile(fileName)
	
	mzLimits <- findMz(cpdID, mode, retrieval=retrieval)
	mz <- mzLimits$mzCenter
	limit.fine <- ppm(mz, ppmFine, p=TRUE)
	if(!useRtLimit)
		rtLimits <- NA
	else
	{
		dbRt <- findRt(cpdID)
		rtLimits <- c(dbRt$RT - rtMargin, dbRt$RT + rtMargin) * 60
	}
	spectra <- findMsMsHR.mass(msRaw, mz, mzCoarse, limit.fine, rtLimits, confirmMode + 1,headerCache
			,fillPrecursorScan, deprofile, peaksCache, cpdID, polarity=polarity, diaWindows = diaWindows)
	# check whether a) spectrum was found and b) enough spectra were found
	if(length(spectra) < (confirmMode + 1))
		sp <- new("RmbSpectraSet", found=FALSE)
	else
		sp <- spectra[[confirmMode + 1]]
	
	#sp@mz <- mzLimits
	sp@id <- as.character(as.integer(cpdID))
	sp@name <- findName(cpdID)
  if(retrieval == "standard")
    sp@smiles <- findSmiles(cpdID)
    ENV <- environment()
	if(retrieval == "unknown"){
        sp@formula <- ""
    } else{
        sp@formula <- findFormula(cpdID, retrieval=retrieval)
    }
	sp@mode <- mode
  
  
  # Overwrite the polarity with a value we generate, so it's consistent.
  # Some mzML files give only -1 as a result for polarity, which is useless for us
  sp@parent@polarity <- .polarity[[sp@mode]]
  for(n in seq_len(length(sp@children)))
  {
    sp@children[[n]]@polarity <- .polarity[[sp@mode]]
  }
	
	# If we had to open the file, we have to close it again
	if(!is.null(fileName))
		mzR::close(msRaw)
	
	return(sp)
}



#' @describeIn findMsMsHR A submethod of find MsMsHR that retrieves basic spectrum data 
#' @export
findMsMsHR.mass <- function(msRaw, mz, limit.coarse, limit.fine, rtLimits = NA, maxCount = NA,
		headerCache = NULL, fillPrecursorScan = FALSE,
		deprofile = getOption("RMassBank")$deprofile, peaksCache = NULL, cpdID = NA,
		polarity = NA,
		diaWindows = getOption("RMassBank")$findMsMsRawSettings$diaWindows)
{
	eic <- findEIC(msRaw, mz, limit.fine, rtLimits, headerCache=headerCache, 
			peaksCache=peaksCache)
	#	if(!is.na(rtLimits))
	#	{  
	#		eic <- subset(eic, rt >= rtLimits[[1]] & rt <= rtLimits[[2]])
	#	}
	if(!is.null(headerCache))
		headerData <- headerCache
	else
		headerData <- as.data.frame(header(msRaw))
	
	if(!is.na(polarity))
	  headerData <- headerData[headerData$polarity == polarity,,drop=FALSE]
	
	
	###If no precursor scan number, fill the number
	if(length(unique(headerData$precursorScanNum)) == 1){
		fillPrecursorScan <- TRUE
	}

	if(fillPrecursorScan == TRUE)
	{
		# reset the precursor scan number. first set to NA, then
		# carry forward the precursor scan number from the last parent scan
		headerData$precursorScanNum <- NA
		headerData[which(headerData$msLevel == 1),"precursorScanNum"] <-
				headerData[which(headerData$msLevel == 1),"acquisitionNum"]
		headerData[,"precursorScanNum"] <- .locf(headerData[,"precursorScanNum"])
		# Clear the actual MS1 precursor scan number again
		headerData[which(headerData$msLevel == 1),"precursorScanNum"] <- 0
		# Remove precursors which are still NA in precursor scan num.
		# This removes a bug when filling precursor if the first scan(s) are MS2 before a
		# MS1 scan appears. The resulting NA values in precursorScanNum are problematic downstream.
		headerData <- headerData[!is.na(headerData$precursorScanNum),]
	}
	# bugfix 201803: PRM scans that were performed before the first full scan (found in some files)
	headerData <- headerData[!((headerData$msLevel == 2) &
                                   (is.na(headerData$precursorScanNum))),,
                                 drop = FALSE]
	# Find MS2 spectra with precursors which are in the allowed 
	# scan filter (coarse limit) range; which to get rid of NAs
	if(!is.null(diaWindows))
	{
	  rmb_log_info("using diaWindows")
	  window <- which((diaWindows$mzMin < mz) & (diaWindows$mzMax >= mz))
	  if(length(window) > 1)
	  {
	    warning("Compound mass lies in two DIA windows - this is still in test phase!")
	    diaWin <- diaWindows[window,]
	    diaWin$internality <- min(abs(mz - diaWin$mzMin), abs(mz - diaWin$mzMax))
	    window <- window[which.max(diaWin$internality)]
	  }
	    
	  precursor <- diaWindows$precursorMz[window]
	  # this allows for minimal numeric errors in precursor mass:
	  findValidPrecursors <- headerData[abs(headerData$precursorMZ - precursor) < 0.1,]
	}
	else
	{
	  findValidPrecursors <- headerData[which(
	    (headerData$precursorMZ > (mz - limit.coarse)) &
	      (headerData$precursorMZ < (mz + limit.coarse))),]
	# Find the precursors for the found spectra
	}
	validPrecursors <- unique(findValidPrecursors$precursorScanNum)
	validPrecursors <- validPrecursors[!is.na(validPrecursors)]
	# check whether the precursors are real: must be within fine limits!
	# previously even "bad" precursors were taken. e.g. 1-benzylpiperazine
	which_OK <- lapply(validPrecursors, function(pscan)
			{
	     # Debugging for figuring out the right way to fix a bug
	      # print("AcquisitionNum")
	      # print(pscan)
	      # i <- which(headerData$acquisitionNum == pscan)
	      # print("Which:")
	      # print(i)
	      # print("SeqNum:")
	      # print(headerData$seqNum[[i]])
	      
				pplist <- as.data.frame(
						mzR::peaks(msRaw, headerData$seqNum[[which(headerData$acquisitionNum == pscan)]]))
				colnames(pplist) <- c("mz","int")
				pplist <- pplist[(pplist$mz >= mz -limit.fine)
								& (pplist$mz <= mz + limit.fine),]
				if(nrow(pplist) > 0)
					return(TRUE)
				return(FALSE)
			})
	validPrecursors <- validPrecursors[which(which_OK==TRUE)]
	if(length(validPrecursors) == 0){
		if(!is.na(cpdID))
			warning(paste0("No precursor was detected for compound, ", cpdID, " with m/z ", mz, ". Please check the mass and retention time window."))
		else
			warning(paste0("No precursor was detected for m/z ", mz, ". Please check the mass and retention time window."))
	}
	# Crop the "EIC" to the valid precursor scans
	eic <- eic[eic$scan %in% validPrecursors,]
	# Order by intensity, descending
	eic <- eic[order(eic$intensity, decreasing=TRUE),]
	if(nrow(eic) == 0)
		return(list(
						new("RmbSpectraSet",
								found=FALSE)))
	if(!is.na(maxCount))
	{
		spectraCount <- min(maxCount, nrow(eic))
		eic <- eic[1:spectraCount,]
	}
	# Construct all spectra groups in decreasing intensity order
	spectra <- lapply(eic$scan, function(masterScan)
	{
		masterHeader <- headerData[headerData$acquisitionNum == masterScan,]
		
		if(is.null(diaWindows))
		{
			childHeaders <- headerData[
			  which(headerData$precursorScanNum == masterScan 
			  & headerData$precursorMZ > (mz - limit.coarse) 
			  & headerData$precursorMZ < (mz + limit.coarse)) , ,
			  drop = FALSE]
		}
		else
		{
			childHeaders <- headerData[which(headerData$precursorScanNum == masterScan), drop = FALSE]
			childHeaders <- childHeaders[window, drop = FALSE]
		}
		
		# Fix 9.10.17: headers now include non-numeric columns, leading to errors in data conversion.
		# Remove non-numeric columns
		headerCols <- colnames(masterHeader)
		headerCols <- headerCols[unlist(lapply(headerCols, function(col) is.numeric(masterHeader[,col])))]
		masterHeader <- masterHeader[,headerCols,drop=FALSE]
		childHeaders <- childHeaders[,headerCols,drop=FALSE]
		
		childScans <- childHeaders$seqNum
		
		msPeaks <- mzR::peaks(msRaw, masterHeader$seqNum)
		# if deprofile option is set: run deprofiling
		deprofile.setting <- deprofile
		if(!is.na(deprofile.setting))
			msPeaks <- deprofile.scan(
					msPeaks, method = deprofile.setting, noise = NA, colnames = FALSE
			)
		colnames(msPeaks) <- c("mz","int")
		
		msmsSpecs <- apply(childHeaders, 1, function(line)
		{
			pks <- mzR::peaks(msRaw, line["seqNum"])
			
			if(!is.na(deprofile.setting))
			{
				pks <- deprofile.scan(
				  pks, method = deprofile.setting, noise = NA
				  , colnames = FALSE)
			}
			pks_mz <- pks[,1]
			pks_intensity <- pks[,2]
			scanWindowLowerLimit <- line["scanWindowLowerLimit"]
			scanWindowUpperLimit <- line["scanWindowUpperLimit"]
			limits <- list(
			  scanWindowLowerLimit=scanWindowLowerLimit,
			  scanWindowUpperLimit=scanWindowUpperLimit
			)
			if(!anyNA(limits)) {
				check_mz <- function(m) {isTRUE(
					m > scanWindowLowerLimit &&
					m < scanWindowUpperLimit
				)}
				in_range <- sapply(pks_mz, check_mz)
				if (!all(in_range)) {
					outliers <- pks[!in_range, ]
					warning(paste('There were',
					  nrow(outliers),
					  'peaks out of mass range.'))
				}
			}
			new("RmbSpectrum2",
			  mz = pks_mz,
			  intensity = pks_intensity,
			  precScanNum = as.integer(line["precursorScanNum"]),
			  precursorMz = line["precursorMZ"],
			  precursorIntensity = line["precursorIntensity"],
			  precursorCharge = as.integer(line["precursorCharge"]),
			  collisionEnergy = line["collisionEnergy"],
			  tic = line["totIonCurrent"],
			  peaksCount = line["peaksCount"],
			  rt = line["retentionTime"],
			  acquisitionNum = as.integer(line["seqNum"]),
			  centroided = TRUE,
			  polarity = as.integer(line["polarity"]),
			  info = lapply(limits, unname)
			)
		})
		msmsSpecs <- as(do.call(c, msmsSpecs), "SimpleList")
		
		# build the new objects
		masterSpec <- new("Spectrum1",
				mz = msPeaks[,"mz"],
				intensity = msPeaks[,"int"],
				polarity = as.integer(masterHeader$polarity),
				peaksCount = as.integer(masterHeader$peaksCount),
				rt = masterHeader$retentionTime,
				acquisitionNum = as.integer(masterHeader$seqNum),
				tic = masterHeader$totIonCurrent,
				centroided = TRUE
				)
				
		spectraSet <- new("RmbSpectraSet",
				parent = masterSpec,
				children = msmsSpecs,
				found = TRUE,
				#complete = NA,
				#empty = NA,
				#formula = character(),
				mz = mz
				#name = character(),
				#annotations = list()
				)
		return(spectraSet)
			})
	names(spectra) <- eic$acquisitionNum
	return(spectra)
}


#' Discontinued: find MS/MS spectrum from open raw file
#' 
#' This interface has been discontinued. \code{\link{findMsMsHR}} now supports the same parameters (use named
#' parameters).
#' 
#' @param msRaw x
#' @param cpdID x
#' @param mode x
#' @param confirmMode x 
#' @param useRtLimit x
#' @param ppmFine x
#' @param mzCoarse x
#' @param fillPrecursorScan x 
#' @param rtMargin x
#' @param deprofile x
#' @param headerCache x
#' @return an error
#' 
#' @author stravsmi
#' @export
findMsMsHR.direct <- function(msRaw, cpdID, mode = "pH", confirmMode = 0, useRtLimit = TRUE, 
			ppmFine = getOption("RMassBank")$findMsMsRawSettings$ppmFine,
			mzCoarse = getOption("RMassBank")$findMsMsRawSettings$mzCoarse,
			fillPrecursorScan = getOption("RMassBank")$findMsMsRawSettings$fillPrecursorScan,
			rtMargin = getOption("RMassBank")$rtMargin,
			deprofile = getOption("RMassBank")$deprofile,
      headerCache = NULL)
{
	stop("Support for this interface has been discontinued. Use findMsMsHR with the same parameters instead (use named parameter msRaw)")
}

#' Read in mz-files using XCMS
#' 
#' Picks peaks from mz-files and returns the pseudospectra that CAMERA creates with the help of XCMS
#'
#' @aliases findMsMsHRperxcms.direct findMsMsHRperxcms
#' @param fileName The path to the mz-file that should be read
#' @param cpdID The compoundID(s) of the compound that has been used for the file
#' @param mode The ionization mode that has been used for the spectrum represented by the peaklist
#' @param findPeaksArgs A list of arguments that will be handed to the xcms-method findPeaks via do.call
#' @param plots A parameter that determines whether the spectra should be plotted or not
#' @param MSe A boolean value that determines whether the spectra were recorded using MSe or not
#' @return The spectra generated from XCMS
#' @seealso \code{\link{msmsWorkflow}} \code{\link{toRMB}}
#' @author Erik Mueller
#' @examples \dontrun{
#' 		fileList <- list.files(system.file("XCMSinput", package = "RMassBank"), "Glucolesquerellin", full.names=TRUE)[3]
#'		loadList(system.file("XCMSinput/compoundList.csv",package="RMassBank"))
#'      psp <- findMsMsHRperxcms(fileList,2184)
#' }
#' @export
findMsMsHRperxcms <- function(fileName, cpdID, mode="pH", findPeaksArgs = NULL, plots = FALSE, MSe = FALSE){
	
	# Find mz
	mzLimits <- findMz(cpdID, mode)
	mz <- mzLimits$mzCenter
	
	
	# If there are more files than cpdIDs
	if(length(fileName) > 1){
		fspectra <- list()
		
		for(i in 1:length(fileName)){
			fspectra[[i]] <- findMsMsHRperxcms.direct(fileName[i], cpdID, mode=mode, findPeaksArgs = findPeaksArgs, plots = plots, MSe = MSe)
		}
		
		spectra <- toRMB(unlist(unlist(fspectra, FALSE),FALSE), cpdID, mode)
		
	} else if(length(cpdID) > 1){ # If there are more cpdIDs than files
	
		spectra <- findMsMsHRperxcms.direct(fileName, cpdID, mode=mode, findPeaksArgs = findPeaksArgs, plots = plots, MSe = MSe)
		
		P <- lapply(1:length(spectra), function(i){
			sp <- toRMB(spectra[[i]], cpdID[i], mode)
			sp@id <- as.character(as.integer(cpdID[i]))
			sp@name <- findName(cpdID[i])
			sp@formula <- findFormula(cpdID[i])
			sp@mode <- mode
			sp@polarity <- .polarity[[sp@mode]]
			return(sp)
		})
		return(P)
		
	} else { # There is a file for every cpdID
		spectra <- toRMB(
		  unlist(findMsMsHRperxcms.direct(fileName, cpdID, mode=mode, findPeaksArgs = NULL, plots = FALSE, MSe = FALSE),FALSE)
		  , cpdID, mode)
	}
	
	sp <- spectra
	
	#sp@mz <- mzLimits
	sp@id <- as.character(as.integer(cpdID))
	sp@name <- findName(cpdID)
	sp@formula <- findFormula(cpdID)
	sp@mode <- mode
	return(sp)
}

#' @describeIn findMsMsHRperxcms A submethod of findMsMsHrperxcms that retrieves basic spectrum data
#' @export
findMsMsHRperxcms.direct <- function(fileName, cpdID, mode="pH", findPeaksArgs = NULL, plots = FALSE, MSe = FALSE) {
	
	requireNamespace("CAMERA",quietly=TRUE)
	requireNamespace("xcms",quietly=TRUE)
	
	##
	## getRT function
	##
	
	getRT <- function(xa) {
		rt <- sapply(xa@pspectra, function(x) {median(xcms::peaks(xa@xcmsSet)[x, "rt"])})
	}
	
	##
	## MSMS
	##
	
	# Read file
	suppressWarnings(xrmsms <- xcms::xcmsRaw(fileName, includeMSn=TRUE))
	
	
	# If file is not MSe, split by collision energy
	if(MSe == FALSE){
		# Also, fake MS1 from the MSn data
		suppressWarnings(xrs <- split(xcms::msn2xcmsRaw(xrmsms), f = xrmsms@msnCollisionEnergy))
	} else{
		# Else, MSn data will already be in MS1
		xrs <- list()
		xrs[[1]] <- xrmsms
	}
	
	# Fake a simplistic xcmsSet
	suppressWarnings(setReplicate <- xcms::xcmsSet(files=fileName, method="MS1"))
	xsmsms <- as.list(replicate(length(xrs),setReplicate))
	
	mzabs <- 0.1	

	# Definitions
	whichmissing <- vector()
	metaspec <- list()
	
	##
	## Retrieval over all supplied cpdIDs
	##
	
	for(ID in 1:length(cpdID)){
	
		# Find all relevant information for the current cpdID
		XCMSspectra <- list()
		RT <- findRt(cpdID[ID])$RT * 60
		parentMass <- findMz(cpdID[ID], mode=mode)$mzCenter
		
		# Is the information in the compound list?
		if(is.na(parentMass)){
		  stop(paste("There was no matching entry to the supplied cpdID", cpdID[ID] ,"\n Please check the cpdIDs and the compoundlist."))
		}
		
		# Go over every collision energy of the MS2
		for(i in 1:length(xrs)){
		
			suppressWarnings(capture.output(xcms::peaks(xsmsms[[i]]) <- do.call(xcms::findPeaks,c(findPeaksArgs, object = xrs[[i]]))))
			
			if (nrow(xcms::peaks(xsmsms[[i]])) == 0) {
			  XCMSspectra[[i]] <- matrix(0,2,7)
			  next
			} else{	
			
				# Get the peaklist
				pl <- xcms::peaks(xsmsms[[i]])[,c("mz", "rt"), drop=FALSE]

				# Find precursor peak within limits
				candidates <- which( pl[,"mz", drop=FALSE] < parentMass + mzabs & pl[,"mz", drop=FALSE] > parentMass - mzabs
								& pl[,"rt", drop=FALSE] < RT * 1.1 & pl[,"rt", drop=FALSE] > RT * 0.9 )
				
				# Annotate and group by FWHM (full width at half maximum)
				capture.output(anmsms <- CAMERA::xsAnnotate(xsmsms[[i]]))
				capture.output(anmsms <- CAMERA::groupFWHM(anmsms))
				
				# If a candidate fulfills the condition, choose the closest and retrieve the index of those pesudospectra
				if(length(candidates) > 0){
					closestCandidate <- which.min(abs(RT - pl[candidates, "rt", drop=FALSE]))
					pspIndex <- which(sapply(anmsms@pspectra, function(x) {candidates[[i]][closestCandidate] %in% x}))
				} else{
				# Else choose the candidate with the closest RT
				  if(RMassBank.env$strictMsMsSpectraSelection){
				    pspIndex <- NULL
				  } else {
					  pspIndex <- which.min(abs(getRT(anmsms) - RT))
				  }
				}
				
				# 2nd best: Spectrum closest to MS1
				# pspIndex <- which.min( abs(getRT(anmsms) - actualRT))
				
				# If the plot parameter was supplied, plot it
				if((plots == TRUE) && (length(pspIndex) > 0)){
					CAMERA::plotPsSpectrum(anmsms, pspIndex, log=TRUE,  mzrange=c(0, findMz(cpdID)[[3]]), maxlabel=10)
				}
				
				# If there is a number of indexes, retrieve the pseudospectra
				if(length(pspIndex) != 0){
					XCMSspectra[[i]] <- CAMERA::getpspectra(anmsms, pspIndex)
				} else {
				# Else note the spectrum as missing
					whichmissing <- c(whichmissing,i)
				}
			}
		}
		
		# If XCMSspectra were found but there are some missing for some collision energies, fill these XCMSspectra
		if((length(XCMSspectra) != 0) && length(whichmissing)){
			for(i in whichmissing){
				XCMSspectra[[i]] <- matrix(0,2,7)
			}
		}

		metaspec[[ID]] <- XCMSspectra
	}
	
	return(metaspec)
}

#' Retrieve spectra from msp files
#'
#' This function is currently used to read msp files
#' containing data that were already processed in order to
#' convert the results to MassBank records.
#'
#' @param fileName vector of character-strings
#' The msp files to be searched for spectra
#' @param cpdIDs vector of integers
#' The IDs of compounds in the compoundlist
#' for which spectra should be retrieved
#' @param mode character, default: "pH"
#' The processing mode that was used to produce the spectrum.
#' Should be one of
#' "pH": ([M+H]+)
#' "pNa": ([M+Na]+)
#' "pM": ([M]+)
#' "mH": ([M-H]-)
#' or "mFA": ([M+FA]-)
#' (see the \code{RMassBank} vignette)
#' @return An \code{RmbSpectraSet} with integrated information from the msp files
#' @export
findMsMsHRperMsp <- function(fileName, cpdIDs, mode="pH"){
  # Find mz
  #mzLimits <- findMz(cpdIDs, mode)
  #mz <- mzLimits$mzCenter
  
  # If there are more files than cpdIDs
  if(length(fileName) > 1){
    fspectra <- list()
    
    for(i in 1:length(fileName)){
      fspectra[[i]] <- findMsMsHRperMsp.direct(fileName[i], cpdIDs, mode=mode)
    }
    
    spectra <- toRMB(unlist(unlist(fspectra, FALSE),FALSE), cpdIDs, mode)
    
  } else if(length(cpdIDs) > 1){ # If there are more cpdIDs than files
    
    spectra <- findMsMsHRperMsp.direct(fileName = fileName, cpdIDs = cpdIDs, mode=mode)
    
    P <- lapply(1:length(spectra), function(i){
      sp <- toRMB(msmsXCMSspecs = spectra[[i]], cpdID = cpdIDs[i], mode = mode)
      sp@id <- as.character(as.integer(cpdIDs[i]))
      sp@name <- findName(cpdIDs[i])
      sp@formula <- findFormula(cpdIDs[i])
      sp@mode <- mode
      
      if(length(sp@children) == 1){
        sp@children[[1]]@rawOK <- rep(x = TRUE, times = sp@children[[1]]@peaksCount)
        sp@children[[1]]@good  <- rep(x = TRUE, times = sp@children[[1]]@peaksCount)
        #sp@children[[1]]@good  <- TRUE
      }
      
      return(sp)
    })
    return(P)
    
  } else { # There is a file for every cpdID
    spectra <- toRMB(msmsXCMSspecs = unlist(findMsMsHRperMsp.direct(fileName = fileName, cpdIDs = cpdIDs, mode=mode),FALSE), cpdID = cpdIDs, mode = mode)
  }
  
  sp <- spectra
  
  #sp@mz <- mzLimits
  sp@id <- as.character(as.integer(cpdIDs))
  sp@name <- findName(cpdIDs)
  sp@formula <- findFormula(cpdIDs)
  sp@mode <- mode
  
  return(sp)
}

.retrieve <- function (x, argument) {
	entry <- x[[argument]]
	if(length(entry) == 0 || entry == "NA")
		return(NA)
	else
		return(entry) 
}

#' @describeIn findMsMsHRperMsp A submethod of findMsMsHrperxcms that retrieves basic spectrum data
findMsMsHRperMsp.direct <- function(fileName, cpdIDs, mode="pH") {
  
  #requireNamespace("CAMERA",quietly=TRUE)
  #requireNamespace("xcms",quietly=TRUE)
  
  ##
  ## MSMS
  ##
  
  # Read file
  suppressWarnings(xrmsms <- read.msp(file = fileName))
  xrmsms <- xrmsms[unlist(lapply(X = xrmsms, FUN = function(spectrum){nrow(spectrum$pspectrum)})) > 0]
  
  ## If file is not MSe, split by collision energy
  #if(MSe == FALSE){
  #  # Also, fake MS1 from the MSn data
  #  suppressWarnings(xrs <- split(xcms::msn2xcmsRaw(xrmsms), f = xrmsms@msnCollisionEnergy))
  #} else{
  #  # Else, MSn data will already be in MS1
  #  xrs <- list()
  #  xrs[[1]] <- xrmsms
  #}
  #xrs <- xrmsms
  
  mzabs <- 0.1	
  
  # Definitions
  whichmissing <- vector()
  metaspec <- list()
  
  mzs <- unlist(lapply(X = xrmsms, FUN = function(x){.retrieve(x, 'PRECURSORMZ')}))
  rts <- unlist(lapply(X = xrmsms, FUN = function(x){.retrieve(x, 'RETENTIONTIME')}))
  precursorTable <- data.frame(stringsAsFactors = FALSE,
    mz = as.numeric(mzs),
    rt = as.numeric(rts)
  )
  precursorTable[, "rt"] <- precursorTable[, "rt"] * 60
  
  ##
  ## Retrieval over all supplied cpdIDs
  ##
  
  for(idIdx in seq_along(cpdIDs)){
    
    # Find all relevant information for the current cpdID
    spectrum <- NULL
    RT <- findRt(cpdIDs[[idIdx]])$RT * 60
    parentMass <- findMz(cpdIDs[[idIdx]], mode=mode)$mzCenter
    
    # Is the information in the compound list?
    if(is.na(parentMass)){
      stop(paste("There was no matching entry to the supplied cpdID", cpdIDs[[idIdx]] ,"\n Please check the cpdIDs and the compoundlist."))
    }
    
    # Go over every collision energy of the MS2
    #for(i in seq_along(xrs)){
      
      
      
      if (nrow(precursorTable) == 0) {
        ## no peaks there
        #spectrum <- matrix(0,2,7)
        next
      } else{	
        ## at least one peak there
        
        # Get the peaklist
        #pl <- xrs[[i]]$pspectrum
        #pl <- data.frame("mz" = pl[, "mz"], "rt" = xrs[[i]]$RETENTIONTIME, stringsAsFactors = F)
        
        maximumParentMass <- parentMass + mzabs
        minimumParentMass <- parentMass - mzabs
        maximumRT <- RT * 1.1
        minimumRT <- RT * 0.9
        
        mzMatch <- 
          precursorTable[,"mz", drop=FALSE] < maximumParentMass & 
          precursorTable[,"mz", drop=FALSE] > minimumParentMass
        rtMatch <- 
          precursorTable[,"rt", drop=FALSE] < maximumRT & 
          precursorTable[,"rt", drop=FALSE] > minimumRT
        
        mzMatch[is.na(mzMatch)] <- TRUE ## RT not given
        if(is.na(RT))
          rtMatch <- TRUE
        
        # Find precursor peak within limits
        candidates <- which( mzMatch & rtMatch )
        
        # Annotate and group by FWHM (full width at half maximum)
        #capture.output(anmsms <- CAMERA::xsAnnotate(xsmsms[[i]]))
        #capture.output(anmsms <- CAMERA::groupFWHM(anmsms))
        
        # If a candidate fulfills the condition, choose the closest and retrieve the index of those pesudospectra
        if(length(candidates) > 0){
          if(is.na(RT)){
            pspIndex <- candidates[[1]]
            
            if(RMassBank.env$verbose.output)
              cat(paste("\n### Info ### Compound ", cpdIDs[[idIdx]], ": RT is not given. ", length(candidates), " candidates in range. Taking the first hit: mz[", minimumParentMass, ", ", maximumParentMass , "] vs mz ", precursorTable[pspIndex,"mz"], ".\n", sep = ""))
          } else {
            closestCandidate <- which.min(abs(RT - precursorTable[candidates, "rt"]))
            pspIndex <- candidates[[closestCandidate]]
            if(RMassBank.env$verbose.output)
              cat(paste("\n### Info ### Compound ", cpdIDs[[idIdx]], ": ", length(candidates), " candidates in range. Taking the closest hit regarding RT (", RT, "): mz[", minimumParentMass, ", ", maximumParentMass , "] x rt[", minimumRT, ", ", maximumRT, "] vs (mz ", precursorTable[pspIndex,"mz"], ", rt ", precursorTable[pspIndex,"rt"], ")\n", sep = ""))
          }
        } else{
          # Else choose the candidate with the closest RT
          if(RMassBank.env$strictMsMsSpectraSelection){
            pspIndex <- NULL
            cat(paste("\n### Warning ### Compound ", cpdIDs[[idIdx]], ": No candidates in range.\n", sep = ""))
          } else {
            pspIndex <- which.min(abs(RT - precursorTable[, "rt"]))
            cat(paste("\n### Warning ### Compound ", cpdIDs[[idIdx]], ": No candidates in range. Taking the closest hit regarding RT (", RT, "): mz[", minimumParentMass, ", ", maximumParentMass , "] x rt[", minimumRT, ", ", maximumRT, "] vs (mz ", precursorTable[pspIndex,"mz"], ", rt ", precursorTable[pspIndex,"rt"], ")\n", sep = ""))
          }
        }
        
        # 2nd best: Spectrum closest to MS1
        # pspIndex <- which.min( abs(getRT(anmsms) - actualRT))
        
        ## If the plot parameter was supplied, plot it
        #if((plots == TRUE) && (length(pspIndex) > 0)){
        #  CAMERA::plotPsSpectrum(anmsms, pspIndex, log=TRUE,  mzrange=c(0, findMz(cpdIDs)[[3]]), maxlabel=10)
        #}
        
        # If there is a number of indexes, retrieve the pseudospectra
        if(length(pspIndex) != 0){
          spectrum <- xrmsms[[pspIndex]]
        } else {
          # Else note the spectrum as missing
          whichmissing <- c(whichmissing,idIdx)
          #spectrum <- matrix(0,2,7)
        }
      }
    #}
    
    # If XCMSspectra were found but there are some missing for some collision energies, fill these XCMSspectra
    #if((length(XCMSspectra) != 0) && length(whichmissing)){
    #  for(i in whichmissing){
    #    XCMSspectra[[idIdx]] <- matrix(0,2,7)
    #  }
    #}
    
    if(is.null(spectrum)){
      metaspec[[idIdx]] <- list(matrix(0,1,7))
    } else {
      mz <- as.numeric(spectrum$pspectrum[, "mz"])
      rt <- as.numeric(.retrieve(spectrum, 'RETENTIONTIME'))
      metaspec[[idIdx]] <- list(data.frame(
        stringsAsFactors = F,
        "mz"      = mz,
        "mzmin"   = mz,
        "mzmax"   = mz,
        "rt"      = rt,
        "rtmin"   = rt,
        "rtmax"   = rt,
        "into"    = as.numeric(spectrum$pspectrum[, "intensity"]),
        "into_parent" = as.numeric(spectrum$INTENSITY)
      ))
    }
  }
  
  return(metaspec)
}

## adapted from the Bioconductor package 'metaMS' (method 'read.msp')
read.msp <- function(file){
  get.text.value <- function(x, field, do.err = TRUE) {
    if(trimws(x) == field) return("")
    woppa <- strsplit(x, field)
    woppa.lengths <- sapply(woppa, length)
    if (all(woppa.lengths == 2)) {
      sapply(woppa, function(y) gsub("^ +", "", y[2]))
    }
    else {
      if (do.err) {
        stop(paste("Invalid field", field, "in", x[woppa.lengths != 2]))
      }
      else {
        NULL
      }
    }
  }
  read.compound <- function(strs) {
    fields.idx <- grep(":", strs)
    fields <- sapply(strsplit(strs[fields.idx], ":"), "[[", 1)
    pk.idx <- which(fields == "Num Peaks")
    if (length(pk.idx) == 0) 
      stop("No spectrum found")
    cmpnd <- lapply(fields.idx[-pk.idx], function(x) get.text.value(strs[x], paste(fields[x], ":", sep = "")))
    names(cmpnd) <- fields[-pk.idx]
    if(!("INTENSITY" %in% names(cmpnd))) cmpnd$"INTENSITY" <- 100
    
    cmpnd$PRECURSORMZ   <- gsub(x = cmpnd$PRECURSORMZ,   pattern = ",", replacement = ".")
    cmpnd$RETENTIONTIME <- gsub(x = cmpnd$RETENTIONTIME, pattern = ",", replacement = ".")
    cmpnd$INTENSITY     <- gsub(x = cmpnd$INTENSITY,     pattern = ",", replacement = ".")
    
    ## minutes to seconds
    #cmpnd$RETENTIONTIME <- as.numeric(cmpnd$RETENTIONTIME) * 60
    
    nlines <- length(strs)
    npeaks <- as.numeric(get.text.value(strs[pk.idx], "Num Peaks:"))
    peaks.idx <- (pk.idx + 1):nlines
    pks <- gsub("^ +", "", unlist(strsplit(strs[peaks.idx], ";")))
    pks <- pks[pks != ""]
    if (length(pks) != npeaks) 
      stop(paste("Not the right number of peaks in compound '", cmpnd$Name, "' (", npeaks, " vs ", length(pks), ") in file '", file, "'", sep = ""))
    pklst <- strsplit(x = pks, split = "\t| ")
    pklst <- lapply(pklst, function(x) x[x != ""])
    cmz <- as.numeric(gsub(x = sapply(pklst, "[[", 1), pattern = ",", replacement = "."))
    cintens <- as.numeric(gsub(x = sapply(pklst, "[[", 2), pattern = ",", replacement = "."))
    finaltab <- matrix(c(cmz, cintens), ncol = 2)
    if (any(table(cmz) > 1)) {
      warning("Duplicate mass in compound ", cmpnd$Name, " (CAS ", cmpnd$CAS, ")... summing up intensities")
      finaltab <- aggregate(finaltab[, 2], by = list(finaltab[, 1]), FUN = sum)
    }
    colnames(finaltab) <- c("mz", "intensity")
    c(cmpnd, list(pspectrum = finaltab))
  }
  huhn <- readLines(con = file)
  starts <- which(regexpr("(Name:)|(NAME:) ", huhn) == 1)
  ends <- c(starts[-1] - 1, length(huhn))
  lapply(1:length(starts), function(i){
    read.compound(huhn[starts[[i]]:ends[[i]]])
  })
}
## new
################################################################################

# Finds the EIC for a mass trace with a window of x ppm.
# (For ppm = 10, this is +5 / -5 ppm from the non-recalibrated mz.)
#' Extract EICs 
#' 
#' Extract EICs from raw data for a determined mass window.
#' 
#' @param msRaw The mzR file handle 
#' @param mz The mass or mass range to extract the EIC for: either a single mass
#' 			(with the range specified by \code{limit} below) or a mass range
#' 			in the form of \code{c(min, max)}. 
#' @param limit If a single mass was given for \code{mz}: the mass window to extract.
#' 			A limit of 0.001 means that the EIC will be returned for \code{[mz - 0.001, mz + 0.001]}.
#' @param rtLimit If given, the retention time limits in form \code{c(rtmin, rtmax)} in seconds.
#' @param headerCache If present, the complete \code{mzR::header(msRaw)}. Passing
#' 			this value is useful if spectra for multiple compounds should be 
#' 			extracted from the same mzML file, since it avoids getting the data
#' 			freshly from \code{msRaw} for every compound.
#' @param peaksCache If present, the complete output of \code{mzR::peaks(msRaw)}. This speeds up the lookup
#' 			if multiple compounds should be searched in the same file.
#' @param floatingRecalibration 
#' 			A fitting function that \code{predict()}s a mass shift based on the retention time. Can be used
#' 			if a lockmass calibration is known (however you have to build the calibration yourself.)
#' @return A \code{[rt, intensity, scan]} matrix (\code{scan} being the scan number.) 
#' @author Michael A. Stravs, Eawag <michael.stravs@@eawag.ch>
#' @seealso findMsMsHR
#' @export
findEIC <- function(msRaw, mz, limit = NULL, rtLimit = NA, headerCache = NULL, floatingRecalibration = NULL,
		peaksCache = NULL, polarity = NA, msLevel = 1, precursor = NULL)
{
	# calculate mz upper and lower limits for "integration"
	if(all(c("mzMin", "mzMax") %in% names(mz)))
		mzlimits <- c(mz$mzMin, mz$mzMax)
	else
		mzlimits <- c(mz - limit, mz + limit)
	# Find peaklists for all MS1 scans
	if(!is.null(headerCache))
		headerData <- as.data.frame(headerCache)
	else
		headerData <- as.data.frame(header(msRaw))
	# Add row numbering because I'm not sure if seqNum or acquisitionNum correspond to anything really
	if(nrow(headerData) > 0)
		headerData$rowNum <- 1:nrow(headerData)
	else
		headerData$rowNum <- integer(0)
	
	# If RT limit is already given, retrieve only candidates in the first place,
	# since this makes everything much faster.
	if(all(!is.na(rtLimit)))
		headerMS1 <- headerData[
				(headerData$msLevel == msLevel) & (headerData$retentionTime >= rtLimit[[1]])
						& (headerData$retentionTime <= rtLimit[[2]])
				,]
	else
		headerMS1 <- headerData[headerData$msLevel == msLevel,]
	# DIA handling:
	if(!is.null(precursor))
	{
	  headerMS1 <- headerMS1[abs(headerMS1$precursorMZ - precursor) < 0.1, ]
	}
	if(!is.na(polarity))
	{
	  if(is.character(polarity))
	    polarity <- .polarity[[polarity]]
	  headerMS1 <- headerMS1[headerMS1$polarity == polarity,]
	}
	
	
	if(is.null(peaksCache))
		pks <- mzR::peaks(msRaw, headerMS1$seqNum)
	else
		pks <- peaksCache[headerMS1$rowNum]
		
	# Sum intensities in the given mass window for each scan
	if(is.null(floatingRecalibration))
	{
		headerMS1$mzMin <- mzlimits[[1]]
		headerMS1$mzMax <- mzlimits[[2]]
	}
	else
	{
		headerMS1$mzMin <- mzlimits[[1]] + predict(floatingRecalibration, headerMS1$retentionTime)
		headerMS1$mzMax <- mzlimits[[2]] + predict(floatingRecalibration, headerMS1$retentionTime)
	}
	intensity <- unlist(lapply(1:nrow(headerMS1), function(row){
						peaktable <- pks[[row]]
						sum(peaktable[
							which((peaktable[,1] >= headerMS1[row,"mzMin"]) & (peaktable[,1] <= headerMS1[row,"mzMax"])),2
						])
						
	}))
	return(data.frame(rt = headerMS1$retentionTime, intensity=intensity, scan=headerMS1$acquisitionNum))
}


#' Generate peaks cache
#' 
#' Generates a peak cache table for use with \code{\link{findMsMsHR}} functions.
#' 
#' @param msRaw the input raw datafile (opened)
#' @param headerCache the cached header, or subset thereof for which peaks should be extracted. Peak extraction goes
#' 		by \code{seqNum}.
#' @return A list of dataframes as from \code{mzR::peaks}.
#' 
#' @author stravsmi
#' @export
makePeaksCache <- function(msRaw, headerCache) 
{
	mzR::peaks(msRaw, headerCache$seqNum)
}

#' Conversion of XCMS-pseudospectra into RMassBank-spectra
#' 
#' Converts a pseudospectrum extracted from XCMS using CAMERA into the msmsWorkspace(at)spectrum-format that RMassBank uses
#'
#' @usage toRMB(msmsXCMSspecs, cpdID, mode, MS1spec)
#' @param msmsXCMSspecs The compoundID of the compound that has been used for the peaklist
#' @param cpdID The compound ID of the substance of the given spectrum
#' @param mode The ionization mode that has been used for the spectrum
#' @param MS1spec The MS1-spectrum from XCMS, which can be optionally supplied
#' @return One list element of the (at)specs-entry from an msmsWorkspace
#' @seealso \code{\link{msmsWorkspace-class}}
#' @author Erik Mueller
#' @examples \dontrun{
#' 		XCMSpspectra <- findmsmsHRperxcms.direct("Glucolesquerellin_2184_1.mzdata", 2184)
#'      wspecs <- toRMB(XCMSpspectra)
#' }
#' @export
toRMB <- function(msmsXCMSspecs = NA, cpdID = NA, mode="pH", MS1spec = NA){

	
	##Basic parameters
	mz <- findMz(cpdID,mode=mode)$mzCenter
	id <- cpdID
	formula <- findFormula(cpdID)
	
	if(length(msmsXCMSspecs) == 0){
		return(new("RmbSpectraSet",found=FALSE))
	}
	
	foundOK <- !any(sapply(msmsXCMSspecs, function(x) all(x == 0)))
	
	
	if(!foundOK){
		return(new("RmbSpectraSet",found=FALSE))
	}
	
	if(suppressWarnings(is.na(msmsXCMSspecs)[1])){
			stop("You need a readable spectrum!")
	}
	
	if(is.na(cpdID)){
			stop("Please supply the compoundID!")
	}
	
	mockAcqnum <- 1
	mockenv <- environment()
	
	msmsSpecs <- lapply(msmsXCMSspecs, function(spec){
		## Mock acquisition num
		mockenv$mockAcqnum <- mockenv$mockAcqnum + 1
		
		## Find peak table
		pks <- matrix(nrow = nrow(spec), ncol = 2)
		colnames(pks) <- c("mz","int")
		pks[,1] <- spec[,"mz"]
		pks[,2] <- spec[,"into"]
		
		## Deprofiling not necessary for XCMS
		
		## New spectrum object
		return(new("RmbSpectrum2",
				mz = pks[,"mz"],
				intensity = pks[,"int"],
				precScanNum = as.integer(1),
				precursorMz = findMz(cpdID, mode=mode)$mzCenter,
				precursorIntensity = ifelse(test = "into_parent" %in% colnames(spec), yes = spec[,"into_parent"], no = 0),
				precursorCharge = as.integer(1),
				collisionEnergy = 0,
				polarity = .polarity[[mode]],
				tic = 0,
				peaksCount = nrow(spec),
				rt = median(spec[,"rt"]),
				acquisitionNum = as.integer(mockenv$mockAcqnum),
				centroided = TRUE
		))
	})
	
	msmsSpecs <- as(do.call(c, msmsSpecs), "SimpleList")			
	
	##Build the new objects
	masterSpec <- new("Spectrum1",
				mz = findMz(cpdID,mode=mode)$mzCenter,
				intensity = ifelse(test = msmsSpecs[[1]]@precursorIntensity != 0, yes = msmsSpecs[[1]]@precursorIntensity, no = 100),
				polarity = as.integer(0),
				peaksCount = as.integer(1),
				rt = msmsSpecs[[1]]@rt,
				acquisitionNum = as.integer(1),
				tic = 0,
				centroided = TRUE
			)
						
	spectraSet <- new("RmbSpectraSet",
				parent = masterSpec,
				children = msmsSpecs,
				found = TRUE,
				#complete = NA,
				#empty = NA,
				#formula = character(),
				mz = mz
				#name = character(),
				#annotations = list()
			)
	
	return(spectraSet)
}

#' Addition of manual peaklists
#' 
#' Adds a manual peaklist in matrix-format
#'
#' @usage addPeaksManually(w, cpdID, handSpec, mode)
#' @param w The msmsWorkspace that the peaklist should be added to.
#' @param cpdID The compoundID of the compound that has been used for the peaklist
#' @param handSpec A peaklist with 2 columns, one with "mz", one with "int" 
#' @param mode The ionization mode that has been used for the spectrum represented by the peaklist
#' @return The \code{msmsWorkspace} with the additional peaklist added to the right spectrum
#' @seealso \code{\link{msmsWorkflow}}
#' @author Erik Mueller
#' @examples \dontrun{
#' 		handSpec <- cbind(mz=c(274.986685367956, 259.012401087427, 95.9493025990907, 96.9573002472772),
#'                                int=c(357,761, 2821, 3446))
#' 		addPeaksManually(w, cpdID, handSpec)
#' }
#' @export
addPeaksManually <- function(w, cpdID = NA, handSpec, mode = "pH"){

			
	if(is.na(cpdID)){
			stop("Please supply the compoundID!")
	}
	
	# For the case that the cpdID turns up for the first time
	# a new spectrumset needs to be created
	if(!(cpdID %in% sapply(w@spectra,function(s) s@id))){
		
		# Create fake MS1 spectrum
		masterSpec <- new("Spectrum1",
			mz = findMz(cpdID,mode=mode)$mzCenter,
			intensity = 100,
			polarity = as.integer(0),
			peaksCount = as.integer(1),
			rt = findRt(cpdID)$RT,
			acquisitionNum = as.integer(1),
			tic = 0,
			centroided = TRUE
		)
		
		# Create fake spectrumset
		spectraSet <- new("RmbSpectraSet",
			parent = masterSpec,
			found = TRUE,
			#complete = NA,
			#empty = NA,
			id = as.character(as.integer(cpdID)),
			formula = findFormula(cpdID),
			mz = findMz(cpdID,mode=mode)$mzCenter,
			name = findName(cpdID),
			mode = mode
			#annotations = list()
		)
		
		w@spectra[[length(w@spectra) + 1]] <- spectraSet
	}
	
	specIndex <- which(cpdID == sapply(w@spectra, function(s) s@id))
	
	# New spectrum object
	w@spectra[[specIndex]]@children[[length(w@spectra[[specIndex]]@children) + 1]] <- new("RmbSpectrum2",
				mz = handSpec[,"mz"],
				intensity = handSpec[,"int"],
				precScanNum = as.integer(1),
				precursorMz = findMz(cpdID,mode=mode)$mzCenter,
				precursorIntensity = 0,
				precursorCharge = as.integer(1),
				collisionEnergy = 0,
				tic = 0,
				peaksCount = nrow(handSpec),
				rt = findRt(cpdID)$RT,
				acquisitionNum = as.integer(length(w@spectra[[specIndex]]@children) + 2),
				centroided = TRUE)
	return(w)
}


createSpecsFromPeaklists <- function(w, cpdIDs, filenames, mode="pH"){
	for(j in 1:length(filenames)){
		w <- addPeaksManually(w,cpdIDs[j],as.matrix(read.csv(filenames[j]), header=TRUE),mode)
	}
	
	return(w)
}


#' MassBank-record Addition
#' 
#' Adds the peaklist of a MassBank-Record to the specs of an msmsWorkspace
#'
#' @aliases addMB
#' @usage addMB(w, cpdID, fileName, mode)
#' @param w The msmsWorkspace that the peaklist should be added to.
#' @param cpdID The compoundID of the compound that has been used for the record
#' @param fileName The path to the record
#' @param mode The ionization mode that has been used to create the record
#' @return The \code{msmsWorkspace} with the additional peaklist from the record
#' @seealso \code{\link{addPeaksManually}}
#' @author Erik Mueller
#' @examples \dontrun{
#' 		addMB("filepath_to_records/RC00001.txt")
#' }
#' @export
addMB <- function(w, cpdID, fileName, mode){
	mb <- parseMassBank(fileName)
	peaklist <- list()
	peaklist[[1]] <- mb@compiled_ok[[1]][["PK$PEAK"]][,1:2]
	w <- addPeaksManually(w, cpdID, peaklist[[1]], mode)
	return(w)
}


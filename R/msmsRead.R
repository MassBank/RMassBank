#' @import R.utils
NULL

#' 
#' Extracts and processes spectra from a specified file list, according to 
#' loaded options and given parameters.
#' 
#' The filenames of the raw LC-MS runs are read from the array \code{files} 
#' in the global enviroment.
#' See the vignette \code{vignette("RMassBank")} for further details about the
#' workflow.
#' 
#' @param w A \code{msmsWorkspace} to work with.
#' @param filetable The path to a .csv-file that contains the columns "Files" and "ID" supplying
#' 			the relationships between files and compound IDs. Either this or the parameter "files" need
#'			to be specified. 
#' @param files A vector or list containing the filenames of the files that are to be read as spectra. 
#'				For the IDs to be inferred from the filenames alone, there need to be exactly 2 underscores.
#' @param cpdids A vector or list containing the compound IDs of the files that are to be read as spectra.
#'				The ordering of this and \code{files} implicitly assigns each ID to the corresponding file.
#'				If this is supplied, then the IDs implicitly named in the filenames are ignored.
#' @param readMethod Several methods are available to get peak lists from the files.
#'        Currently supported are "mzR", "xcms", "MassBank" and "peaklist".
#'        The first two read MS/MS raw data, and differ in the strategy 
#'        used to extract peaks. MassBank will read existing records, 
#'        so that e.g. a recalibration can be performed, and "peaklist" 
#'        just requires a CSV with two columns and the column header "mz", "int".
#' @param mode \code{"pH", "pNa", "pM", "pNH4", "mH", "mM", "mFA"} for different ions 
#' 			([M+H]+, [M+Na]+, [M]+, [M+NH4]+, [M-H]-, [M]-, [M+FA]-).
#' 			For `readMethod == "mzR"`, a vector of `mode` entries is supported. The user 
#' 			should check that they are either all positive or negative. If this isn't the case,
#' 			the recalibration will be incorrect.
#' @param confirmMode Defaults to false (use most intense precursor). Value 1 uses
#' 			the 2nd-most intense precursor for a chosen ion (and its data-dependent scans)
#' 			, etc.
#' @param useRtLimit Whether to enforce the given retention time window.
#' @param Args A list of arguments that will be handed to the xcms-method findPeaks via do.call
#' @param settings Options to be used for processing. Defaults to the options loaded via
#' 			\code{\link{loadRmbSettings}} et al. Refer to there for specific settings.
#' @param progressbar The progress bar callback to use. Only needed for specialized applications.
#' 			Cf. the documentation of \code{\link{progressBarHook}} for usage.
#' @param MSe A boolean value that determines whether the spectra were recorded using MSe or not
#' @param plots A boolean value that determines whether the pseudospectra in XCMS should be plotted
#' @return The \code{msmsWorkspace} with msms-spectra read.
#' @seealso \code{\link{msmsWorkspace-class}}, \code{\link{msmsWorkflow}}
#' @author Michael Stravs, Eawag <michael.stravs@@eawag.ch>
#' @author Erik Mueller, UFZ
#' @export
msmsRead <- function(w, filetable = NULL, files = NULL, cpdids = NULL, 
					readMethod, mode = NULL, confirmMode = FALSE, useRtLimit = TRUE, 
					Args = NULL, settings = getOption("RMassBank"),
                    progressbar = "progressBarHook", MSe = FALSE, plots = FALSE){
	.checkMbSettings()

  
	if(is.null(filetable)){
		##If no filetable is supplied, filenames must be named explicitly
		if(is.null(files))
			stop("Please supply the files")
	  if(is.null(mode))
	    stop("Please supply the mode(s)")
		
		##Assign the filenames to the workspace
		w@files <- unlist(files)
		
		##If no filetable is supplied, cpdids must be delivered explicitly or implicitly within the filenames
		if(is.null(cpdids)){
			splitfn <- strsplit(files,"_")
			splitsfn <- sapply(splitfn, function(x) x[length(x)-1])
			if(suppressWarnings(any(is.na(as.numeric(splitsfn)[1]))))
				stop("Please supply the cpdids corresponding to the files in the filetable or the filenames")
			cpdids <- splitsfn
		}
	} else{
		##If a filetable is supplied read it
		tab <- read.csv(filetable, stringsAsFactors = FALSE)
		# Check if we have absolute or relative paths.
		# If relative, they are assumed to be relative to the filetable path
		
		tab[,"Files"] <- ifelse(
		  isAbsolutePath(tab[,"Files"]),
		  tab[,"Files"],
		  paste(dirname(filetable), tab[,"Files"], sep="/")
		)
		w@files <- tab[,"Files"]
		cpdids <- tab[,"ID"]
            if ("mode" %in% colnames(tab)) {
                mode <- tab[,"mode"]
            }
	}
  
  ##Read the files and cpdids according to the definition
  ##All cases are silently accepted, as long as they can be handled according to one definition
  if(!all(mode %in% knownAdducts())) stop(paste("The ionization mode", mode, "is unknown."))
	
	##If there's more cpdids than filenames or the other way around, then abort
	if(length(w@files) != length(cpdids)){
		stop("There are a different number of cpdids than files")
	}
	
	if(!(readMethod %in% c("mzR","peaklist","xcms","minimal","msp"))){
		stop("The supplied method does not exist")
	}
	
	if(!all(file.exists(w@files))){
		stop("The supplied files ", paste(w@files[!file.exists(w@files)]), " don't exist. Paths in the Filelist were interpreted relative to the location of the Filelist.")
	}

    # na.ids <- which(is.na(sapply(cpdids, findSmiles)))
    
    # if(length(na.ids)){
        # stop("The supplied compound ids ", paste(cpdids[na.ids], collapse=" "), " don't have a corresponding smiles entry. Maybe they are missing from the compound list")
    # }

	##This should work
  if(readMethod == "minimal"){
      ##Edit options
      opt <- getOption("RMassBank")
      opt$recalibrator$MS1 <- "recalibrate.identity"
      opt$recalibrator$MS2 <- "recalibrate.identity"
      opt$add_annotation==FALSE
      options(RMassBank=opt)
      ##Edit analyzemethod
      analyzeMethod <- "intensity"
  }
	
	if(readMethod == "mzR"){
	  
	  # To do: check if we can use this verbatim in xcms method too
	  mode_ <- mode
	  if(length(mode) == 1)
	    mode_ <- rep(mode, length(w@files))
	  if(length(mode_) != length(w@files))
	    stop("Supply either one mode or a vector for one mode per file")
	  
		##Progressbar
		nLen <- length(w@files)
		nProg <- 0
		pb <- do.call(progressbar, list(object=NULL, value=0, min=0, max=nLen))
		
		w@spectra <- as(lapply(seq_along(w@files), function(i) {
							
		          fileName <- w@files[i]
							# Find compound ID
							cpdID <- cpdids[i]


							# Retrieve spectrum data
							spec <- findMsMsHR(fileName = fileName, 
									cpdID = cpdID, mode = mode_[i], confirmMode = confirmMode, useRtLimit = useRtLimit,
									ppmFine = settings$findMsMsRawSettings$ppmFine,
									mzCoarse = settings$findMsMsRawSettings$mzCoarse,
									fillPrecursorScan = settings$findMsMsRawSettings$fillPrecursorScan,
									rtMargin = settings$rtMargin,
									deprofile = settings$deprofile)
							gc()
														
							# Progress:
							nProg <<- nProg + 1
							pb <- do.call(progressbar, list(object=pb, value= nProg))
							
							return(spec)
						} ), "SimpleList")
		names(w@spectra) <- basename(as.character(w@files))
	}
	
	##xcms-readmethod 
	if(readMethod == "xcms"){
		
		##Load libraries
		requireNamespace("xcms",quietly=TRUE)
		requireNamespace("CAMERA",quietly=TRUE)
		
		##Find unique files and cpdIDs
		ufiles <- unique(w@files)
		uIDs <- unique(cpdids)
		nLen <- length(ufiles)
		
		##Progressbar
		nProg <- 0
		pb <- do.call(progressbar, list(object=NULL, value=0, min=0, max=nLen))
		i <- 1
		
		##Routine for the case of multiple cpdIDs per file
		if(length(uIDs) > length(ufiles)){
			w@spectra <- as(unlist(lapply(ufiles, function(currentFile){
						fileIDs <- cpdids[which(w@files == currentFile)]
						spec <- findMsMsHRperxcms(currentFile, fileIDs, mode=mode, findPeaksArgs=Args, plots, MSe = MSe)
						gc()
					
						# Progress:
						nProg <<- nProg + 1
						pb <- do.call(progressbar, list(object=pb, value= nProg))
						
						return(spec)
					}),FALSE),"SimpleList")
		} else {
		  ##Routine for the other cases
		  w@spectra <- as(lapply(uIDs, function(ID){
		    # Find files corresponding to the compoundID
		    currentFile <- w@files[which(cpdids == ID)]
		    
		    # Retrieve spectrum data
		    spec <- findMsMsHRperxcms(currentFile, ID, mode=mode, findPeaksArgs=Args, plots, MSe = MSe)
		    gc()
		    
		    # Progress:
		    nProg <<- nProg + 1
		    pb <- do.call(progressbar, list(object=pb, value= nProg))
		    
		    return(spec)
		  }),"SimpleList")
		  ##If there are more files than unique cpdIDs, only remember the first file for every cpdID
		  w@files <- w@files[sapply(uIDs, function(ID){
		    return(which(cpdids == ID)[1])
		  })]
		}
	}
	
	##Peaklist-readmethod 
	if((readMethod == "peaklist") || (readMethod=="minimal")){
		w <- createSpecsFromPeaklists(w, cpdids, filenames=w@files, mode=mode)
		uIDs <- unique(cpdids)
		files <- list()
		
		for(i in 1:length(uIDs)){
			indices <- sapply(cpdids,function(a){return(uIDs[i] %in% a)})
			files[[i]] <- w@files[indices]
		}
		
		w@files <- sapply(files,function(file){return(file[1])})
		rmb_log_info("Peaks read")
	}
  
  ##MSP-readmethod 
  if(readMethod == "msp"){
    ##Find unique files and cpdIDs
    ufiles <- unique(w@files)
    uIDs <- unique(cpdids)
    nLen <- length(ufiles)
    
    ##Progressbar
    nProg <- 0
    pb <- do.call(progressbar, list(object=NULL, value=0, min=0, max=nLen))
    i <- 1
    
    ##Routine for the case of multiple cpdIDs per file
    if(length(uIDs) > length(ufiles)){
      w@spectra <- as(unlist(lapply(ufiles, function(currentFile){
        fileIDs <- cpdids[which(w@files == currentFile)]
        spec <- findMsMsHRperMsp(fileName = currentFile, cpdIDs = fileIDs, mode=mode)
        gc()
        
        # Progress:
        nProg <<- nProg + 1
        pb <- do.call(progressbar, list(object=pb, value= nProg))
        
        return(spec)
      }),FALSE),"SimpleList")
    } else {
      ##Routine for the other cases
      w@spectra <- as(lapply(uIDs, function(ID){
        # Find files corresponding to the compoundID
        currentFile <- w@files[which(cpdids == ID)]
        
        # Retrieve spectrum data
        spec <- findMsMsHRperMsp(fileName = currentFile, cpdIDs = ID, mode=mode)
        gc()
        
        # Progress:
        nProg <<- nProg + 1
        pb <- do.call(progressbar, list(object=pb, value= nProg))
        
        return(spec)
      }),"SimpleList")
      ##If there are more files than unique cpdIDs, only remember the first file for every cpdID
      w@files <- w@files[sapply(uIDs, function(ID){
        return(which(cpdids == ID)[1])
      })]
    }
  }
	
  ## verbose output
  if(RMassBank.env$verbose.output)
    for(parentIdx in seq_along(w@spectra))
      if(!w@spectra[[parentIdx]]@found)
        cat(paste("### Warning ### No precursor ion was detected for ID '", w@spectra[[parentIdx]]@id, "'\n", sep = ""))
  
  return(w)
}

#' 
#' Extracts and processes spectra from a list of xcms-Objects
#' 
#' The filenames of the raw LC-MS runs are read from the array \code{files} 
#' in the global enviroment.
#' See the vignette \code{vignette("RMassBank")} for further details about the
#' workflow.
#' 
#' @param w A \code{msmsWorkspace} to work with.
#' @param xRAW A list of xcmsRaw objects whose peaks should be detected and added to the workspace.
#'				The relevant data must be in the MS1 data of the xcmsRaw object.  You can coerce the
#'				msn-data in a usable object with the \code{msn2xcmsRaw} function of xcms.
#' @param cpdids A vector or list containing the compound IDs of the files that are to be read as spectra.
#'				The ordering of this and \code{files} implicitly assigns each ID to the corresponding file.
#'				If this is supplied, then the IDs implicitly named in the filenames are ignored.
#' @param mode \code{"pH", "pNa", "pM", "pNH4", "mH", "mM", "mFA"} for different ions 
#' 			([M+H]+, [M+Na]+, [M]+, [M+NH4]+, [M-H]-, [M]-, [M+FA]-).
#' @param findPeaksArgs A list of arguments that will be handed to the xcms-method findPeaks via do.call
#' @param settings Options to be used for processing. Defaults to the options loaded via
#' 			\code{\link{loadRmbSettings}} et al. Refer to there for specific settings.
#' @param progressbar The progress bar callback to use. Only needed for specialized applications.
#' 			Cf. the documentation of \code{\link{progressBarHook}} for usage.
#' @param plots A boolean value that determines whether the pseudospectra in XCMS should be plotted
#' @return The \code{msmsWorkspace} with msms-spectra read.
#' @seealso \code{\link{msmsWorkspace-class}}, \code{\link{msmsWorkflow}}
#' @author Michael Stravs, Eawag <michael.stravs@@eawag.ch>
#' @author Erik Mueller, UFZ
#' @export
msmsRead.RAW <- function(w, xRAW = NULL, cpdids = NULL, mode, findPeaksArgs = NULL, 
							settings = getOption("RMassBank"), progressbar = "progressBarHook", plots = FALSE){
	
	requireNamespace("xcms", quietly=TRUE)
	
	##xRAW will be coerced into a list of length 1 if it is an xcmsRaw-object
	if(class(xRAW) == "xcmsRaw"){
		xRAW <- list(xRAW)
	}
	
	##Error messages
	if((class(xRAW) != "list") || any(sapply(xRAW, function(x) class(x) != "xcmsRaw"))){
		stop("No list of xcmsRaw-objects supplied")
	}
	
	if(is.null(cpdids)){
		stop("No cpdids supplied")
	}
		
	#msnExist <- which(sapply(xRAW,function(x) length(x@msnPrecursorScan) != 0))
	#print(length(msnExist))
	#print(length(xRAW))
	
	#if(length(msnExist) != length(xRAW)){
	#	stop(paste("No msn data in list elements", setdiff(1:length(xRAW),msnExist)))
	#}
	
	requireNamespace("CAMERA",quietly=TRUE)
	
	parentMass <- findMz(cpdids[1], mode=mode)$mzCenter
	if(is.na(parentMass)){
		stop(paste("There was no matching entry to the supplied cpdID", cpdids[1] ,"\n Please check the cpdIDs and the compoundlist."))
	}
		
	RT <- findRt(cpdids[1])$RT * 60
	mzabs <- 0.1
	
	getRT <- function(xa) {
		rt <- sapply(xa@pspectra, function(x) {median(peaks(xa@xcmsSet)[x, "rt"])})
	}
	
	suppressWarnings(setReplicate <- xcms::xcmsSet(files=xRAW[[1]]@filepath, method="MS1"))
	xsmsms <- as.list(replicate(length(xRAW),setReplicate))
	candidates <- list()
	anmsms <- list()
	psp <- list()
	spectra <- list()
	whichmissing <- vector()
	metaspec <- list()
	for(i in 1:length(xRAW)){
		devnull <- suppressWarnings(capture.output(xcms::peaks(xsmsms[[i]]) <- do.call(xcms::findPeaks,c(findPeaksArgs, object = xRAW[[i]]))))
		
		if (nrow(xcms::peaks(xsmsms[[i]])) == 0) { ##If there are no peaks
			spectra[[i]] <- matrix(0,2,7)
			next
		} else{	
			## Get pspec 
			pl <- xcms::peaks(xsmsms[[i]])[,c("mz", "rt"), drop=FALSE]

			## Best: find precursor peak
			candidates[[i]] <- which( pl[,"mz", drop=FALSE] < parentMass + mzabs & pl[,"mz", drop=FALSE] > parentMass - mzabs
							& pl[,"rt", drop=FALSE] < RT * 1.1 & pl[,"rt", drop=FALSE] > RT * 0.9 )
			devnull <- capture.output(anmsms[[i]] <- CAMERA::xsAnnotate(xsmsms[[i]]))
			devnull <- capture.output(anmsms[[i]] <- CAMERA::groupFWHM(anmsms[[i]]))

			if(length(candidates[[i]]) > 0){
				closestCandidate <- which.min (abs( RT - pl[candidates[[i]], "rt", drop=FALSE]))
				psp[[i]] <- which(sapply(anmsms[[i]]@pspectra, function(x) {candidates[[i]][closestCandidate] %in% x}))
			} else{
				psp[[i]] <- which.min( abs(getRT(anmsms[[i]]) - RT) )
			}
			## Now find the pspec for compound

			## 2nd best: Spectrum closest to MS1
			##psp <- which.min( abs(getRT(anmsms) - actualRT))

			## 3rd Best: find pspec closest to RT from spreadsheet
			##psp <- which.min( abs(getRT(anmsms) - RT) )
			if((plots == TRUE) && (length(psp[[i]]) > 0)){
				CAMERA::plotPsSpectrum(anmsms[[i]], psp[[i]], log=TRUE,  mzrange=c(0, findMz(cpdids[1])[[3]]), maxlabel=10)
			}
			if(length(psp[[i]]) != 0){
				spectra[[i]] <- CAMERA::getpspectra(anmsms[[i]], psp[[i]])
			} else {
				whichmissing <- c(whichmissing,i)
			}
		}
	}
	if(length(spectra) != 0){
		for(i in whichmissing){
			spectra[[i]] <- matrix(0,2,7)
		}
	}
	
	sp <- toRMB(spectra,cpdids,"mH")
	sp@id <- as.character(as.integer(cpdids))
	sp@name <- findName(cpdids)
	sp@formula <- findFormula(cpdids)
	sp@mode <- mode

	if(length(w@spectra) != 0){
		IDindex <- sapply(w@spectra,function(s) s@id == cpdids)
		if(length(IDindex)){
			spectraNum <- length(w@spectra[[which(IDindex)]]@children)
			w@spectra[[which(IDindex)]]@children[[spectraNum+1]] <- sp@children[[1]]
		} else {
			w@spectra[[length(w@spectra)+1]] <- sp
		}
	} else{
		w@spectra[[1]] <- sp
	}
	
	if(all(w@files != xRAW[[1]]@filepath)){
		w@files <- c(w@files,xRAW[[1]]@filepath)
	} else{
		for(i in 2:(length(w@files)+1)){
			currentFPath <- paste0(xRAW[[1]]@filepath,"_",i)
			if(all(w@files != currentFPath)){
				w@files <- c(w@files,currentFPath)
				break
			}
		}
	}
	
	return(w)
}


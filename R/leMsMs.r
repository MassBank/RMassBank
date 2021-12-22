
#library(xcms)

#' Backup \code{msmsWorkflow} results
#' 
#' Writes the results from different \code{msmsWorkflow} steps to a file.
#' 
#' @aliases archiveResults
#' @usage archiveResults(w, fileName, settings = getOption("RMassBank"))
#' @param w The \code{msmsWorkspace} to be saved.
#' @param fileName The filename to store the results under.
#' @param settings The settings to be stored into the msmsWorkspace image.
#' @examples 
#' 
#' 		# This doesn't really make a lot of sense,
#' 		# it stores an empty workspace.
#' 		RmbDefaultSettings()
#' 		w <- newMsmsWorkspace()
#' 		archiveResults(w, "narcotics.RData")
#' 
#' @export
archiveResults <- function(w, fileName, settings = getOption("RMassBank"))
{
  # save the settings into the settings slot
  w@settings <- settings
    # save
  save(w, file=fileName)

}


#' RMassBank mass spectrometry pipeline
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
#' @param mode \code{"pH", "pNa", "pM", "mH", "mM", "mFA", "pNH4"} for different ions 
#' 			([M+H]+, [M+Na]+, [M]+, [M-H]-, [M]-, [M+FA]-, [M+NH4]+).
#' @param steps Which steps of the workflow to process. See the vignette 
#' 			\code{vignette("RMassBank")} for details.
#' @param confirmMode Defaults to false (use most intense precursor). Value 1 uses
#' 			the 2nd-most intense precursor for a chosen ion (and its data-dependent scans)
#' 			, etc.
#' @param newRecalibration Whether to generate a new recalibration curve (\code{TRUE}, default) or
#' 			to reuse the currently stored curve (\code{FALSE}, useful e.g. for adduct-processing runs.) 
#' @param useRtLimit Whether to enforce the given retention time window.
#' @param archivename The prefix under which to store the analyzed result files.
#' @param readMethod Several methods are available to get peak lists from the files.
#'        Currently supported are "mzR", "xcms", "MassBank" and "peaklist".
#'        The first two read MS/MS raw data, and differ in the strategy 
#'        used to extract peaks. MassBank will read existing records, 
#'        so that e.g. a recalibration can be performed, and "peaklist" 
#'        just requires a CSV with two columns and the column header "mz", "int".
#' @param findPeaksArgs A list of arguments that will be handed to the xcms-method findPeaks via do.call
#' @param plots A parameter that determines whether the spectra should be plotted or not (This parameter is only used for the xcms-method)
#' @param precursorscan.cf Whether to fill precursor scans. To be used with files which for
#' 		some reasons do not contain precursor scan IDs in the mzML, e.g. AB Sciex converted
#' 		files.
#' @param settings Options to be used for processing. Defaults to the options loaded via
#' 			\code{\link{loadRmbSettings}} et al. Refer to there for specific settings.
#' @param analyzeMethod The "method" parameter to pass to \code{\link{analyzeMsMs}}.
#' @param progressbar The progress bar callback to use. Only needed for specialized applications.
#' 			Cf. the documentation of \code{\link{progressBarHook}} for usage.
#' @param MSe A boolean value that determines whether the spectra were recorded using MSe or not
#' @return The processed \code{msmsWorkspace}.
#' @seealso \code{\link{msmsWorkspace-class}}
#' @author Michael Stravs, Eawag <michael.stravs@@eawag.ch>
#' @export
msmsWorkflow <- function(w, mode="pH", steps=c(1:8), confirmMode = FALSE, newRecalibration = TRUE, 
		useRtLimit = TRUE, archivename=NA, readMethod = "mzR", filetable=NULL, findPeaksArgs = NULL, plots = FALSE,
		precursorscan.cf = FALSE,
		settings = getOption("RMassBank"), analyzeMethod = "formula",
		progressbar = "progressBarHook", MSe = FALSE)
{
    .checkMbSettings()
    if(!any(mode %in% knownAdducts())) stop(paste("The ionization mode", mode, "is unknown."))

  if(!is.na(archivename))
	  w@archivename <- archivename
  
  # Make a progress bar:
  nProg <- 0
  nLen <- length(w@files)
    
    allUnknown <- FALSE
    
    # If all compounds are unknown some specific conditions apply
  allUnknown <- FALSE
    if(all(.listEnvEnv$listEnv$compoundList$Level == "5")){
        allUnknown <- TRUE
        rmb_log_info("All compounds are unknown, the workflow will be adjusted accordingly")
    }
    
    if(readMethod == "minimal"){
	##Edit options
	opt <- getOption("RMassBank")
	opt$recalibrator$MS1 <- "recalibrate.identity"
	opt$recalibrator$MS2 <- "recalibrate.identity"
        opt$add_annotation <- FALSE
        opt$multiplicityFilter <- 1
	options(RMassBank=opt)
        settings <- getOption("RMassBank")
	##Edit analyzemethod
	analyzeMethod <- "intensity"
  }

  # clean rerun functionality:
  # if any step after 3 has been run, rerunning steps 4 or below needs moving back to the parent workspace.
  # However, the recalibration must be preserved, because:
  # if someone runs       
  # w <- msmsWorkflow(w, steps=c(1:4)),
  # then substitutes the recalibration
  # w@rc <- myrecal
  # then runs step 4 again
  # w <- msmsWorkflow(w, steps=c(4), newRecalibration=FALSE)
  # the rc and rc.ms1 must be preserved and not taken from the parent workspace
  if(!all(steps > 4) & !is.null(w@parent))
  {
    rc <- w@rc
    rc.ms1 <- w@rc.ms1
    w <- w@parent
    w@rc <- rc
    w@rc.ms1 <- rc.ms1
  }
  
  # Step 1: acquire all MSMS spectra from files
  if(1 %in% steps)
  {
    rmb_log_info("msmsWorkflow: Step 1. Acquire all MSMS spectra from files")
	w <- msmsRead(w = w, files = w@files, readMethod=readMethod, filetable=filetable, mode=mode, confirmMode = confirmMode, useRtLimit = useRtLimit, 
                        Args = findPeaksArgs, settings = settings, progressbar = progressbar, MSe = MSe)
  }
  # Step 2: first run analysis before recalibration
  if(2 %in% steps)
  {
	  nProg <- 0
	  rmb_log_info("msmsWorkflow: Step 2. First analysis pre recalibration")
        if(allUnknown){
            analyzeMethod <- "intensity"
        }
	  pb <- do.call(progressbar, list(object=NULL, value=0, min=0, max=nLen))
	  w@spectra <- as(lapply(w@spectra, function(spec) {
                        s <- analyzeMsMs(msmsPeaks = spec, mode=mode, detail=TRUE, run="preliminary",
						  filterSettings = settings$filterSettings,
						  spectraList = settings$spectraList, method = analyzeMethod)
				  # Progress:
				  nProg <<- nProg + 1
				  pb <- do.call(progressbar, list(object=pb, value= nProg))

				  return(s)
			  }), "SimpleList")
	## for(f in w@files)
	## w@spectra[[basename(as.character(f))]]@name <- basename(as.character(f))
	  suppressWarnings(do.call(progressbar, list(object=pb, close=TRUE)))
  }
  # Step 3: aggregate all spectra
  if(3 %in% steps)
  {
	rmb_log_info("msmsWorkflow: Step 3. Aggregate all spectra")
        w@aggregated <- aggregateSpectra(spec = w@spectra, addIncomplete=TRUE)
        
        if(RMassBank.env$verbose.output){
          numberOfPeaksThere <- sum(unlist(lapply(X = w@spectra, FUN = function(spec){ sum(unlist(lapply(X = spec@children, FUN = function(child){ child@peaksCount }))) })))
          if(nrow(w@aggregated) < numberOfPeaksThere)
            cat(paste("### Warning ### The aggregation of spectra lead to the removal of ", (numberOfPeaksThere-nrow(w@aggregated)), " / ", numberOfPeaksThere, " peaks\n", sep = ""))
        }
  }
    
    if(allUnknown){
        w@aggregated$noise <- FALSE
        w@aggregated$noise <- FALSE
        w@aggregated$reanalyzed.formula <- NA
        w@aggregated$reanalyzed.mzCalc <- NA
        w@aggregated$reanalyzed.dppm <- NA
        w@aggregated$reanalyzed.formulaCount <- NA
        w@aggregated$reanalyzed.dbe <- NA
        w@aggregated$matchedReanalysis <- NA
        w@aggregated$filterOK <- TRUE
        w@aggregated$problematicPeak <- FALSE
        w@aggregated$formulaMultiplicity <- unlist(sapply(table(w@aggregated$cpdID),function(x) rep(x,x)))
        return(w)
    }
    
    
  # Step 4: recalibrate all m/z values in raw spectra
  if(4 %in% steps)
  {
	rmb_log_info("msmsWorkflow: Step 4. Recalibrate m/z values in raw spectra")
	if(newRecalibration)
	{
		# note: makeRecalibration takes w as argument now, because it needs to get the MS1 spectra from @spectra
		recal <- makeRecalibration(w, mode,
				recalibrateBy = settings$recalibrateBy,
				recalibrateMS1 = settings$recalibrateMS1,
				recalibrator = settings$recalibrator,
				recalibrateMS1Window = settings$recalibrateMS1Window)
		w@rc <- recal$rc
		w@rc.ms1 <- recal$rc.ms1
	}
	w@parent <- w
	w@aggregated <- data.frame()
    spectra <- recalibrateSpectra(mode, w@spectra, w = w,
			recalibrateBy = settings$recalibrateBy,
			recalibrateMS1 = settings$recalibrateMS1)
	w@spectra <- spectra
  }
  # Step 5: re-analysis on recalibrated spectra
  if(5 %in% steps)
  {
	nProg <- 0
	rmb_log_info("msmsWorkflow: Step 5. Reanalyze recalibrated spectra")
	pb <- do.call(progressbar, list(object=NULL, value=0, min=0, max=nLen))
	
	w@spectra <- as(lapply(w@spectra, function(spec) {
						#print(spec$id)
                            if(findLevel(spec@id,TRUE) == "unknown"){
                                analyzeMethod <- "intensity"
                            } else {
                                analyzeMethod <- "formula"
                            }
                            s <- analyzeMsMs(msmsPeaks = spec, mode=mode, detail=TRUE, run="recalibrated",
								filterSettings = settings$filterSettings,
								spectraList = settings$spectraList, method = analyzeMethod)
						# Progress:
						nProg <<- nProg + 1
						pb <- do.call(progressbar, list(object=pb, value= nProg))
						
						return(s)
					}), "SimpleList")
	## for(f in w@files)
	## w@spectra[[basename(as.character(f))]]@name <- basename(as.character(f))
	suppressWarnings(do.call(progressbar, list(object=pb, close=TRUE)))
	
  do.call(progressbar, list(object=pb, close=TRUE))
  }
  # Step 6: aggregate recalibrated results
  if(6 %in% steps)
  {
    rmb_log_info("msmsWorkflow: Step 6. Aggregate recalibrated results")
        w@aggregated <- aggregateSpectra(spec = w@spectra, addIncomplete=TRUE)
        
        if(RMassBank.env$verbose.output){
          numberOfPeaksThere <- sum(unlist(lapply(X = w@spectra, FUN = function(spec){ sum(unlist(lapply(X = spec@children, FUN = function(child){ child@peaksCount }))) })))
          if(nrow(w@aggregated) < numberOfPeaksThere)
            cat(paste("### Warning ### The aggregation of spectra lead to the removal of ", (numberOfPeaksThere-nrow(w@aggregated)), " / ", numberOfPeaksThere, " peaks\n", sep = ""))
        }
        
    if(!is.na(archivename))
        archiveResults(w, paste(archivename, ".RData", sep=''), settings)

    ## clean electronic noise if specified in the settings.ini file
    if (length(settings$electronicNoise) > 0 && settings$electronicNoiseWidth > 0) {  
        spectra <- lapply(w@spectra, cleanElnoise, noise=settings$electronicNoise, width=settings$electronicNoiseWidth)
        w@spectra <- as(spectra, "SimpleList")
        if(RMassBank.env$verbose.output)
            if(sum(w@aggregated$noise) > 0)
                cat(paste("### Warning ### ", sum(w@aggregated$noise), " / ", nrow(w@aggregated), " peaks have been identified as electronic noise\n", sep = ""))
    }
  }  
  # Step 7: reanalyze failpeaks for (mono)oxidation and N2 adduct peaks
  if(7 %in% steps)
  {
	rmb_log_info("msmsWorkflow: Step 7. Reanalyze fail peaks for N2 + O")
    w <- reanalyzeFailpeaks(
			w, custom_additions="N2O", mode=mode,
				filterSettings=settings$filterSettings,
                    progressbar=progressbar)
    if(!is.na(archivename))
      archiveResults(w, paste(archivename, "_RA.RData", sep=''), settings)
        
        if(RMassBank.env$verbose.output){
          isNoFormula <- is.na(w@aggregated$formula) & is.na(w@aggregated$reanalyzed.formula)
          noFormulaCount <- sum(isNoFormula)
          numberOfPeaksThere <- sum(unlist(lapply(X = w@spectra, FUN = function(spec){ sum(unlist(lapply(X = spec@children, FUN = function(child){ child@peaksCount }))) })))
          if(noFormulaCount > 0){
            cat(paste("### Warning ### ", noFormulaCount, " / ", numberOfPeaksThere, " peaks have no molecular formula:\n", sep = ""))
            print(w@aggregated[isNoFormula, c("mzFound","intensity","cpdID")])
          }
        }
    }
    
  # Step 8: heuristic filtering based on peak multiplicity;
  #         creation of failpeak list
  if(8 %in% steps)
  {
	rmb_log_info("msmsWorkflow: Step 8. Peak multiplicity filtering")
    if (is.null(settings$multiplicityFilter)) {
      rmb_log_info("msmsWorkflow: Step 8. Peak multiplicity filtering skipped because multiplicityFilter parameter is not set.")
          w@aggregated <- addProperty(w@aggregated, "formulaMultiplicity", "integer", 1)
          w@aggregated <- addProperty(w@aggregated, "filterOK",            "logical", FALSE)
          w@aggregated$filterOK <- !((is.na(w@aggregated$formulaCount) | w@aggregated$formulaCount==0) & (is.na(w@aggregated$reanalyzed.formulaCount) | w@aggregated$reanalyzed.formulaCount==0))
          w@aggregated <- addProperty(w@aggregated, "problematicPeak",     "logical", FALSE)
    } else {
      # apply heuristic filter      
      w <- filterMultiplicity(
                w, archivename, mode, multiplicityFilter = settings$multiplicityFilter)
            if(RMassBank.env$verbose.output){
              peakDfs <- split(x = w@aggregated, f = list("mzFound"=w@aggregated$mzFound, "cpdID"=w@aggregated$cpdID))
              numberOfPeaks <- length(peakDfs)
              multiplicityNotOkCount <- numberOfPeaks - sum(unlist(lapply(X = peakDfs, FUN = function(x){any(x$filterOK)})))
              if(multiplicityNotOkCount > 0)
                cat(paste("### Warning ### ", multiplicityNotOkCount, " / ", numberOfPeaks, " peaks do not fulfill the multiplicity criterion\n", sep = ""))
            }
            
            
            if(RMassBank.env$verbose.output){
              problematicPeakCount <- sum(w@aggregated$problematicPeak)
              if(problematicPeakCount > 0)
                cat(paste("### Warning ### ", problematicPeakCount, " / ", nrow(w@aggregated), " peaks are problematic\n", sep = ""))
            }
            
      if(!is.na(archivename))
        archiveResults(w, paste(archivename, "_RF.RData", sep=''), settings)   
    }
  }
  rmb_log_info("msmsWorkflow: Done.")
  return(w)
}

#' Analyze MSMS spectra
#' 
#' Analyzes MSMS spectra of a compound by fitting formulas to each subpeak.
#' 
#' The analysis function uses Rcdk. Note
#' that in this step, \emph{satellite peaks} are removed by a simple heuristic
#' rule (refer to the documentation of \code{\link{filterPeakSatellites}} for details.)
#' 
## # @usage analyzeMsMs(msmsPeaks, mode="pH", detail=FALSE, run="preliminary",
## # 			filterSettings = getOption("RMassBank")$filterSettings,
## # 			spectraList = getOption("RMassBank")$spectraList, method="formula")
## # 
## # 		analyzeMsMs.formula(msmsPeaks, mode="pH", detail=FALSE, run="preliminary",
## # 			filterSettings = getOption("RMassBank")$filterSettings,
## # 			spectraList = getOption("RMassBank")$spectraList)
## # 
## # 		analyzeMsMs.intensity(msmsPeaks, mode="pH", detail=FALSE, run="preliminary",
## # 			filterSettings = getOption("RMassBank")$filterSettings,
## # 			spectraList = getOption("RMassBank")$spectraList)
#' 
#' @param msmsPeaks A \code{RmbSpectraSet} object.
#' 	Corresponds to a parent spectrum and children MSMS spectra of one compound (plus some metadata).
#'  The objects are typically generated with \code{\link{findMsMsHR}}, and populate the \code{@@spectrum} slot
#'  in a \code{msmsWorkspace} (refer to the corresponding
#' documentation for the precise format specifications).
#' @param mode Specifies the processing mode, i.e. which molecule species the
#' spectra contain. \code{\var{pH}} (positive H) specifies [M+H]+,
#' \code{\var{pNa}} specifies [M+Na]+, \code{\var{pM}} specifies [M]+,
#' \code{\var{mH}} and \code{\var{mNa}} specify [M-H]- and [M-Na]-,
#' respectively. (I apologize for the naming of \code{\var{pH}} which has
#' absolutely nothing to do with chemical \emph{pH} values.)
#' @param detail Whether detailed return information should be provided
#' (defaults to \code{FALSE}). See below.
#' @param run \code{"preliminary"} or \code{"recalibrated"}. In the
#' \code{preliminary} run, mass tolerance is set to 10 ppm (above m/z 120) and
#' 15 ppm (below m/z 120), the default intensity cutoff is $10^4$ for positive
#' mode (no default cutoff in negative mode), and the column \code{"mz"} from
#' the spectra is used as data source.  In the \code{recalibrated} run, the
#' mass tolerance is set to 5 ppm over the whole mass range, the default cutoff
#' is 0 and the column \code{"mzRecal"} is used as source for the m/z values.
#' Defaults to \code{"preliminary"}.
#' @param filterSettings
#' 		Settings for the filter parameters, by default loaded from the RMassBank settings
#' 		set with e.g. \code{\link{loadRmbSettings}}. Must contain:
#' 		\itemize{
#' 			\item \code{ppmHighMass}, allowed ppm deviation before recalibration
#' 				for high mass range
#' 			\item \code{ppmLowMass}, allowed ppm deviation before recalibration
#' 				for low mass range
#' 			\item \code{massRangeDivision}, division point between high and low mass
#' 				range (before recalibration)
#' 			\item \code{ppmFine}, allowed ppm deviation overall after recalibration
#' 			\item \code{prelimCut}, intensity cutoff for peaks in preliminary run
#' 			\item \code{prelimCutRatio}, relative intensity cutoff for peaks in 
#' 				preliminary run, e.g. 0.01 = 1%
#' 			\item \code{fineCut}, intensity cutoff for peaks in second run
#' 			\item \code{fineCutRatio}, relative intensity cutoff for peaks in 
#' 				second run
#' 			\item \code{specOkLimit}, minimum intensity of base peak for spectrum
#' 				to be accepted for processing
#' 			\item \code{dbeMinLimit}, minimum double bond equivalent for accepted
#' 				molecular subformula.
#' 			\item \code{satelliteMzLimit}, for satellite peak filtering 
#' 				(\code{\link{filterPeakSatellites}}: mass window to use for satellite
#' 				removal
#' 			\item \code{satelliteIntLimit}, the relative intensity below which to 
#' 				discard "satellites". (refer to  \code{\link{filterPeakSatellites}}).
#' 	}
#' @param spectraList The list of MS/MS spectra present in each data block. As also
#' 		defined in the settings file.  
#' @param method Selects which function to actually use for data evaluation. The default
#' 		"formula" runs a full analysis via formula assignment to fragment peaks. The
#' 		alternative setting "intensity" calls a "mock" implementation which circumvents
#' 		formula assignment and filters peaks purely based on intensity cutoffs and the
#' 		satellite filtering. (In this case, the ppm and dbe related settings in filterSettings
#' 		are ignored.)
#' @return The processed \code{RmbSpectraSet} object.
#' Added (or filled, respectively, since the slots are present before) data include
#' \item{list("complete")}{whether all spectra have useful value}
#' \item{list("empty")}{whether there are no useful spectra}
#' \item{list("children")}{
#' 		The processed \code{RmbSpectrum2} objects (in a \code{RmbSpectrum2List}).
#' 		\itemize{
#' 			\item \code{ok} if the spectrum was successfully processed with at least one resulting peak
#' 			\item \code{mz}, \code{intensity}: note that mz/int pairs can be duplicated when multiple matches
#' 				are found for one mz value, therefore the two slots are not necessarily unchanged from before
#'			\item \code{rawOK} (logical) whether the m/z peak passes satellite/low removal
#' 			\item \code{low}, \code{satellite} if \code{TRUE}, the peak failed cutoff (\code{low}) or was removed as \code{satellite}
#' 			\item \code{formula}, \code{mzCalc}, \code{dppm}, \code{dbe} Formula, calculated mass, ppm deviation and dbe assigned to a peak
#' 			\item \code{formulaCount}, \code{dppmBest} Number of formulae matched for this m/z value and ppm deviation of the best match
#' 			\item \code{info} Spectrum identifying information (collision energy, resolution, collision mode) from
#' 				the \code{spectraList} 
#' 			\item All other entries are retained from the original \code{RmbSpectrum2}.
#' 			}
#' }
#' @aliases analyzeMsMs analyzeMsMs.formula analyzeMsMs.intensity
#' @author Michael Stravs
#' @seealso \code{\link{msmsWorkflow}}, \code{\link{filterLowaccResults}},
#' \code{\link{filterPeakSatellites}}, \code{\link{reanalyzeFailpeaks}}
#' @examples
#' 
#' 	\dontrun{analyzed <- analyzeMsMs(spec, "pH", TRUE)}
#' 
#' @export
analyzeMsMs <- function(msmsPeaks, mode="pH", detail=FALSE, run="preliminary",
		filterSettings = getOption("RMassBank")$filterSettings,
		spectraList = getOption("RMassBank")$spectraList, method="formula")
{
	## .RmbSpectraSet <- setClass("RmbSpectraSet",
	##         representation = representation(
	##                 parent = "Spectrum1",
	##                 children = "RmbSpectrum2List",
	##                 # These are done as slots and not as S4 functions, because they are set during the workflow
	##                 # in "checking" steps. It's easier.
	##                 found = "logical",
	##                 complete = "logical",
	##                 empty = "logical",
	##                 formula = "character",
	##                 id = "integer",
	##                 mz = "numeric",
	##                 name = "character",
	##                 annotations = "list"
	##         ),
	##         prototype = prototype(
	##                 parent = new("Spectrum1"),
	##                 children = new("RmbSpectrum2List"),
	##                 found = FALSE,
	##                 complete = NA,
	##                 empty = NA,
	##                 formula = character(),
	##                 id = integer(),
	##                 mz = numeric(),
	##                 name = character(),
	##                 annotations = list()
	##         )
	## );
	.checkMbSettings()
	
	if(msmsPeaks@found == FALSE)
		return(msmsPeaks)
	
	# Check whether the spectra can be fitted to the spectra list correctly!
	if(length(msmsPeaks@children) != length(spectraList))
	{
		warning(paste0(
						"The spectra count of the substance ", msmsPeaks@id, " (", length(msmsPeaks@children), " spectra) doesn't match the provided spectra list (", length(spectraList), " spectra)."
				))
		msmsPeaks@found <- FALSE
		return(msmsPeaks)
		
	}
	
	if(method=="formula")
	{
		r <- analyzeMsMs.formula(msmsPeaks, mode, detail, run, filterSettings)
	}
	else if(method == "intensity")
	{
		r <- analyzeMsMs.intensity(msmsPeaks, mode, detail, run, filterSettings)
	}
	
	# Add the spectrum labels to the spectra here.
	# If there is any better place to do this, please tell me. I hate it.
	# However, the info should be added in msmsWorkflow not in mbWorkflow, because two msmsWorkspaces with different spectraLists can be
	# merged together with all the combine / pack stuff.
	children <- mapply(function(spec, info)
			{
				spec@info <- c(info, spec@info)
				spec
			}, r@children, spectraList, SIMPLIFY=FALSE)
	r@children <- as(children, "SimpleList")
	
	
	#nspectra <- length(spectraList)
	ok <- unlist(lapply(r@children, function(c) c@ok))
	r@complete <- FALSE
	r@empty <- FALSE
	if(all(ok))
		r@complete <- TRUE
	if(all(!ok))
		r@empty <- TRUE
	return(r)
}


#' @describeIn analyzeMsMs Analyze the peaks using formula annotation
#' @export
analyzeMsMs.formula <- function(msmsPeaks, mode="pH", detail=FALSE, run="preliminary",
			filterSettings = getOption("RMassBank")$filterSettings)
{
  cut <- 0
  cut_ratio <- 0
  if(run=="preliminary"){
    filterMode <- "coarse"
	  cut <- filterSettings$prelimCut
    if(is.na(cut))
    {
      adductProperties <- getAdductProperties(mode, msmsPeaks@formula)
      if(adductProperties$charge > 0) cut <- 1e4
      if(adductProperties$charge < 0) cut <- 0
    }
	cut_ratio <- filterSettings$prelimCutRatio
  } else {
    filterMode <- "fine"
	cut <- filterSettings$fineCut
	cut_ratio <- filterSettings$fineCutRatio
    if(is.na(cut)) cut <- 0
  }

  # find whole spectrum of parent peak, so we have reasonable data to feed into
  # MolgenMsMs
  parentSpectrum <- msmsPeaks@parent

  
  # On each spectrum the following function analyzeTandemShot will be applied.
  # It takes the raw peaks matrix as argument (mz, int) and processes the spectrum by
  # filtering out low-intensity (<1e4) and shoulder peaks (deltam/z < 0.5, intensity
  # < 5%) and subsequently matching the peaks to formulas using Rcdk, discarding peaks
  # with insufficient match accuracy or no match.
  analyzeTandemShot <- function(child, childIdx = 0)
  {
  	shot <- getData(child)
  	shot$row <- which(!is.na(shot$mz))
  	
  	
    # Filter out low intensity peaks:
    child@low <- (shot$intensity < cut) | (shot$intensity < max(shot$intensity)*cut_ratio)
    shot <- shot[!child@low,,drop=FALSE]
    shot_full <- shot
    
    # Is there still anything left?
    if(length(which(!child@low))==0)
  	{
  		child@ok <- FALSE
  		
  		if(RMassBank.env$verbose.output)
        cat(paste("\n### Warning ### The spectrum '#", childIdx, "' for ID '", msmsPeaks@id, "' contains only low-intensity peaks\n", sep = ""))
  		
  		return(child)
  	}
      
    # Filter out satellite peaks:
    shot <- filterPeakSatellites(shot, filterSettings)
  	child@satellite <- rep(TRUE, child@peaksCount)
  	child@satellite[which(child@low == TRUE)] <- NA
  	child@satellite[shot$row] <- FALSE
  	
    # Is there still anything left?
    if(nrow(shot)==0)
  	{
  		child@ok <- FALSE
  		
  		if(RMassBank.env$verbose.output)
  		  cat(paste("\n### Warning ### The spectrum '#", childIdx, "' for ID '", msmsPeaks@id, "' contains no peaks after satellite filtering\n", sep = ""))
  		
  		return(child)
  	}
      
    if(max(shot$intensity) < as.numeric(filterSettings$specOkLimit))
  	{
  		child@ok <- FALSE
  		
  		if(RMassBank.env$verbose.output)
  		  cat(paste("\n### Warning ### The spectrum '#", childIdx, "' for ID '", msmsPeaks@id, "' is discarded due to parameter 'specOkLimit'\n", sep = ""))
  		
  		return(child)
  	}
  	
    # Crop to 4 digits (necessary because of the recalibrated values)
  	# this was done for the MOLGEN MSMS type analysis, is not necessary anymore now (23.1.15 MST)
    # shot[,mzColname] <- round(shot[,mzColname], 5)
      
  	# here follows the Rcdk analysis
  	#------------------------------------
  	parentPeaks <- data.frame(mzFound=msmsPeaks@mz, 
  			formula=msmsPeaks@formula,
  			dppm=0,
  			x1=0,x2=0,x3=0)
  	
  	# get the adduct additions
  	adductProperties <- getAdductProperties(mode, msmsPeaks@formula)
  	allowed_additions <- adductProperties$addition
  	mode.charge <- adductProperties$charge
  	
  	# the ppm range is two-sided here.
  	# The range is slightly expanded because dppm calculation of
  	# generate.formula starts from empirical mass, but dppm cal-
  	# culation of the evaluation starts from theoretical mass.
  	# So we don't miss the points on 'the border'.
  	
  	if(run=="preliminary")
  		ppmlimit <- 2 * max(filterSettings$ppmLowMass, filterSettings$ppmHighMass)
  	else
  		ppmlimit <- 2.25 * filterSettings$ppmFine
  	
  	parent_formula <- add.formula(msmsPeaks@formula, allowed_additions)
  	dbe_parent <- dbe(parent_formula)
  	# check whether the formula is valid, i.e. has no negative or zero element numbers.
  	#print(parent_formula)
  	if(!is.valid.formula(parent_formula))
  	{
  		child@ok <- FALSE
  		
  		if(RMassBank.env$verbose.output)
  		  cat(paste("\n### Warning ### The precursor ion formula of spectrum '#", childIdx, "' for ID '", msmsPeaks@id, "' is invalid\n", sep = ""))
  		
  		return(child)
  	}
  
  	limits <- to.limits.rcdk(parent_formula)
	maximal.mass <- get.formula(parent_formula)@mass
  	
  	peakmatrix <- lapply(split(shot,shot$row), function(shot.row)  {
  				# Circumvent bug in rcdk: correct the mass for the charge first, then calculate uncharged formulae
  				# finally back-correct calculated masses for the charge
  				mass <- shot.row[["mz"]]
  				mass.calc <- mass + mode.charge * .emass
				tolerance <- ppm(mass.calc, ppmlimit, p=TRUE)
				if (mass.calc > maximal.mass + tolerance)
  					return(t(c(row=shot.row[["row"]], intensity = shot.row[["intensity"]], mz=mass, formula=NA, mzCalc=NA)))
				peakformula <- tryCatch(
				  suppressWarnings(generate.formula(mass.calc, tolerance,
				                                    limits, charge=0)),
				  error = function(e) list())
				# was a formula found? If not, return empty result
  				# ppm(mass, ppmlimit, p=TRUE),
  				# limits, charge=1),
  				#error= function(e) list())
  			
			if(length(peakformula)==0)
  				return(t(c(row=shot.row[["row"]], intensity = shot.row[["intensity"]], mz=mass,
  										formula=NA, mzCalc=NA)))
  			else
  			{
  				return(t(sapply(peakformula, function(f)
  										{
  											mzCalc <- f@mass - mode.charge * .emass
  											c(row=shot.row[["row"]], intensity = shot.row[["intensity"]], mz=mass,
  													formula=f@string, 
  													mzCalc=mzCalc)
  										})))
  			}
  	})
  	
  	childPeaks <- as.data.frame(do.call(rbind, peakmatrix))
  	
  	presentElements <- unique(unlist(lapply(X = lapply(X = childPeaks$formula, FUN = formulastring.to.list), FUN = names)))
  	## function dbe from FormulaCalculator.R
  	atomDBEs <- sapply(X = presentElements, FUN = dbe)
  	unknownElements <- names(atomDBEs)[sapply(X = atomDBEs, FUN = function(atomDBE){length(atomDBE)==0})]
  	if(length(unknownElements) > 0) stop(paste("Element(s)", paste(unknownElements, collapse = "; "), "cannot be assigned a DBE"))
  	
  	# Reformat the deformatted output correctly (why doesn't R have a better way to do this, e.g. avoid deformatting?)
  
  	childPeaks$row <- as.numeric(as.character(childPeaks$row))
  	childPeaks$intensity <- as.numeric(as.character(childPeaks$intensity))
  	childPeaks$mz <- as.numeric(as.character(childPeaks$mz))
  	childPeaks$formula <- as.character(childPeaks$formula)
  	childPeaks$mzCalc <- as.numeric(as.character(childPeaks$mzCalc))
  	childPeaks$dppm <- (childPeaks$mz / childPeaks$mzCalc - 1) * 1e6
  	childPeaks$dbe <- unlist(lapply(childPeaks$formula, dbe))
  	
  	# childPeaks now contains all the good and unmatched peaks
  	# but not the ones which were cut as satellites or below threshold.
  	
  	## child@mzFound <- rep(NA, child@peaksCount)
  	## child@mzFound[childPeaks$row] <- as.numeric(as.character(childPeaks$mzFound))
  	## 
  	## child@formula <- rep(NA, child@peaksCount)
  	## child@formula[childPeaks$row] <- as.character(childPeaks$formula)
  	## 
  	## child@mzCalc <- rep(NA, child@peaksCount)
  	## child@mzCalc[childPeaks$row] <- as.numeric(as.character(childPeaks$mzCalc))
  	## 
  	## child@dppm<- rep(NA, child@peaksCount)
  	## child@dppm[childPeaks$row] <- (childPeaks$mzFound / childPeaks$mzCalc - 1) * 1e6
  	# delete the NA data out again, because MolgenMsMs doesn't have them
  	# here and they will be re-added later
  	# (this is just left like this for "historical" reasons)
  	#childPeaks <- childPeaks[!is.na(childPeaks$formula),]
  	# check if a peak was recognized (here for the first time,
  	# otherwise the next command would fail)
  
  	if(nrow(childPeaks)==0)
  	{
  		child@ok <- FALSE
  		
  		chdata <- getData(child)
  		chdata$formulaSource <- character(nrow(chdata))
  		child <- setData(child, chdata)
  		
  		
  		if(RMassBank.env$verbose.output)
  		  cat(paste("\n### Warning ### The spectrum '#", childIdx, "' for ID '", msmsPeaks@id, "' is empty\n", sep = ""))
  		
  		return(child)
  	}
    
  	if(all(is.na(childPeaks$formula)))
  	{
  	  child@ok <- FALSE
  	  child@good <- rep(FALSE, length(childPeaks$formula))
  	  
  	  chdata <- getData(child)
  	  chdata$formulaSource <- character(nrow(chdata))
  	  child <- setData(child, chdata)
  	  
  	  
  	  if(RMassBank.env$verbose.output)
  	    cat(paste("\n### Warning ### The spectrum '#", childIdx, "' for ID '", msmsPeaks@id, "' comprises no peaks which could be assiged to a molecular formula\n", sep = ""))
  	  
  	  return(child)
  	}
  	
  	# now apply the rule-based filters to get rid of total junk:
  	# dbe >= -0.5, dbe excess over mother cpd < 3
  	# dbe() has been adapted to return NA for NA input
  	#iff_rcdk_pM_eln$maxvalence <- unlist(lapply(diff_rcdk_pM_eln$formula.rcdk, maxvalence))
  	temp.child.ok <- (childPeaks$dbe >= filterSettings$dbeMinLimit) 
  	# & dbe < dbe_parent + 3)
  	# check if a peak was recognized
  	if(length(which(temp.child.ok)) == 0)
  	{
  		child@ok <- FALSE
  		child@good <- rep(FALSE, length(temp.child.ok))
  		
  		chdata <- getData(child)
  		chdata$formulaSource <- character(nrow(chdata))
  		child <- setData(child, chdata)
  		
  		if(RMassBank.env$verbose.output)
  		  cat(paste("\n### Warning ### The spectrum '#", childIdx, "' for ID '", msmsPeaks@id, "' comprises no peaks which fulfil the dbeMinLimit criterion\n", sep = ""))
  		
  		return(child)
  	}
    #browser()	
  	# find the best ppm value
    bestPpm <- aggregate(
      x   = as.data.frame(childPeaks[!is.na(childPeaks$dppm),"dppm"]),
  		by  = list(childPeaks[!is.na(childPeaks$dppm),"row"]),
      FUN = function(dppm) dppm[[which.min(abs(dppm))]]
  	)
    colnames(bestPpm) <- c("row", "dppmBest")
    childPeaks <- merge(childPeaks, bestPpm, by="row", all.x=TRUE)
  	
  	# Deactivated the following lines because we never actually want to look at the "old" formula count.
  	# To be verified (cf Refiltering, failpeak list and comparable things) 
  
  	## # count formulas found per mass
  	## countFormulasTab <- xtabs( ~formula + mz, data=childPeaks)
  	## countFormulas <- colSums(countFormulasTab)
  	## childPeaks$formulaCount <- countFormulas[as.character(childPeaks$row)]
  	
    # filter results
    childPeaksFilt <- filterLowaccResults(childPeaks, filterMode, filterSettings)
    childPeaksGood <- childPeaksFilt[["TRUE"]]
    childPeaksBad <- childPeaksFilt[["FALSE"]]
  	if(is.null(childPeaksGood)){
  		childPeaksGood <- childPeaks[c(),,drop=FALSE]
      childPeaksGood$good <- logical(0)
    }
  	if(is.null(childPeaksBad))
  		childPeaksBad <- childPeaks[c(),,drop=FALSE]
  	childPeaksUnassigned <- childPeaks[is.na(childPeaks$dppm),,drop=FALSE]
  	childPeaksUnassigned$good <- rep(FALSE, nrow(childPeaksUnassigned))
    
  	# count formulas within new limits
    # (the results of the "old" count stay in childPeaksInt and are returned
    # in $childPeaks)
  	countFormulasTab <- xtabs( ~formula + mz, data=childPeaksGood)
  	countFormulas <- colSums(countFormulasTab)
  	childPeaksGood$formulaCount <- countFormulas[as.character(childPeaksGood$mz)]
  	  
  	childPeaksUnassigned$formulaCount <- rep(0, nrow(childPeaksUnassigned))
  	childPeaksBad$formulaCount <- rep(0, nrow(childPeaksBad))
  	childPeaksBad$good <- rep(FALSE, nrow(childPeaksBad))
      
  	# Now: childPeaksGood (containing the new, recounted peaks with good = TRUE), and childPeaksBad (containing the 
  	# peaks with good=FALSE, i.e. outside filter criteria, with the old formula count even though it is worthless)
  	# are bound together.
  	childPeaksBad <- childPeaksBad[,colnames(childPeaksGood),drop=FALSE]
  	childPeaksUnassigned <- childPeaksUnassigned[,colnames(childPeaksGood),drop=FALSE]
  	childPeaks <- rbind(childPeaksGood, childPeaksBad, childPeaksUnassigned)
  	
  	# Now let's cross fingers. Add a good=NA column to the unmatched peaks and reorder the columns
  	# to match order in childPeaks. After that, setData to the child slot.
  
  	childPeaksOmitted <- getData(child)
  	childPeaksOmitted <- childPeaksOmitted[child@low | child@satellite,,drop=FALSE]
  	childPeaksOmitted$rawOK <- rep(FALSE, nrow(childPeaksOmitted))
  	childPeaksOmitted$good <- rep(FALSE, nrow(childPeaksOmitted))
  	childPeaksOmitted$dppm <- rep(NA, nrow(childPeaksOmitted))
  	childPeaksOmitted$formula <- rep(NA, nrow(childPeaksOmitted))
  	childPeaksOmitted$mzCalc <- rep(NA, nrow(childPeaksOmitted))
  	childPeaksOmitted$dbe <- rep(NA, nrow(childPeaksOmitted))
    childPeaksOmitted$dppmBest <- rep(NA, nrow(childPeaksOmitted))
    childPeaksOmitted$formulaCount <- rep(0, nrow(childPeaksOmitted))
  	childPeaks$satellite <- rep(FALSE, nrow(childPeaks))
  	childPeaks$low <- rep(FALSE, nrow(childPeaks))
  	childPeaks$rawOK <- rep(TRUE, nrow(childPeaks))
	for(col in colnames(child@properties)) {
	  childPeaks[,col] <- child@properties[childPeaks$row, col]
	}
      
  	childPeaks <- childPeaks[,colnames(childPeaksOmitted), drop=FALSE]
  	
  	childPeaksTotal <- rbind(childPeaks, childPeaksOmitted)
	# add the formulaSource column
	childPeaksTotal$formulaSource <- character(nrow(childPeaksTotal))
	childPeaksTotal$formulaSource[!is.na(childPeaksTotal$formula)] <- "analyze"
  	child <- setData(child, childPeaksTotal)
  	child@ok <- TRUE
  	
  	return(child)
  }
  
  # I believe these lines were fixed to remove a warning but in the refactored workflow "mzranges" doesn't exist anymore.
  # Leave here for now 
  ## mzranges <- t(sapply(shots, function(p) {
  ##     if(!is.null(p$childRaw)){
  ##       return(range(p$childRaw[,mzColname]))
  ##     } else {
  ##       return(c(NA,NA))
  ##     }
  ## }))
  ## 
  ## mzmin <- min(mzranges[,1], na.rm=TRUE)
  ## mzmax <- max(mzranges[,2], na.rm=TRUE)
  children <- lapply(seq_along(msmsPeaks@children),
                     function(i) analyzeTandemShot(msmsPeaks@children[[i]], 
                                                   childIdx = i))
  
  ## correct fields in case of invalid data
  children <- lapply(children, function(child){
    if(child@ok)  return(child)
    child@rawOK        <- rep(x = TRUE,  times = child@peaksCount)
    child@good         <- rep(x = FALSE, times = child@peaksCount)
    child@mzCalc       <- as.numeric(  rep(x = NA, times = child@peaksCount))
    child@formula      <- as.character(rep(x = NA, times = child@peaksCount))
    child@dbe          <- as.numeric(  rep(x = NA, times = child@peaksCount))
    child@formulaCount <- as.integer(  rep(x = NA, times = child@peaksCount))
    child@dppm         <- as.numeric(  rep(x = NA, times = child@peaksCount))
    child@dppmBest     <- as.numeric(  rep(x = NA, times = child@peaksCount))
    return(child)
  })
  
## shots <- mapply(function(shot, scan, info)
  ##         {
  ##             shot$scan <- scan
  ##             shot$info <- info
  ##             shot$header <- msmsPeaks$childHeaders[as.character(scan),]
  ##             return(shot)
  ##         }, shots, msmsPeaks$childScans, spectraList, SIMPLIFY=FALSE)
  msmsPeaks@children <- as(children, "SimpleList")
  return(msmsPeaks)
}


#' @describeIn analyzeMsMs Analyze the peaks going only by intensity values
#' @export
analyzeMsMs.intensity <- function(msmsPeaks, mode="pH", detail=FALSE, run="preliminary",
		filterSettings = getOption("RMassBank")$filterSettings)
{
	cut <- 0
	cut_ratio <- 0
	if(run=="preliminary")
	{
		filterMode <- "coarse"
		cut <- filterSettings$prelimCut
		if(is.na(cut))
		{
		  adductProperties <- getAdductProperties(mode, msmsPeaks@formula)
		  if(adductProperties$charge > 0) cut <- 1e4
		  if(adductProperties$charge < 0) cut <- 0
		}
		cut_ratio <- filterSettings$prelimCutRatio
	}
	else
	{
		filterMode <- "fine"
		cut <- filterSettings$fineCut
		cut_ratio <- filterSettings$fineCutRatio
		if(is.na(cut)) cut <- 0
	}
	
	# find whole spectrum of parent peak, so we have reasonable data to feed into
	
	
	# On each spectrum the following function analyzeTandemShot will be applied.
	# It takes the raw peaks matrix as argument (mz, int) and processes the spectrum by
	# filtering out low-intensity (<1e4) and shoulder peaks (deltam/z < 0.5, intensity
	# < 5%) and subsequently matching the peaks to formulas using Rcdk, discarding peaks
	# with insufficient match accuracy or no match.
	analyzeTandemShot <- function(child)
	{
		shot <- getData(child)
		shot$row <- which(!is.na(shot$mz))
		
		# Filter out low intensity peaks:
		child@low <- (shot$intensity < cut) | (shot$intensity < max(shot$intensity)*cut_ratio)
		shot_full <- shot
        shot <- shot[!child@low,,drop=FALSE]
		
		
		# Is there still anything left?
		if(length(which(!child@low))==0)
		{
			child@ok <- FALSE
			return(child)
		}
		
		# Filter out satellite peaks:
		shot <- filterPeakSatellites(shot, filterSettings)
		child@satellite <- rep(TRUE, child@peaksCount)
		child@satellite[which(child@low == TRUE)] <- NA
		child@satellite[shot$row] <- FALSE
		
		# Is there still anything left?
		if(nrow(shot)==0)
		{
			child@ok <- FALSE
			return(child)
		}
		
		if(max(shot$intensity) < as.numeric(filterSettings$specOkLimit))
		{
			child@ok <- FALSE
			return(child)
		}
		

		# here follows the fake analysis
		#------------------------------------
		
		childPeaks <- addProperty(shot_full, "rawOK", "logical", FALSE)
        childPeaks[!(child@low | child@satellite),"rawOK"] <- TRUE
     
        childPeaks <- addProperty(childPeaks, "good", "logical", FALSE)
		childPeaks[childPeaks$rawOK,"good"] <- TRUE

		childPeaks <- addProperty(childPeaks, "mzCalc", "numeric")
		childPeaks[childPeaks$rawOK,"mzCalc"] <- childPeaks[childPeaks$rawOK,"mz"]
		
		childPeaks <- addProperty(childPeaks, "formula", "character")
		childPeaks[childPeaks$rawOK,"formula"] <- ""
		
		childPeaks <- addProperty(childPeaks, "dbe", "numeric")
		childPeaks[childPeaks$rawOK,"dbe"] <- 0
		
		childPeaks <- addProperty(childPeaks, "formulaCount", "integer")
		childPeaks[childPeaks$rawOK,"formulaCount"] <- 1
		
		childPeaks <- addProperty(childPeaks, "dppm", "numeric")
		childPeaks[childPeaks$rawOK,"dppm"] <- 0
		
		childPeaks <- addProperty(childPeaks, "dppmBest", "numeric")
		childPeaks[childPeaks$rawOK,"dppmBest"] <- 0
		
		child <- setData(child, childPeaks)
		child@ok <- TRUE
		return(child)
	}
	children <- lapply(msmsPeaks@children, analyzeTandemShot)
	msmsPeaks@children <- as(children, "SimpleList")
	#browser()

	return(msmsPeaks)

	# Omit all the stuff below for now, I don't believe it is needed. One thing is that spectraList info will have to be added somewhere else.
	## shots <- mapply(function(shot, scan, info)
	##         {
	##             shot$scan <- scan
	##             shot$info <- info
	##             shot$header <- msmsPeaks$childHeaders[as.character(scan),]
	##             return(shot)
	##         }, shots, msmsPeaks$childScans, spectraList, SIMPLIFY=FALSE)
	## 
	## mzranges <- t(sapply(shots, function(p) {return(range(p$childRaw[,mzColname]))}))
	## mzmin <- min(mzranges[,1], na.rm=TRUE)
	## mzmax <- max(mzranges[,2], na.rm=TRUE)
	## 
	## return(list(
	##                 msmsdata=shots,
	##                 mzrange=c(mzmin, mzmax),
	##                 id=msmsPeaks$id,
	##                 mode=mode,
	##                 parentHeader = msmsPeaks$parentHeader,
	##                 parentMs = msmsPeaks$parentPeak,
	##                 formula = msmsPeaks$formula,
	##                 foundOK = TRUE))
}


#' Filter peaks with low accuracy
#' 
#' Filters a peak table (with annotated formulas) for accuracy. Low-accuracy
#' peaks are removed.
#' 
#' In the \code{coarse} mode, mass tolerance is set to 10 ppm (above m/z 120)
#' and 15 ppm (below m/z 120). This is useful for formula assignment before
#' recalibration, where a wide window is desirable to accomodate the high mass
#' deviations at low m/z values, so we get a nice recalibration curve.
#' 
#' In the \code{fine} run, the mass tolerance is set to 5 ppm over the whole
#' mass range. This should be applied after recalibration.
#' 
#' @usage filterLowaccResults(peaks, mode="fine", filterSettings  = getOption("RMassBank")$filterSettings)
#' @param peaks A data frame with at least the columns \code{mzFound} and
#' \code{dppm}.
#' @param mode \code{coarse} or \code{fine}, see below.
#' @param filterSettings Settings for filtering. For details, see documentation of
#' 		\code{\link{analyzeMsMs}}
#' @return A \code{list(TRUE = goodPeakDataframe, FALSE = badPeakDataframe)} is
#' returned: A data frame with all peaks which are "good" is in
#' \code{return[["TRUE"]]}.
#' @author Michael Stravs
#' @seealso \code{\link{analyzeMsMs}}, \code{\link{filterPeakSatellites}}
#' @examples
#' 
#' # from analyzeMsMs:
#' \dontrun{childPeaksFilt <- filterLowaccResults(childPeaksInt, filterMode)}
#' 
#'
filterLowaccResults <- function(peaks, mode="fine", filterSettings  = getOption("RMassBank")$filterSettings)
{
  # Check if filter settings are properly set, otherwise use defaults 
  if(is.null(filterSettings))
  {
	  filterSettings <- list(
			ppmHighMass = 10,
	  		ppmLowMass = 15,
	  		massRangeDivision = 120,
	  		ppmFine = 5)
  }
	  
  peaks$good = NA
  peaks[!is.na(peaks$dppm), "good"] <- TRUE
  
  # coarse mode: to use for determinating the recalibration function
  if(mode=="coarse")
  {
    if(nrow(peaks[which(abs(peaks$dppm) > filterSettings$ppmHighMass),])>0)
    	peaks[which(abs(peaks$dppm) > filterSettings$ppmHighMass), "good"] <- FALSE
	if(nrow(peaks[which(peaks$mz > filterSettings$massRangeDivision & abs(peaks$dppm) > filterSettings$ppmLowMass),])>0)
    	peaks[which(peaks$mz > filterSettings$massRangeDivision & abs(peaks$dppm) > filterSettings$ppmLowMass), "good"] <- FALSE
  }
  # fine mode: for use after recalibration
  else
  {
	if(nrow(peaks[which(abs(peaks$dppm) > filterSettings$ppmFine),]) > 0)
    	peaks[which(abs(peaks$dppm) > filterSettings$ppmFine), "good"] <- FALSE
  }
  return(split(peaks, peaks$good))
}

#' Aggregate analyzed spectra
#' 
#' Groups an array of analyzed spectra and creates aggregated peak tables
#' 
#' \code{\var{addIncomplete}} is relevant for recalibration. For recalibration,
#' we want to use only high-confidence peaks, therefore we set
#' \code{\var{addIncomplete}} to \code{FALSE}. When we want to generate a peak
#' list for actually generating MassBank records, we want to include all peaks
#' into the peak tables.
#' 
#' @usage aggregateSpectra(spec,  addIncomplete=FALSE)
#' @param spec The \code{RmbSpectraSetList} of spectra sets (\code{RmbSpectraSet} objects) to aggregate
#' @param addIncomplete Whether or not the peaks from incomplete files (files
#' for which less than the maximal number of spectra are present)
#' @return 
#' A summary \code{data.frame} with all peaks (possibly multiple rows for one m/z value from a spectrum, see below) with columns:
#' \item{mzFound, intensity}{Mass and intensity of the peak}
#' \item{good}{if the peak passes filter criteria}
#' \item{mzCalc, formula, dbe, dppm}{calculated mass, formula, dbe and ppm deviation of the assigned formula}
#' \item{formulaCount, dppmBest}{Number of matched formulae for this m/z value, and ppm deviation of the best match}
#' \item{scan, cpdID, parentScan}{Scan number of the child and parent spectrum in the raw file, also the compound ID to which the peak belongs}
#' \item{dppmRc}{ppm deviation recalculated from the aggregation function}
#' \item{index}{Aggregate-table peak index, so the table can be subsetted, edited and results reinserted back into this table easily}
#' Further columns are later added by workflow steps 6 (electronic noise culler), 7 and 8.
#' 
#' @author Michael Stravs
#' @seealso \code{\link{msmsWorkflow}}, \code{\link{analyzeMsMs}}
#' @examples
#' 
#' ## As used in the workflow:
#' \dontrun{%
#' 	w@@spectra <- lapply(w@@spectra, function(spec)
#' 		analyzeMsMs(spec, mode="pH", detail=TRUE, run="recalibrated", cut=0, cut_ratio=0 ) )
#' 	w@@aggregate <- aggregateSpectra(w@@spectra)
#' }
#' 
#' @export
aggregateSpectra <- function(spec,  addIncomplete=FALSE)
{
	
	if(addIncomplete)
		aggSpectra <- selectSpectra(spec, "found", "object")
	else
		aggSpectra <- selectSpectra(spec, "complete", "object")
	
	compoundTables <- lapply(aggSpectra, function(s)
			{
				tables.c <- lapply(s@children, function(c)
						{
							table.c <- getData(c)
							table.c <- table.c[table.c$rawOK,,drop=FALSE]
							# remove superfluous columns, since only rawOK peaks are selected anyway
							table.c$rawOK <- NULL
							table.c$low <- NULL
							table.c$satellite <- NULL
							# add scan no
							table.c$scan <- rep(c@acquisitionNum, nrow(table.c))
							return(table.c)
						})
				table.cpd <- do.call(rbind, tables.c)
				
				## complete missing columns if necessary
				## mz intensity  good                                                                   scan cpdID parentScan
				## mz intensity  good    mzCalc    formula    dbe    formulaCount    dppm    dppmBest   scan cpdID parentScan
				columnNames <- c(       "mzCalc", "formula", "dbe", "formulaCount", "dppm", "dppmBest")
				if(all(!(columnNames %in% colnames(table.cpd))))
				  for(columnName in columnNames)
				    table.cpd[, columnName] <- as.numeric(rep(x = NA, times = nrow(table.cpd)))
				
				table.cpd$cpdID <- rep(s@id, nrow(table.cpd))
				table.cpd$name <- rep(s@name, nrow(table.cpd))
				table.cpd$parentScan <- rep(s@parent@acquisitionNum, nrow(table.cpd))
				return(table.cpd)
			})
	#return(compoundTables)
	aggTable <- do.call(rbind, compoundTables)
	if(is.null(aggTable)) aggTable <- data.frame("mz"=numeric(), "intensity"=numeric(), "good"=logical(), "mzCalc"=numeric(), "formula"=character(), "dbe"=numeric(), "formulaCount"=integer(), "dppm"=numeric(), "dppmBest"=numeric(), "scan"=integer(), "cpdID"=integer(), "parentScan"=integer(), stringsAsFactors=FALSE)
	colnames(aggTable)[1] <- "mzFound"

	aggTable <- addProperty(aggTable, "dppmRc", "numeric")
	aggTable <- addProperty(aggTable, "index", "integer")
	if(nrow(aggTable) > 0)
		aggTable$index <- 1:nrow(aggTable)
	
	aggTable[aggTable$good, "dppmRc"] <- (aggTable[aggTable$good, "mzFound"]/aggTable[aggTable$good, "mzCalc"] - 1)*1e6
	
	return(aggTable)
}

#' Remove electronic noise
#' 
#' Removes known electronic noise peaks from a peak table
#' 
#' @usage cleanElnoise(peaks, noise=getOption("RMassBank")$electronicNoise,
#' 		width = getOption("RMassBank")$electronicNoiseWidth)
#' @param peaks An aggregated peak frame as described in \code{\link{aggregateSpectra}}. Columns
#' \code{mzFound}, \code{dppm} and \code{dppmBest} are needed.
#' @param noise A numeric vector of known m/z of electronic noise peaks from the instrument
#' 		Defaults to the	entries in the RMassBank settings.
#' @param width The window for the noise peak in m/z units. Defaults to the entries in 
#' 		the RMassBank settings.
#' @return Extends the aggregate data frame by column \code{noise} (logical), which is \code{TRUE} if the peak is marked as noise.
#' 
#' @author Michael Stravs
#' @seealso \code{\link{msmsWorkflow}}
#' @examples
#' # As used in the workflow:
#' \dontrun{
#' 	    w@@aggregated <- 
#' 		cleanElnoise(w@@aggregated)	
#' }
#' @export
setGeneric("cleanElnoise",	
		signature=c("peaks", "noise", "width"),
		function(peaks, 
				noise = getOption("RMassBank")$electronicNoise,
				width = getOption("RMassBank")$electronicNoiseWidth
				)
				standardGeneric("cleanElnoise"))



.cleanElnoise.df <- function(peaks, noise=getOption("RMassBank")$electronicNoise,
		width = getOption("RMassBank")$electronicNoiseWidth)
{
	peaks <- addProperty(peaks, "noise", "logical", FALSE)
	
	# I don't think this makes sense if using one big table...
	## # use only best peaks
	## p_best <- peaks[is.na(peaks$dppmBest) | (peaks$dppm == peaks$dppmBest),]
  
  # remove known electronic noise
  p_eln <- peaks
  for(noisePeak in noise)
  {
		noiseMatches <- which(!((p_eln$mzFound > noisePeak + width)	| (p_eln$mzFound < noisePeak - width)))
		if(length(noiseMatches) > 0)
			p_eln[noiseMatches, "noise"] <- TRUE
  }
  return(p_eln)
}

#' @export
setMethod("cleanElnoise", signature(peaks="data.frame", noise="numeric", width="numeric"),
		function(peaks, noise, width) .cleanElnoise.df(peaks, noise, width))


#' @export
setMethod("cleanElnoise", signature(peaks="RmbSpectrum2", noise="numeric", width="numeric"),
		function(peaks, noise, width) 
		{
		 se <- peaks
		 data <- getData(se)
		 data$mzFound <- data$mz # this is a bit of a temporary workaround
		 se <- addProperty(se, "noise", "logical", FALSE)
		 data <- cleanElnoise(data, noise, width)
		 se <- setData(se, data)
		 se
		})

#' @export
#' @describeIn selectPeaks A method to filter spectra to the specified peaks
setMethod("cleanElnoise", c("RmbSpectrum2List",noise="numeric", width="numeric"), function(peaks, noise, width)
		{
			s <- lapply(peaks, function(s) cleanElnoise(s, noise, width))
			for(n in seq_len(length(peaks)))
				peaks[[n]] <- s[[n]]
			return(peaks)
		})

setMethod("cleanElnoise", c("RmbSpectraSet", noise="numeric", width="numeric"), function(peaks, noise, width)
		{
      if(length(peaks@children) == 0)
        return(peaks)
			children <- lapply(peaks@children, cleanElnoise, noise=noise, width=width)
			peaks@children <- as(children, "SimpleList")
			return(peaks)
		})
#' Identify intense peaks (in a list of unmatched peaks)
#' 
#' Finds a list of peaks in spectra with a high relative intensity (>10% and
#' 1e4, or >1% and 1e5) to write a list of peaks which must be manually
#' checked.  Peaks orbiting around the parent peak mass (calculated from the
#' compound ID), which are very likely co-isolated substances, are ignored.
#' 
#' 
#' @usage problematicPeaks(peaks_unmatched, peaks_matched, mode = "pH")
#' @param peaks_unmatched Table of unmatched peaks, with at least \code{cpdID,
#' scan, mzFound, int}.
#' @param peaks_matched Table of matched peaks (used for base peak reference),
#' with at least \code{cpdID, scan, int}.
#' @param mode Processing mode (\code{"pH", "pNa"} etc.)
#' @return A filtered table with the potentially problematic peaks, including
#' the precursor mass and MSMS base peak intensity (\code{aMax}) for reference.
#' @author Michael Stravs
#' @seealso \code{\link{msmsWorkflow}}
#' @examples \dontrun{
#' # As used in the workflow: 
#' fp <- problematicPeaks(specs[!specs$filterOK & !specs$noise & 
#' 						((specs$dppm == specs$dppmBest) | (is.na(specs$dppmBest)))
#' 				,,drop=FALSE], peaksMatched(w), mode)
#' }
#' @export
problematicPeaks <- function(sp)
{
  if(!sp@ok)
    return(sp)
  sp <- addProperty(sp, "problematicPeak", "logical", FALSE)
  sp <- addProperty(sp, "aMax", "numeric", 0)
  sp <- addProperty(sp, "mzCenter", "numeric", sp@precursorMz)
  peaks <- getData(sp)
  peaks$`_index` <- seq_len(nrow(peaks))
  peaks$aMax <- max(peaksMatched(peaks)$intensity)
  
  peaks_unmatched <- peaks[!peaks$good,,drop=FALSE]
  #peaks_unmatched$aMax <- aMax
  p_control <- peaks_unmatched[
    ( (peaks_unmatched$intensity > 1e5) &
        (peaks_unmatched$intensity > 0.01*peaks_unmatched$aMax)) 
    | ( (peaks_unmatched$intensity > 1e4) &
          (peaks_unmatched$intensity > 0.1* peaks_unmatched$aMax)) ,]
  
  # use parent m/z to exclude co-isolated peaks
  #p_control$mzCenter <- numeric(nrow(p_control))
  p_control_noMH <- p_control[
		  (p_control$mzFound < p_control$mzCenter - 1) |
		  (p_control$mzFound > p_control$mzCenter + 1),,drop=FALSE]
  if(nrow(p_control_noMH) > 0)
    peaks[peaks$`_index` %in% p_control_noMH$`_index`]$problematicPeak <- TRUE
  sp <- setData(sp, peaks)
  return(sp)
}

#' Generate list of problematic peaks
#' 
#' Generates a list of intense unmatched peaks for further review (the "failpeak list") and exports it if the archive name is given.
#' 
#' @param w \code{msmsWorkspace} to analyze. 
#' @param mode Processing mode (pH etc)
#' @param archivename Base name of the archive to write to (for "abc" the exported failpeaks list will be "abc_Failpeaks.csv").
#' if the compoundlist is complete, "tentative", if at least a formula is present or "unknown"
#' if the only know thing is the m/z
#' @return  Returns the aggregate data.frame with added column "\code{problematic}" (logical) which marks peaks which match the problematic criteria 
#' 
#' @author stravsmi
#' @export
processProblematicPeaks <- function(w, archivename = NA)
{
	w@spectra <- as(lapply(w@spectra, function(cpd) {
	  if(length(cpd@children) == 0)
	    return(cpd)
	  cpd@children <- as(
	    lapply(cpd@children, problematicPeaks),
	    "SimpleList")
	  cpd
	}), "SimpleList")
	
	w@aggregated <- aggregateSpectra(w, addIncomplete = TRUE)
	
	fp <- w@aggregated[w@aggregated$problematicPeak,,drop=FALSE]
	fp$OK <- rep(FALSE, nrow(fp))
	
	# Select the columns for output into the failpeaks file
	fp <- fp[,c("OK", "name", "cpdID", "scan", "mzFound", "formula", 
	            "mzCalc", "dppm", "dbe", "intensity",
	            "formulaCount", "parentScan", "aMax", "mzCenter")]		
	if(!is.na(archivename))
	  write.csv(fp, file=
	              paste(archivename,"_Failpeaks.csv", sep=''), row.names=FALSE)
	
	return(w)
	# TODO: Handle combineMultiplicities
	
	# Select the correct precursor scans. This serves to filter the list
	# for the cases where multiple workspaces were combined after step 7
	# with combineMultiplicities.
	# Note that this has drawbacks. Leaving the "duplicates" in would make it more easy
	# to identify legitimate unformulaed peaks. We might experiment by marking them up
	# somehow. 
	precursors <- unlist(lapply(selectSpectra(w, "found", "object"), function(s) s@parent@acquisitionNum))
	fp <- fp[
			fp$parentScan %in% precursors
			,]
	
}

#' Recalibrate MS/MS spectra
#' 
#' Recalibrates MS/MS spectra by building a recalibration curve of the
#' assigned putative fragments of all spectra in \code{aggregatedSpecs}
#' (measured mass vs. mass of putative associated fragment) and additionally
#' the parent ion peaks. 
#' 
#' Note that the actually used recalibration functions are governed by the
#' general MassBank settings (see \code{\link{recalibrate}}).
#' 
#' If a set of acquired LC-MS runs contains spectra for two different ion types
#' (e.g. [M+H]+ and [M+Na]+) which should both be processed by RMassBank, it is
#' necessary to do this in two separate runs. Since it is likely that one ion type
#' will be the vast majority of spectra (e.g. most in [M+H]+ mode), and only few
#' spectra will be present for other specific adducts (e.g. only few [M+Na]+ spectra),
#' it is possible that too few spectra are present to build a good recalibration curve
#' using only e.g. the [M+Na]+ ions. Therefore we recommend, for one set of LC/MS runs,
#' to build the recalibration curve for one ion type 
#' (\code{msmsWorkflow(mode="pH", steps=c(1:8), newRecalibration=TRUE)})
#' and reuse the same curve for processing different ion types 
#' (\code{msmsWorkflow(mode="pNa", steps=c(1:8), newRecalibration=FALSE)}).
#' This also ensures a consistent recalibration across all spectra of the same batch. 
#' 
#' @usage makeRecalibration(w, mode, 
#'  	recalibrateBy = getOption("RMassBank")$recalibrateBy,
#' 		recalibrateMS1 = getOption("RMassBank")$recalibrateMS1,
#' 		recalibrator = getOption("RMassBank")$recalibrator,
#' 		recalibrateMS1Window = getOption("RMassBank")$recalibrateMS1Window 
#' 		)
#' 
#'  recalibrateSpectra(mode, rawspec = NULL, rc = NULL, rc.ms1=NULL, w = NULL,
#' 		recalibrateBy = getOption("RMassBank")$recalibrateBy,
#' 		recalibrateMS1 = getOption("RMassBank")$recalibrateMS1)
#' 
#'  recalibrateSingleSpec(spectrum, rc, 
#' 		recalibrateBy = getOption("RMassBank")$recalibrateBy)
#' @aliases makeRecalibration recalibrateSpectra recalibrateSingleSpec
#' @param w For \code{makeRecalibration}: to perform the recalibration with. For \code{recalibrateSpectra}: 
#' 			the \code{msmsWorkspace} which contains the recalibration curves (alternatively to specifying \code{rc, rc.ms1}). 
#' @param spectrum For \code{recalibrateSingleSpec}:
#' 			a \code{MSnbase} \code{Spectrum}-derived object, commonly a \code{RmbSpectrum2} for MS2 or \code{Spectrum1} for MS1.
#' @param mode \code{"pH", "pNa", "pM", "mH", "mM", "mFA"} for different ions 
#' 			([M+H]+, [M+Na]+, [M]+, [M-H]-, [M]-, [M+FA]-).
#' @param rawspec For \code{recalibrateSpectra}:an \code{RmbSpectraSetList} of \code{RmbSpectraSet} objects
#' 			, as the \code{w@@spectra} slot from \code{msmsWorkspace} or any object returned by \code{\link{findMsMsHR}}.
#' 			If empty, no spectra are recalibrated, but the recalibration curve is
#' 			returned.  
#' @param rc,rc.ms1 The recalibration curves to be used in the recalibration.
#' @param recalibrateBy Whether recalibration should be done by ppm ("ppm") or by m/z ("mz").
#' @param recalibrateMS1 Whether MS1 spectra should be recalibrated separately ("separate"),
#' 		together with MS2 ("common") or not at all ("none"). Usually taken from settings.
#' @param recalibrator The recalibrator functions to be used.
#' 		 Refer to \code{\link{recalibrate}} for details. Usually taken from settings.
#' @param recalibrateMS1Window Window width to look for MS1 peaks to recalibrate (in ppm).
#' @return \code{makeRecalibration}: a \code{list(rc, rc.ms1)} with recalibration curves
#' 			for the MS2 and MS1 spectra.
#' 
#' 			\code{recalibrateSpectra}: if \code{rawspec} is not \code{NULL}, returns the recalibrated
#' 			spectra as \code{RmbSpectraSetList}. All spectra have their mass recalibrated and evaluation data deleted.
#' 
#' 			\code{recalibrateSingleSpec}: the recalibrated \code{Spectrum} (same object, recalibrated masses,
#' 				 evaluation data like assigned formulae etc. deleted). 
#' 
#' @examples \dontrun{ 
#' 			rcCurve <- recalibrateSpectra(w, "pH")
#' 			w@@spectra <- recalibrateSpectra(mode="pH", rawspec=w@@spectra, w=myWorkspace)
#' 			w@@spectra <- recalibrateSpectra(mode="pH", rawspec=w@@spectra,	rcCurve$rc, rcCurve$rc.ms1)
#' 			}
#' 
#' @author Michael Stravs, Eawag <michael.stravs@@eawag.ch>
#' @export 
makeRecalibration <- function(w, mode, 
		recalibrateBy = getOption("RMassBank")$recalibrateBy,
		recalibrateMS1 = getOption("RMassBank")$recalibrateMS1,
		recalibrator = getOption("RMassBank")$recalibrator,
		recalibrateMS1Window = getOption("RMassBank")$recalibrateMS1Window
		#cropSingleFormula = getOption("RMassBank")$recalibrateCropSingleFormula
		)
{
	if(is.null(w@spectra))
		stop("No spectra present to generate recalibration curve.")

	rcdata <- peaksMatched(w)
	#if(cropSingleFormula)
	  rcdata <- rcdata[!is.na(rcdata$formulaCount) & rcdata$formulaCount == 1, ,drop=FALSE]
	#else
	#  rcdata <- rcdata[!is.na(rcdata$formulaCount), ,drop=FALSE]
	
	rcdata <- rcdata[,c("mzFound", "dppm", "mzCalc")]
	
	if(nrow(rcdata) == 0)
		stop("No peaks matched to generate recalibration curve.")
	
	ms1data <- recalibrate.addMS1data(w@spectra, mode, recalibrateMS1Window)
	ms1data <- ms1data[,c("mzFound", "dppm", "mzCalc")]
  
	if (recalibrateMS1 != "none") {
          ## Add m/z values from MS1 to calibration datapoints
          rcdata <- rbind(rcdata, ms1data)
        }
        
	rcdata$dmz <- rcdata$mzFound - rcdata$mzCalc
	ms1data$dmz <- ms1data$mzFound - ms1data$mzCalc
	
	if(recalibrateBy == "dppm")
	{
		rcdata$recalfield <- rcdata$dppm
		ms1data$recalfield <- ms1data$dppm
	}
	else
	{
		rcdata$recalfield <- rcdata$dmz
		ms1data$recalfield <- ms1data$dmz
	}
	
	# generate recalibration model
	rc <- do.call(recalibrator$MS2, list(rcdata)) 
	if(recalibrateMS1 == "separate")
		rc.ms1 <- do.call(recalibrator$MS1, list(ms1data)) 
	else
		rc.ms1 <- rc
	
	# plot the model
	par(mfrow=c(2,2))
	if(nrow(rcdata)>0)
		plotRecalibration.direct(rcdata = rcdata, rc = rc, rc.ms1 = rc.ms1, title = "MS2", 
				mzrange = range(rcdata$mzFound),
				recalibrateBy = recalibrateBy)
	if(nrow(ms1data)>0)
		plotRecalibration.direct(ms1data, rc, rc.ms1, "MS1",
				range(ms1data$mzFound),
				recalibrateBy)
	# Return the computed recalibration curves
	return(list(rc=rc, rc.ms1=rc.ms1))
}



#' Plot the recalibration graph.
#' 
#' @aliases plotRecalibration plotRecalibration.direct
#' @usage plotRecalibration(w, recalibrateBy = getOption("RMassBank")$recalibrateBy)
#' 
#' 		plotRecalibration.direct(rcdata, rc, rc.ms1, title, mzrange,
#' 		recalibrateBy = getOption("RMassBank")$recalibrateBy)
#' 
#' @param w The workspace to plot the calibration graph from
#' @param rcdata A data frame with columns \code{recalfield} and \code{mzFound}.
#' @param rc Predictor for MS2 data
#' @param rc.ms1 Predictor for MS1 data
#' @param title Prefix for the graph titles
#' @param mzrange m/z value range for the graph
#' @param recalibrateBy Whether recalibration was done by ppm ("ppm") or by m/z ("mz").
#' 		Important only for graph labeling here.
#' 
#' @author Michele Stravs, Eawag <michael.stravs@@eawag.ch>
#' @export
plotRecalibration <- function(w, recalibrateBy = getOption("RMassBank")$recalibrateBy)
{
	spec <- w@aggregated
	if(!is.null(w@parent))
		spec <- w@parent@aggregated

	rcdata <- data.frame(mzFound = w@rc$x, recalfield = w@rc$y)
	ms1data <- data.frame(mzFound = w@rc.ms1$x, recalfield = w@rc.ms1$y)
	
	
	
	par(mfrow=c(2,2))
	if(nrow(rcdata)>0)
		plotRecalibration.direct(rcdata, w@rc, w@rc.ms1, "MS2", 
				range(spec$mzFound[which(spec$good)]),recalibrateBy)
	if(nrow(ms1data)>0)
		plotRecalibration.direct(ms1data, w@rc, w@rc.ms1, "MS1",
				range(ms1data$mzFound),recalibrateBy)
	
}

#' @export 
plotRecalibration.direct <- function(rcdata, rc, rc.ms1, title, mzrange,
		recalibrateBy = getOption("RMassBank")$recalibrateBy
		)
{
	if(recalibrateBy == "dppm")
		ylab.plot <- expression(paste(delta, "ppm"))
	else
		ylab.plot <- expression(paste(delta, "m/z"))	
	
	plot(recalfield ~ mzFound, data=rcdata,
			xlab = "m/z", ylab = ylab.plot, main=paste(title, "scatterplot"))
	RcModelMz <- seq(mzrange[[1]], mzrange[[2]], by=0.2)
	RcModelRecal <- predict(rc, newdata= data.frame(mzFound =RcModelMz))
	RcModelRecalMs1 <- predict(rc.ms1, newdata= data.frame(mzFound =RcModelMz))
	lines(RcModelMz, RcModelRecal, col="blue")
	lines(RcModelMz, RcModelRecalMs1, col="yellow")
	if((length(unique(rcdata$mzFound))>1) & 
			(length(unique(rcdata$recalfield))>1))
	{
		if(requireNamespace("gplots",quietly=TRUE))
		{
			
			gplots::hist2d(rcdata$mzFound, rcdata$recalfield, 
					col=c("white", heat.colors(12)), xlab="m/z", 
					ylab = ylab.plot, main=paste(title, "density"))
			lines(RcModelMz, RcModelRecal, col="blue")
			lines(RcModelMz, RcModelRecalMs1, col="yellow")
		}
		else
		{
			rmb_log_info("Package gplots not installed. The recalibration density plot will not be displayed.")
			rmb_log_info("To install gplots: install.packages('gplots')")
		}
	}
}


#' @export
recalibrateSpectra <- function(mode, rawspec = NULL, rc = NULL, rc.ms1=NULL, w = NULL,
		recalibrateBy = getOption("RMassBank")$recalibrateBy,
		recalibrateMS1 = getOption("RMassBank")$recalibrateMS1)
{
	# Load the recal curves from the workspace if one is specified.
  if(!is.null(w))
  {
	  rc <- w@rc
	  rc.ms1 <- w@rc.ms1
  }
  if(is.null(rc) || is.null(rc.ms1))
	  stop("Please specify the recalibration curves either via workspace (w) or via parameters rc, rc.ms1.")

  # Do the recalibration
  if(!is.null(rawspec))
  {
	  # go through all raw spectra and recalculate m/z values
	  recalibratedSpecs <- lapply(rawspec, function(s)
			  {
				  if(s@found)
				  {
					  # recalculate tandem spectrum peaks
					  recalSpectra <- lapply(s@children, function(p)
							  {
								  recalibrateSingleSpec(p, rc, recalibrateBy)
							  })
					  s@children <- as(recalSpectra, "SimpleList")
					  # recalculate MS1 spectrum if required
					  if(recalibrateMS1 != "none")
					  {
						  s@parent <- recalibrateSingleSpec(s@parent, rc.ms1, recalibrateBy)
					  }
				  }
				  s@empty <- NA
				  s@complete <- NA
				  return(s)
			  } )
	  return(as(recalibratedSpecs, "SimpleList"))
  }
  else # no rawspec passed
	  return(list())
}

#' @export
recalibrateSingleSpec <- function(spectrum, rc, 
		recalibrateBy = getOption("RMassBank")$recalibrateBy)
{
	spectrum.df <- as.data.frame(as(spectrum, "Spectrum"))
	spectrum.df <- spectrum.df[!duplicated(spectrum.df$mz),,drop=FALSE]
	spectrum.df <- spectrum.df[order(spectrum.df$mz),,drop=FALSE]
	
	mzVals <- spectrum.df
	if(nrow(mzVals) > 0)
	{
		# Fix the column names so our
		# prediction functions choose the right
		# rows. 
		colnames(mzVals) <- c("mzFound", "int")
		drecal <- predict(rc, newdata=mzVals)
		if(recalibrateBy == "dppm")
			mzRecal <- mzVals$mzFound / (1 + drecal/1e6)
		else
			mzRecal <- mzVals$mzFound - drecal
		# And rename them back so our "mz" column is
		# called "mz" again
	}
	spectrum.df$mzRaw <- spectrum.df$mz
	spectrum.df$mz <- mzRecal
	
	
	# now comes the part that I don't like too much; this could be improved by using as.data.frame instead of getData and correspondingly
	# also not use setData. For now I leave it like this.
	# The problem is that I am not sure whether the default behaviour of as.RmbSpectrum2 should be clean=TRUE or FALSE,
	# and vice versa, I am not sure if as.data.frame should return only mz/int or the whole table.
	
	if(is(spectrum, "RmbSpectrum2"))
	{
		# this removes all evaluated data that were added in step 2 except for @ok I think
		colnames(spectrum.df) <- c("mz", "intensity", "mzRaw")
		spectrum <- addProperty(spectrum, "mzRaw", "numeric", NA)
		spectrum <- setData(spectrum, spectrum.df, clean=TRUE)
		# It also avoids making a new object when we don't know what class it should be 
	}
	else
	{
		# for Spectrum1 or all others that we don't know
		spectrum@mz <- spectrum.df$mz
		spectrum@intensity <- spectrum.df$i
	}
		
	return(spectrum)
}





#' Filter satellite peaks
#' 
#' Filters satellite peaks in FT spectra which arise from FT artifacts and from
#' conversion to stick mode. A very simple rule is used which holds mostly true
#' for MSMS spectra (and shouldn't be applied to MS1 spectra which contain
#' isotope structures...)
#' 
#' The function cuts off all peaks within 0.5 m/z from every peak, in
#' decreasing intensity order, which are below 5% of the referring peak's
#' intensity.  E.g. for peaks m/z=100, int=100; m/z=100.2, int=2, m/z=100.3,
#' int=6, m/z 150, int=10: The most intense peak (m/z=100) is selected, all
#' neighborhood peaks below 5% are removed (in this case, only the m/z=100.2
#' peak) and the next less intense peak is selected. Here this is the m/z=150
#' peak. All low-intensity neighborhood peaks are removed (nothing). The next
#' less intense peak is selected (m/z=100.3) and again neighborhood peaks are
#' cut away (nothing to cut here. Note that the m/z = 100.2 peak was alredy
#' removed.)
#' 
#' @usage filterPeakSatellites(peaks, filterSettings = getOption("RMassBank")$filterSettings)
#' @param peaks A peak dataframe with at least the columns \code{mz, int}. Note
#' that \code{mz} is used even for the recalibrated spectra, i.e. the
#' desatellited spectrum is identical for both the unrecalibrated and the
#' recalibrated spectra.
#' @param filterSettings The settings used for filtering. Refer to \code{\link{analyzeMsMs}}
#' 		documentation for filter settings.
#' @return Returns the peak table with satellite peaks removed.
#' @note This is a very crude rule, but works remarkably well for our spectra.
#' @author Michael Stravs
#' @seealso \code{\link{analyzeMsMs}}, \code{\link{filterLowaccResults}}
#' @examples
#' 
#' # From the workflow:
#' \dontrun{
#'     # Filter out satellite peaks:
#'     shot <- filterPeakSatellites(shot)
#'     shot_satellite_n <- setdiff(row.names(shot_full), row.names(shot))
#'     shot_satellite <- shot_full[shot_satellite_n,]
#'     # shot_satellite contains the peaks which were eliminated as satellites.
#' }
#' 
#' @export
filterPeakSatellites <- function(peaks, filterSettings = getOption("RMassBank")$filterSettings)
{
	mergePeaks(peaks, 0, filterSettings$satelliteMzLimit, filterSettings$satelliteIntLimit)
}

#' Reanalyze unmatched peaks
#' 
#' Reanalysis of peaks with no matching molecular formula by allowing
#' additional elements (e.g. "N2O").
#' 
#' \code{reanalyzeFailpeaks} examines the \code{unmatchedPeaksC} table in
#' \code{specs} and sends every peak through \code{reanalyzeFailpeak}.
#' 
#' @aliases reanalyzeFailpeaks reanalyzeFailpeak
#' @usage reanalyzeFailpeaks(aggregated, custom_additions, mode, filterSettings =
#' 				getOption("RMassBank")$filterSettings, progressbar = "progressBarHook")
#' reanalyzeFailpeak(custom_additions, mass, cpdID, counter, pb = NULL, mode,
#' 				filterSettings = getOption("RMassBank")$filterSettings)
#' @param aggregated A peake aggregate table (\code{w@@aggregate}) (after processing electronic noise removal!)
#' @param custom_additions The allowed additions, e.g. "N2O".
#' @param mode Processing mode (\code{"pH", "pNa", "mH"} etc.)
#' @param mass (Usually recalibrated) m/z value of the peak.
#' @param cpdID Compound ID of this spectrum.
#' @param counter Current peak index (used exclusively for the progress
#' indicator)
#' @param pb A progressbar object to display progress on, as passed by
#'  \code{reanalyzeFailpeaks} to \code{reanalyzeFailpeak}. No progress 
#' is displayed if NULL.
#' @param progressbar The progress bar callback to use. Only needed for specialized
#'  applications.	Cf. the documentation of \code{\link{progressBarHook}} for usage.
#' @param filterSettings Settings for filtering data. Refer to\code{\link{analyzeMsMs}} for settings.
#' @return The aggregate data frame extended by the columns:
#' #' \item{reanalyzed.???}{If reanalysis (step 7) has already been processed: matching values from the reanalyzed peaks}
#' \item{matchedReanalysis}{Whether reanalysis has matched (\code{TRUE}), not matched(\code{FALSE}) or has not been conducted for the peak(\code{NA}).}
#' 
#' It would be good to merge the analysis functions of \code{analyzeMsMs} with
#' the one used here, to simplify code changes.
#' @author Michael Stravs
#' @seealso \code{\link{analyzeMsMs}}, \code{\link{msmsWorkflow}}
#' @examples
#' 
#' ## As used in the workflow:
#' \dontrun{    
#' 	reanalyzedRcSpecs <- reanalyzeFailpeaks(w@@aggregated, custom_additions="N2O", mode="pH")
#' # A single peak:
#' reanalyzeFailpeak("N2O", 105.0447, 1234, 1, 1, "pH")
#' }
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' @export
reanalyzeFailpeaks <- function(w, custom_additions, mode, filterSettings =
				getOption("RMassBank")$filterSettings, progressbar = "progressBarHook")
{
	nProg <- 0
	nLen <- length(w@spectra)
	pb <- do.call(progressbar, list(object = NULL, value = 0, 
					min = 0, max = nLen))
	spectra <- lapply(w@spectra, function(sp) {
				if(length(sp@children) == 0)
					return(sp)
				
				children <- lapply(sp@children, function(ch) {
				      if(!ch@ok)
				        return(ch)
							peaks <- getData(ch)
							# get the peaks that have no matching formula, but are considered not noise etc.
							mzGood <- peaks$mz[peaks$good]
							peaksReanalyze <- !(peaks$mz %in% mzGood)
							peaksReanalyze <- peaksReanalyze & peaks$rawOK
							peaks.sp <- split(peaks, !peaksReanalyze)
							if (!is.null(peaks.sp[["FALSE"]])) {
								fp <- peaks.sp[["FALSE"]]
								fp <- fp[!duplicated(fp$mz),,drop=FALSE]
								peaks.rean <- lapply(fp$mz, reanalyzeFailpeak, 
										custom_additions = custom_additions, cpdID = sp@id, 
										mode = mode, filterSettings = filterSettings)
								matched <- (unlist(lapply(peaks.rean, nrow))) > 
										0
								df.rean <- do.call(rbind, peaks.rean[matched])
								fp[matched, colnames(df.rean)] <- df.rean
								fp[matched, "formulaSource"] <- rep("reanalysis", 
										sum(matched))
								peaks.sp[["FALSE"]] <- fp
								peaks <- do.call(rbind, peaks.sp)
								ch <- setData(ch, peaks)
							}
							ch
						})
				sp@children <- as(children, "SimpleList")
				nProg <<- nProg + 1
				pb <- do.call(progressbar, list(object = pb, value = nProg))
				return(sp)
			})
	w@spectra <- as(spectra, "SimpleList")
	do.call(progressbar, list(object = pb, close = TRUE))
	return(w)
}


#' @export
reanalyzeFailpeak <- function(mass, custom_additions, cpdID, mode,
		filterSettings = getOption("RMassBank")$filterSettings)
{
	# here follows the Rcdk analysis
	#------------------------------------
	
	# the ppm range is two-sided here.
	# The range is slightly expanded because dppm calculation of
	# generate.formula starts from empirical mass, but dppm cal-
	# culation of the evaluation starts from theoretical mass.
	# So we don't miss the points on 'the border'.
	
	db_formula <- findFormula(cpdID, retrieval=findLevel(cpdID,TRUE))
	
	# get the adduct additions
	adductProperties <- getAdductProperties(mode, db_formula)
	allowed_additions <- adductProperties$addition
	mode.charge <- adductProperties$charge
	
	ppmlimit <- 2.25 * filterSettings$ppmFine
	parent_formula <- add.formula(db_formula, allowed_additions)
	parent_formula <- add.formula(parent_formula, custom_additions)
	dbe_parent <- dbe(parent_formula)
	# check whether the formula is valid, i.e. has no negative or zero element numbers.
	#print(parent_formula)
	limits <- to.limits.rcdk(parent_formula)        
	
	peakformula <- tryCatch(
	  suppressWarnings(generate.formula(mass.calc, ppm(mass.calc, ppmlimit, p=TRUE),
	                                    limits, charge=0)),
	  error = function(e) list())
	
	if(length(peakformula)==0)
		return(data.frame())
	else # if is.list(peakformula)
	# formula found? then return the one with lowest dppm
	{
		# calculate dppm for all formulas
		peakformula <- sapply(peakformula, function(f)
				{
					l <- list(mz=as.numeric(as.character(mass)),
							formula=as.character(f@string),
							mzCalc=as.numeric(as.character(f@mass))
					)
					
					return(unlist(l))
				})
		
		# filter out bad dbe stuff
		peakformula <- as.data.frame(t(peakformula))
		# for some reason completely oblivious to me, the columns in peakformula
		# are still factors, even though i de-factored them by hand.
		# Therefore, convert them again...
		peakformula$mz <- as.numeric(as.character(peakformula$mz))
		peakformula$formula <- as.character(peakformula$formula)
		peakformula$mzCalc <- as.numeric(as.character(peakformula$mzCalc))
		
		peakformula$dppm <- (peakformula$mz / peakformula$mzCalc - 1) * 1e6
		peakformula$formulaCount=nrow(peakformula)
		
		
		# filter out bad dbe and high ppm stuff          
		peakformula$dbe <- unlist(lapply(peakformula$formula, dbe))
		peakformula <- peakformula[(peakformula$dbe >= filterSettings$dbeMinLimit) 
						& (abs(peakformula$dppm) < filterSettings$ppmFine),]
		# is there still something left?
		if(nrow(peakformula) == 0)
			return(data.frame())
		else
		{
			#update formula count to the remaining formulas
			peakformula$good = TRUE
			peakformula$formulaCount=nrow(peakformula)
			return(peakformula[which.min(abs(peakformula$dppm)),])
		}
		
	} # endif is.list(peakformula)
      

    
    }

#' Multiplicity filtering: Removes peaks which occur only once in a n-spectra set.
#' 
#' For every compound, every peak (with annotated formula) is compared 
#' across all spectra. Peaks whose formula occurs only once for all collision energies
#' / spectra types, are discarded. This eliminates "stochastic formula hits" of pure
#' electronic noise peaks efficiently from the spectra. Note that in the author's 
#' experimental setup two spectra were recorded at every collision energy,
#' and therefore every peak-formula should appear
#' at least twice if it is real, even if it is by chance a fragment which appears
#' on only one collision energy setting. The function was not tested in a different
#' setup. Therefore, use with a bit of caution.
#' @usage filterPeaksMultiplicity(peaks, formulacol, recalcBest = TRUE)
#' @param peaks An aggregate peak data.frame containing all peaks to be analyzed; with at least
#' 			the columns \code{cpdID, scan, mzFound} and one column for the formula
#' 			specified with the \code{formulacol} parameter. 
#' @param formulacol Which column the assigned formula is stored in. (Needed to separately process \code{"formula"} and
#' 			\code{"reanalyzed.formula"} multiplicites.)
#' @param recalcBest Whether the best formula for each peak should be re-determined.
#' 			This is necessary for results from the ordinary \code{\link{analyzeMsMs}}
#' 			analysis which allows multiple potential formulas per peak - the old best match
#' 			could potentially have been dropped because of multiplicity filtering. For results
#' 			from \code{\link{reanalyzeFailpeak}} this is not necessary, since only one potential
#' 			formula is assigned in this case.
#' @return The peak table is returned, enriched with columns:
#' 			\itemize{
#' 				\item{\code{formulaMultiplicity}}{The # of occurrences of this formula
#' 					in the spectra of its compounds.}
#' 			}
#' @examples \dontrun{
#' 		peaksFiltered <- filterPeaksMultiplicity(peaksMatched(w), 
#' 			"formula", TRUE)
#' 		peaksOK <- subset(peaksFiltered, formulaMultiplicity > 1)
#' }
#' @author Michael Stravs, EAWAG <michael.stravs@@eawag.ch>
#' @export
filterPeaksMultiplicity <- function(w, recalcBest = TRUE)
{
	
	spectra <- lapply(w@spectra, function(sp)
			{
	      if(!sp@found)
	        return(sp)
				allPeaks <- lapply(sp@children, getData)
				childrenOK <- unlist(lapply(sp@children, function(c) c@ok))
				allPeaks <- allPeaks[childrenOK]
				allPeaks.df <- do.call(rbind, allPeaks)
				# sum the count per formula
				formulaMultiplicity <- table(allPeaks.df$formula[allPeaks.df$good %in% TRUE])
				allPeaks <- lapply(allPeaks, function(peaks)
						{
							peaks$formulaMultiplicity <- formulaMultiplicity[peaks$formula]
							#peaks$fM_factor <- as.factor(peaks$formulaMultiplicity)
							
							if(recalcBest == FALSE)
								return(peaks)
							# prioritize duplicate peaks
							# get unique peaks with their maximum-multiplicity formula attached
							peaksMult <- peaks[(peaks$good %in% TRUE),,drop=FALSE]
							if(nrow(peaksMult) == 0)
							  best_mult <- data.frame(mz = numeric(), bestMultiplicity = integer())
							else {
							  best_mult <- aggregate(as.integer(peaksMult$formulaMultiplicity), 
							                         by=list(peaksMult$mz), 
							                         function(x) max(c(0,x), na.rm=TRUE))
							  colnames(best_mult) <- c("mz", "bestMultiplicity")
							}
							peaks$bestMultiplicity <- NULL
							peaks <- merge(peaks, best_mult, by="mz", all.x = TRUE)
							
							#peaks <- peaks[peaks$formulaMultiplicity==peaks$bestMultiplicity,]
							
							# now we also have to recalculate dppmBest since the "old best" may have been
							# dropped.
							peaks$dppmBest <- NULL
							peaksBestMult <- peaks[peaks$formulaMultiplicity==peaks$bestMultiplicity & (peaks$good %in% TRUE),, drop=FALSE]
							if(nrow(peaksBestMult) == 0)
							  bestPpm <- data.frame(mz = numeric(), dppmBest = integer())
							else {
							  bestPpm <- aggregate(peaksBestMult$dppm, 
							                       list(peaksBestMult$mz),
							                       function(dppm) dppm[which.min(abs(dppm))])
							  colnames(bestPpm) <- c("mz", "dppmBest")  
							}
							peaks$dppmBest <- NULL
							peaks <- merge(peaks, bestPpm, all.x = TRUE)
							# return the multiplicity-culled (not yet filtered) peak table, which has
							# the new additional columns formulaMultiplicity and bestMultiplicity
							peaks
						})
				# if recalcBest is activated, we need to recalculate the multiplicity too, since peaks could have been dropped in
				# by the dppmBest recalculation
				# (note: this is to make the output consistent, so if you count the multiplicities of the peaks by hand you get
				# the same result as is written in the data. Otherwise a peak might show 2 multiplicity but in reality one of the
				# two was dropped. It is debatable if this is really necessary, since we might kick out one or two peaks too much
				# because of this step)
				# 
				if(recalcBest)
				{
					allPeaks.df <- do.call(rbind, allPeaks)
					# sum the count per formula
					#formulaMultiplicity <- table(allPeaks.df$formula[allPeaks.df$good %in% TRUE])
					formulaMultiplicity <- table(allPeaks.df$formula[(allPeaks.df$good %in% TRUE) & (allPeaks.df$dppm == allPeaks.df$dppmBest)])
					allPeaks <- lapply(allPeaks, function(peaks)
							{
								peaks$formulaMultiplicity <- formulaMultiplicity[peaks$formula]
								peaks
							})
					
				}
				# Write the data back to the spectra
				rewritePos <- which(childrenOK)
				for(i_ in seq_along(rewritePos))
				{
				  i <- rewritePos[[i_]]
					sp@children[[i]] <- addProperty(sp@children[[i]], "formulaMultiplicity", "integer", NA)
					sp@children[[i]] <- addProperty(sp@children[[i]], "bestMultiplicity", "integer", NA)
					sp@children[[i]] <- setData(sp@children[[i]], allPeaks[[i_]])
				}
				sp
			})
	w@spectra <- as(spectra, "SimpleList")
	
	w
}


#' filterMultiplicity
#' 
#' Multiplicity filtering: Removes peaks which occur only once in a n-spectra
#' set.
#' 
#' This function executes multiplicity filtering for a set of spectra using the
#' workhorse function \code{\link{filterPeaksMultiplicity}} (see details there)
#' and retrieves problematic filtered peaks (peaks which are of high intensity
#' but were discarded, because either no formula was assigned or it was not
#' present at least 2x), using the workhorse function
#' \code{\link{problematicPeaks}}. The results are returned in a format ready
#' for further processing with \code{\link{mbWorkflow}}.
#' 
#' @usage filterMultiplicity(w, archivename=NA, mode="pH", recalcBest = TRUE,
#' 		multiplicityFilter = getOption("RMassBank")$multiplicityFilter)
#' @param w Workspace containing the data to be processed (aggregate table and \code{RmbSpectraSet} objects)
#' @param archivename The archive name, used for generation of
#' archivename_Failpeaks.csv
#' @param mode Mode of ion analysis
#' @param recalcBest Boolean, whether to recalculate the formula multiplicity 
#' 		after the first multiplicity filtering step. Sometimes, setting this
#' 		to FALSE can be a solution if you have many compounds with e.g. fluorine
#' 		atoms, which often have multiple assigned formulas per peak and might occasionally
#' 		lose peaks because of that. 
#' @param multiplicityFilter Threshold for the multiplicity filter. If set to 1,
#' 		no filtering will apply (minimum 1 occurrence of peak). 2 equals minimum
#' 		2 occurrences etc. 
#' @return A list object with values: 
#' \item{peaksOK}{ Peaks with >1-fold formula multiplicity from the
#' 		"normal" peak analysis.  } 
#' \item{peaksReanOK}{ Peaks with >1-fold formula multiplicity from
#' 		peak reanalysis.  }
#' \item{peaksFiltered}{ All peaks with annotated formula multiplicity from
#' 		first analysis.  } 
#' \item{peaksFilteredReanalysis}{ All peaks with annotated
#' 		formula multiplicity from peak reanalysis.  } 
#' \item{peaksProblematic}{ Peaks with high intensity which do not match 
#' 		inclusion criteria -> possible false negatives. The list will be
#' 		exported into archivename_failpeaks.csv.
#' }
#' @author Michael Stravs
#' @seealso
#' \code{\link{filterPeaksMultiplicity}},\code{\link{problematicPeaks}}
#' @examples
#' \dontrun{
#'     refilteredRcSpecs <- filterMultiplicity(
#' 			w, "myarchive", "pH")
#' }
#' @export
filterMultiplicity <- function(w, archivename=NA, mode="pH", recalcBest = TRUE,
		multiplicityFilter = getOption("RMassBank")$multiplicityFilter)
{
  # Read multiplicity filter setting
  # For backwards compatibility: If the option is not set, define as 2
  # (which was the behaviour before we introduced the option)
  if(is.null(multiplicityFilter))
    multiplicityFilter <- 2
  
	# perform the counting:
    w <- filterPeaksMultiplicity(w, recalcBest)

	# now add the filter criterion to the spectra:
	spectra <- lapply(w@spectra, function(sp)
			{
				mz <- sp@mz
				if(length(sp@children) == 0)
					return(sp)
				children <- lapply(sp@children, function(ch)
						{
				      if(ch@ok == FALSE)
				        return(ch)
							# filterOK TRUE if multiplicity is sufficient
							ch <- addProperty(ch, "filterOK", "logical", NA)
							property(ch, "filterOK") <- (property(ch, "formulaMultiplicity") > (multiplicityFilter - 1)) &
							  (ch@good == TRUE) 
							# filterOK FALSE if the peak is a reanalyzed peak and is a M+H+ satellite
							peaks <- getData(ch)
							reanSat <- which(
									(peaks$formulaSource == "reanalysis") &
									(abs(peaks$mz - mz) < 1) &
									(peaks$filterOK)
									)
							peaks[reanSat, "filterOK"] <- rep(FALSE, length(reanSat))
							ch <- setData(ch, peaks)
							ch
						})
				sp@children <- as(children, "SimpleList")
				sp
			})
	w@spectra <- as(spectra, "SimpleList")
	return(w)
}

#' Return MS1 peaks to be used for recalibration
#' 
#' Returns the precursor peaks for all MS1 spectra in the \code{spec} dataset
#' with annotated formula to be used in recalibration.
#'  
#' For all spectra in \code{spec$specFound}, the precursor ion is extracted from
#' the MS1 precursor spectrum. All found ions are returned in a data frame with a
#' format matching \code{spec$peaksMatched} and therefore suitable for \code{rbind}ing
#' to the \code{spec$peaksMatched} table. However, only minimal information needed for
#' recalibration is returned. 
#' 
#' @usage  recalibrate.addMS1data(spec,mode="pH", recalibrateMS1Window = 
#' 				getOption("RMassBank")$recalibrateMS1Window)
#' @param spec A \code{msmsWorkspace} or \code{RmbSpectraSetList} containing spectra for which MS1 "peaks" should be "constructed". 
#' @param mode \code{"pH", "pNa", "pM", "pNH4",  "mH", "mM", "mFA"} for different ions 
#' 			([M+H]+, [M+Na]+, [M]+,  [M+NH4]+, [M-H]-, [M]-, [M+FA]-).
#' @param recalibrateMS1Window Window width to look for MS1 peaks to recalibrate (in ppm).
#' @return A dataframe with columns \code{mzFound, formula, mzCalc, dppm, dbe, int,
#' 		dppmBest, formulaCount, good, cpdID, scan, parentScan, dppmRc}. However,
#' 		columns \code{dbe, int, formulaCount, good, scan, parentScan} do not contain
#' 		real information and are provided only as fillers.
#' @examples \dontrun{
#' # More or less as used in recalibrateSpectra:
#' 		rcdata <- peaksMatched(w)
#' 		rcdata <- rcdata[rcdata$formulaCount == 1, ,drop=FALSE]
#' 		ms1data <- recalibrate.addMS1data(w, "pH", 15)
#' 		rcdata <- rbind(rcdata, ms1data)
#'  # ... continue constructing recalibration curve with rcdata
#' }
#' @author Michael Stravs, EAWAG <michael.stravs@@eawag.ch>
#' @export
recalibrate.addMS1data <- function(spec,mode="pH", recalibrateMS1Window = 
				getOption("RMassBank")$recalibrateMS1Window)
{
	## which_OK <- lapply(validPrecursors, function(pscan)
	##         {
	##             pplist <- as.data.frame(
	##                     mzR::peaks(msRaw, which(headerData$acquisitionNum == pscan)))
	##             colnames(pplist) <- c("mz","int")
	##             pplist <- subset(pplist, mz >= mzLimits$mzMin & mz <= mzLimits$mzMax)
	##             if(nrow(pplist) > 0)
	##                 return(TRUE)
	##             return(FALSE)
	##         })
	
	specFound <- selectSpectra(spec, "found", "object")
	
	ms1peaks <- lapply(specFound, function(cpd){
	    if(cpd@formula == "")
	      return(NULL)
			mzL <- findMz.formula(cpd@formula,mode,recalibrateMS1Window,0)
			mzCalc <- mzL$mzCenter
			ms1 <- mz(cpd@parent)
           
			mzFound <- ms1[which.min(abs(ms1 - mzL$mzCenter))]
			if(!length(mzFound)){
				return(c(
					mzFound = NA,
					mzCalc = mzCalc,
					dppm = NA
				))
			} else {
				dppmRc <- (mzFound/mzCalc - 1)*1e6
				return(c(
					mzFound = mzFound,
					mzCalc = mzCalc,
					dppm = dppmRc,
                    id=cpd@id
				))
			}
		})
	ms1peaks <- ms1peaks[which(!unlist(lapply(ms1peaks, is.null)))]
	ms1peaks <- as.data.frame(do.call(rbind, ms1peaks), stringsAsFactors=FALSE)
	# convert numbers to numeric
	tonum <- c("mzFound", "dppm", "mzCalc")
	ms1peaks[,tonum] <- as.numeric(unlist(ms1peaks[,tonum]))
	# throw out NA stuff
	ms1peaks <- ms1peaks[!is.na(ms1peaks$mzFound),]
	return(ms1peaks)
}


# Custom recalibration function: You can overwrite the recal function by
# making any function which takes rcdata$recalfield ~ rcdata$mzFound.
# The settings define which recal function is used
# getOption("RMassBank")$recalibrator = list(
#	MS1 = "recalibrate.loess",
#	MS2 = "recalibrate.loess")

#' Predefined recalibration functions.
#' 
#' Predefined fits to use for recalibration: Loess fit and GAM fit.
#' 
#' \code{recalibrate.loess()} provides a Loess fit (\code{recalibrate.loess}) 
#' to a given recalibration parameter.  
#' If MS and MS/MS data should be fit together, recalibrate.loess 
#' provides good default settings for Orbitrap instruments.
#' 
#' \code{recalibrate.identity()} returns a non-recalibration, i.e. a predictor
#' which predicts 0 for all input values. This can be used if the user wants to
#' skip recalibration in the RMassBank workflow.
#' 
#' #' \code{recalibrate.mean()} and \code{recalibrate.linear()} are simple recalibrations
#' which return a constant shift or a linear recalibration. They will be only useful
#' in particular cases.
#' 
#' \code{recalibrate()} itself is only a dummy function and does not do anything.
#' 
#' Alternatively other functions can be defined. Which functions are used for recalibration
#' is specified by the RMassBank options file. (Note: if \code{recalibrateMS1: common}, the
#' \code{recalibrator: MS1} value is irrelevant, since for a common curve generated with
#' the function specified in \code{recalibrator: MS2} will be used.)
#' 
#' @aliases recalibrate.loess recalibrate recalibrate.identity recalibrate.mean recalibrate.linear
#' @usage recalibrate.loess(rcdata)
#' 
#' 		recalibrate.identity(rcdata)
#' 
#' 		recalibrate.mean(rcdata)
#' 
#' 		recalibrate.linear(rcdata)
#' 
#' @param rcdata A data frame with at least the columns \code{recalfield} and
#' 			\code{mzFound}. \code{recalfield} will usually contain delta(ppm) or
#' 			delta(mz) values and is the target parameter for the recalibration.
#' @return Returns a model for recalibration to be used with \code{predict} and the like.
#' @examples \dontrun{
#' rcdata <- subset(spec$peaksMatched, formulaCount==1)
#' ms1data <- recalibrate.addMS1data(spec, mode, 15)
#' rcdata <- rbind(rcdata, ms1data)
#' rcdata$recalfield <- rcdata$dppm
#' rcCurve <- recalibrate.loess(rcdata)
#' # define a spectrum and recalibrate it
#' s <- matrix(c(100,150,200,88.8887,95.0005,222.2223), ncol=2)
#' colnames(s) <- c("mz", "int")
#' recalS <- recalibrateSingleSpec(s, rcCurve)
#' 
#' Alternative: define an custom recalibrator function with different parameters
#' recalibrate.MyOwnLoess <- function(rcdata)
#' {
#' 	return(loess(recalfield ~ mzFound, data=rcdata, family=c("symmetric"),
#' 					degree = 2, span=0.4))
#' }
#' # This can then be specified in the RMassBank settings file:
#' # recalibrateMS1: common
#' # recalibrator:
#' #    MS1: recalibrate.loess
#' #    MS2: recalibrate.MyOwnLoess")
#' # [...]
#' }
#' @author Michael Stravs, EAWAG <michael.stravs@@eawag.ch>
#' @export
recalibrate <- function()
{
	return(NA)
}

#' @export
recalibrate.loess <- function(rcdata)
{
  span <- 0.25
  # ex XCMS (permission by Steffen): heuristically decide on loess vs linear
  mingroups <- nrow(rcdata[!is.na(rcdata$mzFound),])
  if(mingroups < 4)
  {
    warning("recalibrate.loess: Not enough data points, omitting recalibration")
    return(recalibrate.identity(rcdata))
  } else if (mingroups*span < 4) {
    span <- 4/mingroups
    warning("recalibrate.loess: Span too small, resetting to ", round(span, 2))
  }
	return(loess(recalfield ~ mzFound, data=rcdata, family=c("symmetric"),
					degree = 1, span=0.25, surface="direct" ))
}

#' @export 
recalibrate.identity <- function(rcdata)
{
	return(lm(recalfield ~ 0, data=rcdata))
}

#' @export 
recalibrate.mean <- function(rcdata)
{
  return(lm(recalfield ~ 1, data=rcdata))
}

#' @export 
recalibrate.linear <- function(rcdata)
{
  return(lm(recalfield ~ mzFound, data=rcdata))
}

#' Standard progress bar hook.
#' 
#' This function provides a standard implementation for the progress bar in RMassBank.
#' 
#' RMassBank calls the progress bar function in the following three ways:
#' \code{pb <- progressBarHook(object=NULL, value=0, min=0, max=LEN)}
#' to create a new progress bar.
#' \code{pb <- progressBarHook(object=pb, value= VAL)}
#' to set the progress bar to a new value (between the set \code{min} and \code{max})
#' \code{progressBarHook(object=pb, close=TRUE)}
#' to close the progress bar. (The actual calls are performed with \code{do.call}, 
#' e.g. 
#' \code{progressbar <- "progressBarHook"
#' pb <- do.call(progressbar, list(object=pb, value= nProg))
#' }. See the source code for details.)
#' 
#' To substitute the standard progress bar for an alternative implementation (e.g. for
#' use in a GUI), the developer can write his own function which behaves in the same way
#' as \code{progressBarHook}, i.e. takes the same parameters and can be called in the 
#' same way. 
#'  
#' @param object An identifier representing an instance of a progress bar. 
#' @param value The new value to assign to the progress indicator
#' @param min The minimal value of the progress indicator
#' @param max The maximal value of the progress indicator
#' @param close If \code{TRUE}, the progress bar is closed.
#' @return Returns a progress bar instance identifier (i.e. an identifier
#' 		which can be used as \code{object} in subsequent calls.) 
#' 
#' @author Michele Stravs, Eawag <stravsmi@@eawag.ch>
#' @export
progressBarHook <- function(object = NULL, value = 0, min = 0, max = 100, close = FALSE)
{
	if(is.null(object))
	{
		object <- txtProgressBar(min, max, value, style=3, file=stderr())
	}
	if(close)
		close(object)
	else
	{
		setTxtProgressBar(object, value)
		return(object)
	}
}

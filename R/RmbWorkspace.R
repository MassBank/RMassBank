#' @import methods
NULL

setClassUnion("msmsWorkspaceOrNULL", "NULL")

#' Workspace for \code{msmsWorkflow} data
#' 
#' A workspace which stores input and output data for \code{\link{msmsWorkflow}}.
#' 
#' Slots:
#' 
#'  \describe{
#' 	\item{files}{The input file names}
#' 	\item{spectra}{The spectra per compound (\code{RmbSpectraSet}) extracted from the raw files}
#'  \item{aggregated}{A data.frame with an aggregated peak table from all \code{spectra}.
#' 		Further columns are added during processing.}
#'  \item{rc, rc.ms1}{The recalibration curves generated in workflow step 4.}
#'  \item{parent}{For the workflow steps after 4: the parent workspace containing the state (spectra, aggregate)
#' 		before recalibration, such that the workflow can be reprocessed from start.}
#' 	\item{archivename}{The base name of the files the archive is stored to during the workflow.}
#' \item{settings}{The RMassBank settings used during the workflow, if stored with the workspace.}#' 
#' }
#' 
#' Methods: \describe{
#' 	\item{show}{Shows a brief summary of the object and processing progress.}
#' 	}
#' 
## ' @slot files The input file names
## ' @slot spectra The spectra extracted from the raw files
## ' @slot analyzedSpecs The spectra with annotated peaks after workflow step 2.
## ' @slot aggregatedSpecs The \code{analyzedSpec} data regrouped and aggregated, 
## ' 		after workflow step 3.
## ' @slot rc,rc.ms1 The recalibration curves generated in workflow step 4.
## ' @slot recalibratedSpecs The spectra from \code{specs} recalibrated with the curves
## ' 		from \code{rc, rc,ms1}.
## ' @slot analyzedRcSpecs The recalibrated spectra with annotated peaks after 
## ' 		workflow step 5.
## ' @slot aggregatedRcSpecs The \code{analyzedRcSpec} data regrouped and aggregated, 
## ' 		after workflow step 6.
## ' @slot reanalyzedRcSpecs The regrouped and aggregated spectra, with added reanalyzed
## ' 		peaks (after step 7, see \code{\link{reanalyzeFailpeaks}}).
## ' @slot refilteredRcSpecs Final data to use for MassBank record creation after 
## ' 		multiplicity filtering (step 8).
##' 
## ' @method show,msmsWorkflow Shows a brief summary of the object. Currently only the included files.
#' 
#' @seealso \code{\link{msmsWorkflow}}
#' @author Michael Stravs, Eawag <michael.stravs@@eawag.ch>
#' @name msmsWorkspace-class
#' @docType class
#' @exportClass msmsWorkspace
#' @export
.msmsWorkspace <- setClass("msmsWorkspace",
		representation = representation(
				files = "character",
				spectra = "RmbSpectraSetList",
				aggregated = "data.frame",
				parent = "msmsWorkspaceOrNULL",
				rc = "ANY",
				rc.ms1 = "ANY",
				archivename = "character",
				settings = "list"
		),
		contains=c("Versioned"),
		prototype = prototype(
				# Version 2.0.2: the RmbSpectrum2 objects were raised to version 0.1.1.
        # Version 2.0.3: the RmbSpectraSet objects are now versioned and at 0.1.1.
        # Version 2.0.4: polarity into Spectrum1 and RmbSpectrum2 objects
				new("Versioned", versions=c(msmsWorkspace = "2.0.4")),
				parent = NULL
		)
)

setIs("msmsWorkspace", "msmsWorkspaceOrNULL")

#.msmsWorkspace <- setClass("msmsWorkspace",
#		representation = representation(
#				files = "character",
#				specs = "list",
#				analyzedSpecs = "list",
#				aggregatedSpecs = "list",
#				rc = "ANY",
#				rc.ms1 = "ANY",
#				recalibratedSpecs = "list",
#				analyzedRcSpecs = "list",
#				aggregatedRcSpecs = "list",
#				reanalyzedRcSpecs = "list",
#				refilteredRcSpecs = "list",
#				archivename = "character",
#				settings = "list"
#				),
#		contains=c("Versioned"),
#		prototype = prototype(
#				new("Versioned", versions=c(msmsWorkspace = "1.0.1"))
#				)
#		)
		


#' Workspace for \code{mbWorkflow} data
#' 
#' A workspace which stores input and output data for use with \code{mbWorkflow}.
#' 
#' Slots:
#'  \describe{
#' 	\item{spectra, aggregated}{The corresponding
#' 		 input data from \code{\link{msmsWorkspace-class}}}
#'  \item{additionalPeaks}{A list of additional peaks which can be loaded
#' 		using \code{\link{addPeaks}}.}
#'  \item{mbdata, mbdata_archive, mbdata_relisted}{Infolist data: Data for
#' 		annotation of MassBank records, which can be loaded using
#' 		\code{\link{loadInfolists}}.}
#'  \item{compiled, compiled_ok}{
#' 		Compiled tree-structured MassBank records. \code{compiled_ok} contains
#' 		only the compounds with at least one valid spectrum.}
#'  \item{mbfiles}{Compiled MassBank records in text representation.}
#'  \item{molfile}{MOL files with the compound structures.}
#'  \item{ok,problems}{Index lists for internal use which denote which compounds
#' 		have valid spectra.}
#' }
#' 
## ' @method show,mbWorkflow Shows a brief summary of the object. Currently only a stub. 
#' Methods:
#' \describe{
#' 		\item{show}{Shows a brief summary of the object. Currently only a stub.}
#' }
#' 
#' @seealso \code{\link{mbWorkflow}}
#' 
#' @docType class
#' @author Michael Stravs, Eawag <michael.stravs@@eawag.ch>
#' @name mbWorkspace-class
#' @exportClass mbWorkspace
#' 
setOldClass(c("tbl_df", "tbl", "data.frame"))

setClass("mbWorkspace",
		representation(
				# input data:
				spectra = "RmbSpectraSetList",
				aggregated = "data.frame",
				## aggregatedRcSpecs = "list",
				## refilteredRcSpecs = "list",
				additionalPeaks = "data.frame", # ex additional_peaks
				# infolists data:
				mbdata = "list",
				mbdata_archive = "data.frame",
				mbdata_relisted = "list",				
				# output data:
				compiled = "list",
				compiled_ok = "list",
				compiled_notOk = "list",
				mbfiles = "list",
				mbfiles_notOk = "list",
				molfile = "list",
				ok = "integer",
				problems = "integer"
				)
		)
		

#' Create new empty workspace or load saved data for \code{msmsWorkflow}
#' 
#' Creates an empty workspace or loads an existing workspace from disk.
#' 
#' \code{newMsmsWorkspace} creates a new empty workspace for use with \code{msmsWorkflow.}
#' 
#' \code{loadMsmsWorkspace} loads a workspace saved using \code{\link{archiveResults}}. 
#' Note that it also successfully loads data saved with the old RMassBank data format
#' into the new \code{msmsWorkspace} object.   
#' 
#' @aliases newMsmsWorkspace loadMsmsWorkspace
#' @param files If given, the files list to initialize the workspace with.
#' @return A new \code{msmsWorkspace} object
#' @seealso \code{\link{msmsWorkflow}}, \code{\link{msmsWorkspace-class}}
#' 
#' @author Michael Stravs, Eawag <michael.stravs@@eawag.ch>
#' @export
newMsmsWorkspace <- function(files = character(0))
{
	w <- new("msmsWorkspace")
	w@files <- files
	return(w)
}

#' @export
loadMsmsWorkspace <- function(fileName, loadSettings = FALSE)
{
	tempEnv <- new.env()
	load(fileName, envir=tempEnv)
	# Look if there is a msmsWorkspace in the file
	objs <- ls(tempEnv)
	isWs <- unlist(lapply(objs, function(obj) is(tempEnv[[obj]], "msmsWorkspace")))
	whichWs <- match(TRUE, isWs)
	# Found? Then just return it.
	if(!is.na(whichWs))
	{
		w <- tempEnv[[objs[[whichWs]]]]
		# If settings don't have to be loaded, erase them from the object
		if(loadSettings == FALSE)
			w@settings <- list()
	}
	# If there is no msmsWorkspace object in the workspace, this means that the workspace is version 1!
	else
	{
		w <- new("msmsWorkspace")
		dataset <- c("rc",
				"rc.ms1",
				"files",
				"specs",
				"analyzedSpecs",
				"aggregatedSpecs",
				"recalibratedSpecs",
				"analyzedRcSpecs",
				"aggregatedRcSpecs",
				"reanalyzedRcSpecs",
				"refilteredRcSpecs",
				"archivename")
		for(var in dataset)
		{
			if(exists(var, envir=tempEnv))
				slot(w, var, check=FALSE) <- tempEnv[[var]]
			classVersion(w)["msmsWorkspace"] <- "1.0.0"
		}
		# Check if settings exist...
		if((loadSettings == TRUE) && exists("RmbSettings", envir=tempEnv))
			w@settings <- tempEnv$RmbSettings
	}
	# process version updates
	updateClass <- FALSE
	if(!isVersioned(w)) updateClass <- TRUE
	else if(!all(isCurrent(w))) updateClass <- TRUE
	if(updateClass)
	{
		w <- updateObject(w)
	}
	
	# If loadSettings is set: load the settings into RMassBank
	if((loadSettings == TRUE) && (length(w@settings) > 0))
		loadRmbSettings(w@settings)
	else if (loadSettings == TRUE)
		warning("You specified to load RMassBank settings from the saved workspace, but no stored settings were found.")
	return(w)
}

#' Create new workspace for \code{mbWorkflow}
#' 
#' Creates a new workspace for use with \code{\link{mbWorkflow}}.
#' 
#' The workspace input data will be  loaded from the \code{\link{msmsWorkspace-class}}
#' object provided by the parameter \code{w}.
#' 
#' @param w The input \code{msmsWorkspace} to load input data from.
#' @return A new \code{mbWorkflow} object with the loaded input data.
#' 
#' @seealso \code{\link{mbWorkflow}}, \code{\link{msmsWorkspace-class}}
#' @author Michael Stravs, Eawag <michael.stravs@@eawag.ch>
#' @export
newMbWorkspace <- function(w)
{
  #w <- fillback(w)
	mb <- new("mbWorkspace",
			spectra = w@spectra#,
	#		aggregated = w@aggregated
			)
	return(mb)
}


#' @name show,msmsWorkspace-method
#' @aliases show,msmsWorkspace-method
#' @param object The \code{msmsWorkspace} to display.
#' @docType methods
#' @rdname msmsWorkspace-class
#' @export
setMethod("show", "msmsWorkspace",
		function(object) {
			cat("Object of class \"",class(object),"\"\n",sep="")
			cat(" with files: \n")
			sapply(basename(object@files), function(x) cat("  -", x, "\n"))
                        
						progress <- findProgress(object)
						
						if(4 %in% progress)
							o123 <- object@parent
						else
							o123 <- object
						
						## msmsWorkflow: Step 1. Get peaks ?
						if(1 %in% progress){
							numspecs <- 0
							dummy <- sapply(o123@spectra, function(x) cat(" -", x@id, "\t foundOK:", x@found, "\n"))
							cat("Peaks found:\n")
							ids <- vector()
							dummy1 <- sapply(o123@spectra, function(x) {
								cat(" -", x@id, "\t peaks:",
								sapply(x@children, function(s) s@peaksCount), "\n")
								ids <<- c(ids, x@id)
								numspecs <<- numspecs + 1
								return(sapply(x@children, function(s) s@peaksCount))
							})
						}
									
									
                        ## msmsWorkflow: Step 2. First analysis pre recalibration
                        if(2 %in% progress){
							cat("Peaks annotated after Step 2:\n")
							dummy2 <- sapply(o123@spectra, function(x) {
								cat(" -", x@id, "\t peaks:",
								sapply(x@children, function(s) length(which(s@good))), "\n")
								return(sapply(x@children, function(s) length(which(s@good))))
							})
							
							whichok <- which(sapply(dummy1,length) != 0)
							anyok <- whichok[1]
							
							nrowd1 <- nrow(dummy1)
							if(is.null(nrowd1)) 
							nrowd1 <- 1
							
							
							dummy2 <- matrix(unlist(dummy2), ncol = nrowd1, byrow = TRUE)
							dummy1 <- matrix(unlist(dummy1), ncol = nrowd1, byrow = TRUE)
							
							PeakMat <- dummy1-dummy2
							
							tempids <- ids[whichok]
							cat("Peaks without annotation:\n")
							sapply(split(PeakMat, rep(1:nrow(PeakMat), each = ncol(PeakMat))), function(x){
								cat(" -", tempids[1], "\t number of peaks filtered:", x, "\n")
								tempids <<- tempids[-1]
							})
						}
                        ## msmsWorkflow: Step 3. Aggregate all spectra
						if(3 %in% progress){
							cat("Peaks aggregated after Step 3:\n")
							cat("Matched Peaks:\n")
							pm <- peaksMatched(o123)
							dummy <- sapply(unique(pm[,"cpdID"]), function(x) cat(" -", x, "\t peaks:",
																		  sapply(x, function(y){ 
																			  compoundIndex <- which(pm[,"cpdID"] == y)
																			  peaksMatched <- pm[compoundIndex,]
																			  uscans <- unique(peaksMatched[,"scan"])
																			  return(sapply(uscans, function(z){
																				uscantemp <- which(peaksMatched[,"scan"] == z)
																				return(length(unique(peaksMatched[uscantemp,"mzFound"])))
																			  }))
																		  }), "\n"))
							cat("Unmatched Peaks:\n")
							pu <- peaksUnmatched(o123)
							dummy <- sapply(unique(pu[,"cpdID"]), function(x) cat(" -", x, "\t peaks:",
																		  sapply(x, function(y){ 
																			  compoundIndex <- which(pu[,"cpdID"] == y)
																			  peaksUnmatched <- pu[compoundIndex,]
																			  uscans <- unique(peaksUnmatched[,"scan"])
																			  return(sapply(uscans, function(z){
																				uscantemp <- which(peaksUnmatched[,"scan"] == z)
																				return(length(unique(peaksUnmatched[uscantemp,"mzFound"])))
																			  }))
																		  }), "\n"))
						}
						
                        ## msmsWorkflow: Step 4. Recalibrate m/z values in raw spectra
						if(4 %in% progress){
							cat("Peaks successfully recalibrated:\n")
						}
						
                        ## msmsWorkflow: Step 5. Reanalyze recalibrated spectra
						if(5 %in% progress){
							cat("Peaks found after Step 5:\n")

							dummy5 <- sapply(object@spectra, function(x) {
										cat(" -", x@id, "\t peaks:",
												sapply(x@children, function(s) length(which(s@good))), "\n")
										return(sapply(x@children, function(s) length(which(s@good))))
									})
							
							dummy5 <- matrix(unlist(dummy5), ncol = nrowd1, byrow = TRUE)
							
							PeakMat <- dummy1-dummy5
							
							tempids <- ids[whichok]
							cat("Peaks without annotation in reanalyzed recalibrated peaks:\n")
							sapply(split(PeakMat, rep(1:nrow(PeakMat), each = ncol(PeakMat))), function(x){
								cat(" -", tempids[1], "\t number of peaks filtered:", x, "\n")
								tempids <<- tempids[-1]
							})
						}
						
                        ## msmsWorkflow: Step 6. Aggregate recalibrated results
						if(6 %in% progress){
							cat("Peaks found after Step 6:\n")

							cat("Matched Peaks:\n")
							pm <- peaksMatched(object)
							dummy <- sapply(unique(pm[,"cpdID"]), function(x) cat(" -", x, "\t peaks:",
												sapply(x, function(y){ 
															compoundIndex <- which(pm[,"cpdID"] == y)
															peaksMatched <- pm[compoundIndex,]
															uscans <- unique(peaksMatched[,"scan"])
															return(sapply(uscans, function(z){
																				uscantemp <- which(peaksMatched[,"scan"] == z)
																				return(length(unique(peaksMatched[uscantemp,"mzFound"])))
																			}))
														}), "\n"))
							cat("Unmatched Peaks:\n")
							pu <- peaksUnmatched(object)
							dummy <- sapply(unique(pu[,"cpdID"]), function(x) cat(" -", x, "\t peaks:",
												sapply(x, function(y){ 
															compoundIndex <- which(pu[,"cpdID"] == y)
															peaksUnmatched <- pu[compoundIndex,]
															uscans <- unique(peaksUnmatched[,"scan"])
															return(sapply(uscans, function(z){
																				uscantemp <- which(peaksUnmatched[,"scan"] == z)
																				return(length(unique(peaksUnmatched[uscantemp,"mzFound"])))
																			}))
														}), "\n"))
						}
							
                        ## msmsWorkflow: Step 7. Reanalyze fail peaks for N2 + O
						if(7 %in% progress){
							cat("Peaks added in Step 7:\n")
							
							pmr <- object@aggregated[object@aggregated$matchedReanalysis & !is.na(object@aggregated$matchedReanalysis),,drop=FALSE]
							
							dummy <- sapply(unique(pmr[,"cpdID"]), function(x) cat(" -", x, "\t peaks:",
																		  sapply(x, function(y){ 
																			  compoundIndex <- which(pmr[,"cpdID"] == y)
																			  peaksMatchedReanalysis <- pmr[compoundIndex,]
																			  uscans <- unique(peaksMatchedReanalysis[,"scan"])
																			  return(sapply(uscans, function(z){
																				uscantemp <- which(peaksMatchedReanalysis[,"scan"] == z)
																				return(length(unique(peaksMatchedReanalysis[uscantemp,"mzFound"])))
																			  }))
																		  }), "\n"))
						}
						
                        ## msmsWorkflow: Step 8. Peak multiplicity filtering 
                        if(8 %in% progress){
							cat("After Step 8: multiplicity filtering:\n")
							show(table(object@aggregated[object@aggregated$filterOK & is.na(object@aggregated$matchedReanalysis), "cpdID"]))
							show(table(object@aggregated[object@aggregated$filterOK & !is.na(object@aggregated$matchedReanalysis), "cpdID"]))
						}
                      })


#'
#' @name show,mbWorkspace-method 
#' @aliases show,mbWorkspace-method
#' @rdname mbWorkspace-class
#' @param object The \code{mbWorkspace} to display.
#' @docType methods
#' @export
setMethod("show", "mbWorkspace",
		function(object) {
			cat("Object of class \"",class(object),"\"\n",sep="")
                        str(object, max.level=2)
		}) 
		
#' Plots mbWorkspaces
#' 
#' Plots the peaks of one or two \code{mbWorkspace} to compare them.
#' 
#' This functions plots one or two \code{mbWorkspace}s in case the use has used different methods to acquire
#' similar spectra. \code{w1} must always be supplied, while \code{w2} is optional. The wokspaces need to be fully processed
#' for this function to work.
#'
#' @param w1 The \code{mbWorkspace} to be plotted
#' @param w2 Another optional \code{mbWorkspace} be plotted as a reference.
#' @return A logical indicating whether the information was plotted or not
#' @author Erik Mueller
#' @examples
#' 
#' #
#' \dontrun{plotMbWorkspaces(w1,w2)}
#' 
#' @export
plotMbWorkspaces <- function(w1, w2=NULL){
	if(class(w1) != "mbWorkspace"){
		stop("The first supplied argument must be an mbWorkspace")
	}
	if(!is.null(w2)){
		if(class(w2) != "mbWorkspace"){
			stop("If there is a second argument supplied it must be of class mbWorkspace")
		}
		pl2 <- lapply(w2@compiled_ok,function(x){
			lapply(x,function(y) y[['PK$PEAK']][,c("m/z","rel.int.")])
		})
		plot_title <- lapply(w2@compiled_ok,function(x){
			lapply(x,function(y) y[['ACCESSION']])
		})
	}
	
	
	
	pl1 <- lapply(w1@compiled_ok,function(x){
		lapply(x,function(y) y[['PK$PEAK']][,c("m/z","rel.int.")])
	})
	
	plot_title <- lapply(w1@compiled_ok,function(x){
		lapply(x,function(y) y[['ACCESSION']])
	})
	
	
	
	maxpeaks <- c(-1000,1000)
	
	if(length(pl1) != 0){
		for(i in 1:length(pl1)){
			currentCompound <- pl1[[i]]
			
			for(j in 1:length(currentCompound)){
				currentSpectrum <- currentCompound[[j]]
				plot(currentSpectrum[,"m/z"],currentSpectrum[,"rel.int."], type="h",col="green",lwd=3,ylim=maxpeaks,xlab = "mz", ylab="rel.int", main=plot_title[[i]][[j]])
				if(!is.null(w2)){
				  graphics::points(pl2[[i]][[j]][,"m/z"],-pl2[[i]][[j]][,"rel.int."], type="h",col="red",lwd=3)
				}
			}
		}
	}
	
	return(TRUE)
}
		

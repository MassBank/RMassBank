

eicWorkflow <- function(
    w, mode="pH", steps=c(1:8),
    settings = getOption("RMassBank"),
    progressbar = "progressBarHook"
    ) {

  # Step 1: acquire EIC for parent from files
  if(1 %in% steps)
  {
    
    rmb_log_info("eicWorkflow: Step 1. Acquire parent EIC from files")
    spectra_count <- length(w@spectra)
    pb <- do.call(progressbar, 
                  list(object=NULL, value=0, min=0, max=spectra_count))
    
    w@spectra <- w@spectra |> as.list() |>
      purrr::map2(w@files, 
                  extractCpdParentEic, 
                  mode=mode, 
                  ppm=settings$findMsMsRawSettings$ppmFine,
                  progressbar = list(hook=progressbar, pb=pb)
          ) |>
      as("SimpleList")
    do.call(progressbar, list(object=pb, close=TRUE))
  }
  w
  
}


extractCpdParentEic <- function(cpd, file, mode, ppm, progressbar = NULL) {
  
  d <- mzR::openMSfile(file)
  h <- mzR::header(d)
  h <- h[h$polarity == RMassBank:::getAdductPolarity(mode),]
  mz <- RMassBank::findMz(cpd@id, mode = mode, ppm = ppm)
  eic <- RMassBank::findEIC(d, mz, headerCache = h)
  attr(cpd, "eic") <- eic
  if(!is.null(progressbar)) {
    value <- do.call(progressbar$hook, list(object=progressbar$pb, value=NULL))
    do.call(progressbar$hook, list(object=progressbar$pb, value=value+1))
  }
    
  return(cpd)
}

extractCpdFragmentEics <- function(cpd) {
  
}

extractSpectrumFragmentEics <- function(cpd, sp) {
  
}
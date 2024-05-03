

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
  
  # Step 2: acquire EIC for fragments from files
  if(2 %in% steps)
  {
    
    rmb_log_info("eicWorkflow: Step 2. Acquire fragment EIC from files")
    spectra_count <- length(w@spectra)
    pb <- do.call(progressbar, 
                  list(object=NULL, value=0, min=0, max=spectra_count))
    
    w@spectra <- w@spectra |> as.list() |>
      purrr::map2(w@files, 
                  extractCpdFragmentEics, 
                  mode=mode, 
                  ppm=settings$findMsMsRawSettings$ppmFine,
                  precursor_dmz = settings$findMsMsRawSettings$mzCoarse,
                  rt_window = settings$rtMargin,
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


.fixHeaderData  <- function(headerData) {
  
  headerData$precursorScanNum <- NA
  headerData[which(headerData$msLevel == 1),"precursorScanNum"] <-
    headerData[which(headerData$msLevel == 1),"acquisitionNum"]
  headerData[,"precursorScanNum"] <- RMassBank:::.locf(headerData[,"precursorScanNum"])
  # Clear the actual MS1 precursor scan number again
  headerData[which(headerData$msLevel == 1),"precursorScanNum"] <- 0
  # Remove precursors which are still NA in precursor scan num.
  # This removes a bug when filling precursor if the first scan(s) are MS2 before a
  # MS1 scan appears. The resulting NA values in precursorScanNum are problematic downstream.
  headerData <- headerData[!is.na(headerData$precursorScanNum),]
  
  return(headerData)
}

extractCpdFragmentEics <- function(cpd, file, mode, 
                                   ppm, precursor_dmz, rt_window,
                                   progressbar = NULL) {
  
  attr(cpd, "file") <- file

  d <- mzR::openMSfile(file)
  h <- mzR::header(d)
  p <- RMassBank::makePeaksCache(d, h)
  
  # 1.6.22: For Exploris data, we need to fill precursor here.
  h <- .fixHeaderData(h)
  
  if(length(cpd@children) == 0)
    return(cpd)
  
  if(!cpd@found)
    return(cpd)
  
  cpd@children <- cpd@children |>
    as.list() |>
    purrr::map(extractSpectrumFragmentEics,
        header = h, peaksCache = p,
        ppm=ppm, precursor_dmz=precursor_dmz, rt_window=rt_window * 60) |>
    as("SimpleList")
  
  
  
  
  if(!is.null(progressbar)) {
    value <- do.call(progressbar$hook, list(object=progressbar$pb, value=NULL))
    do.call(progressbar$hook, list(object=progressbar$pb, value=value+1))
  }
  
  return(cpd)
}

extractSpectrumFragmentEics <- function(
    sp, header, peaksCache,
    ppm, precursor_dmz, rt_window) {
  hSub <- header |> dplyr::filter(
    abs(precursorMZ - sp@precursorMz) < precursor_dmz,
    msLevel == sp@msLevel,
    polarity == sp@polarity,
    abs(retentionTime - sp@rt) < rt_window,
    collisionEnergy == sp@collisionEnergy
  )
  if(nrow(hSub) == 0)
  {
    attr(sp, "eics") <- list()
    return(sp)
  }
  hSub$msLevel <- 1
  pSub <- peaksCache[hSub$seqNum]
  
  mzs <- sp@mz
  if(!is.null(RMassBank::property(sp, "mzRaw")))
    mzs <- RMassBank::property(sp, "mzRaw")
  
  eics <- lapply(mzs, function(mz) {
    eic <- RMassBank::findEIC(d, mz, RMassBank::ppm(mz, ppm, p = TRUE), headerCache = hSub, peaksCache = pSub)
    eic$precursorScan <- hSub$precursorScanNum
    eic$mz <- mz
    eic
  })
  
  attr(sp, "eics") <- eics
  return(sp)
}
# Peak EIC correlation while filtering out zero rows
.eicScoreCor <- function(prec, eic) {
  eicSum <- colSums(abs(eic))
  .eicScore <- cor(prec, eic[,eicSum > 0,drop=FALSE])
  eicScore <- rep(NA_real_, ncol(eic))
  eicScore[eicSum > 0] <- .eicScore
  eicScore
} 

# Peak EIC dot product score
.eicScoreDot <- function(prec, eic) {
  normX <- sqrt(sum(prec^2))
  normY <- sqrt(colSums(eic^2))
  dot <- t(prec) %*% eic
  score <- dot / (normX * normY)
  score[normY == 0] <- NA
  score
}


f1 <- function(sens, spec)
  2 * sens * spec / (sens + spec)
fBeta <- function(beta) function(sens, spec)
  (1 + beta^2) * (sens * spec) / ((spec * beta^2) + sens)


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
                  extractParentEic, 
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
    
    rmb_log_info("eicWorkflow: Step 3. Acquire fragment EIC from files")
    spectra_count <- length(w@spectra)
    pb <- do.call(progressbar, 
                  list(object=NULL, value=0, min=0, max=spectra_count))
    
    w@spectra <- w@spectra |> as.list() |>
      purrr::map2(w@files, 
                  extractFragmentEics, 
                  mode=mode, 
                  ppm=settings$findMsMsRawSettings$ppmFine,
                  precursor_dmz = settings$findMsMsRawSettings$mzCoarse,
                  rt_window = settings$rtMargin,
                  progressbar = list(hook=progressbar, pb=pb)
      ) |>
      as("SimpleList")
    do.call(progressbar, list(object=pb, close=TRUE))
  }
  
  # Step 3: calculate dot score and correlations
  if(3 %in% steps)
  {
    
    rmb_log_info("eicWorkflow: Step 2. Calculate fragment EIC correlations")
    spectra_count <- length(w@spectra)
    pb <- do.call(progressbar, 
                  list(object=NULL, value=0, min=0, max=spectra_count))
    
    w@spectra <- w@spectra |> as.list() |>
      purrr::map(correlateEics, 
                  progressbar = list(hook=progressbar, pb=pb)
      ) |>
      as("SimpleList")
    do.call(progressbar, list(object=pb, close=TRUE))
  }
  
  
  w
  
}


extractParentEic <- function(cpd, file, mode, ppm, progressbar = NULL) {
  
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




extractFragmentEics <- function (x, ...) {
  UseMethod("extractFragmentEics", x)
}


extractFragmentEics.RmbSpectraSet <- function(cpd, file, mode, 
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
    purrr::map(extractFragmentEics,
        header = h, peaksCache = p,
        ppm=ppm, precursor_dmz=precursor_dmz, rt_window=rt_window * 60) |>
    as("SimpleList")
  
  
  
  
  if(!is.null(progressbar)) {
    value <- do.call(progressbar$hook, list(object=progressbar$pb, value=NULL))
    do.call(progressbar$hook, list(object=progressbar$pb, value=value+1))
  }
  
  return(cpd)
}

extractFragmentEics.RmbSpectrum2 <- function(
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


correlateEics <- function (x, ...) {
  UseMethod("correlateEics", x)
}

correlateEics.RmbSpectraSet <- function(cpd,
                                     progressbar = NULL) {
  # correlate EICs with precursor EIC
  # note: REMOVE the scan itself and the precursor itself, to zero out the one-point correlations!

  if(length(cpd@children) == 0)
    return(cpd)
  
  if(!cpd@found)
    return(cpd)
    
  cpd@children <- cpd@children |>
    as.list() |>
    purrr::map(correlateEics, cpd = cpd) |>
    as("SimpleList")
  
  if(!is.null(progressbar)) {
    value <- do.call(progressbar$hook, list(object=progressbar$pb, value=NULL))
    do.call(progressbar$hook, list(object=progressbar$pb, value=value+1))
  }
  
      
  cpd
}

correlateEics.RmbSpectrum2 <- function(sp, cpd) {

    eic <- attr(sp, "eic") |>
      dplyr::bind_rows(.id = "mzIndex")
    if(nrow(eic) == 0)
    {
      sp <- RMassBank::addProperty(sp, "eicScoreCor", type = "numeric", NA)
      sp <- RMassBank::addProperty(sp, "eicScoreDot", type = "numeric", NA)
      return(sp)
    }
    
    eic <- eic |>
      dplyr::mutate(mzIndex = as.numeric(as.character(mzIndex))) |>
      dplyr::select(-mz) |>
      tidyr::pivot_wider(names_from = c("mzIndex"), values_from = "intensity", names_prefix = "mz_")
      
    eicPrecursor <- attr(cpd, "eic")
    eic <- eic |> dplyr::left_join(
      eicPrecursor |> dplyr::rename(precursorScan = scan),
      by="precursorScan"
      )
      
      
    # filter out "this scan"
    eic <- eic |> dplyr::filter(precursorScan != sp@precScanNum)
    
    eicPrecursor <- eic[,"intensity"] |> as.matrix()
    eic <- eic |> dplyr::select(starts_with("mz_")) |> as.matrix()
    
    eicScoreCor <- .eicScoreCor(eicPrecursor, eic)
    eicScoreDot <- .eicScoreDot(eicPrecursor, eic)
      
    sp <- RMassBank::addProperty(sp, "eicScoreCor", type = "numeric", NA)
    property(sp, "eicScoreCor") <- as.vector(eicScoreCor)
    sp <- RMassBank::addProperty(sp, "eicScoreDot", type = "numeric", NA)
    property(sp, "eicScoreDot") <- as.vector(eicScoreDot)

    sp
}
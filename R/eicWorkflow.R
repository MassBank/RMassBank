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
    
    rmb_log_info("eicWorkflow: Step 2. Acquire fragment EIC from files")
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
    
    rmb_log_info("eicWorkflow: Step 3. Calculate fragment EIC correlations")
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
  
  # Step 3: calculate dot score and correlations
  if(4 %in% steps)
  {
    rmb_log_info("eicWorkflow: Step 4. Calculate correlation thresholds")
    w <- w |>  computeEicThreshold()
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


.getSetpoint <- function(ag, scoreCol, beta) {
  
  fsc <- fBeta(beta)
  setpoint <- roc(ag$good, ag[[scoreCol]], levels=c("FALSE", "TRUE"), direction="<")
  plot(setpoint)
  fscore <- fsc(setpoint$sensitivities, setpoint$specificities)
  plot(fscore)
  maxFscore <- which.max(fscore)
  return(setpoint$thresholds[maxFscore])
}

computeEicThreshold <- function(w, beta = 1.5) {
  
  # aggregate and keep the best result for every peak
  ag <- RMassBank::aggregateSpectra(w)
  ag <- ag |> 
    dplyr::group_by(cpdID, scan, mzFound) |> 
    dplyr::arrange(!good, abs(dppm)) |> 
    dplyr::slice(1) 
  
  # for every formula, keep the best correlation of all children;
  # when no formula present, keep the peak alone
  ag$formula_ <- factor(ag$formula) |> as.numeric()
  ag$formula_[is.na(ag$formula_)] <- -seq_along(ag$formula_[is.na(ag$formula_)])
  
  ag$eicScoreCor0 <- ag$eicScoreCor |> tidyr::replace_na(0)
  ag$eicScoreDot0 <- ag$eicScoreDot |> tidyr::replace_na(0)
  
  score_cols <- c("eicScoreCor", "eicScoreDot", "eicScoreCor0", "eicScoreDot0")
  setpoints <- purrr::map(score_cols, \(col) .getSetpoint(ag, col, beta))
  names(setpoints) <- score_cols
  
  attr(w, "eicScoreFilter") <- setpoints
  w
}

autoReview <- function (x, ...) {
  UseMethod("autoReview", x)
}

autoReview.RmbSpectraSet <- function(cpd, ...) {
  cpd@children <- cpd@children |>
    as.list() |>
    purrr::map(autoReview, ..., cpd=cpd) |>
    as("SimpleList")
  cpd
}

autoReview.RmbSpectrum2 <- function(sp, cpd, specOkLimit) {
  
  d <- getData(sp)
  mx <- max(d$intensity, 0, na.rm = TRUE)
  sp@ok <- sp@ok & (mx > specOkLimit)
  sp
}


autoReview.msmsWorkspace <- function(w, settings = getOption("RMassBank")) {
  
  specOkLimit <- as.numeric(settings$filterSettings$specOkLimit)
                            
  not_empty <- w@spectra |> 
    as.list() |> 
    purrr::map_lgl( ~ length(.x@children) > 0)
  
  glue::glue("{w@files[!not_empty]}: removing empty compound") |>
    purrr::walk(rmb_log_info)
  
  w@spectra <- w@spectra[not_empty]
  w@files <- w@files[not_empty]
  
  # Apply the specified specOK cutoff to reduce excessive work for nothing
  w@spectra <- w@spectra |>
    as.list() |>
    purrr::map(autoReview, specOkLimit = specOkLimit) |>
    as("SimpleList")
  
  # Check to exclude 0-spectra-found cpds
  w@spectra <- w@spectra |>
    as.list() |>
    purrr::map(function(cpd) {
    n_ok <- cpd@children |>
      as.list() |>
      purrr::map_lgl(~ isTRUE(.x@ok)) |>
      sum()
    cpd@found <- n_ok > 0
    cpd
    }) |>
    as("SimpleList")
  
  not_problematic <- w@spectra |> as.list() |>
    purrr::map_lgl(function(cpd) isTRUE(cpd@found))
  glue::glue("{w@files[!not_problematic]}: removing problematic compound") |>
    purrr::walk(rmb_log_info)
  w@spectra <- w@spectra[not_problematic]
  w@files <- w@files[not_problematic]
  
  w
}


loadReview <- function(w, 
                       from = NULL,
                       cutoffs = NULL,
                       cpds = NULL,
                       spectra = NULL,
                       settings = getOption("RMassBank")
) {
  
  if(!is.null(settings$eicCorrMetric))
    metric <- settings$eicCorrMetric
  else {
    rmb_log_warn("Correlation metric not specified in settings, using eicScoreCor")
    metric <- "eicScoreCor"
  }
    
  
  if(is.null(from) & any(is.null(cutoff, cpds, spectra)))
    stop("Either 'from' or all three source files 'cutoffs, cpds, spectra' need to be specified")
  if(!is.null(from) & any(!is.null(cutoff, cpds, spectra)))
    stop("Only either 'from' or all three source files 'cutoffs, cpds, spectra' may be specified")
  if(!is.null(from)) {
    cutoff <- glue::glue("{from}_score_cutoff.csv")
    cpds <- glue::glue("{from}_cpds_ok.csv")
    spectra <- glue::glue("{from}_spec_ok.csv")
  }
  if(!is.numeric(cutoff))
    cutoff <- read_file(cutoff) %>% as.numeric()
  if(!is.data.frame(cpds))
    cpds <- readr::read_csv(cpds)
  if(!is.data.frame(spectra))
    spectra <- readr::read_csv(spectra)
  specOk <- spectra
  cpdOk <- cpds
  
  # Inject review data into compounds,
  # then select and deselect spectra based on the choices from the review files
  for(i in seq_along(w@spectra)) {
    attr(w@spectra[[i]], "specOK") <- specOk %>% dplyr::filter(cpd == i)
    attr(w@spectra[[i]], "threshold") <- cpdOk$threshold[[i]]
  }
  w@spectra <- w@spectra[cpdOk$ok]
  w@files <- w@files[cpdOk$ok]
  
  # Apply specOk after applying cpdOk
  w@spectra <- w@spectra |>
    as.list() |>
    purrr::map(function(cpd) {
      specOK <- attr(cpd, "specOK")
      assert_that(nrow(specOK) == length(cpd@children))
      for(i in seq_along(cpd@children)) {
        cpd@children[[i]]@ok <- specOK$ok[[i]]
        attr(cpd@children[[i]], "threshold") <- specOK$threshold[[i]]
      }
      cpd
    }) |>
    as("SimpleList")
  
  # Warn if for whatever reason the cutoff couldn't be read
  if(!is.na(cutoff)) {
    eicScoreLimit <- cutoff
  } else {
    rmb_log_warn("score cutoff could not be read from review data, using calculated default")
    eicScoreLimit <- attr(w, "eicScoreFilter")[[metric]]
  }
  attr(w, "eicScoreLimit") <- eicScoreLimit
  w
}



eicCorrFilter <- function (x, ...) {
  UseMethod("eicCorrFilter", x)
}

eicCorrFilter.msmsWorkspace <- function(w, settings = getOption("RMassBank")) {
  
  if(!is.null(settings$eicCorrMetric))
    metric <- settings$eicCorrMetric
  else {
    rmb_log_warn("Correlation metric not specified in settings, using eicScoreCor")
    metric <- "eicScoreCor"
  }
  
  if(is.null(attr(w, "eicScoreLimit"))) {
    # this is info not warn, because reviewing by hand is not mandatory part
    # of the process; reviewing and EIC filtering are distinct processes
    eicScoreLimit <- attr(w, "eicScoreFilter")[[metric]]
    rmb_log_info(glue::glue(
      "No review was applied, using calculated correlation cutoff {eicScoreLimit}"
    ))
  } else {
    eicScoreLimit <- attr(w, "eicScoreLimit")
    rmb_log_info(glue::glue(
      "No review was applied, using calculated correlation cutoff {eicScoreLimit}"
    ))
  }
  
  w@spectra <- w@spectra |>
    as.list() |>
    purrr::map(eicCorrFilter, metric=metric, eicScoreLimit=eicScoreLimit) |>
    as("SimpleList")
  w
}

eicCorrFilter.RmbSpectraSet <- function(cpd, ...) {
  
  if(!cpd@found)
    return(cpd)
  cpd@children <- cpd@children |>
    as.list() |>
    purrr::map(eicCorrFilter, cpd=cpd, ...) |>
    as("SimpleList")
  cpd
}

eicCorrFilter.RmbSpectrum2 <- function(sp, cpd, metric, eicScoreLimit) {
  
  if(!sp@ok)
    return(sp)
  # find the valid threshold for this spectrum:
  # can be the global, compound-specific or spectrum-specific threshold
  # Note: the attr() values are 
  # * NULL if no review was applied
  # * NA if the value for this compound/spectrum was not overridden (the most common case)
  # * not NA if a specific value per compound / spectrum was set.
  # eicScoreLimit is always set and comes either from the automatic estimation
  # or from the review data.
  thresholds <- c(
    eicScoreLimit,
    attr(cpd, "threshold"),
    attr(sp, "threshold")
  )
  threshold <- tail(thresholds[!is.na(thresholds)], 1)
  
  d <- getData(sp)
  d <- d |>
    dplyr::group_by(mz) |>
    dplyr::arrange(!good, desc(formulaMultiplicity), abs(dppm)) |>
    dplyr::slice(1)
  d <- d |>
    dplyr::mutate(
      good_ = good,
      good = (.data[[metric]] > threshold) & !satellite & !low
    ) |> dplyr::filter(good)
  sp <- setData(sp, as.data.frame(d))
  sp
}



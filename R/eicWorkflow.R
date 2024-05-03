

eicWorkflow <- function(w, mode="pH", steps=c(1:8),
          settings = getOption("RMassBank")) {

  # Step 1: acquire EIC for parent from files
  if(1 %in% steps)
  {
    rmb_log_info("eicWorkflow: Step 1. Acquire parent EIC from files")
    w@spectra <- w@spectra |> as.list() |>
      purrr::map2(w@files, 
                  extractCpdParentEic, 
                  mode=mode, 
                  ppm=settings$filterSettings$ppmFine
          ) |>
      as("SimpleList")
  }
  w
  
}


extractCpdParentEic <- function(cpd, file, mode, ppm) {
  
  d <- mzR::openMSfile(file)
  h <- header(d)
  h <- h[h$polarity == RMassBank:::getAdductPolarity(mode),]
  mz <- findMz(cpd@id, mode = mode, ppm = ppm)
  eic <- findEIC(d, mz, headerCache = h)
  attr(cpd, "eic") <- eic
  return(cpd)
}

extractCpdFragmentEics <- function(cpd) {
  
}

extractSpectrumFragmentEics <- function(cpd, sp) {
  
}
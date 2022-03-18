library(S4Vectors)
library(testthat)
library(dplyr)

filterSettings <- list(ppmHighMass = 10, ppmLowMass = 15, massRangeDivision = 120,      
                       ppmFine = 5, prelimCut = 10000, prelimCutRatio = 0, fineCut = 0,      
                       fineCutRatio = 0, specOkLimit = 10000, dbeMinLimit = -0.5,     
                       satelliteMzLimit = 0.5, satelliteIntLimit = 0.05)
# spectraList <- list(list(
#   mode = "CID", ce = "100", ces = "100", res = 10000
# ))

msmsPeaks <- new("RmbSpectraSet")
msmsPeaks@formula <- "C12H2"
msmsPeaks@mz <- findMz.formula("C12H2", "m2H_c2")$mzCenter

mz <- c(
  120 + RMassBank:::.emass, # this can only be [C10]-
  60 + RMassBank:::.emass, # this can be [C5]- or [C10]--
  60 + 1.00783 + RMassBank:::.emass, # this can only be [C5H1]-
  60.5, # this is nothing
  143, 
  142)
msmsPeaks@children <- list(new("RmbSpectrum2",
                               mz = mz,
                               intensity = rep_along(mz, 1e6)
)) %>% as("SimpleList")
                               
results <- analyzeMsMs.formula(msmsPeaks, "m2H_c2", 
                         detail=TRUE, 
                         run="preliminary",
                         filterSettings)
data <- getData(results@children[[1]])

# Check that all formulas for all peaks are correct, see expectations above

mz60 <- data %>% filter(between(mz, 60, 60.1))
expect_equal( nrow(mz60), 2 )
expect_setequal(mz60$formula, c("C5", "C10"))
expect_true(all(mz60$dppm < 0.1))

mz605 <- data %>% filter(between(mz, 60.4, 60.6))
expect_equal( nrow(mz605), 1 )
expect_equal(mz605$formula, NA_character_)

mz61 <- data %>% filter(between(mz, 61., 61.1))
expect_equal( nrow(mz61), 1 )
expect_equal(mz61$formula, "C5H")
expect_true(all(mz61$dppm < 0.1))

mz120 <- data %>% filter(between(mz, 120, 120.1))
expect_equal( nrow(mz120), 1 )
expect_equal(mz120$formula, "C10")
expect_true(all(mz120$dppm < 0.1))

mz142 <- data %>% filter(between(mz, 141.9, 142.1))
expect_equal( nrow(mz142), 1 )
expect_equal(mz142$formula, NA_character_)

mz143 <- data %>% filter(between(mz, 142.9, 143.1))
expect_equal( nrow(mz143), 1 )
expect_equal(mz143$formula, NA_character_)

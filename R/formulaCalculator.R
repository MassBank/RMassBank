#' @import rcdk
NULL



#' Silent getformula
#' 
#' This is a local function which silences Rcdk get.mol2formula 
#' so it doesn't print() on every call. (It is not a modified version 
#' of the LGPL'ed get.mol2formula, just a call to it.)
#'   
#' @param molecule see ?get.mol2formula
#' @param charge get.mol2formula
#' @return get.mol2formula
#' 
#' @author stravsmi
#' @noRd 
.get.mol2formula <- function (molecule, charge = 0) 
{
	capture.output(f <- get.mol2formula(molecule, charge))
	return(f)
}

#
# Formula calculation. Addition and subtraction of chemical formulas.
# Examples:


#' Convert formula to Rcdk limits
#' 
#' Converts a molecular formula e.g. C15H20 into an upper limit appropriate for
#' use with Rcdk's \code{\link{generate.formula}} function's \code{element}
#' argument.
#' 
#' This helper function is used to make the upper limits for
#' \code{\link{generate.formula}} when finding subformulas to match to a MS2
#' fragment peak.
#' 
#' @usage to.limits.rcdk(formula)
#' @param formula A molecular formula in string or list representation
#' (\code{"C6H6"} or \code{list(C=6,H=6)}).
#' @param fix Fix the bug https://github.com/CDK-R/cdkr/issues/129 by inserting dummy element.
#' @return An array in the form \code{c( c("C", "0", "12"), c("H", "0", "12"))}
#' (for input of "C12H12").
#' @author Michael Stravs
#' @seealso \code{\link{generate.formula}}, \code{\link{add.formula}}
#' @examples
#' 
#' #
#' to.limits.rcdk("C6H6")
#' to.limits.rcdk(add.formula("C6H12O6", "H"))
#' 
#' @export
to.limits.rcdk <- function(formula, fix = FALSE)
{
    if(!is.list(formula))
        formula <- formulastring.to.list(formula)
    elelist <- lapply(names(formula), function(element) {
                    return(c(element, 0, formula[[element]]))
                })
    if(fix & (length(elelist) == 1)) {
      # we need to try at least two elements to make sure to add a new one
      if("H" %in% names(formula)) {
        elelist <- c(elelist, list(c("C", "0", "0")))
      } else {
        elelist <- c(elelist, list(c("H", "0", "0")))
      }
    }
    return(elelist)
}


#' Check validity of formula
#' 
#' Checks whether the formula is chemically valid, i.e. has no zero-count or
#' negative-count elements.
#' 
#' The check is only meant to identify formulas which have negative elements,
#' which can arise from the subtraction of adducts.  It is \bold{not} a
#' high-level formula "validity" check like e.g. the Rcdk function
#' \code{isvalid.formula} which uses the nitrogen rule or a DBE rule.
#' 
#' @usage is.valid.formula(formula)
#' @param formula A molecular formula in string or list representation
#' (\code{"C6H6"} or \code{list(C=6,H=6)}).
#' @author Michael Stravs
#' @seealso \code{\link{list.to.formula}}, \code{\link{add.formula}},
#' \code{\link{order.formula}}
#' @examples
#' 
#' #
#' is.valid.formula(list(C=0,H=1,Br=2))
#' is.valid.formula("CH2Cl")
#' is.valid.formula("C0H2")
#' 
#' @export
is.valid.formula <- function(formula)
{
  if(!is.list(formula))
    formula <- formulastring.to.list(formula)
  if(length(which(formula <= 0)) > 0)
    return(FALSE)
  else
    return(TRUE)
}

#' Order a chemical formula correctly
#' 
#' Orders a chemical formula in the commonly accepted order (CH followed by
#' alphabetic ordering).
#' 
#' 
#' @usage order.formula(formula, as.formula = TRUE, as.list = FALSE)
#' @param formula A molecular formula in string or list representation
#' (\code{"C6H6"} or \code{list(C=6,H=6)}).
#' @param as.formula If \code{TRUE}, the return value is returned as a string.
#' This is the default.
#' @param as.list If \code{TRUE}, the return value is returned in list
#' representation.
#' @author Michele Stravs
#' @seealso \code{\link{list.to.formula}}, \code{\link{add.formula}},
#' \code{\link{is.valid.formula}}
#' @examples
#' 
#' #
#' order.formula("H4C9")
#' order.formula("C2N5HClBr")
#' 
#' @export
order.formula <- function(formula, as.formula=TRUE, as.list=FALSE)
{
  if(!is.list(formula))
    formula <- formulastring.to.list(formula)
  result <- list()
  if(!is.null(formula$C))
    result$C <- formula$C
  if(!is.null(formula$H))
    result$H <- formula$H
  elements <- setdiff(names(formula), c("C","H"))
  elements <- sort(elements)
  for(element in elements)
    result[[element]] <- formula[[element]]
  if(!as.list | as.formula)
    return(list.to.formula(result))
  else
    return(result)
}

#' Calculate Double Bond Equivalents
#' 
#' Calculates the Ring and Double Bond Equivalents for a chemical formula. The
#' highest valence state of each atom is used, such that the returned DBE
#' should never be below 0.
#' 
#' 
#' @usage dbe(formula)
#' @param formula A molecular formula in text or list representation (e.g.
#' \code{"C6H12O6"} or \code{list(C=6, H=12, O=6)} ).
#' @return Returns the DBE for the given formula.
#' @author Michael Stravs
#' @examples
#' 
#' #
#' 	dbe("C6H12O6")
#' 
#' @export
dbe <- function(formula)
{
  if(!is.list(formula))
  {
	  if(is.na(formula))
		  return(NA)
	  formula <- formulastring.to.list(formula)
  }
  # Valences are set to the "maximum" state. This is done
  # in order to not exclude peaks from high-valence SPN atoms.
  atomDBE <- list(
    "C" = 1,
    "N" = 1.5,
    "O" = 0,
    "Si" = 1,
    "H" = -0.5,
    "F" = -0.5,
    "Cl"= -0.5,
    "Br" = -0.5,
    "S" = 2,
    "Se" = 2,
    "P" = 1.5,
    "I" = -0.5,
    "As" = 2.5,
    "Hg" = 0,
    "Li" = -0.5,
    "Na" = -0.5,
    "K" = -0.5,    
    "B" = 0.5,
    "Al" = 0.5,  
    "Ga" = 0.5,
    "Ca" = 0  
    )
  count <- 1
  for(element in names(formula))
    count <- count + atomDBE[[element]] * formula[[element]]
  return(count)
}


#' Interconvert molecular formula representations
#' 
#' Converts molecular formulas from string to list representation or vice
#' versa.
#' 
#' The function doesn't care about whether your formula makes sense. However,
#' \code{"C3.5O4"} will give \code{list("C" = 3, "O" = 4)} because regular
#' expressions are used for matching (however, \code{list("C" = 3.5, "O" = 4)}
#' gives \code{"C3.5O4"}.) Duplicate elements cause problems; only "strict"
#' molecular formulas ("CH4O", but not "CH3OH") work correctly.
#' 
#' @aliases list.to.formula formulastring.to.list
#' @usage list.to.formula(flist) 
#' 
#' formulastring.to.list(formula)
#' @param flist A molecular formula in list format, e.g. \code{list( "C" = 6,
#' "H" = 12, "O" = 6 )}.
#' @param formula A molecular formula in string format, e.g. \code{"C6H12O6"}.
#' @return \code{list.to.formula} returns a string representation of the
#' formula; \code{formulastring.to.list} returns the list representation.
#' @author Michael Stravs
#' @seealso \code{\link{add.formula}}, \code{\link{order.formula}},
#' \code{\link{is.valid.formula}}
#' @examples
#' 
#' #
#' 	list.to.formula(list("C" = 4, "H" = 12))
#' 	# This is also OK and useful to calculate e.g. adducts or losses.
#' 	list.to.formula(list("C" = 4, "H" = -1))
#' 	formulastring.to.list(list.to.formula(formulastring.to.list("CHIBr")))
#' 
#' @export
formulastring.to.list <- function(formula)
{
  matches.list <- gregexpr("([A-Z][a-z]*)([-0-9]*)", formula, perl=TRUE)
  matches <- regmatches(formula, matches.list )[[1]]
  
  flist <- list()
  for(match in matches)
  {
    match.element <- sub("([A-Z][a-z]*)([-0-9]*)", "\\1", match)
    match.count <- as.integer(sub("([A-Z][a-z]*)([-0-9]*)", "\\2", match))
    if(!is.na(match.count))
      flist[[match.element]] <- match.count
    else
      flist[[match.element]] <- 1
  }
  return(flist)
  # " (.*?) - (-*[0-9]+), "
}

#' @export
list.to.formula <- function(flist)
{
  formula <- ""
  for(element in names(flist))
    formula <- paste(formula, element, flist[[element]], sep='')
  return(formula)
}

#' Calculations on molecular formulas
#' 
#' Add, subtract, and multiply molecular formulas.
#' 
#' Note that the results are not checked for plausibility at any stage, nor
#' reordered.
#' 
#' @aliases add.formula multiply.formula
#' @usage add.formula(f1, f2, as.formula = TRUE, as.list = FALSE)
#' multiply.formula(f1, n, as.formula = TRUE, as.list = FALSE)
#' @param f1,f2 Molecular formulas (in list form or in text form) to calculate
#' with.
#' @param n Multiplier (positive or negative, integer or non-integer.)
#' @param as.formula Return the result as a text formula (e.g.
#' \code{"C6H12O6"}). This is the default
#' @param as.list Return the result in list format (e.g. \code{list(C=6, H=12,
#' O=6)}).
#' @return The resulting formula, as specified above.
#' @author Michael Stravs
#' @seealso \code{\link{formulastring.to.list}}, \code{\link{is.valid.formula}},
#' \code{\link{order.formula}}
#' @examples
#' 
#' ##
#' 
#' add.formula("C6H12O6", "C3H3")
#' add.formula("C6H12O6", "C-3H-3")
#' add.formula("C6H12O6", multiply.formula("C3H3", -1))
#' 
#' @export
add.formula <- function(f1, f2, as.formula = TRUE, as.list=FALSE)
{
  ret = list()
  if(!is.list(f1)) f1 <- formulastring.to.list(f1)
  if(!is.list(f2)) f2 <- formulastring.to.list(f2)
  add <- intersect(names(f1), names(f2))
  for(element in add)
    ret[[element]] <- f1[[element]] + f2[[element]]
  e_f1 <- setdiff(names(f1), names(f2))
  e_f2 <- setdiff(names(f2), names(f1))
  for(element in e_f1)
    ret[[element]] <- f1[[element]]
  for(element in e_f2)
    ret[[element]] <- f2[[element]]
  
  # eliminate all 0-elements
  ret <- ret[which(ret != 0)]
  
  if(as.formula & !as.list)
    return(list.to.formula(ret))
  else
    return(ret)
}

#' @export
multiply.formula <- function(f1, n, as.formula=TRUE, as.list=FALSE)
{
  if(!is.list(f1)) f1 <- formulastring.to.list(f1)
  ret <- lapply(f1, function(element) n*element)
  if(as.formula & !as.list)
    return(list.to.formula(ret))
  else
    return(ret)
}

#' Calculate ppm values
#' 
#' Calculates ppm values for a given mass.
#' 
#' This is a helper function used in RMassBank code.
#' 
#' @param mass The "real" mass
#' @param dppm The mass deviation to calculate
#' @param l Boolean: return limits? Defaults to FALSE.
#' @param p Boolean: return ppm error itself? Defaults to FALSE.
#' @return By default (\code{l=FALSE, p=FALSE}) the function returns the mass plus the 
#' ppm error (for 123.00000 and 10 ppm: 123.00123, or for 123 and -10 ppm: 
#' 122.99877).
#' 
#' For \code{l=TRUE}, the function returns the upper and lower limit (sic!)
#' For \code{p=TRUE}, just the difference itself is returned (0.00123 for 123/10ppm).
#' @examples ppm(100, 10)
#' @author Michael A. Stravs, Eawag <michael.stravs@@eawag.ch>
#' @export
ppm <- function(mass, dppm, l=FALSE, p=FALSE)
{
    if(p) return(mass*dppm*1e-6)
    dmass <- mass * (1 + dppm*1e-6)
    if(l) dmass <- c(dmass, mass * (1 - dppm*1e-6))
    return(dmass)
}

## # auxiliaries
.emass <- 0.0005485799
## pmass <- 1.007276565
## hmass <- 1.007825


split.formula.posneg <- function(f, as.formula = TRUE, as.list=FALSE)
{
	if(!is.list(f)) f <- formulastring.to.list(f)
	pos <- f[which(f > 0)]
	neg <- f[which(f < 0)]
	neg <- multiply.formula(neg, -1, as.list = TRUE)
	if(as.formula & !as.list)
		return(list(pos=list.to.formula(pos), neg=list.to.formula(neg)))
	else
		return(list(pos=pos, neg=neg))
}



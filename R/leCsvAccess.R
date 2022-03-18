
.listEnvEnv <- new.env()
assign("listEnv", NULL, envir=.listEnvEnv)

#' Load compound list for RMassBank
#' 
#' Loads a CSV compound list with compound IDs
#' 
#' The list is loaded into the variable \code{\var{compoundList}} in the environment
#' \code{listEnv} (which defaults to the global environment) and used by
#' the \code{findMz}, \code{findCAS}, ... functions. The CSV file is required to have at least the following columns, which are used for 
#' further processing and must be named correctly (but present in any order): \var{ID, Name, SMILES, RT,
#' CAS}
#' 
#' resetList() clears a currently loaded list.
#' 
#' @aliases loadList resetList
#' @usage loadList(path, listEnv=NULL, check=TRUE)
#' 
#' resetList()
#' @param path Path to the CSV list.
#' @param listEnv The environment to load the list into. By default, the namelist is loaded
#' 		into an environment internally in RMassBank. 
#' @param check A parameter that specifies whether the SMILES-Codes in the list should be checked for readability by rcdk.
#' @return No return value.
#' @author Michael Stravs
#' @seealso \code{\link{findMz}}
#' @examples
#' 
#' ##
#' \dontrun{loadList("mylist.csv")}
#' 
#' @export
loadList <- function(path, listEnv = NULL, check = TRUE)
{
  	if(is.null(listEnv))
  		listEnv <- .listEnvEnv
  	if(!file.exists(path))
  		stop("The supplied file does not exist, please supply a correct path")
          
    # Try out if the file is comma- or semicolon-separated
    compoundList <- read.csv(path, stringsAsFactors=FALSE, check.names=FALSE)
  	n <- colnames(compoundList)
    if(!("ID" %in% n)){ # If no ID column, it must be semicolon separated
        compoundList <- read.csv2(path, stringsAsFactors=FALSE, check.names=FALSE)
        n <- colnames(compoundList)
        if(!("ID" %in% n)){ # ...or there must be something wrong with the column names
             stop("There is no 'ID' column in the compound list")
        }
    }
    
    # Now everything should be fine, at least regarding csv and ssv
    
    # Are there duplicate compound IDs? 
    if(any(duplicated(compoundList$ID))){
        stop("Duplicate compound IDs are present. Please check the compound list.")
    }
    
    # Do all Compoundlist IDs have 4 characters or less?
    if(any(nchar(compoundList$ID) > 4)){
        stop("The maximum number of digits for compound IDs is 4")
    }
    
    # Evaluate if all strictly necessary columns are in the list
    cols <- c('ID', 'Name', 'SMILES', 'RT', 'CAS')
    d <- setdiff(cols, n)
    if(length(d)>0){
        stop(paste("Some columns are missing in the compound list. It needs at least", paste(cols,collapse=", ")))
    }
    
    # Is there a "Level" column?
    if("Level" %in% n){
        newList <- TRUE
    } else {
        newList <- FALSE
    }
    
    # Assign for now...
    assign("listEnv", listEnv, envir=.listEnvEnv) 
    .listEnvEnv$listEnv$compoundList <- compoundList
    # If "level" is in the compound list we have to check several things:
    
    if(newList){
        # a) Are the levels part of the defined levels?
        # b) Are the values ok for every level? (i.e. all necessary values supplied for each line in the compound list?)
        
        # Check a) (and translate to number levels because handling numbers and text would make the if-statements here even more confusing)
        
        compoundList$Level <- sapply(as.character(compoundList$Level),.translateLevel)
        
        # Need this variable for levels 3b, 3c, 3d, 4 and potentially 3 and 3a
        formulaColPresent <- "Formula" %in% n
        
        if(any(c("3b","3c","3d","4") %in% compoundList$Level)){
    
            if(!(formulaColPresent)){
                .listEnvEnv$listEnv$compoundList <- NULL
                stop('The compound list must contain the column "Formula" if a tentative or formula compound (levels 3b, 3c, 3d, 4) is present')
            }
        }
        # Need this variable for level 5
        if("5" %in% compoundList$Level){
            mzCol <- which(c("mz","m/z","mass","exactMass") %in% n)
            if(length(mzCol) > 1){
                .listEnvEnv$listEnv$compoundList <- NULL
                stop("The compound list can only contain one column with one of the following column names: 'mz','m/z','mass','exactMass' if an unknown compound (level 5) is present")
            }
        
            if(mzCol == 2){
                    n[which(n == "m/z")] <- "mz"
            }
            
            if(mzCol == 4){
                    n[which(n == "exactMass")] <- "mass"
            }
            .listEnvEnv$listEnv$mzCol <- c("mz","m/z","mass","exactMass")[mzCol]
            colnames(compoundList) <- n
        }
        # Check b) 
        for(i in 1:length(compoundList$Level)){
            level <- compoundList$Level[i] 
            # ID must have numbers as values
            if(!is.numeric(compoundList[i,"ID"])){
                .listEnvEnv$listEnv$compoundList <- NULL
                stop(paste("Value for 'ID' in line", i, "of the compound list is not a valid compound ID"))
            }
            
            # If level is "1x" or "2x", the SMILES must be supplied and must be correct
            if(level %in% c("0","1","1a","1b","1c","2","2a","2b")){
                currEnvir <- environment()
                
                tryCatch(
                    findMz(compoundList[i,"ID"]),
                    error = function(e){
                        .listEnvEnv$listEnv$compoundList <- NULL
                        stop(paste("'SMILES' value for compound", compoundList[i,"ID"] ,"in line", i, "of the compound list is not a valid SMILES"))
                    }
                )
            }
            
            # If level is "3" or "3a", a valid smiles or formula must be supplied
            if(level %in% c("3","3a")){
  
                if(!is.na(findSmiles(compoundList[i,"ID"]))){
                    tryCatch(
                        findMz(compoundList[i,"ID"]),
                        error = function(e){
                            .listEnvEnv$listEnv$compoundList <- NULL
                            stop(paste("'SMILES' value for compound", compoundList[i,"ID"] ,"in line", i, "of the compound list is not a valid SMILES"))
                        }
                    )
                } else{
                    if(!formulaColPresent){
                        .listEnvEnv$listEnv$compoundList <- NULL
                        stop(paste("Compound", compoundList[i,"ID"], "in line", i, "of the compound list is marked as level 3, but there is no valid SMILES
                            value nor a 'Formula' column present"))
                    }
                    
                    tryCatch(
                        rcdk::get.formula(findFormula(compoundList[i,"ID"], retrieval="tentative")),
                        error = function(e){
                            .listEnvEnv$listEnv$compoundList <- NULL
                            stop(paste("'Formula' value for compound", compoundList[i,"ID"] ,"in line", i, "of the compound list is not a valid Formula"))
                        }
                    )
                }
            }
            
            # If level is "3b","3c","3d" or "4", a valid formula must be supplied
            if(level %in% c("3b","3c","3d","4")){
                
                tryCatch(
                    rcdk::get.formula(findFormula(compoundList[i,"ID"], retrieval="tentative")),
                    error = function(e){
                        .listEnvEnv$listEnv$compoundList <- NULL
                        stop(paste("'Formula' value for compound", compoundList[i,"ID"] ,"in line", i, "of the compound list is not a valid Formula"))
                    }
                )
            }
            # If level is "5", m/z must be supplied
           
            if(level == "5"){
                if(!is.numeric(findMz(compoundList[i,"ID"], retrieval="unknown")$mzCenter)){
                    .listEnvEnv$listEnv$compoundList <- NULL
                    stop(paste(c("mz","m/z","mass","exactMass")[mzCol],"value for compound", compoundList[i,"ID"] ,"in line", i, "of the compound list is not a valid number"))
                }
            }
        }
    }
    
    # If "Level" is not in the compound list it MUST be a standard list, so process just as before:
    if(!newList){
        cols <- c('ID', 'Name', 'SMILES', 'RT', 'CAS')
        d <- setdiff(cols, n)
        if(length(d)>0){
            stop("Some columns are missing in the compound list. Needs at least ID, Name, SMILES, RT, CAS.")
        }
  
        ###
        ###Test if all the IDs work
        ###
        
        if(check){
            wrongID <- vector()
            currEnvir <- environment()
            sapply(compoundList[,"ID"], function(x){
                tryCatch(
                    findMz(x),
                    error = function(e){
                      if(RMassBank.env$verbose.output)
                        cat(paste("### Warning ### Error finding SMILES for ID '", x, "': ", e, sep = ""))
                      
                      currEnvir$wrongID <- c(currEnvir$wrongID, x)
                    }
                )
            })
            if(length(wrongID)){
                .listEnvEnv$listEnv$compoundList <- NULL
                stop(paste("Unable to interpret the SMILES-strings for ID(s)", paste(wrongID, collapse= " "), "\nPlease check these and load the list again."))
            }
        }
        Level <- rep("0",nrow(compoundList))
        .listEnvEnv$listEnv$compoundList <- cbind(compoundList,Level)
    }
    rmb_log_info("Loaded compoundlist successfully")
}

#' @export
resetList <- function()
{
	if(is.null(.listEnvEnv$listEnv))
		return()
	if(!exists("compoundList", where=.listEnvEnv$listEnv))
		return()
	rm(.listEnvEnv$listEnv$compoundList)
	assign("listEnv", NULL, envir=.listEnvEnv)
}

# Function that translates one entry of the level column of the compound list into the number system
.translateLevel <- function(level){
    if(!(level %in% c("0","1","1a","1b","1c","2","2a","2b","3","3a","3b","3c","3d","4","5"))){
        switch(level,
            standard={
                return("1a")
            },
            parent={
                return("1b")
            },
            confirmed={
                return("1c")
            },
            probable={
                return("2")
            },
            probableLibrary={
                return("2a")
            },
            probableDiagnostic={
                return("2b")
            },
            tentative={
                return("3")
            },
            tentativeStructure={
                return("3a")
            },
            tentativeIsomer={
                return("3b")
            },
            tentativeTPClass={
                return("3c")
            },
            tentativeBestMatch={
                return("3d")
            },
            formula={
                return("4")
            },
            unknown={
                return("5")
            },
            exactMass={
                return("5")
            },
            {
                stop(paste(level, "is not a valid level"))
            }
        )
    } else{
        return(level)
    }
}

#' Create Rcdk molecule from SMILES
#' 
#' Generates a Rcdk molecule object from SMILES code, which is fully typed and
#' usable (in contrast to the built-in \code{parse.smiles}).
#' 
#' \bold{NOTE: As of today (2012-03-16), Rcdk discards stereochemistry when
#' loading the SMILES code!} Therefore, do not trust this function blindly,
#' e.g.  don't generate InChI keys from the result. It is, however, useful if
#' you want to compute the mass (or something else) with Rcdk.
#' 
#' @usage getMolecule(smiles)
#' @param smiles The SMILES code of the compound.
#' @return A Rcdk \code{IAtomContainer} reference.
#' @author Michael Stravs
#' @seealso \code{\link{parse.smiles}}
#' @examples
#' 
#' # Lindane:
#' getMolecule("C1(C(C(C(C(C1Cl)Cl)Cl)Cl)Cl)Cl")
#' # Benzene:
#' getMolecule("C1=CC=CC=C1")
#' 
#' @export
getMolecule <- function(smiles)
{
  capture.output(mol <- parse.smiles(smiles)[[1]])
  do.aromaticity(mol)
  convert.implicit.to.explicit(mol)
  do.aromaticity(mol)
  set.atom.types(mol)
  do.isotopes(mol)
 
  return(mol)
}

knownAdducts <- function(){
  return(getAdductInformation("")$mode)
}
getMonoisotopicMass <- function(formula){
  if(!exists("isotopes")) data("isotopes", package = "enviPat")
  
  if(formula == "") return(0)
  
  if(grepl(x = formula, pattern = "-")){
    starts <- gregexpr(text = formula, pattern = "[A-Z]")[[1]]
    subFormulas <- sapply(X = seq_along(starts), FUN = function(startIdx){
      ifelse(
        test = startIdx < length(starts), 
        yes = substr(x = formula, start = starts[[startIdx]], stop = starts[[startIdx + 1]] - 1), 
        no  = substr(x = formula, start = starts[[startIdx]], stop = nchar(formula))
      )
    })
    
    monoisotopicMass <- sum(sapply(X = subFormulas, FUN = function(subFormula){
      ifelse(
        test = grepl(x = subFormula, pattern = "-"), 
        yes = -enviPat::isopattern(isotopes = isotopes, chemforms = gsub(x = subFormula, pattern = "-", replacement = ""), threshold=0.1, charge = FALSE, verbose = FALSE)[[1]][[1,1]], 
        no  =  enviPat::isopattern(isotopes = isotopes, chemforms = subFormula,                                            threshold=0.1, charge = FALSE, verbose = FALSE)[[1]][[1,1]]
      )
    }))
  } else {
    monoisotopicMass <- enviPat::isopattern(isotopes = isotopes, chemforms = formula, threshold=0.1, charge = FALSE, verbose = FALSE)[[1]][[1,1]]
  }
  return(monoisotopicMass)
}

getAdductPolarity <- function(mode) {
  df <- getAdductInformation("")
  charge <- df[df$mode == mode,"charge"]
  ifelse(charge > 0, 1L, 0L)
}

getIonMode <- function(mode) {
  df <- getAdductInformation("")
  charge <- df[df$mode == mode,"charge"]
  ifelse(charge > 0, "POSITIVE", "NEGATIVE")
}

getAdductInformation <- function(formula){
  adductDf <- as.data.frame(rbind(
    
    ## positive: M+X
    c(mode = "pH",       addition = "H",         charge = 1, adductString = "[M+H]+"),
    c(mode = "pLi",      addition = "Li",        charge = 1, adductString = "[M+Li]+"),
    c(mode = "pNa",      addition = "Na",        charge = 1, adductString = "[M+Na]+"),
    c(mode = "pNa_mO3S_mH",      addition = "Na1O-3S-1H-1",        charge = 1, adductString = "[M-O3S-H+Na]+"),
    c(mode = "pK",       addition = "K",         charge = 1, adductString = "[M+K]+"),
    c(mode = "pM",       addition = "",          charge = 1, adductString = "[M]+"),
    c(mode = "pM_mC7H11NO9S2",       addition = "C-7H-11NO-9S-2",          charge = 1, adductString = "[M-C7H11NO9S2]+"),
    c(mode = "pM_mC7H12NO9S2",       addition = "C-7H-12NO-9S-2",          charge = 1, adductString = "[M-C7H12NO9S2]+"),
    c(mode = "pNH4",     addition = "NH4",       charge = 1, adductString = "[M+NH4]+"),
    c(mode = "p2Na_mH",  addition = "Na2H-1",    charge = 1, adductString = "[M+2Na-H]+"),
    c(mode = "pACN_pH",  addition = "C2H4N1",    charge = 1, adductString = "[M+ACN+H]+"),
    c(mode = "pACN_pNa", addition = "C2H3N1Na1", charge = 1, adductString = "[M+ACN+Na]+"),
    c(mode = "pH_mC7H6O",     addition = "C-7H-5O-1", charge = 1, adductString = "[M-C7H6O+H]+"),
    c(mode = "pH_mC18H30O14", addition = "C-18H-29O-14", charge = 1, adductString = "[M-C18H30O14+H]+"),
    c(mode = "pH_mC6H10O5",   addition = "C-6H-9O-5", charge = 1, adductString = "[M-C6H10O5+H]+"),
    c(mode = "pH_mC12H20O9",  addition = "C-12H-19O-9", charge = 1, adductString = "[M-C12H20O9+H]+"),
    c(mode = "pH_mC9H8O4_mH2O",  addition = "C-9H-9O-5", charge = 1, adductString = "[M-C9H8O4-H2O+H]+"),
    c(mode = "pH_mC6H10O5_mH2O", addition = "C-6H-11O-6", charge = 1, adductString = "[M-C6H10O5-H2O+H]+"),
    c(mode = "pH_mC5H8NO4",   addition = "C-5H-7N-1O-4", charge = 1, adductString = "[M-C5H8NO4+H]+"),
    c(mode = "pH_mO3S",       addition = "O-3S-1H1", charge = 1, adductString = "[M-O3S+H]+"),
    c(mode = "pH_mC6H10O8S",  addition = "C-6H-9O-8S-1", charge = 1, adductString = "[M-C6H10O8S+H]+"),
    c(mode = "pH_mC5H10N2O",  addition = "C-5H-9N-2O-1", charge = 1, adductString = "[M-C5H10N2O+H]+"),
    c(mode = "pH_mHO3P",      addition = "O-3P-1", charge = 1, adductString = "[M-HO3P+H]+"),
    c(mode = "pH_mC4H7",      addition = "C-4H-6", charge = 1, adductString = "[M-C4H7+H]+"),
    c(mode = "pH_mC6H10O4",   addition = "C-6H-9O-4", charge = 1, adductString = "[M-C6H10O4+H]+"),
    c(mode = "pH_mC5H8O3",    addition = "C-5H-7O-3", charge = 1, adductString = "[M-C5H8O3+H]+"),
    c(mode = "pH_mCO",        addition = "H1C-1O-1", charge = 1, adductString = "[M-CO+H]+"),
    c(mode = "p_mCO",         addition = "C-1O-1", charge = 1, adductString = "[M-CO]+"),
    c(mode = "pH_mO3",        addition = "H-1O-3", charge = 1, adductString = "[M-O3+H]+"),
    c(mode = "pH_mC3H6",      addition = "C-3H-5", charge = 1, adductString = "[M-C3H6+H]+"),
    c(mode = "pH_mC4H3O5",    addition = "C-4H-2O-5", charge = 1, adductString = "[M-C4H3O5+H]+"),
    c(mode = "pH_mC6H11O6",   addition = "C-6H-10O-6", charge = 1, adductString = "[M-C6H11O6+H]+"),
    c(mode = "pH_mCH4S",     addition = "C-1H-3S-1", charge = 1, adductString = "[M-CH4S+H]+"),
    c(mode = "pH_mC7H12O6",  addition = "C-7H-11O-6", charge = 1, adductString = "[M-C7H12O6+H]+"),
    c(mode = "pH_mCH4O",     addition = "C-1H-3O-1", charge = 1, adductString = "[M-CH4O+H]+"),
    c(mode = "pH_mCH2O2",    addition = "C-1H-1O-2", charge = 1, adductString = "[M-CH2O2+H]+"),
    c(mode = "pH_mC4H8",     addition = "C-4H-7", charge = 1, adductString = "[M-C4H8+H]+"),
    c(mode = "pH_mC3H6O",    addition = "C-3H-5O-1", charge = 1, adductString = "[M-C3H6O+H]+"),
    c(mode = "pH_mC8H18O2",  addition = "C-8H-17O-2", charge = 1, adductString = "[M-C8H18O2+H]+"),
    c(mode = "pH_mC6H14O2",  addition = "C-6H-13O-2", charge = 1, adductString = "[M-C6H14O2+H]+"),
    c(mode = "pH_mC4H12O2",  addition = "C-4H-11O-2", charge = 1, adductString = "[M-C4H12O2+H]+"),
    c(mode = "pH_mH2O",  addition = "H-1O-1",    charge = 1, adductString = "[M-H2O+H]+"),
    c(mode = "pNa_mH2O", addition = "H-2O-1Na1", charge = 1, adductString = "[M-H2O+Na]+"),
    c(mode = "pH_mCO2",  addition = "C-1O-2H1", charge = 1, adductString = "[M-CO2+H]+"),
    c(mode = "pH_mO",  addition = "O-1H1", charge = 1, adductString = "[M-O+H]+"),
    c(mode = "p_mO",  addition = "O-1", charge = 1, adductString = "[M-O]+"),
    c(mode = "p2H",      addition = "H2",        charge = 2, adductString = "[M+2H]2+"),
    c(mode = "pACN_p2H", addition = "C2H5N1",    charge = 2, adductString = "[M+ACN+2H]2+"),
    ## positive: 2M+X
    c(mode = "pM_pH",      addition = add.formula(formula, "H1"),     charge = 1, adductString = "[2M+H]+"),
    c(mode = "pM_pK",      addition = add.formula(formula, "K1"),     charge = 1, adductString = "[2M+K]+"),
    c(mode = "pM_pNa",     addition = add.formula(formula, "Na1"),    charge = 1, adductString = "[2M+Na]+"),
    c(mode = "pM_pNH4",    addition = add.formula(formula, "N1H4"),   charge = 1, adductString = "[2M+NH4]+"),
    c(mode = "pM_pACN_pH", addition = add.formula(formula, "C2H4N1"), charge = 1, adductString = "[2M+ACN+H]+"),
    ## positive: strange positive adducts
    c(mode = "pCOONa",         addition = "C1O2Na1", charge =  1, adductString = "[M+COONa]+"),
    c(mode = "p3H_c1",         addition = "H3",      charge =  1, adductString = "[M+3H]+"),
    c(mode = "pH2O_c1",        addition = "H2O1",    charge =  1, adductString = "[M+H2O]+"),
    c(mode = "pH_m2H2O",       addition = "H-3O-2",  charge =  1, adductString = "[M-2H2O+H]+"),
    c(mode = "pH_pH2O",        addition = "H3O1",    charge =  1, adductString = "[M+H2O+H]+"),
    c(mode = "pH_mNH3",        addition = "N-1H-2",  charge =  1, adductString = "[M-NH3+H]+"),
    c(mode = "p2H_c1",         addition = "H2",      charge =  1, adductString = "[M+2H]+"),
    c(mode = "p_mNH3_c1",      addition = "N-1H-3",  charge =  1, adductString = "[M-NH2-H]+"),
    c(mode = "p_mNH2_pH_c1",   addition = "N-1H-1",  charge =  1, adductString = "[M-NH2+H]+"),
    c(mode = "pM_p2Na_m3H_c1", addition = add.formula(formula, "Na2H-3"), charge =  1, adductString = "[2M+2Na-3H]+"),
    c(mode = "pM_pNa_m2H_c1",  addition = add.formula(formula, "Na1H-2"), charge =  1, adductString = "[2M+Na-2H]+"),
    c(mode = "pM_pNa_mH_c1",   addition = add.formula(formula, "Na1H-1"), charge =  1, adductString = "[2M+Na-H]+"),
    c(mode = "pM_p2Na_m2H_c1", addition = add.formula(formula, "Na2H-2"), charge =  1, adductString = "[2M+2Na-2H]+"),
    c(mode = "pM_pH_m2H2O_c1", addition = add.formula(formula, "H-3O-2"), charge =  1, adductString = "[2M-2H2O+H]+"),
    c(mode = "pM_pH_mH2O",     addition = add.formula(formula, "H-1O-1"), charge =  1, adductString = "[2M-H2O+H]+"),
    c(mode = "pM_pNa_mH2O",    addition = add.formula(formula, "H-2O-1Na1"), charge =  1, adductString = "[2M-H2O+Na]+"),
    c(mode = "pM_m2H_c1",      addition = add.formula(formula, "H-2"),    charge =  1, adductString = "[2M-2H]+"),
    c(mode = "pM_mH_c2",       addition = add.formula(formula, "H-1"),    charge =  1, adductString = "[2M-2H+H]+"),
    c(mode = "pM_pLi",         addition = add.formula(formula, "Li1"),    charge =  1, adductString = "[2M+Li]+"),
    c(mode = "pM_pH_m2O",      addition = add.formula(formula, "O-2H1"),  charge =  1, adductString = "[2M-2O+H]+"),
    c(mode = "pM_pNa_m2O",     addition = add.formula(formula, "O-2Na1"), charge =  1, adductString = "[2M-2O+Na]+"),
    c(mode = "pM_pH_m3O",      addition = add.formula(formula, "O-3H1"),  charge =  1, adductString = "[2M-3O+H]+"),
    c(mode = "pM_pNa_m3O",     addition = add.formula(formula, "O-3Na1"), charge =  1, adductString = "[2M-3O+Na]+"),
    c(mode = "pM_mH_c1",       addition = add.formula(formula, "H-1"),    charge =  1, adductString = "[2M-H]+"),
    c(mode = "p_p2M_m3H",  addition = add.formula(formula, add.formula(formula, "H-3")), charge = 1, adductString = "[3M-3H]+"),
    c(mode = "pH_p2M_m2H2O",  addition = add.formula(formula, add.formula(formula, "H-3O-2")), charge = 1, adductString = "[3M-2H2O+H]+"),
    c(mode = "pNa_p2M_m2H2O",  addition = add.formula(formula, add.formula(formula, "H-4O-2Na1")), charge = 1, adductString = "[3M-2H2O+Na]+"),
    c(mode = "p_p2M_m2H2O",  addition = add.formula(formula, add.formula(formula, "H-4O-2")), charge = 1, adductString = "[3M-2H2O]+"),
    c(mode = "pM_mH_pH",       addition = formula,   charge =  1, adductString = "[2M-H+H]+"),
    c(mode = "pH_c2",          addition = "H1",      charge =  2, adductString = "[M+H]2+"),
    
    
    ## negative: M-X
    c(mode = "mH",      addition = "H-1",    charge = -1, adductString = "[M-H]-"),
    c(mode = "mCl",     addition = "Cl1",    charge = -1, adductString = "[M+Cl]-"),
    c(mode = "mFA",     addition = "C1O2H",  charge = -1, adductString = "[M+HCOOH-H]-"),
    c(mode = "mH_pTFA", addition = "C2F3O2", charge = -1, adductString = "[M+CF3CO2H-H]-"),
    
    c(mode = "mH_mC6H10O5", addition = "C-6H-11O-5", charge = -1, adductString = "[M-C6H10O5-H]-"),
    
    c(mode = "mFA_pH",  addition = "C1O2H2", charge = -1, adductString = "[M+HCOOH]-"),
    c(mode = "mH_mH2O", addition = "H-3O-1", charge = -1, adductString = "[M-H2O-H]-"),
    c(mode = "mCO2",    addition = "C-1O-2", charge = -1, adductString = "[M-CO2]-"),
    c(mode = "mH_mCH3", addition = "C-1H-4", charge = -1, adductString = "[M-CH3-H]-"),
    c(mode = "mH_mCO2", addition = "C-1H-1O-2", charge = -1, adductString = "[M-CO2-H]-"),
    c(mode = "mCH3",    addition = "C-1H-3", charge = -1, adductString = "[M-CH3]-"),
    c(mode = "m2H_pNa", addition = "H-2Na1", charge = -1, adductString = "[M+Na-2H]-"),
    c(mode = "mM",      addition = "",       charge = -1, adductString = "[M]-"),
    c(mode = "m2H",     addition = "H-2",    charge = -1, adductString = "[M-2H]-"), ## in case of positively charged compounds
    c(mode = "m2H_c2",  addition = "H-2",    charge = -2, adductString = "[M-2H]2-"),
    ## negative: 2M-X
    c(mode = "mH_pM",      addition = add.formula(formula, "H-1"),    charge = -1, adductString = "[2M-H]-"),
    c(mode = "mFA_pM",     addition = add.formula(formula, "C1O2H"),  charge = -1, adductString = "[2M+HCOOH-H]-"),
    c(mode = "mH_pM_mH2O", addition = add.formula(formula, "H-3O-1"), charge = -1, adductString = "[2M-H2O-H]-"),
    c(mode = "m2H_pM_pNa", addition = add.formula(formula, "H-2Na1"), charge = -1, adductString = "[2M+Na-2H]-"),
    ## negative: strange adducts
    c(mode = "mpM",            addition = formula,    charge = -1, adductString = "[2M]-"),
    c(mode = "m2H_pHCOOH_pNa", addition = "Na1C1O2",  charge = -1, adductString = "[M+HCOOH+Na-2H]-"),
    c(mode = "mH_p2H",         addition = "H2",       charge = -1, adductString = "[M+3H-H]-"),
    c(mode = "mH_pH",          addition = "H1",       charge = -1, adductString = "[M+2H-H]-"),
    c(mode = "mH_pH2O",        addition = "H1O1",     charge = -1, adductString = "[M+H2O-H]-"),
    c(mode = "m4H_pM_p3Na",    addition = add.formula(formula, "Na3H-4"),    charge = -1, adductString = "[2M+3Na-4H]-"),
    c(mode = "m2H_mNH3_pNa",   addition = add.formula(formula, "Na1N-1H-5"), charge = -1, adductString = "[2M-NH3+Na-2H]-"),
    c(mode = "m3H_pM_p2Na",    addition = add.formula(formula, "Na2H-3"),    charge = -1, adductString = "[2M+2Na-3H]-"),
    c(mode = "m3H_pM",         addition = add.formula(formula, "H-3"),       charge = -1, adductString = "[2M-3H]-"),
    c(mode = "mH_p2M",         addition = add.formula(formula, add.formula(formula, "H-1")), charge = -1, adductString = "[3M-H]-"),
    c(mode = "mAc", addition = "C2O2H3",  charge = -1, adductString = "[M+CH3COO]-"),
    
    ## ???
    c(mode = "",        addition = "",       charge = 0,  adductString = "[M]")
  ), stringsAsFactors = F)
  adductDf$charge <- as.integer(adductDf$charge)
  
  
  
  if(any(any(duplicated(adductDf$mode)), any(duplicated(adductDf$adductString)))) stop("Invalid adduct table")
  
  return(adductDf)
}
getAdductProperties <- function(mode, formula){
  if(grepl(x = mode, pattern = "pM") & is.null(formula))
    stop("Cannot calculate pM adduct without formula")
  else if(is.null(formula)) formula <- ""
  
  adductDf <- getAdductInformation(formula)
  
  if(!(mode %in% adductDf$mode))
    stop("mode = \"", mode, "\" not defined")
  
  mzopt <- as.list(adductDf[adductDf$mode==mode,])
  return(mzopt)
}

#' Find the exact mass +/- a given margin for a given formula or its ions and adducts.
#' 
#' @param formula The molecular formula  in text or list format (see \code{\link{formulastring.to.list}}
#' @param mode \code{"pH", "pNa", "pM", "pNH4", "mH", "mM", "mFA"} for different ions 
#' 			([M+H]+, [M+Na]+, [M]+, [M+NH4]+, [M-H]-, [M]-, [M+FA]-). "" for the uncharged molecule.
#' @param ppm The ppm margin to add/subtract
#' @param deltaMz The absolute mass to add/subtract. Cumulative with \code{ppm}
#' @return A \code{list(mzMin=, mzCenter=, mzMax=)} with the masses.
#' @author Michael A. Stravs, Eawag <michael.stravs@@eawag.ch>
#' @examples findMz.formula("C6H6")
#' @seealso \code{\link{findMz}}
#' @export
findMz.formula <- function(formula, mode="pH", ppm=10, deltaMz=0) 
{
	if (!any(mode %in% knownAdducts())) 
		stop(paste("The ionization mode", mode, "is unknown."))
  mzopt <- getAdductProperties(mode, formula)
	formula <- add.formula(formula, mzopt$addition)
	# Since in special cases we want to use this with negative and zero number of atoms, we account for this case
	# by splitting up the formula into positive and negative atom counts (this eliminates the zeroes.)
	# Note: the previous implementation was incorrect, since 
	formula.split <- split.formula.posneg(formula)
	m <- 0
	if(formula.split$pos != "")
	{
		formula.pos <- get.formula(formula.split$pos, charge = 0)
		m = m + formula.pos@mass
	}
	if(formula.split$neg != "")
	{
		formula.neg <- get.formula(formula.split$neg, charge = 0)
		m = m - formula.neg@mass
	}
	m <- m + get.formula("H", charge = mzopt$charge)@mass - get.formula("H", charge = 0)@mass
	
	# get.formula only takes "charge" into account to add the electrons - not to 
	# divide by z to get m/z. therefore, we do it ourselves
	if(mzopt$charge != 0)
	  m <- m / abs(mzopt$charge)
	# Note: technically there is no m/z for charge=0
	delta <- ppm(m, ppm, l = TRUE)
	return(list(mzMin = delta[[2]] - deltaMz, mzMax = delta[[1]] + 
							deltaMz, mzCenter = m))
}

#' Find compound information
#' 
#' Retrieves compound information from the loaded compound list or calculates
#' it from the SMILES code in the list.
#' 
#' @aliases findMz findSmiles findFormula findRt findCAS findName findLevel
#' @usage  findMz(cpdID, mode = "pH", ppm = 10, deltaMz = 0, retrieval="standard")
#' 
#' findRt(cpdID) 
#' 
#' findSmiles(cpdID) 
#' 
#' findFormula(cpdID, retrieval="standard") 
#' 
#' findCAS(cpdID)
#' 
#' findName(cpdID)
#'
#' findLevel(cpdID, compact=FALSE)
#' @param cpdID The compound ID in the compound list.
#' @param mode Specifies the species of the molecule: An empty string specifies
#' uncharged monoisotopic mass, \code{\var{pH}} (positive H) specifies [M+H]+,
#' \code{\var{pNa}} specifies [M+Na]+, \code{\var{pM}} specifies [M]+,
#' \code{\var{mH}} and \code{\var{mFA}} specify [M-H]- and [M+FA]-,
#' respectively. (I apologize for the naming of \code{\var{pH}} which has
#' absolutely nothing to do with chemical \emph{pH} values.)
#' @param ppm Specifies ppm window (10 ppm will return the range of the
#' molecular mass + and - 10 ppm).
#' @param deltaMz Specifies additional m/z window to add to the range (deltaMz
#' = 0.02 will return the range of the molecular mass +- 0.02 (and additionally
#' +- the set ppm value).
#' @param retrieval A value that determines whether the files should be handled either as "standard",
#' if the compoundlist is complete, "tentative", if at least a formula is present or "unknown"
#' if the only know thing is the m/z
#' @param compact Only for \code{findLevel}, returns the "retrieval" parameter used for many functions 
#' within RMassBank if TRUE
#' @return \code{findMz} will return a \code{list(mzCenter=, mzMin=, mzMax=)}
#' with the molecular weight of the given ion, as calculated from the SMILES
#' code and Rcdk.
#' 
#' \code{findRt}, \code{findSmiles},\code{findCAS},\code{findName} will return
#' the corresponding entry from the compound list. \code{findFormula} returns
#' the molecular formula as determined from the SMILES code.
#' 
#' @author Michael Stravs
#' @seealso \code{\link{findMass}}, \code{\link{loadList}}, \code{\link{findMz.formula}}
#' @examples
#' 
#' \dontrun{%
#' 	findMz(123, "pH", 5)
#' 	findFormula(123)
#' }
#' 
#' @export
findMz <- function(cpdID, mode="pH", ppm=10, deltaMz=0, retrieval="standard",
                   unknownMass = getOption("RMassBank")$unknownMass)
{
  if(is.null(unknownMass))
    unknownMass = "charged"
  if(!(unknownMass %in% c("charged", "neutral")))
     stop("unknownMass definition must be either 'charged' or 'neutral', default is 'charged'")
     
    # In case of unknown: m/z is in table
    if(retrieval == "unknown"){
        mz <- .listEnvEnv$listEnv$compoundList[which(.listEnvEnv$listEnv$compoundList$ID == cpdID),.listEnvEnv$listEnv$mzCol]
        delta <- ppm(mz, ppm, l=TRUE)
        if(unknownMass == "neutral")
          dmass <- findMz.formula("", mode=mode)$mzCenter
        else
          dmass <- 0
        return(list(mzMin=delta[2] - deltaMz + dmass, mzMax=delta[1] + deltaMz + dmass, mzCenter=mz + dmass))
    } 
    
    # In case of tentative: calculate mass from formula
    if(retrieval == "tentative"){
        return(findMz.formula(findFormula(cpdID, "tentative"),mode=mode))
    }
    
    # All other cases: Use smiles to calculate mass
    if(retrieval == "standard"){
        s <- findSmiles(cpdID)
        #if(s=="-") s <- NA
        if(is.na(s)){
            stop("There was no matching SMILES-entry to the supplied cpdID(s) \n  Please check the cpdIDs and the compoundlist.")
            return(list(mzMin=NA,mzMax=NA,mzCenter=NA))
        }
        formula <- .get.mol2formula(getMolecule(s))
        return(findMz.formula(formula@string, mode, ppm, deltaMz))
    }
}

#findMz <- function(...)	findMz.rcdk(...)
#' @export
findRt <- function(cpdID) {
	if(is.null(.listEnvEnv$listEnv))
		stop("Compound list must be loaded first.")
	if(!exists("compoundList", where=.listEnvEnv$listEnv))
		stop("Compound list must be loaded first.")
	if(is.character(cpdID))
		cpdID <- as.numeric(cpdID)
	rt <- as.numeric(.listEnvEnv$listEnv$compoundList[which(.listEnvEnv$listEnv$compoundList$ID == cpdID),"RT"])
	if(!is.null(getOption("RMassBank")$rtShift))
		rt <- rt + getOption("RMassBank")$rtShift
	if(is.na(rt)) return(list(RT=NA))
	else return(list(RT=rt))
}

#' @export
findSmiles <- function(cpdID) {
	if(is.character(cpdID))
		cpdID <- as.numeric(cpdID)
	if(is.null(.listEnvEnv$listEnv))
		stop("Compound list must be loaded first.")
	if(!exists("compoundList", where=.listEnvEnv$listEnv))
		stop("Compound list must be loaded first.")
  if(.listEnvEnv$listEnv$compoundList[match(cpdID, .listEnvEnv$listEnv$compoundList$ID),"SMILES"] == "")
    return(NA)
	return(.listEnvEnv$listEnv$compoundList[match(cpdID, .listEnvEnv$listEnv$compoundList$ID),"SMILES"])
}

#' @export
findFormula <- function(cpdID, retrieval = "standard") {
    
    # In case of tentative: read formula from table
    if(retrieval=="tentative"){
        return(.listEnvEnv$listEnv$compoundList[which(.listEnvEnv$listEnv$compoundList$ID == cpdID),"Formula"])
    }
    
    # Otherwise: Convert smiles to formula
    if(retrieval=="standard"){
        smiles <- findSmiles(cpdID)
        mol <- getMolecule(smiles)
        f <- .get.mol2formula(mol)
        return(f@string)
    }
}

#' @export
findCAS <- function(cpdID) {
	if(is.character(cpdID))
		cpdID <- as.numeric(cpdID)
	if(is.null(.listEnvEnv$listEnv))
		stop("Compound list must be loaded first.")
	if(!exists("compoundList", where=.listEnvEnv$listEnv))
		stop("Compound list must be loaded first.")
	return(.listEnvEnv$listEnv$compoundList[which(.listEnvEnv$listEnv$compoundList$ID == cpdID),"CAS"])
}

#' @export
findName <- function(cpdID) {
	if(is.character(cpdID))
		cpdID <- as.numeric(cpdID)
	if(is.null(.listEnvEnv$listEnv))
		stop("Compound list must be loaded first.")
	if(!exists("compoundList", where=.listEnvEnv$listEnv))
		stop("Compound list must be loaded first.")
	return(.listEnvEnv$listEnv$compoundList[which(.listEnvEnv$listEnv$compoundList$ID == cpdID),"Name"])
}

#' @export
findLevel <- function(cpdID, compact=FALSE) {
	if(is.character(cpdID))
		cpdID <- as.numeric(cpdID)
	if(is.null(.listEnvEnv$listEnv))
		stop("Compound list must be loaded first.")
	if(!exists("compoundList", where=.listEnvEnv$listEnv))
		stop("Compound list must be loaded first.")
    if(compact){
        level <- .listEnvEnv$listEnv$compoundList[which(.listEnvEnv$listEnv$compoundList$ID == cpdID),"Level"]
        if(level %in% c("0","1","1a","1b","1c","2","2a","2b","3","3a")){
            return("standard")
        }
        if(level %in% c("3b","3c","3d","4")){
            return("tentative")
        }
        if(level == "5"){
            return("unknown")
        }
        
        return("tentative")
    }
	return(.listEnvEnv$listEnv$compoundList[which(.listEnvEnv$listEnv$compoundList$ID == cpdID),"Level"])
}

#' Calculate exact mass
#' 
#' Retrieves the exact mass of the uncharged molecule. It works directly from
#' the SMILES and therefore is used in the MassBank workflow
#' (\code{\link{mbWorkflow}}) - there, all properties are calculated from the
#' SMILES code retrieved from the database. (Alternatively, takes also the
#' compound ID as parameter and looks it up.) Calculation relies on Rcdk.
#' 
#' @param cpdID_or_smiles SMILES code or compound ID of the molecule. (Numerics
#' are treated as compound ID).
#' @param retrieval A value that determines whether the files should be handled either as "standard",
#' if the compoundlist is complete, "tentative", if at least a formula is present or "unknown"
#' if the only know thing is the m/z
#' @param mode \code{"pH", "pNa", "pM", "pNH4", "mH", "mM", "mFA"} for different ions 
#' 			([M+H]+, [M+Na]+, [M]+, [M+NH4]+, [M-H]-, [M]-, [M+FA]-).
#' Only needed for retrieval="unknown"
#' @return  Returns the exact mass of the uncharged molecule.
#' @author Michael Stravs
#' @seealso \code{\link{findMz}}
#' @examples
#' 
#' ##
#' findMass("OC[C@@H]1OC(O)[C@@H](O)[C@@@@H](O)[C@@@@H]1O")
#' 
#' @export
findMass <- function(cpdID_or_smiles, retrieval="standard", mode = "pH")
{
    # Must calculate mass manually if no formula is given
    if(retrieval == "unknown"){
        adductProperties <- getAdductProperties(mode, rcdk::get.formula(findFormula(cpdID_or_smiles)))
        allowed_additions <- adductProperties$addition
        mode.charge <- adductProperties$charge
        mass <- getMonoisotopicMass(allowed_additions)
        return(findMz(cpdID_or_smiles, mode=mode, retrieval=retrieval)$mzCenter - mass + mode.charge * .emass)
    }
    
    # Calculate mass from formula
    if(retrieval == "tentative"){
        return(get.formula(findFormula(cpdID_or_smiles, "tentative"))@mass)
    }
    
    # Calculate mass from SMILES
    if(retrieval == "standard"){
        if(!is.numeric(cpdID_or_smiles))
            s <- cpdID_or_smiles
        else
            s <- findSmiles(cpdID_or_smiles)
        mol <- getMolecule(s)
        return(get.exact.mass(mol))
    }
}

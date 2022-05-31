# Script for writing MassBank files

#testtest change
#' Load MassBank compound information lists
#' 
#' Loads MassBank compound information lists (i.e. the lists which were created
#' in the first two steps of the MassBank \code{\link{mbWorkflow}} and
#' subsequently edited by hand.).
#' 
#' \code{resetInfolists} clears the information lists, i.e. it creates a new
#' empty list in \code{mbdata_archive}. \code{loadInfolist} loads a single CSV
#' file, whereas \code{loadInfolists} loads a whole directory.
#' 
#' @aliases loadInfolists loadInfolist resetInfolists
#' @usage loadInfolists(mb, path)
#' 
#'  loadInfolist(mb, fileName)
#' 
#'  resetInfolists(mb)
#' @param path Directory in which the namelists reside. All CSV files in this
#' directory will be loaded.
#' @param fileName A single namelist to be loaded.
#' @param mb The \code{mbWorkspace} to load/reset the lists in.
#' @return The new workspace with loaded/reset lists.
#' @author Michael Stravs
#' @examples
#' 
#' #
#' \dontrun{mb <- resetInfolists(mb)
#' 	mb <- loadInfolist(mb, "my_csv_infolist.csv")}
#' 
#' @export
loadInfolists <- function(mb, path)
{
  archivefiles <- list.files(path, ".csv", full.names=TRUE)
  for(afile in archivefiles)
    mb <- loadInfolist(mb, afile)
  return(mb)
}

# Load an "infolist". This loads a CSV file which should contain the entries
# edited and controlled by hand. All compound infos from fileName are added into the
# global mbdata_archive. Entries with a cpdID which was already present, are substituted
# by new entries from the fileName file.
#' @export
loadInfolist <- function(mb, fileName)
{
  # Prime a new infolist if it doesn't exist
  if(ncol(mb@mbdata_archive) == 0)
    mb <- resetInfolists(mb)
  mbdata_new <- read.csv(fileName, sep=",", stringsAsFactors=FALSE)
  # Legacy check for loading the Uchem format files.
  # Even if dbname_* are not used downstream of here, it's still good to keep them
  # for debugging reasons.
  n <- colnames(mbdata_new)
  cols <- c("id","dbcas","dataused")
  
  # Check if comma-separated or semicolon-separated
  d <- setdiff(cols, n)
  if(length(d)>0){
		mbdata_new <- read.csv2(fileName, stringsAsFactors=FALSE)
		n <- colnames(mbdata_new)
		d2 <- setdiff(cols, n)
		if(length(d2) > 0){
			stop("Some columns are missing in the infolist.")
		}
	}
  if("dbname_d" %in% colnames(mbdata_new))
  {
    colnames(mbdata_new)[[which(colnames(mbdata_new)=="dbname_d")]] <- "dbname"
    # dbname_e will be dropped because of the select= in the subset below.
  }
  if("COMMENT.EAWAG_UCHEM_ID" %in% colnames(mbdata_new))
    colnames(mbdata_new)[[which(colnames(mbdata_new)== "COMMENT.EAWAG_UCHEM_ID")]] <-
      "COMMENT.ID"
  
  # Clear from padding spaces and NAs
  mbdata_new <- as.data.frame(x = t(apply(mbdata_new, 1, function(r) 
    {
    # Substitute empty spaces by real NA values
    r[which(r == "")] <- NA
    # Trim spaces (in all non-NA fields)
    r[which(!is.na(r))] <- sub(pattern = "^ *([^ ]+) *$", replacement = "\\1", x = r[which(!is.na(r))])
    return(r)
  })), stringsAsFactors = FALSE)
  # use only the columns present in mbdata_archive, no other columns added in excel
  colNames <- colnames(mb@mbdata_archive)
  commentColNames <- colnames(mbdata_new)[grepl(x = colnames(mbdata_new), pattern = "^COMMENT\\.(?!CONFIDENCE)(?!ID)", perl = TRUE)]
  colNames <- c(colNames, commentColNames)

  ## The read infolists might not have all required / expected columns
  missingColNames <- colNames[! colNames %in% colnames(mbdata_new)]
  if (length(missingColNames >0)) {
    missingCols <- matrix(NA, ncol=length(missingColNames))
    colnames(missingCols) <- missingColNames
    mbdata_new <- cbind(mbdata_new, missingCols)
  }
    
  mbdata_new <- mbdata_new[, colNames]
  # substitute the old entires with the ones from our files
  # then find the new (previously inexistent) entries, and rbind them to the table
  new_entries <- setdiff(mbdata_new$id, mb@mbdata_archive$id)
  old_entries <- intersect(mbdata_new$id, mb@mbdata_archive$id)
  
  for(colname in colnames(mb@mbdata_archive))
    mb@mbdata_archive[, colname] <- as.character(mb@mbdata_archive[, colname])
  
  for(entry in old_entries)
    mb@mbdata_archive[mb@mbdata_archive$id == entry,] <- mbdata_new[mbdata_new$id == entry,]
  mb@mbdata_archive <- rbind(mb@mbdata_archive, mbdata_new[mbdata_new$id==new_entries,])
  
  for(colname in colnames(mb@mbdata_archive))
    mb@mbdata_archive[, colname] <- as.factor(mb@mbdata_archive[, colname])
  
  return(mb)
}


# Resets the mbdata_archive to an empty version.
#' @export
resetInfolists <- function(mb) 
{    
	mb@mbdata_archive <-
			structure(list(X = integer(0), id = integer(0), dbcas = character(0), 
							dbname = character(0), dataused = character(0), COMMENT.CONFIDENCE = character(0), 
							COMMENT.ID = integer(0), CH.NAME1 = character(0), 
							CH.NAME2 = character(0), CH.NAME3 = character(0), CH.NAME4 = character(0), CH.NAME5 = character(0), CH.COMPOUND_CLASS = character(0), 
							CH.FORMULA = character(0), CH.EXACT_MASS = numeric(0), CH.SMILES = character(0), 
							CH.IUPAC = character(0), CH.LINK.CAS = character(0), CH.LINK.CHEBI = integer(0), 
							CH.LINK.HMDB = character(0), CH.LINK.KEGG = character(0), CH.LINK.LIPIDMAPS = character(0), 
							CH.LINK.PUBCHEM = character(0), CH.LINK.INCHIKEY = character(0), 
							CH.LINK.CHEMSPIDER = integer(0), CH.LINK.COMPTOX = character(0)), .Names = c("X", "id", "dbcas", 
							"dbname", "dataused", "COMMENT.CONFIDENCE", "COMMENT.ID", 
              "CH.NAME1", "CH.NAME2", "CH.NAME3", "CH.NAME4", "CH.NAME5", "CH.COMPOUND_CLASS", "CH.FORMULA", 
							"CH.EXACT_MASS", "CH.SMILES", "CH.IUPAC", "CH.LINK.CAS", "CH.LINK.CHEBI", 
							"CH.LINK.HMDB", "CH.LINK.KEGG", "CH.LINK.LIPIDMAPS", "CH.LINK.PUBCHEM",
							"CH.LINK.INCHIKEY", "CH.LINK.CHEMSPIDER", "CH.LINK.COMPTOX"), row.names = integer(0), class = "data.frame")
	if(getOption("RMassBank")$include_sp_tags)
	{
	  mb@mbdata_archive["SP.SAMPLE"] <- character(0)
	}
	return(mb)
	
}

# The workflow function, i.e. (almost) the only thing you actually need to call.
# See below for explanation of steps.
#' MassBank record creation workflow
#' 
#' Uses data generated by \code{\link{msmsWorkflow}} to create MassBank records.
#' 
#' See the vignette \code{vignette("RMassBank")} for detailed informations about the usage.
#' 
#' Steps:
#' 
#' Step 1: Find which compounds don't have annotation information yet. For these
#' 		 compounds, pull information from several databases (using gatherData).
#' 
#' Step 2: If new compounds were found, then export the infolist.csv and stop the workflow.
#' 		Otherwise, continue.
#' 
#' Step 3: Take the archive data (in table format) and reformat it to MassBank tree format.
#' 
#' Step 4: Compile the spectra. Using the skeletons from the archive data, create
#'   MassBank records per compound and fill them with peak data for each spectrum.
#'   Also, assign accession numbers based on scan mode and relative scan no.
#' 
#' Step 5: Convert the internal tree-like representation of the MassBank data into
#'  flat-text string arrays (basically, into text-file style, but still in memory)
#' 
#' Step 6: For all OK records, generate a corresponding molfile with the structure
#'   of the compound, based on the SMILES entry from the MassBank record. (This molfile
#'   is still in memory only, not yet a physical file)
#' 
#' Step 7: If necessary, generate the appropriate subdirectories, and actually write
#'   the files to disk.
#' 
#' Step 8: Create the list.tsv in the molfiles folder, which is required by MassBank
#'   to attribute substances to their corresponding structure molfiles. 
#' 
#' @param steps Which steps in the workflow to perform.
#' @param infolist_path A path where to store newly downloaded compound informations,
#' 			which should then be manually inspected.
#' @param mb The \code{mbWorkspace} to work in.
#' @param gatherData A variable denoting whether to retrieve information using several online databases \code{gatherData= "online"}
#' or to use the local babel installation \code{gatherData= "babel"}. Note that babel is used either way, if a directory is given 
#' in the settings. This setting will be ignored if retrieval is set to "standard"
#' @param filter If \code{TRUE}, the peaks will be filtered according to the standard processing workflow in RMassBank - 
#' only the best formula for a peak is retained, and only peaks passing multiplicity filtering are retained. If FALSE, it is assumed
#' that the user has already done filtering, and all peaks in the spectrum should be printed in the record (with or without formula.)
#' @return The processed \code{mbWorkspace}.
#' @seealso \code{\link{mbWorkspace-class}}
#' @author Michael A. Stravs, Eawag <michael.stravs@@eawag.ch>
#' @examples \dontrun{
#' 		mb <- newMbWorkspace(w) # w being a msmsWorkspace
#' 		mb <- loadInfolists(mb, "D:/myInfolistPath")
#' 		mb <- mbWorkflow(mb, steps=c(1:3), "newinfos.csv")
#' 		
#' }
#' @export
mbWorkflow <- function(mb, steps=c(1,2,3,4,5,6,7,8), infolist_path="./infolist.csv", gatherData = "online", filter = TRUE)
{
    # Step 1: Find which compounds don't have annotation information yet. For these
    # compounds, pull information from CTS (using gatherData).
    if(1 %in% steps)
    {
        mbdata_ids <- lapply(selectSpectra(mb@spectra, "found", "object"), function(spec) spec@id)
                rmb_log_info("mbWorkflow: Step 1. Gather info from several databases")
      # Which IDs are not in mbdata_archive yet?
      new_ids <- setdiff(as.numeric(unlist(mbdata_ids)), mb@mbdata_archive$id)
      mb@mbdata <- lapply(new_ids, function(id) 
      {
            if(findLevel(id,TRUE) == "standard"){
            if(gatherData == "online"){
                    
                d <- gatherData(id)
            } 
            if(gatherData == "babel"){
                    # message("mbWorkflow: Step 1. Gather info using babel")
                d <- gatherDataBabel(id)
            }
        } else{
                # message("mbWorkflow: Step 1. Gather no info - Unknown structure")
                d <- gatherDataUnknown(id, mb@spectra[[1]]@mode, retrieval=findLevel(id,TRUE))
        }
		rmb_log_info(paste(id, ": ", d$dataused, sep=''))
        return(d)
      })
  }
  # Step 2: If new compounds were found, then export the infolist.csv and stop the workflow.
  # Otherwise, continue!
  if(2 %in% steps)
  {
	rmb_log_info("mbWorkflow: Step 2. Export infolist (if required)")
    if(length(mb@mbdata)>0)
    {
      mbdata_mat <- flatten(mb@mbdata)
      write.csv(as.data.frame(mbdata_mat),infolist_path, na="")
            rmb_log_info(paste("The file", infolist_path, "was generated with new compound information. Please check and edit the table, and add it to your infolist folder."))
      return(mb)
    }
    else
      rmb_log_info("No new data added.")
  }
  # Step 3: Take the archive data (in table format) and reformat it to MassBank tree format.
  if(3 %in% steps)
  {
	rmb_log_info("mbWorkflow: Step 3. Data reformatting")
    mb@mbdata_relisted <- apply(mb@mbdata_archive, 1, readMbdata)
  }
  # Step 4: Compile the spectra! Using the skeletons from the archive data, create
  # MassBank records per compound and fill them with peak data for each spectrum.
  # Also, assign accession numbers based on scan mode and relative scan no.
  if(4 %in% steps)
  {
	  rmb_log_info("mbWorkflow: Step 4. Spectra compilation")
	  mb@compiled <- lapply(
			  selectSpectra(mb@spectra, "found", "object"),
			  function(r) {
				  rmb_log_info(paste("Compiling: ", r@name, sep=""))
				  mbdata <- mb@mbdata_relisted[[which(mb@mbdata_archive$id == as.numeric(r@id))]]
				  if(filter)
            res <- buildRecord(r, mbdata=mbdata, additionalPeaks=mb@additionalPeaks, filter = filterOK & best)
				  else
				    res <- buildRecord(r, mbdata=mbdata, additionalPeaks=mb@additionalPeaks)
          return(res)
			  })
	  # check which compounds have useful spectra
	  mb@ok <- which(!is.na(mb@compiled) & !(lapply(mb@compiled, length)==0))
        #mb@ok <- which(!is.na(mb@compiled) & !(lapply(mb@compiled, length)==0))
	  mb@problems <- which(is.na(mb@compiled))
	  mb@compiled_ok <- mb@compiled[mb@ok]
    mb@compiled_notOk <- mb@compiled[!mb@ok]
  }
  # Step 5: Convert the internal tree-like representation of the MassBank data into
  # flat-text string arrays (basically, into text-file style, but still in memory)
  if(5 %in% steps)
  {
	rmb_log_info("mbWorkflow: [Legacy Step 5. Flattening records] ignored")
    #mb@mbfiles <- lapply(mb@compiled_ok, function(cpd) toMassbank(cpd, mb@additionalPeaks))
    #mb@mbfiles_notOk <- lapply(mb@compiled_notOk, function(c) lapply(c, toMassbank))
  }
  # Step 6: For all OK records, generate a corresponding molfile with the structure
  # of the compound, based on the SMILES entry from the MassBank record. (This molfile
  # is still in memory only, not yet a physical file)
  if(6 %in% steps)
  {
    if(RMassBank.env$export.molfiles){
      rmb_log_info("mbWorkflow: Step 6. Generate molfiles")
      mb@molfile <- lapply(mb@compiled_ok, function(c) createMolfile(as.numeric(c@id)))
    } else
      warning("RMassBank is configured not to export molfiles (RMassBank.env$export.molfiles). Step 6 is therefore ignored.")
    }
  # Step 7: If necessary, generate the appropriate subdirectories, and actually write
  # the files to disk.
  if(7 %in% steps)
  {
	rmb_log_info("mbWorkflow: Step 7. Generate subdirs and export")
        
        ## create folder
        filePath_recData_valid   <- file.path(getOption("RMassBank")$annotations$entry_prefix, "recdata")
        filePath_recData_invalid <- file.path(getOption("RMassBank")$annotations$entry_prefix, "recdata_invalid")
        filePath_molData         <- file.path(getOption("RMassBank")$annotations$entry_prefix, "moldata")
        
        if(!file.exists(filePath_recData_valid)) if(!dir.create(filePath_recData_valid,recursive=TRUE))  stop(paste("Could not create folder", filePath_recData_valid))
        if(RMassBank.env$export.molfiles)
          if(!file.exists(filePath_molData)) if(!dir.create(filePath_molData,recursive=TRUE))  stop(paste("Could not create folder", filePath_molData))
        if(RMassBank.env$export.invalid & length(mb@mbfiles_notOk) > 0)
          if(!file.exists(filePath_recData_invalid)) if(!dir.create(filePath_recData_invalid,recursive=TRUE))  stop(paste("Could not create folder", filePath_recData_invalid))
        
        if(length(mb@molfile) == 0)
            mb@molfile <- as.list(rep(x = NA, times = length(mb@compiled_ok)))
        
        ## export valid spectra
        for(cnt in seq_along(mb@compiled_ok)){
            exportMassbank_recdata(
                mb@compiled_ok[[cnt]], 
                recDataFolder = filePath_recData_valid
            )
            if(RMassBank.env$export.molfiles)
              exportMassbank_moldata(
                mb@compiled_ok[[cnt]], 
                molfile = mb@molfile[[cnt]], 
                molDataFolder = filePath_molData
              )
        }
        
        ## export invalid spectra
            for(cnt in seq_along(mb@compiled_notOk))
                exportMassbank_recdata(
                    compiled = mb@mbfiles_notOk[[cnt]], 
                    recDataFolder = filePath_recData_invalid
                )
  }
  # Step 8: Create the list.tsv in the molfiles folder, which is required by MassBank
  # to attribute substances to their corresponding structure molfiles.
  if(8 %in% steps)
  {
        if(RMassBank.env$export.molfiles){
          rmb_log_info("mbWorkflow: Step 8. Create list.tsv")
          makeMollist(compiled = mb@compiled_ok)
        } else
            warning("RMassBank is configured not to export molfiles (RMassBank.env$export.molfiles). Step 8 is therefore ignored.")
  }
  return(mb)
}


# Calls openbabel and converts the SMILES code string (or retrieves the SMILES code from
# the ID, and then calls openbabel) to create a molfile in text format.
# If fileName is given, the file is directly stored. Otherwise, it is returned as a 
# character array.
#' Create MOL file for a chemical structure
#' 
#' Creates a MOL file (in memory or on disk) for a compound specified by the
#' compound ID or by a SMILES code.
#' 
#' The function invokes OpenBabel (and therefore needs a correctly set
#' OpenBabel path in the RMassBank settings), using the SMILES code retrieved
#' with \code{findSmiles} or using the SMILES code directly. The current
#' implementation of the workflow uses the latter version, reading the SMILES
#' code directly from the MassBank record itself.
#' 
#' @usage createMolfile(id_or_smiles, fileName = FALSE)
#' @param id_or_smiles The compound ID or a SMILES code.
#' @param fileName If the filename is set, the file is written directly to disk
#' using the specified filename. Otherwise, it is returned as a text array.
#' @return A character array containing the MOL/SDF format file, ready to be
#' written to disk.
#' @author Michael Stravs
#' @seealso \code{\link{findSmiles}}
#' @references OpenBabel: \url{http://openbabel.org}
#' @examples
#' 
#' # Benzene:
#' \dontrun{
#' createMolfile("C1=CC=CC=C1")
#' }
#' 
#' @export
createMolfile <- function(id_or_smiles, fileName = FALSE)
{
	.checkMbSettings()
	babeldir <- getOption("RMassBank")$babeldir
    
	if(!is.numeric(id_or_smiles)){
		smiles <- id_or_smiles
    } else{
        if(findLevel(id_or_smiles,TRUE) != "standard"){
            return(c(" ","$$$$"))
        }
		smiles <- findSmiles(id_or_smiles)
    }
    # if no babeldir was set, get the result from cactus.
	if(is.na(babeldir))
	{
		res <- getCactus(smiles, "sdf")
		
		if(any(is.na(res))){
			res <- getPcSDF(smiles)
		}
		if(any(is.na(res))){
			stop("Pubchem and Cactus both seem to be down.")
		}
		if(is.character(fileName))
			writeLines(res, fileName)
	}
	# otherwise use the better-tested OpenBabel toolkit.
	else
	{
		if(!is.character(fileName))
			cmd <- paste(babeldir, "obabel -ismi -osdf -d -b --gen2D", sep='')
		else
			cmd <- paste(babeldir, "obabel -ismi -osdf ", fileName , " -d -b --gen2D", sep='')
		res <- system(cmd, intern=TRUE, input=smiles, ignore.stderr=TRUE)
		# If we wrote to a file, read it back as return value.
		if(is.character(fileName))
			res <- readLines(fileName)
	} 
  #return(c(" ","$$$$"))
	return(res)
}



# Retrieve annotation data for a compound, from the internet service Pubchem
#' Retrieve supplemental annotation data from Pubchem
#' 
#' Retrieves annotation data for a compound from the internet service Pubchem 
#' based on the inchikey generated by babel or Cactus
#' 
#' The data retrieved is the Pubchem CID, a synonym from the Pubchem database,
#' the IUPAC name (using the preferred if available) and a Chebi link
#' 
#' @usage gatherPubChem(key)
#' @param key An Inchi-Key
#' @return Returns a list with 4 slots:
#' \code{PcID} The Pubchem CID
#' \code{Synonym} An arbitrary synonym for the compound name
#' \code{IUPAC} A IUPAC-name (preferred if available)
#' \code{Chebi} The identification number of the chebi database
#' @author Erik Mueller
#' @seealso \code{\link{mbWorkflow}}
#' @references Pubchem REST:
#' \url{https://pubchem.ncbi.nlm.nih.gov/pug_rest/PUG_REST.html}
#' Chebi:
#' \url{http://www.ebi.ac.uk/chebi}
#' @examples
#' 
#' # Gather data for compound ID 131
#' \dontrun{gatherPubChem("QEIXBXXKTUNWDK-UHFFFAOYSA-N")}
#' 
#' @export
gatherPubChem <- function(key){
	
	PubChemData <- list()
	
	##Trycatches are there because pubchem has connection issues 1 in 50 times.
	##Write NA into the respective fields if something goes wrong with the conenction or the data.
	
	##Retrieve Pubchem CID
	tryCatch(
		PubChemData$PcID <- getPcId(key),
		error=function(e){
		PubChemData$PcID <<- NA
	})
	
	##Retrieve a synonym to the name
	tryCatch(
		PubChemData$Synonym <- getPcSynonym(key),
		error=function(e){
		PubChemData$Synonym <<- NA
	})
	
	##Retrieve the IUPAC-name
	tryCatch(
		PubChemData$IUPAC <- getPcIUPAC(key),
		error=function(e){
		PubChemData$IUPAC <<- NA
	})
	
	##Retrieve the Chebi-ID
	tryCatch(
		PubChemData$Chebi <- getPcCHEBI(key),
		error=function(e){
		PubChemData$Chebi <<- NA
	})
	
	return(PubChemData)
}

# Retrieve annotation data for a compound, from the internet services Cactvs, Pubchem, Chemspider and CTS.
#' Retrieve annotation data
#' 
#' Retrieves annotation data for a compound from the internet services CTS, Pubchem, Chemspider and
#' Cactvs, based on the SMILES code and name of the compounds stored in the
#' compound list.
#' 
#' Composes the "upper part" of a MassBank record filled with chemical data
#' about the compound: name, exact mass, structure, CAS no., links to PubChem,
#' KEGG, ChemSpider.  The instrument type is also written into this block (even
#' if not strictly part of the chemical information). Additionally, index
#' fields are added at the start of the record, which will be removed later:
#' \code{id, dbcas, dbname} from the compound list, \code{dataused} to indicate
#' the used identifier for CTS search (\code{smiles} or \code{dbname}).
#' 
#' Additionally, the fields \code{ACCESSION} and \code{RECORD_TITLE} are
#' inserted empty and will be filled later on.
#' 
#' @usage gatherData(id)
#' @aliases gatherData
#' @param id The compound ID.
#' @return Returns a list of type \code{list(id= \var{compoundID}, ...,
#' 'ACCESSION' = '', 'RECORD_TITLE' = '', )} etc. %% ...
#' @author Michael Stravs
#' @seealso \code{\link{mbWorkflow}}
#' @references Chemical Translation Service:
#' \url{http://uranus.fiehnlab.ucdavis.edu:8080/cts/homePage} 
#' cactus Chemical Identifier Resolver: 
#' \url{http://cactus.nci.nih.gov/chemical/structure}
#' MassBank record format:
#' \url{http://www.massbank.jp/manuals/MassBankRecord_en.pdf}
#' Pubchem REST:
#' \url{https://pubchem.ncbi.nlm.nih.gov/pug_rest/PUG_REST.html}
#' Chemspider InChI conversion:
#' \url{https://www.chemspider.com/InChI.asmx}
#' @examples
#' 
#' # Gather data for compound ID 131
#' \dontrun{gatherData(131)}
#' 
#' @export
gatherData <- function(id)
{ 
	##Preamble: Is a babeldir supplied?
	##If yes, use it
	
	.checkMbSettings()
	usebabel=TRUE
	babeldir <- getOption("RMassBank")$babeldir
	
	if(is.na(babeldir)){
		usebabel=FALSE
	}
	
	
	##Get all useful information from the local "database" (from the CSV sheet)
	
	smiles <- findSmiles(id)
	mass <- findMass(smiles)
	dbcas <- findCAS(id)
	dbname <- findName(id)
	if(is.na(dbname)) dbname <- ""
	if(is.na(dbcas)) dbcas <- ""
	iupacName <- dbname
	synonym <- dbname
	formula <- findFormula(id)
	
	##Convert SMILES to InChI key via Cactvs or babel. CTS doesn't "interpret" the SMILES per se,
	##it just matches identical known SMILES, so we need to convert to a "searchable" and
	##standardized format beforehand. Other databases are able to interpret the smiles.
	
	if(usebabel){
		cmdinchikey <- paste0(babeldir, 'obabel -:"',smiles,'" ', '-oinchikey')
		inchikey_split <- system(cmdinchikey, intern=TRUE, input=smiles, ignore.stderr=TRUE)
	} else{
		inchikey <- getCactus(smiles, 'stdinchikey')
		if(!is.na(inchikey)){
			##Split the "InChiKey=" part off the key
			inchikey_split <- strsplit(inchikey, "=", fixed=TRUE)[[1]][[2]]
		} else{
		    inchikey_split <- getPcInchiKey(smiles)
		}
	}
	
	##Use Pubchem to retrieve information
	PcInfo <- gatherPubChem(inchikey_split)
	
	if(!is.null(PcInfo$Synonym) & !is.na(PcInfo$Synonym)){
		synonym <- PcInfo$Synonym
	}
	
	if(!is.null(PcInfo$IUPAC) & !is.na(PcInfo$IUPAC)){
		iupacName <- PcInfo$IUPAC
	}
	
	##Get Chemspider-ID
	csid <- getCSID(inchikey_split)
	
	if(is.na(csid)){
		##Get ChemSpider ID from Cactus if the Chemspider page is down
		csid <- getCactus(inchikey_split, 'chemspider_id')
	}
	
	## ##Get CompTox
	## comptox <- getCompTox(inchikey_split)
	
	## if(is.null(comptox)){
 	  comptox <- NA
	## }
	
	##Use CTS to retrieve information
	CTSinfo <- getCtsRecord(inchikey_split)
		
	if((CTSinfo[1] == "Sorry, we couldn't find any matching results") || is.null(CTSinfo[1]))
	{
		CTSinfo <- NA
	}
	
	##List the names
	if(iupacName == ""){
		warning(paste0("Compound ID ",id,": no IUPAC name could be identified."))
	}

	if(toupper(dbname) == toupper(synonym)){
		synonym <- dbname
	}
	
	if(toupper(dbname) == toupper(iupacName)){
		iupacName <- dbname
	}
	
	if(toupper(synonym) == toupper(iupacName)){
		synonym <- iupacName
	}
	
	names <- as.list(unique(c(dbname, synonym, iupacName)))
	
	##If no name is found, it must be supplied in one way or another
	if(all(sapply(names, function(x) x == ""))){
		stop("RMassBank wasn't able to extract a usable name for this compound from any database. Please supply a name manually.")
	}
	
	# Start to fill the MassBank record.
	# The top 4 entries will not go into the final record; they are used to identify
	# the record and also to facilitate manual editing of the exported record table.
	mbdata <- list()
	mbdata[['id']] <- id
	mbdata[['dbcas']] <- dbcas
	mbdata[['dbname']] <- dbname
	mbdata[['dataused']] <- "smiles"
	mbdata[['ACCESSION']] <- ""
	mbdata[['RECORD_TITLE']] <- ""
	mbdata[['DATE']] <- format(Sys.Date(), "%Y.%m.%d")
	mbdata[['AUTHORS']] <- getOption("RMassBank")$annotations$authors
	mbdata[['LICENSE']] <- getOption("RMassBank")$annotations$license
	mbdata[['COPYRIGHT']] <- getOption("RMassBank")$annotations$copyright
	# Confidence annotation and internal ID annotation.
	# The ID of the compound will be written like:
	# COMMENT: EAWAG_UCHEM_ID 1234
	# if annotations$internal_id_fieldname is set to "EAWAG_UCHEM_ID"
	mbdata[["COMMENT"]] <- list()
  if(findLevel(id) == "0"){
	mbdata[["COMMENT"]][["CONFIDENCE"]] <- getOption("RMassBank")$annotations$confidence_comment
	} else{
        level <- findLevel(id)
        if(level %in% c("1","1a")){
             mbdata[["COMMENT"]][["CONFIDENCE"]] <- "Reference Standard (Level 1)"
        }
        if(level == c("2")){
             mbdata[["COMMENT"]][["CONFIDENCE"]] <- "Probable structure, tentative identification (Level 2)"
        }
        if(level == c("2a")){
             mbdata[["COMMENT"]][["CONFIDENCE"]] <- "Probable structure via library match, tentative identification (Level 2a)"
        }
        if(level == c("2b")){
             mbdata[["COMMENT"]][["CONFIDENCE"]] <- "Probable structure via diagnostic evidence, tentative identification (Level 2b)"
        }
        if(level == c("3")){
             mbdata[["COMMENT"]][["CONFIDENCE"]] <- "Tentative identification only (Level 3)"
        }
        if(level == c("3a")){
             mbdata[["COMMENT"]][["CONFIDENCE"]] <- "Tentative identification: most likely structure (Level 3)"
        }
        if(level == c("3b")){
             mbdata[["COMMENT"]][["CONFIDENCE"]] <- "Tentative identification: isomers possible (Level 3)"
        }
        if(level == c("3c")){
             mbdata[["COMMENT"]][["CONFIDENCE"]] <- "Tentative identification: substance class known (Level 3)"
        }
        if(level == c("3d")){
             mbdata[["COMMENT"]][["CONFIDENCE"]] <- "Tentative identification: best match only (Level 3)"
        }
        if(level == c("4")){
             mbdata[["COMMENT"]][["CONFIDENCE"]] <- "Tentative identification: molecular formula only (Level 4)"
        }
        if(level == c("5")){
             mbdata[["COMMENT"]][["CONFIDENCE"]] <- "Tentative identification: structure and formula unknown (Level 5)"
        }
	}
	
	mbdata[["COMMENT"]][["ID"]] = id
  
  ## add generic COMMENT information
  rowIdx <- which(.listEnvEnv$listEnv$compoundList$ID == id)
  properties      <- colnames(.listEnvEnv$listEnv$compoundList)
  properties2     <- gsub(x = properties, pattern = "^COMMENT ", replacement = "")
  theseProperties <- grepl(x = properties, pattern = "^COMMENT ")
  theseProperties <- theseProperties & (!(unlist(.listEnvEnv$listEnv$compoundList[rowIdx, ]) == "NA" | is.na(unlist(.listEnvEnv$listEnv$compoundList[rowIdx, ]))))
  mbdata[["COMMENT"]][properties2[theseProperties]] <- unlist(.listEnvEnv$listEnv$compoundList[rowIdx, theseProperties])
  
	# here compound info starts
	mbdata[['CH$NAME']] <- names
	# Currently we use a fixed value for Compound Class, since there is no useful
	# convention of what should go there and what shouldn't, and the field is not used
	# in search queries.
	mbdata[['CH$COMPOUND_CLASS']] <- getOption("RMassBank")$annotations$compound_class
	mbdata[['CH$FORMULA']] <- formula
	mbdata[['CH$EXACT_MASS']] <- mass
	mbdata[['CH$SMILES']] <- smiles
	
	if(usebabel){
		cmdinchi <- paste0(babeldir, 'obabel -:"',smiles,'" ', '-oinchi')
		mbdata[['CH$IUPAC']] <- system(cmdinchi, intern=TRUE, input=smiles, ignore.stderr=TRUE)
	} else{
		mbdata[['CH$IUPAC']] <- getCactus(smiles, "stdinchi")
	}
	

	
	# Add all CH$LINK fields present in the compound datasets
	link <- list()
	# CAS
	if(!is.na(CTSinfo[1])){
		if("CAS" %in% CTS.externalIdTypes(CTSinfo))
		{
			# Prefer database CAS if it is also listed in the CTS results.
			# otherwise take the shortest one.
			cas <- CTS.externalIdSubset(CTSinfo,"CAS")
			if(dbcas %in% cas)
				link[["CAS"]] <- dbcas
			else
				link[["CAS"]] <- cas[[which.min(nchar(cas))]]
		} else{
			if(dbcas != ""){
				link[["CAS"]] <- dbcas
			}
		}
	} else{
		if(dbcas != ""){
			link[["CAS"]] <- dbcas
		}
	}
	
	
	# CHEBI
	if(is.na(PcInfo$Chebi[1])){
		if(!is.na(CTSinfo[1])){
			if("ChEBI" %in% CTS.externalIdTypes(CTSinfo))
			{
				# Cut off front "CHEBI:" if present
				chebi <- CTS.externalIdSubset(CTSinfo,"ChEBI")
				chebi <- chebi[[which.min(nchar(chebi))]]
				chebi <- strsplit(chebi,":")[[1]]
				link[["CHEBI"]] <- chebi[[length(chebi)]]
			}
		}
	} else{
		chebi <- PcInfo$Chebi
		chebi <- chebi[[which.min(nchar(chebi))]]
		chebi <- strsplit(chebi,":")[[1]]
		link[["CHEBI"]] <- chebi[[length(chebi)]]
	}
	# HMDB
	if(!is.na(CTSinfo[1])){
		if("Human Metabolome Database" %in% CTS.externalIdTypes(CTSinfo))
			link[["HMDB"]] <- CTS.externalIdSubset(CTSinfo,"HMDB")[[1]]
		# KEGG
		if("KEGG" %in% CTS.externalIdTypes(CTSinfo))
			link[["KEGG"]] <- CTS.externalIdSubset(CTSinfo,"KEGG")[[1]]
		# LipidMAPS
		if("LipidMAPS" %in% CTS.externalIdTypes(CTSinfo))
			link[["LIPIDMAPS"]] <- CTS.externalIdSubset(CTSinfo,"LipidMAPS")[[1]]
	}
	# PubChem CID
	if(is.na(PcInfo$PcID[1])){
		if(!is.na(CTSinfo[1])){
			if("PubChem CID" %in% CTS.externalIdTypes(CTSinfo))
			{
				pc <- CTS.externalIdSubset(CTSinfo,"PubChem CID")
				link[["PUBCHEM"]] <- paste0(min(pc))
			}
		}
	} else{
		link[["PUBCHEM"]] <- PcInfo$PcID[1]
	}
	
	
	if(!is.null(link[["PUBCHEM"]])){
		if(substr(link[["PUBCHEM"]],1,4) != "CID:"){
			link[["PUBCHEM"]] <- paste0("CID:", link[["PUBCHEM"]])
		}
	}
	
	link[["INCHIKEY"]] <- inchikey_split
	link[["COMPTOX"]] <- comptox
	if(length(csid)>0) if(any(!is.na(csid))) link[["CHEMSPIDER"]] <- min(as.numeric(as.character(csid[!is.na(csid)])))
	mbdata[['CH$LINK']] <- link
		
	return(mbdata)  
}

# Retrieve annotation data for a compound, using only babel
#' Retrieve annotation data
#' 
#' Retrieves annotation data for a compound by using babel,
#' based on the SMILES code and name of the compounds stored in the
#' compound list.
#' 
#' Composes the "upper part" of a MassBank record filled with chemical data
#' about the compound: name, exact mass, structure, CAS no..  
#' The instrument type is also written into this block (even
#' if not strictly part of the chemical information). Additionally, index
#' fields are added at the start of the record, which will be removed later:
#' \code{id, dbcas, dbname} from the compound list.
#' 
#' Additionally, the fields \code{ACCESSION} and \code{RECORD_TITLE} are
#' inserted empty and will be filled later on.
#' 
#' This function is an alternative to gatherData, in case CTS is down or if information
#' on one or more of the compounds in the compound list are sparse
#'
#' @usage gatherDataBabel(id)
#' @param id The compound ID.
#' @return Returns a list of type \code{list(id= \var{compoundID}, ...,
#' 'ACCESSION' = '', 'RECORD_TITLE' = '', )} etc. %% ...
#' @author Michael Stravs, Erik Mueller
#' @seealso \code{\link{mbWorkflow}}
#' @references MassBank record format:
#' \url{http://www.massbank.jp/manuals/MassBankRecord_en.pdf}
#' @examples
#' 
#' # Gather data for compound ID 131
#' \dontrun{gatherDataBabel(131)}
#' 
#' @export
gatherDataBabel <- function(id){
		.checkMbSettings()
		babeldir <- getOption("RMassBank")$babeldir
		smiles <- findSmiles(id)
			
		
		# if no babeldir was set, throw an error that says that either CTS or babel have to be used
		if(is.na(babeldir))
		{
			stop("No babeldir supplied; It is currently not possible to convert the information without either babel or CTS")
		} else {
			###Babel conversion
			cmdinchikey <- paste0(babeldir, 'obabel -:"',smiles,'" ', '-oinchikey')
			inchikey <- system(cmdinchikey, intern=TRUE, input=smiles, ignore.stderr=TRUE)
			cmdinchi <- paste0(babeldir, 'obabel -:"',smiles,'" ', '-oinchi')
			inchi <- system(cmdinchi, intern=TRUE, input=smiles, ignore.stderr=TRUE)
			
			##Read from Compoundlist
			smiles <- findSmiles(id)
			mass <- findMass(smiles)
			dbcas <- findCAS(id)
			dbname <- findName(id)
			if(is.na(dbname)) dbname <- ""
			if(is.na(dbcas)) dbcas <- ""
			formula <- findFormula(id)
			
			##Create 
			mbdata <- list()
			mbdata[['id']] <- id
			mbdata[['dbcas']] <- dbcas
			mbdata[['dbname']] <- dbname
			mbdata[['dataused']] <- "smiles"
			mbdata[['ACCESSION']] <- ""
			mbdata[['RECORD_TITLE']] <- ""
			mbdata[['DATE']] <- format(Sys.Date(), "%Y.%m.%d")
			mbdata[['AUTHORS']] <- getOption("RMassBank")$annotations$authors
			mbdata[['LICENSE']] <- getOption("RMassBank")$annotations$license
			mbdata[['COPYRIGHT']] <- getOption("RMassBank")$annotations$copyright
			# Confidence annotation and internal ID annotation.
			# The ID of the compound will be written like:
			# COMMENT: EAWAG_UCHEM_ID 1234
			# if annotations$internal_id_fieldname is set to "EAWAG_UCHEM_ID"
			mbdata[["COMMENT"]] <- list()
			if(findLevel(id) == "0"){
			mbdata[["COMMENT"]][["CONFIDENCE"]] <- getOption("RMassBank")$annotations$confidence_comment
            } else{
                level <- findLevel(id)
                if(level %in% c("1","1a")){
                     mbdata[["COMMENT"]][["CONFIDENCE"]] <- "Reference Standard (Level 1)"
                }
                if(level == c("2")){
                     mbdata[["COMMENT"]][["CONFIDENCE"]] <- "Probable structure, tentative identification (Level 2)"
                }
                if(level == c("2a")){
                     mbdata[["COMMENT"]][["CONFIDENCE"]] <- "Probable structure via library match, tentative identification (Level 2a)"
                }
                if(level == c("2b")){
                     mbdata[["COMMENT"]][["CONFIDENCE"]] <- "Probable structure via diagnostic evidence, tentative identification (Level 2b)"
                }
                if(level == c("3")){
                     mbdata[["COMMENT"]][["CONFIDENCE"]] <- "Tentative identification only (Level 3)"
                }
                if(level == c("3a")){
                     mbdata[["COMMENT"]][["CONFIDENCE"]] <- "Tentative identification: most likely structure (Level 3)"
                }
                if(level == c("3b")){
                     mbdata[["COMMENT"]][["CONFIDENCE"]] <- "Tentative identification: isomers possible (Level 3)"
                }
                if(level == c("3c")){
                     mbdata[["COMMENT"]][["CONFIDENCE"]] <- "Tentative identification: substance class known (Level 3)"
                }
                if(level == c("3d")){
                     mbdata[["COMMENT"]][["CONFIDENCE"]] <- "Tentative identification: best match only (Level 3)"
                }
                if(level == c("4")){
                     mbdata[["COMMENT"]][["CONFIDENCE"]] <- "Tentative identification: molecular formula only (Level 4)"
                }
                if(level == c("5")){
                     mbdata[["COMMENT"]][["CONFIDENCE"]] <- "Tentative identification: structure and formula unknown (Level 5)"
                }
            }
			mbdata[["COMMENT"]][["ID"]] <- id

			# here compound info starts
			mbdata[['CH$NAME']] <- as.list(dbname)
			
			# Currently we use a fixed value for Compound Class, since there is no useful
			# convention of what should go there and what shouldn't, and the field is not used
			# in search queries.
			mbdata[['CH$COMPOUND_CLASS']] <- getOption("RMassBank")$annotations$compound_class
			mbdata[['CH$FORMULA']] <- formula
			mbdata[['CH$EXACT_MASS']] <- mass
			mbdata[['CH$SMILES']] <- smiles
			mbdata[['CH$IUPAC']] <- inchi
			
			link <- list()
			if(dbcas != "")
			link[["CAS"]] <- dbcas
			link[["INCHIKEY"]] <- inchikey
			mbdata[['CH$LINK']] <- link
		}
		return(mbdata)
}

# Retrieve annotation data for a compound, using only babel
#' Retrieve annotation data
#' 
#' Retrieves annotation data for an unknown compound by using basic information present
#'
#' Composes the "upper part" of a MassBank record filled with chemical data
#' about the compound: name, exact mass, structure, CAS no..  
#' The instrument type is also written into this block (even
#' if not strictly part of the chemical information). Additionally, index
#' fields are added at the start of the record, which will be removed later:
#' \code{id, dbcas, dbname} from the compound list.
#' 
#' Additionally, the fields \code{ACCESSION} and \code{RECORD_TITLE} are
#' inserted empty and will be filled later on.
#' 
#' This function is used to generate the data in case a substance is unknown,
#' i.e. not enough information is present to derive anything about formulas or links
#'
#' @usage gatherDataUnknown(id, mode, retrieval)
#' @param id The compound ID.
#' @param mode \code{"pH", "pNa", "pM", "pNH4", "mH", "mM", "mFA"} for different ions 
#' 			([M+H]+, [M+Na]+, [M]+, [M+NH4]+, [M-H]-, [M]-, [M+FA]-).
#' @param retrieval A value that determines whether the files should be handled either as "standard",
#' if the compoundlist is complete, "tentative", if at least a formula is present or "unknown"
#' if the only know thing is the m/z
#' @return Returns a list of type \code{list(id= \var{compoundID}, ...,
#' 'ACCESSION' = '', 'RECORD_TITLE' = '', )} etc. %% ...
#' @author Michael Stravs, Erik Mueller
#' @seealso \code{\link{mbWorkflow}}
#' @references MassBank record format:
#' \url{http://www.massbank.jp/manuals/MassBankRecord_en.pdf}
#' @examples
#' 
#' # Gather data for compound ID 131
#' \dontrun{gatherDataUnknown(131,"pH")}
#' 
#' @export
gatherDataUnknown <- function(id, mode, retrieval){
    .checkMbSettings()
    
    ##Read from Compoundlist
    smiles <- ""
    if(retrieval == "unknown"){
        mass <- findMass(id, "unknown", mode)
        formula <- ""
    }    
    if(retrieval == "tentative"){
        mass <- findMass(id, "tentative", mode)
        formula <- findFormula(id, "tentative")
    }
    dbcas <- NA
    dbname <- findName(id)
    if(is.na(dbname)) dbname <- paste("Unknown ID:",id)
    if(is.na(dbcas)) dbcas <- ""
    

    
    ##Create 
    mbdata <- list()
    mbdata[['id']] <- id
    mbdata[['dbcas']] <- dbcas
    mbdata[['dbname']] <- dbname
    mbdata[['dataused']] <- "none"
    mbdata[['ACCESSION']] <- ""
    mbdata[['RECORD_TITLE']] <- ""
    mbdata[['DATE']] <- format(Sys.Date(), "%Y.%m.%d")
    mbdata[['AUTHORS']] <- getOption("RMassBank")$annotations$authors
    mbdata[['LICENSE']] <- getOption("RMassBank")$annotations$license
    mbdata[['COPYRIGHT']] <- getOption("RMassBank")$annotations$copyright
    # Confidence annotation and internal ID annotation.
    # The ID of the compound will be written like:
    # COMMENT: EAWAG_UCHEM_ID 1234
    # if annotations$internal_id_fieldname is set to "EAWAG_UCHEM_ID"
    mbdata[["COMMENT"]] <- list()
    if(findLevel(id) == "0"){
    mbdata[["COMMENT"]][["CONFIDENCE"]] <- getOption("RMassBank")$annotations$confidence_comment
	} else{
        level <- findLevel(id)
        if(level %in% c("1","1a")){
             mbdata[["COMMENT"]][["CONFIDENCE"]] <- "Reference Standard (Level 1)"
        }
        if(level == c("2")){
             mbdata[["COMMENT"]][["CONFIDENCE"]] <- "Probable structure, tentative identification (Level 2)"
        }
        if(level == c("2a")){
             mbdata[["COMMENT"]][["CONFIDENCE"]] <- "Probable structure via library match, tentative identification (Level 2a)"
        }
        if(level == c("2b")){
             mbdata[["COMMENT"]][["CONFIDENCE"]] <- "Probable structure via diagnostic evidence, tentative identification (Level 2b)"
        }
        if(level == c("3")){
             mbdata[["COMMENT"]][["CONFIDENCE"]] <- "Tentative identification only (Level 3)"
        }
        if(level == c("3a")){
             mbdata[["COMMENT"]][["CONFIDENCE"]] <- "Tentative identification: most likely structure (Level 3)"
        }
        if(level == c("3b")){
             mbdata[["COMMENT"]][["CONFIDENCE"]] <- "Tentative identification: isomers possible (Level 3)"
        }
        if(level == c("3c")){
             mbdata[["COMMENT"]][["CONFIDENCE"]] <- "Tentative identification: substance class known (Level 3)"
        }
        if(level == c("3d")){
             mbdata[["COMMENT"]][["CONFIDENCE"]] <- "Tentative identification: best match only (Level 3)"
        }
        if(level == c("4")){
             mbdata[["COMMENT"]][["CONFIDENCE"]] <- "Tentative identification: molecular formula only (Level 4)"
        }
        if(level == c("5")){
             mbdata[["COMMENT"]][["CONFIDENCE"]] <- "Tentative identification: structure and formula unknown (Level 5)"
        }
    }
    mbdata[["COMMENT"]][["ID"]] <- id

    # here compound info starts
    mbdata[['CH$NAME']] <- as.list(dbname)
    
    # Currently we use a fixed value for Compound Class, since there is no useful
    # convention of what should go there and what shouldn't, and the field is not used
    # in search queries.
    mbdata[['CH$COMPOUND_CLASS']] <- getOption("RMassBank")$annotations$compound_class
    mbdata[['CH$FORMULA']] <- formula
    mbdata[['CH$EXACT_MASS']] <- mass
    mbdata[['CH$SMILES']] <- ""
    mbdata[['CH$IUPAC']] <- ""
    
    link <- list()
    mbdata[['CH$LINK']] <- link

    return(mbdata)
}

# Flatten the internal tree-like representation of MassBank data to a flat table.
# Note that this limits us, in that the fields should be constant over all records!
# Therefore, e.g. the fixed number of 3 names which may be filled.
# If anybody has a cooler solution, I'll be happy to hear from you :)
#
# Note: the records from gatherData have additional information which is discarded, like
# author, copyright etc. They will be re-filled automatically when reading the file.
#' Flatten, or re-read, MassBank header blocks
#' 
#' \code{flatten} converts a list of MassBank compound information sets (as
#' retrieved by \code{\link{gatherData}}) to a flat table, to be exported into
#' an \link[=loadInfolist]{infolist}. \code{readMbdata} reads a single record
#' from an infolist flat table back into a MassBank (half-)entry.
#' 
#' Neither the flattening system itself nor the implementation are particularly
#' fantastic, but since hand-checking of records is a necessary evil, there is
#' currently no alternative (short of coding a complete GUI for this and
#' working directly on the records.)
#' 
#' @aliases flatten readMbdata
#' @usage flatten(mbdata) 
#' 
#' readMbdata(row)
#' @param mbdata A list of MassBank compound information sets as returned from
#' \code{\link{gatherData}}.
#' @param row One row of MassBank compound information retrieved from an
#' infolist.
#' @return \code{flatten} returns a matrix (not a data frame) to be written to
#' CSV.
#' 
#' \code{readMbdata} returns a list of type \code{list(id= \var{compoundID},
#' ..., 'ACCESSION' = '', 'RECORD_TITLE' = '', )} etc.
#' @author Michael Stravs
#' @seealso \code{\link{gatherData}},\code{\link{loadInfolist}}
#' @references MassBank record format:
#' \url{http://www.massbank.jp/manuals/MassBankRecord_en.pdf}
#' @examples \dontrun{
#' 	# Collect some data to flatten
#' 	ids <- c(40,50,60,70)
#'  data <- lapply(ids, gatherData)
#'  # Flatten the data trees to a table
#'  flat.table <- flatten(data)
#'  # reimport the table into a tree
#'  data.reimported <- apply(flat.table, 1, readMbdata)
#' }
#' 
#' @export
#' 
flatten <- function(mbdata)
{
  .checkMbSettings()
  
  colNames     <- names(unlist(mbdata[[1]]))
  commentNames <- colNames[grepl(x = colNames, pattern = "^COMMENT\\.")]
  
  colList <- c(
              "id",
              "dbcas",
              "dbname",
              "dataused",
              commentNames,
              #"COMMENT.CONFIDENCE",
              # Note: The field name of the internal id field is replaced with the real name
              # at "compilation" time. Therefore, functions DOWNSTREAM from compileRecord() 
              # must use the full name including the info from options("RMassBank").
              #"COMMENT.ID",
              "CH$NAME1",
              "CH$NAME2",
              "CH$NAME3",
              "CH$NAME4",
              "CH$NAME5",
              "CH$COMPOUND_CLASS",
              "CH$FORMULA",
              "CH$EXACT_MASS",
              "CH$SMILES",
              "CH$IUPAC",
              "CH$LINK.CAS",
              "CH$LINK.CHEBI",
              "CH$LINK.HMDB",
              "CH$LINK.KEGG",
              "CH$LINK.LIPIDMAPS",
              "CH$LINK.PUBCHEM",
              "CH$LINK.INCHIKEY",
              "CH$LINK.CHEMSPIDER",
	          "CH$LINK.COMPTOX"
	          )
  # make an empty data frame with the right length
  rows <- length(mbdata)
  cols <- length(colList)
  mbframe <- matrix(data=NA, nrow=rows, ncol=cols)
  colnames(mbframe) <- colList
  #browser()
  for(row in 1:rows)
  {
    # fill in all the data into the dataframe: all columns which 
    # a) exist in the target dataframe and b) exist in the (unlisted) MB record
    # are written into the dataframe.
    data <- unlist(mbdata[[row]])
	# bugfix for the case of only one name
	if(!("CH$NAME1" %in% names(data)))
		data[["CH$NAME1"]] <- data[["CH$NAME"]]
    datacols <- intersect(colList, names(data))
    mbframe[row,datacols] <- data[datacols]
  }
  return(mbframe)
  
}

# Read data from a flat-table MassBank record row and feed it into a
# MassBank tree-like record. Also, prime the ACCESSION and RECORD_TITLE fields in the
# correct position in the record.
#' @export
readMbdata <- function(row)
{
  .checkMbSettings()
  
  # Listify the table row. Lists are just cooler to work with :)
  row <- as.list(row)
  
  mbdata <- list()
  # Accession and title are added empty for now, to have them in the right place.
  # Constants are read from the options or generated.
  mbdata[['ACCESSION']] <- ""
  mbdata[['RECORD_TITLE']] <- ""
  mbdata[['DATE']] <- format(Sys.Date(), "%Y.%m.%d")
  mbdata[['AUTHORS']] <- getOption("RMassBank")$annotations$authors
  mbdata[['LICENSE']] <- getOption("RMassBank")$annotations$license
  mbdata[['COPYRIGHT']] <- getOption("RMassBank")$annotations$copyright
  if(getOption("RMassBank")$annotations$publication!="") {
    mbdata[['PUBLICATION']] <- getOption("RMassBank")$annotations$publication
  }
  commentNames <- names(row)[grepl(x = names(row), pattern = "^COMMENT\\.")]
  commentNames <- commentNames[!is.na(row[commentNames])]
  
  # Read all determined fields from the file
  # This is not very flexible, as you can see...
    colList <- c(
              commentNames,
              #"COMMENT.CONFIDENCE",
              #"COMMENT.ID",
              "CH$NAME1",
              "CH$NAME2",
              "CH$NAME3",
              "CH$NAME4",
              "CH$NAME5",
              "CH$COMPOUND_CLASS",
              "CH$FORMULA",
              "CH$EXACT_MASS",
              "CH$SMILES",
              "CH$IUPAC",
              "CH$LINK.CAS",
              "CH$LINK.CHEBI",
              "CH$LINK.HMDB",
              "CH$LINK.KEGG",
              "CH$LINK.LIPIDMAPS",
              "CH$LINK.PUBCHEM",
              "CH$LINK.INCHIKEY",
              "CH$LINK.CHEMSPIDER",
              "CH$LINK.COMPTOX")
  mbdata[["COMMENT"]] = list()
  #mbdata[["COMMENT"]][["CONFIDENCE"]] <- row[["COMMENT.CONFIDENCE"]]
  # Again, our ID field. 
  #mbdata[["COMMENT"]][["ID"]] <- row[["COMMENT.ID"]]
  mbdata[["COMMENT"]][gsub(x = commentNames, pattern = "^COMMENT\\.", replacement = "")] <- row[commentNames]
  
  names = c(row[["CH.NAME1"]], row[["CH.NAME2"]], row[["CH.NAME3"]], row[["CH.NAME4"]], row[["CH.NAME5"]])
  names = names[which(!is.na(names))]
  
  names <- gsub("'", "`", names) 
  mbdata[["CH$NAME"]] = names
  mbdata[["CH$COMPOUND_CLASS"]] = row[["CH.COMPOUND_CLASS"]]
  mbdata[["CH$FORMULA"]] = row[["CH.FORMULA"]]
  mbdata[["CH$EXACT_MASS"]] = row[["CH.EXACT_MASS"]]
  mbdata[["CH$SMILES"]] = row[["CH.SMILES"]]
  mbdata[["CH$IUPAC"]] = row[["CH.IUPAC"]]
  # Add all links and then eliminate the NA values from the tree.
  link = list()
  link[["CAS"]] = row[["CH.LINK.CAS"]]
  link[["CHEBI"]] = row[["CH.LINK.CHEBI"]]
  link[["HMDB"]] = row[["CH.LINK.HMDB"]]
  link[["KEGG"]] = row[["CH.LINK.KEGG"]]
  link[["LIPIDMAPS"]] = row[["CH.LINK.LIPIDMAPS"]]
  link[["PUBCHEM"]] = row[["CH.LINK.PUBCHEM"]]
  link[["INCHIKEY"]] = row[["CH.LINK.INCHIKEY"]]
  link[["CHEMSPIDER"]] = row[["CH.LINK.CHEMSPIDER"]]
  link[["COMPTOX"]] = row[["CH.LINK.COMPTOX"]]
  link[which(is.na(link))] <- NULL
  mbdata[["CH$LINK"]] <- link

    ## SP$SAMPLE
  if(all(nchar(row[["SP.SAMPLE"]]) > 0, row[["SP.SAMPLE"]] != "NA", !is.na(row[["SP.SAMPLE"]]), na.rm = TRUE))
    mbdata[['SP$SAMPLE']] <- row[["SP.SAMPLE"]]


  
  return(mbdata)
  
}

#' Generate peak annotation from peaklist
#' 
#' Generates the PK$ANNOTATION entry from the peaklist obtained. This function is
#' overridable by using the "annotator" option in the settings file.
#' 
#' @param annotation A peak list to be annotated. Contains columns:
#' \code{"cpdID","formula","mzFound" ,"scan","mzCalc","dppm",
#'      "dbe","mz","int","formulaCount","parentScan","fM_factor","dppmBest",
#'     "formulaMultiplicity","intrel","mzSpec"}
#' 
#' @param formulaTag The ion type to be added to annotated formulas ("+" or "-" usually)
#' 
#' @return The annotated peak table. Table \code{colnames()} will be used for the
#' 		titles (preferrably don't use spaces in the column titles; however no format is
#' 		strictly enforced by the MassBank data format.
#' 
#' @examples 
#' \dontrun{
#' annotation <- annotator.default(annotation)
#' }
#' @author Michele Stravs, Eawag <stravsmi@@eawag.ch>
#' @export
annotator.default <- function(annotation, formulaTag)
{
  if(!is.null(formulaTag))
    type <- formulaTag
  else
    type <- ""
  
  annotation <- annotation[!is.na(annotation$formula),,drop=FALSE]
  annotation <- annotation[annotation$formula != "",,drop=FALSE]
  
  annotation$formula <- paste(annotation$formula, rep(type, length(annotation$formula)), sep='')
  # Select the right columns and name them correctly for output.
  annotation <- annotation[,c("mz","formula", "formulaCount", "mzCalc", "dppm")]
  colnames(annotation) <- c("m/z", "tentative_formula", "formula_count", "mass", "error(ppm)")
  
  return(annotation)
}

#' Parse record title
#' 
#' Parses a title for a single MassBank record using the title format
#' specified in the option titleFormat. Internally used, not exported.
#' 
#' If the option is not set, a standard title format is used (for record definition
#' version 1 or 2).
#' 
#' @usage .parseTitleString(mbdata)
#' @param mbdata list
#' The information data block for the record header, as stored in
#' \code{mbdata_relisted} after loading an infolist.
#' @return A string with the title.
#' @author Michael Stravs, Eawag
#' @seealso \code{\link{buildRecord}}
#' @references MassBank record format:
#' \url{http://www.massbank.jp/manuals/MassBankRecord_en.pdf}
#' @examples
#' \dontrun{
#' 		# used in buildRecord()
#' 		title <- .parseTitleString(mbdata)
#' }
#' 
#' 
#' 
.parseTitleString <- function(mbdata)
{
	
	varlist <- getOption("RMassBank")$titleFormat
	
	# Set the standard title format.
	if(is.null(varlist))
	{
		if(getOption("RMassBank")$use_version == 2)
		{
			varlist <- c(
					"{CH$NAME}",
					"{AC$INSTRUMENT_TYPE}",
					"{AC$MASS_SPECTROMETRY: MS_TYPE}",
					"CE: {RECORD_TITLE_CE}",
					"R={AC$MASS_SPECTROMETRY: RESOLUTION}",
					"{MS$FOCUSED_ION: PRECURSOR_TYPE}"
			)
		}
		else
		{
			varlist <- c(
					"{CH$NAME}",
					"{AC$INSTRUMENT_TYPE}",
					"{AC$ANALYTICAL_CONDITION: MS_TYPE}",
					"CE: {RECORD_TITLE_CE}",
					"R={AC$ANALYTICAL_CONDITION: RESOLUTION}",
					"{MS$FOCUSED_ION: PRECURSOR_TYPE}"
			)
		}
	}
  
	
	# Extract a {XXX} argument from each title section.
	# check that every title has one and only one match
	args <- regexec("\\{(.*)\\}", varlist)
	arglist <- regmatches(varlist, args)
	if(any(unlist(lapply(arglist, length)) != 2))
		stop("Title format is incorrectly specified: a section with not exactly 1 parameters")
	
	parsedVars <- lapply(varlist, function(var)
			{
				# Extract the specified parameter inside the {}.
				# I.e. from a string like "R={BLA: BLUB}" return "BLA: BLUB"
				args <- regexec("\\{(.*)\\}", var)
				arg <- regmatches(var, args)[[1]][[2]]
				# Split the parameter by colon if necessary
				splitVar <- strsplit(arg, ": ")[[1]]
				# Read the parameter value from the record
				if(length(splitVar) == 2)
					replaceVar <- mbdata[[splitVar[[1]]]][[splitVar[[2]]]]
				else if(length(splitVar) ==  1)
					replaceVar <- mbdata[[splitVar]]
				else
					stop(paste(
									"Title format is incorrectly specified:", var)
					)
				# Fix problems: NULL returns
				if(is.null(replaceVar))
					replaceVar <- ""
				# Fix problems: Names will have >= 1 match. Take the first
				if(length(replaceVar) > 1)
					replaceVar <- replaceVar[[1]]
                
                # Fix problems: Unknowns might have no name
                if(!length(replaceVar)){
                    replaceVar <- ""
                }
                
				# Substitute the parameter value into the string
				parsedVar <- sub("\\{(.*)\\}", replaceVar, var)	
				return(parsedVar)
			})
	title <- paste(parsedVars, collapse="; ")
	return(title)
}


# This converts the tree-like list (as obtained e.g. from compileRecord())
# into a plain text array, which can then be dumped to a file suitable for 
# MassBank upload.

#' Write MassBank record into character array
#' 
#' Writes a MassBank record in list format to a text array.
#' 
#' The function is a general conversion tool for the MassBank format; i.e. the
#' field names are not fixed. \code{mbdata} must be a named list, and the
#' entries can be as follows: \itemize{
#'  \item A single text line:
#' 
#' \code{'CH\$EXACT_MASS' = '329.1023'}
#' 
#'  is written as
#' 
#'  \code{CH\$EXACT_MASS: 329.1023} 
#' \item A character array:
#' 
#'  \code{'CH\$NAME' = c('2-Aminobenzimidazole', '1H-Benzimidazol-2-amine')} 
#' 
#' is written as
#' 
#' \code{CH\$NAME: 2-Aminobenzimidazole}
#' 
#' \code{CH\$NAME: 1H-Benzimidazol-2-amine}
#' 
#' \item A named list of strings: 
#' 
#' 	\code{'CH\$LINK' = list('CHEBI' = "27822", "KEGG" = "C10901")} 
#' 
#' is written as 
#' 
#' \code{CH\$LINK: CHEBI 27822}
#' 
#' \code{CH\$LINK: KEGG C10901} 
#' 
#' \item A data frame (e.g. the peak table) is written as specified in
#' the MassBank record format (Section 2.6.3): the column names are used as
#' headers for the first line, all data rows are printed space-separated. 
#' }
#' 
#' @usage toMassbank(mbdata)
#' @param mbdata A MassBank record in list format.
#' @return The result is a text array, which is ready to be written to the disk
#' as a file.
#' @note The function iterates over the list item names. \bold{This means that
#' duplicate entries in \code{mbdata} are (partially) discarded!} The correct
#' way to add them is by making a character array (as specified above): Instead
#' of \code{'CH\$NAME' = 'bla', 'CH\$NAME' = 'blub'} specify \code{'CH\$NAME' =
#' c('bla','blub')}.
#' @author Michael Stravs
#' @seealso \code{\link{buildRecord}}, \code{\link{mbWorkflow}}
#' @references MassBank record format:
#' \url{http://www.massbank.jp/manuals/MassBankRecord_en.pdf}
#' @examples
#' \dontrun{
#' # Read just the compound info skeleton from the Internet for some compound ID
#' id <- 35
#' mbdata <- gatherData(id)
#' #' # Export the mbdata blocks to line arrays 
#' # (there is no spectrum information, just the compound info...)
#' mbtext <- toMassbank(mbdata)
#' }
#' 
#' 
#' @export
setGeneric("toMassbank", function(o, ...) standardGeneric("toMassbank"))

#' @export
setMethod("toMassbank", "RmbSpectraSet", function(o, addAnnotation = getOption("RMassBank")$add_annotation)
    {
      lapply(o@children, function(s) toMassbank(s, addAnnotation))
    })

#' @export
setMethod("toMassbank", "RmbSpectrum2", function(o, addAnnotation = getOption("RMassBank")$add_annotation)
    {
      .toMassbank(o, addAnnotation)
    })

.toMassbank <- function (s, addAnnotation = getOption("RMassBank")$add_annotation)
{
  
  peaks <- getData(s)
  # check that peaks were normalized
  if(!("intrel" %in% colnames(peaks)))
  {
    s <- normalize(s, slot="intrel")
    peaks <- getData(s)
  }  
  
  # Keep only peaks with relative intensity >= 1 o/oo, since the MassBank record
  # makes no sense otherwise. Also, keep only the columns needed in the output.
  peaks <- peaks[ peaks$intrel >= 1,,drop=FALSE]	
  
  peaks$mz <- round(peaks$mz, 4)
  # Also format the other values, which are used in the annotation
  peaks$dppm <- round(peaks$dppm, 2)
  peaks$mzCalc <- round(peaks$mzCalc, 4)
  peaks$intensity <- round(peaks$intensity, 1)
  
  # Get polarity from Spectrum2 now!
  formulaTag <- ""
  if(s@polarity == 1) formulaTag <- "+"
  if(s@polarity == 0) formulaTag <- "-"
  # if polarity is -1, leave it unspecified. the "specs" seem to be 1 for +, 0 for - and -1 for ???
  # (when reading mzML I often get -1, when reading mzXML I get 1 and 0 respectively)
  
  annotator <- getOption("RMassBank")$annotator
  if(is.null(annotator))
    annotator <- "annotator.default"
  
  annotation <- do.call(annotator, list(annotation= peaks, formulaTag = formulaTag))
  
  peaks <- peaks[,c("mz", "intensity", "intrel")]
  peaks <- unique(peaks)
  # Name the columns correctly.
  colnames(peaks) <- c("m/z", "int.", "rel.int.")
  peaknum <- nrow(peaks)
  
  mbdata <- s@info
  
  mbdata[["PK$SPLASH"]] <- list(SPLASH = getSplash(peaks[,c("m/z", "int.")]))
  
  # Annotation:
  if(addAnnotation && (nrow(annotation) > 0))
    mbdata[["PK$ANNOTATION"]] <- annotation
  
  # Peak table
  mbdata[["PK$NUM_PEAK"]] <- peaknum
  mbdata[["PK$PEAK"]] <- peaks
  
  # mbf is an array of lines and count is the line counter.
  # Very old-school, but it works. :)
  mbf <- character(0)
  count <- 1
  lapply(names(mbdata), function(entry)
    {
      # If entry is a char line, add it to the file.
      # If it is a named sublist, add each subentry with name
      # If it is an unnamed sublist, add each subentry without name
      # if it is a dataframe, write in PEAKS mode
    
      # Note: this is were I liked "lapply" a little too much. "for" would
      # be more idiomatic, and wouldn't need the <<- assignments.
      
      # Data frame: table mode. A header line and one space-separated line for
      # each data frame row.
      if(is.data.frame(mbdata[[entry]]))
      {
        mbf[[count]] <<- paste(entry,": " ,
                               paste(colnames(mbdata[[entry]]), collapse=" "),
                               sep='')
        count <<- count+1
        for(row in 1:nrow(mbdata[[entry]]))
        {
          mbf[[count]] <<- paste("  ", 
                                 paste(mbdata[[entry]][row,],collapse=" "), 
                                 sep="")
          count <<- count+1
        }
        #browser()
      }
      # List with named items: Write every entry like CH$LINK: CAS 12-345-678
      else if(is.list(mbdata[[entry]]) & !is.null(names(mbdata[[entry]])))
      {
        
        lapply(names(mbdata[[entry]]), function(subentry)
        {
          if(subentry != "SPLASH"){
            mbf[[count]] <<- paste(entry,": ",subentry, " ", mbdata[[entry]][[subentry]], sep='')
          } else {
            mbf[[count]] <<- paste(entry,": ", mbdata[[entry]][[subentry]], sep='')
          }
          #print(mbf)
          count <<- count + 1
        })
      }
      # Array (or list) of unnamed items: Write every entry like CH$NAME: Paracetamol
      # (iterative entry without subindices)
      else if (length(mbdata[[entry]]) > 1 & is.null(names(mbdata[[entry]])))
      {
        lapply(mbdata[[entry]], function(subentry)
        {
          mbf[[count]] <<- paste(entry,": ",subentry, sep='')
          #print(mbf)
          count <<- count + 1
        })   
      }
      # Length is 1: just write the entry like PK$NUM_PEAKS: 131
      else
      {
        mbf[[count]] <<- paste(entry,": ",mbdata[[entry]], sep='')
        count <<- count + 1
      }
    }
    ) # End of lapply block (per child spectrum)
  # Add mandatory EOF marker
  mbf[[count]] <- "//"
  return(mbf)
}

# Exports compiled and massbanked spectra, with their associated molfiles, to physical files.
# "compiled" is still used here, because we need an accessible accession number.
# In the plain text arrays, the accession number is already "hidden".
# compiled: is ONE "compiled" entry, i.e. ONE compound with e.g. 14 spectra.
# files: is a return value from lapply(toMassbank), i.e. contains 14 plain-text arrays
#  (for a 14-spectra method)
# molfile: a molfile from createMolfile

#' Export internally stored MassBank data to files
#' 
#' Exports MassBank recfile data arrays and corresponding molfiles to physical
#' files on hard disk, for one compound.
#' 
#' The data from \code{compiled} is still used here, because it contains the
#' "visible" accession number. In the plain-text format contained in
#' \code{files}, the accession number is not "accessible" anymore since it's in
#' the file.
#' 
#' @usage exportMassbank(compiled, files, molfile)
#' @param compiled \code{RmbSpectraSet}
#' the spectra of one compound for which files should be exported
#' @param molfile A molfile from \code{\link{createMolfile}}
#' @return No return value.
#' @note An improvement would be to write the accession numbers into
#' \code{names(compiled)} and later into \code{names(files)} so \code{compiled}
#' wouldn't be needed here anymore. (The compound ID would have to go into
#' \code{names(molfile)}, since it is also retrieved from \code{compiled}.)
#' @author Michael Stravs
#' @seealso \code{\link{createMolfile}}, \code{\link{toMassbank}},
#' \code{\link{mbWorkflow}}
#' @references MassBank record format:
#' \url{http://www.massbank.jp/manuals/MassBankRecord_en.pdf}
#' @export
exportMassbank <- function(compiled, molfile = NULL)
{
  exportMassbank_recdata(
    compiled,   
    recDataFolder = file.path(getOption("RMassBank")$annotations$entry_prefix, "recdata")
  )
  if(!is.null(molfile)) {
    exportMassbank_moldata(
      compiled,
      molfile,
      molDataFolder = file.path(getOption("RMassBank")$annotations$entry_prefix, "moldata")
    )
    }
}

exportMassbank_recdata <- function(compiled, recDataFolder)
{
  #mb@mbfiles <- lapply(mb@compiled_ok, function(cpd) toMassbank(cpd, mb@additionalPeaks))
  
  files <- toMassbank(compiled)
  names(files) <- lapply(compiled@children, function(c) c@info[["ACCESSION"]] )
  
  molnames <- c()
  for(file in seq_len(length(files)))
  {
    # Read the accession no. from the corresponding "compiled" entry
    filename <- names(files)[[file]]
    # use this accession no. as filename
    filename <- paste(filename, ".txt", sep="")
    filePath <- file.path(recDataFolder,filename)
    write(files[[file]], filePath)
  }
}

exportMassbank_moldata <- function(compiled, molfile, molDataFolder)
{
  # Use internal ID for naming the molfiles
  if(findLevel(compiled@id,TRUE)=="standard"){
    molname <- sprintf("%04d", as.numeric(compiled@id))
    molname <- paste(molname, ".mol", sep="")
    write(molfile, file.path(molDataFolder,molname))
  }
}





# Makes a list.tsv with molfile -> massbank ch$name attribution.

#' Write list.tsv file
#' 
#' Makes a list.tsv file in the "moldata" folder.
#' 
#' Generates the list.tsv file which is needed by MassBank to connect records with
#' their respective molfiles. The first compound name is linked to a mol-file with
#' the compound ID (e.g. 2334.mol for ID 2334).
#' 
#' @param compiled list of \code{RmbSpectraSet}
#' compiled spectra for multiple compounds (one \code{RmbSpectraSet} each).
#' @return No return value.
#' @author Michael A. Stravs, Eawag <michael.stravs@eawag.ch>
#' @export
makeMollist <- function(compiled)
{
  # For every "compiled" entry (here, compiled is not one "compiled" entry but the total
  # list of all compiled spectra), extract the uppermost CH$NAME and the ID (from the
  # first spectrum.) Make the ID into 0000 format.
  
  emptySpectra <- unlist(lapply(compiled, function(cpd) length(cpd@children) == 0))
  compiled <- compiled[!emptySpectra]
  
  tsvlist <- t(sapply(compiled, function(entry)
    {
    name <- entry@children[[1]]@info[["CH$NAME"]][[1]]
    id <- sprintf("%04d", as.numeric(entry@id))
    molfilename <- paste(id,".mol",sep='')
    return(c(name,molfilename))
  }))
    
    IDs <- sapply(compiled, function(entry) return( sprintf("%04d", as.numeric(
                      entry@id))))
    level <- sapply(IDs, findLevel, compact=TRUE)
    validentries <- which(level == "standard")
  # Write the file with the 
    write.table(tsvlist[validentries,], 
              paste(getOption("RMassBank")$annotations$entry_prefix,"/moldata/list.tsv", sep=''),
              quote = FALSE,
              sep="\t",
              row.names=FALSE,
              col.names=FALSE
              )
}


# Load a dataframe or file into additional_peaks (or add additional points in there.)
# The columns cpdID, scan, mzFound, int, OK are mandatory. OK=1 means that the peaks
# will be added into the spectrum. mzFound and int will be taken for the table.
# No annotation will be written.
# Add peaks to the spectra by hand

#' Add additional peaks to spectra
#' 
#' Loads a table with additional peaks to add to the MassBank spectra. Required
#' columns are \code{cpdID, scan, int, mzFound, OK}.
#' 
#' All peaks with OK=1 will be included in the spectra.
#' 
#' @usage addPeaks(mb, filename_or_dataframe)
#' @param mb The \code{mbWorkspace} to load the peaks into.
#' @param filename_or_dataframe Filename of the csv file, or name of the R
#' dataframe containing the peaklist.
#' @return The \code{mbWorkspace} with loaded additional peaks.
#' @author Michael Stravs
#' @seealso \code{\link{mbWorkflow}}
#' @examples
#' 
#' 	\dontrun{addPeaks("myrun_additionalPeaks.csv")}
#' 
#' @export 
addPeaks <- function(mb, filename_or_dataframe)
{
	
	errorvar <- 0
	currEnvir <- environment()
	d <- 1
	
	if(is.data.frame(filename_or_dataframe))
		df <- filename_or_dataframe
	else
	tryCatch(
		df <- read.csv(filename_or_dataframe),
		error=function(e){
		currEnvir$errorvar <- 1
	})
	# I change your heuristic fix to another heuristic fix, because I will have to test for a column name change...
	
	if(!errorvar){
	
		if(ncol(df) < 2)
			df <- read.csv(filename_or_dataframe, sep=";")
		# here: the column int was renamed to intensity, and we need to be able to read old files. sorry.
		if(!("intensity" %in% colnames(df)) & ("int" %in% colnames(df)))
			df$intensity <- df$int
		
		cols <- c("cpdID", "scan", "mzFound", "intensity", "OK")
		n <- colnames(df)
		# Check if comma-separated or semicolon-separated
		d <- setdiff(cols, n)
		if(length(d)>0){
			stop("Some columns are missing in the additional peak list. Needs at least cpdID, scan, mzFound, intensity, OK.")
		}
	}
	
	culled_df <- df[,c("cpdID", "scan", "mzFound", "intensity", "OK")]
	
	
	if(nrow(mb@additionalPeaks) == 0)
		mb@additionalPeaks <- culled_df
	else
		mb@additionalPeaks <- rbind(mb@additionalPeaks, culled_df)
	return(mb)
}



gatherDataMinimal.cpd <- function(cpd){
  
  ##Read from Compoundlist
  if(length(cpd@smiles) == 1) smiles <- cpd@smiles
  else
    smiles <- ""
  
  ##Create 
  mbdata <- list()
  mbdata[['ACCESSION']] <- ""
  mbdata[['RECORD_TITLE']] <- ""
  mbdata[['DATE']] <- format(Sys.Date(), "%Y.%m.%d")
  # Confidence annotation and internal ID annotation.
  # The ID of the compound will be written like:
  # COMMENT: EAWAG_UCHEM_ID 1234
  # if annotations$internal_id_fieldname is set to "EAWAG_UCHEM_ID"
  if(length(cpd@id) > 0)
    mbdata[["COMMENT"]][["ID"]] <- cpd@id
  
  # here compound info starts
  mbdata[['CH$NAME']] <- cpd@name
  
  # Currently we use a fixed value for Compound Class, since there is no useful
  # convention of what should go there and what shouldn't, and the field is not used
  # in search queries.
  mbdata[['CH$FORMULA']] <- cpd@formula
  mbdata[['CH$EXACT_MASS']] <- round(findMz.formula(cpd@formula, "")$mzCenter, 4)
  
  if(cpd@smiles != "")
    mbdata[['CH$SMILES']] <- cpd@smiles
  
  link <- list()
  mbdata[['CH$LINK']] <- link

  return(mbdata)
}



gatherDataMinimal.spectrum <- function(spectrum){
  
  ##Read from Compoundlist
  if(length(cpd@smiles) == 1) smiles <- cpd@smiles
  else
    smiles <- ""
  
  ##Create 
  mbdata <- list()
  mbdata[['ACCESSION']] <- ""
  mbdata[['RECORD_TITLE']] <- ""
  mbdata[['DATE']] <- format(Sys.Date(), "%Y.%m.%d")
  # Confidence annotation and internal ID annotation.
  # The ID of the compound will be written like:
  # COMMENT: EAWAG_UCHEM_ID 1234
  # if annotations$internal_id_fieldname is set to "EAWAG_UCHEM_ID"
  
  # here compound info starts
  mbdata[['CH$NAME']] <- paste("parent", spectrum@precursorMz, "at RT", spectrum@rt, "- CE", spectrum@collisionEnergy) 
  
  # Currently we use a fixed value for Compound Class, since there is no useful
  # convention of what should go there and what shouldn't, and the field is not used
  # in search queries.
  
  return(mbdata)
}



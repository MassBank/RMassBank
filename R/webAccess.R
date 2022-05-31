#' @import XML RCurl rjson httr
NULL
## library(XML)
## library(httr)
## library(jsonlite)


retrieveDataWithRetry <- function(url, timeout, maximumNumberOfRetries = 5, retryDelayInSeconds = 3){
  #data <- getURL(URLencode(url), timeout=8)
  
  data <- NULL
  queryIsSuccessful <- FALSE
  numberOfRetries <- 0
  while(!queryIsSuccessful & numberOfRetries < maximumNumberOfRetries){
    data <- tryCatch(
      expr = {
        #data <- getURL(url = url, timeout = timeout)
        res <- GET(URLencode(url))
        data <- httr::content(res, type="text", encoding="UTF-8")
        
        queryIsSuccessful <- TRUE
        data
      },
      warning=function(w){
        numberOfRetries <<- numberOfRetries + 1
        if(RMassBank.env$verbose.output)
          cat(paste("### Warning ### Web query failed (", numberOfRetries, " / ", maximumNumberOfRetries, ") for url '", url, "' because of warning '", w, "'\n", sep = ""))
        if(numberOfRetries < maximumNumberOfRetries)
          Sys.sleep(time = retryDelayInSeconds)
      },
      error=function(e){
        numberOfRetries <<- numberOfRetries + 1
        if(RMassBank.env$verbose.output)
          cat(paste("### Warning ### Web query failed (", numberOfRetries, " / ", maximumNumberOfRetries, ") for url '", url, "' because of error '", e, "'\n", sep = ""))
        if(numberOfRetries < maximumNumberOfRetries)
          Sys.sleep(time = retryDelayInSeconds)
      }
    )
  }
  
  return(data)
}

#' Retrieve information from Cactus
#' 
#' Retrieves information from the Cactus Chemical Identifier Resolver
#' (PubChem).
#' 
#' It is not necessary to specify in which format the \code{identifier} is.
#' Somehow, cactus does this automatically.
#' 
#' @usage getCactus(identifier, representation)
#' @param identifier Any identifier interpreted by the resolver, e.g. an InChI
#' key or a SMILES code.
#' @param representation The desired representation, as required from the
#' resolver. e.g. \code{stdinchikey}, \code{chemspider_id}, \code{formula}...
#' Refer to the webpage for details.
#' @return The result of the query, in plain text. Can be NA, or one or
#' multiple lines (character array) of results.
#' @note Note that the InChI key is retrieved with a prefix (\code{InChIkey=}),
#' which must be removed for most database searches in other databases (e.g.
#' CTS).
#' @author Michael Stravs
#' @seealso \code{\link{getCtsRecord}}, \code{\link{getPcId}}
#' @references cactus Chemical Identifier Resolver:
#' \url{http://cactus.nci.nih.gov/chemical/structure}
#' @examples
#' 
#' # Benzene:
#' getCactus("C1=CC=CC=C1", "cas")
#' getCactus("C1=CC=CC=C1", "stdinchikey")
#' getCactus("C1=CC=CC=C1", "chemspider_id")
#' 
#' @export 
#' 
#' 
getCactus <- function(identifier,representation){
  identifier <- gsub('#', '%23', identifier)
  ret <- tryCatch(httr::GET(paste("https://cactus.nci.nih.gov/chemical/structure/",
                            URLencode(identifier), "/", representation, sep = "")),
                  error = function(e) NA)
  if (all(is.na(ret)))
    return(NA)
  if (ret["status_code"] == 404)
    return(NA)
  ret <- httr::content(ret)
  return(unlist(strsplit(ret, "\n")))
  
}

#' Search Pubchem CID
#' 
#' Retrieves PubChem CIDs for a search term.
#' 
#' Only the first result is returned currently. \bold{The function should be
#' regarded as experimental and has not thoroughly been tested.}
#' 
#' @usage getPcId(query, from = "inchikey")
#' @param query ID to be converted
#' @param from Type of input ID
#' @return The PubChem CID (in string type).
#' @author Michael Stravs, Erik Mueller
#' @seealso \code{\link{getCtsRecord}}, \code{\link{getCactus}}
#' @references PubChem search: \url{http://pubchem.ncbi.nlm.nih.gov/} 
#' 
#' Pubchem REST:
#' \url{https://pubchem.ncbi.nlm.nih.gov/pug_rest/PUG_REST.html}
#' @examples
#' getPcId("MKXZASYAUGDDCJ-NJAFHUGGSA-N")
#' 
#' @export
getPcId <- function(query, from = "inchikey")
{
	baseURL <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"
	url <- paste(baseURL, from, query, "description", "json", sep="/")
	
	errorvar <- 0
	currEnvir <- environment()
	
	tryCatch(
    {#		data <- getURL(URLencode(url),timeout=8),
	    res <- GET(URLencode(url))
	    data <- httr::content(res, type="text", encoding="UTF-8")
    },
		error=function(e){
		currEnvir$errorvar <- 1
	})
	
	if(errorvar){
		return(NA)
	}
	
	# This happens if the InChI key is not found:
	r <- fromJSON(data)
	
	if(!is.null(r$Fault))
	return(NA)
	
	titleEntry <- which(unlist(lapply(r$InformationList$Information, function(i) !is.null(i$Title))))
	
	titleEntry <- titleEntry[which.min(sapply(titleEntry, function(x)r$InformationList$Information[[x]]$CID))]

	PcID <- r$InformationList$Information[[titleEntry]]$CID
	
	if(is.null(PcID)){
		return(NA)
	} else{
		return(PcID)
	}
}

# The following function is unfinished.
# getPcRecord <- function(pcid)
# {
#   baseUrl <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
#   term <- paste(baseUrl, "esummary.fcgi?db=pccompound&id=", URLencode(as.character(pcid)), 
#                 
#                 sep='')
#   ret <- getURL(term)
#   xml <- xmlParseDoc(ret,asText=TRUE)
#   browser()
# }


# Note: some of the CHEBI codes returned are erroneous. (When the entry in 
# CTS starts with "CHEBI:" instead of just the number, the XML output)
# Also, there is no ChemSpider ID in the XML output, unfortunately.








#' Retrieve information from CTS
#' 
#' Retrieves a complete CTS record from the InChI key.
#' 
#' @usage getCtsRecord(key)
#' 
#' @param key The InChI key. 
#' @return Returns a list with all information from CTS: \code{inchikey, 
#' 	inchicode, formula, exactmass} contain single values. \code{synonyms} contains
#' an unordered list of scored synonyms (\code{type, name, score}, where \code{type}
#' indicates either a normal name or a specific IUPAC name, see below).
#'  \code{externalIds} contains an unordered list of identifiers of the compound in 
#' various databases (\code{name, value}, where \code{name} is the database name and
#' \code{value} the identifier in that database.)
#' 
#' @note Currently, the CTS results are still incomplete; the name scores are all 0,
#' formula and exact mass return zero.
#' @references Chemical Translation Service:
#' \url{https://cts.fiehnlab.ucdavis.edu}
#' 
#' @examples
#' data <- getCtsRecord("UHOVQNZJYSORNB-UHFFFAOYSA-N")
#' # show all synonym "types"
#' types <- unique(unlist(lapply(data$synonyms, function(i) i$type)))
#' \dontrun{print(types)}
#' 
#' @author Michele Stravs, Eawag <stravsmi@@eawag.ch>
#' @export
getCtsRecord <- function(key)
{
	baseURL <- "https://cts.fiehnlab.ucdavis.edu/service/compound/"
	
	errorvar <- 0
	currEnvir <- environment()
	
	##tryCatch a CTS timeout
	##
	tryCatch(
		{
			#data <- getURL(paste0(baseURL,key), timeout=10)
			url <- paste0(baseURL,key)
			res <- GET(URLencode(url))
			data <- httr::content(res, type="text", encoding="UTF-8")
		},
		error=function(e){
			currEnvir$errorvar <- 1
		}
	)
  
	if(errorvar){
		warning("CTS seems to be currently unavailable or incapable of interpreting your request")
		return(NULL)
	}

	r <- fromJSON(data)
	if(length(r) == 1)
		if(r == "You entered an invalid InChIKey")
			return(list())
	return(r)
}

#' Convert a single ID to another using CTS.
#' 
#' @usage getCtsKey(query, from = "Chemical Name", to = "InChIKey")
#' @param query ID to be converted
#' @param from Type of input ID
#' @param to Desired output ID 
#' @return An unordered array with the resulting converted key(s). 
#' 
#' @examples 
#' k <- getCtsKey("benzene", "Chemical Name", "InChIKey")
#' @author Michele Stravs, Eawag <stravsmi@@eawag.ch>
#' @export
getCtsKey <- function(query, from = "Chemical Name", to = "InChIKey")
{
	baseURL <- "https://cts.fiehnlab.ucdavis.edu/service/convert"
	url <- paste(baseURL, from, to, query, sep='/')
	errorvar <- 0
	currEnvir <- environment()
	
	##tryCatch a CTS timeout
	##
	tryCatch(
		{
			#data <- getURL(URLencode(url), timeout=10)
			res <- GET(URLencode(url))
			data <- httr::content(res, type="text", encoding="UTF-8")
		},
		error=function(e){
			currEnvir$errorvar <- 1
		}
	)
  
	if(errorvar){
		warning("CTS seems to be currently unavailable or incapable of interpreting your request")
		return(NULL)
	}

	if(res$status_code != 200){
	  warning(paste("CTS has return code", res$status_code))
	  return(NULL)
	}
	
	r <- fromJSON(data)
	if(length(r) == 0)
		return(NULL)
	else
	{
		# read out the results in simplest form:
		results <- unlist(lapply(r, function(row) row$result))
		return(results)
	}
}

#' Select a subset of external IDs from a CTS record.
#' 
#' @usage CTS.externalIdSubset(data, database)
#' @param data The complete CTS record as retrieved by \code{\link{getCtsRecord}}. 
#' @param database The database for which keys should be returned. 
#' @return Returns an array of all external identifiers stored in the record for the
#' given database.
#' 
#' @examples 
#' 
#' \dontrun{
#' # Return all CAS registry numbers stored for benzene.
#' data <- getCtsRecord("UHOVQNZJYSORNB-UHFFFAOYSA-N")
#' cas <- CTS.externalIdSubset(data, "CAS")
#' } 
#' 
#' @author Michele Stravs, Eawag <stravsmi@@eawag.ch>
#' @export
CTS.externalIdSubset <- function(data, database)
{
	select <- which(unlist(lapply(data$externalIds, function(id)
							{
								id[["name"]] == database
							})))
	keyEntries <- data$externalIds[select]
	keys <- unlist(lapply(keyEntries, function(e) e[["value"]]))
}

#' Find all available databases for a CTS record
#' 
#' @usage CTS.externalIdTypes(data)
#' @param data The complete CTS record as retrieved by \code{\link{getCtsRecord}}.  
#' @return Returns an array of all database names for which there are external 
#' identifiers stored in the record.
#' 
#' @examples 
#' 
#' \dontrun{
#' # Return all databases for which the benzene entry has
#' # links in the CTS record.
#' 
#' data <- getCTS("UHOVQNZJYSORNB-UHFFFAOYSA-N")
#' databases <- CTS.externalIdTypes(data)
#' } 
#' 
#' @author Michele Stravs, Eawag <stravsmi@@eawag.ch>
#' @export
CTS.externalIdTypes <- function(data)
{
	unique(unlist(lapply(data$externalIds, function(id)
							{
								id[["name"]]
							})))
}

.pubChemOnline <- function(){
	baseURL <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"
	url <- paste(baseURL, "inchikey", "QEIXBXXKTUNWDK-UHFFFAOYSA-N", "description", "json", sep="/")
	
	errorvar <- 0
	currEnvir <- environment()
	tryCatch(
		{#ret <- getURL(URLencode(url), timeout=8),
	    res <- GET(URLencode(url))
	    ret <- httr::content(res, type="text", encoding="UTF-8")
	  },
	  error=function(e){
		currEnvir$errorvar <- 1
	})
  
  if(errorvar){
	warning("Pubchem is currently offline")
	return(FALSE)
  } else{
	return(TRUE)
  }
}



getPcCHEBI <- function(query, from = "inchikey")
{
	# Get the JSON-Data from Pubchem
	baseURL <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"
	url <- paste(baseURL, from, query, "synonyms", "json", sep="/")
	errorvar <- 0
	currEnvir <- environment()
	
	tryCatch(
		{#data <- getURL(URLencode(url),timeout=8),
		  res <- GET(URLencode(url))
		  data <- httr::content(res, type="text", encoding="UTF-8")
		},
		error=function(e){
		currEnvir$errorvar <- 1
	})
	
	if(errorvar){
		return(NA)
	}
	
	r <- fromJSON(data)
	
	# This happens if the InChI key is not found:
	if(!is.null(r$Fault))
	return(NA)
	
	# Find the entries which contain Chebi-links
	synonymEntry <- which(unlist(lapply(r$InformationList$Information, function(i) !is.null(i$Synonym))))
	synonymList <- r$InformationList$Information[[synonymEntry]]$Synonym
	matchChebi <- which(grepl("CHEBI:", synonymList, fixed=TRUE))
	
	# It doesn't matter if the db is down or if chebi isn't found, so return NA also
	if(length(matchChebi) == 0){
		return (NA) 
	} else {
		return (sapply(matchChebi, function(x) synonymList[[x]]))
	}
}

#' Retrieve the Chemspider ID for a given compound
#' 
#' Given an InChIKey, this function queries the chemspider web API to retrieve
#' the Chemspider ID of he compound with that InChIkey. 
#'
#' @usage getCSID(query)
#'
#' @param query The InChIKey of the compound 
#' @return Returns the chemspide
#' 
#' @examples 
#' 
#' \dontrun{
#' # Return all CAS registry numbers stored for benzene.
#' data <- getCtsRecord("UHOVQNZJYSORNB-UHFFFAOYSA-N")
#' cas <- CTS.externalIdSubset(data, "CAS")
#' } 
#' 
#' @author Michele Stravs, Eawag <stravsmi@@eawag.ch>
#' @author Erik Mueller, UFZ <erik.mueller@@ufz.de>
#' @export
getCSID <- function(query)
{
	baseURL <- "http://www.chemspider.com/InChI.asmx/InChIKeyToCSID?inchi_key="
	url <- paste0(baseURL, query)
	
	#errorvar <- 0
	#currEnvir <- environment()
	#
	#tryCatch(
	#	data <- getURL(URLencode(url), timeout=8),
	#	error=function(e){
	#	currEnvir$errorvar <- 1
	#})
	#
	#if(errorvar){
	#	warning("Chemspider is currently offline")
	#	return(NA)
	#}
	
	data <- retrieveDataWithRetry(url = URLencode(url), timeout=8)
	if(is.null(data)){
		warning("Chemspider is currently offline")
		return(NA)
	}
	
	xml <- xmlParseDoc(data,asText=TRUE)
	# the returned XML document contains only the root node called "string" which contains the correct CSID
	idNodes <- getNodeSet(xml, "/")
	id <- xmlValue(idNodes[[1]])
	return(id)
} 

##This function returns a sensible name for the compound
getPcSynonym <- function (query, from = "inchikey")
{
	# Get the JSON-Data from Pubchem
	baseURL <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"
	url <- paste(baseURL, from, query, "description", "json", sep="/")
	
	errorvar <- 0
	currEnvir <- environment()
	
	tryCatch(
		{#data <- getURL(URLencode(url),timeout=8)
		  res <- GET(URLencode(url))
		  data <- httr::content(res, type="text", encoding="UTF-8")
		},
		error=function(e){
		currEnvir$errorvar <- 1
	})
	
	if(errorvar){
		return(NA)
	}
	
	r <- fromJSON(data)
	
	# This happens if the InChI key is not found:
	if(!is.null(r$Fault))
	return(NA)
	
	# Find the synonym
	
	titleEntry <- which(unlist(lapply(r$InformationList$Information, function(i) !is.null(i$Title))))
	
	titleEntry <- titleEntry[which.min(sapply(titleEntry, function(x)r$InformationList$Information[[x]]$CID))]
	
	title <- r$InformationList$Information[[titleEntry]]$Title
	
	if(is.null(title)){
		return(NA)
	} else{
		return(title)
	}
} 


##A function to retrieve a IUPAC Name from Pubchem
getPcIUPAC <- function (query, from = "inchikey")
{
	# Get the JSON-Data from Pubchem
	baseURL <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"
	url <- paste(baseURL, from, query, "record", "json", sep="/")
	
	errorvar <- 0
	currEnvir <- environment()
	
	tryCatch(
		{#data <- getURL(URLencode(url),timeout=8)
		  res <- GET(URLencode(url))
		  data <- httr::content(res, type="text", encoding="UTF-8")
		},
		error=function(e){
		currEnvir$errorvar <- 1
	})
	
	if(errorvar){
		return(NA)
	}
	
	r <- fromJSON(data)
	
	# This happens if the InChI key is not found:
	if(!is.null(r$Fault))
	return(NA)
	
	# Find the IUPAC-Names
	if(!is.null(r$PC_Compounds[[1]]$props)){
		IUPACIndex <- which(unlist(lapply(r$PC_Compounds[[1]]$props, function(i) (i$urn$label == "IUPAC Name"))))
		if(length(IUPACIndex) > 0){
			# Retrieve all IUPAC-Names
			IUPACEntries <- lapply(IUPACIndex, function(x) r$PC_Compounds[[1]]$props[[x]])
			if(!is.null(IUPACEntries)){
				# Is there a preferred IUPAC-Name? If yes, retrieve that
				PrefIUPAC <- which(unlist(lapply(IUPACEntries, function(x) x$urn$name == "Preferred")))
			}	else{return(NA)}
		}	else{return(NA)}
	}	else{return(NA)}
	

	if(length(PrefIUPAC) == 1){
		return(IUPACEntries[[PrefIUPAC]]$value$sval)
	} else{
		# Else it doesn't matter which
		return(IUPACEntries[[1]]$value$sval)
	}
} 

getPcInchiKey <- function(query, from = "smiles"){
	# Get the JSON-Data from Pubchem
	baseURL <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"
	url <- paste(baseURL, from, query, "record", "json", sep="/")
	errorvar <- 0
	currEnvir <- environment()
	
	tryCatch(
		{#data <- getURL(URLencode(url),timeout=8)
		  res <- GET(URLencode(url))
		  data <- httr::content(res, type="text", encoding="UTF-8")
		},
		error=function(e){
		currEnvir$errorvar <- 1
	})
	
	if(errorvar){
		return(NA)
	}
	
	r <- fromJSON(data)
	
	# This happens if the InChI key is not found:
	if(!is.null(r$Fault))
	return(NA)
	
	# Find the entries which contain Chebi-links
	if(!is.null(r$PC_Compounds[[1]]$props)){
		INKEYindex <- which(sapply(r$PC_Compounds[[1]]$props, function(x) x$urn$label) == "InChIKey")
		if(length(INKEYindex) > 0){
			return(r$PC_Compounds[[1]]$props[[INKEYindex]]$value$sval)
		}	else{return(NA)}
	}	else{return(NA)}

	
}

getPcSDF <- function(query, from = "smiles"){
	baseURL <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"
	url <- paste(baseURL, from, query, "sdf", sep="/")
	
	errorvar <- 0
	currEnvir <- environment()
	
	tryCatch(
		{#data <- getURL(URLencode(url),timeout=8)
		  res <- GET(URLencode(url))
		  data <- httr::content(res, type="text", encoding="UTF-8")
		},
		error=function(e){
		currEnvir$errorvar <- 1
	})
	
	if(errorvar){
		return(NA)
	}
	
	molEnd <- regexpr(data,pattern="M  END",fixed=TRUE)+5
	data <- c(strsplit(substring(data,1,molEnd),"\n")[[1]],"$$$$")
	return(data)
}


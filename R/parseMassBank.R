#' MassBank-record Parser
#' 
#' Can parse MassBank-records(only V2)
#'
#' @aliases parseMassBank
#' @usage parseMassBank(Files)
#' @param Files array of character-strings
#' Paths to the plaintext-records that should be read
#' @return The \code{mbWorkspace} that the plaintext-record creates.
#' All parsed information will be stored in the 'compiled_ok' slot.
#' @seealso \code{\link{validate}}
#' @author Erik Mueller
#' @examples \dontrun{
#'		paths <- c("filepath_to_records/RC000001.txt",
#'		           "filepath_to_records/RC000002.txt")
#' 		mb <- parseMassBank(paths)
#' }
#' @export
parseMassBank <- function(Files){
	mb <- new("mbWorkspace")
	mb@compiled_ok <- list()
	for (i in seq_along(Files)) {
		fileConnection <- file(Files[i])
		record <- readLines(fileConnection)
		close(fileConnection)
		mb@compiled_ok[[i]] <- list()
		mb@compiled_ok[[i]][['ACCESSION']] <- substring(grep('ACCESSION:',record, value = TRUE, fixed = TRUE),12)
		mb@compiled_ok[[i]][['RECORD_TITLE']] <- substring(grep('RECORD_TITLE:',record, value = TRUE),12)
		mb@compiled_ok[[i]][['DATE']] <- format(as.Date(substring(grep('DATE:',record, value = TRUE, fixed = TRUE),7), format = "%Y.%m.%d"), "%Y.%m.%d")
		mb@compiled_ok[[i]][['AUTHORS']] <- substring(grep('AUTHORS:',record, value = TRUE, fixed = TRUE),10)
		mb@compiled_ok[[i]][['LICENSE']] <- substring(grep('LICENSE:',record, value = TRUE, fixed = TRUE),10)
		mb@compiled_ok[[i]][['COPYRIGHT']] <- substring(grep('COPYRIGHT:',record, value = TRUE, fixed = TRUE),12)
		##publication <- substring(grep('PUBLICATION:',record, fixed = TRUE),14)
		##if(length(publication) > 0){
		#mb@compiled_ok[[i]][['PUBLICATION']] <- publication
		##}
	
		##The list of comments is handled differently
		##in RMassBank, but the flattening should work anyway, if I'm correct(RMassBank uses internal values for comments)
		commentlist <- list()
		commentlist <- as.list(substring(grep('COMMENT:',record, value = TRUE, fixed = TRUE),10))
		mb@compiled_ok[[i]][['COMMENT']] <- list()
		mb@compiled_ok[[i]][['COMMENT']] <- commentlist
		chnames <- list()
		chnames <- as.list(substring(grep('CH$NAME:',record, value = TRUE, fixed = TRUE),10))
		mb@compiled_ok[[i]][['CH$NAME']] <- chnames
		mb@compiled_ok[[i]][['CH$COMPOUND_CLASS']] <- substring(grep('CH$COMPOUND_CLASS:',record, value = TRUE, fixed = TRUE),20)
		mb@compiled_ok[[i]][['CH$FORMULA']] <- substring(grep('CH$FORMULA:',record, value = TRUE, fixed = TRUE),13)
		mb@compiled_ok[[i]][['CH$EXACT_MASS']] <- as.numeric(substring(grep('CH$EXACT_MASS:',record, value = TRUE, fixed = TRUE),16))
		mb@compiled_ok[[i]][['CH$SMILES']] <- substring(grep('CH$SMILES:',record, value = TRUE, fixed = TRUE),12)
		mb@compiled_ok[[i]][['CH$IUPAC']] <- substring(grep('CH$IUPAC:',record, value = TRUE, fixed = TRUE),11)
	
		##Again: Flattening this should be no Problem, although the structure is different -
		##RMassBank names every type of link, but this isn't necessary here since we're only
		##reading, not creating. If that's a problem, I'll change it.
		links <- list()
		links <- as.list(substring(grep('CH$LINK:',record, value = TRUE, fixed = TRUE),10))
		mb@compiled_ok[[i]][['CH$LINK']] <- links
	
		##SP$ will be included later since it's kind of rarely used
	
		mb@compiled_ok[[i]][['AC$INSTRUMENT']] <- substring(grep('AC$INSTRUMENT:',record, value = TRUE, fixed = TRUE),16)
		mb@compiled_ok[[i]][['AC$INSTRUMENT_TYPE']] <- substring(grep('AC$INSTRUMENT_TYPE:',record, value = TRUE, fixed = TRUE),21)
		##Get the Subvalues just like in RMassBank
	
		##RECORD VERSION SPECIFIC READING INCLUDED
		##This could convert Version 1 -> Version 2 if used right,
		##Although I have no idea how well it'd do that
		##I'll have to find the old specifications to do this right, until then it should only kind of work
		##well enough to do some tests
		Version <- 2
		ac_ms <- list()
		ac_ms[['MS_TYPE']] <- substring(grep('AC$MASS_SPECTROMETRY: MS_TYPE',record, value = TRUE, fixed = TRUE),31)
		if(length(ac_ms[['MS_TYPE']]) == 0){
			ac_ms[['MS_TYPE']] <- substring(grep('AC$ANALYTICAL_CONDITION: MS_TYPE',record, value = TRUE, fixed = TRUE),34)
			Version <- 1
		}
		if(Version == 1){
			##This not a real tag anymore(according to the specifications) but RMassBank still writes it...?
			##I'll include it for the case that I'm reading V1-records
		    ac_ms[['ION_MODE']] <- substring(grep('AC$ANALYTICAL_CONDITION: MODE',record, value = TRUE, fixed = TRUE),31)
		    ac_ms[['IONIZATION']] <- substring(grep('AC$MASS_SPECTROMETRY: IONIZATION',record, value = TRUE, fixed = TRUE),34)
			
		} else{
				ac_ms[['ION_MODE']] <- substring(grep('AC$MASS_SPECTROMETRY: ION_MODE',record, value = TRUE, fixed = TRUE),32)
				
				##Some of the following are part of the (optional) specification, but NOT in RMassBank(!)
				##This is just for the sake of completeness
				ac_ms[['COLLISION_ENERGY']] <- substring(grep('AC$MASS_SPECTROMETRY: COLLISION_ENERGY',record, value = TRUE, fixed = TRUE),40)
				ac_ms[['COLLISION_GAS']] <- substring(grep('AC$MASS_SPECTROMETRY: COLLISION_GAS',record, value = TRUE, fixed = TRUE),37)
				ac_ms[['DATE']] <- substring(grep('AC$MASS_SPECTROMETRY: DATE',record, value = TRUE, fixed = TRUE),28)
				ac_ms[['DESOLVATION_GAS_FLOW']] <- substring(grep('AC$MASS_SPECTROMETRY: DESOLVATION_GAS_FLOW',record, value = TRUE, fixed = TRUE),44)
				ac_ms[['DESOLVATION_TEMPERATURE']] <- substring(grep('AC$MASS_SPECTROMETRY: DESOLVATION_TEMPERATURE',record, value = TRUE, fixed = TRUE),47)
				ac_ms[['IONIZATION']] <- substring(grep('AC$MASS_SPECTROMETRY: IONIZATION',record, value = TRUE, fixed = TRUE),34)
				ac_ms[['IONIZATION_ENERGY']] <- substring(grep('AC$MASS_SPECTROMETRY: IONIZATION_ENERGY',record, value = TRUE, fixed = TRUE),41)
				ac_ms[['LASER']] <- substring(grep('AC$MASS_SPECTROMETRY: LASER',record, value = TRUE, fixed = TRUE),29)
				ac_ms[['MATRIX']] <- substring(grep('AC$MASS_SPECTROMETRY: MATRIX',record, value = TRUE, fixed = TRUE),30)
				ac_ms[['MASS_ACCURACY']] <- substring(grep('AC$MASS_SPECTROMETRY: MASS_ACCURACY',record, value = TRUE, fixed = TRUE),37)
				ac_ms[['REAGENT_GAS']] <- substring(grep('AC$MASS_SPECTROMETRY: REAGENT_GAS',record, value = TRUE, fixed = TRUE),35)
				ac_ms[['SCANNING']] <- substring(grep('AC$MASS_SPECTROMETRY: SCANNING',record, value = TRUE, fixed = TRUE),32)
				
				##These are in RMassBank, but not part of the specification?
				##I think I'm misreading something...
				#ac_ms[['FRAGMENTATION_MODE']] <- msmsdata$info$mode
				#ac_ms[['PRECURSOR_TYPE']] <- precursor_types[spec$mode]
				#ac_ms[['RESOLUTION']] <- msmsdata$info$res
				
				ac_lc <- list();
				ac_lc[['CAPILLARY_VOLTAGE']] <- substring(grep('AC$CHROMATOGRAPHY: CAPILLARY_VOLTAGE',record, value = TRUE, fixed = TRUE),36)
				ac_lc[['COLUMN_NAME']] <- substring(grep('AC$CHROMATOGRAPHY: COLUMN_NAME',record, value = TRUE, fixed = TRUE),32)
				ac_lc[['COLUMN_TEMPERATURE']] <- substring(grep('AC$CHROMATOGRAPHY: COLUMN_TEMPERATURE',record, value = TRUE, fixed = TRUE),39)
				ac_lc[['FLOW_GRADIENT']] <- substring(grep('AC$CHROMATOGRAPHY: FLOW_GRADIENT',record, value = TRUE, fixed = TRUE),34)
				ac_lc[['FLOW_RATE']] <- substring(grep('AC$CHROMATOGRAPHY: FLOW_RATE',record, value = TRUE, fixed = TRUE),30)
				ac_lc[['RETENTION_TIME']] <- substring(grep('AC$CHROMATOGRAPHY: RETENTION_TIME',record, value = TRUE, fixed = TRUE),35)
				ac_lc[['SOLVENT A']] <- substring(grep('AC$CHROMATOGRAPHY: SOLVENT A',record, value = TRUE, fixed = TRUE),30)
				ac_lc[['SOLVENT B']] <- substring(grep('AC$CHROMATOGRAPHY: SOLVENT B',record, value = TRUE, fixed = TRUE),30)
				
				ms_fi <- list()
				ms_fi[['BASE_PEAK']] <- as.double(substring(grep('MS$FOCUSED_ION: BASE_PEAK',record, value = TRUE, fixed = TRUE),27))
				ms_fi[['PRECURSOR_M/Z']] <- substring(grep('MS$FOCUSED_ION: PRECURSOR_M/Z',record, value = TRUE, fixed = TRUE),31)
				ms_fi[['PRECURSOR_TYPE']] <- substring(grep('MS$FOCUSED_ION: PRECURSOR_TYPE',record, value = TRUE, fixed = TRUE),32)
				
				if(ac_ms[['MS_TYPE']] == 'MS2'){
					ms_fi[['PRECURSOR_M/Z']] <- as.double(ms_fi[['PRECURSOR_M/Z']])
				}
		}
		namesAcms <- names(ac_ms)
		namesAclc <- names(ac_lc)
		namesMsfi <- names(ms_fi)
		for(k in 1:length(ac_ms)){
			if(length(ac_ms[[namesAcms[k]]]) == 0){
				ac_ms[[namesAcms[k]]] <- NA
			}
		}
		for(k in 1:length(ac_lc)){
			if(length(ac_lc[[namesAclc[k]]]) == 0){
				ac_lc[[namesAclc[k]]] <- NA
			}
		}
		for(k in 1:length(ms_fi)){
			if(length(ms_fi[[namesMsfi[k]]]) == 0){
				ms_fi[[namesMsfi[k]]] <- NA
			}
		}
		mb@compiled_ok[[i]][['AC$MASS_SPECTROMETRY']] <- list()
		mb@compiled_ok[[i]][['AC$MASS_SPECTROMETRY']] <- ac_ms
		mb@compiled_ok[[i]][['AC$CHROMATOGRAPHY']] <- list()
		mb@compiled_ok[[i]][['AC$CHROMATOGRAPHY']] <- ac_lc
		mb@compiled_ok[[i]][['MS$FOCUSED_ION']] <- list()
		mb@compiled_ok[[i]][['MS$FOCUSED_ION']] <- ms_fi
		
		##Can currently only read annotations of the type "m/z num {formula mass error(ppm)}"
		##and'll only read it properly if there is only one annotation
		##the strange conversion of the data.frames is there so RMassBank can actually write it again
		PKannotationStart <- grep('PK$ANNOTATION:',record, fixed = TRUE) + 1
		numpeak <- grep('PK$NUM_PEAK:',record, fixed = TRUE)
			
			if(length(PKannotationStart) > 0 && ReadAnnotation == TRUE){
				if(PKannotationStart < numpeak){
					splitted <- strsplit(record[PKannotationStart:(numpeak-1)]," ")
					PKannotation <- matrix(nrow = numpeak - PKannotationStart, ncol = 5)
					for(k in 1:length(splitted)){
						splitted[[k]] <- splitted[[k]][which(splitted[[k]] != "")]
						PKannotation[k,] <- splitted[[k]]
					}
					PKannotation <- as.data.frame(PKannotation, stringsAsFactors = FALSE)
					PKannotation[] <- lapply(PKannotation, type.convert)
					colnames(PKannotation) <- c("m/z", "num", "{formula", "mass", "error(ppm)}")
					PKannotation$"{formula" <- as.character(PKannotation$"{formula")
				}
				mb@compiled_ok[[i]][['PK$ANNOTATION']] <- PKannotation
			}
		##Extract the peaks and write the data into a data.frame
		PKStart <- grep('PK$PEAK:',record, fixed = TRUE) + 1
		endslash <- tail(grep('//',record, fixed = TRUE),1)
        if(PKStart < endslash){
            splitted <- strsplit(record[PKStart:(endslash-1)]," ")
            PKPeak <- matrix(nrow = endslash - PKStart, ncol = 3)
            for(k in 1:length(splitted)){
                splitted[[k]] <- splitted[[k]][which(splitted[[k]] != "")]
                PKPeak[k,] <- splitted[[k]]
            }
            PKPeak <- as.data.frame(PKPeak, stringsAsFactors = FALSE)
            PKPeak[] <- lapply(PKPeak, type.convert)
            colnames(PKPeak) <- c("m/z", "int", "rel.int.")
        }
	
		mb@compiled_ok[[i]][['PK$PEAK']] <- PKPeak	
	
		namesComp <- names(mb@compiled_ok[[i]])
		for(k in 1:length(mb@compiled_ok[[i]])){
			if(length(mb@compiled_ok[[i]][[namesComp[k]]]) == 0){
				mb@compiled_ok[[i]][[namesComp[k]]] <- NA
			}
		}
		print(paste("Read",Files[i]))
		flush.console()
	}
	return(mb)
}

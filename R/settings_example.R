.checkMbSettings <- function()
{
	o <- getOption("RMassBank", NULL)
	if(is.null(o)){
		stop("Please load your settings before using the RMassBank workflow.")
	}

}


#' RMassBank settings
#'
#' Describes all settings for the RMassBank settings file.
#'
#' \itemize{
#' 		\item{\code{deprofile}}{
#'   	Whether and how to deprofile input raw files. Leave the
#' 			setting empty if your raw files are already in "centroid" mode. If your
#' 			input files are in profile mode, you have the choice between algorithms
#' 			\code{\link{deprofile}.spline, deprofile.fwhm, deprofile.localMax}; refer to
#' 			the individual manpages for more information.}
#' 		\item{\code{rtMargin, rtShift}}{
#'   	The allowed retention time deviation relative to the
#' 			values specified in your compound list (see \code{\link{loadList}}), and the systematic
#' 			shift (due to the use of, e.g., pre-columns or other special equipment.}
#' 		\item{\code{babeldir}}{
#' 			Directory to OpenBabel. Required for creating molfiles for MassBank export.
#' 			If no OpenBabel directory is given, RMassBank will attempt to use the CACTUS webservice
#' 			for SDF generation. It is strongly advised to install OpenBabel; the CACTUS structures
#' 			have explicit hydrogen atoms.
#'			The path should point to the directory where babel.exe (or the Linux "babel" equivalent) lies.
#' 			}
#' 		\item{\code{use_version}}{
#'   	Which MassBank record format to use; version 2 is strongly advised,
#' 			version 1 is considered outdated and should be used only if for some reason you are running
#' 			old servers and an upgrade is not feasible.}
#' 		\item{\code{use_rean_peaks}}{
#'   	Whether to include peaks from reanalysis (see
#' 			\code{\link{reanalyzeFailpeaks}}) in the MassBank records. Boolean, TRUE or FALSE.
#' 			}
#' 		\item{\code{annotations}}{
#' 			A list of constant annotations to use in the MassBank records. The entries
#' 			\code{authors, copyright, license, instrument, instrument_type, compound_class}
#' 			correspond to the MassBank entries \code{AUTHORS, COPYRIGHT, PUBLICATION, LICENSE, AC$INSTRUMENT,
#' 			AC$INSTRUMENT_TYPE, CH$COMPOUND_CLASS}. The entry \code{confidence_comment} is added as
#' 			\code{COMMENT: CONFIDENCE} entry.
#'
#' 			The entry \code{internal_id_fieldname} is used to name
#' 			the MassBank entry which will keep a reference to the internal compound ID used in
#' 			the workflow: for \code{internal_id_fieldname = MYID} and e.g. compound 1234, an
#' 			entry will be added	to the MassBank record with
#' 			\code{COMMENT: MYID 1234}. The internal fieldname should not be left empty!
#'
#' 			The entries \code{lc_gradient, lc_flow, lc_solvent_a, lc_solvent_b, lc_column} correspond
#' 			to the MassBank entries \code{AC$CHROMATOGRAPHY: FLOW_GRADIENT, FLOW_RATE,
#' 			SOLVENT A, SOLVENT B, COLUMN_NAME}.
#'
#' 			\code{ms_type, ionization} correspond to \code{AC$MASS_SPECTROMETRY: MS_TYPE, IONIZATION}.
#'
#' 			\code{entry_prefix} is the two-letter prefix used when building MassBank accession codes.
#'
#' 			Entries under \code{ms_dataprocessing} are added as \code{MS$DATA_PROCESSING:} entries,
#' 			in addition to the default \code{WHOLE: RMassBank}.
#' 			}
#'   	\item{\code{annotator}}{
#'     For advanced users: option to select your own custom annotator.
#'     Check \code{\link{annotator.default}} and the source code for details.}
#' 		\item{\code{spectraList}}{
#'   	This setting describes the experimental annotations for the single
#' 			data-dependent scans. For every data-dependent scan event, a \code{spectraList} entry with
#' 			\code{mode, ces, ce, res} denoting collision mode, collision energy in short and verbose
#' 			notation, and FT resolution.}
#' 		\item{\code{accessionNumberShifts}}{
#'   	This denotes the starting points for accession numbers
#' 			for different ion types. For example, \code{pH: 0, mH: 50} means that [M+H]+ spectra will
#' 			start at \code{XX123401} (\code{XX} being the \code{entry_prefix} and \code{1234} the compound
#' 			id) and [M-H]- will start at \code{XX123451}.}
#' 		\item{\code{electronicNoise, electronicNoiseWidth}}{
#'   	Known electronic noise peaks and the window
#' 			to be used by \code{\link{cleanElnoise}}}
#' 		\item{\code{recalibrateBy}}{
#'   	\code{dppm} or \code{dmz} to recalibrate either by delta ppm or by
#' 			delta mz.}
#' 		\item{\code{recalibrateMS1}}{
#'   	\code{common} or \code{separate} to recalibrate MS1 data points together
#' 			or separately from MS2 data points.}
#' 		\item{\code{recalibrator: MS1, MS2}}{
#'   	The functions to use for recalibration of MS1 and MS2 data points.
#' 			Note that the \code{MS1} setting is only meaningful if \code{recalibrateMS1: separate}, otherwise
#' 			the \code{MS2} setting is used for a common recalibration curve. See \code{\link{recalibrate.loess}}
#' 			for details.}
#'   	\item{\code{multiplicityFilter}}{
#'     Define the multiplicity filtering level. Default is 2, a value of 1
#'     is off (no filtering) and >2 is harsher filtering.}
#'     \item{\code{titleFormat}}{
#'     The title of MassBank records is a mini-summary
#'     of the record, for example "Dinotefuran; LC-ESI-QFT; MS2; CE: 35%; R=35000; [M+H]+".
#'     By default, the first compound name \code{CH$NAME}, instrument type
#'     \code{AC$INSTRUMENT_TYPE}, MS/MS type \code{AC$MASS_SPECTROMETRY: MS_TYPE},
#'     collision energy \code{RECORD_TITLE_CE}, resolution \code{AC$MASS_SPECTROMETRY: RESOLUTION}
#'     and precursor \code{MS$FOCUSED_ION: PRECURSOR_TYPE} are used. If alternative
#'     information is relevant to differentiate acquired spectra, the title should be adjusted.
#'     For example, many TOFs do not have a resolution setting.
#'     See MassBank documentation for more.}
#'   	\item{\code{filterSettings}}{
#' 			A list of settings that affect the MS/MS processing. The entries
#' 			\code{ppmHighMass, ppmLowMass, massRangeDivision} set values for
#'   		pre-processing, prior to recalibration. \code{ppmHighMass} defines the
#'     	ppm error for the high mass range (default 10 ppm for Orbitraps),
#'       \code{ppmLowMass} is the error for the low mass range (default 15 ppm
#'       for Orbitraps) and \code{massRangeDivision} is the m/z value defining
#'       the split between the high and low mass range (default m/z = 120).
#'
#' 			The entry \code{ppmFine} defines the ppm cut-off post recalibration.
#'   		The default value of 5 ppm is recommended for Orbitraps. For other
#'     	instruments this can be interpreted from the recalibration plot.
#'      All ppm limits are one-sided (e.g. this includes values to +5 ppm or -5 ppm
#'      deviation from the exact mass).
#'
#' 			The entries \code{prelimCut, prelimCutRatio} define the intensity cut-off and
#'   		cut-off ratio (in % of the most intense peak) for pre-processing. This affects
#'     	the peak selection for the recalibration only. Careful: the default value
#'       1e4 for Orbitrap LTQ positive mode could remove all peaks for TOF data
#'       and will remove too many peaks for Orbitrap LTQ negative mode spectra!
#'
#' 			The entry \code{specOKLimit} defines the intensity limit to include MS/MS spectra.
#'   		MS/MS spectra must have at least one peak above this limit to proceed through
#'     	the workflow.
#'
#' 			\code{dbeMinLimit} defines the minimum allowable ring and double bond equivalents (DBE)
#'   		allowed for assigned formulas. This assumes maximum valuences for elements with
#'     	multiple valence states. The default is -0.5 (accounting for fragments being ions).
#'
#' 			The entries \code{satelliteMzLimit, satelliteIntLimit} define the cut-off m/z and
#'   		intensity values for satellite peak removal (an artefact of Fourier Transform
#'     	processing). All peaks within the m/z limit (default 0.5) and intensity ratio
#'       (default 0.05 or 5 %) of the respective peak will be removed. Applicable to
#'       Fourier Transform instruments only (e.g. Orbitrap).
#' 			}
#'     \item{\code{filterSettings}}{
#' 			Parameters for adjusting the raw data retrieval.
#'   		The entry \code{ppmFine} defines the ppm error to look for the precursor in
#'     	the MS1 (parent) spectrum. Default is 10 ppm for Orbitrap.
#'
#' 			\code{mzCoarse} defines the error to search for the precursor specification
#'   		in the MS2 spectrum. This is often only saved to 2 decimal places and thus
#'     	can be quite inaccurate. The accuracy also depends on the isolation window used.
#'       The default settings (for e.g. Orbitrap) is 0.5 (Da, or Th for m/z).
#'
#' 			The entry \code{fillPrecursorScan} is largely untested. The default value
#'   		(FALSE) assumes all necessary precursor information is available in the mzML file.
#'     	A setting ot TRUE tries to fill in the precursor data scan number if it is missing.
#'       Only tested on one case study so far - feedback welcome!
#' 			}
#' }
#'
#'
#' @author Michael Stravs, Emma Schymanski
#' @seealso \code{\link{loadRmbSettings}}
#' @rdname RmbSettings
#' @name RmbSettings
NULL

.settingsList <- list(
  # Deprofile input data?
  # NA if input data is already in "centroid" mode,
  # "deprofile.fwhm" or "deprofile.localMax" to convert the input data with the
  # corresponding algorithm. See ?deprofile
  deprofile = NA,
  # Deviation (in minutes) allowed the for retention time
  rtMargin = 0.4,
  # Systematic RT shift
  rtShift = -0.3,
  # Directory to OpenBabel. Required for MassBank export
  babeldir = NA,
  # Which MassBank format should be used? Version 2 is advised.
  use_version = 2,
  # Include reanalyzed peaks?
  use_rean_peaks = TRUE,
  # annotate the spectra files with (putative) molecular formulas for fragments?
  add_annotation = TRUE,
  # Annotations for the spectrum:
  annotations = list(
    authors = "Nomen Nescio, The Unseen University",
    copyright = "Copyright (C) XXX",
    publication = "",
    license = "CC BY",
    instrument = "LTQ Orbitrap XL Thermo Scientific",
    instrument_type = "LC-ESI-ITFT",
    confidence_comment = "standard compound",
    compound_class = "N/A; Environmental Standard",
    internal_id_fieldname = "INTERNAL_ID",
	#
	# HPLC annotations:
	#
    # example: lc_gradient = "90/10 at 0 min, 50/50 at 4 min, 5/95 at 17 min, 5/95 at 25 min, 90/10 at 25.1 min, 90/10 at 30 min""
    lc_gradient = "",
    # example: lc_flow = "200 uL/min",
    lc_flow = "",
    # example: lc_solvent_a = 'water with 0.1% formic acid',
    lc_solvent_a = '',
    lc_solvent_b = '',
    # example: lc_column "XBridge C18 3.5um, 2.1x50mm, Waters",
    lc_column = "",
	#
	# Prefix for MassBank accession IDs
	#
    entry_prefix = "XX",
	  contributor_prefix = "CONTRIBUTOR",
    ms_type = "MS2",
    ionization = "ESI",
    ms_dataprocessing = list(
			"RECALIBRATE" = "loess on assigned fragments and MS1"
			)
    ),
  include_sp_tags = FALSE,
  # List of data-dependent scans in their order (relative to the parent scan)
	# list(mode, ces, ce, res):
	# mode: fragmentation mode
	# ces: "short" format collision energy (for record title)
	# ce: "long" format collision energy (for annotation field)
	# res: FT resolution
  spectraList = list(
    list(mode="CID", ces = "35%", ce = "35 % (nominal)", res = 7500),
    list(mode="HCD", ces = "15%", ce = "15 % (nominal)", res = 7500),
    list(mode="HCD", ces = "30%", ce = "30 % (nominal)", res = 7500),
    list(mode="HCD", ces = "45%",ce = "45 % (nominal)", res = 7500),
    list(mode="HCD", ces = "60%",ce = "60 % (nominal)", res = 7500),
    list(mode="HCD", ces = "75%",ce = "75 % (nominal)", res = 7500),
    list(mode="HCD", ces = "90%",ce = "90 % (nominal)", res = 7500),
    list(mode="HCD", ces = "15%",ce = "15 % (nominal)", res = 15000),
    list(mode="HCD", ces = "30%",ce = "30 % (nominal)", res = 15000),
    list(mode="HCD", ces = "45%",ce = "45 % (nominal)", res = 15000),
    list(mode="HCD", ces = "60%",ce = "60 % (nominal)", res = 15000),
    list(mode="HCD", ces = "75%", ce = "75 % (nominal)", res = 15000),
    list(mode="HCD", ces = "90%", ce = "90 % (nominal)", res = 15000),
    list(mode="CID", ces = "35%", ce = "35 % (nominal)", res = 15000)
  ),
  accessionNumberShifts = list(
    "pH" = 0, # [M+H]+: Accession numbers 1-14
    "pM" = 16, # [M]+: 17-30
    "pNa" = 32, # [M+Na]+: 33-46
    "mH" = 50, # [M-H]-: 51-64
    "mFA" = 66, # [M+FA]-: 67-80
	"mM" = 80 # [M]-: 81-94
    ),
  accessionBuilderType = NULL,
  accessionBuilder = "MSBNK-{contributor_prefix}-{entry_prefix}{compound_id(4)}{scan_id(2)}",
  # Validate accession? Set to FALSE to bypass accession validation
  accessionValidate = TRUE,
  # Known electronic noise peaks in the Orbitrap data
  electronicNoise = c(189.825, 201.725,196.875),
  # Exclusion width of electronic noise peaks (from unmatched peaks, prior to
  # reanalysis)
  electronicNoiseWidth = 0.3,
  # recalibration settings:
  # recalibrate by: dppm or dmz
  recalibrateBy = "dppm",
  # recalibrate MS1:
  # separately ("separate")
  # with common curve ("common")
  # do not recalibrate ("none")
  recalibrateMS1 = "common",
  # Custom recalibration function: You can overwrite the recal function by
  # making any function which takes rcdata$recalfield ~ rcdata$mzFound.
  # The settings define which recal function is used
  recalibrator = list(
	MS1 = "recalibrate.loess",
	MS2 = "recalibrate.loess"),
  # Window width to look for MS1 peaks to recalibrate (in ppm)
  recalibrateMS1Window= 15,

  # Define the multiplicity filtering level
  # Default is 2 (peak occurs at least twice)
  # Set this to 1 if you want to turn this option off.
  # Set this to anything > 2 if you want harder filtering
  multiplicityFilter = 2,
	# Define the title format.
	# You can use all entries from MassBank records as tokens
	# plus the additional token RECORD_TITLE_CE, which is a shortened
	# version of the collision energy specifically for use in the title.
	# Every line is one entry and must have one token in curly brackets
	# e.g. {CH$NAME} or {AC$MASS_SPECTROMETRY: MS_TYPE} plus optionally
	# additional text in front or behind e.g.
	# R={AC$MASS_SPECTROMETRY: RESOLUTION}
	# If this is not specified, it defaults to a title of the format
	# "Dinotefuran; LC-ESI-QFT; MS2; CE: 35%; R=35000; [M+H]+"
  titleFormat = c(
		  "{CH$NAME}",
		  "{AC$INSTRUMENT_TYPE}",
		  "{AC$MASS_SPECTROMETRY: MS_TYPE}",
		  "CE: {RECORD_TITLE_CE}",
		  "R={AC$MASS_SPECTROMETRY: RESOLUTION}",
		  "{MS$FOCUSED_ION: PRECURSOR_TYPE}"
  ),
  # Define filter settings.
  # For Orbitrap, settings of 15 ppm in low mass range, 10 ppm in high
  # mass range, m/z = 120 as mass range division and 5 ppm for recalibrated
  # data overall are recommended.
  filterSettings = list(
		  ppmHighMass = 10,
  		  ppmLowMass = 15,
		  massRangeDivision= 120,
		  ppmFine= 5,
		  prelimCut= 1e4,
		  prelimCutRatio= 0,
		  fineCut= 0,
		  fineCutRatio= 0,
		  specOkLimit= 1e4,
		  dbeMinLimit= -0.5,
		  satelliteMzLimit= 0.5,
		  satelliteIntLimit= 0.05
		  ),
 # Define raw MS retrieval settings.
 # fillPrecursorScan is FALSE for "good" mzML files which have all the info needed.
 # However, for example AB Sciex files will have missing precursor scan information,
 # in which case fillPrecursorScan = TRUE is needed. Try it out.
 findMsMsRawSettings = list(
		ppmFine= 10,
		mzCoarse= 0.5,
		fillPrecursorScan= FALSE),

 # Select how to treat unknown compound masses:
 # "charged" (the default, also if no option set) treats unknown (level 5) compound masses as the m/z,
 # "neutral" treats unknown (level 5) compound masses as the neutral mass and applies [M+H]+ and [M-H]- calculations accordingly.
 unknownMass = "charged",

 # Add the CCTE api key to retrieve information from https://api-ccte.epa.gov/docs
 # Be aware, this is confidential information, so do not share with unauthorized
 # persons
 ccte_api_key = NULL,

 # Add the RSC api key to retrieve information from https://developer.rsc.org/api-reference
 # Be aware, this is confidential information, so do not share with unauthorized
 # persons
 rcs_api_key = NULL
)

# Writes a file with sample settings which the user can adjust with his values.
#' @export
RmbSettingsTemplate <- function(target)
{
  blub <- file.copy(from=system.file("RMB_options.ini", package="RMassBank"), to=target)
}

# Loads settings from a file or from an object.
#' @export
loadRmbSettings <- function(file_or_list)
{
  # If the object exists in R, it is assumed to be the list itself
  # Otherwise, it's assumed to be a file name and loaded.
  # It will be either an INI file in YAML format or an R file to be directly sourced.
  if(is.list(file_or_list))
    options(RMassBank = file_or_list)
  else if(exists(file_or_list, inherits=TRUE))
    options(RMassBank = get(file_or_list))
  else if(file.exists(file_or_list))
  {
	# Check if the file has an INI extension:
	# If yes, load it with YAML
	isIni <- grepl("\\.[iI][nN][iI]$", file_or_list, perl = T)
	isIni <- isIni || grepl("\\.[yY][mM][lL]$", file_or_list, perl = T)
	isR <- grepl("\\.[rR]$", file_or_list, perl = T)
	if(isIni)
	{
		o <- yaml.load_file(file_or_list)
		# Fix the YAML file to suit our needs
		if(is.null(o$deprofile))
			o$deprofile <- NA
		if(is.null(o$babeldir)){
			o$babeldir <- NA
		} else{
			##Check if babeldir exists
			babelcheck <- gsub('\"','',o$babeldir)
			if(substring(babelcheck, nchar(babelcheck)) == "\\"){
				babelexists <- file.exists(substring(babelcheck, 1, nchar(babelcheck)-1))
			} else{
				babelexists <- file.exists(babelcheck)
			}

			if(!babelexists){
				stop("The babeldir does not exist. Please check the babeldir in the settings and adjust it accordingly.")
			}
		}

		for(name in names(o$annotations))
		{
			if(is.null(o$annotations[[name]]))
				o$annotations[[name]] <- ""
		}
		if (!is.null(o$logging_file)) {
			appender_obj <- logger::appender_file(o$logging_file)
			# This implicitly creates a new namespace in the
			# logger package, that is used to treat calls from
			# RMassBank differently
			log_appender(appender_obj, namespace='RMassBank')
		}
		options(RMassBank = o)
	}
	else if (isR)
	{
		ov <- source(file_or_list)
    o <- ov$value
		options(RMassBank = o)
	}
	else
		stop("Options format not recognized. Use YAML (.ini, .yml) or R file (.R) format.")

  }
  else
    stop("The file path supplied for the options does not exist.")

  # Settings are loaded, now check if they are up to date
  o <- getOption("RMassBank")
  curr <- names(.settingsList)
  problem <- length(setdiff(curr, names(o))) > 0
  # Hesch es problem? He?
  if(problem)
  {
    warning("Your settings are outdated. Missing will be replaced by default values.")
    o <- updateSettings(o)
    options(RMassBank = o)
  }
}

#' @export
loadRmbSettingsFromEnv <- function(env = .GlobalEnv)
{
	loadRmbSettings(env$RmbSettings)
}

#' RMassBank settings
#'
#' Load, set and reset settings for RMassBank.
#'
#' \code{RmbSettingsTemplate} creates a template file in which you can adjust the
#' settings as you like. Before using RMassBank, you must then load the
#' settings file using \code{loadRmbSettings}. \code{RmbDefaultSettings} loads
#' the default settings. \code{loadRmbSettingsFromEnv} loads the settings
#' stored in env$RmbSettings, which is useful when reloading archives with
#' saved settings inside.
#'
#' Note: no settings are loaded upon loading MassBank!
#' This is intended, so that one never forgets to load the correct settings.
#'
#' The settings are described in \code{\link{RmbSettings}}.
#'
#' @aliases loadRmbSettings RmbDefaultSettings RmbSettingsTemplate loadRmbSettingsFromEnv
#' @usage loadRmbSettings(file_or_list)
#'
#' loadRmbSettingsFromEnv(env = .GlobalEnv)
#'
#' RmbDefaultSettings()
#'
#' RmbSettingsTemplate(target)
#' @param file_or_list The file (YML or R format) or R \code{list} with the settings to load.
#' @param target The path where the template setting file should be stored.
#' @param env The environment to load the settings from.
#' @return None.
#' @note \bold{The default settings will not work for you unless you have, by
#' chance, installed OpenBabel into the same directory as I have!}
#' @author Michael Stravs
#' @seealso \code{\link{RmbSettings}}
#' @examples
#'
#'  # Create a standard settings file and load it (unedited)
#' 	RmbSettingsTemplate("mysettings.ini")
#'  loadRmbSettings("mysettings.ini")
#'  unlink("mysettings.ini")
#'
#' @export
RmbDefaultSettings <- function()
{
  options("RMassBank" = .settingsList)
}

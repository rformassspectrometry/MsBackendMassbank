##' @param f `character(1)` with the path to an MassBank file.
##'
##' @param msLevel `numeric(1)` with the MS level. Default is 2.
##'
##' @param metaDataBlocks `data.frame` data frame indicating which metadata shall
##'     be read
##'
##' @param ... Additional parameters, currently ignored.
##'
##' @importFrom S4Vectors DataFrame cbind.DataFrame
##'
##' @importFrom IRanges NumericList
##'
##' @author Michael Witting
##'
##' @noRd
.read_massbank <- function(f, msLevel = 2L,
                           metaBlocks = metaDataBlocks(),
                           ...) {

  if (length(f) != 1L)
    stop("Please provide a single mgf file.")

  mb <- scan(file = f, what = "",
              sep = "\n", quote = "",
              allowEscapes = FALSE,
              quiet = TRUE)

  begin <- grep("ACCESSION:", mb) + 1L
  end <- grep("^//$", mb)

  n <- length(begin)
  spec <- vector("list", length = n)
  ac <- vector("list", length = n)
  ch <- vector("list", length = n)
  sp <- vector("list", length = n)
  ms <- vector("list", length = n)
  record <- vector("list", length = n)
  pk <- vector("list", length = n)
  comment <- vector("list", length = n)

  for (i in seq(along = spec)) {
    spec[[i]] <- .extract_mb_spectrum(mb[begin[i]:end[i]])

    if(metaBlocks$read[which(metaBlocks$metadata == "ac")])
      ac[[i]] <- .extract_mb_ac(mb[begin[i]:end[i]])

    if(metaBlocks$read[which(metaBlocks$metadata == "ch")])
      ch[[i]] <- .extract_mb_ch(mb[begin[i]:end[i]])

    if(metaBlocks$read[which(metaBlocks$metadata == "sp")])
      sp[[i]] <- .extract_mb_sp(mb[begin[i]:end[i]])

    if(metaBlocks$read[which(metaBlocks$metadata == "ms")])
      ms[[i]] <- .extract_mb_ms(mb[begin[i]:end[i]])

    if(metaBlocks$read[which(metaBlocks$metadata == "record")])
      record[[i]] <- .extract_mb_record(mb[begin[i]:end[i]])

    if(metaBlocks$read[which(metaBlocks$metadata == "pk")])
      pk[[i]] <- .extract_mb_pk(mb[begin[i]:end[i]])

    if(metaBlocks$read[which(metaBlocks$metadata == "comment")])
      comment[[i]] <- .extract_mb_comment(mb[begin[i]:end[i]])
  }

  res <- DataFrame(do.call(rbind, spec))
  res_ac <- DataFrame(do.call(rbind, ac))
  res_ch <- DataFrame(do.call(rbind, ch))
  res_sp <- DataFrame(do.call(rbind, sp))
  res_ms <- DataFrame(do.call(rbind, ms))
  res_record <- DataFrame(do.call(rbind, record))
  res_pk <- DataFrame(do.call(rbind, pk))
  res_comment <- DataFrame(do.call(rbind, comment))

  # only bind if it contains elements
  if(length(res_ac)) {
    res <- cbind.DataFrame(res, res_ac)
  }

  if(length(res_ch)) {
    res <- cbind.DataFrame(res, res_ch)
  }

  if(length(res_sp)) {
    res <- cbind.DataFrame(res, res_sp)
  }

  if(length(res_ms)) {
    res <- cbind.DataFrame(res, res_ms)
  }

  if(length(res_record)) {
    res <- cbind.DataFrame(res, res_record)
  }

  if(length(res_pk)) {
    res <- cbind.DataFrame(res, res_pk)
  }

  if(length(res_comment)) {
    res <- cbind.DataFrame(res, res_comment)
  }

  for (i in seq_along(res)) {
    if (all(lengths(res[[i]]) == 1))
      res[[i]] <- unlist(res[[i]])
  }

  res$mz <- IRanges::NumericList(res$mz)
  res$intensity <- IRanges::NumericList(res$intensity)
  res$dataOrigin <- f
  res$msLevel <- as.integer(msLevel)
  res

}

##' @param mb `character()` of lines defining a spectrum in mgf
##'     format.
##'
##' @author Michael Witting
##'
##' @noRd
.extract_mb_spectrum <- function(mb) {

  #read the spectrum
  spectrum_start <- grep("PK$PEAK:", mb, fixed = TRUE) + 1
  spectrum_end <- tail(grep("//", mb, fixed = TRUE), 1) - 1

  splitted <- strsplit(mb[spectrum_start:(spectrum_end)]," ")
  spectrum <- matrix(nrow = spectrum_end + 1 - spectrum_start, ncol = 3)

  for(k in 1:length(splitted)){
    splitted[[k]] <- splitted[[k]][which(splitted[[k]] != "")]
    spectrum[k,] <- splitted[[k]]
  }

  # convert to data frame and adjust data type
  spectrum <- as.data.frame(spectrum, stringsAsFactors = FALSE)
  spectrum[] <- lapply(spectrum, type.convert)
  colnames(spectrum) <- c("mz", "intensity", "rel.intensity")

  # isolate spectrum metadata from record
  rtime <- substring(grep("AC$CHROMATOGRAPHY: RETENTION_TIME",
                          mb,
                          value = TRUE,
                          fixed = TRUE),
                     35)
  rtime <- as.numeric(regmatches(rtime,
                                 regexpr("[[:digit:]]+\\.[[:digit:]]+",
                                         rtime)))

  # back up if no RT is supplied
  if(!length(rtime)) rtime <- NA_real_

  precursorMz <- as.numeric(substring(grep("MS$FOCUSED_ION: PRECURSOR_M/Z",
                                           mb,
                                           value = TRUE,
                                           fixed = TRUE),
                                      30))

  precursorIntensity <- as.numeric(substring(grep("MS$FOCUSED_ION: PRECURSOR_INT",
                                                  mb,
                                                  value = TRUE,
                                                  fixed = TRUE),
                                             31))

  title <- substring(grep("RECORD_TITLE:",
                          mb,
                          value = TRUE,
                          fixed = TRUE),
                     15)

  # check data
  if(!length(precursorIntensity))
    precursorIntensity <- 100

  list(rtime = rtime * 60,
       scanIndex = as.integer(1),
       precursorMz = precursorMz,
       precursorIntensity = precursorIntensity,
       precursorCharge = as.integer(0),
       mz = spectrum$mz,
       intensity = spectrum$intensity,
       title = title)

}

##' @param mb `character()` of lines defining a spectrum in mgf
##'     format.
##'
##' @author Michael Witting
##'
##' @noRd
.extract_mb_ac <- function(mb) {

  # create empty list
  ac <- list()

  # analytical chemistry information, MS instrument ----------------------------
  ac$instrument <- substring(grep("AC$INSTRUMENT:", mb, value = TRUE, fixed = TRUE), 16)
  ac$instrument_type <- substring(grep("AC$INSTRUMENT_TYPE:", mb, value = TRUE, fixed = TRUE), 21)

  # analytical chemistry information, MS settings ------------------------------
  ac$ms_ms_type <- substring(grep("AC$MASS_SPECTROMETRY: MS_TYPE", mb, value = TRUE, fixed = TRUE), 31)
  ac$ms_cap_voltage <- substring(grep("AC$MASS_SPECTROMETRY: CAPILLARY_VOLTAGE", mb, value = TRUE, fixed = TRUE), 41)
  ac$ms_ion_mode <- substring(grep("AC$MASS_SPECTROMETRY: ION_MODE", mb, value = TRUE, fixed = TRUE), 32)
  ac$ms_col_energy <- substring(grep("AC$MASS_SPECTROMETRY: COLLISION_ENERGY", mb, value = TRUE, fixed = TRUE), 40)
  ac$ms_col_gas <- substring(grep("AC$MASS_SPECTROMETRY: COLLISION_GAS", mb, value = TRUE, fixed = TRUE), 37)
  ac$ms_desolv_gas_flow <- substring(grep("AC$MASS_SPECTROMETRY: DESOLVATION_GAS_FLOW", mb, value = TRUE, fixed = TRUE), 44)
  ac$ms_desolv_temp <- substring(grep("AC$MASS_SPECTROMETRY: DESOLVATION_TEMPERATURE", mb, value = TRUE, fixed = TRUE), 47)
  ac$ms_frag_mode <- substring(grep("AC$MASS_SPECTROMETRY: FRAGMENTATION_MODE", mb, value = TRUE, fixed = TRUE), 42)
  ac$ms_ionization <- substring(grep("AC$MASS_SPECTROMETRY: IONIZATION", mb, value = TRUE, fixed = TRUE), 34)
  ac$ms_ionization_energy <- substring(grep("AC$MASS_SPECTROMETRY: IONIZATION_ENERGY", mb, value = TRUE, fixed = TRUE), 41)
  ac$ms_laser <- substring(grep("AC$MASS_SPECTROMETRY: LASER", mb, value = TRUE, fixed = TRUE), 29)
  ac$ms_matrix <- substring(grep("AC$MASS_SPECTROMETRY: MATRIX", mb, value = TRUE, fixed = TRUE), 30)
  ac$ms_mass_accuracy <- substring(grep("AC$MASS_SPECTROMETRY: MASS_ACCURACY", mb, value = TRUE, fixed = TRUE), 37)
  ac$ms_mass_range <- substring(grep("AC$MASS_SPECTROMETRY: MASS_RANGE_MZ", mb, value = TRUE, fixed = TRUE), 37)
  ac$ms_reagent_gas <- substring(grep("AC$MASS_SPECTROMETRY: REAGENT_GAS", mb, value = TRUE, fixed = TRUE), 35)
  ac$ms_resolution <- substring(grep("AC$MASS_SPECTROMETRY: RESOLUTION", mb, value = TRUE, fixed = TRUE), 34)
  ac$ms_scan_setting <- substring(grep("AC$MASS_SPECTROMETRY: SCANNING_SETTING", mb, value = TRUE, fixed = TRUE), 40)
  ac$ms_source_temp <- substring(grep("AC$MASS_SPECTROMETRY: SOURCE_TEMPERATURE", mb, value = TRUE, fixed = TRUE), 42)

  # analytical chemistry information, chromatography ---------------------------
  ac$chrom_carrier_gas <- substring(grep("AC$CHROMATOGRAPHY: CARRIER_GAS", mb, value = TRUE, fixed = TRUE), 32)
  ac$chrom_column <- substring(grep("AC$CHROMATOGRAPHY: COLUMN_NAME", mb, value = TRUE, fixed = TRUE), 32)
  ac$chrom_column_temp <- substring(grep("AC$CHROMATOGRAPHY: COLUMN_TEMPERATURE", mb, value = TRUE, fixed = TRUE), 39)
  ac$chrom_column_temp_gradient <- substring(grep("AC$CHROMATOGRAPHY: COLUMN_TEMPERATURE_GRADIENT", mb, value = TRUE, fixed = TRUE), 48)
  ac$chrom_flow_gradient <- substring(grep("AC$CHROMATOGRAPHY: FLOW_GRADIENT", mb, value = TRUE, fixed = TRUE), 34)
  ac$chrom_flow_rate <- substring(grep("AC$CHROMATOGRAPHY: FLOW_RATE", mb, value = TRUE, fixed = TRUE), 30)
  ac$chrom_inj_temp <- substring(grep("AC$CHROMATOGRAPHY: INJECTION_TEMPERATURE", mb, value = TRUE, fixed = TRUE), 42)
  ac$chrom_inj_temp_gradient <- substring(grep("AC$CHROMATOGRAPHY: INJECTION_TEMPERATURE_GRADIENT", mb, value = TRUE, fixed = TRUE), 51)
  ac$chrom_rti_kovats <- substring(grep("AC$CHROMATOGRAPHY: KOVATS_RTI", mb, value = TRUE, fixed = TRUE), 31)
  ac$chrom_rti_lee <- substring(grep("AC$CHROMATOGRAPHY: LEE_RTI", mb, value = TRUE, fixed = TRUE), 28)
  ac$chrom_rti_naps <- substring(grep("AC$CHROMATOGRAPHY: NAPS_RTI", mb, value = TRUE, fixed = TRUE), 29)
  ac$chrom_rti_uoa <- substring(grep("AC$CHROMATOGRAPHY: UOA_RTI", mb, value = TRUE, fixed = TRUE), 28)
  ac$chrom_rti_uoa_pred <- substring(grep("AC$CHROMATOGRAPHY: UOA_PREDICTED_RTI", mb, value = TRUE, fixed = TRUE), 38)
  ac$chrom_rt <- substring(grep("AC$CHROMATOGRAPHY: RETENTION_TIME", mb, value = TRUE, fixed = TRUE), 35)
  ac$chrom_rt_uoa_pred <- substring(grep("AC$CHROMATOGRAPHY: TRAMS_PREDICTED_RETENTION_TIME", mb, value = TRUE, fixed = TRUE), 51)
  ac$chrom_solvent <- as.list(substring(grep("AC$CHROMATOGRAPHY: SOLVENT", mb, value = TRUE, fixed = TRUE), 28))
  ac$chrom_transfer_temp <- substring(grep("AC$CHROMATOGRAPHY: TRANSFERLINE_TEMPERATURE", mb, value = TRUE, fixed = TRUE), 45)

  # # analytical chemistry information, ion mobility
  # # preparation for IMS update of MassBank format
  # ac$ims_instrument_type <- substring(grep("AC$ION_MOBILITY: INSTRUMENT_TYPE", mb, value = TRUE, fixed = TRUE), 34)
  # ac$ims_drift_gas <- substring(grep("AC$ION_MOBILITY: DRIFT_GAS", mb, value = TRUE, fixed = TRUE), 28)
  # ac$ims_drift_time <- substring(grep("AC$ION_MOBILITY: DRIFT_TIME", mb, value = TRUE, fixed = TRUE), 29)
  # ac$ims_ccs <- substring(grep("AC$ION_MOBILITY: CCS", mb, value = TRUE, fixed = TRUE), 22)

  # analytical chemistry information, general ----------------------------------
  ac$general_conc <- substring(grep("AC$GENERAL: CONCENTRATION", mb, value = TRUE, fixed = TRUE), 27)

  ac

}

##' @param mb `character()` of lines defining a spectrum in mgf
##'     format.
##'
##' @author Michael Witting
##'
##' @noRd
.extract_mb_ch <- function(mb) {

  # create empty list
  ch <- list()

  # isolate chemical information
  ch$name <- as.list(substring(grep("CH$NAME:", mb, value = TRUE, fixed = TRUE), 10))
  ch$compound_class <- substring(grep("CH$COMPOUND_CLASS:", mb, value = TRUE, fixed = TRUE), 20)
  ch$formula <- substring(grep("CH$FORMULA:", mb, value = TRUE, fixed = TRUE), 13)
  ch$exact_mass <- as.numeric(substring(grep("CH$EXACT_MASS:", mb, value = TRUE, fixed = TRUE), 16))
  ch$smiles <- substring(grep("CH$SMILES:", mb, value = TRUE, fixed = TRUE), 12)
  ch$iupac <- substring(grep("CH$IUPAC:", mb, value = TRUE, fixed = TRUE), 11)
  ch$link_cas <- substring(grep("CH$LINK: CAS", mb, value = TRUE, fixed = TRUE), 14)
  ch$link_cayman <- substring(grep("CH$LINK: CAYMAN", mb, value = TRUE, fixed = TRUE), 17)
  ch$link_chebi <- substring(grep("CH$LINK: CHEBI", mb, value = TRUE, fixed = TRUE), 16)
  ch$link_chembl <- substring(grep("CH$LINK: CHEMBL", mb, value = TRUE, fixed = TRUE), 17)
  ch$link_chempdb <- substring(grep("CH$LINK: CHEMPDB", mb, value = TRUE, fixed = TRUE), 18)
  ch$link_chemspider <- substring(grep("CH$LINK: CHEMSPIDER", mb, value = TRUE, fixed = TRUE), 21)
  ch$link_comptox <- substring(grep("CH$LINK: COMPTOX", mb, value = TRUE, fixed = TRUE), 18)
  ch$link_hmdb <- substring(grep("CH$LINK: HMDB", mb, value = TRUE, fixed = TRUE), 15)
  ch$link_inchikey <- substring(grep("CH$LINK: INCHIKEY", mb, value = TRUE, fixed = TRUE), 19)
  ch$link_kappaview <- substring(grep("CH$LINK: KAPPAVIEW", mb, value = TRUE, fixed = TRUE), 20)
  ch$link_kegg <- substring(grep("CH$LINK: KEGG", mb, value = TRUE, fixed = TRUE), 15)
  ch$link_knapsack <- substring(grep("CH$LINK: KNAPSACK", mb, value = TRUE, fixed = TRUE), 19)
  ch$link_lipidbank <- substring(grep("CH$LINK: LIPIDBANK", mb, value = TRUE, fixed = TRUE), 20)
  ch$link_lipidmaps <- substring(grep("CH$LINK: LIPIDMAPS", mb, value = TRUE, fixed = TRUE), 20)
  ch$link_nikkaji <- substring(grep("CH$LINK: NIKKAJI", mb, value = TRUE, fixed = TRUE), 18)
  ch$link_pubchem <- substring(grep("CH$LINK: PUBCHEM", mb, value = TRUE, fixed = TRUE), 18)
  ch$link_zinc <- substring(grep("CH$LINK: ZINC", mb, value = TRUE, fixed = TRUE), 15)

  # clean up data
  # TODO replace character(0) with NA_character

  # return result list
  ch

}

##' @param mb `character()` of lines defining a spectrum in mgf
##'     format.
##'
##' @author Michael Witting
##'
##' @noRd
.extract_mb_sp <- function(mb) {

  # create empty list
  sp <- list()

  # species information
  sp$scientific_name <- substring(grep("SP$SCIENTIFIC_NAME:", mb, value = TRUE, fixed = TRUE), 21)
  sp$lineage <- substring(grep("SP$LINEAGE:", mb, value = TRUE, fixed = TRUE), 13)
  sp$link <- substring(grep("SP$LINK:", mb, value = TRUE, fixed = TRUE), 10)
  sp$sample <- substring(grep("SP$SAMPLE:", mb, value = TRUE, fixed = TRUE), 12)

  # clean up data
  # TODO replace character(0) with NA_character

  # return result list
  sp

}

##' @param mb `character()` of lines defining a spectrum in mgf
##'     format.
##'
##' @author Michael Witting
##'
##' @noRd
.extract_mb_ms <- function(mb) {

  # create empty list
  ms <- list()

  # MS information, base peak
  ms$focus_base_peak <- substring(grep("MS$FOCUSED_ION: BASE_PEAK", mb, value = TRUE, fixed = TRUE), 27)

  # MS information, derivative
  ms$focus_derivative_form <- substring(grep("MS$FOCUSED_ION: DERIVATIVE_FORM", mb, value = TRUE, fixed = TRUE), 33)
  ms$focus_derivative_mass <- substring(grep("MS$FOCUSED_ION: DERIVATIVE_MASS", mb, value = TRUE, fixed = TRUE), 33)
  ms$focus_derivative_type <- substring(grep("MS$FOCUSED_ION: DERIVATIVE_TYPE", mb, value = TRUE, fixed = TRUE), 33)

  # MS information, precursor
  ms$focus_ion_type <- substring(grep("MS$FOCUSED_ION: ION_TYPE", mb, value = TRUE, fixed = TRUE), 26)
  ms$focus_precursor_int <- substring(grep("MS$FOCUSED_ION: PRECURSOR_INT", mb, value = TRUE, fixed = TRUE), 31)
  ms$focus_precursor_mz <- substring(grep("MS$FOCUSED_ION: PRECURSOR_MZ", mb, value = TRUE, fixed = TRUE), 30)
  ms$focus_precursor_type <- substring(grep("MS$FOCUSED_ION: PRECURSOR_TYPE", mb, value = TRUE, fixed = TRUE), 32)

  # MS data processing
  ms$data_processing_comment <- substring(grep("MS$DATA_PROCESSING: COMMENT", mb, value = TRUE, fixed = TRUE), 29)
  ms$data_processing_deprofile <- substring(grep("MS$DATA_PROCESSING: DEPROFILE", mb, value = TRUE, fixed = TRUE), 31)
  ms$data_processing_find <- substring(grep("MS$DATA_PROCESSING: FIND_PEAK", mb, value = TRUE, fixed = TRUE), 31)
  ms$data_processing_reanalyze <- substring(grep("MS$DATA_PROCESSING: REANALYZE", mb, value = TRUE, fixed = TRUE), 31)
  ms$data_processing_recalibrate <- substring(grep("MS$DATA_PROCESSING: RECALIBRATE", mb, value = TRUE, fixed = TRUE), 33)
  ms$data_processing_whole <- substring(grep("MS$DATA_PROCESSING: WHOLE", mb, value = TRUE, fixed = TRUE), 27)

  # clean up data
  # TODO replace character(0) with NA_character_

  # TODO type conversion for numeric data

  # return result list
  ms

}

##' @param mb `character()` of lines defining a spectrum in mgf
##'     format.
##'
##' @author Michael Witting
##'
##' @noRd
.extract_mb_record <- function(mb) {

  # create empty list
  recordinfo <- list()

  # mb information
  recordinfo$accession <- substring(grep("ACCESSION:", mb, value = TRUE, fixed = TRUE), 12)
  recordinfo$deprecated <- substring(grep("DEPRECATED:", mb, value = TRUE, fixed = TRUE), 13)
  recordinfo$record_title <- substring(grep("RECORD_TITLE:", mb, value = TRUE, fixed = TRUE), 15)
  recordinfo$date <- substring(grep("DATE:", mb, value = TRUE, fixed = TRUE), 7)
  recordinfo$authors <- substring(grep("AUTHORS:", mb, value = TRUE, fixed = TRUE), 10)
  recordinfo$license <- substring(grep("LICENSE:", mb, value = TRUE, fixed = TRUE), 10)
  recordinfo$copyright <- substring(grep("COPYRIGHT:", mb, value = TRUE, fixed = TRUE), 12)
  recordinfo$publication <- substring(grep("PUBLICATION:", mb, value = TRUE, fixed = TRUE), 14)
  recordinfo$project <- as.list(substring(grep("PROJECT:", mb, value = TRUE, fixed = TRUE), 10))

  # clean up data
  # TODO replace character(0) with NA_character_

  # TODO type conversion for dates

  # return result list
  recordinfo

}

##' @param mb `character()` of lines defining a spectrum in mgf
##'     format.
##'
##' @author Michael Witting
##'
##' @noRd
.extract_mb_pk <- function(mb) {

  # create empty list
  pk <- list()

  # peak data
  pk$splash <- substring(grep("PK$SPLASH:", mb, value = TRUE, fixed = TRUE), 12)
  pk$num <- substring(grep("PK$NUM_PEAK:", mb, value = TRUE, fixed = TRUE), 14)

  # clean up data
  # TODO replace character(0) with NA_character_

  # TODO type conversion for numeric data

  # return result list
  pk

}

##' @param mb `character()` of lines defining a spectrum in mgf
##'     format.
##'
##' @author Michael Witting
##'
##' @noRd
.extract_mb_comment <- function(mb) {

  # comment section
  comment <- as.list(substring(grep("COMMENT:", mb, value = TRUE, fixed = TRUE), 10))

  # return result list
  comment

}



##' It is possible to
##'
##' @title Metadata blocks to be read
##'
##' @return A `data.frame` with metadata blocks.
##'
##' @author Laurent Gatto
##'
##' @importFrom utils read.csv
##'
##' @export
##'
##' @examples
##' metaDataBlocks()
metaDataBlocks <- function() {
  read.csv(dir(system.file("extdata", package = "MsBackendMassbank"),
               pattern = "metadata_blocks.csv",
               full.names = TRUE),
           header = TRUE,
           stringsAsFactors = FALSE)

}

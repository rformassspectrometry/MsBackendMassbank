# all fields are found here
# https://github.com/MassBank/MassBank-web/blob/dev/Documentation/MassBankRecordFormat.md

#'
#'
#'
.import_massbank_ms_ms_spectrum <- function(x) {

  if (!is.character(x) || length(x) != 1)
    stop("'x' has to be of type character with length 1")

  con <- file(x)
  mb_record <- readLines(con)
  close(con)

  # isolate spectrum and basic metadata (default)
  recordinfo <- .isolate_recordinfo(mb_record)
  spectrum <- .isolate_spectrum(mb_record)

  # idea: isolate only additional information if requested
  if(isolate_ch) {
    ch <- .isolate_chemical(mb_record)
  } else {
    ch <- NULL
  }

  # create Spectra object with metadata

}

#'
#'
#'
.isolate_spectrum <- function(mb_record) {

  #read the spectrum
  spectrum_start <- grep("PK$PEAK:", mb_record, fixed = TRUE) + 1
  spectrum_end <- tail(grep("//", mb_record, fixed = TRUE), 1) - 1

  if(spectrum_start < spectrum_end){
    splitted <- strsplit(mb_record[spectrum_start:(spectrum_end)]," ")
    spectrum <- matrix(nrow = spectrum_end + 1 - spectrum_start, ncol = 3)

    for(k in 1:length(splitted)){
      splitted[[k]] <- splitted[[k]][which(splitted[[k]] != "")]
      spectrum[k,] <- splitted[[k]]
    }

    # convert to data frame and adjust data type
    spectrum <- as.data.frame(spectrum, stringsAsFactors = FALSE)
    spectrum[] <- lapply(spectrum, type.convert)
    colnames(spectrum) <- c("m/z", "int", "rel.int.")

    # return spectrum
    return(spectrum)

  } else {

    stop("Could not find PK$PEAK or //. Please check the supplied MB record file")

  }
}

#'
#'
#'
.isolate_chemical <- function(mb_record) {

  # create empty list
  ch <- list()

  # isolate chemical information
  ch$name <- as.list(substring(grep("CH$NAME:", mb_record, value = TRUE, fixed = TRUE), 10))
  ch$compound_class <- substring(grep("CH$COMPOUND_CLASS:", mb_record, value = TRUE, fixed = TRUE), 20)
  ch$formula <- substring(grep("CH$FORMULA:", mb_record, value = TRUE, fixed = TRUE), 13)
  ch$exact_mass <- as.numeric(substring(grep("CH$EXACT_MASS:", mb_record, value = TRUE, fixed = TRUE), 16))
  ch$smiles <- substring(grep("CH$SMILES:", mb_record, value = TRUE, fixed = TRUE), 12)
  ch$iupac <- substring(grep("CH$IUPAC:", mb_record, value = TRUE, fixed = TRUE), 11)
  ch$link_cas <- substring(grep("CH$LINK: CAS", mb_record, value = TRUE, fixed = TRUE), 14)
  ch$link_cayman <- substring(grep("CH$LINK: CAYMAN", mb_record, value = TRUE, fixed = TRUE), 17)
  ch$link_chebi <- substring(grep("CH$LINK: CHEBI", mb_record, value = TRUE, fixed = TRUE), 16)
  ch$link_chembl <- substring(grep("CH$LINK: CHEMBL", mb_record, value = TRUE, fixed = TRUE), 17)
  ch$link_chempdb <- substring(grep("CH$LINK: CHEMPDB", mb_record, value = TRUE, fixed = TRUE), 18)
  ch$link_chemspider <- substring(grep("CH$LINK: CHEMSPIDER", mb_record, value = TRUE, fixed = TRUE), 21)
  ch$link_comptox <- substring(grep("CH$LINK: COMPTOX", mb_record, value = TRUE, fixed = TRUE), 18)
  ch$link_hmdb <- substring(grep("CH$LINK: HMDB", mb_record, value = TRUE, fixed = TRUE), 15)
  ch$link_inchikey <- substring(grep("CH$LINK: INCHIKEY", mb_record, value = TRUE, fixed = TRUE), 19)
  ch$link_kappaview <- substring(grep("CH$LINK: KAPPAVIEW", mb_record, value = TRUE, fixed = TRUE), 20)
  ch$link_kegg <- substring(grep("CH$LINK: KEGG", mb_record, value = TRUE, fixed = TRUE), 15)
  ch$link_knapsack <- substring(grep("CH$LINK: KNAPSACK", mb_record, value = TRUE, fixed = TRUE), 19)
  ch$link_lipidbank <- substring(grep("CH$LINK: LIPIDBANK", mb_record, value = TRUE, fixed = TRUE), 20)
  ch$link_lipidmaps <- substring(grep("CH$LINK: LIPIDMAPS", mb_record, value = TRUE, fixed = TRUE), 20)
  ch$link_nikkaji <- substring(grep("CH$LINK: NIKKAJI", mb_record, value = TRUE, fixed = TRUE), 18)
  ch$link_pubchem <- substring(grep("CH$LINK: PUBCHEM", mb_record, value = TRUE, fixed = TRUE), 18)
  ch$link_zinc <- substring(grep("CH$LINK: ZINC", mb_record, value = TRUE, fixed = TRUE), 15)

  # clean up data
  # TODO replace character(0) with NA_character

  # return result list
  ch
}

#'
#'
#'
.isolate_species <- function(mb_record) {

  # create empty list
  sp <- list()

  # species information
  sp$scientific_name <- substring(grep("SP$SCIENTIFIC_NAME:", mb_record, value = TRUE, fixed = TRUE), 21)
  sp$lineage <- substring(grep("SP$LINEAGE:", mb_record, value = TRUE, fixed = TRUE), 13)
  sp$link <- substring(grep("SP$LINK:", mb_record, value = TRUE, fixed = TRUE), 10)
  sp$sample <- substring(grep("SP$SAMPLE:", mb_record, value = TRUE, fixed = TRUE), 12)

  # clean up data
  # TODO replace character(0) with NA_character

  # return result list
  sp
}

#'
#'
#'
.isolate_analchem <- function(mb_record) {

  # create empty list
  ac <- list()

  # analytical chemistry information, MS instrument
  ac$instrument <- substring(grep("AC$INSTRUMENT:", mb_record, value = TRUE, fixed = TRUE), 16)
  ac$instrument_type <- substring(grep("AC$INSTRUMENT_TYPE:", mb_record, value = TRUE, fixed = TRUE), 21)

  # analytical chemistry information, MS settings
  ac$ms_ms_type <- substring(grep("AC$MASS_SPECTROMETRY: MS_TYPE", mb_record, value = TRUE, fixed = TRUE), 31)
  ac$ms_cap_voltage <- substring(grep("AC$MASS_SPECTROMETRY: CAPILLARY_VOLTAGE", mb_record, value = TRUE, fixed = TRUE), 41)
  ac$ms_ion_mode <- substring(grep("AC$MASS_SPECTROMETRY: ION_MODE", mb_record, value = TRUE, fixed = TRUE), 32)
  ac$ms_col_energy <- substring(grep("AC$MASS_SPECTROMETRY: COLLISION_ENERGY", mb_record, value = TRUE, fixed = TRUE), 40)
  ac$ms_col_gas <- substring(grep("AC$MASS_SPECTROMETRY: COLLISION_GAS", mb_record, value = TRUE, fixed = TRUE), 37)
  ac$ms_desolv_gas_flow <- substring(grep("AC$MASS_SPECTROMETRY: DESOLVATION_GAS_FLOW", mb_record, value = TRUE, fixed = TRUE), 44)
  ac$ms_desolv_temp <- substring(grep("AC$MASS_SPECTROMETRY: DESOLVATION_TEMPERATURE", mb_record, value = TRUE, fixed = TRUE), 47)
  ac$ms_frag_mode <- substring(grep("AC$MASS_SPECTROMETRY: FRAGMENTATION_MODE", mb_record, value = TRUE, fixed = TRUE), 42)
  ac$ms_ionization <- substring(grep("AC$MASS_SPECTROMETRY: IONIZATION", mb_record, value = TRUE, fixed = TRUE), 34)
  ac$ms_ionization_energy <- substring(grep("AC$MASS_SPECTROMETRY: IONIZATION_ENERGY", mb_record, value = TRUE, fixed = TRUE), 41)
  ac$ms_laser <- substring(grep("AC$MASS_SPECTROMETRY: LASER", mb_record, value = TRUE, fixed = TRUE), 29)
  ac$ms_matrix <- substring(grep("AC$MASS_SPECTROMETRY: MATRIX", mb_record, value = TRUE, fixed = TRUE), 30)
  ac$ms_mass_accuracy <- substring(grep("AC$MASS_SPECTROMETRY: MASS_ACCURACY", mb_record, value = TRUE, fixed = TRUE), 37)
  ac$ms_mass_range <- substring(grep("AC$MASS_SPECTROMETRY: MASS_RANGE_MZ", mb_record, value = TRUE, fixed = TRUE), 37)
  ac$ms_reagent_gas <- substring(grep("AC$MASS_SPECTROMETRY: REAGENT_GAS", mb_record, value = TRUE, fixed = TRUE), 35)
  ac$ms_resolution <- substring(grep("AC$MASS_SPECTROMETRY: RESOLUTION", mb_record, value = TRUE, fixed = TRUE), 34)
  ac$ms_scan_setting <- substring(grep("AC$MASS_SPECTROMETRY: SCANNING_SETTING", mb_record, value = TRUE, fixed = TRUE), 40)
  ac$ms_source_temp <- substring(grep("AC$MASS_SPECTROMETRY: SOURCE_TEMPERATURE", mb_record, value = TRUE, fixed = TRUE), 42)

  # analytical chemistry information, chromatography
  ac$chrom_carrier_gas <- substring(grep("AC$CHROMATOGRAPHY: CARRIER_GAS", mb_record, value = TRUE, fixed = TRUE), 32)
  ac$chrom_column <- substring(grep("AC$CHROMATOGRAPHY: COLUMN_NAME", mb_record, value = TRUE, fixed = TRUE), 32)
  ac$chrom_column_temp <- substring(grep("AC$CHROMATOGRAPHY: COLUMN_TEMPERATURE", mb_record, value = TRUE, fixed = TRUE), 39)
  ac$chrom_column_temp_gradient <- substring(grep("AC$CHROMATOGRAPHY: COLUMN_TEMPERATURE_GRADIENT", mb_record, value = TRUE, fixed = TRUE), 48)
  ac$chrom_flow_gradient <- substring(grep("AC$CHROMATOGRAPHY: FLOW_GRADIENT", mb_record, value = TRUE, fixed = TRUE), 34)
  ac$chrom_flow_rate <- substring(grep("AC$CHROMATOGRAPHY: FLOW_RATE", mb_record, value = TRUE, fixed = TRUE), 30)
  ac$chrom_inj_temp <- substring(grep("AC$CHROMATOGRAPHY: INJECTION_TEMPERATURE", mb_record, value = TRUE, fixed = TRUE), 42)
  ac$chrom_inj_temp_gradient <- substring(grep("AC$CHROMATOGRAPHY: INJECTION_TEMPERATURE_GRADIENT", mb_record, value = TRUE, fixed = TRUE), 51)
  ac$chrom_rti_kovats <- substring(grep("AC$CHROMATOGRAPHY: KOVATS_RTI", mb_record, value = TRUE, fixed = TRUE), 31)
  ac$chrom_rti_lee <- substring(grep("AC$CHROMATOGRAPHY: LEE_RTI", mb_record, value = TRUE, fixed = TRUE), 28)
  ac$chrom_rti_naps <- substring(grep("AC$CHROMATOGRAPHY: NAPS_RTI", mb_record, value = TRUE, fixed = TRUE), 29)
  ac$chrom_rti_uoa <- substring(grep("AC$CHROMATOGRAPHY: UOA_RTI", mb_record, value = TRUE, fixed = TRUE), 28)
  ac$chrom_rti_uoa_pred <- substring(grep("AC$CHROMATOGRAPHY: UOA_PREDICTED_RTI", mb_record, value = TRUE, fixed = TRUE), 38)
  ac$chrom_rt <- substring(grep("AC$CHROMATOGRAPHY: RETENTION_TIME", mb_record, value = TRUE, fixed = TRUE), 35)
  ac$chrom_rt_uoa_pred <- substring(grep("AC$CHROMATOGRAPHY: TRAMS_PREDICTED_RETENTION_TIME", mb_record, value = TRUE, fixed = TRUE), 51)
  ac$chrom_solvent <- as.list(substring(grep("AC$CHROMATOGRAPHY: SOLVENT", mb_record, value = TRUE, fixed = TRUE), 28))
  ac$chrom_transfer_temp <- substring(grep("AC$CHROMATOGRAPHY: TRANSFERLINE_TEMPERATURE", mb_record, value = TRUE, fixed = TRUE), 45)

  # # analytical chemistry information, ion mobility
  # # preparation for IMS update of MassBank format
  # ac$ims_instrument_type <- substring(grep("AC$ION_MOBILITY: INSTRUMENT_TYPE", mb_record, value = TRUE, fixed = TRUE), 34)
  # ac$ims_drift_gas <- substring(grep("AC$ION_MOBILITY: DRIFT_GAS", mb_record, value = TRUE, fixed = TRUE), 28)
  # ac$ims_drift_time <- substring(grep("AC$ION_MOBILITY: DRIFT_TIME", mb_record, value = TRUE, fixed = TRUE), 29)
  # ac$ims_ccs <- substring(grep("AC$ION_MOBILITY: CCS", mb_record, value = TRUE, fixed = TRUE), 22)

  # analytical chemistry information, general
  ac$general_conc <- substring(grep("AC$GENERAL: CONCENTRATION", mb_record, value = TRUE, fixed = TRUE), 27)

  # clean up data
  # TODO replace character(0) with NA_character

  # return result list
  ac
}

#'
#'
#'
.isolate_massspec <- function(mb_record) {

  # create empty list
  ms <- list()

  # MS information, base peak
  ms$focus_base_peak <- substring(grep("MS$FOCUSED_ION: BASE_PEAK", mb_record, value = TRUE, fixed = TRUE), 27)

  # MS information, derivative
  ms$focus_derivative_form <- substring(grep("MS$FOCUSED_ION: DERIVATIVE_FORM", mb_record, value = TRUE, fixed = TRUE), 33)
  ms$focus_derivative_mass <- substring(grep("MS$FOCUSED_ION: DERIVATIVE_MASS", mb_record, value = TRUE, fixed = TRUE), 33)
  ms$focus_derivative_type <- substring(grep("MS$FOCUSED_ION: DERIVATIVE_TYPE", mb_record, value = TRUE, fixed = TRUE), 33)

  # MS information, precursor
  ms$focus_ion_type <- substring(grep("MS$FOCUSED_ION: ION_TYPE", mb_record, value = TRUE, fixed = TRUE), 26)
  ms$focus_precursor_int <- substring(grep("MS$FOCUSED_ION: PRECURSOR_INT", mb_record, value = TRUE, fixed = TRUE), 31)
  ms$focus_precursor_mz <- substring(grep("MS$FOCUSED_ION: PRECURSOR_MZ", mb_record, value = TRUE, fixed = TRUE), 30)
  ms$focus_precursor_type <- substring(grep("MS$FOCUSED_ION: PRECURSOR_TYPE", mb_record, value = TRUE, fixed = TRUE), 32)

  # MS data processing
  ms$data_processing_comment <- substring(grep("MS$DATA_PROCESSING: COMMENT", mb_record, value = TRUE, fixed = TRUE), 29)
  ms$data_processing_deprofile <- substring(grep("MS$DATA_PROCESSING: DEPROFILE", mb_record, value = TRUE, fixed = TRUE), 31)
  ms$data_processing_find <- substring(grep("MS$DATA_PROCESSING: FIND_PEAK", mb_record, value = TRUE, fixed = TRUE), 31)
  ms$data_processing_reanalyze <- substring(grep("MS$DATA_PROCESSING: REANALYZE", mb_record, value = TRUE, fixed = TRUE), 31)
  ms$data_processing_recalibrate <- substring(grep("MS$DATA_PROCESSING: RECALIBRATE", mb_record, value = TRUE, fixed = TRUE), 33)
  ms$data_processing_whole <- substring(grep("MS$DATA_PROCESSING: WHOLE", mb_record, value = TRUE, fixed = TRUE), 27)

  # clean up data
  # TODO replace character(0) with NA_character_

  # TODO type conversion for numeric data

  # return result list
  ms
}

#'
#'
#'
.isolate_recordinfo <- function(mb_record) {

  # create empty list
  recordinfo <- list()

  # mb_record information
  recordinfo$accession <- substring(grep("ACCESSION:", mb_record, value = TRUE, fixed = TRUE), 12)
  recordinfo$deprecated <- substring(grep("DEPRECATED:", mb_record, value = TRUE, fixed = TRUE), 13)
  recordinfo$record_title <- substring(grep("RECORD_TITLE:", mb_record, value = TRUE, fixed = TRUE), 15)
  recordinfo$date <- substring(grep("DATE:", mb_record, value = TRUE, fixed = TRUE), 7)
  recordinfo$authors <- substring(grep("AUTHORS:", mb_record, value = TRUE, fixed = TRUE), 10)
  recordinfo$license <- substring(grep("LICENSE:", mb_record, value = TRUE, fixed = TRUE), 10)
  recordinfo$copyright <- substring(grep("COPYRIGHT:", mb_record, value = TRUE, fixed = TRUE), 12)
  recordinfo$publication <- substring(grep("PUBLICATION:", mb_record, value = TRUE, fixed = TRUE), 14)
  recordinfo$project <- as.list(substring(grep("PROJECT:", mb_record, value = TRUE, fixed = TRUE), 10))

  # clean up data
  # TODO replace character(0) with NA_character_

  # TODO type conversion for dates

  # return result list
  recordinfo

}

#'
#'
#'
.isolate_peakinfo <- function(mb_record) {

  # create empty list
  pk <- list()

  # peak data
  pk$splash <- substring(grep("PK$SPLASH:", mb_record, value = TRUE, fixed = TRUE), 12)
  pk$num <- substring(grep("PK$NUM_PEAK:", mb_record, value = TRUE, fixed = TRUE), 14)

  # clean up data
  # TODO replace character(0) with NA_character_

  # TODO type conversion for numeric data

  # return result list
  pk

}

#'
#'
#'
.isolate_comment <- function(mb_record) {

  # comment section
  comment <- as.list(substring(grep("COMMENT:", mb_record, value = TRUE, fixed = TRUE), 10))

  # return result list
  comment

}

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

}

#'
#'
#'
.isolate_chemical <- function(record) {

  # create empty list
  ch <- list()

  # isolate chemical information
  ch$name <- as.list(substring(grep("CH$NAME:", record, value = TRUE, fixed = TRUE), 10))
  ch$compound_class <- substring(grep("CH$COMPOUND_CLASS:", record, value = TRUE, fixed = TRUE), 20)
  ch$formula <- substring(grep("CH$FORMULA:", record, value = TRUE, fixed = TRUE), 13)
  ch$exact_mass <- as.numeric(substring(grep("CH$EXACT_MASS:", record, value = TRUE, fixed = TRUE), 16))
  ch$smiles <- substring(grep("CH$SMILES:", record, value = TRUE, fixed = TRUE), 12)
  ch$iupac <- substring(grep("CH$IUPAC:", record, value = TRUE, fixed = TRUE), 11)
  ch$link_cas <- substring(grep("CH$LINK: CAS", record, value = TRUE, fixed = TRUE), 14)
  ch$link_cayman <- substring(grep("CH$LINK: CAYMAN", record, value = TRUE, fixed = TRUE), 17)
  ch$link_chebi <- substring(grep("CH$LINK: CHEBI", record, value = TRUE, fixed = TRUE), 16)
  ch$link_chembl <- substring(grep("CH$LINK: CHEMBL", record, value = TRUE, fixed = TRUE), 17)
  ch$link_chempdb <- substring(grep("CH$LINK: CHEMPDB", record, value = TRUE, fixed = TRUE), 18)
  ch$link_chemspider <- substring(grep("CH$LINK: CHEMSPIDER", record, value = TRUE, fixed = TRUE), 21)
  ch$link_comptox <- substring(grep("CH$LINK: COMPTOX", record, value = TRUE, fixed = TRUE), 18)
  ch$link_hmdb <- substring(grep("CH$LINK: HMDB", record, value = TRUE, fixed = TRUE), 15)
  ch$link_inchikey <- substring(grep("CH$LINK: INCHIKEY", record, value = TRUE, fixed = TRUE), 19)
  ch$link_kappaview <- substring(grep("CH$LINK: KAPPAVIEW", record, value = TRUE, fixed = TRUE), 20)
  ch$link_kegg <- substring(grep("CH$LINK: KEGG", record, value = TRUE, fixed = TRUE), 15)
  ch$link_knapsack <- substring(grep("CH$LINK: KNAPSACK", record, value = TRUE, fixed = TRUE), 19)
  ch$link_lipidbank <- substring(grep("CH$LINK: LIPIDBANK", record, value = TRUE, fixed = TRUE), 20)
  ch$link_lipidmaps <- substring(grep("CH$LINK: LIPIDMAPS", record, value = TRUE, fixed = TRUE), 20)
  ch$link_nikkaji <- substring(grep("CH$LINK: NIKKAJI", record, value = TRUE, fixed = TRUE), 18)
  ch$link_pubchem <- substring(grep("CH$LINK: PUBCHEM", record, value = TRUE, fixed = TRUE), 18)
  ch$link_zinc <- substring(grep("CH$LINK: ZINC", record, value = TRUE, fixed = TRUE), 15)

  # clean up data
  # TODO replace character(0) with NA_character

  # return result list
  ch
}

#'
#'
#'
.isolate_species <- function(record) {

  # create empty list
  sp <- list()

  # species information
  # TODO fix all the length
  sp$scientific_name <- substring(grep("SP$SCIENTIFIC_NAME:", record, value = TRUE, fixed = TRUE), 21)
  sp$lineage <- substring(grep("SP$LINEAGE:", record, value = TRUE, fixed = TRUE), 13)
  sp$link <- substring(grep("SP$LINK:", record, value = TRUE, fixed = TRUE), 10)
  sp$sample <- substring(grep("SP$SAMPLE:", record, value = TRUE, fixed = TRUE), 12)

  # clean up data
  # TODO replace character(0) with NA_character

  # return result list
  sp
}

#'
#'
#'
.isolate_analchem <- function(record) {

  # create empty list
  ac <- list()

  # analytical chemistry information, MS
  # TODO fix all the length
  ac$instrument <- substring(grep("AC$INSTRUMENT:", record, value = TRUE, fixed = TRUE), 12)
  ac$instrument_type <- substring(grep("AC$INSTRUMENT_TYPE:", record, value = TRUE, fixed = TRUE), 12)
  ac$ms_ms_type <- substring(grep("AC$MASS_SPECTROMETRY: MS_TYPE", record, value = TRUE, fixed = TRUE), 12)
  ac$ms_cap_voltage <- substring(grep("AC$MASS_SPECTROMETRY: CAPILLARY_VOLTAGE", record, value = TRUE, fixed = TRUE), 12)
  ac$ms_ion_mode <- substring(grep("AC$MASS_SPECTROMETRY: ION_MODE", record, value = TRUE, fixed = TRUE), 12)
  ac$ms_col_energy <- substring(grep("AC$MASS_SPECTROMETRY: COLLISION_ENERGY", record, value = TRUE, fixed = TRUE), 12)
  ac$ms_col_gas <- substring(grep("AC$MASS_SPECTROMETRY: COLLISION_GAS", record, value = TRUE, fixed = TRUE), 12)
  ac$ms_desolv_gas_flow <- substring(grep("AC$MASS_SPECTROMETRY: DESOLVATION_GAS_FLOW", record, value = TRUE, fixed = TRUE), 12)
  ac$ms_desolv_temp <- substring(grep("AC$MASS_SPECTROMETRY: DESOLVATION_TEMPERATURE", record, value = TRUE, fixed = TRUE), 12)
  ac$ms_frag_mode <- substring(grep("AC$MASS_SPECTROMETRY: FRAGMENTATION_MODE", record, value = TRUE, fixed = TRUE), 12)
  ac$ms_ionization <- substring(grep("AC$MASS_SPECTROMETRY: IONIZATION", record, value = TRUE, fixed = TRUE), 12)
  ac$ms_ionization_energy <- substring(grep("AC$MASS_SPECTROMETRY: IONIZATION_ENERGY", record, value = TRUE, fixed = TRUE), 12)
  ac$ms_laser <- substring(grep("AC$MASS_SPECTROMETRY: LASER", record, value = TRUE, fixed = TRUE), 12)
  ac$ms_matrix <- substring(grep("AC$MASS_SPECTROMETRY: MATRIX", record, value = TRUE, fixed = TRUE), 12)
  ac$ms_mass_accuracy <- substring(grep("AC$MASS_SPECTROMETRY: MASS_ACCURACY", record, value = TRUE, fixed = TRUE), 12)
  ac$ms_mass_range <- substring(grep("AC$MASS_SPECTROMETRY: MASS_RANGE_MZ", record, value = TRUE, fixed = TRUE), 12)
  ac$ms_reagent_gas <- substring(grep("AC$MASS_SPECTROMETRY: REAGENT_GAS", record, value = TRUE, fixed = TRUE), 12)
  ac$ms_resolution <- substring(grep("AC$MASS_SPECTROMETRY: RESOLUTION", record, value = TRUE, fixed = TRUE), 12)
  ac$ms_scan_setting <- substring(grep("AC$MASS_SPECTROMETRY: SCANNING_SETTING", record, value = TRUE, fixed = TRUE), 12)
  ac$ms_source_temp <- substring(grep("AC$MASS_SPECTROMETRY: SOURCE_TEMPERATURE", record, value = TRUE, fixed = TRUE), 12)

  # analytical chemistry information, chromatography
  # TODO fix all the length
  ac$chrom_carrier_gas <- substring(grep("AC$CHROMATOGRAPHY: CARRIER_GAS", record, value = TRUE, fixed = TRUE), 12)
  ac$chrom_column <- substring(grep("AC$CHROMATOGRAPHY: COLUMN_NAME", record, value = TRUE, fixed = TRUE), 12)
  ac$chrom_column_temp <- substring(grep("AC$CHROMATOGRAPHY: COLUMN_TEMPERATURE", record, value = TRUE, fixed = TRUE), 12)
  ac$chrom_column_temp_gradient <- substring(grep("AC$CHROMATOGRAPHY: COLUMN_TEMPERATURE_GRADIENT", record, value = TRUE, fixed = TRUE), 12)
  ac$chrom_flow_gradient <- substring(grep("AC$CHROMATOGRAPHY: FLOW_GRADIENT", record, value = TRUE, fixed = TRUE), 12)
  ac$chrom_flow_rate <- substring(grep("AC$CHROMATOGRAPHY: FLOW_RATE 0.25 mL/min", record, value = TRUE, fixed = TRUE), 12)
  ac$chrom_inj_temp <- substring(grep("AC$CHROMATOGRAPHY: INJECTION_TEMPERATURE", record, value = TRUE, fixed = TRUE), 12)
  ac$chrom_inj_temp_gradient <- substring(grep("AC$CHROMATOGRAPHY: INJECTION_TEMPERATURE_GRADIENT", record, value = TRUE, fixed = TRUE), 12)
  ac$chrom_rti_kovats <- substring(grep("AC$CHROMATOGRAPHY: KOVATS_RTI", record, value = TRUE, fixed = TRUE), 12)
  ac$chrom_rti_lee <- substring(grep("AC$CHROMATOGRAPHY: LEE_RTI", record, value = TRUE, fixed = TRUE), 12)
  ac$chrom_rti_naps <- substring(grep("AC$CHROMATOGRAPHY: NAPS_RTI", record, value = TRUE, fixed = TRUE), 12)
  ac$chrom_rti_uoa <- substring(grep("AC$CHROMATOGRAPHY: UOA_RTI", record, value = TRUE, fixed = TRUE), 12)
  ac$chrom_rti_uoa_pred <- substring(grep("AC$CHROMATOGRAPHY: UOA_PREDICTED_RTI", record, value = TRUE, fixed = TRUE), 12)
  ac$chrom_rt <- substring(grep("AC$CHROMATOGRAPHY: RETENTION_TIME", record, value = TRUE, fixed = TRUE), 12)
  ac$chrom_rt_uoa_pred <- substring(grep("AC$CHROMATOGRAPHY: TRAMS_PREDICTED_RETENTION_TIME", record, value = TRUE, fixed = TRUE), 12)
  ac$chrom_solvent <- as.list(substring(grep("AC$CHROMATOGRAPHY: SOLVENT", record, value = TRUE, fixed = TRUE), 10))
  ac$chrom_transfer_temp <- substring(grep("AC$CHROMATOGRAPHY: TRANSFERLINE_TEMPERATURE", record, value = TRUE, fixed = TRUE), 12)

  # analytical chemistry information, general
  # TODO fix all the length
  ac$general_conc <- substring(grep("AC$GENERAL: CONCENTRATION", record, value = TRUE, fixed = TRUE), 12)

  # clean up data
  # TODO replace character(0) with NA_character

  # return result list
  ac
}

#'
#'
#'
.isolate_massspec <- function(record) {

  # create empty list
  ms <- list()

  # MS information, base peak
  # TODO fix all the length
  ms$focus_base_peak <- substring(grep("MS$FOCUSED_ION: BASE_PEAK", record, value = TRUE, fixed = TRUE), 12)

  # MS information, derivative
  # TODO fix all the length
  ms$focus_derivative_form <- substring(grep("MS$FOCUSED_ION: DERIVATIVE_FORM", record, value = TRUE, fixed = TRUE), 12)
  ms$focus_derivative_mass <- substring(grep("MS$FOCUSED_ION: DERIVATIVE_MASS", record, value = TRUE, fixed = TRUE), 12)
  ms$focus_derivative_type <- substring(grep("MS$FOCUSED_ION: DERIVATIVE_TYPE", record, value = TRUE, fixed = TRUE), 12)

  # MS information, precursor
  # TODO fix all the length
  ms$focus_ion_type <- substring(grep("MS$FOCUSED_ION: ION_TYPE", record, value = TRUE, fixed = TRUE), 12)
  ms$focus_precursor_int <- substring(grep("MS$FOCUSED_ION: PRECURSOR_INT", record, value = TRUE, fixed = TRUE), 12)
  ms$focus_precursor_mz <- substring(grep("MS$FOCUSED_ION: PRECURSOR_MZ", record, value = TRUE, fixed = TRUE), 12)
  ms$focus_precursor_type <- substring(grep("MS$FOCUSED_ION: PRECURSOR_TYPE", record, value = TRUE, fixed = TRUE), 12)

  # MS data processing
  # TODO fix all the length
  ms$data_processing_comment <- substring(grep("MS$DATA_PROCESSING: COMMENT", record, value = TRUE, fixed = TRUE), 12)
  ms$data_processing_deprofile <- substring(grep("MS$DATA_PROCESSING: DEPROFILE", record, value = TRUE, fixed = TRUE), 12)
  ms$data_processing_find <- substring(grep("MS$DATA_PROCESSING: FIND_PEAK", record, value = TRUE, fixed = TRUE), 12)
  ms$data_processing_reanalyze <- substring(grep("MS$DATA_PROCESSING: REANALYZE", record, value = TRUE, fixed = TRUE), 12)
  ms$data_processing_recalibrate <- substring(grep("MS$DATA_PROCESSING: RECALIBRATE", record, value = TRUE, fixed = TRUE), 12)
  ms$data_processing_whole <- substring(grep("MS$DATA_PROCESSING: WHOLE", record, value = TRUE, fixed = TRUE), 12)

  # clean up data
  # TODO replace character(0) with NA_character_

  # TODO type conversion for numeric data

  # return result list
  ms
}

#'
#'
#'
.isolate_recordinfo <- function(record) {

  # create empty list
  recordinfo <- list()

  # record information
  recordinfo$accession <- substring(grep("ACCESSION:", record, value = TRUE, fixed = TRUE), 12)
  recordinfo$deprecated <- ""
  recordinfo$record_title <- ""
  recordinfo$date <- ""
  recordinfo$authors <- ""
  recordinfo$license <- ""
  recordinfo$copyright <- ""
  recordinfo$publication <- ""
  recordinfo$project <- "" # multiline

}

#'
#'
#'
.isolate_peakinfo <- function(record) {
  # peak data
  pk_splash <- ""
  pk_num <- ""
}

#'
#'
#'
.isolate_comment <- function(record) {

  # comment section
  comment <- "" # multiline
}

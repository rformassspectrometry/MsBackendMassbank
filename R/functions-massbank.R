.import_massbank_ms_ms_spectrum <- function(x) {
  
  if (!is.character(x) || length(x) != 1)
    stop("'x' has to be of type character with length 1")
  
  con <- file(x)
  mb_record <- readLines(con)
  close(con)

}

.isolate_chemical <- function(record) {
  
  # create empty list
  ch <- list()
  
  # chemical information
  # TODO fix all the length
  ch$name <- as.list(substring(grep("CH$NAME:", record, value = TRUE, fixed = TRUE), 10))
  ch$compound_class <- substring(grep("CH$COMPOUND_CLASS:", record, value = TRUE, fixed = TRUE), 13)
  ch$formula <- substring(grep("CH$FORMULA:", record, value = TRUE, fixed = TRUE), 13)
  ch$exact_mass <- as.numeric(substring(grep("CH$EXACT_MASS:", record, value = TRUE, fixed = TRUE), 16))
  ch$smiles <- substring(grep("CH$SMILES:", record, value = TRUE, fixed = TRUE), 12)
  ch_iupac <- substring(grep("CH$IUPAC:", record, value = TRUE, fixed = TRUE), 11)
  #ch$cdk_depict <- ""
  ch$link_cas <- substring(grep("CH$LINK: CAS", record, value = TRUE, fixed = TRUE), 11)
  ch$link_cayman <- substring(grep("CH$LINK: CAYMAN", record, value = TRUE, fixed = TRUE), 11)
  ch$link_chebi <- substring(grep("CH$LINK: CHEBI", record, value = TRUE, fixed = TRUE), 11)
  ch$link_chembl <- substring(grep("CH$LINK: CHEMBL", record, value = TRUE, fixed = TRUE), 11)
  ch$link_chempdb <- substring(grep("CH$LINK: CHEMPDB", record, value = TRUE, fixed = TRUE), 11)
  ch$link_chemspider <- substring(grep("CH$LINK: CHEMSPIDER", record, value = TRUE, fixed = TRUE), 11)
  ch$link_comptox <- substring(grep("CH$LINK: COMPTOX", record, value = TRUE, fixed = TRUE), 11)
  ch$link_hmdb <- substring(grep("CH$LINK: HMDB", record, value = TRUE, fixed = TRUE), 11)
  ch$link_inchikey <- substring(grep("CH$LINK: INCHIKEY", record, value = TRUE, fixed = TRUE), 11)
  ch$link_kappaview <- substring(grep("CH$LINK: KAPPAVIEW", record, value = TRUE, fixed = TRUE), 11)
  ch$link_kegg <- substring(grep("CH$LINK: KEGG", record, value = TRUE, fixed = TRUE), 11)
  ch$link_knapsack <- substring(grep("CH$LINK: KNAPSACK", record, value = TRUE, fixed = TRUE), 11)
  ch$link_lipidbank <- substring(grep("CH$LINK: LIPIDBANK", record, value = TRUE, fixed = TRUE), 11)
  ch$link_lipidmaps <- substring(grep("CH$LINK: LIPIDMAPS", record, value = TRUE, fixed = TRUE), 11)
  ch$link_nikkaji <- substring(grep("CH$LINK: NIKKAJI", record, value = TRUE, fixed = TRUE), 11)
  ch$link_pubchem <- substring(grep("CH$LINK: PUBCHEM", record, value = TRUE, fixed = TRUE), 11)
  ch$link_zinc <- substring(grep("CH$LINK: ZINC", record, value = TRUE, fixed = TRUE), 11)
  
  # clean up data
  # TODO replace character(0) with NA_character
  
  # return result list
  ch
}

.isolate_species <- function(record) {
  
  # create empty list
  sp <- list()
  
  # species information
  # TODO fix all the length
  sp$scientific_name <- substring(grep("SP$SCIENTIFIC_NAME:", record, value = TRUE, fixed = TRUE), 12)
  sp$lineage <- substring(grep("SP$LINEAGE:", record, value = TRUE, fixed = TRUE), 12)
  sp$link <- substring(grep("SP$LINK:", record, value = TRUE, fixed = TRUE), 12)
  sp$sample <- substring(grep("SP$SAMPLE:", record, value = TRUE, fixed = TRUE), 12)
  
  # clean up data
  # TODO replace character(0) with NA_character
  
  # return result list
  sp
}

.isolate_species <- function(record) {
  
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
  ac$chrom_carrier_gas <- ""
  ac$chrom_column <- ""
  ac$chrom_column_temp <- ""
  ac$chrom_column_temp_gradient <- ""
  ac$chrom_flow_gradient <- ""
  ac$chrom_flow_rate <- ""
  ac$chrom_inj_temp <- ""
  ac$chrom_inj_temp_gradient <- ""
  ac$chrom_rti_kovats <- ""
  ac$chrom_rti_lee <- ""
  ac$chrom_rti_naps <- ""
  ac$chrom_rti_uoa <- ""
  ac$chrom_rti_uoa_pred <- ""
  ac$chrom_rt <- ""
  ac$chrom_rt_uoa_pred <- ""
  ac$chrom_solvent <- "" # list!!!
  ac$chrom_transfer_temp <- ""
  
  # analytical chemistry information, general
  # TODO fix all the length
  ac$general_conc <- ""
  
  # clean up data
  # TODO replace character(0) with NA_character
  
  # return result list
  ac
}
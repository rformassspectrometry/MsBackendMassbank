##' @param f `character(1)` with the path to an MassBank file.
##'
##' @param msLevel `numeric(1)` with the MS level. Default is 2.
##'
##' @param metaDataBlocks `data.frame` data frame indicating which metadata to
##'     be read
##'
##' @param nonStop `logical(1)` whether import should be stopped if an
##'     Massbank file does not contain all required fields. Defaults to
##'     `nonStop = FALSE`.
##'
##' @param ... Additional parameters, currently ignored.
##'
##' @importFrom S4Vectors DataFrame cbind.DataFrame
##' @importFrom IRanges NumericList
##'
##' @author Michael Witting
##'
##' @noRd
.read_massbank <- function(f, msLevel = 2L, metaBlocks = metaDataBlocks(),
                           nonStop = FALSE, ...) {
    requireNamespace("MsBackendMassbank", quietly = TRUE)
    if (length(f) != 1L)
        stop("Please provide a single Massbank file.")
    mb <- scan(file = f, what = "", sep = "\n", quote = "",
               allowEscapes = FALSE, quiet = TRUE)

    begin <- grep("ACCESSION:", mb)
    end <- grep("^//$", mb)
    if (!length(begin) || length(begin) != length(end)) {
        if (nonStop) {
            warning("Unexpected file format: ", basename(f))
            return(DataFrame())
        } else
            stop("Unexpected file format: ", basename(f))
    }

    n <- length(begin)
    spec <- vector("list", length = n)
    ac <- vector("list", length = n)
    ch <- vector("list", length = n)
    sp <- vector("list", length = n)
    ms <- vector("list", length = n)
    record <- vector("list", length = n)
    pk <- vector("list", length = n)
    cmt <- vector("list", length = n)

    for (i in seq(along = spec)) {
        mb_sub <- mb[begin[i]:end[i]]
        spec[[i]] <- MsBackendMassbank:::.extract_mb_spectrum(mb_sub)
        if (metaBlocks$read[which(metaBlocks$metadata == "ac")])
            ac[[i]] <- MsBackendMassbank:::.extract_mb_ac(mb_sub)
        if (metaBlocks$read[which(metaBlocks$metadata == "ch")])
            ch[[i]] <- MsBackendMassbank:::.extract_mb_ch(mb_sub)
        if (metaBlocks$read[which(metaBlocks$metadata == "sp")])
            sp[[i]] <- MsBackendMassbank:::.extract_mb_sp(mb_sub)
        if (metaBlocks$read[which(metaBlocks$metadata == "ms")])
            ms[[i]] <- MsBackendMassbank:::.extract_mb_ms(mb_sub)
        if (metaBlocks$read[which(metaBlocks$metadata == "record")])
            record[[i]] <- MsBackendMassbank:::.extract_mb_record(mb_sub)
        if (metaBlocks$read[which(metaBlocks$metadata == "pk")])
            pk[[i]] <- MsBackendMassbank:::.extract_mb_pk(mb_sub)
        if (metaBlocks$read[which(metaBlocks$metadata == "comment")])
            cmt[[i]] <- list(
                comment = MsBackendMassbank:::.extract_mb_comment(mb_sub))
    }

    res <- DataFrame(do.call(rbind, spec))
    res_ac <- DataFrame(do.call(rbind, ac))
    res_ch <- DataFrame(do.call(rbind, ch))
    res_sp <- DataFrame(do.call(rbind, sp))
    res_ms <- DataFrame(do.call(rbind, ms))
    res_record <- DataFrame(do.call(rbind, record))
    res_pk <- DataFrame(do.call(rbind, pk))
    res_cmt <- DataFrame(do.call(rbind, cmt))

    if (length(res_ac))
        res <- cbind.DataFrame(res, res_ac)
    if (length(res_ch))
        res <- cbind.DataFrame(res, res_ch)
    if (length(res_sp))
        res <- cbind.DataFrame(res, res_sp)
    if (length(res_ms))
        res <- cbind.DataFrame(res, res_ms)
    if (length(res_record))
        res <- cbind.DataFrame(res, res_record)
    if (length(res_pk))
        res <- cbind.DataFrame(res, res_pk)
    if (length(res_cmt))
        res <- cbind.DataFrame(res, res_cmt)
    for (i in seq_along(res)) {
        if (all(lengths(res[[i]]) == 1))
            res[[i]] <- unlist(res[[i]])
    }
    res$mz <- IRanges::NumericList(res$mz, compress = FALSE)
    res$intensity <- IRanges::NumericList(res$intensity, compress = FALSE)
    res$dataOrigin <- f
    res$msLevel <- as.integer(msLevel)
    res
}

##' @param mb `character()` of lines defining a spectrum in mgf
##'     format.
##'
##' @param nonStop `logical(1)` not used.
##'
##' @importFrom utils tail type.convert
##'
##' @author Michael Witting
##'
##' @noRd
.extract_mb_spectrum <- function(mb, nonStop = FALSE) {
    spectrum_start <- grep("PK$PEAK:", mb, fixed = TRUE) + 1
    spectrum_end <- tail(grep("//", mb, fixed = TRUE), 1) - 1

    splitted <- strsplit(mb[spectrum_start:(spectrum_end)]," ")
    spectrum <- matrix(nrow = spectrum_end + 1 - spectrum_start, ncol = 3)

    for (k in seq_along(splitted)) {
        splitted[[k]] <- splitted[[k]][which(splitted[[k]] != "")]
        spectrum[k,] <- splitted[[k]]
    }

    ## convert to data frame and adjust data type
    spectrum <- as.data.frame(spectrum, stringsAsFactors = FALSE)
    spectrum[] <- lapply(spectrum, type.convert, as.is = TRUE)
    colnames(spectrum) <- c("mz", "intensity", "rel.intensity")

    ## isolate Accession, name etc...
    meta <- list()

    meta$accession <-
        substring(grep("ACCESSION:", mb, value = TRUE, fixed = TRUE), 12)
    meta$name <- as.list(substring(grep("CH$NAME:", mb, value = TRUE,
                                        fixed = TRUE), 10))
    meta$smiles <- substring(grep("CH$SMILES:", mb, value = TRUE,
                                  fixed = TRUE), 12)
    meta$exactmass <- as.numeric(
        substring(grep("CH$EXACT_MASS:", mb, value = TRUE, fixed = TRUE), 16))
    meta$formula <- substring(grep("CH$FORMULA:", mb, value = TRUE,
                                   fixed = TRUE), 13)
    meta$inchi <- substring(grep("CH$IUPAC:", mb, value = TRUE,
                                 fixed = TRUE), 11)
    meta$cas <- substring(grep("CH$LINK: CAS", mb, value = TRUE,
                               fixed = TRUE), 14)
    meta$inchikey <- substring(grep("CH$LINK: INCHIKEY", mb, value = TRUE,
                                    fixed = TRUE), 19)
    meta$collisionEnergy <- substring(
        grep("AC$MASS_SPECTROMETRY: COLLISION_ENERGY", mb, value = TRUE,
             fixed = TRUE), 40)
    meta$adduct <- substring(grep("MS$FOCUSED_ION: PRECURSOR_TYPE", mb,
                                  value = TRUE, fixed = TRUE), 32)
    meta$rtime_string <- substring(grep("AC$CHROMATOGRAPHY: RETENTION_TIME",
                                        mb, value = TRUE, fixed = TRUE), 35)
    meta$polarity <- substring(grep("AC$MASS_SPECTROMETRY: ION_MODE", mb,
                                    value = TRUE, fixed = TRUE), 32)
    meta$splash <- substring(grep("PK$SPLASH:", mb, value = TRUE,
                                  fixed = TRUE), 12)
    ## clean NA values
    meta <- .cleanParsing(meta)
    ## type conversion
    ce <- NA_real_
    if (length(meta$collisionEnergy)) {
        ce <- as.numeric(regmatches(
            meta$collisionEnergy, regexpr("[[:digit:]]+\\.*[[:digit:]]*",
                                          meta$collisionEnergy)))
        if (!length(ce))
            ce <- NA_real_
    }
    meta$collisionEnergy <- ce

    ## convert rtime
    if(!is.na(meta$rtime_string)) {
        rtime <- as.numeric(regmatches(
            meta$rtime_string, regexpr("[[:digit:]]+\\.*[[:digit:]]*",
                                       meta$rtime_string)))
        if(grepl("min", meta$rtime_string)) rtime <- rtime * 60
    } else
        rtime <- NA_real_

    ## convert polarity
    if (meta$polarity == "POSITIVE") {
        meta$polarity <- 1L
    } else if(meta$polarity == "NEGATIVE") {
        meta$polarity <- 0L
    } else {
        meta$polarity <- NA_integer_
    }

    precursorMz <- as.numeric(substring(
        grep("MS$FOCUSED_ION: PRECURSOR_M/Z", mb, value = TRUE, fixed = TRUE),
        30))
    precursorIntensity <- as.numeric(substring(
        grep("MS$FOCUSED_ION: PRECURSOR_INT", mb, value = TRUE, fixed = TRUE),
        31))

    ## back up if no values are supplied
    if (!length(precursorMz)) precursorMz <- NA_real_
    if (!length(precursorIntensity)) precursorIntensity <- NA_real_

    title <- substring(grep("RECORD_TITLE:", mb, value = TRUE, fixed = TRUE),
                       15)

    ## first core variables, then others
    list(acquistionNum = 1L,
         centroided = TRUE,
         collisionEnergy = meta$collisionEnergy,
         intensity = spectrum$intensity,
         mz = spectrum$mz,
         polarity = meta$polarity,
         precursorCharge = as.integer(0),
         precursorIntensity = precursorIntensity,
         precursorMz = precursorMz,
         rtime = rtime,
         scanIndex = as.integer(1),
         accession = meta$accession,
         name = meta$name,
         smiles = meta$smiles,
         exactmass = meta$exactmass,
         formula = meta$formula,
         inchi = meta$inchi,
         cas = meta$cas,
         inchikey = meta$inchikey,
         adduct = meta$adduct,
         splash = meta$splash,
         title = title)
}

##' @param mb `character()` of lines defining a spectrum in mgf
##'     format.
##'
##' @author Michael Witting
##'
##' @noRd
.extract_mb_ac <- function(mb) {
    ## create empty list
    ac <- list()
    ## analytical chemistry information, MS instrument -------------------------
    ac$instrument <- substring(grep("AC$INSTRUMENT:", mb, value = TRUE,
                                    fixed = TRUE), 16)
    ac$instrument_type <- substring(grep("AC$INSTRUMENT_TYPE:", mb,
                                         value = TRUE, fixed = TRUE), 21)

    ## analytical chemistry information, MS settings ---------------------------
    ac$ms_ms_type <- substring(grep("AC$MASS_SPECTROMETRY: MS_TYPE", mb,
                                    value = TRUE, fixed = TRUE), 31)
    ac$ms_cap_voltage <- substring(grep(
        "AC$MASS_SPECTROMETRY: CAPILLARY_VOLTAGE", mb, value = TRUE,
        fixed = TRUE), 41)
    ac$ms_col_gas <- substring(grep("AC$MASS_SPECTROMETRY: COLLISION_GAS", mb,
                                    value = TRUE, fixed = TRUE), 37)
    ac$ms_desolv_gas_flow <- substring(
        grep("AC$MASS_SPECTROMETRY: DESOLVATION_GAS_FLOW", mb, value = TRUE,
             fixed = TRUE), 44)
    ac$ms_desolv_temp <- substring(
        grep("AC$MASS_SPECTROMETRY: DESOLVATION_TEMPERATURE", mb, value = TRUE,
             fixed = TRUE), 47)
    ac$ms_frag_mode <- substring(
        grep("AC$MASS_SPECTROMETRY: FRAGMENTATION_MODE", mb, value = TRUE,
             fixed = TRUE), 42)
    ac$ms_ionization <- substring(
        grep("AC$MASS_SPECTROMETRY: IONIZATION", mb, value = TRUE,
             fixed = TRUE), 34)
    ac$ms_ionization_energy <- substring(
        grep("AC$MASS_SPECTROMETRY: IONIZATION_ENERGY", mb, value = TRUE,
             fixed = TRUE), 41)
    ac$ms_laser <- substring(grep("AC$MASS_SPECTROMETRY: LASER", mb,
                                  value = TRUE, fixed = TRUE), 29)
    ac$ms_matrix <- substring(grep("AC$MASS_SPECTROMETRY: MATRIX", mb,
                                   value = TRUE, fixed = TRUE), 30)
    ac$ms_mass_accuracy <- substring(grep("AC$MASS_SPECTROMETRY: MASS_ACCURACY",
                                          mb, value = TRUE, fixed = TRUE), 37)
    ac$ms_mass_range <- substring(grep("AC$MASS_SPECTROMETRY: MASS_RANGE_MZ",
                                       mb, value = TRUE, fixed = TRUE), 37)
    ac$ms_reagent_gas <- substring(grep("AC$MASS_SPECTROMETRY: REAGENT_GAS",
                                        mb, value = TRUE, fixed = TRUE), 35)
    ac$ms_resolution <- substring(grep("AC$MASS_SPECTROMETRY: RESOLUTION",
                                       mb, value = TRUE, fixed = TRUE), 34)
    ac$ms_scan_setting <- substring(
        grep("AC$MASS_SPECTROMETRY: SCANNING_SETTING",
             mb, value = TRUE, fixed = TRUE), 40)
    ac$ms_source_temp <- substring(
        grep("AC$MASS_SPECTROMETRY: SOURCE_TEMPERATURE", mb, value = TRUE,
             fixed = TRUE), 42)
    ac$ms_kinetic_energy <- substring(
      grep("AC$MASS_SPECTROMETRY: KINETIC_ENERGY", mb, value = TRUE,
           fixed = TRUE), 38)
    ac$ms_electron_current <- substring(
      grep("AC$MASS_SPECTROMETRY: ELECTRON_CURRENT", mb, value = TRUE,
           fixed = TRUE), 40)
    ac$ms_reaction_time <- substring(
      grep("AC$MASS_SPECTROMETRY: REACTION_TIME", mb, value = TRUE,
           fixed = TRUE), 37)

    ## analytical chemistry information, chromatography ------------------------
    ac$chrom_carrier_gas <- substring(grep("AC$CHROMATOGRAPHY: CARRIER_GAS",
                                           mb, value = TRUE, fixed = TRUE), 32)
    ac$chrom_column <- substring(grep("AC$CHROMATOGRAPHY: COLUMN_NAME", mb,
                                      value = TRUE, fixed = TRUE), 32)
    ac$chrom_column_temp <- substring(
        grep("AC$CHROMATOGRAPHY: COLUMN_TEMPERATURE", mb, value = TRUE,
             fixed = TRUE), 39)
    ac$chrom_column_temp_gradient <- substring(
        grep("AC$CHROMATOGRAPHY: COLUMN_TEMPERATURE_GRADIENT", mb,
             value = TRUE, fixed = TRUE), 48)
    ac$chrom_flow_gradient <- substring(
        grep("AC$CHROMATOGRAPHY: FLOW_GRADIENT", mb,
             value = TRUE, fixed = TRUE), 34)
    ac$chrom_flow_rate <- substring(grep("AC$CHROMATOGRAPHY: FLOW_RATE",
                                         mb, value = TRUE, fixed = TRUE), 30)
    ac$chrom_inj_temp <- substring(
        grep("AC$CHROMATOGRAPHY: INJECTION_TEMPERATURE", mb, value = TRUE,
             fixed = TRUE), 42)
    ac$chrom_inj_temp_gradient <- substring(
        grep("AC$CHROMATOGRAPHY: INJECTION_TEMPERATURE_GRADIENT", mb,
             value = TRUE, fixed = TRUE), 51)
    ac$chrom_rti_kovats <- substring(grep("AC$CHROMATOGRAPHY: KOVATS_RTI",
                                          mb, value = TRUE, fixed = TRUE), 31)
    ac$chrom_rti_lee <- substring(grep("AC$CHROMATOGRAPHY: LEE_RTI",
                                       mb, value = TRUE, fixed = TRUE), 28)
    ac$chrom_rti_naps <- substring(grep("AC$CHROMATOGRAPHY: NAPS_RTI",
                                        mb, value = TRUE, fixed = TRUE), 29)
    ac$chrom_rti_uoa <- substring(grep("AC$CHROMATOGRAPHY: UOA_RTI",
                                       mb, value = TRUE, fixed = TRUE), 28)
    ac$chrom_rti_uoa_pred <- substring(
        grep("AC$CHROMATOGRAPHY: UOA_PREDICTED_RTI", mb,
             value = TRUE, fixed = TRUE), 38)
    ac$chrom_rt <- substring(grep("AC$CHROMATOGRAPHY: RETENTION_TIME",
                                  mb, value = TRUE, fixed = TRUE), 35)
    ac$chrom_rt_uoa_pred <- substring(
        grep("AC$CHROMATOGRAPHY: TRAMS_PREDICTED_RETENTION_TIME",
             mb, value = TRUE, fixed = TRUE), 51)
    ac$chrom_solvent <- as.list(substring(
        grep("AC$CHROMATOGRAPHY: SOLVENT", mb, value = TRUE, fixed = TRUE), 28))
    ac$chrom_transfer_temp <- substring(
        grep("AC$CHROMATOGRAPHY: TRANSFERLINE_TEMPERATURE", mb,
             value = TRUE, fixed = TRUE), 45)

    ## analytical chemistry information, ion mobility
    ## preparation for IMS update of MassBank format
    ac$ims_instrument_type <- substring(
        grep("AC$ION_MOBILITY: INSTRUMENT_TYPE", mb,
             value = TRUE, fixed = TRUE), 34)
    ac$ims_drift_gas <- substring(grep("AC$ION_MOBILITY: DRIFT_GAS",
                                       mb, value = TRUE, fixed = TRUE), 28)
    ac$ims_drift_time <- substring(grep("AC$ION_MOBILITY: DRIFT_TIME",
                                        mb, value = TRUE, fixed = TRUE), 29)
    ac$ims_ccs <- substring(grep("AC$ION_MOBILITY: CCS", mb, value = TRUE,
                                 fixed = TRUE), 22)

    ## analytical chemistry information, general -------------------------------
    ac$general_conc <- substring(grep("AC$GENERAL: CONCENTRATION",
                                      mb, value = TRUE, fixed = TRUE), 27)

    ac <- .cleanParsing(ac)
    ac
}

##' @param mb `character()` of lines defining a spectrum in mgf
##'     format.
##'
##' @author Michael Witting
##'
##' @noRd
.extract_mb_ch <- function(mb) {
    ch <- list()

    ## isolate chemical information
    ch$compound_class <- substring(grep("CH$COMPOUND_CLASS:", mb,
                                        value = TRUE, fixed = TRUE), 20)
    ch$link_cayman <- substring(grep("CH$LINK: CAYMAN", mb, value = TRUE,
                                     fixed = TRUE), 17)
    ch$link_chebi <- substring(grep("CH$LINK: CHEBI", mb, value = TRUE,
                                    fixed = TRUE), 16)
    ch$link_chembl <- substring(grep("CH$LINK: CHEMBL", mb, value = TRUE,
                                     fixed = TRUE), 17)
    ch$link_chempdb <- substring(grep("CH$LINK: CHEMPDB", mb, value = TRUE,
                                      fixed = TRUE), 18)
    ch$link_chemspider <- substring(grep("CH$LINK: CHEMSPIDER", mb,
                                         value = TRUE, fixed = TRUE), 21)
    ch$link_comptox <- substring(grep("CH$LINK: COMPTOX", mb, value = TRUE,
                                      fixed = TRUE), 18)
    ch$link_hmdb <- substring(grep("CH$LINK: HMDB", mb, value = TRUE,
                                   fixed = TRUE), 15)
    ch$link_kappaview <- substring(grep("CH$LINK: KAPPAVIEW", mb,
                                        value = TRUE, fixed = TRUE), 20)
    ch$link_kegg <- substring(grep("CH$LINK: KEGG", mb, value = TRUE,
                                   fixed = TRUE), 15)
    ch$link_knapsack <- substring(grep("CH$LINK: KNAPSACK", mb,
                                       value = TRUE, fixed = TRUE), 19)
    ch$link_lipidbank <- substring(grep("CH$LINK: LIPIDBANK", mb,
                                        value = TRUE, fixed = TRUE), 20)
    ch$link_lipidmaps <- substring(grep("CH$LINK: LIPIDMAPS", mb,
                                        value = TRUE, fixed = TRUE), 20)
    ch$link_nikkaji <- substring(grep("CH$LINK: NIKKAJI", mb,
                                      value = TRUE, fixed = TRUE), 18)
    ch$link_pubchem <- substring(grep("CH$LINK: PUBCHEM", mb,
                                      value = TRUE, fixed = TRUE), 18)
    ch$link_zinc <- substring(grep("CH$LINK: ZINC", mb,
                                   value = TRUE, fixed = TRUE), 15)
    ch <- .cleanParsing(ch)
    ch
}

##' @param mb `character()` of lines defining a spectrum in mgf
##'     format.
##'
##' @author Michael Witting
##'
##' @noRd
.extract_mb_sp <- function(mb) {
    sp <- list()
    ## species information
    sp$scientific_name <- substring(grep("SP$SCIENTIFIC_NAME:", mb,
                                         value = TRUE, fixed = TRUE), 21)
    sp$lineage <- substring(grep("SP$LINEAGE:", mb, value = TRUE,
                                 fixed = TRUE), 13)
    sp$link <- substring(grep("SP$LINK:", mb, value = TRUE,
                              fixed = TRUE), 10)
    sp$sample <- substring(grep("SP$SAMPLE:", mb, value = TRUE,
                                fixed = TRUE), 12)
    sp <- .cleanParsing(sp)
    sp
}

##' @param mb `character()` of lines defining a spectrum in mgf
##'     format.
##'
##' @author Michael Witting
##'
##' @noRd
.extract_mb_ms <- function(mb) {
    ms <- list()
    ## MS information, base peak
    ms$focus_base_peak <- substring(grep("MS$FOCUSED_ION: BASE_PEAK", mb,
                                         value = TRUE, fixed = TRUE), 27)
    ## MS information, derivative
    ms$focus_derivative_form <- substring(
        grep("MS$FOCUSED_ION: DERIVATIVE_FORM", mb, value = TRUE,
             fixed = TRUE), 33)
    ms$focus_derivative_mass <- substring(
        grep("MS$FOCUSED_ION: DERIVATIVE_MASS", mb, value = TRUE,
             fixed = TRUE), 33)
    ms$focus_derivative_type <- substring(
        grep("MS$FOCUSED_ION: DERIVATIVE_TYPE", mb, value = TRUE,
             fixed = TRUE), 33)
    ## MS information, precursor
    ms$focus_ion_type <- substring(grep("MS$FOCUSED_ION: ION_TYPE", mb,
                                        value = TRUE, fixed = TRUE), 26)

    ## MS data processing
    ms$data_processing_comment <- substring(
        grep("MS$DATA_PROCESSING: COMMENT", mb, value = TRUE, fixed = TRUE), 29)
    ms$data_processing_deprofile <- substring(
        grep("MS$DATA_PROCESSING: DEPROFILE", mb, value = TRUE, fixed = TRUE),
        31)
    ms$data_processing_find <- substring(
        grep("MS$DATA_PROCESSING: FIND_PEAK", mb, value = TRUE, fixed = TRUE),
        31)
    ms$data_processing_reanalyze <- substring(
        grep("MS$DATA_PROCESSING: REANALYZE", mb, value = TRUE, fixed = TRUE),
        31)
    ms$data_processing_recalibrate <- substring(
        grep("MS$DATA_PROCESSING: RECALIBRATE", mb, value = TRUE, fixed = TRUE),
        33)
    ms$data_processing_whole <- substring(
        grep("MS$DATA_PROCESSING: WHOLE", mb, value = TRUE, fixed = TRUE), 27)
    ms <- .cleanParsing(ms)
    ## TODO type conversion for numeric data
    ms
}

##' @param mb `character()` of lines defining a spectrum in mgf
##'     format.
##'
##' @author Michael Witting
##'
##' @noRd
.extract_mb_record <- function(mb) {
    recordinfo <- list()
    ## mb information
    recordinfo$deprecated <- substring(grep("DEPRECATED:", mb, value = TRUE,
                                            fixed = TRUE), 13)
    recordinfo$date <- substring(grep("DATE:", mb, value = TRUE,
                                      fixed = TRUE), 7)
    recordinfo$authors <- substring(grep("AUTHORS:", mb, value = TRUE,
                                         fixed = TRUE), 10)
    recordinfo$license <- substring(grep("LICENSE:", mb, value = TRUE,
                                         fixed = TRUE), 10)
    recordinfo$copyright <- substring(grep("COPYRIGHT:", mb, value = TRUE,
                                           fixed = TRUE), 12)
    recordinfo$publication <- substring(grep("PUBLICATION:", mb, value = TRUE,
                                             fixed = TRUE), 14)
    recordinfo$project <- as.list(substring(grep("PROJECT:", mb, value = TRUE,
                                                 fixed = TRUE), 10))
    recordinfo <- .cleanParsing(recordinfo)
    ## TODO type conversion for dates
    recordinfo
}

##' @param mb `character()` of lines defining a spectrum in mgf
##'     format.
##'
##' @author Michael Witting
##'
##' @noRd
.extract_mb_pk <- function(mb) {
    pk <- list()
    pk$pknum <- as.integer(substring(grep("PK$NUM_PEAK:", mb, value = TRUE,
                                          fixed = TRUE), 14))
    pk <- .cleanParsing(pk)
    pk
}

##' @param mb `character()` of lines defining a spectrum in mgf
##'     format.
##'
##' @author Michael Witting
##'
##' @noRd
.extract_mb_comment <- function(mb) {
    comment <- substring(grep("COMMENT:", mb, value = TRUE, fixed = TRUE), 10)
    comment <- .cleanParsing(comment)
    comment
}

##' @title Metadata blocks to be read
##'
##' @description
##'
##' `metaDataBlocks` returns a `data.frame` with the MassBank metadata blocks
##' and whether they should be imported by default from the MassBank text files.
##'
##' @return A `data.frame` with metadata blocks.
##'
##' @author Michael Witting
##'
##' @importFrom utils read.csv
##'
##' @export
##'
##' @examples
##'
##' metaDataBlocks()
metaDataBlocks <- function() {
    read.csv(dir(system.file("extdata", package = "MsBackendMassbank"),
                 pattern = "metadata_blocks.csv",
                 full.names = TRUE),
             header = TRUE, as.is = TRUE,
             stringsAsFactors = FALSE)
}

##' Clean parsing
##'
##' @title Cleaning meta data parsing
##'
##' @param x `List` with entries
##'
##' @return `List` with cleaned entries
##'
##' @noRd
.cleanParsing <- function(x) {
    x[lengths(x) == 0] <- NA_character_
    x
}

#' @description
#'
#' Function to export a `Spectra` object in MassBank format to `con`.
#'
#' @param x `Spectra`
#'
#' @param con output file.
#'
#' @param mapping named `character` vector that maps from `spectraVariables`
#'    (i.e. `names(mapping)`) to the variable name that should be used in the
#'    MGF file.
#'
#' @author Michael Witting
#'
#' @importMethodsFrom Spectra spectraVariables spectraNames spectraData
#'
#' @noRd
.export_massbank <- function(x, con = stdout(),
                             mapping = spectraVariableMapping(
                                 MsBackendMassbank())) {
    if (is.character(con) && file.exists(con)) {
        message("Overwriting ", con, "!")
        unlink(con)
    }
    if (is.character(con)) {
        con <- file(description = con, open = "at")
        on.exit(close(con))
    }
    .cat <- function(..., file = con, sep = " ", append = TRUE) {
        cat(..., file = file, sep = sep, append = append)
    }
    ## iterate over all spectra
    for (i in seq_along(x)) {
        spv <- spectraVariables(x[i])
        spd <- spectraData(x[i], spv[!(spv %in%
                                       c("dataOrigin", "dataStorage"))])
        idx <- match(colnames(spd), names(mapping))
        colnames(spd)[!is.na(idx)] <- mapping[idx[!is.na(idx)]]
        spp <- peaksData(x[i])
        ## here list with stuff in right order
        entries <- .getEntries()
        for (entry in entries) {
            if (entry %in% colnames(spd)) {
                value <- spd[entry][[1]]
                if (entry == "AC$MASS_SPECTROMETRY: ION_MODE") {
                    if (value == 0L) {
                        value <- "NEGATIVE"
                    } else if (value == 1L) {
                        value <- "POSITIVE"
                    }
                } else if (entry == "AC$CHROMATOGRAPHY: RETENTION_TIME") {
                    if (!is.na(value))
                        value <- paste0(value / 60, " min")
                }
                if (!is.na(value)) {
                    if (is.list(value)) {
                        .cat(paste0(entry, " ", unlist(value), collapse = "\n"))
                        .cat("\n")
                    } else {
                        .cat(entry, paste0(value, "\n"))
                    }
                }
            }
        }
        .cat("PK$PEAK: m/z int. rel.int.\n")
        .cat(paste0("  ",
                    peaksData(x[i])[[1]][,1],
                    " ",
                    peaksData(x[i])[[1]][,2],
                    " ",
                    as.integer(peaksData(x[i])[[1]][,2] /
                               max(peaksData(x[i])[[1]][,2]) * 999),
                    collapse = "\n"))
        .cat("\n//\n")
    }
}

.getEntries <- function() {
    c(
        "ACCESSION:",
        "DEPRECATED:",
        "RECORD_TITLE:",
        "DATE:",
        "AUTHORS:",
        "LICENSE:",
        "COPYRIGHT:",
        "PUBLICATION:",
        "PROJECT:",
        "COMMENT:",
        "CH$NAME:",
        "CH$COMPOUND_CLASS:",
        "CH$FORMULA:",
        "CH$EXACT_MASS:",
        "CH$SMILES:",
        "CH$IUPAC:",
        "CH$CDK_DEPICT:",
        "CH$LINK: CAS",
        "CH$LINK: CAYMAN",
        "CH$LINK: CHEBI",
        "CH$LINK: CHEMBL",
        "CH$LINK: CHEMPDB",
        "CH$LINK: CHEMSPIDER",
        "CH$LINK: COMPTOX",
        "CH$LINK: HMDB",
        "CH$LINK: INCHIKEY",
        "CH$LINK: KAPPAVIEW",
        "CH$LINK: KEGG",
        "CH$LINK: KNAPSACK",
        "CH$LINK: LIPIDBANK",
        "CH$LINK: LIPIDMAPS",
        "CH$LINK: NIKKAJI",
        "CH$LINK: PUBCHEM",
        "CH$LINK: ZINC",
        "SP$SCIENTIFIC_NAME:",
        "SP$LINEAGE:",
        "SP$LINK:",
        "SP$SAMPLE:",
        "AC$INSTRUMENT:",
        "AC$INSTRUMENT_TYPE:",
        "AC$MASS_SPECTROMETRY: MS_TYPE",
        "AC$MASS_SPECTROMETRY: ION_MODE",
        "AC$MASS_SPECTROMETRY: CAPILLARY_VOLTAGE",
        "AC$MASS_SPECTROMETRY: COLLISION_ENERGY",
        "AC$MASS_SPECTROMETRY: COLLISION_GAS",
        "AC$MASS_SPECTROMETRY: DESOLVATION_GAS_FLOW",
        "AC$MASS_SPECTROMETRY: DESOLVATION_TEMPERATURE",
        "AC$MASS_SPECTROMETRY: FRAGMENTATION_MODE",
        "AC$MASS_SPECTROMETRY: IONIZATION",
        "AC$MASS_SPECTROMETRY: IONIZATION_ENERGY",
        "AC$MASS_SPECTROMETRY: LASER",
        "AC$MASS_SPECTROMETRY: MATRIX",
        "AC$MASS_SPECTROMETRY: MASS_ACCURACY",
        "AC$MASS_SPECTROMETRY: MASS_RANGE_MZ",
        "AC$MASS_SPECTROMETRY: REAGENT_GAS",
        "AC$MASS_SPECTROMETRY: RESOLUTION",
        "AC$MASS_SPECTROMETRY: SCANNING_SETTING",
        "AC$MASS_SPECTROMETRY: SOURCE_TEMPERATURE",
        "AC$MASS_SPECTROMETRY: KINETIC_ENERGY",
        "AC$MASS_SPECTROMETRY: ELECTRON_CURRENT",
        "AC$MASS_SPECTROMETRY: REACTION_TIME",
        "AC$CHROMATOGRAPHY: CARRIER_GAS",
        "AC$CHROMATOGRAPHY: COLUMN_NAME",
        "AC$CHROMATOGRAPHY: COLUMN_TEMPERATURE",
        "AC$CHROMATOGRAPHY: COLUMN_TEMPERATURE_GRADIENT",
        "AC$CHROMATOGRAPHY: FLOW_GRADIENT",
        "AC$CHROMATOGRAPHY: FLOW_RATE",
        "AC$CHROMATOGRAPHY: INJECTION_TEMPERATURE",
        "AC$CHROMATOGRAPHY: INJECTION_TEMPERATURE_GRADIENT",
        "AC$CHROMATOGRAPHY: INLET_TYPE",
        "AC$CHROMATOGRAPHY: KOVATS_RTI",
        "AC$CHROMATOGRAPHY: LEE_RTI",
        "AC$CHROMATOGRAPHY: NAPS_RTI",
        "AC$CHROMATOGRAPHY: UOA_RTI",
        "AC$CHROMATOGRAPHY: UOA_PREDICTED_RTI",
        "AC$CHROMATOGRAPHY: RETENTION_TIME",
        "AC$CHROMATOGRAPHY: UOA_PREDICTED_RETENTION_TIME",
        "AC$CHROMATOGRAPHY: SOLVENT",
        "AC$CHROMATOGRAPHY: TRANSFERLINE_TEMPERATURE",
        "AC$ION_MOBILITY: INSTRUMENT_TYPE",
        "AC$ION_MOBILITY: DRIFT_GAS",
        "AC$ION_MOBILITY: DRIFT_TIME",
        "AC$ION_MOBILITY: CCS",
        "MS$FOCUSED_ION: BASE_PEAK",
        "MS$FOCUSED_ION: DERIVATIVE_FORM",
        "MS$FOCUSED_ION: DERIVATIVE_MASS",
        "MS$FOCUSED_ION: DERIVATIVE_TYPE",
        "MS$FOCUSED_ION: ION_TYPE",
        "MS$FOCUSED_ION: PRECURSOR_INT",
        "MS$FOCUSED_ION: PRECURSOR_M/Z",
        "MS$FOCUSED_ION: PRECURSOR_TYPE",
        "MS$DATA_PROCESSING: COMMENT",
        "MS$DATA_PROCESSING: DEPROFILE",
        "MS$DATA_PROCESSING: FIND_PEAK",
        "MS$DATA_PROCESSING: REANALYZE",
        "MS$DATA_PROCESSING: RECALIBRATE",
        "MS$DATA_PROCESSING: WHOLE",
        "PK$SPLASH:",
        "PK$ANNOTATION:",
        "PK$NUM_PEAK:")
}

#' @rdname MsBackendMassbankSql
#'
#' @export MsBackendMassbankSql
MsBackendMassbankSql <- function() {
    if (!requireNamespace("DBI", quietly = TRUE))
        stop("'MsBackendMassbankSql' requires package 'DBI'. Please ",
             "install with 'install.packages(\"DBI\")'")
    new("MsBackendMassbankSql")
}

#' @importFrom DBI dbListTables
#'
#' @noRd
.valid_dbcon <- function(x) {
    if (length(x)) {
        if (!inherits(x, "DBIConnection"))
            return("'dbcon' is expected to be a connection to a database")
        tables <- dbListTables(x)
        if (!all(c("msms_spectrum", "msms_spectrum_peak") %in% tables))
            return(paste0("Database lacks some required tables. Is 'dbcon' a",
                          " connection to a MassBank database?"))
    }
    NULL
}

.valid_local_data <- function(x, y) {
    if (nrow(x) && nrow(x) != length(y))
        "Number of rows in local data and number of spectra don't match"
    else NULL
}

#' Returns the spectra data, from the database and eventually filling with
#' *core* spectra variables, if they are not available in the database.
#'
#' The data can be either:
#' - in the database.
#' - in the local data (if new variables were added with $name <-).
#' - core spectra variables - if they are not in the database they have to be
#'   initialized with `NA` and the correct data type.
#'
#' @return a `data.frame` - always, even if only with a single column.
#'
#' @importFrom IRanges NumericList CharacterList
#'
#' @importFrom S4Vectors extractCOLS
#'
#' @importFrom methods as
#'
#' @author Johannes Rainer
#'
#' @noRd
.spectra_data_massbank_sql <- function(x, columns = spectraVariables(x)) {
    local_cols <- intersect(columns, colnames(x@localData))
    db_cols <- intersect(x@spectraVariables, columns)
    db_cols <- db_cols[!db_cols %in% local_cols]
    mz_cols <- intersect(columns, c("mz", "intensity"))
    core_cols <- intersect(columns, names(Spectra:::.SPECTRA_DATA_COLUMNS))
    core_cols <- core_cols[!core_cols %in% c(db_cols, mz_cols, local_cols)]
    res <- NULL
    ## Get data from database
    if (length(db_cols)) {
        res <- DataFrame(.fetch_spectra_data_sql(x, columns = db_cols))
        if (any(colnames(res) == "synonym"))
            res$synonym <- CharacterList(res$synonym, compress = FALSE)
    }
    ## Get m/z and intensity values
    if (length(mz_cols)) {
        pks <- .fetch_peaks_sql(x, columns = mz_cols)
        f <- factor(pks$spectrum_id)
        if (any(mz_cols == "mz")) {
            mzs <- split(pks$mz, f)
            if (length(res))
                res$mz <- NumericList(mzs[as.character(x@spectraIds)],
                                      compress = FALSE)
            else res <- DataFrame(
                     mz = NumericList(mzs[as.character(x@spectraIds)],
                                      compress = FALSE))
        }
        if (any(mz_cols == "intensity")) {
            ints <- split(pks$intensity, f)
            if (length(res))
                res$intensity <- NumericList(
                    ints[as.character(x@spectraIds)], compress = FALSE)
            else res <- DataFrame(
                     intensity = NumericList(
                         ints[as.character(x@spectraIds)], compress = FALSE))
        }
    }
    ## Get local data
    if (length(local_cols)) {
        if (length(res))
            res <- cbind(res, extractCOLS(x@localData, local_cols))
        else res <- extractCOLS(x@localData, local_cols)
    }
    ## Create missing core variables
    if (length(core_cols)) {
        tmp <- DataFrame(lapply(Spectra:::.SPECTRA_DATA_COLUMNS[core_cols],
                                function(z, n) rep(as(NA, z), n), length(x)))
        if (length(res))
            res <- cbind(res, tmp)
        else res <- tmp
        if (any(core_cols == "dataStorage"))
            res$dataStorage <- dataStorage(x)
    }
    if (!all(columns %in% colnames(res)))
        stop("Column(s) ", paste0(columns[!columns %in% names(res)],
                                  collapse = ", "), " not available.",
             call. = FALSE)
    extractCOLS(res, columns)
}

## .spectra_data_massbank_sql <- function(x, columns = spectraVariables(x)) {
##     local_cols <- intersect(columns, colnames(x@localData))
##     db_cols <- intersect(x@spectraVariables, columns)
##     db_cols <- db_cols[!db_cols %in% local_cols]
##     mz_cols <- intersect(columns, c("mz", "intensity"))
##     core_cols <- intersect(columns, names(Spectra:::.SPECTRA_DATA_COLUMNS))
##     core_cols <- core_cols[!core_cols %in% c(db_cols, mz_cols, local_cols)]
##     res <- NULL
##     ## Get data from database
##     if (length(db_cols)) {
##         res <- DataFrame(.fetch_spectra_data_sql(x, columns = db_cols))
##         if (any(colnames(res) == "synonym"))
##             res$synonym <- CharacterList(res$synonym, compress = FALSE)
##     }
##     ## Get m/z and intensity values
##     if (length(mz_cols)) {
##         pks <- .fetch_peaks_sql(x, columns = mz_cols)
##         f <- factor(pks$spectrum_id)
##         if (any(mz_cols == "mz")) {
##             mzs <- split(pks$mz, f)
##             if (length(res))
##                 res$mz <- NumericList(mzs[as.character(x@spectraIds)],
##                                       compress = FALSE)
##             else res <- DataFrame(
##                      mz = NumericList(mzs[as.character(x@spectraIds)],
##                                       compress = FALSE))
##         }
##         if (any(mz_cols == "intensity")) {
##             ints <- split(pks$intensity, f)
##             if (length(res))
##                 res$intensity <- NumericList(
##                     ints[as.character(x@spectraIds)], compress = FALSE)
##             else res <- DataFrame(
##                      intensity = NumericList(
##                          ints[as.character(x@spectraIds)], compress = FALSE))
##         }
##     }
##     ## Get local data
##     if (length(local_cols)) {
##         if (length(res))
##             res <- cbind(res, extractCOLS(x@localData, local_cols))
##         else res <- extractCOLS(x@localData, local_cols)
##     }
##     ## Create missing core variables
##     if (length(core_cols)) {
##         tmp <- DataFrame(lapply(Spectra:::.SPECTRA_DATA_COLUMNS[core_cols],
##                                 function(z, n) rep(as(NA, z), n), length(x)))
##         if (length(res))
##             res <- cbind(res, tmp)
##         else res <- tmp
##         if (any(core_cols == "dataStorage"))
##             res$dataStorage <- dataStorage(x)
##     }
##     if (!all(columns %in% colnames(res)))
##         stop("Column(s) ", paste0(columns[!columns %in% names(res)],
##                                   collapse = ", "), " not available.",
##              call. = FALSE)
##     extractCOLS(res, columns)
## }

#' @importFrom DBI dbSendQuery dbBind dbFetch dbClearResult
#'
#' @noRd
.fetch_peaks_sql <- function(x, columns = c("mz", "intensity")) {
    if (length(x@dbcon)) {
        dbGetQuery(
            x@dbcon,
            paste0("select spectrum_id,", paste(columns, collapse = ","),
                   " from msms_spectrum_peak where spectrum_id in (",
                   paste0("'", unique(x@spectraIds), "'", collapse = ","),")"))
    } else {
        data.frame(spectrum_id = character(), mz = numeric(),
                   intensity = numeric())
    }
}

## .fetch_peaks_sql <- function(x, columns = c("mz", "intensity")) {
##     if (length(x@dbcon)) {
##         qry <- dbSendQuery(
##             x@dbcon, paste0("select spectrum_id,",
##                             paste(columns, collapse = ","),
##                             " from msms_spectrum_peak where spectrum_id = ?"))
##         qry <- dbBind(qry, list(unique(x@spectraIds)))
##         res <- dbFetch(qry)
##         dbClearResult(qry)
##         res
##     } else {
##         data.frame(spectrum_id = character(), mz = numeric(),
##                    intensity = numeric())
##     }
## }

.columns_sql <- c(
    precursorIntensity = "precursor_intensity",
    precursorMz = "precursor_mz_text",
    msLevel = "ms_level",
    compound_id = "msms_spectrum.compound_id"
)

.map_spectraVariables_to_sql <- function(x) {
    for (i in seq_along(.columns_sql))
        x <- sub(names(.columns_sql)[i], .columns_sql[i], x, fixed = TRUE)
    x
}

.map_sql_to_spectraVariables <- function(x) {
    for (i in seq_along(.columns_sql))
        x <- sub(.columns_sql[i], names(.columns_sql[i]), x, fixed = TRUE)
    x
}

#' Simple helper that creates a join query depending on the provided columns.
#'
#' @param x `Spectra`.
#'
#' @param columns `character` with the column names.
#'
#' @noRd
.join_query <- function(x, columns) {
    res <- "msms_spectrum"
    if (any(columns %in% x@.tables$ms_compound))
        res <- paste0(res, " join ms_compound on (msms_spectrum.compound_id",
                      "=ms_compound.compound_id)")
    res
}

.fetch_spectra_data_sql <- function(x, columns = c("spectrum_id")) {
    orig_columns <- columns
    if (any(columns %in% c("compound_name", "synonym"))) {
        columns <- columns[!columns %in% c("compound_name", "synonym")]
        columns <- unique(c(columns, "compound_id"))
    }
    sql_columns <-
        unique(c("spectrum_id", .map_spectraVariables_to_sql(columns)))
    ## That turns out to be faster than dbBind if we use a field in the
    ## database that is unique (such as spectrum_id).
    res <- dbGetQuery(
        x@dbcon,
        paste0("select ", paste(sql_columns, collapse = ","), " from ",
               .join_query(x, sql_columns), " where spectrum_id in (",
               paste0("'", unique(x@spectraIds), "'", collapse = ", ") ,")"))
    idx <- match(x@spectraIds, res$spectrum_id)
    res <- res[idx[!is.na(idx)], , drop = FALSE]
    rownames(res) <- NULL
    if (any(columns == "msLevel")) {
        res$msLevel <- as.integer(sub("MS", "", res$ms_level))
        res$ms_level <- NULL
    }
    if (any(columns == "polarity")) {
        pol <- rep(NA_integer_, nrow(res))
        pol[res$polarity == "POSITIVE"] <- 1L
        pol[res$polarity == "NEGATIVE"] <- 0L
        res$polarity <- pol
    }
    if (any(columns == "publication"))
        res$dataOrigin <- res$publication
    if (any(columns == "precursorIntensity")) {
        res$precursorIntensity <- as.numeric(res$precursor_intensity)
        res$precursor_intensity <- NULL
    }
    ## So far we're not dealing with multiple precursor m/z here!
    if (any(columns == "precursorMz")) {
        suppressWarnings(
            res$precursorMz <- as.numeric(res$precursor_mz_text))
        if (!any(columns == "precursor_mz_text"))
            res$precursor_intensity_text <- NULL
    }
    if (any(columns == "collisionEnergy")) {
        suppressWarnings(
            res$collisionEnergy <- as.numeric(res$collision_energy_test))
        res$collision_energy_text <- NULL
    }
    ## manage synonym and compound_name. Need a second query for that.
    if (any(orig_columns %in% c("synonym", "compound_name"))) {
        cmps <- dbGetQuery(
            x@dbcon,
            paste0("select * from synonym where compound_id in (",
                   paste0("'", unique(res$compound_id), "'",
                          collapse = ","), ")"))
        cmpl <- split(
            cmps$synonym,
            as.factor(cmps$compound_id))[as.character(res$compound_id)]
        if (any(orig_columns == "compound_name"))
            res$compound_name <- vapply(cmpl, function(z) z[1], character(1))
        if (any(orig_columns == "synonym"))
            res$synonym <- cmpl
    }
    res[, orig_columns, drop = FALSE]
}

.compounds_sql <- function(x, id, columns = "*") {
    id <- force(id)
    qry <- dbSendQuery(
        x, paste0("select ", columns, " from ms_compound join synonym on (",
                  "ms_compound.compound_id=synonym.compound_id)",
                  " where ms_compound.compound_id = ?"))
    qry <- dbBind(qry, list(id))
    res <- dbFetch(qry)
    dbClearResult(qry)
    idx <- grep("^compound_id", colnames(res))
    if (length(idx) > 1)
        res <- res[, -idx[-1]]
    res
}

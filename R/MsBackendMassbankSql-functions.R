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
        if (!all(c("spectra_data", "peaks") %in% tables))
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
#' @importFrom IRanges NumericList
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
    }
    ## Get m/z and intensity values
    if (length(mz_cols)) {
        pks <- .fetch_peaks_sql(x, columns = mz_cols)
        f <- factor(pks$spectrum_id, levels = x@spectraIds)
        if (any(mz_cols == "mz")) {
            if (length(res))
                res$mz <- NumericList(split(pks$mz, f), compress = FALSE)
            else res <- DataFrame(mz = NumericList(split(pks$mz, f),
                                                   compress = FALSE))
        }
        if (any(mz_cols == "intensity")) {
            if (length(res))
                res$intensity <- NumericList(split(pks$intensity, f),
                                             compress = FALSE)
            else res <- DataFrame(intensity = NumericList(
                                      split(pks$intensity, f),
                                      compress = FALSE))
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

#' @importFrom DBI dbSendQuery dbBind dbFetch dbClearResult
#'
#' @noRd
.fetch_peaks_sql <- function(x, columns = c("mz", "intensity")) {
    if (length(x@dbcon)) {
        qry <- dbSendQuery(x@dbcon, paste0("select spectrum_id,",
                                           paste(columns, collapse = ","),
                                           " from peaks where spectrum_id = ?"))
        qry <- dbBind(qry, list(x@spectraIds))
        res <- dbFetch(qry)
        dbClearResult(qry)
        res
    } else {
        data.frame(spectrum_id = character(), mz = numeric(),
                   intensity = numeric())
    }
}

.fetch_spectra_data_sql <- function(x, columns = c("spectrum_id")) {
    qry <- dbSendQuery(
        x@dbcon, paste0("select ", paste(columns, collapse = ","),
                        " from spectra_data where spectrum_id = ?"))
    qry <- dbBind(qry, list(x@spectraIds))
    res <- dbFetch(qry)
    dbClearResult(qry)
    if (any(columns == "msLevel"))
        res$msLevel <- as.integer(sub("MS", "", res$msLevel))
    if (any(columns == "polarity")) {
        pol <- rep(NA_integer_, nrow(res))
        pol[res$polarity == "POSITIVE"] <- 1L
        pol[res$polarity == "NEGATIVE"] <- 0L
        res$polarity <- pol
    }
    if (any(columns == "publication"))
        res$dataOrigin <- res$publication
    if (any(columns == "precursorIntensity"))
        res$precursorIntensity <- as.numeric(res$precursorIntensity)
    ## So far we're not dealing with multiple precursor m/z here!
    if (any(columns == "precursorMz"))
        suppressWarnings(
            res$precursorMz <- as.numeric(res$precursorMz))
    res
}

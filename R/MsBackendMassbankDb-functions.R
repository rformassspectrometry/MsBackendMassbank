#' @rdname MsBackendMassbankDb
#'
#' @export MsBackendMassbankDb
MsBackendMassbankDb <- function() {
    if (!requireNamespace("DBI", quietly = TRUE))
        stop("'MsBackendMassbankDb' requires package 'DBI'. Please ",
             "install with 'install.packages(\"DBI\")'")
    new("MsBackendMassbankDb")
}

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

#' Returns the spectra data, from the database and eventually filling with
#' *core* spectra variables, if they are not available in the database.
#'
#' The data can be either:
#' - in the database.
#' - in the local data (if new variables were added with $name <-).
#'
#' @return a `data.frame` - always, even if only with a single column.
#'
#' @author Johannes Rainer
#'
#' @noRd
.spectra_data_massbank <- function(x, columns = spectraVariables(x)) {
    db_cols <- intersect(x@spectraVariables, columns)
    core_cols <- intersect(columns, names(Spectra:::.SPECTRA_DATA_COLUMNS))
    core_cols <- core_cols[!core_cols %in% c(db_cols, "mz", "intensity")]
    res <- data.frame()
    if (length(db_cols)) {
        qry <- dbSendQuery(
            x@dbcon, paste0("select ", paste(db_cols, collapse = ","),
                            " from spectra_data where spectrum_id = ?"))
        qry <- dbBind(qry, list(x@spectraIds))
        res <- dbFetch(qry)
        dbClearResult(qry)
    }
    if (any(columns %in% c("mz", "intensity"))) {
    }
    if (length(core_cols)) {
        tmp <- as.data.frame(
            lapply(Spectra:::.SPECTRA_DATA_COLUMNS[core_cols],
                   function(z, n) rep(as(NA, z), n), n = length(x)))
        if (nrow(res))
            res <- cbind(res, tmp)
        else res <- tmp
        if (any(core_cols == "dataStorage"))
            res$dataStorage <- dataStorage(x)
    }
    res[, columns, drop = FALSE]
}

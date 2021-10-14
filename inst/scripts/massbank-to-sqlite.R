#' @title Create an SQLite MassBank database for Spectra backends
#'
#' @description
#'
#' Extract all data relevant from the MySQL MassBank database and store that
#' into a (minimal) SQLite database to enable use with the
#' `MsBackendMassbankSql` database.
#'
#' @param con connection to the MySQL MassBank database.
#'
#' @param sqlite `SQLiteConnection` with the connection to the SQLite database.
#'
#' @author Johannes Rainer
massbank_to_sqlite <- function(con, sqlite = dbConnect(SQLite(), tempfile())) {
    tbls <- c("ms_compound", "msms_spectrum", "msms_spectrum_peak", "synonym")
    if (!all(tbls %in% dbListTables(con)))
        stop("One or more required tables not found.")
    for (tbl in tbls) {
        tmp <- dbGetQuery(con, paste0("select * from ", tbl))
        dbWriteTable(sqlite, name = tbl, tmp)
    }
    dbExecute(sqlite, "create index compound_id_idx on ms_compound (compound_id)")
    dbExecute(sqlite, "create index msms_id_idx on msms_spectrum_peak (spectrum_id)")
    dbExecute(sqlite, "create index msms_mid_idx on msms_spectrum (spectrum_id)")
    dbExecute(sqlite, "create index msms_cid_idx on msms_spectrum (compound_id)")
}

#' @title MS backend accessing the MassBank MySQL database
#'
#' @aliases MsBackendMassbankSql-class compounds
#'
#' @description
#'
#' The `MsBackendMassbankSql` provides access to mass spectrometry data from
#' [MassBank](https://massbank.eu/MassBank/) by directly accessing its
#' MySQL/MariaDb database. In addition it supports adding new spectra variables
#' or *locally* changing spectra variables provided by MassBank (without
#' changing the original values in the database).
#'
#' Note that `MsBackendMassbankSql` requires a local installation of the
#' MassBank database since direct database access is not supported for the
#' *main* MassBank instance.
#'
#' Also, some of the fields in the MassBank database are not directly compatible
#' with `Spectra`, such as the *collision energy* which is not available as a
#' numeric value. The collision energy as available in MassBank is reported as
#' spectra variable `"collision_energy_text"`. Also, precursor m/z values
#' reported for some spectra can not be converted to a `numeric` and hence `NA`
#' is reported with the spectra variable `precursorMz` for these spectra. The
#' variable `"precursor_mz_text"` can be used to get the *original* precursor
#' m/z reported in MassBank.
#'
#' @param dbcon For `backendInitialize,MsBackendMassbankSql`: SQL database
#'     connection to the MassBank (MariaDb) database.
#'
#' @param columns For `spectraData` accessor: optional `character` with column
#'     names (spectra variables) that should be included in the
#'     returned `DataFrame`. By default, all columns are returned.
#'     For `peaksData` accessor: optional `character` with requested columns in
#'     the individual `matrix` of the returned `list`. Use
#'     `peaksVariables(object)` for supported columns.
#'
#' @param drop For `[`: not considered.
#'
#' @param initial For `tic`: `logical(1)` whether the initially
#'     reported total ion current should be reported, or whether the
#'     total ion current should be (re)calculated on the actual data
#'     (`initial = FALSE`).
#'
#' @param i For `[`: `integer`, `logical` or `character` to subset the object.
#'
#' @param j For `[`: not supported.
#'
#' @param name name of the variable to replace for `<-` methods. See individual
#'     method description or expected data type.
#'
#' @param object Object extending `MsBackendMassbankSql`.
#'
#' @param spectraVariables For `selectSpectraVariables`: `character` with the
#'     names of the spectra variables to which the backend should be subsetted.
#'
#' @param value replacement value for `<-` methods. See individual
#'     method description or expected data type.
#'
#' @param x Object extending `MsBackendMassbankSql`.
#'
#' @param ... Additional arguments.
#'
#'
#' @section Supported Backend functions:
#'
#' The following functions are supported by the `MsBackendMassbankSqlMassbankDb`.
#'
#' - `[`: subset the backend. Only subsetting by element (*row*/`i`) is
#'   allowed
#'
#' - `$`, `$<-`: access or set/add a single spectrum variable (column) in the
#'   backend.
#'
#' - `acquisitionNum`: returns the acquisition number of each
#'   spectrum. Returns an `integer` of length equal to the number of
#'   spectra (with `NA_integer_` if not available).
#'
#' - `peaksData` returns a `list` with the spectras' peak data. The length of
#'   the list is equal to the number of spectra in `object`. Each element of
#'   the list is a `matrix` with columns `"mz"` and `"intensity"`. For an empty
#'   spectrum, a `matrix` with 0 rows and two columns (named `mz` and
#'   `intensity`) is returned. Parameter `columns` allows to select which peaks
#'   variables to return, but supports currently only `"mz"` and `"intensity"`.
#'
#' - `backendInitialize`: initialises the backend by retrieving the IDs of all
#'   spectra in the database. Parameter `dbcon` with the connection to the
#'   MassBank MySQL database is required.
#'
#' - `dataOrigin`: gets a `character` of length equal to the number of spectra
#'   in `object` with the *data origin* of each spectrum. This could e.g. be
#'   the mzML file from which the data was read.
#'
#' - `dataStorage`: returns `"<MassBank>"` for all spectra.
#'
#' - `centroided`, `centroided<-`: gets or sets the centroiding
#'   information of the spectra. `centroided` returns a `logical`
#'   vector of length equal to the number of spectra with `TRUE` if a
#'   spectrum is centroided, `FALSE` if it is in profile mode and `NA`
#'   if it is undefined. See also `isCentroided` for estimating from
#'   the spectrum data whether the spectrum is centroided.  `value`
#'   for `centroided<-` is either a single `logical` or a `logical` of
#'   length equal to the number of spectra in `object`.
#'
#' - `collisionEnergy`, `collisionEnergy<-`: gets or sets the
#'   collision energy for all spectra in `object`. `collisionEnergy`
#'   returns a `numeric` with length equal to the number of spectra
#'   (`NA_real_` if not present/defined), `collisionEnergy<-` takes a
#'   `numeric` of length equal to the number of spectra in `object`. Note that
#'   the collision energy description from MassBank are provided as spectra
#'   variable `"collisionEnergyText"`.
#'
#' - `intensity`: gets the intensity values from the spectra. Returns
#'   a [NumericList()] of `numeric` vectors (intensity values for each
#'   spectrum). The length of the `list` is equal to the number of
#'   `spectra` in `object`.
#'
#' - `ionCount`: returns a `numeric` with the sum of intensities for
#'   each spectrum. If the spectrum is empty (see `isEmpty`),
#'   `NA_real_` is returned.
#'
#' - `isCentroided`: a heuristic approach assessing if the spectra in
#'   `object` are in profile or centroided mode. The function takes
#'   the `qtl` th quantile top peaks, then calculates the difference
#'   between adjacent m/z value and returns `TRUE` if the first
#'   quartile is greater than `k`. (See `Spectra:::.isCentroided` for
#'   the code.)
#'
#' - `isEmpty`: checks whether a spectrum in `object` is empty
#'   (i.e. does not contain any peaks). Returns a `logical` vector of
#'   length equal number of spectra.
#'
#' - `isolationWindowLowerMz`, `isolationWindowLowerMz<-`: gets or sets the
#'   lower m/z boundary of the isolation window.
#'
#' - `isolationWindowTargetMz`, `isolationWindowTargetMz<-`: gets or sets the
#'   target m/z of the isolation window.
#'
#' - `isolationWindowUpperMz`, `isolationWindowUpperMz<-`: gets or sets the
#'   upper m/z boundary of the isolation window.
#'
#' - `isReadOnly`: returns a `logical(1)` whether the backend is *read
#'   only* or does allow also to write/update data.
#'
#' - `length`: returns the number of spectra in the object.
#'
#' - `lengths`: gets the number of peaks (m/z-intensity values) per
#'   spectrum.  Returns an `integer` vector (length equal to the
#'   number of spectra). For empty spectra, `0` is returned.
#'
#' - `msLevel`: gets the spectra's MS level. Returns an `integer`
#'   vector (of length equal to the number of spectra) with the MS
#'   level for each spectrum (or `NA_integer_` if not available).
#'
#' - `mz`: gets the mass-to-charge ratios (m/z) from the
#'   spectra. Returns a [NumericList()] or length equal to the number of
#'   spectra, each element a `numeric` vector with the m/z values of
#'   one spectrum.
#'
#' - `polarity`, `polarity<-`: gets or sets the polarity for each
#'   spectrum.  `polarity` returns an `integer` vector (length equal
#'   to the number of spectra), with `0` and `1` representing negative
#'   and positive polarities, respectively. `polarity<-` expects an
#'   integer vector of length 1 or equal to the number of spectra.
#'
#' - `precursorCharge`, `precursorIntensity`, `precursorMz`,
#'   `precScanNum`, `precAcquisitionNum`: get the charge (`integer`),
#'   intensity (`numeric`), m/z (`numeric`), scan index (`integer`)
#'   and acquisition number (`interger`) of the precursor for MS level
#'   2 and above spectra from the object. Returns a vector of length equal to
#'   the number of spectra in `object`. `NA` are reported for MS1
#'   spectra of if no precursor information is available.
#'
#' - `reset`: restores the backend to its original state, i.e. deletes all
#'   locally modified data and reinitializes the backend to the full data
#'   available in the database.
#'
#' - `rtime`, `rtime<-`: gets or sets the retention times for each
#'   spectrum (in seconds). `rtime` returns a `numeric` vector (length equal to
#'   the number of spectra) with the retention time for each spectrum.
#'   `rtime<-` expects a numeric vector with length equal to the
#'   number of spectra.
#'
#' - `scanIndex`: returns an `integer` vector with the *scan index*
#'   for each spectrum. This represents the relative index of the
#'   spectrum within each file. Note that this can be different to the
#'   `acquisitionNum` of the spectrum which is the index of the
#'   spectrum as reported in the mzML file.
#'
#' - `selectSpectraVariables`: reduces the information within the backend to
#'   the selected spectra variables.
#'
#' - `smoothed`,`smoothed<-`: gets or sets whether a spectrum is
#'   *smoothed*. `smoothed` returns a `logical` vector of length equal
#'   to the number of spectra. `smoothed<-` takes a `logical` vector
#'   of length 1 or equal to the number of spectra in `object`.
#'
#' - `spectraData`: gets general spectrum metadata (annotation, also called
#'   header).  `spectraData` returns a `DataFrame`. Note that replacing the
#'   spectra data with `spectraData<-` is not supported.
#'
#' - `spectraNames`: returns a `character` vector with the names of
#'   the spectra in `object`.
#'
#' - `spectraVariables`: returns a `character` vector with the
#'   available spectra variables (columns, fields or attributes)
#'   available in `object`. This should return **all** spectra variables which
#'   are present in `object`, also `"mz"` and `"intensity"` (which are by
#'   default not returned by the `spectraVariables,Spectra` method).
#'
#' - `tic`: gets the total ion current/count (sum of signal of a
#'   spectrum) for all spectra in `object`. By default, the value
#'   reported in the original raw data file is returned. For an empty
#'   spectrum, `NA_real_` is returned.
#'
#' @section Not supported Backend functions:
#'
#' The following functions are not supported by the `MsBackendMassbankSql` since
#' the original data can not be changed.
#'
#' `backendMerge`, `export`, `filterDataStorage`, `filterPrecursorScan`,
#' `peaksData<-`, `filterAcquisitionNum`, `intensity<-`, `mz<-`, `precScanNum`,
#' `spectraData<-`, `spectraNames<-`.
#'
#' @section Retrieving compound annotations for spectra:
#'
#' While compound annotations are also provided *via* the `spectraVariables` of
#' the backend, it would also be possible to use the `compounds` function on
#' a `Spectra` object (that uses a `MsBackendMassbankSql` backend) to retrieve
#' compound annotations for the specific spectra.
#'
#' @name MsBackendMassbankSql
#'
#' @return See documentation of respective function.
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @exportClass MsBackendMassbankSql
#'
#' @examples
#'
#' ## Create a connection to a database with MassBank data - in the present
#' ## example we connect to a tiny SQLite database bundled in this package
#' ## as public access to the MassBank MySQL is not (yet) supported. See the
#' ## vignette for more information on how to install MassBank locally and
#' ## enable MySQL database connections
#' library(RSQLite)
#' con <- dbConnect(SQLite(), system.file("sql", "minimassbank.sqlite",
#'     package = "MsBackendMassbank"))
#'
#' ## Given that we have the connection to a MassBank databas we can
#' ## initialize the backend:
#' be <- backendInitialize(MsBackendMassbankSql(), dbcon = con)
#' be
#'
#' ## Access MS level
#' msLevel(be)
#' be$msLevel
#'
#' ## Access m/z values
#' be$mz
#'
#' ## Access the full spectra data (including m/z and intensity values)
#' spectraData(be)
#'
#' ## Add a new spectra variable
#' be$new_variable <- "b"
#' be$new_variable
#'
#' ## Subset the backend
#' be_sub <- be[c(3, 1)]
#'
#' spectraNames(be)
#' spectraNames(be_sub)
NULL

setClassUnion("DBIConnectionOrNULL", c("DBIConnection", "NULL"))

#' @importClassesFrom DBI DBIConnection
#'
#' @importClassesFrom Spectra MsBackendCached
#'
#' @importClassesFrom S4Vectors DataFrame
setClass(
    "MsBackendMassbankSql",
    contains = "MsBackendCached",
    slots = c(
        dbcon = "DBIConnectionOrNULL",
        spectraIds = "character",
        .tables = "list"),
    prototype = prototype(
        dbcon = NULL,
        spectraIds = character(),
        .tables = list(),
        readonly = TRUE, version = "0.2"))

#' @importFrom methods .valueClassTest is new validObject
#'
#' @noRd
setValidity("MsBackendMassbankSql", function(object) {
    msg <- .valid_dbcon(object@dbcon)
    if (is.null(msg)) TRUE
    else msg
})

#' @exportMethod backendInitialize
#'
#' @importFrom DBI dbGetQuery
#'
#' @rdname MsBackendMassbankSql
setMethod("backendInitialize", "MsBackendMassbankSql", function(object,
                                                                dbcon, ...) {
    if (missing(dbcon))
        stop("Parameter 'dbcon' is required for 'MsBackendMassbankSql'")
    msg <- .valid_dbcon(dbcon)
    object@dbcon <- dbcon
    if (length(msg))
        stop(msg)

    res <- dbGetQuery(
        dbcon, "select spectrum_id, precursor_mz_text from msms_spectrum")
    object@spectraIds <- as.character(res[, "spectrum_id"])
    object@.tables <- list(
        msms_spectrum = colnames(
            dbGetQuery(dbcon, "select * from msms_spectrum limit 0")),
        ms_compound = colnames(
            dbGetQuery(dbcon, "select * from ms_compound limit 0")),
        synonym = colnames(
            dbGetQuery(dbcon, "select * from synonym")))
    ## Initialize cached backend
    object <- callNextMethod(
        object, nspectra = length(object@spectraIds),
        spectraVariables = c(
            .map_sql_to_spectraVariables(unique(unlist(object@.tables))),
            "precursor_mz_text", "compound_name"))
    suppressWarnings(object@localData$precursorMz <-
                         as.numeric(res[, "precursor_mz_text"]))
    validObject(object)
    object
})

#' @importMethodsFrom Spectra peaksData
#'
#' @importMethodsFrom Spectra peaksVariables
#'
#' @exportMethod peaksData
#'
#' @rdname MsBackendMassbankSql
setMethod(
    "peaksData", "MsBackendMassbankSql",
    function(object, columns = peaksVariables(object)) {
        if (!all(columns %in% c("mz", "intensity")))
            stop("'peaksData' for 'MsBackendMassbankSql' does only support",
                 " columns \"mz\" and \"intensity\"", call. = FALSE)
        pks <- .fetch_peaks_sql(object, columns = columns)
        f <- factor(pks$spectrum_id)        # using levels does not work because we can have duplicates
        pks <- unname(split.data.frame(pks, f)[as.character(object@spectraIds)])
        idx <- seq_along(columns) + 1
        lapply(pks, function(z) {
            if (nrow(z))
                as.matrix(z[, idx, drop = FALSE], rownames.force = FALSE)
            else matrix(NA_real_, ncol = length(columns), nrow = 0,
                        dimnames = list(character(), columns))
        })
    })

#' @exportMethod dataStorage
#'
#' @importMethodsFrom ProtGenerics dataStorage
#'
#' @rdname MsBackendMassbankSql
setMethod("dataStorage", "MsBackendMassbankSql", function(object) {
    rep("<MassBank>", length(object))
})

#' @exportMethod intensity<-
#'
#' @importMethodsFrom ProtGenerics intensity<-
#'
#' @rdname MsBackendMassbankSql
setReplaceMethod("intensity", "MsBackendMassbankSql", function(object, value) {
    stop("Can not replace original intensity values in MassBank.")
})

#' @exportMethod mz<-
#'
#' @importMethodsFrom ProtGenerics mz<-
#'
#' @rdname MsBackendMassbankSql
setReplaceMethod("mz", "MsBackendMassbankSql", function(object, value) {
    stop("Can not replace original data in MassBank.")
})

#' @exportMethod reset
#'
#' @importMethodsFrom Spectra reset
#'
#' @rdname MsBackendMassbankSql
setMethod("reset", "MsBackendMassbankSql", function(object) {
    message("Restoring original data ...", appendLF = FALSE)
    if (is(object@dbcon, "DBIConnection"))
        object <- backendInitialize(object, object@dbcon)
    message("DONE")
    object
})

#' @exportMethod spectraData
#'
#' @rdname MsBackendMassbankSql
setMethod(
    "spectraData", "MsBackendMassbankSql",
    function(object, columns = spectraVariables(object)) {
        .spectra_data_massbank_sql(object, columns = columns)
    })

#' @exportMethod spectraNames
#'
#' @importMethodsFrom ProtGenerics spectraNames
#'
#' @rdname MsBackendMassbankSql
setMethod("spectraNames", "MsBackendMassbankSql", function(object) {
    object@spectraIds
})

#' @exportMethod spectraNames<-
#'
#' @importMethodsFrom ProtGenerics spectraNames<-
#'
#' @rdname MsBackendMassbankSql
setReplaceMethod("spectraNames", "MsBackendMassbankSql",
                 function(object, value) {
                     stop(class(object)[1],
                          " does not support replacing spectra names (IDs).")
})

#' @exportMethod tic
#'
#' @importMethodsFrom ProtGenerics tic
#'
#' @importFrom Spectra intensity
#'
#' @importFrom MsCoreUtils vapply1d
#'
#' @rdname MsBackendMassbankSql
setMethod("tic", "MsBackendMassbankSql", function(object, initial = TRUE) {
    if (initial) {
        if (any(colnames(object@localData) == "totIonCurrent"))
            object@localData[, "totIonCurrent"]
        else rep(NA_real_, times = length(object))
    } else vapply1d(intensity(object), sum, na.rm = TRUE)
})

#' @exportMethod [
#'
#' @importFrom MsCoreUtils i2index
#'
#' @importFrom methods slot<-
#'
#' @importFrom S4Vectors extractROWS
#'
#' @rdname MsBackendMassbankSql
setMethod("[", "MsBackendMassbankSql", function(x, i, j, ..., drop = FALSE) {
    if (missing(i))
        return(x)
    i <- i2index(i, length(x), x@spectraIds)
    slot(x, "spectraIds", check = FALSE) <- x@spectraIds[i]
    x <- callNextMethod(x, i = i)
    x
})

#' @rdname MsBackendMassbankSql
#'
#' @importFrom ProtGenerics compounds
#'
#' @export
setMethod("compounds", "Spectra", function(object, ...) {
    compounds(object@backend, ...)
})

#' @rdname MsBackendMassbankSql
#'
#' @importFrom IRanges CharacterList
#'
#' @export
setMethod("compounds", "MsBackendMassbankSql", function(object, ...) {
    if (!length(object))
        return(DataFrame())
    if (!any(spectraVariables(object) == "compound_id"))
        stop("Spectra variable 'compound_id' not present.")
    res <- .compounds_sql(object@dbcon, object$compound_id)
    syns <- res$synonym
    syns <- split(
        syns, f = factor(res$compound_id, levels = unique(res$compound_id)))
    res <- DataFrame(unique(res[, colnames(res) != "synonym"]))
    res$synonym <- CharacterList(syns, compress = FALSE)
    res$name <- vapply(syns, function(z) z[1], character(1))
    res[match(object$compound_id, res$compound_id), ]
})

#' @rdname MsBackendMassbankSql
#'
#' @export
setReplaceMethod("$", "MsBackendMassbankSql", function(x, name, value) {
    if (name %in% c("spectrum_id"))
        stop("Spectra IDs can not be changes.", call. = FALSE)
    callNextMethod()
})


#' @rdname MsBackendMassbankSql
#'
#' @importMethodsFrom Spectra precScanNum
#' @export
setMethod("precScanNum", "MsBackendMassbankSql", function(object) {
    message("precursor scan numbers not available")
    rep(NA_integer_, length(object))
})

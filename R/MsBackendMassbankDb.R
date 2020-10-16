#' @title MS backend accessing the MassBank MySQL database
#'
#' @aliases MsBackendMassbankDb-class
#'
#' @description
#'
#' The `MsBackendMassbankDb` provides access to mass spectrometry data from
#' [MassBank](https://massbank.eu/MassBank/) by directly accessing its
#' MySQL/MariaDb database.
#'
#' Note that this requires a local installation or the MassBank database since
#' direct database access is not supported for the *main* MassBank instance.
#'
#' @param acquisitionNum for `filterPrecursorScan`: `integer` with the
#'     acquisition number of the spectra to which the object should be
#'     subsetted.
#'
#' @param columns For `spectraData` accessor: optional `character` with column
#'     names (spectra variables) that should be included in the
#'     returned `DataFrame`. By default, all columns are returned.
#'
#' @param data For `backendInitialize`: `DataFrame` with spectrum
#'     metadata/data. This parameter can be empty for `MsBackendMzR` backends
#'     but needs to be provided for `MsBackendDataFrame` backends.
#'
#' @param dataOrigin For `filterDataOrigin`: `character` to define which
#'     spectra to keep.
#'     For `filterAcquisitionNum`: optionally specify if filtering should occurr
#'     only for spectra of selected `dataOrigin`.
#'
#' @param drop For `[`: not considered.
#'
#' @param f `factor` defining the grouping to split `x`. See [split()].
#'
#' @param file For `filterFile`: index or name of the file(s) to which the data
#'     should be subsetted. For `export`: `character` of length 1 or equal to
#'     the number of spectra.
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
#' @param msLevel `integer` defining the MS level of the spectra to which the
#'     function should be applied. For `filterMsLevel`: the MS level to which
#'     `object` should be subsetted.
#'
#' @param mz For `filterIsolationWindow`: `numeric(1)` with the m/z value to
#'     filter the object. For `filterPrecursorMz`: `numeric(2)` with the lower
#'     and upper m/z boundary.
#'
#' @param n for `filterAcquisitionNum`: `integer` with the acquisition numbers
#'     to filter for.
#'
#' @param name For `$` and `$<-`: the name of the spectra variable to return
#'     or set.
#'
#' @param object Object extending `MsBackendMassbankDb`.
#'
#' @param polarity For `filterPolarity`: `integer` specifying the polarity to
#'     to subset `object`.
#'
#' @param rt for `filterRt`: `numeric(2)` defining the retention time range to
#'     be used to subset/filter `object`.
#'
#' @param spectraVariables For `selectSpectraVariables`: `character` with the
#'     names of the spectra variables to which the backend should be subsetted.
#'
#' @param use.names For `lengths`: whether spectrum names should be used.
#'
#' @param value replacement value for `<-` methods. See individual
#'     method description or expected data type.
#'
#' @param x Object extending `MsBackendMassbankDb`.
#'
#' @param ... Additional arguments.
#'
#'
#' @section Supported Backend functions:
#'
#' The following functions are supported by the `MsBackendMassbankDbMassbankDb`.
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
#'   `intensity`) is returned.
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
#' - `dropNaSpectraVariables`: removes spectra variables (i.e. columns in the
#'   object's `spectraData` that contain only missing values (`NA`). Note that
#'   while columns with only `NA`s are removed, a `spectraData` call after
#'   `dropNaSpectraVariables` might still show columns containing `NA` values
#'   for *core* spectra variables.
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
#'   `numeric` of length equal to the number of spectra in `object`.
#'
#' - `filterAcquisitionNum`: filters the object keeping only spectra matching
#'   the provided acquisition numbers (argument `n`). If `dataOrigin` or
#'   `dataStorage` is also provided, `object` is subsetted to the spectra with
#'   an acquisition number equal to `n` **in spectra with matching dataOrigin
#'   or dataStorage values** retaining all other spectra.
#'
#' - `filterDataOrigin`: filters the object retaining spectra matching the
#'   provided `dataOrigin`. Parameter `dataOrigin` has to be of type
#'   `character` and needs to match exactly the data origin value of the
#'   spectra to subset.
#'   `filterDataOrigin` should return the data ordered by the provided
#'   `dataOrigin` parameter, i.e. if `dataOrigin = c("2", "1")` was provided,
#'   the spectra in the resulting object should be ordered accordingly (first
#'   spectra from data origin `"2"` and then from `"1"`).
#'
#' - `filterEmptySpectra`: removes empty spectra (i.e. spectra without peaks).
#'
#' - `filterFile`: retains data of files matching the file index or file name
#'    provided with parameter `file`.
#'
#' - `filterIsolationWindow`: retains spectra that contain `mz` in their
#'   isolation window m/z range (i.e. with an `isolationWindowLowerMz` `<=` `mz`
#'   and `isolationWindowUpperMz` `>=` `mz`.
#'
#' - `filterMsLevel`: retains spectra of MS level `msLevel`.
#'
#' - `filterPolarity`: retains spectra of polarity `polarity`.
#'
#' - `filterPrecursorMz`: retains spectra with a precursor m/z within the
#'   provided m/z range.
#'
#' - `filterPrecursorScan`: retains parent (e.g. MS1) and children scans (e.g.
#'    MS2) of acquisition number `acquisitionNum`.
#'
#' - `filterRt`: retains spectra of MS level `msLevel` with retention times
#'    within (`>=`) `rt[1]` and (`<=`) `rt[2]`.
#'
#' - `intensity`: gets the intensity values from the spectra. Returns
#'   a [NumericList()] of `numeric` vectors (intensity values for each
#'   spectrum). The length of the `list` is equal to the number of
#'   `spectra` in `object`.
#'
#' - `intensity<-`: replaces the intensity values. `value` has to be a `list`
#'   (or [NumericList()]) of length equal to the number of spectra and the
#'   number of values within each list element identical to the number of
#'   peaks in each spectrum (i.e. the `lengths(x)`). Note that just
#'   writeable backends support this method.
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
#' - `mz<-`: replaces the m/z values. `value` has to be a `list` of length equal
#'   to the number of spectra and the number of values within each list element
#'   identical to the number of peaks in each spectrum (i.e. the
#'   `lengths(x)`). Note that just writeable backends support this method.
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
#' - `peaksData<-` replaces the peak data (m/z and intensity values) of the
#'   backend. This method expects a `list` of `matrix` objects with columns
#'   `"mz"` and `"intensity"` that has the same length as the number of
#'   spectra in the backend. Note that just writeable backends support this
#'   method.
#'
#' - `reset` a backend (if supported). This method will be called on the backend
#'   by the `reset,Spectra` method that is supposed to restore the data to its
#'   original state (see `reset,Spectra` for more details). The function
#'   returns the *reset* backend. The default implementation for `MsBackendMassbankDb`
#'   returns the backend as-is.
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
#' - `spectraData`, `spectraData<-`: gets or sets general spectrum
#'   metadata (annotation, also called header).  `spectraData` returns
#'   a `DataFrame`, `spectraData<-` expects a `DataFrame` with the same number
#'   of rows as there are spectra in `object`. Note that `spectraData` has to
#'   return the full data, i.e. also the m/z and intensity values (as a `list`
#'   or `SimpleList` in columns `"mz"` and `"intensity"`.
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
#' - `split`: splits the backend into a `list` of backends (depending on
#'   parameter `f`). The default method for `MsBackendMassbankDb` uses [split.default()],
#'   thus backends extending `MsBackendMassbankDb` don't necessarily need to implement
#'   this method.
#'
#' - `tic`: gets the total ion current/count (sum of signal of a
#'   spectrum) for all spectra in `object`. By default, the value
#'   reported in the original raw data file is returned. For an empty
#'   spectrum, `NA_real_` is returned.
#'
#' @section Not supported Backend functions:
#'
#' The following functions are not supported by the `MsBackendMassbankDb` since
#' the original data can not be changed.
#'
#' `backendMerge`, `export`, `filterDataStorage`.
#'
#' @name MsBackendMassbankDb
#'
#' @return See documentation of respective function.
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @exportClass MsBackendMassbankDb
NULL

setClassUnion("DBIConnectionOrNULL", c("DBIConnection", "NULL"))

#' @importClassesFrom DBI DBIConnection
#'
#' @importClassesFrom S4Vectors DataFrame
setClass(
    "MsBackendMassbankDb",
    contains = "MsBackend",
    slots = c(
        dbcon = "DBIConnectionOrNULL",
        spectraIds = "character",
        spectraVariables = "character",
        localData = "DataFrame"),
    prototype = prototype(
        dbcon = NULL,
        spectraIds = character(),
        spectraVariables = character(),
        localData = DataFrame(),
        readonly = TRUE, version = "0.1"))

#' @importFrom methods .valueClassTest is new validObject
#'
#' @noRd
setValidity("MsBackendMassbankDb", function(object) {
    msg <- .valid_dbcon(object@dbcon)
    if (is.null(msg)) TRUE
    else msg
})

#' @exportMethod backendInitialize
#'
#' @importFrom DBI dbGetQuery
#'
#' @rdname MsBackendMassbankDb
setMethod("backendInitialize", "MsBackendMassbankDb", function(object, dbcon,
                                                               ...) {
    if (missing(dbcon))
        stop("Parameter 'dbcon' is required for 'MsBackendMassbankDb'")
    msg <- .valid_dbcon(dbcon)
    object@dbcon <- dbcon
    if (length(msg))
        stop(msg)
    res <- dbGetQuery(dbcon, "select spectrum_id from spectra_data")
    object@spectraIds <- res[, 1]
    res <- dbGetQuery(dbcon, "select * from spectra_data limit 1")
    object@spectraVariables <- colnames(res)
    validObject(object)
    object
})

## #' @exportMethod acquisitionNum
## #'
## #' @importMethodsFrom ProtGenerics acquisitionNum
## #'
## #' @rdname MsBackendMassbankDb
## setMethod("acquisitionNum", "MsBackendMassbankDb", function(object) {
##     stop("Not implemented for ", class(object), ".")
## })

## #' @exportMethod peaksData
## #'
## #' @rdname MsBackendMassbankDb
## setMethod("peaksData", "MsBackendMassbankDb", function(object) {
##     stop("Not implemented for ", class(object), ".")
## })

## #' @exportMethod centroided
## #'
## #' @aliases centroided<-,MsBackendMassbankDb-method
## #'
## #' @importMethodsFrom ProtGenerics centroided
## #'
## #' @rdname MsBackendMassbankDb
## setMethod("centroided", "MsBackendMassbankDb", function(object) {
##     stop("Not implemented for ", class(object), ".")
## })

## #' @exportMethod centroided<-
## #'
## #' @importMethodsFrom ProtGenerics centroided<-
## #'
## #' @rdname MsBackendMassbankDb
## setReplaceMethod("centroided", "MsBackendMassbankDb", function(object, value) {
##     stop("Not implemented for ", class(object), ".")
## })

## #' @exportMethod collisionEnergy
## #'
## #' @importMethodsFrom ProtGenerics collisionEnergy
## #'
## #' @rdname MsBackendMassbankDb
## setMethod("collisionEnergy", "MsBackendMassbankDb", function(object) {
##     stop("Not implemented for ", class(object), ".")
## })

## #' @exportMethod collisionEnergy<-
## #'
## #' @importMethodsFrom ProtGenerics collisionEnergy<-
## #'
## #' @rdname MsBackendMassbankDb
## setReplaceMethod("collisionEnergy", "MsBackendMassbankDb", function(object, value) {
##     stop("Not implemented for ", class(object), ".")
## })

## #' @exportMethod dataOrigin
## #'
## #' @importMethodsFrom ProtGenerics dataOrigin
## #'
## #' @rdname MsBackendMassbankDb
## setMethod("dataOrigin", "MsBackendMassbankDb", function(object) {
##     stop("Not implemented for ", class(object), ".")
## })

## #' @exportMethod dataOrigin<-
## #'
## #' @importMethodsFrom ProtGenerics dataOrigin<-
## #'
## #' @rdname MsBackendMassbankDb
## setReplaceMethod("dataOrigin", "MsBackendMassbankDb", function(object, value) {
##     stop("Not implemented for ", class(object), ".")
## })

#' @exportMethod dataStorage
#'
#' @importMethodsFrom ProtGenerics dataStorage
#'
#' @rdname MsBackendMassbankDb
setMethod("dataStorage", "MsBackendMassbankDb", function(object) {
    rep("<MassBank>", length(object))
})

## #' @exportMethod dataStorage<-
## #'
## #' @importMethodsFrom ProtGenerics dataStorage<-
## #'
## #' @rdname MsBackendMassbankDb
## setReplaceMethod("dataStorage", "MsBackendMassbankDb", function(object, value) {
##     stop("Method 'dataStorage' is not implemented for ", class(object), ".")
## })

## #' @exportMethod dropNaSpectraVariables
## #'
## #' @rdname MsBackendMassbankDb
## setMethod("dropNaSpectraVariables", "MsBackendMassbankDb", function(object) {
##     svs <- spectraVariables(object)
##     svs <- svs[!(svs %in% c("mz", "intensity"))]
##     spd <- spectraData(object, columns = svs)
##     keep <- !vapply1l(spd, function(z) all(is.na(z)))
##     selectSpectraVariables(object, c(svs[keep], "mz", "intensity"))
## })

## #' @exportMethod filterAcquisitionNum
## #'
## #' @importMethodsFrom ProtGenerics filterAcquisitionNum
## #'
## #' @rdname MsBackendMassbankDb
## setMethod("filterAcquisitionNum", "MsBackendMassbankDb", function(object, n, file, ...) {
##     stop("Not implemented for ", class(object), ".")
## })

## #' @exportMethod filterDataOrigin
## #'
## #' @importMethodsFrom ProtGenerics filterDataOrigin
## #'
## #' @rdname MsBackendMassbankDb
## setMethod("filterDataOrigin", "MsBackendMassbankDb", function(object, dataOrigin, ...) {
##     stop("Not implemented for ", class(object), ".")
## })

## #' @exportMethod filterDataStorage
## #'
## #' @importMethodsFrom ProtGenerics filterDataStorage
## #'
## #' @rdname MsBackendMassbankDb
## setMethod("filterDataStorage", "MsBackendMassbankDb", function(object, dataStorage, ...) {
##     stop("Not implemented for ", class(object), ".")
## })

## #' @exportMethod filterEmptySpectra
## #'
## #' @importMethodsFrom ProtGenerics filterEmptySpectra
## #'
## #' @rdname MsBackendMassbankDb
## setMethod("filterEmptySpectra", "MsBackendMassbankDb", function(object, ...) {
##     stop("Not implemented for ", class(object), ".")
## })

## #' @exportMethod filterIsolationWindow
## #'
## #' @importMethodsFrom ProtGenerics filterIsolationWindow
## #'
## #' @rdname MsBackendMassbankDb
## setMethod("filterIsolationWindow", "MsBackendMassbankDb", function(object, mz, ...) {
##     stop("Not implemented for ", class(object), ".")
## })

## #' @exportMethod filterMsLevel
## #'
## #' @importMethodsFrom ProtGenerics filterMsLevel
## #'
## #' @rdname MsBackendMassbankDb
## setMethod("filterMsLevel", "MsBackendMassbankDb", function(object, msLevel) {
##     stop("Not implemented for ", class(object), ".")
## })

## #' @exportMethod filterPolarity
## #'
## #' @importMethodsFrom ProtGenerics filterPolarity
## #'
## #' @rdname MsBackendMassbankDb
## setMethod("filterPolarity", "MsBackendMassbankDb", function(object, polarity) {
##     stop("Not implemented for ", class(object), ".")
## })

## #' @exportMethod filterPrecursorMz
## #'
## #' @importMethodsFrom ProtGenerics filterPrecursorMz
## #'
## #' @rdname MsBackendMassbankDb
## setMethod("filterPrecursorMz", "MsBackendMassbankDb", function(object, mz) {
##     stop("Not implemented for ", class(object), ".")
## })

## #' @exportMethod filterPrecursorScan
## #'
## #' @importMethodsFrom ProtGenerics filterPrecursorScan
## #'
## #' @rdname MsBackendMassbankDb
## setMethod("filterPrecursorScan", "MsBackendMassbankDb", function(object,
##                                                        acquisitionNum, ...) {
##     stop("Not implemented for ", class(object), ".")
## })

## #' @exportMethod filterRt
## #'
## #' @importMethodsFrom ProtGenerics filterRt
## #'
## #' @rdname MsBackendMassbankDb
## setMethod("filterRt", "MsBackendMassbankDb", function(object, rt, msLevel, ...) {
##     stop("Not implemented for ", class(object), ".")
## })

## #' @exportMethod intensity
## #'
## #' @importMethodsFrom ProtGenerics intensity
## #'
## #' @rdname MsBackendMassbankDb
## setMethod("intensity", "MsBackendMassbankDb", function(object) {
##     stop("Not implemented for ", class(object), ".")
## })

## #' @exportMethod intensity<-
## #'
## #' @importMethodsFrom ProtGenerics intensity<-
## #'
## #' @rdname MsBackendMassbankDb
## setReplaceMethod("intensity", "MsBackendMassbankDb", function(object, value) {
##     stop("Not implemented for ", class(object), ".")
## })

## #' @exportMethod ionCount
## #'
## #' @importMethodsFrom ProtGenerics ionCount
## #'
## #' @rdname MsBackendMassbankDb
## setMethod("ionCount", "MsBackendMassbankDb", function(object) {
##     stop("Not implemented for ", class(object), ".")
## })

## #' @exportMethod isCentroided
## #'
## #' @importMethodsFrom ProtGenerics isCentroided
## #'
## #' @rdname MsBackendMassbankDb
## setMethod("isCentroided", "MsBackendMassbankDb", function(object, ...) {
##     stop("Not implemented for ", class(object), ".")
## })

## #' @exportMethod isEmpty
## #'
## #' @rdname MsBackendMassbankDb
## #'
## #' @importMethodsFrom S4Vectors isEmpty
## setMethod("isEmpty", "MsBackendMassbankDb", function(x) {
##     stop("Not implemented for ", class(x), ".")
## })

## #' @exportMethod isolationWindowLowerMz
## #'
## #' @importMethodsFrom ProtGenerics isolationWindowLowerMz
## #'
## #' @rdname MsBackendMassbankDb
## setMethod("isolationWindowLowerMz", "MsBackendMassbankDb", function(object) {
##     stop("Not implemented for ", class(object), ".")
## })

## #' @exportMethod isolationWindowLowerMz<-
## #'
## #' @importMethodsFrom ProtGenerics isolationWindowLowerMz<-
## #'
## #' @rdname MsBackendMassbankDb
## setReplaceMethod("isolationWindowLowerMz", "MsBackendMassbankDb", function(object,
##                                                                  value) {
##     stop("Not implemented for ", class(object), ".")
## })

## #' @exportMethod isolationWindowTargetMz
## #'
## #' @importMethodsFrom ProtGenerics isolationWindowTargetMz
## #'
## #' @rdname MsBackendMassbankDb
## setMethod("isolationWindowTargetMz", "MsBackendMassbankDb", function(object) {
##     stop("Not implemented for ", class(object), ".")
## })

## #' @exportMethod isolationWindowTargetMz<-
## #'
## #' @importMethodsFrom ProtGenerics isolationWindowTargetMz<-
## #'
## #' @rdname MsBackendMassbankDb
## setReplaceMethod("isolationWindowTargetMz", "MsBackendMassbankDb", function(object,
##                                                                   value) {
##     stop("Not implemented for ", class(object), ".")
## })

## #' @exportMethod isolationWindowUpperMz
## #'
## #' @importMethodsFrom ProtGenerics isolationWindowUpperMz
## #'
## #' @rdname MsBackendMassbankDb
## setMethod("isolationWindowUpperMz", "MsBackendMassbankDb", function(object) {
##     stop("Not implemented for ", class(object), ".")
## })

## #' @exportMethod isolationWindowUpperMz<-
## #'
## #' @importMethodsFrom ProtGenerics isolationWindowUpperMz<-
## #'
## #' @rdname MsBackendMassbankDb
## setReplaceMethod("isolationWindowUpperMz", "MsBackendMassbankDb", function(object,
##                                                                  value) {
##     stop("Not implemented for ", class(object), ".")
## })

## #' @exportMethod isReadOnly
## #'
## #' @rdname MsBackendMassbankDb
## setMethod("isReadOnly", "MsBackendMassbankDb", function(object) {
##     object@readonly
## })

#' @exportMethod length
#'
#' @rdname MsBackendMassbankDb
setMethod("length", "MsBackendMassbankDb", function(x) {
    length(x@spectraIds)
})

## #' @exportMethod msLevel
## #'
## #' @importMethodsFrom ProtGenerics msLevel
## #'
## #' @rdname MsBackendMassbankDb
## setMethod("msLevel", "MsBackendMassbankDb", function(object) {
##     stop("Not implemented for ", class(object), ".")
## })

## #' @exportMethod mz
## #'
## #' @importMethodsFrom ProtGenerics mz
## #'
## #' @rdname MsBackendMassbankDb
## setMethod("mz", "MsBackendMassbankDb", function(object) {
##     stop("Not implemented for ", class(object), ".")
## })

## #' @exportMethod mz<-
## #'
## #' @importMethodsFrom ProtGenerics mz<-
## #'
## #' @rdname MsBackendMassbankDb
## setReplaceMethod("mz", "MsBackendMassbankDb", function(object, value) {
##     stop("Not implemented for ", class(object), ".")
## })

## #' @rdname MsBackendMassbankDb
## setMethod("lengths", "MsBackendMassbankDb", function(x, use.names = FALSE) {
##     stop("Not implemented for ", class(x), ".")
## })

## #' @exportMethod polarity
## #'
## #' @importMethodsFrom ProtGenerics polarity
## #'
## #' @rdname MsBackendMassbankDb
## setMethod("polarity", "MsBackendMassbankDb", function(object) {
##     stop("Not implemented for ", class(object), ".")
## })

## #' @exportMethod polarity<-
## #'
## #' @importMethodsFrom ProtGenerics polarity<-
## #'
## #' @rdname MsBackendMassbankDb
## setReplaceMethod("polarity", "MsBackendMassbankDb", function(object, value) {
##     stop("Not implemented for ", class(object), ".")
## })

## #' @exportMethod precScanNum
## #'
## #' @importMethodsFrom ProtGenerics precScanNum
## #'
## #' @rdname MsBackendMassbankDb
## setMethod("precScanNum", "MsBackendMassbankDb", function(object) {
##     stop("Not implemented for ", class(object), ".")
## })

## #' @exportMethod precursorCharge
## #'
## #' @importMethodsFrom ProtGenerics precursorCharge
## #'
## #' @rdname MsBackendMassbankDb
## setMethod("precursorCharge", "MsBackendMassbankDb", function(object) {
##     stop("Not implemented for ", class(object), ".")
## })

## #' @exportMethod precursorIntensity
## #'
## #' @importMethodsFrom ProtGenerics precursorIntensity
## #'
## #' @rdname MsBackendMassbankDb
## setMethod("precursorIntensity", "MsBackendMassbankDb", function(object) {
##     stop("Not implemented for ", class(object), ".")
## })

## #' @exportMethod precursorMz
## #'
## #' @importMethodsFrom ProtGenerics precursorMz
## #'
## #' @rdname MsBackendMassbankDb
## setMethod("precursorMz", "MsBackendMassbankDb", function(object) {
##     stop("Not implemented for ", class(object), ".")
## })

## #' @exportMethod peaksData<-
## #'
## #' @rdname MsBackendMassbankDb
## setReplaceMethod("peaksData", "MsBackendMassbankDb", function(object, value) {
##     stop("Not implemented for ", class(object), ".")
## })

## #' @exportMethod reset
## #'
## #' @rdname MsBackendMassbankDb
## setMethod("reset", "MsBackendMassbankDb", function(object) {
##     object
## })

## #' @exportMethod rtime
## #'
## #' @importMethodsFrom ProtGenerics rtime
## #'
## #' @rdname MsBackendMassbankDb
## setMethod("rtime", "MsBackendMassbankDb", function(object) {
##     stop("Not implemented for ", class(object), ".")
## })

## #' @exportMethod rtime<-
## #'
## #' @importMethodsFrom ProtGenerics rtime<-
## #'
## #' @rdname MsBackendMassbankDb
## setReplaceMethod("rtime", "MsBackendMassbankDb", function(object, value) {
##     stop("Not implemented for ", class(object), ".")
## })

## #' @exportMethod scanIndex
## #'
## #' @importMethodsFrom ProtGenerics scanIndex
## #'
## #' @rdname MsBackendMassbankDb
## setMethod("scanIndex", "MsBackendMassbankDb", function(object) {
##     stop("Not implemented for ", class(object), ".")
## })

## #' @exportMethod selectSpectraVariables
## #'
## #' @rdname MsBackendMassbankDb
## setMethod(
##     "selectSpectraVariables", "MsBackendMassbankDb",
##     function(object, spectraVariables = spectraVariables(object)) {
##         stop("Not implemented for ", class(object), ".")
##     })

## #' @exportMethod smoothed
## #'
## #' @importMethodsFrom ProtGenerics smoothed
## #'
## #' @rdname MsBackendMassbankDb
## setMethod("smoothed", "MsBackendMassbankDb", function(object) {
##     stop("Not implemented for ", class(object), ".")
## })

## #' @exportMethod smoothed<-
## #'
## #' @aliases smoothed<-,MsBackendMassbankDb-method
## #'
## #' @importMethodsFrom ProtGenerics smoothed<-
## #'
## #' @rdname MsBackendMassbankDb
## setReplaceMethod("smoothed", "MsBackendMassbankDb", function(object, value) {
##     stop("Not implemented for ", class(object), ".")
## })

## #' @exportMethod spectraData
## #'
## #' @rdname MsBackendMassbankDb
## setMethod(
##     "spectraData", "MsBackendMassbankDb",
##     function(object, columns = spectraVariables(object)) {
##         stop("Not implemented for ", class(object), ".")
##     })

## #' @exportMethod spectraData<-
## #'
## #' @rdname MsBackendMassbankDb
## setReplaceMethod("spectraData", "MsBackendMassbankDb", function(object, value) {
##     stop("Not implemented for ", class(object), ".")
## })

## #' @exportMethod spectraNames
## #'
## #' @importMethodsFrom ProtGenerics spectraNames
## #'
## #' @rdname MsBackendMassbankDb
## setMethod("spectraNames", "MsBackendMassbankDb", function(object) {
##     stop("Not implemented for ", class(object), ".")
## })

## #' @exportMethod spectraNames<-
## #'
## #' @importMethodsFrom ProtGenerics spectraNames<-
## #'
## #' @rdname MsBackendMassbankDb
## setReplaceMethod("spectraNames", "MsBackendMassbankDb", function(object, value) {
##     stop("Not implemented for ", class(object), ".")
## })

#' @exportMethod spectraVariables
#'
#' @importMethodsFrom ProtGenerics spectraVariables
#'
#' @rdname MsBackendMassbankDb
setMethod("spectraVariables", "MsBackendMassbankDb", function(object) {
    unique(c(names(Spectra:::.SPECTRA_DATA_COLUMNS), object@spectraVariables))
})

## #' @exportMethod split
## #'
## #' @importMethodsFrom S4Vectors split
## #'
## #' @rdname MsBackendMassbankDb
## setMethod("split", "MsBackendMassbankDb", function(x, f, drop = FALSE, ...) {
##     split.default(x, f, drop = drop, ...)
## })

## #' @exportMethod tic
## #'
## #' @importMethodsFrom ProtGenerics tic
## #'
## #' @rdname MsBackendMassbankDb
## setMethod("tic", "MsBackendMassbankDb", function(object, initial = TRUE) {
##     stop("Not implemented for ", class(object), ".")
## })

## #' @exportMethod [
## #'
## #' @rdname MsBackendMassbankDb
## setMethod("[", "MsBackendMassbankDb", function(x, i, j, ..., drop = FALSE) {
##     stop("Not implemented for ", class(x), ".")
## })

## #' @exportMethod $
## #'
## #' @rdname MsBackendMassbankDb
## setMethod("$", "MsBackendMassbankDb", function(x, name) {
##     stop("Not implemented for ", class(x), ".")
## })

## #' @exportMethod $<-
## #'
## #' @rdname MsBackendMassbankDb
## setReplaceMethod("$", "MsBackendMassbankDb", function(x, name, value) {
##     stop("Not implemented for ", class(x), ".")
## })

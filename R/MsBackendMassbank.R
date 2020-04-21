
#' @include hidden_aliases.R
NULL

#' @title MS data backend for mgf files
#'
#' @aliases MsBackendMassbank-class
#'
#' @description
#'
#' The `MsBackendMassbank` class supports import of MS/MS spectra data from
#' MS/MS spectrum data from [Massbank](https://github.com/MassBank/MassBank-data)
#' files. After initial import, the full MS data is kept in
#' memory. `MsBackendMassbank` extends the [MsBackendDataFrame()] backend
#' directly and supports thus the [applyProcessing()] function to make
#' data manipulations persistent. The backend does however not
#' support export to mgf files yet.
#'
#' New objects are created with the `MsBackendMassbank` function. The
#' `backendInitialize` method has to be subsequently called to
#' initialize the object and import MS/MS data from (one or more) mgf
#' files.  Optional parameter `nonStop` allows to specify whether the
#' import returns with an error if one of the xml files lacks required
#' data, such as `mz` and `intensity` values (default `nonStop =
#' FALSE`), or whether only affected file(s) is(are) skipped and a
#' warning is shown (`nonStop = TRUE`). Note that any other error
#' (such as xml import error) will abort import regardless of
#' parameter `nonStop`.
#'
#' @param object Instance of `MsBackendMassbank` class.
#'
#' @param files `character` with the (full) file name(s) of the mgf file(s)
#'     from which MS/MS data should be imported.
#'
#' @param metaBlocks `data.frame` data frame indicating which metadata shall
#'     be read Default is [metaDataBlocks()].
#'
#' @param nonStop `logical(1)` whether import should be stopped if an
#'     xml file does not contain all required fields. Defaults to
#'     `nonStop = FALSE`.
#'
#' @param BPPARAM Parameter object defining the parallel processing
#'     setup to import data in parallel. Defaults to `BPPARAM =
#'     bpparam()`. See [bpparam()] for more information.
#'
#' @param ... Currently ignored.
#'
#' @author Michael Witting
#'
#' @importClassesFrom Spectra MsBackendDataFrame
#'
#' @exportClass MsBackendMassbank
#'
#' @name MsBackendMassbank
#'
#' @examples
#'
#' ## Create an MsBackendMassbank backend and import data from test xml files.
#' fls <- dir(system.file("extdata", package = "MsBackendMassbank"),
#'     full.names = TRUE, pattern = "txt$")
#' be <- backendInitialize(MsBackendMassbank(), fls)
#' be
#'
#' be$msLevel
#' be$intensity
#' be$mz
NULL

setClass("MsBackendMassbank",
         contains = "MsBackendDataFrame",
         prototype = prototype(spectraData = DataFrame(),
                               readonly = FALSE,
                               version = "0.1"))

#' @importMethodsFrom Spectra backendInitialize asDataFrame<- $<- $
#'
#' @importFrom BiocParallel bpparam
#'
#' @importMethodsFrom BiocParallel bplapply
#'
#' @importFrom methods validObject
#'
#' @exportMethod backendInitialize
#'
#' @rdname MsBackendMassbank
setMethod("backendInitialize", signature = "MsBackendMassbank",
          function(object, files, metaBlocks = metaDataBlocks(),
                   nonStop = FALSE, ..., BPPARAM = bpparam()) {
            if (missing(files) || !length(files))
              stop("Parameter 'files' is mandatory for ", class(object))
            if (!is.character(files))
              stop("Parameter 'files' is expected to be a character vector",
                   " with the files names from where data should be",
                   " imported")
            files <- normalizePath(files)
            if (any(!file.exists(files)))
              stop("file(s) ",
                   paste(files[!file.exists(files)], collapse = ", "),
                   " not found")
            ## Import data and rbind.
            message("Start data import from ", length(files), " files ... ",
                    appendLF = FALSE)
            res <- bplapply(files, FUN = .read_massbank, metaBlocks = metaBlocks,
                            nonStop = nonStop, BPPARAM = BPPARAM)
            message("done")
            res <- do.call(rbind, res)
            if (nonStop && length(files) > nrow(res))
              warning("Import failed for ", length(files) - nrow(res),
                      " files")
            asDataFrame(object) <- res
            object$dataStorage <- "<memory>"
            object$centroided <- TRUE
            validObject(object)
            object
          })

#' @rdname MsBackendMassbank
#'
#' @importFrom methods new
#'
#' @export MsBackendMassbank
MsBackendMassbank <- function() {
  new("MsBackendMassbank")
}

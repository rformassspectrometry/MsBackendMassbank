##' @param f `character(1)` with the path to an MassBank file.
##'
##' @param msLevel `numeric(1)` with the MS level. Default is 2.
##'
##' @param ... Additional parameters, currently ignored.
##'
##' @importFrom S4Vectors DataFrame
##'
##' @importFrom IRanges NumericList
##'
##' @author Michael Witting
##'
##' @noRd
.read_massbank <- function(f, msLevel = 2L, ...) {

  if (length(f) != 1L)
    stop("Please provide a single mgf file.")

  mb <- scan(file = f, what = "",
              sep = "\n", quote = "",
              allowEscapes = FALSE,
              quiet = TRUE)

  begin <- grep("ACCESSION:", mb) + 1L
  end <- grep("//", mb)
  n <- length(begin)
  sp <- vector("list", length = n)
  meta <- vector("list", length = n)

  for (i in seq(along = sp)) {
    sp[[i]] <- .extract_mb_spectrum(mb[begin[i]:end[i]])
    meta[[i]] <- .extract_mb_metadata(mb[begin[i]:end[i]])
  }

  res <- DataFrame(do.call(rbind, sp))
  res_meta <- DataFrame(do.call(rbind, meta))

  print(res)
  print(res_meta)

  res <- cbind.DataFrame(res, res_meta)

  for (i in seq_along(res)) {
    if (all(lengths(res[[i]]) == 1))
      res[[i]] <- unlist(res[[i]])
  }

  res$mz <- IRanges::NumericList(res$mz)
  res$intensity <- IRanges::NumericList(res$intensity)
  res$dataOrigin <- f
  res$msLevel <- as.integer(msLevel)
  res

}




##' @param mb `character()` of lines defining a spectrum in mgf
##'     format.
##'
##' @author Michael Witting
##'
##' @noRd
.extract_mb_spectrum <- function(mb) {

  #read the spectrum
  spectrum_start <- grep("PK$PEAK:", mb, fixed = TRUE) + 1
  spectrum_end <- tail(grep("//", mb, fixed = TRUE), 1) - 1

  splitted <- strsplit(mb[spectrum_start:(spectrum_end)]," ")
  spectrum <- matrix(nrow = spectrum_end + 1 - spectrum_start, ncol = 3)

  for(k in 1:length(splitted)){
    splitted[[k]] <- splitted[[k]][which(splitted[[k]] != "")]
    spectrum[k,] <- splitted[[k]]
  }

  # convert to data frame and adjust data type
  spectrum <- as.data.frame(spectrum, stringsAsFactors = FALSE)
  spectrum[] <- lapply(spectrum, type.convert)
  colnames(spectrum) <- c("mz", "intensity", "rel.intensity")

  # isolate spectrum metadata from record
  rtime <- substring(grep("AC$CHROMATOGRAPHY: RETENTION_TIME",
                          mb,
                          value = TRUE,
                          fixed = TRUE),
                     35)
  rtime <- as.numeric(regmatches(rtime,
                                 regexpr("[[:digit:]]+\\.[[:digit:]]+",
                                         rtime)))

  precursorMz <- as.numeric(substring(grep("MS$FOCUSED_ION: PRECURSOR_M/Z",
                                           mb,
                                           value = TRUE,
                                           fixed = TRUE),
                                      30))

  precursorIntensity <- as.numeric(substring(grep("MS$FOCUSED_ION: PRECURSOR_INT",
                                                  mb,
                                                  value = TRUE,
                                                  fixed = TRUE),
                                             31))

  title <- substring(grep("RECORD_TITLE:",
                          mb,
                          value = TRUE,
                          fixed = TRUE),
                     15)

  # check data
  if(!length(precursorIntensity))
    precursorIntensity <- 100

  list(rtime = rtime * 60,
       scanIndex = as.integer(1),
       precursorMz = precursorMz,
       precursorIntensity = precursorIntensity,
       precursorCharge = as.integer(0),
       mz = spectrum$mz,
       intensity = spectrum$intensity,
       title = title)

}
##' @param mb `character()` of lines defining a spectrum in mgf
##'     format.
##'
##' @author Michael Witting
##'
##' @noRd
.extract_mb_metadata <- function(mb) {

  # create empty list
  ac <- list()

  # analytical chemistry information, MS instrument ----------------------------
  ac$instrument <- substring(grep("AC$INSTRUMENT:", mb, value = TRUE, fixed = TRUE), 16)
  ac$instrument_type <- substring(grep("AC$INSTRUMENT_TYPE:", mb, value = TRUE, fixed = TRUE), 21)

  ac

}



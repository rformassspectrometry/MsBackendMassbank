test_that("MsBackendMassbankSql class works", {
    obj <- new("MsBackendMassbankSql", dbcon = dbc)
    expect_true(validObject(obj))
    obj <- new("MsBackendMassbankSql")
    expect_true(validObject(obj))

    show(obj)
})

test_that("backendInitialize,MsBackendMassbankSql works", {
    expect_error(backendInitialize(MsBackendMassbankSql()), "required")

    be <- backendInitialize(MsBackendMassbankSql(), dbcon = dbc)
    expect_true(length(be@spectraIds) > 0)
    expect_true(length(be@spectraVariables) > 0)
    show(be)
})

test_that("length,MsBackendMassbankSql works", {
    expect_equal(length(MsBackendMassbankSql()), 0L)
    be <- backendInitialize(MsBackendMassbankSql(), dbc)
    expect_equal(length(be), 70L)
})

test_that("dataStorage,MsBackendMassbankSql works", {
    be <- MsBackendMassbankSql()
    expect_equal(dataStorage(be), character())

    be <- backendInitialize(MsBackendMassbankSql(), dbc)
    expect_equal(dataStorage(be), rep("<MassBank>", length(be)))
})

test_that("spectraVariables,MsBackendMassbankSql works", {
    be <- MsBackendMassbankSql()
    res <- spectraVariables(be)
    expect_equal(res, names(Spectra:::.SPECTRA_DATA_COLUMNS))

    be <- backendInitialize(MsBackendMassbankSql(), dbc)
    res <- spectraVariables(be)
    expect_true(length(res) > length(names(Spectra:::.SPECTRA_DATA_COLUMNS)))

    be$a <- 1
    res_2 <- spectraVariables(be)
    expect_true(length(res_2) == (length(res) + 1))
})

test_that("spectraData,MsBackendMassbankSql works", {
    be <- MsBackendMassbankSql()
    res <- spectraData(be)
    expect_true(is(res, "DataFrame"))
    expect_true(nrow(res) == 0)

    be <- backendInitialize(be, dbc)
    res <- spectraData(be)
    expect_true(is(res, "DataFrame"))
    expect_true(nrow(res) > 0)

    res <- spectraData(be, columns = "authors")
    expect_true(is(res, "DataFrame"))
    expect_true(nrow(res) > 0)
    expect_equal(colnames(res), "authors")
})

test_that("$,MsBackendMassbankSql works", {
    be <- MsBackendMassbankSql()
    res <- be$rtime
    expect_true(is.numeric(res))
    expect_true(length(res) == 0)

    res <- be$mz
    expect_true(is(res, "NumericList"))

    expect_error(be$other, "not available")
})

test_that("$<-,MsBackendMassbankSql works", {
    be <- backendInitialize(MsBackendMassbankSql(), dbc)
    ## Test adding new column
    be$new_col <- "a"
    expect_equal(be$new_col, rep("a", length(be)))

    ## Test replacing column
    be$new_col <- seq_len(length(be))
    expect_equal(be$new_col, seq_len(length(be)))

    be$rtime <- seq_len(length(be))
    expect_equal(be$rtime, seq_len(length(be)))

    be$authors <- "a"
    expect_equal(be$authors, rep("a", length(be)))

    ## Test errors with m/z etc.
    expect_error(be$rtime <- c(1, 2), "length 1 or")
    expect_error(be$mz <- 1:3, "Replacing m/z and intensity")
})

test_that("acquisitionNum,MsBackendMassbankSql works", {
    be <- MsBackendMassbankSql()
    res <- acquisitionNum(be)
    expect_equal(res, integer())

    be <- backendInitialize(be, dbc)
    res <- acquisitionNum(be)
    expect_equal(res, rep(NA_integer_, length(be)))
})

test_that("centroided,centroided<-,MsBackendMassbankSql works", {
    be <- MsBackendMassbankSql()
    res <- centroided(be)
    expect_equal(res, logical())

    be <- backendInitialize(be, dbc)
    res <- centroided(be)
    expect_equal(res, rep(NA, length(be)))

    centroided(be) <- FALSE
    res <- centroided(be)
    expect_equal(res, rep(FALSE, length(be)))

    expect_error(centroided(be) <- 3, "logical")
})

test_that("collisionEnergy,collisionEnergy<-,MsBackendMassbankSql works", {
    be <- MsBackendMassbankSql()
    res <- collisionEnergy(be)
    expect_equal(res, numeric())

    be <- backendInitialize(be, dbc)
    res <- collisionEnergy(be)
    expect_true(is.numeric(res))

    be$collisionEnergy <- 1.2
    res <- collisionEnergy(be)
    expect_equal(res, rep(1.2, length(be)))

    expect_error(collisionEnergy(be) <- "a", "numeric")
})

test_that("peaksData,MsBackendMassbankSql works", {
    be <- MsBackendMassbankSql()
    res <- peaksData(be)
    expect_true(is.list(res))

    be <- backendInitialize(be, dbc)
    res <- peaksData(be)
    expect_true(is.list(res))
    expect_true(length(res) == length(be))
    expect_true(is.matrix(res[[1]]))
    expect_true(is.matrix(res[[2]]))
    expect_true(is.matrix(res[[3]]))

    expect_equal(lapply(res, function(z) unname(z[, 1])),
                 unname(as.list(be$mz)))

    ## duplicated spectra.
    be2 <- be[c(3, 1, 1, 2)]
    res2 <- peaksData(be2)
    expect_equal(lapply(res2, function(z) unname(z[, 1])),
                 unname(as.list(be2$mz)))
    expect_equal(res2[[1]], res[[3]])
    expect_equal(res2[[2]], res[[1]])
    expect_equal(res2[[3]], res[[1]])
    expect_equal(res2[[4]], res[[2]])
})

test_that("dataOrigin, dataOrigin<-,MsBackendMassbankSql works", {
    be <- MsBackendMassbankSql()
    res <- dataOrigin(be)
    expect_equal(res, character())

    be <- backendInitialize(be, dbc)
    res <- dataOrigin(be)
    expect_true(is.character(res))

    dataOrigin(be) <- "b"
    res <- dataOrigin(be)
    expect_equal(res, rep("b", length(be)))

    expect_error(dataOrigin(be) <- 1, "character")
})

test_that("selectSpectraVariables,MsBackendMassbankSql works", {
    be <- backendInitialize(MsBackendMassbankSql(), dbc)
    res <- selectSpectraVariables(be, spectraVariables = spectraVariables(be))
    expect_equal(spectraVariables(be), spectraVariables(res))

    ## errors
    expect_error(selectSpectraVariables(be, c("rtime", "other")), "available")

    res <- selectSpectraVariables(be, c("rtime", "mz", "intensity", "authors"))
    expect_equal(sort(spectraVariables(res)),
                 sort(unique(c(names(Spectra:::.SPECTRA_DATA_COLUMNS), "rtime",
                               "mz", "intensity", "authors"))))
    res$new_col <- "b"
    expect_equal(sort(spectraVariables(res)),
                 sort(unique(c(names(Spectra:::.SPECTRA_DATA_COLUMNS), "rtime",
                               "mz", "intensity", "authors", "new_col"))))
    res_2 <- selectSpectraVariables(res, spectraVariables(res))
    expect_equal(spectraVariables(res), spectraVariables(res_2))
    res <- selectSpectraVariables(res, c("rtime"))
    expect_equal(
        spectraVariables(res),
        unique(c(names(Spectra:::.SPECTRA_DATA_COLUMNS)), "rtime"))
})

test_that("[,MsBackendMassbankSql works", {
    be <- backendInitialize(MsBackendMassbankSql(), dbc)

    res <- be[2:3]
    expect_true(length(res) == 2)
    expect_equal(res@spectraIds, be@spectraIds[2:3])
    expect_equal(res@spectraIds, res$spectrum_id)
    expect_equal(be$mz[2:3], res$mz)

    res <- be[c(3, 1)]
    expect_true(length(res) == 2)
    expect_equal(res@spectraIds, be@spectraIds[c(3, 1)])
    expect_equal(be$mz[c(3, 1)], res$mz)
    expect_equal(be$splash[c(3, 1)], res$splash)

    res <- be[c(3, 1, 2, 3, 1)]
    expect_equal(mz(res)[[1]], mz(be)[[3]])
    expect_equal(mz(res)[[2]], mz(be)[[1]])
    expect_equal(mz(res)[[3]], mz(be)[[2]])
    expect_equal(mz(res)[[4]], mz(be)[[3]])
    expect_equal(mz(res)[[5]], mz(be)[[1]])
    expect_equal(intensity(res)[[1]], intensity(be)[[3]])
    expect_equal(intensity(res)[[2]], intensity(be)[[1]])
    expect_equal(intensity(res)[[3]], intensity(be)[[2]])
    expect_equal(intensity(res)[[4]], intensity(be)[[3]])
    expect_equal(intensity(res)[[5]], intensity(be)[[1]])
    expect_equal(res$spectrum_id[1], be$spectrum_id[3])
    expect_equal(res$spectrum_id[2], be$spectrum_id[1])
    expect_equal(res$spectrum_id[3], be$spectrum_id[2])
    expect_equal(res$spectrum_id[4], be$spectrum_id[3])
    expect_equal(res$spectrum_id[5], be$spectrum_id[1])
    expect_equal(res$spectrum_name[1], be$spectrum_name[3])
    expect_equal(res$spectrum_name[2], be$spectrum_name[1])
    expect_equal(res$spectrum_name[3], be$spectrum_name[2])
    expect_equal(res$spectrum_name[4], be$spectrum_name[3])
    expect_equal(res$spectrum_name[5], be$spectrum_name[1])
})

test_that("lengths,MsBackendMassbankSql works", {
    be <- MsBackendMassbankSql()
    expect_equal(lengths(be), integer())

    be <- backendInitialize(be, dbc)
    res <- lengths(be)
    expect_true(length(res) > 0)
    expect_true(all(res > 0))
})

test_that("intensity,intensity<-,MsBackendMassbankSql works", {
    be <- MsBackendMassbankSql()
    res <- intensity(be)
    expect_equal(res, IRanges::NumericList(compress = FALSE))

    be <- backendInitialize(be, dbc)
    res <- intensity(be)
    expect_true(is(res, "NumericList"))
    expect_true(length(res) == length(be))
    expect_true(all(lengths(res) > 0))

    expect_error(intensity(be) <- 3, "Can not")

    ## duplicated spectra
    be2 <- be[c(4, 1, 9, 1, 4, 1)]
    res2 <- intensity(be2)
    expect_equal(res[[4]], res2[[1]])
    expect_equal(res[[1]], res2[[2]])
    expect_equal(res[[9]], res2[[3]])
    expect_equal(res[[1]], res2[[4]])
    expect_equal(res[[4]], res2[[5]])
    expect_equal(res[[1]], res2[[6]])
})

test_that("mz,mz<-,MsBackendMassbankSql works", {
    be <- MsBackendMassbankSql()
    res <- mz(be)
    expect_equal(res, IRanges::NumericList(compress = FALSE))

    be <- backendInitialize(be, dbc)
    res <- mz(be)
    expect_true(is(res, "NumericList"))
    expect_true(length(res) == length(be))
    expect_true(all(lengths(res) > 0))

    expect_error(mz(be) <- 3, "Can not")

    ## duplicated spectra.
    be2 <- be[c(4, 1, 9, 1, 4, 1)]
    res2 <- mz(be2)
    expect_equal(res[[4]], res2[[1]])
    expect_equal(res[[1]], res2[[2]])
    expect_equal(res[[9]], res2[[3]])
    expect_equal(res[[1]], res2[[4]])
    expect_equal(res[[4]], res2[[5]])
    expect_equal(res[[1]], res2[[6]])
})

test_that("ionCount,MsBackendMassbankSql works", {
    be <- MsBackendMassbankSql()
    res <- ionCount(be)
    expect_equal(res, numeric())

    be <- backendInitialize(be, dbc)
    res <- ionCount(be)
    expect_true(is.numeric(res))
    expect_equal(length(res), length(be))
    expect_true(all(res > 0))
})

test_that("isEmpty,MsBackendMassbankSql works", {
    be <- MsBackendMassbankSql()
    res <- isEmpty(be)
    expect_equal(res, logical())

    be <- backendInitialize(be, dbc)
    res <- isEmpty(be)
    expect_true(all(!res))
})

test_that("isolationWindowLowerMz,&<-,MsBackendMassbankSql works", {
    be <- MsBackendMassbankSql()
    res <- isolationWindowLowerMz(be)
    expect_equal(res, numeric())

    be <- backendInitialize(be, dbc)
    res <- isolationWindowLowerMz(be)
    expect_equal(res, rep(NA_real_, length(be)))

    isolationWindowLowerMz(be) <- seq_len(length(be))
    res <- isolationWindowLowerMz(be)
    expect_equal(res, seq_len(length(be)))

    expect_error(isolationWindowLowerMz(be) <- 1:2, "length 1")
    expect_error(isolationWindowLowerMz(be) <- "a", "numeric")
})

test_that("isolationWindowTargetMz,&<-,MsBackendMassbankSql works", {
    be <- MsBackendMassbankSql()
    res <- isolationWindowTargetMz(be)
    expect_equal(res, numeric())

    be <- backendInitialize(be, dbc)
    res <- isolationWindowTargetMz(be)
    expect_equal(res, rep(NA_real_, length(be)))

    isolationWindowTargetMz(be) <- seq_len(length(be))
    res <- isolationWindowTargetMz(be)
    expect_equal(res, seq_len(length(be)))

    expect_error(isolationWindowTargetMz(be) <- 1:2, "length 1")
    expect_error(isolationWindowTargetMz(be) <- "a", "numeric")
})

test_that("isolationWindowUpperMz,&<-,MsBackendMassbankSql works", {
    be <- MsBackendMassbankSql()
    res <- isolationWindowUpperMz(be)
    expect_equal(res, numeric())

    be <- backendInitialize(be, dbc)
    res <- isolationWindowUpperMz(be)
    expect_equal(res, rep(NA_real_, length(be)))

    isolationWindowUpperMz(be) <- seq_len(length(be))
    res <- isolationWindowUpperMz(be)
    expect_equal(res, seq_len(length(be)))

    expect_error(isolationWindowUpperMz(be) <- 1:2, "length 1")
    expect_error(isolationWindowUpperMz(be) <- "a", "numeric")
})

test_that("polarity,&<-,MsBackendMassbankSql works", {
    be <- MsBackendMassbankSql()
    res <- polarity(be)
    expect_equal(res, integer())

    be <- backendInitialize(be, dbc)
    res <- polarity(be)
    expect_true(all(res %in% c(0L, 1L)))

    polarity(be) <- rep(c(1, 0), length(be)/2)
    res <- polarity(be)
    expect_equal(res, rep(c(1, 0), length(be)/2))

    expect_error(polarity(be) <- 1:2, "length 1")
    expect_error(polarity(be) <- "a", "integer")
})

test_that("precursorCharge,MsBackendMassbankSql works", {
    be <- MsBackendMassbankSql()
    res <- precursorCharge(be)
    expect_equal(res, integer())

    be <- backendInitialize(be, dbc)
    res <- precursorCharge(be)
    expect_equal(res, rep(NA_integer_, length(be)))

    be$precursorCharge <- 1L
    res <- precursorCharge(be)
    expect_equal(res, rep(1L, length(be)))
})

test_that("precursorIntensity,MsBackendMassbankSql works", {
    be <- MsBackendMassbankSql()
    res <- precursorIntensity(be)
    expect_equal(res, numeric())

    be <- backendInitialize(be, dbc)
    res <- precursorIntensity(be)
    expect_true(is.numeric(res))

    nmbrs <- abs(rnorm(length(be)))
    be$precursorIntensity <- nmbrs
    res <- precursorIntensity(be)
    expect_equal(res, nmbrs)
})

test_that("precursorMz,MsBackendMassbankSql works", {
    be <- MsBackendMassbankSql()
    res <- precursorMz(be)
    expect_equal(res, numeric())

    be <- backendInitialize(be, dbc)
    res <- precursorMz(be)
    expect_true(is.numeric(res))

    be$precursorMz <- 12.21
    res <- precursorMz(be)
    expect_equal(res, rep(12.21, length(be)))
})

test_that("reset,MsBackendMassbankSql works", {
    be <- MsBackendMassbankSql()
    res <- reset(be)
    expect_equal(res, be)

    be <- backendInitialize(be, dbc)
    res <- reset(be)
    expect_equal(res, be)

    be$new_col <- "c"
    be <- be[c(3, 1)]
    res <- reset(be)
    expect_true(length(res) == 70)
    expect_true(!any(spectraVariables(res) == "new_col"))
})

test_that("rtime,rtime<-,MsBackendMassbankSql works", {
    be <- MsBackendMassbankSql()
    res <- rtime(be)
    expect_equal(res, numeric())

    be <- backendInitialize(be, dbc)
    res <- rtime(be)
    expect_equal(res, rep(NA_real_, length(be)))

    rtime(be) <- 1.4
    res <- rtime(be)
    expect_equal(res, rep(1.4, length(be)))

    expect_error(rtime(be) <- 1:2, "length 1")
    expect_error(rtime(be) <- "a", "numeric")
})

test_that("scanIndex,MsBackendMassbankSql works", {
    be <- MsBackendMassbankSql()
    res <- scanIndex(be)
    expect_equal(res, integer())

    be <- backendInitialize(be, dbc)
    res <- scanIndex(be)
    expect_equal(res, rep(NA_integer_, length(be)))

    be$scanIndex <- seq_len(length(be))
    res <- scanIndex(be)
    expect_equal(res, seq_len(length(be)))
})

test_that("smoothed,smoothed<-,MsBackendMassbankSql works", {
    be <- MsBackendMassbankSql()
    res <- smoothed(be)
    expect_equal(res, logical())

    be <- backendInitialize(be, dbc)
    res <- smoothed(be)
    expect_equal(res, rep(NA, length(be)))

    smoothed(be) <- rep(c(TRUE, FALSE), length(be)/2)
    res <- smoothed(be)
    expect_equal(res, rep(c(TRUE, FALSE), length(be)/2))

    expect_error(smoothed(be) <- c(TRUE, FALSE), "length 1")
    expect_error(smoothed(be) <- "a", "logical")
})

test_that("spectraData<-,MsBackendMassbankSql works", {
    be <- MsBackendMassbankSql()
    expect_error(spectraData(be) <- spectraData(be), "not support")
})

test_that("spectraNames,spectraNames<-,MsBackendMassbankSql works", {
    be <- MsBackendMassbankSql()
    res <- spectraNames(be)
    expect_equal(res, character())

    be <- backendInitialize(be, dbc)
    res <- spectraNames(be)
    expect_true(is.character(res))
    expect_true(length(res) == length(be))
    expect_equal(res, be@spectraIds)

    expect_error(spectraNames(be) <- "a", "not support")
})

test_that("tic,MsBackendMassbankSql works", {
    be <- MsBackendMassbankSql()
    res <- tic(be)
    expect_equal(res, numeric())

    be <- backendInitialize(be, dbc)
    res <- tic(be)
    expect_true(is.numeric(res))
    expect_true(all(is.na(res)))

    be$totIonCurrent <- 1.12
    res <- tic(be)
    expect_equal(res, rep(1.12, length(be)))

    res <- tic(be, initial = FALSE)
    expect_true(is.numeric(res))
    expect_true(length(res) == length(be))
    expect_true(all(res > 0))
})

test_that("msLevel,MsBackendMassbankSql works", {
    be <- MsBackendMassbankSql()
    res <- msLevel(be)
    expect_equal(res, integer())

    be <- backendInitialize(be, dbc)
    res <- msLevel(be)
    expect_true(is.integer(res))
    expect_true(length(res) == length(be))
    expect_true(all(res %in% c(NA_integer_, 2L)))
})

test_that("compounds,MsBackendMassbankSql works", {
    be <- MsBackendMassbankSql()
    res <- compounds(be)
    expect_true(is(res, "DataFrame"))
    expect_true(nrow(res) == 0)

    be <- backendInitialize(be, dbc)
    be <- be[c(3, 1, 3)]
    res <- compounds(be)
    expect_equal(res$compound_id, be$compound_id)
    expect_equal(res$formula[1], res$formula[3])
})

test_that("compounds,Spectra works", {
    sps <- Spectra(backendInitialize(MsBackendMassbankSql(), dbc))
    res <- compounds(sps)
})

test_that("precursorMz,Spectra works", {
    sps <- Spectra(backendInitialize(MsBackendMassbankSql(), dbc))
    precursors_cached <- precursorMz(sps)
    sps@backend@localData$precursorMz <- NULL
    precursors_direct <- precursorMz(sps)
    expect_equal(precursors_cached, precursors_direct)
})

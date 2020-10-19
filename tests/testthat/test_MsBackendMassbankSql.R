test_that("MsBackendMassbankSql class works", {
    obj <- new("MsBackendMassbankSql", dbcon = dbc)
    expect_true(validObject(obj))
    obj <- new("MsBackendMassbankSql")
    expect_true(validObject(obj))
})

test_that("backendInitialize,MsBackendMassbankSql works", {
    expect_error(backendInitialize(MsBackendMassbankSql()), "required")

    be <- backendInitialize(MsBackendMassbankSql(), dbcon = dbc)
    expect_true(length(be@spectraIds) > 0)
    expect_true(length(be@spectraVariables) > 0)
})

test_that("length,MsBackendMassbankSql works", {
    expect_equal(length(MsBackendMassbankSql()), 0L)
    be <- backendInitialize(MsBackendMassbankSql(), dbc)
    expect_equal(length(be), 3L)
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

    be@localData <- DataFrame(a = 1:3)
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
    expect_equal(be$new_col, c("a", "a", "a"))

    ## Test replacing column
    be$new_col <- 1:3
    expect_equal(be$new_col, 1:3)

    be$rtime <- 1:3
    expect_equal(be$rtime, 1:3)

    be$authors <- "a"
    expect_equal(be$authors, c("a", "a", "a"))

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
    expect_equal(res, rep(NA_integer_, 3))
})

test_that("centroided,centroided<-,MsBackendMassbankSql works", {
    be <- MsBackendMassbankSql()
    res <- centroided(be)
    expect_equal(res, logical())

    be <- backendInitialize(be, dbc)
    res <- centroided(be)
    expect_equal(res, rep(NA, 3))

    centroided(be) <- FALSE
    res <- centroided(be)
    expect_equal(res, rep(FALSE, 3))

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
    expect_equal(res, rep(1.2, 3))

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
    expect_equal(res, rep("b", 3))

    expect_error(dataOrigin(be) <- 1, "character")
})

test_that("selectSpectraVariables,MsBackendMassbankSql works", {
    be <- backendInitialize(MsBackendMassbankSql(), dbc)
    res <- selectSpectraVariables(be, spectraVariables = spectraVariables(be))
    expect_equal(spectraVariables(be), spectraVariables(res))

    ## errors
    expect_error(selectSpectraVariables(be, c("rtime", "other")), "available")

    res <- selectSpectraVariables(be, c("rtime", "mz", "intensity", "authors"))
    expect_equal(spectraVariables(res),
                 c("rtime", "mz", "intensity", "authors"))
    res$new_col <- "b"
    expect_equal(spectraVariables(res),
                 c("rtime", "mz", "intensity", "new_col", "authors"))
    res <- selectSpectraVariables(res, c("rtime"))
    expect_equal(spectraVariables(res), "rtime")
})

test_that("[,MsBackendMassbankSql works", {
    be <- backendInitialize(MsBackendMassbankSql(), dbc)

    res <- be[2:3]
    expect_true(length(res) == 2)
    expect_equal(res@spectraIds, be@spectraIds[2:3])
    expect_equal(be$mz[2:3], res$mz)

    res <- be[c(3, 1)]
    expect_true(length(res) == 2)
    expect_equal(res@spectraIds, be@spectraIds[c(3, 1)])
    expect_equal(be$mz[c(3, 1)], res$mz)
    expect_equal(be$splash[c(3, 1)], res$splash)
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
    expect_equal(res, IRanges::NumericList())

    be <- backendInitialize(be, dbc)
    res <- intensity(be)
    expect_true(is(res, "NumericList"))
    expect_true(length(res) == length(be))
    expect_true(all(lengths(res) > 0))

    expect_error(intensity(be) <- 3, "Can not")
})

test_that("mz,mz<-,MsBackendMassbankSql works", {
    be <- MsBackendMassbankSql()
    res <- mz(be)
    expect_equal(res, IRanges::NumericList())

    be <- backendInitialize(be, dbc)
    res <- mz(be)
    expect_true(is(res, "NumericList"))
    expect_true(length(res) == length(be))
    expect_true(all(lengths(res) > 0))

    expect_error(mz(be) <- 3, "Can not")
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

    isolationWindowLowerMz(be) <- c(1, 2, 3)
    res <- isolationWindowLowerMz(be)
    expect_equal(res, 1:3)

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

    isolationWindowTargetMz(be) <- c(1, 2, 3)
    res <- isolationWindowTargetMz(be)
    expect_equal(res, 1:3)

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

    isolationWindowUpperMz(be) <- c(1, 2, 3)
    res <- isolationWindowUpperMz(be)
    expect_equal(res, 1:3)

    expect_error(isolationWindowUpperMz(be) <- 1:2, "length 1")
    expect_error(isolationWindowUpperMz(be) <- "a", "numeric")
})

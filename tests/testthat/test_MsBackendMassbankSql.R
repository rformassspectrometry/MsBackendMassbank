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

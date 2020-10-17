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

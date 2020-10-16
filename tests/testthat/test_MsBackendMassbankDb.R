test_that("MsBackendMassbankDb class works", {
    obj <- new("MsBackendMassbankDb", dbcon = dbc)
    expect_true(validObject(obj))
    obj <- new("MsBackendMassbankDb")
    expect_true(validObject(obj))
})

test_that("backendInitialize,MsBackendMassbankDb works", {
    expect_error(backendInitialize(MsBackendMassbankDb()), "required")

    be <- backendInitialize(MsBackendMassbankDb(), dbcon = dbc)
    expect_true(length(be@spectraIds) > 0)
    expect_true(length(be@spectraVariables) > 0)
})

test_that("length,MsBackendMassbankDb works", {
    expect_equal(length(MsBackendMassbankDb()), 0L)
    be <- backendInitialize(MsBackendMassbankDb(), dbc)
    expect_equal(length(be), 3L)
})

test_that("dataStorage,MsBackendMassbankDb works", {
    be <- MsBackendMassbankDb()
    expect_equal(dataStorage(be), character())

    be <- backendInitialize(MsBackendMassbankDb(), dbc)
    expect_equal(dataStorage(be), rep("<MassBank>", length(be)))
})

test_that("spectraVariables,MsBackendMassbankDb works", {
    be <- MsBackendMassbankDb()
    res <- spectraVariables(be)
    expect_equal(res, names(Spectra:::.SPECTRA_DATA_COLUMNS))

    be <- backendInitialize(MsBackendMassbankDb(), dbc)
    res <- spectraVariables(be)
    expect_true(length(res) > length(names(Spectra:::.SPECTRA_DATA_COLUMNS)))
})

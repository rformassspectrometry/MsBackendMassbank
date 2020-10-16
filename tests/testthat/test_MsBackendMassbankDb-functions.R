test_that(".valid_dbcon works", {
    expect_true(length(.valid_dbcon(dbc)) == 0)
    expect_true(length(.valid_dbcon(4)) == 1)
})

test_that("MsBackendMassbankDb works", {
    res <- MsBackendMassbankDb()
    expect_true(validObject(res))
    expect_true(is(res, "MsBackendMassbankDb"))
})

test_that(".spectra_data_massbank works", {
    be <- MsBackendMassbankDb()

    res <- .spectra_data_massbank(be)
})

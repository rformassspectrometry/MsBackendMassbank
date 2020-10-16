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

test_that(".fetch_peaks_sql works", {
    be <- MsBackendMassbankDb()
    res <- .fetch_peaks_sql(be)
    expect_true(is.data.frame(res))
    expect_equal(colnames(res), c("spectrum_id", "mz", "intensity"))
    expect_true(nrow(res) == 0)

    be <- backendInitialize(MsBackendMassbankDb(), dbc)
    res <- .fetch_peaks_sql(be)
    expect_true(is.data.frame(res))
    expect_equal(colnames(res), c("spectrum_id", "mz", "intensity"))
    expect_true(nrow(res) > 0)
})

test_that(".fetch_spectra_data_sql works", {
    be <- backendInitialize(MsBackendMassbankDb(), dbc)
    res <- .fetch_spectra_data_sql(be)
    expect_true(nrow(res) > 0)
    expect_true(colnames(res) == "spectrum_id")

    res <- .fetch_spectra_data_sql(be, columns = "msLevel")
    expect_true(all(res$msLevel == 2L))

    res <- .fetch_spectra_data_sql(be, columns = "polarity")
    expect_true(all(res$polarity == 1L))
})

test_that(".spectra_data_massbank_sql works", {
    be <- MsBackendMassbankDb()
    res <- .spectra_data_massbank_sql(be)
    expect_true(is(res, "DataFrame"))
    expect_true(nrow(res) == 0L)
    expect_equal(colnames(res), names(Spectra:::.SPECTRA_DATA_COLUMNS))

    be <- backendInitialize(MsBackendMassbankDb(), dbc)
    ## Full data.
    expect_error(.spectra_data_massbank_sql(be, columns = "other"), "available")
    res <- .spectra_data_massbank_sql(be)
    expect_true(is(res, "DataFrame"))
    expect_equal(colnames(res), spectraVariables(be))
    expect_true(is(res$mz, "NumericList"))
    expect_true(is(res$intensity, "NumericList"))
    expect_true(all(is.na(res$rtime)))

    ## A single column
    res <- .spectra_data_massbank_sql(be, columns = "rtime")
    expect_true(is(res, "DataFrame"))
    expect_true(nrow(res) > 0)
    expect_equal(colnames(res), "rtime")
    expect_true(all(is.na(res$rtime)))

    res <- .spectra_data_massbank_sql(be, "msLevel")
    expect_true(is(res, "DataFrame"))
    expect_true(nrow(res) > 0)
    expect_equal(colnames(res), "msLevel")
    expect_true(all(res$msLevel == 2L))

    res <- .spectra_data_massbank_sql(be, "mz")
    expect_true(is(res, "DataFrame"))
    expect_true(nrow(res) > 0)
    expect_equal(colnames(res), "mz")
    expect_true(is(res$mz, "NumericList"))

    ## A combination of core and db
    res <- .spectra_data_massbank_sql(be, columns = c("rtime", "polarity"))
    expect_true(is(res, "DataFrame"))
    expect_true(nrow(res) > 0)
    expect_equal(colnames(res), c("rtime", "polarity"))
    expect_true(all(res$polarity == 1L))
    expect_true(all(is.na(res$rtime)))

    res <- .spectra_data_massbank_sql(
        be, columns = c("mz", "polarity", "intensity"))
    expect_true(is(res, "DataFrame"))
    expect_true(nrow(res) > 0)
    expect_equal(colnames(res), c("mz", "polarity", "intensity"))
    expect_true(all(res$polarity == 1L))
    expect_true(is(res$mz, "NumericList"))
    expect_true(is(res$intensity, "NumericList"))

    ## With local data.
    stop("implement local data")
    ## - validation method that checks that nrow @localData and @spectrumIds
    ##   have same length.
    ## - spectraVariables adds also columns from there (unique one).
    ## - spectraData takes columns FIRST from there (not from other sources).
)

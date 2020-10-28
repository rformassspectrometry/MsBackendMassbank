test_that(".valid_dbcon works", {
    expect_true(length(.valid_dbcon(dbc)) == 0)
    expect_true(length(.valid_dbcon(4)) == 1)
})

test_that(".valid_local_data works", {
    expect_true(length(.valid_local_data(DataFrame(), 3)) == 0)
    expect_true(length(.valid_local_data(data.frame(a = 4), 1:3)) == 1)
})

test_that("MsBackendMassbankSql works", {
    res <- MsBackendMassbankSql()
    expect_true(validObject(res))
    expect_true(is(res, "MsBackendMassbankSql"))
})

test_that(".fetch_peaks_sql works", {
    be <- MsBackendMassbankSql()
    res <- .fetch_peaks_sql(be)
    expect_true(is.data.frame(res))
    expect_equal(colnames(res), c("spectrum_id", "mz", "intensity"))
    expect_true(nrow(res) == 0)

    be <- backendInitialize(MsBackendMassbankSql(), dbc)
    res <- .fetch_peaks_sql(be)
    expect_true(is.data.frame(res))
    expect_equal(colnames(res), c("spectrum_id", "mz", "intensity"))
    expect_true(nrow(res) > 0)
})

test_that(".fetch_spectra_data_sql works", {
    be <- backendInitialize(MsBackendMassbankSql(), dbc)
    res <- .fetch_spectra_data_sql(be)
    expect_true(nrow(res) > 0)
    expect_true(colnames(res) == "spectrum_id")

    res <- .fetch_spectra_data_sql(be, columns = "msLevel")
    expect_true(all(res$msLevel %in% c(NA_integer_, 2L)))

    res <- .fetch_spectra_data_sql(be, columns = "polarity")
    expect_true(all(res$polarity %in% c(1L, 0L)))
})

test_that(".spectra_data_massbank_sql works", {
    be <- MsBackendMassbankSql()
    res <- .spectra_data_massbank_sql(be)
    expect_true(is(res, "DataFrame"))
    expect_true(nrow(res) == 0L)
    expect_equal(colnames(res), names(Spectra:::.SPECTRA_DATA_COLUMNS))

    be <- backendInitialize(MsBackendMassbankSql(), dbc)
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
    expect_true(all(res$msLevel %in% c(NA, 2L)))

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
    expect_true(all(res$polarity %in% c(0L, 1L)))
    expect_true(all(is.na(res$rtime)))

    res <- .spectra_data_massbank_sql(
        be, columns = c("mz", "polarity", "intensity"))
    expect_true(is(res, "DataFrame"))
    expect_true(nrow(res) > 0)
    expect_equal(colnames(res), c("mz", "polarity", "intensity"))
    expect_true(all(res$polarity %in% c(0L, 1L)))
    expect_true(is(res$mz, "NumericList"))
    expect_true(is(res$intensity, "NumericList"))

    ## With local data.
    be@localData <- DataFrame(rtime = seq_len(length(be)))
    res <- .spectra_data_massbank_sql(
        be, columns = c("mz", "rtime", "spectrum_id"))
    expect_true(is(res, "DataFrame"))
    expect_true(nrow(res) > 0)
    expect_equal(res$rtime, seq_len(length(be)))

    be@localData$authors <- rep("A", length(be))
    res <- .spectra_data_massbank_sql(be, columns = c("rtime", "authors"))
    expect_true(is(res, "DataFrame"))
    expect_true(nrow(res) > 0)
    expect_equal(res$authors, rep("A", length(be)))
})

test_that(".compounds_sql works", {
    res <- .compounds_sql(dbc, id = c(2, 4, 6))
    expect_true(all(c(2, 4, 6) %in% res$compound_id))
})

test_that(".map_spectraVariables_to_sql works", {
    res <- .map_spectraVariables_to_sql("msLevel")
    expect_equal(res, "ms_level")
    res <- .map_spectraVariables_to_sql(c("mz", "collisionEnergy", "ms_level"))
    expect_equal(res, c("mz", "collisionEnergy", "ms_level"))
})

test_that(".map_sql_to_spectraVariables works", {
    res <- .map_sql_to_spectraVariables("ms_level")
    expect_equal(res, "msLevel")
    res <- .map_sql_to_spectraVariables(c("precursor_mz_text", "intensity"))
    expect_equal(res, c("precursorMz", "intensity"))
})

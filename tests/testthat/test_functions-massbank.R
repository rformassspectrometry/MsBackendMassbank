fls <- system.file("extdata", "MassBankRecords.txt",
                   package = "MsBackendMassbank")
mb <- scan(fls, what = "", sep = "\n", quote = "",
           allowEscapes = FALSE, quiet = TRUE)
mb <- mb[seq_len(53)]

test_that(".read_massbank works", {
    expect_error(.read_massbank(c(fls, fls)), "single")

    fail <- dir(system.file("extdata", package = "MsBackendMassbank"),
                 pattern = "metadata_blocks.csv",
                full.names = TRUE)
    expect_error(.read_massbank(fail), "Unexpected file format")
    res <- expect_warning(.read_massbank(fail, nonStop = TRUE), "Unexpected")
    expect_equal(res, DataFrame())

    res <- .read_massbank(fls)
    expect_true(is(res, "DataFrame"))
    expect_true(nrow(res) == 3L)
    expect_true(all(res$msLevel == 2L))

    mblock <- metaDataBlocks()
    mblock[1, 2] <- TRUE
    res_2 <- .read_massbank(fls, metaBlocks = mblock)
    expect_true(ncol(res) < ncol(res_2))

    mblock[2, 2] <- TRUE
    res_3 <- .read_massbank(fls, metaBlocks = mblock)
    expect_true(ncol(res_2) < ncol(res_3))

    mblock[3, 2] <- TRUE
    res_4 <- .read_massbank(fls, metaBlocks = mblock)
    expect_true(ncol(res_3) < ncol(res_4))

    mblock[4, 2] <- TRUE
    res_5 <- .read_massbank(fls, metaBlocks = mblock)
    expect_true(ncol(res_4) < ncol(res_5))

    mblock[5, 2] <- TRUE
    res_6 <- .read_massbank(fls, metaBlocks = mblock)
    expect_true(ncol(res_5) < ncol(res_6))

    mblock[6, 2] <- TRUE
    res_7 <- .read_massbank(fls, metaBlocks = mblock)
    expect_true(ncol(res_6) < ncol(res_7))

    mblock[7, 2] <- TRUE
    res_8 <- .read_massbank(fls, metaBlocks = mblock)
    expect_true(ncol(res_7) < ncol(res_8))
})

test_that(".extract_mb_spectrum works", {
    res <- .extract_mb_spectrum(mb)
    expect_true(is.list(res))
    expect_true(is.integer(res$polarity))
    expect_true(is.numeric(res$rtime))
    expect_true(is.character(res$title))
    expect_true(is.list(res$name))

    fl <- system.file("extdata", "RP000501_mod.txt",
                      package = "MsBackendMassbank")
    tmp <- scan(fl, what = "", sep = "\n", quote = "",
           allowEscapes = FALSE, quiet = TRUE)
    res <- .extract_mb_spectrum(tmp)
    expect_identical(res$collisionEnergy, NA_real_)
    ## collisionenery completely missing
    tmp <- tmp[-grep("COLLISION_ENERGY", tmp)]
    res <- .extract_mb_spectrum(tmp)
    expect_identical(res$collisionEnergy, NA_real_)
    tmp <- tmp[-grep("RETENTION_TIME", tmp)]
    res <- .extract_mb_spectrum(tmp)
    expect_identical(res$rtime, NA_real_)
    idx <- grep("ION_MODE", tmp)
    tmp[idx] <- sub("POSITIVE", "NEGATIVE", tmp[idx])
    res <- .extract_mb_spectrum(tmp)
    expect_identical(res$polarity, 0L)
    tmp[idx] <- sub("NEGATIVE", "OTHER", tmp[idx])
    res <- .extract_mb_spectrum(tmp)
    expect_identical(res$polarity, NA_integer_)
    tmp <- tmp[-grep("PRECURSOR_M/Z", tmp)]
    res <- .extract_mb_spectrum(tmp)
    expect_identical(res$precursorMz, NA_real_)
})

test_that(".extract_mb_ac works", {
    res <- .extract_mb_ac(mb)
    expect_true(is.list(res))
    expect_true(length(res) == 43L)
    expect_true(is.character(res$instrument))
})

test_that(".extract_mb_ch works", {
    res <- .extract_mb_ch(mb)
    expect_true(is.list(res))
    expect_true(length(res) == 16L)
})

test_that(".extract_mb_sp works", {
    res <- .extract_mb_sp(mb)
    expect_true(is.list(res))
    expect_true(length(res) == 4L)
})

test_that(".extract_mb_ms works", {
    res <- .extract_mb_ms(mb)
    expect_true(is.list(res))
    expect_true(length(res) == 11L)
})

test_that(".extract_mb_record works", {
    res <- .extract_mb_record(mb)
    expect_true(is.list(res))
    expect_true(length(res) == 7L)
})

test_that(".extract_mb_pk works", {
    res <- .extract_mb_pk(mb)
    expect_true(is.list(res))
    expect_true(is.integer(res$pknum))
})

test_that(".extract_mb_comment works", {
    res <- .extract_mb_comment(mb)
    expect_true(is.character(res))
})

test_that("metaDataBlocks works", {
    res <- metaDataBlocks()
    expect_true(is.data.frame(res))
    expect_true(all(res$read == FALSE))
})

test_that(".cleanParsing works", {
    tmp <- list(3, "b", "", integer())
    res <- .cleanParsing(tmp)
    expect_true(is.na(res[[4]]))
})

test_that(".export_massbank works", {
    fls <- system.file("extdata","MassBankRecords.txt",
                       package = "MsBackendMassbank")
    be <- backendInitialize(MsBackendMassbank(), fls)
    sps <- Spectra(be)

    tmpf <- tempfile()
    .export_massbank(sps, tmpf)

    be2 <- backendInitialize(MsBackendMassbank(), tmpf)
    expect_equal(be$msLevel, be2$msLevel)
    expect_equal(be$title, be2$title)
    expect_equal(be$formula, be2$formula)
    expect_equal(be$name, be2$name)
    expect_equal(be$mz, be2$mz)
    expect_equal(be$intensity, be2$intensity)
    expect_equal(be$accession, be2$accession)
    expect_equal(be$polarity, be2$polarity)

    ## All with negative polarity...
    sps$polarity <- 0L
    expect_message(.export_massbank(sps, tmpf), "Overwriting")
    be2 <- backendInitialize(MsBackendMassbank(), tmpf)
    expect_equal(be$msLevel, be2$msLevel)
    expect_true(all(be2$polarity == 0L))
})

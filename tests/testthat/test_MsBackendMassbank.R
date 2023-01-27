test_that("backendInitialize,MsBackendMassbank works", {
    ## Import a single file with multiple record
    fls <- system.file("extdata","MassBankRecords.txt",
                       package = "MsBackendMassbank")
    be <- MsBackendMassbank()

    res_single <- backendInitialize(be, fls)

    expect_identical(length(res_single), 3L)
    expect_identical(res_single$dataStorage, rep("<memory>", 3))
    expect_identical(res_single$dataOrigin, rep(normalizePath(fls[1]), 3))
    expect_identical(res_single$msLevel, rep(2L, 3))
    expect_identical(length(spectraVariables(res_single)), 31L)


    ## Import a single file with multiple record, only spectrum
    metaDataBlocks <- data.frame(
        metadata = c("ac", "ch", "sp", "ms", "record", "pk", "comment"),
        read = rep(TRUE, 7),
        stringsAsFactors = FALSE)

    res_single <- backendInitialize(be, fls, metaBlocks = metaDataBlocks)

    expect_identical(length(res_single), 3L)
    expect_identical(res_single$dataStorage, rep("<memory>", 3))
    expect_identical(res_single$dataOrigin, rep(normalizePath(fls[1]), 3))
    expect_identical(res_single$msLevel, rep(2L, 3))
    expect_identical(length(spectraVariables(res_single)), 114L)

    ## Import multiple files, single entries
    fls <- dir(system.file("extdata", package = "MsBackendMassbank"),
               full.names = TRUE, pattern = "^RP.*txt$")

    res_multiple <- backendInitialize(be, fls)

    expect_true(length(res_multiple) == 6)
    expect_true(all(res_multiple$msLevel == 2L))
    expect_identical(res_multiple$dataOrigin, normalizePath(fls))
    expect_true(is.integer(res_multiple@spectraData$msLevel))

    ## Import a single file without RT information
    fls <- system.file("extdata", "BSU00001.txt", package = "MsBackendMassbank")
    be <- MsBackendMassbank()

    res_single <- backendInitialize(be, fls)

    expect_identical(res_single$rtime, c(rtime = NA_real_))

    ## TODO: Import with failing file.
    ## TODO: Import with failing file and nonStop = TRUE

    ## errors
    expect_error(backendInitialize(be), "'files' is mandatory")
    expect_error(backendInitialize(be, 4), "expected to be a character")
    expect_error(backendInitialize(be, "a"), "a not found")
})

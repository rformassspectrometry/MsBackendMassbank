test_that("backendInitialize,MsBackendMgf works", {

  # get test files
  fls <- dir(system.file("extdata", package = "MsBackendMassbank"),
             full.names = TRUE, pattern = "txt$")
  be <- MsBackendMassbank()

  ## Import a single file with multiple record
  res1 <- backendInitialize(be, fls[1])

  n1 <- length(res1) ## 3
  expect_identical(length(res1), n1)
  expect_identical(res1$dataStorage, rep("<memory>", n1))
  expect_identical(res1$dataOrigin, rep(normalizePath(fls[1]), n1))
  expect_identical(res1$msLevel, rep(2L, n1))

  ## Import multiple files.
  res_all <- backendInitialize(be, fls)
  expect_true(length(res_all) == 9)
  #expect_identical(res_all[1]$mz, res1[1]$mz)
  expect_true(all(res_all$msLevel == 2L))
  expect_identical(res_all$dataOrigin, c(rep(normalizePath(fls)[1], 3),
                                         normalizePath(fls[-1])))
  expect_true(is.integer(res_all@spectraData$msLevel))

  # Import a single file with multiple record, only spectrum
  metaDataBlocks <- data.frame(metadata = c("ac", "ch", "sp", "ms", "record", "pk", "comment"),
                               read = rep(FALSE, 7),
                               stringsAsFactors = FALSE)

  be <- MsBackendMassbank()

  ## Import a single file with multiple record
  res1 <- backendInitialize(be, fls[1], metaBlocks = metaDataBlocks)

  n1 <- length(res1) ## 3
  expect_identical(length(res1), n1)
  expect_identical(res1$dataStorage, rep("<memory>", n1))
  expect_identical(res1$dataOrigin, rep(normalizePath(fls[1]), n1))
  expect_identical(res1$msLevel, rep(2L, n1))

  ## TODO: Import with failing file.
  ## TODO: Import with failing file and nonStop = TRUE

  ## errors
  expect_error(backendInitialize(be), "'files' is mandatory")
  expect_error(backendInitialize(be, 4), "expected to be a character")
  expect_error(backendInitialize(be, "a"), "a not found")
})

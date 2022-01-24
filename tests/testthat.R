library(testthat)
library(MsBackendMassbank)

library(RSQLite)
dbc <- dbConnect(SQLite(), system.file("sql", "minimassbank.sqlite",
                                       package = "MsBackendMassbank"))

test_check("MsBackendMassbank")

## Run additional tests from Spectra:
test_suite <- system.file("test_backends", "test_MsBackend",
                          package = "Spectra")

## Test MsBackendMassbank:
fls <- dir(system.file("extdata", package = "MsBackendMassbank"),
           full.names = TRUE, pattern = "^RP.*txt$")
be <- MsBackendMassbank()
be <- backendInitialize(be, fls)

res <- test_file(paste0(test_suite, "/test_spectra_variables.R"),
                 reporter = check_reporter(), stop_on_failure = TRUE)

## Test MsBackendMassbankSql
be <- MsBackendMassbankSql()
be <- backendInitialize(be, dbc)
res <- test_file(paste0(test_suite, "/test_spectra_variables.R"),
                 reporter = check_reporter(), stop_on_failure = TRUE)

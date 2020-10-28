library(testthat)
library(MsBackendMassbank)

library(RSQLite)
dbc <- dbConnect(SQLite(), system.file("sql/minimassbank.sql",
                                       package = "MsBackendMassbank"))

test_check("MsBackendMassbank")

# load required libraries
library(MsBackendMassBank)

# load all functions, including non-exported fuctions
devtools::load_all()

context("parsing of MassBank records")

test_that("correct isolation of chemical information", {

  # get test file with chemical information example
  ch_test <- system.file("testfiles/ch_example.txt", package = "MsBackendMassbank")

  # load test file
  con <- file(ch_test)
  mb_record <- readLines(con)
  close(con)

  # isolate chemical information
  ch <- .isolate_chemical(mb_record)

  # test different chemical information
  expect_equal(length(ch$name), 2)
  expect_equal(ch$name[[1]], "Testcompound 1")
  expect_equal(ch$name[[2]], "Testcompound 2")
  expect_equal(ch$compound_class, "Testclass")
  expect_equal(ch$formula, "C6H12O6")
  expect_equal(ch$exact_mass, 180)
  expect_equal(ch$smiles, "CCCCCCCC")
  expect_equal(ch$iupac, "InChI")
  expect_equal(ch$link_cas, "100-100-100")
  expect_equal(ch$link_cayman, "Testcayman")
  expect_equal(ch$link_chebi, "Testchebi")
  expect_equal(ch$link_chembl, "Testchembl")
  expect_equal(ch$link_chempdb, "Testchempdb")
  expect_equal(ch$link_chemspider, "Testchemspider")
  expect_equal(ch$link_comptox, "Testcomptox")
  expect_equal(ch$link_hmdb, "Testhmdb")
  expect_equal(ch$link_inchikey, "Testinchi")
  expect_equal(ch$link_kappaview, "Testkappaview")
  expect_equal(ch$link_kegg, "Testkegg")
  expect_equal(ch$link_knapsack, "Testknapsack")
  expect_equal(ch$link_lipidbank, "Testlipidbank")
  expect_equal(ch$link_lipidmaps, "Testlipidmaps")
  expect_equal(ch$link_nikkaji, "Testnikkaji")
  expect_equal(ch$link_pubchem, "Testpubchem")
  expect_equal(ch$link_zinc, "Testzinc")

})


test_that("correct isolation of species information", {

  # get test file with chemical information example
  sp_test <- system.file("testfiles/sp_example.txt", package = "MsBackendMassbank")

  # load test file
  con <- file(sp_test)
  mb_record <- readLines(con)
  close(con)

  # isolate chemical information
  sp <- .isolate_species(mb_record)

  # test different chemical information
  expect_equal(sp$scientific_name, "Mus musculus")
  expect_equal(sp$lineage, "cellular organisms; Eukaryota; Fungi/Metazoa group; Metazoa; Eumetazoa; Bilateria; Coelomata; Deuterostomia; Chordata; Craniata; Vertebrata; Gnathostomata; Teleostomi; Euteleostomi; Sarcopterygii; Tetrapoda; Amniota; Mammalia; Theria; Eutheria; Euarchontoglires; Glires; Rodentia; Sciurognathi; Muroidea; Muridae; Murinae; Mus")
  expect_equal(sp$link, "NCBI-TAXONOMY 10090")
  expect_equal(sp$sample, "Liver extracts")


})

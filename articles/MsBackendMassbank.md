# Description and usage of MsBackendMassbank

**Package**:
*[MsBackendMassbank](https://bioconductor.org/packages/3.23/MsBackendMassbank)*  
**Authors**: RforMassSpectrometry Package Maintainer \[cre\], Michael
Witting \[aut\] (ORCID: <https://orcid.org/0000-0002-1462-4426>),
Johannes Rainer \[aut\] (ORCID:
<https://orcid.org/0000-0002-6977-7147>), Michael Stravs \[ctb\]  
**Compiled**: Tue Jan 20 08:55:22 2026

## Introduction

The `Spectra` package provides a central infrastructure for the handling
of Mass Spectrometry (MS) data. The package supports interchangeable use
of different *backends* to import MS data from a variety of sources
(such as mzML files). The `MsBackendMassbank` package allows import and
handling MS/MS spectrum data from
[Massbank](https://massbank.eu/MassBank/). This vignette illustrates the
usage of the `MsBackendMassbank` package to include MassBank data into
MS data analysis workflow with the `Spectra` package in R.

## Installation

The package can be installed with the `BiocManager` package. To install
`BiocManager` use `install.packages("BiocManager")` and, after that,
`BiocManager::install("MsBackendMassbank")` to install this package.

## Importing MS/MS data from MassBank files

MassBank files (as provided by the [Massbank github
repository](https://github.com/MassBank/MassBank-data)) store normally
one library spectrum per file, typically centroided and of MS level 2.
In our short example below, we load data from a file containing multiple
library spectra per file or from files with each a single spectrum
provided with this package. Below we first load all required packages
and define the paths to the Massbank files.

``` r

library(Spectra)
library(MsBackendMassbank)
fls <- dir(system.file("extdata", package = "MsBackendMassbank"),
           full.names = TRUE, pattern = "txt$")
fls
```

    ##  [1] "/__w/_temp/Library/MsBackendMassbank/extdata/BSU00001.txt"          
    ##  [2] "/__w/_temp/Library/MsBackendMassbank/extdata/MassBankRecords.txt"   
    ##  [3] "/__w/_temp/Library/MsBackendMassbank/extdata/MSBNK-UFZ-UA000101.txt"
    ##  [4] "/__w/_temp/Library/MsBackendMassbank/extdata/multi_precursor_mz.txt"
    ##  [5] "/__w/_temp/Library/MsBackendMassbank/extdata/RP000501_mod.txt"      
    ##  [6] "/__w/_temp/Library/MsBackendMassbank/extdata/RP000501.txt"          
    ##  [7] "/__w/_temp/Library/MsBackendMassbank/extdata/RP000502.txt"          
    ##  [8] "/__w/_temp/Library/MsBackendMassbank/extdata/RP000503.txt"          
    ##  [9] "/__w/_temp/Library/MsBackendMassbank/extdata/RP000511.txt"          
    ## [10] "/__w/_temp/Library/MsBackendMassbank/extdata/RP000512.txt"          
    ## [11] "/__w/_temp/Library/MsBackendMassbank/extdata/RP000513.txt"

MS data can be accessed and analyzed through `Spectra` objects. Below we
create a `Spectra` with the data from these mgf files. To this end we
provide the file names and specify to use a
[`MsBackendMassbank()`](https://rformassspectrometry.github.io/MsBackendMassbank/reference/MsBackendMassbank.md)
backend as *source* to enable data import. First we import from a single
file with multiple library spectra.

``` r

sps <- Spectra(fls[1],
               source = MsBackendMassbank(),
               backend = MsBackendDataFrame(),
               nonStop = TRUE)
```

With that we have now full access to all imported spectra variables that
we list below.

``` r

spectraVariables(sps)
```

    ##  [1] "msLevel"                 "rtime"                  
    ##  [3] "acquisitionNum"          "scanIndex"              
    ##  [5] "dataStorage"             "dataOrigin"             
    ##  [7] "centroided"              "smoothed"               
    ##  [9] "polarity"                "precScanNum"            
    ## [11] "precursorMz"             "precursorIntensity"     
    ## [13] "precursorCharge"         "collisionEnergy"        
    ## [15] "isolationWindowLowerMz"  "isolationWindowTargetMz"
    ## [17] "isolationWindowUpperMz"  "acquistionNum"          
    ## [19] "accession"               "name"                   
    ## [21] "smiles"                  "exactmass"              
    ## [23] "formula"                 "inchi"                  
    ## [25] "cas"                     "inchikey"               
    ## [27] "adduct"                  "splash"                 
    ## [29] "title"

The same is possible with multiple files, each containing a library
spectrum.

``` r

sps <- Spectra(fls[-1],
               source = MsBackendMassbank(),
               backend = MsBackendDataFrame(),
               nonStop = TRUE)
```

    ## Warning in .local(object, ...): Import failed for some files

``` r

spectraVariables(sps)
```

    ##  [1] "msLevel"                 "rtime"                  
    ##  [3] "acquisitionNum"          "scanIndex"              
    ##  [5] "dataStorage"             "dataOrigin"             
    ##  [7] "centroided"              "smoothed"               
    ##  [9] "polarity"                "precScanNum"            
    ## [11] "precursorMz"             "precursorIntensity"     
    ## [13] "precursorCharge"         "collisionEnergy"        
    ## [15] "isolationWindowLowerMz"  "isolationWindowTargetMz"
    ## [17] "isolationWindowUpperMz"  "acquistionNum"          
    ## [19] "accession"               "name"                   
    ## [21] "smiles"                  "exactmass"              
    ## [23] "formula"                 "inchi"                  
    ## [25] "cas"                     "inchikey"               
    ## [27] "adduct"                  "splash"                 
    ## [29] "title"

By default the complete metadata is read together with the spectra. This
can increase loading time. The different metadata blocks can be skipped
which reduces import time. This requires to define an additional
`data.frame` indicating what shall be read.

``` r

# create data frame to indicate with metadata blocks shall be read.
metaDataBlocks <- data.frame(metadata = c("ac", "ch", "sp", "ms",
                                          "record", "pk", "comment"),
                             read = rep(TRUE, 7))

sps <- Spectra(fls[-1],
               source = MsBackendMassbank(),
               backeend = MsBackendDataFrame(),
               metaBlock = metaDataBlocks,
               nonStop = TRUE)
```

    ## Warning in .local(object, ...): Import failed for some files

``` r

# all spectraVariables possible in MassBank are read
spectraVariables(sps)
```

    ##   [1] "msLevel"                     "rtime"                      
    ##   [3] "acquisitionNum"              "scanIndex"                  
    ##   [5] "dataStorage"                 "dataOrigin"                 
    ##   [7] "centroided"                  "smoothed"                   
    ##   [9] "polarity"                    "precScanNum"                
    ##  [11] "precursorMz"                 "precursorIntensity"         
    ##  [13] "precursorCharge"             "collisionEnergy"            
    ##  [15] "isolationWindowLowerMz"      "isolationWindowTargetMz"    
    ##  [17] "isolationWindowUpperMz"      "acquistionNum"              
    ##  [19] "accession"                   "name"                       
    ##  [21] "smiles"                      "exactmass"                  
    ##  [23] "formula"                     "inchi"                      
    ##  [25] "cas"                         "inchikey"                   
    ##  [27] "adduct"                      "splash"                     
    ##  [29] "title"                       "instrument"                 
    ##  [31] "instrument_type"             "ms_ms_type"                 
    ##  [33] "ms_cap_voltage"              "ms_col_gas"                 
    ##  [35] "ms_desolv_gas_flow"          "ms_desolv_temp"             
    ##  [37] "ms_frag_mode"                "ms_ionization"              
    ##  [39] "ms_ionization_energy"        "ms_laser"                   
    ##  [41] "ms_matrix"                   "ms_mass_accuracy"           
    ##  [43] "ms_mass_range"               "ms_reagent_gas"             
    ##  [45] "ms_resolution"               "ms_scan_setting"            
    ##  [47] "ms_source_temp"              "ms_kinetic_energy"          
    ##  [49] "ms_electron_current"         "ms_reaction_time"           
    ##  [51] "chrom_carrier_gas"           "chrom_column"               
    ##  [53] "chrom_column_temp"           "chrom_column_temp_gradient" 
    ##  [55] "chrom_flow_gradient"         "chrom_flow_rate"            
    ##  [57] "chrom_inj_temp"              "chrom_inj_temp_gradient"    
    ##  [59] "chrom_rti_kovats"            "chrom_rti_lee"              
    ##  [61] "chrom_rti_naps"              "chrom_rti_uoa"              
    ##  [63] "chrom_rti_uoa_pred"          "chrom_rt"                   
    ##  [65] "chrom_rt_uoa_pred"           "chrom_solvent"              
    ##  [67] "chrom_transfer_temp"         "ims_instrument_type"        
    ##  [69] "ims_drift_gas"               "ims_drift_time"             
    ##  [71] "ims_ccs"                     "general_conc"               
    ##  [73] "compound_class"              "link_cayman"                
    ##  [75] "link_chebi"                  "link_chembl"                
    ##  [77] "link_chempdb"                "link_chemspider"            
    ##  [79] "link_comptox"                "link_hmdb"                  
    ##  [81] "link_kappaview"              "link_kegg"                  
    ##  [83] "link_knapsack"               "link_lipidbank"             
    ##  [85] "link_lipidmaps"              "link_nikkaji"               
    ##  [87] "link_pubchem"                "link_zinc"                  
    ##  [89] "scientific_name"             "lineage"                    
    ##  [91] "link"                        "sample"                     
    ##  [93] "focus_base_peak"             "focus_derivative_form"      
    ##  [95] "focus_derivative_mass"       "focus_derivative_type"      
    ##  [97] "focus_ion_type"              "data_processing_comment"    
    ##  [99] "data_processing_deprofile"   "data_processing_find"       
    ## [101] "data_processing_reanalyze"   "data_processing_recalibrate"
    ## [103] "data_processing_whole"       "deprecated"                 
    ## [105] "date"                        "authors"                    
    ## [107] "license"                     "copyright"                  
    ## [109] "publication"                 "project"                    
    ## [111] "pknum"                       "comment"

``` r

# all NA columns can be dropped
spectraVariables(dropNaSpectraVariables(sps))
```

    ##  [1] "msLevel"                     "rtime"                      
    ##  [3] "acquisitionNum"              "scanIndex"                  
    ##  [5] "dataStorage"                 "dataOrigin"                 
    ##  [7] "centroided"                  "smoothed"                   
    ##  [9] "polarity"                    "precScanNum"                
    ## [11] "precursorMz"                 "precursorIntensity"         
    ## [13] "precursorCharge"             "collisionEnergy"            
    ## [15] "isolationWindowLowerMz"      "isolationWindowTargetMz"    
    ## [17] "isolationWindowUpperMz"      "acquistionNum"              
    ## [19] "accession"                   "name"                       
    ## [21] "smiles"                      "exactmass"                  
    ## [23] "formula"                     "inchi"                      
    ## [25] "cas"                         "inchikey"                   
    ## [27] "adduct"                      "splash"                     
    ## [29] "title"                       "instrument"                 
    ## [31] "instrument_type"             "ms_ms_type"                 
    ## [33] "ms_frag_mode"                "ms_ionization"              
    ## [35] "ms_resolution"               "chrom_column"               
    ## [37] "chrom_flow_gradient"         "chrom_flow_rate"            
    ## [39] "chrom_rt"                    "chrom_solvent"              
    ## [41] "compound_class"              "link_chebi"                 
    ## [43] "link_chemspider"             "link_comptox"               
    ## [45] "link_kegg"                   "link_pubchem"               
    ## [47] "focus_base_peak"             "data_processing_reanalyze"  
    ## [49] "data_processing_recalibrate" "data_processing_whole"      
    ## [51] "date"                        "authors"                    
    ## [53] "license"                     "copyright"                  
    ## [55] "publication"                 "pknum"                      
    ## [57] "comment"

Besides default spectra variables, such as `msLevel`, `rtime`,
`precursorMz`, we also have additional spectra variables such as the
`title` of each spectrum in the mgf file.

``` r

sps$rtime
```

    ##  [1] 142.14 142.14 142.14     NA 142.14 142.14 142.14 142.14 143.94 143.94
    ## [11] 143.94

``` r

sps$title
```

    ##  [1] "L-Tryptophan; LC-ESI-QTOF; MS2; CE: 10; R=; [M+H]+"
    ##  [2] "L-Tryptophan; LC-ESI-QTOF; MS2; CE: 20; R=; [M+H]+"
    ##  [3] "L-Tryptophan; LC-ESI-QTOF; MS2; CE: 40; R=; [M+H]+"
    ##  [4] "Carbazole; ESI-ITFT; MS2; CE: 35%; R=30000; [M+H]+"
    ##  [5] "L-Tryptophan; LC-ESI-QTOF; MS2; CE: 10; R=; [M+H]+"
    ##  [6] "L-Tryptophan; LC-ESI-QTOF; MS2; CE: 10; R=; [M+H]+"
    ##  [7] "L-Tryptophan; LC-ESI-QTOF; MS2; CE: 20; R=; [M+H]+"
    ##  [8] "L-Tryptophan; LC-ESI-QTOF; MS2; CE: 40; R=; [M+H]+"
    ##  [9] "L-Tryptophan; LC-ESI-QTOF; MS2; CE: 10; R=; [M-H]-"
    ## [10] "L-Tryptophan; LC-ESI-QTOF; MS2; CE: 20; R=; [M-H]-"
    ## [11] "L-Tryptophan; LC-ESI-QTOF; MS2; CE: 40; R=; [M-H]-"

In addition we can also access the m/z and intensity values of each
spectrum.

``` r

mz(sps)
```

    ## NumericList of length 11
    ## [[1]] 74.0233 132.0807 144.0805 146.0598 ... 170.0597 188.0699 205.0965
    ## [[2]] 74.0232 77.0381 86.0027 91.0539 ... 160.0947 170.0596 171.0625 188.07
    ## [[3]] 53.0019 53.0383 63.0225 65.0381 ... 158.0817 159.0921 160.0755 170.06
    ## [[4]] 115.0167 168.0809
    ## [[5]] 74.0233 132.0807 144.0805 146.0598 ... 170.0597 188.0699 205.0965
    ## [[6]] 74.0233 132.0807 144.0805 146.0598 ... 170.0597 188.0699 205.0965
    ## [[7]] 74.0232 77.0381 86.0027 91.0539 ... 160.0947 170.0596 171.0625 188.07
    ## [[8]] 53.0019 53.0383 63.0225 65.0381 ... 158.0817 159.0921 160.0755 170.06
    ## [[9]] 72.0095 116.0517 117.0554 159.0935 186.0558 203.0826
    ## [[10]] 72.0094 74.0253 116.0511 117.0539 ... 162.0307 186.0548 203.0818
    ## ...
    ## <1 more element>

``` r

intensity(sps)
```

    ## NumericList of length 11
    ## [[1]] 646 980 2114 20052 1248 7628 2036 494048 75708
    ## [[2]] 10186 142 142 750 138 490 126 ... 11254 14266 1478 1600 16504 1446 109762
    ## [[3]] 324 184 138 3770 500 800 7214 3238 ... 2802 206 898 162 166 814 250 1840
    ## [[4]] 650.7 14157.3
    ## [[5]] 646 980 2114 20052 1248 7628 2036 494048 75708
    ## [[6]] 646 980 2114 20052 1248 7628 2036 494048 75708
    ## [[7]] 10186 142 142 750 138 490 126 ... 11254 14266 1478 1600 16504 1446 109762
    ## [[8]] 324 184 138 3770 500 800 7214 3238 ... 2802 206 898 162 166 814 250 1840
    ## [[9]] 150 200 32 232 80 12162
    ## [[10]] 460 3630 658 74 28 190 44 142 50 52 634
    ## ...
    ## <1 more element>

When importing a large number of mgf files, setting `nonStop = TRUE`
prevents the call to stop whenever problematic mgf files are
encountered.

``` r

sps <- Spectra(fls, source = MsBackendMassbank(), nonStop = TRUE)
```

## Accessing the MassBank MySQL database

An alternative to the import of the MassBank data from individual text
files (which can take a considerable amount of time) is to directly
access the MS/MS data in the MassBank MySQL database. For demonstration
purposes we are using here a tiny subset of the MassBank data which is
stored as a SQLite database within this package.

### Pre-requisites

At present it is not possible to directly connect to the main MassBank
*production* MySQL server, thus, to use the `MsBackendMassbankSql`
backend it is required to install the database locally. The MySQL
database dump for each MassBank release can be downloaded from
[here](https://rformassspectrometry.github.io/MsBackendMassbank/articles/).
This dump could be imported to a local MySQL server.

### Direct access to the MassBank database

To use the `MsBackendMassbankSql` it is required to first connect to a
*MassBank* database. Below we show the R code which could be used for
that - but the actual settings (user name, password, database name, or
host) will depend on where and how the MassBank database was installed.

``` r

library(RMariaDB)
con <- dbConnect(MariaDB(), host = "localhost", user = "massbank",
                 dbname = "MassBank")
```

To illustrate the general functionality of this backend we use a tiny
subset of the MassBank (release 2020.10) which is provided as an small
SQLite database within this package. Below we connect to this database.

``` r

library(RSQLite)
con <- dbConnect(SQLite(), system.file("sql", "minimassbank.sqlite",
                                       package = "MsBackendMassbank"))
```

We next *initialize* the `MsBackendMassbankSql` backend which supports
direct access to the MassBank in a SQL database and create a `Spectra`
object from that.

``` r

mb <- Spectra(con, source = MsBackendMassbankSql())
mb
```

    ## MSn data (Spectra) with 70 spectra in a MsBackendMassbankSql backend:
    ##       msLevel precursorMz  polarity
    ##     <integer>   <numeric> <integer>
    ## 1           2         506         0
    ## 2          NA          NA         1
    ## 3          NA          NA         0
    ## 4          NA          NA         1
    ## 5          NA          NA         0
    ## ...       ...         ...       ...
    ## 66          2     185.028         0
    ## 67          2     455.290         1
    ## 68          2     253.051         0
    ## 69          2     358.238         1
    ## 70          2     256.170         1
    ##  ... 41 more variables/columns.
    ##  Use  'spectraVariables' to list all of them.

We can now use this `Spectra` object to access and use the MassBank data
for our analysis. Note that the `Spectra` object itself does not contain
any data from MassBank. Any data will be fetched on demand from the
database backend.

To get a listing of all available annotations for each spectrum (the
so-called *spectra variables*) we can use the `spectraVariables`
function.

``` r

spectraVariables(mb)
```

    ##  [1] "msLevel"                 "rtime"                  
    ##  [3] "acquisitionNum"          "scanIndex"              
    ##  [5] "dataStorage"             "dataOrigin"             
    ##  [7] "centroided"              "smoothed"               
    ##  [9] "polarity"                "precScanNum"            
    ## [11] "precursorMz"             "precursorIntensity"     
    ## [13] "precursorCharge"         "collisionEnergy"        
    ## [15] "isolationWindowLowerMz"  "isolationWindowTargetMz"
    ## [17] "isolationWindowUpperMz"  "spectrum_id"            
    ## [19] "spectrum_name"           "date"                   
    ## [21] "authors"                 "license"                
    ## [23] "copyright"               "publication"            
    ## [25] "splash"                  "compound_id"            
    ## [27] "adduct"                  "ionization"             
    ## [29] "ionization_voltage"      "fragmentation_mode"     
    ## [31] "instrument"              "instrument_type"        
    ## [33] "formula"                 "exactmass"              
    ## [35] "smiles"                  "inchi"                  
    ## [37] "inchikey"                "cas"                    
    ## [39] "pubchem"                 "synonym"                
    ## [41] "precursor_mz_text"       "compound_name"

Through the `MsBackendMassbankSql` we can thus access spectra
information as well as its annotation.

We can access *core* spectra variables, such as the MS level with the
corresponding function `msLevel`.

``` r

head(msLevel(mb))
```

    ## [1]  2 NA NA NA NA  2

Spectra variables can also be accessed with `$` and the name of the
variable. Thus, MS levels can also be accessed with `$msLevel`:

``` r

head(mb$msLevel)
```

    ## [1]  2 NA NA NA NA  2

In addition to spectra variables, we can also get the actual peaks
(i.e. m/z and intensity values) with the `mz` and `intensity` functions:

``` r

mz(mb)
```

    ## NumericList of length 70
    ## [[1]] 146.760803 158.863541 174.988785 ... 470.057434 487.989319 585.88446
    ## [[2]] 22.99 23.07 23.19 38.98 53.15 60.08 ... 391.18 413.22 414.22 429.16 495.2
    ## [[3]] 108.099566 152.005378 153.01824 ... 341.031964 409.06729 499.103447
    ## [[4]] 137.0227 138.025605 139.034602 ... 304.271886 316.303106 352.234228
    ## [[5]] 112.977759 205.053323 208.041128 ... 449.137152 671.1778 673.200866
    ## [[6]] 78.893929 183.011719 193.957916 ... 328.091949 408.010193 426.021942
    ## [[7]] 167.03389 168.040758 333.079965 ... 357.073039 365.040698 373.073473
    ## [[8]] 195.09167 196.095033
    ## [[9]] 108.096149 152.004656 153.01824 ... 359.997563 409.065314 491.043775
    ## [[10]] 1278.12 1279.11 1279.18 1306.21 ... 2064.29 2091.32 2092.3 2136.33
    ## ...
    ## <60 more elements>

Note that not all spectra from the database were generated using the
same instrumentation. Below we list the number of spectra for each type
of instrument.

``` r

table(mb$instrument_type)
```

    ## 
    ## LC-ESI-ITFT  LC-ESI-QFT   MALDI-TOF 
    ##           3          50          17

We next subset the data to all spectra from ions generated by electro
spray ionization (ESI).

``` r

mb <- mb[mb$ionization == "ESI"]
length(mb)
```

    ## [1] 50

As a simple example to illustrate the `Spectra` functionality we next
calculate spectra similarity between one spectrum against all other
spectra in the database. To this end we use the `compareSpectra`
function with the normalized dot product as similarity function and
allowing 20 ppm difference in m/z between matching peaks

``` r

library(MsCoreUtils)
sims <- compareSpectra(mb[11], mb[-11], FUN = ndotproduct, ppm = 40)
max(sims)
```

    ## [1] 0.7507467

We plot next a mirror plot for the two best matching spectra.

``` r

plotSpectraMirror(mb[11], mb[(which.max(sims) + 1)], ppm = 40)
```

![](MsBackendMassbank_files/figure-html/unnamed-chunk-10-1.png)

We can also retrieve the *compound* information for these two best
matching spectra. Note that this `compounds` function works only with
the `MsBackendMassbankSql` backend as it retrieves the corresponding
information from the database’s compound annotation table.

``` r

mb_match <- mb[c(11, which.max(sims) + 1)]
compounds(mb_match)
```

    ## DataFrame with 2 rows and 10 columns
    ##   compound_id     formula exactmass                 smiles
    ##     <integer> <character> <numeric>            <character>
    ## 1          31    C12H10O2   186.068 COC1=C(C=O)C2=CC=CC=..
    ## 2          45    C12H10O2   186.068 COC1=CC=C(C=O)C2=CC=..
    ##                    inchi               inchikey         cas     pubchem
    ##              <character>            <character> <character> <character>
    ## 1 InChI=1S/C12H10O2/c1.. YIQGLTKAOHRZOL-UHFFF..   1/12/5392   CID:79352
    ## 2 InChI=1S/C12H10O2/c1.. MVXMNHYVCLMLDD-UHFFF..  15971-29-6   CID:85217
    ##                                         synonym                   name
    ##                                 <CharacterList>            <character>
    ## 1 2-Methoxy-1-naphthal..,2-methoxynaphthalene.. 2-Methoxy-1-naphthal..
    ## 2 4-Methoxy-1-Naphthal..,4-methoxynaphthalene.. 4-Methoxy-1-Naphthal..

Note that the `MsBackendMassbankSql` backend does not support parallel
processing because the database connection within the backend can not be
shared across parallel processes. Any function on a `Spectra` object
that uses a `MsBackendMassbankSql` will thus (silently) disable any
parallel processing, even if the user might have passed one along to the
function using the `BPPARAM` parameter. In general, the `backendBpparam`
function can be used on any `Spectra` object to test whether its backend
supports the provided parallel processing setup (which might be helpful
for developers).

## Session information

``` r

sessionInfo()
```

    ## R Under development (unstable) (2026-01-18 r89306)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.3 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ## [1] MsCoreUtils_1.23.2       RSQLite_2.4.5            MsBackendMassbank_1.19.1
    ## [4] Spectra_1.21.1           BiocParallel_1.45.0      S4Vectors_0.49.0        
    ## [7] BiocGenerics_0.57.0      generics_0.1.4           BiocStyle_2.39.0        
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] bit_4.6.0              jsonlite_2.0.0         compiler_4.6.0        
    ##  [4] BiocManager_1.30.27    blob_1.3.0             parallel_4.6.0        
    ##  [7] cluster_2.1.8.1        jquerylib_0.1.4        systemfonts_1.3.1     
    ## [10] IRanges_2.45.0         textshaping_1.0.4      yaml_2.3.12           
    ## [13] fastmap_1.2.0          R6_2.6.1               ProtGenerics_1.39.2   
    ## [16] knitr_1.51             htmlwidgets_1.6.4      MASS_7.3-65           
    ## [19] bookdown_0.46          desc_1.4.3             DBI_1.2.3             
    ## [22] bslib_0.9.0            rlang_1.1.7            cachem_1.1.0          
    ## [25] xfun_0.56              fs_1.6.6               sass_0.4.10           
    ## [28] bit64_4.6.0-1          otel_0.2.0             memoise_2.0.1         
    ## [31] cli_3.6.5              pkgdown_2.2.0.9000     digest_0.6.39         
    ## [34] MetaboCoreUtils_1.19.1 lifecycle_1.0.5        clue_0.3-66           
    ## [37] vctrs_0.7.0            data.table_1.18.0      evaluate_1.0.5        
    ## [40] codetools_0.2-20       ragg_1.5.0             rmarkdown_2.30        
    ## [43] pkgconfig_2.0.3        tools_4.6.0            htmltools_0.5.9

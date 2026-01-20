# MS data backend for mgf files

The `MsBackendMassbank` class supports import of MS/MS spectra data from
MS/MS spectrum data from
[Massbank](https://github.com/MassBank/MassBank-data) files. After
initial import, the full MS data is kept in memory. `MsBackendMassbank`
extends the
[`Spectra::MsBackendDataFrame()`](https://rdrr.io/pkg/Spectra/man/MsBackend.html)
backend directly and supports thus the
[`Spectra::applyProcessing()`](https://rdrr.io/pkg/ProtGenerics/man/processingQueue.html)
function to make data manipulations persistent.

New objects are created with the `MsBackendMassbank` function. The
`backendInitialize` method has to be subsequently called to initialize
the object and import MS/MS data from (one or more) MassBank files.
Optional parameter `nonStop` allows to specify whether the import
returns with an error if one of the text files lacks required data, such
as `mz` and `intensity` values (default `nonStop = FALSE`), or whether
only affected file(s) is(are) skipped and a warning is shown
(`nonStop = TRUE`). Note that any other error will abort import
regardless of parameter `nonStop`.

## Usage

``` r
# S4 method for class 'MsBackendMassbank'
backendInitialize(
  object,
  files,
  metaBlocks = metaDataBlocks(),
  nonStop = FALSE,
  ...,
  BPPARAM = bpparam()
)

MsBackendMassbank()

# S4 method for class 'MsBackendMassbank'
spectraVariableMapping(object, format = c("Massbank"))

# S4 method for class 'MsBackendMassbank'
export(
  object,
  x,
  file = tempfile(),
  mapping = spectraVariableMapping(MsBackendMassbank()),
  ...
)
```

## Arguments

- object:

  Instance of `MsBackendMassbank` class.

- files:

  `character` with the (full) file name(s) of the MassBank file(s) from
  which MS/MS data should be imported.

- metaBlocks:

  `data.frame` indicating which metadata shall be imported. Default is
  [`metaDataBlocks()`](https://rformassspectrometry.github.io/MsBackendMassbank/reference/metaDataBlocks.md).

- nonStop:

  `logical(1)` whether import should be stopped if an xml file does not
  contain all required fields. Defaults to `nonStop = FALSE`.

- ...:

  Currently ignored.

- BPPARAM:

  Parameter object defining the parallel processing setup to import data
  in parallel. Defaults to `BPPARAM = bpparam()`. See
  [`BiocParallel::bpparam()`](https://rdrr.io/pkg/BiocParallel/man/register.html)
  for more information.

- format:

  for `spectraVariableMapping`: `character(1)` defining the format to be
  used. Currently only `format = "Massbank"` is supported.

- x:

  [`Spectra::Spectra()`](https://rdrr.io/pkg/Spectra/man/Spectra.html)
  object that should be exported.

- file:

  for `export`: `character(1)` defining the output file.

- mapping:

  for `export`: named `character` vector allowing to specify how fields
  from the Massbank file should be renamed. Names are supposed to be the
  spectra variable name and values of the vector the field names in the
  Massbank file. See output of
  `spectraVariableMapping(MsBackendMassbank())` for the expected format.

## Value

`backendInitialize` and `MsBackendMassbank` return an instance of
`MsBackendMassbank-class`.

## Author

Michael Witting

## Examples

``` r

## Create an MsBackendMassbank backend and import data from a test file.
fls <- dir(system.file("extdata", package = "MsBackendMassbank"),
    full.names = TRUE, pattern = "txt$")
be <- backendInitialize(MsBackendMassbank(), fls)
#> Start data import from 11 files ... 
#> done
#> Merging results ...
#> done
be
#> MsBackendMassbank with 12 spectra
#>       msLevel     rtime scanIndex
#>     <integer> <numeric> <integer>
#> 1           2        NA         1
#> 2           2        NA         1
#> 3           2    142.14         1
#> 4           2    142.14         1
#> 5           2    142.14         1
#> ...       ...       ...       ...
#> 8           2    142.14         1
#> 9           2    142.14         1
#> 10          2    143.94         1
#> 11          2    143.94         1
#> 12          2    143.94         1
#>  ... 28 more variables/columns.

be$msLevel
#>  [1] 2 2 2 2 2 2 2 2 2 2 2 2
be$intensity
#> NumericList of length 12
#> [[1]] 12461 2208 2394 40390 2816 3122 2233 ... 23807 6937 2914 6059 1871 9233
#> [[2]] 650.7 14157.3
#> [[3]] 646 980 2114 20052 1248 7628 2036 494048 75708
#> [[4]] 10186 142 142 750 138 490 126 ... 11254 14266 1478 1600 16504 1446 109762
#> [[5]] 324 184 138 3770 500 800 7214 3238 ... 2802 206 898 162 166 814 250 1840
#> [[6]] 646 980 2114 20052 1248 7628 2036 494048 75708
#> [[7]] 646 980 2114 20052 1248 7628 2036 494048 75708
#> [[8]] 10186 142 142 750 138 490 126 ... 11254 14266 1478 1600 16504 1446 109762
#> [[9]] 324 184 138 3770 500 800 7214 3238 ... 2802 206 898 162 166 814 250 1840
#> [[10]] 150 200 32 232 80 12162
#> ...
#> <2 more elements>
be$mz
#> NumericList of length 12
#> [[1]] 84.1 105.1 107.1 114.1 115.1 119.1 ... 393.3 396.3 410.3 411.3 414.3
#> [[2]] 115.0167 168.0809
#> [[3]] 74.0233 132.0807 144.0805 146.0598 ... 170.0597 188.0699 205.0965
#> [[4]] 74.0232 77.0381 86.0027 91.0539 ... 160.0947 170.0596 171.0625 188.07
#> [[5]] 53.0019 53.0383 63.0225 65.0381 ... 158.0817 159.0921 160.0755 170.06
#> [[6]] 74.0233 132.0807 144.0805 146.0598 ... 170.0597 188.0699 205.0965
#> [[7]] 74.0233 132.0807 144.0805 146.0598 ... 170.0597 188.0699 205.0965
#> [[8]] 74.0232 77.0381 86.0027 91.0539 ... 160.0947 170.0596 171.0625 188.07
#> [[9]] 53.0019 53.0383 63.0225 65.0381 ... 158.0817 159.0921 160.0755 170.06
#> [[10]] 72.0095 116.0517 117.0554 159.0935 186.0558 203.0826
#> ...
#> <2 more elements>

## Initializing a backend reading additional metadata columns/information
mb <- metaDataBlocks()
mb
#>   metadata  read
#> 1       ac FALSE
#> 2       ch FALSE
#> 3       sp FALSE
#> 4       ms FALSE
#> 5   record FALSE
#> 6       pk FALSE
#> 7  comment FALSE
mb[1, 2] <- TRUE

be <- backendInitialize(MsBackendMassbank(), fls, metaBlocks = mb)
#> Start data import from 11 files ... 
#> done
#> Merging results ...
#> done
spectraVariables(be)
#>  [1] "msLevel"                    "rtime"                     
#>  [3] "acquisitionNum"             "scanIndex"                 
#>  [5] "mz"                         "intensity"                 
#>  [7] "dataStorage"                "dataOrigin"                
#>  [9] "centroided"                 "smoothed"                  
#> [11] "polarity"                   "precScanNum"               
#> [13] "precursorMz"                "precursorIntensity"        
#> [15] "precursorCharge"            "collisionEnergy"           
#> [17] "isolationWindowLowerMz"     "isolationWindowTargetMz"   
#> [19] "isolationWindowUpperMz"     "acquistionNum"             
#> [21] "accession"                  "name"                      
#> [23] "smiles"                     "exactmass"                 
#> [25] "formula"                    "inchi"                     
#> [27] "cas"                        "inchikey"                  
#> [29] "adduct"                     "splash"                    
#> [31] "title"                      "instrument"                
#> [33] "instrument_type"            "ms_ms_type"                
#> [35] "ms_cap_voltage"             "ms_col_gas"                
#> [37] "ms_desolv_gas_flow"         "ms_desolv_temp"            
#> [39] "ms_frag_mode"               "ms_ionization"             
#> [41] "ms_ionization_energy"       "ms_laser"                  
#> [43] "ms_matrix"                  "ms_mass_accuracy"          
#> [45] "ms_mass_range"              "ms_reagent_gas"            
#> [47] "ms_resolution"              "ms_scan_setting"           
#> [49] "ms_source_temp"             "ms_kinetic_energy"         
#> [51] "ms_electron_current"        "ms_reaction_time"          
#> [53] "chrom_carrier_gas"          "chrom_column"              
#> [55] "chrom_column_temp"          "chrom_column_temp_gradient"
#> [57] "chrom_flow_gradient"        "chrom_flow_rate"           
#> [59] "chrom_inj_temp"             "chrom_inj_temp_gradient"   
#> [61] "chrom_rti_kovats"           "chrom_rti_lee"             
#> [63] "chrom_rti_naps"             "chrom_rti_uoa"             
#> [65] "chrom_rti_uoa_pred"         "chrom_rt"                  
#> [67] "chrom_rt_uoa_pred"          "chrom_solvent"             
#> [69] "chrom_transfer_temp"        "ims_instrument_type"       
#> [71] "ims_drift_gas"              "ims_drift_time"            
#> [73] "ims_ccs"                    "general_conc"              
be$instrument
#>  [1] "Bruker maXis ESI-QTOF"                  
#>  [2] "LTQ Orbitrap XL Thermo Scientific"      
#>  [3] "maXis plus UHR-ToF-MS, Bruker Daltonics"
#>  [4] "maXis plus UHR-ToF-MS, Bruker Daltonics"
#>  [5] "maXis plus UHR-ToF-MS, Bruker Daltonics"
#>  [6] "maXis plus UHR-ToF-MS, Bruker Daltonics"
#>  [7] "maXis plus UHR-ToF-MS, Bruker Daltonics"
#>  [8] "maXis plus UHR-ToF-MS, Bruker Daltonics"
#>  [9] "maXis plus UHR-ToF-MS, Bruker Daltonics"
#> [10] "maXis plus UHR-ToF-MS, Bruker Daltonics"
#> [11] "maXis plus UHR-ToF-MS, Bruker Daltonics"
#> [12] "maXis plus UHR-ToF-MS, Bruker Daltonics"
```

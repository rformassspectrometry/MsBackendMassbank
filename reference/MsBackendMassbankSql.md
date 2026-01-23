# MS backend accessing the MassBank MySQL database

The `MsBackendMassbankSql` provides access to mass spectrometry data
from [MassBank](https://massbank.eu/MassBank/) by directly accessing its
MySQL/MariaDb database. In addition it supports adding new spectra
variables or *locally* changing spectra variables provided by MassBank
(without changing the original values in the database).

Note that `MsBackendMassbankSql` requires a local installation of the
MassBank database since direct database access is not supported for the
*main* MassBank instance.

Also, some of the fields in the MassBank database are not directly
compatible with `Spectra`, such as the *collision energy* which is not
available as a numeric value. The collision energy as available in
MassBank is reported as spectra variable `"collision_energy_text"`.
Also, precursor m/z values reported for some spectra can not be
converted to a `numeric` and hence `NA` is reported with the spectra
variable `precursorMz` for these spectra. The variable
`"precursor_mz_text"` can be used to get the *original* precursor m/z
reported in MassBank.

Finally, `MsBackendMassbankSql` does **not support** parallel processing
because the database connection stored within the object can not be
shared acrcoss parallel processes. All functions on `Spectra` objects
with a `MsBackendMassbankSql` will (silently) disable parallel
processing even if the user provides a dedicated parallel processing
setup with the `BPPARAM` parameter.

## Usage

``` r
MsBackendMassbankSql()

# S4 method for class 'MsBackendMassbankSql'
backendInitialize(object, dbcon, ...)

# S4 method for class 'MsBackendMassbankSql'
peaksData(object, columns = peaksVariables(object))

# S4 method for class 'MsBackendMassbankSql'
dataStorage(object)

# S4 method for class 'MsBackendMassbankSql'
intensity(object) <- value

# S4 method for class 'MsBackendMassbankSql'
mz(object) <- value

# S4 method for class 'MsBackendMassbankSql'
reset(object)

# S4 method for class 'MsBackendMassbankSql'
spectraData(object, columns = spectraVariables(object))

# S4 method for class 'MsBackendMassbankSql'
spectraNames(object)

# S4 method for class 'MsBackendMassbankSql'
spectraNames(object) <- value

# S4 method for class 'MsBackendMassbankSql'
tic(object, initial = TRUE)

# S4 method for class 'MsBackendMassbankSql'
x[i, j, ..., drop = FALSE]

# S4 method for class 'MsBackendMassbankSql,ANY'
extractByIndex(object, i)

# S4 method for class 'Spectra'
compounds(object, ...)

# S4 method for class 'MsBackendMassbankSql'
compounds(object, ...)

# S4 method for class 'MsBackendMassbankSql'
x$name <- value

# S4 method for class 'MsBackendMassbankSql'
precScanNum(object)

# S4 method for class 'MsBackendMassbankSql'
backendBpparam(object, BPPARAM = bpparam())
```

## Arguments

- object:

  Object extending `MsBackendMassbankSql`.

- dbcon:

  For `backendInitialize,MsBackendMassbankSql`: SQL database connection
  to the MassBank (MariaDb) database.

- ...:

  Additional arguments.

- columns:

  For `spectraData()` accessor: optional `character` with column names
  (spectra variables) that should be included in the returned
  `DataFrame`. By default, all columns are returned. For `peaksData`
  accessor: optional `character` with requested columns in the
  individual `matrix` of the returned `list`. Use
  `peaksVariables(object)` for supported columns.

- value:

  replacement value for `<-` methods. See individual method description
  or expected data type.

- initial:

  For `tic`: `logical(1)` whether the initially reported total ion
  current should be reported, or whether the total ion current should be
  (re)calculated on the actual data (`initial = FALSE`).

- x:

  Object extending `MsBackendMassbankSql`.

- i:

  For `[`: `integer`, `logical` or `character` to subset the object.

- j:

  For `[`: not supported.

- drop:

  For `[`: not considered.

- name:

  name of the variable to replace for `<-` methods. See individual
  method description or expected data type.

- BPPARAM:

  for `backendBpparam()`: `BiocParallel` parallel processing setup. See
  [`BiocParallel::bpparam()`](https://rdrr.io/pkg/BiocParallel/man/register.html)
  for more information.

- spectraVariables:

  For `selectSpectraVariables()`: `character` with the names of the
  spectra variables to which the backend should be subsetted.

## Value

See documentation of respective function.

## Supported Backend functions

The following functions are supported by the `MsBackendMassbankSql`.

- `[`: subset the backend. Only subsetting by element (*row*/`i`) is
  allowed

- `$`, `$<-`: access or set/add a single spectrum variable (column) in
  the backend.

- `acquisitionNum()`: returns the acquisition number of each spectrum.
  Returns an `integer` of length equal to the number of spectra (with
  `NA_integer_` if not available).

- `peaksData()` returns a `list` with the spectras' peak data. The
  length of the list is equal to the number of spectra in `object`. Each
  element of the list is a `matrix` with columns `"mz"` and
  `"intensity"`. For an empty spectrum, a `matrix` with 0 rows and two
  columns (named `mz` and `intensity`) is returned. Parameter `columns`
  allows to select which peaks variables to return, but supports
  currently only `"mz"` and `"intensity"`.

- `backendBpparam()`: whether the backend supports parallel processing.
  Takes a `MsBackendMassbankSql` and a parallel processing setup (see
  [`BiocParallel::bpparam()`](https://rdrr.io/pkg/BiocParallel/man/register.html)
  for details) as input and **always** returns a
  [`BiocParallel::SerialParam()`](https://rdrr.io/pkg/BiocParallel/man/SerialParam-class.html).
  This function can be used to test whether a provided parallel
  processing setup is supported by the backend and returns the supported
  setup.

- `backendInitialize()`: initialises the backend by retrieving the IDs
  of all spectra in the database. Parameter `dbcon` with the connection
  to the MassBank MySQL database is required.

- `dataOrigin()`: gets a `character` of length equal to the number of
  spectra in `object` with the *data origin* of each spectrum. This
  could e.g. be the mzML file from which the data was read.

- `dataStorage()`: returns `"<MassBank>"` for all spectra.

- `centroided()`, `centroided<-`: gets or sets the centroiding
  information of the spectra. `centroided()` returns a `logical` vector
  of length equal to the number of spectra with `TRUE` if a spectrum is
  centroided, `FALSE` if it is in profile mode and `NA` if it is
  undefined. See also `isCentroided()` for estimating from the spectrum
  data whether the spectrum is centroided. `value` for `centroided<-` is
  either a single `logical` or a `logical` of length equal to the number
  of spectra in `object`.

- `collisionEnergy()`, `collisionEnergy<-`: gets or sets the collision
  energy for all spectra in `object`. `collisionEnergy` returns a
  `numeric` with length equal to the number of spectra (`NA_real_` if
  not present/defined), `collisionEnergy<-` takes a `numeric` of length
  equal to the number of spectra in `object`. Note that the collision
  energy description from MassBank are provided as spectra variable
  `"collisionEnergyText"`.

- `intensity()`: gets the intensity values from the spectra. Returns a
  [`IRanges::NumericList()`](https://rdrr.io/pkg/IRanges/man/AtomicList-class.html)
  of `numeric` vectors (intensity values for each spectrum). The length
  of the `list` is equal to the number of `spectra` in `object`.

- `ionCount()`: returns a `numeric` with the sum of intensities for each
  spectrum. If the spectrum is empty (see `isEmpty()`), `NA_real_` is
  returned.

- `isCentroided()`: a heuristic approach assessing if the spectra in
  `object` are in profile or centroided mode. The function takes the
  `qtl` th quantile top peaks, then calculates the difference between
  adjacent m/z value and returns `TRUE` if the first quartile is greater
  than `k`. (See `Spectra:::.isCentroided` for the code.)

- `isEmpty()`: checks whether a spectrum in `object` is empty (i.e. does
  not contain any peaks). Returns a `logical` vector of length equal
  number of spectra.

- `isolationWindowLowerMz()`, `isolationWindowLowerMz<-`: gets or sets
  the lower m/z boundary of the isolation window.

- `isolationWindowTargetMz()`, `isolationWindowTargetMz<-`: gets or sets
  the target m/z of the isolation window.

- `isolationWindowUpperMz()`, `isolationWindowUpperMz<-`: gets or sets
  the upper m/z boundary of the isolation window.

- `isReadOnly()`: returns a `logical(1)` whether the backend is *read
  only* or does allow also to write/update data.

- [`length()`](https://rdrr.io/r/base/length.html): returns the number
  of spectra in the object.

- [`lengths()`](https://rdrr.io/r/base/lengths.html): gets the number of
  peaks (m/z-intensity values) per spectrum. Returns an `integer` vector
  (length equal to the number of spectra). For empty spectra, `0` is
  returned.

- `msLevel()`: gets the spectra's MS level. Returns an `integer` vector
  (of length equal to the number of spectra) with the MS level for each
  spectrum (or `NA_integer_` if not available).

- `mz()`: gets the mass-to-charge ratios (m/z) from the spectra. Returns
  a
  [`IRanges::NumericList()`](https://rdrr.io/pkg/IRanges/man/AtomicList-class.html)
  or length equal to the number of spectra, each element a `numeric`
  vector with the m/z values of one spectrum.

- `polarity()`, `polarity<-`: gets or sets the polarity for each
  spectrum. `polarity` returns an `integer` vector (length equal to the
  number of spectra), with `0` and `1` representing negative and
  positive polarities, respectively. `polarity<-` expects an integer
  vector of length 1 or equal to the number of spectra.

- `precursorCharge90`, `precursorIntensity()`, `precursorMz()`,
  `precScanNum()`, `precAcquisitionNum()`: get the charge (`integer`),
  intensity (`numeric`), m/z (`numeric`), scan index (`integer`) and
  acquisition number (`interger`) of the precursor for MS level 2 and
  above spectra from the object. Returns a vector of length equal to the
  number of spectra in `object`. `NA` are reported for MS1 spectra of if
  no precursor information is available.

- `reset()`: restores the backend to its original state, i.e. deletes
  all locally modified data and reinitializes the backend to the full
  data available in the database.

- `rtime()`, `rtime<-`: gets or sets the retention times for each
  spectrum (in seconds). `rtime` returns a `numeric` vector (length
  equal to the number of spectra) with the retention time for each
  spectrum. `rtime<-` expects a numeric vector with length equal to the
  number of spectra.

- `scanIndex()`: returns an `integer` vector with the *scan index* for
  each spectrum. This represents the relative index of the spectrum
  within each file. Note that this can be different to the
  `acquisitionNum` of the spectrum which is the index of the spectrum as
  reported in the mzML file.

- `selectSpectraVariables()`: reduces the information within the backend
  to the selected spectra variables.

- `smoothed()`,`smoothed<-`: gets or sets whether a spectrum is
  *smoothed*. `smoothed` returns a `logical` vector of length equal to
  the number of spectra. `smoothed<-` takes a `logical` vector of length
  1 or equal to the number of spectra in `object`.

- `spectraData()`: gets general spectrum metadata (annotation, also
  called header). `spectraData` returns a `DataFrame`. Note that
  replacing the spectra data with `spectraData<-` is not supported.

- `spectraNames()`: returns a `character` vector with the names of the
  spectra in `object`.

- `spectraVariables()`: returns a `character` vector with the available
  spectra variables (columns, fields or attributes) available in
  `object`. This should return **all** spectra variables which are
  present in `object`, also `"mz"` and `"intensity"` (which are by
  default not returned by the `spectraVariables,Spectra` method).

- `tic()`: gets the total ion current/count (sum of signal of a
  spectrum) for all spectra in `object`. By default, the value reported
  in the original raw data file is returned. For an empty spectrum,
  `NA_real_` is returned.

## Not supported Backend functions

The following functions are not supported by the `MsBackendMassbankSql`
since the original data can not be changed.

`backendMerge()`, `export()`, `filterDataStorage()`,
`filterPrecursorScan()`, `peaksData<-`, `filterAcquisitionNum()`,
`intensity<-`, `mz<-`, `precScanNum()`, `spectraData<-`,
`spectraNames<-`.

## Retrieving compound annotations for spectra

While compound annotations are also provided *via* the
`spectraVariables()` of the backend, it would also be possible to use
the `compounds` function on a `Spectra` object (that uses a
`MsBackendMassbankSql` backend) to retrieve compound annotations for the
specific spectra.

## Author

Johannes Rainer

## Examples

``` r

## Create a connection to a database with MassBank data - in the present
## example we connect to a tiny SQLite database bundled in this package
## as public access to the MassBank MySQL is not (yet) supported. See the
## vignette for more information on how to install MassBank locally and
## enable MySQL database connections
library(RSQLite)
con <- dbConnect(SQLite(), system.file("sql", "minimassbank.sqlite",
    package = "MsBackendMassbank"))

## Given that we have the connection to a MassBank databas we can
## initialize the backend:
be <- backendInitialize(MsBackendMassbankSql(), dbcon = con)
be
#> MsBackendMassbankSql with 70 spectra
#>       msLevel precursorMz  polarity
#>     <integer>   <numeric> <integer>
#> 1           2         506         0
#> 2          NA          NA         1
#> 3          NA          NA         0
#> 4          NA          NA         1
#> 5          NA          NA         0
#> ...       ...         ...       ...
#> 66          2     185.028         0
#> 67          2     455.290         1
#> 68          2     253.051         0
#> 69          2     358.238         1
#> 70          2     256.170         1
#>  ... 41 more variables/columns.
#>  Use  'spectraVariables' to list all of them.

## Access MS level
msLevel(be)
#>  [1]  2 NA NA NA NA  2 NA NA NA NA NA NA NA NA NA NA NA NA NA NA  2  2  2  2  2
#> [26]  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2
#> [51]  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2
be$msLevel
#>  [1]  2 NA NA NA NA  2 NA NA NA NA NA NA NA NA NA NA NA NA NA NA  2  2  2  2  2
#> [26]  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2
#> [51]  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2

## Access m/z values
be$mz
#> NumericList of length 70
#> [[1]] 146.760803 158.863541 174.988785 ... 470.057434 487.989319 585.88446
#> [[2]] 22.99 23.07 23.19 38.98 53.15 60.08 ... 391.18 413.22 414.22 429.16 495.2
#> [[3]] 108.099566 152.005378 153.01824 ... 341.031964 409.06729 499.103447
#> [[4]] 137.0227 138.025605 139.034602 ... 304.271886 316.303106 352.234228
#> [[5]] 112.977759 205.053323 208.041128 ... 449.137152 671.1778 673.200866
#> [[6]] 78.893929 183.011719 193.957916 ... 328.091949 408.010193 426.021942
#> [[7]] 167.03389 168.040758 333.079965 ... 357.073039 365.040698 373.073473
#> [[8]] 195.09167 196.095033
#> [[9]] 108.096149 152.004656 153.01824 ... 359.997563 409.065314 491.043775
#> [[10]] 1278.12 1279.11 1279.18 1306.21 ... 2064.29 2091.32 2092.3 2136.33
#> ...
#> <60 more elements>

## Access the full spectra data (including m/z and intensity values)
spectraData(be)
#> DataFrame with 70 rows and 44 columns
#>       msLevel     rtime acquisitionNum scanIndex                             mz
#>     <integer> <numeric>      <integer> <integer>                  <NumericList>
#> 1           2        NA             NA        NA    146.761,158.864,174.989,...
#> 2          NA        NA             NA        NA          22.99,23.07,23.19,...
#> 3          NA        NA             NA        NA    108.100,152.005,153.018,...
#> 4          NA        NA             NA        NA    137.023,138.026,139.035,...
#> 5          NA        NA             NA        NA    112.978,205.053,208.041,...
#> ...       ...       ...            ...       ...                            ...
#> 66          2        NA             NA        NA                        185.028
#> 67          2        NA             NA        NA    58.0653,77.0386,79.0543,...
#> 68          2        NA             NA        NA    132.021,133.029,135.009,...
#> 69          2        NA             NA        NA    53.0388,54.0340,55.0544,...
#> 70          2        NA             NA        NA  88.0757,151.0540,152.0620,...
#>                                intensity dataStorage  dataOrigin centroided
#>                            <NumericList> <character> <character>  <logical>
#> 1     0.980184,938.145447, 56.718353,...  <MassBank>          NA         NA
#> 2                     3173,2266, 880,...  <MassBank>          NA         NA
#> 3                 2809,119550, 91935,...  <MassBank>          NA         NA
#> 4               157703, 23566,  5190,...  <MassBank>          NA         NA
#> 5                     5686, 927, 756,...  <MassBank>          NA         NA
#> ...                                  ...         ...         ...        ...
#> 66                              34967548  <MassBank>          NA         NA
#> 67           9143284, 341391,5278218,...  <MassBank>          NA         NA
#> 68         82825.8,290504.6,301494.3,...  <MassBank>          NA         NA
#> 69         1568742,19224644, 4676170,...  <MassBank>          NA         NA
#> 70          678846,  403851,11111366,...  <MassBank>          NA         NA
#>      smoothed  polarity precScanNum precursorMz precursorIntensity
#>     <logical> <integer>   <integer>   <numeric>          <numeric>
#> 1          NA         0          NA         506                408
#> 2          NA         1          NA          NA                 NA
#> 3          NA         0          NA          NA                 NA
#> 4          NA         1          NA          NA                 NA
#> 5          NA         0          NA          NA                 NA
#> ...       ...       ...         ...         ...                ...
#> 66         NA         0          NA     185.028            185.028
#> 67         NA         1          NA     455.290            455.290
#> 68         NA         0          NA     253.051            253.051
#> 69         NA         1          NA     358.238            358.237
#> 70         NA         1          NA     256.170            256.169
#>     precursorCharge collisionEnergy isolationWindowLowerMz
#>           <integer>       <numeric>              <numeric>
#> 1                NA              NA                     NA
#> 2                NA              NA                     NA
#> 3                NA              NA                     NA
#> 4                NA              NA                     NA
#> 5                NA              NA                     NA
#> ...             ...             ...                    ...
#> 66               NA              NA                     NA
#> 67               NA              NA                     NA
#> 68               NA              NA                     NA
#> 69               NA              NA                     NA
#> 70               NA              NA                     NA
#>     isolationWindowTargetMz isolationWindowUpperMz spectrum_id
#>                   <numeric>              <numeric> <character>
#> 1                        NA                     NA    MCH00020
#> 2                        NA                     NA    MCH00002
#> 3                        NA                     NA    MCH00006
#> 4                        NA                     NA    MCH00009
#> 5                        NA                     NA    MCH00014
#> ...                     ...                    ...         ...
#> 66                       NA                     NA    SM834451
#> 67                       NA                     NA    SM855703
#> 68                       NA                     NA    SM882851
#> 69                       NA                     NA    SM849703
#> 70                       NA                     NA    SM858801
#>              spectrum_name                   date                authors
#>                <character>            <character>            <character>
#> 1   Adenosine 5'-triphos.. 2016.01.19 (Created .. Yoshikuni K, Tajiri ..
#> 2   Gentisic acid; MALDI.. 2016.01.19 (Created .. Wada Y, Osaka Medica..
#> 3   Gentisic acid; MALDI.. 2016.01.19 (Created .. Wada Y, Osaka Medica..
#> 4   5-Methoxysalicylic a.. 2016.01.19 (Created .. Wada Y, Osaka Medica..
#> 5   Sinapic acid; MALDI-.. 2016.01.19 (Created .. Wada Y, Osaka Medica..
#> ...                    ...                    ...                    ...
#> 66  m-Xylene-4-sulfonic ..             2016.12.12 Krauss M, Schymanski..
#> 67  Verapamil; LC-ESI-QF..             2016.12.12 Krauss M, Schymanski..
#> 68  Daidzein; LC-ESI-QFT..             2016.12.12 Krauss M, Schymanski..
#> 69  Oxybutynin; LC-ESI-Q..             2016.12.12 Krauss M, Schymanski..
#> 70  Diphenhydramine; LC-.. 2019.11.20 (Created .. Krauss M, Schymanski..
#>         license              copyright            publication
#>     <character>            <character>            <character>
#> 1      CC BY-SA                     NA                     NA
#> 2      CC BY-SA                     NA                     NA
#> 3      CC BY-SA                     NA                     NA
#> 4      CC BY-SA                     NA                     NA
#> 5      CC BY-SA                     NA                     NA
#> ...         ...                    ...                    ...
#> 66        CC BY Copyright (C) 2016 U.. Schymanski, E. L.; R..
#> 67        CC BY Copyright (C) 2016 U.. Schymanski, E. L.; R..
#> 68        CC BY Copyright (C) 2016 U.. Schymanski, E. L.; R..
#> 69        CC BY Copyright (C) 2016 U.. Schymanski, E. L.; R..
#> 70        CC BY Copyright (C) 2016 U.. Schymanski, E. L.; R..
#>                     splash compound_id      adduct  ionization
#>                <character>   <integer> <character> <character>
#> 1   splash10-0a4i-000090..           1      [M-H]-          NA
#> 2   splash10-004r-191300..           2          NA          NA
#> 3   splash10-0udi-090200..           3          NA          NA
#> 4   splash10-0fri-093000..           4          NA          NA
#> 5   splash10-006t-019080..           5          NA          NA
#> ...                    ...         ...         ...         ...
#> 66  splash10-000i-090000..          66      [M-H]-         ESI
#> 67  splash10-066r-090030..          67      [M+H]+         ESI
#> 68  splash10-0udi-009000..          68      [M-H]-         ESI
#> 69  splash10-0a4i-860900..          69      [M+H]+         ESI
#> 70  splash10-014i-090000..          70      [M+H]+         ESI
#>     ionization_voltage fragmentation_mode             instrument
#>            <character>        <character>            <character>
#> 1                   NA                 NA LTQ Orbitrap XL, The..
#> 2                   NA                 NA Voyager DE-PRO, Appl..
#> 3                   NA                 NA Voyager DE-PRO, Appl..
#> 4                   NA                 NA Voyager DE-PRO, Appl..
#> 5                   NA                 NA Voyager DE-PRO, Appl..
#> ...                ...                ...                    ...
#> 66                  NA                HCD Q Exactive Plus Orbi..
#> 67                  NA                HCD Q Exactive Plus Orbi..
#> 68                  NA                HCD Q Exactive Plus Orbi..
#> 69                  NA                HCD Q Exactive Plus Orbi..
#> 70                  NA                HCD Q Exactive Plus Orbi..
#>     instrument_type       formula exactmass                 smiles
#>         <character>   <character> <numeric>            <character>
#> 1       LC-ESI-ITFT C10H16N5O13P3   506.996 Nc(n3)c(n2)c(nc3)n(c..
#> 2         MALDI-TOF        C7H6O4   154.027 Oc(c1)cc(C(O)=O)c(O)c1
#> 3         MALDI-TOF        C7H6O4   154.027 Oc(c1)cc(C(O)=O)c(O)c1
#> 4         MALDI-TOF        C8H8O4   168.042 COc(c1)cc(C(O)=O)c(O..
#> 5         MALDI-TOF      C11H12O5   224.068 COC1=CC(=CC(=C1O)OC)..
#> ...             ...           ...       ...                    ...
#> 66       LC-ESI-QFT      C8H10O3S   186.035 Cc1ccc(c(c1)C)S(=O)(..
#> 67       LC-ESI-QFT    C27H38N2O4   454.283 COc1ccc(CCN(C)CCCC(C..
#> 68       LC-ESI-QFT      C15H10O4   254.058 OC1=CC=C(C=C1)C1=COC..
#> 69       LC-ESI-QFT     C22H31NO3   357.230 CCN(CC)CC#CCOC(=O)C(..
#> 70       LC-ESI-QFT      C17H21NO   255.162 CN(C)CCOC(c1ccccc1)c..
#>                      inchi               inchikey         cas     pubchem
#>                <character>            <character> <character> <character>
#> 1   InChI=1S/C10H16N5O13.. ZKHQWZAMYRWXGA-KQYNX..     56-65-5    CID:5957
#> 2   InChI=1S/C7H6O4/c8-4.. WXTMDXOMEHJXQO-UHFFF..          NA    CID:3469
#> 3   InChI=1S/C7H6O4/c8-4.. WXTMDXOMEHJXQO-UHFFF..          NA    CID:3469
#> 4   InChI=1S/C8H8O4/c1-1.. IZZIWIAOVZOBLF-UHFFF..          NA   CID:75787
#> 5   InChI=1S/C11H12O5/c1.. PCMORTLOPMLEFB-ONEGZ..    530-59-6  CID:637775
#> ...                    ...                    ...         ...         ...
#> 66  InChI=1S/C8H10O3S/c1.. CHZLVSBMXZSPNN-UHFFF..     88-61-9    CID:6938
#> 67  InChI=1S/C27H38N2O4/.. SGTNSNPWRIOYBX-UHFFF..     52-53-9    CID:2520
#> 68  InChI=1S/C15H10O4/c1.. ZQSIJRDFPHDXIC-UHFFF..    486-66-8 CID:5281708
#> 69  InChI=1S/C22H31NO3/c.. XIQVNETUBQGFHX-UHFFF..   5633-20-5    CID:4634
#> 70  InChI=1S/C17H21NO/c1.. ZZVUWRFHKOJYTH-UHFFF..     58-73-1    CID:3100
#>                                                        synonym
#>                                                <CharacterList>
#> 1                                   Adenosine 5'-triphos..,ATP
#> 2                         2,5-Dihydroxybenzoic..,Gentisic acid
#> 3                     2,5-Dihydroxybenzoic..,DHB,Gentisic acid
#> 4                      5-Methoxysalicylic a..,DHBs (Super DHB)
#> 5   Sinapic acid,3,5-dimethoxy-4-hydr..,3-(4-hydroxy-3,5-dim..
#> ...                                                        ...
#> 66               m-Xylene-4-sulfonic ..,2,4-dimethylbenzenes..
#> 67                            Verapamil,2-(3,4-dimethoxyphen..
#> 68                             Daidzein,7-hydroxy-3-(4-hydro..
#> 69                           Oxybutynin,4-(diethylamino)but-..
#> 70                      Diphenhydramine,2-benzhydryloxy-N,N-..
#>     precursor_mz_text          compound_name
#>           <character>            <character>
#> 1          506.000000 Adenosine 5'-triphos..
#> 2                  NA 2,5-Dihydroxybenzoic..
#> 3                  NA 2,5-Dihydroxybenzoic..
#> 4                  NA 5-Methoxysalicylic a..
#> 5                  NA           Sinapic acid
#> ...               ...                    ...
#> 66           185.0278 m-Xylene-4-sulfonic ..
#> 67           455.2904              Verapamil
#> 68           253.0506               Daidzein
#> 69           358.2377             Oxybutynin
#> 70           256.1696        Diphenhydramine

## Add a new spectra variable
be$new_variable <- "b"
be$new_variable
#>  [1] "b" "b" "b" "b" "b" "b" "b" "b" "b" "b" "b" "b" "b" "b" "b" "b" "b" "b" "b"
#> [20] "b" "b" "b" "b" "b" "b" "b" "b" "b" "b" "b" "b" "b" "b" "b" "b" "b" "b" "b"
#> [39] "b" "b" "b" "b" "b" "b" "b" "b" "b" "b" "b" "b" "b" "b" "b" "b" "b" "b" "b"
#> [58] "b" "b" "b" "b" "b" "b" "b" "b" "b" "b" "b" "b" "b"

## Subset the backend
be_sub <- be[c(3, 1)]

spectraNames(be)
#>  [1] "MCH00020" "MCH00002" "MCH00006" "MCH00009" "MCH00014" "MCH00019"
#>  [7] "MCH00012" "MCH00003" "MCH00016" "MCH00001" "MCH00018" "MCH00007"
#> [13] "MCH00013" "MCH00015" "MCH00008" "MCH00004" "MCH00005" "MCH00011"
#> [19] "MCH00017" "MCH00010" "SM864952" "SM823851" "SM868004" "SM862204"
#> [25] "SM834302" "SM810352" "SM847101" "SM861203" "SM839603" "SM829401"
#> [31] "SM814501" "SM858902" "SM850104" "SM854403" "SM836003" "SM821703"
#> [37] "SM880602" "SM868101" "SM851552" "SM804503" "SM848503" "SM842251"
#> [43] "SM817703" "SM847302" "SM805801" "SM820101" "SM816753" "SM805001"
#> [49] "SM838051" "SM812904" "SM840401" "SM812701" "SM861901" "SM874401"
#> [55] "SM838502" "SM857001" "SM858203" "SM847702" "SM883903" "SM818403"
#> [61] "SM871101" "SM842156" "SM841901" "SM854306" "SM878551" "SM834451"
#> [67] "SM855703" "SM882851" "SM849703" "SM858801"
spectraNames(be_sub)
#> [1] "MCH00006" "MCH00020"
```

# Metadata blocks to be read

`metaDataBlocks()` allows to define the metadata *blocks* to imported
from the MassBank record files.

## Usage

``` r
metaDataBlocks(
  ac = FALSE,
  ch = FALSE,
  sp = FALSE,
  ms = FALSE,
  record = FALSE,
  pk = FALSE,
  comment = FALSE
)
```

## Arguments

- ac:

  `logical(1)`: read and parse the `"AC$"` entries. These include
  information on the mass spectrometry instrument, ionization applied,
  fragmentation mode etc.

- ch:

  `logical(1)`: read and parse the `"CH$"` entries with compound related
  information/annotation, such as IDs to external databases.

- sp:

  `logical(1)`: read and parse the `"SP$"` entries with sample related
  information.

- ms:

  `logical(1)`: read and parse the `"MS$"` entries with mass
  spectrometry related information and data processing applied.

- record:

  `logical(1)`: read and parse *record* related information such as the
  authors, the date, license etc.

- pk:

  `logical(1)`: read the number of peaks.

- comment:

  `logical(1)`: read optional comments.

## Value

A `data.frame` with information which metadata blocks should be mported.

## Author

Michael Witting

## Examples

``` r

metaDataBlocks()
#>   metadata  read
#> 1       ac FALSE
#> 2       ch FALSE
#> 3       sp FALSE
#> 4       ms FALSE
#> 5   record FALSE
#> 6       pk FALSE
#> 7  comment FALSE
```

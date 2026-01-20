# Metadata blocks to be read

`metaDataBlocks` returns a `data.frame` with the MassBank metadata
blocks and whether they should be imported by default from the MassBank
text files.

## Usage

``` r
metaDataBlocks()
```

## Value

A `data.frame` with metadata blocks.

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

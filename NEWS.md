# MsBackendMassbank 1.1

## Changes in 1.1.2

- Fix bug in `show,MsBackendMassbankSql`.

## Changes in 1.1.1

- Cache precursor m/z to allow faster queries for spectral matching.

# MsBackendMassbank 0.99

## Changes in 0.99.4

- Fix an issue in which `selectSpectraVariables,MsBackendMassbankSql` would fail
  if subsetted to a single variable/column in `@localData` [issue
  #29](https://github.com/rformassspectrometry/MsBackendMassbank/issues/29)

## Changes in 0.99.3

- Drop names on the `peaksData`.

## Changes in 0.99.2

- Minor updates and changes as requested during package review process.

## Changes in 0.99.1

- Directly call package internal functions using `:::` if used within `bplapply`
  to avoid errors on Windows.

## Changes in 0.99.0

- Prepare for Bioconductor submission.
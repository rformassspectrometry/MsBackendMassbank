# MsBackendMassbank 1.7

## Changes in 1.7.2

- Avoid export of retention time data if not available.

## Changes in 1.7.1

- Run the full test suite from the `Spectra` package to validate the `MsBackend`
  implementations.

# MsBackendMassbank 1.3

## Changes in 1.3.5

- Add parameter `columns` to `peaksData`.

## Changes in 1.3.4

- Import and use `spectraVariableMapping` method from `Spectra`.

## Changes in 1.3.3

- Map comment to spectra variable.

## Changes in 1.3.2

- Use in addition unit tests from the `Spectra` package.
- Add `filterPrecursorMzValues` and `filterPrecursorMzRange`.

## Changes in 1.3.1

- `MsBackendMassbankSql` extends `Spectra::MsBackendCached` to re-use the
  general caching mechanism provided by that backend.

# MsBackendMassbank 1.1

## Changes in 1.1.4

- Fix wrong database column name for collision energy.

## Changes in 1.1.3

- Change SQL queries to increase performance.

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
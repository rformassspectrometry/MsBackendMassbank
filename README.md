# Mass Spectrometry Data Backend for MassBank Files

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check-bioc](https://github.com/RforMassSpectrometry/MsBackendMassbank/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/RforMassSpectrometry/MsBackendMassbank/actions?query=workflow%3AR-CMD-check-bioc)
[![codecov.io](http://codecov.io/github/RforMassSpectrometry/MsBackendMassbank/coverage.svg?branch=master)](http://codecov.io/github/RforMassSpectrometry/MsBackendMassbank?branch=master)
[![license](https://img.shields.io/badge/license-Artistic--2.0-brightgreen.svg)](https://opensource.org/licenses/Artistic-2.0)

<img
src="https://raw.githubusercontent.com/rformassspectrometry/stickers/master/MsBackendMassbank/MsBackendMassbank.png"
height="150">

The `MsBackendMassbank` package provides functionality to import and handle
MS/MS spectrum data from [Massbank](https://github.com/MassBank/MassBank-data)
files or to directly access a MassBank SQL database.  The package defines the
`MsBackendMassbank` and `MsBackendMassbankSql` backends which can be used to
import and use MS2 spectrum data from mgf files respectively MySQL databases
with the [Spectra](https://github.com/rformassspectrometry/Spectra) R package.

For more information see the package
[homepage](https://github.com/RforMassSpectrometry/MsBackendMassbank).

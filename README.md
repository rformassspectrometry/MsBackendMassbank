# Mass Spectrometry Data Backend for MassBank Files

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check-bioc](https://github.com/RforMassSpectrometry/MsBackendMassbank/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/RforMassSpectrometry/MsBackendMassbank/actions?query=workflow%3AR-CMD-check-bioc)
[![codecov](https://codecov.io/gh/rformassspectrometry/MsBackendMassbank/branch/main/graph/badge.svg?token=OZ4Z5VN50J)](https://codecov.io/gh/rformassspectrometry/MsBackendMassbank)
[![license](https://img.shields.io/badge/license-Artistic--2.0-brightgreen.svg)](https://opensource.org/licenses/Artistic-2.0)
[![years in bioc](http://bioconductor.org/shields/years-in-bioc/MsBackendMassbank.svg)](https://bioconductor.org/packages/release/bioc/html/MsBackendMassbank.html)
[![Ranking by downloads](http://bioconductor.org/shields/downloads/release/MsBackendMassbank.svg)](https://bioconductor.org/packages/stats/bioc/MsBackendMassbank/)
[![build release](http://bioconductor.org/shields/build/release/bioc/MsBackendMassbank.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/MsBackendMassbank/)
[![build devel](http://bioconductor.org/shields/build/devel/bioc/MsBackendMassbank.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/MsBackendMassbank/)

The `MsBackendMassbank` package provides functionality to import and handle
MS/MS spectrum data from [Massbank](https://github.com/MassBank/MassBank-data)
files or to directly access a MassBank SQL database.  The package defines the
`MsBackendMassbank` and `MsBackendMassbankSql` backends which can be used to
import and use MS2 spectrum data from mgf files respectively MySQL databases
with the [Spectra](https://github.com/rformassspectrometry/Spectra) R package.

For more information see the package
[homepage](https://rformassspectrometry.github.io/MsBackendMassbank).


# Installation

The package can be installed with

```r
install.packages("BiocManager")
BiocManager::install("MsBackendMassbank")
```


# Contributions

Contributions are highly welcome and should follow the [contribution
guidelines](https://rformassspectrometry.github.io/RforMassSpectrometry/articles/RforMassSpectrometry.html#contributions).
Also, please check the coding style guidelines in the [RforMassSpectrometry
vignette](https://rformassspectrometry.github.io/RforMassSpectrometry/articles/RforMassSpectrometry.html).

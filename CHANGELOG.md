# Changelog

Release changelog are available here : https://gitlab.inria.fr/croco-ocean/croco_pytools/-/releases

## [1.0.1] - 2024-05-30

### Added
- In prepro/ :
  - In tools_fort_routines : solve issue with large netCDF input file processing bry files using using netcdf4
  - Add `s_rho` variables writing in `make_ini.py` and `make_bry.py`

### Fixed
- **Major In prepro/** : 
  - wrong vertical integration to compute barotropic transport when using `conserv=1` in `prepro/make_bry.py` . It induce a wrong velocity fields and strong discontinuity from month to month, especially problematic at depth when using the MONTHLY putput format.
  - See issue [#3](https://gitlab.inria.fr/croco-ocean/croco_pytools/-/issues/3) ;
check merge request and documentation if interested

### Changed

### Deprecated

### Removed

### Other


## [1.0.0] - 2024-04-22 : First croco_pytools release
### Added

### Fixed

### Changed

### Deprecated

### Removed

### Other

# Changelog

Release changelog are available here : https://gitlab.inria.fr/croco-ocean/croco_pytools/-/releases

## [x.x.x] - xxxx-xx-xx

### Added

### Fixed

- Nesting with AGRIF: correct issue #20: child grid corners indices are now correctly written in the croco_grd file and AGRIF_FixedGrids.in

- prepro/make_tides: some numpy meshgrid was lost at a previous commit, correct to have make_tides work again ( solve #31 )


### Changed

- prepro/make_bry: add Morig Dorig in user changes instead of considering 01 01 (solve issue #28)

### Deprecated

### Removed

### Other

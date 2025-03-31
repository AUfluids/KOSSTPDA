# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.0.0] - 2025-03-28

### Added
- Progressive Data Augmentation (PDA) framework for turbulence modeling
- Separation factor correction for improved flow separation prediction
- Secondary flow prediction capabilities
- Support for both incompressible and compressible flows
- Configurable model coefficients through dictionary
- Decay control functionality for far-field conditions
- Comprehensive test cases for validation

### Changed
- Updated model structure to use base class architecture
- Improved code organization and documentation
- Enhanced error handling and boundary conditions
- Optimized performance for 2D and 3D flows
- Updated coefficient initialization and validation

### Fixed
- Improved stability in high adverse pressure gradient regions
- Enhanced secondary flow prediction accuracy
- Fixed boundary condition handling
- Corrected tensor operations for better numerical stability

### Removed
- Legacy model implementations
- Deprecated coefficient sets
- Outdated documentation

## [1.0.0] - 2023-12-15

### Added
- Initial release of kOmegaSSTPDA model
- Basic separation and secondary flow corrections
- Support for incompressible flows
- Standard test cases

### Changed
- Initial implementation of model structure
- Basic documentation

### Fixed
- Initial bug fixes and stability improvements

### Removed
- None (initial release) 
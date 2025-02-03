# KOSSTPDA
Progressive data-augmented k-omega SST model as proposed by Amarloo and Rincón (2023) for OpenFOAM.
Developed by Fluid Physics & Turbulence research group at Aarhus University.

## Table of Contents
- [Overview](#overview)
- [Features](#features)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Model Configuration](#model-configuration)
- [Test Results](#test-results)
- [License](#license)
- [Citations](#how-to-cite)
- [Authors](#authors)
- [References](#references)
- [Disclaimer](#disclaimer)

## Overview
Progressive augmentation of the kOmegaSST turbulence model in OpenFOAM. Standard optimised values are given by default in the model.
More information is available in the publications listed at the end of this file.
This correction enhances the performance of kOmegaSST turbulence model in capturing:

### Features
- **Separation Flow**: Optimised for separation after bumps or cases with high adverse pressure gradient, including:
  - Periodic hills
  - Curved backward-facing step
  - Converging-diverging channel
  - Parametric bumps

- **Secondary Flow**: Specifically the Prandtl's second kind of secondary flow (corner flow) induced by heterogenities of the Reynolds stresses.

### Implementation Details
- Optimised for 2D cases
- Tested and validated on 3D flows
- Configurable through five user-modifiable coefficients
- Default optimised values provided

## Requirements
- Compatible with OpenFOAM-v2412 and previous ESI versions

## Installation
1. Clone the repository:
     ```
     git clone https://github.com/AUfluids/KOSSTPDA.git
     ```

2. Make the installation script executable:
     ```
     cd KOSSTPDA
     chmod a+x Allwmake
     ```

3. Compile the model:
     ```
     ./Allwmake
     ```

## Usage
1. Add required library to `controlDict`:
   - For incompressible flow:
     ```
     libs ( "libPDAIncompressibleTurbulenceModels" );
     ```
   - For compressible flow:
     ```
     libs ( "libPDAcompressibleTurbulenceModels" );
     ```

2. Specify in `turbulentProperties`:
   ```
   RASModel kOmegaSSTPDA;
   ```

3. Add to `system/fvSchemes`:
   ```
   divSchemes
   {
       div(dev(((2*k)*bijDelta)))    Gauss linear;
   }
   ```

> **Note**: It is recommended to first run your case with standard kOmegaSST before switching to kOmegaSSTPDA.

## Model Configuration
### Mode Selection
By default, the model activates both secondary and separation effects. If desired, one can change the models as follows:
   ```
   separationMode 4; // default: 4 - off: 0 | ModelI: 1 | ModelII: 2 | ModelIII: 3 | ModelIV: 4
   secondaryMode 2; // default: 2 - off: 0 | ModelI: 1 | ModelII: 2
   ```

If you use 0, the extra effects are deactivated, and the standard kOmegaSST is used.
For info about the differences within these models, users are referred to the publications corresponding to the development of the each model (can be found at the end of the document)
   
### Optional Stability Settings
In case of stability and convergence issues, we also suggest the following setting for the new model to be tested.
Otherwise, these coefficients are automatically assigned with values corresponding to the models.

   ```
   // Separation Flow coefficients
   separationMode 1;
   separationLambda1 1;
   separationLambda2 1;
   C0 -1;
   C1 0;
   C2 0;
   // Secondary Flow coefficients
   secondaryMode 1;
   A0 -1;
   A1 0;
   A2 0;
   ```

## Test Results
### Separation Effect
Results for curved backward-facing step (Reb = 13700, Models I and III):
![Separation Effect Contours](https://github.com/AUfluids/KOSSTPDA/blob/main/testCases/CBFS_Reb13700/contours_comparisonCBFS.png)
![Separation Effect Comparison](https://github.com/AUfluids/KOSSTPDA/blob/main/testCases/CBFS_Reb13700/quantitative_comparison_CBFS.png)

### Secondary Effect
Results for duct flow (AR = 1, Reb = 3500, Model II):
![Secondary Effect U](https://github.com/AUfluids/KOSSTPDA/blob/main/testCases/ductFlowAR1Reb3500/SD_u.png)
![Secondary Effect Profiles](https://github.com/AUfluids/KOSSTPDA/blob/main/testCases/ductFlowAR1Reb3500/SD_profiles.png)

## Target platform
The code is known to work with OpenFOAM-v2412 and previous ESI versions.

## Authors
Mario Javier Rincón <mjrp@mpe.au.dk>
Ali Amarloo <amarloo@mpe.au.dk>

## References
For more details about the separation and secondary flow effects, refer to the publications at: 
 - [Progressive augmentation of turbulence models for **flow separation** by multi-case computational fluid dynamics driven surrogate optimization](https://doi.org/10.1063/5.0174470)
 - [Progressive augmentation of Reynolds stress tensor models for **secondary flow** prediction by computational fluid dynamics driven surrogate optimisation](https://doi.org/10.1016/j.ijheatfluidflow.2023.109242)


## How to cite
Please, cite this library using the following two publications: 

Amarloo and Rincón (2023):

    @article{amarloo2023progressive,
      title={Progressive augmentation of turbulence models for flow separation by multi-case computational fluid dynamics driven surrogate optimization},
      author={Amarloo, Ali and Rinc{\'o}n, Mario Javier and Reclari, Martino and Abkar, Mahdi},
      journal={Physics of Fluids},
      volume={35},
      number={12},
      year={2023},
      publisher={AIP Publishing}
    }

Rincón and Amarloo (2023)

         @article{rincon2023progressive,
         title = {Progressive augmentation of Reynolds stress tensor models for secondary flow prediction by computational fluid dynamics driven surrogate optimisation},
         journal = {International Journal of Heat and Fluid Flow},
         volume = {104},
         pages = {109242},
         year = {2023},
         issn = {0142-727X},
         doi = {https://doi.org/10.1016/j.ijheatfluidflow.2023.109242},
         author = {Mario Javier Rincón and Ali Amarloo and Martino Reclari and Xiang I.A. Yang and Mahdi Abkar},
         }


## Disclaimer
This offering is not approved or endorsed by OpenCFD Limited, the producer of the OpenFOAM software and owner of the OPENFOAM® and OpenCFD® trade marks.

Detailed information on the OpenFOAM trademark can be found at

http://www.openfoam.com/legal/trademark-policy.php
http://www.openfoam.com/legal/trademark-guidelines.php
For further information on OpenCFD and OpenFOAM, please refer to

http://www.openfoam.com

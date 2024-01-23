# KOSSTPDA
Progressive data-augmented k-omega SST model
as proposed by Amarloo and Rincón (2023) for OpenFOAM.
Developed by Fluid Physics & Turbulence research group at Aarhus University.


## License
    This program is free software: you can redistribute and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

## Description
Progressive augmentation of the kOmegaSST turbulence model in OpenFOAM.
This correction enhances the performance of kOmegaSST turbulence model in 
capturing the separation flow and secondary flow:
- **Separation** Flow: specifically optimized for separation after bumps or for cases with high adverse pressure gradient
  including periodic hills, curved backward-facing step, converging-diverging channel, and parametric bumps.
- **Secondary** Flow: specifically the Prandtl's second kind of secondary flow (aka corner flow) induced by Reynolds stress tensor.


 
The implementation of the augmented model has been optimized for 2D cases but 
its applicability has been tested on 3D flows, yielding similar results.
Five coefficients can be modified by the user to change the model's performance. 
Standard optimised values are given by default in the model.
More information is available in the publication listed at the of this file.

## Target platform
The code is known to work with OpenFOAM v2312 and previous ESI versions.

## Authors
Ali Amarloo <amarloo@mpe.au.dk>
Mario Javier Rincón <mjrp@mpe.au.dk>

## Instructions

1. Download the source code using git:

         git clone https://github.com/AUfluids/KOSSTPDA.git

2. Enter the directory where the source code has been extracted, and compile it by typing: 

         wmake

3. Add the following line to the _controlDict_ of your case:

         libs ( "libPDAIncompressibleTurbulenceModels" ) ;

4. Specify the following in _turbulentProperties_.

         RASModel kOmegaSSTPDA;
   
(We suggest first running your case with the standard kOmegaSST and then changing the model to kOmegaSSTPDA, to avoid possible errors)

NOTE: You might have to define the bijDelta term in system/fvSchemes file, here you have an example:

         divSchemes
         {
                  div(dev(((2*k)*bijDelta)))          Gauss linear;
         }

5. (Optional) By default, the model activates both secondary and separation effects. If desired, one can change the models as follows: 

         separationMode      4; //optional - default:4 - off:0 | ModelI:1 | ModelII:2 | ModelIII:3 | ModelIV:4
         secondaryMode       2; //optional - default:2 - off:0 | ModelI:1 | ModelII:2
   
If you use 0, the extra effects are deactivated, and the standard kOmegaSST is used.
For info about the differences within these models, users are referred to the publications corresponding to the development of the each model (can be found at the end of the document)
   
6. (Optional) In case of stability and convergence issues, we also suggest the following setting for the new model to be tested.
   Otherwise, these coefficients are automatically assigned with values corresponding to the models. 

           //Separation Flow coefficients
            separationLambda1   1;              //optional - default taken from separationMode 4
            separationLambda2   1;              //optional - default taken from separationMode 4
            C0                  -1;             //optional - default taken from separationMode 4
            C1                  0;              //optional - default taken from separationMode 4
            C2                  0;              //optional - default taken from separationMode 4
            //Secodnary Flow coefficients
            A0                  -1;             //optional - default taken from secondaryMode 2
            A1                  0;              //optional - default taken from secondaryMode 2
            A2                  0;              //optional - default taken from secondaryMode 2



You can also check the test cases of CBFS_Reb13700 and ductFlowAR1Reb3500 in the folder testCases.

## Test results of the separation effect

Results of using KOSSTPDA for curved backward-facing step with bulk Reynolds number of 13700 **I** and **III**:
![alt text](https://github.com/AUfluids/KOSSTPDA/blob/main/testCases/CBFS_Reb13700/contours_comparisonCBFS.png)
![alt text](https://github.com/AUfluids/KOSSTPDA/blob/main/testCases/CBFS_Reb13700/quantitative_comparison_CBFS.png)

For more details about the separation effect, refer to the publication at: 
[Progressive augmentation of Reynolds stress tensor models for secondary flow prediction by computational fluid dynamics driven surrogate optimisation](https://doi.org/10.1016/j.ijheatfluidflow.2023.109242)
or [arXiv](https://arxiv.org/abs/2308.12720)


## Test results of the secondary effect

Results for a duct flow of aspect ratio 1 at bulk Reynolds number 3500 for Model **II**:
![alt text](https://github.com/AUfluids/KOSSTPDA/blob/main/testCases/ductFlowAR1Reb3500/SD_u.png)
![alt text](https://github.com/AUfluids/KOSSTPDA/blob/main/testCases/ductFlowAR1Reb3500/SD_profiles.png)

For more details about the secondary effect, refer to the publication at: 
[Progressive augmentation of Reynolds stress tensor models for secondary flow prediction by computational fluid dynamics driven surrogate optimisation](https://doi.org/10.1016/j.ijheatfluidflow.2023.109242)
or [arXiv](https://arxiv.org/abs/2308.12720)

## How to cite
Please, cite this library using the following publications: 

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

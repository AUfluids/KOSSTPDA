#!/bin/bash
#
#-------------------------------------------------#
# Contributor: Mario Javier Rincon                #
# Updated on:  25 March 2024                      #
#-------------------------------------------------#
# Topic:       KOSSTPDA                           #
# OpenFOAM:    v2212                              #
#-------------------------------------------------#
# Institution: Aarhus University / Kamstrup       #
# Email:       mjrp@mpe.au.dk                     #
#-------------------------------------------------#

# Source OpenFOAM
# source /lib/openfoam/openfoam2212/etc/bashrc

# processors=4

foamListTimes -rm
rm -r postProcessing
rm -r logs
mkdir logs

# decomposePar -force
#------------------------------------------------------------------------------
# SIMULATION
#------------------------------------------------------------------------------
# echo -e "   - Running simpleFoam"
simpleFoam > logs/simulation
#------------------------------------------------------------------------------
# echo -e "   - Running post processing"
# reconstructPar -latestTime

# rm -r processor*
python postProcess.py
foamLog logs/simulation
python plotResiduals.py

# rhoSimpleFoam -postProcess -latestTime -funcs '(wallShearStress writeCellCentres airfoilValues airfoilForces)'

#------------------------------------------------------------------------------
# echo -e "   - Simulation done :)"
#------------------------------------------------------------------------------
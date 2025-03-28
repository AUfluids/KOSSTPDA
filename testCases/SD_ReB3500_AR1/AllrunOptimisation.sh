#!/bin/bash
#
#-------------------------------------------------#
# Contributor: Mario Javier Rincon                #
# Updated on:  10 December 2020                   #
#-------------------------------------------------#
# Topic:       Optimisation Flow Meters           #
# OpenFOAM:    v2006                              #
#-------------------------------------------------#
# Institution: Aarhus University / Kamstrup       #
# Email:       mjrp@eng.au.dk                     #
#-------------------------------------------------#
#
source /lib/openfoam/openfoam2212/etc/bashrc
#------------------------------------------------------------------------------
# Source OpenFOAM run functions
# .$WM_PROJECT_DIR/bin/tools/RunFunctions

#------------------------------------------------------------------------------
#cd ${0%/*} || exit 1 # Run from this directory, otherwise exit
# clear

touch STARTFILE

rm -r logs

foamListTimes -rm

mkdir logs
#------------------------------------------------------------------------------
echo -e "   - Running simpleFoam"
simpleFoam > logs/simulation

#------------------------------------------------------------------------------
echo -e "   - Running post processing"
postProcess -func vorticity
simpleFoam -postProcess -func R
# python saveObjFunc_SD.py
# python plot_SD.py
python postProcessMulti.py

foamLog logs/simulation

python plotResiduals.py

#------------------------------------------------------------------------------
echo -e "   - Simulation done :)"

touch DONEFILE
#------------------------------------------------------------------------------

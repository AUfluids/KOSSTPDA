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

for case in CF_590 PH_2800 SD_ReB3500_AR1  
do
    cd "${case}"

    pwd > log

	touch STARTFILE

	#------------------------------------------------------------------------------
	echo -e "   - Running simulation on $case"
	source AllrunOptimisation.sh&

	touch DONESIMULATION

	pwd
	
	cd ..

done

wait

touch DONESIMULATIONS

# for case in PH CBFS CF5200
# do
#     cd "${case}"

# 	#------------------------------------------------------------------------------
# 	echo -e "   - Running post processing on $case"
# 	# simpleFoam -postProcess -func wallShearStress -latestTime

# 	python postProcess.py

# 	foamLog logs/simulation >> dictionaries

# 	python plotResiduals.py

# 	# #------------------------------------------------------------------------------
# 	echo -e "   - Simulation $case done :)"
# 	touch DONEFILE

# 	cd ..
# done

#------------------------------------------------------------------------------
echo -e "   - Simulations done :)"

touch DONEFILE

# #------------------------------------------------------------------------------
# wait; echo -e "   - Running post processing on $case" &
# wait; simpleFoam -postProcess -func wallShearStress -latestTime &

# wait; python postProcess.py &

# wait; foamLog logs/simulation &

# wait; python plotResiduals.py &

# #------------------------------------------------------------------------------
# wait; echo -e "   - Simulation $case done :)" &
# wait; touch DONEFILE &

# cd ..

# done



#------------------------------------------------------------------------------

This README file is in development



:::::::::::::::::::::::::::::::::
1) Get the code and needed inputs
:::::::::::::::::::::::::::::::::

git clone https://github.com/knuenz/NRQCDglobalfit
cd NRQCDglobalfit
mkdir Figures

I use ROOT 5.34/17, but all newer versions should also work.

Compile with 'make' to test if it compiles

The script runNRQCDglobalfit.sh is the global steering macro, that contains all settings and steers all executables.
In this script, change the 'storagedir' to a directory where you can store several GBs.

There are some needed inputs that are not on git. They are stored on /afs/hephy.at/user/k/knuenz/public/forNRQCD/InputsNotOnGit.
Copy the subdirectories into your storagedir.

If you want to update the code with the newest version in git, go to NRQCDglobalfit and do
git pull



::::::::::::::::::::::::::::::::::::::::::
2) Some general information about the code
::::::::::::::::::::::::::::::::::::::::::

a) Class NRQCDglobalfitObject (see NRQCDglobalfitObject.cc/.h files in src/ and interface/):

This class defines an object that contains all information about one measurement point:
-Central value, all uncertainties (stat, syst, global, polarization 'uncertainty'...)
-pT/y of the bin, experiment, state, MeasurementID (cross section, lamth, ...)
It also contains the information about the theoretical model, calculated for this specific result and kinematic bin:
-SDCs for all channels, as well as the SDCs for all states that feed down to this state
-Polaritation for all channels, as well as the SDCs for all states that feed down to this state
-Uncertainties for all model ingredients
The NRQCDglobalfitObjects are saved as text files, and can be read back into the framework from these text files


b) Global definitions:

The global variables are defined in interface/GlobalVar.h
Important variables:
-If you want to include or not include the 3PJ, this can be changed by changing the variables nColorChannels_S, nColorChannels_P, ColorChannelNameTexS, ColorChannelNameTexP
-Including or not including the singlet can happen by changing the ColorSingletME
-Including the theoretical uncertainties by changing nModelSystematicScales (0... no theoretical uncertainties, 3... include theoretical uncertainties (3: for the three color channels CS, 1S0, 3S1))

::::::::::::::::::::::::::::::::::::::::
3) Convert data inputs too needed format
::::::::::::::::::::::::::::::::::::::::

runNRQCDglobalfit.sh:
-set run_ConvertDataInput=1
-define DataID (=output)

Input files: HEPDATA/...
Output files: storagedir/DataID/${DataID}/ConvertedData....txt

There is one output text file for each measurement, containing one NRQCDglobalfitObject, that contains the data-part of the object, not yet the model.



:::::::::::::::::::::::::::::::::::::::::
4) Convert model inputs too needed format
:::::::::::::::::::::::::::::::::::::::::

a) ScaleBKmodel


runNRQCDglobalfit.sh:
-set run_ScaleBKmodel=1
-define OriginalModelID (=input, output)
-define OriginalModelID

Input files: storagedir/OriginalModelID/${OriginalModelID}/...txt
Output files: storagedir/OriginalModelID/${OriginalModelID}/TGraphs_scaled.root

The purpose of this step is to take the original theory inputs, include all the assumptions of our study (e.g. extrapolation towards high/low pT and throughout all rapidity regions needed, scaling pT/M), and convert the new theory ingredients as TGraphs.
The input files are the BK model ingredients (direcly the values sent by the authors)
The output file contains TGraphs for all scaled SDCs and polarizations


b) ConvertBKmodelToTTree

runNRQCDglobalfit.sh:
-set run_ConvertBKmodelToTTree=1
-define OriginalModelID (=input)
-define ModelID (=output)

src/ConvertBKmodelToTTree.cc:
n_nTuple=5e7 -> for debugging/testing you can set it to 5e5
ignoreState: This array of bools defines for which states this is run (speeds up the running if you are only interested in the models of a subset of states. Default in git: run only psi2S, ups3S)

Input files: storagedir/OriginalModelID/${OriginalModelID}/TGraphs_scaled.root
Output file: storagedir/ModelID/${ModelID}/OriginalModelTTree.root

The purpose of this step is to generate, for each state and each color channel, a number of n_nTuple events that are distributed in the kinematic region where the scaled model ingredients are defined.
The events are weighed according to the values of the SDCs at a given pT/rap of the event. A costh and phi (HX) are generated as well, according to the lamth, lamph, lamtp of the input model.


c) ConvertModelInput

runNRQCDglobalfit.sh:
-set run_ConvertModelInput=1
-define ModelID (=input)
-define ModelID (=output)

Input files: storagedir/ModelID/${ModelID}/OriginalModelTTree.root
Output files: storagedir/ModelID/${ModelID}/ModelIngredients.root and ModelIngredients_consts_star.txt

The purpose of this step: Starting from the OriginalModelTTree.root events, looping over all color channels and possible feed-down decay cascades (including e.g. 3S->2P->2S->1P->1S decays), saving a TTree with all the decayed kinematic properties (changed pT/rap) and the transferred polarization of the decay (as new costh and phi values).
In the default git case (run only on psi2S and ups3S) the output file is identical to the input file, as there is no feed down into the psi2S and ups3S (chib3P->ups3S decays are switched off, but can be switched on in interface/GlobalVars.h, changing the branchings)
Also, a text file is created that contains the estimated values of the SDCs of each color channel, at a given pT/M* and rap* defined in interface/GlobalVars.h to be 6 and 0, respectively. These 'consts_star' values are needed to define the fractions (of the individual color octed components wrt the inclusive color octet contribution, at pT/M* and rap*) used in the fit.



::::::::::::::::::::::::::::::::
5) Combine data points and model
::::::::::::::::::::::::::::::::

runNRQCDglobalfit.sh:
-set run_CombineDataModel=1
-define DataID (=input)
-define ModelID (=input)
-define ModelSystScaleID (=input)
-define DataModelCombID (=output)

Input files: storagedir/DataID/${DataID}/...txt, storagedir/ModelID/${ModelID}/ModelIngredients.root and storagedir/ModelSystScaleID/${ModelSystScaleID}/....root
Output files: storagedir/DataModelCombID/${DataModelCombID}/ConvertedDataModel...txt

This executable loops over all data points, and calculates
-The SDCs for each color channel and feed down state, for the given kinematic cell
-The same for the polarization parameters
-The theoretical uncertainties for each of the above mentioned values

It then saves the final NRQCDglobalfitObjects into text files, containing both the data and model information.

Remark: The inputs for the theoretical uncertainties are stored in storagedir/ModelSystScaleID. These TGraphs are calculated 'externally' to the framework in src/CalcTheoryUncertaintyGraphs.cc, steered by runTheoryUncertainty.sh. These files are not yet ready to be used within the framework, and some changes are needed. But the current TGraph files are copied in 1)



:::::::::::::::::
6) The global fit
:::::::::::::::::

runNRQCDglobalfit.sh:
-set run_SamplePPD=1
-define ModelID (=input)
-define DataModelCombID (=input)
-define JobID (=output)

Input files: storagedir/DataModelCombID/${DataModelCombID}/ConvertedDataModel...txt, storagedir/ModelID/${ModelID}/ModelIngredients_consts_star.txt
Output files: storagedir/JobID/${JobID}/results.root and two text files containing information about the chi2, ndf, and which are the free parameters of the fit (for the plotting macros)

In short: The code runs a Markov chain. For each set of the fit parameters (the fractions), the code loops over all data points, calculates the corresponding model (e.g. model cross section for this data point), compares it to the data central value and uncertainties, and decides if this set of fit parameters values are kept. If so, the output TTree is filled (very similar to the polarization analysis approach)

There are several settings that can be set:

Defining which data to be used in the fit:
pTMin
rapMin
rapMax
useSstatesOnly 
usePstatesOnly
useCharmoniumOnly
useBottomoniumOnly
useOnlyState
useCrossSectionOnly
usePolarizationOnly

Which minimizer to use
Minimizer		#0...MH, 1...Minuit

Metropolis Hastings settings:
nBurnIn
nSample

Define if the nuisance parameters are floating in the fit:
SampleNp
SampleNp_consts_star



:::::::::::::::::::::::::::::::::::::::::
7) Interpretation of the fit and plotting
:::::::::::::::::::::::::::::::::::::::::

a) InterpretPPD

runNRQCDglobalfit.sh:
-set run_InterpretPPD=1
-define JobID (=input)
-define JobID (=output)

Input files: storagedir/JobID/${JobID}/*
Output files: storagedir/JobID/${JobID}/results*.txt

This code runs over the output TTree of the fit, calculating the MPV and uncertainties on the fit parameters, making some plots of the PPD projections


b) PlotCompareDataModel

runNRQCDglobalfit.sh:
-set run_PlotCompareDataModel=1
-define JobID (=input)
-define JobID (=output)

Input files: storagedir/JobID/${JobID}/*
Output files: storagedir/JobID/${JobID}/Figures/* (copied to Figures/${JobID}/)

This code makes the results plots and prediction plots.
In case of the results plots, set
PredictionPlot=false
PlotVSpToverM=false
PredictionDataPlot=false

In case of the prediction plots, set
PredictionPlot=true
PlotVSpToverM=true
PredictionDataPlot=false
For the prediction plots, you need to generate ToyData up to high pT values, sucht that the plotting can be done... (to be explained in person)

The plotting code is very complex, and can not reasonably be explained in this README. One tip: the plotting takes a long time (if in the fit there were many samplings used). For debugging, to speed up, you have tu ensure in line 681 that LoopThroughPPD=false. This speeds up the plotting.





::::::::::::::::::
X) Additional Code
::::::::::::::::::

a) ToyMC

There is the possibility to generate pseudo data corresponding to a certain model, for certain kinematic bins, experiments, measurementIDs, with certain uncertainties, ...
The settings are in interface/ToyData.h, the source code is in src/GenerateToyData.cc

This code generates text files containing the data information of the NRQCDglobalfitObjects. These have to be complemented by the model information with the CombineDataModel step.


b) run_FitPtDists
Macro used to do the initial pT/M scaling study


c) run_PlotPPD
Plotting certain projections of the PPD, e.g. comparing ellipses of results with/without polarization constraint


d) run_PlotPPDderivative
Plotting certain projections of the PPD, e.g. ratios of matrix elements


e) run_PlotPPD_vs_pTmin
Plotting the resutls of the pT scan vs pTmin






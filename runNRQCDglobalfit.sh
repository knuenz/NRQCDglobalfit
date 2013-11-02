#!/bin/sh

storagedir=/scratch/knuenz/NRQCD/NRQCDglobalfit

#store here the NRQCDglobalfitObjects containing data measurements
DataID=October9_DataForPtFits #September27_Psi2SUps3Sonly 

### This is the location of the original model (nTuple provided by Sergey or BK model txt file)
#OriginalModelID=GWWZmodel
OriginalModelID=BKmodel70
### store here the ModelIngredients.root file and consts_star file
ModelID=October31_UNCORR_BKmodel70_OriginalAllStates_pTstar_over_m6_UpsAssumption #HP original
#ModelID=October27_BKmodel70_ScaledAllStates_pTstar_over_m6_UpsAssumption #HP scaled

### store here the NRQCD objects combining data and model predictions
ModelSystScaleID=NLOmLO

DataModelCombID=August24_Psi2Sonly_BKmodel_NoAbsDef_NoLamphLamtp
#DataModelCombID=October27_BKmodel70_ScaledAllStates_pTstar_over_m6_UpsAssumption_AddTheoryUncertainty_AddNonModelDataPointsForPlots

### store here the output TTree of the likelihood sampling, and all Figures of the results:
for JobID in November1_HPsequence_3_cs_pol_pTmin9_psi2S_wo3PJ;do
#for JobID in BKtruth;do

run_ConvertDataInput=0

#Generate model Trees
run_ConvertNTupleToTTree=0
run_ScaleBKmodel=0
run_ScaleGWWZmodel=0
run_ConvertBKmodelToTTree=0
run_ConvertModelInput=0
run_CombineDataModel=0
run_GenerateToyData=0

#Fit and plot
run_SamplePPD=0
run_InterpretPPD=0
run_PlotCompareDataModel=1

#individual macros:::
run_FitPtDists=0
run_PlotPPD=0

##################################
########## SETTINGS ##############
##################################

### ConvertNTupleToTTree

### ConvertDataInput

### ConvertModelInput
useToyModel=false
calcOnlyConstsStat=true

### CombineDataModel

### GenerateToyData

### SamplePPD
Minimizer=0		#0...MH, 1...Minuit
nBurnIn=5000
nSample=50000
SampleNp=true
SampleNp_consts_star=true

### InterpretPPD
nSigma=1
MPValgo=4 		#1...mean,2...gauss,3...gauss-loop with chi2<2, 4...(new)...mode

### PlotCompareDataModel

##################################
########## SELECTION #############
##################################
# To be added in SamplePPD and the plotting

pTMin=9
pTMax=100
rapMin=-10
rapMax=10
useSstatesOnly=false
usePstatesOnly=false
useCharmoniumOnly=false
useBottomoniumOnly=false
useOnlyState=3 #switch off by setting it to an int > N_STATES

##################################
############# CODE ###############
##################################

touch src/*
make

cp ConvertDataInput ConvertDataInput_${DataID}
cp FitPtDists FitPtDists_${DataID}
cp ConvertNTupleToTTree ConvertNTupleToTTree_${ModelID}
cp ScaleBKmodel ScaleBKmodel_${ModelID}
cp ScaleGWWZmodel ScaleGWWZmodel_${ModelID}
cp ConvertBKmodelToTTree ConvertBKmodelToTTree_${ModelID}
cp ConvertModelInput ConvertModelInput_${ModelID}
cp CombineDataModel CombineDataModel_${DataModelCombID}
cp GenerateToyData GenerateToyData_${DataModelCombID}
cp SamplePPD SamplePPD_${JobID}
cp InterpretPPD InterpretPPD_${JobID}
cp PlotCompareDataModel PlotCompareDataModel_${JobID}
cp PlotPPD PlotPPD_${JobID}

if [ ${run_ConvertDataInput} -eq 1 ]
then
./ConvertDataInput_${DataID} ${DataID}=DataID ${storagedir}=storagedir
fi
if [ ${run_FitPtDists} -eq 1 ]
then
./FitPtDists_${DataID} ${DataID}=DataID ${storagedir}=storagedir
mkdir Figures/${DataID}
cp -r ${storagedir}/DataID/${DataID}/Figures/* Figures/${DataID}/
fi
if [ ${run_ConvertNTupleToTTree} -eq 1 ]
then
./ConvertNTupleToTTree_${ModelID} ${OriginalModelID}=OriginalModelID ${ModelID}=ModelID ${storagedir}=storagedir
fi
if [ ${run_ScaleBKmodel} -eq 1 ]
then
./ScaleBKmodel_${ModelID} ${OriginalModelID}=OriginalModelID ${storagedir}=storagedir
fi
if [ ${run_ScaleGWWZmodel} -eq 1 ]
then
./ScaleGWWZmodel_${ModelID} ${OriginalModelID}=OriginalModelID ${storagedir}=storagedir
fi
if [ ${run_ConvertBKmodelToTTree} -eq 1 ]
then
./ConvertBKmodelToTTree_${ModelID} ${OriginalModelID}=OriginalModelID ${ModelID}=ModelID ${storagedir}=storagedir
fi
if [ ${run_ConvertModelInput} -eq 1 ]
then
./ConvertModelInput_${ModelID} ${ModelID}=ModelID useToyModel=${useToyModel} ${storagedir}=storagedir calcOnlyConstsStat=${calcOnlyConstsStat}
fi
if [ ${run_CombineDataModel} -eq 1 ]
then
./CombineDataModel_${DataModelCombID} ${ModelSystScaleID}=ModelSystScaleID ${ModelID}=ModelID ${DataID}=DataID ${DataModelCombID}=DataModelCombID ${storagedir}=storagedir
fi
if [ ${run_GenerateToyData} -eq 1 ]
then
./GenerateToyData_${DataModelCombID} ${ModelID}=ModelID ${DataID}=DataID ${DataModelCombID}=DataModelCombID ${storagedir}=storagedir
cp interface/ToyData.h ${storagedir}/DataModelCombID/${DataModelCombID}/ToyData.h
fi
if [ ${run_SamplePPD} -eq 1 ]
then
./SamplePPD_${JobID} ${nBurnIn}nBurnIn ${nSample}nSample ${ModelID}=ModelID  ${JobID}=JobID ${DataModelCombID}=DataModelCombID ${pTMin}pTMin ${pTMax}pTMax ${rapMin}rapMin ${rapMax}rapMax ${useOnlyState}useOnlyState ${Minimizer}Minimizer useSstatesOnly=${useSstatesOnly} usePstatesOnly=${usePstatesOnly} useCharmoniumOnly=${useCharmoniumOnly} useBottomoniumOnly=${useBottomoniumOnly} ${storagedir}=storagedir SampleNp=${SampleNp} SampleNp_consts_star=${SampleNp_consts_star} ${MPValgo}MPValgo ${nSigma}nSigma
fi
if [ ${run_InterpretPPD} -eq 1 ]
then
./InterpretPPD_${JobID} ${JobID}=JobID ${MPValgo}MPValgo ${nSigma}nSigma ${storagedir}=storagedir
mkdir Figures/${JobID}
cp -r ${storagedir}/JobID/${JobID}/Figures/* Figures/${JobID}/
cp -r ${storagedir}/JobID/${JobID}/*.txt Figures/${JobID}/
fi
if [ ${run_PlotCompareDataModel} -eq 1 ]
then
./PlotCompareDataModel_${JobID} ${JobID}=JobID ${ModelID}=ModelID ${nSigma}nSigma ${MPValgo}MPValgo ${storagedir}=storagedir ${DataModelCombID}=DataModelCombID ${pTMin}pTMin ${pTMax}pTMax ${rapMin}rapMin ${rapMax}rapMax ${useOnlyState}useOnlyState useSstatesOnly=${useSstatesOnly} usePstatesOnly=${usePstatesOnly} useCharmoniumOnly=${useCharmoniumOnly} useBottomoniumOnly=${useBottomoniumOnly} SampleNp=${SampleNp} SampleNp_consts_star=${SampleNp_consts_star}
mkdir Figures/${JobID}
cp -r ${storagedir}/JobID/${JobID}/Figures/* Figures/${JobID}/
cp -r ${storagedir}/JobID/${JobID}/*.txt Figures/${JobID}/
fi
if [ ${run_PlotPPD} -eq 1 ]
then
./PlotPPD_${JobID} ${JobID}=JobID ${storagedir}=storagedir
mkdir Figures/PlotPPD/
mkdir Figures/PlotPPD/${JobID}
fi



rm ConvertDataInput_${DataID}
rm ConvertNTupleToTTree_${ModelID}
rm ConvertBKmodelToTTree_${ModelID}
rm ConvertModelInput_${ModelID}
rm CombineDataModel_${DataModelCombID}
rm GenerateToyData_${DataModelCombID}
rm SamplePPD_${JobID}
rm InterpretPPD_${JobID}
rm PlotCompareDataModel_${JobID}
rm PlotPPD_${JobID}

done

#!/bin/sh

storagedir=/scratch/knuenz/NRQCD/NRQCDglobalfit

#store here the NRQCDglobalfitObjects containing data measurements
#DataID=October9_DataForPtFits #September27_Psi2SUps3Sonly 
DataID=December7_addJpsiData_addSomPolUncertainties 
#DataID=February7_addNewPrelCMSPsi2SUps3SmergedRapData


### This is the location of the original model (nTuple provided by Sergey or BK model txt file)
#OriginalModelID=GWWZmodel
OriginalModelID=BKmodel70
### store here the ModelIngredients.root file and consts_star file
########################ModelID=October31_UNCORR_BKmodel70_OriginalAllStates_pTstar_over_m6_UpsAssumption #HP original
########################ModelID=October27_BKmodel70_ScaledAllStates_pTstar_over_m6_UpsAssumption #HP scaled
ModelID=December2_BKmodel70_EpScaled_newRapNorm_pTstar_over_m6
#ModelID=March1_BKmodel70_EpScaled_newRapNorm_pTstar_over_m6_inclPJ_ModelOnly2States

### store here the NRQCD objects combining data and model predictions
ModelSystScaleID=NLOmLO_Nov30

##########################DataModelCombID=August24_Psi2Sonly_BKmodel_NoAbsDef_NoLamphLamtp
##########################DataModelCombID=October27_BKmodel70_ScaledAllStates_pTstar_over_m6_UpsAssumption_AddTheoryUncertainty_AddNonModelDataPointsForPlots_NEW_THuncert
##########################DataModelCombID=February7_allNew_newRapNorm_newTHunc_realChiPol_addJpsiData_addNewPrelCMSPsi2SUps3SmergedRapData
DataModelCombID=December7_allNew_newRapNorm_newTHunc_realChiPol_addJpsiData_addSomePolUncertainties
#DataModelCombID=March1_allNew_newRapNorm_noTHunc_realChiPol_addJpsiData_addSomePolUncertainties_inclPJ_ModelOnly2States
#DataModelCombID=March3_ToyData_THunc

##### November30: new TH uncertainty, new P-wave CCs, Ep BK model, more data incl chi, new iExpYear, allScaled

for pTMin in 0;do
#for pTMin in 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 22 25 27 30;do #new psi2S scan
#for pTMin in 4 5 6 7 8 9 10 11 12 13 14 15 16;do #psi2S scan
#for pTMin in 10 11 12 13 14 16 18 20 22 24 26 28 30 32 34 36 38 40 43 46 49 54;do #ups3S scan
#for pTMin in 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19;do  #allCharm scan
#for pTMin in 14;do
pTMax=200

### store here the output TTree of the likelihood sampling, and all Figures of the results:
#for JobID in pTMinPlots_allCharmScan4t18_Chi2LDMEs;do
#for JobID in pTMinPlots_Ups3S_March3_pToverM;do
#for JobID in BKtruth_NewPwaveCCs;do

#for JobID in ScanPsi2S_PrelData_pTmin${pTMin}_February7;do
#for JobID in ScanUps3S_PrelData_pTmin${pTMin}_February7;do

#for JobID in BestFitPsi2S_pTmin${pTMin}_January22_1e4BI_1e6steps;do
#for JobID in BestFitUps3S_pTmin${pTMin}_January22_1e4BI_1e6steps;do

#for JobID in BestFit_Psi2SpTmin12_Ups3SpTmin30_March1_1e4BI_5e5steps_DebugTHuncert;do
#for JobID in BestFit_inclPJ_Psi2SpTmin12_Ups3SpTmin30_March2_1e4BI_5e4steps_noTHuncert;do
for JobID in Test_JobID;do


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
run_PlotCompareDataModel=0

#individual macros:::
run_FitPtDists=0
run_PlotPPD=0
run_PlotPPDderivative=0
run_PlotPPD_vs_pTmin=0

reCompile=0 #dangerous in loops...

##################################
########## SELECTION #############
##################################
# To be added in SamplePPD and the plotting

rapMin=-10
rapMax=10
useSstatesOnly=false
usePstatesOnly=false
useCharmoniumOnly=false
useBottomoniumOnly=false
useOnlyState=999 #switch off by setting it to an int > N_STATES
useCrossSectionOnly=false
usePolarizationOnly=false

##################################
########## SETTINGS ##############
##################################

### ConvertNTupleToTTree

### ConvertDataInput

### ConvertModelInput
useToyModel=false
calcOnlyConstsStat=false

### CombineDataModel

### GenerateToyData

### PlotCompareDataModel
PredictionPlot=false
PlotVSpToverM=false
PredictionDataPlot=false

### SamplePPD
Minimizer=0		#0...MH, 1...Minuit
#nBurnIn=5000
#nSample=50000
nBurnIn=5000
nSample=50000
SampleNp=true
SampleNp_consts_star=true

### InterpretPPD
nSigma=1
MPValgo=4 		#1...mean,2...gauss,3...gauss-loop with chi2<2, 4...(new)...mode


##################################
############# CODE ###############
##################################

if [ ${reCompile} -eq 1 ]
then
touch src/*
fi
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
cp PlotCompareDataModel PlotCompareDataModel_${JobID}_PredictionPlot${PredictionPlot}_PredictionDataPlot${PredictionDataPlot}
cp PlotPPD PlotPPD_${JobID}
cp PlotPPDderivative PlotPPDderivative_${JobID}
cp PlotPPD_vs_pTmin PlotPPD_vs_pTmin_${JobID}


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
./GenerateToyData_${DataModelCombID} ${ModelSystScaleID}=ModelSystScaleID ${ModelID}=ModelID ${DataID}=DataID ${DataModelCombID}=DataModelCombID ${storagedir}=storagedir
cp interface/ToyData.h ${storagedir}/DataModelCombID/${DataModelCombID}/ToyData.h
fi
if [ ${run_SamplePPD} -eq 1 ]
then
./SamplePPD_${JobID} ${nBurnIn}nBurnIn ${nSample}nSample ${ModelID}=ModelID  ${JobID}=JobID ${DataModelCombID}=DataModelCombID ${pTMin}pTMin ${pTMax}pTMax ${rapMin}rapMin ${rapMax}rapMax ${useOnlyState}useOnlyState useCrossSectionOnly=${useCrossSectionOnly} usePolarizationOnly=${usePolarizationOnly} ${Minimizer}Minimizer useSstatesOnly=${useSstatesOnly} usePstatesOnly=${usePstatesOnly} useCharmoniumOnly=${useCharmoniumOnly} useBottomoniumOnly=${useBottomoniumOnly} ${storagedir}=storagedir SampleNp=${SampleNp} SampleNp_consts_star=${SampleNp_consts_star} ${MPValgo}MPValgo ${nSigma}nSigma
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
./PlotCompareDataModel_${JobID}_PredictionPlot${PredictionPlot}_PredictionDataPlot${PredictionDataPlot} PredictionDataPlot=${PredictionDataPlot} PlotVSpToverM=${PlotVSpToverM} PredictionPlot=${PredictionPlot} ${JobID}=JobID ${ModelID}=ModelID ${nSigma}nSigma ${MPValgo}MPValgo ${storagedir}=storagedir ${DataModelCombID}=DataModelCombID ${pTMin}pTMin ${pTMax}pTMax ${rapMin}rapMin ${rapMax}rapMax ${useOnlyState}useOnlyState useCrossSectionOnly=${useCrossSectionOnly} usePolarizationOnly=${usePolarizationOnly} useSstatesOnly=${useSstatesOnly} usePstatesOnly=${usePstatesOnly} useCharmoniumOnly=${useCharmoniumOnly} useBottomoniumOnly=${useBottomoniumOnly} SampleNp=${SampleNp} SampleNp_consts_star=${SampleNp_consts_star}
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

if [ ${run_PlotPPDderivative} -eq 1 ]
then
echo 'run PlotPPDderivative'
./PlotPPDderivative_${JobID} ${JobID}=JobID ${storagedir}=storagedir ${useOnlyState}useOnlyState ${MPValgo}MPValgo ${nSigma}nSigma
mkdir Figures/PlotPPDderivative/
mkdir Figures/PlotPPDderivative/${JobID}
fi

if [ ${run_PlotPPD_vs_pTmin} -eq 1 ]
then
echo 'run PlotPPD_vs_pTmin'
./PlotPPD_vs_pTmin_${JobID} ${JobID}=JobID PlotVSpToverM=${PlotVSpToverM} ${storagedir}=storagedir ${useOnlyState}useOnlyState ${MPValgo}MPValgo ${nSigma}nSigma
mkdir Figures/PlotPPD_vs_pTmin/
mkdir Figures/PlotPPD_vs_pTmin/${JobID}
fi



rm ConvertDataInput_${DataID}
rm ConvertNTupleToTTree_${ModelID}
rm ConvertBKmodelToTTree_${ModelID}
rm ConvertModelInput_${ModelID}
rm CombineDataModel_${DataModelCombID}
rm GenerateToyData_${DataModelCombID}
rm SamplePPD_${JobID}
rm InterpretPPD_${JobID}
rm PlotCompareDataModel_${JobID}_PredictionPlot${PredictionPlot}_PredictionDataPlot${PredictionDataPlot}
rm PlotPPD_${JobID}
rm PlotPPDderivative_${JobID}
rm PlotPPD_vs_pTmin_${JobID}
rm ScaleBKmodel_${ModelID}
rm ScaleGWWZmodel_${ModelID}
rm FitPtDists_${DataID}

done
done
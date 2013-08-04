#!/bin/sh

storagedir=/scratch/knuenz/NRQCD/NRQCDglobalfit

OriginalNTupleID=BaranovSmallFile #This is the location of the original model nTuple (provided by Sergey)
DataID=August4_Psi2Sonly #store here the NRQCDglobalfitObjects containing data measurements
ModelID=July28_ToyAddPol #store here the ModelIngredients.root file and consts_star file
DataModelCombID=July29_ToyDataToyModel_ToyAddPol_polCorr #store here the NRQCD objects combining data and model predictions
for JobID in August4_psi2Sonly_SampleNp_MH;do #store here the output TTree of the likelihood sampling, and all Figures of the results

run_ConvertDataInput=1
run_ConvertNTupleToTTree=0
run_ConvertModelInput=0
run_CombineDataModel=0
run_GenerateToyData=0
run_SamplePPD=0
run_InterpretPPD=0
run_PlotCompareDataModel=0

##################################
########## SETTINGS ##############
##################################

### ConvertNTupleToTTree

### ConvertDataInput

### ConvertModelInput
useToyModel=true

### CombineDataModel

### GenerateToyData

### SamplePPD
nBurnIn=1000
nSample=10000

### InterpretPPD
nSigma=1
MPValgo=3 		#1...mean,2...gauss,3...gauss-loop with chi2<2

### PlotCompareDataModel

##################################
########## SELECTION #############
##################################
# To be added in SamplePPD and the plotting

pTMin=0
pTMax=100
rapMin=-10
rapMax=10
useSstatesOnly=false
usePstatesOnly=false
useCharmoniumOnly=true
useBottomoniumOnly=false
useOnlyState=3 #switch off by setting it to an int > N_STATES

##################################
############# CODE ###############
##################################

touch src/*
make

if [ ${run_ConvertDataInput} -eq 1 ]
then
./ConvertDataInput ${DataID}=DataID ${storagedir}=storagedir
fi
if [ ${run_ConvertNTupleToTTree} -eq 1 ]
then
./ConvertNTupleToTTree ${OriginalNTupleID}=OriginalNTupleID ${ModelID}=ModelID ${storagedir}=storagedir
fi
if [ ${run_ConvertModelInput} -eq 1 ]
then
./ConvertModelInput ${ModelID}=ModelID useToyModel=${useToyModel} ${storagedir}=storagedir
fi
if [ ${run_CombineDataModel} -eq 1 ]
then
./CombineDataModel ${ModelID}=ModelID ${DataID}=DataID ${DataModelCombID}=DataModelCombID ${storagedir}=storagedir
fi
if [ ${run_GenerateToyData} -eq 1 ]
then
./GenerateToyData ${ModelID}=ModelID ${DataID}=DataID ${DataModelCombID}=DataModelCombID ${storagedir}=storagedir
cp interface/ToyData.h ${storagedir}/DataModelCombID/${DataModelCombID}/ToyData.h
fi
if [ ${run_SamplePPD} -eq 1 ]
then
./SamplePPD ${nBurnIn}nBurnIn ${nSample}nSample ${ModelID}=ModelID  ${JobID}=JobID ${DataModelCombID}=DataModelCombID ${pTMin}pTMin ${pTMax}pTMax ${rapMin}rapMin ${rapMax}rapMax ${useOnlyState}useOnlyState useSstatesOnly=${useSstatesOnly} usePstatesOnly=${usePstatesOnly} useCharmoniumOnly=${useCharmoniumOnly} useBottomoniumOnly=${useBottomoniumOnly} ${storagedir}=storagedir
fi
if [ ${run_InterpretPPD} -eq 1 ]
then
./InterpretPPD ${JobID}=JobID ${MPValgo}MPValgo ${nSigma}nSigma ${storagedir}=storagedir
fi
if [ ${run_PlotCompareDataModel} -eq 1 ]
then
./PlotCompareDataModel ${JobID}=JobID ${ModelID}=ModelID ${nSigma}nSigma ${storagedir}=storagedir ${DataModelCombID}=DataModelCombID ${pTMin}pTMin ${pTMax}pTMax ${rapMin}rapMin ${rapMax}rapMax ${useOnlyState}useOnlyState useSstatesOnly=${useSstatesOnly} usePstatesOnly=${usePstatesOnly} useCharmoniumOnly=${useCharmoniumOnly} useBottomoniumOnly=${useBottomoniumOnly}
fi

done

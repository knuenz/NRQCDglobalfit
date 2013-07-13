#!/bin/sh

OriginalNTupleID=BaranovSmallFile #This is the location of the original model nTuple (provided by Sergey)
DataID=July13_Data #store here the NRQCDglobalfitObjects containing data and model prediction
ModelID=July13_Model #store here the ModelIngredients.root file
JobID=July13_Sample #store here the output TTree of the likelihood sampling

run_ConvertDataInput=0
run_ConvertNTupleToTTree=0
run_ConvertModelInput=0
run_CombineDataModel=1
run_SamplePPD=0

##################################
########## SETTINGS ##############
##################################

### ConvertNTupleToTTree

### ConvertDataInput

### ConvertModelInput
useToyModel=false

### CombineDataModel

### SamplePPD
nSample=10000
nBurnIn=2000

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
useCharmoniumOnly=false
useBottomoniumOnly=false
useOnlyState=999 #switch off by setting it to an int > N_STATES

##################################
############# CODE ###############
##################################

touch src/*
make

if [ ${run_ConvertDataInput} -eq 1 ]
then
./ConvertDataInput ${DataID}=DataID
fi
if [ ${run_ConvertNTupleToTTree} -eq 1 ]
then
./ConvertNTupleToTTree ${OriginalNTupleID}=OriginalNTupleID ${ModelID}=ModelID
fi
if [ ${run_ConvertModelInput} -eq 1 ]
then
./ConvertModelInput ${ModelID}=ModelID useToyModel=${useToyModel}
fi
if [ ${run_CombineDataModel} -eq 1 ]
then
./CombineDataModel ${ModelID}=ModelID ${DataID}=DataID
fi
if [ ${run_SamplePPD} -eq 1 ]
then
./SamplePPD ${JobID}=JobID ${DataID}=DataID ${pTMin}pTMin ${pTMax}pTMax ${rapMin}rapMin ${rapMax}rapMax ${useOnlyState}useOnlyState useSstatesOnly=${useSstatesOnly} usePstatesOnly=${usePstatesOnly} useCharmoniumOnly=${useCharmoniumOnly} useBottomoniumOnly=${useBottomoniumOnly}
fi

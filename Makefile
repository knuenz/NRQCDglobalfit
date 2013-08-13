CXX=$(shell root-config --cxx --cflags)
LIBS=$(shell root-config --libs)

CXX+=-I/afs/cern.ch/cms/sw/slc5_ia32_gcc434/lcg/roofit/5.26.00-cms5/include
LIBS+=-L/afs/cern.ch/cms/sw/slc5_ia32_gcc434/lcg/roofit/5.26.00-cms5/lib
#CXX+=-I/Users/valentinknuenz/usr/local/root/roofit/roofit/inc
#LIBS+=-L/Users/valentinknuenz/usr/local/root/roofit/roofit/inc

#CXX+=-I/Users/valentinknuenz/usr/local/root/lib
#LIBS+=-L/Users/valentinknuenz/usr/local/root/lib


#DATA=

%.o : %.cc
	$(CXX) -c $<


all: SamplePPD ConvertDataInput ConvertModelInput CombineDataModel ConvertNTupleToTTree GenerateToyData InterpretPPD PlotCompareDataModel ConvertBKmodelToTTree

SamplePPD: src/SamplePPD.cc
	$(CXX) $^ -o $@ $(LIBS) -lRooFit -lRooFitCore -lFoam -lMinuit  -lMinuit2

ConvertDataInput: src/ConvertDataInput.cc
	$(CXX) $^ -o $@ $(LIBS) -lRooFit -lRooFitCore -lFoam -lMinuit

ConvertModelInput: src/ConvertModelInput.cc
	$(CXX) $^ -o $@ $(LIBS) -lRooFit -lRooFitCore -lFoam -lMinuit

CombineDataModel: src/CombineDataModel.cc
	$(CXX) $^ -o $@ $(LIBS) -lRooFit -lRooFitCore -lFoam -lMinuit

ConvertNTupleToTTree: src/ConvertNTupleToTTree.cc
	$(CXX) $^ -o $@ $(LIBS) -lRooFit -lRooFitCore -lFoam -lMinuit

GenerateToyData: src/GenerateToyData.cc
	$(CXX) $^ -o $@ $(LIBS) -lRooFit -lRooFitCore -lFoam -lMinuit

InterpretPPD: src/InterpretPPD.cc
	$(CXX) $^ -o $@ $(LIBS) -lRooFit -lRooFitCore -lFoam -lMinuit

PlotCompareDataModel: src/PlotCompareDataModel.cc
	$(CXX) $^ -o $@ $(LIBS) -lRooFit -lRooFitCore -lFoam -lMinuit

ConvertBKmodelToTTree: src/ConvertBKmodelToTTree.cc
	$(CXX) $^ -o $@ $(LIBS) -lRooFit -lRooFitCore -lFoam -lMinuit
					
								
clean: 
	rm SamplePPD ConvertDataInput ConvertModelInput CombineDataModel ConvertNTupleToTTree GenerateToyData InterpretPPD PlotCompareDataModel ConvertBKmodelToTTree*.o
	

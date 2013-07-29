/*
 * PlotCompareDataModel.cc
 *
 *  Created on: Jul 27, 2013
 *      Author: valentinknuenz
 */



#include "../interface/NRQCDglobalfitObject.h"
#include "../src/NRQCDglobalfitObject.cc"
#include "../interface/GlobalVar.h"

#include <iostream>
#include <sstream>
#include <iomanip>
#include <iterator>
#include <vector>
#include <string>

//rootincludes
#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TSystem.h"
#include "TBranch.h"
#include "TLine.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TGraphAsymmErrors.h"
#include "TGraph.h"
#include "TROOT.h"


using namespace NRQCDvars;


void plotComp(int iState, int iMeasurementID, int iExperiment, int iRap, char jobdirname[200], char rapchar[200], TGraphAsymmErrors *data_Graph, TGraph *model_Graph, TGraph *model_Graph_low, TGraph *model_Graph_high, bool plotInclusiveFeedDown, bool plotIndividualFeedDown, bool plotDirectColorChannels, const int StatesCont_c, const int ColorChannels_c, TGraph *model_Graph_FeedDowns[500], TGraph *model_Graph_ColorChannels[500], vector<int> StatesContributing);
vector<double> addPolarizations(vector<vector<double> > lamMatrix, vector<double> lamVecContributions);

int main(int argc, char** argv) {

  	Char_t *JobID = "Default";
  	Char_t *ModelID = "Default";
  	Char_t *DataModelCombID = "Default";
	Char_t *storagedir = "Default"; //Storage Directory

  	double 	pTMin=0;
  	double 	pTMax=1000.;
  	double 	rapMin=-10.;
  	double 	rapMax=10.;
  	bool 	useSstatesOnly=false;
  	bool 	usePstatesOnly=false;
  	bool 	useCharmoniumOnly=false;
  	bool 	useBottomoniumOnly=false;
  	int 	useOnlyState=999;

  	int 	MPValgo=-1;
  	double 	nSigma=-1;


  	for( int i=0;i < argc; ++i ) {
		if(std::string(argv[i]).find("storagedir") != std::string::npos) {char* storagedirchar = argv[i]; char* storagedirchar2 = strtok (storagedirchar, "="); storagedir = storagedirchar2; cout<<"storagedir = "<<storagedir<<endl;}
	    if(std::string(argv[i]).find("nSigma") != std::string::npos) { char* nSigmachar = argv[i]; char* nSigmachar2 = strtok (nSigmachar, "n"); nSigma = atof(nSigmachar2); cout<<"nSigma = "<<nSigma<<endl; }
	    if(std::string(argv[i]).find("MPValgo") != std::string::npos) { char* MPValgochar = argv[i]; char* MPValgochar2 = strtok (MPValgochar, "M"); MPValgo = atoi(MPValgochar2); cout<<"MPValgo = "<<MPValgo<<endl; }
		if(std::string(argv[i]).find("JobID") != std::string::npos) { char* JobIDchar = argv[i]; char* JobIDchar2 = strtok (JobIDchar, "="); JobID = JobIDchar2; cout<< "JobID = " << JobID << endl; }
		if(std::string(argv[i]).find("ModelID") != std::string::npos) {char* ModelIDchar = argv[i]; char* ModelIDchar2 = strtok (ModelIDchar, "="); ModelID = ModelIDchar2; cout<<"ModelID = "<<ModelID<<endl;}
		if(std::string(argv[i]).find("DataModelCombID") != std::string::npos) {char* DataModelCombIDchar = argv[i]; char* DataModelCombIDchar2 = strtok (DataModelCombIDchar, "="); DataModelCombID = DataModelCombIDchar2; cout<<"DataModelCombID = "<<DataModelCombID<<endl;}
	    if(std::string(argv[i]).find("pTMin") != std::string::npos) {
	    	char* pTMinchar = argv[i];
	    	char* pTMinchar2 = strtok (pTMinchar, "p");
	    	pTMin = atof(pTMinchar2);
	    	cout<<"pTMin = "<<pTMin<<endl;
	    }
	    if(std::string(argv[i]).find("pTMax") != std::string::npos) {
	    	char* pTMaxchar = argv[i];
	    	char* pTMaxchar2 = strtok (pTMaxchar, "p");
	    	pTMax = atof(pTMaxchar2);
	    	cout<<"pTMax = "<<pTMax<<endl;
	    }
	    if(std::string(argv[i]).find("rapMin") != std::string::npos) {
	    	char* rapMinchar = argv[i];
	    	char* rapMinchar2 = strtok (rapMinchar, "r");
	    	rapMin = atof(rapMinchar2);
	    	cout<<"rapMin = "<<rapMin<<endl;
	    }
	    if(std::string(argv[i]).find("rapMax") != std::string::npos) {
	    	char* rapMaxchar = argv[i];
	    	char* rapMaxchar2 = strtok (rapMaxchar, "r");
	    	rapMax = atof(rapMaxchar2);
	    	cout<<"rapMax = "<<rapMax<<endl;
	    }
	    if(std::string(argv[i]).find("useOnlyState") != std::string::npos) {
	    	char* useOnlyStatechar = argv[i];
	    	char* useOnlyStatechar2 = strtok (useOnlyStatechar, "u");
	    	useOnlyState = atoi(useOnlyStatechar2);
	    	cout<<"useOnlyState = "<<useOnlyState<<endl;
	    }
	    if(std::string(argv[i]).find("useSstatesOnly=true") != std::string::npos) {
	    	useSstatesOnly=true;
	    	cout<<"useSstatesOnly=true"<<endl;
	    }
	    if(std::string(argv[i]).find("usePstatesOnly=true") != std::string::npos) {
	    	usePstatesOnly=true;
	    	cout<<"usePstatesOnly=true"<<endl;
	    }
	    if(std::string(argv[i]).find("useCharmoniumOnly=true") != std::string::npos) {
	    	useCharmoniumOnly=true;
	    	cout<<"useCharmoniumOnly=true"<<endl;
	    }
	    if(std::string(argv[i]).find("useBottomoniumOnly=true") != std::string::npos) {
	    	useBottomoniumOnly=true;
	    	cout<<"useBottomoniumOnly=true"<<endl;
	    }

  	}


	char inname[200];
	char outname[200];
	char jobdirname[200];
	char modeldirname[200];
	char datamodeldirname[200];
	sprintf(jobdirname,"%s/JobID", storagedir);
	gSystem->mkdir(jobdirname);
	sprintf(jobdirname,"%s/%s",jobdirname,JobID);
	gSystem->mkdir(jobdirname);
	sprintf(modeldirname,"%s/ModelID", storagedir);
	gSystem->mkdir(modeldirname);
	sprintf(modeldirname,"%s/%s",modeldirname,ModelID);
	gSystem->mkdir(modeldirname);
	sprintf(datamodeldirname,"%s/DataModelCombID", storagedir);
	gSystem->mkdir(datamodeldirname);
	sprintf(datamodeldirname,"%s/%s",datamodeldirname,DataModelCombID);
	gSystem->mkdir(datamodeldirname);

	sprintf(inname,"%s/results.txt",jobdirname);
	cout<<"read results from "<<inname<<endl;

  	dmatrix fi_MPV;
  	dmatrix fi_errlow;
 	dmatrix fi_errhigh;
  	dmatrix Oi_MPV;
  	dmatrix Oi_errlow;
 	dmatrix Oi_errhigh;

    ifstream in;
    in.open(inname);//, std::ofstream::app);

	in >> fi_MPV;
	in >> fi_errlow;
	in >> fi_errhigh;
	in >> Oi_MPV;
	in >> Oi_errlow;
	in >> Oi_errhigh;

    in.close();


    //dmatrix consts_star;
	//ifstream instar;
	//instar.open(inname);
	//instar >> consts_star;
	//instar.close();
    //
	//cout << consts_star << endl;



    cout<<"fi_MPV:"<<endl;
	cout<<fi_MPV<<endl;
    cout<<"fi_errlow:"<<endl;
	cout<<fi_errlow<<endl;
    cout<<"fi_errhigh:"<<endl;
	cout<<fi_errhigh<<endl;
    cout<<"Oi_MPV:"<<endl;
	cout<<Oi_MPV<<endl;
    cout<<"Oi_errlow:"<<endl;
	cout<<Oi_errlow<<endl;
    cout<<"Oi_errhigh:"<<endl;
	cout<<Oi_errhigh<<endl;



    cout<<"Inizialize Op_plus:"<<endl;

	dmatrix Op_plus(NRQCDvars::nStates);
	vector<double> Op_S_plus (NRQCDvars::nColorChannels_S,0);
	vector<double> Op_P_plus (NRQCDvars::nColorChannels_P,0);

	for (int i=0; i<NRQCDvars::nStates; i++){
		bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
		if(isSstate){
			for (int j=0; j<NRQCDvars::nColorChannels_S; j++){
				Op_S_plus.at(j)=Oi_MPV[i][j]+Oi_errhigh[i][j];
			}
			Op_plus.at(i)=Op_S_plus;
		}
		else{
			for (int j=0; j<NRQCDvars::nColorChannels_P; j++){
				Op_P_plus.at(j)=Oi_MPV[i][j]+Oi_errhigh[i][j];
			}
			Op_plus.at(i)=Op_P_plus;
		}
	}

    cout<<"Inizialize Op_minus:"<<endl;

    dmatrix Op_minus(NRQCDvars::nStates);
	vector<double> Op_S_minus (NRQCDvars::nColorChannels_S,0);
	vector<double> Op_P_minus (NRQCDvars::nColorChannels_P,0);

	for (int i=0; i<NRQCDvars::nStates; i++){
		bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
		if(isSstate){
			for (int j=0; j<NRQCDvars::nColorChannels_S; j++){
				Op_S_minus.at(j)=Oi_MPV[i][j]-Oi_errlow[i][j];
			}
			Op_minus.at(i)=Op_S_minus;
		}
		else{
			for (int j=0; j<NRQCDvars::nColorChannels_P; j++){
				Op_P_minus.at(j)=Oi_MPV[i][j]-Oi_errlow[i][j];
			}
			Op_minus.at(i)=Op_P_minus;
		}
	}

    cout<<"Op_plus:"<<endl;
	cout<<Op_plus<<endl;
    cout<<"Op_minus:"<<endl;
	cout<<Op_minus<<endl;

    cout<<"Inizialize Np's:"<<endl;

    dmatrix Np_BR; //Branching ratios, [nDauthgers][nMothers]
	dmatrix Np_US; //Uncertainty scales, [0=Data, 1=Model][nScales]

	dvector Np_BR_0 (NRQCDvars::nStates,0.);
	for(int i=0; i < NRQCDvars::nStates; i++) Np_BR.push_back(Np_BR_0);

	dvector Np_US_0 (NRQCDvars::nDataSystematicScales, 0.);
	dvector Np_US_1 (NRQCDvars::nModelSystematicScales, 0.);
	Np_US.push_back(Np_US_0);
	Np_US.push_back(Np_US_1);

	for(int j=0; j < NRQCDvars::nDataSystematicScales; j++){
		Np_US[0][j]=0.;
	}
	for(int j=0; j < NRQCDvars::nModelSystematicScales; j++){
		Np_US[1][j]=0.;
	}

	for(int i=0; i < nStates; i++){
		for(int j=0; j < nStates; j++){
			if(NRQCDvars::FeedDownBranchingRatio[i][j]>0){
				Np_BR[i][j]=NRQCDvars::FeedDownBranchingRatio[i][j]*0.01;
			}
			else Np_BR[i][j]=0;
		}
	}



    cout<<"Loop through data to plot:"<<endl;

	bool DataSelected=false;

	for(int iState=0; iState < nStates; iState++){
		for(int iMeasurementID=0; iMeasurementID < NRQCDvars::nMeasurementIDs; iMeasurementID++){
			for(int iExperiment=0; iExperiment < NRQCDvars::nExperiments; iExperiment++){


				NRQCDglobalfitObject *DataModelObject[NRQCDvars::nMaxRapBins][NRQCDvars::nMaxPtBins];
				bool DataPresentAndSelected[NRQCDvars::nMaxRapBins][NRQCDvars::nMaxPtBins];

				bool definedRap=false;
				bool isAbsRap=false;
				double rapMinObject, rapMaxObject;

				int StatesCont;//Number of states contributing to result, including direct state

				for(int iRap = 0; iRap < NRQCDvars::nMaxRapBins; iRap++){
			    	cout<<"plot model for iState="<<StateName[iState]<<", iMeasurementID="<<MeasurementIDName[iMeasurementID]<<", iExperiment="<<ExpName[iExperiment]<<", iRap="<<iRap/*", iP="<<iP*/<<endl;

					definedRap=false;
					int nPtBinsSel=0;
				    for(int iP = 0; iP < NRQCDvars::nMaxPtBins; iP++){

				    	DataPresentAndSelected[iRap][iP]=false;

				    	//if(iState!=0 || iMeasurementID!=0 || iExperiment!=0 || iRap!=0 || iP!=0) continue;

						sprintf(inname,"%s/ConvertedDataModel_%s_%s_%s_rap%d_pT%d.txt",datamodeldirname, StateName[iState],
								MeasurementIDName[iMeasurementID],  ExpName[iExperiment], iRap, iP);

						ifstream inData;
						inData.open(inname);

						if( inData.is_open() ){//Measurement present -> calculate model components:: Modified by Joao: more correct from ios point of view
							//if(NRQCDvars::debug){
								cout << "Read in iState=" << StateName[iState] << ", iMeasurementID=" << MeasurementIDName[iMeasurementID];
								cout << ", iExperiment=" << ExpName[iExperiment] << ", iRap=" << iRap << ", iP=" << iP << endl;
							//}
							DataModelObject[iRap][iP] = new NRQCDglobalfitObject();
							inData >> *DataModelObject[iRap][iP];

							if(NRQCDvars::debug) DataModelObject[iRap][iP]->Dump(nStates, true, true);
							inData.close();

							if(
									DataModelObject[iRap][iP]->getpTMin() >= pTMin
								&& DataModelObject[iRap][iP]->getpTMax() <= pTMax
								&& DataModelObject[iRap][iP]->getyMin() >= rapMin
								&& DataModelObject[iRap][iP]->getyMax() <= rapMax
							) DataSelected=true;


							if(useSstatesOnly && StateQuantumID[DataModelObject[iRap][iP]->getState()] != NRQCDvars::quID_S)
								DataSelected=false;
							if(usePstatesOnly && StateQuantumID[DataModelObject[iRap][iP]->getState()] == NRQCDvars::quID_S)
								DataSelected=false;
							if(useCharmoniumOnly && DataModelObject[iRap][iP]->getState() > 3)
								DataSelected=false;
							if(useBottomoniumOnly && DataModelObject[iRap][iP]->getState() < 4)
								DataSelected=false;
							if(useOnlyState < nStates &&  DataModelObject[iRap][iP]->getState() != useOnlyState)
								DataSelected=false;


							if(DataSelected){

								DataPresentAndSelected[iRap][iP]=true;
								nPtBinsSel++;
								if(!definedRap){
									isAbsRap=DataModelObject[iRap][iP]->getisAbsRap();
									rapMinObject=DataModelObject[iRap][iP]->getyMin();
									rapMaxObject=DataModelObject[iRap][iP]->getyMax();
									definedRap=true;

								}
							}
							DataSelected=false;
						}
				    }//pT loop

				    dvector data_centralval;
				    dvector data_errlow_centralval;
				    dvector data_errhigh_centralval;
				    dvector data_pTmean;
				    dvector data_errlow_pT;
				    dvector data_errhigh_pT;

				    dvector model_centralval;
				    dvector model_errlow_centralval;
				    dvector model_errhigh_centralval;

				    dvector model_inclusive_FeedDown;

				    dmatrix model_directProduction; //(nColorChannels, iP)
				    dmatrix model_FeedDown; //(StatesCont, iP)

				    //TODO: code unbinned models
				    //TODO: code submodels -> color channels and feed-down

				    // -> getDirectProduction: add to dvector(4) return dvector CS(nColorC) + dmatrix of lambda(nColorC, nLambda)
				    // -> getPromptProduction: add to dvector(4) return dvector CS(nStates) + dmatrix of lambda(nStates, nLambda)
				    // -> getPromptCorrProduction: add return polCorrFactor, to be used to correct the data -> show data-points polarization-corrected

				    // dcube (nStates, nColorChannels, 4) directProductionCube
				    // dmatrix (nStates, 4) promptProductionCube
				    // double polCorrFactor

				    //Plot: all feed-downs (CC inclusive) individually + inclusive feed-down + all CCs from direct state

		    		vector<int> StatesContributing_;
		    		vector<double>  polCorrFactorVec;

				    for(int iP = 0; iP < NRQCDvars::nMaxPtBins; iP++){

				    	if(DataPresentAndSelected[iRap][iP]){


							dvector ObjectLikelihoodVec;
							double ModelPrediction;
							double errModelPrediction;
							//TODO: Improve model uncertainty (draw from PPD or so)

							dcube directProductionCube;
							dmatrix promptProductionMatrix;
							double polCorrFactor;

							dcube Dummy_directProductionCube;
							dmatrix Dummy_promptProductionMatrix;
							double Dummy_polCorrFactor;

							ObjectLikelihoodVec=DataModelObject[iRap][iP]->getObjectLikelihood(Oi_MPV, Np_BR, Np_US, true, directProductionCube, promptProductionMatrix, polCorrFactor);
							ModelPrediction=ObjectLikelihoodVec[1];
							model_centralval.push_back(ModelPrediction);

							//cout<<"directProductionCube:"<<endl;
							//cout<<directProductionCube<<endl;
							//cout<<"promptProductionMatrix:"<<endl;
							//cout<<promptProductionMatrix<<endl;
							//cout<<"polCorrFactor:"<<endl;
							//cout<<polCorrFactor<<endl;

							polCorrFactorVec.push_back(polCorrFactor);

							ObjectLikelihoodVec=DataModelObject[iRap][iP]->getObjectLikelihood(Op_minus, Np_BR, Np_US, false, Dummy_directProductionCube, Dummy_promptProductionMatrix, Dummy_polCorrFactor);
							errModelPrediction=ObjectLikelihoodVec[1];
							model_errlow_centralval.push_back(fabs(ModelPrediction-errModelPrediction));

							ObjectLikelihoodVec=DataModelObject[iRap][iP]->getObjectLikelihood(Op_plus, Np_BR, Np_US, false, Dummy_directProductionCube, Dummy_promptProductionMatrix, Dummy_polCorrFactor);
							errModelPrediction=ObjectLikelihoodVec[1];
							model_errhigh_centralval.push_back(fabs(ModelPrediction-errModelPrediction));


				    		data_centralval.push_back(DataModelObject[iRap][iP]->getCentralValue()/polCorrFactor);//correct data to match predicted polariztion
				    		data_errlow_centralval.push_back(DataModelObject[iRap][iP]->getErrTotNeg()/polCorrFactor);
				    		data_errhigh_centralval.push_back(DataModelObject[iRap][iP]->getErrTotPos()/polCorrFactor);
				    		data_pTmean.push_back(DataModelObject[iRap][iP]->getpTMean());
				    		data_errlow_pT.push_back(DataModelObject[iRap][iP]->getpTMean()-DataModelObject[iRap][iP]->getpTMin());
				    		data_errhigh_pT.push_back(DataModelObject[iRap][iP]->getpTMax()-DataModelObject[iRap][iP]->getpTMean());

				    		//cout<<"Calculate Inclusive FeedDown"<<endl;
				    		//Calculate Inclusive FeedDown
				    		dvector FeedDownContVec;
				    		//dvector FeedDownLamthVec;
				    		//dvector FeedDownLamphVec;
				    		//dvector FeedDownLamtpVec;
				    		dvector FeedDownLamVec(3,0);
				    		dmatrix FeedDownLamMatrix;
				    		vector<int> StatesContributing;
				    		double FeedDownCont=0;
				    		double BufferCrossSection;
				    		for(int i=0; i<promptProductionMatrix.size();i++){
				    			dvector promptProductionVector=promptProductionMatrix.at(i);
				    			BufferCrossSection = promptProductionVector.at(0);
				    			if(BufferCrossSection>1e-100) StatesContributing.push_back(i);
				    			if(i!=iState){
									FeedDownCont+=BufferCrossSection;
									FeedDownContVec.push_back(BufferCrossSection);
									FeedDownLamVec.at(0)=promptProductionVector.at(1);
									FeedDownLamVec.at(1)=promptProductionVector.at(2);
									FeedDownLamVec.at(2)=promptProductionVector.at(3);
									FeedDownLamMatrix.push_back(FeedDownLamVec);
									//FeedDownLamthVec.push_back(promptProductionVector.at(1));
									//FeedDownLamphVec.push_back(promptProductionVector.at(2));
									//FeedDownLamtpVec.push_back(promptProductionVector.at(3));
				    			}
				    		}



				    		StatesCont=StatesContributing.size();
						    StatesContributing_=StatesContributing;

							 double model_inclusive_FeedDown_=0;

				    		if(StatesCont>2){

								dvector FeedDownInclusiveLambdas;
								 FeedDownInclusiveLambdas = addPolarizations(FeedDownLamMatrix, FeedDownContVec);



								if(iMeasurementID==0) model_inclusive_FeedDown_=FeedDownCont;
								else if(iMeasurementID>0 && iMeasurementID<4){
									model_inclusive_FeedDown_=FeedDownInclusiveLambdas[iMeasurementID-1];
								}

				    		}

				    		//cout<<"Calculate individual FeedDown contributions"<<endl;
				    		//Calculate individual FeedDown contributions
						    dvector Buff_model_FeedDownVec;
						    dvector model_FeedDownVec;
						    for(int i=0; i<StatesCont;i++){
						    	if(i==0){
						    		if(StatesCont>2) model_FeedDownVec.push_back(model_inclusive_FeedDown_);
						    		if(StatesCont==2){
										Buff_model_FeedDownVec=promptProductionMatrix.at(StatesContributing[1]);
										model_FeedDownVec.push_back(Buff_model_FeedDownVec[iMeasurementID]);
						    		}
						    	}
						    	else{
									Buff_model_FeedDownVec=promptProductionMatrix.at(StatesContributing[i]);
									model_FeedDownVec.push_back(Buff_model_FeedDownVec[iMeasurementID]);
						    	}
						    }
						    model_FeedDown.push_back(model_FeedDownVec);


				    		//cout<<"Calculate individual color channel contributions of the directly produced state"<<endl;
				    		//Calculate individual color channel contributions of the directly produced state

						    dmatrix Buff_model_directProductionMatrix=directProductionCube.at(0);
						    dvector model_directProductionVec;
						    model_directProductionVec=Buff_model_directProductionMatrix.at(iMeasurementID);
						    model_directProduction.push_back(model_directProductionVec);


				    	}

				    }


				    if(nPtBinsSel>0){
						cout << "Plot iState=" << StateName[iState] << ", iMeasurementID=" << MeasurementIDName[iMeasurementID] << " from "<< ExpName[iExperiment] << endl;
				    	cout<<"rap "<<iRap+1<<endl;

						const int nPtBinsSel_=nPtBinsSel;
						double d_data_centralval[nPtBinsSel_];
						double d_data_errlow_centralval[nPtBinsSel_];
						double d_data_errhigh_centralval[nPtBinsSel_];
						double d_data_pTmean[nPtBinsSel_];
						double d_data_errlow_pT[nPtBinsSel_];
						double d_data_errhigh_pT[nPtBinsSel_];
						double d_model_centralval[nPtBinsSel_];
						double d_model_errlow_centralval[nPtBinsSel_];
						double d_model_errhigh_centralval[nPtBinsSel_];
						double d_model_errlow_absolute_centralval[nPtBinsSel_];
						double d_model_errhigh_absolute_centralval[nPtBinsSel_];

						for(int iP = 0; iP < nPtBinsSel; iP++){

							d_data_centralval[iP] =         	data_centralval[iP];
							d_data_errlow_centralval[iP] =  	data_errlow_centralval[iP];
							d_data_errhigh_centralval[iP] = 	data_errhigh_centralval[iP];
							d_data_pTmean[iP] =             	data_pTmean[iP];
							d_data_errlow_pT[iP] =          	data_errlow_pT[iP];
							d_data_errhigh_pT[iP] =         	data_errhigh_pT[iP];
							d_model_centralval[iP] =        	model_centralval[iP];
							d_model_errlow_centralval[iP] = 	model_errlow_centralval[iP];
							d_model_errhigh_centralval[iP] =	model_errhigh_centralval[iP];
							d_model_errlow_absolute_centralval[iP] = 	model_centralval[iP]-model_errlow_centralval[iP];
							d_model_errhigh_absolute_centralval[iP] =	model_centralval[iP]+model_errhigh_centralval[iP];


						}

						TGraphAsymmErrors *data_Graph = new TGraphAsymmErrors(nPtBinsSel,d_data_pTmean,d_data_centralval,d_data_errlow_pT,d_data_errhigh_pT,d_data_errlow_centralval,d_data_errhigh_centralval);
						TGraph *model_Graph = new TGraphAsymmErrors(nPtBinsSel,d_data_pTmean,d_model_centralval);
						TGraph *model_Graph_low = new TGraphAsymmErrors(nPtBinsSel,d_data_pTmean,d_model_errlow_absolute_centralval);
						TGraph *model_Graph_high = new TGraphAsymmErrors(nPtBinsSel,d_data_pTmean,d_model_errhigh_absolute_centralval);

						//Make TGraphs for model contributions
						const int StatesCont_c=StatesCont;
						int nColorChannels_state;
						bool isSstate=(StateQuantumID[iState] > NRQCDvars::quID_S)?false:true;
						if(isSstate) nColorChannels_state=NRQCDvars::nColorChannels_S;
						else nColorChannels_state=NRQCDvars::nColorChannels_P;
						const int ColorChannels_c=nColorChannels_state;

						TGraph *model_Graph_FeedDowns[StatesCont_c];//fist element: inclusive FeedDown

						//cout<<model_FeedDown<<endl;
						//cout<<"StatesCont_c "<<StatesCont_c<<endl;
						for(int i=0;i<StatesCont_c;i++){
							for(int iP = 0; iP < nPtBinsSel; iP++){
								d_model_centralval[iP] =  model_FeedDown[iP][i];
								if(iMeasurementID==0) d_model_centralval[iP]/=polCorrFactorVec[iP];
							}
							model_Graph_FeedDowns[i] = new TGraphAsymmErrors(nPtBinsSel,d_data_pTmean,d_model_centralval);
							//model_Graph_FeedDowns[i]->Print();
						}


						TGraph *model_Graph_ColorChannels[ColorChannels_c];

						for(int i=0;i<ColorChannels_c;i++){
							for(int iP = 0; iP < nPtBinsSel; iP++){
								d_model_centralval[iP] =  model_directProduction[iP][i];
								if(iMeasurementID==0) d_model_centralval[iP]/=polCorrFactorVec[iP];
							}
							model_Graph_ColorChannels[i] = new TGraphAsymmErrors(nPtBinsSel,d_data_pTmean,d_model_centralval);
						}


						bool plotInclusiveFeedDown=true;
						if(StatesCont_c<3) plotInclusiveFeedDown=false;
						bool plotIndividualFeedDown=true;
						if(StatesCont_c<2) plotInclusiveFeedDown=false;
						bool plotDirectColorChannels=true;
						if(ColorChannels_c<2) plotDirectColorChannels=false;


						char rapchar[200];
						if(isAbsRap) sprintf(rapchar,"%1.1f < |y| < %1.1f",rapMinObject, rapMaxObject);
						else sprintf(rapchar,"%1.1f < y < %1.1f",rapMinObject, rapMaxObject);
						cout<<"plotComp"<<endl;
						plotComp( iState,  iMeasurementID,  iExperiment,  iRap,  jobdirname, rapchar, data_Graph, model_Graph, model_Graph_low, model_Graph_high, plotInclusiveFeedDown, plotIndividualFeedDown, plotDirectColorChannels, StatesCont_c, ColorChannels_c, model_Graph_FeedDowns, model_Graph_ColorChannels, StatesContributing_);

				}


				}//rap loop



			}
		}
	}











		return 0;

  	}



void plotComp(int iState, int iMeasurementID, int iExperiment, int iRap, char jobdirname[200], char rapchar[200], TGraphAsymmErrors *data_Graph, TGraph *model_Graph, TGraph *model_Graph_low, TGraph *model_Graph_high, bool plotInclusiveFeedDown, bool plotIndividualFeedDown, bool plotDirectColorChannels, const int StatesCont_c, const int ColorChannels_c, TGraph *model_Graph_FeedDowns[500], TGraph *model_Graph_ColorChannels[500], vector<int> StatesContributing){
	gROOT->Reset();
	gStyle->SetOptStat(11);
	gStyle->SetOptFit(101);
	gStyle->SetFillColor(kWhite);

	gStyle->SetPadBottomMargin(0.12);
	gStyle->SetPadLeftMargin(0.12);
	gStyle->SetPadRightMargin(0.02);
	gStyle->SetPadTopMargin(0.02);

	bool isSstate=(StateQuantumID[iState] > NRQCDvars::quID_S)?false:true;

	TCanvas *plotCanvas = new TCanvas("plotCanvas","plotCanvas",1200,800);
	plotCanvas->SetFillColor(kWhite);
	plotCanvas->SetFrameBorderMode(0);


	double x_min=0.;
	double x_max=160;
	double y_min;
	double y_max;
	if(iMeasurementID==0) { y_min=5e-2; y_max =1e2; }
	if(iMeasurementID==1) { y_min=-1.3; y_max =1.3; }
	if(iMeasurementID==2) { y_min=-0.7; y_max =0.7; }
	if(iMeasurementID==3) { y_min=-0.7; y_max =0.7; }


	TH1F *axishist = new TH1F;
	axishist = plotCanvas->DrawFrame(x_min,y_min,x_max,y_max);

	axishist->SetTitle(0);
	axishist->GetXaxis()->SetTitle("#it{p}_{T} [GeV]");
	axishist->GetYaxis()->SetTitle(Form("%s",NRQCDvars::MeasurementIDNameTex[iMeasurementID]));
	axishist->GetYaxis()->SetTitleSize(0.05);
	axishist->GetYaxis()->SetTitleOffset(1.);



	int colorData=1;
	int colorModel=600;
	double linewidthModel=2;
	double linewidthModelErrors=1;

	int colorCC[6]={632, 418, 616, 800, 0, 0};
	int colorFD[6]={632, 418, 616, 800, 0, 0};
	int linestyleCC=1;
	int linestyleFD=3;
	double linewidthCC=1;
	double linewidthFD=1;




	data_Graph->SetMarkerStyle(20);
	data_Graph->SetMarkerSize(1.5);
	data_Graph->SetMarkerColor(colorData);
	data_Graph->SetLineColor(colorData);

	model_Graph->SetLineColor(colorModel);
	model_Graph->SetLineStyle(1);
	model_Graph->SetLineWidth(linewidthModel);
	model_Graph->SetFillColor(colorModel);

	model_Graph_low->SetLineColor(colorModel);
	model_Graph_low->SetLineStyle(2);
	model_Graph_low->SetLineWidth(linewidthModelErrors);
	model_Graph_low->SetFillColor(colorModel);

	model_Graph_high->SetLineColor(colorModel);
	model_Graph_high->SetLineStyle(2);
	model_Graph_high->SetLineWidth(linewidthModelErrors);
	model_Graph_high->SetFillColor(colorModel);




	int nLegendEntries=2;
	if(plotDirectColorChannels) nLegendEntries+=ColorChannels_c;
	if(plotIndividualFeedDown) nLegendEntries+=StatesCont_c;
	if(plotInclusiveFeedDown) nLegendEntries+=1;





	//TLine* extreme0 = new TLine( pTmin, 0, pTmax, 0);
	//extreme0->SetLineWidth( 1 );
	//extreme0->SetLineStyle( 2 );
	//extreme0->SetLineColor( kBlack );
	//extreme0->Draw( "same" );

	//data_Graph->Draw("ap");
	model_Graph->Draw("Csame");
	model_Graph_low->Draw("Csame");
	model_Graph_high->Draw("Csame");
	data_Graph->Draw("psame");


	double delta_y_legend=0.05*nLegendEntries;
	double max_y_legend;
	double min_x_legend=0.65;
	double max_x_legend=0.95;
	if(iMeasurementID==0) max_y_legend=0.95;
	if(iMeasurementID>0 && iMeasurementID<4) max_y_legend=0.95;

	TLegend* legend;
	legend = new TLegend(min_x_legend,max_y_legend-delta_y_legend,max_x_legend,max_y_legend);
	legend->SetFillColor(0);
	legend->SetTextSize(0.03);
	legend->SetBorderSize(0);
	legend->AddEntry(data_Graph,Form("%s %s, %s", ExpName[iExperiment], StateNameTex[iState], rapchar),"lp");
	legend->AddEntry(model_Graph,"NRQCD inclusive model","l");


	if(plotDirectColorChannels){
		//cout<<"ColorChannels_c "<<ColorChannels_c<<endl;
		for(int i=0;i<ColorChannels_c;i++){
			//model_Graph_ColorChannels[i]->Print();
			model_Graph_ColorChannels[i]->SetLineColor(colorCC[i]);
			model_Graph_ColorChannels[i]->SetLineStyle(linestyleCC);
			model_Graph_ColorChannels[i]->SetLineWidth(linewidthCC);
			model_Graph_ColorChannels[i]->Draw("lsame");
			if(isSstate) legend->AddEntry(model_Graph_ColorChannels[i],Form("Direct production, O_{%s}^{%s}", ColorChannelNameTexS[i], StateNameTex[iState]),"l");
			else legend->AddEntry(model_Graph_ColorChannels[i],Form("Direct production, O_{%s}^{%s}", ColorChannelNameTexP[i], StateNameTex[iState]),"l");
		}
	}

	if(plotIndividualFeedDown){
		//cout<<"StatesCont_c "<<StatesCont_c<<endl;
		for(int i=1;i<StatesCont_c;i++){
			//cout<<"StateNameTex[StatesContributing[i]] "<<StateNameTex[StatesContributing[i]]<<endl;
			model_Graph_FeedDowns[i]->SetLineColor(colorFD[i]);
			model_Graph_FeedDowns[i]->SetLineStyle(linestyleFD);
			model_Graph_FeedDowns[i]->SetLineWidth(linewidthFD);
			model_Graph_FeedDowns[i]->Draw("lsame");
			legend->AddEntry(model_Graph_FeedDowns[i],Form("Feed-down from %s", StateNameTex[StatesContributing[i]]),"l");
		}
	}

	if(plotInclusiveFeedDown){
		//cout<<"draw inclusive feed-down"<<endl;
		// inclusive feed-down
		//cout<<"model_Graph_FeedDowns: "<<0<<endl;
		//model_Graph_FeedDowns[0]->Print();
		model_Graph_FeedDowns[0]->SetLineColor(colorFD[0]);
		model_Graph_FeedDowns[0]->SetLineStyle(linestyleFD);
		model_Graph_FeedDowns[0]->SetLineWidth(linewidthFD*2.);
		model_Graph_FeedDowns[0]->Draw("lsame");
		legend->AddEntry(model_Graph_FeedDowns[0],"Inclusive feed-down","l");
	}


	legend->Draw("same");

	////draw latex
	//double left=0.155, top=0.92, textSize=0.045;
	//TLatex *latex=new TLatex();
	//latex->SetTextFont(42);
	//latex->SetNDC(kTRUE);
	//latex->SetTextSize(textSize);
	//double stepLatex=textSize*1.3;
	//
	//latex->DrawLatex(left,top, "pp  #sqrt{s} = 7 TeV");

	if(iMeasurementID==0) plotCanvas->SetLogy(true);

	char savename[1000];
	sprintf(savename,"%s/Figures/DataModelComp_%s_%s_%s_rap%d.pdf",jobdirname, ExpName[iExperiment], StateName[iState], MeasurementIDName[iMeasurementID], iRap+1);
	cout<<"saved here: "<<savename<<endl;
	plotCanvas->SaveAs(savename);



}

vector<double> addPolarizations(vector<vector<double> > lamMatrix, vector<double> lamVecContributions){

	vector<double> lamSumVec (3,0); //return this

	int nContributions=lamVecContributions.size();

	double ContributionIntegral=0;
	for(int i=0; i<nContributions; i++){
		ContributionIntegral+=lamVecContributions[i];
	}
	if(ContributionIntegral==0) {cerr << "ContributionIntegral is 0 -> exit!"<<endl; exit(2);}
	vector<double> lamVecContributionFraction (nContributions,0);
	for(int i=0; i<nContributions; i++){
		lamVecContributionFraction[i]=lamVecContributions[i]/ContributionIntegral;

	}

	vector<double> LamthVec (nContributions);
	vector<double> LamphVec (nContributions);
	vector<double> LamtpVec (nContributions);

	for(int i=0; i<nContributions; i++){
		if(lamVecContributionFraction[i]>0){
			LamthVec[i]=lamMatrix[i][0];
			LamphVec[i]=lamMatrix[i][1];
			LamtpVec[i]=lamMatrix[i][2];

		}
	}

	double Lamth_numerator=0;
	double Lamth_denominator=0;
	for(int i=0; i<nContributions; i++){
		if(lamVecContributionFraction[i]>0){
			Lamth_numerator+=lamVecContributionFraction[i]*LamthVec[i]/(3+LamthVec[i]);
			Lamth_denominator+=lamVecContributionFraction[i]/(3+LamthVec[i]);
		}
	}

	double Lamph_numerator=0;
	for(int i=0; i<nContributions; i++){
		if(lamVecContributionFraction[i]>0) Lamph_numerator+=lamVecContributionFraction[i]*LamphVec[i]/(3+LamthVec[i]);
	}

	double Lamtp_numerator=0;
	for(int i=0; i<nContributions; i++){
		if(lamVecContributionFraction[i]>0) Lamtp_numerator+=lamVecContributionFraction[i]*LamtpVec[i]/(3+LamthVec[i]);
	}

	lamSumVec[0]=Lamth_numerator/Lamth_denominator;
	lamSumVec[1]=Lamph_numerator/Lamth_denominator;
	lamSumVec[2]=Lamtp_numerator/Lamth_denominator;

	//cout<<"lamSumVec[2] "<<lamSumVec[2]<<endl;
	//cout<<"Lamtp_numerator "<<Lamth_denominator<<endl;
	//cout<<"lamMatrix: "<<endl;
	//cout<<lamMatrix<<endl;
	//cout<<"lamVecContributionFraction: "<<endl;
	//cout<<lamVecContributionFraction<<endl;
	//cout<<"LamtpVec: "<<endl;
	//cout<<LamtpVec<<endl;



	return lamSumVec;
}

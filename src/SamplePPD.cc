/*
 * CalcPPD.cc
 *
 *  Created on: Jun 14, 2013
 *      Author: valentinknuenz
 */
#include "../interface/NRQCDglobalfitObject.h"
#include "../src/NRQCDglobalfitObject.cc"
#include "../interface/GlobalVar.h"
#include "../interface/MinuitLikelihoodFunction.h"

#include <iostream>
#include <sstream>
#include <iomanip>
#include <iterator>
#include <vector>
#include <string>
#include <cmath>

//rootincludes
#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TRandom3.h"
#include "TSystem.h"
//#include "TMatrixD.h"
#include "TFitterMinuit.h"
#include "Minuit2/FCNBase.h"

#include "TPad.h"
#include "TLine.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TH2.h"


using namespace NRQCDvars;

void setKernel( const int& kernelType );
void setFractionDimension( const int& nDim );
void genFractionValues1D(dvector &fraction, dvector &candidate, dvector &width );
void genFractionValues2D(dvector &fraction, dvector &candidate, dvector &width );
void genFractionValues3D(dvector &fraction, dvector &candidate, dvector &width );
void genFractionValues4D(dvector &fraction, dvector &candidate, dvector &width );
void genFractionValues5D(dvector &fraction, dvector &candidate, dvector &width );
void genFractionValues6D(dvector &fraction, dvector &candidate, dvector &width );
double MetropolisKernel(const double& candidate, const double& proposalWidth );
double MetropolisHastingsKernel(const double& candidate, const double& proposalWidth );
void transformFractionsToOps(dmatrix &Op, dmatrix &Fractions, dmatrix consts_star);
void (*getFractionValues)(dvector &fraction, dvector &candidate, dvector &width );
double (*kernelFunction)(const double& par0, const double& par1 );
void FindMPV(TH1* PosteriorDist , double& MPV , double& MPVerrorLow, double& MPVerrorHigh, int MPValgo, int nSigma);


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
	int 	Minimizer;
  	std::string filename("data");
  	std::string properties;

  	int 	nBurnIn=-1;
  	int 	nSample=-1;
  	bool 	SampleNp=false;
  	bool 	SampleNp_consts_star=false;
  	int 	MPValgo=-1;
  	int 	nSigma=-1;

  	setKernel(NRQCDvars::MetropolisHastings ); // modified by Joao: Set environment for sampling of fractions
  	//setKernel(NRQCDvars::Metropolis ); // modified by Joao: Set environment for sampling of fractions
	//int Minimizer=NRQCDvars::Minuit;

  	for( int i=0;i < argc; ++i ) {
	    if(std::string(argv[i]).find("SampleNp=true") != std::string::npos) {SampleNp=true; cout<<"SampleNp"<<endl;}
	    if(std::string(argv[i]).find("SampleNp_consts_star=true") != std::string::npos) {SampleNp_consts_star=true; cout<<"SampleNp_consts_star"<<endl;}
		if(std::string(argv[i]).find("storagedir") != std::string::npos) {char* storagedirchar = argv[i]; char* storagedirchar2 = strtok (storagedirchar, "="); storagedir = storagedirchar2; cout<<"storagedir = "<<storagedir<<endl;}
	    if(std::string(argv[i]).find("nSigma") != std::string::npos) { char* nSigmachar = argv[i]; char* nSigmachar2 = strtok (nSigmachar, "n"); nSigma = atoi(nSigmachar2); cout<<"nSigma = "<<nSigma<<endl; }
	    if(std::string(argv[i]).find("MPValgo") != std::string::npos) { char* MPValgochar = argv[i]; char* MPValgochar2 = strtok (MPValgochar, "M"); MPValgo = atoi(MPValgochar2); cout<<"MPValgo = "<<MPValgo<<endl; }
	    if(std::string(argv[i]).find("nBurnIn") != std::string::npos) {
	    	char* nBurnInchar = argv[i];
	    	char* nBurnInchar2 = strtok (nBurnInchar, "n");
	    	nBurnIn = atoi(nBurnInchar2);
	    	cout<<"nBurnIn = "<<nBurnIn<<endl;
			properties.append("_nBurnIn");
			properties.append(nBurnInchar2);
	    }
	    if(std::string(argv[i]).find("nSample") != std::string::npos) {
	    	char* nSamplechar = argv[i];
	    	char* nSamplechar2 = strtok (nSamplechar, "n");
	    	nSample = atoi(nSamplechar2);
	    	cout<<"nSample = "<<nSample<<endl;
			properties.append("_nSample");
			properties.append(nSamplechar2);
	    }
		if(std::string(argv[i]).find("JobID") != std::string::npos) {
			char* JobIDchar = argv[i];
			char* JobIDchar2 = strtok (JobIDchar, "=");
			JobID = JobIDchar2;
			cout<< "JobID = " << JobID << endl;
			properties.append("_JobID");
			properties.append(JobIDchar2);
		}
		if(std::string(argv[i]).find("DataModelCombID") != std::string::npos) {
			char* DataModelCombIDchar = argv[i];
			char* DataModelCombIDchar2 = strtok (DataModelCombIDchar, "=");
			DataModelCombID = DataModelCombIDchar2;
			cout<<"DataModelCombID = "<<DataModelCombID<<endl;
			properties.append("_DataModelCombID");
			properties.append(DataModelCombID);
		}
		if(std::string(argv[i]).find("ModelID") != std::string::npos) {
			char* ModelIDchar = argv[i];
			char* ModelIDchar2 = strtok (ModelIDchar, "=");
			ModelID = ModelIDchar2;
			cout<<"ModelID = "<<ModelID<<endl;
			properties.append("_ModelID");
			properties.append(ModelID);
		}
	    if(std::string(argv[i]).find("pTMin") != std::string::npos) {
	    	char* pTMinchar = argv[i];
	    	char* pTMinchar2 = strtok (pTMinchar, "p");
	    	pTMin = atof(pTMinchar2);
	    	cout<<"pTMin = "<<pTMin<<endl;
			properties.append("_pTMin");
			properties.append(pTMinchar2);
	    }
	    if(std::string(argv[i]).find("pTMax") != std::string::npos) {
	    	char* pTMaxchar = argv[i];
	    	char* pTMaxchar2 = strtok (pTMaxchar, "p");
	    	pTMax = atof(pTMaxchar2);
	    	cout<<"pTMax = "<<pTMax<<endl;
			properties.append("_pTMax");
			properties.append(pTMaxchar2);
	    }
	    if(std::string(argv[i]).find("rapMin") != std::string::npos) {
	    	char* rapMinchar = argv[i];
	    	char* rapMinchar2 = strtok (rapMinchar, "r");
	    	rapMin = atof(rapMinchar2);
	    	cout<<"rapMin = "<<rapMin<<endl;
			properties.append("_rapMin");
			properties.append(rapMinchar2);
	    }
	    if(std::string(argv[i]).find("rapMax") != std::string::npos) {
	    	char* rapMaxchar = argv[i];
	    	char* rapMaxchar2 = strtok (rapMaxchar, "r");
	    	rapMax = atof(rapMaxchar2);
	    	cout<<"rapMax = "<<rapMax<<endl;
			properties.append("_rapMax");
			properties.append(rapMaxchar2);
	    }
	    if(std::string(argv[i]).find("useOnlyState") != std::string::npos) {
	    	char* useOnlyStatechar = argv[i];
	    	char* useOnlyStatechar2 = strtok (useOnlyStatechar, "u");
	    	useOnlyState = atoi(useOnlyStatechar2);
	    	cout<<"useOnlyState = "<<useOnlyState<<endl;
			properties.append("_useOnlyState");
			properties.append(useOnlyStatechar2);
	    }
	    if(std::string(argv[i]).find("Minimizer") != std::string::npos) {
	    	char* Minimizerchar = argv[i];
	    	char* Minimizerchar2 = strtok (Minimizerchar, "M");
	    	Minimizer = atoi(Minimizerchar2);
	    	cout<<"Minimizer = "<<Minimizer<<endl;
			properties.append("_Minimizer");
			properties.append(Minimizerchar2);
	    }
	    if(std::string(argv[i]).find("useSstatesOnly=true") != std::string::npos) {
	    	useSstatesOnly=true;
	    	cout<<"useSstatesOnly=true"<<endl;
			properties.append("_useSstatesOnly");
	    }
	    if(std::string(argv[i]).find("usePstatesOnly=true") != std::string::npos) {
	    	usePstatesOnly=true;
	    	cout<<"usePstatesOnly=true"<<endl;
			properties.append("_usePstatesOnly");
	    }
	    if(std::string(argv[i]).find("useCharmoniumOnly=true") != std::string::npos) {
	    	useCharmoniumOnly=true;
	    	cout<<"useCharmoniumOnly=true"<<endl;
			properties.append("_useCharmoniumOnly");
	    }
	    if(std::string(argv[i]).find("useBottomoniumOnly=true") != std::string::npos) {
	    	useBottomoniumOnly=true;
	    	cout<<"useBottomoniumOnly=true"<<endl;
			properties.append("_useBottomoniumOnly");
	    }
	}
  	filename.append(properties);
	filename.append(".dat");

	delete gRandom;
	gRandom = new TRandom3(NRQCDvars::randomSeed);  // better random generator

	if(NRQCDvars::debug) cout << "Start SamplePPD" << endl;

	cout << "Read in measurements ( data + model )" << endl;

	bool DataSelected=false; // conect this bool with all the selection criteria of the measurements

	vector< NRQCDglobalfitObject > DataModelObject;

	char inname[2000];
	char modeldirname[2000];
	char datamodeldirname[2000];
	sprintf(modeldirname,"%s/ModelID",storagedir);
	gSystem->mkdir(modeldirname);
	sprintf(modeldirname,"%s/%s",modeldirname,ModelID);
	gSystem->mkdir(modeldirname);
	sprintf(datamodeldirname,"%s/DataModelCombID",storagedir);
	gSystem->mkdir(datamodeldirname);
	sprintf(datamodeldirname,"%s/%s",datamodeldirname,DataModelCombID);
	gSystem->mkdir(datamodeldirname);



	bool DataFromExperimentExists[NRQCDvars::nExperiments];

	for(int iExperiment=0; iExperiment < NRQCDvars::nExperiments; iExperiment++){
		DataFromExperimentExists[iExperiment]=false;
	}

	for(int iState=0; iState < nStates; iState++){
		for(int iMeasurementID=0; iMeasurementID < NRQCDvars::nMeasurementIDs; iMeasurementID++){
			for(int iExperiment=0; iExperiment < NRQCDvars::nExperiments; iExperiment++){
				for(int iRap = 0; iRap < NRQCDvars::nMaxRapBins; iRap++){
				    for(int iP = 0; iP < NRQCDvars::nMaxPtBins; iP++){

				    	//if(iState!=0 || iMeasurementID!=0 || iExperiment!=0 || iRap!=0 || iP!=0) continue;

						sprintf(inname,"%s/ConvertedDataModel_%s_%s_%s_rap%d_pT%d.txt",datamodeldirname, StateName[iState],
								MeasurementIDName[iMeasurementID],  ExpName[iExperiment], iRap, iP);

						ifstream in;
						in.open(inname);

						if( in.is_open() ){//Measurement present -> calculate model components:: Modified by Joao: more correct from ios point of view
							if(NRQCDvars::debug){
								cout << "Read in iState=" << StateName[iState] << ", iMeasurementID=" << MeasurementIDName[iMeasurementID];
								cout << ", iExperiment=" << ExpName[iExperiment] << ", iRap=" << iRap << ", iP=" << iP << endl;
							}
							NRQCDglobalfitObject readDataModelObject;
							in >> readDataModelObject;

							if(NRQCDvars::debug) readDataModelObject.Dump(nStates, true, true);
							in.close();

							if(
								   readDataModelObject.getpTMin() >= pTMin
								&& readDataModelObject.getpTMax() <= pTMax
								&& readDataModelObject.getyMin() >= rapMin
								&& readDataModelObject.getyMax() <= rapMax
							) DataSelected=true;

							if(useSstatesOnly && StateQuantumID[readDataModelObject.getState()] != NRQCDvars::quID_S)
								DataSelected=false;
							if(usePstatesOnly && StateQuantumID[readDataModelObject.getState()] == NRQCDvars::quID_S)
								DataSelected=false;
							if(useCharmoniumOnly && readDataModelObject.getState() > 3)
								DataSelected=false;
							if(useBottomoniumOnly && readDataModelObject.getState() < 4)
								DataSelected=false;
							if(useOnlyState < nStates &&  readDataModelObject.getState() != useOnlyState)
								DataSelected=false;

							if(DataSelected){
								DataModelObject.push_back(readDataModelObject);
								DataFromExperimentExists[iExperiment]=true;

								//readDataModelObject.Dump(NRQCDvars::nStates, true, true);

							}

							DataSelected=false;
						}
				    }
				}
			}
		}
	}

	int nDataPoints=DataModelObject.size();
	cout << "Number of considered data points: " << nDataPoints << endl;





	dcube directProductionCube;
	dmatrix promptProductionMatrix;
	double polCorrFactor;





	//Nuisance parameters Np

	dmatrix Np_BR; //Branching ratios, [nDauthgers][nMothers]
	dmatrix Np_US; //Uncertainty scales, [0=Data, 1=Model][nScales]
	dmatrix Np_BR_PreviousCandidates;
	dmatrix Np_US_PreviousCandidates;
	dmatrix Np_BR_SampleWidths;
	dmatrix Np_US_SampleWidths;
	dmatrix Np_BR_Expectation;
	dmatrix Np_US_Expectation;
	dmatrix Np_BR_ExpectationUncertainty;
	dmatrix Np_US_ExpectationUncertainty;

	dmatrix errlow_Np_BR;
	dmatrix errhigh_Np_BR;
	dmatrix errlow_Np_US;
	dmatrix errhigh_Np_US;

	//const_star matrix varied by uncertainty (Nuisance parameter in the sampling)
	dmatrix consts_star_var(NRQCDvars::nStates);
	dmatrix consts_star_var_var(NRQCDvars::nStates);
	dmatrix consts_star_var_PreviousCandidates(NRQCDvars::nStates);
	dmatrix consts_star_var_SampleWidths(NRQCDvars::nStates);
	dmatrix consts_star_var_Expectation(NRQCDvars::nStates);
	dmatrix consts_star_var_ExpectationUncertainty(NRQCDvars::nStates);

	dmatrix errlow_consts_star_var;
	dmatrix errhigh_consts_star_var;

	vector<double> consts_star_var_S (NRQCDvars::nColorChannels_S,0);
	vector<double> consts_star_var_P (NRQCDvars::nColorChannels_P,0);

	cout<<"Initialize Np's"<<endl;// (starting values of MH chain)

	dvector Np_BR_0 (NRQCDvars::nStates,0.);
	for(int i=0; i < NRQCDvars::nStates; i++) Np_BR.push_back(Np_BR_0);
	for(int i=0; i < NRQCDvars::nStates; i++) Np_BR_SampleWidths.push_back(Np_BR_0);
	for(int i=0; i < NRQCDvars::nStates; i++) Np_BR_ExpectationUncertainty.push_back(Np_BR_0);

	bool isBRzero[NRQCDvars::nStates][NRQCDvars::nStates];

	for(int i=0; i < NRQCDvars::nStates; i++){
		for(int j=0; j < NRQCDvars::nStates; j++){
			if(NRQCDvars::FeedDownBranchingRatio[i][j]>0){
				isBRzero[i][j]=false;
				Np_BR[i][j]=NRQCDvars::FeedDownBranchingRatio[i][j]*0.01;
				Np_BR_SampleWidths[i][j]=NRQCDvars::errFeedDownBranchingRatio[i][j]*0.01*NRQCDvars::proposalWidthBurnIn_Np_relToUnc;
				Np_BR_ExpectationUncertainty[i][j]=NRQCDvars::errFeedDownBranchingRatio[i][j]*0.01;
			}
			else{
				isBRzero[i][j]=true;
				Np_BR[i][j]=0;
				Np_BR_SampleWidths[i][j]=0;
				Np_BR_ExpectationUncertainty[i][j]=0;
			}
		}
	}
	Np_BR_PreviousCandidates=Np_BR;
	Np_BR_Expectation=Np_BR_PreviousCandidates;

	dvector Np_US_0 (NRQCDvars::nDataSystematicScales, 0.);
	dvector Np_US_1 (NRQCDvars::nModelSystematicScales, 0.);
	Np_US.push_back(Np_US_0);
	Np_US.push_back(Np_US_1);
	Np_US_PreviousCandidates.push_back(Np_US_0);
	Np_US_PreviousCandidates.push_back(Np_US_1);
	dvector Np_US_0_ (NRQCDvars::nDataSystematicScales, 1.*NRQCDvars::proposalWidthBurnIn_Np_relToUnc);
	dvector Np_US_1_ (NRQCDvars::nModelSystematicScales, 1.*NRQCDvars::proposalWidthBurnIn_Np_relToUnc);
	Np_US_SampleWidths.push_back(Np_US_0_);
	Np_US_SampleWidths.push_back(Np_US_1_);
	dvector Np_US_0__ (NRQCDvars::nDataSystematicScales, 1.);
	dvector Np_US_1__ (NRQCDvars::nModelSystematicScales, 1.);
	Np_US_ExpectationUncertainty.push_back(Np_US_0__);
	Np_US_ExpectationUncertainty.push_back(Np_US_1__);
	Np_US_Expectation=Np_US_PreviousCandidates;

	errlow_Np_BR=Np_BR_ExpectationUncertainty;
	errhigh_Np_BR=Np_BR_ExpectationUncertainty;
	errlow_Np_US=Np_US_ExpectationUncertainty;
	errhigh_Np_US=Np_US_ExpectationUncertainty;



















	int acceptedSampling;
	double loglikelihood=0;
	dvector ObjectLikelihoodVec;



	int nSampledPoints=nBurnIn+nSample;
	int nStep = nSampledPoints/100;  // visualize progress of the parameter sampling
	int nStep_ = 1;

	if(nSampledPoints<10) nStep=nSampledPoints/1;


  	//General S and P vectors, can be used at any time and place
	vector<double> S_vector (NRQCDvars::nColorChannels_S,0);
	vector<double> P_vector (NRQCDvars::nColorChannels_P,0);

	//Observable parameters Op (Matrix elements)
	dmatrix Op(NRQCDvars::nStates);
	vector<double> Op_S (NRQCDvars::nColorChannels_S,0);
	vector<double> Op_P (NRQCDvars::nColorChannels_P,0);

	dmatrix Op_var(NRQCDvars::nStates);
	dmatrix Op_plus(NRQCDvars::nStates);
	dmatrix Op_minus(NRQCDvars::nStates);
	dmatrix err_Op(NRQCDvars::nStates);
	dmatrix errhigh_Op(NRQCDvars::nStates);
	dmatrix errlow_Op(NRQCDvars::nStates);

	dmatrix Fractions(NRQCDvars::nStates);
	vector<double> Fractions_S (NRQCDvars::nColorChannels_S,0);//f0...R, fi: i going from 1 to n=nColorChannels, fn=1-sum(fi_i-(n-1))
	vector<double> Fractions_P (NRQCDvars::nColorChannels_P,0);

	dmatrix Fractions_var(NRQCDvars::nStates);
	dmatrix Fractions_plus(NRQCDvars::nStates);
	dmatrix Fractions_minus(NRQCDvars::nStates);
	dmatrix err_Fractions(NRQCDvars::nStates);
	dmatrix errhigh_Fractions(NRQCDvars::nStates);
	dmatrix errlow_Fractions(NRQCDvars::nStates);

	dmatrix Candidates(NRQCDvars::nStates);
	vector<double> Candidates_S (NRQCDvars::nColorChannels_S,0);
	vector<double> Candidates_P (NRQCDvars::nColorChannels_P,0);

	dmatrix SampleWidths(NRQCDvars::nStates);
	vector<double> SampleWidths_S (NRQCDvars::nColorChannels_S,0);
	vector<double> SampleWidths_P (NRQCDvars::nColorChannels_P,0);

	dmatrix PreviousCandidates (NRQCDvars::nStates);
	vector<double> PreviousCandidates_S (NRQCDvars::nColorChannels_S,0);
	vector<double> PreviousCandidates_P (NRQCDvars::nColorChannels_P,0);

	dmatrix NullMatrix(NRQCDvars::nStates);
	vector<double> NullMatrix_S (NRQCDvars::nColorChannels_S,0);
	vector<double> NullMatrix_P (NRQCDvars::nColorChannels_P,0);

	dmatrix consts_star;
	dmatrix err_consts_star;



	sprintf(inname,"%s/ModelIngredients_consts_star.txt",modeldirname);
	cout<<"read in dmatrix consts_star from "<<inname<<endl;
	ifstream instar;
	instar.open(inname);
	instar >> consts_star;
	instar >> err_consts_star;
	instar.close();

	cout << consts_star << endl;
	cout << err_consts_star << endl;

	consts_star_var=consts_star;
	consts_star_var_PreviousCandidates=consts_star_var;
	consts_star_var_Expectation=consts_star_var_PreviousCandidates;

	cout<<"Initialize consts_star_var's"<<endl;// (starting values of MH chain set further below)

	for (int i=0; i<NRQCDvars::nStates; i++){
		bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
		if(isSstate){
			for (int j=0; j<NRQCDvars::nColorChannels_S; j++){
				consts_star_var_S.at(j)=err_consts_star[i][j]*NRQCDvars::proposalWidthBurnIn_Np_relToUnc;
			}
			consts_star_var_SampleWidths.at(i)=consts_star_var_S;
		}
		else{
			for (int j=0; j<NRQCDvars::nColorChannels_P; j++){
				consts_star_var_P.at(j)=err_consts_star[i][j]*NRQCDvars::proposalWidthBurnIn_Np_relToUnc;
			}
			consts_star_var_SampleWidths.at(i)=consts_star_var_P;
		}
	}

	for (int i=0; i<NRQCDvars::nStates; i++){
		bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
		if(isSstate){
			for (int j=0; j<NRQCDvars::nColorChannels_S; j++){
				consts_star_var_S.at(j)=err_consts_star[i][j];
			}
			consts_star_var_ExpectationUncertainty.at(i)=consts_star_var_S;
		}
		else{
			for (int j=0; j<NRQCDvars::nColorChannels_P; j++){
				consts_star_var_P.at(j)=err_consts_star[i][j];
			}
			consts_star_var_ExpectationUncertainty.at(i)=consts_star_var_P;
		}
	}

	errlow_consts_star_var=consts_star_var_ExpectationUncertainty;
	errhigh_consts_star_var=consts_star_var_ExpectationUncertainty;

	double loglikelihood_PreviousCandidate = -1.e30;  // intial (arbitrary) values
	double loglikelihood_Candidate = -1.e30;  // intial (arbitrary) values

	bool BurnIn=true;


	cout<<"set NullMatrix:"<<endl;

	for (int i=0; i<NRQCDvars::nStates; i++){
		bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
		if(isSstate){
			for (int j=0; j<NRQCDvars::nColorChannels_S; j++){
				NullMatrix_S.at(j)=0.;
			}
			NullMatrix.at(i)=NullMatrix_S;
		}
		else{
			for (int j=0; j<NRQCDvars::nColorChannels_P; j++){
				NullMatrix_P.at(j)=0.;
			}
			NullMatrix.at(i)=NullMatrix_P;
		}
	}


	cout<<"set starting point of Candidates:"<<endl;

	double RstartingVal=26.;

	for (int i=0; i<NRQCDvars::nStates; i++){
		bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
		if(isSstate){
			for (int j=0; j<NRQCDvars::nColorChannels_S; j++){
				if(j==0) Candidates_S.at(j)=RstartingVal;
				else Candidates_S.at(j)=1./double(NRQCDvars::nColorChannels_S-1);

				//if(j==1) Candidates_S.at(j)=0.;
				//if(j==2) Candidates_S.at(j)=1.;
				//if(j==3) Candidates_S.at(j)=0.;
			}
			Candidates.at(i)=Candidates_S;
		}
		else{
			for (int j=0; j<NRQCDvars::nColorChannels_P; j++){
				if(j==0) Candidates_P.at(j)=RstartingVal;
				else Candidates_P.at(j)=1./double(NRQCDvars::nColorChannels_P-1);
			}
			Candidates.at(i)=Candidates_P;
		}
	}
	PreviousCandidates=Candidates;
	Fractions=Candidates;

	cout<<Candidates<<endl;

	consts_star_var=consts_star;

	Op=NullMatrix;
	transformFractionsToOps(Op, Fractions, consts_star);



	cout<<"set widths"<<endl;

	for (int i=0; i<NRQCDvars::nStates; i++){
		bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
		if(isSstate){
			for (int j=0; j<NRQCDvars::nColorChannels_S; j++){
				if(j==0) SampleWidths_S[j]=NRQCDvars::proposalWidthBurnIn_R;
				else SampleWidths_S.at(j)=NRQCDvars::proposalWidthBurnIn_f;
			}
			SampleWidths.at(i)=SampleWidths_S;
		}
		else{
			for (int j=0; j<NRQCDvars::nColorChannels_P; j++){
				if(j==0) SampleWidths_P[j]=NRQCDvars::proposalWidthBurnIn_R;
				else SampleWidths_P.at(j)=NRQCDvars::proposalWidthBurnIn_f;
			}
			SampleWidths.at(i)=SampleWidths_P;
		}
	}




//////////////////////////////////////////////////////////////////////////////
/////// Check influence of Op and Nps on likelihood //////////////////////////
//////////////////////////////////////////////////////////////////////////////

	cout<<"Check influence of Op and Nps on likelihood"<<endl;
	cout<<"-->> Decide which Op's and Np's to fix in MH and Minuit"<<endl;

	double FreeParam_LikelihoodBaseline, FreeParam_LikelihoodVariation;


	imatrix FreeParam_Fractions(NRQCDvars::nStates);
	ivector S_ivector (NRQCDvars::nColorChannels_S,0);
	ivector P_ivector (NRQCDvars::nColorChannels_P,0);

	for (int i=0; i<NRQCDvars::nStates; i++){
		bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
		if(isSstate){
			FreeParam_Fractions.at(i)=S_ivector;
		}
		else{
			FreeParam_Fractions.at(i)=P_ivector;
		}
	}

	for (int i=0; i<NRQCDvars::nStates; i++){
		int nColorChannels_state;
		bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
		if(isSstate) nColorChannels_state=NRQCDvars::nColorChannels_S;
		else nColorChannels_state=NRQCDvars::nColorChannels_P;
		for (int j=0; j<nColorChannels_state; j++){

			transformFractionsToOps(Op, Fractions, consts_star);
			loglikelihood=0;
			for(vector< NRQCDglobalfitObject >::iterator state = DataModelObject.begin(); state != DataModelObject.end(); ++state){
				ObjectLikelihoodVec=state->getObjectLikelihood(Op, Np_BR, Np_US, false, directProductionCube, promptProductionMatrix, polCorrFactor);
				loglikelihood+=ObjectLikelihoodVec[0];
			}
			FreeParam_LikelihoodBaseline=loglikelihood;

			Fractions[i][j]=Fractions[i][j]+10;
			transformFractionsToOps(Op, Fractions, consts_star);
			loglikelihood=0;
			for(vector< NRQCDglobalfitObject >::iterator state = DataModelObject.begin(); state != DataModelObject.end(); ++state){
				ObjectLikelihoodVec=state->getObjectLikelihood(Op, Np_BR, Np_US, false, directProductionCube, promptProductionMatrix, polCorrFactor);
				loglikelihood+=ObjectLikelihoodVec[0];
			}
			Fractions[i][j]=Fractions[i][j]-10;
			FreeParam_LikelihoodVariation=loglikelihood;

			if(fabs(FreeParam_LikelihoodBaseline-FreeParam_LikelihoodVariation)>1e-20) FreeParam_Fractions[i][j]=1;
		}
	}


	cout<<"FreeParam_Fractions"<<endl;
	cout<<FreeParam_Fractions<<endl;

	ivector FreeParam_Fractions_States(NRQCDvars::nStates);
	for (int i=0; i<NRQCDvars::nStates; i++){
		int nColorChannels_state;
		bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
		if(isSstate) nColorChannels_state=NRQCDvars::nColorChannels_S;
		else nColorChannels_state=NRQCDvars::nColorChannels_P;
		bool FreeParam_Fractions_State=false;
		for (int j=0; j<nColorChannels_state; j++){
			if(FreeParam_Fractions[i][j]==1) FreeParam_Fractions_State=true;
		}
		if(FreeParam_Fractions_State) FreeParam_Fractions_States.at(i)=1;
	}


	imatrix FreeParam_Np_BR;
	ivector FreeParam_Np_BR_0 (NRQCDvars::nStates,0);
	for(int i=0; i < NRQCDvars::nStates; i++) FreeParam_Np_BR.push_back(FreeParam_Np_BR_0);
	double FreeParam_nSigmaVar=5;

	for (int i=0; i<NRQCDvars::nStates; i++){
		for (int j=0; j<NRQCDvars::nStates; j++){
			if(NRQCDvars::FeedDownBranchingRatio[i][j]>0){

				loglikelihood=0;
				for(vector< NRQCDglobalfitObject >::iterator state = DataModelObject.begin(); state != DataModelObject.end(); ++state){
					ObjectLikelihoodVec=state->getObjectLikelihood(Op, Np_BR, Np_US, false, directProductionCube, promptProductionMatrix, polCorrFactor);
					loglikelihood+=ObjectLikelihoodVec[0];
				}
				FreeParam_LikelihoodBaseline=loglikelihood;

				Np_BR[i][j]=NRQCDvars::FeedDownBranchingRatio[i][j]*0.01+FreeParam_nSigmaVar*NRQCDvars::errFeedDownBranchingRatio[i][j]*0.01;
				loglikelihood=0;
				for(vector< NRQCDglobalfitObject >::iterator state = DataModelObject.begin(); state != DataModelObject.end(); ++state){
					ObjectLikelihoodVec=state->getObjectLikelihood(Op, Np_BR, Np_US, false, directProductionCube, promptProductionMatrix, polCorrFactor);
					loglikelihood+=ObjectLikelihoodVec[0];
				}
				Np_BR[i][j]=NRQCDvars::FeedDownBranchingRatio[i][j]*0.01-FreeParam_nSigmaVar*NRQCDvars::errFeedDownBranchingRatio[i][j]*0.01;
				FreeParam_LikelihoodVariation=loglikelihood;

				cout<<"FreeParam_Np_BR "<<i<<" "<<j<<" FreeParam_LikelihoodBaseline "<<FreeParam_LikelihoodBaseline<<endl;
				cout<<"FreeParam_Np_BR "<<i<<" "<<j<<" FreeParam_LikelihoodVariation "<<FreeParam_LikelihoodVariation<<endl;

				if(fabs(FreeParam_LikelihoodBaseline-FreeParam_LikelihoodVariation)>1e-20) FreeParam_Np_BR[i][j]=1;

			}
		}
	}


	cout<<"FreeParam_Np_BR"<<endl;
	cout<<FreeParam_Np_BR<<endl;


	imatrix FreeParam_Np_US;
	ivector FreeParam_Np_US_0_ (NRQCDvars::nDataSystematicScales, 0);
	ivector FreeParam_Np_US_1_ (NRQCDvars::nModelSystematicScales, 0);
	FreeParam_Np_US.push_back(FreeParam_Np_US_0_);
	FreeParam_Np_US.push_back(FreeParam_Np_US_1_);
	FreeParam_nSigmaVar=5;

	for(int j=0; j < NRQCDvars::nDataSystematicScales; j++){
		loglikelihood=0;
		for(vector< NRQCDglobalfitObject >::iterator state = DataModelObject.begin(); state != DataModelObject.end(); ++state){
			ObjectLikelihoodVec=state->getObjectLikelihood(Op, Np_BR, Np_US, false, directProductionCube, promptProductionMatrix, polCorrFactor);
			loglikelihood+=ObjectLikelihoodVec[0];
		}
		FreeParam_LikelihoodBaseline=loglikelihood;

		Np_US[0][j]=Np_US[0][j]+Np_US_ExpectationUncertainty[0][j]*FreeParam_nSigmaVar;
		loglikelihood=0;
		for(vector< NRQCDglobalfitObject >::iterator state = DataModelObject.begin(); state != DataModelObject.end(); ++state){
			ObjectLikelihoodVec=state->getObjectLikelihood(Op, Np_BR, Np_US, false, directProductionCube, promptProductionMatrix, polCorrFactor);
			loglikelihood+=ObjectLikelihoodVec[0];
		}
		Np_US[0][j]=Np_US[0][j]-Np_US_ExpectationUncertainty[0][j]*FreeParam_nSigmaVar;
		FreeParam_LikelihoodVariation=loglikelihood;
		cout<<"FreeParam_Np_US 0"<<" "<<j<<" FreeParam_LikelihoodBaseline "<<FreeParam_LikelihoodBaseline<<endl;
		cout<<"FreeParam_Np_US 0"<<" "<<j<<" FreeParam_LikelihoodVariation "<<FreeParam_LikelihoodVariation<<endl;

		if(fabs(FreeParam_LikelihoodBaseline-FreeParam_LikelihoodVariation)>1e-20) FreeParam_Np_US[0][j]=1;
	}
	for(int j=0; j < NRQCDvars::nModelSystematicScales; j++){
		loglikelihood=0;
		for(vector< NRQCDglobalfitObject >::iterator state = DataModelObject.begin(); state != DataModelObject.end(); ++state){
			ObjectLikelihoodVec=state->getObjectLikelihood(Op, Np_BR, Np_US, false, directProductionCube, promptProductionMatrix, polCorrFactor);
			loglikelihood+=ObjectLikelihoodVec[0];
		}
		FreeParam_LikelihoodBaseline=loglikelihood;

		Np_US[1][j]=Np_US[1][j]+Np_US_ExpectationUncertainty[1][j]*FreeParam_nSigmaVar;
		loglikelihood=0;
		for(vector< NRQCDglobalfitObject >::iterator state = DataModelObject.begin(); state != DataModelObject.end(); ++state){
			ObjectLikelihoodVec=state->getObjectLikelihood(Op, Np_BR, Np_US, false, directProductionCube, promptProductionMatrix, polCorrFactor);
			loglikelihood+=ObjectLikelihoodVec[0];
		}
		Np_US[1][j]=Np_US[1][j]-Np_US_ExpectationUncertainty[1][j]*FreeParam_nSigmaVar;
		FreeParam_LikelihoodVariation=loglikelihood;
		cout<<"FreeParam_Np_US 0"<<" "<<j<<" FreeParam_LikelihoodBaseline "<<FreeParam_LikelihoodBaseline<<endl;
		cout<<"FreeParam_Np_US 0"<<" "<<j<<" FreeParam_LikelihoodVariation "<<FreeParam_LikelihoodVariation<<endl;

		if(fabs(FreeParam_LikelihoodBaseline-FreeParam_LikelihoodVariation)>1e-20) FreeParam_Np_US[1][j]=1;
	}


	cout<<"FreeParam_Np_US"<<endl;
	cout<<FreeParam_Np_US<<endl;


	imatrix FreeParam_consts_star(NRQCDvars::nStates);
	FreeParam_nSigmaVar=5;

	for (int i=0; i<NRQCDvars::nStates; i++){
		bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
		if(isSstate){
			FreeParam_consts_star.at(i)=S_ivector;
		}
		else{
			FreeParam_consts_star.at(i)=P_ivector;
		}
	}

	for (int i=0; i<NRQCDvars::nStates; i++){
		int nColorChannels_state;
		bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
		if(isSstate) nColorChannels_state=NRQCDvars::nColorChannels_S;
		else nColorChannels_state=NRQCDvars::nColorChannels_P;
		for (int j=0; j<nColorChannels_state; j++){

			transformFractionsToOps(Op, Fractions, consts_star);
			loglikelihood=0;
			for(vector< NRQCDglobalfitObject >::iterator state = DataModelObject.begin(); state != DataModelObject.end(); ++state){
				ObjectLikelihoodVec=state->getObjectLikelihood(Op, Np_BR, Np_US, false, directProductionCube, promptProductionMatrix, polCorrFactor);
				loglikelihood+=ObjectLikelihoodVec[0];
			}
			FreeParam_LikelihoodBaseline=loglikelihood;

			consts_star[i][j]=consts_star[i][j]+FreeParam_nSigmaVar*err_consts_star[i][j];
			transformFractionsToOps(Op, Fractions, consts_star);
			loglikelihood=0;
			for(vector< NRQCDglobalfitObject >::iterator state = DataModelObject.begin(); state != DataModelObject.end(); ++state){
				ObjectLikelihoodVec=state->getObjectLikelihood(Op, Np_BR, Np_US, false, directProductionCube, promptProductionMatrix, polCorrFactor);
				loglikelihood+=ObjectLikelihoodVec[0];
			}
			consts_star[i][j]=consts_star[i][j]-FreeParam_nSigmaVar*err_consts_star[i][j];
			FreeParam_LikelihoodVariation=loglikelihood;
			cout<<"FreeParam_consts_star "<<i<<" "<<j<<" FreeParam_LikelihoodBaseline "<<FreeParam_LikelihoodBaseline<<endl;
			cout<<"FreeParam_consts_star "<<i<<" "<<j<<" FreeParam_LikelihoodVariation "<<FreeParam_LikelihoodVariation<<endl;

			if(fabs(FreeParam_LikelihoodBaseline-FreeParam_LikelihoodVariation)>1e-20) FreeParam_consts_star[i][j]=1;
		}
	}


	cout<<"FreeParam_consts_star"<<endl;
	cout<<FreeParam_consts_star<<endl;


////////////////////////////////////////////////
/////// DEFINE OUTPUT //////////////////////////
////////////////////////////////////////////////

	char outname[200];
	char jobdirname[200];
	sprintf(jobdirname,"%s/JobID", storagedir);
	gSystem->mkdir(jobdirname);
	sprintf(jobdirname,"%s/%s",jobdirname,JobID);
	gSystem->mkdir(jobdirname);
	sprintf(outname,"%s/results.root",jobdirname);

	TFile *ResultsFile = new TFile(outname, "RECREATE");




	if(Minimizer==NRQCDvars::Minuit){//Use Minuit to minimize the likelihood

		   double amin, edm, errdef;
		   int nvpar, nparx;

		   double Test=2;
		   double TestWidth=0.1;

			dvector Np_expected;
			dvector Np_uncertainty;

			double dNp_expected;
			double dNp_uncertainty;

			if(SampleNp){

				for(int j=0; j < NRQCDvars::nDataSystematicScales; j++){
					dNp_expected=0.;
					dNp_uncertainty=1.;
					Np_expected.push_back(dNp_expected);
					Np_uncertainty.push_back(dNp_uncertainty);
				}
				for(int j=0; j < NRQCDvars::nModelSystematicScales; j++){
					dNp_expected=0.;
					dNp_uncertainty=1.;
					Np_expected.push_back(dNp_expected);
					Np_uncertainty.push_back(dNp_uncertainty);
				}
				for(int i=0; i < nStates; i++){
					for(int j=0; j < nStates; j++){
						if(NRQCDvars::FeedDownBranchingRatio[i][j]>0){
							dNp_expected=NRQCDvars::FeedDownBranchingRatio[i][j]*0.01;
							dNp_uncertainty=NRQCDvars::errFeedDownBranchingRatio[i][j]*0.01;
							Np_expected.push_back(dNp_expected);
							Np_uncertainty.push_back(dNp_uncertainty);
						}
					}
				}


			}

			if(SampleNp_consts_star){

				for (int i=0; i<NRQCDvars::nStates; i++){
					int nColorChannels_state;
					bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
					if(isSstate) nColorChannels_state=NRQCDvars::nColorChannels_S;
					else nColorChannels_state=NRQCDvars::nColorChannels_P;
					for (int j=0; j<nColorChannels_state; j++){
							dNp_expected=consts_star[i][j];
							dNp_uncertainty=err_consts_star[i][j];
							Np_expected.push_back(dNp_expected);
							Np_uncertainty.push_back(dNp_uncertainty);
					}
				}

			}


		   //void Gauss3DFunc(int& /*npar*/, double* /*gin*/, double& fval, double* par, int /*iflag*/){}
		   FcnPVLogL* fcnx = new FcnPVLogL(SampleNp, SampleNp_consts_star, DataModelObject, consts_star, err_consts_star, Np_expected, Np_uncertainty);

		   TFitterMinuit *minuitx = new TFitterMinuit();
		   minuitx->SetMinuitFCN(fcnx);

		   	int nFixedPars=0;
		   	double start_value, verr, vlow, vhigh;
			char iParName[200];
		   	int iPar=0;
			for (int i=0; i<NRQCDvars::nStates; i++){
				int nColorChannels_state;
				bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
				if(isSstate) nColorChannels_state=NRQCDvars::nColorChannels_S;
				else nColorChannels_state=NRQCDvars::nColorChannels_P;
				cout<<"state "<<i<<", nColorChannels_state "<<nColorChannels_state<<endl;
				for (int j=0; j<nColorChannels_state; j++){

					//define starting values
					if(j==0){
						start_value=26.;
						vlow=-1;
						vhigh=100;
						verr=1e-1;
					}
					else{
						start_value=1./double(nColorChannels_state-1);
						vlow=-5;
						vhigh=5;
						verr=1e-2;
					}


					sprintf(iParName,"FitPar_state%d_f%d",i,j);
					minuitx->SetParameter(iPar,iParName,start_value, verr, vlow, vhigh);
					cout<<"set parameter "<<iParName<<endl;

					bool fixStateFractions=false;
					if(j==nColorChannels_state-1)
						fixStateFractions=true;
					if(FreeParam_Fractions[i][j]!=1)
						fixStateFractions=true;

					if(fixStateFractions) {
						cout<<"...and fixed it"<<endl;
						   minuitx->FixParameter(iPar);
						   nFixedPars++;
					}

					iPar++;

				}
			}

			const int nPOI=iPar;

			cout<<"nPOI "<<nPOI<<endl;
			cout<<"nFixedPars "<<nFixedPars<<endl;

			//add nuisance parameters

		   	int nFixedPars_Nuis=0;
			char iParName_Nuis[200];
		   	int iPar_Nuis=0;
			bool fix_NuisPar=false;
			double nSigmaRegion_Nuis=5;
			double verr_Nuis=1e-3;

			if(SampleNp){

				for(int j=0; j < NRQCDvars::nDataSystematicScales; j++){
					dNp_expected=0.;
					dNp_uncertainty=1.;
					start_value=dNp_expected; vlow=dNp_expected-nSigmaRegion_Nuis*dNp_uncertainty; vhigh=dNp_expected+nSigmaRegion_Nuis*dNp_uncertainty; verr=verr_Nuis;
					sprintf(iParName_Nuis,"NuisPar_LuminosityScale_Experiment%d",j);
					minuitx->SetParameter(nPOI+iPar_Nuis,iParName_Nuis,start_value, verr, vlow, vhigh);
					cout<<"set nuisance parameter "<<iParName_Nuis<<endl;
					if(FreeParam_Np_US[0][j]!=1){
						cout<<"...and fixed it"<<endl;
						minuitx->FixParameter(nPOI+iPar_Nuis);
						nFixedPars_Nuis++;
					}
					iPar_Nuis++;
				}
				for(int j=0; j < NRQCDvars::nModelSystematicScales; j++){
					dNp_expected=0.;
					dNp_uncertainty=1.;
					start_value=dNp_expected; vlow=dNp_expected-nSigmaRegion_Nuis*dNp_uncertainty; vhigh=dNp_expected+nSigmaRegion_Nuis*dNp_uncertainty; verr=verr_Nuis;
					sprintf(iParName_Nuis,"NuisPar_ModelSystematicScale%d",j);
					minuitx->SetParameter(nPOI+iPar_Nuis,iParName_Nuis,start_value, verr, vlow, vhigh);
					cout<<"set nuisance parameter "<<iParName_Nuis<<endl;
					if(FreeParam_Np_US[1][j]!=1){
						cout<<"...and fixed it"<<endl;
						minuitx->FixParameter(nPOI+iPar_Nuis);
						nFixedPars_Nuis++;
					}
					iPar_Nuis++;
				}
				for(int i=0; i < nStates; i++){
					for(int j=0; j < nStates; j++){
						if(NRQCDvars::FeedDownBranchingRatio[i][j]>0){
							dNp_expected=NRQCDvars::FeedDownBranchingRatio[i][j]*0.01;
							dNp_uncertainty=NRQCDvars::errFeedDownBranchingRatio[i][j]*0.01;
							start_value=dNp_expected; vlow=dNp_expected-nSigmaRegion_Nuis*dNp_uncertainty; vhigh=dNp_expected+nSigmaRegion_Nuis*dNp_uncertainty; verr=verr_Nuis;
							sprintf(iParName_Nuis,"NuisPar_BR%dto%d",j,i);
							minuitx->SetParameter(nPOI+iPar_Nuis,iParName_Nuis,start_value, verr, vlow, vhigh);
							cout<<"set nuisance parameter "<<iParName_Nuis<<endl;
							if(FreeParam_Np_BR[i][j]!=1){
								cout<<"...and fixed it"<<endl;
								minuitx->FixParameter(nPOI+iPar_Nuis);
								nFixedPars_Nuis++;
							}
							iPar_Nuis++;
						}
					}
				}



			}


			if(SampleNp_consts_star){

				for (int i=0; i<NRQCDvars::nStates; i++){
					int nColorChannels_state;
					bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
					if(isSstate) nColorChannels_state=NRQCDvars::nColorChannels_S;
					else nColorChannels_state=NRQCDvars::nColorChannels_P;
					for (int j=0; j<nColorChannels_state; j++){
							dNp_expected=consts_star[i][j];
							dNp_uncertainty=err_consts_star[i][j];
							start_value=dNp_expected; vlow=dNp_expected-nSigmaRegion_Nuis*dNp_uncertainty; vhigh=dNp_expected+nSigmaRegion_Nuis*dNp_uncertainty; verr=verr_Nuis;
							sprintf(iParName_Nuis,"NuisPar_consts_star_state%d_CC%d",i,j);
							minuitx->SetParameter(nPOI+iPar_Nuis,iParName_Nuis,start_value, verr, vlow, vhigh);
							cout<<"set nuisance parameter "<<iParName_Nuis<<endl;
							if(FreeParam_consts_star[i][j]!=1){
								cout<<"...and fixed it"<<endl;
								minuitx->FixParameter(nPOI+iPar_Nuis);
								nFixedPars_Nuis++;
							}
							iPar_Nuis++;
					}
				}

			}


			const int nNUIS=iPar_Nuis;

			cout<<"nNUIS "<<nNUIS<<endl;
			cout<<"nFixedPars_Nuis "<<nFixedPars_Nuis<<endl;

			const int nPar = nPOI+nNUIS;




		   minuitx->SetPrintLevel(3);
		   minuitx->CreateMinimizer();
		   int nFits=10;
		   for(int iFits=0;iFits<nFits;iFits++){
			   cout<<"Fit number "<<iFits<<endl;
			   minuitx->Minimize();
			   minuitx->PrintResults(3,1.);
			   minuitx->GetStats(amin, edm, errdef, nvpar, nparx);
			   cout<<"GetStats"<<endl;
			   cout<<"amin = "<<amin<<endl;
			   cout<<"edm = "<<edm<<endl;
			   cout<<"errdef = "<<errdef<<endl;
			   cout<<"nvpar = "<<nvpar<<endl;
			   cout<<"nparx = "<<nparx<<endl;
			   if(edm>0 && edm<1e-8) break;
		   }


		   Double_t* arglist;

		   bool exeMINOS=false;
		   if(exeMINOS){
		   cout<<"STARTING MINOS"<<endl;
		   minuitx->ExecuteCommand("MINOS",arglist,nvpar);
		   minuitx->PrintResults(3,1.);
		   }



		   double resultPar[nPar];
		   double err_resultPar[nPar];

		   double eplus_resultPar[nPar];
		   double eminus_resultPar[nPar];
		   double eparab_resultPar[nPar];
		   double globcc_resultPar[nPar];

		  cout<<"MINUIT result:"<<endl;
		  for(int i=0;i<nPar;i++){
			   resultPar[i]=minuitx->GetParameter(i);
			   err_resultPar[i]=minuitx->GetParError(i);
			   minuitx->GetErrors(i, eplus_resultPar[i], eminus_resultPar[i], eparab_resultPar[i], globcc_resultPar[i]);
			   cout<<minuitx->GetParName(i)<<": "<<resultPar[i]<<" +- "<<err_resultPar[i]<<endl;
			   if(!exeMINOS){
				   eplus_resultPar[i]=err_resultPar[i];
				   eminus_resultPar[i]=err_resultPar[i];
			   }
		  }




		  cout<<"Write result to file:"<<endl;
			//Normalize fractions such that they sum up to 1
			double sum_fi[NRQCDvars::nStates];
			iPar=0;
			for (int i=0; i<NRQCDvars::nStates; i++){
				int nColorChannels_state;
				bool isSstate=(NRQCDvars::StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
				if(isSstate) nColorChannels_state=NRQCDvars::nColorChannels_S;
				else nColorChannels_state=NRQCDvars::nColorChannels_P;
				//cout<<"state "<<i<<", nColorChannels_state "<<nColorChannels_state<<endl;
				sum_fi[i]=0;
				iPar++;
				for (int j=1; j<nColorChannels_state-1; j++){
					sum_fi[i]+=resultPar[iPar];
					iPar++;
				}
				iPar++;
				cout<<"sum_fi["<<i<<"]="<<sum_fi[i]<<endl;
			}


			iPar=0;
			for (int i=0; i<NRQCDvars::nStates; i++){
				bool isSstate=(NRQCDvars::StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
				if(isSstate){
					for (int j=0; j<NRQCDvars::nColorChannels_S; j++){
						if(j<NRQCDvars::nColorChannels_S-1) Fractions_S.at(j)=resultPar[iPar];
						else Fractions_S.at(j)=1-sum_fi[i];
						iPar++;
					}
					Fractions.at(i)=Fractions_S;
				}
				else{
					for (int j=0; j<NRQCDvars::nColorChannels_P; j++){
						if(j<NRQCDvars::nColorChannels_P-1) Fractions_P.at(j)=resultPar[iPar];
						else Fractions_P.at(j)=1-sum_fi[i];
						iPar++;
					}
					Fractions.at(i)=Fractions_P;
				}
			}

			iPar=0;
			for (int i=0; i<NRQCDvars::nStates; i++){
				bool isSstate=(NRQCDvars::StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
				if(isSstate){
					for (int j=0; j<NRQCDvars::nColorChannels_S; j++){
						if(j<NRQCDvars::nColorChannels_S-1) S_vector.at(j)=err_resultPar[iPar];
						else S_vector.at(j)=0;
						iPar++;
					}
					err_Fractions.at(i)=S_vector;
				}
				else{
					for (int j=0; j<NRQCDvars::nColorChannels_P; j++){
						if(j<NRQCDvars::nColorChannels_P-1) P_vector.at(j)=err_resultPar[iPar];
						else P_vector.at(j)=0;
						iPar++;
					}
					err_Fractions.at(i)=P_vector;
				}
			}

			iPar=0;
			for (int i=0; i<NRQCDvars::nStates; i++){
				bool isSstate=(NRQCDvars::StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
				if(isSstate){
					for (int j=0; j<NRQCDvars::nColorChannels_S; j++){
						if(j<NRQCDvars::nColorChannels_S-1) S_vector.at(j)=eplus_resultPar[iPar];
						else S_vector.at(j)=0;
						iPar++;
					}
					errhigh_Fractions.at(i)=S_vector;
				}
				else{
					for (int j=0; j<NRQCDvars::nColorChannels_P; j++){
						if(j<NRQCDvars::nColorChannels_P-1) P_vector.at(j)=eplus_resultPar[iPar];
						else P_vector.at(j)=0;
						iPar++;
					}
					errhigh_Fractions.at(i)=P_vector;
				}
			}

			iPar=0;
			for (int i=0; i<NRQCDvars::nStates; i++){
				bool isSstate=(NRQCDvars::StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
				if(isSstate){
					for (int j=0; j<NRQCDvars::nColorChannels_S; j++){
						if(j<NRQCDvars::nColorChannels_S-1) S_vector.at(j)=fabs(eminus_resultPar[iPar]);
						else S_vector.at(j)=0;
						iPar++;
					}
					errlow_Fractions.at(i)=S_vector;
				}
				else{
					for (int j=0; j<NRQCDvars::nColorChannels_P; j++){
						if(j<NRQCDvars::nColorChannels_P-1) P_vector.at(j)=fabs(eminus_resultPar[iPar]);
						else P_vector.at(j)=0;
						iPar++;
					}
					errlow_Fractions.at(i)=P_vector;
				}
			}


			for (int i=0; i<NRQCDvars::nStates; i++){
				bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
				if(isSstate){
					for (int j=0; j<NRQCDvars::nColorChannels_S; j++){
						S_vector.at(j)=Fractions[i][j]+errhigh_Fractions[i][j];
					}
					Fractions_plus.at(i)=S_vector;
				}
				else{
					for (int j=0; j<NRQCDvars::nColorChannels_P; j++){
						P_vector.at(j)=Fractions[i][j]+errhigh_Fractions[i][j];
					}
					Fractions_plus.at(i)=P_vector;
				}
			}

			for (int i=0; i<NRQCDvars::nStates; i++){
				bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
				if(isSstate){
					for (int j=0; j<NRQCDvars::nColorChannels_S; j++){
						S_vector.at(j)=Fractions[i][j]-errlow_Fractions[i][j];
					}
					Fractions_minus.at(i)=S_vector;
				}
				else{
					for (int j=0; j<NRQCDvars::nColorChannels_P; j++){
						P_vector.at(j)=Fractions[i][j]-errlow_Fractions[i][j];
					}
					Fractions_minus.at(i)=P_vector;
				}
			}





			for (int i=0; i<NRQCDvars::nStates; i++){
				bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
				if(isSstate){
					for (int j=0; j<NRQCDvars::nColorChannels_S; j++){
						S_vector.at(j)=fabs(Op[i][j]-Op_minus[i][j]);
					}
					errlow_Op.at(i)=S_vector;
				}
				else{
					for (int j=0; j<NRQCDvars::nColorChannels_P; j++){
						P_vector.at(j)=fabs(Op[i][j]-Op_minus[i][j]);
					}
					errlow_Op.at(i)=P_vector;
				}
			}

			for (int i=0; i<NRQCDvars::nStates; i++){
				bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
				if(isSstate){
					for (int j=0; j<NRQCDvars::nColorChannels_S; j++){
						S_vector.at(j)=fabs(Op_plus[i][j]-Op[i][j]);
					}
					errhigh_Op.at(i)=S_vector;
				}
				else{
					for (int j=0; j<NRQCDvars::nColorChannels_P; j++){
						P_vector.at(j)=fabs(Op_plus[i][j]-Op[i][j]);
					}
					errhigh_Op.at(i)=P_vector;
				}
			}



			cout<<"FitResult for Fractions:"<<endl;
			cout<<Fractions<<endl;
			cout<<"FitResult for errFractions:"<<endl;
			cout<<err_Fractions<<endl;


			cout<<"FitResult for Op:"<<endl;
			cout<<Op<<endl;









			if(SampleNp){

				iPar=nPOI;
				cout<<"FitResults NP_US"<<endl;
				for(int j=0; j < NRQCDvars::nDataSystematicScales; j++){
					if(FreeParam_Np_US[0][j]==1){
						Np_US[0][j]=resultPar[iPar];
						errlow_Np_US[0][j]=err_resultPar[iPar];
						errhigh_Np_US[0][j]=err_resultPar[iPar];
					}
					else Np_US[0][j]=Np_US_Expectation[0][j];
					cout<<"iPar "<<iPar<<endl;
					iPar++;
				}
				for(int j=0; j < NRQCDvars::nModelSystematicScales; j++){
					if(FreeParam_Np_US[1][j]==1){
						Np_US[1][j]=resultPar[iPar];
						errlow_Np_US[1][j]=err_resultPar[iPar];
						errhigh_Np_US[1][j]=err_resultPar[iPar];
					}
					else Np_US[1][j]=Np_US_Expectation[1][j];
					cout<<"iPar "<<iPar<<endl;
					iPar++;
				}
				cout<<"FitResults NP_BR"<<endl;
				for(int i=0; i < nStates; i++){
					for(int j=0; j < nStates; j++){
						if(NRQCDvars::FeedDownBranchingRatio[i][j]>0){
							if(FreeParam_Np_BR[i][j]==1){
								Np_BR[i][j]=resultPar[iPar];
								errlow_Np_BR[i][j]=err_resultPar[iPar];
								errhigh_Np_BR[i][j]=err_resultPar[iPar];
							}
							else Np_BR[i][j]=Np_BR_Expectation[i][j];
							cout<<"iPar "<<iPar<<endl;
							iPar++;
						}
					}
				}
			}
			else{
				Np_US=Np_US_Expectation;
				Np_BR=Np_BR_Expectation;
			}


			if(SampleNp_consts_star){

				cout<<"FitResults NP_consts_star"<<endl;
				for (int i=0; i<NRQCDvars::nStates; i++){
					bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
					if(isSstate){
						for (int j=0; j<NRQCDvars::nColorChannels_S; j++){
							if(FreeParam_consts_star[i][j]==1){
								errlow_consts_star_var[i][j]=err_resultPar[iPar];
								errhigh_consts_star_var[i][j]=err_resultPar[iPar];
								consts_star_var[i][j]=resultPar[iPar];
							}
							else consts_star_var[i][j]=consts_star_var_Expectation[i][j];
							S_vector.at(j)=consts_star_var[i][j];
							cout<<"iPar "<<iPar<<endl;
							iPar++;
						}
						consts_star_var.at(i)=S_vector;
					}
					else{
						for (int j=0; j<NRQCDvars::nColorChannels_P; j++){
							if(FreeParam_consts_star[i][j]==1){
								consts_star_var[i][j]=resultPar[iPar];
								errlow_consts_star_var[i][j]=err_resultPar[iPar];
								errhigh_consts_star_var[i][j]=err_resultPar[iPar];
							}
							else consts_star_var[i][j]=consts_star_var_Expectation[i][j];
							P_vector.at(j)=consts_star_var[i][j];
							cout<<"iPar "<<iPar<<endl;
							iPar++;
						}
						consts_star_var.at(i)=P_vector;
					}
				}

			}
			else{
				consts_star_var=consts_star_var_Expectation;
			}





			//TODO: separately for each state, vary fi within uncertainties, ensure sum=1, calc Oi hist and define errlow_Oi and errhigh_Oi

			transformFractionsToOps(Op, Fractions, consts_star);
			transformFractionsToOps(Op_plus, Fractions_plus, consts_star);
			transformFractionsToOps(Op_minus, Fractions_minus, consts_star);

			TH1D *h_Oi[NRQCDvars::nStates][NRQCDvars::nColorChannels];
			int h_Oi_nBins=100;
			char hist_name[200];
			double nSigmaOp=1;

			for (int i=0; i<NRQCDvars::nStates; i++){
				int nColorChannels_state;
				bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
				if(isSstate) nColorChannels_state=NRQCDvars::nColorChannels_S;
				else nColorChannels_state=NRQCDvars::nColorChannels_P;
				for (int j=0; j<nColorChannels_state; j++){

					sprintf(hist_name,"h_Oi_state%d_CC%d",i,j);
					h_Oi[i][j]= new TH1D( hist_name, hist_name, h_Oi_nBins, Op[i][j]-nSigmaOp*errlow_Op[i][j], Op[i][j]+nSigmaOp*errhigh_Op[i][j]);
					//if(i==3&&j==1) {
					//	cout<<"Op[i][j] "<<Op[i][j]<<endl;
					//	cout<<"Op_plus[i][j] "<<Op_plus[i][j]<<endl;
					//	cout<<"Op_minus[i][j] "<<Op_minus[i][j]<<endl;
					//}
				}
			}



			int nRand=10000;
			for(int iRand=0;iRand<nRand;iRand++){



				Fractions_var=Fractions;
				consts_star_var_var=consts_star_var;
				Op_var=Op;

				for (int i=0; i<NRQCDvars::nStates; i++){
					int nColorChannels_state;
					bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
					if(isSstate) nColorChannels_state=NRQCDvars::nColorChannels_S;
					else nColorChannels_state=NRQCDvars::nColorChannels_P;
					double buffSum=0;
					for (int j=0; j<nColorChannels_state-1; j++){
						if(FreeParam_Fractions[i][j]!=1) continue;
						Fractions_var[i][j]=gRandom->Gaus(Fractions[i][j], err_Fractions[i][j]);
						if(j>0) buffSum+=Fractions_var[i][j];
					}
					Fractions_var[i][nColorChannels_state-1]=1-buffSum;
					for (int j=0; j<nColorChannels_state; j++){
						if(FreeParam_Fractions[i][j]!=1) continue;
						consts_star_var_var[i][j]=gRandom->Gaus(consts_star_var[i][j], (errlow_consts_star_var[i][j]+errhigh_consts_star_var[i][j])/2.);
					}
				}

				transformFractionsToOps(Op_var, Fractions_var, consts_star_var_var);


				for (int i=0; i<NRQCDvars::nStates; i++){
					int nColorChannels_state;
					bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
					if(isSstate) nColorChannels_state=NRQCDvars::nColorChannels_S;
					else nColorChannels_state=NRQCDvars::nColorChannels_P;
					for (int j=0; j<nColorChannels_state; j++){
						if(FreeParam_Fractions[i][j]!=1) continue;
						h_Oi[i][j]->Fill(Op_var[i][j]);
						//if(i==3&&j==1) cout<<Op_var[i][j]<<endl;

					}
				}

				//cout<<"Fractions_var"<<endl;
				//cout<<Fractions_var<<endl;
				//cout<<"consts_star_var_var"<<endl;
				//cout<<consts_star_var_var<<endl;
				//cout<<"Op_var"<<endl;
				//cout<<Op_var<<endl;



			}


			double buff_MPV[NRQCDvars::nStates][NRQCDvars::nColorChannels], buff_errlow[NRQCDvars::nStates][NRQCDvars::nColorChannels], buff_errhigh[NRQCDvars::nStates][NRQCDvars::nColorChannels];
			for (int i=0; i<NRQCDvars::nStates; i++){
				int nColorChannels_state;
				bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
				if(isSstate) nColorChannels_state=NRQCDvars::nColorChannels_S;
				else nColorChannels_state=NRQCDvars::nColorChannels_P;
				for (int j=0; j<nColorChannels_state; j++){

					if(FreeParam_Fractions[i][j]!=1){
						buff_errlow[i][j]=0;
						buff_errhigh[i][j]=0;
						continue;
					}

					if(i==3&&j==2) h_Oi[i][j]->SaveAs("tmpHist.root");
					cout<<i<<" "<<j<<endl;
					cout<<h_Oi[i][j]->GetMean()<<endl;
					cout<<h_Oi[i][j]->GetRMS()<<endl;
					FindMPV(h_Oi[i][j], buff_MPV[i][j], buff_errlow[i][j], buff_errhigh[i][j], MPValgo, nSigma);
					cout<<buff_MPV[i][j]<<endl;
					cout<<buff_errlow[i][j]<<endl;
					cout<<buff_errhigh[i][j]<<endl;
					buff_errlow[i][j]+=Op[i][j]-buff_MPV[i][j];
					buff_errhigh[i][j]-=Op[i][j]-buff_MPV[i][j];
					buff_errlow[i][j]=fabs(buff_errlow[i][j]);
					buff_errhigh[i][j]=fabs(buff_errhigh[i][j]);
					cout<<buff_errlow[i][j]<<endl;
					cout<<buff_errhigh[i][j]<<endl;

				}
			}


			for (int i=0; i<NRQCDvars::nStates; i++){
				bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
				if(isSstate){
					for (int j=0; j<NRQCDvars::nColorChannels_S; j++){
						S_vector.at(j)=buff_errlow[i][j];
					}
					errlow_Op.at(i)=S_vector;
				}
				else{
					for (int j=0; j<NRQCDvars::nColorChannels_P; j++){
						P_vector.at(j)=buff_errlow[i][j];
					}
					errlow_Op.at(i)=P_vector;
				}
			}
			for (int i=0; i<NRQCDvars::nStates; i++){
				bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
				if(isSstate){
					for (int j=0; j<NRQCDvars::nColorChannels_S; j++){
						S_vector.at(j)=buff_errhigh[i][j];
					}
					errhigh_Op.at(i)=S_vector;
				}
				else{
					for (int j=0; j<NRQCDvars::nColorChannels_P; j++){
						P_vector.at(j)=buff_errhigh[i][j];
					}
					errhigh_Op.at(i)=P_vector;
				}
			}




			sprintf(outname,"%s/results.txt",jobdirname);
			cout<<"save results to "<<outname<<endl;

		    ofstream out;
		    out.open(outname);//, std::ofstream::app);

			out << Fractions;
			out << errlow_Fractions;
			out << errhigh_Fractions;
			out << Op;
			out << errlow_Op;
			out << errhigh_Op;

		    out.close();


		    dmatrix Np_US_Altered;
		    dmatrix errlow_Np_US_Altered;
		    dmatrix errhigh_Np_US_Altered;
		    Np_US_Altered.push_back(Np_US[0]);
		    errlow_Np_US_Altered.push_back(errlow_Np_US[0]);
		    errhigh_Np_US_Altered.push_back(errhigh_Np_US[0]);
		  	if(NRQCDvars::nModelSystematicScales!=0){
			    Np_US_Altered.push_back(Np_US[1]);
		  		errlow_Np_US_Altered.push_back(errlow_Np_US[1]);
		  		errhigh_Np_US_Altered.push_back(errhigh_Np_US[1]);
		  	}


			sprintf(outname,"%s/results_Np.txt",jobdirname);
			cout<<"save results to "<<outname<<endl;

		    out.open(outname);//, std::ofstream::app);

			out << Np_BR;
			out << errlow_Np_BR;
			out << errhigh_Np_BR;
			out << consts_star_var;
			out << errlow_consts_star_var;
			out << errhigh_consts_star_var;
			out << Np_US_Altered;
			out << errlow_Np_US_Altered;
			out << errhigh_Np_US_Altered;

		    out.close();



			bool DrawLikelihood=true;

			if(DrawLikelihood){

				gSystem->mkdir(Form("%s/Figures",jobdirname));

				int nParToPlot;

				double iScanStart=1.;
				double iScanDelta=0.5;

				double jScanStart=0.5;
				double jScanDelta=0.25;

				bool Plot1D=true;
				bool Plot2D=false;

				double yTitleOffset=1.5;

				if(Plot1D){


					const int nPlots=3;
					int nParToPlotArray[nPlots]={8,9,10};

					for(int k=0;k<nPlots;k++){

						nParToPlot=nParToPlotArray[k];
						cout<<"Plotting 1D likelihood of parameter "<<minuitx->GetParName(nParToPlot)<<endl;
						dvector myPar(nPar,0);

						char xTitle[200];
						if(k==0) sprintf(xTitle,"R^{#psi(2S)}");
						if(k==1) sprintf(xTitle,"f_{%s}^{#psi(2S)}",NRQCDvars::ColorChannelNameTexS[1]);
						if(k==2) sprintf(xTitle,"f_{%s}^{#psi(2S)}",NRQCDvars::ColorChannelNameTexS[2]);

						//HowToFind the parameter to vary
						int iParChoice;
						int jParChoice;
						iParChoice=3;
						if(k==0) jParChoice=0;
						if(k==1) jParChoice=1;
						if(k==2) jParChoice=2;

						char drawchar[200];
						sprintf(drawchar,"e");

						int nSigma=3;

						iScanStart=resultPar[nParToPlot];
						iScanDelta=nSigma*err_resultPar[nParToPlot];

						int chi2CheckBins1D=50;
						TH1D* chi2_par0_1D   = new TH1D( "chi2_par0_1D", "chi2_par0_1D", chi2CheckBins1D,  iScanStart-iScanDelta, iScanStart+iScanDelta);
						int iLikeScan=0;
						for(int iScan=1;iScan<chi2CheckBins1D+1;iScan++){
							iLikeScan++;
							cout<<iLikeScan<<" scans out of "<<chi2CheckBins1D<<endl;


							iPar=0;
							for (int i=0; i<NRQCDvars::nStates; i++){
								int nColorChannels_state;
								bool isSstate=(NRQCDvars::StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
								if(isSstate) nColorChannels_state=NRQCDvars::nColorChannels_S;
								else nColorChannels_state=NRQCDvars::nColorChannels_P;
								for (int j=0; j<nColorChannels_state; j++){
									myPar[iPar]=Fractions[i][j];
									if(i==iParChoice && j==jParChoice){
										myPar[iPar]=iScanStart-iScanDelta+2.*iScanDelta*iScan/double(chi2CheckBins1D);
									}
									iPar++;
								}
							}

						for (int i=0; i<nNUIS; i++){
							myPar[iPar+i]=resultPar[iPar+i];
						}


						double chi2 = (*fcnx)(myPar);
						cout<<"chi2 "<<chi2<<endl;
						if(chi2<10*amin){
							chi2_par0_1D->SetBinContent(iScan,chi2);
							chi2_par0_1D->SetBinError(iScan,1e-10);
						}
						else cout<<"chi2> "<<10*amin<<" -> entry not plotted"<<endl;

						}

						chi2_par0_1D->Print("all");

						cout<<minuitx->GetParName(nParToPlot)<<" done..."<<endl;

						chi2_par0_1D->SetStats(0);
						chi2_par0_1D->SetTitle(0);

						chi2_par0_1D->Print();

						TCanvas *chi2CheckCanvas2D = new TCanvas("chi2CheckCanvas2D","chi2CheckCanvas2D",1000,800);
						chi2CheckCanvas2D->SetFillColor(kWhite);
						chi2CheckCanvas2D->GetFrame()->SetFillColor(kWhite);
						chi2CheckCanvas2D->GetFrame()->SetBorderSize(0);
						chi2CheckCanvas2D->cd(1); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);chi2_par0_1D->GetYaxis()->SetTitle("#chi^{2}");	chi2_par0_1D->GetXaxis()->SetTitle(xTitle);	chi2_par0_1D->GetYaxis()->SetTitleOffset(yTitleOffset); chi2_par0_1D->SetMarkerStyle(20); chi2_par0_1D->Draw(drawchar);

						TLine* centralLine = new TLine( iScanStart, chi2_par0_1D->GetMinimum(), iScanStart ,chi2_par0_1D->GetMaximum());
						centralLine->SetLineWidth( 1 );
						centralLine->SetLineStyle( 1 );
						centralLine->SetLineColor( kGreen+2 );
						centralLine->Draw( "same" );

						TLine* minus1Sig = new TLine( iScanStart-err_resultPar[nParToPlot], chi2_par0_1D->GetMinimum(), iScanStart-err_resultPar[nParToPlot] ,chi2_par0_1D->GetMaximum());
						minus1Sig->SetLineWidth( 1 );
						minus1Sig->SetLineStyle( 2 );
						minus1Sig->SetLineColor( kRed );
						minus1Sig->Draw( "same" );

						TLine* plus1Sig = new TLine( iScanStart+err_resultPar[nParToPlot], chi2_par0_1D->GetMinimum(), iScanStart+err_resultPar[nParToPlot] ,chi2_par0_1D->GetMaximum());
						plus1Sig->SetLineWidth( 1 );
						plus1Sig->SetLineStyle( 2 );
						plus1Sig->SetLineColor( kRed );
						plus1Sig->Draw( "same" );

						sprintf(outname,"%s/Figures/chi2_1DforPar_%s.pdf",jobdirname,minuitx->GetParName(nParToPlot));
						chi2CheckCanvas2D->SaveAs(outname);
						delete chi2CheckCanvas2D;

					}

					/*
					nParToPlot=9;
					cout<<"Plotting 1D likelihood of parameter "<<minuitx->GetParName(nParToPlot)<<endl;
					dvector myPar(nPar,0);


					int nSigma=3;

					iScanStart=resultPar[nParToPlot];
					iScanDelta=nSigma*err_resultPar[nParToPlot];

					int chi2CheckBins1D=50;
					TH1D* chi2_par1_1D   = new TH1D( "chi2_par1_1D", "chi2_par1_1D", chi2CheckBins1D,  iScanStart-iScanDelta, iScanStart+iScanDelta);
					int iLikeScan=0;
					for(int iScan=1;iScan<chi2CheckBins1D+1;iScan++){
							iLikeScan++;
							cout<<iLikeScan<<" scans out of "<<chi2CheckBins1D<<endl;

						iPar=0;
						for (int i=0; i<NRQCDvars::nStates; i++){
							int nColorChannels_state;
							bool isSstate=(NRQCDvars::StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
							if(isSstate) nColorChannels_state=NRQCDvars::nColorChannels_S;
							else nColorChannels_state=NRQCDvars::nColorChannels_P;
							for (int j=0; j<nColorChannels_state; j++){
								myPar[iPar]=Fractions[i][j];
								if(i==3 && j==1){
									myPar[iPar]=iScanStart-iScanDelta+2.*iScanDelta*iScan/double(chi2CheckBins1D);
								}
								iPar++;

						}
						for (int i=0; i<nNUIS; i++){
							myPar[iPar+i]=resultPar[iPar+i];
						}


					double chi2 = (*fcnx)(myPar);
					chi2_par1_1D->SetBinContent(iScan,chi2);
					}}
					cout<<"chi2_par0par1_2D done..."<<endl;

					chi2_par1_1D->SetStats(0);
					chi2_par1_1D->SetTitle(0);

					TCanvas *chi2CheckCanvas2D = new TCanvas("chi2CheckCanvas2D","chi2CheckCanvas2D",1000,800);
					chi2CheckCanvas2D->SetFillColor(kWhite);
					chi2CheckCanvas2D->GetFrame()->SetFillColor(kWhite);
					chi2CheckCanvas2D->GetFrame()->SetBorderSize(0);
					chi2CheckCanvas2D->cd(1); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);chi2_par1_1D->GetYaxis()->SetTitle("likelihood");	chi2_par1_1D->GetXaxis()->SetTitle("f_{0}^{#psi(2S)}");	chi2_par1_1D->GetYaxis()->SetTitleOffset(yTitleOffset); chi2_par1_1D->Draw("c");

					TLine* centralLine = new TLine( iScanStart, chi2_par1_1D->GetMinimum(), iScanStart ,chi2_par1_1D->GetMaximum());
					centralLine->SetLineWidth( 1 );
					centralLine->SetLineStyle( 1 );
					centralLine->SetLineColor( kGreen+2 );
					centralLine->Draw( "same" );

					TLine* minus1Sig = new TLine( iScanStart-err_resultPar[nParToPlot], chi2_par1_1D->GetMinimum(), iScanStart-err_resultPar[nParToPlot] ,chi2_par1_1D->GetMaximum());
					minus1Sig->SetLineWidth( 1 );
					minus1Sig->SetLineStyle( 2 );
					minus1Sig->SetLineColor( kRed );
					minus1Sig->Draw( "same" );

					TLine* plus1Sig = new TLine( iScanStart+err_resultPar[nParToPlot], chi2_par1_1D->GetMinimum(), iScanStart+err_resultPar[nParToPlot] ,chi2_par1_1D->GetMaximum());
					plus1Sig->SetLineWidth( 1 );
					plus1Sig->SetLineStyle( 2 );
					plus1Sig->SetLineColor( kRed );
					plus1Sig->Draw( "same" );


					sprintf(outname,"%s/Figures/likelihood1D_par1.pdf",jobdirname);
					chi2CheckCanvas2D->SaveAs(outname);
					delete chi2CheckCanvas2D;

*/
				}



				if(Plot2D){
					cout<<"Plotting 2D likelihood"<<endl;
					dvector myPar(nPar,0);

					int chi2CheckBins2D=30;
					TH2D* chi2_par0par1_2D   = new TH2D( "chi2_par0par1_2D", "chi2_par0par1_2D", chi2CheckBins2D,  iScanStart-iScanDelta, iScanStart+iScanDelta, chi2CheckBins2D,  jScanStart-jScanDelta, jScanStart+jScanDelta);
					int iLikeScan=0;
					for(int iScan=1;iScan<chi2CheckBins2D+1;iScan++){
						for(int jScan=1;jScan<chi2CheckBins2D+1;jScan++){
							iLikeScan++;
							cout<<iLikeScan<<" scans out of "<<chi2CheckBins2D*chi2CheckBins2D<<endl;

						iPar=0;
						for (int i=0; i<NRQCDvars::nStates; i++){
							int nColorChannels_state;
							bool isSstate=(NRQCDvars::StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
							if(isSstate) nColorChannels_state=NRQCDvars::nColorChannels_S;
							else nColorChannels_state=NRQCDvars::nColorChannels_P;
							for (int j=0; j<nColorChannels_state; j++){
								myPar[iPar]=Fractions[i][j];
								if(i==3 && j==0){
									myPar[iPar]=iScanStart-iScanDelta+2.*iScanDelta*iScan/double(chi2CheckBins2D);
								}
								else if(i==3 && j==1){
									myPar[iPar]=jScanStart-jScanDelta+2.*jScanDelta*jScan/double(chi2CheckBins2D);
								}
								//cout<<"state "<<i<<" color "<<j<<" iPar "<<iPar<<" myPar "<<myPar[iPar]<<endl;

								iPar++;
							}
						}
						for (int i=0; i<nNUIS; i++){
							myPar[iPar+i]=resultPar[iPar+i];
						}


					double chi2 = (*fcnx)(myPar);
					/*if(chi2<10*amin)*/chi2_par0par1_2D->SetBinContent(iScan,jScan,chi2);
					}}
					cout<<"chi2_par0par1_2D done..."<<endl;

					chi2_par0par1_2D->SetStats(0);

					TCanvas *chi2CheckCanvas2D = new TCanvas("chi2CheckCanvas2D","chi2CheckCanvas2D",1000,800);
					chi2CheckCanvas2D->SetFillColor(kWhite);
					chi2CheckCanvas2D->GetFrame()->SetFillColor(kWhite);
					chi2CheckCanvas2D->GetFrame()->SetBorderSize(0);
					chi2CheckCanvas2D->cd(1); gPad->SetLeftMargin(0.2); gPad->SetFillColor(kWhite);chi2_par0par1_2D->GetYaxis()->SetTitle("f1");	chi2_par0par1_2D->GetXaxis()->SetTitle("f0");	chi2_par0par1_2D->Draw("colz");
					sprintf(outname,"%s/Fig.pdf",jobdirname);
					chi2CheckCanvas2D->SaveAs(outname);
					delete chi2CheckCanvas2D;
				}


			}




	}

	if(Minimizer==NRQCDvars::MH){//Use Metropolis-Hastings to minimize the likelihood

		TTree*  outputTreeAllSamplings = new TTree("AllSamplings","AllSamplings"); // tree filled in all samplings after burnin
		//TTree*  outputTreeAccSamplings = new TTree("AccSamplings","AccSamplings"); // tree filled only in accepted samplings after burnin


		char branch_name[200];
		char branch_address[200];

		// Make BRANCHES for outputTreeAccSamplings:::
		// Make BRANCHES for fi
		for (int i=0; i<NRQCDvars::nStates; i++){
			int nColorChannels_state;
			bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
			if(isSstate) nColorChannels_state=NRQCDvars::nColorChannels_S;
			else nColorChannels_state=NRQCDvars::nColorChannels_P;
			for (int j=0; j<nColorChannels_state; j++){
					sprintf(branch_name,"state%d_f%d",i,j);
					sprintf(branch_address,"%s/D",branch_name);
					outputTreeAllSamplings->Branch(branch_name,     &Candidates[i][j],     branch_address);
				}
			}

		// Make BRANCHES for Op
		for (int i=0; i<NRQCDvars::nStates; i++){
			int nColorChannels_state;
			bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
			if(isSstate) nColorChannels_state=NRQCDvars::nColorChannels_S;
			else nColorChannels_state=NRQCDvars::nColorChannels_P;
			for (int j=0; j<nColorChannels_state; j++){
					sprintf(branch_name,"state%d_Op%d",i,j);
					sprintf(branch_address,"%s/D",branch_name);
					outputTreeAllSamplings->Branch(branch_name,     &Op[i][j],     branch_address);
				}
			}


		int iAccSampling;
		sprintf(branch_name,"iAccSampling");
		sprintf(branch_address,"%s/I",branch_name);
		outputTreeAllSamplings->Branch(branch_name,     &iAccSampling,     branch_address);

		int nSampledPointsTotal;
		sprintf(branch_name,"nSampledPointsTotal");
		sprintf(branch_address,"%s/I",branch_name);
		outputTreeAllSamplings->Branch(branch_name,     &nSampledPointsTotal,     branch_address);
		int BurnInInt;
		sprintf(branch_name,"BurnInInt");
		sprintf(branch_address,"%s/I",branch_name);
		outputTreeAllSamplings->Branch(branch_name,     &BurnInInt,     branch_address);



		// Make BRANCHES for outputTreeAllSamplings
		outputTreeAllSamplings->Branch("loglikelihood", &loglikelihood, "loglikelihood/D");
		outputTreeAllSamplings->Branch("acceptedSampling", &acceptedSampling, "acceptedSampling/I");
		double MH_av_eff;
		outputTreeAllSamplings->Branch("MH_av_eff",     &MH_av_eff,     "MH_av_eff/D");

		for(int i=0; i < NRQCDvars::nStates; i++){
			for(int j=0; j < NRQCDvars::nStates; j++){
				if(NRQCDvars::FeedDownBranchingRatio[i][j]>0){
					sprintf(branch_name,"Np_BR_Daughter%d_Mother%d",i,j);
					sprintf(branch_address,"%s/D",branch_name);
					outputTreeAllSamplings->Branch(branch_name, &Np_BR[i][j], branch_address);
				}
			}
		}

		char Np_US_name[200];
		char Np_US_address[200];
		for(int j=0; j < NRQCDvars::nDataSystematicScales; j++){
			sprintf(branch_name,"Np_US_DataSystematicScale%d",j);
			sprintf(branch_address,"%s/D",branch_name);
			outputTreeAllSamplings->Branch(branch_name, &Np_US[0][j], branch_address);
		}
		for(int j=0; j < NRQCDvars::nModelSystematicScales; j++){
			sprintf(branch_name,"Np_US_ModelSystematicScale%d",j);
			sprintf(branch_address,"%s/D",branch_name);
			outputTreeAllSamplings->Branch(branch_name, &Np_US[1][j], branch_address);
		}


		// Make BRANCHES for consts_star_var
		for (int i=0; i<NRQCDvars::nStates; i++){
			int nColorChannels_state;
			bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
			if(isSstate) nColorChannels_state=NRQCDvars::nColorChannels_S;
			else nColorChannels_state=NRQCDvars::nColorChannels_P;
			for (int j=0; j<nColorChannels_state; j++){
					sprintf(branch_name,"state%d_const_star%d",i,j);
					sprintf(branch_address,"%s/D",branch_name);
					outputTreeAllSamplings->Branch(branch_name,     &consts_star_var[i][j],     branch_address);
				}
			}


////////////////////////////////////////////////
/////// END DEFINE OUTPUT //////////////////////
////////////////////////////////////////////////

		cout<<"Start sampling ("<<nBurnIn<<" + "<<nSample<<" points)"<<endl;
		cout<<"Progress:"<<endl;

		int nMH_rejectedInARow=0;
		int iTotalSinceLastStep=0;
		int iAcceptedSinceLastStep=0;

		BurnInInt=1;
		iAccSampling=0;
		nSampledPointsTotal=1;
		for(int iSampledPoint = 1; iSampledPoint <= nSampledPoints; nSampledPointsTotal++){ // Sampling loop

			iAccSampling=iSampledPoint;
			iTotalSinceLastStep++;

			//if(nSampledPointsTotal>2) break;

			if(iSampledPoint==nBurnIn+1 && BurnIn){

				BurnIn=false;
				cout<<"BurnIn period finished"<<endl;
				double fi_burnin_in[NRQCDvars::nStates][NRQCDvars::nColorChannels];
				char fi_name_in[200];
				TH1D *h_BurnIn[NRQCDvars::nStates][NRQCDvars::nColorChannels];
				int nBins_fi_burnin=100;

				TTree *copy_outputTreeAllSamplings=(TTree*)outputTreeAllSamplings->CopyTree("acceptedSampling<1000");

				for (int i=0; i<NRQCDvars::nStates; i++){
					int nColorChannels_state;
					bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
					if(isSstate) nColorChannels_state=NRQCDvars::nColorChannels_S;
					else nColorChannels_state=NRQCDvars::nColorChannels_P;
					for (int j=0; j<nColorChannels_state; j++){
						sprintf(fi_name_in,"state%d_f%d",i,j);
						copy_outputTreeAllSamplings->SetBranchAddress( fi_name_in,  &fi_burnin_in[i][j] );
						h_BurnIn[i][j] = new TH1D (fi_name_in, fi_name_in, nBins_fi_burnin, outputTreeAllSamplings->GetMinimum(fi_name_in), outputTreeAllSamplings->GetMaximum(fi_name_in));
					}
				}

				int n_events = int( copy_outputTreeAllSamplings->GetEntries() );
				int nBinsh_pT=50;
				int nBinsh_rap=50;

				// loop over  events in the burnin ntuple
				for ( int i_event = 1; i_event <= n_events; i_event++ ) {

					copy_outputTreeAllSamplings->GetEvent( i_event-1 );

					//fill histos
					for (int i=0; i<NRQCDvars::nStates; i++){
						int nColorChannels_state;
						bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
						if(isSstate) nColorChannels_state=NRQCDvars::nColorChannels_S;
						else nColorChannels_state=NRQCDvars::nColorChannels_P;
						for (int j=0; j<nColorChannels_state; j++){
							h_BurnIn[i][j] -> Fill(fi_burnin_in[i][j]);
						}
					}
				}

				bool changeProposalWidthsAfterBurnIn=true;

				for (int i=0; i<NRQCDvars::nStates; i++){
					bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
					if(isSstate){
						for (int j=0; j<NRQCDvars::nColorChannels_S; j++){
							const int nQu=2;
							double qu[nQu];
							double q[nQu];
							double coverage=0.683;
							q[0]=(1.-coverage)/2.;
							q[1]=1-q[0];
							h_BurnIn[i][j]->GetQuantiles(nQu,qu,q);
							double newWidth=(qu[1]-qu[0])/2.;
							if(changeProposalWidthsAfterBurnIn){
								if(j==0 && newWidth<NRQCDvars::proposalWidthBurnIn_R) SampleWidths_S[j]=newWidth;
								if(j>1 && newWidth<NRQCDvars::proposalWidthBurnIn_f) SampleWidths_S[j]=newWidth;
								SampleWidths_S[j]=newWidth/2.;
							}
							//SampleWidths_S[j]=h_BurnIn[i][j]->GetRMS();
							cout<<"Proposal width for state "<<i<<", ColorChannel "<<j<<" = "<<SampleWidths_S[j]<<endl;
						}
						SampleWidths.at(i)=SampleWidths_S;
					}
					else{
						for (int j=0; j<NRQCDvars::nColorChannels_P; j++){
							const int nQu=2;
							double qu[nQu];
							double q[nQu];
							double coverage=0.683;
							q[0]=(1.-coverage)/2.;
							q[1]=1-q[0];
							h_BurnIn[i][j]->GetQuantiles(nQu,qu,q);
							double newWidth=(qu[1]-qu[0])/2.;
							if(changeProposalWidthsAfterBurnIn){
								if(j==0 && newWidth<NRQCDvars::proposalWidthBurnIn_R) SampleWidths_P[j]=newWidth;
								if(j>1 && newWidth<NRQCDvars::proposalWidthBurnIn_f) SampleWidths_P[j]=newWidth;
								SampleWidths_P[j]=newWidth/2.;
							}

							//SampleWidths_P[j]=h_BurnIn[i][j]->GetRMS();
							cout<<"Proposal width for state "<<i<<", ColorChannel "<<j<<" = "<<SampleWidths_P[j]<<endl;
						}
						SampleWidths.at(i)=SampleWidths_P;
					}
				}
				//BurnInSamplings->Write();
				delete copy_outputTreeAllSamplings;

				BurnInInt=0;

			}





			// Fill Nuisance parameter matrices Np_BR and Np_US. Np_US[0]: Data-related, Np_US[1]: Model-related uncertainty scales
			//cout << "Fill Nuisance parameter matrices" <<endl;

			 //DataFromExperimentExists
			//candidate[i]=kernelFunction(candidate[i], width[i]);

			if(NRQCDvars::debug) cout<<"Sampling Nps"<<endl;

			for(int j=0; j < NRQCDvars::nDataSystematicScales; j++){
				if(SampleNp && FreeParam_Np_US[0][j]==1) Np_US[0][j]=kernelFunction(Np_US_PreviousCandidates[0][j], Np_US_SampleWidths[0][j]);
				else Np_US[0][j]=0.;
			}
			for(int j=0; j < NRQCDvars::nModelSystematicScales; j++){
				if(SampleNp && FreeParam_Np_US[1][j]==1) Np_US[1][j]=kernelFunction(Np_US_PreviousCandidates[1][j], Np_US_SampleWidths[1][j]);
				else Np_US[1][j]=0.;
			}

			for(int i=0; i < nStates; i++){
				for(int j=0; j < nStates; j++){
					if(NRQCDvars::FeedDownBranchingRatio[i][j]>0  && FreeParam_Np_BR[i][j]==1){
						if(SampleNp) Np_BR[i][j]=kernelFunction(Np_BR_PreviousCandidates[i][j], Np_BR_SampleWidths[i][j]);
						else Np_BR[i][j]=NRQCDvars::FeedDownBranchingRatio[i][j]*0.01;
					}
					else Np_BR[i][j]=0;
				}
			}

			for (int i=0; i<NRQCDvars::nStates; i++){
				bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
				if(isSstate){
					for (int j=0; j<NRQCDvars::nColorChannels_S; j++){
						if(SampleNp_consts_star && FreeParam_consts_star[i][j]==1) consts_star_var_S.at(j)=kernelFunction(consts_star_var_PreviousCandidates[i][j], consts_star_var_SampleWidths[i][j]);
						else consts_star_var_S.at(j)=consts_star[i][j];
					}
					consts_star_var.at(i)=consts_star_var_S;
				}
				else{
					for (int j=0; j<NRQCDvars::nColorChannels_P; j++){
						if(SampleNp_consts_star && FreeParam_consts_star[i][j]==1) consts_star_var_P.at(j)=kernelFunction(consts_star_var_PreviousCandidates[i][j], consts_star_var_SampleWidths[i][j]);
						else consts_star_var_P.at(j)=consts_star[i][j];
					}
					consts_star_var.at(i)=consts_star_var_P;
				}
			}



			if(NRQCDvars::debug) cout<<"Sampling Fractions"<<endl;

			///// Set Fractions -> calculate Matrix Elements Op

			if(NRQCDvars::debug) cout<< "fill Op" << endl;


			//cout << "getFractionValues" <<endl;

			for (int i=0; i<NRQCDvars::nStates; i++){

				//cout << "i" <<i<<endl;
				bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
				if(isSstate){
					//cout << "Sstate" <<endl;
					//cout << "setFractionDimension:getFractionValues" <<NRQCDvars::nColorChannels_S-1<<"D"<<endl;
					setFractionDimension(NRQCDvars::nColorChannels_S-1);
					PreviousCandidates_S=PreviousCandidates.at(i);
					SampleWidths_S=SampleWidths.at(i);
					if(FreeParam_Fractions_States[i]==1) getFractionValues(Fractions_S, PreviousCandidates_S, SampleWidths_S);
					else Fractions_S=NullMatrix_S;
					//cout << "sampled fractions:" <<endl;
					//double SumFractions=0;
					//for (int j=0; j<Fractions_S.size(); j++){
					//	cout << "Fractions_S "<< Fractions_S[j] <<endl;
					//	if(j>0) SumFractions+=Fractions_S[j];
					//}
					//cout << "Sum f[1]-f[n] "<< SumFractions <<endl;
					Fractions.at(i)=Fractions_S;
				}
				else{
					//cout << "Pstate" <<endl;
					//cout << "setFractionDimension:getFractionValues" <<NRQCDvars::nColorChannels_P-1<<"D"<<endl;
					setFractionDimension(NRQCDvars::nColorChannels_P-1);
					PreviousCandidates_P=PreviousCandidates.at(i);
					SampleWidths_P=SampleWidths.at(i);
					if(FreeParam_Fractions_States[i]==1) getFractionValues(Fractions_P, PreviousCandidates_P, SampleWidths_P);
					else Fractions_P=NullMatrix_P;
					//cout << "sampled fractions:" <<endl;
					//double SumFractions=0;
					//for (int j=0; j<Fractions_P.size(); j++){
					//	cout << "Fractions_P "<< Fractions_P[j] <<endl;
					//	if(j>0) SumFractions+=Fractions_P[j];
					//}
					//cout << "Sum f[1]-f[n] "<< SumFractions <<endl;
					Fractions.at(i)=Fractions_P;
				}
			}

			Candidates=Fractions;

			// Relate Op to R, fi -> getObjectLikelihood
			//cout << "transformFractionsToOps" <<endl;


			transformFractionsToOps(Op, Candidates, consts_star_var);


			//cout<<Op;
			//
			//cout<<"Op[0].size() "<<Op[0].size()<<endl;
			//cout<<"Op[1].size() "<<Op[1].size()<<endl;

			if(NRQCDvars::debug) cout << "getObjectLikelihood" << endl;

			if(NRQCDvars::debug) cout<<"getObjectLikelihood"<<endl;

			//cout<<"Candidates"<<endl;
			//cout<<Candidates<<endl;
			//cout<<"Op"<<endl;
			//cout<<Op<<endl;

			loglikelihood=0;
			for(vector< NRQCDglobalfitObject >::iterator state = DataModelObject.begin(); state != DataModelObject.end(); ++state){
				//cout << "getObjectLikelihood" <<endl;
				//state->Dump(NRQCDvars::nStates, true, true);
				//int nOp=Op[0].size();
				//cout<<"nOp "<<nOp<<endl;

				ObjectLikelihoodVec=state->getObjectLikelihood(Op, Np_BR, Np_US, false, directProductionCube, promptProductionMatrix, polCorrFactor);
				loglikelihood+=ObjectLikelihoodVec[0];
				//cout << "likelihood: " << loglikelihood << endl;
				//if(loglikelihood>1e99) state->Dump(NRQCDvars::nStates, true, true);
				//cout << "likelihood: " << likelihood << endl;
				if(NRQCDvars::debug) {
					cout << "iDataPoint: " << distance(DataModelObject.begin(), state) << endl;
					cout << "likelihood: " << loglikelihood << endl;
				}
			}


			//cout << "Sum_likelihood " << likelihood << endl;




			if(NRQCDvars::debug) cout<<"add chi-square terms to constrain Nps"<<endl;

			if(NRQCDvars::debug) cout<<"for Np_US"<<endl;


			for(int j=0; j < NRQCDvars::nDataSystematicScales; j++){
				if(SampleNp && FreeParam_Np_US[0][j]) loglikelihood+=((Np_US[0][j]-Np_US_Expectation[0][j])/Np_US_ExpectationUncertainty[0][j])*((Np_US[0][j]-Np_US_Expectation[0][j])/Np_US_ExpectationUncertainty[0][j]);
			}
			for(int j=0; j < NRQCDvars::nModelSystematicScales; j++){
				if(SampleNp && FreeParam_Np_US[1][j]) loglikelihood+=((Np_US[1][j]-Np_US_Expectation[1][j])/Np_US_ExpectationUncertainty[1][j])*((Np_US[1][j]-Np_US_Expectation[1][j])/Np_US_ExpectationUncertainty[1][j]);
				else Np_US[1][j]=0.;
			}

			if(NRQCDvars::debug) cout<<"for Np_BR"<<endl;
			for(int i=0; i < nStates; i++){
				for(int j=0; j < nStates; j++){
					if(NRQCDvars::FeedDownBranchingRatio[i][j]>0){
						if(SampleNp && FreeParam_Np_BR[i][j]) loglikelihood+=((Np_BR[i][j]-Np_BR_Expectation[i][j])/Np_BR_ExpectationUncertainty[i][j])*((Np_BR[i][j]-Np_BR_Expectation[i][j])/Np_BR_ExpectationUncertainty[i][j]);
					}
				}
			}

			if(NRQCDvars::debug) cout<<"for consts_star"<<endl;
			for (int i=0; i<NRQCDvars::nStates; i++){
				bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
				int nColorChannels_state;
				if(isSstate) nColorChannels_state=NRQCDvars::nColorChannels_S;
				else nColorChannels_state=NRQCDvars::nColorChannels_P;
				for (int j=0; j<nColorChannels_state; j++){
					if(SampleNp_consts_star && FreeParam_consts_star[i][j]) loglikelihood+=((consts_star_var[i][j]-consts_star_var_Expectation[i][j])/consts_star_var_ExpectationUncertainty[i][j])*((consts_star_var[i][j]-consts_star_var_Expectation[i][j])/consts_star_var_ExpectationUncertainty[i][j]);
				}
			}

			loglikelihood*=-0.5;

				//cout << "Sum_likelihood*0.5 " << likelihood << endl;

			if(NRQCDvars::debug) cout<<"likeliohood ratio criterion"<<endl;

			 loglikelihood_Candidate=loglikelihood;
				//cout << "loglikelihood_Candidate: " << loglikelihood_Candidate << endl;
				//cout << "loglikelihood_PreviousCandidate: " << loglikelihood_PreviousCandidate << endl;

			 double loglikelihood_difference = loglikelihood_Candidate - loglikelihood_PreviousCandidate;
			 if(  loglikelihood_difference > 0.  ||  loglikelihood_difference > log( gRandom->Uniform(1.) )   ) {
				 PreviousCandidates=Candidates;
				 Np_BR_PreviousCandidates=Np_BR;
				 Np_US_PreviousCandidates=Np_US;
				 consts_star_var_PreviousCandidates=consts_star_var;
				 loglikelihood_PreviousCandidate=loglikelihood_Candidate;
	//	    	 if(!BurnIn) data->Fill();
				 iSampledPoint++;
				 acceptedSampling=1;
				 iAcceptedSinceLastStep++;
				 nMH_rejectedInARow=0;
				 //cout<<"MH accepts event: YES"<<endl;
				// cout<<"loglikelihood_Candidate = "<<loglikelihood_Candidate<<endl;
				//	cout << "PreviousCandidates:"<<endl;
				//	cout << PreviousCandidates;
				//	cout << "Candidates:"<<endl;
				//	cout << Candidates;

				 //Fill BurnIn Tree:
					//for (int i=0; i<NRQCDvars::nStates; i++){
					//	bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
					//	if(isSstate){
					//		for (int j=0; j<NRQCDvars::nColorChannels_S; j++){
					//			fi_burnin[i][j]=Candidates[i][j];
					//		}
					//	}
					//	else{
					//		for (int j=0; j<NRQCDvars::nColorChannels_P; j++){
					//			fi_burnin[i][j]=Candidates[i][j];
					//		}
					//	}
					//}

				 MH_av_eff=double(iSampledPoint-1)/double(nSampledPointsTotal);




				if (iSampledPoint%nStep == 0) {
					double MH_step_eff=double(iAcceptedSinceLastStep)/double(iTotalSinceLastStep)*100;
					iTotalSinceLastStep=0;
					iAcceptedSinceLastStep=0;
					cout << nStep_ << " %, MH efficiency so far: "<<MH_av_eff*100.<<"%, since last step: "<<MH_step_eff<<"%"<<endl;
					++nStep_;
					cout << PreviousCandidates;
				}

			 }
			else{
			 //cout<<"MH accepts event: NO"<<endl;
			 acceptedSampling=0;
			 nMH_rejectedInARow++;
			 MH_av_eff=double(iSampledPoint)/double(nSampledPointsTotal);

			}
			//	cout << "PreviousCandidates:"<<endl;
			//	cout << PreviousCandidates;
			//	cout << "Candidates:"<<endl;
			//	cout << Candidates;
			//}



				if(iSampledPoint<10){
					cout << "Some output for the beginning of the sampling -> accepted candidate below? " <<acceptedSampling<<endl;
					cout << Candidates;
				}
				if(iSampledPoint==10 && acceptedSampling){
					cout << "MH Eff for the first 10 accepted samplings: " <<double(iSampledPoint)/double(nSampledPointsTotal)*100<<"%"<<endl;
				}




			 if(nMH_rejectedInARow==100) {cout<<"I'd be worried, 100 points rejected in a row... iSampledPoints: "<<iSampledPoint<<", MH_av_eff = "<<MH_av_eff*100<<"%"<<endl;}
			 if(nMH_rejectedInARow==1000) {cout<<"I'd be even more worried, 1000 points rejected in a row.. iSampledPoints: "<<iSampledPoint<<", MH_av_eff = "<<MH_av_eff*100<<"%"<<endl;}
			 if(nMH_rejectedInARow==10000) {cout<<"Sorry - 10000 points rejected in a row is enough -> EXIT"<<endl; exit(0);}

			outputTreeAllSamplings->Fill();

			if(NRQCDvars::debug) cout << "Sum loglikelihood: " << loglikelihood << endl;
		}

		outputTreeAllSamplings->Write();

	}



    ResultsFile->Close();


    ivector int_Sample_Np(1,0);
    ivector int_Sample_Np_consts_star(1,0);
    if(SampleNp) int_Sample_Np.at(0)=1;
    if(SampleNp_consts_star) int_Sample_Np_consts_star.at(0)=1;
    ivector int_Minimizer(1,0);
    int_Minimizer.at(0)=Minimizer;

	sprintf(outname,"%s/FreeParam.txt",jobdirname);
	cout<<"save FreeParam to "<<outname<<endl;

	ofstream out;
    out.open(outname);

    cout << "FreeParam_Fractions"<<endl;
    out << FreeParam_Fractions;
    cout << "FreeParam_Np_BR"<<endl;
	out << FreeParam_Np_BR;
    cout << "FreeParam_consts_star"<<endl;
	out << FreeParam_consts_star;
    cout << "FreeParam_Fractions_States"<<endl;
    out << FreeParam_Fractions_States;
    cout << "int_Sample_Np"<<endl;
	out << int_Sample_Np;
    cout << "int_Sample_Np_consts_star"<<endl;
	out << int_Sample_Np_consts_star;
    cout << "int_Minimizer"<<endl;
	out << int_Minimizer;
    cout << "FreeParam_Np_US"<<endl;
	out << FreeParam_Np_US;

    out.close();


	delete gRandom;
	//delete TestObject;

	return 0;

}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Sets the kernel function

void setKernel( const int& kernelType ){
	cout<<"setKernel"<<endl;

	switch(kernelType){
   case NRQCDvars::Metropolis:
       kernelFunction=MetropolisKernel;
       break;
   case NRQCDvars::MetropolisHastings:
       kernelFunction=MetropolisHastingsKernel;
       break;
   default:
       cerr << "Error: Unsupported kernel type! Execution stop" << endl;
       exit(1);
   }
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Sets the fractions generating function depending on dimensionality

void setFractionDimension( const int& nDim ){//nDim number of color octet channels
   switch(nDim){
   case 1:
       getFractionValues=genFractionValues1D;
       break;
   case 2:
       getFractionValues=genFractionValues2D;
       break;
   case 3:
       getFractionValues=genFractionValues3D;
       break;
   case 4:
       getFractionValues=genFractionValues4D;
       break;
   case 5:
       getFractionValues=genFractionValues5D;
       break;
   case 6:
       getFractionValues=genFractionValues6D;
       break;
   default:
       cerr << "Error: Fractions for " << nDim << "dimensions are not available! Execution stop." << endl;
       exit(1);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                //
//                                              Function to gen the fractions in 1D case (1 CO channel)                                                         //
//                                                                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void genFractionValues1D(dvector &fraction, dvector &candidate, dvector &width ){
	//basis vectors for random sampling: basis vectors in hyperplane are orthonormal

	//Generate new candidates
	for(int i=0; i < 1; ++i){
		candidate[i]=kernelFunction(candidate[i], width[i]);
	}

	fraction[0]=candidate[0];

	fraction[1]=1.;

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                //
//                                              Function to gen the fractions in 2D case (2 CO channels)                                                          //
//                                                                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void genFractionValues2D(dvector &fraction, dvector &candidate, dvector &width ){
	//basis vectors for random sampling: basis vectors in hyperplane are orthonormal
	static const double a0=1./2., a1=1./TMath::Sqrt(2.), a2=1./TMath::Sqrt(6.);


	//Generate new candidates
	for(int i=0; i < 2; ++i){
		candidate[i]=kernelFunction(candidate[i], width[i]);
	}

	fraction[0]=candidate[0];



	//TODO: Valentin changed definition:
	//fraction[1]=1./2.*(1+candidate[1]);
	//fraction[2]=1-fraction[1];
	fraction[1]=candidate[1];
	fraction[2]=1-fraction[1];

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                //
//                                              Function to gen the fractions in 3D case (3 CO channels)                                                          //
//                                                                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void genFractionValues3D(dvector &fraction, dvector &candidate, dvector &width ){
	//basis vectors for random sampling: basis vectors in hyperplane are orthonormal
	static const double a0=1./3., a1=1./TMath::Sqrt(2.), a2=1./TMath::Sqrt(6.);

	//Generate new candidates
	for(int i=0; i < 3; ++i){
		candidate[i]=kernelFunction(candidate[i], width[i]);
	}

	fraction[0]=candidate[0];

	//TODO: Valentin changed definition from this:
	//fraction[1]=a0-(a1*candidate[1]+a2*candidate[2]);
	//fraction[2]=a0+(a1*candidate[1]-a2*candidate[2]);
	//fraction[3]=a0+2.*a2*candidate[2];

	//to this:
	candidate[3]=kernelFunction(candidate[3], width[3]);
	double sum=candidate[1]+candidate[2]+candidate[3];
	fraction[1]=candidate[1]/sum;
	fraction[2]=candidate[2]/sum;
	fraction[3]=candidate[3]/sum;

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                //
//                                              Function to gen the fractions in 4D case (4 CO channels)                                                          //
//                                                                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void genFractionValues4D(dvector &fraction, dvector &candidate, dvector &width ){
	//basis vectors for random sampling: basis vectors in hyperplane are orthonormal
	static const double a0=0.25, a1=1./TMath::Sqrt(2.), a2=1./TMath::Sqrt(2.), a3=0.5;

	//Generate new candidates
	for(int i=0; i < 4; ++i){
		candidate[i]=kernelFunction(candidate[i], width[i]);
	}

	fraction[0]=candidate[0];

	fraction[1]=a0-(a1*candidate[1]-a3*candidate[3]);
	fraction[2]=a0+(a1*candidate[1]+a3*candidate[3]);
	fraction[3]=a0+(a2*candidate[2]-a3*candidate[3]);
	fraction[4]=a0-(a2*candidate[2]+a3*candidate[3]);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                //
//                                              Function to gen the fractions in 5D case (5 CO channels)                                                          //
//                                                                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void genFractionValues5D(dvector &fraction, dvector &candidate, dvector &width ){
	//basis vectors for random sampling: basis vectors in hyperplane are orthonormal
	static const double a0=0.2, a1=1./TMath::Sqrt(2.), a2=1./TMath::Sqrt(6.), a3=1./TMath::Sqrt(2.), a4=1./TMath::Sqrt(30.);

	//Generate new candidates
	for(int i=0; i < 5; ++i){
		candidate[i]=kernelFunction(candidate[i], width[i]);
	}

	fraction[0]=candidate[0];

	fraction[1]=a0-(a1*candidate[1]+a2*candidate[2])+2.*a4*candidate[4];
	fraction[2]=a0+(a1*candidate[1]-a2*candidate[2])+2.*a4*candidate[4];
	fraction[3]=a0+2.*(a2*candidate[2]+a4*candidate[4]);
	fraction[4]=a0+(a3*candidate[3]-3.*a4*candidate[4]);
	fraction[5]=a0-(a3*candidate[3]+3.*a4*candidate[4]);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                //
//                                              Function to gen the fractions in 6D case (6 CO channels)                                                          //
//                                                                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void genFractionValues6D(dvector &fraction, dvector &candidate, dvector &width ){
	//basis vectors for random sampling: basis vectors in hyperplane are orthonormal
	static const double a0=1./6., a1=1./TMath::Sqrt(2.), a2=1./TMath::Sqrt(6.), a3=1./TMath::Sqrt(2.),
			            a4=1./TMath::Sqrt(6.), a5=1./TMath::Sqrt(6.);

	//Generate new candidates
	for(int i=0; i < 6; ++i){
		candidate[i]=kernelFunction(candidate[i], width[i]);
	}

	fraction[0]=candidate[0];

	fraction[1]=a0-(a1*candidate[1]+a2*candidate[2])+a5*candidate[5];
	fraction[2]=a0+(a1*candidate[1]-a2*candidate[2])+a5*candidate[5];
	fraction[3]=a0+2.*a2*candidate[2]+a5*candidate[5];
	fraction[4]=a0+2.*a4*candidate[4]-a5*candidate[5];
	fraction[5]=a0+a3*candidate[3]-(a4*candidate[4]+a5*candidate[5]);
	fraction[6]=a0-(a3*candidate[3]+a4*candidate[4]+a5*candidate[5]);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                //
//                                              Kernel function for the Metropolis case                                                           //
//                                                                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double MetropolisKernel(const double& candidate, const double& proposalWidth ){
	double new_candidate;
		new_candidate = gRandom->Uniform( candidate-0.5*proposalWidth, candidate+0.5*proposalWidth );
	return new_candidate;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                //
//                                              Kernel function for the Metropolis-Hastings case                                                  //
//                                                                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double MetropolisHastingsKernel(const double& candidate, const double& proposalWidth ){
	double new_candidate;
		new_candidate = gRandom->Gaus ( candidate, proposalWidth );
	return new_candidate;
}


void transformFractionsToOps(dmatrix &Op, dmatrix &Fractions, dmatrix consts_star){


	int iMax = Fractions.size();
	for(int i=0; i < iMax; i++){
		//cout<<"i "<<i<<endl;
		dvector Op_state;
		Op_state.push_back(NRQCDvars::ColorSingletME[i]);
		for(dvector::iterator j = Fractions[i].begin()+1; j != Fractions[i].end(); ++j){
			int k=j-Fractions[i].begin();
			//cout<<"k "<<k<<endl;
			Op_state.push_back(Fractions[i][k]*Fractions[i][0]*consts_star[i][0]*NRQCDvars::ColorSingletME[i]/consts_star[i][k]);
		}
		Op.at(i)=Op_state;
	}
}

void FindMPV(TH1* PosteriorDist , double& MPV , double& MPVerrorLow, double& MPVerrorHigh, int MPValgo, int nSigma){

	//PosteriorDist->Print();

	bool illPPD=false;
	int iFilledBinsPDD=0;
	for(int i=0;i<PosteriorDist->GetNbinsX();i++){
		if(PosteriorDist->GetBinContent(i)>0) iFilledBinsPDD++;
	}
	if(iFilledBinsPDD<2) illPPD=true;

	if(MPValgo==1){
		MPV=PosteriorDist->GetMean();
		MPVerrorLow=PosteriorDist->GetRMS();
		MPVerrorHigh=PosteriorDist->GetRMS();
	}

	if(MPValgo==2||MPValgo==3){

		int nBins = PosteriorDist->GetNbinsX();
		int maxbin_PosteriorDist = PosteriorDist->GetMaximumBin();
		double PosteriorDist_initial = PosteriorDist->GetBinCenter(maxbin_PosteriorDist);
		double err_PosteriorDist_initial=PosteriorDist->GetRMS();
		double PosteriorDist_par [3];

		TF1 *gauss;

		int nMaxFits=1;
		if(MPValgo==3) nMaxFits=20;
		for(int iFits=0;iFits<nMaxFits;iFits++){
			gauss = new TF1("f1", "gaus", PosteriorDist_initial-err_PosteriorDist_initial, PosteriorDist_initial+err_PosteriorDist_initial);
			gauss->SetParameters(PosteriorDist_initial,err_PosteriorDist_initial);
			PosteriorDist->Fit(gauss, "R");
			gauss->GetParameters(PosteriorDist_par);
			double ndof = 2*err_PosteriorDist_initial/PosteriorDist->GetBinWidth(1)-3;
			cout<<"chi2/ndf = "<<gauss->GetChisquare()/ndof<<endl;
			PosteriorDist_initial=PosteriorDist_par[1];
			err_PosteriorDist_initial=err_PosteriorDist_initial/2;
			if(gauss->GetChisquare()/ndof<5 && gauss->GetChisquare()/ndof<5 >0) break;
			if(iFits==nMaxFits-1) illPPD=true;
		}
		MPV=PosteriorDist_par[1];

		double OneSigmaCL;
		if(nSigma==1) OneSigmaCL=0.682689492137;
		if(nSigma==2) OneSigmaCL=0.954499736104;
		if(nSigma==3) OneSigmaCL=0.997300203937;
		double fullInt=PosteriorDist->Integral(1,nBins);
		//cout<<(1-OneSigmaCL)/2.<<endl;

		for(int i = 1; i < nBins+1; i++){
			//	cout<<i<<" "<<PosteriorDist->Integral(1,i)/fullInt<<endl;
			if(PosteriorDist->Integral(1,i)/fullInt > (1-OneSigmaCL)/2.) {MPVerrorLow=MPV-PosteriorDist->GetBinCenter(i-1); break;}
		}
		for(int i = 1; i < nBins+1; i++){
			//	cout<<i<<" "<<PosteriorDist->Integral(nBins+1-i,nBins)/fullInt<<endl;
			if(PosteriorDist->Integral(nBins+1-i,nBins)/fullInt > (1-OneSigmaCL)/2.) {MPVerrorHigh=PosteriorDist->GetBinCenter(nBins-i)-MPV; break;}
		}

	}

	if(illPPD){
		MPV=PosteriorDist->GetMean();
		MPVerrorLow=1e-10;
		MPVerrorHigh=1e-10;
	}


	return;

}


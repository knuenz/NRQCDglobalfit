/*
 * CalcPPD.cc
 *
 *  Created on: Jun 14, 2013
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
#include "TRandom3.h"
#include "TSystem.h"
//#include "TMatrixD.h"

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
  	std::string filename("data");
  	std::string properties;

  	int 	nBurnIn=-1;
  	int 	nSample=-1;
  	bool 	SampleNp=false;
  	bool 	SampleNp_consts_star=false;


  	for( int i=0;i < argc; ++i ) {
		if(std::string(argv[i]).find("storagedir") != std::string::npos) {char* storagedirchar = argv[i]; char* storagedirchar2 = strtok (storagedirchar, "="); storagedir = storagedirchar2; cout<<"storagedir = "<<storagedir<<endl;}
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

							if(DataSelected) DataModelObject.push_back(readDataModelObject);

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
	int acceptedSampling;

	dvector Np_BR_0 (NRQCDvars::nStates,0.);
	for(int i=0; i < NRQCDvars::nStates; i++) Np_BR.push_back(Np_BR_0);

	dvector Np_US_0 (NRQCDvars::nDataSystematicScales, 0.);
	dvector Np_US_1 (NRQCDvars::nModelSystematicScales, 0.);
	Np_US.push_back(Np_US_0);
	Np_US.push_back(Np_US_1);
	double loglikelihood=0;
	dvector ObjectLikelihoodVec;



	int nSampledPoints=nBurnIn+nSample;
	int nStep = nSampledPoints/100;  // visualize progress of the parameter sampling
	int nStep_ = 1;

	if(nSampledPoints<10) nStep=nSampledPoints/1;

  	setKernel(NRQCDvars::MetropolisHastings ); // modified by Joao: Set environment for sampling of fractions

	//Observable parameters Op (Matrix elements)
	dmatrix Op(NRQCDvars::nStates);
	vector<double> Op_S (NRQCDvars::nColorChannels_S,0);
	vector<double> Op_P (NRQCDvars::nColorChannels_P,0);

	dmatrix Fractions(NRQCDvars::nStates);
	vector<double> Fractions_S (NRQCDvars::nColorChannels_S,0);//f0...R, fi: i going from 1 to n=nColorChannels, fn=1-sum(fi_i-(n-1))
	vector<double> Fractions_P (NRQCDvars::nColorChannels_P,0);

	dmatrix Candidates(NRQCDvars::nStates);
	vector<double> Candidates_S (NRQCDvars::nColorChannels_S,0);
	vector<double> Candidates_P (NRQCDvars::nColorChannels_P,0);

	dmatrix SampleWidths(NRQCDvars::nStates);
	vector<double> SampleWidths_S (NRQCDvars::nColorChannels_S,0);
	vector<double> SampleWidths_P (NRQCDvars::nColorChannels_P,0);

	dmatrix PreviousCandidates (NRQCDvars::nStates);
	vector<double> PreviousCandidates_S (NRQCDvars::nColorChannels_S,0);
	vector<double> PreviousCandidates_P (NRQCDvars::nColorChannels_P,0);

	dmatrix consts_star;
	dmatrix err_consts_star;

	//const_star matrix varied by uncertainty (Nuisance parameter inthe sampling)
	dmatrix consts_star_var(NRQCDvars::nStates);
	vector<double> consts_star_var_S (NRQCDvars::nColorChannels_S,0);
	vector<double> consts_star_var_P (NRQCDvars::nColorChannels_P,0);


	sprintf(inname,"%s/ModelIngredients_consts_star.txt",modeldirname);
	cout<<"read in dmatrix consts_star from "<<inname<<endl;
	ifstream instar;
	instar.open(inname);
	instar >> consts_star;
	instar >> err_consts_star;
	instar.close();

	cout << consts_star << endl;
	cout << err_consts_star << endl;

	double loglikelihood_PreviousCandidate = -1.e30;  // intial (arbitrary) values
	double loglikelihood_Candidate = -1.e30;  // intial (arbitrary) values

	bool BurnIn=true;

	cout<<"set starting point of Candidates:"<<endl;

	for (int i=0; i<NRQCDvars::nStates; i++){
		bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
		if(isSstate){
			for (int j=0; j<NRQCDvars::nColorChannels_S; j++){
				if(j==0) Candidates_S.at(j)=1.;
				else Candidates_S.at(j)=1./double(NRQCDvars::nColorChannels_S-1);
			}
			Candidates.at(i)=Candidates_S;
		}
		else{
			for (int j=0; j<NRQCDvars::nColorChannels_P; j++){
				if(j==0) Candidates_P.at(j)=1.;
				else Candidates_P.at(j)=1./double(NRQCDvars::nColorChannels_P-1);
			}
			Candidates.at(i)=Candidates_P;
		}
	}
	PreviousCandidates=Candidates;

	cout<<Candidates<<endl;

	cout<<"Initialize Op's"<<endl;

	for (int i=0; i<NRQCDvars::nStates; i++){
		bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
		if(isSstate){
			for (int j=0; j<NRQCDvars::nColorChannels_S; j++){
				Op_S.at(j)=0.;
			}
			Op.at(i)=Op_S;
		}
		else{
			for (int j=0; j<NRQCDvars::nColorChannels_P; j++){
				Op_P.at(j)=0.;
			}
			Op.at(i)=Op_P;
		}
	}

	cout<<"Initialize consts_star_var's"<<endl;

	for (int i=0; i<NRQCDvars::nStates; i++){
		bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
		if(isSstate){
			for (int j=0; j<NRQCDvars::nColorChannels_S; j++){
				consts_star_var_S.at(j)=0.;
			}
			consts_star_var.at(i)=consts_star_var_S;
		}
		else{
			for (int j=0; j<NRQCDvars::nColorChannels_P; j++){
				consts_star_var_P.at(j)=0.;
			}
			consts_star_var.at(i)=consts_star_var_P;
		}
	}

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

  	TTree*  outputTreeAllSamplings = new TTree("AllSamplings","AllSamplings"); // tree filled in all samplings after burnin
  	TTree*  outputTreeAccSamplings = new TTree("AccSamplings","AccSamplings"); // tree filled only in accepted samplings after burnin


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
				outputTreeAccSamplings->Branch(branch_name,     &Candidates[i][j],     branch_address);
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
				outputTreeAccSamplings->Branch(branch_name,     &Op[i][j],     branch_address);
			}
		}


	int iAccSampling;
	sprintf(branch_name,"iAccSampling");
	sprintf(branch_address,"%s/I",branch_name);
	outputTreeAccSamplings->Branch(branch_name,     &iAccSampling,     branch_address);



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
				outputTreeAccSamplings->Branch(branch_name,     &consts_star_var[i][j],     branch_address);
			}
		}

	//Make TTree to be filled with BurnIn samplings

	char BurnInSamplingsName[200];
	sprintf(BurnInSamplingsName,"BurnInSamplingsName");
	TTree *BurnInSamplings = new TTree(BurnInSamplingsName, BurnInSamplingsName);


	for (int i=0; i<NRQCDvars::nStates; i++){
		bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
		if(isSstate){
			for (int j=0; j<NRQCDvars::nColorChannels_S; j++){
				sprintf(branch_name,"state%d_f%d_burnin",i,j);
				sprintf(branch_address,"%s/D",branch_name);
				BurnInSamplings->Branch(branch_name,     &Candidates[i][j],     branch_address);
			}
		}
		else{
			for (int j=0; j<NRQCDvars::nColorChannels_P; j++){
				sprintf(branch_name,"state%d_f%d_burnin",i,j);
				sprintf(branch_address,"%s/D",branch_name);
				BurnInSamplings->Branch(branch_name,     &Candidates[i][j],     branch_address);
			}
		}
	}

	sprintf(branch_name,"iAccSampling");
	sprintf(branch_address,"%s/I",branch_name);
	BurnInSamplings->Branch(branch_name,     &iAccSampling,     branch_address);


////////////////////////////////////////////////
/////// END DEFINE OUTPUT //////////////////////
////////////////////////////////////////////////

	cout<<"Start sampling ("<<nBurnIn<<" + "<<nSample<<" points)"<<endl;
	cout<<"Progress:"<<endl;

	int nMH_rejectedInARow=0;
	int iTotalSinceLastStep=0;
	int iAcceptedSinceLastStep=0;

	iAccSampling=0;
	int nSampledPointsTotal=1;
	for(int iSampledPoint = 1; iSampledPoint <= nSampledPoints; nSampledPointsTotal++){ // Sampling loop

		//if(nSampledPointsTotal>2) break;

		//cout << "iSampledPoint " << iSampledPoint <<endl;

		iAccSampling=iSampledPoint;
		iTotalSinceLastStep++;

		if(iSampledPoint==nBurnIn && BurnIn){
			BurnIn=false;
			cout<<"BurnIn period finished"<<endl;
			double fi_burnin_in[NRQCDvars::nStates][NRQCDvars::nColorChannels];
			char fi_name_in[200];
			TH1D *h_BurnIn[NRQCDvars::nStates][NRQCDvars::nColorChannels];
			int nBins_fi_burnin=100;

			for (int i=0; i<NRQCDvars::nStates; i++){
				int nColorChannels_state;
				bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
				if(isSstate) nColorChannels_state=NRQCDvars::nColorChannels_S;
				else nColorChannels_state=NRQCDvars::nColorChannels_P;
				for (int j=0; j<nColorChannels_state; j++){
					sprintf(fi_name_in,"state%d_f%d_burnin",i,j);
					BurnInSamplings->SetBranchAddress( fi_name_in,  &fi_burnin_in[i][j] );
					h_BurnIn[i][j] = new TH1D (fi_name_in, fi_name_in, nBins_fi_burnin, BurnInSamplings->GetMinimum(fi_name_in), BurnInSamplings->GetMaximum(fi_name_in));
				}
			}

			int n_events = int( BurnInSamplings->GetEntries() );
			int nBinsh_pT=50;
			int nBinsh_rap=50;

			// loop over  events in the burnin ntuple
			for ( int i_event = 1; i_event <= n_events; i_event++ ) {

				BurnInSamplings->GetEvent( i_event-1 );

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
						if(j==0 && newWidth<NRQCDvars::proposalWidthBurnIn_R) SampleWidths_S[j]=newWidth;
						if(j>1 && newWidth<NRQCDvars::proposalWidthBurnIn_f) SampleWidths_S[j]=newWidth;
						SampleWidths_S[j]=newWidth/2.;

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
						if(j==0 && newWidth<NRQCDvars::proposalWidthBurnIn_R) SampleWidths_P[j]=newWidth;
						if(j>1 && newWidth<NRQCDvars::proposalWidthBurnIn_f) SampleWidths_P[j]=newWidth;
						SampleWidths_P[j]=newWidth/2.;

						//SampleWidths_P[j]=h_BurnIn[i][j]->GetRMS();
						cout<<"Proposal width for state "<<i<<", ColorChannel "<<j<<" = "<<SampleWidths_P[j]<<endl;
					}
					SampleWidths.at(i)=SampleWidths_P;
				}
			}
			BurnInSamplings->Write();
		}

		// Fill Nuisance parameter matrices Np_BR and Np_US. Np_US[0]: Data-related, Np_US[1]: Model-related uncertainty scales
		//cout << "Fill Nuisance parameter matrices" <<endl;

		for(int j=0; j < NRQCDvars::nDataSystematicScales; j++){
			if(SampleNp) Np_US[0][j]=gRandom->Gaus(0,1); // Luminosity scaling
			else Np_US[0][j]=0.;
		}
		for(int j=0; j < NRQCDvars::nModelSystematicScales; j++){
			if(SampleNp) Np_US[1][j]=gRandom->Gaus(0,1);
			else Np_US[1][j]=0.;
		}

		for(int i=0; i < nStates; i++){
			for(int j=0; j < nStates; j++){
				if(NRQCDvars::FeedDownBranchingRatio[i][j]>0){
					if(SampleNp) Np_BR[i][j]=gRandom->Gaus( NRQCDvars::FeedDownBranchingRatio[i][j]  , NRQCDvars::errFeedDownBranchingRatio[i][j])*0.01;
					else Np_BR[i][j]=NRQCDvars::FeedDownBranchingRatio[i][j]*0.01;
				}
				else Np_BR[i][j]=0;
			}
		}

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
				getFractionValues(Fractions_S, PreviousCandidates_S, SampleWidths_S);
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
			  	getFractionValues(Fractions_P, PreviousCandidates_P, SampleWidths_P);
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

		//cout<<"varying consts_star by uncertainty"<<endl;

		for (int i=0; i<NRQCDvars::nStates; i++){
			bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
			if(isSstate){
				for (int j=0; j<NRQCDvars::nColorChannels_S; j++){
					if(SampleNp_consts_star) consts_star_var_S.at(j)=gRandom->Gaus(consts_star[i][j],err_consts_star[i][j]);
					else consts_star_var_S.at(j)=consts_star[i][j];
				}
				consts_star_var.at(i)=consts_star_var_S;
			}
			else{
				for (int j=0; j<NRQCDvars::nColorChannels_P; j++){
					if(SampleNp_consts_star) consts_star_var_P.at(j)=gRandom->Gaus(consts_star[i][j], err_consts_star[i][j]);
					else consts_star_var_P.at(j)=consts_star[i][j];
				}
				consts_star_var.at(i)=consts_star_var_P;
			}
		}

		transformFractionsToOps(Op, Candidates, consts_star_var);


		//cout<<Op;
        //
		//cout<<"Op[0].size() "<<Op[0].size()<<endl;
		//cout<<"Op[1].size() "<<Op[1].size()<<endl;

		if(NRQCDvars::debug) cout << "getObjectLikelihood" << endl;

		loglikelihood=0;
		for(vector< NRQCDglobalfitObject >::iterator state = DataModelObject.begin(); state != DataModelObject.end(); ++state){
			//cout << "getObjectLikelihood" <<endl;
			//state->Dump(NRQCDvars::nStates, true, true);
			//int nOp=Op[0].size();
			//cout<<"nOp "<<nOp<<endl;

			ObjectLikelihoodVec=state->getObjectLikelihood(Op, Np_BR, Np_US, false, directProductionCube, promptProductionMatrix, polCorrFactor);
			loglikelihood+=ObjectLikelihoodVec[0];
			//cout << "likelihood: " << likelihood << endl;
			if(NRQCDvars::debug) {
				cout << "iDataPoint: " << distance(DataModelObject.begin(), state) << endl;
				cout << "likelihood: " << loglikelihood << endl;
			}
		}


		//cout << "Sum_likelihood " << likelihood << endl;

		loglikelihood*=-0.5;

			//cout << "Sum_likelihood*0.5 " << likelihood << endl;

		 loglikelihood_Candidate=loglikelihood;
			//cout << "loglikelihood_Candidate: " << loglikelihood_Candidate << endl;
			//cout << "loglikelihood_PreviousCandidate: " << loglikelihood_PreviousCandidate << endl;

	     double loglikelihood_difference = loglikelihood_Candidate - loglikelihood_PreviousCandidate;
	     if(  loglikelihood_difference > 0.  ||  loglikelihood_difference > log( gRandom->Uniform(1.) )   ) {
	    	 PreviousCandidates=Candidates;
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

	    	 MH_av_eff=double(iSampledPoint)/double(nSampledPointsTotal);

	    		if(BurnIn) BurnInSamplings->Fill();
	    		else outputTreeAccSamplings->Fill();



	 		if (iSampledPoint%nStep == 0) {
	 			double MH_step_eff=double(iAcceptedSinceLastStep+1)/double(iTotalSinceLastStep)*100;
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



	     if(nMH_rejectedInARow==100) cout<<"I'd be worried, 100 points rejected in a row... keep you posted!"<<endl;
	     if(nMH_rejectedInARow==1000) cout<<"I'd be even more worried, 1000 points rejected in a row... keep you posted!"<<endl;
	     if(nMH_rejectedInARow==10000) {cout<<"Sorry - 10000 points rejected in a row is enough -> EXIT"<<endl; exit(0);}

		    if(!BurnIn) outputTreeAllSamplings->Fill();

		if(NRQCDvars::debug) cout << "Sum loglikelihood: " << loglikelihood << endl;
	}

    outputTreeAllSamplings->Write();
    outputTreeAccSamplings->Write();
    ResultsFile->Close();

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

	fraction[1]=a0-(a1*candidate[1]+a2*candidate[2]);
	fraction[2]=a0+(a1*candidate[1]-a2*candidate[2]);
	fraction[3]=a0+2.*a2*candidate[2];
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

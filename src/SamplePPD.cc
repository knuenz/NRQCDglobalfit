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
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TRandom3.h"
#include "TSystem.h"
//#include "TMatrixD.h"

using namespace NRQCDvars;


int main(int argc, char** argv) {

  	Char_t *JobID = "Default";
  	Char_t *DataID = "Default";

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

  	for( int i=0;i < argc; ++i ) {
		if(std::string(argv[i]).find("JobID") != std::string::npos) {
			char* JobIDchar = argv[i];
			char* JobIDchar2 = strtok (JobIDchar, "=");
			JobID = JobIDchar2;
			cout<< "_JobID = " << JobID << endl;
			properties.append("_JobID");
			properties.append(JobIDchar2);
		}
		if(std::string(argv[i]).find("DataID") != std::string::npos) {
			char* DataIDchar = argv[i];
			char* DataIDchar2 = strtok (DataIDchar, "=");
			DataID = DataIDchar2;
			cout<<"DataID = "<<DataID<<endl;
			properties.append("_DataID");
			properties.append(DataID);
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
	    	useOnlyState = atof(useOnlyStatechar2);
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

	for(int iState=0; iState < nStates; iState++){
		for(int iMeasurementID=0; iMeasurementID < NRQCDvars::nMeasurementIDs; iMeasurementID++){
			for(int iExperiment=0; iExperiment < NRQCDvars::nExperiments; iExperiment++){
				for(int iRap = 0; iRap < NRQCDvars::nMaxRapBins; iRap++){
				    for(int iP = 0; iP < NRQCDvars::nMaxPtBins; iP++){

				    	//if(iState!=4 || iMeasurementID!=0 || iExperiment!=2 || iRap!=0 || iP!=0) continue;

						char inname[2000];
						char datadirname[2000];
						sprintf(datadirname,"DataID");
						gSystem->mkdir(datadirname);
						sprintf(datadirname,"%s/%s",datadirname,DataID);
						gSystem->mkdir(datadirname);
						sprintf(inname,"%s/ConvertedDataModel_%s_%s_%s_rap%d_pT%d.txt",datadirname, StateName[iState],
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
							readDataModelObject.initFractions(nDim); //TODO (?): dimension of fractions vector to be read from state

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

							if(!DataSelected){
								readDataModelObject.Dump(nStates, true, true);
							}
							else DataModelObject.push_back(readDataModelObject);
							DataSelected=false;
						}
				    }
				}
			}
		}
	}

	int nDataPoints=DataModelObject.size();
	cout << "Number of considered data points: " << nDataPoints << endl;


	//Observable parameters Op (Matrix elements)
	dmatrix Op(nStates);

	//Nuisance parameters Np

	dmatrix Np_BR; //Branching ratios, [nDauthgers][nMothers]
	dmatrix Np_US; //Uncertainty scales, [0=Data, 1=Model][nScales]

	dvector Np_BR_0 (nStates,0.);
	for(int i=0; i < nStates; i++) Np_BR.push_back(Np_BR_0);

	dvector Np_US_0 (NRQCDvars::nDataSystematicScales, 0.);
	dvector Np_US_1 (NRQCDvars::nModelSystematicScales, 0.);
	Np_US.push_back(Np_US_0);
	Np_US.push_back(Np_US_1);
	double likelihood=0;

  	setFractionsSetup(NRQCDvars::nFractionDim, NRQCDvars::MetropolisHastings ); // modified by Joao: Set environment for sampling of fractions

  	TTree*  data = new TTree("data","data"); // data tree

  	data->Branch("likelihood", &likelihood, "likelihood/D");
  	int fractionDim=nFractionDim;
  	data->Branch("nFractionDim", &fractionDim, "nFractionDim/I");
//  	data->Branch("fraction", fraction, "fraction[nFractionDim}/D"); //TODO: How do we keep the fractions in the TTree?
//  	data->Branch("nSampledPoints", &nSampledPoints, "nSampledPoints/I");~
//  	data->Branch("nStates", &nStates, "nStates/I");
//  	data->Branch("nColorChannels", &nColorChannels, "nColorChannels/I");
  	data->Branch("Np_BR", &Np_BR);
  	data->Branch("Np_US", &Np_US);
  	data->Branch("Op", &Op);

	int nSampledPoints=100;

	int nStep = nSampledPoints/10;  // visualize progress of the parameter sampling
	int nStep_ = 1;

	cout<<"Start sampling"<<endl;
	cout<<"Progress:"<<endl;

	for(int iSampledPoint = 1; iSampledPoint <= nSampledPoints; iSampledPoint++){ // Sampling loop

		if (iSampledPoint%nStep == 0) {
			cout << nStep_*10 << " % " <<endl; ++nStep_;
		}


		///// This code is to be run in the loop over all samplings:::

		// Fill Nuisance parameter matrices Np_BR and Np_US. Np_US[0]: Data-related, Np_US[1]: Model-related uncertainty scales

		Np_US[0][0]=gRandom->Gaus( 0.,1. ); // Luminosity scaling

		for(int j=0; j < NRQCDvars::nModelSystematicScales; j++){
			Np_US[1][j]=gRandom->Gaus(0,1);
		}

		for(int i=0; i < nStates; i++){
			for(int j=0; j < nStates; j++){
				Np_BR[i][j]=gRandom->Gaus( NRQCDvars::FeedDownBranchingRatio[i][j]  , NRQCDvars::errFeedDownBranchingRatio[i][j])*0.01;
			}
		}

		///// Set toyMatrixElements

		if(NRQCDvars::debug) cout<< "fill Op" << endl;

		vector<double> Op_state(NRQCDvars::nColorChannels,0);
		for (int i=0; i<NRQCDvars::nStates; i++){
			for (int j=0; j<NRQCDvars::nColorChannels; j++){
				if(j==0) Op_state[j]=gRandom->Gaus(NRQCDvars::ColorSingletME[i], NRQCDvars::errColorSingletME[i]); //Just to set it to some value for debugging purposes
				else Op_state[j]=gRandom->Uniform(0,1);
			}
			Op[i]=Op_state;
		}

		if(NRQCDvars::debug) cout << "getObjectLikelihood" << endl;

		likelihood=0;
		for(vector< NRQCDglobalfitObject >::iterator state = DataModelObject.begin(); state != DataModelObject.end(); ++state){
			getFractionValues(state->fraction, state->sampleValues, state->sampleWidths);
			likelihood+= state->getObjectLikelihood(Op, Np_BR, Np_US, state->fraction);
			if(NRQCDvars::debug) {
				cout << "iDataPoint: " << distance(DataModelObject.begin(), state) << endl;
				cout << "likelihood: " << likelihood << endl;
			}
		}
		likelihood*=-0.5;
		data->Fill();
		if(NRQCDvars::debug) cout << "Sum likelihood: " << likelihood << endl;
	}

	delete gRandom;
	//delete TestObject;

	return 0;

}

void (*getFractionValues)(double* fraction, double* candidate, double* width );
double (*kernelFunction)(double& par0, double& par1 );

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Sets the quantities to calculate with operator fractions

void setFractionsSetup(const int& nDim, const int& kernelType ){
	switch(kernelType){                     // Set kernel function
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
	switch(nDim){                          // Set the fractions generating function depending on dimensionality
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
//                                              Function to gen the fractions in 3D case                                                          //
//                                                                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void genFractionValues3D(double* fraction, double* candidate, double* width ){
	//basis vectors for random sampling: basis vectors in hyperplane are orthonormal
	static const double a0=1./3., a1=1./TMath::Sqrt(2.), a2=1./TMath::Sqrt(6.);

	//Generate new candidates
	for(int i=1; i < 3; ++i){
		candidate[i]=kernelFunction(candidate[i], width[i]);
	}
	fraction[0]=a0-(a1*candidate[1]+a2*candidate[2]);
	fraction[1]=a0+(a1*candidate[1]-a2*candidate[2]);
	fraction[2]=a0+2.*a2*candidate[2];
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                //
//                                              Function to gen the fractions in 4D case                                                          //
//                                                                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void genFractionValues4D(double* fraction, double* candidate, double* width ){
	//basis vectors for random sampling: basis vectors in hyperplane are orthonormal
	static const double a0=0.25, a1=1./TMath::Sqrt(2.), a2=1./TMath::Sqrt(2.), a3=0.5;

	//Generate new candidates
	for(int i=1; i < 4; ++i){
		candidate[i]=kernelFunction(candidate[i], width[i]);
	}
	fraction[0]=a0-(a1*candidate[1]-a3*candidate[3]);
	fraction[1]=a0+(a1*candidate[1]+a3*candidate[3]);
	fraction[2]=a0+(a2*candidate[2]-a3*candidate[3]);
	fraction[3]=a0-(a2*candidate[2]+a3*candidate[3]);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                //
//                                              Function to gen the fractions in 5D case                                                          //
//                                                                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void genFractionValues5D(double* fraction, double* candidate, double* width ){
	//basis vectors for random sampling: basis vectors in hyperplane are orthonormal
	static const double a0=0.2, a1=1./TMath::Sqrt(2.), a2=1./TMath::Sqrt(6.), a3=1./TMath::Sqrt(2.), a4=1./TMath::Sqrt(30.);

	//Generate new candidates
	for(int i=1; i < 5; ++i){
		candidate[i]=kernelFunction(candidate[i], width[i]);
	}
	fraction[0]=a0-(a1*candidate[1]+a2*candidate[2])+2.*a4*candidate[4];
	fraction[1]=a0+(a1*candidate[1]-a2*candidate[2])+2.*a4*candidate[4];
	fraction[2]=a0+2.*(a2*candidate[2]+a4*candidate[4]);
	fraction[3]=a0+(a3*candidate[3]-3.*a4*candidate[4]);
	fraction[4]=a0-(a3*candidate[3]+3.*a4*candidate[4]);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                //
//                                              Function to gen the fractions in 6D case                                                          //
//                                                                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void genFractionValues6D(double* fraction, double* candidate, double* width ){
	//basis vectors for random sampling: basis vectors in hyperplane are orthonormal
	static const double a0=1./6., a1=1./TMath::Sqrt(2.), a2=1./TMath::Sqrt(6.), a3=1./TMath::Sqrt(2.),
			            a4=1./TMath::Sqrt(6.), a5=1./TMath::Sqrt(6.);

	//Generate new candidates
	for(int i=1; i < 6; ++i){
		candidate[i]=kernelFunction(candidate[i], width[i]);
	}
	fraction[0]=a0-(a1*candidate[1]+a2*candidate[2])+a5*candidate[5];
	fraction[1]=a0+(a1*candidate[1]-a2*candidate[2])+a5*candidate[5];
	fraction[2]=a0+2.*a2*candidate[2]+a5*candidate[5];
	fraction[3]=a0+2.*a4*candidate[4]-a5*candidate[5];
	fraction[4]=a0+a3*candidate[3]-(a4*candidate[4]+a5*candidate[5]);
	fraction[5]=a0-(a3*candidate[3]+a4*candidate[4]+a5*candidate[5]);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                //
//                                              Kernel function for the Metropolis case                                                           //
//                                                                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double MetropolisKernel(const double& candidate, const double& proposalWidth ){
	if ( proposalWidth < 0. ) {
	   do {                      // Uniform proposal pdf with "large" sigma
		   candidate = gRandom->Uniform( candidate-0.5*NRQCDvars::proposalWidthBurnIn, candidate+0.5*NRQCDvars::proposalWidthBurnIn );
	   }
	   while ( candidate < 0 );
	}
	else {                      // Uniform proposal pdf with possibly smaller sigmas
		candidate = gRandom->Uniform( candidate-0.5*proposalWidth, candidate+0.5*proposalWidth );; // Gaussian proposal pdf with possibly smaller sigmas
	}
	return candidate;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                //
//                                              Kernel function for the Metropolis-Hastings case                                                  //
//                                                                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double MetropolisHastingsKernel(double& candidate, double& proposalWidth ){
	if ( proposalWidth < 0. ) {
	   do {
		   candidate = gRandom->Gaus( candidate, NRQCDvars::proposalWidthBurnIn ); // Gaussian proposal pdf with "large" sigma
	   }
	   while ( candidate < 0 );
	}
	else {
		candidate = gRandom->Gaus ( candidate, proposalWidth ); // Gaussian proposal pdf with possibly smaller sigmas
	}
	return candidate;
}

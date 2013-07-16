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

  	int 	nBurnIn=-1;
  	int 	nSample=-1;


  	for( int i=0;i < argc; ++i ) {
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

	for(int iState=0; iState < nStates; iState++){
		for(int iMeasurementID=0; iMeasurementID < NRQCDvars::nMeasurementIDs; iMeasurementID++){
			for(int iExperiment=0; iExperiment < NRQCDvars::nExperiments; iExperiment++){
				for(int iRap = 0; iRap < NRQCDvars::nMaxRapBins; iRap++){
				    for(int iP = 0; iP < NRQCDvars::nMaxPtBins; iP++){

				    	if(iState!=0 || iMeasurementID!=0 || iExperiment!=0 || iRap!=0 || iP!=0) continue;

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
							//readDataModelObject.initFractions(nDim); //TODO (?): dimension of fractions vector to be read from state

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


	//TODO: open model file, get sigma_star, c_star, close model file



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
//  	data->Branch("Op", &Op);


	int nSampledPoints=nBurnIn+nSample;
	int nStep = nSampledPoints/10;  // visualize progress of the parameter sampling
	int nStep_ = 1;


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

	dmatrix consts_star(NRQCDvars::nStates);
	vector<double> consts_star_S (NRQCDvars::nColorChannels_S,0);//0...sigma_star_CS, i going from 1 to n=nColorChannels, representing c_i star for a given state and pT_star
	vector<double> consts_star_P (NRQCDvars::nColorChannels_P,0);


	double loglikelihood_PreviousCandidate = -1.e30;  // intial (arbitrary) values
	double loglikelihood_Candidate = -1.e30;  // intial (arbitrary) values

	// Set all matrices to a starting value

	cout<<"set starting point of Candidates"<<endl;

	for (int i=0; i<NRQCDvars::nStates; i++){
		bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
		if(isSstate){
			for (int j=0; j<NRQCDvars::nColorChannels_S; j++){
				consts_star_S.at(j)=1.;//TODO: define starting values for these constants
			}
			consts_star.at(i)=consts_star_S;
		}
		else{
			for (int j=0; j<NRQCDvars::nColorChannels_P; j++){
				consts_star_P.at(j)=1.;
			}
			consts_star.at(i)=consts_star_P;
		}
	}

	//TODO: set proposal width matrices for burn in
	bool BurnIn=true;

	cout<<"set starting point of Candidates"<<endl;

	for (int i=0; i<NRQCDvars::nStates; i++){
		bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
		if(isSstate){
			for (int j=0; j<NRQCDvars::nColorChannels_S; j++){
				Candidates_S.at(j)=0.;//TODO: define starting values for Candidates
			}
			Candidates.at(i)=Candidates_S;
		}
		else{
			for (int j=0; j<NRQCDvars::nColorChannels_P; j++){
				Candidates_P.at(j)=0.;
			}
			Candidates.at(i)=Candidates_P;
		}
	}

	cout<<"set starting point of Op's"<<endl;

	for (int i=0; i<NRQCDvars::nStates; i++){
		bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
		if(isSstate){
			for (int j=0; j<NRQCDvars::nColorChannels_S; j++){
				Op_S.at(j)=0.;//TODO: define starting values for Candidates
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

	//cout << "Candidates:"<<endl;
	//cout << Candidates;

	cout<<"set widths"<<endl;

	for (int i=0; i<NRQCDvars::nStates; i++){
		bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
		if(isSstate){
			for (int j=0; j<NRQCDvars::nColorChannels_S; j++){
				SampleWidths_S[j]=-1.;
			}
			SampleWidths.at(i)=SampleWidths_S;
		}
		else{
			for (int j=0; j<NRQCDvars::nColorChannels_P; j++){
				SampleWidths_P[j]=-1.;
			}
			SampleWidths.at(i)=SampleWidths_P;
		}
	}

	int nSampledPointsTotal=1;
	for(int iSampledPoint = 1; iSampledPoint <= nSampledPoints; nSampledPointsTotal++){ // Sampling loop

		//if(nSampledPointsTotal>2) break;

		cout << "iSampledPoint " << iSampledPoint <<endl;

	//	if (iSampledPoint%nStep == 0) {
	//		cout << nStep_*10 << " % " <<endl; ++nStep_;
	//	}
    //
	//	if(iSampledPoint==nBurnIn){
	//		BurnIn=false;
	//		cout<<"BurnIn period finished"<<endl;
    //
	//		for (int i=0; i<NRQCDvars::nStates; i++){
	//			bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
	//			if(isSstate){
	//				for (int j=0; j<NRQCDvars::nColorChannels_S; j++){
	//					//TODO: implement histos, take rms as estimate of proposal width
	//					SampleWidths_S[j]=NRQCDvars::proposalWidthBurnIn;
	//				}
	//				SampleWidths.at(i)=SampleWidths_S;
	//			}
	//			else{
	//				for (int j=0; j<NRQCDvars::nColorChannels_P; j++){
	//					SampleWidths_P[j]=NRQCDvars::proposalWidthBurnIn;
	//				}
	//				SampleWidths.at(i)=SampleWidths_P;
	//			}
	//		}
    //
	//		break;
	//	}
    //
	//	///// This code is to be run in the loop over all samplings:::
    //
	//	// Fill Nuisance parameter matrices Np_BR and Np_US. Np_US[0]: Data-related, Np_US[1]: Model-related uncertainty scales
		//cout << "Fill Nuisance parameter matrices" <<endl;

		Np_US[0][0]=gRandom->Gaus( 0.,1. ); // Luminosity scaling

		for(int j=0; j < NRQCDvars::nModelSystematicScales; j++){
			Np_US[1][j]=gRandom->Gaus(0,1);
		}

		for(int i=0; i < nStates; i++){
			for(int j=0; j < nStates; j++){
				Np_BR[i][j]=gRandom->Gaus( NRQCDvars::FeedDownBranchingRatio[i][j]  , NRQCDvars::errFeedDownBranchingRatio[i][j])*0.01;
			}
		}

		///// Set Fractions -> calculate Matrix Elements Op

		if(NRQCDvars::debug) cout<< "fill Op" << endl;

		//cout << "PreviousCandidates=Candidates;" <<endl;

		PreviousCandidates=Candidates;

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



		// Relate Op to R, fi -> getObjectLikelihood
		//cout << "transformFractionsToOps" <<endl;

		transformFractionsToOps(Op, Fractions, consts_star);


		//cout<<Op;
        //
		//cout<<"Op[0].size() "<<Op[0].size()<<endl;
		//cout<<"Op[1].size() "<<Op[1].size()<<endl;

		if(NRQCDvars::debug) cout << "getObjectLikelihood" << endl;

		likelihood=0;
		for(vector< NRQCDglobalfitObject >::iterator state = DataModelObject.begin(); state != DataModelObject.end(); ++state){
			//cout << "getObjectLikelihood" <<endl;
			//state->Dump(NRQCDvars::nStates, true, true);
			//int nOp=Op[0].size();
			//cout<<"nOp "<<nOp<<endl;

			likelihood+= state->getObjectLikelihood(Op, Np_BR, Np_US);
			//cout << "likelihood: " << likelihood << endl;
			if(NRQCDvars::debug) {
				cout << "iDataPoint: " << distance(DataModelObject.begin(), state) << endl;
				cout << "likelihood: " << likelihood << endl;
			}
		}

		//TODO: calculate log-likelihood, make likelihood criterion decision

		//cout << "Sum_likelihood " << likelihood << endl;

		 likelihood*=0.5;

			//cout << "Sum_likelihood*0.5 " << likelihood << endl;

		 loglikelihood_Candidate=log(likelihood);
			cout << "loglikelihood_Candidate: " << loglikelihood_Candidate << endl;
			cout << "loglikelihood_PreviousCandidate: " << loglikelihood_PreviousCandidate << endl;

	     double loglikelihood_difference = loglikelihood_Candidate - loglikelihood_PreviousCandidate;
	     if(  loglikelihood_difference > 0.  ||  log( gRandom->Uniform(1.) ) < loglikelihood_difference  ) {
	    	 PreviousCandidates=Candidates;
	    	 loglikelihood_PreviousCandidate=loglikelihood_Candidate;
//	    	 if(!BurnIn) data->Fill();
	    	 iSampledPoint++;
	    	 cout<<"MH accepts event: YES"<<endl;
	    	// cout<<"loglikelihood_Candidate = "<<loglikelihood_Candidate<<endl;
	    	//	cout << "PreviousCandidates:"<<endl;
	    	//	cout << PreviousCandidates;
	    	//	cout << "Candidates:"<<endl;
	    	//	cout << Candidates;

	     }
	    else{
	     cout<<"MH accepts event: NO"<<endl;
	    }
	    //	cout << "PreviousCandidates:"<<endl;
	    //	cout << PreviousCandidates;
	    //	cout << "Candidates:"<<endl;
	    //	cout << Candidates;
	    //}

		if(NRQCDvars::debug) cout << "Sum likelihood: " << likelihood << endl;
	}

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

	fraction[1]=1./2.*(1+candidate[1]);
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
	if ( proposalWidth < 0. ) {// Uniform proposal pdf with "large" sigma
		new_candidate = gRandom->Uniform( candidate-0.5*NRQCDvars::proposalWidthBurnIn, candidate+0.5*NRQCDvars::proposalWidthBurnIn );
	}
	else {                      // Uniform proposal pdf with possibly smaller sigmas
		new_candidate = gRandom->Uniform( candidate-0.5*proposalWidth, candidate+0.5*proposalWidth );; // Gaussian proposal pdf with possibly smaller sigmas
	}
	return new_candidate;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                //
//                                              Kernel function for the Metropolis-Hastings case                                                  //
//                                                                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double MetropolisHastingsKernel(const double& candidate, const double& proposalWidth ){
	double new_candidate;
	if ( proposalWidth < 0. ) {
		new_candidate = gRandom->Gaus( candidate, NRQCDvars::proposalWidthBurnIn ); // Gaussian proposal pdf with "large" sigma
	}
	else {
		new_candidate = gRandom->Gaus ( candidate, proposalWidth ); // Gaussian proposal pdf with possibly smaller sigmas
	}
	//cout<<"candidate: "<<candidate<<endl;
	//cout<<"new_candidate: "<<new_candidate<<endl;
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
			Op_state.push_back(Fractions[i][k]*Fractions[i][0]*consts_star[i][0]/consts_star[i][k]);
		}
		Op.at(i)=Op_state;
	}
}

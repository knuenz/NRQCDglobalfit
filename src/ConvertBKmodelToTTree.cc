/*
 * ConvertBKmodelToTTree.cc
 *
 *  Created on: Aug 10, 2013
 *      Author: valentinknuenz
 */


#include "../interface/NRQCDglobalfitObject.h"
#include "../src/NRQCDglobalfitObject.cc"
#include "../interface/GlobalVar.h"

#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>
#include <fstream>

//rootincludes
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TF2.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TSystem.h"
#include "Riostream.h"
#include "TString.h"
#include "TROOT.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TRotation.h"
#include "TGraph.h"


using namespace NRQCDvars;

dvector func_lam_gen(int iMother, int iColorChannel);
double func_pT_gen(double* x, double* par);

int main(int argc, char** argv) {

  	Char_t *OriginalModelID = "Default";
  	Char_t *ModelID = "Default";
	Char_t *storagedir = "Default"; //Storage Directory

  	for( int i=0;i < argc; ++i ) {
		if(std::string(argv[i]).find("OriginalModelID") != std::string::npos) {char* OriginalModelIDchar = argv[i]; char* OriginalModelIDchar2 = strtok (OriginalModelIDchar, "="); OriginalModelID = OriginalModelIDchar2; cout<<"OriginalModelID = "<<OriginalModelID<<endl;}
		if(std::string(argv[i]).find("ModelID") != std::string::npos) {char* ModelIDchar = argv[i]; char* ModelIDchar2 = strtok (ModelIDchar, "="); ModelID = ModelIDchar2; cout<<"ModelID = "<<ModelID<<endl;}
		if(std::string(argv[i]).find("storagedir") != std::string::npos) {char* storagedirchar = argv[i]; char* storagedirchar2 = strtok (storagedirchar, "="); storagedir = storagedirchar2; cout<<"storagedir = "<<storagedir<<endl;}
  	}

	gRandom = new TRandom3(NRQCDvars::randomSeed);  // better random generator

	char outname[2000];
	char inname[2000];
	char premodeldirname[2000];
	char modeldirname[2000];
	char originalmodeldirname[2000];


	sprintf(premodeldirname,"%s/ModelID", storagedir);
	gSystem->mkdir(premodeldirname);
	sprintf(modeldirname,"%s/%s",premodeldirname,ModelID);
	gSystem->mkdir(modeldirname);

	sprintf(premodeldirname,"%s/OriginalModelID",storagedir);
	gSystem->mkdir(premodeldirname);
	sprintf(originalmodeldirname,"%s/%s",premodeldirname,OriginalModelID);
	gSystem->mkdir(originalmodeldirname);


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Convert BK input model
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	const int nRapIntervals=3;
	bool isAbsRap[nRapIntervals]={true, true, false};
	double RapIntervalBordersMin[nRapIntervals]={0, 0.6, 2.};
	double RapIntervalBordersMax[nRapIntervals]={0.6, 1.2, 4.5};


		sprintf(inname,"%s/TGraphs_scaled.root",originalmodeldirname);
		TFile* inputFile = new TFile(inname, "READ");
		cout<<"opened TGraph file"<<endl;

		TGraph *SDC_Graph[NRQCDvars::nStates][nRapIntervals][NRQCDvars::nColorChannels];
		TGraph *Lamth_Graph[NRQCDvars::nStates][nRapIntervals][NRQCDvars::nColorChannels];
		TGraph *Lamph_Graph[NRQCDvars::nStates][nRapIntervals][NRQCDvars::nColorChannels];
		TGraph *Lamtp_Graph[NRQCDvars::nStates][nRapIntervals][NRQCDvars::nColorChannels];

		char nameSDCgraph[200];
		char nameLamthgraph[200];
		char nameLamphgraph[200];
		char nameLamtpgraph[200];

		double model_pTMin_state[NRQCDvars::nStates], model_pTMax_state[NRQCDvars::nStates];
		double model_pTMin_staterap[NRQCDvars::nStates][nRapIntervals], model_pTMax_staterap[NRQCDvars::nStates][nRapIntervals];

		bool ModelGiven[NRQCDvars::nStates];

		for(int i=0;i<NRQCDvars::nStates;i++){
			model_pTMin_state[i]=10000;
			model_pTMax_state[i]=0;
			ModelGiven[i]=false;

			for(int j=0;j<nRapIntervals;j++){
				int nColorChannels_state;
				bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
				if(isSstate) nColorChannels_state=NRQCDvars::nColorChannels_S;
				else nColorChannels_state=NRQCDvars::nColorChannels_P;
				for (int k=0; k<nColorChannels_state; k++){

					char TGid[200];
					sprintf(TGid,"original");

					sprintf(nameSDCgraph,"SDC_Graph_%s_state%d_rap%d_CC%d",TGid,i,j,k);
					SDC_Graph[i][j][k]=(TGraph*)inputFile->Get(nameSDCgraph);
					sprintf(nameSDCgraph,"Lamth_Graph_%s_state%d_rap%d_CC%d",TGid,i,j,k);
					Lamth_Graph[i][j][k]=(TGraph*)inputFile->Get(nameSDCgraph);
					sprintf(nameSDCgraph,"Lamph_Graph_%s_state%d_rap%d_CC%d",TGid,i,j,k);
					Lamph_Graph[i][j][k]=(TGraph*)inputFile->Get(nameSDCgraph);
					sprintf(nameSDCgraph,"Lamtp_Graph_%s_state%d_rap%d_CC%d",TGid,i,j,k);
					Lamtp_Graph[i][j][k]=(TGraph*)inputFile->Get(nameSDCgraph);

			                cout<<"opened TGraphs for state"<<i<<" rap"<<j<<" CC"<<k<<endl;

					if(inputFile->Get(nameSDCgraph)!=NULL)
			                        ModelGiven[i]=true;
					else continue;

					double buff_pT, buff_val;
					SDC_Graph[i][j][k]->GetPoint(0, buff_pT, buff_val);
					if(buff_pT<model_pTMin_state[i])
						model_pTMin_state[i]=buff_pT;
					model_pTMin_staterap[i][j]=buff_pT;
					SDC_Graph[i][j][k]->GetPoint(SDC_Graph[i][j][k]->GetN()-1, buff_pT, buff_val);
					if(buff_pT>model_pTMax_state[i])
						model_pTMax_state[i]=buff_pT;
					model_pTMax_staterap[i][j]=buff_pT;

			}
		}
	}




	sprintf(outname,"%s/OriginalModelTTree.root",modeldirname);
	TFile* ModelIngredientsFile = new TFile(outname, "RECREATE");

	double model_pTMin, model_pTMax;
	double model_rapMin, model_rapMax;
	double model_costhMin, model_costhMax;
	double model_phiMin, model_phiMax;

	//model_pTMin=3;
	//model_pTMax=40;
	model_rapMin=-1.2;
	model_rapMax=4.5;
	model_costhMin=-1;
	model_costhMax=1;
	model_phiMin=-180;
	model_phiMax=180;


	int n_nTuple=1000000;

	TTree* nTupleModel[NRQCDvars::nStates][NRQCDvars::nColorChannels];
	char nTupleModelName[1000];

	for (int iMother=0; iMother<NRQCDvars::nStates; iMother++){
		if(!ModelGiven[iMother]) continue;
		if(iMother!=0 && iMother!=3&& iMother!=4 && iMother!=7 && iMother!=10)
		cout<<"Converting original model nTuple for nState = "<<iMother<<endl;
		int nColorChannels_state;
		bool isSstate=(StateQuantumID[iMother] > NRQCDvars::quID_S)?false:true;
		if(isSstate) nColorChannels_state=NRQCDvars::nColorChannels_S;
		else nColorChannels_state=NRQCDvars::nColorChannels_P;
		for (int iColorChannel=0; iColorChannel<nColorChannels_state; iColorChannel++){
			cout<<"		Color channel = "<<iColorChannel<<endl;


			sprintf(nTupleModelName,"nTupleModel_Mother%d_Daughter%d_ColorChannel%d_Cascade0",iMother, iMother, iColorChannel);
			nTupleModel[iMother][iColorChannel] = new TTree(nTupleModelName, nTupleModelName);

			Double_t model_pT;  nTupleModel[iMother][iColorChannel]->Branch("model_pT",     &model_pT,     "model_pT/D");
			Double_t model_rap;    nTupleModel[iMother][iColorChannel]->Branch("model_rap",       &model_rap,       "model_rap/D"  );
			Double_t model_costh;  nTupleModel[iMother][iColorChannel]->Branch("model_costh",      &model_costh,     "model_costh/D");
			Double_t model_phi;  nTupleModel[iMother][iColorChannel]->Branch("model_phi",      &model_phi,     "model_phi/D");
			Double_t model_weight; nTupleModel[iMother][iColorChannel]->Branch("weight",      &model_weight,     "weight/D");

			int nStep = n_nTuple/10;  // visualize progress of the parameter sampling
			int nStep_ = 0;
			bool genAngDist2D=true;

			//if(iColorChannel==3) n_nTuple=100000;

			for (int k=0; k<n_nTuple; k++){



				model_pT=gRandom->Uniform(model_pTMin_state[iMother],model_pTMax_state[iMother]);
				model_rap=gRandom->Uniform(model_rapMin,model_rapMax);
				model_costh=gRandom->Uniform(model_costhMin,model_costhMax);
				model_phi=gRandom->Uniform(model_phiMin,model_phiMax);

				int rapIndex = -1;
				for (int j=0; j<nRapIntervals; j++){
					if(isAbsRap[j]){
						if(fabs(model_rap)>RapIntervalBordersMin[j] && fabs(model_rap)<RapIntervalBordersMax[j]) {rapIndex=j; break;}
					}
					else{
						if(model_rap>RapIntervalBordersMin[j] && model_rap<RapIntervalBordersMax[j]) {rapIndex=j; break;}
					}
					}


				//cout<<"model_rap "<<model_rap<<endl;
				//cout<<"rapIndex "<<rapIndex<<endl;
				//cout<<"model_pT "<<model_pT<<endl;
				bool FillEvent=true;
				if(rapIndex==-1) FillEvent=false;

				if(model_pT<model_pTMin_staterap[iMother][rapIndex] || model_pT>model_pTMax_staterap[iMother][rapIndex]) FillEvent=false;

				if(FillEvent){


					model_weight=SDC_Graph[iMother][rapIndex][iColorChannel]->Eval(model_pT);
					if(log(model_weight)>8 && TMath::Abs(model_rap)<1.2&&model_pT>12&&model_pT<13.5) cout<<"model_weight pT"<<model_weight<<endl;
					//Define angular distribution:
					//cout<<"model_weight pT"<<model_weight<<endl;

					//dvector model_lam = func_lam_gen(iMother, iColorChannel);
					double model_lamth, model_lamph, model_lamtp;
					model_lamth=Lamth_Graph[iMother][rapIndex][iColorChannel]->Eval(model_pT);
					model_lamph=Lamph_Graph[iMother][rapIndex][iColorChannel]->Eval(model_pT);
					model_lamtp=Lamtp_Graph[iMother][rapIndex][iColorChannel]->Eval(model_pT);

					if(model_lamth>-4.&&model_lamth<-2.){
						double polDecision = gRandom->Uniform(-1,1);
						if(polDecision>0) model_lamth=-2.;
						else model_lamth=-4.;
					}

					//cout<<"model_lamth"<<model_lamth<<endl;
					//cout<<"model_lamph"<<model_lamph<<endl;
					//cout<<"model_lamtp"<<model_lamtp<<endl;

					//cout<<"model_costh"<<model_costh<<endl;
					//cout<<"model_phi"<<model_phi<<endl;

					double polNormFactor;
					if(genAngDist2D){
						TF2 *fcosthphi;
						fcosthphi = new TF2( "fcosthphi", "[0]*(1.+[1]*x[0]*x[0]+[2]*(1.-x[0]*x[0])*cos(2.*x[1]*0.0174532925)+[3]*2.*x[0]*sqrt(1.-x[0]*x[0])*cos(x[1]*0.0174532925))", -1., 1., -180., 180. );
						fcosthphi->SetParameters(3./(3.+model_lamth),model_lamth, model_lamph, model_lamtp);
						model_weight*=fcosthphi->Eval(model_costh,model_phi);
						polNormFactor=fcosthphi->Integral(-1., 1., -180., 180. );
						model_weight*=2.*360./fcosthphi->Integral(-1., 1., -180., 180. );
						if(2.*360./fcosthphi->Integral(-1., 1., -180., 180. )-3>1e-2) cout<<"model_weight pol norm not 3 but "<<2.*360./fcosthphi->Integral(-1., 1., -180., 180. )<<endl;
						if(log(model_weight)>8 && TMath::Abs(model_rap)<1.2&&model_pT>12&&model_pT<13.5){
							cout<<"model_weight pT "<<SDC_Graph[iMother][rapIndex][iColorChannel]->Eval(model_pT)<<endl;
							cout<<"model_lamth "<<model_lamth<<endl;
							cout<<"model_weight pol function "<<fcosthphi->Eval(model_costh,model_phi)<<endl;
							cout<<"model_weight pol norm "<<2.*360./fcosthphi->Integral(-1., 1., -180., 180. )<<endl;
							cout<<"resulting weight "<<model_weight<<endl;
						}

						delete fcosthphi;
					}
					if(!genAngDist2D){
						TF1 *fcosth;
						fcosth = new TF1( "fcosth", "[0]*(1.+[1]*x[0]*x[0])", -1., 1.);
						fcosth->SetParameters(1.,model_lamth);
						model_weight*=fcosth->Eval(model_costh);
						model_weight*=2./fcosth->Integral(-1., 1.);

						delete fcosth;
					}
					//cout<<"model_weight pol"<<model_weight<<endl;

					nTupleModel[iMother][iColorChannel]->Fill();

			 		if (k%nStep == 0) {
			 			cout << nStep_*10 <<"% (nState = "<<iMother<<", Color channel = "<<iColorChannel<<")"<<endl;
			 			++nStep_;
			 		}

				}
				else k--;




			}


			bool interpretModelIntegratedInRapidityInterval=true;

			double deltay=model_rapMax-model_rapMin;

			double phasespaceFactor=0;//30.*1.2+20.*1.2; //sum_i ( deltaPt_i*deltaRap_i ) with i...rap bin
			for(int j=0;j<nRapIntervals;j++){
				double deltaPt=model_pTMax_staterap[iMother][j]-model_pTMin_staterap[iMother][j];
				double deltaY=RapIntervalBordersMax[j]-RapIntervalBordersMin[j];
				if(isAbsRap[j]) deltaY*=2;

				if(interpretModelIntegratedInRapidityInterval) phasespaceFactor+=deltaPt;
				else phasespaceFactor+=deltaPt*deltaY;
			}

			double globalWeight=1./double(n_nTuple)*phasespaceFactor;
			nTupleModel[iMother][iColorChannel]->SetWeight(globalWeight);

			cout<<"globalWeight "<<globalWeight<<endl;

			cout<<"Generated "<<n_nTuple<<" events"<<endl;
			cout<<"nTupleModel[iMother][iColorChannel]->GetEntries() "<<nTupleModel[iMother][iColorChannel]->GetEntries()<<endl;


			nTupleModel[iMother][iColorChannel]->Write();


		}
	}



	ModelIngredientsFile->Close();

  	return 0;

  	}



double func_pT_gen(double* x, double* par) {


	double funcval=1.;

	int nState=par[0];
	int nColorChannel=par[1];

	//if(nColorChannel==0){
	//	double beta=-1.52;
	//	funcval=pow(x[0],beta);
	//}
    //
	//if(nColorChannel==1){
	//	double alpha=-1e-5;
	//	double gamma=1e-2;
	//	funcval=alpha*x[0]+gamma;
	//}
    //
	//if(nColorChannel==2){
	//	double gamma=0.0075;
	//	funcval=gamma;
	//}
    //
	//if(nColorChannel==3){
	//	double alpha=1e-5;
	//	double gamma=5e-3;
	//	funcval=alpha*x[0]+gamma;
	//}


	if(nColorChannel==0){
		double beta;
		double pTsq;
		double constfact = 9;
		beta=3;
		pTsq=10;//NRQCDvars::mass[nState]*NRQCDvars::mass[nState];
		funcval = constfact * x[0] * pow( 1. + 1./(beta - 2.) * x[0]*x[0] / pTsq, -beta  );
	}
	if(nColorChannel==1){
		double beta=-2;
		double constfact = 1.01;
		funcval= constfact * pow(x[0],beta);
	}
	if(nColorChannel==2){
		double beta=-1;
		double constfact = 0.05;
		funcval= constfact * pow(x[0],beta);
	}





	return funcval;

}

dvector func_lam_gen(int iMother, int iColorChannel){

	double model_lamth, model_lamph, model_lamtp;

	bool isSstate=(StateQuantumID[iMother] > NRQCDvars::quID_S)?false:true;

	if(isSstate){
		if(iColorChannel==0) {model_lamth=-0.5; model_lamph=0; model_lamtp=0;}
		if(iColorChannel==1) {model_lamth=+0.5; model_lamph=0; model_lamtp=0;}
		if(iColorChannel==2) {model_lamth=+1; model_lamph=0; model_lamtp=0;}
	}
	else{
		if(iColorChannel==0) {model_lamth=-1./3.; model_lamph=0; model_lamtp=0;}
		if(iColorChannel==1) {model_lamth=-1./3.; model_lamph=0; model_lamtp=0;}
	}

	dvector model_lam(3,0);
	model_lam.at(0)=model_lamth;
	model_lam.at(1)=model_lamph;
	model_lam.at(2)=model_lamtp;

	return model_lam;

}

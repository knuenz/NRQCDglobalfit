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

// Matrix elements:

	const int nStatesGiven=2;//Given models in any form
	int StatesGiven[nStatesGiven]={0, 3};
	const int nColorChannelsGiven=4;

	const int nCalcModelForStates=1;//Model is only generated for these states
	int CalcModelForStates[nCalcModelForStates]={3};

	double fittedLDMEs[nStatesGiven][nColorChannelsGiven]={
			{1.32, 0.0304, 0.00168, -0.00908},
			{0.76, -0.00247, 0.00280, 0.00168}
	};

	double errfittedLDMEs[nStatesGiven][nColorChannelsGiven]={
			{0., 0.0035, 0.00046, 0.00161},
			{0., 0.00370, 0.00049, 0.00185}
	};

	const int nHelicityChannels=3;	//00, 11, 1m1

	const int nRapIntervals=3;
	bool isAbsRap[nRapIntervals]={true, true, false};
	double RapIntervalBordersMin[nRapIntervals]={0, 0.6, 2.};
	double RapIntervalBordersMax[nRapIntervals]={0.6, 1.2, 4.5};
	const int npTBinsPerRap[nRapIntervals][nHelicityChannels]={
			{9, 9, 9},
			{7, 7, 7},
			{8, 10, 11}
	};
	const int maxpTBinsPerRap=11;

	//Scale charmonium models to bottomonium with mass ratio
	bool ScaleCharmModelToBottom=true;
	int ScaleFromCharmState=3;
	const int nScaleToBottomStates=1;
	int ScaleToBottomStates[nScaleToBottomStates]={10};
	double MassRatioCtB[nScaleToBottomStates];

	for(int iCtB=0; iCtB<nScaleToBottomStates; iCtB++){
		MassRatioCtB[iCtB]=NRQCDvars::mass[ScaleToBottomStates[iCtB]]/NRQCDvars::mass[ScaleFromCharmState];
	}




// CMS rapidities




	double pTmean[nStatesGiven][nRapIntervals][nHelicityChannels][maxpTBinsPerRap];//={
	//		{
	//				{10.00, 12.50, 15.00, 17.50, 20.00, 25.00, 30.00, 35.00, 40.00},
	//				{10.00, 12.50, 15.00, 17.50, 20.00, 25.00, 30.00, 00.00, 00.00}
	//		},
	//		{
	//				{10.00, 12.50, 15.00, 17.50, 20.00, 25.00, 30.00, 35.00, 40.00},
	//				{10.00, 12.50, 15.00, 17.50, 20.00, 25.00, 30.00, 00.00, 00.00}
	//		}
	//};

	double pTMin[nRapIntervals]={10, 10, 3};
	double pTMax[nRapIntervals]={40, 30, 10};

	double model[nStatesGiven][nRapIntervals][nColorChannelsGiven][nHelicityChannels][maxpTBinsPerRap];
	double total_model[nStatesGiven][nRapIntervals][nHelicityChannels][maxpTBinsPerRap];

	cout<<"read in file"<<endl;

	for(int j=0;j<nRapIntervals;j++){

		cout<<"rap"<<j<<endl;
		sprintf(inname,"%s/BKmodel_Psi_rap%d.txt",originalmodeldirname,j+1);
		FILE *fIn;
		fIn = fopen(inname, "read");
		cout<<inname<<endl;
		Char_t line[1000];

		for(int l=0;l<nHelicityChannels;l++){
			cout<<"helicity"<<l<<endl;
			for(int m=0;m<npTBinsPerRap[j][l];m++){

				double buffer;
				fgets(line, sizeof(line), fIn); //comment
				cout<<line<<endl;
				sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &pTmean[0][j][l][m], &model[0][j][0][l][m], &model[0][j][1][l][m], &model[0][j][2][l][m], &model[0][j][3][l][m], &total_model[0][j][l][m], &model[1][j][0][l][m], &model[1][j][1][l][m], &model[1][j][2][l][m], &model[1][j][3][l][m], &total_model[1][j][l][m]);

				pTmean[1][j][l][m]=pTmean[0][j][l][m];

				cout<<"ratio sum_cc/total (1S) = "<<(model[0][j][0][l][m]+model[0][j][1][l][m]+model[0][j][2][l][m]+model[0][j][3][l][m])/total_model[0][j][l][m]<<endl;
				cout<<"ratio sum_cc/total (2S) = "<<(model[1][j][0][l][m]+model[1][j][1][l][m]+model[1][j][2][l][m]+model[1][j][3][l][m])/total_model[1][j][l][m]<<endl;

			}
			fgets(line, sizeof(line), fIn); //comment
		}
	}




	TGraph *model_Graph[nStatesGiven][nRapIntervals][nColorChannelsGiven][nHelicityChannels];

	for(int i=0;i<nStatesGiven;i++){
		for(int j=0;j<nRapIntervals;j++){
			for(int k=0;k<nColorChannelsGiven;k++){
				for(int l=0;l<nHelicityChannels;l++){

					double pTmean_graph[npTBinsPerRap[j][l]];
					double model_graph[npTBinsPerRap[j][l]];
					for(int m=0;m<npTBinsPerRap[j][l];m++){
						pTmean_graph[m]=pTmean[i][j][l][m];
						model_graph[m]=model[i][j][k][l][m];
					}

					model_Graph[i][j][k][l] = new TGraph(npTBinsPerRap[j][l],pTmean_graph,model_graph);

					cout<<"model_Graph[state"<<i<<"][rap"<<j<<"][CC"<<k<<"][HC"<<l<<"]"<<endl;
					model_Graph[i][j][k][l]->Print();

				}
			}
		}
	}


	double scaleFactor[nColorChannelsGiven];
	int countAverage;
		for(int k=0;k<nColorChannelsGiven;k++){
			cout<<"CC "<<k<<endl;
			scaleFactor[k]=0;
			countAverage=0;
			for(int l=0;l<nHelicityChannels;l++){
				cout<<"HC "<<l<<endl;
				for(int j=0;j<nRapIntervals-1;j++){
					cout<<"rap "<<j<<endl;

				double meanRapRatio_1ov0=0;
				int nMissedBecauseOfZeroDemon=0;
				for(int m=0;m<npTBinsPerRap[j][l];m++){
					double RapRatio_1ov0;
					double buffx[nStatesGiven], buffy[nStatesGiven];
					model_Graph[0][j][k][l]->GetPoint(m,buffx[0], buffy[0]);
					model_Graph[1][j][k][l]->GetPoint(m,buffx[1], buffy[1]);
					if(buffy[0]!=0){
					RapRatio_1ov0=buffy[1]/buffy[0];
					meanRapRatio_1ov0+=RapRatio_1ov0;
					}
					else nMissedBecauseOfZeroDemon++;
					//cout<<"RapRatio_1ov0 "<<RapRatio_1ov0<<endl;
				}
				if(npTBinsPerRap[j][l]-nMissedBecauseOfZeroDemon!=0){
					meanRapRatio_1ov0/=(npTBinsPerRap[j][l]-nMissedBecauseOfZeroDemon);
				}
				else meanRapRatio_1ov0=0;
				cout<<"----------------------------"<<endl;
				cout<<"meanRapRatio_1ov0 "<<meanRapRatio_1ov0<<endl;
				cout<<"----------------------------"<<endl;

				if(meanRapRatio_1ov0!=0){
				scaleFactor[k]+=meanRapRatio_1ov0;
				countAverage++;
				}
			}
		}
			scaleFactor[k]/=double(countAverage);
			cout<<"(sigma) scaleFactor[k] "<<scaleFactor[k]<<endl;
			cout<<"(SDC) scaleFactor[k]*fittedLDMEs[0][k]/fittedLDMEs[1][k] "<<scaleFactor[k]*fittedLDMEs[0][k]/fittedLDMEs[1][k]<<endl;
	}

		//scaling of sigmas (from Jpsi to psi2S) does *not* depend on rapidity, helicity channel, but only on the color channel
		// SDC2S/SDC1S=1 -> scaleFactor[k] is equal to fittedLDMEs[1][k]/fittedLDMEs[0][k]

		cout<<"Rescaling psi2S models for LHCb rapidities"<<endl;

		for(int i=1;i<2;i++){
			for(int j=2;j<3;j++){
				for(int k=0;k<nColorChannelsGiven;k++){
					for(int l=0;l<nHelicityChannels;l++){

						double pTmean_graph[npTBinsPerRap[j][l]];
						double model_graph[npTBinsPerRap[j][l]];
						for(int m=0;m<npTBinsPerRap[j][l];m++){
							pTmean_graph[m]=pTmean[i][j][l][m];
							model_graph[m]=model[0][j][k][l][m]*scaleFactor[k];
							model[1][j][k][l][m]=model[0][j][k][l][m]*scaleFactor[k];
						}

						model_Graph[i][j][k][l] = new TGraph(npTBinsPerRap[j][l],pTmean_graph,model_graph);

						cout<<"model_Graph[state"<<i<<"][rap"<<j<<"][CC"<<k<<"][HC"<<l<<"]"<<endl;
						model_Graph[i][j][k][l]->Print();

					}
				}
			}
		}



	TGraph *SDC_Graph[nStatesGiven][nRapIntervals][nColorChannelsGiven];
	TGraph *Lamth_Graph[nStatesGiven][nRapIntervals][nColorChannelsGiven];
	TGraph *Lamph_Graph[nStatesGiven][nRapIntervals][nColorChannelsGiven];
	TGraph *Lamtp_Graph[nStatesGiven][nRapIntervals][nColorChannelsGiven];

	bool UseChaoPolDenominatorDef=false;

	int nBinsFinalModels=100;
	for(int i=0;i<nStatesGiven;i++){
		for(int j=0;j<nRapIntervals;j++){
			for(int k=0;k<nColorChannelsGiven;k++){

				cout<<"rap "<<j<<" CC "<<k<<" state "<<i<<endl;
				double pTmean_graph[nBinsFinalModels];
				double SDC_graph[nBinsFinalModels];
				double Lamth_graph[nBinsFinalModels];
				double Lamph_graph[nBinsFinalModels];
				double Lamtp_graph[nBinsFinalModels];
				for(int m=0;m<nBinsFinalModels;m++){
					double deltaPt=(pTMax[j]-pTMin[j])/double(nBinsFinalModels);
					pTmean_graph[m]=pTMin[j]+m*deltaPt;
					SDC_graph[m]=model_Graph[i][j][k][0]->Eval(pTmean_graph[m])/fittedLDMEs[i][k]+2*model_Graph[i][j][k][1]->Eval(pTmean_graph[m])/fittedLDMEs[i][k];
					double polDenominator=(model_Graph[i][j][k][1]->Eval(pTmean_graph[m])/fittedLDMEs[i][k]+model_Graph[i][j][k][0]->Eval(pTmean_graph[m])/fittedLDMEs[i][k]);
					if(UseChaoPolDenominatorDef && model_Graph[i][j][k][1]->Eval(pTmean_graph[m]) * model_Graph[i][j][k][0]->Eval(pTmean_graph[m]) < 0) polDenominator=fabs(polDenominator);
					Lamth_graph[m]=(model_Graph[i][j][k][1]->Eval(pTmean_graph[m])/fittedLDMEs[i][k]-model_Graph[i][j][k][0]->Eval(pTmean_graph[m])/fittedLDMEs[i][k])/polDenominator;
					Lamph_graph[m]=(model_Graph[i][j][k][2]->Eval(pTmean_graph[m])/fittedLDMEs[i][k])/polDenominator;
					Lamtp_graph[m]=0;
				}

				SDC_Graph[i][j][k]= new TGraph(nBinsFinalModels,pTmean_graph,SDC_graph);
				Lamth_Graph[i][j][k]= new TGraph(nBinsFinalModels,pTmean_graph,Lamth_graph);
				Lamph_Graph[i][j][k]= new TGraph(nBinsFinalModels,pTmean_graph,Lamph_graph);
				Lamtp_Graph[i][j][k]= new TGraph(nBinsFinalModels,pTmean_graph,Lamtp_graph);

				cout<<"SDC_Graph[i][j][k]"<<endl;
				SDC_Graph[i][j][k]->Print();
				cout<<"Lamth_Graph[i][j][k]"<<endl;
				Lamth_Graph[i][j][k]->Print();
				cout<<"Lamph_Graph[i][j][k]"<<endl;
				Lamph_Graph[i][j][k]->Print();
				cout<<"Lamtp_Graph[i][j][k]"<<endl;
				Lamtp_Graph[i][j][k]->Print();
			}
		}
	}

	double evalPt[6]={10.5, 11.48, 12.72, 14.23, 16.3, 21.92};
	for(int u=0;u<6;u++){
		double crossSection=0;
		for(int k=0;k<nColorChannelsGiven;k++){
			crossSection+=SDC_Graph[1][1][k]->Eval(evalPt[u])*fittedLDMEs[1][k];
		}
		cout<<"crossSection at "<<evalPt[u]<<" = "<<crossSection<<endl;
	}



	//x[0]=10.5, y[0]=5.43118
	//x[1]=11.48, y[1]=3.89092
	//x[2]=12.72, y[2]=2.33473
	//x[3]=14.23, y[3]=1.45768
	//x[4]=16.3, y[4]=0.73896
	//x[5]=21.92, y[5]=0.168861


	for(int i=0;i<nStatesGiven;i++){
		cout<<"state "<<i<<endl;
		for(int k=0;k<nColorChannelsGiven;k++){
			cout<<"CC "<<k<<endl;

			double RapRatio_1ov0;
			double RapRatio_2ov1;
			RapRatio_1ov0=SDC_Graph[i][1][k]->Eval(10.)/SDC_Graph[i][0][k]->Eval(10.);
			RapRatio_2ov1=SDC_Graph[i][2][k]->Eval(10.)/SDC_Graph[i][1][k]->Eval(10.);
			cout<<"RapRatio_1ov0 "<<RapRatio_1ov0<<" RapRatio_2ov1 "<<RapRatio_2ov1<<endl;

			double RapDiffLamth_1m0;
			double RapDiffLamth_2m1;
			RapDiffLamth_1m0=Lamth_Graph[i][1][k]->Eval(10.)-Lamth_Graph[i][0][k]->Eval(10.);
			RapDiffLamth_2m1=Lamth_Graph[i][2][k]->Eval(10.)-Lamth_Graph[i][1][k]->Eval(10.);
			cout<<"RapDiffLamth_1m0 "<<RapDiffLamth_1m0<<" RapDiffLamth_2m1 "<<RapDiffLamth_2m1<<endl;
			double RapDiffLamph_1m0;
			double RapDiffLamph_2m1;
			RapDiffLamph_1m0=Lamph_Graph[i][1][k]->Eval(10.)-Lamph_Graph[i][0][k]->Eval(10.);
			RapDiffLamph_2m1=Lamph_Graph[i][2][k]->Eval(10.)-Lamph_Graph[i][1][k]->Eval(10.);
			cout<<"RapDiffLamph_1m0 "<<RapDiffLamph_1m0<<" RapDiffLamph_2m1 "<<RapDiffLamph_2m1<<endl;

		}
	}





	sprintf(outname,"%s/OriginalModelTTree.root",modeldirname);
	TFile* ModelIngredientsFile = new TFile(outname, "RECREATE");

	double model_pTMin, model_pTMax;
	double model_rapMin, model_rapMax;
	double model_costhMin, model_costhMax;
	double model_phiMin, model_phiMax;

	model_pTMin=3;
	model_pTMax=40;
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
		cout<<"Converting original model nTuple for nState = "<<iMother<<endl;
		bool CalcState=false;
		for (int jMother=0; jMother<nCalcModelForStates; jMother++){
			if(CalcModelForStates[jMother]==iMother) CalcState=true;
		}
		if(!CalcState) continue;
		int StatesGivenIndex=-1;
		for (int jMother=0; jMother<nStatesGiven; jMother++){
			if(StatesGiven[jMother]==iMother) StatesGivenIndex=jMother;

		}
		if(StatesGivenIndex==-1) continue;

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



				model_pT=gRandom->Uniform(model_pTMin,model_pTMax);
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
				if(model_pT<pTMin[rapIndex] || model_pT>pTMax[rapIndex]) FillEvent=false;

				if(FillEvent){


					model_weight=SDC_Graph[StatesGivenIndex][rapIndex][iColorChannel]->Eval(model_pT);
					if(log(model_weight)>8 && TMath::Abs(model_rap)<1.2&&model_pT>12&&model_pT<13.5) cout<<"model_weight pT"<<model_weight<<endl;
					//Define angular distribution:
					//cout<<"model_weight pT"<<model_weight<<endl;

					//dvector model_lam = func_lam_gen(iMother, iColorChannel);
					double model_lamth, model_lamph, model_lamtp;
					model_lamth=Lamth_Graph[StatesGivenIndex][rapIndex][iColorChannel]->Eval(model_pT);
					model_lamph=Lamph_Graph[StatesGivenIndex][rapIndex][iColorChannel]->Eval(model_pT);
					model_lamtp=Lamtp_Graph[StatesGivenIndex][rapIndex][iColorChannel]->Eval(model_pT);

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
							cout<<"model_weight pT "<<SDC_Graph[StatesGivenIndex][rapIndex][iColorChannel]->Eval(model_pT)<<endl;
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

			double deltapT=model_pTMax-model_pTMin;
			double deltay=model_rapMax-model_rapMin;

			double phasespaceFactor=0;//30.*1.2+20.*1.2; //sum_i ( deltaPt_i*deltaRap_i ) with i...rap bin
			for(int j=0;j<nRapIntervals;j++){
				double deltaPt=pTMax[j]-pTMin[j];
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



	if(ScaleCharmModelToBottom){
	int iMother=ScaleFromCharmState;
		for (int iScaleToMother=0; iScaleToMother<NRQCDvars::nStates; iScaleToMother++){
			bool ScaleToState=false;
			cout<<"Scaling original model from nState = "<<iMother<<" to nState = "<<iScaleToMother<<endl;
			int ScaleStateIndex=-1;
			for (int jMother=0; jMother<nScaleToBottomStates; jMother++){
				if(ScaleToBottomStates[jMother]==iScaleToMother) {ScaleToState=true; ScaleStateIndex=jMother;}
			}
			if(!ScaleToState) continue;
			int StatesGivenIndex=-1;
			for (int jMother=0; jMother<nStatesGiven; jMother++){
				if(StatesGiven[jMother]==iMother) StatesGivenIndex=jMother;
			}
			if(StatesGivenIndex==-1) continue;
			cout<<"MassRatioCtB[ScaleStateIndex] "<<MassRatioCtB[ScaleStateIndex]<<endl;

			int nColorChannels_state;
			bool isSstate=(StateQuantumID[iMother] > NRQCDvars::quID_S)?false:true;
			if(isSstate) nColorChannels_state=NRQCDvars::nColorChannels_S;
			else nColorChannels_state=NRQCDvars::nColorChannels_P;
			for (int iColorChannel=0; iColorChannel<nColorChannels_state; iColorChannel++){
				cout<<"		Color channel = "<<iColorChannel<<endl;


				sprintf(nTupleModelName,"nTupleModel_Mother%d_Daughter%d_ColorChannel%d_Cascade0",iScaleToMother, iScaleToMother, iColorChannel);
				nTupleModel[iScaleToMother][iColorChannel] = new TTree(nTupleModelName, nTupleModelName);


				Double_t model_pT;  nTupleModel[iScaleToMother][iColorChannel]->Branch("model_pT",     &model_pT,     "model_pT/D");
				Double_t model_rap;    nTupleModel[iScaleToMother][iColorChannel]->Branch("model_rap",       &model_rap,       "model_rap/D"  );
				Double_t model_costh;  nTupleModel[iScaleToMother][iColorChannel]->Branch("model_costh",      &model_costh,     "model_costh/D");
				Double_t model_phi;  nTupleModel[iScaleToMother][iColorChannel]->Branch("model_phi",      &model_phi,     "model_phi/D");
				Double_t model_weight; nTupleModel[iScaleToMother][iColorChannel]->Branch("weight",      &model_weight,     "weight/D");

				int nStep = n_nTuple/10;  // visualize progress of the parameter sampling
				int nStep_ = 0;
				bool genAngDist2D=true;

				//if(iColorChannel==3) n_nTuple=100000;


				for (int k=0; k<n_nTuple; k++){



					model_pT=gRandom->Uniform(model_pTMin,model_pTMax);
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
					if(model_pT<pTMin[rapIndex] || model_pT>pTMax[rapIndex]) FillEvent=false;

					if(FillEvent){



						model_weight=SDC_Graph[StatesGivenIndex][rapIndex][iColorChannel]->Eval(model_pT);
						//if(log(model_weight)>8 && TMath::Abs(model_rap)<1.2&&model_pT>12&&model_pT<13.5) cout<<"model_weight pT"<<model_weight<<endl;
						//Define angular distribution:
						//cout<<"model_weight pT"<<model_weight<<endl;

						//dvector model_lam = func_lam_gen(iMother, iColorChannel);
						double model_lamth, model_lamph, model_lamtp;
						model_lamth=Lamth_Graph[StatesGivenIndex][rapIndex][iColorChannel]->Eval(model_pT);
						model_lamph=Lamph_Graph[StatesGivenIndex][rapIndex][iColorChannel]->Eval(model_pT);
						model_lamtp=Lamtp_Graph[StatesGivenIndex][rapIndex][iColorChannel]->Eval(model_pT);

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
								cout<<"model_weight pT "<<SDC_Graph[StatesGivenIndex][rapIndex][iColorChannel]->Eval(model_pT)<<endl;
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

						model_pT*=MassRatioCtB[ScaleStateIndex];

						nTupleModel[iScaleToMother][iColorChannel]->Fill();

						if (k%nStep == 0) {
							cout << nStep_*10 <<"% (nState = "<<iScaleToMother<<", Color channel = "<<iColorChannel<<")"<<endl;
							++nStep_;
						}

					}
					else k--;




				}


				bool interpretModelIntegratedInRapidityInterval=true;

				double deltapT=(model_pTMax-model_pTMin)*MassRatioCtB[ScaleStateIndex];
				double deltay=model_rapMax-model_rapMin;

				double phasespaceFactor=0;//30.*1.2+20.*1.2; //sum_i ( deltaPt_i*deltaRap_i ) with i...rap bin
				for(int j=0;j<nRapIntervals;j++){
					double deltaPt=(pTMax[j]-pTMin[j])*MassRatioCtB[ScaleStateIndex];
					double deltaY=RapIntervalBordersMax[j]-RapIntervalBordersMin[j];
					if(isAbsRap[j]) deltaY*=2;

					if(interpretModelIntegratedInRapidityInterval) phasespaceFactor+=deltaPt;
					else phasespaceFactor+=deltaPt*deltaY;
				}

				double globalWeight=1./double(n_nTuple)*phasespaceFactor;
				nTupleModel[iScaleToMother][iColorChannel]->SetWeight(globalWeight);

				cout<<"globalWeight "<<globalWeight<<endl;

				cout<<"Generated "<<n_nTuple<<" events"<<endl;
				cout<<"nTupleModel[iScaleToMother][iColorChannel]->GetEntries() "<<nTupleModel[iScaleToMother][iColorChannel]->GetEntries()<<endl;


				nTupleModel[iScaleToMother][iColorChannel]->Write();

			}
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

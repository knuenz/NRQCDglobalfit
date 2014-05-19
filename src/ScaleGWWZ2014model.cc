/*
 * ScaleGWWZ2014model.cc
 *
 *  Created on: May 16, 2014
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
#include "TCanvas.h"
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
double exponential(double* x, double* par);
Double_t paramMassRapParabola(Double_t *x, Double_t *par);

int main(int argc, char** argv) {

  	Char_t *OriginalModelID = "Default";
	Char_t *storagedir = "Default"; //Storage Directory
	bool useTF1inputs=false;

  	for( int i=0;i < argc; ++i ) {
		if(std::string(argv[i]).find("OriginalModelID") != std::string::npos) {char* OriginalModelIDchar = argv[i]; char* OriginalModelIDchar2 = strtok (OriginalModelIDchar, "="); OriginalModelID = OriginalModelIDchar2; cout<<"OriginalModelID = "<<OriginalModelID<<endl;}
		if(std::string(argv[i]).find("storagedir") != std::string::npos) {char* storagedirchar = argv[i]; char* storagedirchar2 = strtok (storagedirchar, "="); storagedir = storagedirchar2; cout<<"storagedir = "<<storagedir<<endl;}
	    if(std::string(argv[i]).find("useTF1inputs=true") != std::string::npos) {useTF1inputs=true;cout<<"useTF1inputs=true"<<endl;}
  	}

	gRandom = new TRandom3(NRQCDvars::randomSeed);  // better random generator

	char outname[2000];
	char inname[2000];
	char premodeldirname[2000];
	char originalmodeldirname[2000];

	sprintf(premodeldirname,"%s/OriginalModelID",storagedir);
	gSystem->mkdir(premodeldirname);
	sprintf(originalmodeldirname,"%s/%s",premodeldirname,OriginalModelID);
	gSystem->mkdir(originalmodeldirname);


	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Convert GWWZ2014 input model
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//STEP1: read in the original inputs (lp and tp SDC components)

		const int nStatesGiven=7;//Given models in any form
		int StatesGiven[nStatesGiven]={0, 1, 2, 3, 4, 7, 10};
		const int nColorChannelsGiven=4;

		double fittedLDMEs[nStatesGiven][nColorChannelsGiven]={//if inputs are given as SDC*LDME, the actual LDMEs have to be put here. If the inputs are given as SDC already, the LDME correction is 1
				{1.16, 1., 1., 1.},
				{0.1075, 1., 1., 1.},
				{0.1075, 1., 1., 1.},
				{0.76, 1., 1., 1.},
				{1., 1., 1., 1.},
				{1., 1., 1., 1.},
				{1., 1., 1., 1.}
		};


		const int nHelicityChannels=2;	//00, 11, 1m1

		const int nRapIntervals=3;
		bool isAbsRap[nRapIntervals]={true, true, true};
		double RapIntervalBordersMin[nRapIntervals]={0, 0.6, 1.2};
		double RapIntervalBordersMax[nRapIntervals]={0.6, 1.2, 1.5};
		double AverageRap[nRapIntervals]={0.3, 0.9, 1.35};


		double deltaRapOriginal[nRapIntervals];
		for(int iRap=0;iRap<nRapIntervals;iRap++){
			deltaRapOriginal[iRap]=RapIntervalBordersMax[iRap]-RapIntervalBordersMin[iRap];
			if(isAbsRap[iRap]) deltaRapOriginal[iRap]*=2;
		}

		bool interpretOriginalModelAsIntegratedInRap=true;

		bool isDummyRap[nStatesGiven][nRapIntervals]={
				{false, false, false},
				{false, false, false},
				{false, false, false},
				{false, false, false},
				{false, false, true},
				{false, false, true},
				{false, false, true}
		};
		bool isDummyColorChannel[nStatesGiven][nColorChannelsGiven]={
				{false, false, false, false},
				{false, true, false, true},
				{false, true, false, true},
				{false, false, false, false},
				{false, false, false, false},
				{false, false, false, false},
				{false, false, false, false},
		};
		double dummyVal=999;


		const int maxpTBinsPerRap=40;

		double pTmean[nStatesGiven][nRapIntervals][nColorChannelsGiven][nHelicityChannels][maxpTBinsPerRap];//={

		double pTMin[nStatesGiven][nRapIntervals]={
				{3, 3, 3},
				{3, 3, 3},
				{3, 3, 3},
				{3, 3, 3},
				{3, 3, 3},
				{3, 3, 3},
				{3, 3, 3}
		};
		double pTMax[nStatesGiven][nRapIntervals]={
				{75, 75, 75},
				{75, 75, 75},
				{75, 75, 75},
				{75, 75, 75},
				{100, 100, 100},
				{110, 110, 110},
				{110, 110, 110}
		};

		cout<<"read in file"<<endl;

		//Objects used for TF1 input
		TF1* SDC_funcInp[nStatesGiven][nRapIntervals][nColorChannelsGiven][nHelicityChannels];
		//TGraph *SDC_graphInp[nStatesGiven][nRapIntervals][nColorChannelsGiven][nHelicityChannels];

		//Objects used for txt input
		TGraph *model_Graph[nStatesGiven][nRapIntervals][nColorChannelsGiven][nHelicityChannels];


		if(!useTF1inputs){

			const int npTBinsPerRap[nStatesGiven][nRapIntervals]={
					{24, 24, 24},
					{24, 24, 24},
					{24, 24, 24},
					{24, 24, 24},
					{40, 40, 40},
					{32, 32, 32},
					{32, 32, 32}
			};

			double model[nStatesGiven][nRapIntervals][nColorChannelsGiven][nHelicityChannels][maxpTBinsPerRap];
			double err_model[nStatesGiven][nRapIntervals][nColorChannelsGiven][nHelicityChannels][maxpTBinsPerRap];

			char in_rapChar[200];
			char in_stateChar[200];
			char in_helcolChar[200];

			for(int i=0;i<nStatesGiven;i++){
				cout<<"State"<<StatesGiven[i]<<endl;

				if(StatesGiven[i]==0) sprintf(in_stateChar,"Jpsi/direct");
				if(StatesGiven[i]==1) sprintf(in_stateChar,"Jpsi/chic_feeddown");
				if(StatesGiven[i]==2) sprintf(in_stateChar,"Jpsi/chic_feeddown");
				if(StatesGiven[i]==3) sprintf(in_stateChar,"Psi2s/direct");
				if(StatesGiven[i]==4) sprintf(in_stateChar,"upsilon1s/direct");
				if(StatesGiven[i]==7) sprintf(in_stateChar,"upsilon2s/direct");
				if(StatesGiven[i]==10) sprintf(in_stateChar,"upsilon3s/direct");

				for(int j=0;j<nRapIntervals;j++){
					cout<<"rap"<<j<<endl;

					if(j==0) sprintf(in_rapChar,"y_0-0.6");
					if(j==1) sprintf(in_rapChar,"y_0.6-1.2");
					if(j==2) sprintf(in_rapChar,"y_1.2-1.5");

					for(int l=0;l<nHelicityChannels;l++){
						cout<<"helicity"<<l<<endl;


						for(int k=0;k<nColorChannelsGiven;k++){
							cout<<"color channel"<<k<<endl;

							if(l==0 && k==0) sprintf(in_helcolChar,"lp/cs");
							if(l==1 && k==0) sprintf(in_helcolChar,"tp/cs");

							if(l==0 && k==1) sprintf(in_helcolChar,"pt/1s08");
							if(l==1 && k==1) sprintf(in_helcolChar,"pt/1s08");

							if(l==0 && k==2) sprintf(in_helcolChar,"lp/3s18");
							if(l==1 && k==2) sprintf(in_helcolChar,"tp/3s18");

							if(l==0 && k==3) sprintf(in_helcolChar,"lp/3pj8");
							if(l==1 && k==3) sprintf(in_helcolChar,"tp/3pj8");

							if(StatesGiven[i]==1 || StatesGiven[i]==2){
								if(l==0 && k==0) sprintf(in_helcolChar,"lp/3pj1");
								if(l==1 && k==0) sprintf(in_helcolChar,"tp/3pj1");
							}

							if(StatesGiven[i]<4) sprintf(inname,"%s/data-jpsi/%s/%s/%s/data",originalmodeldirname,in_rapChar,in_stateChar,in_helcolChar);
							else sprintf(inname,"%s/data-ups/%s/%s/%s/data",originalmodeldirname,in_rapChar,in_stateChar,in_helcolChar);

							FILE *fIn;
							fIn = fopen(inname, "read");
							cout<<inname<<endl;
							Char_t line[1000];

							bool readIn=true;
							if(isDummyColorChannel[i][k]){
								cout<<"isDummyColorChannel"<<endl;
								readIn=false;
							}
							if(isDummyRap[i][j]){
								cout<<"isDummyRap"<<endl;
								readIn=false;
							}

							for(int m=0;m<npTBinsPerRap[i][j];m++){


								if(!readIn){
									pTmean[i][j][k][l][m]=m;
									model[i][j][k][l][m]=dummyVal;
									err_model[i][j][k][l][m]=dummyVal;
								}

								if(readIn){
									double buffer;
									fgets(line, sizeof(line), fIn); //comment
									cout<<line<<endl;
									sscanf(line, "%lf %lf %lf", &pTmean[i][j][k][l][m], &model[i][j][k][l][m], &err_model[i][j][k][l][m]);
								}

								if(k==1){
									model[i][j][k][l][m]/=3.;
									err_model[i][j][k][l][m]/=3.;
								}

								cout<<pTmean[i][j][k][l][m]<<" "<<model[i][j][k][l][m]<<" "<<err_model[i][j][k][l][m]<<" "<<endl;

							}
							if(readIn) fgets(line, sizeof(line), fIn); //comment
						}
					}
				}
			}


			for(int i=0;i<nStatesGiven;i++){
				for(int j=0;j<nRapIntervals;j++){
					for(int k=0;k<nColorChannelsGiven;k++){
						for(int l=0;l<nHelicityChannels;l++){

							double pTmean_graph[npTBinsPerRap[i][j]];
							double model_graph[npTBinsPerRap[i][j]];
							for(int m=0;m<npTBinsPerRap[i][j];m++){
								pTmean_graph[m]=pTmean[i][j][k][l][m];
								model_graph[m]=model[i][j][k][l][m];
							}

							model_Graph[i][j][k][l] = new TGraph(npTBinsPerRap[i][j],pTmean_graph,model_graph);

							cout<<"model_Graph[state"<<i<<"][rap"<<j<<"][CC"<<k<<"][HC"<<l<<"]"<<endl;
							model_Graph[i][j][k][l]->Print();

						}
					}
				}
			}


		}



		else if(useTF1inputs){


			char in_rapChar[200];
			char in_stateChar[200];
			char in_helcolChar[200];
			char in_colorChar[200];
			char in_name_func[200];
			char in_name_graph[200];
			char in_name_part[200];

			for(int i=0;i<nStatesGiven;i++){
				cout<<"State"<<StatesGiven[i]<<endl;

				if(StatesGiven[i]==0) sprintf(in_stateChar,"JPsi");
				if(StatesGiven[i]==1) sprintf(in_stateChar,"JPsi");
				if(StatesGiven[i]==2) sprintf(in_stateChar,"JPsi");
				if(StatesGiven[i]==3) sprintf(in_stateChar,"JPsi");
				if(StatesGiven[i]==4) sprintf(in_stateChar,"Ups1S");
				if(StatesGiven[i]==7) sprintf(in_stateChar,"Ups2S");
				if(StatesGiven[i]==10) sprintf(in_stateChar,"Ups3S");

				sprintf(inname,"%s/data-TF1/param_GWWZ_longTransvComp_direct%s_pT.root",originalmodeldirname,in_stateChar);

				TFile *infile = new TFile(inname,"READ");
				cout<<inname<<endl;

				for(int j=0;j<nRapIntervals;j++){
					cout<<"rap"<<j<<endl;

					for(int l=0;l<nHelicityChannels;l++){
						cout<<"helicity"<<l<<endl;

						if(l==0) sprintf(in_helcolChar,"LP");
						if(l==1) sprintf(in_helcolChar,"TP");

						for(int k=0;k<nColorChannelsGiven;k++){
							cout<<"color channel"<<k<<endl;

							//set the LDMEs to 1, given that the TF1 inputs do not have to be corrected by the LDMEs
							fittedLDMEs[i][k]=1.;

							if(k==0) sprintf(in_colorChar,"cs");
							if(k==1) sprintf(in_colorChar,"1s08");
							if(k==2) sprintf(in_colorChar,"3s18");
							if(k==3) sprintf(in_colorChar,"3pj8");

							sprintf(in_name_part,"%s_%s_%s_rap%d",in_stateChar, in_colorChar, in_helcolChar, j);

							sprintf(in_name_func,"func_%s",in_name_part);
							sprintf(in_name_graph,"graph_%s",in_name_part);

							bool readIn=true;
							if(isDummyColorChannel[i][k]){
								cout<<"isDummyColorChannel"<<endl;
								readIn=false;
							}
							if(isDummyRap[i][j]){
								cout<<"isDummyRap"<<endl;
								readIn=false;
							}

							if(readIn){
								cout<<"readIn "<<in_name_graph<<endl;
								SDC_funcInp[i][j][k][l] = (TF1*)infile->Get(in_name_func);
								model_Graph[i][j][k][l] = (TGraph*)infile->Get(in_name_graph);
							}
							else{
								cout<<"!readIn"<<endl;
								SDC_funcInp[i][j][k][l] = new TF1("gaus", "gaus", 0., 1.);
								const int nBinsBuff=1000;
								double buffX[nBinsBuff] = {NULL};
								model_Graph[i][j][k][l] = new TGraph(nBinsBuff,buffX,buffX);
							}

							cout<<"model_Graph[i][j][k][l]->GetN() "<<model_Graph[i][j][k][l]->GetN()<<endl;

						}
					}
				}

				infile->Close();
				delete infile;
				infile = NULL;

			}


		}



		//STEP2: convert lp and tp SDC components to SDCs and polarizations

		TGraph *SDC_Graph[nStatesGiven][nRapIntervals][nColorChannelsGiven];
		TGraph *Lamth_Graph[nStatesGiven][nRapIntervals][nColorChannelsGiven];
		TGraph *Lamph_Graph[nStatesGiven][nRapIntervals][nColorChannelsGiven];
		TGraph *Lamtp_Graph[nStatesGiven][nRapIntervals][nColorChannelsGiven];

		TGraph *SDC_Graph_original[nStatesGiven][nRapIntervals][nColorChannelsGiven];
		TGraph *Lamth_Graph_original[nStatesGiven][nRapIntervals][nColorChannelsGiven];
		TGraph *Lamph_Graph_original[nStatesGiven][nRapIntervals][nColorChannelsGiven];
		TGraph *Lamtp_Graph_original[nStatesGiven][nRapIntervals][nColorChannelsGiven];

	bool UseChaoPolDenominatorDef=false;
	char graphName[200];

	int nBinsFinalModels=500;
	for(int i=0;i<nStatesGiven;i++){
		for(int j=0;j<nRapIntervals;j++){
			for(int k=0;k<nColorChannelsGiven;k++){




				cout<<"originals: rap "<<j<<" CC "<<k<<" state "<<i<<endl;

				int nBinsOriginal=model_Graph[i][j][k][0]->GetN();
				cout<<"nBinsOriginal "<<nBinsOriginal<<endl;

				double pTmean_graph_original[nBinsOriginal];
				double SDC_graph_original[nBinsOriginal];
				double Lamth_graph_original[nBinsOriginal];
				double Lamph_graph_original[nBinsOriginal];
				double Lamtp_graph_original[nBinsOriginal];
				for(int m=0;m<nBinsOriginal;m++){
					double buff_sigma00;
					double buff_sigma11;
					double buff_sigma1m1;
					model_Graph[i][j][k][0]->GetPoint(m,pTmean_graph_original[m],buff_sigma00);
					model_Graph[i][j][k][1]->GetPoint(m,pTmean_graph_original[m],buff_sigma11);
					//model_Graph[i][j][k][2]->GetPoint(m,pTmean_graph_original[m],buff_sigma1m1);
					//if(pTmean_graph_original[m]<45) model_Graph[i][j][k][2]->GetPoint(m,pTmean_graph_original[m],buff_sigma1m1);
					//else buff_sigma1m1=0.;

					SDC_graph_original[m]=buff_sigma00/fittedLDMEs[i][k]+2*buff_sigma11/fittedLDMEs[i][k];
					double polDenominator=(buff_sigma11/fittedLDMEs[i][k]+buff_sigma00/fittedLDMEs[i][k]);
					if(UseChaoPolDenominatorDef && buff_sigma11 * buff_sigma00 < 0) polDenominator=fabs(polDenominator);
					Lamth_graph_original[m]=(buff_sigma11/fittedLDMEs[i][k]-buff_sigma00/fittedLDMEs[i][k])/polDenominator;
					Lamph_graph_original[m]=0.;//(buff_sigma1m1/fittedLDMEs[i][k])/polDenominator;
					Lamtp_graph_original[m]=0;

					if(isDummyRap[i][j]){
						SDC_graph_original[m]=dummyVal;
						Lamth_graph_original[m]=dummyVal;
						Lamph_graph_original[m]=dummyVal;
						Lamtp_graph_original[m]=dummyVal;
					}
				}





				SDC_Graph_original[i][j][k]= new TGraph(nBinsOriginal,pTmean_graph_original,SDC_graph_original);
				Lamth_Graph_original[i][j][k]= new TGraph(nBinsOriginal,pTmean_graph_original,Lamth_graph_original);
				Lamph_Graph_original[i][j][k]= new TGraph(nBinsOriginal,pTmean_graph_original,Lamph_graph_original);
				Lamtp_Graph_original[i][j][k]= new TGraph(nBinsOriginal,pTmean_graph_original,Lamtp_graph_original);

				sprintf(graphName,"SDC_Graph_original_state%d_rap%d_CC%d",StatesGiven[i],j,k);
				SDC_Graph_original[i][j][k]->SetName(graphName);
				sprintf(graphName,"Lamth_Graph_original_state%d_rap%d_CC%d",StatesGiven[i],j,k);
				Lamth_Graph_original[i][j][k]->SetName(graphName);
				sprintf(graphName,"Lamph_Graph_original_state%d_rap%d_CC%d",StatesGiven[i],j,k);
				Lamph_Graph_original[i][j][k]->SetName(graphName);
				sprintf(graphName,"Lamtp_Graph_original_state%d_rap%d_CC%d",StatesGiven[i],j,k);
				Lamtp_Graph_original[i][j][k]->SetName(graphName);




				cout<<"fine-binned: rap "<<j<<" CC "<<k<<" state "<<i<<endl;
				double pTmean_graph[nBinsFinalModels];
				double SDC_graph[nBinsFinalModels];
				double Lamth_graph[nBinsFinalModels];
				double Lamph_graph[nBinsFinalModels];
				double Lamtp_graph[nBinsFinalModels];
				for(int m=0;m<nBinsFinalModels;m++){
					double deltaPt=(pTMax[i][j]-pTMin[i][j])/double(nBinsFinalModels);
					pTmean_graph[m]=pTMin[i][j]+m*deltaPt;
					//SDC_graph[m]=SDC_Graph_original[i][j][k]->Eval(pTmean_graph[m],0,"S");

					double rapCorrFactor=deltaRapOriginal[j];
					if(!interpretOriginalModelAsIntegratedInRap) rapCorrFactor=1.;

					SDC_graph[m]=SDC_Graph_original[i][j][k]->Eval(pTmean_graph[m])/rapCorrFactor;
					Lamth_graph[m]=Lamth_Graph_original[i][j][k]->Eval(pTmean_graph[m]);
					Lamph_graph[m]=0;//Lamph_Graph_original[i][j][k]->Eval(pTmean_graph[m]);
					Lamtp_graph[m]=0;

					if(useTF1inputs){
						SDC_graph[m]=(SDC_funcInp[i][j][k][0]->Eval(pTmean_graph[m])+2*SDC_funcInp[i][j][k][1]->Eval(pTmean_graph[m]))/rapCorrFactor;
						Lamth_graph[m]=(SDC_funcInp[i][j][k][1]->Eval(pTmean_graph[m])-SDC_funcInp[i][j][k][0]->Eval(pTmean_graph[m])) / (SDC_funcInp[i][j][k][1]->Eval(pTmean_graph[m])+SDC_funcInp[i][j][k][0]->Eval(pTmean_graph[m]));
						Lamph_graph[m]=0;
						Lamtp_graph[m]=0;

					}

				}

				SDC_Graph[i][j][k]= new TGraph(nBinsFinalModels,pTmean_graph,SDC_graph);
				Lamth_Graph[i][j][k]= new TGraph(nBinsFinalModels,pTmean_graph,Lamth_graph);
				Lamph_Graph[i][j][k]= new TGraph(nBinsFinalModels,pTmean_graph,Lamph_graph);
				Lamtp_Graph[i][j][k]= new TGraph(nBinsFinalModels,pTmean_graph,Lamtp_graph);

				sprintf(graphName,"SDC_Graph_state%d_rap%d_CC%d",StatesGiven[i],j,k);
				SDC_Graph[i][j][k]->SetName(graphName);
				sprintf(graphName,"Lamth_Graph_state%d_rap%d_CC%d",StatesGiven[i],j,k);
				Lamth_Graph[i][j][k]->SetName(graphName);
				sprintf(graphName,"Lamph_Graph_state%d_rap%d_CC%d",StatesGiven[i],j,k);
				Lamph_Graph[i][j][k]->SetName(graphName);
				sprintf(graphName,"Lamtp_Graph_state%d_rap%d_CC%d",StatesGiven[i],j,k);
				Lamtp_Graph[i][j][k]->SetName(graphName);

				//cout<<"SDC_Graph[i][j][k]"<<endl;
				//SDC_Graph[i][j][k]->Print();
				//cout<<"Lamth_Graph[i][j][k]"<<endl;
				//Lamth_Graph[i][j][k]->Print();
				//cout<<"Lamph_Graph[i][j][k]"<<endl;
				//Lamph_Graph[i][j][k]->Print();
				//cout<<"Lamtp_Graph[i][j][k]"<<endl;
				//Lamtp_Graph[i][j][k]->Print();
			}
		}
	}




	//STEP3: Extrapolate SDCs and polarizations to desired pT and y region


	const int nEpRapIntervals=9;
	//bool isAbsEpRap[nEpRapIntervals]={true, true, true, true, true, true, true, true, true};
	double EpRapIntervalBordersMin[nEpRapIntervals]={0, 0.6, 1.2, 1.5, 2., 2.5, 3., 3.5, 4.};
	double EpRapIntervalBordersMax[nEpRapIntervals]={0.6, 1.2, 1.5, 2., 2.5, 3., 3.5, 4., 4.5};
	double AverageEpRap[nEpRapIntervals]={0.3, 0.9, 1.35, 1.75, 2.25, 2.75, 3.25, 3.75, 4.25};

	double pTMinEp[nStatesGiven][nEpRapIntervals]={
			{2.875, 2.875, 2.875, 2.875, 2.875, 2.875, 2.875, 2.875, 2.875},
			{2.875, 2.875, 2.875, 2.875, 2.875, 2.875, 2.875, 2.875, 2.875},
			{2.875, 2.875, 2.875, 2.875, 2.875, 2.875, 2.875, 2.875, 2.875},
			{2.875, 2.875, 2.875, 2.875, 2.875, 2.875, 2.875, 2.875, 2.875},
			{2.875, 2.875, 2.875, 2.875, 2.875, 2.875, 2.875, 2.875, 2.875},
			{2.875, 2.875, 2.875, 2.875, 2.875, 2.875, 2.875, 2.875, 2.875},
			{2.875, 2.875, 2.875, 2.875, 2.875, 2.875, 2.875, 2.875, 2.875}
	};
	double pTMaxEp[nStatesGiven][nEpRapIntervals]={
			{80., 80., 80., 80., 80., 80., 80., 80., 80.},
			{80., 80., 80., 80., 80., 80., 80., 80., 80.},
			{80., 80., 80., 80., 80., 80., 80., 80., 80.},
			{80., 80., 80., 80., 80., 80., 80., 80., 80.},
			{110., 110., 110., 110., 110., 110., 110., 110., 110.},
			{110., 110., 110., 110., 110., 110., 110., 110., 110.},
			{110., 110., 110., 110., 110., 110., 110., 110., 110.}
	};

	bool useOriginalRap[nStatesGiven][nEpRapIntervals]={
			{true, true, true, false, false, false, false, false, false},
			{true, true, true, false, false, false, false, false, false},
			{true, true, false, false, false, false, false, false, false},
			{true, true, false, false, false, false, false, false, false},
			{true, true, false, false, false, false, false, false, false}
	};

	double deltaRap[nEpRapIntervals];
	double NormRapXToRap0Factor[nEpRapIntervals][nColorChannelsGiven];
	double NormRapXToRap0Factor_dycorr[nEpRapIntervals][nColorChannelsGiven];
	int normToRap=0;
	double normAtPt=10.;

	double LamthDiffToRap2[nEpRapIntervals][nColorChannelsGiven];
	int diffToRap=2;
	double diffAtPt=10.;

	for(int i=0;i<1;i++){
		for(int j=0;j<nRapIntervals;j++){

			deltaRap[j]=RapIntervalBordersMax[j]-RapIntervalBordersMin[j];
			if(isAbsRap[j]) deltaRap[j]*=2;

			for(int k=0;k<nColorChannelsGiven;k++){

				NormRapXToRap0Factor[j][k]=SDC_Graph[i][j][k]->Eval(normAtPt)/SDC_Graph[i][normToRap][k]->Eval(normAtPt);

				LamthDiffToRap2[j][k]=Lamth_Graph[i][j][k]->Eval(diffAtPt)-Lamth_Graph[i][diffToRap][k]->Eval(diffAtPt);

				cout<<"LamthDiffToRap2["<<j<<"]["<<k<<"] = "<<LamthDiffToRap2[j][k]<<endl;

				cout<<"NormRapXToRap0Factor["<<j<<"]["<<k<<"] = "<<NormRapXToRap0Factor[j][k]<<endl;

			}
		}
	}


	//Get rapidity-dependence of normalization


	double rapDepNorm[nRapIntervals];

	TGraph *g_rapDepNorm[nColorChannelsGiven];
	char name[200];
	sprintf(name, "fParabola");
	TF1* fParabola[nColorChannelsGiven];
	int color[nColorChannelsGiven]={kGray, kGreen+2, kRed, 1};

	TCanvas *tmpC = new TCanvas("tmpC","tmpC",1000,800);
	tmpC->SetFillColor(kWhite);
	tmpC->SetFrameBorderMode(0);

	TGraph *g_rapDepNorm_XC;
	const int nRapIntervals_XC=9;
	double rapDepNorm_XC[nRapIntervals_XC]={4.35e-1/4.31780303030303014e-01, 4.31e-1/4.31780303030303014e-01, 3.84e-1/4.31780303030303014e-01, 3.49e-1/4.31780303030303014e-01, 2.94e-1/4.31780303030303014e-01, 2.56e-1/4.31780303030303014e-01, 2.17e-1/4.31780303030303014e-01, 1.28e-1/4.31780303030303014e-01, 0.715e-1/4.31780303030303014e-01};
	double rap_XC[nRapIntervals_XC]={0.375, 0.6, 1.4, 1.875, 2.25, 2.75, 3.25, 3.75, 4.25};
	g_rapDepNorm_XC = new TGraph(nRapIntervals_XC, rap_XC, rapDepNorm_XC);
	g_rapDepNorm_XC->SetMarkerColor(kMagenta);
	g_rapDepNorm_XC->SetMarkerStyle(24);

	for(int i=0;i<1;i++){
		for(int k=0;k<nColorChannelsGiven;k++){
			for(int j=0;j<nRapIntervals;j++){
				rapDepNorm[j]=SDC_Graph[i][j][k]->Eval(normAtPt)/SDC_Graph[i][normToRap][k]->Eval(normAtPt);;
				if(j==2&&k==0) rapDepNorm[j]+=0.035;
				if(j==2&&k==2) rapDepNorm[j]+=0.0375;
			}
			g_rapDepNorm[k] = new TGraph(nRapIntervals, AverageRap, rapDepNorm);
			g_rapDepNorm[k]->Print();
			fParabola[k] = new TF1(name, paramMassRapParabola, 0, 5, 3.);

			double a=1;
			double b=1;
			double c=1;
			fParabola[k]->SetParameter(0,a);
			fParabola[k]->SetParameter(1,b);
			fParabola[k]->SetParameter(2,c);
			g_rapDepNorm[k]->Fit(fParabola[k], "0", "", 0., 5.);

			double a_fit=fParabola[k]->GetParameter(0);
			double b_fit=fParabola[k]->GetParameter(1);
			double c_fit=fParabola[k]->GetParameter(2);

			cout<<"double a_fit = "<<a_fit<<";"<<endl;
			cout<<"double b_fit = "<<b_fit<<";"<<endl;
			cout<<"double c_fit = "<<c_fit<<";"<<endl;

			fParabola[k]->SetLineColor(color[k]);
			fParabola[k]->SetLineWidth(1.);
			g_rapDepNorm[k]->SetMarkerStyle(20);
			g_rapDepNorm[k]->SetTitle(0);
			g_rapDepNorm[k]->GetYaxis()->SetRangeUser(0,1.2);
			g_rapDepNorm[k]->GetXaxis()->SetLimits(0.,5.);
			g_rapDepNorm[k]->GetXaxis()->SetTitle("#it{y}");
			g_rapDepNorm[k]->GetYaxis()->SetTitle("normalized to #bar{y}=0.3");

			if(k==0) g_rapDepNorm[k]->Draw("ap");
			else g_rapDepNorm[k]->Draw("psame");
			fParabola[k]->Draw("lsame");

		}
	}

	g_rapDepNorm_XC->Draw("psame");

	char tmpN[500];
	sprintf(tmpN,"%s/RapNormBK.pdf",originalmodeldirname);
	tmpC->SaveAs(tmpN);

	// Fit with some function






	TGraph *EpSDC_Graph[nStatesGiven][nEpRapIntervals][nColorChannelsGiven];
	TGraph *EpLamth_Graph[nStatesGiven][nEpRapIntervals][nColorChannelsGiven];
	TGraph *EpLamph_Graph[nStatesGiven][nEpRapIntervals][nColorChannelsGiven];
	TGraph *EpLamtp_Graph[nStatesGiven][nEpRapIntervals][nColorChannelsGiven];


	for(int i=0;i<nStatesGiven;i++){
		for(int j=0;j<nEpRapIntervals;j++){
			for(int k=0;k<nColorChannelsGiven;k++){

				cout<<"Ep: rap "<<j<<" CC "<<k<<" state "<<i<<endl;

				cout<<"NormRapXToRap0Factor["<<j<<"]["<<k<<"] = "<<NormRapXToRap0Factor[j][k]<<endl;

				double pTmean_graph_Ep[nBinsFinalModels];
				double SDC_graph_Ep[nBinsFinalModels];
				double Lamth_graph_Ep[nBinsFinalModels];
				double Lamph_graph_Ep[nBinsFinalModels];
				double Lamtp_graph_Ep[nBinsFinalModels];
				for(int m=0;m<nBinsFinalModels;m++){

					double deltaPt=(pTMaxEp[i][j]-pTMinEp[i][j])/double(nBinsFinalModels);
					pTmean_graph_Ep[m]=pTMinEp[i][j]+m*deltaPt;


					if(useOriginalRap[i][j]){
						SDC_graph_Ep[m]=SDC_Graph[i][j][k]->Eval(pTmean_graph_Ep[m]);
					}
					else{
						SDC_graph_Ep[m]=SDC_Graph[i][0][k]->Eval(pTmean_graph_Ep[m])*fParabola[k]->Eval(AverageEpRap[j]);

					}

					if(k==1){
						Lamth_graph_Ep[m]=0;
					}
					else{

						if(useOriginalRap[i][j]){
							Lamth_graph_Ep[m]=Lamth_Graph[i][j][k]->Eval(pTmean_graph_Ep[m]);
						}
						else{
							for(int jj=j;jj>-1;jj--){
								if(useOriginalRap[i][jj]){
									Lamth_graph_Ep[m]=Lamth_Graph[i][jj][k]->Eval(pTmean_graph_Ep[m]);
									break;
								}
							}
						}
					}




					Lamph_graph_Ep[m]=0;
					Lamtp_graph_Ep[m]=0;

				}

				EpSDC_Graph[i][j][k]= new TGraph(nBinsFinalModels,pTmean_graph_Ep,SDC_graph_Ep);
				EpLamth_Graph[i][j][k]= new TGraph(nBinsFinalModels,pTmean_graph_Ep,Lamth_graph_Ep);
				EpLamph_Graph[i][j][k]= new TGraph(nBinsFinalModels,pTmean_graph_Ep,Lamph_graph_Ep);
				EpLamtp_Graph[i][j][k]= new TGraph(nBinsFinalModels,pTmean_graph_Ep,Lamtp_graph_Ep);

				sprintf(graphName,"EpSDC_Graph_state%d_rap%d_CC%d",StatesGiven[i],j,k);
				EpSDC_Graph[i][j][k]->SetName(graphName);
				sprintf(graphName,"EpLamth_Graph_state%d_rap%d_CC%d",StatesGiven[i],j,k);
				EpLamth_Graph[i][j][k]->SetName(graphName);
				sprintf(graphName,"EpLamph_Graph_state%d_rap%d_CC%d",StatesGiven[i],j,k);
				EpLamph_Graph[i][j][k]->SetName(graphName);
				sprintf(graphName,"EpLamtp_Graph_state%d_rap%d_CC%d",StatesGiven[i],j,k);
				EpLamtp_Graph[i][j][k]->SetName(graphName);

			}
		}
	}





	/*

 1) scale the "horizontal" variable, pT, such that p2/p1 = m2/m1 (i.e. taking into account the average rapidity to calculate p for each pT point)

2) make the ratio of psiprime over "scaled Jpsi" (for example, or the inverse): check if this ratio is more or less flat, which means that the "shape scaling" was correct, so that we can afterwards apply a "normalization scaling". I would even use this test to see if the scaling p2/p1 = m2/m1 is better or worse that the pT2/pT1 = m2/m1. My idea that the good one is the one with p is based on the analogy with the feeddowns, where I have explicitly calculated that p is the correct variable. But clearly here it is not the same thing. Our aim is to *find* a scaling rule, even if it is not the most obvious one. So whatever works best, we should use it.

3) If the ratio is reasonably flat, take its constant value as the
multiplicative factor defining the "normalization scaling".

4) Now we have a general rule to scale from direct Jspi to any mass: first the shape scaling (p2/p1 = m2/m1 or pT2/pT1 = m2/m1), then a normalization scaling equal to the "maximum" scaling determined from the psiprime/Jpsi_scaled ratio, multiplied by
(m_X - m_Jspi)/(m_psiprime - m_Jspi) (or something similar: is this the linear interpolation? you know better, after doing a lot of interpolations of sidebands).

5) What if the shape scaling does not work? Not very good, but it is true that whatever scaling we may found to work and that is general enough would be good, even if it is not factorized into "shape" and "normalization". We should think abut other options. In the worst case, we can determine the "maximum scaling" simply as a pT-dependent scaling equal to the ratio psiprime/Jspi, possibly as a function of pT/M (not of pT, otherwise we know that the scaling would not work for sure for the Upsilons).

6) As you say, the normalization scaling is less crucial than the shape scaling. Physically we do not care too much about this. Our considerations (dominance of this or that) will always be based on the product SDC * LDME. So I think that the linear interpolation for the normalization part is not a catastrophe. On the other hand it is better than a pure arbitrary scaling. It would be good that the LDMEs resulting from the fit were not completely crazy at least in order of magnitude.

7) Concerning the colour singlet part, I think that the problem does not really exist for the charmonium case, where CS is negligible. We can do the intentional approximation to neglect it, also in the Upsilon case, where we will use curves *extra*polated from the charmonium case. Actually, concerning the S states we can simply determine the CS scaling just as done for the colour octet states, and apply it to scale the CS contributions from charmonium to bottomonium (not obvious that the normalizations will be correct, but we may try, and compare with the literature: I think that the CS Upsilon curves can be found, even in NRQCD papers). The only thing that we cannot do is to interpolate (charmonium) or extrapolate (bottomonium) the CS SDCs to determine those of the chi's, because of the different quantum numbers. But maybe we could do this anyway, assuming that PDF and mass dependence count more than differences in the participating Feynman processes. Provided that we state the approximation and what we obtain makes sense, we can do it.

:::rather than a linear interpolation, an exponential one would be better and perhaps more physical,

NEW:

1) scale all with pT'=pT*M/3
-> Fit psi' and Ups3S

2) normalization scaling
-Make TGraphs with many bins, linearly interpolating the data.
-Rescale data to some base mass (Jpsi) and build ratio with Jpsi distributions.
-Do this for all Psi(nS) and Ups(nS) data
-If reasonably flat, fit constant for each state, plot constants as function of M(quarkonium)/M(Jpsi)
-> Do we see a clear trend with the mass ratio?


	 */

	//STEP4: Scale the extrapolated SDCs and polarizations (pT and normalization scaling)

	//Current output 'scaled': all charmonium states scaled from Jpsi, all bottomonium states scaled from Ups1S
	//Exception: chic1,2 singlet scaled from chic1,2 singlet

	TGraph *SDC_Graph_scaled[NRQCDvars::nStates][nEpRapIntervals][nColorChannelsGiven];
	TGraph *Lamth_Graph_scaled[NRQCDvars::nStates][nEpRapIntervals][nColorChannelsGiven];
	TGraph *Lamph_Graph_scaled[NRQCDvars::nStates][nEpRapIntervals][nColorChannelsGiven];
	TGraph *Lamtp_Graph_scaled[NRQCDvars::nStates][nEpRapIntervals][nColorChannelsGiven];

	const int const_nBinsFinalModels=nBinsFinalModels;
	double buff_array1[const_nBinsFinalModels];//={NULL};
	double buff_array2[const_nBinsFinalModels];//={NULL};

	for(int i=0;i<NRQCDvars::nStates;i++){
		for(int j=0;j<nEpRapIntervals;j++){
			for(int k=0;k<nColorChannelsGiven;k++){

				SDC_Graph_scaled[i][j][k]= new TGraph(nBinsFinalModels,buff_array1,buff_array2);
				Lamth_Graph_scaled[i][j][k]= new TGraph(nBinsFinalModels,buff_array1,buff_array2);
				Lamph_Graph_scaled[i][j][k]= new TGraph(nBinsFinalModels,buff_array1,buff_array2);
				Lamtp_Graph_scaled[i][j][k]= new TGraph(nBinsFinalModels,buff_array1,buff_array2);
				sprintf(graphName,"SDC_Graph_scaled_state%d_rap%d_CC%d",i,j,k);
				SDC_Graph_scaled[i][j][k]->SetName(graphName);
				sprintf(graphName,"Lamth_Graph_scaled_state%d_rap%d_CC%d",i,j,k);
				Lamth_Graph_scaled[i][j][k]->SetName(graphName);
				sprintf(graphName,"Lamph_Graph_scaled_state%d_rap%d_CC%d",i,j,k);
				Lamph_Graph_scaled[i][j][k]->SetName(graphName);
				sprintf(graphName,"Lamtp_Graph_scaled_state%d_rap%d_CC%d",i,j,k);
				Lamtp_Graph_scaled[i][j][k]->SetName(graphName);

			}
		}
	}

	//Scale charmonium models to bottomonium with mass ratio
	//bool ScaleCharmModelToBottom=true;
	//const int nScaleToBottomStates=1;
	//int ScaleToBottomStates[nScaleToBottomStates]={10};

   //Luping Information:
   //For the charmonium case, the mass of J/\psi  is 2x1.5 = 3.0 GeV
   //For psi(2S), it is also 1.5*2=3.0 GeV;
   //For Upsilon(1S), it is 4.75*2=9.5 GeV;
   //For Upsilon(2S), it is 5.11*2=10.22 GeV;
   //For Upsilon(3S), it is 5.2*2=10.4 GeV;

	int ScaleFromStateCC=0;
	int ScaleFromStateBB=4;
	double ScaleFromMassCC=3;
	double ScaleFromMassBB=9.5;
	double MassRatioScale[NRQCDvars::nStates];
	double NormalizationScale[NRQCDvars::nStates];

	// x-axis scale:::

	int ScaleFromStateCC_asStatesGiven=-1;
	int ScaleFromStateBB_asStatesGiven=-1;
	for(int iScaleState=0; iScaleState<nStatesGiven; iScaleState++){
		if(StatesGiven[iScaleState]==ScaleFromStateCC) ScaleFromStateCC_asStatesGiven=iScaleState;
		if(StatesGiven[iScaleState]==ScaleFromStateBB) ScaleFromStateBB_asStatesGiven=iScaleState;
	}

	for(int iScaleState=0; iScaleState<NRQCDvars::nStates; iScaleState++){
		if(iScaleState<4) MassRatioScale[iScaleState]=NRQCDvars::mass[iScaleState]/ScaleFromMassCC;
		else MassRatioScale[iScaleState]=NRQCDvars::mass[iScaleState]/ScaleFromMassBB;
		cout<<"MassRatioScale["<<StateName[iScaleState]<<"] = "<<MassRatioScale[iScaleState]<<endl;
	}



	// Normalization scale:::


	double expnorm = 4.03556;
	double exponent = -0.940807;

	double massmin=1;
	double massmax=13;

	sprintf(name, "expo_fit");
	TF1* expo_fit  = new TF1(name, exponential, massmin, massmax, 2);
	expo_fit->FixParameter(0,expnorm);
	expo_fit->FixParameter(1,exponent);


	for(int iScaleState=0; iScaleState<NRQCDvars::nStates; iScaleState++){
		if(iScaleState<4) NormalizationScale[iScaleState]=expo_fit->Eval(NRQCDvars::mass[iScaleState])/expo_fit->Eval(ScaleFromMassCC);
		else NormalizationScale[iScaleState]=expo_fit->Eval(NRQCDvars::mass[iScaleState])/expo_fit->Eval(ScaleFromMassBB);
		cout<<"NormalizationScale["<<StateName[iScaleState]<<"] = "<<NormalizationScale[iScaleState]<<endl;
	}

	/////////////////////////

	//inputs reminder:
	// SDC_Graph[nStatesGiven][nRapIntervals][nColorChannelsGiven];
	// Lamth_Graph[nStatesGiven][nRapIntervals][nColorChannelsGiven];
	// Lamph_Graph[nStatesGiven][nRapIntervals][nColorChannelsGiven];
	// Lamtp_Graph[nStatesGiven][nRapIntervals][nColorChannelsGiven];



	//TODO: open file
	char outfilename[200];
	sprintf(outfilename,"%s/TGraphs_scaled.root",originalmodeldirname);cout<<outfilename<<endl;
	TFile *outfile = new TFile(outfilename,"RECREATE");



	for(int i=0;i<nStatesGiven;i++){
		for(int j=0;j<nRapIntervals;j++){
			for(int k=0;k<nColorChannelsGiven;k++){

				SDC_Graph[i][j][k]->Write();
				Lamth_Graph[i][j][k]->Write();
				Lamph_Graph[i][j][k]->Write();
				Lamtp_Graph[i][j][k]->Write();

				SDC_Graph_original[i][j][k]->Write();
				Lamth_Graph_original[i][j][k]->Write();
				Lamph_Graph_original[i][j][k]->Write();
				Lamtp_Graph_original[i][j][k]->Write();


			}
		}
	}

	for(int i=0;i<nStatesGiven;i++){
		for(int j=0;j<nEpRapIntervals;j++){
			for(int k=0;k<nColorChannelsGiven;k++){

				EpSDC_Graph[i][j][k]->Write();
				EpLamth_Graph[i][j][k]->Write();
				EpLamph_Graph[i][j][k]->Write();
				EpLamtp_Graph[i][j][k]->Write();

			}
		}
	}




	bool Scale_pT=true;
	double buff_pT, buff_val;
	double buff_pT_scaled, buff_val_scaled;
	double buff_p, buff_p_scaled;

	for(int i=0;i<NRQCDvars::nStates;i++){
		for(int j=0;j<nEpRapIntervals;j++){
			int nColorChannels_state;
			bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
			if(isSstate) nColorChannels_state=NRQCDvars::nColorChannels_S;
			else nColorChannels_state=NRQCDvars::nColorChannels_P;
			for (int k=0; k<nColorChannels_state; k++){


				for(int l=0;l<nBinsFinalModels;l++){

					int ScaleFromState_asStatesGiven;
					int ScaleFromState;
					if(i<4) {
						ScaleFromState_asStatesGiven=ScaleFromStateCC_asStatesGiven;
						ScaleFromState=ScaleFromStateCC;
					}
					else{
						ScaleFromState_asStatesGiven=ScaleFromStateBB_asStatesGiven;
						ScaleFromState=ScaleFromStateBB;
					}

					if( (i==1||i==2) && k==0 ){
						ScaleFromState_asStatesGiven=i;
					}

					int kPrime=k;

					EpSDC_Graph[ScaleFromState_asStatesGiven][j][kPrime]->GetPoint(l,buff_pT, buff_val);

					//actually scale pT:
					if(Scale_pT){
						buff_pT_scaled=buff_pT*MassRatioScale[i];
						buff_val_scaled=buff_val*NormalizationScale[i];
					}
					else{

						buff_p=TMath::Sqrt(NRQCDvars::mass[ScaleFromState]*NRQCDvars::mass[ScaleFromState]*TMath::SinH(AverageRap[j])*TMath::SinH(AverageRap[j])+buff_pT*buff_pT*TMath::CosH(AverageRap[j])*TMath::CosH(AverageRap[j]));
						buff_p_scaled=buff_p*MassRatioScale[i];
						buff_pT_scaled=TMath::Sqrt((buff_p_scaled*buff_p_scaled-NRQCDvars::mass[i]*NRQCDvars::mass[i]*TMath::SinH(AverageRap[j])*TMath::SinH(AverageRap[j]))/(TMath::CosH(AverageRap[j])*TMath::CosH(AverageRap[j])));
						buff_val_scaled=buff_val*NormalizationScale[i];
					}

					//buff_pT_scaled=buff_pT;

					//Dummy model for low pT Upsilon polarization
					//if(l==0&&i==10) buff_pT_scaled=10;

					SDC_Graph_scaled[i][j][k]->SetPoint(l,buff_pT_scaled, buff_val_scaled);

					EpLamth_Graph[ScaleFromState_asStatesGiven][j][kPrime]->GetPoint(l,buff_pT, buff_val);
					Lamth_Graph_scaled[i][j][k]->SetPoint(l,buff_pT_scaled, buff_val);
					EpLamph_Graph[ScaleFromState_asStatesGiven][j][kPrime]->GetPoint(l,buff_pT, buff_val);
					Lamph_Graph_scaled[i][j][k]->SetPoint(l,buff_pT_scaled, buff_val);
					EpLamtp_Graph[ScaleFromState_asStatesGiven][j][kPrime]->GetPoint(l,buff_pT, buff_val);
					Lamtp_Graph_scaled[i][j][k]->SetPoint(l,buff_pT_scaled, buff_val);


				}

				//TODO: write graph to file
				outfile->cd();
				//plotGraph->Draw("P");
				SDC_Graph_scaled[i][j][k]->Write();
				Lamth_Graph_scaled[i][j][k]->Write();
				Lamph_Graph_scaled[i][j][k]->Write();
				Lamtp_Graph_scaled[i][j][k]->Write();

			}
		}
	}

	outfile->Write();
	outfile->Close();
	delete outfile;
	outfile = NULL;


	char outfilenameTo[200];
	if(useTF1inputs) {sprintf(outfilenameTo,"%s/TGraphs_scaled_fromTF1.root",originalmodeldirname);cout<<outfilenameTo<<endl;}
	if(!useTF1inputs) {sprintf(outfilenameTo,"%s/TGraphs_scaled_noTF1.root",originalmodeldirname);cout<<outfilenameTo<<endl;}


	gSystem->CopyFile(outfilename,outfilenameTo,kTRUE);


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

double exponential(double* x, double* par) {

	double funcval;

	funcval=par[0]*exp(par[1]*x[0]);

	return funcval;

}

Double_t paramMassRapParabola(Double_t *x, Double_t *par){
	Double_t result=par[0]+x[0]*par[1]+x[0]*x[0]*par[2];
	return result;
}

/*
 * ScaleBKmodel.cc
 *
 *  Created on: Oct 1, 2013
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

  	for( int i=0;i < argc; ++i ) {
		if(std::string(argv[i]).find("OriginalModelID") != std::string::npos) {char* OriginalModelIDchar = argv[i]; char* OriginalModelIDchar2 = strtok (OriginalModelIDchar, "="); OriginalModelID = OriginalModelIDchar2; cout<<"OriginalModelID = "<<OriginalModelID<<endl;}
		if(std::string(argv[i]).find("storagedir") != std::string::npos) {char* storagedirchar = argv[i]; char* storagedirchar2 = strtok (storagedirchar, "="); storagedir = storagedirchar2; cout<<"storagedir = "<<storagedir<<endl;}
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
// Convert BK input model
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Matrix elements:

	const int nStatesGiven=2;//Given models in any form
	int StatesGiven[nStatesGiven]={0, 3};
	const int nColorChannelsGiven=4;

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
	double RapIntervalBordersMin[nRapIntervals]={0, 0.6, 2.5};
	double RapIntervalBordersMax[nRapIntervals]={0.6, 1.2, 4.};
	double AverageRap[nRapIntervals]={0.3, 0.9, 3.25};
	const int npTBinsPerRap[nRapIntervals][nHelicityChannels]={
			{12, 12, 12},
			{7, 7, 7},
			{8, 10, 11}
	};
	const int maxpTBinsPerRap=12;

	double deltaRapOriginal[nRapIntervals];
	for(int iRap=0;iRap<nRapIntervals;iRap++){
		deltaRapOriginal[iRap]=RapIntervalBordersMax[iRap]-RapIntervalBordersMin[iRap];
		if(isAbsRap[iRap]) deltaRapOriginal[iRap]*=2;
	}

	bool interpretOriginalModelAsIntegratedInRap=true;

	bool isDummyRap[nRapIntervals]={false, false, false};



// CMS rapidities




	double pTmean[nStatesGiven][nRapIntervals][nHelicityChannels][maxpTBinsPerRap];//={

	double pTMin[nRapIntervals]={10, 10, 3};
	double pTMax[nRapIntervals]={70, 30, 10};

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

				//cout<<"ratio sum_cc/total (1S) = "<<(model[0][j][0][l][m]+model[0][j][1][l][m]+model[0][j][2][l][m]+model[0][j][3][l][m])/total_model[0][j][l][m]<<endl;
				//cout<<"ratio sum_cc/total (2S) = "<<(model[1][j][0][l][m]+model[1][j][1][l][m]+model[1][j][2][l][m]+model[1][j][3][l][m])/total_model[1][j][l][m]<<endl;

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

					//cout<<"model_Graph[state"<<i<<"][rap"<<j<<"][CC"<<k<<"][HC"<<l<<"]"<<endl;
					//model_Graph[i][j][k][l]->Print();

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
				//cout<<"----------------------------"<<endl;
				//cout<<"meanRapRatio_1ov0 "<<meanRapRatio_1ov0<<endl;
				//cout<<"----------------------------"<<endl;

				if(meanRapRatio_1ov0!=0){
				scaleFactor[k]+=meanRapRatio_1ov0;
				countAverage++;
				}
			}
		}
			scaleFactor[k]/=double(countAverage);
			cout<<"(sigma) scaleFactor[k] "<<scaleFactor[k]<<endl;
			//cout<<"(SDC) scaleFactor[k]*fittedLDMEs[0][k]/fittedLDMEs[1][k] "<<scaleFactor[k]*fittedLDMEs[0][k]/fittedLDMEs[1][k]<<endl;
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

						//cout<<"model_Graph[state"<<i<<"][rap"<<j<<"][CC"<<k<<"][HC"<<l<<"]"<<endl;
						//model_Graph[i][j][k][l]->Print();

					}
				}
			}
		}



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

				int nBinsOriginal=npTBinsPerRap[j][0];
				if(j==2) nBinsOriginal=nBinsFinalModels;

				double pTmean_graph_original[nBinsOriginal];
				double SDC_graph_original[nBinsOriginal];
				double Lamth_graph_original[nBinsOriginal];
				double Lamph_graph_original[nBinsOriginal];
				double Lamtp_graph_original[nBinsOriginal];
				if(j<2){
					for(int m=0;m<nBinsOriginal;m++){
						double buff_sigma00;
						double buff_sigma11;
						double buff_sigma1m1;
						model_Graph[i][j][k][0]->GetPoint(m,pTmean_graph_original[m],buff_sigma00);
						model_Graph[i][j][k][1]->GetPoint(m,pTmean_graph_original[m],buff_sigma11);
						model_Graph[i][j][k][2]->GetPoint(m,pTmean_graph_original[m],buff_sigma1m1);
						//if(pTmean_graph_original[m]<45) model_Graph[i][j][k][2]->GetPoint(m,pTmean_graph_original[m],buff_sigma1m1);
						//else buff_sigma1m1=0.;

						SDC_graph_original[m]=buff_sigma00/fittedLDMEs[i][k]+2*buff_sigma11/fittedLDMEs[i][k];
						double polDenominator=(buff_sigma11/fittedLDMEs[i][k]+buff_sigma00/fittedLDMEs[i][k]);
						if(UseChaoPolDenominatorDef && buff_sigma11 * buff_sigma00 < 0) polDenominator=fabs(polDenominator);
						Lamth_graph_original[m]=(buff_sigma11/fittedLDMEs[i][k]-buff_sigma00/fittedLDMEs[i][k])/polDenominator;
						Lamph_graph_original[m]=(buff_sigma1m1/fittedLDMEs[i][k])/polDenominator;
						Lamtp_graph_original[m]=0;

						double dummyVal=999;
						if(isDummyRap[j]){
							SDC_graph_original[m]=dummyVal;
							Lamth_graph_original[m]=dummyVal;
							Lamph_graph_original[m]=dummyVal;
							Lamtp_graph_original[m]=dummyVal;
						}
					}
				}
				else{
					for(int m=0;m<nBinsOriginal;m++){

						double deltaPt=(pTMax[j]-pTMin[j])/double(nBinsOriginal);
						pTmean_graph_original[m]=pTMin[j]+m*deltaPt;

						double buff_sigma00;
						double buff_sigma11;
						double buff_sigma1m1;

						buff_sigma00 = model_Graph[i][j][k][0]->Eval(pTmean_graph_original[m]);
						buff_sigma11 = model_Graph[i][j][k][1]->Eval(pTmean_graph_original[m]);
						buff_sigma1m1 = model_Graph[i][j][k][2]->Eval(pTmean_graph_original[m]);

						//if(pTmean_graph_original[m]<45) model_Graph[i][j][k][2]->GetPoint(m,pTmean_graph_original[m],buff_sigma1m1);
						//else buff_sigma1m1=0.;

						SDC_graph_original[m]=buff_sigma00/fittedLDMEs[i][k]+2*buff_sigma11/fittedLDMEs[i][k];
						double polDenominator=(buff_sigma11/fittedLDMEs[i][k]+buff_sigma00/fittedLDMEs[i][k]);
						if(UseChaoPolDenominatorDef && buff_sigma11 * buff_sigma00 < 0) polDenominator=fabs(polDenominator);
						Lamth_graph_original[m]=(buff_sigma11/fittedLDMEs[i][k]-buff_sigma00/fittedLDMEs[i][k])/polDenominator;
						Lamph_graph_original[m]=(buff_sigma1m1/fittedLDMEs[i][k])/polDenominator;
						Lamtp_graph_original[m]=0;

						double dummyVal=999;
						if(isDummyRap[j]){
							SDC_graph_original[m]=dummyVal;
							Lamth_graph_original[m]=dummyVal;
							Lamph_graph_original[m]=dummyVal;
							Lamtp_graph_original[m]=dummyVal;
						}
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
					double deltaPt=(pTMax[j]-pTMin[j])/double(nBinsFinalModels);
					pTmean_graph[m]=pTMin[j]+m*deltaPt;
					//SDC_graph[m]=SDC_Graph_original[i][j][k]->Eval(pTmean_graph[m],0,"S");

					double rapCorrFactor=deltaRapOriginal[j];
					if(!interpretOriginalModelAsIntegratedInRap) rapCorrFactor=1.;

					SDC_graph[m]=SDC_Graph_original[i][j][k]->Eval(pTmean_graph[m])/rapCorrFactor;
					Lamth_graph[m]=Lamth_Graph_original[i][j][k]->Eval(pTmean_graph[m]);
					Lamph_graph[m]=0;//Lamph_Graph_original[i][j][k]->Eval(pTmean_graph[m]);
					Lamtp_graph[m]=0;
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




	SDC_Graph_original[0][0][3]->Print();

	double evalPt[6]={10.5, 11.48, 12.72, 14.23, 16.3, 21.92};
	for(int u=0;u<6;u++){
		double crossSection=0;
		for(int k=0;k<nColorChannelsGiven;k++){
			crossSection+=SDC_Graph[1][1][k]->Eval(evalPt[u])*fittedLDMEs[1][k];
		}
		//cout<<"crossSection at "<<evalPt[u]<<" = "<<crossSection<<endl;
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
			//cout<<"RapRatio_1ov0 "<<RapRatio_1ov0<<" RapRatio_2ov1 "<<RapRatio_2ov1<<endl;

			double RapDiffLamth_1m0;
			double RapDiffLamth_2m1;
			RapDiffLamth_1m0=Lamth_Graph[i][1][k]->Eval(10.)-Lamth_Graph[i][0][k]->Eval(10.);
			RapDiffLamth_2m1=Lamth_Graph[i][2][k]->Eval(10.)-Lamth_Graph[i][1][k]->Eval(10.);
			//cout<<"RapDiffLamth_1m0 "<<RapDiffLamth_1m0<<" RapDiffLamth_2m1 "<<RapDiffLamth_2m1<<endl;
			double RapDiffLamph_1m0;
			double RapDiffLamph_2m1;
			RapDiffLamph_1m0=Lamph_Graph[i][1][k]->Eval(10.)-Lamph_Graph[i][0][k]->Eval(10.);
			RapDiffLamph_2m1=Lamph_Graph[i][2][k]->Eval(10.)-Lamph_Graph[i][1][k]->Eval(10.);
			//cout<<"RapDiffLamph_1m0 "<<RapDiffLamph_1m0<<" RapDiffLamph_2m1 "<<RapDiffLamph_2m1<<endl;

		}
	}


	double SDC_StateRatio_3ov0;
	for(int j=0;j<nRapIntervals;j++){
		for(int k=0;k<nColorChannelsGiven;k++){
			SDC_StateRatio_3ov0=SDC_Graph[1][j][k]->Eval(10.)/SDC_Graph[0][j][k]->Eval(10.);
			cout<<"SDC_StateRatio_3ov0 at 10 GeV"<<SDC_StateRatio_3ov0<<endl;
			SDC_StateRatio_3ov0=SDC_Graph[1][j][k]->Eval(20.)/SDC_Graph[0][j][k]->Eval(20.);
			cout<<"SDC_StateRatio_3ov0 at 20 GeV"<<SDC_StateRatio_3ov0<<endl;
			SDC_StateRatio_3ov0=SDC_Graph[1][j][k]->Eval(30.)/SDC_Graph[0][j][k]->Eval(30.);
			cout<<"SDC_StateRatio_3ov0 at 30 GeV"<<SDC_StateRatio_3ov0<<endl;
			SDC_StateRatio_3ov0=SDC_Graph[1][j][k]->Eval(40.)/SDC_Graph[0][j][k]->Eval(40.);
			cout<<"SDC_StateRatio_3ov0 at 40 GeV"<<SDC_StateRatio_3ov0<<endl;
		}
	}




/*
	SDC extrapolation:
	1) normalize the rap2 curve to match the rap0 curve at 10 GeV. The same can be done for rap1, and (similar) for rap2
	2) We want to fill the hole between y=1.2 and y=2.5, to be able to add more data. For this, I will interpolate linearly the normalization factor between the rap1 and rap2 curves. With this, we will be able to add a rapidity bin of 1.2<|y|<2.5

	Lamth extrapolation:
	S0: everything is 0 :)
	S1 CO rap0: Shift rap2 to match rap0
	S1 CO rap1: Shift rap2 to match rap1
	S1 CO rap2: Shift rap1 to match rap2
	S1 CS rap0: Shift rap2 to match rap0
	S1 CS rap1: Shift rap2 to match rap1
	S1 CS rap2: Shift rap1 to match rap2
*/

/*

	define finalRap-stuff
	start from these:
	SDC_Graph[i][j][k]
	Lamth_Graph[i][j][k]
	Lamph_Graph[i][j][k]
	Lamtp_Graph[i][j][k]
	Take Care Of integrated Rap-factor

*/


	const int nEpRapIntervals=9;
	//bool isAbsEpRap[nEpRapIntervals]={true, true, true, true, true, true, true, true, true};
	double EpRapIntervalBordersMin[nEpRapIntervals]={0, 0.6, 1.2, 1.6, 2., 2.5, 3., 3.5, 4.};
	double EpRapIntervalBordersMax[nEpRapIntervals]={0.6, 1.2, 1.6, 2., 2.5, 3., 3.5, 4., 4.5};
	double AverageEpRap[nEpRapIntervals]={0.3, 0.9, 1.4, 1.8, 2.25, 2.75, 3.25, 3.75, 4.25};

	double pTMinEp[nEpRapIntervals]={2.875, 2.875, 2.875, 2.875, 2.875, 2.875, 2.875, 2.875, 2.875};
	double pTMaxEp[nEpRapIntervals]={70, 70, 70, 70, 70, 70, 70, 70, 70};

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

					double deltaPt=(pTMaxEp[j]-pTMinEp[j])/double(nBinsFinalModels);
					pTmean_graph_Ep[m]=pTMinEp[j]+m*deltaPt;

					if(j==0){

						if(pTmean_graph_Ep[m]>=normAtPt) SDC_graph_Ep[m]=SDC_Graph[i][j][k]->Eval(pTmean_graph_Ep[m]);
						else SDC_graph_Ep[m]=SDC_Graph[i][2][k]->Eval(pTmean_graph_Ep[m])/NormRapXToRap0Factor[2][k];

					}
					else{

						SDC_graph_Ep[m]=EpSDC_Graph[i][0][k]->Eval(pTmean_graph_Ep[m])*fParabola[k]->Eval(AverageEpRap[j]);

					}



					if(k==1){
						Lamth_graph_Ep[m]=0;//Lamth_Graph[i][j][k]->Eval(pTmean_graph_Ep[m]);
					}
					else{

						double Lamth_j0;
						double Lamth_j1;
						double Lamth_j2;

						if(pTmean_graph_Ep[m]>=diffAtPt){
							Lamth_j0=Lamth_Graph[i][0][k]->Eval(pTmean_graph_Ep[m]);
							Lamth_j1=Lamth_Graph[i][1][k]->Eval(pTmean_graph_Ep[m]);
							if(pTmean_graph_Ep[m]>=30) Lamth_j1=Lamth_Graph[i][0][k]->Eval(pTmean_graph_Ep[m])+Lamth_Graph[i][1][k]->Eval(30)-Lamth_Graph[i][0][k]->Eval(30);
							Lamth_j2=Lamth_j1-LamthDiffToRap2[1][k];
						}
						else{
							Lamth_j0=Lamth_Graph[i][diffToRap][k]->Eval(pTmean_graph_Ep[m])+LamthDiffToRap2[0][k];
							Lamth_j1=Lamth_Graph[i][diffToRap][k]->Eval(pTmean_graph_Ep[m])+LamthDiffToRap2[1][k];
							Lamth_j2=Lamth_Graph[i][diffToRap][k]->Eval(pTmean_graph_Ep[m])+LamthDiffToRap2[2][k];
						}

						if(k==3 && pTmean_graph_Ep[m] < 8.71){
							Lamth_j0=Lamth_Graph[i][diffToRap][k]->Eval(pTmean_graph_Ep[m]);
							Lamth_j1=Lamth_Graph[i][diffToRap][k]->Eval(pTmean_graph_Ep[m]);
							Lamth_j2=Lamth_Graph[i][diffToRap][k]->Eval(pTmean_graph_Ep[m]);
						}

						Lamth_graph_Ep[m]=Lamth_j1+(Lamth_j2-Lamth_j1)*(AverageEpRap[j]-AverageRap[1])/(AverageRap[2]-AverageRap[1]);


						if(j==0) Lamth_graph_Ep[m]=Lamth_j0;
						if(j==1) Lamth_graph_Ep[m]=Lamth_j1;


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

	int ScaleFromState=0;
	double ScaleFromMass=3;
	bool doScaleFromMass=true;
	double MassRatioScale[NRQCDvars::nStates];
	double NormalizationScale[NRQCDvars::nStates];

	// x-axis scale:::

	int ScaleFromState_asStatesGiven=-1;
	for(int iScaleState=0; iScaleState<nStatesGiven; iScaleState++){
		if(StatesGiven[iScaleState]==ScaleFromState) ScaleFromState_asStatesGiven=iScaleState;
	}

	for(int iScaleState=0; iScaleState<NRQCDvars::nStates; iScaleState++){
		if(!doScaleFromMass) MassRatioScale[iScaleState]=NRQCDvars::mass[iScaleState]/NRQCDvars::mass[ScaleFromState];
		else MassRatioScale[iScaleState]=NRQCDvars::mass[iScaleState]/ScaleFromMass;
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
		if(!doScaleFromMass) NormalizationScale[iScaleState]=expo_fit->Eval(NRQCDvars::mass[iScaleState])/expo_fit->Eval(NRQCDvars::mass[ScaleFromState]);
		else NormalizationScale[iScaleState]=expo_fit->Eval(NRQCDvars::mass[iScaleState])/expo_fit->Eval(ScaleFromMass);
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

				int kPrime=k;
				//if(!isSstate&&k==1) kPrime=2;

				for(int l=0;l<nBinsFinalModels;l++){

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

	//TODO: close, write file
	outfile->Write();
	outfile->Close();
	delete outfile;
	outfile = NULL;













/*

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

				cout<<"globalWeight "<<globalWeight<<endl;

				cout<<"Generated "<<n_nTuple<<" events"<<endl;
				cout<<"nTupleModel[iScaleToMother][iColorChannel]->GetEntries() "<<nTupleModel[iScaleToMother][iColorChannel]->GetEntries()<<endl;



			}
		}

*/


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

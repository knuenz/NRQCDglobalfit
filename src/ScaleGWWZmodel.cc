/*
 * ScaleGWWZmodel.cc
 *
 *  Created on: Oct 27, 2013
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
double exponential(double* x, double* par);

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


	const int nHelicityChannels=2;	//00, 11, 1m1

	const int nRapIntervals=3;
	bool isAbsRap[nRapIntervals]={true, true, true};
	double RapIntervalBordersMin[nRapIntervals]={0, 0.6, 1.2};
	double RapIntervalBordersMax[nRapIntervals]={0.6, 1.2, 1.5};
	double AverageRap[nRapIntervals]={0.3, 0.9, 1.35};
	const int npTBinsPerRap[nRapIntervals][nHelicityChannels]={
			{24, 24},
			{24, 24},
			{24, 24}
	};
	const int maxpTBinsPerRap=24;


	bool isDummyRap[nRapIntervals]={false, false, false};



// CMS rapidities




	double pTmean[nStatesGiven][nRapIntervals][nColorChannelsGiven][nHelicityChannels][maxpTBinsPerRap];//={
	//		{
	//				{10.00, 12.50, 15.00, 17.50, 20.00, 25.00, 30.00, 35.00, 40.00},
	//				{10.00, 12.50, 15.00, 17.50, 20.00, 25.00, 30.00, 00.00, 00.00}
	//		},
	//		{
	//				{10.00, 12.50, 15.00, 17.50, 20.00, 25.00, 30.00, 35.00, 40.00},
	//				{10.00, 12.50, 15.00, 17.50, 20.00, 25.00, 30.00, 00.00, 00.00}
	//		}
	//};

	double pTMin[nRapIntervals]={3, 3, 3};
	double pTMax[nRapIntervals]={75, 75, 75};

	double model[nStatesGiven][nRapIntervals][nColorChannelsGiven][nHelicityChannels][maxpTBinsPerRap];
	double err_model[nStatesGiven][nRapIntervals][nColorChannelsGiven][nHelicityChannels][maxpTBinsPerRap];

	cout<<"read in file"<<endl;

	char in_rapChar[200];
	char in_stateChar[200];
	char in_helcolChar[200];

	for(int i=0;i<nStatesGiven;i++){
		cout<<"State"<<StatesGiven[i]<<endl;

	    if(i==0) sprintf(in_stateChar,"Jpsi/direct");
	    if(i==1) sprintf(in_stateChar,"Psi2s/direct");

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

					sprintf(inname,"%s/Jpsi_at_CMS/%s/%s/%s/data",originalmodeldirname,in_rapChar,in_stateChar,in_helcolChar);

					FILE *fIn;
					fIn = fopen(inname, "read");
					cout<<inname<<endl;
					Char_t line[1000];

					for(int m=0;m<npTBinsPerRap[j][l];m++){

						double buffer;
						fgets(line, sizeof(line), fIn); //comment
						cout<<line<<endl;
						sscanf(line, "%lf %lf %lf", &pTmean[i][j][k][l][m], &model[i][j][k][l][m], &err_model[i][j][k][l][m]);

						if(k==1){
							model[i][j][k][l][m]/=3.;
							err_model[i][j][k][l][m]/=3.;
						}

						cout<<pTmean[i][j][k][l][m]<<" "<<model[i][j][k][l][m]<<" "<<err_model[i][j][k][l][m]<<" "<<endl;

					}
					fgets(line, sizeof(line), fIn); //comment
				}
			}
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
						pTmean_graph[m]=pTmean[i][j][k][l][m];
						model_graph[m]=model[i][j][k][l][m];
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

	char graphName[200];

	int nBinsFinalModels=500;
	for(int i=0;i<nStatesGiven;i++){
		for(int j=0;j<nRapIntervals;j++){
			for(int k=0;k<nColorChannelsGiven;k++){




				cout<<"originals: rap "<<j<<" CC "<<k<<" state "<<i<<endl;
				double pTmean_graph_original[npTBinsPerRap[j][0]];
				double SDC_graph_original[npTBinsPerRap[j][0]];
				double Lamth_graph_original[npTBinsPerRap[j][0]];
				double Lamph_graph_original[npTBinsPerRap[j][0]];
				double Lamtp_graph_original[npTBinsPerRap[j][0]];
				for(int m=0;m<npTBinsPerRap[j][0];m++){
					double buff_sigma00;
					double buff_sigma11;
					model_Graph[i][j][k][0]->GetPoint(m,pTmean_graph_original[m],buff_sigma00);
					model_Graph[i][j][k][1]->GetPoint(m,pTmean_graph_original[m],buff_sigma11);

					//if(k==1) cout<<"pol asymmetry "<<buff_sigma00-buff_sigma11<<endl;

					double dummyVal=999;

					SDC_graph_original[m]=buff_sigma00+2*buff_sigma11;
					double polDenominator=(buff_sigma11+buff_sigma00);
					Lamth_graph_original[m]=(buff_sigma11-buff_sigma00)/polDenominator;
					Lamph_graph_original[m]=dummyVal;
					Lamtp_graph_original[m]=dummyVal;

					if(isDummyRap[j]){
						SDC_graph_original[m]=dummyVal;
						Lamth_graph_original[m]=dummyVal;
						Lamph_graph_original[m]=dummyVal;
						Lamtp_graph_original[m]=dummyVal;
					}
				}

				SDC_Graph_original[i][j][k]= new TGraph(npTBinsPerRap[j][0],pTmean_graph_original,SDC_graph_original);
				Lamth_Graph_original[i][j][k]= new TGraph(npTBinsPerRap[j][0],pTmean_graph_original,Lamth_graph_original);
				Lamph_Graph_original[i][j][k]= new TGraph(npTBinsPerRap[j][0],pTmean_graph_original,Lamph_graph_original);
				Lamtp_Graph_original[i][j][k]= new TGraph(npTBinsPerRap[j][0],pTmean_graph_original,Lamtp_graph_original);

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
					SDC_graph[m]=SDC_Graph_original[i][j][k]->Eval(pTmean_graph[m]);
					Lamth_graph[m]=Lamth_Graph_original[i][j][k]->Eval(pTmean_graph[m]);
					Lamph_graph[m]=Lamph_Graph_original[i][j][k]->Eval(pTmean_graph[m]);
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


	TGraph *SDC_Graph_scaled[NRQCDvars::nStates][nRapIntervals][nColorChannelsGiven];
	TGraph *Lamth_Graph_scaled[NRQCDvars::nStates][nRapIntervals][nColorChannelsGiven];
	TGraph *Lamph_Graph_scaled[NRQCDvars::nStates][nRapIntervals][nColorChannelsGiven];
	TGraph *Lamtp_Graph_scaled[NRQCDvars::nStates][nRapIntervals][nColorChannelsGiven];

	const int const_nBinsFinalModels=nBinsFinalModels;
	double buff_array1[const_nBinsFinalModels];//={NULL};
	double buff_array2[const_nBinsFinalModels];//={NULL};

	for(int i=0;i<NRQCDvars::nStates;i++){
		for(int j=0;j<nRapIntervals;j++){
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

	char name[200];
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




	bool Scale_pT=true;
	double buff_pT, buff_val;
	double buff_pT_scaled, buff_val_scaled;
	double buff_p, buff_p_scaled;

	for(int i=0;i<NRQCDvars::nStates;i++){
		for(int j=0;j<nRapIntervals;j++){
			int nColorChannels_state;
			bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
			if(isSstate) nColorChannels_state=NRQCDvars::nColorChannels_S;
			else nColorChannels_state=NRQCDvars::nColorChannels_P;
			for (int k=0; k<nColorChannels_state; k++){

				int kPrime=k;
				if(!isSstate&&k==1) kPrime=2;

				for(int l=0;l<nBinsFinalModels;l++){

					SDC_Graph[ScaleFromState_asStatesGiven][j][kPrime]->GetPoint(l,buff_pT, buff_val);

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
					if(l==0&&i==10) buff_pT_scaled=10;

					SDC_Graph_scaled[i][j][k]->SetPoint(l,buff_pT_scaled, buff_val_scaled);

					Lamth_Graph[ScaleFromState_asStatesGiven][j][kPrime]->GetPoint(l,buff_pT, buff_val);
					Lamth_Graph_scaled[i][j][k]->SetPoint(l,buff_pT_scaled, buff_val);
					Lamph_Graph[ScaleFromState_asStatesGiven][j][kPrime]->GetPoint(l,buff_pT, buff_val);
					Lamph_Graph_scaled[i][j][k]->SetPoint(l,buff_pT_scaled, buff_val);
					Lamtp_Graph[ScaleFromState_asStatesGiven][j][kPrime]->GetPoint(l,buff_pT, buff_val);
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


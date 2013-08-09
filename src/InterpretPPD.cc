/*
 * InterpretPPD.cc
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

//#include "TMatrixD.h"

using namespace NRQCDvars;

void PlotPosterior(char xAxisTitle[200], char filename[200], TH1* histo, double MPV, double MPVerrLow, double MPVerrHigh);
void FindMPV(TH1* PosteriorDist , double& MPV , double& MPVerrorLow, double& MPVerrorHigh, int MPValgo, int nSigma);

dmatrix transformDmatrixToSPformat(dmatrix dmatrixFull);

int main(int argc, char** argv) {

  	Char_t *JobID = "Default";
	Char_t *storagedir = "Default"; //Storage Directory

  	int 	MPValgo=-1;
  	double 	nSigma=-1;


  	for( int i=0;i < argc; ++i ) {
		if(std::string(argv[i]).find("storagedir") != std::string::npos) {char* storagedirchar = argv[i]; char* storagedirchar2 = strtok (storagedirchar, "="); storagedir = storagedirchar2; cout<<"storagedir = "<<storagedir<<endl;}
	    if(std::string(argv[i]).find("nSigma") != std::string::npos) { char* nSigmachar = argv[i]; char* nSigmachar2 = strtok (nSigmachar, "n"); nSigma = atof(nSigmachar2); cout<<"nSigma = "<<nSigma<<endl; }
	    if(std::string(argv[i]).find("MPValgo") != std::string::npos) { char* MPValgochar = argv[i]; char* MPValgochar2 = strtok (MPValgochar, "M"); MPValgo = atoi(MPValgochar2); cout<<"MPValgo = "<<MPValgo<<endl; }
		if(std::string(argv[i]).find("JobID") != std::string::npos) { char* JobIDchar = argv[i]; char* JobIDchar2 = strtok (JobIDchar, "="); JobID = JobIDchar2; cout<< "JobID = " << JobID << endl; }

  	}


	char inname[200];
	char outname[200];
	char jobdirname[200];
	sprintf(jobdirname,"%s/JobID", storagedir);
	gSystem->mkdir(jobdirname);
	sprintf(jobdirname,"%s/%s",jobdirname,JobID);
	gSystem->mkdir(jobdirname);

	imatrix FreeParam_Fractions;
	ivector FreeParam_Fractions_States;
	imatrix FreeParam_Np_BR;
	imatrix FreeParam_Np_US;
	imatrix FreeParam_consts_star;
	ivector int_Sample_Np;
	ivector int_Sample_Np_consts_star;
	ivector int_Minimizer;


	sprintf(inname,"%s/FreeParam.txt",jobdirname);
	cout<<"read in FreeParam from "<<inname<<endl;

	ifstream in;
	in.open(inname);

    in >> FreeParam_Fractions;
    cout<<"FreeParam_Fractions"<<endl;
    cout<<FreeParam_Fractions<<endl;
	in >> FreeParam_Np_BR;
    cout<<"FreeParam_Np_BR"<<endl;
    cout<<FreeParam_Np_BR<<endl;
	in >> FreeParam_consts_star;
    cout<<"FreeParam_consts_star"<<endl;
    cout<<FreeParam_consts_star<<endl;
    in >> FreeParam_Fractions_States;
    cout<<"FreeParam_Fractions_States"<<endl;
    cout<<FreeParam_Fractions_States<<endl;
	in >> int_Sample_Np;
    cout<<"int_Sample_Np"<<endl;
    cout<<int_Sample_Np<<endl;
	in >> int_Sample_Np_consts_star;
    cout<<"int_Sample_Np_consts_star"<<endl;
    cout<<int_Sample_Np_consts_star<<endl;
	in >> int_Minimizer;
    cout<<"int_Minimizer"<<endl;
    cout<<int_Minimizer<<endl;
	in >> FreeParam_Np_US;
    cout<<"FreeParam_Np_US"<<endl;
    cout<<FreeParam_Np_US<<endl;

    in.close();

    bool SampleNp=false;
    bool SampleNp_consts_star=false;
    if(int_Sample_Np[0]==1) SampleNp=true;
    if(int_Sample_Np_consts_star[0]==1) SampleNp_consts_star=true;

	sprintf(inname,"%s/results.root",jobdirname);
	TFile *ResultsFile = new TFile(inname, "READ");

  	TTree*  outputTreeAllSamplings = (TTree*) ResultsFile->Get("AllSamplings"); // tree filled in all samplings after burnin


  	double DummyVal=0.;

  	char branch_name[200];
  	char hist_name[200];
  	char hist_var_name[200];
  	char projectchar[200];
  	char selectchar[200];

  	double expandMinMaxBy=0.01;


  	double buff_fi_MPV;
  	double buff_fi_errlow;
  	double buff_fi_errhign;

  	//double fi_sampling[NRQCDvars::nStates][NRQCDvars::nColorChannels];
  	int nBins_h_fi=100;
  	dmatrix fi_MPV(NRQCDvars::nStates);
  	dvector fi_MPV_state(NRQCDvars::nColorChannels);
  	dmatrix fi_errlow(NRQCDvars::nStates);
  	dvector fi_errlow_state(NRQCDvars::nColorChannels);
 	dmatrix fi_errhigh(NRQCDvars::nStates);
  	dvector fi_errhigh_state(NRQCDvars::nColorChannels);
  	double h_fi_min[NRQCDvars::nStates][NRQCDvars::nColorChannels];
  	double h_fi_max[NRQCDvars::nStates][NRQCDvars::nColorChannels];
	TH1D* h_fi[NRQCDvars::nStates][NRQCDvars::nColorChannels];

	for (int i=0; i<NRQCDvars::nStates; i++){
		int nColorChannels_state;
		bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
		if(isSstate) nColorChannels_state=NRQCDvars::nColorChannels_S;
		else nColorChannels_state=NRQCDvars::nColorChannels_P;
		for (int j=0; j<nColorChannels_state; j++){
				sprintf(branch_name,"state%d_f%d",i,j);
				//outputTreeAccSamplings->SetBranchAddress( branch_name,  &fi_sampling[i][j] );
				h_fi_min[i][j]=outputTreeAllSamplings->GetMinimum(branch_name)-expandMinMaxBy*outputTreeAllSamplings->GetMinimum(branch_name);
				h_fi_max[i][j]=outputTreeAllSamplings->GetMaximum(branch_name)+expandMinMaxBy*outputTreeAllSamplings->GetMaximum(branch_name);
				sprintf(hist_name,"h_state%d_f%d",i,j);
				h_fi[i][j] = new TH1D( hist_name, hist_name, nBins_h_fi, h_fi_min[i][j], h_fi_max[i][j] );
				if(isSstate) sprintf(hist_var_name,"f_{%s}^{%s}", ColorChannelNameTexS[j], StateNameTex[i]);
				else sprintf(hist_var_name,"f_{%s}^{%s}", ColorChannelNameTexP[j], StateNameTex[i]);
				if(j==0) sprintf(hist_var_name,"R^{%s}", StateNameTex[i]);
				h_fi[i][j] -> SetXTitle(hist_var_name);
				sprintf(projectchar,"state%d_f%d>>%s",i,j,hist_name);
				sprintf(selectchar,"acceptedSampling==1 && BurnInInt==0");
				cout<<projectchar<<" , "<<selectchar<<endl;
				outputTreeAllSamplings->Draw(projectchar, selectchar);
				h_fi[i][j]->SetTitle(0);
				if(FreeParam_Fractions[i][j]==1) FindMPV(h_fi[i][j], buff_fi_MPV, buff_fi_errlow, buff_fi_errhign, MPValgo, nSigma);
				else{buff_fi_MPV=DummyVal; buff_fi_errlow=DummyVal; buff_fi_errhign=DummyVal; }
				fi_MPV_state.at(j)=buff_fi_MPV;
				fi_errlow_state.at(j)=buff_fi_errlow;
				fi_errhigh_state.at(j)=buff_fi_errhign;
				char filename[200];
				sprintf(filename,"%s/Figures",jobdirname);
				gSystem->mkdir(filename);
				sprintf(filename,"%s/PPD_%s.pdf",filename, branch_name);
				if(FreeParam_Fractions[i][j]==1) PlotPosterior(hist_var_name, filename, h_fi[i][j], buff_fi_MPV, buff_fi_errlow, buff_fi_errhign);
			}
		fi_MPV.at(i)=fi_MPV_state;
		fi_errlow.at(i)=fi_errlow_state;
		fi_errhigh.at(i)=fi_errhigh_state;
		}

	dmatrix fi_MPV_dmatrix = transformDmatrixToSPformat(fi_MPV);
	dmatrix fi_errlow_dmatrix = transformDmatrixToSPformat(fi_errlow);
	dmatrix fi_errhigh_dmatrix = transformDmatrixToSPformat(fi_errhigh);





  	double buff_Opi_MPV;
  	double buff_Opi_errlow;
  	double buff_Opi_errhign;

  	//double Opi_sampling[NRQCDvars::nStates][NRQCDvars::nColorChannels];
  	int nBins_h_Opi=100;
  	dmatrix Opi_MPV(NRQCDvars::nStates);
  	dvector Opi_MPV_state(NRQCDvars::nColorChannels);
  	dmatrix Opi_errlow(NRQCDvars::nStates);
  	dvector Opi_errlow_state(NRQCDvars::nColorChannels);
 	dmatrix Opi_errhigh(NRQCDvars::nStates);
  	dvector Opi_errhigh_state(NRQCDvars::nColorChannels);
  	double h_Opi_min[NRQCDvars::nStates][NRQCDvars::nColorChannels];
  	double h_Opi_max[NRQCDvars::nStates][NRQCDvars::nColorChannels];
	TH1D* h_Opi[NRQCDvars::nStates][NRQCDvars::nColorChannels];

	for (int i=0; i<NRQCDvars::nStates; i++){
		int nColorChannels_state;
		bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
		if(isSstate) nColorChannels_state=NRQCDvars::nColorChannels_S;
		else nColorChannels_state=NRQCDvars::nColorChannels_P;
		for (int j=0; j<nColorChannels_state; j++){
				sprintf(branch_name,"state%d_Op%d",i,j);
				//outputTreeAccSamplings->SetBranchAddress( branch_name,  &Opi_sampling[i][j] );
				h_Opi_min[i][j]=outputTreeAllSamplings->GetMinimum(branch_name)-expandMinMaxBy*outputTreeAllSamplings->GetMinimum(branch_name);
				h_Opi_max[i][j]=outputTreeAllSamplings->GetMaximum(branch_name)+expandMinMaxBy*outputTreeAllSamplings->GetMaximum(branch_name);
				sprintf(hist_name,"h_state%d_Op%d",i,j);
				h_Opi[i][j] = new TH1D( hist_name, hist_name, nBins_h_Opi, h_Opi_min[i][j], h_Opi_max[i][j] );
				if(isSstate) sprintf(hist_var_name,"O_{%s}^{%s}", ColorChannelNameTexS[j], StateNameTex[i]);
				else sprintf(hist_var_name,"O_{%s}^{%s}", ColorChannelNameTexP[j], StateNameTex[i]);
				h_Opi[i][j] -> SetXTitle(hist_var_name);
				sprintf(projectchar,"state%d_Op%d>>%s",i,j,hist_name);
				sprintf(selectchar,"acceptedSampling==1 && BurnInInt==0");
				outputTreeAllSamplings->Draw(projectchar, selectchar);
				h_Opi[i][j]->SetTitle(0);
				if(FreeParam_Fractions[i][j]==1) FindMPV(h_Opi[i][j], buff_Opi_MPV, buff_Opi_errlow, buff_Opi_errhign, MPValgo, nSigma);
				else{buff_Opi_MPV=DummyVal; buff_Opi_errlow=DummyVal; buff_Opi_errhign=DummyVal; }
				Opi_MPV_state.at(j)=buff_Opi_MPV;
				Opi_errlow_state.at(j)=buff_Opi_errlow;
				Opi_errhigh_state.at(j)=buff_Opi_errhign;
				char filename[200];
				sprintf(filename,"%s/Figures",jobdirname);
				gSystem->mkdir(filename);
				sprintf(filename,"%s/PPD_%s.pdf",filename, branch_name);
				if(FreeParam_Fractions[i][j]==1) PlotPosterior(hist_var_name, filename, h_Opi[i][j], buff_Opi_MPV, buff_Opi_errlow, buff_Opi_errhign);
			}
		Opi_MPV.at(i)=Opi_MPV_state;
		Opi_errlow.at(i)=Opi_errlow_state;
		Opi_errhigh.at(i)=Opi_errhigh_state;
		}

	dmatrix Opi_MPV_dmatrix = transformDmatrixToSPformat(Opi_MPV);
	dmatrix Opi_errlow_dmatrix = transformDmatrixToSPformat(Opi_errlow);
	dmatrix Opi_errhigh_dmatrix = transformDmatrixToSPformat(Opi_errhigh);







	sprintf(outname,"%s/results.txt",jobdirname);
	cout<<"save results to "<<outname<<endl;

    ofstream out;
    out.open(outname);//, std::ofstream::app);

	out << fi_MPV_dmatrix;
	out << fi_errlow_dmatrix;
	out << fi_errhigh_dmatrix;
	out << Opi_MPV_dmatrix;
	out << Opi_errlow_dmatrix;
	out << Opi_errhigh_dmatrix;

    out.close();

    cout<<"fi_MPV:"<<endl;
	cout<<fi_MPV_dmatrix<<endl;
    cout<<"fi_errlow:"<<endl;
	cout<<fi_errlow_dmatrix<<endl;
    cout<<"fi_errhigh:"<<endl;
	cout<<fi_errhigh_dmatrix<<endl;

    cout<<"Opi_MPV:"<<endl;
	cout<<Opi_MPV_dmatrix<<endl;
    cout<<"Opi_errlow:"<<endl;
	cout<<Opi_errlow_dmatrix<<endl;
    cout<<"Opi_errhigh:"<<endl;
	cout<<Opi_errhigh_dmatrix<<endl;








  	double buff_consts_star_var_MPV;
  	double buff_consts_star_var_errlow;
  	double buff_consts_star_var_errhign;

  	int nBins_h_consts_star_var=100;
  	dmatrix consts_star_var_MPV(NRQCDvars::nStates);
  	dvector consts_star_var_MPV_state(NRQCDvars::nColorChannels);
  	dmatrix consts_star_var_errlow(NRQCDvars::nStates);
  	dvector consts_star_var_errlow_state(NRQCDvars::nColorChannels);
 	dmatrix consts_star_var_errhigh(NRQCDvars::nStates);
  	dvector consts_star_var_errhigh_state(NRQCDvars::nColorChannels);
  	double h_consts_star_var_min[NRQCDvars::nStates][NRQCDvars::nColorChannels];
  	double h_consts_star_var_max[NRQCDvars::nStates][NRQCDvars::nColorChannels];
	TH1D* h_consts_star_var[NRQCDvars::nStates][NRQCDvars::nColorChannels];

	for (int i=0; i<NRQCDvars::nStates; i++){
		int nColorChannels_state;
		bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
		if(isSstate) nColorChannels_state=NRQCDvars::nColorChannels_S;
		else nColorChannels_state=NRQCDvars::nColorChannels_P;
		for (int j=0; j<nColorChannels_state; j++){
				sprintf(branch_name,"state%d_const_star%d",i,j);
				//outputTreeAccSamplings->SetBranchAddress( branch_name,  &consts_star_var_sampling[i][j] );
				h_consts_star_var_min[i][j]=outputTreeAllSamplings->GetMinimum(branch_name)-expandMinMaxBy*outputTreeAllSamplings->GetMinimum(branch_name);
				h_consts_star_var_max[i][j]=outputTreeAllSamplings->GetMaximum(branch_name)+expandMinMaxBy*outputTreeAllSamplings->GetMaximum(branch_name);
				sprintf(hist_name,"h_state%d_consts_star%d",i,j);
				h_consts_star_var[i][j] = new TH1D( hist_name, hist_name, nBins_h_consts_star_var, h_consts_star_var_min[i][j], h_consts_star_var_max[i][j] );
				if(isSstate) sprintf(hist_var_name,"c*_{%s}^{%s}", ColorChannelNameTexS[j], StateNameTex[i]);
				else sprintf(hist_var_name,"c*_{%s}^{%s}", ColorChannelNameTexP[j], StateNameTex[i]);
				if(j==0) sprintf(hist_var_name,"#sigma_{s}*^{%s}", StateNameTex[i]);
				h_consts_star_var[i][j] -> SetXTitle(hist_var_name);
				sprintf(projectchar,"%s>>%s",branch_name,hist_name);
				sprintf(selectchar,"acceptedSampling==1 && BurnInInt==0");
				outputTreeAllSamplings->Draw(projectchar, selectchar);
				h_consts_star_var[i][j]->SetTitle(0);
				if(FreeParam_consts_star[i][j]==1 && SampleNp_consts_star) FindMPV(h_consts_star_var[i][j], buff_consts_star_var_MPV, buff_consts_star_var_errlow, buff_consts_star_var_errhign, MPValgo, nSigma);
				else{buff_consts_star_var_MPV=DummyVal; buff_consts_star_var_errlow=DummyVal; buff_consts_star_var_errhign=DummyVal; }
				consts_star_var_MPV_state.at(j)=buff_consts_star_var_MPV;
				consts_star_var_errlow_state.at(j)=buff_consts_star_var_errlow;
				consts_star_var_errhigh_state.at(j)=buff_consts_star_var_errhign;
				char filename[200];
				sprintf(filename,"%s/Figures",jobdirname);
				gSystem->mkdir(filename);
				sprintf(filename,"%s/Np_%s.pdf",filename, branch_name);
				if(FreeParam_consts_star[i][j]==1 && SampleNp_consts_star) PlotPosterior(hist_var_name, filename, h_consts_star_var[i][j], buff_consts_star_var_MPV, buff_consts_star_var_errlow, buff_consts_star_var_errhign);
			}
		consts_star_var_MPV.at(i)=consts_star_var_MPV_state;
		consts_star_var_errlow.at(i)=consts_star_var_errlow_state;
		consts_star_var_errhigh.at(i)=consts_star_var_errhigh_state;
		}



  	double buff_Np_BR_MPV;
  	double buff_Np_BR_errlow;
  	double buff_Np_BR_errhign;

  	int nBins_h_Np_BR=100;
  	dmatrix Np_BR_MPV;
  	dvector Np_BR_MPV_0 (NRQCDvars::nStates,0.);
	for(int i=0; i < NRQCDvars::nStates; i++) Np_BR_MPV.push_back(Np_BR_MPV_0);
  	dmatrix Np_BR_errlow;
  	dvector Np_BR_errlow_0 (NRQCDvars::nStates,0.);
	for(int i=0; i < NRQCDvars::nStates; i++) Np_BR_errlow.push_back(Np_BR_errlow_0);
 	dmatrix Np_BR_errhigh;
  	dvector Np_BR_errhigh_0 (NRQCDvars::nStates,0.);
	for(int i=0; i < NRQCDvars::nStates; i++) Np_BR_errhigh.push_back(Np_BR_errhigh_0);
  	double h_Np_BR_min[NRQCDvars::nStates][NRQCDvars::nStates];
  	double h_Np_BR_max[NRQCDvars::nStates][NRQCDvars::nStates];
	TH1D* h_Np_BR[NRQCDvars::nStates][NRQCDvars::nStates];

	for (int i=0; i<NRQCDvars::nStates; i++){
		for (int j=0; j<NRQCDvars::nStates; j++){
			if(NRQCDvars::FeedDownBranchingRatio[i][j]>0){
				sprintf(branch_name,"Np_BR_Daughter%d_Mother%d",i,j);
				h_Np_BR_min[i][j]=outputTreeAllSamplings->GetMinimum(branch_name)-expandMinMaxBy*outputTreeAllSamplings->GetMinimum(branch_name);
				h_Np_BR_max[i][j]=outputTreeAllSamplings->GetMaximum(branch_name)+expandMinMaxBy*outputTreeAllSamplings->GetMaximum(branch_name);
				sprintf(hist_name,"h_Np_BR_Daughter%d_Mother%d",i,j);
				h_Np_BR[i][j] = new TH1D( hist_name, hist_name, nBins_h_Np_BR, h_Np_BR_min[i][j], h_Np_BR_max[i][j] );
				sprintf(hist_var_name,"BR( %s -> %s )", StateNameTex[j], StateNameTex[i]);
				h_Np_BR[i][j] -> SetXTitle(hist_var_name);
				sprintf(projectchar,"%s>>%s",branch_name,hist_name);
				sprintf(selectchar,"acceptedSampling==1 && BurnInInt==0");
				outputTreeAllSamplings->Draw(projectchar, selectchar);
				h_Np_BR[i][j]->SetTitle(0);
				if(FreeParam_Np_BR[i][j]==1 && SampleNp) FindMPV(h_Np_BR[i][j], buff_Np_BR_MPV, buff_Np_BR_errlow, buff_Np_BR_errhign, MPValgo, nSigma);
				else{buff_Np_BR_MPV=DummyVal; buff_Np_BR_errlow=DummyVal; buff_Np_BR_errhign=DummyVal; }
				Np_BR_MPV[i][j]=buff_Np_BR_MPV;
				Np_BR_errlow[i][j]=buff_Np_BR_errlow;
				Np_BR_errhigh[i][j]=buff_Np_BR_errhign;
				char filename[200];
				sprintf(filename,"%s/Figures",jobdirname);
				gSystem->mkdir(filename);
				sprintf(filename,"%s/%s.pdf",filename, branch_name);
				if(FreeParam_Np_BR[i][j]==1 && SampleNp) PlotPosterior(hist_var_name, filename, h_Np_BR[i][j], buff_Np_BR_MPV, buff_Np_BR_errlow, buff_Np_BR_errhign);
				cout<<branch_name<<endl;
				cout<<buff_Np_BR_MPV<<endl;
			}
			else{
				buff_Np_BR_MPV=DummyVal; buff_Np_BR_errlow=DummyVal; buff_Np_BR_errhign=DummyVal;
				Np_BR_MPV[i][j]=buff_Np_BR_MPV;
				Np_BR_errlow[i][j]=buff_Np_BR_errlow;
				Np_BR_errhigh[i][j]=buff_Np_BR_errhign;
				sprintf(branch_name,"Np_BR_Daughter%d_Mother%d",i,j);
				cout<<branch_name<<endl;
				cout<<buff_Np_BR_MPV<<endl;
			}
		}
	}



  	double buff_Np_US_MPV;
  	double buff_Np_US_errlow;
  	double buff_Np_US_errhign;

  	int nBins_h_Np_US=100;

  	int nScales=2;
  	if(NRQCDvars::nModelSystematicScales==0) nScales=1;
  	dmatrix Np_US_MPV(nScales);
  	dmatrix Np_US_errlow(nScales);
 	dmatrix Np_US_errhigh(nScales);
  	double h_Np_US_min[2][max(NRQCDvars::nDataSystematicScales,NRQCDvars::nModelSystematicScales)];
  	double h_Np_US_max[2][max(NRQCDvars::nDataSystematicScales,NRQCDvars::nModelSystematicScales)];
	TH1D* h_Np_US[2][max(NRQCDvars::nDataSystematicScales,NRQCDvars::nModelSystematicScales)];


  	dvector Np_US_MPV_Data_vec(NRQCDvars::nDataSystematicScales);
  	dvector Np_US_errlow_Data_state(NRQCDvars::nDataSystematicScales);
  	dvector Np_US_errhigh_Data_state(NRQCDvars::nDataSystematicScales);

	int i;
	i=0;
	for(int j=0; j < NRQCDvars::nDataSystematicScales; j++){
		sprintf(branch_name,"Np_US_DataSystematicScale%d",j);
		h_Np_US_min[i][j]=outputTreeAllSamplings->GetMinimum(branch_name)-expandMinMaxBy*outputTreeAllSamplings->GetMinimum(branch_name);
		h_Np_US_max[i][j]=outputTreeAllSamplings->GetMaximum(branch_name)+expandMinMaxBy*outputTreeAllSamplings->GetMaximum(branch_name);
		sprintf(hist_name,"h_Np_US_DataSystematicScale%d",j);
		h_Np_US[i][j] = new TH1D( hist_name, hist_name, nBins_h_Np_US, h_Np_US_min[i][j], h_Np_US_max[i][j] );
		sprintf(hist_var_name,"LumiCorrFactor(%s)", ExpName[j]);
		h_Np_US[i][j] -> SetXTitle(hist_var_name);
		sprintf(projectchar,"%s>>%s",branch_name,hist_name);
		sprintf(selectchar,"acceptedSampling==1 && BurnInInt==0");
		outputTreeAllSamplings->Draw(projectchar, selectchar);
		h_Np_US[i][j]->SetTitle(0);
		if(FreeParam_Np_US[i][j]==1 && SampleNp) FindMPV(h_Np_US[i][j], buff_Np_US_MPV, buff_Np_US_errlow, buff_Np_US_errhign, MPValgo, nSigma);
		else{buff_Np_US_MPV=DummyVal; buff_Np_US_errlow=DummyVal; buff_Np_US_errhign=DummyVal; }
		Np_US_MPV_Data_vec.at(j)=buff_Np_US_MPV;
		Np_US_errlow_Data_state.at(j)=buff_Np_US_errlow;
		Np_US_errhigh_Data_state.at(j)=buff_Np_US_errhign;
		char filename[200];
		sprintf(filename,"%s/Figures",jobdirname);
		gSystem->mkdir(filename);
		sprintf(filename,"%s/%s.pdf",filename, branch_name);
		if(FreeParam_Np_US[i][j]==1 && SampleNp) PlotPosterior(hist_var_name, filename, h_Np_US[i][j], buff_Np_US_MPV, buff_Np_US_errlow, buff_Np_US_errhign);
	}
	Np_US_MPV.at(i)=Np_US_MPV_Data_vec;
	Np_US_errlow.at(i)=Np_US_errlow_Data_state;
	Np_US_errhigh.at(i)=Np_US_errhigh_Data_state;

  	dvector Np_US_MPV_Model_vec(NRQCDvars::nModelSystematicScales);
  	dvector Np_US_errlow_Model_state(NRQCDvars::nModelSystematicScales);
  	dvector Np_US_errhigh_Model_state(NRQCDvars::nModelSystematicScales);

	i=1;
	for(int j=0; j < NRQCDvars::nModelSystematicScales; j++){
		sprintf(branch_name,"Np_US_ModelSystematicScale%d",j);
		h_Np_US_min[i][j]=outputTreeAllSamplings->GetMinimum(branch_name)-expandMinMaxBy*outputTreeAllSamplings->GetMinimum(branch_name);
		h_Np_US_max[i][j]=outputTreeAllSamplings->GetMaximum(branch_name)+expandMinMaxBy*outputTreeAllSamplings->GetMaximum(branch_name);
		sprintf(hist_name,"h_Np_US_ModelSystematicScale%d",j);
		h_Np_US[i][j] = new TH1D( hist_name, hist_name, nBins_h_Np_US, h_Np_US_min[i][j], h_Np_US_max[i][j] );
		sprintf(hist_var_name,"ModelSystScale%d", j);
		h_Np_US[i][j] -> SetXTitle(hist_var_name);
		sprintf(projectchar,"%s>>%s",branch_name,hist_name);
		sprintf(selectchar,"acceptedSampling==1 && BurnInInt==0");
		outputTreeAllSamplings->Draw(projectchar, selectchar);
		h_Np_US[i][j]->SetTitle(0);
		if(FreeParam_Np_US[i][j]==1 && SampleNp) FindMPV(h_Np_US[i][j], buff_Np_US_MPV, buff_Np_US_errlow, buff_Np_US_errhign, MPValgo, nSigma);
		else{buff_Np_US_MPV=DummyVal; buff_Np_US_errlow=DummyVal; buff_Np_US_errhign=DummyVal; }
		Np_US_MPV_Model_vec.at(j)=buff_Np_US_MPV;
		Np_US_errlow_Model_state.at(j)=buff_Np_US_errlow;
		Np_US_errhigh_Model_state.at(j)=buff_Np_US_errhign;
		char filename[200];
		sprintf(filename,"%s/Figures",jobdirname);
		gSystem->mkdir(filename);
		sprintf(filename,"%s/%s.pdf",filename, branch_name);
		if(FreeParam_Np_US[i][j]==1 && SampleNp) PlotPosterior(hist_var_name, filename, h_Np_US[i][j], buff_Np_US_MPV, buff_Np_US_errlow, buff_Np_US_errhign);
	}

	if(NRQCDvars::nModelSystematicScales!=0) {
		Np_US_MPV.at(i)=Np_US_MPV_Model_vec;
		Np_US_errlow.at(i)=Np_US_errlow_Model_state;
		Np_US_errhigh.at(i)=Np_US_errhigh_Model_state;
	}



	sprintf(outname,"%s/results_Np.txt",jobdirname);
	cout<<"save results to "<<outname<<endl;

    out.open(outname);//, std::ofstream::app);

    cout<<"Np_BR:"<<endl;
	out << Np_BR_MPV;
    cout<<"errlow_Np_BR:"<<endl;
	out << Np_BR_errlow;
    cout<<"errlow_Np_BR:"<<endl;
	out << Np_BR_errhigh;
    cout<<"consts_star_var:"<<endl;
	out << consts_star_var_MPV;
    cout<<"errlow_consts_star_var:"<<endl;
	out << consts_star_var_errlow;
    cout<<"errhigh_consts_star_var:"<<endl;
	out << consts_star_var_errhigh;
    cout<<"Np_US:"<<endl;
	out << Np_US_MPV;
    cout<<"errlow_Np_US:"<<endl;
	out << Np_US_errlow;
    cout<<"errhigh_Np_US:"<<endl;
	out << Np_US_errhigh;

    out.close();

	//sprintf(outname,"%s/results3.txt",jobdirname);
	//cout<<"save test results to "<<outname<<endl;
    //
    //ofstream out2;
    //out2.open(outname);//, std::ofstream::app);
    //
	//out2 << fi_MPV_dmatrix;
    //
    //out2.close();
    //
    //cout << fi_MPV_dmatrix << endl;



		return 0;

  	}







void FindMPV(TH1* PosteriorDist , double& MPV , double& MPVerrorLow, double& MPVerrorHigh, int MPValgo, int nSigma){

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
			if(gauss->GetChisquare()/ndof<5) break;
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

void PlotPosterior(char xAxisTitle[200], char filename[200], TH1* histo, double MPV, double MPVerrLow, double MPVerrHigh){

	gStyle->SetPalette(1,0);
	gStyle->SetPadBottomMargin(0.12);
	gStyle->SetPadLeftMargin(0.13);
	gStyle->SetPadRightMargin(0.15);

	gStyle->SetTickLength(-0.02, "xyz");
	gStyle->SetLabelOffset(0.02, "x");
	gStyle->SetLabelOffset(0.02, "y");
	gStyle->SetTitleOffset(1.3, "x");
	gStyle->SetTitleOffset(1.4, "y");
	gStyle->SetTitleFillColor(kWhite);

	TLegend* plotLegend=new TLegend(0.775,0.8,0.95,0.9);
	plotLegend->SetFillColor(kWhite);
	plotLegend->SetTextFont(72);
	plotLegend->SetTextSize(0.0175);
	plotLegend->SetBorderSize(1);

	TCanvas *plotCanvas = new TCanvas("plotCanvas","plotCanvas",1000,800);

	plotCanvas->SetFillColor(kWhite);
	//	plotCanvas->SetGrid();
	plotCanvas->GetFrame()->SetFillColor(kWhite);
	plotCanvas->GetFrame()->SetBorderSize(0);
	plotCanvas->SetRightMargin(0.05) ;

	histo->SetStats(kFALSE);
	histo->SetLineColor(kBlack);
	histo->SetYTitle("Posterior Probability");
	histo->SetXTitle(xAxisTitle);
	histo->GetYaxis()->SetTitleOffset(1.5);
	histo->Draw();


	int maxbin_PosteriorDist = histo->GetMaximumBin();
	double maxval = histo->GetBinContent(maxbin_PosteriorDist);

	TLine* MeanLine = new TLine( histo->GetMean(), 0, histo->GetMean(), maxval );
	MeanLine->SetLineWidth( 1 );
	MeanLine->SetLineStyle( 1 );
	MeanLine->SetLineColor( kGreen+2 );
	MeanLine->Draw( "same" );
	TLine* MPVLine = new TLine( MPV, 0, MPV, maxval );
	MPVLine->SetLineWidth( 2.5 );
	MPVLine->SetLineStyle( 1 );
	MPVLine->SetLineColor( kRed );
	MPVLine->Draw( "same" );
	TLine* MPVerrLowLine = new TLine( MPV-MPVerrLow, 0, MPV-MPVerrLow, maxval );
	MPVerrLowLine->SetLineWidth( 2.5 );
	MPVerrLowLine->SetLineStyle( 2 );
	MPVerrLowLine->SetLineColor( kRed );
	MPVerrLowLine->Draw( "same" );
	TLine* MPVerrHighLine = new TLine( MPV+MPVerrHigh, 0, MPV+MPVerrHigh, maxval );
	MPVerrHighLine->SetLineWidth( 2.5 );
	MPVerrHighLine->SetLineStyle( 2 );
	MPVerrHighLine->SetLineColor( kRed );
	MPVerrHighLine->Draw( "same" );
	plotLegend->AddEntry(MeanLine,"Mean","l");
	plotLegend->AddEntry(MPVLine,"MPV","l");
	plotLegend->AddEntry(MPVerrLowLine,"1#sigma high/low","l");

	plotLegend->Draw();
	plotCanvas->SaveAs(filename);
	plotCanvas->Close();

	delete plotCanvas;

}

dmatrix transformDmatrixToSPformat(dmatrix dmatrixFull){

	dmatrix transformed_dmatrix(NRQCDvars::nStates);
	vector<double> transformed_dmatrix_S (NRQCDvars::nColorChannels_S,0);
	vector<double> transformed_dmatrix_P (NRQCDvars::nColorChannels_P,0);

	for (int i=0; i<NRQCDvars::nStates; i++){
		bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
		if(isSstate){
			for (int j=0; j<NRQCDvars::nColorChannels_S; j++){
				transformed_dmatrix_S.at(j)=dmatrixFull[i][j];
			}
			transformed_dmatrix.at(i)=transformed_dmatrix_S;
		}
		else{
			for (int j=0; j<NRQCDvars::nColorChannels_P; j++){
				transformed_dmatrix_P.at(j)=dmatrixFull[i][j];
			}
			transformed_dmatrix.at(i)=transformed_dmatrix_P;
		}
	}

	return transformed_dmatrix;

}

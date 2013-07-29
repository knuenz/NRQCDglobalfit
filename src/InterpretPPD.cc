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
	sprintf(inname,"%s/results.root",jobdirname);

	TFile *ResultsFile = new TFile(inname, "READ");

  	TTree*  outputTreeAllSamplings = (TTree*) ResultsFile->Get("AllSamplings"); // tree filled in all samplings after burnin
  	TTree*  outputTreeAccSamplings = (TTree*) ResultsFile->Get("AccSamplings"); // tree filled only in accepted samplings after burnin


  	char branch_name[200];
  	char hist_name[200];
  	char hist_var_name[200];
  	char projectchar[200];

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
				h_fi_min[i][j]=outputTreeAccSamplings->GetMinimum(branch_name)-expandMinMaxBy*outputTreeAccSamplings->GetMinimum(branch_name);
				h_fi_max[i][j]=outputTreeAccSamplings->GetMaximum(branch_name)+expandMinMaxBy*outputTreeAccSamplings->GetMaximum(branch_name);
				sprintf(hist_name,"h_state%d_f%d",i,j);
				h_fi[i][j] = new TH1D( hist_name, hist_name, nBins_h_fi, h_fi_min[i][j], h_fi_max[i][j] );
				if(isSstate) sprintf(hist_var_name,"f_{%s}^{%s}", ColorChannelNameTexS[j], StateNameTex[i]);
				else sprintf(hist_var_name,"f_{%s}^{%s}", ColorChannelNameTexP[j], StateNameTex[i]);
				if(j==0) sprintf(hist_var_name,"R^{%s}", StateNameTex[i]);
				h_fi[i][j] -> SetXTitle(hist_var_name);
				sprintf(projectchar,"state%d_f%d>>%s",i,j,hist_name);
				outputTreeAccSamplings->Draw(projectchar);
				h_fi[i][j]->SetTitle(0);
				FindMPV(h_fi[i][j], buff_fi_MPV, buff_fi_errlow, buff_fi_errhign, MPValgo, nSigma);
				fi_MPV_state.at(j)=buff_fi_MPV;
				fi_errlow_state.at(j)=buff_fi_errlow;
				fi_errhigh_state.at(j)=buff_fi_errhign;
				char filename[200];
				sprintf(filename,"%s/Figures",jobdirname);
				gSystem->mkdir(filename);
				sprintf(filename,"%s/PPD_%s.pdf",filename, branch_name);
				PlotPosterior(hist_var_name, filename, h_fi[i][j], buff_fi_MPV, buff_fi_errlow, buff_fi_errhign);
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
				h_Opi_min[i][j]=outputTreeAccSamplings->GetMinimum(branch_name)-expandMinMaxBy*outputTreeAccSamplings->GetMinimum(branch_name);
				h_Opi_max[i][j]=outputTreeAccSamplings->GetMaximum(branch_name)+expandMinMaxBy*outputTreeAccSamplings->GetMaximum(branch_name);
				sprintf(hist_name,"h_state%d_Op%d",i,j);
				h_Opi[i][j] = new TH1D( hist_name, hist_name, nBins_h_Opi, h_Opi_min[i][j], h_Opi_max[i][j] );
				if(isSstate) sprintf(hist_var_name,"O_{%s}^{%s}", ColorChannelNameTexS[j], StateNameTex[i]);
				else sprintf(hist_var_name,"O_{%s}^{%s}", ColorChannelNameTexP[j], StateNameTex[i]);
				h_Opi[i][j] -> SetXTitle(hist_var_name);
				sprintf(projectchar,"state%d_Op%d>>%s",i,j,hist_name);
				outputTreeAccSamplings->Draw(projectchar);
				h_Opi[i][j]->SetTitle(0);
				FindMPV(h_Opi[i][j], buff_Opi_MPV, buff_Opi_errlow, buff_Opi_errhign, MPValgo, nSigma);
				Opi_MPV_state.at(j)=buff_Opi_MPV;
				Opi_errlow_state.at(j)=buff_Opi_errlow;
				Opi_errhigh_state.at(j)=buff_Opi_errhign;
				char filename[200];
				sprintf(filename,"%s/Figures",jobdirname);
				gSystem->mkdir(filename);
				sprintf(filename,"%s/PPD_%s.pdf",filename, branch_name);
				PlotPosterior(hist_var_name, filename, h_Opi[i][j], buff_Opi_MPV, buff_Opi_errlow, buff_Opi_errhign);
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

/*
 * PlotPPD.cc
 *
 *  Created on: Oct 15, 2013
 *      Author: valentinknuenz
 */



/*
 * FitPtDists.cc
 *
 *  Created on: Oct 10, 2013
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

//rootincludes
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2.h"
#include "TF1.h"
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
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLatex.h"


using namespace NRQCDvars;

inline void setContourHistogram ( TH2D *h );
inline double contourHeight2D ( TH2D *h, double confidenceLevel );

int main(int argc, char** argv) {

  	Char_t *JobID = "Default";
	Char_t *storagedir = "Default"; //Storage Directory

  	for( int i=0;i < argc; ++i ) {
		if(std::string(argv[i]).find("JobID") != std::string::npos) {char* JobIDchar = argv[i]; char* JobIDchar2 = strtok (JobIDchar, "="); JobID = JobIDchar2; cout<<"JobID = "<<JobID<<endl;}
		if(std::string(argv[i]).find("storagedir") != std::string::npos) {char* storagedirchar = argv[i]; char* storagedirchar2 = strtok (storagedirchar, "="); storagedir = storagedirchar2; cout<<"storagedir = "<<storagedir<<endl;}
  	}


  	gROOT->SetBatch();

	delete gRandom;
	gRandom = new TRandom3(NRQCDvars::randomSeed);  // better random generator


	// TODO: add loop for nSystematics (now nSystematics will be 0)

	char outname[200];
	char inname[200];
	char plotdirname[200];
	char jobdirname[200];
	sprintf(plotdirname,"Figures/PlotPPD/%s",JobID);
	gSystem->mkdir(plotdirname);
	sprintf(jobdirname,"%s/JobID/%s",storagedir, JobID);

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


	char JobID2[200];
	char JobID3[200];
	sprintf(JobID2, "October16_BKmodel_ScaledAllStates_psi2S_1S03S1_pTstar_over_m6_MH_pTmin13_POLonly");
	sprintf(JobID3, "October16_BKmodel_ScaledAllStates_psi2S_1S03S1_pTstar_over_m6_MH_pTmin13_CSonly");

	sprintf(inname,"%s/JobID/%s/results.root",storagedir, JobID);
	TFile *ResultsFile = new TFile(inname, "READ");
  	TTree*  outputTreeAllSamplings = (TTree*) ResultsFile->Get("AllSamplings"); // tree filled in all samplings after burnin

	sprintf(inname,"%s/JobID/%s/results.root",storagedir, JobID2);
	TFile *ResultsFile2 = new TFile(inname, "READ");
  	TTree*  outputTreeAllSamplings2 = (TTree*) ResultsFile2->Get("AllSamplings"); // tree filled in all samplings after burnin
	sprintf(inname,"%s/JobID/%s/results.root",storagedir, JobID3);
	TFile *ResultsFile3 = new TFile(inname, "READ");
  	TTree*  outputTreeAllSamplings3 = (TTree*) ResultsFile3->Get("AllSamplings"); // tree filled in all samplings after burnin



	////////////////////////////////////////
	/// Plot 2D contours in fraction space
	////////////////////////////////////////


  	double expandMinMaxBy=0.01;
  	char projectchar[200];
  	char selectchar[200];

		for (int i=0; i<NRQCDvars::nStates; i++){
			if(FreeParam_Fractions_States[i]!=1) continue;

			bool plot2DLDMEComb[NRQCDvars::nColorChannels][NRQCDvars::nColorChannels];
			TH2D* h_2Df[NRQCDvars::nColorChannels][NRQCDvars::nColorChannels];
			TH2D* h_2Df2[NRQCDvars::nColorChannels][NRQCDvars::nColorChannels];
			TH2D* h_2Df3[NRQCDvars::nColorChannels][NRQCDvars::nColorChannels];
			int nBins2D=20;
			double f1Min, f1Max, f2Min, f2Max;
			char branch_name_f1[200];
			char branch_name_f2[200];
			char hist_var_name_f1[200];
			char hist_var_name_f2[200];


			int nColorChannels_state;
			bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
			if(isSstate) nColorChannels_state=NRQCDvars::nColorChannels_S;
			else nColorChannels_state=NRQCDvars::nColorChannels_P;
			for (int j=1; j<nColorChannels_state; j++){
				for (int k=1; k<nColorChannels_state; k++){
					plot2DLDMEComb[j][k]=false;
				}
			}
			for (int j=0; j<nColorChannels_state; j++){
				for (int k=0; k<nColorChannels_state; k++){
					if(j==k) continue;
					plot2DLDMEComb[j][k]=true;
					if(plot2DLDMEComb[k][j]) plot2DLDMEComb[j][k]=false;
					if(plot2DLDMEComb[j][k]){
						cout<<"j"<<j<<"k"<<k<<endl;
						sprintf(branch_name_f1,"state%d_f%d",i,j);
						sprintf(branch_name_f2,"state%d_f%d",i,k);
						//outputTreeAccSamplings->SetBranchAddress( branch_name,  &Opi_sampling[i][j] );
						f1Min=outputTreeAllSamplings->GetMinimum(branch_name_f1)-expandMinMaxBy*outputTreeAllSamplings->GetMinimum(branch_name_f1);
						f1Max=outputTreeAllSamplings->GetMaximum(branch_name_f1)+expandMinMaxBy*outputTreeAllSamplings->GetMaximum(branch_name_f1);
						f2Min=outputTreeAllSamplings->GetMinimum(branch_name_f2)-expandMinMaxBy*outputTreeAllSamplings->GetMinimum(branch_name_f2);
						f2Max=outputTreeAllSamplings->GetMaximum(branch_name_f2)+expandMinMaxBy*outputTreeAllSamplings->GetMaximum(branch_name_f2);

						f1Min=-0.5;
						f1Max=1.5;
						f2Min=-3;
						f2Max=2;

						if(j==1){
							f1Min=-0.5;
							f1Max=1.5;
						}
						if(j==2){
							f1Min=-3.;
							f1Max=3.;
						}
						if(j==3){
							f1Min=-3.;
							f1Max=3.;
						}

						if(k==1){
							f2Min=-0.5;
							f2Max=1.5;
						}
						if(k==2){
							f2Min=-3.;
							f2Max=3.;
						}
						if(k==3){
							f2Min=-3.;
							f2Max=3.;
						}

						if(isSstate){
							sprintf(hist_var_name_f1,"f_{%s}^{%s}", ColorChannelNameTexS[j], StateNameTex[i]);
							sprintf(hist_var_name_f2,"f_{%s}^{%s}", ColorChannelNameTexS[k], StateNameTex[i]);
						}
						else{
							sprintf(hist_var_name_f1,"f_{%s}^{%s}", ColorChannelNameTexP[j], StateNameTex[i]);
							sprintf(hist_var_name_f2,"f_{%s}^{%s}", ColorChannelNameTexP[k], StateNameTex[i]);
						}
						if(j==0) sprintf(hist_var_name_f1,"R^{%s}", StateNameTex[i]);
						if(k==0) sprintf(hist_var_name_f2,"R^{%s}", StateNameTex[i]);

						char h_2Df_name[200];

						sprintf(h_2Df_name,"h_2Df_state%d_f%d_vs_f%d",i,j,k);
						sprintf(projectchar,"%s:%s>>%s",branch_name_f2, branch_name_f1, h_2Df_name);
						sprintf(selectchar,"acceptedSampling==1 && BurnInInt==0");
						h_2Df[j][k] 	= new TH2D( h_2Df_name, h_2Df_name, nBins2D, f1Min, f1Max, nBins2D, f2Min, f2Max );
						outputTreeAllSamplings->Draw(projectchar, selectchar);

						sprintf(h_2Df_name,"h_2Df2_state%d_f%d_vs_f%d",i,j,k);
						sprintf(projectchar,"%s:%s>>%s",branch_name_f2, branch_name_f1, h_2Df_name);
						sprintf(selectchar,"acceptedSampling==1 && BurnInInt==0");
						h_2Df2[j][k] 	= new TH2D( h_2Df_name, h_2Df_name, nBins2D, f1Min, f1Max, nBins2D, f2Min, f2Max );
						outputTreeAllSamplings2->Draw(projectchar, selectchar);

						sprintf(h_2Df_name,"h_2Df3_state%d_f%d_vs_f%d",i,j,k);
						sprintf(projectchar,"%s:%s>>%s",branch_name_f2, branch_name_f1, h_2Df_name);
						sprintf(selectchar,"acceptedSampling==1 && BurnInInt==0");
						h_2Df3[j][k] 	= new TH2D( h_2Df_name, h_2Df_name, nBins2D, f1Min, f1Max, nBins2D, f2Min, f2Max );
						outputTreeAllSamplings3->Draw(projectchar, selectchar);



						  TCanvas *c1 = new TCanvas("c1", "c1", 10, 28, 650,571);
						  c1->Range(-237.541,-66.47556,187.377,434.8609);
						  c1->SetFillColor(0);
						  c1->SetBorderMode(0);
						  c1->SetBorderSize(0);
						  c1->SetLeftMargin(0.215);
						  c1->SetRightMargin(0.03);
						  c1->SetTopMargin(0.01841621);
						  c1->SetBottomMargin(0.16);
						  c1->SetFrameBorderMode(0);


						  TH2D* h_2Df_axis;
						  h_2Df_axis 	= new TH2D( Form("%s_axis",h_2Df_name), Form("%s_axis",h_2Df_name), nBins2D, f1Min, f1Max, nBins2D, f2Min, f2Max );

						  char DrawContourStyle[200];
						  sprintf(DrawContourStyle,"cont2,same");
						  int LineWidth=4;
						  int LineStyle=2;

						  h_2Df_axis->GetXaxis()->SetTitle(hist_var_name_f1);
						  h_2Df_axis->GetXaxis()->SetLabelOffset(0.028);
						  h_2Df_axis->GetXaxis()->SetTitleSize(0.05);
						  h_2Df_axis->GetXaxis()->SetTickLength(-0.03);
						  h_2Df_axis->GetXaxis()->SetTitleOffset(1.4);
						  h_2Df_axis->GetYaxis()->SetTitle(hist_var_name_f2);
						  h_2Df_axis->GetYaxis()->SetLabelOffset(0.032);
						  h_2Df_axis->GetYaxis()->SetTitleSize(0.05);
						  h_2Df_axis->GetYaxis()->SetTickLength(-0.03);
						  h_2Df_axis->GetYaxis()->SetTitleOffset(1.95);
						  h_2Df_axis->SetTitle(0);
						  h_2Df_axis->SetStats(0);
						  h_2Df_axis->Draw("");

						  setContourHistogram ( h_2Df[j][k] );
						  h_2Df[j][k]->SetLineColor( kGreen+2 );
						  h_2Df[j][k]->SetLineWidth( LineWidth );
						  h_2Df[j][k]->SetLineStyle( LineStyle  );
						  h_2Df[j][k]->SetTitle(0);
						  h_2Df[j][k]->SetStats(0);
						  h_2Df[j][k]->Draw( DrawContourStyle );

						  setContourHistogram ( h_2Df2[j][k] );
						  h_2Df2[j][k]->SetLineColor( kBlue );
						  h_2Df2[j][k]->SetLineWidth( LineWidth );
						  h_2Df2[j][k]->SetLineStyle( LineStyle  );
						  h_2Df2[j][k]->SetTitle(0);
						  h_2Df2[j][k]->SetStats(0);
						  h_2Df2[j][k]->Draw( DrawContourStyle );

						  setContourHistogram ( h_2Df3[j][k] );
						  h_2Df3[j][k]->SetLineColor( kRed );
						  h_2Df3[j][k]->SetLineWidth( LineWidth );
						  h_2Df3[j][k]->SetLineStyle( LineStyle  );
						  h_2Df3[j][k]->SetTitle(0);
						  h_2Df3[j][k]->SetStats(0);
						  h_2Df3[j][k]->Draw( DrawContourStyle );


							TLegend* legend;
							legend = new TLegend(0.25,0.20,0.50,0.40);
							legend->SetFillColor(0);
							legend->SetTextSize(0.03);
							legend->SetBorderSize(0);
							legend->AddEntry(h_2Df3[j][k],Form("Cross section constraint"),"l");
							legend->AddEntry(h_2Df2[j][k],Form("Polarization constraint"),"l");
							legend->AddEntry(h_2Df[j][k],Form("Both constraints"),"l");

							legend->Draw("same");

						  char filename[200];
						  sprintf(filename,"Figures/PlotPPD/%s/2D_PPD_contours_comp_state%d_f%d_vs_f%d.pdf",JobID, i,j,k);
						  c1->SaveAs(filename);
						  c1->Close();

						  delete c1;

					}
				}
			}
		}



		////////////////////////////////////////
		/// Plot 1D contours in fraction space
		////////////////////////////////////////


		for (int i=0; i<NRQCDvars::nStates; i++){
			if(FreeParam_Fractions_States[i]!=1) continue;

			TH1D* h_1Df[NRQCDvars::nColorChannels];
			TH1D* h_1Df2[NRQCDvars::nColorChannels];
			TH1D* h_1Df3[NRQCDvars::nColorChannels];
			int nBins1D=50;
			double f1Min, f1Max;
			char branch_name_f1[200];
			char hist_var_name_f1[200];


			int nColorChannels_state;
			bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
			if(isSstate) nColorChannels_state=NRQCDvars::nColorChannels_S;
			else nColorChannels_state=NRQCDvars::nColorChannels_P;
			for (int j=0; j<nColorChannels_state; j++){
						sprintf(branch_name_f1,"state%d_f%d",i,j);
						//outputTreeAccSamplings->SetBranchAddress( branch_name,  &Opi_sampling[i][j] );
						f1Min=outputTreeAllSamplings->GetMinimum(branch_name_f1)-expandMinMaxBy*outputTreeAllSamplings->GetMinimum(branch_name_f1);
						f1Max=outputTreeAllSamplings->GetMaximum(branch_name_f1)+expandMinMaxBy*outputTreeAllSamplings->GetMaximum(branch_name_f1);

						if(j==1){
							f1Min=-0.5;
							f1Max=1.5;
						}
						if(j==2){
							f1Min=-0.5;
							f1Max=1.5;
						}
						if(j==3){
							f1Min=-3.;
							f1Max=3.;
						}

						if(isSstate){
							sprintf(hist_var_name_f1,"f_{%s}^{%s}", ColorChannelNameTexS[j], StateNameTex[i]);
						}
						else{
							sprintf(hist_var_name_f1,"f_{%s}^{%s}", ColorChannelNameTexP[j], StateNameTex[i]);
						}
						if(j==0) sprintf(hist_var_name_f1,"R^{%s}", StateNameTex[i]);

						char h_1Df_name[200];

						sprintf(h_1Df_name,"h_1Df_state%d_f%d",i,j);
						sprintf(projectchar,"%s>>%s", branch_name_f1, h_1Df_name);
						sprintf(selectchar,"acceptedSampling==1 && BurnInInt==0");
						h_1Df[j] 	= new TH1D( h_1Df_name, h_1Df_name, nBins1D, f1Min, f1Max );
						outputTreeAllSamplings->Draw(projectchar, selectchar);

						sprintf(h_1Df_name,"h_1Df2_state%d_f%d",i,j);
						sprintf(projectchar,"%s>>%s", branch_name_f1, h_1Df_name);
						sprintf(selectchar,"acceptedSampling==1 && BurnInInt==0");
						h_1Df2[j] 	= new TH1D( h_1Df_name, h_1Df_name, nBins1D, f1Min, f1Max );
						outputTreeAllSamplings2->Draw(projectchar, selectchar);

						sprintf(h_1Df_name,"h_1Df3_state%d_f%d",i,j);
						sprintf(projectchar,"%s>>%s", branch_name_f1, h_1Df_name);
						sprintf(selectchar,"acceptedSampling==1 && BurnInInt==0");
						h_1Df3[j] 	= new TH1D( h_1Df_name, h_1Df_name, nBins1D, f1Min, f1Max );
						outputTreeAllSamplings3->Draw(projectchar, selectchar);

						h_1Df[j]->Scale(1./h_1Df[j]->GetSumOfWeights());
						h_1Df2[j]->Scale(1./h_1Df2[j]->GetSumOfWeights());
						h_1Df3[j]->Scale(1./h_1Df3[j]->GetSumOfWeights());

						  TCanvas *c1 = new TCanvas("c1", "c1", 10, 28, 650,571);
						  c1->Range(-237.541,-66.47556,187.377,434.8609);
						  c1->SetFillColor(0);
						  c1->SetBorderMode(0);
						  c1->SetBorderSize(0);
						  c1->SetLeftMargin(0.215);
						  c1->SetRightMargin(0.03);
						  c1->SetTopMargin(0.01841621);
						  c1->SetBottomMargin(0.16);
						  c1->SetFrameBorderMode(0);


						  TH1D* h_1Df_axis;
						  h_1Df_axis 	= new TH1D( Form("%s_axis",h_1Df_name), Form("%s_axis",h_1Df_name), nBins1D, f1Min, f1Max );

						 double maximum=h_1Df[j]->GetMaximum();
						 if(h_1Df[j]->GetMaximum() < h_1Df2[j]->GetMaximum()) maximum=h_1Df2[j]->GetMaximum();
						 if(h_1Df2[j]->GetMaximum() < h_1Df3[j]->GetMaximum()) maximum=h_1Df3[j]->GetMaximum();
						 maximum*=1.3;

						 h_1Df_axis->SetMaximum(maximum);
						 h_1Df_axis->SetMinimum(1e-3);

						  char DrawContourStyle[200];
						  sprintf(DrawContourStyle,"c,same");
						  int LineWidth=2;
						  int LineStyle=1;

						  h_1Df_axis->GetXaxis()->SetTitle(hist_var_name_f1);
						  h_1Df_axis->GetXaxis()->SetLabelOffset(0.028);
						  h_1Df_axis->GetXaxis()->SetTitleSize(0.05);
						  h_1Df_axis->GetXaxis()->SetTickLength(-0.03);
						  h_1Df_axis->GetXaxis()->SetTitleOffset(1.4);
						  h_1Df_axis->GetYaxis()->SetTitle("posterior density");
						  h_1Df_axis->GetYaxis()->SetLabelOffset(0.032);
						  h_1Df_axis->GetYaxis()->SetTitleSize(0.05);
						  h_1Df_axis->GetYaxis()->SetTickLength(-0.03);
						  h_1Df_axis->GetYaxis()->SetTitleOffset(1.95);
						  h_1Df_axis->SetTitle(0);
						  h_1Df_axis->SetStats(0);
						  h_1Df_axis->Draw("");

						  h_1Df[j]->SetLineColor( kGreen+2 );
						  h_1Df[j]->SetLineWidth( LineWidth );
						  h_1Df[j]->SetLineStyle( LineStyle  );
						  h_1Df[j]->SetTitle(0);
						  h_1Df[j]->SetStats(0);
						  h_1Df[j]->Draw(DrawContourStyle);

						  h_1Df2[j]->SetLineColor( kBlue );
						  h_1Df2[j]->SetLineWidth( LineWidth );
						  h_1Df2[j]->SetLineStyle( LineStyle  );
						  h_1Df2[j]->SetTitle(0);
						  h_1Df2[j]->SetStats(0);
						  h_1Df2[j]->Draw(DrawContourStyle);

						  h_1Df3[j]->SetLineColor( kRed );
						  h_1Df3[j]->SetLineWidth( LineWidth );
						  h_1Df3[j]->SetLineStyle( LineStyle  );
						  h_1Df3[j]->SetTitle(0);
						  h_1Df3[j]->SetStats(0);
						  h_1Df3[j]->Draw(DrawContourStyle);


							TLegend* legend;
							legend = new TLegend(0.25,0.75,0.475,0.95);
							if(j==2) legend = new TLegend(0.625,0.75,0.955,0.95);
							legend->SetFillColor(0);
							legend->SetTextSize(0.03);
							legend->SetBorderSize(0);
							legend->AddEntry(h_1Df3[j],Form("Cross section constraint"),"l");
							legend->AddEntry(h_1Df2[j],Form("Polarization constraint"),"l");
							legend->AddEntry(h_1Df[j],Form("Both constraints"),"l");

							legend->Draw("same");

						  char filename[200];
						  sprintf(filename,"Figures/PlotPPD/%s/1D_PPD_comp_state%d_f%d.pdf",JobID, i,j);
						  c1->SaveAs(filename);
						  c1->Close();

						  delete c1;

			}
		}



return 0;

}



inline double contourHeight2D ( TH2D *h, double confidenceLevel ) {
  int Nx = h->GetXaxis()->GetNbins();
  int Ny = h->GetYaxis()->GetNbins();

  double totSum = h->GetSum();
  double targetSum = confidenceLevel * totSum;
  double maxHeight = h->GetMaximum();
  double step = 0.001*maxHeight;

  double tempHeight = 0.;
  double tempSum = totSum;

  while ( tempSum > targetSum && tempHeight < maxHeight ) {
    tempHeight += step;
    tempSum = 0.;
    for ( int ix = 0; ix < Nx; ix++ ) {
      for ( int iy = 0; iy < Ny; iy++ ) {
      	double binContent = h->GetBinContent(ix,iy);
        if ( binContent > tempHeight ) tempSum += binContent;
      }
    }
  }
  return tempHeight;
}

// function to set the 99% and 68% C.L. contours of a 2D histogram
inline void setContourHistogram ( TH2D *h ) {
	  double cont0 = contourHeight2D( h, 0.683 );
	  //double cont0 = contourHeight2D( h, 0.953 );
  h->SetContour(1);
  h->SetContourLevel(0,cont0);

}


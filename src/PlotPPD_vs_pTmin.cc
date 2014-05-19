/*
 * PlotPPD_vs_pTmin.cc
 *
 *  Created on: Nov 16, 2013
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
#include "TFrame.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TSystem.h"
#include "Riostream.h"
#include "TString.h"
#include "TROOT.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLine.h"
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
void FindMPV(TH1* PosteriorDist , double& MPV , double& MPVerrorLow, double& MPVerrorHigh, int MPValgo, int nSigma);
vector<double> addPolarizations(vector<vector<double> > lamMatrix, vector<double> lamVecContributions);
void PlotPosterior(char xAxisTitle[200], char filename[200], TH1* histo, double MPV, double MPVerrLow, double MPVerrHigh, int nSigma);
double func_pT_ratio_NLO(double* x, double* par);


int main(int argc, char** argv) {


	cout<<"start PlotPPD"<<endl;

  	Char_t *JobID = "Default";
	Char_t *storagedir = "Default"; //Storage Directory
	int useOnlyState=999;
  	int 	MPValgo=-1;
  	double 	nSigma=-1;
	bool 	PlotVSpToverM=false;

  	for( int i=0;i < argc; ++i ) {
  		cout<<"i arcgc "<<i<<endl;
		if(std::string(argv[i]).find("JobID") != std::string::npos) {char* JobIDchar = argv[i]; char* JobIDchar2 = strtok (JobIDchar, "="); JobID = JobIDchar2; cout<<"JobID = "<<JobID<<endl;}
		if(std::string(argv[i]).find("storagedir") != std::string::npos) {char* storagedirchar = argv[i]; char* storagedirchar2 = strtok (storagedirchar, "="); storagedir = storagedirchar2; cout<<"storagedir = "<<storagedir<<endl;}
	    if(std::string(argv[i]).find("nSigma") != std::string::npos) { char* nSigmachar = argv[i]; char* nSigmachar2 = strtok (nSigmachar, "n"); nSigma = atof(nSigmachar2); cout<<"nSigma = "<<nSigma<<endl; }
	    if(std::string(argv[i]).find("MPValgo") != std::string::npos) { char* MPValgochar = argv[i]; char* MPValgochar2 = strtok (MPValgochar, "M"); MPValgo = atoi(MPValgochar2); cout<<"MPValgo = "<<MPValgo<<endl; }
	    if(std::string(argv[i]).find("useOnlyState") != std::string::npos) {
	    	char* useOnlyStatechar = argv[i];
	    	char* useOnlyStatechar2 = strtok (useOnlyStatechar, "u");
	    	useOnlyState = atoi(useOnlyStatechar2);
	    	cout<<"useOnlyState = "<<useOnlyState<<endl;
	    }
	    if(std::string(argv[i]).find("PlotVSpToverM=true") != std::string::npos) {
	    	PlotVSpToverM=true;
	    	cout<<"PlotVSpToverM=true"<<endl;
	    }
  	}



  	nSigma=2;

  	gROOT->SetBatch();

	delete gRandom;
	gRandom = new TRandom3(NRQCDvars::randomSeed);  // better random generator


	// TODO: add loop for nSystematics (now nSystematics will be 0)

	char outname[200];
	char inname[200];
	char plotdirname[200];
	char jobdirname[200];
	sprintf(plotdirname,"Figures/PlotPPD_vs_pTmin/%s",JobID);
	gSystem->mkdir(plotdirname);
	sprintf(jobdirname,"%s/JobID/%s",storagedir, JobID);
	gSystem->mkdir(jobdirname);



	//allCharm
	//const int nSeq=15;
	//double pTmin[nSeq]={4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18};
	//Psi2S
	//const int nSeq=13;
	//double pTmin[nSeq]={4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 18};
	//Psi2S prel CMS data
	//const int nSeq=21;
	//double pTmin[nSeq]={4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 25, 27, 30};
	//Ups3S
	const int nSeq=19;
	double pTmin[nSeq]={10, 11, 12, 13, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 43};

	int ScanState=10;

	char JobIDs[200];
	TTree*  outputTreeAllSamplings[nSeq];
	double  chi2Prob[nSeq];
	double  chi2ndf[nSeq];
	double errpTmin[nSeq];
	double zero[nSeq];

	for(int iSeq=0;iSeq<nSeq;iSeq++){

		if(ScanState==3) sprintf(JobIDs,"ScanPsi2S_pTmin%1.0f_December5",pTmin[iSeq]);
		if(ScanState==10) sprintf(JobIDs,"ScanUps3S_pTmin%1.0f_December5",pTmin[iSeq]);

		//sprintf(JobIDs,"ScanPsi1S_pTmin%1.0f_December7",pTmin[iSeq]);
		//if(iSeq<5) sprintf(JobIDs,"ScanPsi1S_pTmin%1.0f_December7_v2",pTmin[iSeq]);

		sprintf(inname,"%s/JobID/%s/results.root",storagedir, JobIDs);
		TFile *ResultsFile = new TFile(inname, "READ");
	  	outputTreeAllSamplings[iSeq] = (TTree*) ResultsFile->Get("AllSamplings"); // tree filled in all samplings after burnin

		sprintf(inname,"%s/JobID/%s/chi2ndf.txt",storagedir, JobIDs);
		cout<<"read chi2 from "<<inname<<endl;

	  	double chi2Min;
	  	double reduced_chi2Min;
	  	int nDataPoints;
	  	int nOps;
	  	int ndf;

		ifstream in;
	    in.open(inname);//, std::ofstream::app);

		in >> chi2Min;
		in >> nDataPoints;
		in >> nOps;
		in >> ndf;
		in >> reduced_chi2Min;

	    in.close();

	    cout<<"chi2Min:"<<chi2Min<<endl;
	    cout<<"nDataPoints:"<<nDataPoints<<endl;
	    cout<<"nOps:"<<nOps<<endl;
	    cout<<"ndf:"<<ndf<<endl;
	    cout<<"reduced_chi2Min:"<<reduced_chi2Min<<endl;

	    chi2Prob[iSeq]=TMath::Prob(chi2Min, ndf);
	    chi2ndf[iSeq]=reduced_chi2Min;

	    errpTmin[iSeq]=0.25;
	    zero[iSeq]=0.;
		if(PlotVSpToverM) pTmin[iSeq]/=NRQCDvars::mass[ScanState];


	}

	char PlotID[200];
	char hist_var_name[200];
	char branch_name[200];
	char branch_name_plot[200];
	char branch_name_plotD[200];
	char hist_var_name_plot[200];
	int i;
	int j;
	double minY;
	double maxY;
	bool plotLogy;
	bool plotDashedLine;
	double plotDashedLineAt;
	bool plotRatio;
	int LegendPos;
	bool plotChiStyle;
	int ColorSetting;
	int k;

	int nPlots=7;


	for(int iPlot=0;iPlot<nPlots;iPlot++){

		gStyle->SetPalette(1,0);
		gStyle->SetPadBottomMargin(0.12);
		gStyle->SetPadLeftMargin(0.13);
		gStyle->SetPadRightMargin(0.15);
		gStyle->SetPadTopMargin(0.1);

		gStyle->SetTickLength(-0.02, "xyz");
		gStyle->SetLabelOffset(0.02, "x");
		gStyle->SetLabelOffset(0.02, "y");
		gStyle->SetTitleOffset(1.3, "x");
		gStyle->SetTitleOffset(1.4, "y");
		gStyle->SetTitleFillColor(kWhite);



		if(iPlot==0){
			//continue;
			i=ScanState;
			j=1;
			minY=-0.25;
			maxY=1.25;
			if(ScanState==1 || ScanState==2){
				minY=-0.5;
				maxY=1.5;
			}
			if(ScanState==0){
				minY=-0.5;
				maxY=1.5;
			}
			plotLogy=false;
			plotChiStyle=false;
			plotRatio=false;
			sprintf(PlotID,"FRAC_S0");
			sprintf(hist_var_name,"Relative fraction of %s contribution", ColorChannelNameTexS[j]);
			sprintf(branch_name_plot,"state%d_f%d",i,j);
			LegendPos=3;
			ColorSetting=0;
			plotDashedLineAt=0.;
			plotDashedLine=false;
		}
		if(iPlot==1){
			//continue;
			i=ScanState;
			j=2;
			minY=-0.25;
			maxY=1.25;
			if(ScanState==1 || ScanState==2){
				minY=-0.5;
				maxY=1.5;
			}
			if(ScanState==0){
				minY=-0.5;
				maxY=1.5;
			}
			plotLogy=false;
			plotChiStyle=false;
			plotRatio=false;
			sprintf(PlotID,"FRAC_S1");
			sprintf(hist_var_name,"Relative fraction of %s contribution", ColorChannelNameTexS[j]);
			sprintf(branch_name_plot,"state%d_f%d",i,j);
			LegendPos=2;
			ColorSetting=1;
			plotDashedLineAt=0.;
			plotDashedLine=false;
		}
		if(iPlot==2){
			//continue;
			i=ScanState;
			j=1;
			if(ScanState==3 || ScanState==10){
				minY=-0.01;
				maxY=0.05;
			}
			if(ScanState==0){
				minY=-0.01;
				maxY=0.11;
			}
			if(ScanState==1 || ScanState==2){
				minY=-0.02;
				maxY=0.25;
			}
			plotLogy=false;
			plotChiStyle=false;
			plotRatio=false;
			sprintf(PlotID,"LDME_S0");
			sprintf(hist_var_name,"%s long distance matrix element [GeV^{3}]", ColorChannelNameTexS[j]);
			sprintf(branch_name_plot,"state%d_Op%d",i,j);
			LegendPos=1;
			ColorSetting=0;
			plotDashedLineAt=0.;
			plotDashedLine=false;
		}
		if(iPlot==3){
			//continue;
			i=ScanState;
			j=2;
			minY=-0.6;
			maxY=3.;
			if(ScanState==1 || ScanState==2){
				minY=-1.5;
				maxY=7.;
			}
			plotLogy=false;
			plotChiStyle=false;
			plotRatio=false;
			sprintf(PlotID,"LDME_S1");
			sprintf(hist_var_name,"%s long distance matrix element [GeV^{3}]", ColorChannelNameTexS[j]);
			sprintf(branch_name_plot,"state%d_Op%d",i,j);
			LegendPos=2;
			ColorSetting=1;
			plotDashedLineAt=0.;
			plotDashedLine=false;
		}
		if(iPlot==4){
			//continue;
			i=ScanState;
			j=2;
			minY=1e-26;
			maxY=5.;
			plotLogy=true;
			plotChiStyle=true;
			plotRatio=false;
			sprintf(PlotID,"CHI2PROB");
			sprintf(hist_var_name,"#chi^{2} probability", ColorChannelNameTexS[j]);
			sprintf(branch_name_plot,"DUMMY");
			LegendPos=2;
			ColorSetting=1;
			plotDashedLineAt=1.;
			plotDashedLine=true;
		}
		if(iPlot==5){
			//continue;
			i=ScanState;
			j=2;
			minY=-0.0009999;
			maxY=8.;
			if(ScanState==1 || ScanState==2 || ScanState==0){
				minY=-0.0009999;
				maxY=14.;
			}
			plotLogy=false;
			plotChiStyle=true;
			plotRatio=false;
			sprintf(PlotID,"CHI2NDF");
			sprintf(hist_var_name,"#chi^{2} / ndf", ColorChannelNameTexS[j]);
			sprintf(branch_name_plot,"DUMMY");
			LegendPos=2;
			ColorSetting=1;
			plotDashedLineAt=1.;
			plotDashedLine=true;
		}
		if(iPlot==6){
			continue;
			i=ScanState;
			j=2;
			k=1;
			minY=5e-4;
			maxY=1.;
			plotLogy=true;
			plotChiStyle=false;
			plotRatio=true;
			sprintf(PlotID,"LDME_S1_ov_LDME_S0");
			sprintf(hist_var_name,"%s / %s long distance matrix elements", ColorChannelNameTexS[j], ColorChannelNameTexS[k]);
			sprintf(branch_name_plot,"state%d_Op%d",i,j);
			sprintf(branch_name_plotD,"state%d_Op%d",i,k);
			LegendPos=2;
			ColorSetting=2;
			plotDashedLineAt=1.;
			plotDashedLine=false;
		}

		cout<<"branch name to plot: "<<branch_name_plot<<endl;


		double buff_MPV[nSeq];
		double buff_errlow_1sig[nSeq];
		double buff_errhigh_1sig[nSeq];
		double buff_errlow_2sig[nSeq];
		double buff_errhigh_2sig[nSeq];
		double buff_errlow_3sig[nSeq];
		double buff_errhigh_3sig[nSeq];


		for(int iSeq=0;iSeq<nSeq;iSeq++){
			if(plotChiStyle) continue;
			cout<<"iSeq: "<<iSeq<<endl;



			double plotObservable;
			double plotObservableD;

			int BurnInInt;
			sprintf (branch_name,"BurnInInt");
			outputTreeAllSamplings[iSeq]->SetBranchAddress( branch_name,  &BurnInInt );
			int acceptedSampling;
			sprintf (branch_name,"acceptedSampling");
			outputTreeAllSamplings[iSeq]->SetBranchAddress( branch_name,  &acceptedSampling );


			outputTreeAllSamplings[iSeq]->SetBranchAddress( branch_name_plot,  &plotObservable );
			if(plotRatio) 	outputTreeAllSamplings[iSeq]->SetBranchAddress( branch_name_plotD,  &plotObservableD );

			sprintf(outname,"%s/plotTree.root", jobdirname);
			TFile *DummyFile = new TFile(outname, "RECREATE");

			TTree*  plotTree;
			double plotObservableOut;
			plotTree = new TTree("plotTree", "plotTree");
			plotTree->Branch("plotObservableOut",     &plotObservableOut,     "plotObservableOut/D");

			int n_events = int( outputTreeAllSamplings[iSeq]->GetEntries() );
			cout<<n_events<<" events: Looping through nTuple of PPD"<<endl;

			// loop over  events in the model ntuple
			for ( int i_event = 1; i_event <= n_events; i_event++ ) {


				outputTreeAllSamplings[iSeq]->GetEvent( i_event-1 );

				if(acceptedSampling<0 || BurnInInt==1) continue;

				plotObservableOut=plotObservable;
				if(plotRatio) plotObservableOut=plotObservable/plotObservableD;

				bool fillTree=true;

				double maxOut=1e3;
				if(iSeq==2) maxOut=100;
				if(plotRatio){
					if(TMath::Abs(plotObservableOut)>maxOut) fillTree=false;
				}

				//plotObservableOut=TMath::Abs(plotObservableOut);

				if(iPlot==3) plotObservableOut*=1000.;

				if(fillTree) plotTree->Fill();


			}

			//plotTree->Print();


			int nBins_h=100;
			char hist_name[200];
			char projectchar[200];
			char selectchar[200];
			TH1D* h_plotObservable;
			double h_plotObservable_min;
			double h_plotObservable_max;

			double expandMinMaxBy=0.01;

			h_plotObservable_min=plotTree->GetMinimum("plotObservableOut")-expandMinMaxBy*TMath::Abs(plotTree->GetMinimum("plotObservableOut"));
			h_plotObservable_max=plotTree->GetMaximum("plotObservableOut")+expandMinMaxBy*TMath::Abs(plotTree->GetMaximum("plotObservableOut"));

			sprintf(hist_name,"h_plotObservable");
			h_plotObservable = new TH1D( hist_name, hist_name, nBins_h,h_plotObservable_min, h_plotObservable_max );
			sprintf(projectchar,"plotObservableOut>>%s",hist_name);
			plotTree->Draw(projectchar);

			FindMPV(h_plotObservable, buff_MPV[iSeq], buff_errlow_1sig[iSeq], buff_errhigh_1sig[iSeq], MPValgo, 1);
			FindMPV(h_plotObservable, buff_MPV[iSeq], buff_errlow_2sig[iSeq], buff_errhigh_2sig[iSeq], MPValgo, 2);
			FindMPV(h_plotObservable, buff_MPV[iSeq], buff_errlow_3sig[iSeq], buff_errhigh_3sig[iSeq], MPValgo, 3);

			//if(iPlot==2){
			//	h_plotObservable->Print("all");
			//}

			char filename[200];
			sprintf(filename,"Figures/PlotPPD_vs_pTmin/%s",JobID);
			gSystem->mkdir(filename);
			sprintf(filename,"%s/PPD_%s_pTmin%d.pdf",filename, PlotID, iSeq);

			cout<<filename<<endl;

			h_plotObservable -> SetTitle("");
			h_plotObservable -> SetXTitle(hist_var_name);

			PlotPosterior(hist_var_name, filename, h_plotObservable, buff_MPV[iSeq], buff_errlow_1sig[iSeq], buff_errhigh_1sig[iSeq], 1);



			delete plotTree;
			delete h_plotObservable;

			DummyFile->Close();
			delete DummyFile;


		}


		TGraphAsymmErrors *g_plotObservable_vs_pTMin = new TGraphAsymmErrors(nSeq,pTmin,buff_MPV, zero, zero, zero, zero);
		TGraphAsymmErrors *g_plotObservable_vs_pTMin_1sig = new TGraphAsymmErrors(nSeq,pTmin,buff_MPV, errpTmin, errpTmin, buff_errlow_1sig, buff_errhigh_1sig);
		TGraphAsymmErrors *g_plotObservable_vs_pTMin_2sig = new TGraphAsymmErrors(nSeq,pTmin,buff_MPV, errpTmin, errpTmin, buff_errlow_2sig, buff_errhigh_2sig);
		TGraphAsymmErrors *g_plotObservable_vs_pTMin_3sig = new TGraphAsymmErrors(nSeq,pTmin,buff_MPV, errpTmin, errpTmin, buff_errlow_3sig, buff_errhigh_3sig);

		if(plotChiStyle&&iPlot==4) g_plotObservable_vs_pTMin = new TGraphAsymmErrors(nSeq,pTmin,chi2Prob, zero, zero, zero, zero);
		if(plotChiStyle&&iPlot==5) g_plotObservable_vs_pTMin = new TGraphAsymmErrors(nSeq,pTmin,chi2ndf, zero, zero, zero, zero);

		char graph_var_name[200];
		sprintf(graph_var_name,"#it{p}^{min}_{T} [GeV]");
		if(PlotVSpToverM) sprintf(graph_var_name,"#it{p}^{min}_{T} / M");
		g_plotObservable_vs_pTMin -> GetXaxis() -> SetTitle(graph_var_name);
		g_plotObservable_vs_pTMin -> GetYaxis() -> SetTitle(hist_var_name);

		g_plotObservable_vs_pTMin -> GetXaxis() -> SetTitleOffset(1.5);
		g_plotObservable_vs_pTMin -> GetYaxis() -> SetTitleOffset(1.75);

		//if(iPlot==4) g_plotObservable_vs_pTMin -> GetYaxis() -> SetNdivisions(8);

		//g_plotObservable_vs_pTMin -> GetYaxis() -> SetLimits(minY, maxY);
		g_plotObservable_vs_pTMin -> SetMinimum(minY);
		g_plotObservable_vs_pTMin -> SetMaximum(maxY);
		g_plotObservable_vs_pTMin -> SetTitle("");

		g_plotObservable_vs_pTMin_1sig->Print();

		////////////////////////////////////////
		/// Plot vs pTmin
		////////////////////////////////////////






			gROOT->Reset();
			gStyle->SetOptStat(11);
			gStyle->SetOptFit(101);
			gStyle->SetFillColor(kWhite);

			gStyle->SetPadBottomMargin(0.12);
			gStyle->SetPadLeftMargin(0.12);
			gStyle->SetPadRightMargin(0.02);
			gStyle->SetPadTopMargin(0.02);






			TCanvas *plotCanvaspTmin = new TCanvas("plotCanvaspTmin","plotCanvaspTmin",1000,800);
			plotCanvaspTmin->SetFillColor(kWhite);
			plotCanvaspTmin->SetFrameBorderMode(0);
			plotCanvaspTmin->SetRightMargin(0.02);
			plotCanvaspTmin->SetLeftMargin(0.125);
			plotCanvaspTmin->SetTopMargin(0.05);


			int LatexPos=4;
			TLegend* legend;
			if(LegendPos==1) legend = new TLegend(0.16,0.7510989,0.31,0.925008242);
			if(LegendPos==2) legend = new TLegend(0.8,0.7510989,0.95,0.925008242);
			if(LegendPos==3) legend = new TLegend(0.8,0.1310989,0.95,0.30008242);
			if(LegendPos==4) legend = new TLegend(0.16,0.1310989,0.31,0.30008242);
			legend->SetFillColor(0);
			legend->SetTextSize(0.025);
			legend->SetBorderSize(0);


			char DrawChar[200];
			char legendentry[200];


			const int nColorSettings=3;

			int FillStyle_sig[nColorSettings]={1001,1001, 1001};
			int LineStyle_Central[nColorSettings]={1,1, 1};
			int LineColor_Central[nColorSettings]={416+2,632+2, 600};
			int LineColor_1Sig[nColorSettings]={416-3,632+0, 600-7};
			int LineColor_2Sig[nColorSettings]={416-7,632-7, 600-9};
			int LineColor_3Sig[nColorSettings]={416-10,632-10, 600-10};
			double linewidthCC=2.;


			g_plotObservable_vs_pTMin->SetMarkerStyle(20);
			g_plotObservable_vs_pTMin->SetMarkerColor(1);
			g_plotObservable_vs_pTMin->SetMarkerSize(1.25);

			g_plotObservable_vs_pTMin->SetLineStyle(LineStyle_Central[ColorSetting]);
			g_plotObservable_vs_pTMin->SetLineColor(LineColor_Central[ColorSetting]);
			g_plotObservable_vs_pTMin->SetLineWidth(linewidthCC);

			g_plotObservable_vs_pTMin_1sig->SetFillStyle(FillStyle_sig[ColorSetting]);
			g_plotObservable_vs_pTMin_1sig->SetFillColor(LineColor_1Sig[ColorSetting]);
			g_plotObservable_vs_pTMin_1sig->SetLineColor(LineColor_1Sig[ColorSetting]);
			g_plotObservable_vs_pTMin_2sig->SetFillStyle(FillStyle_sig[ColorSetting]);
			g_plotObservable_vs_pTMin_2sig->SetFillColor(LineColor_2Sig[ColorSetting]);
			g_plotObservable_vs_pTMin_2sig->SetLineColor(LineColor_2Sig[ColorSetting]);
			g_plotObservable_vs_pTMin_3sig->SetFillStyle(FillStyle_sig[ColorSetting]);
			g_plotObservable_vs_pTMin_3sig->SetFillColor(LineColor_3Sig[ColorSetting]);
			g_plotObservable_vs_pTMin_3sig->SetLineColor(LineColor_3Sig[ColorSetting]);


			double expandBy=0.1;
			double pTminAxis=pTmin[0]-(pTmin[nSeq-1]-pTmin[0])*expandBy;
			double pTmaxAxis=pTmin[nSeq-1]+(pTmin[nSeq-1]-pTmin[0])*expandBy;

			if(PlotVSpToverM){
				pTminAxis=0.5001;
				pTmaxAxis=5.4999;
			}

			g_plotObservable_vs_pTMin->GetXaxis()->SetLimits(pTminAxis, pTmaxAxis);

			//g_plotObservable_vs_pTMin_3sig->Draw("a");
			if(!plotChiStyle){
				g_plotObservable_vs_pTMin->Draw("al");
				g_plotObservable_vs_pTMin_3sig->Draw("3same");
				g_plotObservable_vs_pTMin_2sig->Draw("3same");
				g_plotObservable_vs_pTMin_1sig->Draw("3same");
				g_plotObservable_vs_pTMin->Draw("lsame");
			}

			if(plotChiStyle){
				g_plotObservable_vs_pTMin->Draw("alp");
			}

			if(plotDashedLine){
				TLine* Dashedline = new TLine( pTminAxis, plotDashedLineAt, pTmaxAxis, plotDashedLineAt);
				Dashedline->SetLineWidth( 1 );
				Dashedline->SetLineStyle( 2 );
				Dashedline->SetLineColor( kBlack );
				Dashedline->Draw( "same" );
			}

			sprintf(legendentry,"68.3%% CL");
			legend->AddEntry(g_plotObservable_vs_pTMin_1sig,legendentry,"f");
			sprintf(legendentry,"95.5%% CL");
			legend->AddEntry(g_plotObservable_vs_pTMin_2sig,legendentry,"f");
			sprintf(legendentry,"99.7%% CL");
			legend->AddEntry(g_plotObservable_vs_pTMin_3sig,legendentry,"f");

			int iStatebuff;
			sprintf(DrawChar,"psame");
			char rapchar[200];



			if(!plotChiStyle) legend->Draw("same");

			//draw latex
			double left=0.175, top=0.92, textSize=0.045;
			TLatex *latex=new TLatex();
			latex->SetTextFont(42);
			latex->SetNDC(kTRUE);
			latex->SetTextSize(textSize);
			double stepLatex=textSize*1.3;
			latex->SetTextFont(62);

			if(iPlot==3){
				left=0.125;
				top=0.96;
				textSize=0.035;
				latex->SetTextSize(textSize);
				latex->DrawLatex(left,top, Form("x 10^{-3}"));
			}


			char saveDirName[500];
			char saveName[500];
			sprintf(saveDirName,"Figures/PlotPPD_vs_pTmin/%s",JobID);
			sprintf(saveName,"%s/Plot_vs_pTmin_%s.pdf",saveDirName, PlotID);
			cout<<saveName<<endl;
			if(plotLogy) plotCanvaspTmin->SetLogy(true);
			else plotCanvaspTmin->SetLogy(false);
			plotCanvaspTmin->SaveAs(saveName);
			sprintf(saveName,"%s/PNG_Plot_vs_pTmin_%s.png",saveDirName, PlotID);
			cout<<saveName<<endl;
			plotCanvaspTmin->SaveAs(saveName);


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

vector<double> addPolarizations(vector<vector<double> > lamMatrix, vector<double> lamVecContributions){

	vector<double> lamSumVec (3,0); //return this

	int nContributions=lamVecContributions.size();

	double ContributionIntegral=0;
	for(int i=0; i<nContributions; i++){
		ContributionIntegral+=lamVecContributions[i];
	}
	if(ContributionIntegral==0) {cerr << "ContributionIntegral is 0 -> exit!"<<endl; exit(2);}
	vector<double> lamVecContributionFraction (nContributions,0);
	for(int i=0; i<nContributions; i++){
		lamVecContributionFraction[i]=lamVecContributions[i]/ContributionIntegral;

	}

	vector<double> LamthVec (nContributions);
	vector<double> LamphVec (nContributions);
	vector<double> LamtpVec (nContributions);

	for(int i=0; i<nContributions; i++){
		LamthVec[i]=lamMatrix[i][0];
		LamphVec[i]=lamMatrix[i][1];
		LamtpVec[i]=lamMatrix[i][2];
	}

	double Lamth_numerator=0;
	double Lamth_denominator=0;
	for(int i=0; i<nContributions; i++){
		Lamth_numerator+=lamVecContributionFraction[i]*LamthVec[i]/(3+LamthVec[i]);
		Lamth_denominator+=lamVecContributionFraction[i]/(3+LamthVec[i]);
	}

	double Lamph_numerator=0;
	for(int i=0; i<nContributions; i++){
		Lamph_numerator+=lamVecContributionFraction[i]*LamphVec[i]/(3+LamthVec[i]);
	}

	double Lamtp_numerator=0;
	for(int i=0; i<nContributions; i++){
		Lamtp_numerator+=lamVecContributionFraction[i]*LamtpVec[i]/(3+LamthVec[i]);
	}

	lamSumVec[0]=Lamth_numerator/Lamth_denominator;
	lamSumVec[1]=Lamph_numerator/Lamth_denominator;
	lamSumVec[2]=Lamtp_numerator/Lamth_denominator;


	return lamSumVec;
}


void FindMPV(TH1* PosteriorDist , double& MPV , double& MPVerrorLow, double& MPVerrorHigh, int MPValgo, int nSigma){

	//PosteriorDist->Print();

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

	if(MPValgo==2||MPValgo==3||MPValgo==4){

		int nBins = PosteriorDist->GetNbinsX();
		int maxbin_PosteriorDist = PosteriorDist->GetMaximumBin();
		double PosteriorDist_initial = PosteriorDist->GetBinCenter(maxbin_PosteriorDist);
		double err_PosteriorDist_initial=PosteriorDist->GetRMS();
		double PosteriorDist_par [3];

		TF1 *gauss;

		if(MPValgo==4) MPV=PosteriorDist_initial;
				else{
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
						double chi2max=4.;
						if(gauss->GetChisquare()/ndof<chi2max) {cout<<"chi2 < "<<chi2max<<" -> good enough"<<endl; break;}
						if(iFits==nMaxFits-1) illPPD=true;
					}
					MPV=PosteriorDist_par[1];
				}

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



void PlotPosterior(char xAxisTitle[200], char filename[200], TH1* histo, double MPV, double MPVerrLow, double MPVerrHigh, int nSigma){


	TLegend* plotLegend=new TLegend(0.745,0.8,0.92,0.9);
	plotLegend->SetFillColor(kWhite);
	plotLegend->SetTextFont(72);
	plotLegend->SetTextSize(0.0175);
	plotLegend->SetBorderSize(1);

	TCanvas *plotCanvas = new TCanvas("plotCanvas","plotCanvas",1000,800);

	plotCanvas->SetFillColor(kWhite);
	//	plotCanvas->SetGrid();
	plotCanvas->GetFrame()->SetFillColor(kWhite);
	plotCanvas->GetFrame()->SetBorderSize(0);
	plotCanvas->SetRightMargin(0.08) ;

	histo->SetStats(kFALSE);
	histo->SetLineColor(kBlack);
	histo->SetYTitle("Posterior density");
	histo->SetXTitle(xAxisTitle);
	histo->GetYaxis()->SetTitleOffset(2.);

	histo->Scale(1./histo->Integral(0,histo->GetNbinsX()+1));
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
	plotLegend->AddEntry(MPVerrLowLine,Form("%d#sigma high/low",nSigma),"l");

	plotLegend->Draw();
	plotCanvas->SaveAs(filename);
	plotCanvas->Close();

	delete plotCanvas;

}

double func_pT_ratio_NLO(double* x, double* par) {


	double funcval;

	double norm=par[0];
	double beta=par[1];
	double gamma1=par[2];
	double gamma2=par[3];

	//norm * [(gamma1 + pT^2)/(gamma2 + pT^2)]^(-beta)

	double pTFit;


		pTFit=x[0];

	funcval = norm*pow( (gamma1 + pTFit*pTFit)/(gamma2 + pTFit*pTFit), -beta );


	return funcval;

}


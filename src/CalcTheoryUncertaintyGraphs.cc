/*
 * CalcTheoryUncertaintyGraphs.cc
 *
 *  Created on: Mar 27, 2014
 *      Author: valentinknuenz
 */



#ifndef __CINT__
#endif

#include <stdio.h>
#include <iostream.h>
#include <fstream.h>
#include "TROOT.h"
#include "iomanip.h"
#include "TStyle.h"
#include "TF1.h"
#include "TROOT.h"
#include "TGraphAsymmErrors.h"
#include "TFile.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TMath.h"

using namespace std;


double exponential(double* x, double* par);
Double_t paramMassRapParabola(Double_t *x, Double_t *par);
Double_t paramLog(Double_t *x, Double_t *par);

void CalcTheoryUncertaintyGraphs(int iModelN, int iModelD, int iMeasID=0, int nRapModel=0, int iRapModel=0, int iState=0, int iFrame=0, double pTMin=0, double pTMax=0, double ResMin=0, double ResMax=0, int LegendPos=1, int LatexPos=1){ //nState=1,2,3:Upsi1S,2S,3S
	gROOT->Reset();
	gStyle->SetOptStat(11);
	gStyle->SetOptFit(101);
	gStyle->SetFillColor(kWhite);

	gStyle->SetPadBottomMargin(0.12);
	gStyle->SetPadLeftMargin(0.12);
	gStyle->SetPadRightMargin(0.02);
	gStyle->SetPadTopMargin(0.02);

	//storagedir=/scratch/knuenz/NRQCD/NRQCDglobalfit

	bool scaleToOne=false;//THuncertainty:false

	char ModelIDN[500];
	char ModelIDD[500];

	if(iModelN==0)
		sprintf(ModelIDN,"NLO_BK_Ep");
	if(iModelN==1)
		sprintf(ModelIDN,"NLO_GWWZ");
	if(iModelN==2)
		sprintf(ModelIDN,"LO_BK");

	if(iModelD==0)
		sprintf(ModelIDD,"NLO_BK_Ep");
	if(iModelD==1)
		sprintf(ModelIDD,"NLO_GWWZ");
	if(iModelD==2)
		sprintf(ModelIDD,"LO_BK");


	char ModelIDNModelIDD[500];
	sprintf(ModelIDNModelIDD,"%s/%s", ModelIDN, ModelIDD);
	if(iModelN==0&&iModelD==1) sprintf(ModelIDNModelIDD,"NLO NRQCD BK/GWWZ");
	if(iModelN==2&&iModelD==0) sprintf(ModelIDNModelIDD,"LO/NLO NRQCD");


	bool debug=false;
	enum {charm, bottom};
	enum {Metropolis, MetropolisHastings}; // modified by Joao: choice of kernel
	enum {MH, Minuit}; //choice of Minimizer
	enum {PSI_1S, CHIC1_1P, CHIC2_1P, PSI_2S, UPS_1S, CHIB1_1P, CHIB2_1P, UPS_2S,
		  CHIB1_2P, CHIB2_2P, UPS_3S, CHIB1_3P, CHIB2_3P};
	const int nStates = (CHIB2_3P-PSI_1S)+1;// modified by Joao: always use last state - first state to set this number
	enum {CrossSection, Lamth, Lamph, Lamtp, CrossSectionRatio, FeedDownFraction};
	const int nMeasurementIDs = (FeedDownFraction-CrossSection)+1; // modified by Joao: : always use last measur. - first measur. to set this number
	enum {CMS, LHCb, ATLAS, ALICE, CDF};
	const int nExperiments = (CDF-CMS)+1; // modified by Joao: always use last exp. - first exp. to set this number
	bool isAbsRapExp[nExperiments] = {true, false, true, false, true};
	Char_t *StateName[nStates] = {"PSI_1S", "CHIC1_1P", "CHIC2_1P", "PSI_2S", "UPS_1S", "CHIB1_1P", "CHIB2_1P", "UPS_2S", "CHIB1_2P", "CHIB2_2P", "UPS_3S", "CHIB1_3P", "CHIB2_3P"};
	Char_t *StateNameTex[nStates] = {"#psi(1S)", "#chi_{c1}(1P)", "#chi_{c2}(1P)", "#psi(2S)", "#Upsilon(1S)", "#chi_{b1}(1P)", "#chi_{b2}(1P)", "#Upsilon(2S)", "#chi_{b1}(2P)", "#chi_{b2}(2P)", "#Upsilon(3S)", "#chi_{b1}(3P)", "#chi_{b2}(3P)"};
	Char_t *MeasurementIDName[nMeasurementIDs] = {"CrossSection", "Lamth", "Lamph", "Lamtp", "CrossSectionRatio", "FeedDownFraction"};
	Char_t *MeasurementIDNameTex[nMeasurementIDs] = {"#sigma [nB]", "#lambda_{#vartheta}^{HX}", "#lambda_{#varphi}^{HX}", "#lambda_{#vartheta#varphi}^{HX}", "CrossSectionRatio", "FeedDownFraction"};
	Char_t *ExpName[nExperiments] = {"CMS", "LHCb", "ATLAS", "ALICE", "CDF"};
	Double_t mass[nStates] = {3.096916, 3.51066, 3.55620, 3.686109, 9.46030, 9.89278, 9.91221, 10.02326, 10.25546,
			                   10.26865, 10.3552, 10.530, 10.530};

	enum {quID_S, quID_P1, quID_P2}; //Definition of QuantumID (S, P1, P2 states)
	int StateQuantumID[nStates]={quID_S, quID_P1, quID_P2, quID_S, quID_S, quID_P1, quID_P2, quID_S, quID_P1, quID_P2, quID_S, quID_P1, quID_P2};

	const int nModelSystematicScales = 0;
	const int nDataSystematicScales = nExperiments;

	const int nColorChannels_S=3;//includes CS
	const int nColorChannels_P=3;//includes CS
	Char_t *ColorChannelNameTexS[nColorChannels_S] = {"^{3}S_{1}^{(1)}", "^{1}S_{0}^{(8)}", "^{3}S_{1}^{(8)}"};
	Char_t *ColorChannelNameTexP[nColorChannels_P] = {"^{3}P_{J}^{(1)}", "^{1}S_{0}^{(8)}", "^{3}S_{1}^{(8)}"};


	//const int nColorChannels_S=3;//includes CS
	//const int nColorChannels_P=2;//includes CS
	//Char_t *ColorChannelNameTexS[nColorChannels_S] = {"^{3}S_{1}^{(1)}", "^{1}S_{0}^{(8)}", "^{3}S_{1}^{(8)}"};
	//Char_t *ColorChannelNameTexP[nColorChannels_P] = {"^{3}P_{J}^{(1)}", "^{3}S_{1}^{(8)}"};

	//const int nColorChannels_S=2;//includes CS
	//const int nColorChannels_P=2;//includes CS
	//Char_t *ColorChannelNameTexS[nColorChannels_S] = {"^{3}S_{1}^{(1)}", "^{1}S_{0}^{(8)}"};
	//Char_t *ColorChannelNameTexP[nColorChannels_P] = {"^{3}P_{J}^{(1)}", "^{3}S_{1}^{(8)}"};





	//const int nColorChannels_S=4;//includes CS
	//const int nColorChannels_P=2;//includes CS
	//Char_t *ColorChannelNameTexS[nColorChannels_S] = {"^{1}S_{0}^{(8)}, LO", "^{3}S_{1}^{(8)}, LO", "^{1}S_{0}^{(8)}, NLO", "^{3}S_{1}^{(8)}, NLO"};
	//Char_t *ColorChannelNameTexP[nColorChannels_P] = {"^{3}P_{J}^{(1)}", "^{3}S_{1}^{(8)}"};

	//const int nColorChannels_S=2;//includes CS
	//const int nColorChannels_P=2;//includes CS
	//Char_t *ColorChannelNameTexS[nColorChannels_S] = {"^{3}S_{1}^{(8)} / ^{1}S_{0}^{(8)}, LO", "^{3}S_{1}^{(8)} / ^{1}S_{0}^{(8)}, NLO"};
	//Char_t *ColorChannelNameTexP[nColorChannels_P] = {"^{3}P_{J}^{(1)}", "^{3}S_{1}^{(8)}"};

	const int nColorChannels = std::max(nColorChannels_S, nColorChannels_P);
	//const int nColorChannels4States[nStates]={nColorChannels_S, nColorChannels_P, nColorChannels_P, nColorChannels_S, nColorChannels_S, nColorChannels_P, nColorChannels_P, nColorChannels_S, nColorChannels_P, nColorChannels_P, nColorChannels_S, nColorChannels_P, nColorChannels_P};







	char nameModelN[500], nameModelD[500], StateNameOldConvention[200];

	if(iState==0) sprintf(StateNameOldConvention,"Psi1S");
	if(iState==3) sprintf(StateNameOldConvention,"Psi2S");
	if(iState==4) sprintf(StateNameOldConvention,"1SUps");
	if(iState==7) sprintf(StateNameOldConvention,"2SUps");
	if(iState==10) sprintf(StateNameOldConvention,"3SUps");

	if(iModelN==0) sprintf(nameModelN,"TGraph/CrossSections/TGraphs_BK70_Ep.root");
	if(iModelN==1) sprintf(nameModelN,"TGraph/CrossSections/TGraphs_GWWZ.root");
	if(iModelN==2) sprintf(nameModelN,"TGraph/CrossSections/TGraphs_BK70_LO.root");

	if(iModelD==0) sprintf(nameModelD,"TGraph/CrossSections/TGraphs_BK70_Ep.root");
	if(iModelD==1) sprintf(nameModelD,"TGraph/CrossSections/TGraphs_GWWZ.root");
	if(iModelD==2) sprintf(nameModelD,"TGraph/CrossSections/TGraphs_BK70_LO.root");


	cout<<"model file Numerator "<<nameModelN<<endl;
	cout<<"model file Denominator "<<nameModelD<<endl;

	TFile *fileModelN   = new TFile(nameModelN,"R");
	TFile *fileModelD   = new TFile(nameModelD,"R");

	char FrameChar[200];
	if(iFrame==0) sprintf(FrameChar, "HX");
	if(iFrame==1) sprintf(FrameChar, "CS");


	const int nStateMax=5;
	const int nRapModelMax=3;

	const int nColorChannelsConstp1=5;
	double LineWidthCC[nColorChannelsConstp1]={1.5, 1.5, 1.5, 1.5, 3};
	double LineStyleCC[nRapModelMax]={1, 2, 3};
	int ColorCC[nColorChannelsConstp1]={418, 632, 616, 807, 600};

	TGraphAsymmErrors* graphModelN[nRapModel+iRapModel][nColorChannels+1];
	TGraphAsymmErrors* graphModelD[nRapModel+iRapModel][nColorChannels+1];
	TGraph* graphModelR[nRapModel+iRapModel][nColorChannels+1];
	TGraph* graphModelDiff[nRapModel+iRapModel][nColorChannels+1];
	TGraph* graphModelDiffNeg[nRapModel+iRapModel][nColorChannels+1];

	double pTmin=pTMin;
	double pTmax=pTMax;

	/////////////



	double rapminModel[nRapModelMax]={0, 0.6, 2.};
	double rapmaxModel[nRapModelMax]={0.6, 1.2, 4.5};

	double rapScaleFactor[nRapModelMax]={1./1.2, 1./1.2, 1./2.5};


	char yAxisTitle[200];
	if(iMeasID==0) sprintf(yAxisTitle,"ratio of short distance coefficients");
	if(iMeasID==1) sprintf(yAxisTitle,"ratio of #lambda_{#vartheta}");

	//if(iMeasID==0) sprintf(yAxisTitle,"short distance coefficients, scaled");

	double buffxD, buffyD, buffxN, buffyN;

	double BorderBeforeEp=10;

	for(int iRap=iRapModel;iRap<nRapModel+iRapModel;iRap++){
		for(int iCC=0;iCC<nColorChannels;iCC++){
			char graphName[500];
			if(iMeasID==0) sprintf(graphName,"SDC_Graph_original_state%d_rap%d_CC%d", iState, iRap, iCC);
			if(iMeasID==1) sprintf(graphName,"Lamth_Graph_original_state%d_rap%d_CC%d", iState, iRap, iCC);
			if(iMeasID==2) sprintf(graphName,"Lamph_Graph_original_state%d_rap%d_CC%d", iState, iRap, iCC);
			if(iMeasID==3) sprintf(graphName,"Lamtp_Graph_original_state%d_rap%d_CC%d", iState, iRap, iCC);
			cout<<graphName<<endl;
			graphModelD[iRap][iCC] = (TGraphAsymmErrors*)fileModelD->Get(graphName);

			if(iMeasID==0) sprintf(graphName,"EpSDC_Graph_state%d_rap%d_CC%d", iState, iRap, iCC);
			if(iMeasID==1) sprintf(graphName,"EpLamth_Graph_state%d_rap%d_CC%d", iState, iRap, iCC);
			if(iMeasID==2) sprintf(graphName,"EpLamph_Graph_state%d_rap%d_CC%d", iState, iRap, iCC);
			if(iMeasID==3) sprintf(graphName,"EpLamtp_Graph_state%d_rap%d_CC%d", iState, iRap, iCC);
			cout<<graphName<<endl;
			graphModelN[iRap][iCC] = (TGraphAsymmErrors*)fileModelN->Get(graphName);

			//if(iCC==0){
			//graphModelN[iRap][iCC] = (TGraphAsymmErrors*)fileModelD->Get("SDC_Graph_original_state0_rap0_CC2");
			//graphModelD[iRap][iCC] = (TGraphAsymmErrors*)fileModelD->Get("SDC_Graph_original_state0_rap0_CC1");
			//}
			//if(iCC==1){
			//graphModelN[iRap][iCC] = (TGraphAsymmErrors*)fileModelN->Get("SDC_Graph_original_state0_rap0_CC2");
			//graphModelD[iRap][iCC] = (TGraphAsymmErrors*)fileModelN->Get("SDC_Graph_original_state0_rap0_CC1");
			//}

			//if(iCC==0){
			//	graphModelN[iRap][iCC] = (TGraphAsymmErrors*)fileModelN->Get("SDC_Graph_original_state0_rap0_CC1");
			//	graphModelD[iRap][iCC] = (TGraphAsymmErrors*)fileModelD->Get("SDC_Graph_original_state0_rap0_CC1");
			//}
			//if(iCC==1){
			//	graphModelN[iRap][iCC] = (TGraphAsymmErrors*)fileModelN->Get("SDC_Graph_original_state0_rap0_CC2");
			//	graphModelD[iRap][iCC] = (TGraphAsymmErrors*)fileModelD->Get("SDC_Graph_original_state0_rap0_CC2");
			//}
            //
			//if(iCC==2){
			//	graphModelN[iRap][iCC] = (TGraphAsymmErrors*)fileModelD->Get("SDC_Graph_original_state0_rap0_CC1");
			//	graphModelD[iRap][iCC] = (TGraphAsymmErrors*)fileModelD->Get("SDC_Graph_original_state0_rap0_CC1");
			//}
			//if(iCC==3){
			//	graphModelN[iRap][iCC] = (TGraphAsymmErrors*)fileModelD->Get("SDC_Graph_original_state0_rap0_CC2");
			//	graphModelD[iRap][iCC] = (TGraphAsymmErrors*)fileModelD->Get("SDC_Graph_original_state0_rap0_CC2");
			//}


			//Make ratio NNLO*/NLO in CS JPL paper (Fig. 5) to define the CS SDC THuncertainty:

			const int n_rCS=10;
			double rCSpT_Ups3S[n_rCS]={5., 10., 15., 20., 25., 30., 35., 40., 45., 50.};
			double rCSpT[n_rCS];
			for(int rCSi=0;rCSi<n_rCS;rCSi++){
				rCSpT[rCSi]=rCSpT_Ups3S[rCSi]/10.355*3;
			}
			double rCSd[n_rCS]={2.35, 3.65, 4.8, 5.6, 6.0, 6.6, 7.0, 7.35, 7.5, 7.75};
			TGraph *g_rCS = new TGraph(n_rCS, rCSpT, rCSd);
			g_rCS->Print();

			// Fit with some function

			char name[200];
			sprintf(name, "fParabola");
			TF1* fParabola = new TF1(name, paramLog, 1., 71., 3.);
			double a=0;
			double b=0;
			double c=1;
			fParabola->SetParameter(0,a);
			fParabola->SetParameter(1,b);
			fParabola->SetParameter(2,c);
			g_rCS->Fit(fParabola, "0", "", 1., 71.);

			fParabola->SetLineColor(kBlue);
			fParabola->SetLineWidth(1.);
			g_rCS->SetMarkerStyle(20);
			g_rCS->SetTitle(0);
			g_rCS->GetYaxis()->SetRangeUser(0,18);
			g_rCS->GetXaxis()->SetLimits(0.1,75.);
			g_rCS->GetXaxis()->SetTitle("p_{T}");
			g_rCS->GetYaxis()->SetTitle("r_{CS}");

			TCanvas *tmpC = new TCanvas("tmpC","tmpC",1000,800);
			tmpC->SetFillColor(kWhite);
			tmpC->SetFrameBorderMode(0);
			g_rCS->Draw("ap");
			fParabola->Draw("lsame");
			char tmpN[500];
			sprintf(tmpN,"Figures/Plot_HP_Data_vs_NRQCD/tmp_r_CS_fit.pdf");
			tmpC->SaveAs(tmpN);


			const int nBinsRatio=500;
			double buffArray[nBinsRatio];
			double RatioNorm=graphModelD[iRap][iCC]->Eval(pTMax)/graphModelN[iRap][iCC]->Eval(pTMax);

			double addRatioNorm=0.308;

			double RelativeTHuncertaintyAtBorder;
			RelativeTHuncertaintyAtBorder=graphModelN[iRap][iCC]->Eval(BorderBeforeEp)-graphModelD[iRap][iCC]->Eval(BorderBeforeEp);
			if(iCC==1&&iMeasID>0) RelativeTHuncertaintyAtBorder=0.;//1S0: no change from LO to NLO
			if(iCC==0&&iMeasID==1) RelativeTHuncertaintyAtBorder=graphModelN[iRap][iCC]->Eval(BorderBeforeEp)*0.2;//CS: constrain reduced by NNLO*

			RelativeTHuncertaintyAtBorder/=graphModelN[iRap][iCC]->Eval(BorderBeforeEp);

			if(!scaleToOne) RatioNorm=1;
			graphModelR[iRap][iCC] = new TGraph(nBinsRatio, buffArray, buffArray);
			graphModelDiff[iRap][iCC] = new TGraph(nBinsRatio, buffArray, buffArray);
			graphModelDiffNeg[iRap][iCC] = new TGraph(nBinsRatio, buffArray, buffArray);
			for(int i=0;i<nBinsRatio;i++){
				double pTRatio=pTMin+i*(pTMax-pTMin)/nBinsRatio;
				double Ratio=graphModelN[iRap][iCC]->Eval(pTRatio)/graphModelD[iRap][iCC]->Eval(pTRatio)*RatioNorm;
				double Diff=graphModelN[iRap][iCC]->Eval(pTRatio)-graphModelD[iRap][iCC]->Eval(pTRatio);
				if(iCC==1&&iMeasID>0) Ratio=1.;//1S0: no change from LO to NLO
				if(iCC==1&&iMeasID>0) Diff=0.;//1S0: no change from LO to NLO
				if(iCC==0&&iMeasID==1) Diff=graphModelN[iRap][iCC]->Eval(pTRatio)*0.2;//CS: constrain reduced by NNLO*

				//for Extrapolated NLO calculations, define THuncertainty such that the relative uncertainty is the same as at 10 GeV:
				if(pTRatio<BorderBeforeEp){
					Diff=RelativeTHuncertaintyAtBorder*graphModelN[iRap][iCC]->Eval(pTRatio);
					if(iCC==1&&iMeasID>0) Diff=0.;//1S0: no change from LO to NLO

				}

				double DiffNeg=-Diff;

				if(iCC==0&&iMeasID==0) {
					double rCS=fParabola->Eval(pTRatio);
					Diff=graphModelN[iRap][iCC]->Eval(pTRatio)*(rCS-1);
					//cout<<"pT "<<pTRatio<<endl;
					//cout<<"rCS "<<rCS<<endl;
					//cout<<"Diff "<<Diff<<endl;
					DiffNeg=-1.*(graphModelN[iRap][iCC]->Eval(pTRatio)-graphModelD[iRap][iCC]->Eval(pTRatio));
				}

				//Diff/=graphModelN[iRap][iCC]->Eval(pTRatio);

				//Ratio=graphModelN[iRap][iCC]->Eval(pTRatio)*RatioNorm;
				//if(iCC==1 || iCC==3) Ratio*=addRatioNorm;
				graphModelR[iRap][iCC]->SetPoint(i, pTRatio, Ratio);
				graphModelDiff[iRap][iCC]->SetPoint(i, pTRatio, Diff);
				graphModelDiffNeg[iRap][iCC]->SetPoint(i, pTRatio, DiffNeg);
			}

			cout<<"Calculated ratio of iCC"<<iCC<<endl;
			//graphModelR[iRap][iCC]->Print();

		}//CCModel
	}//RapModel

	for(int iRap=iRapModel;iRap<nRapModel+iRapModel;iRap++){
		for(int iCC=0;iCC<nColorChannels;iCC++){

			graphModelR[iRap][iCC]->GetYaxis()->SetRangeUser(ResMin, ResMax);
			graphModelR[iRap][iCC]->GetXaxis()->SetLimits(pTmin, pTmax);
			graphModelR[iRap][iCC]->SetTitle("");
			graphModelR[iRap][iCC]->GetXaxis()->SetTitle("#it{p}_{T} [GeV]");
			graphModelR[iRap][iCC]->GetYaxis()->SetTitle(yAxisTitle);
			graphModelR[iRap][iCC]->GetYaxis()->SetTitleSize(0.05);
			graphModelR[iRap][iCC]->GetYaxis()->SetTitleOffset(1.);
			graphModelR[iRap][iCC]->SetLineColor(ColorCC[iCC]);
			graphModelR[iRap][iCC]->SetLineWidth(LineWidthCC[iCC]);
			graphModelR[iRap][iCC]->SetLineStyle(LineStyleCC[iRap]);

		}//CCModel
	}//rapModel





		TCanvas *plotCanvas = new TCanvas("plotCanvas","plotCanvas",1000,800);
		plotCanvas->SetFillColor(kWhite);
		plotCanvas->SetFrameBorderMode(0);

		TLegend* legend;
		if(LegendPos==1) legend = new TLegend(0.15,0.6010989,0.35,0.95008242);
		if(LegendPos==2) legend = new TLegend(0.65,0.6010989,0.95,0.95008242);
		if(LegendPos==3) legend = new TLegend(0.65,0.1310989,0.95,0.48008242);
		if(LegendPos==4) legend = new TLegend(0.15,0.1310989,0.35,0.48008242);
		legend->SetFillColor(0);
		legend->SetTextSize(0.025);
		legend->SetBorderSize(0);

		char DrawChar[200];
		char legendentry[200];


		bool IsAxisDrawnAlready=false;

		for(int iRap=iRapModel;iRap<nRapModel+iRapModel;iRap++){
			for(int iCC=0;iCC<nColorChannels;iCC++){

					if(IsAxisDrawnAlready)
					sprintf(DrawChar,"l");
					else{
						sprintf(DrawChar,"al");
						IsAxisDrawnAlready=true;
					}
					graphModelR[iRap][iCC]->Draw(DrawChar);
					sprintf(legendentry,"%s, %s", ColorChannelNameTexS[iCC], ModelIDNModelIDD);
					//sprintf(legendentry,"%s", ColorChannelNameTexS[iCC]);
					legend->AddEntry(graphModelR[iRap][iCC],legendentry,"l");

				}


			}





		TLine* extreme0 = new TLine( pTmin, 1., pTmax, 1.);
		extreme0->SetLineWidth( 1 );
		extreme0->SetLineStyle( 2 );
		extreme0->SetLineColor( kBlack );
		extreme0->Draw( "same" );

		legend->Draw("same");

		//draw latex
		double left=0.175, top=0.92, textSize=0.045;
		TLatex *latex=new TLatex();
		latex->SetTextFont(42);
		latex->SetNDC(kTRUE);
		latex->SetTextSize(textSize);
		double stepLatex=textSize*1.3;
		latex->SetTextFont(62);

		if(LatexPos==1 || LatexPos==4) left=0.175;
		if(LatexPos==2 || LatexPos==3) left=0.78;

		if(LatexPos==1 || LatexPos==2) top=0.90;
		if(LatexPos==3 || LatexPos==4) top=0.27;

		textSize=0.06;
		latex->SetTextSize(textSize);
		latex->DrawLatex(left,top, Form("%s",StateNameTex[iState]));

		if(LatexPos==1 || LatexPos==2) top=0.84;
		if(LatexPos==3 || LatexPos==4) top=0.21;
		textSize=0.045;
		latex->SetTextSize(textSize);
		if(rapminModel[iRapModel]<1e-5) latex->DrawLatex(left,top, Form("|y| < %1.1f", rapmaxModel[iRapModel]));
		else latex->DrawLatex(left,top, Form("%1.1f < |y| < %1.1f", rapminModel[iRapModel], rapmaxModel[iRapModel]));

		if(iMeasID>0){
			if(LatexPos==1 || LatexPos==2) top=0.79;
			if(LatexPos==3 || LatexPos==4) top=0.16;
			textSize=0.045;
			latex->SetTextSize(textSize);
			latex->DrawLatex(left,top, Form("%s frame",FrameChar));
		}





		char MeasurementIDNameFileName[200];
		sprintf(MeasurementIDNameFileName,"%s",MeasurementIDName[iMeasID]);
		if(iMeasID>0)sprintf(MeasurementIDNameFileName,"%s_%s",MeasurementIDName[iMeasID], FrameChar);


		char saveName[500];
		sprintf(saveName,"Figures/Plot_HP_Data_vs_NRQCD/NRQCD_Ratios_%s_over_%s_%s_%s_nRapModel%d_iRapModel%d.pdf", ModelIDN, ModelIDD, StateName[iState], MeasurementIDNameFileName, nRapModel, iRapModel);
		cout<<saveName<<endl;
		//plotCanvas->SetLogy(true);
		plotCanvas->SaveAs(saveName);




		///// THIS PART IS HERE TO SAVE TGRAPHS FOR THE THEORETICAL UNCERTAINTY NLO-LO:::

		//scale graphModelR[iRap][iCC] and save TGraphs, for all states, one rap, all iCCs

		double massToScale=3;

		for(int iCC=0;iCC<nColorChannels;iCC++){

			double massScale[nStates];
			double NormalizationScale[nStates];

			//open file
			char outfilename[200];
			sprintf(outfilename,"TGraph/ModelSystematicScales/TGraphs_ModelSystematicScale%d_iMeasurementID%d.root", iCC, iMeasID);
			cout<<outfilename<<endl;
			TFile *outfile = new TFile(outfilename,"RECREATE");


			for(int iState=0;iState<nStates;iState++){


				massScale[iState]=mass[iState]/massToScale;

				double expnorm = 4.03556;
				double exponent = -0.940807;

				double massmin=1;
				double massmax=13;

				char name[200];
				sprintf(name, "expo_fit");
				TF1* expo_fit  = new TF1(name, exponential, massmin, massmax, 2);
				expo_fit->FixParameter(0,expnorm);
				expo_fit->FixParameter(1,exponent);

				NormalizationScale[iState]=expo_fit->Eval(mass[iState])/expo_fit->Eval(massToScale);


				TGraph* graphModelDiffScaledPos;
				TGraph* graphModelDiffScaledNeg;
				const int nPt=graphModelR[0][iCC]->GetN();

				double pTmean_graph[nPt];
				double pTmeanScaled_graph[nPt];
				double Diff_graph_Pos[nPt];
				double Diff_graph_Neg[nPt];
				for(int m=0;m<nPt;m++){
					graphModelDiff[0][iCC]->GetPoint(m,pTmean_graph[m],Diff_graph_Pos[m]);
					graphModelDiffNeg[0][iCC]->GetPoint(m,pTmean_graph[m],Diff_graph_Neg[m]);
					//apply normalization scale, and rap phase space correction (for cross section only)
					if(iMeasID==0){
						Diff_graph_Pos[m]*=NormalizationScale[iState]*rapScaleFactor[0];
						Diff_graph_Neg[m]*=NormalizationScale[iState]*rapScaleFactor[0];
					}
					if(iMeasID>1){
						Diff_graph_Pos[m]=0.;
						Diff_graph_Neg[m]=0.;
					}
					//Diff_graph_Neg[m]=-Diff_graph_Pos[m];

					//apply pT scale
					pTmeanScaled_graph[m]=pTmean_graph[m]*massScale[iState];
					}

				graphModelDiffScaledPos = new TGraph(nPt,pTmeanScaled_graph,Diff_graph_Pos);
				graphModelDiffScaledNeg = new TGraph(nPt,pTmeanScaled_graph,Diff_graph_Neg);

				char ModelSystematicScaleName[200];
				sprintf(ModelSystematicScaleName,"ModelSystematicScale_iState%d_Pos",iState);
				graphModelDiffScaledPos->SetName(ModelSystematicScaleName);
				sprintf(ModelSystematicScaleName,"ModelSystematicScale_iState%d_Neg",iState);
				graphModelDiffScaledNeg->SetName(ModelSystematicScaleName);

				double extrapolateTo=2.85;
				if(iMeasID>0){

					int nGraph=graphModelDiffScaledPos->GetN();
					int jTG=0;
					for(int j=0;j<nGraph;j++){
						double buffx, buffy;
						graphModelDiffScaledPos->GetPoint(jTG, buffx, buffy);
						if(buffx<BorderBeforeEp){
							graphModelDiffScaledPos->RemovePoint(0);
							graphModelDiffScaledNeg->RemovePoint(0);
							jTG--;
						}

						jTG++;
					}

					if(iCC==2){
						graphModelDiffScaledPos->RemovePoint(0);
						graphModelDiffScaledNeg->RemovePoint(0);
						graphModelDiffScaledPos->RemovePoint(0);
						graphModelDiffScaledNeg->RemovePoint(0);
						graphModelDiffScaledPos->RemovePoint(0);
						graphModelDiffScaledNeg->RemovePoint(0);
					}

					graphModelDiffScaledPos->SetPoint(0,extrapolateTo,graphModelDiffScaledPos->Eval(extrapolateTo));
					graphModelDiffScaledNeg->SetPoint(0,extrapolateTo,graphModelDiffScaledNeg->Eval(extrapolateTo));
				}
				outfile->cd();
				graphModelDiffScaledPos->Write();
				graphModelDiffScaledNeg->Write();


			}

			//close, write file
			outfile->Write();
			outfile->Close();
			delete outfile;
			outfile = NULL;

		}

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

Double_t paramLog(Double_t *x, Double_t *par){
	Double_t result=par[0]+TMath::Log(x[0]-par[1])*par[2];
	return result;
}


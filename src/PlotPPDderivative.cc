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

	cout<<"start PlotPPDderivative"<<endl;

  	Char_t *JobID = "Default";
	Char_t *storagedir = "Default"; //Storage Directory
	int useOnlyState=999;
  	int 	MPValgo=-1;
  	double 	nSigma=-1;

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
	sprintf(plotdirname,"Figures/PlotPPDderivative/%s",JobID);
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


	//char JobID2[200];
	//char JobID3[200];
	//sprintf(JobID2, "October16_BKmodel_ScaledAllStates_psi2S_1S03S1_pTstar_over_m6_MH_pTmin13_POLonly");
	//sprintf(JobID3, "October16_BKmodel_ScaledAllStates_psi2S_1S03S1_pTstar_over_m6_MH_pTmin13_CSonly");

	sprintf(inname,"%s/JobID/%s/results.root",storagedir, JobID);
	TFile *ResultsFile = new TFile(inname, "READ");
  	TTree*  outputTreeAllSamplings = (TTree*) ResultsFile->Get("AllSamplings"); // tree filled in all samplings after burnin

	//sprintf(inname,"%s/JobID/%s/results.root",storagedir, JobID2);
	//TFile *ResultsFile2 = new TFile(inname, "READ");
  	//TTree*  outputTreeAllSamplings2 = (TTree*) ResultsFile2->Get("AllSamplings"); // tree filled in all samplings after burnin
	//sprintf(inname,"%s/JobID/%s/results.root",storagedir, JobID3);
	//TFile *ResultsFile3 = new TFile(inname, "READ");
  	//TTree*  outputTreeAllSamplings3 = (TTree*) ResultsFile3->Get("AllSamplings"); // tree filled in all samplings after burnin



	////////////////////////////////////////
	/// Plot 1D ratio of Ops
	////////////////////////////////////////


			char branch_name[200];
			char branch_name_O1[200];
			char branch_name_O2[200];
			char hist_var_name_f1[200];
			char hist_var_name_f2[200];


			int i=useOnlyState;

			int j=2;//Numerator in ratio
			int k=1;//Denominator in ratio


			cout<<"j"<<j<<"k"<<k<<endl;
			sprintf(branch_name_O1,"state%d_Op%d",i,j);
			sprintf(branch_name_O2,"state%d_Op%d",i,k);

			double O1;
			double O2;
			double OpRatio;

	      	int BurnInInt;
			sprintf (branch_name,"BurnInInt");
			outputTreeAllSamplings->SetBranchAddress( branch_name,  &BurnInInt );
	      	int acceptedSampling;
			sprintf (branch_name,"acceptedSampling");
			outputTreeAllSamplings->SetBranchAddress( branch_name,  &acceptedSampling );


			outputTreeAllSamplings->SetBranchAddress( branch_name_O1,  &O1 );
			outputTreeAllSamplings->SetBranchAddress( branch_name_O2,  &O2 );

			sprintf(outname,"%s/RatioTree.root", jobdirname);
	    	TFile *DummyFile = new TFile(outname, "RECREATE");

	      	TTree*  RatioTree;
	      	RatioTree = new TTree("RatioTree", "RatioTree");
	      	RatioTree->Branch("OpRatio",     &OpRatio,     "OpRatio/D");

			int n_events = int( outputTreeAllSamplings->GetEntries() );
			cout<<n_events<<" events: Looping through nTuple of PPD"<<endl;

			// loop over  events in the model ntuple
			for ( int i_event = 1; i_event <= n_events; i_event++ ) {


				outputTreeAllSamplings->GetEvent( i_event-1 );

				if(acceptedSampling!=1 || BurnInInt==1) continue;
				//if(O1<0 || O2<0) continue;

				OpRatio=O1/O2;
				RatioTree->Fill();


			}

			RatioTree->Print();


			int nBins_h=200;
			char hist_name[200];
			char projectchar[200];
			char selectchar[200];
			TH1D* h_Ratio;
			double h_Ratio_min;
			double h_Ratio_max;

		  	double expandMinMaxBy=0.01;
		  	double buff_MPV;
		  	double buff_errlow;
		  	double buff_errhigh;

		  	h_Ratio_min=RatioTree->GetMinimum("OpRatio")-expandMinMaxBy*RatioTree->GetMinimum("OpRatio");
		  	h_Ratio_max=RatioTree->GetMaximum("OpRatio")+expandMinMaxBy*RatioTree->GetMaximum("OpRatio");

		  	double minmax_min=-5e-2;
		  	double minmax_max=3e-1;
		  	//if(h_Ratio_max>minmax_max) h_Ratio_max=minmax_max;
		  	//if(h_Ratio_min<minmax_min) h_Ratio_min=minmax_min;

		  	sprintf(hist_name,"h_Ratio");
		  	h_Ratio = new TH1D( hist_name, hist_name, nBins_h,h_Ratio_min, h_Ratio_max );
			sprintf(projectchar,"OpRatio>>%s",hist_name);
			RatioTree->Draw(projectchar);


			FindMPV(h_Ratio, buff_MPV, buff_errlow, buff_errhigh, MPValgo, nSigma);

			h_Ratio->Print("all");

			char filename[200];
			sprintf(filename,"Figures/PlotPPDderivative/%s",JobID);
			gSystem->mkdir(filename);
			sprintf(filename,"%s/PPD_LDME_ratio_state%d_Op%d_ov_Op%d.pdf",filename, i, j, k);

			cout<<filename<<endl;

			char hist_var_name[200];
			sprintf(hist_var_name,"O(%s) / O(%s)", ColorChannelNameTexS[j], ColorChannelNameTexS[k]);
			h_Ratio -> SetXTitle(hist_var_name);

			PlotPosterior(hist_var_name, filename, h_Ratio, buff_MPV, buff_errlow, buff_errhigh, nSigma);

			double buff_MPV_lowborder, buff_MPV_highborder;
			double buff_MPV_inv, buff_MPV_inv_lowborder, buff_MPV_inv_highborder;

			buff_MPV_lowborder=buff_MPV-buff_errlow;
			buff_MPV_highborder=buff_MPV+buff_errhigh;


			buff_MPV_inv=1./buff_MPV;
			buff_MPV_inv_lowborder=1./buff_MPV_highborder;
			buff_MPV_inv_highborder=1./buff_MPV_lowborder;

			//if(buff_MPV_inv_highborder<0) buff_MPV_inv_highborder=1e99;

			cout<<StateName[i]<<", "<<nSigma<<"sigma CL:"<<endl;
			cout<<ColorChannelNameTexS[j]<<" / "<<ColorChannelNameTexS[k]<<" = (MPV) "<<buff_MPV<<" = (CL) ["<<buff_MPV_lowborder<<", "<<buff_MPV_highborder<<"]"<<endl;
			cout<<ColorChannelNameTexS[k]<<" / "<<ColorChannelNameTexS[j]<<" = (MPV) "<<buff_MPV_inv<<" = (CL) ["<<buff_MPV_inv_lowborder<<", "<<buff_MPV_inv_highborder<<"]"<<endl;







			gROOT->Reset();
			gStyle->SetOptStat(11);
			gStyle->SetOptFit(101);
			gStyle->SetFillColor(kWhite);

			gStyle->SetPadBottomMargin(0.12);
			gStyle->SetPadLeftMargin(0.12);
			gStyle->SetPadRightMargin(0.02);
			gStyle->SetPadTopMargin(0.02);





			/////////// Pol-predictions plots



			double MassScale;
			double MassToScale=3.;

			MassScale=NRQCDvars::mass[i]/MassToScale;

			bool functionOfPT=true;

			bool plotAllData=true;
			bool plotModel=true;

			bool plotBK=false;
			if(i==3 && plotModel) plotBK=true;
			bool plotGong=false;
			if(i==10 && plotModel) plotGong=true;

			int nxVar=13;
			double xVarminPred=1.;
			double xVarmaxPred=15;

			if(functionOfPT){
				xVarminPred*=NRQCDvars::mass[i];
				xVarmaxPred*=NRQCDvars::mass[i];
			}

			xVarminPred=10;
			xVarmaxPred=120;

			//if(!functionOfPT){
			//	xVarminPred=1;
			//	xVarmaxPred=35;
			//}
			//else{
			//	xVarminPred=10;
			//	xVarmaxPred=150;
			//}

			double xVarMean[nxVar];
			double zero[nxVar];
			double PredLamth[nxVar];
			double PredLamth_lowborder[nxVar];
			double PredLamth_highborder[nxVar];
			double PredLamth_errlow[nxVar];
			double PredLamth_errhigh[nxVar];

			double PredLamth_errlow_1sig[nxVar];
			double PredLamth_errhigh_1sig[nxVar];
			double PredLamth_errlow_2sig[nxVar];
			double PredLamth_errhigh_2sig[nxVar];
			double PredLamth_errlow_3sig[nxVar];
			double PredLamth_errhigh_3sig[nxVar];

			for(int ixVar=0;ixVar<nxVar;ixVar++){

				double pTMean;
				double pTMeanRaw;

				xVarMean[ixVar]=(xVarmaxPred-xVarminPred)/(nxVar-1)*ixVar+xVarminPred;

				if(functionOfPT) pTMean=xVarMean[ixVar];
				else pTMean=xVarMean[ixVar]*NRQCDvars::mass[i];

				pTMeanRaw=pTMean/MassScale;

				char name[200];
				sprintf(name,"fPtDist");
				TF1* fPtDist = new TF1(name, func_pT_ratio_NLO, xVarminPred, xVarmaxPred, 4);
				//
				////norm * [(gamma1 + pT^2)/(gamma2 + pT^2)]^(-beta)
				//

				// Valid for pT>45GeV
				//  1  p0           2.55717e+02   3.98514e+00   4.29487e-04   3.36349e-05
				//  2  p1           4.00000e+00   4.25270e-01   1.59050e-05  -9.81018e-04
				//  3  p2           1.11429e+03   5.33368e+01   1.95649e-03  -7.94055e-06
				//  4  p3           6.48915e+02   4.54545e+01   1.74703e-03   8.96707e-06

				// Valid for pT<45GeV
				// 1  p0           2.10170e+02   4.04053e+00   1.49025e-03   1.55868e-06
				// 2  p1           4.00000e+00   2.72604e-01   2.94104e-05   3.78117e-01
				// 3  p2           4.74459e+02   2.08775e+01   2.40173e-03  -1.30859e-06
				// 4  p3           1.94259e+02   1.27278e+01   1.72845e-03   2.04339e-06

				double func_pT_ratio_NLO_pars_lowpT[4]={2.10170e+02, 4.00000e+00, 4.74459e+02, 1.94259e+02};
				double func_pT_ratio_NLO_pars_highpT[4]={2.55717e+02, 4.00000e+00, 1.11429e+03, 6.48915e+02};

				double norm;
				double beta;
				double gamma1;
				double gamma2;

				double lowPtBoder=45;

				if(pTMeanRaw<lowPtBoder){
					norm=func_pT_ratio_NLO_pars_lowpT[0];
					beta=func_pT_ratio_NLO_pars_lowpT[1];
					gamma1=func_pT_ratio_NLO_pars_lowpT[2];
					gamma2=func_pT_ratio_NLO_pars_lowpT[3];
				}
				else{
					norm=func_pT_ratio_NLO_pars_highpT[0];
					beta=func_pT_ratio_NLO_pars_highpT[1];
					gamma1=func_pT_ratio_NLO_pars_highpT[2];
					gamma2=func_pT_ratio_NLO_pars_highpT[3];
				}

				fPtDist->SetParameter(0,norm);
				fPtDist->SetParameter(1,beta);
				fPtDist->SetParameter(2,gamma1);
				fPtDist->SetParameter(3,gamma2);

				//SDC is defined!


				//sprintf(outname,"%s/LamthTree.root", jobdirname);
		    	//TFile *DummyFile = new TFile(outname, "RECREATE");

		    	double LamthPred;
		      	TTree*  LamthTree;
		      	LamthTree = new TTree("LamthTree", "LamthTree");
		      	LamthTree->Branch("LamthPred",     &LamthPred,     "LamthPred/D");

				sprintf (branch_name,"OpRatio");
				RatioTree->SetBranchAddress( branch_name,  &OpRatio );

				n_events = int( RatioTree->GetEntries() );
				cout<<n_events<<" events: Looping through nTuple of PPD"<<endl;

				// loop over  events in the model ntuple
				for ( int i_event = 1; i_event <= n_events; i_event++ ) {


					RatioTree->GetEvent( i_event-1 );



					double SDC_R;

					SDC_R=fPtDist->Eval(pTMeanRaw);


					double LDME_R=OpRatio;

					//cout<<"SDC_R "<<SDC_R<<endl;
					//cout<<"LDME_R "<<LDME_R<<endl;

					vector<vector<double> > lamMatrix;
					vector<double> lamMatrixCont1;
					vector<double> lamMatrixCont2;

					lamMatrixCont1.push_back(0.);
					lamMatrixCont1.push_back(0.);
					lamMatrixCont1.push_back(0.);
					lamMatrixCont2.push_back(1.);
					lamMatrixCont2.push_back(0.);
					lamMatrixCont2.push_back(0.);

					lamMatrix.push_back(lamMatrixCont1);
					lamMatrix.push_back(lamMatrixCont2);


					vector<double> lamVecContributions;

					lamVecContributions.push_back(1.);
					lamVecContributions.push_back(SDC_R*LDME_R);

					vector<double> LamthPredVec = addPolarizations(lamMatrix, lamVecContributions);

					//cout<<lamMatrix<<endl;
					//cout<<lamVecContributions<<endl;
					//cout<<LamthPredVec<<endl;

					LamthPred=LamthPredVec[0];

					//cout<<LamthPred<<endl;

					LamthTree->Fill();


				}

				LamthTree->Print();



				TH1D* h_Lamth;
				double h_Lamth_min;
				double h_Lamth_max;


			  	h_Lamth_min=LamthTree->GetMinimum("LamthPred")-expandMinMaxBy*LamthTree->GetMinimum("LamthPred");
			  	h_Lamth_max=LamthTree->GetMaximum("LamthPred")+expandMinMaxBy*LamthTree->GetMaximum("LamthPred");

			  	double minmax_min=-1.;
			  	double minmax_max=1.;
			  	h_Lamth_max=minmax_max;
			  	h_Lamth_min=minmax_min;

			  	sprintf(hist_name,"h_Lamth");
			  	h_Lamth = new TH1D( hist_name, hist_name, nBins_h,h_Lamth_min, h_Lamth_max );
				sprintf(projectchar,"LamthPred>>%s",hist_name);
				LamthTree->Draw(projectchar);


				FindMPV(h_Lamth, buff_MPV, buff_errlow, buff_errhigh, MPValgo, nSigma);

				//h_Lamth->Print("all");

				char filename[200];
				sprintf(filename,"Figures/PlotPPDderivative/%s",JobID);
				gSystem->mkdir(filename);
				sprintf(filename,"%s/PPD_LamthPred_state%d_xVar%d.pdf",filename,i, ixVar);

				cout<<filename<<endl;

				char hist_var_name[200];
				sprintf(hist_var_name,"#lambda_{#vartheta}^{HX}");
				h_Lamth->SetXTitle(hist_var_name);
				h_Lamth->SetTitle("");

				PlotPosterior(hist_var_name, filename, h_Lamth, buff_MPV, buff_errlow, buff_errhigh, nSigma);

				buff_MPV_lowborder=buff_MPV-buff_errlow;
				buff_MPV_highborder=buff_MPV+buff_errhigh;

				PredLamth_lowborder[ixVar]=buff_MPV_lowborder;
				PredLamth_highborder[ixVar]=buff_MPV_highborder;

				PredLamth_errlow[ixVar]=buff_errlow;
				PredLamth_errhigh[ixVar]=buff_errhigh;

				PredLamth[ixVar]=buff_MPV;

				zero[ixVar]=0;

				FindMPV(h_Lamth, buff_MPV, buff_errlow, buff_errhigh, MPValgo, 1);
				PredLamth_errlow_1sig[ixVar]=buff_errlow;
				PredLamth_errhigh_1sig[ixVar]=buff_errhigh;
				FindMPV(h_Lamth, buff_MPV, buff_errlow, buff_errhigh, MPValgo, 2);
				PredLamth_errlow_2sig[ixVar]=buff_errlow;
				PredLamth_errhigh_2sig[ixVar]=buff_errhigh;
				FindMPV(h_Lamth, buff_MPV, buff_errlow, buff_errhigh, MPValgo, 3);
				PredLamth_errlow_3sig[ixVar]=buff_errlow;
				PredLamth_errhigh_3sig[ixVar]=buff_errhigh;

				LamthTree->Write();

			}





			TGraph *LamthPredGraph = new TGraph(nxVar,xVarMean,PredLamth);
			TGraph *LamthPredGraph_lowborder = new TGraph(nxVar,xVarMean,PredLamth_lowborder);
			TGraph *LamthPredGraph_highborder = new TGraph(nxVar,xVarMean,PredLamth_highborder);

			LamthPredGraph->Print();
			LamthPredGraph_lowborder->Print();
			LamthPredGraph_highborder->Print();



			TGraphAsymmErrors *LamthPredGraph1sig = new TGraphAsymmErrors(nxVar,xVarMean,PredLamth, zero, zero, PredLamth_errlow_1sig, PredLamth_errhigh_1sig);
			TGraphAsymmErrors *LamthPredGraph2sig = new TGraphAsymmErrors(nxVar,xVarMean,PredLamth, zero, zero, PredLamth_errlow_2sig, PredLamth_errhigh_2sig);
			TGraphAsymmErrors *LamthPredGraph3sig = new TGraphAsymmErrors(nxVar,xVarMean,PredLamth, zero, zero, PredLamth_errlow_3sig, PredLamth_errhigh_3sig);

			//getData


			char TGdir[200];
			char nameData[200];
			sprintf(TGdir, "/afs/hephy.at/scratch/k/knuenz/NRQCDglobalfit");

			sprintf(nameData,"%s/TGraph/CMS/TGraphResults_Psi1S.root", TGdir);
			TFile *fileData_Psi1S   = new TFile(nameData,"R");
			sprintf(nameData,"%s/TGraph/CMS/TGraphResults_Psi2S.root", TGdir);
			TFile *fileData_Psi2S   = new TFile(nameData,"R");
			sprintf(nameData,"%s/TGraph/CMS/TGraphResults_1SUps.root", TGdir);
			TFile *fileData_Ups1S   = new TFile(nameData,"R");
			sprintf(nameData,"%s/TGraph/CMS/TGraphResults_2SUps.root", TGdir);
			TFile *fileData_Ups2S   = new TFile(nameData,"R");
			sprintf(nameData,"%s/TGraph/CMS/TGraphResults_3SUps.root", TGdir);
			TFile *fileData_Ups3S   = new TFile(nameData,"R");

			sprintf(nameData,"%s/TGraph/NRQCD/TGraphResults_3SUps.root", TGdir);
			TFile *fileGong   = new TFile(nameData,"R");

			TGraphAsymmErrors* graphData[2][13];


			int MarkerStyleData[2]={20, 24};
			double MarkerSizeData=1.5;
			int ColorData[13]={1, 0,0, 600, 418, 0,0, 616, 0,0, 632, 0,0};

			for(int iRap=0;iRap<2;iRap++){
					char graphName[500];
					sprintf(graphName,"lth_HX_rap%d", iRap+1);
					graphData[iRap][0] = (TGraphAsymmErrors*)fileData_Psi1S->Get(graphName);
					graphData[iRap][3] = (TGraphAsymmErrors*)fileData_Psi2S->Get(graphName);
					graphData[iRap][4] = (TGraphAsymmErrors*)fileData_Ups1S->Get(graphName);
					graphData[iRap][7] = (TGraphAsymmErrors*)fileData_Ups2S->Get(graphName);
					graphData[iRap][10] = (TGraphAsymmErrors*)fileData_Ups3S->Get(graphName);

					if(!plotAllData){
						graphData[iRap][10]->RemovePoint(0);
						graphData[iRap][10]->RemovePoint(0);
						graphData[iRap][10]->RemovePoint(0);
						graphData[iRap][10]->RemovePoint(0);
						graphData[iRap][10]->RemovePoint(0);
						graphData[iRap][10]->RemovePoint(0);
						graphData[iRap][10]->RemovePoint(0);
						graphData[iRap][10]->RemovePoint(0);
						graphData[iRap][10]->RemovePoint(0);
					}

					int iStatebuff;

					iStatebuff=0;
					graphData[iRap][iStatebuff]->SetMarkerStyle(MarkerStyleData[iRap]);
					graphData[iRap][iStatebuff]->SetMarkerSize(MarkerSizeData);
					graphData[iRap][iStatebuff]->SetMarkerColor(ColorData[iStatebuff]);
					graphData[iRap][iStatebuff]->SetLineColor(ColorData[iStatebuff]);
					for(int m=0;m<graphData[iRap][iStatebuff]->GetN();m++){
						double buffx, buffy;
						graphData[iRap][iStatebuff]->GetPoint(m, buffx, buffy);
						buffx=buffx/NRQCDvars::mass[iStatebuff];
						if(!functionOfPT) graphData[iRap][iStatebuff]->SetPoint(m, buffx, buffy);
						if(!functionOfPT) graphData[iRap][iStatebuff]->SetPointEXhigh(m, 	graphData[iRap][iStatebuff]->GetErrorXhigh(m)/NRQCDvars::mass[iStatebuff]);
						if(!functionOfPT) graphData[iRap][iStatebuff]->SetPointEXlow(m, 	graphData[iRap][iStatebuff]->GetErrorXlow(m)/NRQCDvars::mass[iStatebuff]);
					}

					iStatebuff=3;
					graphData[iRap][iStatebuff]->SetMarkerStyle(MarkerStyleData[iRap]);
					graphData[iRap][iStatebuff]->SetMarkerSize(MarkerSizeData);
					graphData[iRap][iStatebuff]->SetMarkerColor(ColorData[iStatebuff]);
					graphData[iRap][iStatebuff]->SetLineColor(ColorData[iStatebuff]);
					for(int m=0;m<graphData[iRap][iStatebuff]->GetN();m++){
						double buffx, buffy;
						graphData[iRap][iStatebuff]->GetPoint(m, buffx, buffy);
						buffx=buffx/NRQCDvars::mass[iStatebuff];
						if(!functionOfPT) graphData[iRap][iStatebuff]->SetPoint(m, buffx, buffy);
						if(!functionOfPT) graphData[iRap][iStatebuff]->SetPointEXhigh(m, 	graphData[iRap][iStatebuff]->GetErrorXhigh(m)/NRQCDvars::mass[iStatebuff]);
						if(!functionOfPT) graphData[iRap][iStatebuff]->SetPointEXlow(m, 	graphData[iRap][iStatebuff]->GetErrorXlow(m)/NRQCDvars::mass[iStatebuff]);
					}

					iStatebuff=4;
					graphData[iRap][iStatebuff]->SetMarkerStyle(MarkerStyleData[iRap]);
					graphData[iRap][iStatebuff]->SetMarkerSize(MarkerSizeData);
					graphData[iRap][iStatebuff]->SetMarkerColor(ColorData[iStatebuff]);
					graphData[iRap][iStatebuff]->SetLineColor(ColorData[iStatebuff]);
					for(int m=0;m<graphData[iRap][iStatebuff]->GetN();m++){
						double buffx, buffy;
						graphData[iRap][iStatebuff]->GetPoint(m, buffx, buffy);
						buffx=buffx/NRQCDvars::mass[iStatebuff];
						if(!functionOfPT) graphData[iRap][iStatebuff]->SetPoint(m, buffx, buffy);
						if(!functionOfPT) graphData[iRap][iStatebuff]->SetPointEXhigh(m, 	graphData[iRap][iStatebuff]->GetErrorXhigh(m)/NRQCDvars::mass[iStatebuff]);
						if(!functionOfPT) graphData[iRap][iStatebuff]->SetPointEXlow(m, 	graphData[iRap][iStatebuff]->GetErrorXlow(m)/NRQCDvars::mass[iStatebuff]);
					}

					iStatebuff=7;
					graphData[iRap][iStatebuff]->SetMarkerStyle(MarkerStyleData[iRap]);
					graphData[iRap][iStatebuff]->SetMarkerSize(MarkerSizeData);
					graphData[iRap][iStatebuff]->SetMarkerColor(ColorData[iStatebuff]);
					graphData[iRap][iStatebuff]->SetLineColor(ColorData[iStatebuff]);
					for(int m=0;m<graphData[iRap][iStatebuff]->GetN();m++){
						double buffx, buffy;
						graphData[iRap][iStatebuff]->GetPoint(m, buffx, buffy);
						buffx=buffx/NRQCDvars::mass[iStatebuff];
						if(!functionOfPT) graphData[iRap][iStatebuff]->SetPoint(m, buffx, buffy);
						if(!functionOfPT) graphData[iRap][iStatebuff]->SetPointEXhigh(m, 	graphData[iRap][iStatebuff]->GetErrorXhigh(m)/NRQCDvars::mass[iStatebuff]);
						if(!functionOfPT) graphData[iRap][iStatebuff]->SetPointEXlow(m, 	graphData[iRap][iStatebuff]->GetErrorXlow(m)/NRQCDvars::mass[iStatebuff]);
					}

					iStatebuff=10;
					graphData[iRap][iStatebuff]->SetMarkerStyle(MarkerStyleData[iRap]);
					graphData[iRap][iStatebuff]->SetMarkerSize(MarkerSizeData);
					graphData[iRap][iStatebuff]->SetMarkerColor(ColorData[iStatebuff]);
					graphData[iRap][iStatebuff]->SetLineColor(ColorData[iStatebuff]);
					for(int m=0;m<graphData[iRap][iStatebuff]->GetN();m++){
						double buffx, buffy;
						graphData[iRap][iStatebuff]->GetPoint(m, buffx, buffy);
						buffx=buffx/NRQCDvars::mass[iStatebuff];
						if(!functionOfPT) graphData[iRap][iStatebuff]->SetPoint(m, buffx, buffy);
						if(!functionOfPT) graphData[iRap][iStatebuff]->SetPointEXhigh(m, 	graphData[iRap][iStatebuff]->GetErrorXhigh(m)/NRQCDvars::mass[iStatebuff]);
						if(!functionOfPT) graphData[iRap][iStatebuff]->SetPointEXlow(m, 	graphData[iRap][iStatebuff]->GetErrorXlow(m)/NRQCDvars::mass[iStatebuff]);
					}


			}//rapData



			//BK model lamth for |y|<0.6, but practically identical for 0.6<|y|<1.2

			 const int BK_npT=12;
			 double BK_pT[BK_npT]={10.000000, 12.500000, 15.000000, 17.500000, 20.000000, 25.000000, 30.000000, 35.000000, 40.000000, 50.000000, 60.000000, 70.000000};
			 double BK_model_low[BK_npT]={3.1748E-01, 4.8479E-01, 5.9863E-01, 6.7580E-01, 7.2974E-01, 7.9861E-01, 8.3635E-01, 8.7283E-01, 8.9537E-01, 9.3761E-01, 9.6689E-01, 9.5996E-01};
			 double BK_model[BK_npT]={6.5545E-01, 7.6398E-01, 8.2642E-01, 8.6443E-01, 8.8990E-01, 9.2014E-01, 9.3643E-01, 9.5046E-01, 9.5920E-01, 9.7393E-01, 9.9995E-01, 1.0000E+00};
			 double BK_model_high[BK_npT]={9.6956E-01, 9.7603E-01, 9.7966E-01, 9.8224E-01, 9.8509E-01, 9.8861E-01, 9.9186E-01, 9.9197E-01, 9.9316E-01, 9.9338E-01, 1.0000E+00, 1.0000E+00};

			 double BK_model_xVar[BK_npT];
			 double BK_model_errlow[BK_npT];
			 double BK_model_errhigh[BK_npT];

			 for(int BKi=0;BKi<BK_npT; BKi++){
				 BK_model_xVar[BKi]=BK_pT[BKi];
				 if(!functionOfPT) BK_model_xVar[BKi]=BK_pT[BKi]/NRQCDvars::mass[i];
				 BK_model_errlow[BKi]=BK_model[BKi]-BK_model_low[BKi];
				 BK_model_errhigh[BKi]=BK_model_high[BKi]-BK_model[BKi];
			 }

			 TGraphAsymmErrors *BK_model_graph = new TGraphAsymmErrors(BK_npT,BK_model_xVar,BK_model, zero, zero, BK_model_errlow, BK_model_errhigh);

			 cout<<"BK_model_graph:"<<endl;
			 BK_model_graph->Print();

				TGraphAsymmErrors* graphGong[3][2];


				graphGong[1][0] = (TGraphAsymmErrors*)fileGong->Get("lth_HX_rap1_central");
				graphGong[0][0] = (TGraphAsymmErrors*)fileGong->Get("lth_HX_rap1_down");
				graphGong[2][0] = (TGraphAsymmErrors*)fileGong->Get("lth_HX_rap1_up");
				graphGong[1][1] = (TGraphAsymmErrors*)fileGong->Get("lth_HX_rap2_central");
				graphGong[0][1] = (TGraphAsymmErrors*)fileGong->Get("lth_HX_rap2_down");
				graphGong[2][1] = (TGraphAsymmErrors*)fileGong->Get("lth_HX_rap2_up");

				cout<<"lth_HX_rap1_central"<<endl;
				graphGong[1][0]->Print();
				cout<<"lth_HX_rap2_central"<<endl;
				graphGong[1][1]->Print();

				cout<<"lth_HX_rap1_downl"<<endl;
				graphGong[0][0]->Print();
				cout<<"lth_HX_rap2_downl"<<endl;
				graphGong[0][1]->Print();

				cout<<"lth_HX_rap1_up"<<endl;
				graphGong[2][0]->Print();
				cout<<"lth_HX_rap2_up"<<endl;
				graphGong[2][1]->Print();

				const int Gong_npT=27;
				double Gong_model_low[Gong_npT];
				double Gong_model[Gong_npT];
				double Gong_model_high[Gong_npT];

				 double Gong_pT[Gong_npT];
				 double Gong_model_xVar[Gong_npT];
				 double Gong_model_errlow[Gong_npT];
				 double Gong_model_errhigh[Gong_npT];

				 int nGongToRemove=0;
				 for(int Gongi=0;Gongi<Gong_npT; Gongi++){
					 double buffx1, buffx2, buffy1, buffy2;
					 graphGong[1][0]->GetPoint(Gongi, buffx1, buffy1);
					 graphGong[1][1]->GetPoint(Gongi, buffx2, buffy2);
					 Gong_pT[Gongi]=(buffx1+buffx2)/2;
					 Gong_model_xVar[Gongi]=Gong_pT[Gongi];
					 if(Gong_model_xVar[Gongi]<xVarminPred) nGongToRemove++;
					 if(!functionOfPT) Gong_model_xVar[Gongi]=Gong_pT[Gongi]/NRQCDvars::mass[i];

					 Gong_model[Gongi]=(buffy1+buffy2)/2.;

					 graphGong[0][0]->GetPoint(Gongi, buffx1, buffy1);
					 graphGong[0][1]->GetPoint(Gongi, buffx2, buffy2);

					 Gong_model_low[Gongi]=buffy1;
					 if(buffy2<buffy1) Gong_model_low[Gongi]=buffy2;


					 graphGong[2][0]->GetPoint(Gongi, buffx1, buffy1);
					 graphGong[2][1]->GetPoint(Gongi, buffx2, buffy2);

					 Gong_model_high[Gongi]=buffy1;
					 if(buffy2>buffy1) Gong_model_high[Gongi]=buffy2;

					 Gong_model_errlow[Gongi]=Gong_model[Gongi]-Gong_model_low[Gongi];
					 Gong_model_errhigh[Gongi]=Gong_model_high[Gongi]-Gong_model[Gongi];

				 }

				 TGraphAsymmErrors *Gong_model_graph = new TGraphAsymmErrors(Gong_npT,Gong_model_xVar,Gong_model, zero, zero, Gong_model_errlow, Gong_model_errhigh);

				 //for(int iRM=0;iRM<nGongToRemove;iRM++){
				//	 Gong_model_graph->RemovePoint(0);
				 //}

				 cout<<"Gong_model_graph"<<endl;
				 Gong_model_graph->Print();

			double ResMin=-1.3;
			double ResMax=1.3;
			double xVarmin=xVarminPred;
			double xVarmax=xVarmaxPred;
			double fracEmptyPlot=0.04;


			int modelStyle=1;
			int modelWidth=0.;
			int modelColor=1;
			int errmodelStyle=2;
			int errmodelWidth=1;
			int errmodelColor=600;



			LamthPredGraph->GetYaxis()->SetRangeUser(ResMin, ResMax);
			LamthPredGraph->GetXaxis()->SetLimits(xVarmin-(xVarmax-xVarmin)*fracEmptyPlot, xVarmax+(xVarmax-xVarmin)*fracEmptyPlot);
			LamthPredGraph->SetTitle("");
			if(functionOfPT) LamthPredGraph->GetXaxis()->SetTitle("#it{p}_{T}");
			else LamthPredGraph->GetXaxis()->SetTitle("#it{p}_{T}/M");
			LamthPredGraph->GetYaxis()->SetTitleSize(0.05);
			LamthPredGraph->GetYaxis()->SetTitleOffset(1.);
			LamthPredGraph->SetLineStyle(modelStyle);
			LamthPredGraph->SetLineWidth(modelWidth);
			LamthPredGraph->SetLineColor(modelColor);

			LamthPredGraph_lowborder->SetLineStyle(errmodelStyle);
			LamthPredGraph_lowborder->SetLineWidth(errmodelWidth);
			LamthPredGraph_lowborder->SetLineColor(errmodelColor);
			LamthPredGraph_highborder->SetLineStyle(errmodelStyle);
			LamthPredGraph_highborder->SetLineWidth(errmodelWidth);
			LamthPredGraph_highborder->SetLineColor(errmodelColor);

			int FillStyle_sig=1001;
			int Color_1Sig=416;
			int Color_2Sig=400;
			int Color_3Sig=423;
			int Color_ColorBK=632-7;
			int Color_ColorGong=600-6;

			LamthPredGraph1sig->SetFillStyle(FillStyle_sig);
			LamthPredGraph1sig->SetFillColor(Color_1Sig);
			LamthPredGraph1sig->SetLineColor(Color_1Sig);
			LamthPredGraph2sig->SetFillStyle(FillStyle_sig);
			LamthPredGraph2sig->SetFillColor(Color_2Sig);
			LamthPredGraph2sig->SetLineColor(Color_2Sig);
			LamthPredGraph3sig->SetFillStyle(FillStyle_sig);
			LamthPredGraph3sig->SetFillColor(Color_3Sig);
			LamthPredGraph3sig->SetLineColor(Color_3Sig);


			BK_model_graph->SetFillStyle(FillStyle_sig);
			BK_model_graph->SetFillColor(Color_ColorBK);
			BK_model_graph->SetLineColor(Color_ColorBK);

			Gong_model_graph->SetFillStyle(FillStyle_sig);
			Gong_model_graph->SetFillColor(Color_ColorGong);
			Gong_model_graph->SetLineColor(Color_ColorGong);



			TCanvas *plotCanvas = new TCanvas("plotCanvas","plotCanvas",1000,800);
			plotCanvas->SetFillColor(kWhite);
			plotCanvas->SetFrameBorderMode(0);
			plotCanvas->SetRightMargin(0.02);


			int LegendPos=4;
			int LatexPos=4;
			TLegend* legend;
			if(LegendPos==1) legend = new TLegend(0.15,0.7510989,0.35,0.95008242);
			if(LegendPos==2) legend = new TLegend(0.65,0.7510989,0.95,0.95008242);
			if(LegendPos==3) legend = new TLegend(0.6875,0.1310989,0.95,0.48008242);
			if(LegendPos==4) legend = new TLegend(0.15,0.1310989,0.425,0.28008242);
			legend->SetFillColor(0);
			legend->SetTextSize(0.025);
			legend->SetBorderSize(0);

			TLegend* legendData;
			LegendPos=3;
			if(LegendPos==1) legendData = new TLegend(0.15,0.7510989,0.35,0.95008242);
			if(LegendPos==2) legendData = new TLegend(0.65,0.7510989,0.95,0.95008242);
			if(LegendPos==3) legendData = new TLegend(0.6375,0.1310989,0.95,0.42008242);
			if(LegendPos==3 && i==3) legendData = new TLegend(0.6375,0.1310989,0.95,0.36008242);
			if(LegendPos==4) legendData = new TLegend(0.15,0.1310989,0.35,0.33008242);
			legendData->SetFillColor(0);
			legendData->SetFillStyle(0);
			legendData->SetTextSize(0.025);
			legendData->SetBorderSize(0);

			TLegend* legendBK;
			LegendPos=3;
			legendBK = new TLegend(0.6375,0.85, 0.95,0.89);
			legendBK->SetFillColor(0);
			legendBK->SetTextSize(0.025);
			legendBK->SetBorderSize(0);

			TLegend* legendGong;
			LegendPos=3;
			legendGong = new TLegend(0.6375,0.88, 0.95,0.92);
			legendGong->SetFillColor(0);
			legendGong->SetTextSize(0.025);
			legendGong->SetBorderSize(0);

			char DrawChar[200];
			char legendentry[200];



			sprintf(DrawChar,"al");
			LamthPredGraph->Draw(DrawChar);
			LamthPredGraph3sig->Draw("3same");
			LamthPredGraph2sig->Draw("3same");
			LamthPredGraph1sig->Draw("3same");


			if(plotBK) BK_model_graph->Draw("3same");
			if(plotGong) Gong_model_graph->Draw("3same");

			//LamthPredGraph->Draw(DrawChar);

			//LamthPredGraph_lowborder->Draw(DrawChar);
			//LamthPredGraph_highborder->Draw(DrawChar);

			//sprintf(legendentry,"%1.1f < |y| < %1.1f",rapminData[iRap], rapmaxData[iRap]);
			//sprintf(legendentry,"Central value");
			//legend->AddEntry(LamthPredGraph,legendentry,"lp");
			sprintf(legendentry,"68.3%% CL model");
			legend->AddEntry(LamthPredGraph1sig,legendentry,"f");
			sprintf(legendentry,"95.5%% CL model");
			legend->AddEntry(LamthPredGraph2sig,legendentry,"f");
			sprintf(legendentry,"99.7%% CL model");
			legend->AddEntry(LamthPredGraph3sig,legendentry,"f");

			sprintf(legendentry,"BK, PRL 108 (2012) 172002");
			if(plotBK) legendBK->AddEntry(BK_model_graph,legendentry,"f");
			sprintf(legendentry,"GWWZ, arXiv:1305.0748");
			if(plotGong) legendGong->AddEntry(Gong_model_graph,legendentry,"f");


			int iStatebuff;
			sprintf(DrawChar,"psame");
			char rapchar[200];

			if(useOnlyState==3){
			for(int iRap=0;iRap<2;iRap++){
				if(iRap==0) sprintf(rapchar,"|#it{y}| < 0.6");
				if(iRap==1) sprintf(rapchar,"0.6 < |#it{y}| < 1.2");
					iStatebuff=0;
					if(plotAllData) graphData[iRap][iStatebuff]->Draw(DrawChar);
					sprintf(legendentry,"%s, CMS %s", NRQCDvars::StateNameTex[iStatebuff], rapchar);
					legendData->AddEntry(graphData[iRap][iStatebuff],legendentry,"lp");
			}//rapData
			for(int iRap=0;iRap<2;iRap++){
				if(iRap==0) sprintf(rapchar,"|#it{y}| < 0.6");
				if(iRap==1) sprintf(rapchar,"0.6 < |#it{y}| < 1.2");
					iStatebuff=3;
					graphData[iRap][iStatebuff]->Draw(DrawChar);
					sprintf(legendentry,"%s, CMS %s", NRQCDvars::StateNameTex[iStatebuff], rapchar);
					legendData->AddEntry(graphData[iRap][iStatebuff],legendentry,"lp");
			}//rapData
			}


			if(useOnlyState==10){
			for(int iRap=0;iRap<2;iRap++){
			if(iRap==0) sprintf(rapchar,"|#it{y}| < 0.6");
			if(iRap==1) sprintf(rapchar,"0.6 < |#it{y}| < 1.2");
					iStatebuff=4;
					if(plotAllData) graphData[iRap][iStatebuff]->Draw(DrawChar);
					sprintf(legendentry,"%s, CMS %s", NRQCDvars::StateNameTex[iStatebuff], rapchar);
					legendData->AddEntry(graphData[iRap][iStatebuff],legendentry,"lp");
			}//rapData
			for(int iRap=0;iRap<2;iRap++){
			if(iRap==0) sprintf(rapchar,"|#it{y}| < 0.6");
			if(iRap==1) sprintf(rapchar,"0.6 < |#it{y}| < 1.2");
					iStatebuff=7;
					if(plotAllData) graphData[iRap][iStatebuff]->Draw(DrawChar);
					sprintf(legendentry,"%s, CMS %s", NRQCDvars::StateNameTex[iStatebuff], rapchar);
					legendData->AddEntry(graphData[iRap][iStatebuff],legendentry,"lp");
			}//rapData
			for(int iRap=0;iRap<2;iRap++){
					if(iRap==0) sprintf(rapchar,"|#it{y}| < 0.6");
					if(iRap==1) sprintf(rapchar,"0.6 < |#it{y}| < 1.2");
					iStatebuff=10;
					graphData[iRap][iStatebuff]->Draw(DrawChar);
					sprintf(legendentry,"%s, CMS %s", NRQCDvars::StateNameTex[iStatebuff], rapchar);
					legendData->AddEntry(graphData[iRap][iStatebuff],legendentry,"lp");
			}//rapData
			}


			TLine* extreme0 = new TLine( xVarmin, 0, xVarmax, 0);
			extreme0->SetLineWidth( 1 );
			extreme0->SetLineStyle( 2 );
			extreme0->SetLineColor( kBlack );
			extreme0->Draw( "same" );

			legend->Draw("same");
			legendData->Draw("same");
			if(plotBK) legendBK->Draw("same");
			if(plotGong) legendGong->Draw("same");

			//draw latex
			double left=0.175, top=0.92, textSize=0.045;
			TLatex *latex=new TLatex();
			//latex->SetTextFont(42);
			latex->SetNDC(kTRUE);
			latex->SetTextSize(textSize);
			double stepLatex=textSize*1.3;
			//latex->SetTextFont(62);

			if(LatexPos==1 || LatexPos==4) left=0.175;
			if(LatexPos==2 || LatexPos==3) left=0.78;

			if(LatexPos==1 || LatexPos==2) top=0.90;
			if(LatexPos==3 || LatexPos==4) top=0.27;

			textSize=0.06;
			latex->SetTextSize(textSize);
			//latex->DrawLatex(left,top, Form("%s",StateNameTex[i]));

			left=0.16;
			if(LatexPos==1 || LatexPos==2) top=0.84;
			if(LatexPos==3 || LatexPos==4) top=0.30;
			textSize=0.0285;
			latex->SetTextSize(textSize);
			latex->DrawLatex(left,top, Form("|#it{y}| < 1.2"));


			left=0.0075; top=0.54;
			textSize=0.07;
			latex->SetTextFont(42);
			latex->SetTextSize(textSize);
			latex->DrawLatex(left,top, "#lambda_{#vartheta}^{HX}");


			char saveName[500];
			sprintf(saveName,"Figures/PlotPPDderivative/%s",JobID);
			//if(plotAllData && !functionOfPT) sprintf(saveName,"%s/Lamth_Pred_state%d_withAllData_pTovM.pdf",saveName, i);
			//if(!plotAllData && !functionOfPT) sprintf(saveName,"%s/Lamth_Pred_state%d_withData_pTofM.pdf",saveName, i);


			if(!plotModel && !plotAllData) sprintf(saveName,"%s/Lamth_Pred_state%d_withFittedData_pT.pdf",saveName, i);
			if(!plotModel && plotAllData) sprintf(saveName,"%s/Lamth_Pred_state%d_withAllData_pT.pdf",saveName, i);
			if(plotModel && plotAllData) sprintf(saveName,"%s/Lamth_Pred_state%d_withAllDataAndModel_pT.pdf",saveName, i);
			cout<<saveName<<endl;
			plotCanvas->SaveAs(saveName);
			if(!plotModel && !plotAllData) sprintf(saveName,"%s/Lamth_Pred_state%d_withFittedData_pT.png",saveName, i);
			if(!plotModel && plotAllData) sprintf(saveName,"%s/Lamth_Pred_state%d_withAllData_pT.png",saveName, i);
			if(plotModel && plotAllData) sprintf(saveName,"%s/Lamth_Pred_state%d_withAllDataAndModel_pT.png",saveName, i);
			cout<<saveName<<endl;
			plotCanvas->SaveAs(saveName);









			DummyFile->Close();
			delete DummyFile;



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


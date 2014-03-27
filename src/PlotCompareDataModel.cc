/*
 * PlotCompareDataModel.cc
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
#include "TBox.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TLatex.h"
#include "TLine.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TGraphAsymmErrors.h"
#include "TGraph.h"
#include "TROOT.h"


using namespace NRQCDvars;


void plotComp(int iState, int iMeasurementID, int iExperiment, int iRap, char jobdirname[200], char rapchar[200], TGraphAsymmErrors *data_Graph, TGraph *model_Graph, TGraph *model_Graph_low, TGraph *model_Graph_high, bool plotInclusiveFeedDown, bool plotIndividualFeedDown, bool plotDirectColorChannels, const int StatesCont_c, const int ColorChannels_c, TGraph *model_Graph_FeedDowns[500], TGraph *model_Graph_ColorChannels[500], /*TGraphAsymmErrors *model_Graph_ColorChannels[3][500],*/ TGraph *data_Graph_nopolcorr, vector<int> StatesContributing, double chi2Min, double chi2Prob, int ndf, bool HPbool, double pTMinModel, double pTMaxModel, bool longrapchar, TGraphAsymmErrors *model_Graph_ColorChannels_Bands[3][500], TGraph *model_Graph_low_Bands[3], TGraph *model_Graph_high_Bands[3], TGraphAsymmErrors *model_Graph_Bands[3], TGraph *model_CS_THunc[3], TGraphAsymmErrors *model_CS_THuncBandbool, bool plotDirectInclusive, int nDataPointsInRange, bool PredictionPlot, bool PredictionDataPlot, bool smoothPol, bool PlotVSpToverM);
vector<double> addPolarizations(vector<vector<double> > lamMatrix, vector<double> lamVecContributions);
void FindMPV(TH1* PosteriorDist , double& MPV , double& MPVerrorLow, double& MPVerrorHigh, int MPValgo, int nSigma);
TGraph* smoothPolFunc(TGraph *graphToSmooth, int monotonicInt, int c_sign);
TGraphAsymmErrors* CloneCentralValueFromGraph(TGraphAsymmErrors *graphToAdapt, TGraph *graphToAdaptFrom);
TGraphAsymmErrors* smoothPolFuncErrors(TGraphAsymmErrors *graphToAdapt, int polOrder);


Double_t paramParabola(Double_t *x, Double_t *par);

int main(int argc, char** argv) {


	bool HPbool=true;


  	Char_t *JobID = "Default";
  	Char_t *ModelID = "Default";
  	Char_t *DataModelCombID = "Default";
	Char_t *storagedir = "Default"; //Storage Directory

  	double 	pTMin=0;
  	double 	pTMax=1000.;
  	double 	rapMin=-10.;
  	double 	rapMax=10.;
  	bool 	useSstatesOnly=false;
  	bool 	usePstatesOnly=false;
  	bool 	useCharmoniumOnly=false;
  	bool 	useBottomoniumOnly=false;
  	bool 	useCrossSectionOnly=false;
  	bool 	usePolarizationOnly=false;
  	int 	useOnlyState=999;
  	bool	PredictionPlot=false;
  	bool	PredictionDataPlot=false;
	bool 	PlotVSpToverM=false;

  	int 	MPValgo=-1;
  	int 	nSigma=-1;


  	for( int i=0;i < argc; ++i ) {
		if(std::string(argv[i]).find("storagedir") != std::string::npos) {char* storagedirchar = argv[i]; char* storagedirchar2 = strtok (storagedirchar, "="); storagedir = storagedirchar2; cout<<"storagedir = "<<storagedir<<endl;}
	    if(std::string(argv[i]).find("nSigma") != std::string::npos) { char* nSigmachar = argv[i]; char* nSigmachar2 = strtok (nSigmachar, "n"); nSigma = atoi(nSigmachar2); cout<<"nSigma = "<<nSigma<<endl; }
	    if(std::string(argv[i]).find("MPValgo") != std::string::npos) { char* MPValgochar = argv[i]; char* MPValgochar2 = strtok (MPValgochar, "M"); MPValgo = atoi(MPValgochar2); cout<<"MPValgo = "<<MPValgo<<endl; }
		if(std::string(argv[i]).find("JobID") != std::string::npos) { char* JobIDchar = argv[i]; char* JobIDchar2 = strtok (JobIDchar, "="); JobID = JobIDchar2; cout<< "JobID = " << JobID << endl; }
		if(std::string(argv[i]).find("ModelID") != std::string::npos) {char* ModelIDchar = argv[i]; char* ModelIDchar2 = strtok (ModelIDchar, "="); ModelID = ModelIDchar2; cout<<"ModelID = "<<ModelID<<endl;}
		if(std::string(argv[i]).find("DataModelCombID") != std::string::npos) {char* DataModelCombIDchar = argv[i]; char* DataModelCombIDchar2 = strtok (DataModelCombIDchar, "="); DataModelCombID = DataModelCombIDchar2; cout<<"DataModelCombID = "<<DataModelCombID<<endl;}
	    if(std::string(argv[i]).find("pTMin") != std::string::npos) {
	    	char* pTMinchar = argv[i];
	    	char* pTMinchar2 = strtok (pTMinchar, "p");
	    	pTMin = atof(pTMinchar2);
	    	cout<<"pTMin = "<<pTMin<<endl;
	    }
	    if(std::string(argv[i]).find("pTMax") != std::string::npos) {
	    	char* pTMaxchar = argv[i];
	    	char* pTMaxchar2 = strtok (pTMaxchar, "p");
	    	pTMax = atof(pTMaxchar2);
	    	cout<<"pTMax = "<<pTMax<<endl;
	    }
	    if(std::string(argv[i]).find("rapMin") != std::string::npos) {
	    	char* rapMinchar = argv[i];
	    	char* rapMinchar2 = strtok (rapMinchar, "r");
	    	rapMin = atof(rapMinchar2);
	    	cout<<"rapMin = "<<rapMin<<endl;
	    }
	    if(std::string(argv[i]).find("rapMax") != std::string::npos) {
	    	char* rapMaxchar = argv[i];
	    	char* rapMaxchar2 = strtok (rapMaxchar, "r");
	    	rapMax = atof(rapMaxchar2);
	    	cout<<"rapMax = "<<rapMax<<endl;
	    }
	    if(std::string(argv[i]).find("useOnlyState") != std::string::npos) {
	    	char* useOnlyStatechar = argv[i];
	    	char* useOnlyStatechar2 = strtok (useOnlyStatechar, "u");
	    	useOnlyState = atoi(useOnlyStatechar2);
	    	cout<<"useOnlyState = "<<useOnlyState<<endl;
	    }
	    if(std::string(argv[i]).find("useSstatesOnly=true") != std::string::npos) {
	    	useSstatesOnly=true;
	    	cout<<"useSstatesOnly=true"<<endl;
	    }
	    if(std::string(argv[i]).find("usePstatesOnly=true") != std::string::npos) {
	    	usePstatesOnly=true;
	    	cout<<"usePstatesOnly=true"<<endl;
	    }
	    if(std::string(argv[i]).find("useCharmoniumOnly=true") != std::string::npos) {
	    	useCharmoniumOnly=true;
	    	cout<<"useCharmoniumOnly=true"<<endl;
	    }
	    if(std::string(argv[i]).find("useBottomoniumOnly=true") != std::string::npos) {
	    	useBottomoniumOnly=true;
	    	cout<<"useBottomoniumOnly=true"<<endl;
	    }
	    if(std::string(argv[i]).find("useCrossSectionOnly=true") != std::string::npos) {
	    	useCrossSectionOnly=true;
	    	cout<<"useCrossSectionOnly=true"<<endl;
	    }
	    if(std::string(argv[i]).find("usePolarizationOnly=true") != std::string::npos) {
	    	usePolarizationOnly=true;
	    	cout<<"usePolarizationOnly=true"<<endl;
	    }
	    if(std::string(argv[i]).find("PredictionPlot=true") != std::string::npos) {
	    	PredictionPlot=true;
	    	cout<<"PredictionPlot=true"<<endl;
	    }
	    if(std::string(argv[i]).find("PredictionDataPlot=true") != std::string::npos) {
	    	PredictionDataPlot=true;
	    	cout<<"PredictionDataPlot=true"<<endl;
	    }
	    if(std::string(argv[i]).find("PlotVSpToverM=true") != std::string::npos) {
	    	PlotVSpToverM=true;
	    	cout<<"PlotVSpToverM=true"<<endl;
	    }

  	}


	char inname[200];
	char outname[200];
	char jobdirname[200];
	char modeldirname[200];
	char datamodeldirname[200];
	sprintf(jobdirname,"%s/JobID", storagedir);
	gSystem->mkdir(jobdirname);
	sprintf(jobdirname,"%s/%s",jobdirname,JobID);
	gSystem->mkdir(jobdirname);
	sprintf(modeldirname,"%s/ModelID", storagedir);
	gSystem->mkdir(modeldirname);
	sprintf(modeldirname,"%s/%s",modeldirname,ModelID);
	gSystem->mkdir(modeldirname);
	sprintf(datamodeldirname,"%s/DataModelCombID", storagedir);
	gSystem->mkdir(datamodeldirname);
	sprintf(datamodeldirname,"%s/%s",datamodeldirname,DataModelCombID);
	gSystem->mkdir(datamodeldirname);

	gSystem->mkdir(Form("%s/Figures",jobdirname));

	sprintf(inname,"%s/results.txt",jobdirname);
	cout<<"read results from "<<inname<<endl;

  	dmatrix fi_MPV;
  	dmatrix fi_errlow;
 	dmatrix fi_errhigh;
  	dmatrix Oi_MPV;
  	dmatrix Oi_errlow;
 	dmatrix Oi_errhigh;

    ifstream in;
    in.open(inname);//, std::ofstream::app);

	in >> fi_MPV;
	in >> fi_errlow;
	in >> fi_errhigh;
	in >> Oi_MPV;
	in >> Oi_errlow;
	in >> Oi_errhigh;

    in.close();

    cout<<"fi_MPV:"<<endl;
	cout<<fi_MPV<<endl;
    cout<<"fi_errlow:"<<endl;
	cout<<fi_errlow<<endl;
    cout<<"fi_errhigh:"<<endl;
	cout<<fi_errhigh<<endl;
    cout<<"Oi_MPV:"<<endl;
	cout<<Oi_MPV<<endl;
    cout<<"Oi_errlow:"<<endl;
	cout<<Oi_errlow<<endl;
    cout<<"Oi_errhigh:"<<endl;
	cout<<Oi_errhigh<<endl;

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
    bool MinimizerMH=true;
    bool MinimizerMinuit=false;
    if(int_Minimizer[0]==1) {MinimizerMH=false; MinimizerMinuit=true;}



  	dmatrix Np_BR_MPV;
  	dmatrix Np_BR_errlow;
 	dmatrix Np_BR_errhigh;
  	dmatrix Np_US_MPV_ill;
  	dmatrix Np_US_errlow_ill;
 	dmatrix Np_US_errhigh_ill;
  	dmatrix consts_star_var_MPV;
  	dmatrix consts_star_var_errlow;
 	dmatrix consts_star_var_errhigh;

	sprintf(inname,"%s/results_Np.txt",jobdirname);
	cout<<"get results_Np from "<<inname<<endl;

    in.open(inname);//, std::ofstream::app);

    cout<<"Np_BR:"<<endl;
	in >> Np_BR_MPV;
    cout<<Np_BR_MPV<<endl;
    cout<<"errlow_Np_BR:"<<endl;
	in >> Np_BR_errlow;
    cout<<Np_BR_errlow<<endl;
    cout<<"errlow_Np_BR:"<<endl;
	in >> Np_BR_errhigh;
    cout<<Np_BR_errhigh<<endl;
    cout<<"consts_star_var:"<<endl;
	in >> consts_star_var_MPV;
    cout<<consts_star_var_MPV<<endl;
    cout<<"errlow_consts_star_var:"<<endl;
	in >> consts_star_var_errlow;
    cout<<consts_star_var_errlow<<endl;
    cout<<"errhigh_consts_star_var:"<<endl;
	in >> consts_star_var_errhigh;
    cout<<consts_star_var_errhigh<<endl;
    cout<<"Np_US:"<<endl;
	in >> Np_US_MPV_ill;
    cout<<Np_US_MPV_ill<<endl;
    cout<<"errlow_Np_US:"<<endl;
	in >> Np_US_errlow_ill;
    cout<<Np_US_errlow_ill<<endl;
    cout<<"errhigh_Np_US:"<<endl;
	in >> Np_US_errhigh_ill;
    cout<<Np_US_errhigh_ill<<endl;

    in.close();


	sprintf(inname,"%s/chi2ndf.txt",jobdirname);
	cout<<"read chi2 from "<<inname<<endl;

  	double chi2Min;
  	double reduced_chi2Min;
  	int nDataPoints;
  	int nOps;
  	int ndf;

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

    double chi2Prob;
    chi2Prob=TMath::Prob(chi2Min, ndf);


	dmatrix Np_US_MPV; //Uncertainty scales, [0=Data, 1=Model][nScales]
	dmatrix Np_US_errlow; //Uncertainty scales, [0=Data, 1=Model][nScales]
	dmatrix Np_US_errhigh; //Uncertainty scales, [0=Data, 1=Model][nScales]

	dvector Np_US_0 (NRQCDvars::nDataSystematicScales, 0.);
	dvector Np_US_1 (NRQCDvars::nModelSystematicScales, 0.);
	Np_US_MPV.push_back(Np_US_0);
	Np_US_MPV.push_back(Np_US_1);
	Np_US_errlow.push_back(Np_US_0);
	Np_US_errlow.push_back(Np_US_1);
	Np_US_errhigh.push_back(Np_US_0);
	Np_US_errhigh.push_back(Np_US_1);

	for(int j=0; j < NRQCDvars::nDataSystematicScales; j++){
		Np_US_MPV[0][j]=Np_US_MPV_ill[0][j];
		Np_US_errlow[0][j]=Np_US_errlow_ill[0][j];
		Np_US_errhigh[0][j]=Np_US_errhigh_ill[0][j];
	}
	for(int j=0; j < NRQCDvars::nModelSystematicScales; j++){
		Np_US_MPV[1][j]=Np_US_MPV_ill[1][j];
		Np_US_errlow[1][j]=Np_US_errlow_ill[1][j];
		Np_US_errhigh[1][j]=Np_US_errhigh_ill[1][j];
	}


    cout<<"Np_US_MPV[0].size() "<<Np_US_MPV[0].size() <<endl;
    cout<<"Np_US_MPV[1].size() "<<Np_US_MPV[1].size() <<endl;





    cout<<"Inizialize Op_plus:"<<endl;

	dmatrix Op_plus(NRQCDvars::nStates);
	vector<double> Op_S_plus (NRQCDvars::nColorChannels_S,0);
	vector<double> Op_P_plus (NRQCDvars::nColorChannels_P,0);

	for (int i=0; i<NRQCDvars::nStates; i++){
		bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
		if(isSstate){
			for (int j=0; j<NRQCDvars::nColorChannels_S; j++){
				Op_S_plus.at(j)=Oi_MPV[i][j]+Oi_errhigh[i][j];
			}
			Op_plus.at(i)=Op_S_plus;
		}
		else{
			for (int j=0; j<NRQCDvars::nColorChannels_P; j++){
				Op_P_plus.at(j)=Oi_MPV[i][j]+Oi_errhigh[i][j];
			}
			Op_plus.at(i)=Op_P_plus;
		}
	}

    cout<<"Inizialize Op_minus:"<<endl;

    dmatrix Op_minus(NRQCDvars::nStates);
	vector<double> Op_S_minus (NRQCDvars::nColorChannels_S,0);
	vector<double> Op_P_minus (NRQCDvars::nColorChannels_P,0);

	for (int i=0; i<NRQCDvars::nStates; i++){
		bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
		if(isSstate){
			for (int j=0; j<NRQCDvars::nColorChannels_S; j++){
				Op_S_minus.at(j)=Oi_MPV[i][j]-Oi_errlow[i][j];
			}
			Op_minus.at(i)=Op_S_minus;
		}
		else{
			for (int j=0; j<NRQCDvars::nColorChannels_P; j++){
				Op_P_minus.at(j)=Oi_MPV[i][j]-Oi_errlow[i][j];
			}
			Op_minus.at(i)=Op_P_minus;
		}
	}

    cout<<"Op_plus:"<<endl;
	cout<<Op_plus<<endl;
    cout<<"Op_minus:"<<endl;
	cout<<Op_minus<<endl;

    //cout<<"Inizialize Np's:"<<endl;
    //
    //dmatrix Np_BR; //Branching ratios, [nDauthgers][nMothers]
	//dmatrix Np_US; //Uncertainty scales, [0=Data, 1=Model][nScales]
    //
	//dvector Np_BR_0 (NRQCDvars::nStates,0.);
	//for(int i=0; i < NRQCDvars::nStates; i++) Np_BR.push_back(Np_BR_0);
    //
	//dvector Np_US_0 (NRQCDvars::nDataSystematicScales, 0.);
	//dvector Np_US_1 (NRQCDvars::nModelSystematicScales, 0.);
	//Np_US.push_back(Np_US_0);
	//Np_US.push_back(Np_US_1);
    //
	//for(int j=0; j < NRQCDvars::nDataSystematicScales; j++){
	//	Np_US[0][j]=0.;
	//}
	//for(int j=0; j < NRQCDvars::nModelSystematicScales; j++){
	//	Np_US[1][j]=0.;
	//}
    //
	//for(int i=0; i < nStates; i++){
	//	for(int j=0; j < nStates; j++){
	//		if(NRQCDvars::FeedDownBranchingRatio[i][j]>0){
	//			Np_BR[i][j]=NRQCDvars::FeedDownBranchingRatio[i][j]*0.01;
	//		}
	//		else Np_BR[i][j]=0;
	//	}
	//}

	double pTMinModel=pTMin;
	double pTMaxModel=pTMax;


    cout<<"Loop through data to plot:"<<endl;

	bool DataSelected=false;
	bool warnInconsistency[nStates][NRQCDvars::nMeasurementIDs][NRQCDvars::nExperiments][NRQCDvars::nMaxRapBins][NRQCDvars::nMaxPtBins];

	for(int iState=0; iState < nStates; iState++){
		for(int iMeasurementID=0; iMeasurementID < NRQCDvars::nMeasurementIDs; iMeasurementID++){
			for(int iExperiment=0; iExperiment < NRQCDvars::nExperiments; iExperiment++){

		    	//cout<<"Plot iState="<<StateName[iState]<<", iMeasurementID="<<MeasurementIDName[iMeasurementID]<<", iExperiment="<<ExpName[iExperiment]<<endl;

				NRQCDglobalfitObject *DataModelObject[NRQCDvars::nMaxRapBins][NRQCDvars::nMaxPtBins];
				bool DataPresentAndSelected[NRQCDvars::nMaxRapBins][NRQCDvars::nMaxPtBins];

				bool definedRap=false;
				bool isAbsRap=false;
				double rapMinObject, rapMaxObject;

				int StatesCont;//Number of states contributing to result, including direct state

				for(int iRap = 0; iRap < NRQCDvars::nMaxRapBins; iRap++){
			    	//cout<<"plot model for iState="<<StateName[iState]<<", iMeasurementID="<<MeasurementIDName[iMeasurementID]<<", iExperiment="<<ExpName[iExperiment]<<", iRap="<<iRap/*", iP="<<iP*/<<endl;

					definedRap=false;
					int nPtBinsSel=0;
				    for(int iP = 0; iP < NRQCDvars::nMaxPtBins; iP++){

				    	warnInconsistency[iState][iMeasurementID][iExperiment][iRap][iP]=false;

				    	DataPresentAndSelected[iRap][iP]=false;

				    	//if(iState!=0) continue;
				    	//if(iExperiment!=0) continue;
				    	//if(iMeasurementID!=0) continue;
				    	//if(iRap!=0) continue;
				    	//if(iP>0) continue;
				    	//if(iState<1) continue;



				    	if(iState==3) {pTMin=12; pTMinModel=12;}//CHANGE_BACK
				    	if(iState==10) {pTMin=30; pTMinModel=30;}//CHANGE_BACK
				    	if(iState!=3 && iState!=10) continue;//CHANGE_BACK

				    	if(iState==3 && iExperiment==3) {pTMin=8; pTMinModel=8;}//CHANGE_BACK
				    	if(iState==10 && iMeasurementID==1) {pTMin=16; pTMinModel=16;}//CHANGE_BACK



				    	if(PredictionPlot){
					    	if(iState==3) {pTMin=9; pTMinModel=9; pTMax=10000; pTMaxModel=10000;}//CHANGE_BACK
					    	if(iState==10) {pTMin=22.5; pTMinModel=22.5; pTMax=20000; pTMaxModel=20000;}//CHANGE_BACK
					    	if(iState!=3 && iState!=10) continue;//CHANGE_BACK
					    	if(PredictionDataPlot){
						    	if(iState==3) {pTMin=12; pTMinModel=12; pTMax=10000; pTMaxModel=10000;}//CHANGE_BACK
						    	if(iState==10) {pTMin=30; pTMinModel=30; pTMax=20000; pTMaxModel=20000;}//CHANGE_BACK
						    	rapMax=1.2;
					    	}

				    	}




				    	//if(iState!=0 || iMeasurementID!=0 || iExperiment!=0 || iRap!=0 || iP!=0) continue;

						sprintf(inname,"%s/ConvertedDataModel_%s_%s_%s_rap%d_pT%d.txt",datamodeldirname, StateName[iState],
								MeasurementIDName[iMeasurementID],  ExpName[iExperiment], iRap, iP);

						ifstream inData;
						inData.open(inname);

						//cout<<inname<<endl;

						if( inData.is_open() ){//Measurement present -> calculate model components:: Modified by Joao: more correct from ios point of view
							//if(NRQCDvars::debug){
								//cout << "Read in iState=" << StateName[iState] << ", iMeasurementID=" << MeasurementIDName[iMeasurementID];
								//cout << ", iExperiment=" << ExpName[iExperiment] << ", iRap=" << iRap << ", iP=" << iP << endl;
							//}
							DataModelObject[iRap][iP] = new NRQCDglobalfitObject();
							inData >> *DataModelObject[iRap][iP];

							if(NRQCDvars::debug) DataModelObject[iRap][iP]->Dump(nStates, true, true);
							inData.close();

							if(
									DataModelObject[iRap][iP]->getpTMin() >= pTMin
								&& DataModelObject[iRap][iP]->getpTMax() <= pTMax
								&& DataModelObject[iRap][iP]->getyMin() >= rapMin
								&& DataModelObject[iRap][iP]->getyMax() <= rapMax
							) DataSelected=true;

							if(
								HPbool
								&& DataModelObject[iRap][iP]->getyMin() >= rapMin
								&& DataModelObject[iRap][iP]->getyMax() <= rapMax
							) DataSelected=true;
							if(
								PredictionDataPlot &&
							(DataModelObject[iRap][iP]->getpTMin() < pTMin
							|| DataModelObject[iRap][iP]->getpTMax() > pTMax)
							) DataSelected=false;

							//if(DataModelObject[iRap][iP]->getpTMin() >= 29 && iMeasurementID==1 && pTMin==35) DataSelected=true;

							if(useSstatesOnly && StateQuantumID[DataModelObject[iRap][iP]->getState()] != NRQCDvars::quID_S)
								DataSelected=false;
							if(usePstatesOnly && StateQuantumID[DataModelObject[iRap][iP]->getState()] == NRQCDvars::quID_S)
								DataSelected=false;
							if(useCharmoniumOnly && DataModelObject[iRap][iP]->getState() > 3)
								DataSelected=false;
							if(useBottomoniumOnly && DataModelObject[iRap][iP]->getState() < 4)
								DataSelected=false;
							if(useOnlyState < nStates &&  DataModelObject[iRap][iP]->getState() != useOnlyState)
								DataSelected=false;


							if(DataSelected){
								//cout<<"Data selected"<<endl;
								DataPresentAndSelected[iRap][iP]=true;
								nPtBinsSel++;
								if(!definedRap){
									isAbsRap=DataModelObject[iRap][iP]->getisAbsRap();
									rapMinObject=DataModelObject[iRap][iP]->getyMin();
									rapMaxObject=DataModelObject[iRap][iP]->getyMax();
									definedRap=true;

								}
							}
							DataSelected=false;


						}
				    }//pT loop

				    //cout<<"nPtBinsSel "<<nPtBinsSel<<endl;

			    	if(PlotVSpToverM){
			    		pTMinModel/=NRQCDvars::mass[iState];
			    		pTMaxModel/=NRQCDvars::mass[iState];
			    		pTMin/=NRQCDvars::mass[iState];
			    		pTMax/=NRQCDvars::mass[iState];
			    	}

				    dvector data_centralval;
				    dvector data_errlow_centralval;
				    dvector data_errhigh_centralval;
				    dvector data_pTmean;
				    dvector data_errlow_pT;
				    dvector data_errhigh_pT;

				    dvector model_centralval;
				    dvector model_errlow_centralval;
				    dvector model_errhigh_centralval;
				    dvector model_errlow_centralval_1Sig;
				    dvector model_errhigh_centralval_1Sig;
				    dvector model_errlow_centralval_2Sig;
				    dvector model_errhigh_centralval_2Sig;
				    dvector model_errlow_centralval_3Sig;
				    dvector model_errhigh_centralval_3Sig;

				    dvector model_inclusive_FeedDown;

				    dmatrix model_directProduction; //(nColorChannels, iP)
				    dmatrix model_FeedDown; //(StatesCont, iP)

				    dmatrix model_directProduction_Bands;
				    dmatrix model_directProduction_Bands_errlow_1sig;
				    dmatrix model_directProduction_Bands_errhigh_1sig;
				    dmatrix model_directProduction_Bands_errlow_2sig;
				    dmatrix model_directProduction_Bands_errhigh_2sig;
				    dmatrix model_directProduction_Bands_errlow_3sig;
				    dmatrix model_directProduction_Bands_errhigh_3sig;

				    dmatrix model_directProduction_THunc_CS;

					//PAPER: implement the restriction of the PPD to SDC*LDME > 0... (difficult for P, as it changes sign -> what to do? -> Pietro mail)

				    //TODO: code unbinned models
				    //TODO: code submodels -> color channels and feed-down

				    // -> getDirectProduction: add to dvector(4) return dvector CS(nColorC) + dmatrix of lambda(nColorC, nLambda)
				    // -> getPromptProduction: add to dvector(4) return dvector CS(nStates) + dmatrix of lambda(nStates, nLambda)
				    // -> getPromptCorrProduction: add return polCorrFactor, to be used to correct the data -> show data-points polarization-corrected

				    // dcube (nStates, nColorChannels, 4) directProductionCube
				    // dmatrix (nStates, 4) promptProductionCube
				    // double polCorrFactor

				    //Plot: all feed-downs (CC inclusive) individually + inclusive feed-down + all CCs from direct state

		    		vector<int> StatesContributing_;
		    		vector<double>  polCorrFactorVec;
		    		vector<double>  lumiCorrFactorVec;
		    		int globalStatesCont;
		    		int nDataPointsInRange=0;


				    for(int iP = 0; iP < NRQCDvars::nMaxPtBins; iP++){

				    	if(DataPresentAndSelected[iRap][iP]){

				    		cout<<"Plot state: "<<NRQCDvars::StateName[iState]<<" Measurement: "<<NRQCDvars::MeasurementIDName[iMeasurementID]<<" Experiment: "<<NRQCDvars::ExpName[iExperiment]<<" iRap: "<<iRap<<" iP "<<iP<<endl;

							dvector ObjectLikelihoodVec;
							double ModelPrediction;
							double errModelPrediction;
							double polCorrFactor;
							double lumiCorrFactor;

							lumiCorrFactor=1.+DataModelObject[iRap][iP]->getErrGlobal()/DataModelObject[iRap][iP]->getCentralValue()*Np_US_MPV[0][iExperiment];
							lumiCorrFactorVec.push_back(lumiCorrFactor);

							dcube directProductionCube;
							dmatrix promptProductionMatrix;


							bool LoopThroughPPD=false;
							if(MinimizerMH) LoopThroughPPD=true;

							//if(iMeasurementID==0) LoopThroughPPD=false;
							//LoopThroughPPD=false;
							//if(iMeasurementID==1)  LoopThroughPPD=true;
							//if(iExperiment>5 && iState==10)  LoopThroughPPD=true;
							//if(iMeasurementID==1 && iState==10)  LoopThroughPPD=true;

							if(
									!PlotVSpToverM && (DataModelObject[iRap][iP]->getpTMin()<pTMinModel || DataModelObject[iRap][iP]->getpTMax()>pTMaxModel)
									|| PlotVSpToverM && (DataModelObject[iRap][iP]->getpTMin()<pTMinModel*NRQCDvars::mass[iState] || DataModelObject[iRap][iP]->getpTMax()>pTMaxModel*NRQCDvars::mass[iState])
									) LoopThroughPPD=false;
							else nDataPointsInRange++;

							//Find out number of contributing feed-down states:
							dcube global_directProductionCube;
							dmatrix global_promptProductionMatrix;
							double global_polCorrFactor;

							ObjectLikelihoodVec=DataModelObject[iRap][iP]->getObjectLikelihood(Oi_MPV, Np_BR_MPV, Np_US_MPV, true, global_directProductionCube, global_promptProductionMatrix, global_polCorrFactor);

				    		dvector globalFeedDownContVec;
				    		dvector globalFeedDownLamVec(3,0);
				    		dmatrix globalFeedDownLamMatrix;
				    		vector<int> globalStatesContributing;
				    		double globalBufferCrossSection;
				    		for(int i=0; i<global_promptProductionMatrix.size();i++){
				    			dvector globalpromptProductionVector=global_promptProductionMatrix.at(i);
				    			globalBufferCrossSection = globalpromptProductionVector.at(0);
				    			if(globalBufferCrossSection>1e-100) globalStatesContributing.push_back(i);
				    		}

				    		globalStatesCont=globalStatesContributing.size();
				    		StatesContributing_=globalStatesContributing;
				    		cout<<"globalStatesCont "<<globalStatesCont<<endl;

				    		for(int i=0; i<globalStatesContributing.size();i++){
				    			cout<<"state "<<globalStatesContributing[i]<<endl;
				    		}

						    if(!LoopThroughPPD){

						    	cout<<"!LoopThroughPPD"<<endl;
								//cout<<"Oi_MPV"<<endl;
								//cout<<Oi_MPV<<endl;
								//cout<<"Np_BR_MPV"<<endl;
								//cout<<Np_BR_MPV<<endl;
								//cout<<"Np_US_MPV"<<endl;
								//cout<<Np_US_MPV<<endl;

								ObjectLikelihoodVec=DataModelObject[iRap][iP]->getObjectLikelihood(Oi_MPV, Np_BR_MPV, Np_US_MPV, true, directProductionCube, promptProductionMatrix, polCorrFactor);
								ModelPrediction=ObjectLikelihoodVec[1];
								model_centralval.push_back(ModelPrediction);

								//cout<<"directProductionCube"<<endl;
								//cout<<directProductionCube<<endl;
								//cout<<"promptProductionMatrix"<<endl;
								//cout<<promptProductionMatrix<<endl;
								//cout<<"polCorrFactor"<<endl;
								//cout<<polCorrFactor<<endl;

								polCorrFactorVec.push_back(polCorrFactor);

								dcube Dummy_directProductionCube;
								dmatrix Dummy_promptProductionMatrix;
								double Dummy_polCorrFactor;

								ObjectLikelihoodVec=DataModelObject[iRap][iP]->getObjectLikelihood(Op_minus, Np_BR_MPV, Np_US_MPV, false, Dummy_directProductionCube, Dummy_promptProductionMatrix, Dummy_polCorrFactor);
								errModelPrediction=ObjectLikelihoodVec[1];
								model_errlow_centralval.push_back(fabs(ModelPrediction-errModelPrediction));

								model_errlow_centralval_1Sig.push_back(fabs(ModelPrediction-errModelPrediction)*1.);
								model_errlow_centralval_2Sig.push_back(fabs(ModelPrediction-errModelPrediction)*2.);
								model_errlow_centralval_3Sig.push_back(fabs(ModelPrediction-errModelPrediction)*3.);

								ObjectLikelihoodVec=DataModelObject[iRap][iP]->getObjectLikelihood(Op_plus, Np_BR_MPV, Np_US_MPV, false, Dummy_directProductionCube, Dummy_promptProductionMatrix, Dummy_polCorrFactor);
								errModelPrediction=ObjectLikelihoodVec[1];
								model_errhigh_centralval.push_back(fabs(ModelPrediction-errModelPrediction));

								model_errhigh_centralval_1Sig.push_back(fabs(ModelPrediction-errModelPrediction)*1.);
								model_errhigh_centralval_2Sig.push_back(fabs(ModelPrediction-errModelPrediction)*2.);
								model_errhigh_centralval_3Sig.push_back(fabs(ModelPrediction-errModelPrediction)*3.);


								dvector model_directProduction_THunc_CSVec;

								for(int iTHunc=0;iTHunc<3;iTHunc++){

									if(iTHunc==0) model_directProduction_THunc_CSVec.push_back(ModelPrediction*.9);
									if(iTHunc==1) model_directProduction_THunc_CSVec.push_back(ModelPrediction*1.1);
									if(iTHunc==2) model_directProduction_THunc_CSVec.push_back(ModelPrediction*1.2);

								}

								model_directProduction_THunc_CS.push_back(model_directProduction_THunc_CSVec);


								//PAPER: also here add to the vectors (to be initialized) the 1, 2, 3 sigma errors of the color channel predictions (as in the 'Loop'-case)



					    		//cout<<"Calculate Inclusive FeedDown"<<endl;
					    		//Calculate Inclusive FeedDown
					    		dvector FeedDownContVec;
					    		//dvector FeedDownLamthVec;
					    		//dvector FeedDownLamphVec;
					    		//dvector FeedDownLamtpVec;
					    		dvector FeedDownLamVec(3,0);
					    		dmatrix FeedDownLamMatrix;
					    		vector<int> StatesContributing;
					    		double FeedDownCont=0;
					    		double BufferCrossSection;
					    		for(int i=0; i<promptProductionMatrix.size();i++){
					    			dvector promptProductionVector=promptProductionMatrix.at(i);
					    			BufferCrossSection = promptProductionVector.at(0);
					    			if(BufferCrossSection>1e-100) StatesContributing.push_back(i);
					    			if(i!=iState){
										FeedDownCont+=BufferCrossSection;
										FeedDownContVec.push_back(BufferCrossSection);
										FeedDownLamVec.at(0)=promptProductionVector.at(1);
										FeedDownLamVec.at(1)=promptProductionVector.at(2);
										FeedDownLamVec.at(2)=promptProductionVector.at(3);
										FeedDownLamMatrix.push_back(FeedDownLamVec);
										//FeedDownLamthVec.push_back(promptProductionVector.at(1));
										//FeedDownLamphVec.push_back(promptProductionVector.at(2));
										//FeedDownLamtpVec.push_back(promptProductionVector.at(3));
					    			}
					    		}



					    		StatesCont=StatesContributing.size();
							    StatesContributing_=StatesContributing;

								 double model_inclusive_FeedDown_=0;

					    		if(StatesCont>2){

									dvector FeedDownInclusiveLambdas;
									 FeedDownInclusiveLambdas = addPolarizations(FeedDownLamMatrix, FeedDownContVec);



									if(iMeasurementID==0) model_inclusive_FeedDown_=FeedDownCont;
									else if(iMeasurementID>0 && iMeasurementID<4){
										model_inclusive_FeedDown_=FeedDownInclusiveLambdas[iMeasurementID-1];
									}

					    		}


					    		//cout<<"Calculate individual FeedDown contributions"<<endl;
					    		//Calculate individual FeedDown contributions
							    dvector Buff_model_FeedDownVec;
							    dvector model_FeedDownVec;
							    for(int i=0; i<StatesCont+1;i++){
							    	if(i==0){
							    		if(StatesCont>2) model_FeedDownVec.push_back(model_inclusive_FeedDown_);
							    		if(StatesCont==2){
											Buff_model_FeedDownVec=promptProductionMatrix.at(StatesContributing[1]);
											model_FeedDownVec.push_back(Buff_model_FeedDownVec[iMeasurementID]);
							    		}
							    	}
							    	else if(i>0&&i<StatesCont){
										Buff_model_FeedDownVec=promptProductionMatrix.at(StatesContributing[i]);
										model_FeedDownVec.push_back(Buff_model_FeedDownVec[iMeasurementID]);
							    	}
							    	else if(i==StatesCont){
										Buff_model_FeedDownVec=promptProductionMatrix.at(StatesContributing[0]);
										model_FeedDownVec.push_back(Buff_model_FeedDownVec[iMeasurementID]);
							    	}
							    }
							    model_FeedDown.push_back(model_FeedDownVec);


					    		//cout<<"Calculate individual color channel contributions of the directly produced state"<<endl;
					    		//Calculate individual color channel contributions of the directly produced state

							    dmatrix Buff_model_directProductionMatrix=directProductionCube.at(0);
							    dvector model_directProductionVec;
							    model_directProductionVec=Buff_model_directProductionMatrix.at(iMeasurementID);
							    model_directProduction.push_back(model_directProductionVec);

								//cout<<"model_directProductionVec: "<<endl;
								//cout<<model_directProductionVec<<endl;

									int nColorChannels_state;
									bool isSstate=(StateQuantumID[iState] > NRQCDvars::quID_S)?false:true;
									if(isSstate) nColorChannels_state=NRQCDvars::nColorChannels_S;
									else nColorChannels_state=NRQCDvars::nColorChannels_P;

									dvector model_directProduction_BandsVec;
								    dvector model_directProduction_BandsVec_errlow_1sig;
								    dvector model_directProduction_BandsVec_errhigh_1sig;
								    dvector model_directProduction_BandsVec_errlow_2sig;
								    dvector model_directProduction_BandsVec_errhigh_2sig;
								    dvector model_directProduction_BandsVec_errlow_3sig;
								    dvector model_directProduction_BandsVec_errhigh_3sig;

									int nSigmaBand;
									double buff_errlow, buff_errhigh;
									for (int j=0; j<nColorChannels_state; j++){

										model_directProduction_BandsVec.push_back(model_directProductionVec[j]);

										nSigmaBand=1; if(iMeasurementID>0) nSigmaBand=1.;
										model_directProduction_BandsVec_errlow_1sig.push_back( nSigmaBand*Oi_errlow[iState][j]/Oi_MPV[iState][j]*model_directProductionVec[j] );
										model_directProduction_BandsVec_errhigh_1sig.push_back(nSigmaBand*Oi_errhigh[iState][j]/Oi_MPV[iState][j]*model_directProductionVec[j]);
										nSigmaBand=2; if(iMeasurementID>0) nSigmaBand=2.;
										model_directProduction_BandsVec_errlow_2sig.push_back( nSigmaBand*Oi_errlow[iState][j]/Oi_MPV[iState][j]*model_directProductionVec[j] );
										model_directProduction_BandsVec_errhigh_2sig.push_back(nSigmaBand*Oi_errhigh[iState][j]/Oi_MPV[iState][j]*model_directProductionVec[j]);
										nSigmaBand=3; if(iMeasurementID>0) nSigmaBand=3.;
										model_directProduction_BandsVec_errlow_3sig.push_back( nSigmaBand*Oi_errlow[iState][j]/Oi_MPV[iState][j]*model_directProductionVec[j] );
										model_directProduction_BandsVec_errhigh_3sig.push_back(nSigmaBand*Oi_errhigh[iState][j]/Oi_MPV[iState][j]*model_directProductionVec[j]);

									}

									cout<<"model_directProduction_BandsVec"<<endl;
									cout<<model_directProduction_BandsVec<<endl;

									model_directProduction_Bands.push_back(model_directProduction_BandsVec);
									model_directProduction_Bands_errlow_1sig.push_back(model_directProduction_BandsVec_errlow_1sig);
									model_directProduction_Bands_errhigh_1sig.push_back(model_directProduction_BandsVec_errhigh_1sig);
									model_directProduction_Bands_errlow_2sig.push_back(model_directProduction_BandsVec_errlow_2sig);
									model_directProduction_Bands_errhigh_2sig.push_back(model_directProduction_BandsVec_errhigh_2sig);
									model_directProduction_Bands_errlow_3sig.push_back(model_directProduction_BandsVec_errlow_3sig);
									model_directProduction_Bands_errhigh_3sig.push_back(model_directProduction_BandsVec_errhigh_3sig);



						    }

						    else {

						    	cout<<"LOOP"<<endl;
						    	cout<<"LoopThroughPPD"<<endl;
								dcube MH_directProductionCube;
								dmatrix MH_promptProductionMatrix;

						    	sprintf(inname,"%s/results.root",jobdirname);
						    	TFile *ResultsFile = new TFile(inname, "READ");

						      	TTree*  outputTreeAllSamplings = (TTree*) ResultsFile->Get("AllSamplings");


						      	double PPD_Op[NRQCDvars::nStates][NRQCDvars::nColorChannels];
						      	double PPD_Np_US[2][max(NRQCDvars::nDataSystematicScales,NRQCDvars::nModelSystematicScales)];
						      	double PPD_Np_BR[NRQCDvars::nStates][NRQCDvars::nStates];
						      	double PPD_consts_star[NRQCDvars::nStates][NRQCDvars::nColorChannels];


						      	char branch_name[200];

						      	int BurnInInt;
								sprintf (branch_name,"BurnInInt");
								outputTreeAllSamplings->SetBranchAddress( branch_name,  &BurnInInt );
						      	int acceptedSampling;
								sprintf (branch_name,"acceptedSampling");
								outputTreeAllSamplings->SetBranchAddress( branch_name,  &acceptedSampling );


								for (int i=0; i<NRQCDvars::nStates; i++){
									int nColorChannels_state;
									bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
									if(isSstate) nColorChannels_state=NRQCDvars::nColorChannels_S;
									else nColorChannels_state=NRQCDvars::nColorChannels_P;
									for (int j=0; j<nColorChannels_state; j++){

										sprintf (branch_name,"state%d_Op%d",i,j);
										outputTreeAllSamplings->SetBranchAddress( branch_name,  &PPD_Op[i][j] );
										sprintf (branch_name,"state%d_const_star%d",i,j);
										outputTreeAllSamplings->SetBranchAddress( branch_name,  &PPD_consts_star[i][j] );

									}
								}

								for(int i=0; i < NRQCDvars::nStates; i++){
									for(int j=0; j < NRQCDvars::nStates; j++){
										if(NRQCDvars::FeedDownBranchingRatio[i][j]>0){

											sprintf (branch_name,"Np_BR_Daughter%d_Mother%d",i,j);
											outputTreeAllSamplings->SetBranchAddress( branch_name,  &PPD_Np_BR[i][j] );

										}
									}
								}

								for(int j=0; j < NRQCDvars::nDataSystematicScales; j++){
									sprintf (branch_name,"Np_US_DataSystematicScale%d",j);
									outputTreeAllSamplings->SetBranchAddress( branch_name,  &PPD_Np_US[0][j] );
								}
								for(int j=0; j < NRQCDvars::nModelSystematicScales; j++){
									sprintf (branch_name,"Np_US_ModelSystematicScale%d",j);
									outputTreeAllSamplings->SetBranchAddress( branch_name,  &PPD_Np_US[1][j] );
								}


								sprintf(outname,"%s/ModelPrediction_exp%d_meas%d_state%d_pt%d_rap%d.root",jobdirname, iExperiment, iMeasurementID, iState, iP, iRap);
						    	TFile *DummyFile = new TFile(outname, "RECREATE");

						      	TTree*  ModelPredictionTree;
						      	ModelPredictionTree = new TTree("ModelPredictionTree", "ModelPredictionTree");
								ModelPredictionTree->Branch("ModelPrediction",     &ModelPrediction,     "ModelPrediction/D");
								ModelPredictionTree->Branch("polCorrFactor",     &polCorrFactor,     "polCorrFactor/D");

								char branch_name_D[200];
								double CCprediction[NRQCDvars::nColorChannels];

								int nColorChannels_state;
								bool isSstate=(StateQuantumID[iState] > NRQCDvars::quID_S)?false:true;
								if(isSstate) nColorChannels_state=NRQCDvars::nColorChannels_S;
								else nColorChannels_state=NRQCDvars::nColorChannels_P;
								for (int j=0; j<nColorChannels_state; j++){

									sprintf (branch_name,"state%d_CC%d",iState,j);
									sprintf (branch_name_D,"%s/D",branch_name);
									if(FreeParam_Fractions_States.at(iState)==1)  ModelPredictionTree->Branch(branch_name,     &CCprediction[j],     branch_name_D);

								}

								double FeedDownContributions[globalStatesCont+1];
							    for(int iFD=0; iFD<globalStatesCont+1;iFD++){
									sprintf (branch_name,"FD%d",iFD);
									sprintf (branch_name_D,"%s/D",branch_name);
									ModelPredictionTree->Branch(branch_name,     &FeedDownContributions[iFD],     branch_name_D);
							    }


								dmatrix Op;
								vector<double> S_vector (NRQCDvars::nColorChannels_S,0);
								vector<double> P_vector (NRQCDvars::nColorChannels_P,0);
								Op=Oi_MPV;

								dmatrix Np_BR; //Branching ratios, [nDauthgers][nMothers]
								dvector Np_BR_0 (NRQCDvars::nStates,0.);
								for(int i=0; i < NRQCDvars::nStates; i++) Np_BR.push_back(Np_BR_0);

								dmatrix Np_US; //Uncertainty scales, [0=Data, 1=Model][nScales]
								dvector Np_US_0 (NRQCDvars::nDataSystematicScales, 0.);
								dvector Np_US_1 (NRQCDvars::nModelSystematicScales, 0.);
								Np_US.push_back(Np_US_0);
								Np_US.push_back(Np_US_1);

								dmatrix consts_star_var;
								consts_star_var=Oi_MPV;

								double MPmin=1e100;
								double MPmax=-1e100;

								int n_events = int( outputTreeAllSamplings->GetEntries() );
								cout<<n_events<<" events: Looping through nTuple of PPD"<<endl;
								int nFilled=0;

								// loop over  events in the model ntuple
								for ( int i_event = 1; i_event <= n_events; i_event++ ) {

									//cout<<i_event<<endl;

									outputTreeAllSamplings->GetEvent( i_event-1 );

									if(acceptedSampling!=1 || BurnInInt==1) continue;

									//cout<<"set Op"<<endl;

									for (int i=0; i<NRQCDvars::nStates; i++){
										bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
										if(isSstate){
											for (int j=0; j<NRQCDvars::nColorChannels_S; j++){
												S_vector.at(j)=PPD_Op[i][j];
											}
											Op.at(i)=S_vector;
										}
										else{
											for (int j=0; j<NRQCDvars::nColorChannels_P; j++){
												P_vector.at(j)=PPD_Op[i][j];
											}
											Op.at(i)=P_vector;
										}
									}

									//cout<<"set consts_star_var"<<endl;
									for (int i=0; i<NRQCDvars::nStates; i++){
										bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
										if(isSstate){
											for (int j=0; j<NRQCDvars::nColorChannels_S; j++){
												S_vector.at(j)=PPD_consts_star[i][j];
											}
											consts_star_var.at(i)=S_vector;
										}
										else{
											for (int j=0; j<NRQCDvars::nColorChannels_P; j++){
												P_vector.at(j)=PPD_consts_star[i][j];
											}
											consts_star_var.at(i)=P_vector;
										}
									}

									//cout<<"set Np_BR"<<endl;
									for(int i=0; i < NRQCDvars::nStates; i++){
										for(int j=0; j < NRQCDvars::nStates; j++){
											if(NRQCDvars::FeedDownBranchingRatio[i][j]>0){
												Np_BR[i][j]=PPD_Np_BR[i][j];
											}
											else Np_BR[i][j]=0;
										}
									}

									//cout<<"set Np_US"<<endl;
									for(int j=0; j < NRQCDvars::nDataSystematicScales; j++){
										Np_US[0][j]=PPD_Np_US[0][j];
										//Np_US[0][j]=20.;
									}
									for(int j=0; j < NRQCDvars::nModelSystematicScales; j++){
										Np_US[1][j]=PPD_Np_US[1][j];
									}

									dcube directProductionCube_lonely;
									dmatrix directProductionMatrix_lonely;
									ObjectLikelihoodVec=DataModelObject[iRap][iP]->getObjectLikelihood(Op, Np_BR, Np_US, true, directProductionCube_lonely, MH_promptProductionMatrix, polCorrFactor);
									ModelPrediction=ObjectLikelihoodVec[1];

									//cout<<"Np_US: "<<endl;
									//cout<<Np_US<<endl;

									//if(iMeasurementID==0){
										//cout<<"directProductionMatrix_lonely: "<<endl;
										//dvector directProductionVector = DataModelObject[iRap][iP]->getDirectProduction(DataModelObject[iRap][iP]->getState(), Op, Np_BR, Np_US, true, directProductionMatrix_lonely);
										//cout<<directProductionMatrix_lonely<<endl;
										//cout<<"directProductionVector: "<<endl;
										//cout<<directProductionVector<<endl;

									    dvector model_directProductionVec_lonely;
									    dmatrix Buff_model_directProductionMatrix_lonely=directProductionCube_lonely.at(0);

										//cout<<"Buff_model_directProductionMatrix_lonely"<<endl;
										//cout<<Buff_model_directProductionMatrix_lonely<<endl;

									    model_directProductionVec_lonely=Buff_model_directProductionMatrix_lonely.at(iMeasurementID);
									    //model_directProduction_lonely.push_back(model_directProductionVec_lonely);
										//cout<<"model_directProductionVec_lonely: "<<endl;
										//cout<<model_directProductionVec_lonely<<endl;


										//for (int i=0; i<NRQCDvars::nStates; i++){
											int nColorChannels_state;
											bool isSstate=(StateQuantumID[iState] > NRQCDvars::quID_S)?false:true;
											if(isSstate) nColorChannels_state=NRQCDvars::nColorChannels_S;
											else nColorChannels_state=NRQCDvars::nColorChannels_P;
											for (int j=0; j<nColorChannels_state; j++){

												CCprediction[j]=model_directProductionVec_lonely[j];
												//continue;
												//if(j==1) cout<<"CCprediction[1] "<<CCprediction[j]<<endl;

											}
										//}

									//}


											//cout<<"calc FD"<<endl;

								    		dvector FeedDownContVec;
								    		dvector FeedDownLamVec(3,0);
								    		dmatrix FeedDownLamMatrix;
								    		double FeedDownCont=0;
								    		double BufferCrossSection;
								    		for(int i=0; i<globalStatesCont;i++){
								    			dvector promptProductionVector=MH_promptProductionMatrix.at(globalStatesContributing[i]);
								    			BufferCrossSection = promptProductionVector.at(0);
								    			if(globalStatesContributing[i]!=iState){
													FeedDownCont+=BufferCrossSection;
													FeedDownContVec.push_back(BufferCrossSection);
													FeedDownLamVec.at(0)=promptProductionVector.at(1);
													FeedDownLamVec.at(1)=promptProductionVector.at(2);
													FeedDownLamVec.at(2)=promptProductionVector.at(3);
													FeedDownLamMatrix.push_back(FeedDownLamVec);
													//cout<<"globalStatesContributing in: "<<BufferCrossSection<<endl;
								    			}
								    		}

											//if(iP==8) cout<<"FeedDownCont "<<FeedDownCont<<endl;

											double model_inclusive_FeedDown_=0;

								    		if(globalStatesCont>2){
												//cout<<"calc inclusive FD"<<endl;

												if(iMeasurementID==0) model_inclusive_FeedDown_=FeedDownCont;
												else if(iMeasurementID>0 && iMeasurementID<4){
													dvector FeedDownInclusiveLambdas;
													FeedDownInclusiveLambdas = addPolarizations(FeedDownLamMatrix, FeedDownContVec);
													model_inclusive_FeedDown_=FeedDownInclusiveLambdas[iMeasurementID-1];
												}

								    		}

											//if(iP==8) cout<<"model_inclusive_FeedDown_ "<<model_inclusive_FeedDown_<<endl;

								    		//Calculate individual FeedDown contributions
										    dvector Buff_model_FeedDownVec;
										    dvector model_FeedDownVec;
										    for(int i=0; i<globalStatesCont+1;i++){
												//cout<<"calc individual FD"<<endl;
										    	if(i==0){
										    		if(globalStatesCont>2) FeedDownContributions[i]=model_inclusive_FeedDown_;
										    		if(globalStatesCont==2){
														Buff_model_FeedDownVec=MH_promptProductionMatrix.at(globalStatesContributing[1]);
														FeedDownContributions[i]=Buff_model_FeedDownVec[iMeasurementID];
										    		}
										    	}
										    	else if(i>0&&i<globalStatesCont){
													Buff_model_FeedDownVec=MH_promptProductionMatrix.at(globalStatesContributing[i]);
													FeedDownContributions[i]=Buff_model_FeedDownVec[iMeasurementID];
										    	}
										    	else if(i==globalStatesCont){
													Buff_model_FeedDownVec=MH_promptProductionMatrix.at(globalStatesContributing[0]);
													FeedDownContributions[i]=Buff_model_FeedDownVec[iMeasurementID];
										    	}
												//cout<<"globalStatesContributing out: "<<i<<endl;
												//cout<<"FeedDownContributions[i] "<<FeedDownContributions[i]<<endl;

										    }

											//if(iExperiment==0 && iMeasurementID==0 && iRap==0 && iP==12){
											//	cout<<"ModelPrediction "<<ModelPrediction<<endl;
											//}

										    if(ModelPrediction<MPmin){
										    	if(TMath::Abs(ModelPrediction/MPmin)<1e-5 && nFilled>1){
										    		cout<<"ignored sampling because ModelPrediction ("<<ModelPrediction<<") seems to be out of line - too low..."<<endl;
										    		continue;
										    	}
										    	MPmin=ModelPrediction;
										    	cout<<"new MPmin = "<<ModelPrediction<<" (nFilled="<<nFilled<<")"<<endl;
										    }
										    if(ModelPrediction>MPmax){
										    	if(TMath::Abs(ModelPrediction/MPmax)>1e5 && nFilled>1){
										    		cout<<"ignored sampling because ModelPrediction ("<<ModelPrediction<<") seems to be out of line - too large..."<<endl;
										    		continue;
										    	}
										    	MPmax=ModelPrediction;
										    	cout<<"new MPmax = "<<ModelPrediction<<" (nFilled="<<nFilled<<")"<<endl;
										    }

										    nFilled++;
									ModelPredictionTree->Fill();

									//PAPER: Add the color channel 'prediction'-values to the ModelPredictionTree


								}

								//ModelPredictionTree->Print();

								//PAPER: add here the histos of the individual color channels, get FindMPV for 1, 2, 3 sigmas

								int nBins_h=100;
								const int nHists=2;
								char hist_name[200];
								char projectchar[200];
								char selectchar[200];
								TH1D* h_ModelPrediction[nHists];
								double h_ModelPrediction_min[nHists];
								double h_ModelPrediction_max[nHists];

							  	double expandMinMaxBy=0.01;
							  	double buff_MPV;
							  	double buff_errlow;
							  	double buff_errhigh;

							  	h_ModelPrediction_min[0]=ModelPredictionTree->GetMinimum("ModelPrediction")-expandMinMaxBy*TMath::Abs(ModelPredictionTree->GetMinimum("ModelPrediction"));
							  	h_ModelPrediction_max[0]=ModelPredictionTree->GetMaximum("ModelPrediction")+expandMinMaxBy*TMath::Abs(ModelPredictionTree->GetMaximum("ModelPrediction"));
							  	//cout<<"h_ModelPrediction_min[0] "<<h_ModelPrediction_min[0]<<" h_ModelPrediction_max[0] "<<h_ModelPrediction_max[0]<<endl;
							  	sprintf(hist_name,"h_ModelPrediction");
								h_ModelPrediction[0] = new TH1D( hist_name, hist_name, nBins_h, h_ModelPrediction_min[0], h_ModelPrediction_max[0] );
								sprintf(projectchar,"ModelPrediction>>h_ModelPrediction");
								ModelPredictionTree->Draw(projectchar);

							  	h_ModelPrediction_min[1]=ModelPredictionTree->GetMinimum("polCorrFactor")-expandMinMaxBy*TMath::Abs(ModelPredictionTree->GetMinimum("polCorrFactor"));
							  	h_ModelPrediction_max[1]=ModelPredictionTree->GetMaximum("polCorrFactor")+expandMinMaxBy*TMath::Abs(ModelPredictionTree->GetMaximum("polCorrFactor"));
								sprintf(hist_name,"h_polCorrFactor");
								h_ModelPrediction[1] = new TH1D( hist_name, hist_name, nBins_h, h_ModelPrediction_min[1], h_ModelPrediction_max[1] );
								sprintf(projectchar,"polCorrFactor>>h_polCorrFactor");
								ModelPredictionTree->Draw(projectchar);


								FindMPV(h_ModelPrediction[0], buff_MPV, buff_errlow, buff_errhigh, MPValgo, nSigma);
								cout<<"nSigma "<<nSigma<<endl;
								cout<<"MPValgo "<<MPValgo<<endl;
								cout<<"buff_MPV "<<buff_MPV<<endl;
								cout<<"buff_errlow "<<buff_errlow<<endl;
								cout<<"buff_errhigh "<<buff_errhigh<<endl;

								double buff_errlow_Bands[3];
								double buff_errhigh_Bands[3];

								FindMPV(h_ModelPrediction[0], buff_MPV, buff_errlow_Bands[0], buff_errhigh_Bands[0], MPValgo, 1);
								FindMPV(h_ModelPrediction[0], buff_MPV, buff_errlow_Bands[1], buff_errhigh_Bands[1], MPValgo, 2);
								FindMPV(h_ModelPrediction[0], buff_MPV, buff_errlow_Bands[2], buff_errhigh_Bands[2], MPValgo, 3);

								if(iExperiment==0 && iMeasurementID==0 && iRap==0 && iP==12){
									h_ModelPrediction[0]->SaveAs("tmp_Modelpred.root");
									h_ModelPrediction[1]->SaveAs("tmp_polcorrpred.root");
									h_ModelPrediction[0]->Print();
									cout<<h_ModelPrediction[0]->GetMean()<<endl;
									cout<<h_ModelPrediction[0]->GetRMS()<<endl;
									cout<<buff_MPV<<endl;
									cout<<buff_errlow<<endl;
									cout<<buff_errhigh<<endl;
								}

								if(iExperiment==0 && iMeasurementID==0 && iRap==0 && iP==12){
									ModelPredictionTree->SaveAs("tmp_ModelpredTree_CS.root");
								}
								if(iExperiment==0 && iMeasurementID==1 && iRap==0 && iP==12){
									ModelPredictionTree->SaveAs("tmp_ModelpredTree_pol.root");
								}

								//ObjectLikelihoodVec=DataModelObject[iRap][iP]->getObjectLikelihood(Oi_MPV, Np_BR_MPV, Np_US_MPV, true, directProductionCube, promptProductionMatrix, polCorrFactor);
								//ModelPrediction=ObjectLikelihoodVec[1];
                                //
								//cout<<"MPV directProductionCube"<<endl;
								//cout<<directProductionCube<<endl;
								//cout<<"MPV promptProductionMatrix"<<endl;
								//cout<<promptProductionMatrix<<endl;
								//cout<<"MPV polCorrFactor"<<endl;
								//cout<<polCorrFactor<<endl;

								//PAPER: probably use MPV as central model, given that we now display the bands, and not the central values anymore
								//bool useHistModeAsCentralModel=true;

								//if(useHistModeAsCentralModel){
									model_centralval.push_back(buff_MPV);
									model_errlow_centralval.push_back(fabs(buff_errlow));
									model_errhigh_centralval.push_back(fabs(buff_errhigh));
									cout<<"model_centralval "<<buff_MPV<<endl;
									cout<<"model_errlow_centralval "<<fabs(buff_errlow)<<endl;
									cout<<"model_errhigh_centralval "<<fabs(buff_errhigh)<<endl;
									FindMPV(h_ModelPrediction[1], buff_MPV, buff_errlow, buff_errhigh, MPValgo, nSigma);
									polCorrFactorVec.push_back(buff_MPV);
									model_errlow_centralval_1Sig.push_back(fabs(buff_errlow_Bands[0]));
									model_errhigh_centralval_1Sig.push_back(fabs(buff_errhigh_Bands[0]));
									model_errlow_centralval_2Sig.push_back(fabs(buff_errlow_Bands[1]));
									model_errhigh_centralval_2Sig.push_back(fabs(buff_errhigh_Bands[1]));
									model_errlow_centralval_3Sig.push_back(fabs(buff_errlow_Bands[2]));
									model_errhigh_centralval_3Sig.push_back(fabs(buff_errhigh_Bands[2]));

								//}
								//else{
								//	model_centralval.push_back(ModelPrediction);
								//	//model_errlow_centralval.push_back(fabs(buff_errlow)+ModelPrediction-buff_MPV);
								//	//model_errhigh_centralval.push_back(fabs(buff_errhigh)-ModelPrediction+buff_MPV);
								//	model_errlow_centralval.push_back(fabs(buff_errlow));
								//	model_errhigh_centralval.push_back(fabs(buff_errhigh));
								//	polCorrFactorVec.push_back(polCorrFactor);
								//	model_errlow_centralval_1Sig.push_back(fabs(buff_errlow_Bands[0]));
								//	model_errhigh_centralval_1Sig.push_back(fabs(buff_errhigh_Bands[0]));
								//	model_errlow_centralval_2Sig.push_back(fabs(buff_errlow_Bands[1]));
								//	model_errhigh_centralval_2Sig.push_back(fabs(buff_errhigh_Bands[1]));
								//	model_errlow_centralval_3Sig.push_back(fabs(buff_errlow_Bands[2]));
								//	model_errhigh_centralval_3Sig.push_back(fabs(buff_errhigh_Bands[2]));
                                //
								//	cout<<"model_centralval "<<ModelPrediction<<endl;
								//	cout<<"model_errlow_centralval "<<fabs(buff_errlow)<<endl;
								//	cout<<"model_errhigh_centralval "<<fabs(buff_errhigh)<<endl;
								//}



								cout<<"iMeasurementID "<<iMeasurementID<<endl;
								cout<<"iP "<<iP<<endl;

								//if(iMeasurementID==0){
									if(isSstate) nColorChannels_state=NRQCDvars::nColorChannels_S;
									else nColorChannels_state=NRQCDvars::nColorChannels_P;


									double evalTHunc[3]={-1., 0., 1.};
									dvector model_directProduction_THunc_CSVec;

									for(int iTHunc=0;iTHunc<3;iTHunc++){

										for(int j=0; j < NRQCDvars::nDataSystematicScales; j++){
											Np_US[0][j]=0.;
										}
										for(int j=0; j < NRQCDvars::nModelSystematicScales; j++){
											Np_US[1][j]=evalTHunc[iTHunc];
										}

										dcube model_directProduction_THunc_CS_lonely;
										dmatrix model_directProduction_THunc_CSMatrix_lonely;
										ObjectLikelihoodVec=DataModelObject[iRap][iP]->getObjectLikelihood(Oi_MPV, Np_BR, Np_US, true, model_directProduction_THunc_CS_lonely, model_directProduction_THunc_CSMatrix_lonely, polCorrFactor);
										dvector model_directProduction_THunc_CSVec_lonely;
										dmatrix Buff_model_directProductionCSMatrix_lonely=model_directProduction_THunc_CS_lonely.at(0);
										model_directProduction_THunc_CSVec_lonely=Buff_model_directProductionCSMatrix_lonely.at(iMeasurementID);

										int CC=0;
										model_directProduction_THunc_CSVec.push_back(model_directProduction_THunc_CSVec_lonely[CC]);

											cout<<"THunc"<<iTHunc<<" dddd: "<<model_directProduction_THunc_CSVec_lonely[CC]<<endl;
									}

									model_directProduction_THunc_CS.push_back(model_directProduction_THunc_CSVec);





									dvector model_directProduction_BandsVec;
								    dvector model_directProduction_BandsVec_errlow_1sig;
								    dvector model_directProduction_BandsVec_errhigh_1sig;
								    dvector model_directProduction_BandsVec_errlow_2sig;
								    dvector model_directProduction_BandsVec_errhigh_2sig;
								    dvector model_directProduction_BandsVec_errlow_3sig;
								    dvector model_directProduction_BandsVec_errhigh_3sig;

									for (int j=0; j<nColorChannels_state; j++){

										sprintf (branch_name,"state%d_CC%d",iState,j);
										cout<<branch_name<<endl;

										int nBins_h_CC=100;
										char hist_name_CC[200];
										char projectchar_CC[200];
										char selectchar_CC[200];
										TH1D* h_ModelPrediction_CC;
										double h_ModelPrediction_min_CC;
										double h_ModelPrediction_max_CC;

										double expandMinMaxBy=0.01;
										double buff_MPV;
										double buff_errlow;
										double buff_errhigh;

										h_ModelPrediction_min_CC=ModelPredictionTree->GetMinimum(branch_name)-expandMinMaxBy*TMath::Abs(ModelPredictionTree->GetMinimum(branch_name));
										h_ModelPrediction_max_CC=ModelPredictionTree->GetMaximum(branch_name)+expandMinMaxBy*TMath::Abs(ModelPredictionTree->GetMaximum(branch_name));
										sprintf(hist_name_CC,"h_ModelPrediction_CC");
										h_ModelPrediction_CC = new TH1D( hist_name_CC, hist_name_CC, nBins_h_CC, h_ModelPrediction_min_CC, h_ModelPrediction_max_CC );
										sprintf(projectchar_CC,"%s>>h_ModelPrediction_CC",branch_name);
										ModelPredictionTree->Draw(projectchar_CC);

										if(j==0){
											cout<<"h_ModelPrediction_min_CC "<<h_ModelPrediction_min_CC<<endl;
											cout<<"ModelPredictionTree->GetMinimum(branch_name); "<<ModelPredictionTree->GetMinimum(branch_name)<<endl;
											h_ModelPrediction_CC->Print("all");
										}
										//if(j==1) h_ModelPrediction_CC->Print("all");
										FindMPV(h_ModelPrediction_CC, buff_MPV, buff_errlow, buff_errhigh, MPValgo, 1);
										cout<<"ColorChannelsCalc buff_MPV "<<buff_MPV<<endl;
										cout<<"ColorChannelsCalc buff_errlow "<<buff_errlow<<endl;
										cout<<"ColorChannelsCalc buff_errhigh "<<buff_errhigh<<endl;
										model_directProduction_BandsVec.push_back(buff_MPV);
										model_directProduction_BandsVec_errlow_1sig.push_back(buff_errlow);
										model_directProduction_BandsVec_errhigh_1sig.push_back(buff_errhigh);
										FindMPV(h_ModelPrediction_CC, buff_MPV, buff_errlow, buff_errhigh, MPValgo, 2);
										model_directProduction_BandsVec_errlow_2sig.push_back(buff_errlow);
										model_directProduction_BandsVec_errhigh_2sig.push_back(buff_errhigh);
										FindMPV(h_ModelPrediction_CC, buff_MPV, buff_errlow, buff_errhigh, MPValgo, 3);
										model_directProduction_BandsVec_errlow_3sig.push_back(buff_errlow);
										model_directProduction_BandsVec_errhigh_3sig.push_back(buff_errhigh);





									}

									cout<<"ColorChannelsCalc model_directProduction_BandsVec"<<endl;
									cout<<model_directProduction_BandsVec<<endl;
									cout<<"ColorChannelsCalc model_directProduction_Bands_errlow_1sig"<<endl;
									cout<<model_directProduction_Bands_errlow_1sig<<endl;

									model_directProduction_Bands.push_back(model_directProduction_BandsVec);
									model_directProduction_Bands_errlow_1sig.push_back(model_directProduction_BandsVec_errlow_1sig);
									model_directProduction_Bands_errhigh_1sig.push_back(model_directProduction_BandsVec_errhigh_1sig);
									model_directProduction_Bands_errlow_2sig.push_back(model_directProduction_BandsVec_errlow_2sig);
									model_directProduction_Bands_errhigh_2sig.push_back(model_directProduction_BandsVec_errhigh_2sig);
									model_directProduction_Bands_errlow_3sig.push_back(model_directProduction_BandsVec_errlow_3sig);
									model_directProduction_Bands_errhigh_3sig.push_back(model_directProduction_BandsVec_errhigh_3sig);

								//}




								    dvector model_FeedDownVec;
								    for(int i=0; i<globalStatesCont+1;i++){

										sprintf (branch_name,"FD%d",i);

										int nBins_h_FD=100;
										char hist_name_FD[200];
										char projectchar_FD[200];
										char selectchar_FD[200];
										TH1D* h_ModelPrediction_FD;
										double h_ModelPrediction_min_FD;
										double h_ModelPrediction_max_FD;

										double expandMinMaxBy=0.01;
										double buff_MPV;
										double buff_errlow;
										double buff_errhigh;

										h_ModelPrediction_min_FD=ModelPredictionTree->GetMinimum(branch_name)-expandMinMaxBy*TMath::Abs(ModelPredictionTree->GetMinimum(branch_name));
										h_ModelPrediction_max_FD=ModelPredictionTree->GetMaximum(branch_name)+expandMinMaxBy*TMath::Abs(ModelPredictionTree->GetMaximum(branch_name));
										sprintf(hist_name_FD,"h_ModelPrediction_FD");
										h_ModelPrediction_FD = new TH1D( hist_name_FD, hist_name_FD, nBins_h_FD, h_ModelPrediction_min_FD, h_ModelPrediction_max_FD );
										sprintf(projectchar_FD,"%s>>h_ModelPrediction_FD",branch_name);
										ModelPredictionTree->Draw(projectchar_FD);

										if(iP==8) h_ModelPrediction_FD->Print("all");
										FindMPV(h_ModelPrediction_FD, buff_MPV, buff_errlow, buff_errhigh, MPValgo, 1);
										cout<<"FD_MPV "<<buff_MPV<<endl;
										cout<<"FD_errlow "<<buff_errlow<<endl;
										cout<<"FD_errhigh "<<buff_errhigh<<endl;

										model_FeedDownVec.push_back(buff_MPV);

								    }

								    cout<<"model_FeedDownVec:"<<endl;
								    cout<<model_FeedDownVec<<endl;

								    model_FeedDown.push_back(model_FeedDownVec);









								ModelPredictionTree->Write();
								DummyFile->Close();

								delete DummyFile;
								delete ResultsFile;

								//gSystem->Unlink(outname);
								//cout<<"end loop"<<endl;

						    }

							//PAPER: store here the MPV and errors of the color channels in objects still to create

				    		data_centralval.push_back(DataModelObject[iRap][iP]->getCentralValue());//correct data to match predicted polariztion
				    		data_errlow_centralval.push_back(DataModelObject[iRap][iP]->getErrTotNeg());
				    		data_errhigh_centralval.push_back(DataModelObject[iRap][iP]->getErrTotPos());


				    		if(!PlotVSpToverM){
								data_pTmean.push_back(DataModelObject[iRap][iP]->getpTMean());
								data_errlow_pT.push_back(DataModelObject[iRap][iP]->getpTMean()-DataModelObject[iRap][iP]->getpTMin());
								data_errhigh_pT.push_back(DataModelObject[iRap][iP]->getpTMax()-DataModelObject[iRap][iP]->getpTMean());
				    		}
				    		else{
								data_pTmean.push_back(DataModelObject[iRap][iP]->getpTMean()/NRQCDvars::mass[iState]);
								data_errlow_pT.push_back(DataModelObject[iRap][iP]->getpTMean()/NRQCDvars::mass[iState]-DataModelObject[iRap][iP]->getpTMin()/NRQCDvars::mass[iState]);
								data_errhigh_pT.push_back(DataModelObject[iRap][iP]->getpTMax()/NRQCDvars::mass[iState]-DataModelObject[iRap][iP]->getpTMean()/NRQCDvars::mass[iState]);
				    		}




							//TODO: after ifs: calc ObjetLikelihood with 'central' inputs, to be used for component calculations
							ObjectLikelihoodVec=DataModelObject[iRap][iP]->getObjectLikelihood(Oi_MPV, Np_BR_MPV, Np_US_MPV, true, directProductionCube, promptProductionMatrix, polCorrFactor);


						    }


				    }


				    //cout<<"model_directProduction_THunc_CS"<<endl;
				    //cout<<model_directProduction_THunc_CS<<endl;

					const int StatesCont_c=globalStatesCont;
					int nColorChannels_state;
					bool isSstate=(StateQuantumID[iState] > NRQCDvars::quID_S)?false:true;
					if(isSstate) nColorChannels_state=NRQCDvars::nColorChannels_S;
					else nColorChannels_state=NRQCDvars::nColorChannels_P;
					const int ColorChannels_c=nColorChannels_state;


				    if(nPtBinsSel>0){
						cout << "Plot iState=" << StateName[iState] << ", iMeasurementID=" << MeasurementIDName[iMeasurementID] << " from "<< ExpName[iExperiment] << endl;
				    	cout<<"rap "<<iRap+1<<endl;

					    cout<<"model_FeedDown all:"<<endl;
					    cout<<model_FeedDown<<endl;

						const int nPtBinsSel_=nPtBinsSel;
						double d_data_centralval[nPtBinsSel_];
						double d_data_errlow_centralval[nPtBinsSel_];
						double d_data_errhigh_centralval[nPtBinsSel_];
						double d_data_pTmean[nPtBinsSel_];
						double d_data_errlow_pT[nPtBinsSel_];
						double d_data_errhigh_pT[nPtBinsSel_];
						double d_model_centralval[nPtBinsSel_];
						double d_data_centralval_nopolcorr[nPtBinsSel_];
						double d_model_errlow_centralval[nPtBinsSel_];
						double d_model_errhigh_centralval[nPtBinsSel_];
						double d_model_errlow_absolute_centralval[nPtBinsSel_];
						double d_model_errhigh_absolute_centralval[nPtBinsSel_];

						double d_model_errlow_centralval_1Sig[nPtBinsSel_];
						double d_model_errhigh_centralval_1Sig[nPtBinsSel_];
						double d_model_errlow_absolute_centralval_1Sig[nPtBinsSel_];
						double d_model_errhigh_absolute_centralval_1Sig[nPtBinsSel_];
						double d_model_errlow_centralval_2Sig[nPtBinsSel_];
						double d_model_errhigh_centralval_2Sig[nPtBinsSel_];
						double d_model_errlow_absolute_centralval_2Sig[nPtBinsSel_];
						double d_model_errhigh_absolute_centralval_2Sig[nPtBinsSel_];
						double d_model_errlow_centralval_3Sig[nPtBinsSel_];
						double d_model_errhigh_centralval_3Sig[nPtBinsSel_];
						double d_model_errlow_absolute_centralval_3Sig[nPtBinsSel_];
						double d_model_errhigh_absolute_centralval_3Sig[nPtBinsSel_];

						double d_model_CS_THunc_p1[nPtBinsSel_];
						double d_model_CS_THunc_m1[nPtBinsSel_];
						double d_model_CS_THunc_central[nPtBinsSel_];

						double d_model_CS_THunc_delta_p1_m1_diffhalf[nPtBinsSel_];
						double d_model_CS_THunc_delta_p1_m1_mean[nPtBinsSel_];
						double d_zero[nPtBinsSel_];

						for(int iP = 0; iP < nPtBinsSel; iP++){

							bool correctForPol=true;
							bool correctForLumi=true;

							//if(iP<3) correctForPol=false;

							if(data_pTmean[iP]<pTMinModel || data_pTmean[iP]>pTMaxModel) correctForPol=false;

							if(!correctForPol) polCorrFactorVec[iP]=1.;
							else if(correctForPol && iMeasurementID==0) cout<<"polCorrFactor (* to data and all models, incl submodels) for pT bin "<<iP<<" = "<<polCorrFactorVec[iP]<<endl;
							if(!correctForLumi) lumiCorrFactorVec[iP]=1.;
							else if(correctForLumi && iMeasurementID==0) cout<<"lumiCorrFactor (* to data and total models only) for pT bin "<<iP<<" = "<<lumiCorrFactorVec[iP]<<endl;

							d_data_centralval[iP] =         	data_centralval[iP];
							d_data_errlow_centralval[iP] =  	data_errlow_centralval[iP];
							d_data_errhigh_centralval[iP] = 	data_errhigh_centralval[iP];
							d_data_pTmean[iP] =             	data_pTmean[iP];
							d_data_errlow_pT[iP] =          	data_errlow_pT[iP];
							d_data_errhigh_pT[iP] =         	data_errhigh_pT[iP];
							d_model_centralval[iP] =        	model_centralval[iP];
							d_data_centralval_nopolcorr[iP] =         	data_centralval[iP];
							d_model_errlow_centralval[iP] = 	model_errlow_centralval[iP];
							d_model_errhigh_centralval[iP] =	model_errhigh_centralval[iP];
							d_model_errlow_absolute_centralval[iP] = 	(model_centralval[iP]-model_errlow_centralval[iP]);
							d_model_errhigh_absolute_centralval[iP] =	(model_centralval[iP]+model_errhigh_centralval[iP]);

							d_model_errlow_centralval_1Sig[iP] = 	model_errlow_centralval_1Sig[iP];
							d_model_errhigh_centralval_1Sig[iP] =	model_errhigh_centralval_1Sig[iP];
							d_model_errlow_absolute_centralval_1Sig[iP] = 	(model_centralval[iP]-model_errlow_centralval_1Sig[iP]);
							d_model_errhigh_absolute_centralval_1Sig[iP] =	(model_centralval[iP]+model_errhigh_centralval_1Sig[iP]);

							d_model_errlow_centralval_2Sig[iP] = 	model_errlow_centralval_2Sig[iP];
							d_model_errhigh_centralval_2Sig[iP] =	model_errhigh_centralval_2Sig[iP];
							d_model_errlow_absolute_centralval_2Sig[iP] = 	(model_centralval[iP]-model_errlow_centralval_2Sig[iP]);
							d_model_errhigh_absolute_centralval_2Sig[iP] =	(model_centralval[iP]+model_errhigh_centralval_2Sig[iP]);

							d_model_errlow_centralval_3Sig[iP] = 	model_errlow_centralval_3Sig[iP];
							d_model_errhigh_centralval_3Sig[iP] =	model_errhigh_centralval_3Sig[iP];
							d_model_errlow_absolute_centralval_3Sig[iP] = 	(model_centralval[iP]-model_errlow_centralval_3Sig[iP]);
							d_model_errhigh_absolute_centralval_3Sig[iP] =	(model_centralval[iP]+model_errhigh_centralval_3Sig[iP]);

							d_model_CS_THunc_m1[iP]=model_directProduction_THunc_CS[iP][0];
							d_model_CS_THunc_central[iP]=model_directProduction_THunc_CS[iP][1];
							d_model_CS_THunc_p1[iP]=model_directProduction_THunc_CS[iP][2];

							d_model_CS_THunc_delta_p1_m1_diffhalf[iP]=TMath::Abs(model_directProduction_THunc_CS[iP][2]-model_directProduction_THunc_CS[iP][0])/2.;
							d_model_CS_THunc_delta_p1_m1_mean[iP]=(model_directProduction_THunc_CS[iP][0]+model_directProduction_THunc_CS[iP][2])/2.;
							d_zero[iP] =         	0.;

							if(iMeasurementID==0){
								d_data_centralval[iP] *=         	polCorrFactorVec[iP];//polCorrect the data (according to polarization prediction of model)
								d_data_errlow_centralval[iP] *=  	polCorrFactorVec[iP];
								d_data_errhigh_centralval[iP] *= 	polCorrFactorVec[iP];
								d_model_centralval[iP] *=        	polCorrFactorVec[iP];//un-polCorrect the model, as it was corrected in the fit (getCorrPromptProduction...)
								d_model_errlow_centralval[iP] *= 	polCorrFactorVec[iP];
								d_model_errhigh_centralval[iP] *=	polCorrFactorVec[iP];
								d_model_errlow_absolute_centralval[iP] *= 	polCorrFactorVec[iP];
								d_model_errhigh_absolute_centralval[iP] *=	polCorrFactorVec[iP];

								d_model_errlow_centralval_1Sig[iP] *= 	polCorrFactorVec[iP];
								d_model_errhigh_centralval_1Sig[iP] *=	polCorrFactorVec[iP];
								d_model_errlow_absolute_centralval_1Sig[iP] *= 	polCorrFactorVec[iP];
								d_model_errhigh_absolute_centralval_1Sig[iP] *=	polCorrFactorVec[iP];

								d_model_errlow_centralval_2Sig[iP] *= 	polCorrFactorVec[iP];
								d_model_errhigh_centralval_2Sig[iP] *=	polCorrFactorVec[iP];
								d_model_errlow_absolute_centralval_2Sig[iP] *= 	polCorrFactorVec[iP];
								d_model_errhigh_absolute_centralval_2Sig[iP] *=	polCorrFactorVec[iP];

								d_model_errlow_centralval_3Sig[iP] *= 	polCorrFactorVec[iP];
								d_model_errhigh_centralval_3Sig[iP] *=	polCorrFactorVec[iP];
								d_model_errlow_absolute_centralval_3Sig[iP] *= 	polCorrFactorVec[iP];
								d_model_errhigh_absolute_centralval_3Sig[iP] *=	polCorrFactorVec[iP];

								//d_model_CS_THunc_m1[iP] *=	polCorrFactorVec[iP];
								//d_model_CS_THunc_central[iP] *=	polCorrFactorVec[iP];
								//d_model_CS_THunc_p1[iP] *=	polCorrFactorVec[iP];
                                //
								//for(int i=0;i<ColorChannels_c;i++){
								//	model_directProduction_Bands[iP][i] *=	polCorrFactorVec[iP];
								//	model_directProduction_Bands_errlow_1sig[iP][i] *=	polCorrFactorVec[iP];
								//	model_directProduction_Bands_errlow_2sig[iP][i] *=	polCorrFactorVec[iP];
								//	model_directProduction_Bands_errlow_3sig[iP][i] *=	polCorrFactorVec[iP];
								//	model_directProduction_Bands_errhigh_1sig[iP][i] *=	polCorrFactorVec[iP];
								//	model_directProduction_Bands_errhigh_2sig[iP][i] *=	polCorrFactorVec[iP];
								//	model_directProduction_Bands_errhigh_3sig[iP][i] *=	polCorrFactorVec[iP];
								//}
								//for(int i=0;i<StatesCont_c+1;i++){
								//	model_FeedDown[iP][i] *=	polCorrFactorVec[iP];
								//}

								d_data_centralval[iP] *= lumiCorrFactorVec[iP];//lumiCorrect the data (according to MPV of corresponding Np)
								d_data_errlow_centralval[iP] *=  	lumiCorrFactorVec[iP];
								d_data_errhigh_centralval[iP] *= 	lumiCorrFactorVec[iP];
								d_model_centralval[iP] *=        	lumiCorrFactorVec[iP];//un-lumiCorrect the model, as it was corrected in the fit (getCorrPromptProduction...)
								d_model_errlow_centralval[iP] *= 	lumiCorrFactorVec[iP];
								d_model_errhigh_centralval[iP] *=	lumiCorrFactorVec[iP];
								d_model_errlow_absolute_centralval[iP] *= 	lumiCorrFactorVec[iP];
								d_model_errhigh_absolute_centralval[iP] *=	lumiCorrFactorVec[iP];

								d_model_errlow_centralval_1Sig[iP] *= 	lumiCorrFactorVec[iP];
								d_model_errhigh_centralval_1Sig[iP] *=	lumiCorrFactorVec[iP];
								d_model_errlow_absolute_centralval_1Sig[iP] *= 	lumiCorrFactorVec[iP];
								d_model_errhigh_absolute_centralval_1Sig[iP] *=	lumiCorrFactorVec[iP];

								d_model_errlow_centralval_2Sig[iP] *= 	lumiCorrFactorVec[iP];
								d_model_errhigh_centralval_2Sig[iP] *=	lumiCorrFactorVec[iP];
								d_model_errlow_absolute_centralval_2Sig[iP] *= 	lumiCorrFactorVec[iP];
								d_model_errhigh_absolute_centralval_2Sig[iP] *=	lumiCorrFactorVec[iP];

								d_model_errlow_centralval_3Sig[iP] *= 	lumiCorrFactorVec[iP];
								d_model_errhigh_centralval_3Sig[iP] *=	lumiCorrFactorVec[iP];
								d_model_errlow_absolute_centralval_3Sig[iP] *= 	lumiCorrFactorVec[iP];
								d_model_errhigh_absolute_centralval_3Sig[iP] *=	lumiCorrFactorVec[iP];


							}

						}





						bool plotHorizontalErrorBars=false;

						TGraphAsymmErrors *data_Graph;
						data_Graph = new TGraphAsymmErrors(nPtBinsSel,d_data_pTmean,d_data_centralval,d_data_errlow_pT,d_data_errhigh_pT,d_data_errlow_centralval,d_data_errhigh_centralval);
						if(!plotHorizontalErrorBars) data_Graph = new TGraphAsymmErrors(nPtBinsSel,d_data_pTmean,d_data_centralval,d_zero,d_zero,d_data_errlow_centralval,d_data_errhigh_centralval);
						TGraph *model_Graph = new TGraph(nPtBinsSel,d_data_pTmean,d_model_centralval);
						TGraph *model_Graph_low = new TGraph(nPtBinsSel,d_data_pTmean,d_model_errlow_absolute_centralval);
						TGraph *model_Graph_high = new TGraph(nPtBinsSel,d_data_pTmean,d_model_errhigh_absolute_centralval);
						TGraph *data_Graph_nopolcorr = new TGraph(nPtBinsSel,d_data_pTmean,d_data_centralval_nopolcorr);

						TGraph *model_Graph_low_Bands[3];
						TGraph *model_Graph_high_Bands[3];

						model_Graph_low_Bands[0] = new TGraph(nPtBinsSel,d_data_pTmean,d_model_errlow_absolute_centralval_1Sig);
						model_Graph_high_Bands[0] = new TGraph(nPtBinsSel,d_data_pTmean,d_model_errhigh_absolute_centralval_1Sig);
						model_Graph_low_Bands[1] = new TGraph(nPtBinsSel,d_data_pTmean,d_model_errlow_absolute_centralval_2Sig);
						model_Graph_high_Bands[1] = new TGraph(nPtBinsSel,d_data_pTmean,d_model_errhigh_absolute_centralval_2Sig);
						model_Graph_low_Bands[2] = new TGraph(nPtBinsSel,d_data_pTmean,d_model_errlow_absolute_centralval_3Sig);
						model_Graph_high_Bands[2] = new TGraph(nPtBinsSel,d_data_pTmean,d_model_errhigh_absolute_centralval_3Sig);

						TGraphAsymmErrors *model_Graph_Bands[3];
						model_Graph_Bands[0] = new TGraphAsymmErrors(nPtBinsSel,d_data_pTmean,d_model_centralval,d_data_errlow_pT,d_data_errhigh_pT,d_model_errlow_centralval_1Sig,d_model_errhigh_centralval_1Sig);
						model_Graph_Bands[1] = new TGraphAsymmErrors(nPtBinsSel,d_data_pTmean,d_model_centralval,d_data_errlow_pT,d_data_errhigh_pT,d_model_errlow_centralval_2Sig,d_model_errhigh_centralval_2Sig);
						model_Graph_Bands[2] = new TGraphAsymmErrors(nPtBinsSel,d_data_pTmean,d_model_centralval,d_data_errlow_pT,d_data_errhigh_pT,d_model_errlow_centralval_3Sig,d_model_errhigh_centralval_3Sig);

						TGraph *model_CS_THunc[3];
						model_CS_THunc[0] = new TGraph(nPtBinsSel,d_data_pTmean,d_model_CS_THunc_m1);
						model_CS_THunc[1] = new TGraph(nPtBinsSel,d_data_pTmean,d_model_CS_THunc_central);
						model_CS_THunc[2] = new TGraph(nPtBinsSel,d_data_pTmean,d_model_CS_THunc_p1);

						TGraphAsymmErrors *model_CS_THuncBand = new TGraphAsymmErrors(nPtBinsSel,d_data_pTmean,d_model_CS_THunc_delta_p1_m1_mean,d_data_errlow_pT,d_data_errhigh_pT,d_model_CS_THunc_delta_p1_m1_diffhalf,d_model_CS_THunc_delta_p1_m1_diffhalf);

						cout<<"Model Graph for MeasurementID "<<iMeasurementID<<", Experiment "<<iExperiment<<", State "<<iState<<endl;
						model_Graph->Print();
						model_Graph_low->Print();
						model_Graph_high->Print();

						//Make TGraphs for model contributions

						TGraph *model_Graph_FeedDowns[StatesCont_c+1];//fist element: inclusive FeedDown

						//cout<<model_FeedDown<<endl;
						//cout<<"StatesCont_c "<<StatesCont_c<<endl;
						for(int i=0;i<StatesCont_c+1;i++){
							for(int iP = 0; iP < nPtBinsSel; iP++){
								d_model_centralval[iP] =  model_FeedDown[iP][i];
								if(iMeasurementID==0) d_model_centralval[iP];
							}
							model_Graph_FeedDowns[i] = new TGraphAsymmErrors(nPtBinsSel,d_data_pTmean,d_model_centralval);

							cout<<"model_Graph_FeedDowns["<<i<<"]:"<<endl;
							model_Graph_FeedDowns[i]->Print();
						}


						//PAPER: decide here if you want to use the MPV of the color channels, or the prediction from the Op_MPVs, add uncertainties (TGraphAsymmErrors)
						TGraph *model_Graph_ColorChannels[ColorChannels_c];

						for(int i=0;i<ColorChannels_c;i++){
							for(int iP = 0; iP < nPtBinsSel; iP++){
								//d_model_centralval[iP] =  model_directProduction[iP][i];
								d_model_centralval[iP] =  model_directProduction_Bands[iP][i];
								//if(iMeasurementID==0) d_model_centralval[iP];
							}
							model_Graph_ColorChannels[i] = new TGraphAsymmErrors(nPtBinsSel,d_data_pTmean,d_model_centralval);
							cout<<"graph of CC"<<i<<endl;
							model_Graph_ColorChannels[i]->Print();
						}


						TGraphAsymmErrors *model_Graph_ColorChannels_Bands[3][500];//[ColorChannels_c];
						double d_model_centralval_Bands[nPtBinsSel_];
						double d_model_errlow_Bands[nPtBinsSel_];
						double d_model_errhigh_Bands[nPtBinsSel_];

						cout<<"model_directProduction_Bands"<<endl;
						cout<<model_directProduction_Bands<<endl;

						for(int iSig=0;iSig<3;iSig++){
							for(int i=0;i<ColorChannels_c;i++){
								for(int iP = 0; iP < nPtBinsSel; iP++){
									d_model_centralval_Bands[iP] =  model_directProduction_Bands[iP][i];
									if(iSig==0) d_model_errlow_Bands[iP] =   model_directProduction_Bands_errlow_1sig[iP][i];
									if(iSig==1) d_model_errlow_Bands[iP] =   model_directProduction_Bands_errlow_2sig[iP][i];
									if(iSig==2) d_model_errlow_Bands[iP] =   model_directProduction_Bands_errlow_3sig[iP][i];
									if(iSig==0) d_model_errhigh_Bands[iP] =  model_directProduction_Bands_errhigh_1sig[iP][i];
									if(iSig==1) d_model_errhigh_Bands[iP] =  model_directProduction_Bands_errhigh_2sig[iP][i];
									if(iSig==2) d_model_errhigh_Bands[iP] =  model_directProduction_Bands_errhigh_3sig[iP][i];
								}
								model_Graph_ColorChannels_Bands[iSig][i] = new TGraphAsymmErrors(nPtBinsSel,d_data_pTmean,d_model_centralval_Bands,d_data_errlow_pT,d_data_errhigh_pT,d_model_errlow_Bands,d_model_errhigh_Bands);
								cout<<"model_Graph_ColorChannels_Bands CC: "<<i<<", CL:"<<iSig+1<<endl;
								model_Graph_ColorChannels_Bands[iSig][i]->Print();
								if(iSig==0&&i==0&&iMeasurementID==0) cout<<"CHECK_IT_OUT"<<endl;
							}
						}




						bool checkSums=true;
						if(checkSums){

						double checkTotal;
						double checkDirectColors[ColorChannels_c];
						double checkFDconts[StatesCont_c+1];

						double checkFD=0;
						double checkDirect=0;


						double checkBuffx, checkBuffy;


						for(int iP = 0; iP < nPtBinsSel; iP++){

							checkFD=0;
							checkDirect=0;

							model_Graph->GetPoint(iP, checkBuffx, checkBuffy);
							checkTotal=checkBuffy;
							model_Graph->GetPoint(iP, checkBuffx, checkBuffy);
							checkTotal=checkBuffy;

							cout<<"checkTotal = "<<checkTotal<<endl;

							for(int i=0;i<ColorChannels_c;i++){
								model_Graph_ColorChannels[i]->GetPoint(iP, checkBuffx, checkBuffy);
								checkDirectColors[i]=checkBuffy;
								checkDirect+=checkDirectColors[i];
								cout<<"checkDirectColors["<<i<<"] = "<<checkDirectColors[i]<<endl;
							}
							cout<<"checkDirect = "<<checkDirect<<endl;
							for(int i=0;i<StatesCont_c+1;i++){
								model_Graph_FeedDowns[i]->GetPoint(iP, checkBuffx, checkBuffy);
								checkFDconts[i]=checkBuffy;
								if(i>0&&i<StatesCont_c) checkFD+=checkFDconts[i];
								cout<<"checkFDconts["<<i<<"] = "<<checkFDconts[i]<<endl;
							}
							cout<<"checkFD = "<<checkFD<<endl;

							double directRatio=(checkFDconts[StatesCont_c])/checkDirect;
							double FDRatio=(checkFDconts[0])/checkFD;
							double TotalRatio=(checkTotal)/(checkFDconts[0]+checkFDconts[StatesCont_c]);

							cout<<"pT data = "<<checkBuffx<<endl;
							cout<<"Direct from FD (in plot) - sumColors = "<<(checkFDconts[StatesCont_c]-checkDirect)<<", ratio (plot wrt sum) = "<<(checkFDconts[StatesCont_c])/checkDirect <<endl;
							cout<<"FD inclusive (in plot) - sumFD = "<<(checkFDconts[0]-checkFD)<<", ratio (plot wrt sum) = "<<(checkFDconts[0])/checkFD  <<endl;
							cout<<"Total (in plot) - FD+direct (in plot) = "<<(checkTotal-(checkFDconts[0]+checkFDconts[StatesCont_c]))<<", ratio (plot wrt sum) = "<<(checkTotal)/(checkFDconts[0]+checkFDconts[StatesCont_c])  <<endl;

							if(TMath::Abs(directRatio-1)>0.25 || TMath::Abs(FDRatio-1)>0.25 || TMath::Abs(TotalRatio-1)>0.25)
								warnInconsistency[iState][iMeasurementID][iExperiment][iRap][iP]=true;

						}


						}



/*					//PAPER...:
 				    dmatrix model_directProduction_Bands;
				    dmatrix model_directProduction_Bands_errlow_1sig;
				    dmatrix model_directProduction_Bands_errhigh_1sig;
				    dmatrix model_directProduction_Bands_errlow_2sig;
				    dmatrix model_directProduction_Bands_errhigh_2sig;
				    dmatrix model_directProduction_Bands_errlow_3sig;
				    dmatrix model_directProduction_Bands_errhigh_3sig;

 */
						bool plotInclusiveFeedDown=true;
						if(StatesCont_c<3) plotInclusiveFeedDown=false;
						bool plotIndividualFeedDown=true;
						if(StatesCont_c<2) plotInclusiveFeedDown=false;
						bool plotDirectInclusive=false;
						if(StatesCont_c>1) plotDirectInclusive=true;
						bool plotDirectColorChannels=true;
						if(ColorChannels_c<2) plotDirectColorChannels=false;

						cout<<"StatesCont_c "<<StatesCont_c<<endl;

						char rapchar[200];
						char rapmiddlechar[200];

						if(isAbsRap) sprintf(rapmiddlechar,"|#it{y}|");
						else sprintf(rapmiddlechar,"#it{y}");

						sprintf(rapchar,"%1.1f < %s < %1.1f",rapMinObject, rapmiddlechar, rapMaxObject);
						if(TMath::Abs(rapMinObject-0.75)<1e-3) sprintf(rapchar,"%1.2f < %s < %1.1f",rapMinObject, rapmiddlechar, rapMaxObject);
						if(TMath::Abs(rapMaxObject-0.75)<1e-3) sprintf(rapchar,"%1.1f < %s < %1.2f",rapMinObject, rapmiddlechar, rapMaxObject);

						bool longrapchar=true;

						if(isAbsRap&&rapMinObject<1e-3){
							sprintf(rapchar,"%s < %1.1f", rapmiddlechar, rapMaxObject);
							if(TMath::Abs(rapMaxObject-0.75)<1e-3) sprintf(rapchar,"%s < %1.2f", rapmiddlechar, rapMaxObject);
							longrapchar=false;
						}

						bool smoothPol=false;
						for(int iSmoothPlot=0;iSmoothPlot<2;iSmoothPlot++){
							if(iMeasurementID!=1 && iSmoothPlot>0) continue;
							if(PredictionDataPlot && iSmoothPlot>0) continue;
							if(iSmoothPlot==0) smoothPol=false;
							else smoothPol=true;
							plotComp( iState,  iMeasurementID,  iExperiment,  iRap,  jobdirname, rapchar, data_Graph, model_Graph, model_Graph_low, model_Graph_high, plotInclusiveFeedDown, plotIndividualFeedDown, plotDirectColorChannels, StatesCont_c, ColorChannels_c, model_Graph_FeedDowns, model_Graph_ColorChannels, data_Graph_nopolcorr, StatesContributing_, chi2Min, chi2Prob, ndf, HPbool, pTMinModel, pTMaxModel, longrapchar, model_Graph_ColorChannels_Bands, model_Graph_low_Bands, model_Graph_high_Bands, model_Graph_Bands, model_CS_THunc, model_CS_THuncBand, plotDirectInclusive, nDataPointsInRange, PredictionPlot, PredictionDataPlot, smoothPol, PlotVSpToverM);

						}

				}


				}//rap loop



			}
		}
	}




	for(int iState=0; iState < nStates; iState++){
		for(int iMeasurementID=0; iMeasurementID < NRQCDvars::nMeasurementIDs; iMeasurementID++){
			for(int iExperiment=0; iExperiment < NRQCDvars::nExperiments; iExperiment++){
				for(int iRap = 0; iRap < NRQCDvars::nMaxRapBins; iRap++){
				    for(int iP = 0; iP < NRQCDvars::nMaxPtBins; iP++){

				    	if(warnInconsistency[iState][iMeasurementID][iExperiment][iRap][iP]){
				    		cout<<" "<<endl;
				    		cout<<"-------------------------------------------------------------------------------------------------------------------------------------------------"<<endl;
				    		cout<<"WARNING in "<<StateName[iState]<<", iMeasurementID="<<MeasurementIDName[iMeasurementID]<<", iExperiment="<<ExpName[iExperiment]<<", iRap="<<iRap<<", iP="<<iP<<" (sum of individual components does not add up -> check output)"<<endl;
				    		cout<<"-------------------------------------------------------------------------------------------------------------------------------------------------"<<endl;
				    	}


				    }
				}
			}
		}
	}
	cout<<" "<<endl;





		return 0;

  	}



void plotComp(int iState, int iMeasurementID, int iExperiment, int iRap, char jobdirname[200], char rapchar[200], TGraphAsymmErrors *data_Graph, TGraph *model_Graph, TGraph *model_Graph_low, TGraph *model_Graph_high, bool plotInclusiveFeedDown, bool plotIndividualFeedDown, bool plotDirectColorChannels, const int StatesCont_c, const int ColorChannels_c, TGraph *model_Graph_FeedDowns[500], TGraph *model_Graph_ColorChannels[500], /*TGraphAsymmErrors *model_Graph_ColorChannels[3][500],*/ TGraph *data_Graph_nopolcorr, vector<int> StatesContributing, double chi2Min, double chi2Prob, int ndf, bool HPbool, double pTMinModel, double pTMaxModel, bool longrapchar, TGraphAsymmErrors *model_Graph_ColorChannels_Bands[3][500], TGraph *model_Graph_low_Bands[3], TGraph *model_Graph_high_Bands[3], TGraphAsymmErrors *model_Graph_Bands[3], TGraph *model_CS_THunc[3], TGraphAsymmErrors *model_CS_THuncBand, bool plotDirectInclusive, int nDataPointsInRange, bool PredictionPlot, bool PredictionDataPlot, bool smoothPol, bool PlotVSpToverM){
	gROOT->Reset();
	gStyle->SetOptStat(11);
	gStyle->SetOptFit(101);
	gStyle->SetFillColor(kWhite);

	gStyle->SetPadBottomMargin(0.13);
	gStyle->SetPadLeftMargin(0.12);
	gStyle->SetPadRightMargin(0.02);
	gStyle->SetPadTopMargin(0.02);

	bool narrowCanvas=false;
	if(PredictionPlot) narrowCanvas=true;
	double narrFac=0.8;
	if(narrowCanvas){
		gStyle->SetPadBottomMargin(0.13);
		gStyle->SetPadLeftMargin(0.15);
		gStyle->SetPadRightMargin(0.02);
		gStyle->SetPadTopMargin(0.02);
	}

	bool isSstate=(StateQuantumID[iState] > NRQCDvars::quID_S)?false:true;

	bool singlePointModel=false;
	if(nDataPointsInRange<2) singlePointModel=true;
	cout<<"nDataPointsInRange = "<<nDataPointsInRange<<endl;

	//0
	//bool plotColorChannelBands[4]={true, true, true,     false};
	//bool plotColorChannel[4]={true, true, true,     false};
	//bool plotTotalModelBands=true;

	//1
	//bool plotColorChannelBands[4]={true, false, false,     false};
	//bool plotTotalModelBands=true;

	//2
	//bool plotColorChannelBands[4]={false, true, false,     false};
	//bool plotTotalModelBands=false;

	//3
	//bool plotColorChannelBands[4]={false, false, true,     false};
	//bool plotTotalModelBands=false;

	//000

	bool plotData=true;
	bool plotModel=true;

	bool plotColorChannelBands[4]={false, true, true,     true};
	bool plotColorChannel[4]={true, true, true,     true};
	bool plotTotalModelBands=false;
	bool plotTotalModel=true;
	bool plotCSTHuncLines=true;
	bool plotCSTHuncBand=false;//always false (not fully implemented)

	if(!isSstate){
		plotCSTHuncLines=false;
		plotCSTHuncBand=false;
		plotColorChannel[0]=false;
		plotColorChannelBands[0]=false;
	}

	if(PredictionPlot){
		plotCSTHuncLines=false;
		plotCSTHuncBand=false;
		plotColorChannel[0]=true;
		plotColorChannel[1]=true;
		plotColorChannel[2]=true;
		plotColorChannelBands[0]=false;
		plotColorChannelBands[1]=false;
		plotColorChannelBands[2]=false;
		plotTotalModel=true;
		plotTotalModelBands=true;
		plotData=false;
		if(PredictionDataPlot){
			plotCSTHuncLines=false;
			plotCSTHuncBand=false;
			plotColorChannel[0]=false;
			plotColorChannel[1]=false;
			plotColorChannel[2]=false;
			plotColorChannelBands[0]=false;
			plotColorChannelBands[1]=false;
			plotColorChannelBands[2]=false;
			plotTotalModel=false;
			plotTotalModelBands=false;
			plotModel=false;
			plotData=true;
		}
	}

	//if(iMeasurementID==1){
	//	plotColorChannelBands[0]=true;
	//	plotColorChannelBands[2]=true;
	//}

	///////////////////////

	const int nTHunc=3;
	bool plotSigmaBands[3]={true, true, true};
	int ColorChannelSequence[4]={0, 2, 1,     3};
	bool plotCSTHuncLine[nTHunc]={true, true, true};

	char plotOption_ColorChannelCentrals[200];
	sprintf(plotOption_ColorChannelCentrals,"lsame");
	if(iMeasurementID!=0) sprintf(plotOption_ColorChannelCentrals,"lsame");

	char plotOptionBands[200];
	sprintf(plotOptionBands,"3same");
	if(iMeasurementID!=0) sprintf(plotOptionBands,"E3same");


	TCanvas *plotCanvas = new TCanvas("plotCanvas","plotCanvas",1200,800);


	if(narrowCanvas) plotCanvas = new TCanvas("plotCanvas","plotCanvas",800,800);

	plotCanvas->SetFillColor(kWhite);
	plotCanvas->SetFrameBorderMode(0);


	double x_min=0.;
	double x_max=80;
	double y_min;
	double y_max;

	if(iState==0){
		y_min=1.001e-5; y_max =4.999e3;
		x_min=0.001; x_max=89.999;
	}
	if(iState==1){
		y_min=1.001e-3; y_max =0.999e1;
		x_min=10.001; x_max=44.999;
	}
	if(iState==2){
		y_min=1.001e-3; y_max =0.999e1;
		x_min=10.001; x_max=44.999;
	}
	if(iState==3){
		y_min=2.001e-3; y_max =2.999e2;
		x_min=0.001; x_max=49.999;
		if(PredictionPlot){
			y_min=1.001e-6; y_max =9.999e0;
			x_min=10.001; x_max=99.999;
			y_min=5.001e-9; y_max =9.999e0;
		}

	}
	if(iState==10){
		y_min=1.001e-6; y_max=1.4999;
		x_min=0.001; x_max =134.999;
		if(PredictionPlot){
			y_min=5.001e-9; y_max =4.999e-2;
			x_min=20.001; x_max =259.999;
			y_min=5.001e-9; y_max =9.999e0;
		}
	}

	//LHCb specials
	if(iExperiment==3 || iExperiment==4 || iExperiment==5){
		if(iState==0 || iState==3){
			x_min=0.001; x_max =24.999;
			y_min=2.001e-3; y_max =2.999e2;
		}
		if(iState==10){
			x_min=8.001; x_max =21.999;
			y_min=5.001e-3; y_max=0.999;
		}
	}

	//Polarization specials
	if(iMeasurementID==1) {
		y_min=-1.399; y_max =1.399;
		if(iState==0){
			x_min=10.001; x_max=99.999;
		}
		if(iState==3){
			x_min=8.001; x_max =71.999;
			if(PredictionPlot){
				x_min=10.001; x_max =99.999;
			}
		}
		if(iState==10){
			x_min=8.001; x_max =71.999;
			if(PredictionPlot){
				x_min=20.001; x_max =259.999;
			}
		}
	}
	if(iMeasurementID==2) { y_min=-0.699; y_max =0.699; }
	if(iMeasurementID==3) { y_min=-0.699; y_max =0.699; }

	if(PlotVSpToverM){
		x_min/=NRQCDvars::mass[iState]; x_max/=NRQCDvars::mass[iState];
		if(PredictionPlot || PredictionDataPlot){
			x_min=2.001; x_max=27.9999;
		}
	}



	bool logY=true;
	if(iMeasurementID!=0) logY=false;

	TH1F *axishist = new TH1F;

	axishist = new TH1F("axishist", "axishist", 100, x_min, x_max);


	axishist->SetTitle(0);
	axishist = plotCanvas->DrawFrame(x_min,y_min,x_max,y_max);
	axishist->GetXaxis()->SetTitle("#it{p}_{T} [GeV]");
	if(PlotVSpToverM) axishist->GetXaxis()->SetTitle("#it{p}_{T} / M");
	if(iMeasurementID==0) axishist->GetYaxis()->SetTitle(Form("%s",NRQCDvars::MeasurementIDNameTex[iMeasurementID]));
	else axishist->GetYaxis()->SetTitle("");
	axishist->GetYaxis()->SetTitleSize(0.05);
	axishist->GetXaxis()->SetTitleSize(0.05);
	axishist->GetYaxis()->SetTitleOffset(1.05);
	axishist->GetXaxis()->SetTitleOffset(1.2);
	axishist->GetXaxis()->SetTickLength(-0.02);
	axishist->GetYaxis()->SetTickLength(-0.02);
	//axishist->GetXaxis()->SetTicks("-");
	//axishist->GetYaxis()->SetTicks("-");
	//axishist->GetYaxis()->SetRangeUser(y_min, y_max);
	axishist->GetYaxis()->SetLabelOffset(0.025);
	axishist->GetXaxis()->SetLabelOffset(0.025);

	if(narrowCanvas){
		axishist->GetYaxis()->SetTitleSize(0.05*narrFac);
		axishist->GetXaxis()->SetTitleSize(0.05*narrFac);
		axishist->GetYaxis()->SetTitleOffset(1.775);
		axishist->GetXaxis()->SetTitleOffset(1.5);
	}
	//axishist->Draw();

	//char xaxisOption[200];
	//char yaxisOption[200];
	//sprintf(xaxisOption,"-US");
	//sprintf(yaxisOption,"+US");
	//if(logY) sprintf(yaxisOption,"+GUS");
    //
	//int axisDivisionsX=511;
	//int axisDivisionsY=511;
	//if(logY) axisDivisionsY=50510;
    //
    //TGaxis *axisy = new TGaxis(x_min,y_min,x_min,y_max,y_min,y_max,axisDivisionsY,yaxisOption);
    //axisy->SetTickSize(0.02);
    //axisy->Draw("same");
    //TGaxis *axisx = new TGaxis(x_min,y_min,x_max,y_min,x_min,x_max,axisDivisionsX,xaxisOption);
    //axisx->SetTickSize(0.02);
    //axisx->Draw("same");


	//Remove model below bin-center of first data bin:
	double errbuffx, buffx, buffy;
	for(int m=0;m<data_Graph->GetN();m++){
		data_Graph->GetPoint(m,buffx, buffy);
		errbuffx=data_Graph->GetErrorXlow(m);
		if(buffx-errbuffx>pTMinModel-1e-3)
			{pTMinModel=buffx*0.9999;
			break;}
		if(m==data_Graph->GetN()-1) pTMinModel=999;
	}
	for(int m=data_Graph->GetN()-1;m>-1;m--){
		data_Graph->GetPoint(m,buffx, buffy);
		errbuffx=data_Graph->GetErrorXlow(m);
		if(buffx+errbuffx<pTMaxModel+1e-3)
			{pTMaxModel=buffx*1.0001;
			break;}
		if(m==0) pTMaxModel=0;
	}

	cout<<"final pTMinModel = "<<pTMinModel<<endl;
	cout<<"final pTMaxModel = "<<pTMaxModel<<endl;

	if(iState==10&&iMeasurementID==1) pTMaxModel=50.;

	/////




	if(singlePointModel){
		cout<<"setting smoothPol=false for state "<<iState<<endl;
		smoothPol=false;
	}

	//smoothPol=false;

	if(smoothPol){


		cout<<"smoothPol for state "<<iState<<endl;

		int monotonicInt=0;
		int c_sign=0;
		int polOrder=0;

		for(int iColorChannel=0;iColorChannel<ColorChannels_c;iColorChannel++){
			if(iColorChannel==1) continue;
			monotonicInt=0;
			if(iColorChannel==0) {monotonicInt=-1; c_sign=+1;}
			if(iColorChannel==2) {monotonicInt=+1; c_sign=-1;}

			//cout<<"lll "<<iState<<endl;
			//model_Graph_ColorChannels[iColorChannel]->Print();
			model_Graph_ColorChannels[iColorChannel]=smoothPolFunc(model_Graph_ColorChannels[iColorChannel], monotonicInt, c_sign);
			//model_Graph_ColorChannels[iColorChannel]->Print();


			for(int iSig=0;iSig<3;iSig++){
				model_Graph_ColorChannels_Bands[iSig][iColorChannel] = CloneCentralValueFromGraph(model_Graph_ColorChannels_Bands[iSig][iColorChannel], model_Graph_ColorChannels[iColorChannel]);
				polOrder=1;
				model_Graph_ColorChannels_Bands[iSig][iColorChannel] = smoothPolFuncErrors(model_Graph_ColorChannels_Bands[iSig][iColorChannel], polOrder);

				if(iState==10){//extrapolate the ups3S polarization S1 uncertainty bands to 50 GeV, for the first paper
					double buffx, buffy;
					int buffN=model_Graph_ColorChannels_Bands[iSig][iColorChannel]->GetN();
					model_Graph_ColorChannels_Bands[iSig][iColorChannel]->GetPoint(buffN-1, buffx, buffy);
					model_Graph_ColorChannels_Bands[iSig][iColorChannel]->SetPoint(buffN-1, 50., buffy);
					//cout<<"lll "<<iSig<<endl;
					//model_Graph_ColorChannels_Bands[iSig][iColorChannel]->Print();
				}
			}

		}

		for(int i=0;i<StatesCont_c+1;i++){
			monotonicInt=0;
			c_sign=0;
			model_Graph_FeedDowns[i]=smoothPolFunc(model_Graph_FeedDowns[i], monotonicInt, c_sign);
		}

		monotonicInt=0;
		if(PredictionPlot&&iState==10) monotonicInt=+1;
		c_sign=0;
		model_Graph=smoothPolFunc(model_Graph, monotonicInt, c_sign);
		for(int iSig=0;iSig<3;iSig++){
			model_Graph_Bands[iSig] = CloneCentralValueFromGraph(model_Graph_Bands[iSig], model_Graph);
			polOrder=2;
			model_Graph_Bands[iSig] = smoothPolFuncErrors(model_Graph_Bands[iSig], polOrder);
		}

		for(int iTHunc=0;iTHunc<nTHunc;iTHunc++){
			monotonicInt=-1;
			c_sign=+1;
			model_CS_THunc[iTHunc]=smoothPolFunc(model_CS_THunc[iTHunc], monotonicInt, c_sign);
		}



	}


	cout<<"after smoothPol"<<endl;



	int colorData=1;
	int colorModel=600;
	int colorData_nopolcorr=418;
	double linewidthModel=2;
	double linewidthModelErrors=1;

	int colorCC[6]={435, 1, 1, 1, 1, 1};//{632, 418, 616, 800, 0, 0};
	int colorCC_neg[6]={632, 632, 632, 632, 632, 632};

	int colorFD[NRQCDvars::nStates]={0, 616, 800, 632, 0, 0, 0, 0, 0, 0, 0, 0};
	int colorFDDirect=810;
	int colorFDInclusive=434;
	int linestyleFDDirect=7;
	int linestyleFDInclusive=7;
	double linewidthFDDirect=1.;
	double linewidthFDInclusive=1.;

	int linestyleCC[6]={4, 5, 7, 9, 10, 6};
	int linestyleCC_neg[6]={4, 5, 7, 9, 10, 6};
	int linestyleFD=3;
	int MarkerStyle[6]={21, 24, 25, 26, 27, 28};//{632, 418, 616, 800, 0, 0};
	double linewidthCC=2;
	double linewidthFD=1;

	int styleModel[3]={2, 4, 3};

	//int FillStyle_sig[6]={1001, 1001, 1001, 1, 1, 1};
	int FillStyle_sig[6]={1001, 1001, 1001, 1001, 1, 1};
	int LineStyle_Central[6]={1, 1, 1, 1, 1, 1};
	int LineColor_Central[6]={923, 416+2,632+2, 616+2, 1, 1};
	int LineColor_nSig[3][6]={
			{16, 416-3,632+0, 616, 1, 1},
			{17, 416-7,632-7, 616-7, 1, 1},
			{18, 416-10,632-10, 616-10, 1, 1}
	};

	int TotalModel_FillStyle_sig=1001;
	int TotalModel_LineStyle_Central=1;
	int TotalModel_LineColor_Central=600+0;
	int TotalModel_LineColor_nSig[3]={600-7, 600-9, 600-10};

	//TotalModel_LineColor_Central=14;
	//TotalModel_LineColor_nSig[0]=16;
	//TotalModel_LineColor_nSig[1]=17;
	//TotalModel_LineColor_nSig[2]=18;


	int MarkerStyleExp[9]={20, 20, 20, 33, 33, 33, 25, 25, 25};
	int MarkerStylePolRap[3]={20, 24, 31};
	double MarkerSizeExp[9]={1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.25, 1.25, 1.25};

	double LineStyle_THunc[nTHunc]={2, 4, 3};
	double LineColor_THunc[nTHunc]={923, 923, 923};
	double linewidthTHunc=1;

	int FillStyle_THunc=1001;
	int FillColor_THunc=18;


	if(HPbool){

		colorData=1;
		colorModel=600;
		linewidthModel=2;
		linewidthModelErrors=1;
		linewidthCC=2;


		colorCC[0]=923;
		colorCC[1]=418;
		colorCC[2]=634;
		colorCC[3]=618;
		colorCC[4]=1;
		colorCC[5]=1;

		colorCC_neg[0]=923;
		colorCC_neg[1]=418;
		colorCC_neg[2]=634;
		colorCC_neg[3]=435;
		colorCC_neg[4]=618;
		colorCC_neg[5]=1;

		linestyleCC[0]=1;
		linestyleCC[1]=1;
		linestyleCC[2]=1;
		linestyleCC[3]=1;
		linestyleCC[4]=1;
		linestyleCC[5]=1;
		linestyleCC_neg[0]=2;
		linestyleCC_neg[1]=2;
		linestyleCC_neg[2]=2;
		linestyleCC_neg[3]=2;
		linestyleCC_neg[4]=2;
		linestyleCC_neg[5]=2;

	}

	cout<<"start setting up cosmetics"<<endl;

	data_Graph->SetMarkerStyle(MarkerStyleExp[iExperiment]);
	data_Graph->SetMarkerSize(MarkerSizeExp[iExperiment]);
	data_Graph->SetMarkerColor(colorData);
	data_Graph->SetLineColor(colorData);

	if(iMeasurementID>0){
		data_Graph->SetMarkerStyle(MarkerStylePolRap[iRap]);
	}

	data_Graph_nopolcorr->SetMarkerStyle(24);
	data_Graph_nopolcorr->SetMarkerSize(1.);
	data_Graph_nopolcorr->SetMarkerColor(colorData_nopolcorr);
	data_Graph_nopolcorr->SetLineColor(colorData_nopolcorr);

	model_Graph->SetLineColor(colorModel);
	model_Graph->SetLineStyle(1);
	model_Graph->SetLineWidth(linewidthModel);
	model_Graph->SetFillColor(colorModel);

	model_Graph->SetMarkerStyle(MarkerStyle[0]);
	model_Graph->SetMarkerSize(1.5);
	model_Graph->SetMarkerColor(colorModel);

	//cout<<"Model Graph for MeasurementID "<<iMeasurementID<<", Experiment "<<iExperiment<<", State "<<iState<<endl;
	//model_Graph->Print();

	model_Graph_low->SetLineColor(colorModel);
	model_Graph_low->SetLineStyle(styleModel[0]);
	model_Graph_low->SetLineWidth(linewidthModelErrors);
	model_Graph_low->SetFillColor(colorModel);

	model_Graph_high->SetLineColor(colorModel);
	model_Graph_high->SetLineStyle(styleModel[0]);
	model_Graph_high->SetLineWidth(linewidthModelErrors);
	model_Graph_high->SetFillColor(colorModel);

	for(int iSig=0;iSig<3;iSig++){
		model_Graph_low_Bands[iSig]->SetLineColor(colorModel);
		model_Graph_low_Bands[iSig]->SetLineStyle(styleModel[iSig]);
		model_Graph_low_Bands[iSig]->SetLineWidth(linewidthModelErrors);
		model_Graph_low_Bands[iSig]->SetFillColor(colorModel);
		model_Graph_high_Bands[iSig]->SetLineColor(colorModel);
		model_Graph_high_Bands[iSig]->SetLineStyle(styleModel[iSig]);
		model_Graph_high_Bands[iSig]->SetLineWidth(linewidthModelErrors);
		model_Graph_high_Bands[iSig]->SetFillColor(colorModel);
	}

	model_Graph->SetLineColor(TotalModel_LineColor_Central);
	model_Graph->SetLineStyle(TotalModel_LineStyle_Central);
	model_Graph->SetLineWidth(linewidthModel);
	model_Graph->SetFillColor(TotalModel_LineColor_Central);

	for(int iSig=0;iSig<3;iSig++){
		model_Graph_Bands[iSig]->SetFillStyle(TotalModel_FillStyle_sig);
		model_Graph_Bands[iSig]->SetFillColor(TotalModel_LineColor_nSig[iSig]);
		model_Graph_Bands[iSig]->SetLineColor(TotalModel_LineColor_nSig[iSig]);
	}

	cout<<"start setting up cosmetics for CCs"<<endl;
	for(int iColorChannelSequence=0;iColorChannelSequence<ColorChannels_c;iColorChannelSequence++){
		int i=ColorChannelSequence[iColorChannelSequence];
		for(int iSig=0;iSig<3;iSig++){

			model_Graph_ColorChannels[i]->SetLineStyle(LineStyle_Central[i]);
			model_Graph_ColorChannels[i]->SetLineColor(LineColor_Central[i]);
			model_Graph_ColorChannels[i]->SetLineWidth(linewidthCC);

			model_Graph_ColorChannels_Bands[iSig][i]->SetFillStyle(FillStyle_sig[i]);
			model_Graph_ColorChannels_Bands[iSig][i]->SetFillColor(LineColor_nSig[iSig][i]);
			model_Graph_ColorChannels_Bands[iSig][i]->SetLineColor(LineColor_nSig[iSig][i]);

		}
	}




	cout<<"start setting up cosmetics for THuncert"<<endl;

	for(int iTHunc=0;iTHunc<nTHunc;iTHunc++){
		model_CS_THunc[iTHunc]->SetLineStyle(LineStyle_THunc[iTHunc]);
		model_CS_THunc[iTHunc]->SetLineColor(LineColor_THunc[iTHunc]);
		model_CS_THunc[iTHunc]->SetLineWidth(linewidthTHunc);
	}

	model_CS_THuncBand->SetFillStyle(FillStyle_THunc);
	model_CS_THuncBand->SetFillColor(FillColor_THunc);
	model_CS_THuncBand->SetLineColor(FillColor_THunc);

	cout<<"pTMinModel = "<<pTMinModel<<endl;
	cout<<"model_Graph->GetN() = "<<model_Graph->GetN()<<endl;

	//if(iExperiment==3 && pTMinModel>9.99){
	//	plotModel=false;
	//	plotCSTHuncLines=false;
	//	plotCSTHuncBand=false;
    //
	//}

	if(HPbool){
		int nGraph=model_Graph->GetN();
		int jTG=0;
		for(int j=0;j<nGraph;j++){
			double buffx, buffy;
			model_Graph->GetPoint(jTG, buffx, buffy);
			//cout<<"pTbuff "<<buffx<<endl;
			if(buffx<pTMinModel || buffx>pTMaxModel){
				model_Graph->RemovePoint(jTG);
				//cout<<"remove point "<<j<<endl;
				jTG--;
			}
			jTG++;
		}
		nGraph=model_Graph_low->GetN();
		jTG=0;
		for(int j=0;j<nGraph;j++){
			double buffx, buffy;
			model_Graph_low->GetPoint(jTG, buffx, buffy);
			if(buffx<pTMinModel || buffx>pTMaxModel){
				model_Graph_low->RemovePoint(jTG);
				model_Graph_high->RemovePoint(jTG);
				jTG--;
			}

			jTG++;
		}
		nGraph=model_Graph_low_Bands[0]->GetN();
		jTG=0;
		for(int j=0;j<nGraph;j++){
			double buffx, buffy;
			model_Graph_low_Bands[0]->GetPoint(jTG, buffx, buffy);
			if(buffx<pTMinModel || buffx>pTMaxModel){
				for(int iSig=0;iSig<3;iSig++){
					model_Graph_low_Bands[iSig]->RemovePoint(jTG);
					model_Graph_high_Bands[iSig]->RemovePoint(jTG);
				}
				jTG--;
			}

			jTG++;
		}
		nGraph=model_Graph_Bands[0]->GetN();
		jTG=0;
		for(int j=0;j<nGraph;j++){
			double buffx, buffy;
			model_Graph_Bands[0]->GetPoint(jTG, buffx, buffy);
			if(buffx<pTMinModel || buffx>pTMaxModel){
				for(int iSig=0;iSig<3;iSig++){
					model_Graph_Bands[iSig]->RemovePoint(jTG);
				}
				jTG--;
			}

			jTG++;
		}
		nGraph=model_CS_THunc[0]->GetN();
		jTG=0;
		for(int j=0;j<nGraph;j++){
			double buffx, buffy;
			model_CS_THunc[0]->GetPoint(jTG, buffx, buffy);
			if(buffx<pTMinModel || buffx>pTMaxModel){
				for(int iTHunc=0;iTHunc<nTHunc;iTHunc++){
					model_CS_THunc[iTHunc]->RemovePoint(jTG);
				}
				model_CS_THuncBand->RemovePoint(jTG);
				jTG--;
			}

			jTG++;
		}
		nGraph=model_Graph_FeedDowns[0]->GetN();
		jTG=0;
		for(int j=0;j<nGraph;j++){
			double buffx, buffy;
			model_Graph_FeedDowns[0]->GetPoint(jTG, buffx, buffy);
			if(buffx<pTMinModel || buffx>pTMaxModel){
				for(int i=0;i<StatesCont_c+1;i++){
					model_Graph_FeedDowns[i]->RemovePoint(jTG);
				}
				jTG--;
			}

			jTG++;
		}

	}



	int nLegendEntries=1;
	if(plotDirectColorChannels) nLegendEntries+=ColorChannels_c;
	if(plotIndividualFeedDown) nLegendEntries+=StatesCont_c;
	if(plotInclusiveFeedDown) nLegendEntries+=1;
	if(plotDirectInclusive) nLegendEntries+=1;

	if(iMeasurementID==0) nLegendEntries++;

	if(iMeasurementID==0&&HPbool) nLegendEntries--;


	double buff_pT;
	double buff_errpT_low;
	double buff_errpT_high;
	double buff_res;
	double buff_errres;



	if(model_Graph->GetN()<2){
		double buffModel_pT, buffModel_res;
		double buffData_pT, buffData_err_pTlow, buffData_err_pThigh, buffData_res;
		model_Graph->GetPoint(0, buffModel_pT, buffModel_res);
		cout<<"buffModel_pT "<<buffModel_pT<<endl;
		cout<<"buffModel_pT "<<buffModel_pT<<endl;

		for(int iData=0; iData<data_Graph->GetN();iData++){
			data_Graph->GetPoint(iData, buffData_pT, buffData_res);
			buffData_err_pTlow=data_Graph->GetErrorXlow(iData);
			buffData_err_pThigh=data_Graph->GetErrorXhigh(iData);
			cout<<"buffData_pT "<<buffData_pT<<endl;
			cout<<"buffData_err_pTlow "<<buffData_err_pTlow<<endl;
			cout<<"buffData_err_pThigh "<<buffData_err_pThigh<<endl;
			if(buffData_pT-buffData_err_pTlow<buffModel_pT && buffData_pT+buffData_err_pThigh>buffModel_pT){
				buff_pT=buffData_pT;
				buff_errpT_low=buffData_err_pTlow;
				buff_errpT_high=buffData_err_pThigh;
			}
		}

		cout<<"buff_pT "<<buff_pT<<endl;
		cout<<"buff_errpT_low "<<buff_errpT_low<<endl;
		cout<<"buff_errpT_high "<<buff_errpT_high<<endl;

	}


	if(iMeasurementID==0) if(!HPbool) data_Graph_nopolcorr->Draw("psame");



	double delta_y_legend=0.085*nLegendEntries;
	double max_y_legend;
	double min_x_legend=0.60;
	double max_x_legend=0.95;
	min_x_legend=0.7;
	if(iMeasurementID==0) max_y_legend=0.85;
	if(iMeasurementID>0 && iMeasurementID<4) max_y_legend=0.85;

	if(PredictionPlot) min_x_legend=0.8125;
	if(narrowCanvas){
		min_x_legend=0.845;
		delta_y_legend*=narrFac;
	}

	TLegend* legend;
	legend = new TLegend(min_x_legend,max_y_legend-delta_y_legend,max_x_legend,max_y_legend);
	if(iState==0) legend = new TLegend(min_x_legend,max_y_legend-delta_y_legend/1.55,max_x_legend,max_y_legend);
	if(iState==1) legend = new TLegend(min_x_legend,max_y_legend-delta_y_legend/1.25,max_x_legend,max_y_legend);
	if(iState==2) legend = new TLegend(min_x_legend,max_y_legend-delta_y_legend/1.25,max_x_legend,max_y_legend);
	legend->SetFillColor(0);
	legend->SetTextSize(0.04);
	legend->SetBorderSize(0);
	legend->SetFillStyle(0);

	if(narrowCanvas) legend->SetTextSize(0.04*narrFac);

	if(plotData) legend->AddEntry(data_Graph,Form("%s", ExpNameTex[iExperiment]),"lp");
	if(iMeasurementID==0) if(!HPbool) legend->AddEntry(data_Graph_nopolcorr,Form("%s %s (#vec{#lambda}=0)", ExpNameTex[iExperiment], StateNameTex[iState]),"p");

	if(plotModel){
		if(model_Graph->GetN()>1) legend->AddEntry(model_Graph,"Total","l");
		else legend->AddEntry(model_Graph,"Total","l");
	}

	TGraph* model_Graph_ColorChannels_neg[ColorChannels_c];
	TGraph* model_Graph_ColorChannels_pos[ColorChannels_c];
	TGraphAsymmErrors* model_Graph_ColorChannels_Bands_neg[3][ColorChannels_c];
	TGraphAsymmErrors* model_Graph_ColorChannels_Bands_pos[3][ColorChannels_c];

	bool PlotNegChannels[ColorChannels_c];
	bool PlotPosChannels[ColorChannels_c];
	if(plotDirectColorChannels&&plotModel){
		//cout<<"ColorChannels_c "<<ColorChannels_c<<endl;
		for(int iColorChannelSequence=0;iColorChannelSequence<ColorChannels_c;iColorChannelSequence++){
			int i=ColorChannelSequence[iColorChannelSequence];
			PlotNegChannels[i]=false;
			PlotPosChannels[i]=false;
			if(iMeasurementID!=0) PlotPosChannels[i]=true;

			//PAPER: change plotting of color channels: 3 error bands (could be tricky: pos/neg contributions)

			//model_Graph_ColorChannels[i]->Print();
			model_Graph_ColorChannels[i]->SetLineColor(colorCC[i]);
			model_Graph_ColorChannels[i]->SetLineStyle(linestyleCC[i]);
			model_Graph_ColorChannels[i]->SetLineWidth(linewidthCC);

			model_Graph_ColorChannels[i]->SetMarkerStyle(MarkerStyle[i+1]);
			model_Graph_ColorChannels[i]->SetMarkerSize(1.5);
			model_Graph_ColorChannels[i]->SetMarkerColor(colorCC[i]);



			if(HPbool){
				int nGraph=model_Graph_ColorChannels[i]->GetN();
				int jTG=0;
				for(int j=0;j<nGraph;j++){
					double buffx, buffy;
					model_Graph_ColorChannels[i]->GetPoint(jTG, buffx, buffy);
					if(buffx<pTMinModel || buffx>pTMaxModel){
						model_Graph_ColorChannels[i]->RemovePoint(jTG);
						jTG--;
					}

					jTG++;
				}
				nGraph=model_Graph_ColorChannels_Bands[0][i]->GetN();
				jTG=0;
				for(int j=0;j<nGraph;j++){
					double buffx, buffy;
					model_Graph_ColorChannels_Bands[0][i]->GetPoint(jTG, buffx, buffy);
					if(buffx<pTMinModel || buffx>pTMaxModel){
						for(int iSig=0;iSig<3;iSig++){
							model_Graph_ColorChannels_Bands[iSig][i]->RemovePoint(jTG);
						}
						jTG--;
					}

					jTG++;
				}
			}

			//Check if plotting of neg component is necessary
			for(int j=0;j<model_Graph_ColorChannels[i]->GetN();j++){
				double buffx, buffy;
				model_Graph_ColorChannels[i]->GetPoint(j, buffx, buffy);
				if(buffy<0 && iMeasurementID==0) PlotNegChannels[i]=true;
				if(buffy>0 && iMeasurementID==0) PlotPosChannels[i]=true;
			}
			if(PlotNegChannels[i]){
				//cout<<"model_Graph_ColorChannels_original["<<i<<"]:"<<endl;
				//model_Graph_ColorChannels[i]->Print();


				double set_buffx[model_Graph_ColorChannels[i]->GetN()];
				double set_buffy[model_Graph_ColorChannels[i]->GetN()];

				for(int j=0;j<model_Graph_ColorChannels[i]->GetN();j++){
					double buffx, buffy;
					model_Graph_ColorChannels[i]->GetPoint(j, buffx, buffy);
					set_buffx[j]=buffx;
					set_buffy[j]=buffy;
				}
				model_Graph_ColorChannels_pos[i] = new TGraph(model_Graph_ColorChannels[i]->GetN(), set_buffx, set_buffy);
				model_Graph_ColorChannels_pos[i]->SetLineColor(model_Graph_ColorChannels[i]->GetLineColor());
				model_Graph_ColorChannels_pos[i]->SetLineStyle(linestyleCC[i]);
				model_Graph_ColorChannels_pos[i]->SetLineWidth(linewidthCC);
				model_Graph_ColorChannels_neg[i] = new TGraph(model_Graph_ColorChannels[i]->GetN(), set_buffx, set_buffy);
				model_Graph_ColorChannels_neg[i]->SetLineColor(model_Graph_ColorChannels[i]->GetLineColor());
				model_Graph_ColorChannels_neg[i]->SetLineStyle(linestyleCC_neg[i]);
				model_Graph_ColorChannels_neg[i]->SetLineWidth(linewidthCC);


				for(int iSig=0;iSig<3;iSig++){
					double set_buffx_Bands[model_Graph_ColorChannels_Bands[iSig][i]->GetN()];
					double set_buffy_Bands[model_Graph_ColorChannels_Bands[iSig][i]->GetN()];
					double set_bufferrlowx_Bands[model_Graph_ColorChannels_Bands[iSig][i]->GetN()];
					double set_bufferrlowy_Bands[model_Graph_ColorChannels_Bands[iSig][i]->GetN()];
					double set_bufferrhighx_Bands[model_Graph_ColorChannels_Bands[iSig][i]->GetN()];
					double set_bufferrhighy_Bands[model_Graph_ColorChannels_Bands[iSig][i]->GetN()];

					for(int j=0;j<model_Graph_ColorChannels[i]->GetN();j++){
						double buffx, buffy;
						model_Graph_ColorChannels_Bands[iSig][i]->GetPoint(j, buffx, buffy);
						set_buffx[j]=buffx;
						set_buffy[j]=buffy;
						set_bufferrhighx_Bands[j]=model_Graph_ColorChannels_Bands[iSig][i]->GetErrorXhigh(j);
						set_bufferrhighy_Bands[j]=model_Graph_ColorChannels_Bands[iSig][i]->GetErrorYhigh(j);
						set_bufferrlowx_Bands[j]=model_Graph_ColorChannels_Bands[iSig][i]->GetErrorXlow(j);
						set_bufferrlowy_Bands[j]=model_Graph_ColorChannels_Bands[iSig][i]->GetErrorYlow(j);
					}
					model_Graph_ColorChannels_Bands_pos[iSig][i] = new TGraphAsymmErrors(model_Graph_ColorChannels[i]->GetN(), set_buffx, set_buffy, set_bufferrlowx_Bands, set_bufferrhighx_Bands, set_bufferrlowy_Bands, set_bufferrhighy_Bands);
					model_Graph_ColorChannels_Bands_pos[iSig][i]->SetFillStyle(model_Graph_ColorChannels_Bands[iSig][i]->GetFillStyle());
					model_Graph_ColorChannels_Bands_pos[iSig][i]->SetFillColor(model_Graph_ColorChannels_Bands[iSig][i]->GetFillColor());
					model_Graph_ColorChannels_Bands_pos[iSig][i]->SetLineColor(model_Graph_ColorChannels_Bands[iSig][i]->GetLineColor());
					model_Graph_ColorChannels_Bands_neg[iSig][i] = new TGraphAsymmErrors(model_Graph_ColorChannels[i]->GetN(), set_buffx, set_buffy, set_bufferrlowx_Bands, set_bufferrhighx_Bands, set_bufferrlowy_Bands, set_bufferrhighy_Bands);
					model_Graph_ColorChannels_Bands_neg[iSig][i]->SetFillStyle(model_Graph_ColorChannels_Bands[iSig][i]->GetFillStyle());
					model_Graph_ColorChannels_Bands_neg[iSig][i]->SetFillColor(model_Graph_ColorChannels_Bands[iSig][i]->GetFillColor());
					model_Graph_ColorChannels_Bands_neg[iSig][i]->SetLineColor(model_Graph_ColorChannels_Bands[iSig][i]->GetLineColor());
				}

				bool isContributionPositive[model_Graph_ColorChannels[i]->GetN()];
				for(int j=0;j<model_Graph_ColorChannels[i]->GetN();j++){
					double buffx, buffy;
					model_Graph_ColorChannels[i]->GetPoint(j, buffx, buffy);
					if(buffy<0) isContributionPositive[j]=false;
					else isContributionPositive[j]=true;
				}
				for(int j=0;j<model_Graph_ColorChannels[i]->GetN();j++){
					double buffx, buffy;
					model_Graph_ColorChannels[i]->GetPoint(j, buffx, buffy);
					model_Graph_ColorChannels_neg[i]->SetPoint(j,buffx, -buffy);
					for(int iSig=0;iSig<3;iSig++){
						model_Graph_ColorChannels_Bands_neg[iSig][i]->SetPoint(j,buffx, -buffy);
					}
				}
				int removedPoint_pos=0;
				int removedPoint_neg=0;
				for(int j=0;j<model_Graph_ColorChannels[i]->GetN();j++){
					cout<<"j "<<j<<endl;
					if(isContributionPositive[j]){
						cout<<"remove neg point j-removedPoint_neg "<<j-removedPoint_neg<<endl;
						model_Graph_ColorChannels_neg[i]->RemovePoint(j-removedPoint_neg);
						for(int iSig=0;iSig<3;iSig++){
							model_Graph_ColorChannels_Bands_neg[iSig][i]->RemovePoint(j-removedPoint_neg);
						}
						removedPoint_neg++;
					}
					else{
						cout<<"remove pos point j-removedPoint_pos "<<j-removedPoint_pos<<endl;
						model_Graph_ColorChannels_pos[i]->RemovePoint(j-removedPoint_pos);
						for(int iSig=0;iSig<3;iSig++){
							model_Graph_ColorChannels_Bands_pos[iSig][i]->RemovePoint(j-removedPoint_pos);
						}
						removedPoint_pos++;
					}
				}



			}


			if(!PlotNegChannels[i]){
				if(plotColorChannel[i]) model_Graph_ColorChannels[i]->Draw(plotOption_ColorChannelCentrals);
				if(singlePointModel){
					double buffModel_pT, buffModel_res;
					model_Graph_ColorChannels[i]->GetPoint(0, buffModel_pT, buffModel_res);
					TLine* modelReplacement = new TLine( buff_pT-buff_errpT_low, buffModel_res, buff_pT+buff_errpT_high, buffModel_res);
					modelReplacement->SetLineColor(model_Graph_ColorChannels[i]->GetLineColor());
					modelReplacement->SetLineStyle(model_Graph_ColorChannels[i]->GetLineStyle());
					modelReplacement->SetLineWidth(model_Graph_ColorChannels[i]->GetLineWidth());
					if(plotColorChannel[i] && buffModel_res > y_min && buffModel_res < y_max) modelReplacement->Draw( "same" );
					cout<<i<<" modelReplacement plot1"<<endl; modelReplacement->Print();
				}

				//cout<<"model_Graph_ColorChannels["<<i<<"]:"<<endl;
				//model_Graph_ColorChannels[i]->Print();
				//if(model_Graph_ColorChannels[i]->GetN()>1){
					//if(isSstate) legend->AddEntry(model_Graph_ColorChannels[i],Form("Direct production, O_{%s}^{%s}", ColorChannelNameTexS[i], StateNameTex[iState]),"l");
					//else legend->AddEntry(model_Graph_ColorChannels[i],Form("Direct production, O_{%s}^{%s}", ColorChannelNameTexP[i], StateNameTex[iState]),"l");
				if(plotColorChannel[i]) {
					if(isSstate) legend->AddEntry(model_Graph_ColorChannels[i],Form("%s", ColorChannelNameTexS[i]),"l");
					else legend->AddEntry(model_Graph_ColorChannels[i],Form("%s", ColorChannelNameTexP[i]),"l");
				}
				//	}
				//else{
				//	//if(isSstate) legend->AddEntry(model_Graph_ColorChannels[i],Form("Direct production, O_{%s}^{%s}", ColorChannelNameTexS[i], StateNameTex[iState]),"p");
				//	//else legend->AddEntry(model_Graph_ColorChannels[i],Form("Direct production, O_{%s}^{%s}", ColorChannelNameTexP[i], StateNameTex[iState]),"p");
				//	if(isSstate) legend->AddEntry(model_Graph_ColorChannels[i],Form("%s", ColorChannelNameTexS[i]),"l");
				//	else legend->AddEntry(model_Graph_ColorChannels[i],Form("%s", ColorChannelNameTexP[i]),"p");
				//}

					//if(iMeasurementID==0){
						for(int iSig=2;iSig>-1;iSig--){
							//if(iSig==0&&i==0) cout<<"CHECK_IT_OUT_AGAIN"<<endl;
							//model_Graph_ColorChannels_Bands[iSig][i]->Print();
							if(plotColorChannelBands[i] && plotSigmaBands[iSig]) model_Graph_ColorChannels_Bands[iSig][i]->Draw(plotOptionBands);
						}
						if(plotColorChannel[i]) model_Graph_ColorChannels[i]->Draw(plotOption_ColorChannelCentrals);
						if(singlePointModel){
							double buffModel_pT, buffModel_res, buffModel_reserrlow, buffModel_reserrhigh;
							for(int iSig=2;iSig>-1;iSig--){
								model_Graph_ColorChannels_Bands[iSig][i]->GetPoint(0, buffModel_pT, buffModel_res);
								buffModel_reserrhigh=model_Graph_ColorChannels_Bands[iSig][i]->GetErrorYhigh(0);
								buffModel_reserrlow=model_Graph_ColorChannels_Bands[iSig][i]->GetErrorYlow(0);
								double y_min_box, y_max_box;
								y_min_box=buffModel_res-buffModel_reserrlow;
								y_max_box=buffModel_res+buffModel_reserrhigh;
								if(y_min_box<y_min) y_min_box=y_min;
								if(y_max_box>y_max) y_max_box=y_max;
								TBox* modelReplacementbox = new TBox( buff_pT-buff_errpT_low, y_min_box, buff_pT+buff_errpT_high, y_max_box);
								modelReplacementbox->SetFillStyle(model_Graph_ColorChannels_Bands[iSig][i]->GetFillStyle());
								modelReplacementbox->SetFillColor(model_Graph_ColorChannels_Bands[iSig][i]->GetFillColor());
								modelReplacementbox->SetLineColor(model_Graph_ColorChannels_Bands[iSig][i]->GetLineColor());
								modelReplacementbox->SetLineWidth(0.);
								if(plotColorChannelBands[i] && plotSigmaBands[iSig]) modelReplacementbox->Draw( "same" );

							}
							model_Graph_ColorChannels[i]->GetPoint(0, buffModel_pT, buffModel_res);
							TLine* modelReplacement = new TLine( buff_pT-buff_errpT_low, buffModel_res, buff_pT+buff_errpT_high, buffModel_res);
							modelReplacement->SetLineColor(model_Graph_ColorChannels[i]->GetLineColor());
							modelReplacement->SetLineStyle(model_Graph_ColorChannels[i]->GetLineStyle());
							modelReplacement->SetLineWidth(model_Graph_ColorChannels[i]->GetLineWidth());
							if(plotColorChannel[i] && buffModel_res > y_min && buffModel_res < y_max) modelReplacement->Draw( "same" );
							cout<<i<<" modelReplacement plot2"<<endl; modelReplacement->Print();

						}

					//}

			}
			if(PlotNegChannels[i]){
				if(PlotPosChannels[i]){
					for(int iSig=2;iSig>-1;iSig--){
						if(plotColorChannelBands[i] && plotSigmaBands[iSig]) model_Graph_ColorChannels_Bands_pos[iSig][i]->Draw(plotOptionBands);
					}
					if(plotColorChannel[i]) model_Graph_ColorChannels_pos[i]->Draw("lsame");
					if(singlePointModel){
						double buffModel_pT, buffModel_res, buffModel_reserrlow, buffModel_reserrhigh;
						for(int iSig=2;iSig>-1;iSig--){
							model_Graph_ColorChannels_Bands_pos[iSig][i]->GetPoint(0, buffModel_pT, buffModel_res);
							buffModel_reserrhigh=model_Graph_ColorChannels_Bands_pos[iSig][i]->GetErrorYhigh(0);
							buffModel_reserrlow=model_Graph_ColorChannels_Bands_pos[iSig][i]->GetErrorYlow(0);
							double y_min_box, y_max_box;
							y_min_box=buffModel_res-buffModel_reserrlow;
							y_max_box=buffModel_res+buffModel_reserrhigh;
							if(y_min_box<y_min) y_min_box=y_min;
							if(y_max_box>y_max) y_max_box=y_max;
							TBox* modelReplacementbox = new TBox( buff_pT-buff_errpT_low, y_min_box, buff_pT+buff_errpT_high, y_max_box);
							modelReplacementbox->SetFillStyle(model_Graph_ColorChannels_Bands_pos[iSig][i]->GetFillStyle());
							modelReplacementbox->SetFillColor(model_Graph_ColorChannels_Bands_pos[iSig][i]->GetFillColor());
							modelReplacementbox->SetLineColor(model_Graph_ColorChannels_Bands_pos[iSig][i]->GetLineColor());
							modelReplacementbox->SetLineWidth(0.);
							if(plotColorChannelBands[i] && plotSigmaBands[iSig]) modelReplacementbox->Draw( "same" );

						}
						model_Graph_ColorChannels_pos[i]->GetPoint(0, buffModel_pT, buffModel_res);
						TLine* modelReplacement = new TLine( buff_pT-buff_errpT_low, buffModel_res, buff_pT+buff_errpT_high, buffModel_res);
						modelReplacement->SetLineColor(model_Graph_ColorChannels_pos[i]->GetLineColor());
						modelReplacement->SetLineStyle(model_Graph_ColorChannels_pos[i]->GetLineStyle());
						modelReplacement->SetLineWidth(model_Graph_ColorChannels_pos[i]->GetLineWidth());
						if(plotColorChannel[i] && buffModel_res > y_min && buffModel_res < y_max) modelReplacement->Draw( "same" );

					}

					//cout<<"model_Graph_ColorChannels_pos["<<i<<"]:"<<endl;
					//[i]->Print();
					//if(model_Graph_ColorChannels_pos[i]->GetN()>1){
						//if(isSstate) legend->AddEntry(model_Graph_ColorChannels_pos[i],Form("Direct production, O_{%s}^{%s}", ColorChannelNameTexS[i], StateNameTex[iState]),"l");
						//else legend->AddEntry(model_Graph_ColorChannels_pos[i],Form("Direct production, O_{%s}^{%s}", ColorChannelNameTexP[i], StateNameTex[iState]),"l");
						if(isSstate) legend->AddEntry(model_Graph_ColorChannels_pos[i],Form("%s", ColorChannelNameTexS[i]),"l");
						else legend->AddEntry(model_Graph_ColorChannels_pos[i],Form("%s", ColorChannelNameTexP[i]),"l");
					//}
					//else{
					//	//if(isSstate) legend->AddEntry(model_Graph_ColorChannels_pos[i],Form("Direct production, O_{%s}^{%s}", ColorChannelNameTexS[i], StateNameTex[iState]),"p");
					//	//else legend->AddEntry(model_Graph_ColorChannels_pos[i],Form("Direct production, O_{%s}^{%s}", ColorChannelNameTexP[i], StateNameTex[iState]),"p");
					//	if(isSstate) legend->AddEntry(model_Graph_ColorChannels_pos[i],Form("%s", ColorChannelNameTexS[i]),"p");
					//	else legend->AddEntry(model_Graph_ColorChannels_pos[i],Form("%s", ColorChannelNameTexP[i]),"p");
					//}

					}
				for(int iSig=2;iSig>-1;iSig--){
					if(plotColorChannelBands[i] && plotSigmaBands[iSig]) model_Graph_ColorChannels_Bands_neg[iSig][i]->Draw(plotOptionBands);
				}
				if(plotColorChannel[i]) model_Graph_ColorChannels_neg[i]->Draw("lsame");
				if(singlePointModel){
					double buffModel_pT, buffModel_res, buffModel_reserrlow, buffModel_reserrhigh;
					for(int iSig=2;iSig>-1;iSig--){
						model_Graph_ColorChannels_Bands_neg[iSig][i]->GetPoint(0, buffModel_pT, buffModel_res);
						buffModel_reserrhigh=model_Graph_ColorChannels_Bands_neg[iSig][i]->GetErrorYhigh(0);
						buffModel_reserrlow=model_Graph_ColorChannels_Bands_neg[iSig][i]->GetErrorYlow(0);
						double y_min_box, y_max_box;
						y_min_box=buffModel_res-buffModel_reserrlow;
						y_max_box=buffModel_res+buffModel_reserrhigh;
						if(y_min_box<y_min) y_min_box=y_min;
						if(y_max_box>y_max) y_max_box=y_max;
						TBox* modelReplacementbox = new TBox( buff_pT-buff_errpT_low, y_min_box, buff_pT+buff_errpT_high, y_max_box);
						modelReplacementbox->SetFillStyle(model_Graph_ColorChannels_Bands_neg[iSig][i]->GetFillStyle());
						modelReplacementbox->SetFillColor(model_Graph_ColorChannels_Bands_neg[iSig][i]->GetFillColor());
						modelReplacementbox->SetLineColor(model_Graph_ColorChannels_Bands_neg[iSig][i]->GetLineColor());
						modelReplacementbox->SetLineWidth(0.);
						if(plotColorChannelBands[i] && plotSigmaBands[iSig]) modelReplacementbox->Draw( "same" );

					}
					model_Graph_ColorChannels_neg[i]->GetPoint(0, buffModel_pT, buffModel_res);
					TLine* modelReplacement = new TLine( buff_pT-buff_errpT_low, buffModel_res, buff_pT+buff_errpT_high, buffModel_res);
					modelReplacement->SetLineColor(model_Graph_ColorChannels_neg[i]->GetLineColor());
					modelReplacement->SetLineStyle(model_Graph_ColorChannels_neg[i]->GetLineStyle());
					modelReplacement->SetLineWidth(model_Graph_ColorChannels_neg[i]->GetLineWidth());
					if(plotColorChannel[i] && buffModel_res > y_min && buffModel_res < y_max) modelReplacement->Draw( "same" );

				}
				//cout<<"model_Graph_ColorChannels_neg["<<i<<"]:"<<endl;
				//model_Graph_ColorChannels_neg[i]->Print();
				//if(isSstate) legend->AddEntry(model_Graph_ColorChannels_neg[i],Form("(neg.) Direct production, O_{%s}^{%s}", ColorChannelNameTexS[i], StateNameTex[iState]),"l");
				//else legend->AddEntry(model_Graph_ColorChannels_neg[i],Form("(neg.) Direct production, O_{%s}^{%s}", ColorChannelNameTexP[i], StateNameTex[iState]),"l");
				if(isSstate) legend->AddEntry(model_Graph_ColorChannels_neg[i],Form("%s (neg.)", ColorChannelNameTexS[i]),"l");
				else legend->AddEntry(model_Graph_ColorChannels_neg[i],Form("%s (neg.)", ColorChannelNameTexP[i]),"l");

			}







		}
	}


	if(plotModel){
		for(int iSig=2;iSig>-1;iSig--){
			if(plotTotalModelBands && plotSigmaBands[iSig]) model_Graph_Bands[iSig]->Draw(plotOptionBands);
			if(singlePointModel){
				double buffModel_pT, buffModel_res, buffModel_reserrlow, buffModel_reserrhigh;
				for(int iSig=2;iSig>-1;iSig--){
					model_Graph_Bands[iSig]->GetPoint(0, buffModel_pT, buffModel_res);
					buffModel_reserrhigh=model_Graph_Bands[iSig]->GetErrorYhigh(0);
					buffModel_reserrlow=model_Graph_Bands[iSig]->GetErrorYlow(0);
					double y_min_box, y_max_box;
					y_min_box=buffModel_res-buffModel_reserrlow;
					y_max_box=buffModel_res+buffModel_reserrhigh;
					if(y_min_box<y_min) y_min_box=y_min;
					if(y_max_box>y_max) y_max_box=y_max;
					TBox* modelReplacementbox = new TBox( buff_pT-buff_errpT_low, y_min_box, buff_pT+buff_errpT_high, y_max_box);
					modelReplacementbox->SetFillStyle(model_Graph_Bands[iSig]->GetFillStyle());
					modelReplacementbox->SetFillColor(model_Graph_Bands[iSig]->GetFillColor());
					modelReplacementbox->SetLineColor(model_Graph_Bands[iSig]->GetLineColor());
					modelReplacementbox->SetLineWidth(0.);
					if(plotTotalModelBands && plotSigmaBands[iSig]) modelReplacementbox->Draw( "same" );

				}
			}
		}
		if(plotTotalModel) model_Graph->Draw("lsame");
		if(singlePointModel){
			double buffModel_pT, buffModel_res;
			model_Graph->GetPoint(0, buffModel_pT, buffModel_res);
			TLine* modelReplacement = new TLine( buff_pT-buff_errpT_low, buffModel_res, buff_pT+buff_errpT_high, buffModel_res);
			modelReplacement->SetLineColor(colorModel);
			modelReplacement->SetLineStyle(1);
			modelReplacement->SetLineWidth(linewidthModel);
			if(plotTotalModel && buffModel_res > y_min && buffModel_res < y_max) modelReplacement->Draw( "same" );
		}
/*		for(int iSig=0;iSig<3;iSig++){
			model_Graph_low_Bands[iSig]->Draw("lsame");
			if(singlePointModel){
				double buffModel_pT, buffModel_res;
				model_Graph_low_Bands[iSig]->GetPoint(0, buffModel_pT, buffModel_res);
				TLine* modelReplacement = new TLine( buff_pT-buff_errpT_low, buffModel_res, buff_pT+buff_errpT_high, buffModel_res);
				modelReplacement->SetLineColor(model_Graph_low_Bands[iSig]->GetLineColor());
				modelReplacement->SetLineStyle(model_Graph_low_Bands[iSig]->GetLineStyle());
				modelReplacement->SetLineWidth(model_Graph_low_Bands[iSig]->GetLineWidth());
				modelReplacement->Draw( "same" );
			}
		}
		for(int iSig=0;iSig<3;iSig++){
			model_Graph_high_Bands[iSig]->Draw("lsame");
			if(singlePointModel){
				double buffModel_pT, buffModel_res;
				model_Graph_high_Bands[iSig]->GetPoint(0, buffModel_pT, buffModel_res);
				TLine* modelReplacement = new TLine( buff_pT-buff_errpT_low, buffModel_res, buff_pT+buff_errpT_high, buffModel_res);
				modelReplacement->SetLineColor(model_Graph_high_Bands[iSig]->GetLineColor());
				modelReplacement->SetLineStyle(model_Graph_high_Bands[iSig]->GetLineStyle());
				modelReplacement->SetLineWidth(model_Graph_high_Bands[iSig]->GetLineWidth());
				modelReplacement->Draw( "same" );
			}
		}
*/
	}

	if(plotCSTHuncBand){
		model_CS_THuncBand->Draw(plotOptionBands);
		if(singlePointModel){
			double buffModel_pT, buffModel_res, buffModel_reserrlow, buffModel_reserrhigh;
			model_CS_THuncBand->GetPoint(0, buffModel_pT, buffModel_res);
			buffModel_reserrhigh=model_CS_THuncBand->GetErrorYhigh(0);
			buffModel_reserrlow=model_CS_THuncBand->GetErrorYlow(0);
			double y_min_box, y_max_box;
			y_min_box=buffModel_res-buffModel_reserrlow;
			y_max_box=buffModel_res+buffModel_reserrhigh;
			if(y_min_box<y_min) y_min_box=y_min;
			if(y_max_box>y_max) y_max_box=y_max;
			TBox* modelReplacementbox = new TBox( buff_pT-buff_errpT_low, y_min_box, buff_pT+buff_errpT_high, y_max_box);
			modelReplacementbox->SetFillStyle(model_CS_THuncBand->GetFillStyle());
			modelReplacementbox->SetFillColor(model_CS_THuncBand->GetFillColor());
			modelReplacementbox->SetLineColor(model_CS_THuncBand->GetLineColor());
			modelReplacementbox->SetLineWidth(0.);
			modelReplacementbox->Draw( "same" );

		}

	}

	if(plotCSTHuncLines){

		for(int iTHunc=0;iTHunc<nTHunc;iTHunc++){
			model_CS_THunc[iTHunc]->SetLineStyle(LineStyle_THunc[iTHunc]);
			model_CS_THunc[iTHunc]->SetLineColor(LineColor_THunc[iTHunc]);
			model_CS_THunc[iTHunc]->SetLineWidth(linewidthTHunc);

			if(plotCSTHuncLine[iTHunc]) model_CS_THunc[iTHunc]->Draw("lsame");
			if(singlePointModel){
				double buffModel_pT, buffModel_res;
				model_CS_THunc[iTHunc]->GetPoint(0, buffModel_pT, buffModel_res);
				TLine* modelReplacement = new TLine( buff_pT-buff_errpT_low, buffModel_res, buff_pT+buff_errpT_high, buffModel_res);
				modelReplacement->SetLineColor(model_CS_THunc[iTHunc]->GetLineColor());
				modelReplacement->SetLineStyle(model_CS_THunc[iTHunc]->GetLineStyle());
				modelReplacement->SetLineWidth(model_CS_THunc[iTHunc]->GetLineWidth());
				if(plotCSTHuncLine[iTHunc] && buffModel_res > y_min && buffModel_res < y_max) modelReplacement->Draw( "same" );

			}

		}



	}

	cout<<"PlotFeeddown:::"<<endl;


	if(plotDirectInclusive){
		cout<<"draw direct inclusive"<<endl;
		// inclusive feed-down
		//cout<<"model_Graph_FeedDowns: "<<0<<endl;
		//model_Graph_FeedDowns[0]->Print();
		model_Graph_FeedDowns[StatesCont_c]->SetLineColor(colorFDDirect);
		model_Graph_FeedDowns[StatesCont_c]->SetLineStyle(linestyleFDDirect);
		model_Graph_FeedDowns[StatesCont_c]->SetLineWidth(linewidthFDDirect);
		model_Graph_FeedDowns[StatesCont_c]->Draw("lsame");
		legend->AddEntry(model_Graph_FeedDowns[StatesCont_c],Form("Direct %s",StateNameTex[iState]),"l");
		if(singlePointModel){
			double buffModel_pT, buffModel_res;
			model_Graph_FeedDowns[StatesCont_c]->GetPoint(0, buffModel_pT, buffModel_res);
			TLine* modelReplacement = new TLine( buff_pT-buff_errpT_low, buffModel_res, buff_pT+buff_errpT_high, buffModel_res);
			modelReplacement->SetLineColor(model_Graph_FeedDowns[StatesCont_c]->GetLineColor());
			modelReplacement->SetLineStyle(model_Graph_FeedDowns[StatesCont_c]->GetLineStyle());
			modelReplacement->SetLineWidth(model_Graph_FeedDowns[StatesCont_c]->GetLineWidth());
			if(buffModel_res > y_min && buffModel_res < y_max) modelReplacement->Draw( "same" );			cout<<"ModelReplacement direct component"<<endl; modelReplacement->Print();
		}
	}

	if(plotIndividualFeedDown){
		//cout<<"StatesCont_c "<<StatesCont_c<<endl;
		for(int i=1;i<StatesCont_c;i++){
			cout<<Form("Draw feed-down %s #rightarrow %s", StateNameTex[StatesContributing[i]], StateNameTex[iState])<<endl;
			model_Graph_FeedDowns[i]->SetLineColor(colorFD[StatesContributing[i]]);
			model_Graph_FeedDowns[i]->SetLineStyle(linestyleFD);
			model_Graph_FeedDowns[i]->SetLineWidth(linewidthFD);
			model_Graph_FeedDowns[i]->Draw("lsame");
			legend->AddEntry(model_Graph_FeedDowns[i],Form("%s #rightarrow %s", StateNameTex[StatesContributing[i]], StateNameTex[iState]),"l");
			if(singlePointModel){
				double buffModel_pT, buffModel_res;
				model_Graph_FeedDowns[i]->GetPoint(0, buffModel_pT, buffModel_res);
				TLine* modelReplacement = new TLine( buff_pT-buff_errpT_low, buffModel_res, buff_pT+buff_errpT_high, buffModel_res);
				modelReplacement->SetLineColor(model_Graph_FeedDowns[i]->GetLineColor());
				modelReplacement->SetLineStyle(model_Graph_FeedDowns[i]->GetLineStyle());
				modelReplacement->SetLineWidth(model_Graph_FeedDowns[i]->GetLineWidth());
				if(buffModel_res > y_min && buffModel_res < y_max) modelReplacement->Draw( "same" );
				cout<<i<<". contribution of FD, modelReplacement"<<endl; modelReplacement->Print();
			}
		}
	}

	if(plotInclusiveFeedDown){
		cout<<"draw inclusive feed-down"<<endl;
		// inclusive feed-down
		//cout<<"model_Graph_FeedDowns: "<<0<<endl;
		//model_Graph_FeedDowns[0]->Print();
		model_Graph_FeedDowns[0]->SetLineColor(colorFDInclusive);
		model_Graph_FeedDowns[0]->SetLineStyle(linestyleFDInclusive);
		model_Graph_FeedDowns[0]->SetLineWidth(linewidthFDInclusive);
		model_Graph_FeedDowns[0]->Draw("lsame");
		legend->AddEntry(model_Graph_FeedDowns[0],Form("X #rightarrow %s",StateNameTex[iState]),"l");
		if(singlePointModel){
			double buffModel_pT, buffModel_res;
			model_Graph_FeedDowns[0]->GetPoint(0, buffModel_pT, buffModel_res);
			TLine* modelReplacement = new TLine( buff_pT-buff_errpT_low, buffModel_res, buff_pT+buff_errpT_high, buffModel_res);
			modelReplacement->SetLineColor(model_Graph_FeedDowns[0]->GetLineColor());
			modelReplacement->SetLineStyle(model_Graph_FeedDowns[0]->GetLineStyle());
			modelReplacement->SetLineWidth(model_Graph_FeedDowns[0]->GetLineWidth());
			if(buffModel_res > y_min && buffModel_res < y_max) modelReplacement->Draw( "same" );
			cout<<"ModelReplacement inclusive FD"<<endl; modelReplacement->Print();
		}
	}

	if(plotDirectColorChannels&&plotModel){
	//Replot central color channels, to assure that they are on top of the total model bands... (!!!)
		for(int iColorChannelSequence=0;iColorChannelSequence<ColorChannels_c;iColorChannelSequence++){
			int i=ColorChannelSequence[iColorChannelSequence];
			if(!PlotNegChannels[i]){
				if(plotColorChannel[i]) model_Graph_ColorChannels[i]->Draw(plotOption_ColorChannelCentrals);
				if(singlePointModel){
					double buffModel_pT, buffModel_res;
					model_Graph_ColorChannels[i]->GetPoint(0, buffModel_pT, buffModel_res);
					TLine* modelReplacement = new TLine( buff_pT-buff_errpT_low, buffModel_res, buff_pT+buff_errpT_high, buffModel_res);
					modelReplacement->SetLineColor(model_Graph_ColorChannels[i]->GetLineColor());
					modelReplacement->SetLineStyle(model_Graph_ColorChannels[i]->GetLineStyle());
					modelReplacement->SetLineWidth(model_Graph_ColorChannels[i]->GetLineWidth());
					if(plotColorChannel[i] && buffModel_res > y_min && buffModel_res < y_max) modelReplacement->Draw( "same" );
				}
			}
			if(PlotNegChannels[i]){
				if(PlotPosChannels[i]){
					if(plotColorChannel[i]) model_Graph_ColorChannels_pos[i]->Draw("lsame");
					if(singlePointModel){
						double buffModel_pT, buffModel_res, buffModel_reserrlow, buffModel_reserrhigh;
						model_Graph_ColorChannels_pos[i]->GetPoint(0, buffModel_pT, buffModel_res);
						TLine* modelReplacement = new TLine( buff_pT-buff_errpT_low, buffModel_res, buff_pT+buff_errpT_high, buffModel_res);
						modelReplacement->SetLineColor(model_Graph_ColorChannels_pos[i]->GetLineColor());
						modelReplacement->SetLineStyle(model_Graph_ColorChannels_pos[i]->GetLineStyle());
						modelReplacement->SetLineWidth(model_Graph_ColorChannels_pos[i]->GetLineWidth());
						if(plotColorChannel[i] && buffModel_res > y_min && buffModel_res < y_max) modelReplacement->Draw( "same" );

					}
				}
				model_Graph_ColorChannels_neg[i]->Draw("lsame");
				if(singlePointModel){
					double buffModel_pT, buffModel_res, buffModel_reserrlow, buffModel_reserrhigh;
					model_Graph_ColorChannels_neg[i]->GetPoint(0, buffModel_pT, buffModel_res);
					TLine* modelReplacement = new TLine( buff_pT-buff_errpT_low, buffModel_res, buff_pT+buff_errpT_high, buffModel_res);
					modelReplacement->SetLineColor(model_Graph_ColorChannels_neg[i]->GetLineColor());
					modelReplacement->SetLineStyle(model_Graph_ColorChannels_neg[i]->GetLineStyle());
					modelReplacement->SetLineWidth(model_Graph_ColorChannels_neg[i]->GetLineWidth());
					if(plotColorChannel[i] && buffModel_res > y_min && buffModel_res < y_max) modelReplacement->Draw( "same" );

				}
			}
		}
	}

	if(plotData) data_Graph->Draw("psame");

	axishist->Draw("AXISsame");
    //axisx->Draw("same");
    //axisy->Draw("same");

	if(plotModel) legend->Draw("same");

	////draw latex
	double left=0.715, top=0.1875, textSize=0.03625;
	TLatex *latex=new TLatex();
	latex->SetTextFont(42);
	latex->SetNDC(kTRUE);
	latex->SetTextSize(textSize);
	//double stepLatex=textSize*1.3;

	if(narrowCanvas) latex->SetTextSize(textSize*narrFac);

	bool plotChi2latex=false;
	if(chi2Min>1e-2 && chi2Min/double(ndf)<998) plotChi2latex=true;
	if(PredictionPlot) plotChi2latex=false;

	if(plotChi2latex) latex->DrawLatex(left,top, Form("#chi^{2} / ndf = %1.1f / %d", chi2Min, ndf));
	top=0.25;

	if(plotChi2latex && chi2Prob>0.) latex->DrawLatex(left,top, Form("P(#chi^{2}, ndf) = %1.2G", chi2Prob));
	if(plotChi2latex && chi2Prob<1e-320) latex->DrawLatex(left,top, Form("P(#chi^{2}, ndf) < 1E-320"));

	if(logY) plotCanvas->SetLogy(true);

	left=0.705; top=0.895;
	if(longrapchar) left=0.66125;
	textSize=0.05;
	latex->SetTextFont(42);
	latex->SetTextSize(textSize);
	if(narrowCanvas){
		left-=0.0675;
		latex->SetTextSize(textSize*narrFac);
	}
	latex->DrawLatex(left,top, Form("%s, %s", StateNameTex[iState], rapchar));

	if(!plotModel){
		left=0.7875; top=0.825;
		textSize=0.05;
		latex->SetTextFont(42);
		latex->SetTextSize(textSize);
		if(narrowCanvas){
			left=0.705;
			if(longrapchar) left=0.66125;
			left-=0.0675;
			latex->SetTextSize(textSize*narrFac);
		}
		latex->DrawLatex(left,top, Form("%s", ExpNameTex[iExperiment]));
	}

	left=0.01; top=0.55;
	textSize=0.08;
	latex->SetTextFont(42);
	latex->SetTextSize(textSize);
	if(narrowCanvas) latex->SetTextSize(textSize*narrFac);
	if(iMeasurementID==1) latex->DrawLatex(left,top, "#lambda_{#vartheta}^{#scale[0.7]{HX}}");
	if(iMeasurementID==2) latex->DrawLatex(left,top, "#lambda_{#varphi}^{#scale[0.7]{HX}}");
	if(iMeasurementID==3) latex->DrawLatex(left,top, "#lambda_{#vartheta#varphi}^{#scale[0.7]{HX}}");

	char savename[1000];
	char beginsavename[1000];
	sprintf(beginsavename,"DataModelComp");
	if(PredictionPlot){
		sprintf(beginsavename,"Prediction");
		if(PredictionDataPlot) sprintf(beginsavename,"PredictionData");
	}

	sprintf(savename,"%s/Figures/%s_%s_%s_%s_rap%d.pdf",jobdirname, beginsavename, StateName[iState], MeasurementIDName[iMeasurementID], ExpName[iExperiment], iRap+1);
	if(iMeasurementID==1 && smoothPol) sprintf(savename,"%s/Figures/%s_%s_%s_%s_rap%d_smoothPol.pdf",jobdirname, beginsavename, StateName[iState], MeasurementIDName[iMeasurementID], ExpName[iExperiment], iRap+1);
	cout<<"saved here...: "<<savename<<endl;
	plotCanvas->SaveAs(savename);
	//sprintf(savename,"%s/Figures/PNG_DataModelComp_%s_%s_%s_rap%d.png",jobdirname, StateName[iState], MeasurementIDName[iMeasurementID], ExpName[iExperiment], iRap+1);
	//cout<<"... and saved here: "<<savename<<endl;
	//plotCanvas->SaveAs(savename);



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

TGraph* smoothPolFunc(TGraph *graphToSmooth, int monotonicInt, int c_sign){

	double graphToSmooth_pTmin;
	double graphToSmooth_pTmax;
	double buffy;

	graphToSmooth->GetPoint(0, graphToSmooth_pTmin, buffy);
	graphToSmooth->GetPoint(graphToSmooth->GetN()-1, graphToSmooth_pTmax, buffy);

	char name[200];
	sprintf(name, "fParabola");
	TF1* fParabola;
	fParabola = new TF1(name, paramParabola, graphToSmooth_pTmin , graphToSmooth_pTmax, 4);

	double a=1;
	double b=1;
	double c=1;
	double c_sign_d=c_sign;
	fParabola->SetParameter(0,a);
	fParabola->SetParameter(1,b);
	fParabola->SetParameter(2,c);
	fParabola->FixParameter(3,c_sign_d);
	graphToSmooth->Fit(fParabola, "0", "", graphToSmooth_pTmin , graphToSmooth_pTmax);


	int nBins=500;
	const int buff_nPsmooth=nBins;
	double buff_pT_smooth[buff_nPsmooth];
	double buff_res_smooth[buff_nPsmooth];

	for(int iBin=0;iBin<buff_nPsmooth;iBin++){

		buff_pT_smooth[iBin]=graphToSmooth_pTmin+iBin*(graphToSmooth_pTmax-graphToSmooth_pTmin)/double(buff_nPsmooth);
		buff_res_smooth[iBin]=fParabola->Eval(buff_pT_smooth[iBin]);

		if(monotonicInt==-1 && iBin>0){
			if(buff_res_smooth[iBin]>buff_res_smooth[iBin-1]) buff_res_smooth[iBin]=buff_res_smooth[iBin-1];
		}
		if(monotonicInt==+1 && iBin>0){
			if(buff_res_smooth[iBin]<buff_res_smooth[iBin-1]) buff_res_smooth[iBin]=buff_res_smooth[iBin-1];
		}

	}

	TGraph *smoothedGraph = new TGraph(buff_nPsmooth, buff_pT_smooth, buff_res_smooth);

	return smoothedGraph;

}

TGraphAsymmErrors* CloneCentralValueFromGraph(TGraphAsymmErrors *graphToAdapt, TGraph *graphToAdaptFrom){

	const int nClone=graphToAdapt->GetN();

	double buff_pT_clone[nClone];
	double buff_res_clone[nClone];
	double buff_pT_ErrLow_clone[nClone];
	double buff_pT_ErrHigh_clone[nClone];
	double buff_res_ErrLow_clone[nClone];
	double buff_res_ErrHigh_clone[nClone];

	for(int m=0;m<nClone;m++){

		graphToAdapt->GetPoint(m, buff_pT_clone[m], buff_res_clone[m]);
		buff_pT_ErrHigh_clone[m]=graphToAdapt->GetErrorXhigh(m);
		buff_pT_ErrLow_clone[m]=graphToAdapt->GetErrorXlow(m);
		buff_res_ErrHigh_clone[m]=graphToAdapt->GetErrorYhigh(m);
		buff_res_ErrLow_clone[m]=graphToAdapt->GetErrorYlow(m);

	}

	TGraph *ErrGraphLow = new TGraph(nClone, buff_pT_clone, buff_res_ErrLow_clone);
	TGraph *ErrGraphHigh = new TGraph(nClone, buff_pT_clone, buff_res_ErrHigh_clone);

	const int nCloneTo=graphToAdaptFrom->GetN();

	double buff_pT_cloneTo[nCloneTo];
	double buff_res_cloneTo[nCloneTo];
	double buff_pT_ErrLow_cloneTo[nCloneTo];
	double buff_pT_ErrHigh_cloneTo[nCloneTo];
	double buff_res_ErrLow_cloneTo[nCloneTo];
	double buff_res_ErrHigh_cloneTo[nCloneTo];

	for(int m=0;m<nCloneTo;m++){

		graphToAdaptFrom->GetPoint(m, buff_pT_cloneTo[m], buff_res_cloneTo[m]);
		buff_pT_ErrHigh_cloneTo[m]=0;
		buff_pT_ErrLow_cloneTo[m]=0;
		buff_res_ErrHigh_cloneTo[m]=ErrGraphHigh->Eval(buff_pT_cloneTo[m]);
		buff_res_ErrLow_cloneTo[m]=ErrGraphLow->Eval(buff_pT_cloneTo[m]);


	}

	TGraphAsymmErrors* ClonedGraphAsymmErrors = new TGraphAsymmErrors(nCloneTo, buff_pT_cloneTo,buff_res_cloneTo,buff_pT_ErrLow_cloneTo,buff_pT_ErrHigh_cloneTo,buff_res_ErrLow_cloneTo,buff_res_ErrHigh_cloneTo);

	return ClonedGraphAsymmErrors;

}

TGraphAsymmErrors* smoothPolFuncErrors(TGraphAsymmErrors *graphToAdapt, int polOrder){

	double graphToSmooth_pTmin;
	double graphToSmooth_pTmax;
	double buffy;

	graphToAdapt->GetPoint(0, graphToSmooth_pTmin, buffy);
	graphToAdapt->GetPoint(graphToAdapt->GetN()-1, graphToSmooth_pTmax, buffy);

	char name[200];
	sprintf(name, "fParabola");
	TF1* fParabola;
	fParabola = new TF1(name, paramParabola, graphToSmooth_pTmin , graphToSmooth_pTmax, 4);

	double a=1;
	double b=1;
	double c=1;

	if(polOrder==1) c=0;

	double c_sign_d=0;
	fParabola->SetParameter(0,a);
	fParabola->SetParameter(1,b);
	fParabola->SetParameter(2,c);
	if(polOrder==1) fParabola->FixParameter(2,c);
	fParabola->FixParameter(3,c_sign_d);



	const int nClone=graphToAdapt->GetN();
	double buff_pT_clone[nClone];
	double buff_res_clone[nClone];
	double buff_pT_ErrLow_clone[nClone];
	double buff_pT_ErrHigh_clone[nClone];
	double buff_res_ErrLow_clone[nClone];
	double buff_res_ErrHigh_clone[nClone];

	for(int m=0;m<nClone;m++){

		graphToAdapt->GetPoint(m, buff_pT_clone[m], buff_res_clone[m]);
		buff_pT_ErrHigh_clone[m]=graphToAdapt->GetErrorXhigh(m);
		buff_pT_ErrLow_clone[m]=graphToAdapt->GetErrorXlow(m);
		buff_res_ErrHigh_clone[m]=graphToAdapt->GetErrorYhigh(m);
		buff_res_ErrLow_clone[m]=graphToAdapt->GetErrorYlow(m);

	}

	TGraph *ErrGraphLow = new TGraph(nClone, buff_pT_clone, buff_res_ErrLow_clone);
	TGraph *ErrGraphHigh = new TGraph(nClone, buff_pT_clone, buff_res_ErrHigh_clone);

	ErrGraphLow->Fit(fParabola, "0", "", graphToSmooth_pTmin , graphToSmooth_pTmax);
	for(int m=0;m<nClone;m++){
		buff_res_ErrLow_clone[m]=fParabola->Eval(buff_pT_clone[m]);
	}

	fParabola->SetParameter(0,a);
	fParabola->SetParameter(1,b);
	fParabola->SetParameter(2,c);
	if(polOrder==1) fParabola->FixParameter(2,c);
	fParabola->FixParameter(3,c_sign_d);
	ErrGraphHigh->Fit(fParabola, "0", "", graphToSmooth_pTmin , graphToSmooth_pTmax);
	for(int m=0;m<nClone;m++){
		buff_res_ErrHigh_clone[m]=fParabola->Eval(buff_pT_clone[m]);
	}

	TGraphAsymmErrors* smoothedGraphErrors = new TGraphAsymmErrors(nClone, buff_pT_clone,buff_res_clone,buff_pT_ErrLow_clone,buff_pT_ErrHigh_clone,buff_res_ErrLow_clone,buff_res_ErrHigh_clone);

	return smoothedGraphErrors;

}




Double_t paramParabola(Double_t *x, Double_t *par){

	double quadraticCoef=par[2];
	if(par[3]>0.5) quadraticCoef=TMath::Abs(par[2]);
	if(par[3]<0.5) quadraticCoef=(-1)*TMath::Abs(par[2]);


	Double_t result=par[0]+x[0]*par[1]+x[0]*x[0]*quadraticCoef;
	return result;
}

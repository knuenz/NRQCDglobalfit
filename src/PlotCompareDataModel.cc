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


void plotComp(int iState, int iMeasurementID, int iExperiment, int iRap, char jobdirname[200], char rapchar[200], TGraphAsymmErrors *data_Graph, TGraph *model_Graph, TGraph *model_Graph_low, TGraph *model_Graph_high, bool plotInclusiveFeedDown, bool plotIndividualFeedDown, bool plotDirectColorChannels, const int StatesCont_c, const int ColorChannels_c, TGraph *model_Graph_FeedDowns[500], TGraph *model_Graph_ColorChannels[500], /*TGraphAsymmErrors *model_Graph_ColorChannels[3][500],*/ TGraph *data_Graph_nopolcorr, vector<int> StatesContributing, double chi2Min, double chi2Prob, int ndf, bool HPbool, double pTMinModel, bool longrapchar, TGraphAsymmErrors *model_Graph_ColorChannels_Bands[3][500], TGraph *model_Graph_low_Bands[3], TGraph *model_Graph_high_Bands[3], TGraphAsymmErrors *model_Graph_Bands[3], TGraph *model_CS_THunc[3], TGraphAsymmErrors *model_CS_THuncBand);
vector<double> addPolarizations(vector<vector<double> > lamMatrix, vector<double> lamVecContributions);
void FindMPV(TH1* PosteriorDist , double& MPV , double& MPVerrorLow, double& MPVerrorHigh, int MPValgo, int nSigma);

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
  	int 	useOnlyState=999;

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

    double chi2Prob=TMath::Prob(chi2Min, ndf);


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


    cout<<"Loop through data to plot:"<<endl;

	bool DataSelected=false;

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

				    	DataPresentAndSelected[iRap][iP]=false;

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
								&& DataModelObject[iRap][iP]->getpTMax() <= pTMax
								&& DataModelObject[iRap][iP]->getyMin() >= rapMin
								&& DataModelObject[iRap][iP]->getyMax() <= rapMax
							) DataSelected=true;

							if(DataModelObject[iRap][iP]->getpTMin() >= 29 && iMeasurementID==1 && pTMin==35) DataSelected=true;

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

				    for(int iP = 0; iP < NRQCDvars::nMaxPtBins; iP++){

				    	if(DataPresentAndSelected[iRap][iP]){

				    		cout<<"Plot state: "<<NRQCDvars::StateName[iState]<<" Measurement: "<<NRQCDvars::MeasurementIDName[iMeasurementID]<<" Experiment: "<<NRQCDvars::ExpName[iExperiment]<<" iRap: "<<iRap<<" iP "<<iP<<endl;

							dvector ObjectLikelihoodVec;
							double ModelPrediction;
							double errModelPrediction;
							double polCorrFactor;
							//TODO: Improve model uncertainty (draw from PPD or so)

							dcube directProductionCube;
							dmatrix promptProductionMatrix;


							bool LoopThroughPPD=false;
							if(MinimizerMH) LoopThroughPPD=true;

							//if(iMeasurementID==0) LoopThroughPPD=false;
							//LoopThroughPPD=false;

						    if(!LoopThroughPPD){


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



						    }

						    else {

						    	cout<<"LOOP"<<endl;
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


								int n_events = int( outputTreeAllSamplings->GetEntries() );
								cout<<n_events<<" events: Looping through nTuple of PPD"<<endl;

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
												//if(iMeasurementID==1 && j==0) cout<<"CCprediction[0] "<<CCprediction[j]<<endl;

											}
										//}

									//}

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

							  	h_ModelPrediction_min[0]=ModelPredictionTree->GetMinimum("ModelPrediction")-expandMinMaxBy*ModelPredictionTree->GetMinimum("ModelPrediction");
							  	h_ModelPrediction_max[0]=ModelPredictionTree->GetMaximum("ModelPrediction")+expandMinMaxBy*ModelPredictionTree->GetMaximum("ModelPrediction");
							  	//cout<<"h_ModelPrediction_min[0] "<<h_ModelPrediction_min[0]<<" h_ModelPrediction_max[0] "<<h_ModelPrediction_max[0]<<endl;
							  	sprintf(hist_name,"h_ModelPrediction");
								h_ModelPrediction[0] = new TH1D( hist_name, hist_name, nBins_h, h_ModelPrediction_min[0], h_ModelPrediction_max[0] );
								sprintf(projectchar,"ModelPrediction>>h_ModelPrediction");
								ModelPredictionTree->Draw(projectchar);

							  	h_ModelPrediction_min[1]=ModelPredictionTree->GetMinimum("polCorrFactor")-expandMinMaxBy*ModelPredictionTree->GetMinimum("polCorrFactor");
							  	h_ModelPrediction_max[1]=ModelPredictionTree->GetMaximum("polCorrFactor")+expandMinMaxBy*ModelPredictionTree->GetMaximum("polCorrFactor");
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

								if(iExperiment==0 && iMeasurementID==0 && iRap==0 && iP==3){
									h_ModelPrediction[0]->SaveAs("tmp_Modelpred.root");
									h_ModelPrediction[1]->SaveAs("tmp_polcorrpred.root");
									h_ModelPrediction[0]->Print();
									cout<<h_ModelPrediction[0]->GetMean()<<endl;
									cout<<h_ModelPrediction[0]->GetRMS()<<endl;
									cout<<buff_MPV<<endl;
									cout<<buff_errlow<<endl;
									cout<<buff_errhigh<<endl;
								}

								ObjectLikelihoodVec=DataModelObject[iRap][iP]->getObjectLikelihood(Oi_MPV, Np_BR_MPV, Np_US_MPV, true, directProductionCube, promptProductionMatrix, polCorrFactor);
								ModelPrediction=ObjectLikelihoodVec[1];

								cout<<"MPV directProductionCube"<<endl;
								cout<<directProductionCube<<endl;
								cout<<"MPV promptProductionMatrix"<<endl;
								cout<<promptProductionMatrix<<endl;
								cout<<"MPV polCorrFactor"<<endl;
								cout<<polCorrFactor<<endl;

								//PAPER: probably use MPV as central model, given that we now display the bands, and not the central values anymore
								bool useHistModeAsCentralModel=true;

								if(useHistModeAsCentralModel){
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

								}
								else{
									model_centralval.push_back(ModelPrediction);
									//model_errlow_centralval.push_back(fabs(buff_errlow)+ModelPrediction-buff_MPV);
									//model_errhigh_centralval.push_back(fabs(buff_errhigh)-ModelPrediction+buff_MPV);
									model_errlow_centralval.push_back(fabs(buff_errlow));
									model_errhigh_centralval.push_back(fabs(buff_errhigh));
									polCorrFactorVec.push_back(polCorrFactor);
									model_errlow_centralval_1Sig.push_back(fabs(buff_errlow_Bands[0]));
									model_errhigh_centralval_1Sig.push_back(fabs(buff_errhigh_Bands[0]));
									model_errlow_centralval_2Sig.push_back(fabs(buff_errlow_Bands[1]));
									model_errhigh_centralval_2Sig.push_back(fabs(buff_errhigh_Bands[1]));
									model_errlow_centralval_3Sig.push_back(fabs(buff_errlow_Bands[2]));
									model_errhigh_centralval_3Sig.push_back(fabs(buff_errhigh_Bands[2]));

									cout<<"model_centralval "<<ModelPrediction<<endl;
									cout<<"model_errlow_centralval "<<fabs(buff_errlow)<<endl;
									cout<<"model_errhigh_centralval "<<fabs(buff_errhigh)<<endl;
								}



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

										h_ModelPrediction_min_CC=ModelPredictionTree->GetMinimum(branch_name)-expandMinMaxBy*ModelPredictionTree->GetMinimum(branch_name);
										h_ModelPrediction_max_CC=ModelPredictionTree->GetMaximum(branch_name)+expandMinMaxBy*ModelPredictionTree->GetMaximum(branch_name);
										sprintf(hist_name_CC,"h_ModelPrediction_CC");
										h_ModelPrediction_CC = new TH1D( hist_name_CC, hist_name_CC, nBins_h_CC, h_ModelPrediction_min_CC, h_ModelPrediction_max_CC );
										sprintf(projectchar_CC,"%s>>h_ModelPrediction_CC",branch_name);
										ModelPredictionTree->Draw(projectchar_CC);

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
				    		data_pTmean.push_back(DataModelObject[iRap][iP]->getpTMean());
				    		data_errlow_pT.push_back(DataModelObject[iRap][iP]->getpTMean()-DataModelObject[iRap][iP]->getpTMin());
				    		data_errhigh_pT.push_back(DataModelObject[iRap][iP]->getpTMax()-DataModelObject[iRap][iP]->getpTMean());

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
						    for(int i=0; i<StatesCont;i++){
						    	if(i==0){
						    		if(StatesCont>2) model_FeedDownVec.push_back(model_inclusive_FeedDown_);
						    		if(StatesCont==2){
										Buff_model_FeedDownVec=promptProductionMatrix.at(StatesContributing[1]);
										model_FeedDownVec.push_back(Buff_model_FeedDownVec[iMeasurementID]);
						    		}
						    	}
						    	else{
									Buff_model_FeedDownVec=promptProductionMatrix.at(StatesContributing[i]);
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

							if(!LoopThroughPPD){
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

									nSigmaBand=1; if(iMeasurementID>0) nSigmaBand=0.;
									model_directProduction_BandsVec_errlow_1sig.push_back( nSigmaBand*Oi_errlow[iState][j]/Oi_MPV[iState][j]*model_directProductionVec[j] );
									model_directProduction_BandsVec_errhigh_1sig.push_back(nSigmaBand*Oi_errhigh[iState][j]/Oi_MPV[iState][j]*model_directProductionVec[j]);
									nSigmaBand=2; if(iMeasurementID>0) nSigmaBand=0.;
									model_directProduction_BandsVec_errlow_2sig.push_back( nSigmaBand*Oi_errlow[iState][j]/Oi_MPV[iState][j]*model_directProductionVec[j] );
									model_directProduction_BandsVec_errhigh_2sig.push_back(nSigmaBand*Oi_errhigh[iState][j]/Oi_MPV[iState][j]*model_directProductionVec[j]);
									nSigmaBand=3; if(iMeasurementID>0) nSigmaBand=0.;
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

							//if(iMeasurementID>0){
							//	int nColorChannels_state;
							//	bool isSstate=(StateQuantumID[iState] > NRQCDvars::quID_S)?false:true;
							//	if(isSstate) nColorChannels_state=NRQCDvars::nColorChannels_S;
							//	else nColorChannels_state=NRQCDvars::nColorChannels_P;
                            //
							//	dvector model_directProduction_BandsVec;
							//    dvector model_directProduction_BandsVec_errlow_1sig;
							//    dvector model_directProduction_BandsVec_errhigh_1sig;
							//    dvector model_directProduction_BandsVec_errlow_2sig;
							//    dvector model_directProduction_BandsVec_errhigh_2sig;
							//    dvector model_directProduction_BandsVec_errlow_3sig;
							//    dvector model_directProduction_BandsVec_errhigh_3sig;
                            //
							//	int nSigmaBand;
							//	double DummyPush=0.999;
							//	for (int j=0; j<nColorChannels_state; j++){
                            //
							//		model_directProduction_BandsVec.push_back(DummyPush);
                            //
							//		model_directProduction_BandsVec_errlow_1sig.push_back( DummyPush );
							//		model_directProduction_BandsVec_errhigh_1sig.push_back(DummyPush);
							//		model_directProduction_BandsVec_errlow_2sig.push_back( DummyPush );
							//		model_directProduction_BandsVec_errhigh_2sig.push_back(DummyPush);
							//		model_directProduction_BandsVec_errlow_3sig.push_back( DummyPush );
							//		model_directProduction_BandsVec_errhigh_3sig.push_back(DummyPush);
                            //
							//	}
                            //
							//	//cout<<"model_directProduction_BandsVec"<<endl;
							//	//cout<<model_directProduction_BandsVec<<endl;
                            //
							//	model_directProduction_Bands.push_back(model_directProduction_BandsVec);
							//	model_directProduction_Bands_errlow_1sig.push_back(model_directProduction_BandsVec_errlow_1sig);
							//	model_directProduction_Bands_errhigh_1sig.push_back(model_directProduction_BandsVec_errhigh_1sig);
							//	model_directProduction_Bands_errlow_2sig.push_back(model_directProduction_BandsVec_errlow_2sig);
							//	model_directProduction_Bands_errhigh_2sig.push_back(model_directProduction_BandsVec_errhigh_2sig);
							//	model_directProduction_Bands_errlow_3sig.push_back(model_directProduction_BandsVec_errlow_3sig);
							//	model_directProduction_Bands_errhigh_3sig.push_back(model_directProduction_BandsVec_errhigh_3sig);
                            //
							//}

							//cout<<"model_directProduction_Bands"<<endl;
							//cout<<model_directProduction_Bands<<endl;

							//TODO: after ifs: calc ObjetLikelihood with 'central' inputs, to be used for component calculations
							ObjectLikelihoodVec=DataModelObject[iRap][iP]->getObjectLikelihood(Oi_MPV, Np_BR_MPV, Np_US_MPV, true, directProductionCube, promptProductionMatrix, polCorrFactor);


						    }


				    }

				    //cout<<"model_directProduction_THunc_CS"<<endl;
				    //cout<<model_directProduction_THunc_CS<<endl;

					const int StatesCont_c=StatesCont;
					int nColorChannels_state;
					bool isSstate=(StateQuantumID[iState] > NRQCDvars::quID_S)?false:true;
					if(isSstate) nColorChannels_state=NRQCDvars::nColorChannels_S;
					else nColorChannels_state=NRQCDvars::nColorChannels_P;
					const int ColorChannels_c=nColorChannels_state;


				    if(nPtBinsSel>0){
						cout << "Plot iState=" << StateName[iState] << ", iMeasurementID=" << MeasurementIDName[iMeasurementID] << " from "<< ExpName[iExperiment] << endl;
				    	cout<<"rap "<<iRap+1<<endl;

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

						for(int iP = 0; iP < nPtBinsSel; iP++){

							bool correctForPol=true;

							//if(iP<3) correctForPol=false;

							if(!correctForPol) polCorrFactorVec[iP]=1.;

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

								d_model_CS_THunc_m1[iP] *=	polCorrFactorVec[iP];
								d_model_CS_THunc_central[iP] *=	polCorrFactorVec[iP];
								d_model_CS_THunc_p1[iP] *=	polCorrFactorVec[iP];

								for(int i=0;i<ColorChannels_c;i++){
									model_directProduction_Bands[iP][i] *=	polCorrFactorVec[iP];
									model_directProduction_Bands_errlow_1sig[iP][i] *=	polCorrFactorVec[iP];
									model_directProduction_Bands_errlow_2sig[iP][i] *=	polCorrFactorVec[iP];
									model_directProduction_Bands_errlow_3sig[iP][i] *=	polCorrFactorVec[iP];
									model_directProduction_Bands_errhigh_1sig[iP][i] *=	polCorrFactorVec[iP];
									model_directProduction_Bands_errhigh_2sig[iP][i] *=	polCorrFactorVec[iP];
									model_directProduction_Bands_errhigh_3sig[iP][i] *=	polCorrFactorVec[iP];
								}
								for(int i=0;i<StatesCont_c;i++){
									model_FeedDown[iP][i] *=	polCorrFactorVec[iP];
								}

							}

						}

						TGraphAsymmErrors *data_Graph = new TGraphAsymmErrors(nPtBinsSel,d_data_pTmean,d_data_centralval,d_data_errlow_pT,d_data_errhigh_pT,d_data_errlow_centralval,d_data_errhigh_centralval);
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

						model_Graph->Print();
						model_Graph_low->Print();
						model_Graph_high->Print();

						//Make TGraphs for model contributions

						TGraph *model_Graph_FeedDowns[StatesCont_c];//fist element: inclusive FeedDown

						//cout<<model_FeedDown<<endl;
						//cout<<"StatesCont_c "<<StatesCont_c<<endl;
						for(int i=0;i<StatesCont_c;i++){
							for(int iP = 0; iP < nPtBinsSel; iP++){
								d_model_centralval[iP] =  model_FeedDown[iP][i];
								if(iMeasurementID==0) d_model_centralval[iP];
							}
							model_Graph_FeedDowns[i] = new TGraphAsymmErrors(nPtBinsSel,d_data_pTmean,d_model_centralval);
							//model_Graph_FeedDowns[i]->Print();
						}


						//PAPER: decide here if you want to use the MPV of the color channels, or the prediction from the Op_MPVs, add uncertainties (TGraphAsymmErrors)
						TGraph *model_Graph_ColorChannels[ColorChannels_c];

						for(int i=0;i<ColorChannels_c;i++){
							for(int iP = 0; iP < nPtBinsSel; iP++){
								d_model_centralval[iP] =  model_directProduction[iP][i];
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
						bool plotDirectColorChannels=true;
						if(ColorChannels_c<2) plotDirectColorChannels=false;


						char rapchar[200];
						if(isAbsRap) sprintf(rapchar,"%1.1f < |#it{y}| < %1.1f",rapMinObject, rapMaxObject);
						else sprintf(rapchar,"%1.1f < #it{y} < %1.1f",rapMinObject, rapMaxObject);

						bool longrapchar=true;

						if(isAbsRap&&rapMinObject<1e-3){
							sprintf(rapchar,"|#it{y}| < %1.1f", rapMaxObject);
							longrapchar=false;
						}

						plotComp( iState,  iMeasurementID,  iExperiment,  iRap,  jobdirname, rapchar, data_Graph, model_Graph, model_Graph_low, model_Graph_high, plotInclusiveFeedDown, plotIndividualFeedDown, plotDirectColorChannels, StatesCont_c, ColorChannels_c, model_Graph_FeedDowns, model_Graph_ColorChannels, data_Graph_nopolcorr, StatesContributing_, chi2Min, chi2Prob, ndf, HPbool, pTMinModel, longrapchar, model_Graph_ColorChannels_Bands, model_Graph_low_Bands, model_Graph_high_Bands, model_Graph_Bands, model_CS_THunc, model_CS_THuncBand);

				}


				}//rap loop



			}
		}
	}











		return 0;

  	}



void plotComp(int iState, int iMeasurementID, int iExperiment, int iRap, char jobdirname[200], char rapchar[200], TGraphAsymmErrors *data_Graph, TGraph *model_Graph, TGraph *model_Graph_low, TGraph *model_Graph_high, bool plotInclusiveFeedDown, bool plotIndividualFeedDown, bool plotDirectColorChannels, const int StatesCont_c, const int ColorChannels_c, TGraph *model_Graph_FeedDowns[500], TGraph *model_Graph_ColorChannels[500], /*TGraphAsymmErrors *model_Graph_ColorChannels[3][500],*/ TGraph *data_Graph_nopolcorr, vector<int> StatesContributing, double chi2Min, double chi2Prob, int ndf, bool HPbool, double pTMinModel, bool longrapchar, TGraphAsymmErrors *model_Graph_ColorChannels_Bands[3][500], TGraph *model_Graph_low_Bands[3], TGraph *model_Graph_high_Bands[3], TGraphAsymmErrors *model_Graph_Bands[3], TGraph *model_CS_THunc[3], TGraphAsymmErrors *model_CS_THuncBand){
	gROOT->Reset();
	gStyle->SetOptStat(11);
	gStyle->SetOptFit(101);
	gStyle->SetFillColor(kWhite);

	gStyle->SetPadBottomMargin(0.12);
	gStyle->SetPadLeftMargin(0.12);
	gStyle->SetPadRightMargin(0.02);
	gStyle->SetPadTopMargin(0.02);

	bool isSstate=(StateQuantumID[iState] > NRQCDvars::quID_S)?false:true;

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
	bool plotColorChannelBands[4]={false, true, true,     false};
	bool plotColorChannel[4]={true, true, true,     false};
	bool plotTotalModelBands=false;
	bool plotTotalModel=true;
	bool plotCSTHuncLines=true;
	bool plotCSTHuncBand=false;

	///////////////////////

	const int nTHunc=3;
	bool plotSigmaBands[3]={true, true, true};
	int ColorChannelSequence[4]={0, 2, 1,     999};
	bool plotCSTHuncLine[nTHunc]={true, true, true};

	TCanvas *plotCanvas = new TCanvas("plotCanvas","plotCanvas",1200,800);
	plotCanvas->SetFillColor(kWhite);
	plotCanvas->SetFrameBorderMode(0);


	double x_min=0.;
	double x_max=65;
	double y_min;
	double y_max;
	if(iMeasurementID==0) { y_min=5e-2; y_max =1e1; }
	if(iMeasurementID==1) { y_min=-1.3; y_max =1.3; }
	if(iMeasurementID==2) { y_min=-0.7; y_max =0.7; }
	if(iMeasurementID==3) { y_min=-0.7; y_max =0.7; }

	if(iMeasurementID==0) { y_min=1e-3; y_max =1.5e1; }
	if(iMeasurementID==0 && iExperiment==1) { y_min=1e-1; y_max =1.5e2; }
	if(iMeasurementID==0) { x_min=8; x_max=50; }
	if(iMeasurementID==0 && iExperiment==1) { x_min=2; x_max=15; }
	if(iMeasurementID!=0) { x_min=8; x_max=65; }

	//Ups3S

	if(iState==10){

		if(iMeasurementID==0) { y_min=1e-6; y_max =1e-1; x_min=20; x_max =135;  }
		if(iMeasurementID==1) { y_min=-1.3; y_max =1.3; x_min=8; x_max =70; }
		if(iMeasurementID==2) { y_min=-0.7; y_max =0.7; x_min=8; x_max =70;  }
		if(iMeasurementID==3) { y_min=-0.7; y_max =0.7; x_min=8; x_max =70;  }

	}


	TH1F *axishist = new TH1F;
	axishist = plotCanvas->DrawFrame(x_min,y_min,x_max,y_max);

	axishist->SetTitle(0);
	axishist->GetXaxis()->SetTitle("#it{p}_{T} [GeV]");
	if(iMeasurementID==0) axishist->GetYaxis()->SetTitle(Form("%s",NRQCDvars::MeasurementIDNameTex[iMeasurementID]));
	else axishist->GetYaxis()->SetTitle("");
	axishist->GetYaxis()->SetTitleSize(0.05);
	axishist->GetXaxis()->SetTitleSize(0.05);
	axishist->GetYaxis()->SetTitleOffset(0.8);
	axishist->GetXaxis()->SetTitleOffset(1.1);


	bool logY=true;

	int colorData=1;
	int colorModel=600;
	int colorData_nopolcorr=418;
	double linewidthModel=2;
	double linewidthModelErrors=1;

	int colorCC[6]={435, 1, 1, 1, 1, 1};//{632, 418, 616, 800, 0, 0};
	int colorCC_neg[6]={632, 632, 632, 632, 632, 632};
	int colorFD[6]={418, 616, 800, 632, 0, 0};
	int linestyleCC[6]={4, 5, 7, 9, 10, 6};
	int linestyleCC_neg[6]={4, 5, 7, 9, 10, 6};
	int linestyleFD=3;
	int MarkerStyle[6]={21, 24, 25, 26, 27, 28};//{632, 418, 616, 800, 0, 0};
	double linewidthCC=2;
	double linewidthFD=1;

	int styleModel[3]={2, 4, 3};

	//int FillStyle_sig[6]={1001, 1001, 1001, 1, 1, 1};
	int FillStyle_sig[6]={1001, 1001, 1001, 1, 1, 1};
	int LineStyle_Central[6]={1, 1, 1, 1, 1, 1};
	int LineColor_Central[6]={923, 416+2,632+2, 1, 1, 1};
	int LineColor_nSig[3][6]={
			{16, 416-3,632+0, 1, 1, 1},
			{17, 416-7,632-7, 1, 1, 1},
			{18, 416-10,632-10, 1, 1, 1}
	};

	int TotalModel_FillStyle_sig=1001;
	int TotalModel_LineStyle_Central=1;
	int TotalModel_LineColor_Central=600+0;
	int TotalModel_LineColor_nSig[3]={600-7, 600-9, 600-10};

	//TotalModel_LineColor_Central=14;
	//TotalModel_LineColor_nSig[0]=16;
	//TotalModel_LineColor_nSig[1]=17;
	//TotalModel_LineColor_nSig[2]=18;


	int MarkerStyleExp[3]={20, 30, 25};
	int MarkerStylePolRap[2]={20, 24};
	double MarkerSizeExp[3]={1.5, 1.8, 1.5};

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
		colorCC[3]=1;
		colorCC[4]=1;
		colorCC[5]=1;

		colorCC_neg[0]=923;
		colorCC_neg[1]=418;
		colorCC_neg[2]=634;
		colorCC_neg[3]=435;
		colorCC_neg[4]=1;
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

	bool PlotModel=true;
	if(iExperiment==1 && pTMinModel>9.99) PlotModel=false;

	if(HPbool){
		int nGraph=model_Graph->GetN();
		int jTG=0;
		for(int j=0;j<nGraph;j++){
			cout<<"j = "<<j<<endl;
			cout<<"jTG = "<<jTG<<endl;
			model_Graph->Print();
			double buffx, buffy;
			model_Graph->GetPoint(jTG, buffx, buffy);
			cout<<"working bin at pT = "<<buffx<<endl;
			if(buffx<pTMinModel){
				if(iExperiment==1){
					cout<<"remove bin at pT = "<<buffx<<endl;
				}
				model_Graph->RemovePoint(0);
				jTG--;
			}

			jTG++;
		}
		nGraph=model_Graph_low->GetN();
		jTG=0;
		for(int j=0;j<nGraph;j++){
			double buffx, buffy;
			model_Graph_low->GetPoint(jTG, buffx, buffy);
			if(buffx<pTMinModel){
				model_Graph_low->RemovePoint(jTG);
				model_Graph_high->RemovePoint(jTG);
				for(int iSig=0;iSig<3;iSig++){
					model_Graph_low_Bands[iSig]->RemovePoint(jTG);
					model_Graph_high_Bands[iSig]->RemovePoint(jTG);
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
			if(buffx<pTMinModel){
				for(int iTHunc=0;iTHunc<nTHunc;iTHunc++){
					model_CS_THunc[iTHunc]->RemovePoint(jTG);
				}
				model_CS_THuncBand->RemovePoint(jTG);
				jTG--;
			}

			jTG++;
		}

	}



	int nLegendEntries=1;
	if(plotDirectColorChannels) nLegendEntries+=ColorChannels_c;
	if(plotIndividualFeedDown) nLegendEntries+=StatesCont_c;
	if(plotInclusiveFeedDown) nLegendEntries+=1;

	if(iMeasurementID==0) nLegendEntries++;

	if(iMeasurementID==0&&HPbool) nLegendEntries--;


	double buff_pT;
	double buff_errpT_low;
	double buff_errpT_high;
	double buff_res;
	double buff_errres;

	bool singlePointModel=false;
	if(model_Graph->GetN()<2){
		singlePointModel=true;

		double buffModel_pT, buffModel_res;
		double buffData_pT, buffData_err_pTlow, buffData_err_pThigh, buffData_res;
		model_Graph->GetPoint(0, buffModel_pT, buffModel_res);

		for(int iData=0; iData<data_Graph->GetN();iData++){
			data_Graph->GetPoint(iData, buffData_pT, buffData_res);
			buffData_err_pTlow=data_Graph->GetErrorXlow(iData);
			buffData_err_pThigh=data_Graph->GetErrorXhigh(iData);
			if(buffData_pT-buffData_err_pTlow<buffModel_pT && buffData_pT+buffData_err_pThigh>buffModel_pT){
				buff_pT=buffData_pT;
				buff_errpT_low=buffData_err_pTlow;
				buff_errpT_high=buffData_err_pThigh;
			}
		}
	}


	if(iMeasurementID==0) if(!HPbool) data_Graph_nopolcorr->Draw("psame");



	double delta_y_legend=0.085*nLegendEntries;
	double max_y_legend;
	double min_x_legend=0.65;
	double max_x_legend=0.95;
	min_x_legend=0.775;
	if(iMeasurementID==0) max_y_legend=0.85;
	if(iMeasurementID>0 && iMeasurementID<4) max_y_legend=0.85;

	TLegend* legend;
	legend = new TLegend(min_x_legend,max_y_legend-delta_y_legend,max_x_legend,max_y_legend);
	legend->SetFillColor(0);
	legend->SetTextSize(0.04);
	legend->SetBorderSize(0);
	legend->AddEntry(data_Graph,Form("%s", ExpName[iExperiment]),"lp");
	if(iMeasurementID==0) if(!HPbool) legend->AddEntry(data_Graph_nopolcorr,Form("%s %s (#vec{#lambda}=0)", ExpName[iExperiment], StateNameTex[iState]),"p");

	if(PlotModel){
		if(model_Graph->GetN()>1) legend->AddEntry(model_Graph,"Total","l");
		else legend->AddEntry(model_Graph,"Total","l");
	}

	TGraph* model_Graph_ColorChannels_neg[ColorChannels_c];
	TGraph* model_Graph_ColorChannels_pos[ColorChannels_c];
	TGraphAsymmErrors* model_Graph_ColorChannels_Bands_neg[3][ColorChannels_c];
	TGraphAsymmErrors* model_Graph_ColorChannels_Bands_pos[3][ColorChannels_c];

	bool PlotNegChannels[ColorChannels_c];
	bool PlotPosChannels[ColorChannels_c];
	if(plotDirectColorChannels&&PlotModel){
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
					if(buffx<pTMinModel){
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
					if(buffx<pTMinModel){
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
				cout<<"model_Graph_ColorChannels_original["<<i<<"]:"<<endl;
				model_Graph_ColorChannels[i]->Print();


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
				if(plotColorChannel[i]) model_Graph_ColorChannels[i]->Draw("lsame");
				if(singlePointModel){
					double buffModel_pT, buffModel_res;
					model_Graph_ColorChannels[i]->GetPoint(0, buffModel_pT, buffModel_res);
					TLine* modelReplacement = new TLine( buff_pT-buff_errpT_low, buffModel_res, buff_pT+buff_errpT_high, buffModel_res);
					modelReplacement->SetLineColor(model_Graph_ColorChannels[i]->GetLineColor());
					modelReplacement->SetLineStyle(model_Graph_ColorChannels[i]->GetLineStyle());
					modelReplacement->SetLineWidth(model_Graph_ColorChannels[i]->GetLineWidth());
					if(plotColorChannel[i]) modelReplacement->Draw( "same" );
				}

				cout<<"model_Graph_ColorChannels["<<i<<"]:"<<endl;
				model_Graph_ColorChannels[i]->Print();
				//if(model_Graph_ColorChannels[i]->GetN()>1){
					//if(isSstate) legend->AddEntry(model_Graph_ColorChannels[i],Form("Direct production, O_{%s}^{%s}", ColorChannelNameTexS[i], StateNameTex[iState]),"l");
					//else legend->AddEntry(model_Graph_ColorChannels[i],Form("Direct production, O_{%s}^{%s}", ColorChannelNameTexP[i], StateNameTex[iState]),"l");
					if(isSstate) legend->AddEntry(model_Graph_ColorChannels[i],Form("%s", ColorChannelNameTexS[i]),"l");
					else legend->AddEntry(model_Graph_ColorChannels[i],Form("%s", ColorChannelNameTexP[i]),"l");
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
							model_Graph_ColorChannels_Bands[iSig][i]->Print();
							if(plotColorChannelBands[i] && plotSigmaBands[iSig]) model_Graph_ColorChannels_Bands[iSig][i]->Draw("3same");
						}
						if(plotColorChannel[i]) model_Graph_ColorChannels[i]->Draw("lsame");
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
							TLine* modelReplacement = new TLine( buff_pT-buff_errpT_low, buffModel_res, buff_pT+buff_errpT_high, buffModel_res);
							modelReplacement->SetLineColor(model_Graph_ColorChannels[i]->GetLineColor());
							modelReplacement->SetLineStyle(model_Graph_ColorChannels[i]->GetLineStyle());
							modelReplacement->SetLineWidth(model_Graph_ColorChannels[i]->GetLineWidth());
							if(plotColorChannel[i]) modelReplacement->Draw( "same" );

						}

					//}

			}
			if(PlotNegChannels[i]){
				if(PlotPosChannels[i]){
					for(int iSig=2;iSig>-1;iSig--){
						if(plotColorChannelBands[i] && plotSigmaBands[iSig]) model_Graph_ColorChannels_Bands_pos[iSig][i]->Draw("3same");
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
						TLine* modelReplacement = new TLine( buff_pT-buff_errpT_low, buffModel_res, buff_pT+buff_errpT_high, buffModel_res);
						modelReplacement->SetLineColor(model_Graph_ColorChannels_pos[i]->GetLineColor());
						modelReplacement->SetLineStyle(model_Graph_ColorChannels_pos[i]->GetLineStyle());
						modelReplacement->SetLineWidth(model_Graph_ColorChannels_pos[i]->GetLineWidth());
						if(plotColorChannel[i]) modelReplacement->Draw( "same" );

					}

					cout<<"model_Graph_ColorChannels_pos["<<i<<"]:"<<endl;
					model_Graph_ColorChannels_pos[i]->Print();
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
					if(plotColorChannelBands[i] && plotSigmaBands[iSig]) model_Graph_ColorChannels_Bands_neg[iSig][i]->Draw("3same");
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
					TLine* modelReplacement = new TLine( buff_pT-buff_errpT_low, buffModel_res, buff_pT+buff_errpT_high, buffModel_res);
					modelReplacement->SetLineColor(model_Graph_ColorChannels_neg[i]->GetLineColor());
					modelReplacement->SetLineStyle(model_Graph_ColorChannels_neg[i]->GetLineStyle());
					modelReplacement->SetLineWidth(model_Graph_ColorChannels_neg[i]->GetLineWidth());
					if(plotColorChannel[i]) modelReplacement->Draw( "same" );

				}
				cout<<"model_Graph_ColorChannels_neg["<<i<<"]:"<<endl;
				model_Graph_ColorChannels_neg[i]->Print();
				//if(isSstate) legend->AddEntry(model_Graph_ColorChannels_neg[i],Form("(neg.) Direct production, O_{%s}^{%s}", ColorChannelNameTexS[i], StateNameTex[iState]),"l");
				//else legend->AddEntry(model_Graph_ColorChannels_neg[i],Form("(neg.) Direct production, O_{%s}^{%s}", ColorChannelNameTexP[i], StateNameTex[iState]),"l");
				if(isSstate) legend->AddEntry(model_Graph_ColorChannels_neg[i],Form("%s (neg.)", ColorChannelNameTexS[i]),"l");
				else legend->AddEntry(model_Graph_ColorChannels_neg[i],Form("%s (neg.)", ColorChannelNameTexP[i]),"l");

			}







		}
	}

	if(plotIndividualFeedDown){
		//cout<<"StatesCont_c "<<StatesCont_c<<endl;
		for(int i=1;i<StatesCont_c;i++){
			//cout<<"StateNameTex[StatesContributing[i]] "<<StateNameTex[StatesContributing[i]]<<endl;
			model_Graph_FeedDowns[i]->SetLineColor(colorFD[i]);
			model_Graph_FeedDowns[i]->SetLineStyle(linestyleFD);
			model_Graph_FeedDowns[i]->SetLineWidth(linewidthFD);
			model_Graph_FeedDowns[i]->Draw("lsame");
			legend->AddEntry(model_Graph_FeedDowns[i],Form("Feed-down from %s", StateNameTex[StatesContributing[i]]),"l");
		}
	}

	if(plotInclusiveFeedDown){
		//cout<<"draw inclusive feed-down"<<endl;
		// inclusive feed-down
		//cout<<"model_Graph_FeedDowns: "<<0<<endl;
		//model_Graph_FeedDowns[0]->Print();
		model_Graph_FeedDowns[0]->SetLineColor(colorFD[0]);
		model_Graph_FeedDowns[0]->SetLineStyle(linestyleFD);
		model_Graph_FeedDowns[0]->SetLineWidth(linewidthFD*2.);
		model_Graph_FeedDowns[0]->Draw("lsame");
		legend->AddEntry(model_Graph_FeedDowns[0],"Inclusive feed-down","l");
	}

	if(PlotModel){
		for(int iSig=2;iSig>-1;iSig--){
			if(plotTotalModelBands && plotSigmaBands[iSig]) model_Graph_Bands[iSig]->Draw("3same");
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
			if(plotTotalModel) modelReplacement->Draw( "same" );
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
		model_CS_THuncBand->Draw("3same");
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
				if(plotCSTHuncLine[iTHunc]) modelReplacement->Draw( "same" );

			}

		}



	}

	if(plotDirectColorChannels&&PlotModel){
	//Replot central color channels, to assure that they are on top of the total model bands... (!!!)
		for(int iColorChannelSequence=0;iColorChannelSequence<ColorChannels_c;iColorChannelSequence++){
			int i=ColorChannelSequence[iColorChannelSequence];
			if(!PlotNegChannels[i]){
				if(plotColorChannel[i]) model_Graph_ColorChannels[i]->Draw("lsame");
				if(singlePointModel){
					double buffModel_pT, buffModel_res;
					model_Graph_ColorChannels[i]->GetPoint(0, buffModel_pT, buffModel_res);
					TLine* modelReplacement = new TLine( buff_pT-buff_errpT_low, buffModel_res, buff_pT+buff_errpT_high, buffModel_res);
					modelReplacement->SetLineColor(model_Graph_ColorChannels[i]->GetLineColor());
					modelReplacement->SetLineStyle(model_Graph_ColorChannels[i]->GetLineStyle());
					modelReplacement->SetLineWidth(model_Graph_ColorChannels[i]->GetLineWidth());
					if(plotColorChannel[i]) modelReplacement->Draw( "same" );
				}
			}
			if(PlotNegChannels[i]){
				if(PlotPosChannels[i]){
					if(plotColorChannel[i]) model_Graph_ColorChannels_pos[i]->Draw("lsame");
					if(singlePointModel){
						double buffModel_pT, buffModel_res, buffModel_reserrlow, buffModel_reserrhigh;
						TLine* modelReplacement = new TLine( buff_pT-buff_errpT_low, buffModel_res, buff_pT+buff_errpT_high, buffModel_res);
						modelReplacement->SetLineColor(model_Graph_ColorChannels_pos[i]->GetLineColor());
						modelReplacement->SetLineStyle(model_Graph_ColorChannels_pos[i]->GetLineStyle());
						modelReplacement->SetLineWidth(model_Graph_ColorChannels_pos[i]->GetLineWidth());
						if(plotColorChannel[i]) modelReplacement->Draw( "same" );

					}
				}
				model_Graph_ColorChannels_neg[i]->Draw("lsame");
				if(singlePointModel){
					double buffModel_pT, buffModel_res, buffModel_reserrlow, buffModel_reserrhigh;
					TLine* modelReplacement = new TLine( buff_pT-buff_errpT_low, buffModel_res, buff_pT+buff_errpT_high, buffModel_res);
					modelReplacement->SetLineColor(model_Graph_ColorChannels_neg[i]->GetLineColor());
					modelReplacement->SetLineStyle(model_Graph_ColorChannels_neg[i]->GetLineStyle());
					modelReplacement->SetLineWidth(model_Graph_ColorChannels_neg[i]->GetLineWidth());
					if(plotColorChannel[i]) modelReplacement->Draw( "same" );

				}
			}
		}
	}

	data_Graph->Draw("psame");

	axishist->Draw("AXISsame");

	if(PlotModel) legend->Draw("same");

	////draw latex
	double left=0.725, top=0.1875, textSize=0.03625;
	TLatex *latex=new TLatex();
	latex->SetTextFont(42);
	latex->SetNDC(kTRUE);
	latex->SetTextSize(textSize);
	//double stepLatex=textSize*1.3;

	if(chi2Min>1e-2 && chi2Min<998) latex->DrawLatex(left,top, Form("#chi^{2} / ndf = %1.1f / %d", chi2Min, ndf));
	top=0.25;
	if(chi2Min>1e-2 && chi2Min<998) latex->DrawLatex(left,top, Form("P(#chi^{2}, ndf) = %1.2G", chi2Prob));

	if(logY) if(iMeasurementID==0) plotCanvas->SetLogy(true);

	left=0.725; top=0.895;
	if(longrapchar) left=0.68125;
	textSize=0.05;
	latex->SetTextFont(42);
	latex->SetTextSize(textSize);
	latex->DrawLatex(left,top, Form("%s, %s", StateNameTex[iState], rapchar));

	if(!PlotModel){
		left=0.7875; top=0.825;
		textSize=0.05;
		latex->SetTextFont(42);
		latex->SetTextSize(textSize);
		latex->DrawLatex(left,top, Form("%s", ExpName[iExperiment]));

	}

	left=0.03; top=0.55;
	textSize=0.08;
	latex->SetTextFont(42);
	latex->SetTextSize(textSize);
	if(iMeasurementID==1) latex->DrawLatex(left,top, "#lambda_{#vartheta}^{#scale[0.7]{HX}}");
	if(iMeasurementID==2) latex->DrawLatex(left,top, "#lambda_{#varphi}^{#scale[0.7]{HX}}");
	if(iMeasurementID==3) latex->DrawLatex(left,top, "#lambda_{#vartheta#varphi}^{#scale[0.7]{HX}}");

	char savename[1000];
	sprintf(savename,"%s/Figures/DataModelComp_%s_%s_%s_rap%d.pdf",jobdirname, ExpName[iExperiment], StateName[iState], MeasurementIDName[iMeasurementID], iRap+1);
	cout<<"saved here...: "<<savename<<endl;
	plotCanvas->SaveAs(savename);
	sprintf(savename,"%s/Figures/PNG_DataModelComp_%s_%s_%s_rap%d.png",jobdirname, ExpName[iExperiment], StateName[iState], MeasurementIDName[iMeasurementID], iRap+1);
	cout<<"... and saved here: "<<savename<<endl;
	plotCanvas->SaveAs(savename);



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


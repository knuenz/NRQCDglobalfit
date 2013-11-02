/*
 * CombineDataModel.cc
 *
 *  Created on: Jun 30, 2013
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


using namespace NRQCDvars;


vector<double> PolarizationTransfer(vector<double> lamVecOriginal, vector<int> DecayChain);
vector<double> addPolarizations(vector<vector<double> > lamMatrix, vector<double> lamVecContributions);
vector<double> GeneralPolarizationTransformation(vector<double> lamVecOriginal, double t1, double t2, double t3);

int main(int argc, char** argv) {

  	Char_t *ModelID = "Default";
  	Char_t *DataID = "Default";
  	Char_t *DataModelCombID = "Default";
  	Char_t *ModelSystScaleID = "Default";

	Char_t *storagedir = "Default"; //Storage Directory

  	for( int i=0;i < argc; ++i ) {
		if(std::string(argv[i]).find("ModelID") != std::string::npos) {char* ModelIDchar = argv[i]; char* ModelIDchar2 = strtok (ModelIDchar, "="); ModelID = ModelIDchar2; cout<<"ModelID = "<<ModelID<<endl;}
		if(std::string(argv[i]).find("DataID") != std::string::npos) {char* DataIDchar = argv[i]; char* DataIDchar2 = strtok (DataIDchar, "="); DataID = DataIDchar2; cout<<"DataID = "<<DataID<<endl;}
		if(std::string(argv[i]).find("DataModelCombID") != std::string::npos) {char* DataModelCombIDchar = argv[i]; char* DataModelCombIDchar2 = strtok (DataModelCombIDchar, "="); DataModelCombID = DataModelCombIDchar2; cout<<"DataModelCombID = "<<DataModelCombID<<endl;}
		if(std::string(argv[i]).find("ModelSystScaleID") != std::string::npos) {char* ModelSystScaleIDchar = argv[i]; char* ModelSystScaleIDchar2 = strtok (ModelSystScaleIDchar, "="); ModelSystScaleID = ModelSystScaleIDchar2; cout<<"ModelSystScaleID = "<<ModelSystScaleID<<endl;}
		if(std::string(argv[i]).find("storagedir") != std::string::npos) {char* storagedirchar = argv[i]; char* storagedirchar2 = strtok (storagedirchar, "="); storagedir = storagedirchar2; cout<<"storagedir = "<<storagedir<<endl;}
  	}


  	gROOT->SetBatch();

	delete gRandom;
	gRandom = new TRandom3(NRQCDvars::randomSeed);  // better random generator


	// TODO: add loop for nSystematics (now nSystematics will be 0)

	char outname[200];
	char inname[200];
	char datadirname[200];
	char modeldirname[200];
	char datamodeldirname[200];
	char modelsystematicscalesdirname[200];
	sprintf(modeldirname,"%s/ModelID", storagedir);
	gSystem->mkdir(modeldirname);
	sprintf(modeldirname,"%s/%s",modeldirname,ModelID);
	gSystem->mkdir(modeldirname);
	sprintf(datadirname,"%s/DataID", storagedir);
	gSystem->mkdir(datadirname);
	sprintf(datadirname,"%s/%s",datadirname,DataID);
	gSystem->mkdir(datadirname);
	sprintf(datamodeldirname,"%s/DataModelCombID", storagedir);
	gSystem->mkdir(datamodeldirname);
	sprintf(datamodeldirname,"%s/%s",datamodeldirname,DataModelCombID);
	gSystem->mkdir(datamodeldirname);
	sprintf(modelsystematicscalesdirname,"%s/ModelSystScaleID", storagedir);
	gSystem->mkdir(modelsystematicscalesdirname);
	sprintf(modelsystematicscalesdirname,"%s/%s",modelsystematicscalesdirname,ModelSystScaleID);
	gSystem->mkdir(modelsystematicscalesdirname);


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Read in NRQCDglobalfitObjects, Model input: nTuples of previous step
// Loop through NRQCDglobalfitObjects, calculate predictions specific for the MeasurementID, pT/rapidity bins and nState
// Fill NRQCDglobalfitObjects with model and model uncertainties, save into output file
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//SDC: number of events in pT rap bin, weighted
	//Polarization: project lambdas, for each decayChannel of iMother-> iDaughter do all transformations (simple function: from Pietro), sum up contributions with Branchingratiofractions.


  	sprintf(inname,"%s/ModelIngredients.root",modeldirname);
	TFile* ModelIngredientsFile = new TFile(inname, "READ");




	for(int iState=0;iState<NRQCDvars::nStates;iState++){
		for(int iMeasurementID=0;iMeasurementID<NRQCDvars::nMeasurementIDs;iMeasurementID++){
			for(int iExperiment=0;iExperiment<NRQCDvars::nExperiments;iExperiment++){
//				for(int iRap = 0; iRap < NRQCDvars::nMaxRapBins; iRap++){
//				    for(int iP = 0; iP < NRQCDvars::nMaxPtBins; iP++){

				    	//if(iState!=3 || iMeasurementID!=0 || iExperiment!=4 /*|| iRap!=0 || iP!=0*/) continue;
				    	cout<<"calculate and set model for iState="<<StateName[iState]<<", iMeasurementID="<<MeasurementIDName[iMeasurementID]<<", iExperiment="<<ExpName[iExperiment]/*<<", iRap="<<iRap<<", iP="<<iP*/<<endl;


				    	NRQCDglobalfitObject *DataObject[NRQCDvars::nMaxRapBins][NRQCDvars::nMaxPtBins];
				    	NRQCDglobalfitObject *ModelObject[NRQCDvars::nMaxRapBins][NRQCDvars::nMaxPtBins];
				    	bool doesDataExist[NRQCDvars::nMaxRapBins][NRQCDvars::nMaxPtBins];
						bool doesAnyDataExist=false;

						//cout<<"Reading in all DataObjects (for all rap/pT bins)"<<endl;
						for(int iRap = 0; iRap < NRQCDvars::nMaxRapBins; iRap++){
						    for(int iP = 0; iP < NRQCDvars::nMaxPtBins; iP++){

								//get data input file
								sprintf(inname,"%s/ConvertedData_%s_%s_%s_rap%d_pT%d.txt",datadirname, StateName[iState],  MeasurementIDName[iMeasurementID],  ExpName[iExperiment], iRap, iP);
								ifstream in;
								in.open(inname);

								doesDataExist[iRap][iP]=false;
								//read NRQCDglobalfitObject members of data file
								if(in!=NULL){
									doesAnyDataExist=true;
									doesDataExist[iRap][iP]=true;
									DataObject[iRap][iP] = new NRQCDglobalfitObject();
									in >> *DataObject[iRap][iP];
									in.close();

									//Clone DataObject -> ModelObject and set Data members correctly
									ModelObject[iRap][iP] = new NRQCDglobalfitObject();
									ModelObject[iRap][iP]->setData(DataObject[iRap][iP]->getState(),
											DataObject[iRap][iP]->getStateDenom(),
											DataObject[iRap][iP]->getCentralValue(),
											DataObject[iRap][iP]->getErrStatPos(),
											DataObject[iRap][iP]->getErrStatNeg(),
											DataObject[iRap][iP]->getErrSystPos(),
											DataObject[iRap][iP]->getErrSystNeg(),
											DataObject[iRap][iP]->getErrTotPos(),
											DataObject[iRap][iP]->getErrTotNeg(),
											DataObject[iRap][iP]->getErrGlobalPos(),
											DataObject[iRap][iP]->getErrGlobalNeg(),
											DataObject[iRap][iP]->getyMin(),
											DataObject[iRap][iP]->getyMax(),
											DataObject[iRap][iP]->getyMean(),
											DataObject[iRap][iP]->getpTMin(),
											DataObject[iRap][iP]->getpTMax(),
											DataObject[iRap][iP]->getpTMean(),
											DataObject[iRap][iP]->getPolCorrParams(),
											DataObject[iRap][iP]->getMeasurementID(),
											DataObject[iRap][iP]->getExperiment(),
											DataObject[iRap][iP]->getObjectID(),
											DataObject[iRap][iP]->getisDataValid(),
											DataObject[iRap][iP]->getisAbsRap());
								}

						    }
						}
						//cout<<"End reading in all DataObjects (for all rap/pT bins)"<<endl;

						if(!doesAnyDataExist) {cout<<"Measurement does not exist"<<endl; continue;}

							double Dummy_SDC=999.;
							double Dummy_Lam=0.0999;

							for(int iMother=0;iMother<NRQCDvars::nStates;iMother++){

								vector<double>setModel_ShortDistanceCoef[NRQCDvars::nMaxRapBins][NRQCDvars::nMaxPtBins];
								vector<double>setModel_OctetCompLamth[NRQCDvars::nMaxRapBins][NRQCDvars::nMaxPtBins];
								vector<double>setModel_OctetCompLamph[NRQCDvars::nMaxRapBins][NRQCDvars::nMaxPtBins];
								vector<double>setModel_OctetCompLamtp[NRQCDvars::nMaxRapBins][NRQCDvars::nMaxPtBins];
								vector<double>setModel_ShortDistanceCoef_PointByPointSyst[NRQCDvars::nMaxRapBins][NRQCDvars::nMaxPtBins];
								vector<double>setModel_OctetCompLamth_PointByPointSyst[NRQCDvars::nMaxRapBins][NRQCDvars::nMaxPtBins];
								vector<double>setModel_OctetCompLamph_PointByPointSyst[NRQCDvars::nMaxRapBins][NRQCDvars::nMaxPtBins];
								vector<double>setModel_OctetCompLamtp_PointByPointSyst[NRQCDvars::nMaxRapBins][NRQCDvars::nMaxPtBins];
								vector<vector<double> >setModel_ShortDistanceCoef_globalSystPos[NRQCDvars::nMaxRapBins][NRQCDvars::nMaxPtBins];
								vector<vector<double> >setModel_ShortDistanceCoef_globalSystNeg[NRQCDvars::nMaxRapBins][NRQCDvars::nMaxPtBins];
								vector<vector<double> >setModel_OctetCompLamth_globalSystPos[NRQCDvars::nMaxRapBins][NRQCDvars::nMaxPtBins];
								vector<vector<double> >setModel_OctetCompLamth_globalSystNeg[NRQCDvars::nMaxRapBins][NRQCDvars::nMaxPtBins];
								vector<vector<double> >setModel_OctetCompLamph_globalSystPos[NRQCDvars::nMaxRapBins][NRQCDvars::nMaxPtBins];
								vector<vector<double> >setModel_OctetCompLamph_globalSystNeg[NRQCDvars::nMaxRapBins][NRQCDvars::nMaxPtBins];
								vector<vector<double> >setModel_OctetCompLamtp_globalSystPos[NRQCDvars::nMaxRapBins][NRQCDvars::nMaxPtBins];
								vector<vector<double> >setModel_OctetCompLamtp_globalSystNeg[NRQCDvars::nMaxRapBins][NRQCDvars::nMaxPtBins];

								cout<<"iMother "<<iMother<<endl;
								int nColorChannels_state;
								bool isSstate=(StateQuantumID[iMother] > NRQCDvars::quID_S)?false:true;
								if(isSstate) nColorChannels_state=NRQCDvars::nColorChannels_S;
								else nColorChannels_state=NRQCDvars::nColorChannels_P;
								for (int iColorChannel=0; iColorChannel<nColorChannels_state; iColorChannel++){

									bool ContributionExists=false;
									cout<<"iColorChannel "<<iColorChannel<<endl;

									//Calculation of short distance coefficients for each feed-down component + direct production
									double Buffer_ShortDistanceCoef[NRQCDvars::nMaxRapBins][NRQCDvars::nMaxPtBins];
									vector<vector<double> > lamVecTransformedCascades[NRQCDvars::nMaxRapBins][NRQCDvars::nMaxPtBins];
									vector<double> lamVecTransformedCascadesContributions[NRQCDvars::nMaxRapBins][NRQCDvars::nMaxPtBins];
									//cout<<"..."<<endl;

									for (int iCascade=0; iCascade<NRQCDvars::nMaxCascades; iCascade++){//Loop over all different decay chains linking iMother->iDaughter
										//cout<<"iCascade "<<iCascade<<endl;
										char nTupleModelName[1000];
										sprintf(nTupleModelName,"nTupleModel_Mother%d_Daughter%d_ColorChannel%d_Cascade%d",iMother, iState, iColorChannel, iCascade);
										TTree* nTupleModel_Transformed = (TTree*)ModelIngredientsFile->Get(nTupleModelName);
										//cout<<"nTupleModelName "<<nTupleModelName<<endl;
										//cout<<"ContributionExists "<<ContributionExists<<endl;
										if(nTupleModel_Transformed!=NULL){
											ContributionExists=true;
											//cout<<"opened file"<<endl;
										//SDC:
											//char nTupleModelSelection[1000];
											//sprintf(nTupleModelSelection,"model_pT > %f && model_pT < %f && model_rap > %f && model_rap < %f", ModelObject->getpTMin(), ModelObject->getpTMax(), ModelObject->getyMin(), ModelObject->getyMax());
											//Buffer_ShortDistanceCoef+=nTupleModel_Transformed->GetEntries(nTupleModelSelection)*nTupleModel_Transformed->GetWeight();
											//POL:
											const double l_min = -1, l_max =  1, l_step_1D = 0.02;

											TH1D* h_costh2[NRQCDvars::nMaxRapBins][NRQCDvars::nMaxPtBins];
											TH1D* h_cos2ph[NRQCDvars::nMaxRapBins][NRQCDvars::nMaxPtBins];
											TH1D* h_sin2thcosph[NRQCDvars::nMaxRapBins][NRQCDvars::nMaxPtBins];


											//nTupleModel_Transformed->Draw("model_costh*model_costh>>h_costh2",nTupleModelSelection);
											//nTupleModel_Transformed->Draw("TMath::Cos(2*model_phi*TMath::Pi()/180)>>h_cos2ph",nTupleModelSelection);
											//nTupleModel_Transformed->Draw("TMath::Sin(2*TMath::ACos(model_costh))*TMath::Cos(model_phi*TMath::Pi()/180)>>h_sin2thcosph",nTupleModelSelection);


											// SetBranches of original nTuple
											double model_pT_mother; nTupleModel_Transformed->SetBranchAddress( "model_pT",  &model_pT_mother );
											double model_rap_mother; nTupleModel_Transformed->SetBranchAddress( "model_rap",  &model_rap_mother );
											double model_costh_mother; nTupleModel_Transformed->SetBranchAddress( "model_costh",  &model_costh_mother );
											double model_phi_mother; nTupleModel_Transformed->SetBranchAddress( "model_phi",  &model_phi_mother );
											double weight_mother; nTupleModel_Transformed->SetBranchAddress( "weight",  &weight_mother );

											for(int iRap = 0; iRap < NRQCDvars::nMaxRapBins; iRap++){
											    for(int iP = 0; iP < NRQCDvars::nMaxPtBins; iP++){
											    	if(!doesDataExist[iRap][iP]) continue;
											    	Buffer_ShortDistanceCoef[iRap][iP]=0;
											    	char histname[200];
											    	sprintf(histname,"h_costh2_iRap%d_iP%d",iRap,iP);
													h_costh2[iRap][iP] = new TH1D( histname, "", int((l_max-l_min)/l_step_1D), l_min, l_max );
											    	sprintf(histname,"h_cos2ph_iRap%d_iP%d",iRap,iP);
													h_cos2ph[iRap][iP] = new TH1D( histname, "", int((l_max-l_min)/l_step_1D), l_min, l_max );
											    	sprintf(histname,"h_sin2thcosph_iRap%d_iP%d",iRap,iP);
													h_sin2thcosph[iRap][iP] = new TH1D( histname, "", int((l_max-l_min)/l_step_1D), l_min, l_max );
											    }
											}

											int n_events = int( nTupleModel_Transformed->GetEntries() );
											cout<<n_events<<" events: Looping through nTuple of iMother="<<iMother<<", iColorChannel="<<iColorChannel<<", iCascade="<<iCascade<<endl;

											TH1F* h_weight_pT3 = new TH1F("h_weight_pT3","h_weight_pT3",100,-4000,4000);
											TH1F* h_weight_pT4 = new TH1F("h_weight_pT4","h_weight_pT4",100,-4000,4000);
											TH1F* h_weight_pT5 = new TH1F("h_weight_pT5","h_weight_pT5",100,-4000,4000);
											TH1F* h_weight_pT6 = new TH1F("h_weight_pT6","h_weight_pT6",100,-4000,4000);

											int i_eventPerBin[NRQCDvars::nMaxRapBins][NRQCDvars::nMaxPtBins];
											for(int iRap = 0; iRap < NRQCDvars::nMaxRapBins; iRap++){
											    for(int iP = 0; iP < NRQCDvars::nMaxPtBins; iP++){
											    	i_eventPerBin[iRap][iP]=0;
											    }
											}
											// loop over  events in the model ntuple
											for ( int i_event = 1; i_event <= n_events; i_event++ ) {

												nTupleModel_Transformed->GetEvent( i_event-1 );

												//cout<<"i_event="<<i_event<<endl;


												int pT_index=-1;
												int rap_index=-1;
												for(int iRap = 0; iRap < NRQCDvars::nMaxRapBins; iRap++){
												    for(int iP = 0; iP < NRQCDvars::nMaxPtBins; iP++){
												    	if(!doesDataExist[iRap][iP]){
												    		continue;
												    	}
												    	else if(!DataObject[iRap][iP]->getisAbsRap() && (model_pT_mother > ModelObject[iRap][iP]->getpTMin() && model_pT_mother < ModelObject[iRap][iP]->getpTMax() && model_rap_mother > ModelObject[iRap][iP]->getyMin() && model_rap_mother < ModelObject[iRap][iP]->getyMax())){
															rap_index=iRap;
															pT_index=iP;
														}
												    	else if(DataObject[iRap][iP]->getisAbsRap() && (model_pT_mother > ModelObject[iRap][iP]->getpTMin() && model_pT_mother < ModelObject[iRap][iP]->getpTMax() && fabs(model_rap_mother) > ModelObject[iRap][iP]->getyMin() && fabs(model_rap_mother) < ModelObject[iRap][iP]->getyMax())){
															rap_index=iRap;
															pT_index=iP;
														}
														if(pT_index>-1 || rap_index>-1) continue;
												    }
													if(pT_index>-1 || rap_index>-1) continue;
												}

												if(pT_index<0 || rap_index<0) continue;

												if(rap_index==0){
													if(pT_index==3){
														h_weight_pT3->Fill(weight_mother);
													}
													if(pT_index==4){
														h_weight_pT4->Fill(weight_mother);
													}
													if(pT_index==5){
														h_weight_pT5->Fill(weight_mother);
													}
													if(pT_index==6){
														h_weight_pT6->Fill(weight_mother);
													}
												}

												i_eventPerBin[rap_index][pT_index]++;
												Buffer_ShortDistanceCoef[rap_index][pT_index]+=weight_mother;
												h_costh2[rap_index][pT_index]->Fill(model_costh_mother*model_costh_mother, weight_mother);
												h_cos2ph[rap_index][pT_index]->Fill(TMath::Cos(2*model_phi_mother*TMath::Pi()/180), weight_mother);
												h_sin2thcosph[rap_index][pT_index]->Fill(TMath::Sin(2*TMath::ACos(model_costh_mother))*TMath::Cos(model_phi_mother*TMath::Pi()/180), weight_mother);

											}

									    	if(iExperiment==0&&iMeasurementID==0&&iColorChannel==3){
												h_weight_pT3->SaveAs("tmp/h_weight_pT3.root");
												h_weight_pT4->SaveAs("tmp/h_weight_pT4.root");
												h_weight_pT5->SaveAs("tmp/h_weight_pT5.root");
												h_weight_pT6->SaveAs("tmp/h_weight_pT6.root");
									    	}
											//cout<<"Calculating SDC and vecLam of iMother="<<iMother<<", iColorChannel="<<iColorChannel<<", iCascade="<<iCascade<<endl;

											for(int iRap = 0; iRap < NRQCDvars::nMaxRapBins; iRap++){
											    for(int iP = 0; iP < NRQCDvars::nMaxPtBins; iP++){
											    	if(!doesDataExist[iRap][iP]) continue;
											    	if(i_eventPerBin[iRap][iP]<1) continue;


											    	if(iExperiment==0&&iMeasurementID==0&&iColorChannel==3){
											    		cout<<"Buffer_ShortDistanceCoef["<<iRap<<"]["<<iP<<"]"<<Buffer_ShortDistanceCoef[iRap][iP]<<endl;
											    	}

											    	Buffer_ShortDistanceCoef[iRap][iP]*=nTupleModel_Transformed->GetWeight();



											    	double deltaPt=DataObject[iRap][iP]->getpTMax()-DataObject[iRap][iP]->getpTMin();
											    	double deltay=DataObject[iRap][iP]->getyMax()-DataObject[iRap][iP]->getyMin();
											    	if(DataObject[iRap][iP]->getisAbsRap()) deltay*=2;
											    	Buffer_ShortDistanceCoef[iRap][iP]/=deltaPt*deltay;


											vector<double> lamVecRaw(3,0);
											double costh2=h_costh2[iRap][iP]->GetMean();
											lamVecRaw[0] = (1. - 3. * costh2 ) / ( costh2 - 3./5. );
											double cos2ph=h_cos2ph[iRap][iP]->GetMean();
											lamVecRaw[1] = cos2ph * (3. + lamVecRaw[0]);
											double sin2thcosph=h_sin2thcosph[iRap][iP]->GetMean();
											lamVecRaw[2] = sin2thcosph * 5./4. * (3. + lamVecRaw[0]);

											vector<int> DecayChain;
											if(iMother!=iState){
												char DecayChainIDName[1000];
												sprintf(DecayChainIDName,"DecayChainID_Mother%d_Daughter%d_ColorChannel%d_Cascade%d",iMother, iState, iColorChannel, iCascade);
												TGraph* DecayChainID = (TGraph*)ModelIngredientsFile->Get(DecayChainIDName);
												int nCascadeDecays=DecayChainID->GetN();
												for (int iCascadeDecay=0; iCascadeDecay<nCascadeDecays; iCascadeDecay++){//Loop over all different decay chains linking iMother->iDaughter
													double BufferX, BufferY;
													DecayChainID->GetPoint(iCascadeDecay, BufferX, BufferY);
													DecayChain.push_back(int(BufferY));
												}
												//if(iColorChannel==0) cout<<"PolarizationTransfer of "<<DecayChainIDName<<endl;
											}
											else{
												DecayChain.push_back(iState);
												//if(iColorChannel==0) cout<<"PolarizationTransfer of direct state"<<iState<<endl;
											}

											//vector<int> DecayChain2;
											//DecayChain2.push_back(7);
											//DecayChain2.push_back(6);
											//DecayChain2.push_back(4);

											vector<double> lamVecTransformed=PolarizationTransfer(lamVecRaw,DecayChain);

											lamVecTransformedCascades[iRap][iP].push_back(lamVecTransformed);
											lamVecTransformedCascadesContributions[iRap][iP].push_back(Buffer_ShortDistanceCoef[iRap][iP]);

											delete h_costh2[iRap][iP];
											delete h_cos2ph[iRap][iP];
											delete h_sin2thcosph[iRap][iP];

											    }
											}

										}
									}
									//cout<<"here...?"<<endl;
									//cout<<"ContributionExists "<<ContributionExists<<endl;

									//Calculation of the polarization components of each feed-down component + direct production
									//Strategy:
										//Get lambdas of direct production (via projections of costh/phi of original nTuple in kinematic region)
										//Applying Pietro's polarization transfer formulas, calculate lambdas for each cascade of iMother->iDaughter -> sum up the lambdas with given weights


									for(int iRap = 0; iRap < NRQCDvars::nMaxRapBins; iRap++){
									    for(int iP = 0; iP < NRQCDvars::nMaxPtBins; iP++){
									    	if(!doesDataExist[iRap][iP]) continue;
									    	cout<<"iRap "<<iRap<<endl;
									    	cout<<"iP "<<iP<<endl;

									    	//cout<<"here?"<<endl;

									vector<double> lamVecSumCascades;
									if(ContributionExists){

										vector<double> BuffVec=lamVecTransformedCascadesContributions[iRap][iP];
										int nContributions=BuffVec.size();
										double ContributionIntegral=0;
										for(int i=0; i<nContributions; i++){
											ContributionIntegral+=BuffVec[i];
										}
										if(ContributionIntegral==0) {doesDataExist[iRap][iP]=false; continue;}
										else lamVecSumCascades = addPolarizations(lamVecTransformedCascades[iRap][iP], lamVecTransformedCascadesContributions[iRap][iP]);
										//// These are the only values that need to be set to fully define the model:
										setModel_ShortDistanceCoef[iRap][iP].push_back(Buffer_ShortDistanceCoef[iRap][iP]);
										setModel_OctetCompLamth[iRap][iP].push_back(lamVecSumCascades[0]);
										setModel_OctetCompLamph[iRap][iP].push_back(lamVecSumCascades[1]);
										setModel_OctetCompLamtp[iRap][iP].push_back(lamVecSumCascades[2]);
										/////////////////////////////////////////////////////////////////////////////
									}
									else{
										//cout<<"outthere?"<<endl;
										setModel_ShortDistanceCoef[iRap][iP].push_back(Dummy_SDC);
										setModel_OctetCompLamth[iRap][iP].push_back(Dummy_Lam);
										setModel_OctetCompLamph[iRap][iP].push_back(Dummy_Lam);
										setModel_OctetCompLamtp[iRap][iP].push_back(Dummy_Lam);

									}

									//cout<<"here?"<<endl;
									setModel_ShortDistanceCoef_PointByPointSyst[iRap][iP].push_back(Dummy_SDC);//Not implemented so far -> dummy value
									setModel_OctetCompLamth_PointByPointSyst[iRap][iP].push_back(Dummy_Lam);//Not implemented so far -> dummy value
									setModel_OctetCompLamph_PointByPointSyst[iRap][iP].push_back(Dummy_Lam);//Not implemented so far -> dummy value
									setModel_OctetCompLamtp_PointByPointSyst[iRap][iP].push_back(Dummy_Lam);//Not implemented so far -> dummy value

									//cout<<"here?"<<endl;
									vector<double> setModel_Buff_ShortDistanceCoef_globalSystPos;
									vector<double> setModel_Buff_ShortDistanceCoef_globalSystNeg;
									vector<double> setModel_Buff_OctetCompLamth_globalSystPos;
									vector<double> setModel_Buff_OctetCompLamth_globalSystNeg;
									vector<double> setModel_Buff_OctetCompLamph_globalSystPos;
									vector<double> setModel_Buff_OctetCompLamph_globalSystNeg;
									vector<double> setModel_Buff_OctetCompLamtp_globalSystPos;
									vector<double> setModel_Buff_OctetCompLamtp_globalSystNeg;



									// Get here the ModelSystematicsScales-TGraph-files of SDC, Lamth, Lamph, Lamtp
									// Get TGraphs pos and neg for state iState, and for the specific iColorChannel (we are inside the loop)
									// -> in the iModelSystematicScale loop: get Pos and Neg value (use TGraph_iColorChannel for iModelScale, for each iMeasuID)
									// Do this by integration over pT range of measurement and /=deltapT

									// Open TGraph files containing information about the ModelSystematicScales

									//cout<<"open files and get graphs"<<endl;
									TFile* ModelSystematicsFile;
									char ModelSystematicScaleName[200];

									TGraph* ModelSystematicScale_SDC_Pos;
									TGraph* ModelSystematicScale_SDC_Neg;
									TGraph* ModelSystematicScale_Lamth_Pos;
									TGraph* ModelSystematicScale_Lamth_Neg;
									TGraph* ModelSystematicScale_Lamph_Pos;
									TGraph* ModelSystematicScale_Lamph_Neg;
									TGraph* ModelSystematicScale_Lamtp_Pos;
									TGraph* ModelSystematicScale_Lamtp_Neg;

									if(ContributionExists){

										sprintf(inname,"%s/TGraphs_ModelSystematicScale%d_iMeasurementID%d.root",modelsystematicscalesdirname, iColorChannel, 0);
										ModelSystematicsFile = new TFile(inname, "READ");
										sprintf(ModelSystematicScaleName,"ModelSystematicScale_iState%d_Pos",iMother);
										ModelSystematicScale_SDC_Pos = (TGraph*) ModelSystematicsFile->Get(ModelSystematicScaleName);
										sprintf(ModelSystematicScaleName,"ModelSystematicScale_iState%d_Neg",iMother);
										ModelSystematicScale_SDC_Neg = (TGraph*) ModelSystematicsFile->Get(ModelSystematicScaleName);
										ModelSystematicsFile->Close();

										sprintf(inname,"%s/TGraphs_ModelSystematicScale%d_iMeasurementID%d.root",modelsystematicscalesdirname, iColorChannel, 1);
										ModelSystematicsFile = new TFile(inname, "READ");
										sprintf(ModelSystematicScaleName,"ModelSystematicScale_iState%d_Pos",iMother);
										ModelSystematicScale_Lamth_Pos = (TGraph*) ModelSystematicsFile->Get(ModelSystematicScaleName);
										sprintf(ModelSystematicScaleName,"ModelSystematicScale_iState%d_Neg",iMother);
										ModelSystematicScale_Lamth_Neg = (TGraph*) ModelSystematicsFile->Get(ModelSystematicScaleName);
										ModelSystematicsFile->Close();

										sprintf(inname,"%s/TGraphs_ModelSystematicScale%d_iMeasurementID%d.root",modelsystematicscalesdirname, iColorChannel, 2);
										ModelSystematicsFile = new TFile(inname, "READ");
										sprintf(ModelSystematicScaleName,"ModelSystematicScale_iState%d_Pos",iMother);
										ModelSystematicScale_Lamph_Pos = (TGraph*) ModelSystematicsFile->Get(ModelSystematicScaleName);
										sprintf(ModelSystematicScaleName,"ModelSystematicScale_iState%d_Neg",iMother);
										ModelSystematicScale_Lamph_Neg = (TGraph*) ModelSystematicsFile->Get(ModelSystematicScaleName);
										ModelSystematicsFile->Close();

										sprintf(inname,"%s/TGraphs_ModelSystematicScale%d_iMeasurementID%d.root",modelsystematicscalesdirname, iColorChannel, 3);
										ModelSystematicsFile = new TFile(inname, "READ");
										sprintf(ModelSystematicScaleName,"ModelSystematicScale_iState%d_Pos",iMother);
										ModelSystematicScale_Lamtp_Pos = (TGraph*) ModelSystematicsFile->Get(ModelSystematicScaleName);
										sprintf(ModelSystematicScaleName,"ModelSystematicScale_iState%d_Neg",iMother);
										ModelSystematicScale_Lamtp_Neg = (TGraph*) ModelSystematicsFile->Get(ModelSystematicScaleName);
										ModelSystematicsFile->Close();

										//cout<<"...finished"<<endl;


									}

									double setModel_SDC_globalSystPos_val=Dummy_SDC;
									double setModel_SDC_globalSystNeg_val=Dummy_SDC;
									double setModel_Lamth_globalSystPos_val=Dummy_SDC;
									double setModel_Lamth_globalSystNeg_val=Dummy_SDC;
									double setModel_Lamph_globalSystPos_val=Dummy_SDC;
									double setModel_Lamph_globalSystNeg_val=Dummy_SDC;
									double setModel_Lamtp_globalSystPos_val=Dummy_SDC;
									double setModel_Lamtp_globalSystNeg_val=Dummy_SDC;


									for (int iModelSystematicScale=0; iModelSystematicScale<NRQCDvars::nModelSystematicScales; iModelSystematicScale++){
										//cout<<"iModelSystematicScale = "<<iModelSystematicScale<<endl;

										if(ContributionExists){

											if(iModelSystematicScale==iColorChannel){

												double buffx, buffy;
												double Data_pTMin=DataObject[iRap][iP]->getpTMin();
												double Data_pTMax=DataObject[iRap][iP]->getpTMax();
												int nTGbins;

												nTGbins=0;
												setModel_SDC_globalSystPos_val=0;
												for(int iTGbin=0;iTGbin<ModelSystematicScale_SDC_Pos->GetN();iTGbin++){
													ModelSystematicScale_SDC_Pos->GetPoint(iTGbin, buffx, buffy);
													if(buffx>Data_pTMin&&buffx<Data_pTMax){
														setModel_SDC_globalSystPos_val+=buffy;
														nTGbins++;
													}

												}
												setModel_SDC_globalSystPos_val/=nTGbins;
												if(nTGbins<1) setModel_SDC_globalSystPos_val=0;
												nTGbins=0;
												setModel_SDC_globalSystNeg_val=0;
												for(int iTGbin=0;iTGbin<ModelSystematicScale_SDC_Neg->GetN();iTGbin++){
													ModelSystematicScale_SDC_Neg->GetPoint(iTGbin, buffx, buffy);
													if(buffx>Data_pTMin&&buffx<Data_pTMax){
														setModel_SDC_globalSystNeg_val+=buffy;
														nTGbins++;
													}

												}
												setModel_SDC_globalSystNeg_val/=nTGbins;
												if(nTGbins<1) setModel_SDC_globalSystNeg_val=0;

												nTGbins=0;
												setModel_Lamth_globalSystPos_val=0;
												for(int iTGbin=0;iTGbin<ModelSystematicScale_Lamth_Pos->GetN();iTGbin++){
													ModelSystematicScale_Lamth_Pos->GetPoint(iTGbin, buffx, buffy);
													if(buffx>Data_pTMin&&buffx<Data_pTMax){
														setModel_Lamth_globalSystPos_val+=buffy;
														nTGbins++;
													}

												}
												setModel_Lamth_globalSystPos_val/=nTGbins;
												if(nTGbins<1) setModel_Lamth_globalSystPos_val=0;
												nTGbins=0;
												setModel_Lamth_globalSystNeg_val=0;
												for(int iTGbin=0;iTGbin<ModelSystematicScale_Lamth_Neg->GetN();iTGbin++){
													ModelSystematicScale_Lamth_Neg->GetPoint(iTGbin, buffx, buffy);
													if(buffx>Data_pTMin&&buffx<Data_pTMax){
														setModel_Lamth_globalSystNeg_val+=buffy;
														nTGbins++;
													}

												}
												setModel_Lamth_globalSystNeg_val/=nTGbins;
												if(nTGbins<1) setModel_Lamth_globalSystNeg_val=0;

												nTGbins=0;
												setModel_Lamph_globalSystPos_val=0;
												for(int iTGbin=0;iTGbin<ModelSystematicScale_Lamph_Pos->GetN();iTGbin++){
													ModelSystematicScale_Lamph_Pos->GetPoint(iTGbin, buffx, buffy);
													if(buffx>Data_pTMin&&buffx<Data_pTMax){
														setModel_Lamph_globalSystPos_val+=buffy;
														nTGbins++;
													}

												}
												setModel_Lamph_globalSystPos_val/=nTGbins;
												if(nTGbins<1) setModel_Lamph_globalSystPos_val=0;
												nTGbins=0;
												setModel_Lamph_globalSystNeg_val=0;
												for(int iTGbin=0;iTGbin<ModelSystematicScale_Lamph_Neg->GetN();iTGbin++){
													ModelSystematicScale_Lamph_Neg->GetPoint(iTGbin, buffx, buffy);
													if(buffx>Data_pTMin&&buffx<Data_pTMax){
														setModel_Lamph_globalSystNeg_val+=buffy;
														nTGbins++;
													}

												}
												setModel_Lamph_globalSystNeg_val/=nTGbins;
												if(nTGbins<1) setModel_Lamph_globalSystNeg_val=0;

												nTGbins=0;
												setModel_Lamtp_globalSystPos_val=0;
												for(int iTGbin=0;iTGbin<ModelSystematicScale_Lamtp_Pos->GetN();iTGbin++){
													ModelSystematicScale_Lamtp_Pos->GetPoint(iTGbin, buffx, buffy);
													if(buffx>Data_pTMin&&buffx<Data_pTMax){
														setModel_Lamtp_globalSystPos_val+=buffy;
														nTGbins++;
													}

												}
												setModel_Lamtp_globalSystPos_val/=nTGbins;
												if(nTGbins<1) setModel_Lamtp_globalSystPos_val=0;
												nTGbins=0;
												setModel_Lamtp_globalSystNeg_val=0;
												for(int iTGbin=0;iTGbin<ModelSystematicScale_Lamtp_Neg->GetN();iTGbin++){
													ModelSystematicScale_Lamtp_Neg->GetPoint(iTGbin, buffx, buffy);
													if(buffx>Data_pTMin&&buffx<Data_pTMax){
														setModel_Lamtp_globalSystNeg_val+=buffy;
														nTGbins++;
													}

												}
												setModel_Lamtp_globalSystNeg_val/=nTGbins;
												if(nTGbins<1) setModel_Lamtp_globalSystNeg_val=0;


												//Here the setModel_-_globalSyst-_val's are converted from differences into lower and upper 'absolute' values around the 'best' model
												setModel_SDC_globalSystPos_val=Buffer_ShortDistanceCoef[iRap][iP]+setModel_SDC_globalSystPos_val;
												setModel_SDC_globalSystNeg_val=Buffer_ShortDistanceCoef[iRap][iP]+setModel_SDC_globalSystNeg_val;
												setModel_Lamth_globalSystPos_val=lamVecSumCascades[0]+setModel_Lamth_globalSystPos_val;
												setModel_Lamth_globalSystNeg_val=lamVecSumCascades[0]+setModel_Lamth_globalSystNeg_val;
												setModel_Lamph_globalSystPos_val=lamVecSumCascades[1]+setModel_Lamph_globalSystPos_val;
												setModel_Lamph_globalSystNeg_val=lamVecSumCascades[1]+setModel_Lamph_globalSystNeg_val;
												setModel_Lamtp_globalSystPos_val=lamVecSumCascades[2]+setModel_Lamtp_globalSystPos_val;
												setModel_Lamtp_globalSystNeg_val=lamVecSumCascades[2]+setModel_Lamtp_globalSystNeg_val;

												//cout<<"Buffer_ShortDistanceCoef "<<Buffer_ShortDistanceCoef[iRap][iP]<<endl;
												//cout<<"setModel_SDC_globalSystPos_val "<<setModel_SDC_globalSystPos_val<<endl;
												//cout<<"setModel_SDC_globalSystNeg_val "<<setModel_SDC_globalSystNeg_val<<endl;
												//cout<<"lamVecSumCascades[0] "<<lamVecSumCascades[0]<<endl;
												//cout<<"setModel_Lamth_globalSystPos_val "<<setModel_Lamth_globalSystPos_val<<endl;
												//cout<<"setModel_Lamth_globalSystNeg_val "<<setModel_Lamth_globalSystNeg_val<<endl;

											}
											else{
												setModel_SDC_globalSystPos_val=Buffer_ShortDistanceCoef[iRap][iP];
												setModel_SDC_globalSystNeg_val=Buffer_ShortDistanceCoef[iRap][iP];
												setModel_Lamth_globalSystPos_val=lamVecSumCascades[0];
												setModel_Lamth_globalSystNeg_val=lamVecSumCascades[0];
												setModel_Lamph_globalSystPos_val=lamVecSumCascades[1];
												setModel_Lamph_globalSystNeg_val=lamVecSumCascades[1];
												setModel_Lamtp_globalSystPos_val=lamVecSumCascades[2];
												setModel_Lamtp_globalSystNeg_val=lamVecSumCascades[2];
											}
										}




										setModel_Buff_ShortDistanceCoef_globalSystPos.push_back(setModel_SDC_globalSystPos_val);//Not used so far, but implemented -> dummy value
										setModel_Buff_OctetCompLamth_globalSystPos.push_back(setModel_Lamth_globalSystPos_val);//Not used so far, but implemented -> dummy value
										setModel_Buff_OctetCompLamph_globalSystPos.push_back(setModel_Lamph_globalSystPos_val);//Not used so far, but implemented -> dummy value
										setModel_Buff_OctetCompLamtp_globalSystPos.push_back(setModel_Lamtp_globalSystPos_val);//Not used so far, but implemented -> dummy value
										setModel_Buff_ShortDistanceCoef_globalSystNeg.push_back(setModel_SDC_globalSystNeg_val);//Not used so far, but implemented -> dummy value
										setModel_Buff_OctetCompLamth_globalSystNeg.push_back(setModel_Lamth_globalSystNeg_val);//Not used so far, but implemented -> dummy value
										setModel_Buff_OctetCompLamph_globalSystNeg.push_back(setModel_Lamph_globalSystNeg_val);//Not used so far, but implemented -> dummy value
										setModel_Buff_OctetCompLamtp_globalSystNeg.push_back(setModel_Lamtp_globalSystNeg_val);//Not used so far, but implemented -> dummy value

									}
									if(NRQCDvars::nModelSystematicScales==0){//fill with dummy values, because >> << operators can not cope with non-existing dcubes
										for (int iModelSystematicScale=0; iModelSystematicScale<2; iModelSystematicScale++){
											//cout<<"iModelSystematicScale?"<<iModelSystematicScale<<endl;
											setModel_Buff_ShortDistanceCoef_globalSystPos.push_back(Dummy_SDC);//Not used so far, but implemented -> dummy value
											setModel_Buff_OctetCompLamth_globalSystPos.push_back(Dummy_Lam);//Not used so far, but implemented -> dummy value
											setModel_Buff_OctetCompLamph_globalSystPos.push_back(Dummy_Lam);//Not used so far, but implemented -> dummy value
											setModel_Buff_OctetCompLamtp_globalSystPos.push_back(Dummy_Lam);//Not used so far, but implemented -> dummy value
											setModel_Buff_ShortDistanceCoef_globalSystNeg.push_back(Dummy_SDC);//Not used so far, but implemented -> dummy value
											setModel_Buff_OctetCompLamth_globalSystNeg.push_back(Dummy_Lam);//Not used so far, but implemented -> dummy value
											setModel_Buff_OctetCompLamph_globalSystNeg.push_back(Dummy_Lam);//Not used so far, but implemented -> dummy value
											setModel_Buff_OctetCompLamtp_globalSystNeg.push_back(Dummy_Lam);//Not used so far, but implemented -> dummy value
										}
									}

									setModel_ShortDistanceCoef_globalSystPos[iRap][iP].push_back(setModel_Buff_ShortDistanceCoef_globalSystPos);
									setModel_OctetCompLamth_globalSystPos[iRap][iP].push_back(setModel_Buff_OctetCompLamth_globalSystPos);
									setModel_OctetCompLamph_globalSystPos[iRap][iP].push_back(setModel_Buff_OctetCompLamph_globalSystPos);
									setModel_OctetCompLamtp_globalSystPos[iRap][iP].push_back(setModel_Buff_OctetCompLamtp_globalSystPos);
									setModel_ShortDistanceCoef_globalSystNeg[iRap][iP].push_back(setModel_Buff_ShortDistanceCoef_globalSystNeg);
									setModel_OctetCompLamth_globalSystNeg[iRap][iP].push_back(setModel_Buff_OctetCompLamth_globalSystNeg);
									setModel_OctetCompLamph_globalSystNeg[iRap][iP].push_back(setModel_Buff_OctetCompLamph_globalSystNeg);
									setModel_OctetCompLamtp_globalSystNeg[iRap][iP].push_back(setModel_Buff_OctetCompLamtp_globalSystNeg);
									//cout<<"here?"<<endl;

									    }
									}

								}

								for(int iRap = 0; iRap < NRQCDvars::nMaxRapBins; iRap++){
								    for(int iP = 0; iP < NRQCDvars::nMaxPtBins; iP++){
								    	if(!doesDataExist[iRap][iP]) continue;

								    	//cout<<"ModelObject->setModel(iMother, ..."<<endl;
								ModelObject[iRap][iP]->setModel(iMother,
										  setModel_ShortDistanceCoef[iRap][iP], setModel_OctetCompLamth[iRap][iP], setModel_OctetCompLamph[iRap][iP], setModel_OctetCompLamtp[iRap][iP],
										  setModel_ShortDistanceCoef_PointByPointSyst[iRap][iP], setModel_OctetCompLamth_PointByPointSyst[iRap][iP], setModel_OctetCompLamph_PointByPointSyst[iRap][iP], setModel_OctetCompLamtp_PointByPointSyst[iRap][iP],
										  setModel_ShortDistanceCoef_globalSystPos[iRap][iP], setModel_ShortDistanceCoef_globalSystNeg[iRap][iP],
										  setModel_OctetCompLamth_globalSystPos[iRap][iP], setModel_OctetCompLamth_globalSystNeg[iRap][iP],
										  setModel_OctetCompLamph_globalSystPos[iRap][iP], setModel_OctetCompLamph_globalSystNeg[iRap][iP],
										  setModel_OctetCompLamtp_globalSystPos[iRap][iP], setModel_OctetCompLamtp_globalSystNeg[iRap][iP],
										  true);
								    }
								}

							}

							for(int iRap = 0; iRap < NRQCDvars::nMaxRapBins; iRap++){
							    for(int iP = 0; iP < NRQCDvars::nMaxPtBins; iP++){
							    	if(!doesDataExist[iRap][iP]) continue;

							if(NRQCDvars::debug) ModelObject[iRap][iP]->Dump(NRQCDvars::nStates,true,true);

							//write updated NRQCDglobalfitObject containing data+model components to output file
							sprintf(outname,"%s/ConvertedDataModel_%s_%s_%s_rap%d_pT%d.txt",datamodeldirname, StateName[iState],  MeasurementIDName[iMeasurementID],  ExpName[iExperiment], iRap, iP);
							ofstream out;
							out.open(outname);
							out << *ModelObject[iRap][iP];
							out.close();

							    }
							}
//						}//end if(in!=NULL)

//				    }//iP
//				}//iRap
			}//iExperiment
		}//iMeasurementID
	}//iState


	ModelIngredientsFile->Close();

return 0;

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

vector<double> GeneralPolarizationTransformation(vector<double> lamVecOriginal, double t1, double t2, double t3){


	vector<double> lamVecTransformed (3,0);

	lamVecTransformed[0]=t1*lamVecOriginal[0]/(t2+t3*lamVecOriginal[0]);
	lamVecTransformed[1]=t1*lamVecOriginal[1]/(t2+t3*lamVecOriginal[0]);
	lamVecTransformed[2]=t1*lamVecOriginal[2]/(t2+t3*lamVecOriginal[0]);

	return lamVecTransformed;

}

vector<double> PolarizationTransfer(vector<double> lamVecOriginal, vector<int> DecayChain){//DecayChain contains the state ID for the decay chain, e.g. nMother, nIntermediateState, nDaughter

	vector<double> lamVecTransfer (3,0);
	bool debugPolTransfer=false;
	if(debugPolTransfer){
		cout<<"lamVecOriginal[0] "<<lamVecOriginal[0]<<endl;
		cout<<"lamVecOriginal[1] "<<lamVecOriginal[1]<<endl;
		cout<<"lamVecOriginal[2] "<<lamVecOriginal[2]<<endl;
	}
	int nCascadeDecays=DecayChain.size();

	if(debugPolTransfer){
		for (int iCascadeDecay=0; iCascadeDecay<nCascadeDecays; iCascadeDecay++){//Loop over all different decay chains linking iMother->iDaughter
			cout<<"DecayChain["<<iCascadeDecay<<"]: "<<DecayChain[iCascadeDecay]<<endl;
		}
		for (int iCascadeDecay=0; iCascadeDecay<nCascadeDecays; iCascadeDecay++){//Loop over all different decay chains linking iMother->iDaughter
			cout<<"NRQCDvars::StateQuantumID[DecayChain["<<iCascadeDecay<<"]]: "<<NRQCDvars::StateQuantumID[DecayChain[iCascadeDecay]]<<endl;
		}
	}
	//If nCascadeDecays==1: No Decay -> direct component -> no polarization transfer
	if(nCascadeDecays==1){
		if(debugPolTransfer) cout<<"return lamVecOriginal; (nCascadeDecays==1)"<<endl;
		return lamVecOriginal;
	}

	//Remove S -> S links in the decay chain
	vector<int> AlteredDecayChain;
	bool ignoreThisDecay;
	for (int iCascadeDecay=0; iCascadeDecay<nCascadeDecays; iCascadeDecay++){//Loop over all different decay chains linking iMother->iDaughter
		ignoreThisDecay=false;
		if(iCascadeDecay<nCascadeDecays-1){
			if(NRQCDvars::StateQuantumID[DecayChain[iCascadeDecay]]==NRQCDvars::quID_S && NRQCDvars::StateQuantumID[DecayChain[iCascadeDecay+1]]==NRQCDvars::quID_S) ignoreThisDecay=true;
		}
		if(!ignoreThisDecay) AlteredDecayChain.push_back(DecayChain[iCascadeDecay]);
	}

	int nAlteredCascadeDecays=AlteredDecayChain.size();
	//If nAlteredCascadeDecays==1: polarization transfer identical to input (only S states in chain)
	if(nAlteredCascadeDecays==1){
		if(debugPolTransfer) cout<<"return lamVecOriginal; (nAlteredCascadeDecays==1)"<<endl;
		return lamVecOriginal;
	}


	if(debugPolTransfer) {
		for (int iCascadeDecay=0; iCascadeDecay<nAlteredCascadeDecays; iCascadeDecay++){//Loop over all different decay chains linking iMother->iDaughter
			cout<<"NRQCDvars::StateQuantumID[AlteredDecayChain["<<iCascadeDecay<<"]]: "<<NRQCDvars::StateQuantumID[AlteredDecayChain[iCascadeDecay]]<<endl;
		}
	}

	vector<double> lamVecStart (3,0);

	//If first decay in altered decay chain is a S->P decay -> apply the transformation from the dilepton decay (given by Sergey) to the radiative decay
	bool applyTransformationFromDileptonToRadiativeDecay=false;
	if(NRQCDvars::StateQuantumID[AlteredDecayChain[0]]==NRQCDvars::quID_S && NRQCDvars::StateQuantumID[AlteredDecayChain[1]]==NRQCDvars::quID_P1){
		//chi_1:  vec{lambda}(psi->chi_1 gamma) = - vec{lambda}(psi->l+l-) / [2 + lambda_theta(psi->l+l-)]
		lamVecStart = GeneralPolarizationTransformation(lamVecOriginal, -1., 2., 1.);
		if(debugPolTransfer) cout<<"Transform S -> ll to S -> P1 lambdas"<<endl;
	}
	else if(NRQCDvars::StateQuantumID[AlteredDecayChain[0]]==NRQCDvars::quID_S && NRQCDvars::StateQuantumID[AlteredDecayChain[1]]==NRQCDvars::quID_P2){
		//chi_2:  vec{lambda}(psi->chi_2 gamma) = vec{lambda}(psi->l+l-) / [10 + 3 lambda_theta(psi->l+l-)]
		lamVecStart = GeneralPolarizationTransformation(lamVecOriginal, 1., 10., 3.);
		if(debugPolTransfer) cout<<"Transform S -> ll to S -> P2 lambdas"<<endl;
	}
	else if(NRQCDvars::StateQuantumID[AlteredDecayChain[0]]==NRQCDvars::quID_P1 || NRQCDvars::StateQuantumID[AlteredDecayChain[0]]==NRQCDvars::quID_P2){
		// If first state in chain is a P state -> transformation not necessary
		lamVecStart = GeneralPolarizationTransformation(lamVecOriginal, 1., 1., 0.);
		if(debugPolTransfer) cout<<"Use PJ -> S gamma (input) lambdas"<<endl;
	}
	else throw;

	if(debugPolTransfer){
		cout<<"lamVecStart[0] "<<lamVecStart[0]<<endl;
		cout<<"lamVecStart[1] "<<lamVecStart[1]<<endl;
		cout<<"lamVecStart[2] "<<lamVecStart[2]<<endl;
	}

	if(nAlteredCascadeDecays==2){
		//In this case, we either have a
		//S -> P or P -> S decay
		//S -> P: return lamVecStart
		//P -> S: return lamVecStart
		if(debugPolTransfer) cout<<"return lamVecOriginal; (nAlteredCascadeDecays==2)"<<endl;
		return lamVecOriginal;
	}

	//Now, for the more complicated decay chains (at least two P-S cross decays):
	//Loop over altered decay chain and apply all transformations
	vector<vector<double> > lamVecTransformationMatrix;
	lamVecTransformationMatrix.push_back(lamVecStart);
	double t1, t2, t3;
	int TransformationID;
	for (int iCascadeDecay=1; iCascadeDecay<nAlteredCascadeDecays-1; iCascadeDecay++){//Loop over all different decay chains linking iMother->iDaughter

		if(NRQCDvars::StateQuantumID[AlteredDecayChain[iCascadeDecay-1]]==NRQCDvars::quID_S
		   && NRQCDvars::StateQuantumID[AlteredDecayChain[iCascadeDecay]]==NRQCDvars::quID_P1
		   && NRQCDvars::StateQuantumID[AlteredDecayChain[iCascadeDecay+1]]==NRQCDvars::quID_S) {t1=-1.; t2=2.; t3=1.; 	if(debugPolTransfer) cout<<"Transformation 3a)"<<endl;}
		if(NRQCDvars::StateQuantumID[AlteredDecayChain[iCascadeDecay-1]]==NRQCDvars::quID_S
		   && NRQCDvars::StateQuantumID[AlteredDecayChain[iCascadeDecay]]==NRQCDvars::quID_P2
		   && NRQCDvars::StateQuantumID[AlteredDecayChain[iCascadeDecay+1]]==NRQCDvars::quID_S) {t1=21.; t2=6.; t3=-5.; if(debugPolTransfer) cout<<"Transformation 3b)"<<endl;}
		if(NRQCDvars::StateQuantumID[AlteredDecayChain[iCascadeDecay-1]]==NRQCDvars::quID_P1
		   && NRQCDvars::StateQuantumID[AlteredDecayChain[iCascadeDecay]]==NRQCDvars::quID_S
		   && NRQCDvars::StateQuantumID[AlteredDecayChain[iCascadeDecay+1]]==NRQCDvars::quID_P1) {t1=-1.; t2=2.; t3=1.; if(debugPolTransfer) cout<<"Transformation 4a)"<<endl;}
		if(NRQCDvars::StateQuantumID[AlteredDecayChain[iCascadeDecay-1]]==NRQCDvars::quID_P1
		   && NRQCDvars::StateQuantumID[AlteredDecayChain[iCascadeDecay]]==NRQCDvars::quID_S
		   && NRQCDvars::StateQuantumID[AlteredDecayChain[iCascadeDecay+1]]==NRQCDvars::quID_P2) {t1=1.; t2=10.; t3=3.; if(debugPolTransfer) cout<<"Transformation 4b)"<<endl;}
		if(NRQCDvars::StateQuantumID[AlteredDecayChain[iCascadeDecay-1]]==NRQCDvars::quID_P2
		   && NRQCDvars::StateQuantumID[AlteredDecayChain[iCascadeDecay]]==NRQCDvars::quID_S
		   && NRQCDvars::StateQuantumID[AlteredDecayChain[iCascadeDecay+1]]==NRQCDvars::quID_P1) {t1=-1.; t2=2.; t3=1.; if(debugPolTransfer) cout<<"Transformation 4c)"<<endl;}
		if(NRQCDvars::StateQuantumID[AlteredDecayChain[iCascadeDecay-1]]==NRQCDvars::quID_P2
		   && NRQCDvars::StateQuantumID[AlteredDecayChain[iCascadeDecay]]==NRQCDvars::quID_S
		   && NRQCDvars::StateQuantumID[AlteredDecayChain[iCascadeDecay+1]]==NRQCDvars::quID_P2) {t1=1.; t2=10.; t3=3.; if(debugPolTransfer) cout<<"Transformation 4d)"<<endl;}

		vector<double> lamVecTransformed=GeneralPolarizationTransformation(lamVecTransformationMatrix[iCascadeDecay-1], t1, t2, t3);

		if(debugPolTransfer){
			cout<<"lamVecTransformed[0] "<<lamVecTransformed[0]<<endl;
			cout<<"lamVecTransformed[1] "<<lamVecTransformed[1]<<endl;
			cout<<"lamVecTransformed[2] "<<lamVecTransformed[2]<<endl;
		}

		lamVecTransformationMatrix.push_back(lamVecTransformed);
	}

	lamVecTransfer=lamVecTransformationMatrix[nAlteredCascadeDecays-2];

	if(debugPolTransfer){
		cout<<"lamVecTransfer[0] "<<lamVecTransfer[0]<<endl;
		cout<<"lamVecTransfer[1] "<<lamVecTransfer[1]<<endl;
		cout<<"lamVecTransfer[2] "<<lamVecTransfer[2]<<endl;
	}


	return lamVecTransfer;
}

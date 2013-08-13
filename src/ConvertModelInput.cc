/*
 * ConvertModelInput.cc
 *
 *  Created on: Jun 18, 2013
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
#include "TGraph.h"
#include "TGraph2D.h"
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
vector<double> Transform_kinematics_MotherToDaughter(int MotherState, int DaughterState, vector<double> Kinematics_Mother);


int main(int argc, char** argv) {

  	Char_t *ModelID = "Default";
	Char_t *storagedir = "Default"; //Storage Directory
  	bool useToyModel=false;
  	bool calcOnlyConstsStat=false;

  	for( int i=0;i < argc; ++i ) {
		if(std::string(argv[i]).find("ModelID") != std::string::npos) {char* ModelIDchar = argv[i]; char* ModelIDchar2 = strtok (ModelIDchar, "="); ModelID = ModelIDchar2; cout<<"ModelID = "<<ModelID<<endl;}
		if(std::string(argv[i]).find("useToyModel=true") != std::string::npos) {  useToyModel=true, cout<<"useToyModel"<<endl;	}
		if(std::string(argv[i]).find("calcOnlyConstsStat=true") != std::string::npos) {  calcOnlyConstsStat=true, cout<<"calcOnlyConstsStat"<<endl;	}
		if(std::string(argv[i]).find("storagedir") != std::string::npos) {char* storagedirchar = argv[i]; char* storagedirchar2 = strtok (storagedirchar, "="); storagedir = storagedirchar2; cout<<"storagedir = "<<storagedir<<endl;}
		}




	delete gRandom;
	gRandom = new TRandom3(NRQCDvars::randomSeed);  // better random generator

	double model_pTMin, model_pTMax;
	double model_rapMin, model_rapMax;
	double model_costhMin, model_costhMax;
	double model_phiMin, model_phiMax;

	model_pTMin=8;
	model_pTMax=100;
	model_rapMin=-5;
	model_rapMax=5;
	model_costhMin=-1;
	model_costhMax=1;
	model_phiMin=-180;
	model_phiMax=180;

	int nBins_pT=100;
	int nBins_rap=100;


	// TODO: add loop for nSystematics (now nSystematics will be 0)

	char outname[2000];
	char inname[2000];
	char predirname[2000];
	char dirname[2000];
	sprintf(predirname,"%s/ModelID", storagedir);
	gSystem->mkdir(predirname);
	sprintf(dirname,"%s/%s",predirname,ModelID);
	gSystem->mkdir(dirname);


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Read input model
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	sprintf(inname,"%s/OriginalModelTTree.root",dirname);
	TFile* OriginalModelTTreeFile = new TFile(inname, "READ");


		TTree* nTupleModel[NRQCDvars::nStates][NRQCDvars::nColorChannels];
		char nTupleModelName[1000];

		bool nTupleModelGiven[NRQCDvars::nStates];

		if(!useToyModel){


			//Read in TTrees from ModelID/JobID folder
			for (int iMother=0; iMother<NRQCDvars::nStates; iMother++){
				nTupleModelGiven[iMother]=true;
				int nColorChannels_state;
				bool isSstate=(StateQuantumID[iMother] > NRQCDvars::quID_S)?false:true;
				if(isSstate) nColorChannels_state=NRQCDvars::nColorChannels_S;
				else nColorChannels_state=NRQCDvars::nColorChannels_P;
				/*if(NRQCDvars::debug) */cout<<"Read in original TTree model for nState = "<<iMother<<endl;
				for (int iColorChannel=0; iColorChannel<nColorChannels_state; iColorChannel++){
					if(NRQCDvars::debug) cout<<"		Color channel = "<<iColorChannel<<endl;
					sprintf(nTupleModelName,"nTupleModel_Mother%d_Daughter%d_ColorChannel%d_Cascade0",iMother, iMother, iColorChannel);
					if(OriginalModelTTreeFile->Get(nTupleModelName)!=NULL){
						nTupleModel[iMother][iColorChannel]=(TTree*)OriginalModelTTreeFile->Get(nTupleModelName);
					}
					else{
						nTupleModelGiven[iMother]=false;
					}
					//cout<<"nTupleModel[iMother][iColorChannel]->GetEntries() "<<nTupleModel[iMother][iColorChannel]->GetEntries();
				}
			}


		}




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Alternatively, let's generate some toy data
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


		if(useToyModel){
			//double LuminosityPerModelEvent[NRQCDvars::nStates][NRQCDvars::nColorChannels];
			int n_nTuple=1000000;
			for (int iMother=0; iMother<NRQCDvars::nStates; iMother++){
				nTupleModelGiven[iMother]=true;
				int nColorChannels_state;
				bool isSstate=(StateQuantumID[iMother] > NRQCDvars::quID_S)?false:true;
				if(isSstate) nColorChannels_state=NRQCDvars::nColorChannels_S;
				else nColorChannels_state=NRQCDvars::nColorChannels_P;
				cout<<"Generating toy model for nState = "<<iMother<<endl;
				for (int iColorChannel=0; iColorChannel<nColorChannels_state; iColorChannel++){
					cout<<"		Color channel = "<<iColorChannel<<endl;

					//LuminosityPerModelEvent[iMother][iColorChannel]=100.; // Sets weight of event, to take care of the normalization of the individual states and color channels

					sprintf(nTupleModelName,"nTupleModel_Mother%d_Daughter%d_ColorChannel%d_Cascade0",iMother, iMother, iColorChannel);
					nTupleModel[iMother][iColorChannel] = new TTree(nTupleModelName, nTupleModelName);

					double model_pT;  nTupleModel[iMother][iColorChannel]->Branch("model_pT",     &model_pT,     "model_pT/D");
					double model_rap;    nTupleModel[iMother][iColorChannel]->Branch("model_rap",       &model_rap,       "model_rap/D"  );
					double model_costh;  nTupleModel[iMother][iColorChannel]->Branch("model_costh",      &model_costh,     "model_costh/D");
					double model_phi;  nTupleModel[iMother][iColorChannel]->Branch("model_phi",      &model_phi,     "model_phi/D");
					double model_weight; nTupleModel[iMother][iColorChannel]->Branch("weight",      &model_weight,     "weight/D");

					int nStep = n_nTuple/10;  // visualize progress of the parameter sampling
					int nStep_ = 0;
					bool genAngDist2D=false;

					for (int k=0; k<n_nTuple; k++){

				 		if (k%nStep == 0) {
				 			cout << nStep_*10 <<"% (nState = "<<iMother<<", Color channel = "<<iColorChannel<<")"<<endl;;
				 			++nStep_;
				 		}

						TF1* pT_distr;

						model_pT=gRandom->Uniform(model_pTMin,model_pTMax);
						model_rap=gRandom->Uniform(model_rapMin,model_rapMax);
						model_costh=gRandom->Uniform(model_costhMin,model_costhMax);
						model_phi=gRandom->Uniform(model_phiMin,model_phiMax);

						//Define pT distribution:
						pT_distr = new TF1("pT_distr",func_pT_gen,model_pTMin,model_pTMax,2);
						pT_distr->SetParameter(0,iMother);
						pT_distr->SetParameter(1,iColorChannel);
						model_weight=pT_distr->Eval(model_pT);

						//Define angular distribution:


						dvector model_lam = func_lam_gen(iMother, iColorChannel);
						double model_lamth, model_lamph, model_lamtp;
						model_lamth=model_lam[0];
						model_lamph=model_lam[1];
						model_lamtp=model_lam[2];

						double polNormFactor;
						if(genAngDist2D){
							TF2 *fcosthphi;
							fcosthphi = new TF2( "fcosthphi", "[0]*(1.+[1]*x[0]*x[0]+[2]*(1.-x[0]*x[0])*cos(2.*x[1]*0.0174532925)+[3]*2.*x[0]*sqrt(1.-x[0]*x[0])*cos(x[1]*0.0174532925))", -1., 1., -180., 180. );
							fcosthphi->SetParameters(1.,model_lamth, model_lamph, model_lamtp);
							model_weight*=fcosthphi->Eval(model_costh,model_phi);
							polNormFactor=fcosthphi->Integral(-1., 1., -180., 180. );
							model_weight*=2.*360./fcosthphi->Integral(-1., 1., -180., 180. );

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

						nTupleModel[iMother][iColorChannel]->Fill();

						delete pT_distr;

					}


					if(NRQCDvars::debug) cout<<"Number of generated events: "<<nTupleModel[iMother][iColorChannel]->GetEntries()<<endl;
					double globalweight=1.;
					nTupleModel[iMother][iColorChannel]->SetWeight(globalweight);

				}
			}
		}




		sprintf(outname,"%s/ModelIngredients.root",dirname);
		TFile* ModelIngredientsFile;

		if(!useToyModel){
			gSystem->Unlink(outname);
			gSystem->CopyFile(inname,outname,kTRUE);
			ModelIngredientsFile = new TFile(outname, "UPDATE");
		}
		if(useToyModel){
			ModelIngredientsFile = new TFile(outname, "RECREATE");
			for (int iMother=0; iMother<NRQCDvars::nStates; iMother++){
				int nColorChannels_state;
				bool isSstate=(StateQuantumID[iMother] > NRQCDvars::quID_S)?false:true;
				if(isSstate) nColorChannels_state=NRQCDvars::nColorChannels_S;
				else nColorChannels_state=NRQCDvars::nColorChannels_P;
				/*if(NRQCDvars::debug) */cout<<"Save original TTree model for nState = "<<iMother<<endl;
				for (int iColorChannel=0; iColorChannel<nColorChannels_state; iColorChannel++){
					nTupleModel[iMother][iColorChannel]->Write();
				}
			}
		}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Transformation of the kinematics for each mother->daughter channel and color channel -> Save these TTree's[nMothers][nDaughters][nColorChannels], (+nSystematics) in a root file
// Loop through [nMothers][nDaughters][nColorChannels]
// Transfrom momenta to momenta of the daughter
// Transform lambdas to the lambdas of the daughter
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		if(NRQCDvars::debug) cout<<"Transform Mother->Daughter kinematics"<<endl;

		TTree* nTupleModel_Transformed[NRQCDvars::nStates][NRQCDvars::nStates][NRQCDvars::nColorChannels][NRQCDvars::nMaxCascades];//[iMother][iDaughter][iColorChannel]
		TGraph* DecayChainID[NRQCDvars::nStates][NRQCDvars::nStates][NRQCDvars::nColorChannels][NRQCDvars::nMaxCascades];//[iMother][iDaughter][iColorChannel]
		char DecayChainIDName[1000];

		for (int iMother=0; iMother<NRQCDvars::nStates; iMother++){
			int nColorChannels_state;
			bool isSstate=(StateQuantumID[iMother] > NRQCDvars::quID_S)?false:true;
			if(isSstate) nColorChannels_state=NRQCDvars::nColorChannels_S;
			else nColorChannels_state=NRQCDvars::nColorChannels_P;
			for (int iDaughter=0; iDaughter<iMother; iDaughter++){

				if(!nTupleModelGiven[iMother] || !nTupleModelGiven[iDaughter]) continue;

					cout<<"transformation Mother "<<StateName[iMother]<<" to Daughter "<<StateName[iDaughter]<<endl;

				vector<vector<int> > Buffer_DecayChains;//Includes the nState integers of each decay chain from iMother -> iDaughter; includes iDaughter and iMother
				vector<double> Buffer_DecayChainBranchingRatios;
				vector<double> Buffer_DecayChainBranchingRatioFractions;

				// Fill ListFeedDownStates, ListFeedDownStateBranchingFractions
				for(int i=iDaughter; i<iMother+1; i++){
					vector<int> Buffer_DecayChain;
					// Loop over the states decaying into iDaughter
					if(NRQCDvars::FeedDownBranchingRatio[iDaughter][i]>0){
						if(i==iMother) {vector<int> Buffer_DecayChain; Buffer_DecayChain.push_back(iMother);  Buffer_DecayChain.push_back(iDaughter);  Buffer_DecayChains.push_back(Buffer_DecayChain); Buffer_DecayChainBranchingRatios.push_back(0.01*NRQCDvars::FeedDownBranchingRatio[iDaughter][i]);}
						// Add Cascade decays to the prompt components (maximum considered: cascade with 5 decays)
						for(int j=i+1; j<iMother+1; j++){
							if(NRQCDvars::FeedDownBranchingRatio[i][j]>0){
								if(j==iMother) {vector<int> Buffer_DecayChain; Buffer_DecayChain.push_back(iMother); Buffer_DecayChain.push_back(i); Buffer_DecayChain.push_back(iDaughter);  Buffer_DecayChains.push_back(Buffer_DecayChain); Buffer_DecayChainBranchingRatios.push_back(0.01*NRQCDvars::FeedDownBranchingRatio[iDaughter][i]*0.01*NRQCDvars::FeedDownBranchingRatio[i][j]);}
								for(int k=j+1; k<iMother+1; k++){
									if(NRQCDvars::FeedDownBranchingRatio[j][k]>0){
										if(k==iMother) {vector<int> Buffer_DecayChain; Buffer_DecayChain.push_back(iMother); Buffer_DecayChain.push_back(j); Buffer_DecayChain.push_back(i); Buffer_DecayChain.push_back(iDaughter);  Buffer_DecayChains.push_back(Buffer_DecayChain); Buffer_DecayChainBranchingRatios.push_back(0.01*NRQCDvars::FeedDownBranchingRatio[iDaughter][i]*0.01*NRQCDvars::FeedDownBranchingRatio[i][j]*0.01*NRQCDvars::FeedDownBranchingRatio[j][k]);}
										for(int l=k+1; l<iMother+1; l++){
											if(NRQCDvars::FeedDownBranchingRatio[k][l]>0){
												if(l==iMother) {vector<int> Buffer_DecayChain; Buffer_DecayChain.push_back(iMother); Buffer_DecayChain.push_back(k); Buffer_DecayChain.push_back(j); Buffer_DecayChain.push_back(i); Buffer_DecayChain.push_back(iDaughter);  Buffer_DecayChains.push_back(Buffer_DecayChain); Buffer_DecayChainBranchingRatios.push_back(0.01*NRQCDvars::FeedDownBranchingRatio[iDaughter][i]*0.01*NRQCDvars::FeedDownBranchingRatio[i][j]*0.01*NRQCDvars::FeedDownBranchingRatio[j][k]*0.01*NRQCDvars::FeedDownBranchingRatio[k][l]);}
												for(int m=l+1; m<iMother+1; m++){
													if(NRQCDvars::FeedDownBranchingRatio[l][m]>0){
														if(m==iMother) {vector<int> Buffer_DecayChain; Buffer_DecayChain.push_back(iMother); Buffer_DecayChain.push_back(l); Buffer_DecayChain.push_back(k); Buffer_DecayChain.push_back(j); Buffer_DecayChain.push_back(i); Buffer_DecayChain.push_back(iDaughter); Buffer_DecayChains.push_back(Buffer_DecayChain); Buffer_DecayChainBranchingRatios.push_back(0.01*NRQCDvars::FeedDownBranchingRatio[iDaughter][i]*0.01*NRQCDvars::FeedDownBranchingRatio[i][j]*0.01*NRQCDvars::FeedDownBranchingRatio[j][k]*0.01*NRQCDvars::FeedDownBranchingRatio[k][l]*0.01*NRQCDvars::FeedDownBranchingRatio[l][m]);}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}

				if(Buffer_DecayChains.size()<1) {cout<<"No contributing decay chain"<<endl; continue;} // continue, if there is not a single deca chain linkin iMother -> iDaughter

				double FullBR=0;
				for(int iDecayChain=0;iDecayChain<Buffer_DecayChainBranchingRatios.size();iDecayChain++){
					FullBR+=Buffer_DecayChainBranchingRatios[iDecayChain];
				}
				for(int iDecayChain=0;iDecayChain<Buffer_DecayChainBranchingRatios.size();iDecayChain++){
					Buffer_DecayChainBranchingRatioFractions.push_back(Buffer_DecayChainBranchingRatios[iDecayChain]/FullBR);
				}


				for (int iColorChannel=0; iColorChannel<nColorChannels_state; iColorChannel++){

					// find all decay chains nDecayChains from iMother -> iDaughter
					// Transfrom mother momenta to momenta of the daughters, for each decay chain
					// weight the events
					// Save the transformed nTuples and the projections (TH2's) in an external file


					for (int iCascade=0; iCascade<Buffer_DecayChains.size(); iCascade++){//Loop over all different decay chains linking iMother->iDaughter

						int nCascadeDecays = Buffer_DecayChains[iCascade].size()-1;

						// SetBranches of original nTuple
						double model_pT_mother; nTupleModel[iMother][iColorChannel]->SetBranchAddress( "model_pT",  &model_pT_mother );
						double model_rap_mother; nTupleModel[iMother][iColorChannel]->SetBranchAddress( "model_rap",  &model_rap_mother );
						double model_costh_mother; nTupleModel[iMother][iColorChannel]->SetBranchAddress( "model_costh",  &model_costh_mother );
						double model_phi_mother; nTupleModel[iMother][iColorChannel]->SetBranchAddress( "model_phi",  &model_phi_mother );
						double weight_mother; nTupleModel[iMother][iColorChannel]->SetBranchAddress( "weight",  &weight_mother );

						sprintf(nTupleModelName,"nTupleModel_Mother%d_Daughter%d_ColorChannel%d_Cascade%d",iMother, iDaughter, iColorChannel, iCascade);
						nTupleModel_Transformed[iMother][iDaughter][iColorChannel][iCascade] = new TTree(nTupleModelName, nTupleModelName);

						// SetBranches of transformed nTuple
						double model_pT_daughter; nTupleModel_Transformed[iMother][iDaughter][iColorChannel][iCascade]->Branch("model_pT",     &model_pT_daughter,     "model_pT/D");
						double model_rap_daughter; nTupleModel_Transformed[iMother][iDaughter][iColorChannel][iCascade]->Branch("model_rap",       &model_rap_daughter,       "model_rap/D"  );
						double model_costh_daughter; nTupleModel_Transformed[iMother][iDaughter][iColorChannel][iCascade]->Branch("model_costh",      &model_costh_daughter,     "model_costh/D");
						double model_phi_daughter; nTupleModel_Transformed[iMother][iDaughter][iColorChannel][iCascade]->Branch("model_phi",      &model_phi_daughter,     "model_phi/D");
						double weight_daughter; nTupleModel_Transformed[iMother][iDaughter][iColorChannel][iCascade]->Branch("weight",      &weight_daughter,     "weight/D");

						sprintf(DecayChainIDName,"DecayChainID_Mother%d_Daughter%d_ColorChannel%d_Cascade%d",iMother, iDaughter, iColorChannel, iCascade);

						int DecayChainIterator[nCascadeDecays+1], DecayChainState[nCascadeDecays+1];
						for (int iCascadeDecay=0; iCascadeDecay<nCascadeDecays+1; iCascadeDecay++){//Loop over all different decay chains linking iMother->iDaughter
							DecayChainIterator[iCascadeDecay]=iCascadeDecay;
							DecayChainState[iCascadeDecay]=Buffer_DecayChains[iCascade][iCascadeDecay];
						}
						DecayChainID[iMother][iDaughter][iColorChannel][iCascade] = new TGraph(nCascadeDecays+1, DecayChainIterator, DecayChainState);
						DecayChainID[iMother][iDaughter][iColorChannel][iCascade] -> SetName(DecayChainIDName);

						int n_events = int( nTupleModel[iMother][iColorChannel]->GetEntries() );

						// loop over  events in the model ntuple
						for ( int i_event = 1; i_event <= n_events; i_event++ ) {

							  nTupleModel[iMother][iColorChannel]->GetEvent( i_event-1 );

							  //In here, transform the kinematics of mother -> X -> daughter
							  vector<double> Kinematics_Mother;
							  vector<vector<double> > Kinematics_TransformationMatrix;

							  Kinematics_Mother.push_back(model_pT_mother);
							  Kinematics_Mother.push_back(model_rap_mother);
							  Kinematics_Mother.push_back(model_costh_mother);
							  Kinematics_Mother.push_back(model_phi_mother);

							  Kinematics_TransformationMatrix.push_back(Kinematics_Mother);

							  for (int iCascadeDecay=0; iCascadeDecay<nCascadeDecays; iCascadeDecay++){//Loop over all different decay chains linking iMother->iDaughter
								  //if(i_event==1&&iColorChannel==0)cout<<"transform decay "<<StateName[Buffer_DecayChains[iCascade][iCascadeDecay]]<<" -> "<<StateName[Buffer_DecayChains[iCascade][iCascadeDecay+1]]<<endl;
								  vector<double> Kinematics_Daughter = Transform_kinematics_MotherToDaughter(Buffer_DecayChains[iCascade][iCascadeDecay], Buffer_DecayChains[iCascade][iCascadeDecay+1], Kinematics_TransformationMatrix[iCascadeDecay]);
								  Kinematics_TransformationMatrix.push_back(Kinematics_Daughter);
							  }


							  model_pT_daughter=Kinematics_TransformationMatrix[nCascadeDecays][0];
							  model_rap_daughter=Kinematics_TransformationMatrix[nCascadeDecays][1];
							  //Assigning the costh/phi values of the mother to the daughter. The transformed polarization due to the decay mechanism is calculated in CombineDataModel, for a given kinematic interval (with Pietro's polarization transfer formulas)
							  model_costh_daughter=model_costh_mother;
							  model_phi_daughter=model_phi_mother;
							  weight_daughter=weight_mother;

							  //if(i_event==1&&iColorChannel==0)cout<<"weight_daughter "<<weight_daughter<<endl;

							  nTupleModel_Transformed[iMother][iDaughter][iColorChannel][iCascade]->Fill();

						}

						double globalweight_daughter=nTupleModel[iMother][iColorChannel]->GetWeight()*Buffer_DecayChainBranchingRatioFractions[iCascade]; // Luminosity/event times product of Branching ratios of the full decay chain
						nTupleModel_Transformed[iMother][iDaughter][iColorChannel][iCascade]->SetWeight(globalweight_daughter);
						nTupleModel_Transformed[iMother][iDaughter][iColorChannel][iCascade]->Write();
						DecayChainID[iMother][iDaughter][iColorChannel][iCascade]->Write();

					}


				}
			}
		}


		if(NRQCDvars::debug) cout<<"Calc and save dmatrix consts_star"<<endl;


		cout<<"Calc and save dmatrix consts_star"<<endl;

		double consts_star_vals[NRQCDvars::nStates][NRQCDvars::nColorChannels];
		double err_consts_star_vals[NRQCDvars::nStates][NRQCDvars::nColorChannels];

		for (int iMother=0; iMother<NRQCDvars::nStates; iMother++){
			cout<<"iMother "<<iMother<<endl;
			int nColorChannels_state;
			bool isSstate=(StateQuantumID[iMother] > NRQCDvars::quID_S)?false:true;
			if(isSstate) nColorChannels_state=NRQCDvars::nColorChannels_S;
			else nColorChannels_state=NRQCDvars::nColorChannels_P;
			for (int iColorChannel=0; iColorChannel<nColorChannels_state; iColorChannel++){
				cout<<"iColorChannel "<<iColorChannel<<endl;

				if(nTupleModelGiven[iMother]){

					sprintf(nTupleModelName,"nTupleModel_Mother%d_Daughter%d_ColorChannel%d_Cascade0",iMother, iMother, iColorChannel);

					TTree* nTupleModelDirectProd = (TTree*) ModelIngredientsFile->Get(nTupleModelName);

					//TODO: calc here the const_star values consts_star_vals[iMother][iColorChannel]

					//Loop over Tree, define asymmetric area around NRQCDvars::pT_star and NRQCDvars::rap_star
					//Fill weighted histogram, check if mean is compatible with NRQCDvars::pT_star and NRQCDvars::rap_star
					//Calc consts_star_vals[iMother][iColorChannel], check if uncertainty is below threshold
					//set consts_star_vals[iMother][iColorChannel]

					double model_pT_mother; nTupleModelDirectProd->SetBranchAddress( "model_pT",  &model_pT_mother );
					double model_rap_mother; nTupleModelDirectProd->SetBranchAddress( "model_rap",  &model_rap_mother );
					double model_costh_mother; nTupleModelDirectProd->SetBranchAddress( "model_costh",  &model_costh_mother );
					double model_phi_mother; nTupleModelDirectProd->SetBranchAddress( "model_phi",  &model_phi_mother );
					double weight_mother; nTupleModelDirectProd->SetBranchAddress( "weight",  &weight_mother );

					int n_events = int( nTupleModelDirectProd->GetEntries() );



					int nBinsh_pT=200;//TODO: fine-tune
					int nBinsh_rap=15;//TODO: fine-tune
					TH2D *h_pTrap_DirectProd = new TH2D ("h_pTrap_DirectProd", "h_pTrap_DirectProd", nBinsh_pT, model_pTMin, model_pTMax, nBinsh_rap, model_rapMin, model_rapMax);
					h_pTrap_DirectProd->Sumw2();
					int gStarBin = h_pTrap_DirectProd->FindBin(NRQCDvars::pT_star, NRQCDvars::rap_star);
					int StarBin_pT, StarBin_rap, StarBinz;
					h_pTrap_DirectProd->GetBinXYZ(gStarBin, StarBin_pT, StarBin_rap, StarBinz);
					double BinWidth_pT=(model_pTMax-model_pTMin)/double(nBinsh_pT);
					double BinWidth_rap=(model_rapMax-model_rapMin)/double(nBinsh_rap);

					double StarBin_rap_binCenter=model_rapMin+BinWidth_rap*(StarBin_rap-0.5);
					double StarBin_rap_binLowerEdge=StarBin_rap_binCenter-BinWidth_rap*0.5;
					double StarBin_rap_binUpperEdge=StarBin_rap_binCenter+BinWidth_rap*0.5;

					TH1D *h_pT_DirectProd_MeanPt[nBinsh_pT];
					char histName_MeanPt[200];
					int nBinsh_pT_MeanPt=100;
					double MeanPt_binCenter[nBinsh_pT];
					double MeanPt_rap_binLowerEdge[nBinsh_pT];
					double MeanPt_rap_binUpperEdge[nBinsh_pT];
					for (Int_t i=0; i<nBinsh_pT; i++) {
						MeanPt_binCenter[i]=model_pTMin+BinWidth_pT*(i+0.5);
						MeanPt_rap_binLowerEdge[i]=MeanPt_binCenter[i]-BinWidth_pT*0.5;
						MeanPt_rap_binUpperEdge[i]=MeanPt_binCenter[i]+BinWidth_pT*0.5;
						sprintf(histName_MeanPt,"histName_MeanPt%d",i);
						h_pT_DirectProd_MeanPt[i] = new TH1D (histName_MeanPt, histName_MeanPt, nBinsh_pT_MeanPt, MeanPt_rap_binLowerEdge[i], MeanPt_rap_binUpperEdge[i]);

					}

					// loop over  events in the model ntuple
					for ( int i_event = 1; i_event <= n_events; i_event++ ) {

						nTupleModelDirectProd->GetEvent( i_event-1 );
						h_pTrap_DirectProd->Fill(model_pT_mother, model_rap_mother, weight_mother);

						if(model_rap_mother>StarBin_rap_binLowerEdge&&model_rap_mother<StarBin_rap_binUpperEdge){
							for (Int_t i=0; i<nBinsh_pT; i++) {
								if(model_pT_mother>MeanPt_rap_binLowerEdge[i]&&model_pT_mother<MeanPt_rap_binUpperEdge[i]){
									h_pT_DirectProd_MeanPt[i]->Fill(model_pT_mother, weight_mother);
								}
							}
						}

					}

					h_pTrap_DirectProd->Scale(nTupleModelDirectProd->GetWeight());


					TGraph2D *g_pTrap_DirectProd = new TGraph2D(nBinsh_pT*nBinsh_rap);
					int k=0;
					double g_pT, g_rap, z;
					for (Int_t i=1; i<nBinsh_pT+1; i++) {
						g_pT=model_pTMin+BinWidth_pT*(i-0.5);
						for (Int_t j=1; j<nBinsh_rap+1; j++) {
							g_rap=model_rapMin+BinWidth_rap*(j-0.5);
							z = h_pTrap_DirectProd->GetBinContent(i,j);
							g_pTrap_DirectProd->SetPoint(k,g_pT,g_rap,z);
							k++;
						}
					}


					TGraph *g1D_pTrap_DirectProd;

					double pTmean_graph[nBinsh_pT];
					double model_graph[nBinsh_pT];
					for(int m=1;m<nBinsh_pT+1;m++){
						pTmean_graph[m-1]=h_pT_DirectProd_MeanPt[m-1]->GetMean();
						model_graph[m-1] = h_pTrap_DirectProd->GetBinContent(m,StarBin_rap);
					}

					g1D_pTrap_DirectProd = new TGraph(nBinsh_pT,pTmean_graph,model_graph);

					double DirectProd_previousBin=h_pTrap_DirectProd->GetBinContent(StarBin_pT-1, StarBin_rap);
					double DirectProd_nextBin=h_pTrap_DirectProd->GetBinContent(StarBin_pT+1, StarBin_rap);


					g1D_pTrap_DirectProd->Print();

					double g_consts_star_vals=g_pTrap_DirectProd->Interpolate(NRQCDvars::pT_star, NRQCDvars::rap_star)/(BinWidth_pT*BinWidth_rap);
					double h_consts_star_vals=h_pTrap_DirectProd->GetBinContent(StarBin_pT, StarBin_rap)/(BinWidth_pT*BinWidth_rap);
					double g1D_consts_star_vals=g1D_pTrap_DirectProd->Eval(NRQCDvars::pT_star)/(BinWidth_pT*BinWidth_rap);

					cout<<"g_consts_star_vals "<<g_consts_star_vals<<endl;
					cout<<"g1D_consts_star_vals "<<g1D_consts_star_vals<<endl;
					cout<<"h_consts_star_vals "<<h_consts_star_vals<<endl;
					cout<<"relative difference g/h in % "<<(g_consts_star_vals/h_consts_star_vals-1)*100<<endl;
					cout<<"relative difference g1D/h in % "<<(g1D_consts_star_vals/h_consts_star_vals-1)*100<<endl;

					h_pTrap_DirectProd->SaveAs(Form("tmp_%d.root",iColorChannel));

					//TODO: check uncertainties, propagate them

					consts_star_vals[iMother][iColorChannel]=(g1D_consts_star_vals+h_consts_star_vals)/2.;
					err_consts_star_vals[iMother][iColorChannel]=h_pTrap_DirectProd->GetBinError(StarBin_pT, StarBin_rap)/(BinWidth_pT*BinWidth_rap);

					delete g_pTrap_DirectProd;
					delete h_pTrap_DirectProd;

				}
				else{
					consts_star_vals[iMother][iColorChannel]=0.;
					err_consts_star_vals[iMother][iColorChannel]=0.;

				}
			}
		}

		// define dmatrix consts_star:::


		dmatrix consts_star(NRQCDvars::nStates);
		vector<double> consts_star_S (NRQCDvars::nColorChannels_S,0);//0...sigma_star_CS, i going from 1 to n=nColorChannels, representing c_i star for a given state and pT_star
		vector<double> consts_star_P (NRQCDvars::nColorChannels_P,0);

		dmatrix err_consts_star(NRQCDvars::nStates);
		vector<double> err_consts_star_S (NRQCDvars::nColorChannels_S,0);//0...sigma_star_CS, i going from 1 to n=nColorChannels, representing c_i star for a given state and pT_star
		vector<double> err_consts_star_P (NRQCDvars::nColorChannels_P,0);

		for (int i=0; i<NRQCDvars::nStates; i++){
			bool isSstate=(StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
			if(isSstate){
				for (int j=0; j<NRQCDvars::nColorChannels_S; j++){
					consts_star_S.at(j)=consts_star_vals[i][j];
					err_consts_star_S.at(j)=err_consts_star_vals[i][j];
				}
				consts_star.at(i)=consts_star_S;
				err_consts_star.at(i)=err_consts_star_S;
			}
			else{
				for (int j=0; j<NRQCDvars::nColorChannels_P; j++){
					consts_star_P.at(j)=consts_star_vals[i][j];
					err_consts_star_P.at(j)=err_consts_star_vals[i][j];
				}
				consts_star.at(i)=consts_star_P;
				err_consts_star.at(i)=err_consts_star_P;
			}
		}

		sprintf(outname,"%s/ModelIngredients_consts_star.txt",dirname);

	    ofstream out;
	    out.open(outname);//, std::ofstream::app);

    	out << consts_star;
    	out << err_consts_star;

	    out.close();

	    cout << consts_star << endl;
	    cout << err_consts_star << endl;


		if(NRQCDvars::debug) cout<<"Save ModelIngredients.root"<<endl;
		ModelIngredientsFile->Close();
		OriginalModelTTreeFile->Close();


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

vector<double> Transform_kinematics_MotherToDaughter(int MotherState, int DaughterState, vector<double> Kinematics_Mother) {


	  double model_pT_mother=Kinematics_Mother[0];
	  double model_rap_mother=Kinematics_Mother[1];
	  double model_costh_mother=Kinematics_Mother[2];
	  double model_phi_mother=Kinematics_Mother[3];



	  // Actual transformation of pT, rap, costh, phi

	  // pT, rap, costh, phi -> pT', rap'
	  // costh', phi': draw uniform in -1,1, -180,180


//////// Pietros Transformation Start

	    ////// INPUT PARAMETERS ///////////////////
	    //////////////////////////////////////////

	    double M  = mass[MotherState];
	    double m  = mass[DaughterState];
	    double mX = 0.0;

	    double PT = model_pT_mother;
	    double Y  = model_rap_mother;

	    double cosTheta = model_costh_mother;
	    double Phi = model_phi_mother*TMath::Pi()/180.; // radiants

	    //////////////////////////////////////////


	    double MT = sqrt( M*M + PT*PT );
	    double PL1 =   0.5 *MT * exp( Y);
	    double PL2 = - 0.5 *MT * exp(-Y);
	    double PL = PL1 + PL2;

	    // mother's 4-vector:

	    TLorentzVector mother;
	    mother.SetXYZM( PT, 0., PL, M );

	    // daughter's 4-momentum in mother's rest frame, wrt to axes of mother's helicity frame:

	    double p_mrf = 0.5 * sqrt( pow( M*M + m*m - mX*mX, 2. ) - 4.* M*M*m*m ) / M;

	    double sinTheta = sqrt( 1. - cosTheta*cosTheta );

	    TLorentzVector daughter_mrf;
	    daughter_mrf.SetXYZM( p_mrf * sinTheta * cos(Phi),
	                          p_mrf * sinTheta * sin(Phi),
	                          p_mrf * cosTheta,
	                          m );


	    // calculate daughter's 4-momentum in mother's rest frame, wrt the "xyz" laboratory axes:
	    // rotate from helicity system of axes to the "xyz" system of axes

	    TVector3 mother_direction = mother.Vect().Unit();

	    // transform (rotation) psi momentum components from polarization axis system
	    // to the system with x,y,z axes as in the laboratory
	    // It is a rotation about the Y axis (perpendicular to the production plane)

	    TVector3 oldZaxis = mother_direction;
	    TVector3 oldYaxis(0.,1.,0.);
	    TVector3 oldXaxis = oldYaxis.Cross(oldZaxis);

	    TRotation rotation;
	    rotation.RotateAxes(oldXaxis, oldYaxis, oldZaxis);
	                     // transforms coordinates from the "old" frame to the "xyz" frame

	    daughter_mrf.Transform(rotation);
	                     // daughter in mother's rest frame
	                     // wrt to the xyz axes


	    // boost daughter from the mother rest frame to the lab:

	    TVector3 mother_to_lab = mother.BoostVector();
	    TLorentzVector daughter = daughter_mrf;
	    daughter.Boost(mother_to_lab);


	    ////// OUTPUT VALUES //////////////////////
	    //////////////////////////////////////////

	    double pT = daughter.Perp();
	    double pL = daughter.Pz();
	    double y  = daughter.Rapidity();

//////// Pietros Transformation End


	  double model_pT_daughter;
	  double model_rap_daughter;
	  double model_costh_daughter;
	  double model_phi_daughter;


///////////////////////////// Actual transformation:
	  model_pT_daughter=pT;
	  model_rap_daughter=y;
	  //Generate uniform costh/phi distributions, assuming no polarization (neglecting the polarization term in the transformation formulas -> pT-rap of the following decay in same cascade can be wrong by ~1-2%)
	  model_costh_daughter=gRandom->Uniform(-1,1);
	  model_phi_daughter=gRandom->Uniform(-180,180);
//////////////////////////////////////////////////////

	  vector<double> Kinematics_Daughter;
	  Kinematics_Daughter.push_back(model_pT_daughter);
	  Kinematics_Daughter.push_back(model_rap_daughter);
	  Kinematics_Daughter.push_back(model_costh_daughter);
	  Kinematics_Daughter.push_back(model_phi_daughter);

	  return Kinematics_Daughter;
}

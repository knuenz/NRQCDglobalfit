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


double func_pT_gen(double* x, double* par);
vector<double> Transform_kinematics_MotherToDaughter(int MotherState, int DaughterState, vector<double> Kinematics_Mother);


int main(int argc, char** argv) {

  	Char_t *ModelID = "Default";
  	bool useToyModel=false;

  	for( int i=0;i < argc; ++i ) {
		if(std::string(argv[i]).find("ModelID") != std::string::npos) {char* ModelIDchar = argv[i]; char* ModelIDchar2 = strtok (ModelIDchar, "="); ModelID = ModelIDchar2; cout<<"ModelID = "<<ModelID<<endl;}
		if(std::string(argv[i]).find("useToyModel=true") != std::string::npos) {  useToyModel=true, cout<<"useToyModel"<<endl;	}
		}




	delete gRandom;
	gRandom = new TRandom3(NRQCDvars::randomSeed);  // better random generator

	double model_pTMin, model_pTMax;
	double model_rapMin, model_rapMax;
	double model_costhMin, model_costhMax;
	double model_phiMin, model_phiMax;

	model_pTMin=0;
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
	sprintf(predirname,"ModelID");
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

		if(!useToyModel){


			//Read in TTrees from ModelID/JobID folder
			for (int iMother=0; iMother<NRQCDvars::nStates; iMother++){
				/*if(NRQCDvars::debug) */cout<<"Read in original TTree model for nState = "<<iMother<<endl;
				for (int iColorChannel=0; iColorChannel<NRQCDvars::nColorChannels; iColorChannel++){
					if(NRQCDvars::debug) cout<<"		Color channel = "<<iColorChannel<<endl;
					sprintf(nTupleModelName,"nTupleModel_Mother%d_Daughter%d_ColorChannel%d_Cascade0",iMother, iMother, iColorChannel);
					nTupleModel[iMother][iColorChannel]=(TTree*)OriginalModelTTreeFile->Get(nTupleModelName);
					//cout<<"nTupleModel[iMother][iColorChannel]->GetEntries() "<<nTupleModel[iMother][iColorChannel]->GetEntries();
				}
			}


		}




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Alternatively, let's generate some toy data
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


		if(useToyModel){
			//double LuminosityPerModelEvent[NRQCDvars::nStates][NRQCDvars::nColorChannels];
			int n_nTuple=1000;
			for (int iMother=0; iMother<NRQCDvars::nStates; iMother++){
				/*if(NRQCDvars::debug) */cout<<"Generating toy model for nState = "<<iMother<<endl;
				for (int iColorChannel=0; iColorChannel<NRQCDvars::nColorChannels; iColorChannel++){
					if(NRQCDvars::debug) cout<<"		Color channel = "<<iColorChannel<<endl;

					//LuminosityPerModelEvent[iMother][iColorChannel]=100.; // Sets weight of event, to take care of the normalization of the individual states and color channels

					sprintf(nTupleModelName,"nTupleModel_Mother%d_Daughter%d_ColorChannel%d_Cascade0",iMother, iMother, iColorChannel);
					nTupleModel[iMother][iColorChannel] = new TTree(nTupleModelName, nTupleModelName);

					double model_pT;  nTupleModel[iMother][iColorChannel]->Branch("model_pT",     &model_pT,     "model_pT/D");
					double model_rap;    nTupleModel[iMother][iColorChannel]->Branch("model_rap",       &model_rap,       "model_rap/D"  );
					double model_costh;  nTupleModel[iMother][iColorChannel]->Branch("model_costh",      &model_costh,     "model_costh/D");
					double model_phi;  nTupleModel[iMother][iColorChannel]->Branch("model_phi",      &model_phi,     "model_phi/D");
					double model_weight; nTupleModel[iMother][iColorChannel]->Branch("weight",      &model_weight,     "weight/D");

					for (int k=0; k<n_nTuple; k++){

						//if(iColorChannel<2){
						//TF1* pT_distr = new TF1("pT_distr",func_pT_gen,model_pTMin,model_pTMax,1);
						//pT_distr->SetParameter(0,iMother);
						//model_pT=pT_distr->GetRandom();
						//}
						//else
							model_pT=gRandom->Uniform(model_pTMin,model_pTMax);


						model_rap=gRandom->Uniform(model_rapMin,model_rapMax);
						model_costh=gRandom->Uniform(model_costhMin,model_costhMax);
						model_phi=gRandom->Uniform(model_phiMin,model_phiMax);
						model_weight=1.;

						nTupleModel[iMother][iColorChannel]->Fill();

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
				/*if(NRQCDvars::debug) */cout<<"Save original TTree model for nState = "<<iMother<<endl;
				for (int iColorChannel=0; iColorChannel<NRQCDvars::nColorChannels; iColorChannel++){
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
			for (int iDaughter=0; iDaughter<iMother; iDaughter++){

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

				if(Buffer_DecayChains.size()<1) continue; // continue, if there is not a single deca chain linkin iMother -> iDaughter

				double FullBR=0;
				for(int iDecayChain=0;iDecayChain<Buffer_DecayChainBranchingRatios.size();iDecayChain++){
					FullBR+=Buffer_DecayChainBranchingRatios[iDecayChain];
				}
				for(int iDecayChain=0;iDecayChain<Buffer_DecayChainBranchingRatios.size();iDecayChain++){
					Buffer_DecayChainBranchingRatioFractions.push_back(Buffer_DecayChainBranchingRatios[iDecayChain]/FullBR);
				}


				for (int iColorChannel=0; iColorChannel<NRQCDvars::nColorChannels; iColorChannel++){

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

						double globalweight_daughter=nTupleModel[iMother][iColorChannel]->GetWeight()*Buffer_DecayChainBranchingRatios[iCascade]; // Luminosity/event times product of Branching ratios of the full decay chain
						nTupleModel_Transformed[iMother][iDaughter][iColorChannel][iCascade]->SetWeight(globalweight_daughter);
						nTupleModel_Transformed[iMother][iDaughter][iColorChannel][iCascade]->Write();
						DecayChainID[iMother][iDaughter][iColorChannel][iCascade]->Write();

					}


				}
			}
		}

		if(NRQCDvars::debug) cout<<"Save ModelIngredients.root"<<endl;
		ModelIngredientsFile->Close();
		OriginalModelTTreeFile->Close();


return 0;

}


double func_pT_gen(double* x, double* par) {

	double beta;
	double pTsq;

	//if(par[0]==0) {beta = 3.69; pTsq = 12.0;}  // Psi(1S)
	//if(par[0]==3) {beta = 3.71; pTsq = 19.5;}  // Psi(2S)
	//if(par[0]==4) {beta = 3.46; pTsq = 47.3;}  // Upsi(1S)
	//if(par[0]==7) {beta = 3.27; pTsq = 65.7;}  // Upsi(2S)
	//if(par[0]==10) {beta = 3.05; pTsq = 80.5;}  // Upsi(3S)

	int nState=par[0];
	beta=3;
	pTsq=NRQCDvars::mass[nState]*NRQCDvars::mass[nState];

	return x[0] * pow( 1. + 1./(beta - 2.) * x[0]*x[0] / pTsq, -beta  );

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

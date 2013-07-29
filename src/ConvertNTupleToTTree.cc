/*
 * ConvertNTupleToTTree.cc
 *
 *  Created on: Jul 12, 2013
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
#include <fstream>

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


int main(int argc, char** argv) {

  	Char_t *OriginalNTupleID = "Default";
  	Char_t *ModelID = "Default";
	Char_t *storagedir = "Default"; //Storage Directory

  	for( int i=0;i < argc; ++i ) {
		if(std::string(argv[i]).find("OriginalNTupleID") != std::string::npos) {char* OriginalNTupleIDchar = argv[i]; char* OriginalNTupleIDchar2 = strtok (OriginalNTupleIDchar, "="); OriginalNTupleID = OriginalNTupleIDchar2; cout<<"OriginalNTupleID = "<<OriginalNTupleID<<endl;}
		if(std::string(argv[i]).find("ModelID") != std::string::npos) {char* ModelIDchar = argv[i]; char* ModelIDchar2 = strtok (ModelIDchar, "="); ModelID = ModelIDchar2; cout<<"ModelID = "<<ModelID<<endl;}
		if(std::string(argv[i]).find("storagedir") != std::string::npos) {char* storagedirchar = argv[i]; char* storagedirchar2 = strtok (storagedirchar, "="); storagedir = storagedirchar2; cout<<"storagedir = "<<storagedir<<endl;}
  	}


	char outname[2000];
	char inname[2000];
	char predirname[2000];
	char dirname[2000];


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Read input model
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	sprintf(predirname,"%s/ModelID", storagedir);
	gSystem->mkdir(predirname);
	sprintf(dirname,"%s/%s",predirname,ModelID);
	gSystem->mkdir(dirname);

	sprintf(outname,"%s/OriginalModelTTree.root",dirname);
	TFile* ModelIngredientsFile = new TFile(outname, "RECREATE");


	TTree* nTupleModel[NRQCDvars::nStates][NRQCDvars::nColorChannels];
	char nTupleModelName[1000];

	for (int iMother=0; iMother<NRQCDvars::nStates; iMother++){
		cout<<"Converting original model nTuple for nState = "<<iMother<<endl;
		int nColorChannels_state;
		bool isSstate=(StateQuantumID[iMother] > NRQCDvars::quID_S)?false:true;
		if(isSstate) nColorChannels_state=NRQCDvars::nColorChannels_S;
		else nColorChannels_state=NRQCDvars::nColorChannels_P;
		for (int iColorChannel=0; iColorChannel<nColorChannels_state; iColorChannel++){
			cout<<"		Color channel = "<<iColorChannel<<endl;

			sprintf(predirname,"OriginalNTupleID");
			gSystem->mkdir(predirname);
			sprintf(dirname,"%s/%s",predirname,OriginalNTupleID);
			gSystem->mkdir(dirname);
			sprintf(inname,"%s/ups1s_hx.outpt",dirname);
			FILE *fIn;
			fIn = fopen(inname, "read");
			cout<<inname<<endl;
			cout<<"debug"<<endl;

			sprintf(nTupleModelName,"nTupleModel_Mother%d_Daughter%d_ColorChannel%d_Cascade0",iMother, iMother, iColorChannel);
			nTupleModel[iMother][iColorChannel] = new TTree(nTupleModelName, nTupleModelName);

			Double_t model_pT;  nTupleModel[iMother][iColorChannel]->Branch("model_pT",     &model_pT,     "model_pT/D");
			Double_t model_rap;    nTupleModel[iMother][iColorChannel]->Branch("model_rap",       &model_rap,       "model_rap/D"  );
			Double_t model_costh;  nTupleModel[iMother][iColorChannel]->Branch("model_costh",      &model_costh,     "model_costh/D");
			Double_t model_phi;  nTupleModel[iMother][iColorChannel]->Branch("model_phi",      &model_phi,     "model_phi/D");
			Double_t weight; nTupleModel[iMother][iColorChannel]->Branch("weight",      &weight,     "weight/D");

			int iEvents=0;
			while ( !(feof(fIn)) ){
				if(feof(fIn)) {cout<<"break"<<endl; break;}
				iEvents++;

				Char_t line[1000];
				fgets(line, sizeof(line), fIn); //comment
				Float_t in_model_pT, in_model_rap, in_model_costh, in_model_phi, in_weight;
				sscanf(line, "%f %f %f %f %f", &in_model_pT, &in_model_rap, &in_model_costh, &in_model_phi, &in_weight);

				model_pT=in_model_pT;
				model_rap=in_model_rap;
				model_costh=in_model_costh;
				model_phi=in_model_phi;
				weight=in_weight;

				nTupleModel[iMother][iColorChannel]->Fill();

				if(iEvents<10) cout<<model_pT<<" "<<model_rap<<" "<<model_costh<<" "<<model_phi<<" "<<weight<<" "<<endl;

			}

			cout<<"Read in "<<iEvents<<" events"<<endl;
			cout<<"nTupleModel[iMother][iColorChannel]->GetEntries() "<<nTupleModel[iMother][iColorChannel]->GetEntries()<<endl;


			nTupleModel[iMother][iColorChannel]->Write();
			fclose(fIn);


		}
	}


	ModelIngredientsFile->Close();

  	return 0;

  	}



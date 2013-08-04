/*
 * MinuitLikelihoodFunction.h
 *
 *  Created on: Jul 31, 2013
 *      Author: valentinknuenz
 */

//#ifndef MINUITLIKELIHOODFUNCTION_H_
//#define MINUITLIKELIHOODFUNCTION_H_

#ifndef MN_FcnPVLogL_H_
#define MN_FcnPVLogL_H_
#include "Minuit2/FCNBase.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TGraphAsymmErrors.h"
#include <iostream>
#include <string>
#include <sstream>

#include "TRandom3.h"
//#include "../interface/NRQCDglobalfitObject.h"
//#include "../src/NRQCDglobalfitObject.cc"
#include "../interface/GlobalVar.h"


using namespace std;

#if !defined(__CINT__) && !defined(__MAKECINT__)

#include <vector>
#endif

class FcnPVLogL : public ROOT::Minuit2::FCNBase {
public:
//   FcnPVLogL() {}
FcnPVLogL( bool SampleNp_, bool SampleNp_consts_star_, vector<NRQCDglobalfitObject> DataModelObject_, dmatrix consts_star_, dmatrix err_consts_star_, dvector Np_expected_, dvector Np_uncertainty_) : DataModelObject(DataModelObject_), consts_star(consts_star_), err_consts_star(err_consts_star_), SampleNp(SampleNp_), SampleNp_consts_star(SampleNp_consts_star_), Np_expected(Np_expected_), Np_uncertainty(Np_uncertainty_), theErrorDef(1.) {}
 ~FcnPVLogL() {}
 double Up() const {return theErrorDef;}
 double operator()(const std::vector<double>&) const;
private:

 bool SampleNp;
 bool SampleNp_consts_star;
 vector<NRQCDglobalfitObject> DataModelObject;
 dmatrix consts_star;
 dmatrix err_consts_star;
 dvector Np_uncertainty;
 dvector Np_expected;


 double theErrorDef;

};
#endif //MN_FcnPVLogL_H_

double
FcnPVLogL::operator() (const std::vector<double>& pars) const
{






	double loglikelihood;//return this

	delete gRandom;
	gRandom = new TRandom3(0);  // better random generator



//////////////////////////
// INIZIALIZE vars

	dvector ObjectLikelihoodVec(2,0);
	dcube directProductionCube;
	dmatrix promptProductionMatrix;
	double polCorrFactor;

	dmatrix Np_BR; //Branching ratios, [nDauthgers][nMothers]
	dmatrix Np_US; //Uncertainty scales, [0=Data, 1=Model][nScales]

	dvector Np_BR_0 (NRQCDvars::nStates,0.);
	for(int i=0; i < NRQCDvars::nStates; i++) Np_BR.push_back(Np_BR_0);

	dvector Np_US_0 (NRQCDvars::nDataSystematicScales, 0.);
	dvector Np_US_1 (NRQCDvars::nModelSystematicScales, 0.);
	Np_US.push_back(Np_US_0);
	Np_US.push_back(Np_US_1);

	//Observable parameters Op (Matrix elements)
	dmatrix Op(NRQCDvars::nStates);
	vector<double> Op_S (NRQCDvars::nColorChannels_S,0);
	vector<double> Op_P (NRQCDvars::nColorChannels_P,0);

	dmatrix Fractions(NRQCDvars::nStates);
	vector<double> Fractions_S (NRQCDvars::nColorChannels_S,0);//f0...R, fi: i going from 1 to n=nColorChannels, fn=1-sum(fi_i-(n-1))
	vector<double> Fractions_P (NRQCDvars::nColorChannels_P,0);

	//const_star matrix varied by uncertainty (Nuisance parameter inthe sampling)
	dmatrix consts_star_var(NRQCDvars::nStates);
	vector<double> consts_star_var_S (NRQCDvars::nColorChannels_S,0);
	vector<double> consts_star_var_P (NRQCDvars::nColorChannels_P,0);



//////////////////////////
// CONNECT pars with Fractions

	int nPOI=0;
	//Normalize fractions such that they sum up to 1
	double sum_fi[NRQCDvars::nStates];
	int iPar=0;
	for (int i=0; i<NRQCDvars::nStates; i++){
		int nColorChannels_state;
		bool isSstate=(NRQCDvars::StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
		if(isSstate) nColorChannels_state=NRQCDvars::nColorChannels_S;
		else nColorChannels_state=NRQCDvars::nColorChannels_P;
		//cout<<"state "<<i<<", nColorChannels_state "<<nColorChannels_state<<endl;
		sum_fi[i]=0.;
		iPar++;
		for (int j=1; j<nColorChannels_state-1; j++){
			sum_fi[i]+=pars[iPar];
			iPar++;
		}
		iPar++;
		//cout<<"sum_fi["<<i<<"]="<<sum_fi[i]<<endl;
	}
	nPOI=iPar;

	iPar=0;
	for (int i=0; i<NRQCDvars::nStates; i++){
		bool isSstate=(NRQCDvars::StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
		if(isSstate){
			for (int j=0; j<NRQCDvars::nColorChannels_S; j++){
				if(j<NRQCDvars::nColorChannels_S-1) Fractions_S.at(j)=pars[iPar];
				else Fractions_S.at(j)=1-sum_fi[i];
				iPar++;
			}
			Fractions.at(i)=Fractions_S;
		}
		else{
			for (int j=0; j<NRQCDvars::nColorChannels_P; j++){
				if(j<NRQCDvars::nColorChannels_P-1) Fractions_P.at(j)=pars[iPar];
				else Fractions_P.at(j)=1-sum_fi[i];
				iPar++;
			}
			Fractions.at(i)=Fractions_P;
		}
	}

	//cout<<Fractions<<endl;


//////////////////////////
// VARY Nps

	int iPar_Nuis=0;
	for(int j=0; j < NRQCDvars::nDataSystematicScales; j++){
		if(SampleNp) {Np_US[0][j]=pars[nPOI+iPar_Nuis]; iPar_Nuis++;} // Luminosity scaling
		else Np_US[0][j]=0.;
	}
	for(int j=0; j < NRQCDvars::nModelSystematicScales; j++){
		if(SampleNp) {Np_US[1][j]=pars[nPOI+iPar_Nuis]; iPar_Nuis++;} // Luminosity scaling
		else Np_US[1][j]=0.;
	}
	for(int i=0; i < NRQCDvars::nStates; i++){
		for(int j=0; j < NRQCDvars::nStates; j++){
			if(NRQCDvars::FeedDownBranchingRatio[i][j]>0){
				if(SampleNp) {Np_BR[i][j]=pars[nPOI+iPar_Nuis]; iPar_Nuis++;}
				else Np_BR[i][j]=NRQCDvars::FeedDownBranchingRatio[i][j]*0.01;
			}
			else Np_BR[i][j]=0;
		}
	}

	for (int i=0; i<NRQCDvars::nStates; i++){
		bool isSstate=(NRQCDvars::StateQuantumID[i] > NRQCDvars::quID_S)?false:true;
		if(isSstate){
			for (int j=0; j<NRQCDvars::nColorChannels_S; j++){
				if(SampleNp_consts_star) {consts_star_var_S.at(j)=pars[nPOI+iPar_Nuis]; iPar_Nuis++;}
				else consts_star_var_S.at(j)=consts_star[i][j];
			}
			consts_star_var.at(i)=consts_star_var_S;
		}
		else{
			for (int j=0; j<NRQCDvars::nColorChannels_P; j++){
				if(SampleNp_consts_star) {consts_star_var_P.at(j)=pars[nPOI+iPar_Nuis]; iPar_Nuis++;}
				else consts_star_var_P.at(j)=consts_star[i][j];
			}
			consts_star_var.at(i)=consts_star_var_P;
		}
	}

	const int nPar_Nuis=iPar_Nuis;



	// Transform Fractions to Op's

		int iMax = Fractions.size();
		for(int i=0; i < iMax; i++){
			//cout<<"i "<<i<<endl;
			dvector Op_state;
			Op_state.push_back(NRQCDvars::ColorSingletME[i]);
			for(dvector::iterator j = Fractions[i].begin()+1; j != Fractions[i].end(); ++j){
				int k=j-Fractions[i].begin();
				//cout<<"k "<<k<<endl;
				Op_state.push_back(Fractions[i][k]*Fractions[i][0]*consts_star_var[i][0]*NRQCDvars::ColorSingletME[i]/consts_star_var[i][k]);
			}
			Op.at(i)=Op_state;
		}



//////////////////////////
// CALC likelihood


	if(NRQCDvars::debug) cout << "getObjectLikelihood for Minuit" << endl;

	loglikelihood=0;
//	for(vector< NRQCDglobalfitObject >::iterator state = *DataModelObject.begin(); state != *DataModelObject.end(); ++state){

	int nData=DataModelObject.size();
	for(int iData=0;iData<nData;iData++){
		//cout << "getObjectLikelihood" <<endl;
		//state->Dump(NRQCDvars::nStates, true, true);
		//int nOp=Op[0].size();
		//cout<<"nOp "<<nOp<<endl;

		NRQCDglobalfitObject state=DataModelObject.at(iData);
		ObjectLikelihoodVec=state.getObjectLikelihood(Op, Np_BR, Np_US, false, directProductionCube, promptProductionMatrix, polCorrFactor);
		loglikelihood+=ObjectLikelihoodVec[0];
		//cout << "likelihood: " << likelihood << endl;
	}


	//cout << "Sum_likelihood " << likelihood << endl;

	loglikelihood*=0.5;

	//add constraints of nuisance parameters to likelihood
	for(int i=0;i<nPar_Nuis;i++){

		loglikelihood+=((pars[nPOI+i]-Np_expected[i])/Np_uncertainty[i])*((pars[nPOI+i]-Np_expected[i])/Np_uncertainty[i]);

	}


	return loglikelihood;
}


//#endif /* MINUITLIKELIHOODFUNCTION_H_ */

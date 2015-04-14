#include "../interface/NRQCDglobalfitObject.h"
#include "../interface/GlobalVar.h"

#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include "stdio.h"

#include <TMath.h>



//#define mydebug


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//NRQCDglobalfitObject::NRQCDglobalfitObject(){}

//void NRQCDglobalfitObject::initFractions(const int& nDim){
//	fraction = newCVector<double>(nDim);
//  	sampleValues = newCVector<double>(NRQCDvars::nFractionDim-1);       // TODO (?): nFractionDim from object
//  	--sampleValues;                                                     // index starts at 1
//  	sampleWidths = newCVector<double>(NRQCDvars::nFractionDim-1);       //
//  	--sampleWidths;                                                     // index starts at 1
//}

//NRQCDglobalfitObject::NRQCDglobalfitObject(const int& numOperatorStates){
//	this->setFractionsBasis(const int& numOperatorStates);
//}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//NRQCDglobalfitObject::NRQCDglobalfitObject(): Object_nState( N_STATES ), Object_CentralValue(0), Object_ErrStatPos(0), Object_ErrStatNeg(0),
//                                              Object_ErrSystPos(0), Object_ErrSystNeg(0), Object_ErrTotPos(0), Object_ErrTotNeg(0),
//                                              Object_ErrGlobalPos(0), Object_ErrGlobalNeg(0), Object_yMin(0), Object_yMax(0), Object_yMean(0),
//                                              Object_pTMin(0), Object_pTMax(0), Object_pTMean(0)
//{
//	  Object_ShortDistanceCoef = this->dmatrix_allocate(1, 1);
//	  Object_OctetCompLamth = this->dmatrix_allocate(1, 1);
//	  Object_OctetCompLamph = this->dmatrix_allocate(1, 1);
//	  Object_OctetCompLamtp = this->dmatrix_allocate(1, 1);
//}
//
//NRQCDglobalfitObject::NRQCDglobalfitObject(const int nrows, const int ncolumns): Object_nState( N_STATES ), Object_CentralValue(0), Object_ErrStatPos(0),
//		                                                                         Object_ErrStatNeg(0), Object_ErrSystPos(0), Object_ErrSystNeg(0),
//		                                                                         Object_ErrTotPos(0), Object_ErrTotNeg(0), Object_ErrGlobalPos(0),
//		                                                                         Object_ErrGlobalNeg(0), Object_yMin(0), Object_yMax(0),
//		                                                                         Object_yMean(0), Object_pTMin(0), Object_pTMax(0), Object_pTMean(0)
//{
//	  Object_ShortDistanceCoef = this->dmatrix_allocate(nrows, ncolumns);
//	  Object_OctetCompLamth = this->dmatrix_allocate(nrows, ncolumns);
//	  Object_OctetCompLamph = this->dmatrix_allocate(nrows, ncolumns);
//	  Object_OctetCompLamtp = this->dmatrix_allocate(nrows, ncolumns);
//}
//
//NRQCDglobalfitObject::NRQCDglobalfitObject(const NRQCDglobalfitObject &other)
//{
//	// Measurement members
//	  Object_nState = other.Object_nState;
//	  Object_CentralValue = other.Object_CentralValue;
//	  Object_ErrStatPos = other.Object_ErrStatPos;
//	  Object_ErrStatNeg = other.Object_ErrStatNeg;
//	  Object_ErrSystPos = other.Object_ErrSystPos;
//	  Object_ErrSystNeg = other.Object_ErrSystNeg;
//	  Object_ErrTotPos = other.Object_ErrTotPos;
//	  Object_ErrTotNeg = other.Object_ErrTotNeg;
//	  Object_ErrGlobalPos = other.Object_ErrGlobalPos;
//	  Object_ErrGlobalNeg = other.Object_ErrGlobalNeg;
//	  Object_yMin = other.Object_yMin;
//	  Object_yMax = other.Object_yMax;
//	  Object_yMean = other.Object_yMean;
//	  Object_pTMin = other.Object_pTMin;
//	  Object_pTMax = other.Object_pTMax;
//	  Object_pTMean = other.Object_pTMean;
//	  Object_PolCorrParams = other.Object_PolCorrParams;
//	  Object_ObjectID = other.Object_ObjectID;
//
//	//Model members
//	  dmatrix_copy(other.Object_ShortDistanceCoef, Object_ShortDistanceCoef);
//	  dmatrix_copy(other.Object_OctetCompLamth, Object_OctetCompLamth);
//	  dmatrix_copy(other.Object_OctetCompLamph, Object_OctetCompLamph);
//	  dmatrix_copy(other.Object_OctetCompLamtp, Object_OctetCompLamtp);
//}

dvector NRQCDglobalfitObject::getDirectProduction(int nState, dmatrix  &Op, dmatrix &Np_BR, dmatrix &Np_US, bool returnMPDetails, dmatrix &directProductionMatrix) {
//	cout<<"NRQCDglobalfitObject::getDirectProduction of state "<<nState_<<endl;


	dvector DirectProductionVec;//return this

	int nOp=Op[nState].size();
	int nModelSystematicScales_ = Np_US[1].size();
	//cout<<"nOp "<<nOp<<endl;
	//cout<<"nModelSystematicScales_ "<<nModelSystematicScales_<<endl;

	//calculate direct production cross section (sum of individual octet contributions)
	dvector ShortDistanceCoef_ = getShortDistanceCoef(nState);
	dvector ShortDistanceCoef_corrected (nOp,0.);

	//Apply the corrections to each ColorChannel-element: global systematics

	dmatrix  ShortDistanceCoef_globalSystPos = getShortDistanceCoef_globalSystPos(nState);
	dmatrix  ShortDistanceCoef_globalSystNeg = getShortDistanceCoef_globalSystNeg(nState);

	for(int i=0; i<nOp ; i++){
		double sigma_Central=ShortDistanceCoef_[i];
		ShortDistanceCoef_corrected[i]=sigma_Central;
		for(int j=0; j<nModelSystematicScales_ ; j++){
			double sigma_Pos=ShortDistanceCoef_globalSystPos[i][j]; //if(j==0) sigma_Pos=sigma_Central+(sigma_Pos-sigma_Central)*10;
			double sigma_Neg=ShortDistanceCoef_globalSystNeg[i][j];
			double RandomUS=Np_US[1][j];
			double GlobalSystCorrectionFactor=1+(sigma_Pos-sigma_Neg)/(2.*sigma_Central)*RandomUS+((sigma_Pos+sigma_Neg)/(2.*sigma_Central)-1)*RandomUS*RandomUS;

			//old bugged version
			//if(RandomUS>0) GlobalSystCorrectionFactor=sigma_Pos/sigma_Central*RandomUS;
			//if(RandomUS<0) GlobalSystCorrectionFactor=sigma_Neg/sigma_Central*(-1.*RandomUS);

			//new debugged version
			if(RandomUS>0) GlobalSystCorrectionFactor=1+(sigma_Pos/sigma_Central-1)*RandomUS;
			if(RandomUS<0) GlobalSystCorrectionFactor=1+(sigma_Neg/sigma_Central-1)*(-1.*RandomUS);


			if(i==j) ShortDistanceCoef_corrected[i]*=GlobalSystCorrectionFactor;

			//if(i==j){
			//	cout<<"CC"<<j<<endl;
			//	cout<<"sigma_Central "<<sigma_Central<<endl;
			//	cout<<"sigma_Pos "<<ShortDistanceCoef_globalSystPos[i][j]<<endl;
			//	cout<<"sigma_Neg "<<ShortDistanceCoef_globalSystNeg[i][j]<<endl;
			//	cout<<"RandomUS "<<RandomUS<<endl;
			//	cout<<"GlobalSystCorrectionFactor "<<GlobalSystCorrectionFactor<<endl;
			//	cout<<"sigma_Corr "<<ShortDistanceCoef_corrected[i]<<endl;
			//	cout<<" "<<endl;
			//}
}
		//cout<<"ShortDistanceCoef_corrected["<<i<<"] "<<ShortDistanceCoef_corrected[i]<<endl;
	}

	double DirectCrossSect=0;
	dvector DirectCrossSectCont;
	for (int i=0 ; i<nOp ; i++){
		DirectCrossSectCont.push_back(ShortDistanceCoef_corrected[i]*Op[nState][i]);
		DirectCrossSect+=DirectCrossSectCont[i];
		//cout<<"DirectCrossSect[i] "<<DirectCrossSect<<endl;
	}

	//calculate fractions of individual octet contributions
	dvector DirectOctetContFraction;
	for (int i=0 ; i<nOp ; i++){
		DirectOctetContFraction.push_back(DirectCrossSectCont[i]/DirectCrossSect);
	}




	//calculate the polarization of the direct production
	dvector OctetCompLamth_ = getOctetCompLamth(nState);
	dvector OctetCompLamph_ = getOctetCompLamph(nState);
	dvector OctetCompLamtp_ = getOctetCompLamtp(nState);
	dvector OctetCompLamth_corrected (nOp,0.);
	dvector OctetCompLamph_corrected (nOp,0.);
	dvector OctetCompLamtp_corrected (nOp,0.);

	//Apply the corrections to each ColorChannel-element: global systematics

	dmatrix  OctetCompLamth_globalSystPos = getOctetCompLamth_globalSystPos(nState);
	dmatrix  OctetCompLamth_globalSystNeg = getOctetCompLamth_globalSystNeg(nState);
	dmatrix  OctetCompLamph_globalSystPos = getOctetCompLamph_globalSystPos(nState);
	dmatrix  OctetCompLamph_globalSystNeg = getOctetCompLamph_globalSystNeg(nState);
	dmatrix  OctetCompLamtp_globalSystPos = getOctetCompLamtp_globalSystPos(nState);
	dmatrix  OctetCompLamtp_globalSystNeg = getOctetCompLamtp_globalSystNeg(nState);

	for(int i=0; i<nOp ; i++){
		double Lamth_Central=OctetCompLamth_[i];
		double Lamph_Central=OctetCompLamph_[i];
		double Lamtp_Central=OctetCompLamtp_[i];
		OctetCompLamth_corrected[i]=Lamth_Central;
		OctetCompLamph_corrected[i]=Lamph_Central;
		OctetCompLamtp_corrected[i]=Lamtp_Central;
		for(int j=0 ; j<nModelSystematicScales_ ; j++){

			double Lamth_Pos=OctetCompLamth_globalSystPos[i][j];
			double Lamth_Neg=OctetCompLamth_globalSystNeg[i][j];
			double Lamph_Pos=OctetCompLamph_globalSystPos[i][j];
			double Lamph_Neg=OctetCompLamph_globalSystNeg[i][j];
			double Lamtp_Pos=OctetCompLamtp_globalSystPos[i][j];
			double Lamtp_Neg=OctetCompLamtp_globalSystNeg[i][j];
			double RandomUS=Np_US[1][j];

			//cout<<"i "<<i<<endl;
				//cout<<"j "<<j<<endl;

			//if(i==j && j==0) {
			//cout<<"Lamth_Pos "<<Lamth_Pos<<endl;
			//cout<<"Lamth_Neg "<<Lamth_Neg<<endl;
			//}

			//Old lamth correction - did not work
			double GlobalLamthCorrectionTerm=1+(Lamth_Pos-Lamth_Neg)/(2.)*RandomUS+((Lamth_Pos+Lamth_Neg)/(2.)-Lamth_Central)*RandomUS*RandomUS;
			double GlobalLamphCorrectionTerm=1+(Lamph_Pos-Lamph_Neg)/(2.)*RandomUS+((Lamph_Pos+Lamph_Neg)/(2.)-Lamph_Central)*RandomUS*RandomUS;
			double GlobalLamtpCorrectionTerm=1+(Lamtp_Pos-Lamtp_Neg)/(2.)*RandomUS+((Lamtp_Pos+Lamtp_Neg)/(2.)-Lamtp_Central)*RandomUS*RandomUS;

			//OctetCompLamth_corrected[i]*=GlobalLamthCorrectionTerm;
			//OctetCompLamph_corrected[i]*=GlobalLamphCorrectionTerm;
			//OctetCompLamtp_corrected[i]*=GlobalLamtpCorrectionTerm;

			//New lamth correction: lam'=kx+d, assuming symmetric LampPos and LamNeg
			if(i==j) OctetCompLamth_corrected[i] = (Lamth_Pos-Lamth_Neg)/(2.) * RandomUS - (Lamth_Pos-Lamth_Neg)/(2.)+Lamth_Pos;
			if(i==j) OctetCompLamph_corrected[i] = (Lamph_Pos-Lamph_Neg)/(2.) * RandomUS - (Lamph_Pos-Lamph_Neg)/(2.)+Lamph_Pos;
			if(i==j) OctetCompLamtp_corrected[i] = (Lamtp_Pos-Lamtp_Neg)/(2.) * RandomUS - (Lamtp_Pos-Lamtp_Neg)/(2.)+Lamtp_Pos;

			//if(i==j && j==0) {
			//cout<<"Lamth_Central "<<Lamth_Central<<endl;
			//cout<<"OctetCompLamth_corrected "<<OctetCompLamth_corrected[i]<<endl;
			//}
		}

	}



#ifdef mydebug
	for(int i=0 ; i<nOp ; i++){
		cout<<"DirectOctetContFraction[ "<<i<<"] = "<<DirectOctetContFraction[i]<<endl;
		cout<<"OctetCompLamth_corrected[ "<<i<<"] = "<<OctetCompLamth_corrected[i]<<endl;
		cout<<"OctetCompLamph_corrected[ "<<i<<"] = "<<OctetCompLamph_corrected[i]<<endl;
		cout<<"OctetCompLamtp_corrected[ "<<i<<"] = "<<OctetCompLamtp_corrected[i]<<endl;
	}
#endif


	double DirectLamth_numerator=0;
	double DirectLamth_denominator=0;
	for(int i=0 ; i<nOp ; i++){
		DirectLamth_numerator+=DirectOctetContFraction[i]*OctetCompLamth_corrected[i]/(3+OctetCompLamth_corrected[i]);
		DirectLamth_denominator+=DirectOctetContFraction[i]/(3+OctetCompLamth_corrected[i]);
	}

	double DirectLamph_numerator=0;
	for(int i=0 ; i<nOp ; i++){
		DirectLamph_numerator+=DirectOctetContFraction[i]*OctetCompLamph_corrected[i]/(3+OctetCompLamth_corrected[i]);
	}

	double DirectLamtp_numerator=0;
	for(int i=0 ; i<nOp ; i++){
		DirectLamtp_numerator+=DirectOctetContFraction[i]*OctetCompLamtp_corrected[i]/(3+OctetCompLamth_corrected[i]);
	}

	double DirectLamth=DirectLamth_numerator/DirectLamth_denominator;
	double DirectLamph=DirectLamph_numerator/DirectLamth_denominator;
	double DirectLamtp=DirectLamtp_numerator/DirectLamth_denominator;


	//set the return vector
	DirectProductionVec.push_back(DirectCrossSect);
	DirectProductionVec.push_back(DirectLamth);
	DirectProductionVec.push_back(DirectLamph);
	DirectProductionVec.push_back(DirectLamtp);

#ifdef mydebug
	cout<<"DirectCrossSect "<<DirectCrossSect<<endl;
	cout<<"DirectLamth "<<DirectLamth<<endl;
	cout<<"DirectLamph "<<DirectLamph<<endl;
	cout<<"DirectLamtp "<<DirectLamtp<<endl;
#endif


	if(returnMPDetails){
		directProductionMatrix.push_back(DirectCrossSectCont);
		directProductionMatrix.push_back(OctetCompLamth_corrected);
		directProductionMatrix.push_back(OctetCompLamph_corrected);
		directProductionMatrix.push_back(OctetCompLamtp_corrected);
	}

	//cout<<"OctetCompLamth_corrected[0]"<<endl;
	//cout<<OctetCompLamth_corrected[0]<<endl;

	return DirectProductionVec;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

dvector NRQCDglobalfitObject::getPromptProduction(int state, dmatrix  &Op, dmatrix &Np_BR, dmatrix &Np_US, bool returnMPDetails, dcube &directProductionCube, dmatrix &promptProductionMatrix)  {
//	cout<<"NRQCDglobalfitObject::getPromptProduction"<<endl;

	dvector PromptProductionVec;//return this


	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	//   Make a vector ListFeedDownStates (NULL) containing a list of int nFeedDownStates for which the directProduction needs to be calculated
	//   Make a vector ListFeedDownStateBranchingFractions (NULL) containing a list of double BranchingRatio for each nFeedDownStates, used as weight for the directProduction results
	//   :::
	//   Go through the list of FeedDownBranchings[state-1][i], starting from index i=state.
	//   For each non-zero element i
	//	  	-> fill the vector ListFeedDownStates with i+1
	//	  	-> fill ListFeedDownStateBranchingFractions with the BR identified with the element i itself
	//	  	-> go through FeedDownBranchings[i][j], starting from index j=i+1.
	//   		-> For each non-zero element j
	//	  		-> check if the element of FeedDownStates[state-1][j] is also non-zero. If it is zero, fill ListFeedDownStates with j+1.
	//	  		-> add the corresponding BR of element j to the j'th entry of ListFeedDownStateBranchingFractions as FeedDownBranchings[state-1][i]*FeedDownBranchings[i][j]
	//	  		-> go through FeedDownBranchings[j][k], starting from index k=j+1.
	//	  			-> repeat previous loop (until nCascade reached)
	//
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	// Initialize the arrays ListFeedDownStates, ListFeedDownStateBranchingFractions (and set entries 0)
	//int ListFeedDownStates[states];
	double *ListFeedDownStateBranchingFractions = newCVector<double>(NRQCDvars::nStates);
//	for(int i=0 ; i<nStates ; i++){
//		//ListFeedDownStates[i]=0;
//		ListFeedDownStateBranchingFractions[i]=0.;
//	}

	// Fill ListFeedDownStates, ListFeedDownStateBranchingFractions
	for(int i=state; i < NRQCDvars::nStates; i++){
		// Loop over the states decaying into state
		if(Np_BR[state][i]>0){
			ListFeedDownStateBranchingFractions[i]+=Np_BR[state][i];
			#ifdef mydebug
			cout<<"Branching "<<i<<"->"<<state<<" to state: Np_BR = "<<Np_BR[state][i]<<endl;
			#endif
			// Add Cascade decays to the prompt components (maximum considered: cascade with 5 decays)
			for(int j=i+1; j < NRQCDvars::nStates; j++){
				if(Np_BR[i][j]>0){
					ListFeedDownStateBranchingFractions[j]+=Np_BR[state][i]*Np_BR[i][j];
					#ifdef mydebug
					cout<<"Branching "<<j<<"->"<<state<<" to state: Np_BR = "<<Np_BR[state][i]*Np_BR[i][j]<<" intermediate state: "<<i<<endl;
					#endif
					for(int k=j+1; k < NRQCDvars::nStates; k++){
						if(Np_BR[j][k]>0){
							ListFeedDownStateBranchingFractions[k]+=Np_BR[state][i]*Np_BR[i][j]*Np_BR[j][k];
							#ifdef mydebug
							cout<<"Branching "<<k<<"->"<<state<<" to state: Np_BR = "<<Np_BR[state][i]*Np_BR[i][j]*Np_BR[j][k]<<" intermediate states: "<<j<<", "<<i<<endl;
							#endif
							for(int l=k+1; l < NRQCDvars::nStates; l++){
								if(Np_BR[k][l]>0){
									ListFeedDownStateBranchingFractions[l]+=Np_BR[state][i]*Np_BR[i][j]*Np_BR[j][k]*Np_BR[k][l];
									#ifdef mydebug
									cout<<"Branching "<<l<<"->"<<state<<" to state: Np_BR = "<<Np_BR[state][i]*Np_BR[i][j]*Np_BR[j][k]*Np_BR[k][l]<<" intermediate states: "<<k<<", "<<j<<", "<<i<<endl;
									#endif
									for(int m=l+1; m < NRQCDvars::nStates; m++){
										if(Np_BR[l][m]>0){
											ListFeedDownStateBranchingFractions[k]+=Np_BR[state][i]*Np_BR[i][j]*Np_BR[j][k]*Np_BR[k][l]*Np_BR[l][m];
											#ifdef mydebug
											cout<<"Branching "<<m<<"->"<<state<<" to state: Np_BR = "<<Np_BR[state][i]*Np_BR[i][j]*Np_BR[j][k]*Np_BR[k][l]*Np_BR[l][m]<<" intermediate states: "<<l<<", "<<k<<", "<<j<<", "<<i<<endl;
											#endif
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


	ListFeedDownStateBranchingFractions[state]=1; //The directly produced contribution has 'BR'=1

	// Calculate directProduction of necessary states
	dmatrix  BufferProductionVec;
	dmatrix  BufferPromptProductionMatrix;
	dvector BufferNullVec (4,0);

	for(int iFeedDownState=0; iFeedDownState<state; iFeedDownState++){
		BufferProductionVec.push_back(BufferNullVec);
		if(returnMPDetails) BufferPromptProductionMatrix.push_back(BufferNullVec);
	}

	for(int iFeedDownState=state; iFeedDownState < NRQCDvars::nStates; iFeedDownState++){
		if(ListFeedDownStateBranchingFractions[iFeedDownState]>0){
			dmatrix directProductionMatrix;
			dvector directProductionVector = getDirectProduction(iFeedDownState, Op, Np_BR, Np_US, returnMPDetails, directProductionMatrix);
			if(returnMPDetails) directProductionCube.push_back(directProductionMatrix);
			if(directProductionVector[0]>1e-100){
				BufferProductionVec.push_back(directProductionVector);
				if(returnMPDetails) BufferPromptProductionMatrix.push_back(directProductionVector);
			}
			else{
				if(returnMPDetails) BufferPromptProductionMatrix.push_back(BufferNullVec);
				return BufferNullVec;
			}
		}
		else{
			BufferProductionVec.push_back(BufferNullVec);
			if(returnMPDetails) BufferPromptProductionMatrix.push_back(BufferNullVec);
		}
	}


	double PromptCrossSect=0;

	// Calculate PromptCrossSect
	dvector PromptCrossSectCont (NRQCDvars::nStates, 0);
	for(int iFeedDownState=state; iFeedDownState < NRQCDvars::nStates; iFeedDownState++){
		if(ListFeedDownStateBranchingFractions[iFeedDownState]>0){
			PromptCrossSectCont[iFeedDownState]=BufferProductionVec[iFeedDownState][0]*ListFeedDownStateBranchingFractions[iFeedDownState];
			PromptCrossSect+=PromptCrossSectCont[iFeedDownState];
			if(returnMPDetails) BufferPromptProductionMatrix[iFeedDownState][0]=PromptCrossSectCont[iFeedDownState];
		}
		//else PromptCrossSectCont[iFeedDownState]=0;
	}

	//calculate fractions of individual feed-down contributions
	dvector PromptContFraction (NRQCDvars::nStates,0);
	for(int iFeedDownState=state; iFeedDownState < NRQCDvars::nStates; iFeedDownState++){
		PromptContFraction[iFeedDownState]=PromptCrossSectCont[iFeedDownState]/PromptCrossSect;
	}

#ifdef mydebug
	for(int iFeedDownState=0; iFeedDownState<NRQCDvars::nStates; iFeedDownState++){
		cout<<"BranchingRatio from state "<<iFeedDownState<<": "<<ListFeedDownStateBranchingFractions[iFeedDownState]<<endl;
	}
	for(int iFeedDownState=0; iFeedDownState<NRQCDvars::nStates; iFeedDownState++){
		cout<<"PromptCrossSectCont[ "<<iFeedDownState<<"]: "<<PromptCrossSectCont[iFeedDownState]<<endl;
	}
	for(int iFeedDownState=0; iFeedDownState<NRQCDvars::nStates; iFeedDownState++){
		cout<<"PromptContFraction[ "<<iFeedDownState<<"]: "<<PromptContFraction[iFeedDownState]<<endl;
	}
#endif


	//calculate the polarization of the prompt production
	dvector FeedDownContLamth_ (NRQCDvars::nStates);
	dvector FeedDownContLamph_ (NRQCDvars::nStates);
	dvector FeedDownContLamtp_ (NRQCDvars::nStates);

	for(int iFeedDownState=state; iFeedDownState < NRQCDvars::nStates; iFeedDownState++){
		if(PromptCrossSectCont[iFeedDownState]>0){
			FeedDownContLamth_[iFeedDownState]=BufferProductionVec[iFeedDownState][1];
			FeedDownContLamph_[iFeedDownState]=BufferProductionVec[iFeedDownState][2];
			FeedDownContLamtp_[iFeedDownState]=BufferProductionVec[iFeedDownState][3];
		}
	}

	double PromptLamth_numerator=0;
	double PromptLamth_denominator=0;
	for(int iFeedDownState=state; iFeedDownState < NRQCDvars::nStates; iFeedDownState++){
		if(PromptCrossSectCont[iFeedDownState]>0){
			PromptLamth_numerator+=PromptContFraction[iFeedDownState]*FeedDownContLamth_[iFeedDownState]/(3+FeedDownContLamth_[iFeedDownState]);
			PromptLamth_denominator+=PromptContFraction[iFeedDownState]/(3+FeedDownContLamth_[iFeedDownState]);
		}
	}

	double PromptLamph_numerator=0;
	for(int iFeedDownState=state; iFeedDownState < NRQCDvars::nStates; iFeedDownState++){
		if(PromptCrossSectCont[iFeedDownState]>0)
			PromptLamph_numerator+=PromptContFraction[iFeedDownState]*FeedDownContLamph_[iFeedDownState]/(3+FeedDownContLamth_[iFeedDownState]);
	}

	double PromptLamtp_numerator=0;
	for(int iFeedDownState=state; iFeedDownState < NRQCDvars::nStates; iFeedDownState++){
		if(PromptCrossSectCont[iFeedDownState]>0)
			PromptLamtp_numerator+=PromptContFraction[iFeedDownState]*FeedDownContLamtp_[iFeedDownState]/(3+FeedDownContLamth_[iFeedDownState]);
	}


	double PromptLamth=PromptLamth_numerator/PromptLamth_denominator;
	double PromptLamph=PromptLamph_numerator/PromptLamth_denominator;
	double PromptLamtp=PromptLamtp_numerator/PromptLamth_denominator;






	//set the return vector
	PromptProductionVec.push_back(PromptCrossSect);
	PromptProductionVec.push_back(PromptLamth);
	PromptProductionVec.push_back(PromptLamph);
	PromptProductionVec.push_back(PromptLamtp);

#ifdef mydebug
	cout<<"PromptCrossSect "<<PromptCrossSect<<endl;
	cout<<"PromptLamth "<<PromptLamth<<endl;
	cout<<"PromptLamph "<<PromptLamph<<endl;
	cout<<"PromptLamtp "<<PromptLamtp<<endl;
#endif

	if(returnMPDetails) promptProductionMatrix=BufferPromptProductionMatrix;

	return PromptProductionVec;
}





////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double NRQCDglobalfitObject::getCorrPromptCrossSect(dvector &PredPromptCrossSect, dmatrix &Np_BR, dmatrix &Np_US, bool returnMPDetails, double &polCorrFactor)  {
// Correction to the prompt cross section due to effects concerning the measurement (global Luminosity uncertainty, polarization dependence)
//	cout<<"NRQCDglobalfitObject::getCorrPromptCrossSect"<<endl;

	double CorrPromptCrossSect=PredPromptCrossSect[0];

	// Luminosity scaling

	double ErrGlobalScale = 1.+getErrGlobal()/getCentralValue()*Np_US[0][getExperiment()];
	CorrPromptCrossSect/=ErrGlobalScale;

#ifdef mydebug
	cout<<"getMeasurementID(): "<<getMeasurementID()<<endl;
	cout<<"getState(): "<<getState()<<endl;
	cout<<"getpTMin(): "<<getpTMin()<<endl;
	cout<<"getyMin(): "<<getyMin()<<endl;
	cout<<"getErrGlobal(): "<<getErrGlobal()<<endl;
	cout<<"getCentralValue(): "<<getCentralValue()<<endl;
	cout<<"getExperiment(): "<<getExperiment()<<endl;
	cout<<"gNp_US[0][getExperiment()]: "<<Np_US[0][getExperiment()]<<endl;
	cout<<"Luminosity correction factor: "<<ErrGlobalScale<<endl;
	cout<<"Luminosity corrected PromptCrossSect: "<<CorrPromptCrossSect<<endl;
#endif

	//Polarization correction
	dvector PolCorrParams = getPolCorrParams();
	double sigma_LongHX=PolCorrParams[0];
	double sigma_TransHX=PolCorrParams[1];
	double sigma_Unpolarized=getCentralValue();
	double PredPromptLamth=PredPromptCrossSect[1];

	//cout<<"PredPromptLamth "<<PredPromptLamth<<endl;

	double PolarizationCorrectionFactor=1+(sigma_TransHX-sigma_LongHX)/(2.*sigma_Unpolarized)*PredPromptLamth+
			                              ((sigma_TransHX+sigma_LongHX)/(2.*sigma_Unpolarized)-1)*PredPromptLamth*PredPromptLamth;

	CorrPromptCrossSect/=PolarizationCorrectionFactor;
	//PolarizationCorrectionFactor: sigma_unpolarized*PolarizationCorrectionFactor(lambda)=correct polarization


#ifdef mydebug
	cout<<"Polarization correction factor: "<<PolarizationCorrectionFactor<<endl;
	cout<<"Polarization corrected PromptCrossSect: "<<CorrPromptCrossSect<<endl;
#endif

	if(returnMPDetails) polCorrFactor=PolarizationCorrectionFactor;

	return CorrPromptCrossSect;

}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

vector<double> NRQCDglobalfitObject::getObjectLikelihood(dmatrix  &Op, dmatrix &Np_BR, dmatrix &Np_US, bool returnMPDetails, dcube &directProductionCube, dmatrix &promptProductionMatrix, double &polCorrFactor)  {
	//Returns per-object likelihood [0] and ModelPrediction [1]
//	cout<<"NRQCDglobalfitObject::getObjectLikelihood"<<endl;

    // dcube (nStates, nColorChannels, 4) directProductionCube
    // dmatrix (nStates, 4) promptProductionMatrix
    // double polCorrFactor

	vector<double> ObjectLikelihoodVec(2);

	dvector PromptProduction = getPromptProduction(getState(), Op, Np_BR, Np_US, returnMPDetails, directProductionCube, promptProductionMatrix);

	if(PromptProduction[0]<1e-100) {//cout<<"Penalty on likelihood (DirectProduction[0] of a contribution < 0)"<<endl;
		ObjectLikelihoodVec.at(0)=1e100;
		ObjectLikelihoodVec.at(1)=1e100;
		return ObjectLikelihoodVec;} // Penalty on likelihood if PromptProduction[0] = 0 -> true if any DirectProduction[0]<0


	double ModelPrediction;
	double Data_PointByPoint_Untertainty=getErrTot();
	double Model_PointByPoint_Untertainty=0;
	double Combined_DataModel_PointByPoint_Untertainty=TMath::Sqrt(Model_PointByPoint_Untertainty*Model_PointByPoint_Untertainty+
			                                                       Data_PointByPoint_Untertainty*Data_PointByPoint_Untertainty);

	dvector PromptProductionRatioDenom;
	dcube directProductionCubeDenom;
	//dmatrix directProductionMatrix_;
    //dvector model_directProductionVec_;
	//dmatrix directProductionMatrixDenom_;
    //dvector model_directProductionVecDenom_;


//	int MeasurementID=getMeasurementID();
	switch(Object_MeasurementID){ //modified by Joao
	case 0:
		ModelPrediction=getCorrPromptCrossSect(PromptProduction, Np_BR, Np_US, returnMPDetails, polCorrFactor);
		break;
	case 1:
		ModelPrediction=PromptProduction[Object_MeasurementID];
		break;
	case 2:
		ModelPrediction=PromptProduction[Object_MeasurementID];
		break;
	case 3:
		ModelPrediction=PromptProduction[Object_MeasurementID];
		break;
	case 4:
		PromptProductionRatioDenom = getPromptProduction(getStateDenom(), Op, Np_BR, Np_US, returnMPDetails, directProductionCubeDenom, promptProductionMatrix);
		ModelPrediction=getCorrPromptCrossSect(PromptProduction, Np_BR, Np_US, returnMPDetails, polCorrFactor)/getCorrPromptCrossSect(PromptProductionRatioDenom, Np_BR, Np_US, returnMPDetails, polCorrFactor);



		//directProductionMatrix_=directProductionCube.at(0);
	    //model_directProductionVec_=directProductionMatrix_.at(0);
		//directProductionMatrixDenom_=directProductionCubeDenom.at(0);
	    //model_directProductionVecDenom_=directProductionMatrixDenom_.at(0);

		// That would be the solution:
		//directProductionCube[0][0][0]=directProductionCube[0][0][0]/directProductionCubeDenom[0][0][0];
		//similar for promptProductionMatrix... make a different object for the denominator and divide


	    break;
	case 5:
		PromptProductionRatioDenom = getPromptProduction(getStateDenom(), Op, Np_BR, Np_US, returnMPDetails, directProductionCubeDenom, promptProductionMatrix);//Prompt Jpsi
		ModelPrediction=getCorrPromptCrossSect(PromptProduction, Np_BR, Np_US, returnMPDetails, polCorrFactor)/getCorrPromptCrossSect(PromptProductionRatioDenom, Np_BR, Np_US, returnMPDetails, polCorrFactor);
		break;
	default:
		cerr << "Error: Unknown MeasurementID. Execution stop!" << endl;
		exit(1);
	}


    double ObjectLikelihood;
    double BufferLikelihood=(ModelPrediction-Object_CentralValue)/Combined_DataModel_PointByPoint_Untertainty;
    ObjectLikelihood=BufferLikelihood*BufferLikelihood;

#ifdef mydebug
	cout << "ModelPrediction: " << ModelPrediction << endl;
	cout << "Object_CentralValue: " << Object_CentralValue << endl;
	cout << "Combined_DataModel_PointByPoint_Untertainty: " << Combined_DataModel_PointByPoint_Untertainty << endl;
	cout << "ObjectLikelihood: " << ObjectLikelihood << endl;
#endif

	ObjectLikelihoodVec.at(0)=ObjectLikelihood;
	ObjectLikelihoodVec.at(1)=ModelPrediction;

    return ObjectLikelihoodVec;

}


void NRQCDglobalfitObject::Dump(const int nStates, bool DumpData, bool DumpModel) {


	if(DumpData){
	  cout << endl;
	  cout << "DATA members of NRQCDglobalfitObject:::"    << endl;
	  cout << "ObjectID()       : " << getObjectID()       << endl;
	  cout << "State()          : " << getState()          << endl;
	  cout << "StateRatioDenom(): " << getStateDenom()     << endl;
	  cout << "MeasurementID()  : " << getMeasurementID()  << endl;
	  cout << "Experiment()     : " << getExperiment()     << endl;
	  cout << "CentralValue()   : " << getCentralValue()   << endl;
	  cout << "ErrStatPos()     : " << getErrStatPos()     << endl;
	  cout << "ErrStatNeg()     : " << getErrStatNeg()     << endl;
	  cout << "ErrStat()        : " << getErrStat()        << endl;
	  cout << "ErrSystPos()     : " << getErrSystPos()     << endl;
	  cout << "ErrSystNeg()     : " << getErrSystNeg()     << endl;
	  cout << "ErrSyst()        : " << getErrSyst()        << endl;
	  cout << "ErrTotPos()      : " << getErrTotPos()      << endl;
	  cout << "ErrTotNeg()      : " << getErrTotNeg()      << endl;
	  cout << "ErrTot()         : " << getErrTot()         << endl;
	  cout << "ErrGlobalPos()   : " << getErrGlobalPos()   << endl;
	  cout << "ErrGlobalNeg()   : " << getErrGlobalNeg()   << endl;
	  cout << "ErrGlobal()      : " << getErrGlobal()      << endl;
	  cout << "yMin()           : " << getyMin()           << endl;
	  cout << "yMax()           : " << getyMax()           << endl;
	  cout << "yMean()          : " << getyMean()          << endl;
	  cout << "pTMin()          : " << getpTMin()          << endl;
	  cout << "pTMax()          : " << getpTMax()          << endl;
	  cout << "pTMean()         : " << getpTMean()         << endl;
	  dvector BufferPolCorrParams=getPolCorrParams();
	  cout<<"PolCorrParams = ";
	  for(int j=0;j<BufferPolCorrParams.size();j++){
		  cout<<"   "<<BufferPolCorrParams[j];
	  }
	  cout<<endl;
	  if(getisDataValid())
	  cout<< "isDataValid()     : DATA is valid"           << endl;
	  else
	  cout<< "isDataValid()     : DATA is NOT valid"           << endl;
	  if(getisAbsRap())
	  cout<< "sAbsRap()     : rapidity region to be interpreted in absolute rapidity"           << endl;
	  else
	  cout<< "sAbsRap()     : rapidity region NOT to be interpreted in absolute rapidity"           << endl;


	}

	if(DumpModel){
	  cout << endl;
	  cout << "MODEL members of NRQCDglobalfitObject:::" << endl;
	  for(int i=0;i<nStates;i++){
		  int nCC=getShortDistanceCoef(getState()).size();
		  cout<< "nMotherState = "<<i<<endl;
		  dvector BufferShortDistanceCoef=getShortDistanceCoef(i);
		  dvector BufferOctetCompLamth=getOctetCompLamth(i);
		  dvector BufferOctetCompLamph=getOctetCompLamph(i);
		  dvector BufferOctetCompLamtp=getOctetCompLamtp(i);

		  dvector BufferShortDistanceCoef_PointByPointSyst=getShortDistanceCoef_PointByPointSyst(i);
		  dvector BufferOctetCompLamth_PointByPointSyst=getOctetCompLamth_PointByPointSyst(i);
		  dvector BufferOctetCompLamph_PointByPointSyst=getOctetCompLamph_PointByPointSyst(i);
		  dvector BufferOctetCompLamtp_PointByPointSyst=getOctetCompLamtp_PointByPointSyst(i);

		  dmatrix  BufferShortDistanceCoef_globalSystPos=getShortDistanceCoef_globalSystPos(i);
		  dmatrix  BufferShortDistanceCoef_globalSystNeg=getShortDistanceCoef_globalSystNeg(i);
		  dmatrix  BufferOctetCompLamth_globalSystPos=getOctetCompLamth_globalSystPos(i);
		  dmatrix  BufferOctetCompLamth_globalSystNeg=getOctetCompLamth_globalSystNeg(i);
		  dmatrix  BufferOctetCompLamph_globalSystPos=getOctetCompLamph_globalSystPos(i);
		  dmatrix  BufferOctetCompLamph_globalSystNeg=getOctetCompLamph_globalSystNeg(i);
		  dmatrix  BufferOctetCompLamtp_globalSystPos=getOctetCompLamtp_globalSystPos(i);
		  dmatrix  BufferOctetCompLamtp_globalSystNeg=getOctetCompLamtp_globalSystNeg(i);

		  cout<<"ShortDistanceCoef = "; for(int j=0;j<nCC;j++)cout<<"   "<<BufferShortDistanceCoef[j]; cout<<endl;
		  cout<<"OctetCompLamth = "; for(int j=0;j<nCC;j++) cout<<"   "<<BufferOctetCompLamth[j]; cout<<endl;
		  cout<<"OctetCompLamph = "; for(int j=0;j<nCC;j++) cout<<"   "<<BufferOctetCompLamph[j]; cout<<endl;
		  cout<<"OctetCompLamtp = "; for(int j=0;j<nCC;j++) cout<<"   "<<BufferOctetCompLamtp[j]; cout<<endl;

		  cout<<"ShortDistanceCoef_PointByPointSyst = "; for(int j=0;j<nCC;j++)cout<<"   "<<BufferShortDistanceCoef_PointByPointSyst[j]; cout<<endl;
		  cout<<"OctetCompLamth_PointByPointSyst = "; for(int j=0;j<nCC;j++) cout<<"   "<<BufferOctetCompLamth_PointByPointSyst[j]; cout<<endl;
		  cout<<"OctetCompLamph_PointByPointSyst = "; for(int j=0;j<nCC;j++) cout<<"   "<<BufferOctetCompLamph_PointByPointSyst[j]; cout<<endl;
		  cout<<"OctetCompLamtp_PointByPointSyst = "; for(int j=0;j<nCC;j++) cout<<"   "<<BufferOctetCompLamtp_PointByPointSyst[j]; cout<<endl;

		  int nMSS=BufferShortDistanceCoef_globalSystPos[0].size();
		  for(int i=0;i<nMSS;i++){
		  cout<<"ShortDistanceCoef_globalSystPos(nMSS="<<i<<") = "; for(int j=0;j<nCC;j++)cout<<"   "<<BufferShortDistanceCoef_globalSystPos[j][i]; cout<<endl;
		  cout<<"OctetCompLamth_globalSystPos(nMSS="<<i<<") = "; for(int j=0;j<nCC;j++) cout<<"   "<<BufferOctetCompLamth_globalSystPos[j][i]; cout<<endl;
		  cout<<"OctetCompLamph_globalSystPos(nMSS="<<i<<") = "; for(int j=0;j<nCC;j++) cout<<"   "<<BufferOctetCompLamph_globalSystPos[j][i]; cout<<endl;
		  cout<<"OctetCompLamtp_globalSystPos(nMSS="<<i<<") = "; for(int j=0;j<nCC;j++) cout<<"   "<<BufferOctetCompLamtp_globalSystPos[j][i]; cout<<endl;
		  }



		  cout<<endl;

	  }
	  if(getisModelValid())
		  cout<<"isModelValid()   : MODEL is valid"           << endl;
	  else
		  cout<<"isModelValid()   : MODEL is NOT valid"           << endl;

	}
	  cout << endl;

}


void NRQCDglobalfitObject::setInvalidModel(int Input_nStates, int Input_nColorChannels, int Input_nModelSystematicScales) {

	double Dummy_SDC=999.;
	double Dummy_Lam=0.0999;

	//Dummy_SDC=1.;
	//Dummy_Lam=0.;

	for(int iMotherState=0;iMotherState<Input_nStates;iMotherState++){

		dvector setModel_ShortDistanceCoef;
		dvector setModel_OctetCompLamth;
		dvector setModel_OctetCompLamph;
		dvector setModel_OctetCompLamtp;
		dvector setModel_ShortDistanceCoef_PointByPointSyst;
		dvector setModel_OctetCompLamth_PointByPointSyst;
		dvector setModel_OctetCompLamph_PointByPointSyst;
		dvector setModel_OctetCompLamtp_PointByPointSyst;
		dmatrix setModel_ShortDistanceCoef_globalSystPos;
		dmatrix setModel_ShortDistanceCoef_globalSystNeg;
		dmatrix setModel_OctetCompLamth_globalSystPos;
		dmatrix setModel_OctetCompLamth_globalSystNeg;
		dmatrix setModel_OctetCompLamph_globalSystPos;
		dmatrix setModel_OctetCompLamph_globalSystNeg;
		dmatrix setModel_OctetCompLamtp_globalSystPos;
		dmatrix setModel_OctetCompLamtp_globalSystNeg;


		for (int iColorChannel=0; iColorChannel<Input_nColorChannels; iColorChannel++){

			setModel_ShortDistanceCoef.push_back(Dummy_SDC);
			setModel_OctetCompLamth.push_back(Dummy_Lam);
			setModel_OctetCompLamph.push_back(Dummy_Lam);
			setModel_OctetCompLamtp.push_back(Dummy_Lam);
			setModel_ShortDistanceCoef_PointByPointSyst.push_back(Dummy_SDC);
			setModel_OctetCompLamth_PointByPointSyst.push_back(Dummy_Lam);
			setModel_OctetCompLamph_PointByPointSyst.push_back(Dummy_Lam);
			setModel_OctetCompLamtp_PointByPointSyst.push_back(Dummy_Lam);

			dvector setModel_Buff_ShortDistanceCoef_globalSystPos;
			dvector setModel_Buff_ShortDistanceCoef_globalSystNeg;
			dvector setModel_Buff_OctetCompLamth_globalSystPos;
			dvector setModel_Buff_OctetCompLamth_globalSystNeg;
			dvector setModel_Buff_OctetCompLamph_globalSystPos;
			dvector setModel_Buff_OctetCompLamph_globalSystNeg;
			dvector setModel_Buff_OctetCompLamtp_globalSystPos;
			dvector setModel_Buff_OctetCompLamtp_globalSystNeg;

			for (int iModelSystematicScale=0; iModelSystematicScale<Input_nModelSystematicScales; iModelSystematicScale++){

				setModel_Buff_ShortDistanceCoef_globalSystPos.push_back(Dummy_SDC);
				setModel_Buff_OctetCompLamth_globalSystPos.push_back(Dummy_Lam);
				setModel_Buff_OctetCompLamph_globalSystPos.push_back(Dummy_Lam);
				setModel_Buff_OctetCompLamtp_globalSystPos.push_back(Dummy_Lam);
				setModel_Buff_ShortDistanceCoef_globalSystNeg.push_back(Dummy_SDC);
				setModel_Buff_OctetCompLamth_globalSystNeg.push_back(Dummy_Lam);
				setModel_Buff_OctetCompLamph_globalSystNeg.push_back(Dummy_Lam);
				setModel_Buff_OctetCompLamtp_globalSystNeg.push_back(Dummy_Lam);

			}
			if(Input_nModelSystematicScales==0){//fill with dummy values, because >> << operators can not cope with non-existing dcubes
				for (int iModelSystematicScale=0; iModelSystematicScale<2; iModelSystematicScale++){
					setModel_Buff_ShortDistanceCoef_globalSystPos.push_back(Dummy_SDC);
					setModel_Buff_OctetCompLamth_globalSystPos.push_back(Dummy_Lam);
					setModel_Buff_OctetCompLamph_globalSystPos.push_back(Dummy_Lam);
					setModel_Buff_OctetCompLamtp_globalSystPos.push_back(Dummy_Lam);
					setModel_Buff_ShortDistanceCoef_globalSystNeg.push_back(Dummy_SDC);
					setModel_Buff_OctetCompLamth_globalSystNeg.push_back(Dummy_Lam);
					setModel_Buff_OctetCompLamph_globalSystNeg.push_back(Dummy_Lam);
					setModel_Buff_OctetCompLamtp_globalSystNeg.push_back(Dummy_Lam);
				}
			}

			setModel_ShortDistanceCoef_globalSystPos.push_back(setModel_Buff_ShortDistanceCoef_globalSystPos);
			setModel_OctetCompLamth_globalSystPos.push_back(setModel_Buff_OctetCompLamth_globalSystPos);
			setModel_OctetCompLamph_globalSystPos.push_back(setModel_Buff_OctetCompLamph_globalSystPos);
			setModel_OctetCompLamtp_globalSystPos.push_back(setModel_Buff_OctetCompLamtp_globalSystPos);
			setModel_ShortDistanceCoef_globalSystNeg.push_back(setModel_Buff_ShortDistanceCoef_globalSystNeg);
			setModel_OctetCompLamth_globalSystNeg.push_back(setModel_Buff_OctetCompLamth_globalSystNeg);
			setModel_OctetCompLamph_globalSystNeg.push_back(setModel_Buff_OctetCompLamph_globalSystNeg);
			setModel_OctetCompLamtp_globalSystNeg.push_back(setModel_Buff_OctetCompLamtp_globalSystNeg);

		}


		setModel(iMotherState,
				  setModel_ShortDistanceCoef, setModel_OctetCompLamth, setModel_OctetCompLamph, setModel_OctetCompLamtp,
				  setModel_ShortDistanceCoef_PointByPointSyst, setModel_OctetCompLamth_PointByPointSyst, setModel_OctetCompLamph_PointByPointSyst, setModel_OctetCompLamtp_PointByPointSyst,
				  setModel_ShortDistanceCoef_globalSystPos, setModel_ShortDistanceCoef_globalSystNeg,
				  setModel_OctetCompLamth_globalSystPos, setModel_OctetCompLamth_globalSystNeg,
				  setModel_OctetCompLamph_globalSystPos, setModel_OctetCompLamph_globalSystNeg,
				  setModel_OctetCompLamtp_globalSystPos, setModel_OctetCompLamtp_globalSystNeg,
				  true);

	}


}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                //
//                                             Sets Operator matrix elements basis                                                                //
//                               First vector f[0] is always the vector perpendicular to he hyperplane                                            //
//                                                                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//void NRQCDglobalfitObject::setFractionsBasis(const int& numOperatorStates){
//	f = newCMatrix(numOperatorStates, numOperatorStates);
//	double component;
//	ifstream fFile("./BasisVectors.dat", ifstream::in);
//	if( !fFile.is_open() ){
//		cerr << "Error: Basis vector files for operators not opened! Execution stop." << endl;
//		exit(1);
//	}
//	char* searchString[32];
//	sprintf(searchString, "Dimension%d", numOperatorStates);
//	string line;
//	do{                      // search for the correct block in the basis operator file
//		getline(fFile, line);
//		size_t found = line.find(searchString);
//	}
//	while( found == string::npos );
//	if( fFile.eof() ){
//		cerr << "Error: Basis vector block not found in operator basis file. Execution stop." << endl;
//		exit(1);
//	}
//	for(int i=0; i < numOperatorStates; ++i){ // now read the basis vectors onto the matrix (each row is a basis vector)
//		getline(fFile, line);
//		istringstream stream(line);
//		int j=0;
//	    while ( stream >> number, !stream.fail() ){
//	    	f[i][j]=number;
//	    	++j;
//	    }
//	    if( j >= numOperatorStates ){
//			cerr << "Error: Basis vector with wrong size in operator basis file. Execution stop." << endl;
//			exit(1);
//	    }
//	}
//}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                //
//                                              Copies a matrix                                                                                   //
//                                                                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void NRQCDglobalfitObject::dmatrix_copy(dmatrix& original, dmatrix& copy){
	for(dmatrix::iterator iti = original.begin(); iti != original.end(); ++iti){
	  int i=0;
	  for(dvector::iterator itj = iti->begin(); itj != iti->end(); ++itj){
		  int j=0;
		  copy[i][j] = *itj;
		  ++j;
	  }
	  ++i;
	}
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                //
//                                              Allocates memory to a matrix                                                                      //
//                                                                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

dmatrix* NRQCDglobalfitObject::dmatrix_allocate(const int nrow, const int ncolumn, const int value){
	dmatrix* alloc_matrix = new dmatrix;
	for (int i = 0; i < nrow; i++) {
	    dvector row; // Create an empty row
	    for (int j = 0; j < ncolumn; j++) {
	        row.push_back(value); // Add an element (column) to the row
	    }
	    alloc_matrix->push_back(row); // Add the row to the main vector
	}
	return alloc_matrix;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                //
//                                             Adds rows to a dmatrix                                                                             //
//                                                                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void NRQCDglobalfitObject::addRows(dmatrix& matrix, const int nrows, const int value){
	dvector newRow(matrix[0]);
	for(dvector::iterator it = newRow.begin(); it != newRow.end(); ++it){ //initialize all elements of the new row to value
		*it = value;
	}
	for(int i=0; i < nrows; ++i){
		matrix.push_back(newRow);
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                //
//                                              Adds columns to a dmatrix                                                                         //
//                                                                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void NRQCDglobalfitObject::addColumns(dmatrix& matrix, const int ncols, int value){
	for(dmatrix::iterator it = matrix.begin(); it != matrix.end(); ++it){ //initialize all elements of the new column to value
		for (int j = 0; j < ncols; ++j) {
			it->push_back(value);
		}
	}
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////                          I/O part                      //////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                //
//                                              Overloaded << ostream operator dvector                                                            //
//                                                                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

ostream& operator<< (ostream &out, dvector& my_vector){
	for (dvector::iterator it = my_vector.begin(); it != my_vector.end(); ++it) {
		out << *it << " ";
	}
	out << endl;
	return out;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                //
//                                              Overloaded << ostream operator dmatrix                                                            //
//                                                                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

ostream& operator<< (ostream &out, dmatrix& my_matrix){
	for (dmatrix::iterator iti = my_matrix.begin(); iti != my_matrix.end(); ++iti) {
		out << *iti;
	}
	out << endl;
	return out;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                //
//                                              Overloaded << ostream operator dcube                                                            //
//                                                                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

ostream& operator<< (ostream &out, dcube& my_cube){
	for (dcube::iterator iti = my_cube.begin(); iti != my_cube.end(); ++iti) {
		out << *iti;
	}
	out << endl;
	return out;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                //
//                                              Overloaded >> istream operator dvector                                                            //
//                                                                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

istream& operator>> (istream &in, dvector& my_vector){
    string line;
    double number;

    if( my_vector.size() != 0 ){
		cerr << "WARNING: dvector not empty: data will be appended!" << endl;
	}
    getline(in, line);
    if( line.empty() ) {
    	in.setstate(ios::failbit); // no data to read sets failbit!
    	return in;
    }
    istringstream stream(line);
    while ( stream >> number, !stream.fail() )
    	my_vector.push_back(number);
    in.clear();
	return in;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                //
//                                              Overloaded >> istream operator dmatrix                                                            //
//                                                                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

istream& operator>> (istream &in, dmatrix& my_matrix){
	dvector row;

	if( my_matrix.size() != 0 ){
		cerr << "READ WARNING: dmatrix not empty: data will be appended!" << endl;
	}
	if( in >> row , !in.fail() ){
		my_matrix.push_back(row);
		row.erase(row.begin(), row.end());
	}
	else{	// no data to read sets failbit (dvector version) and returns the stream!
		return in;
	}
	while( in >> row , !in.fail() ){
		my_matrix.push_back(row);
		row.erase(row.begin(), row.end());
    }
	in.clear();
	return in;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                //
//                                              Overloaded >> istream operator dcube                                                              //
//                                                                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

istream& operator>> (istream &in, dcube& my_cube){
	dmatrix matrix;

	if( my_cube.size() != 0 ){
		cerr << "READ WARNING: dcube not empty: data will be appended!" << endl;
	}
	if( in >> matrix , !in.fail() ){
		my_cube.push_back(matrix);
		for(dmatrix::iterator it=matrix.begin(); it != matrix.end(); ++it){
			it->erase(it->begin(),it->end());
		}
		matrix.erase(matrix.begin(), matrix.end());
	}
	else{    // no data to read sets failbit (dmatrix version) and returns the stream!
		return in;
	}
	while( in >> matrix , !in.fail() ){
		my_cube.push_back(matrix);
		for(dmatrix::iterator it=matrix.begin(); it != matrix.end(); ++it){
			it->erase(it->begin(),it->end());
		}
		matrix.erase(matrix.begin(), matrix.end());
    }
	in.clear();
	return in;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                //
//                                              Overloaded << ostream operator fvector                                                            //
//                                                                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

ostream& operator<< (ostream &out, fvector& my_vector){
	for (fvector::iterator it = my_vector.begin(); it != my_vector.end(); ++it) {
		out << *it << " ";
	}
	out << endl;
	return out;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                //
//                                              Overloaded << ostream operator fmatrix                                                            //
//                                                                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

ostream& operator<< (ostream &out, fmatrix& my_matrix){
	for (fmatrix::iterator iti = my_matrix.begin(); iti != my_matrix.end(); ++iti) {
		out << *iti;
	}
	out << endl;
	return out;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                //
//                                              Overloaded << ostream operator fcube                                                            //
//                                                                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

ostream& operator<< (ostream &out, fcube& my_cube){
	for (fcube::iterator iti = my_cube.begin(); iti != my_cube.end(); ++iti) {
		out << *iti;
	}
	out << endl;
	return out;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                //
//                                              Overloaded >> istream operator fvector                                                            //
//                                                                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

istream& operator>> (istream &in, fvector& my_vector){
    string line;
    float number;

    if( my_vector.size() != 0 ){
		cerr << "WARNING: fvector not empty: data will be appended!" << endl;
	}
    getline(in, line);
    if( line.empty() ) {
    	in.setstate(ios::failbit); // no data to read sets failbit!
    	return in;
    }
    istringstream stream(line);
    while ( stream >> number, !stream.fail() )
    	my_vector.push_back(number);
    in.clear();
	return in;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                //
//                                              Overloaded >> istream operator fmatrix                                                            //
//                                                                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

istream& operator>> (istream &in, fmatrix& my_matrix){
	fvector row;

	if( my_matrix.size() != 0 ){
		cerr << "READ WARNING: fmatrix not empty: data will be appended!" << endl;
	}
	if( in >> row , !in.fail() ){
		my_matrix.push_back(row);
		row.erase(row.begin(), row.end());
	}
	else{	// no data to read sets failbit (fvector version) and returns the stream!
		return in;
	}
	while( in >> row , !in.fail() ){
		my_matrix.push_back(row);
		row.erase(row.begin(), row.end());
    }
	in.clear();
	return in;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                //
//                                              Overloaded >> istream operator fcube                                                              //
//                                                                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

istream& operator>> (istream &in, fcube& my_cube){
	fmatrix matrix;

	if( my_cube.size() != 0 ){
		cerr << "READ WARNING: fcube not empty: data will be appended!" << endl;
	}
	if( in >> matrix , !in.fail() ){
		my_cube.push_back(matrix);
		for(fmatrix::iterator it=matrix.begin(); it != matrix.end(); ++it){
			it->erase(it->begin(),it->end());
		}
		matrix.erase(matrix.begin(), matrix.end());
	}
	else{    // no data to read sets failbit (fmatrix version) and returns the stream!
		return in;
	}
	while( in >> matrix , !in.fail() ){
		my_cube.push_back(matrix);
		for(fmatrix::iterator it=matrix.begin(); it != matrix.end(); ++it){
			it->erase(it->begin(),it->end());
		}
		matrix.erase(matrix.begin(), matrix.end());
    }
	in.clear();
	return in;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                //
//                                              Overloaded << ostream operator ivector                                                            //
//                                                                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

ostream& operator<< (ostream &out, ivector& my_vector){
	for (ivector::iterator it = my_vector.begin(); it != my_vector.end(); ++it) {
		out << *it << " ";
	}
	out << endl;
	return out;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                //
//                                              Overloaded << ostream operator imatrix                                                            //
//                                                                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

ostream& operator<< (ostream &out, imatrix& my_matrix){
	for (imatrix::iterator iti = my_matrix.begin(); iti != my_matrix.end(); ++iti) {
		out << *iti;
	}
	out << endl;
	return out;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                //
//                                              Overloaded << ostream operator icube                                                            //
//                                                                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

ostream& operator<< (ostream &out, icube& my_cube){
	for (icube::iterator iti = my_cube.begin(); iti != my_cube.end(); ++iti) {
		out << *iti;
	}
	out << endl;
	return out;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                //
//                                              Overloaded >> istream operator ivector                                                            //
//                                                                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

istream& operator>> (istream &in, ivector& my_vector){
    string line;
    int number;

    if( my_vector.size() != 0 ){
		cerr << "WARNING: ivector not empty: data will be appended!" << endl;
	}
    getline(in, line);
    if( line.empty() ) {
    	in.setstate(ios::failbit); // no data to read sets failbit!
    	return in;
    }
    istringstream stream(line);
    while ( stream >> number, !stream.fail() )
    	my_vector.push_back(number);
    in.clear();
	return in;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                //
//                                              Overloaded >> istream operator imatrix                                                            //
//                                                                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

istream& operator>> (istream &in, imatrix& my_matrix){
	ivector row;

	if( my_matrix.size() != 0 ){
		cerr << "READ WARNING: imatrix not empty: data will be appended!" << endl;
	}
	if( in >> row , !in.fail() ){
		my_matrix.push_back(row);
		row.erase(row.begin(), row.end());
	}
	else{	// no data to read sets failbit (ivector version) and returns the stream!
		return in;
	}
	while( in >> row , !in.fail() ){
		my_matrix.push_back(row);
		row.erase(row.begin(), row.end());
    }
	in.clear();
	return in;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                //
//                                              Overloaded >> istream operator fcube                                                              //
//                                                                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

istream& operator>> (istream &in, icube& my_cube){
	imatrix matrix;

	if( my_cube.size() != 0 ){
		cerr << "READ WARNING: dcube not empty: data will be appended!" << endl;
	}
	if( in >> matrix , !in.fail() ){
		my_cube.push_back(matrix);
		for(imatrix::iterator it=matrix.begin(); it != matrix.end(); ++it){
			it->erase(it->begin(),it->end());
		}
		matrix.erase(matrix.begin(), matrix.end());
	}
	else{    // no data to read sets failbit (imatrix version) and returns the stream!
		return in;
	}
	while( in >> matrix , !in.fail() ){
		my_cube.push_back(matrix);
		for(imatrix::iterator it=matrix.begin(); it != matrix.end(); ++it){
			it->erase(it->begin(),it->end());
		}
		matrix.erase(matrix.begin(), matrix.end());
    }
	in.clear();
	return in;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                //
//                                              Overloaded << ostream operator for NRQCDglobalfitObject                                           //
//                                                                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

ostream& operator<< (ostream &out, NRQCDglobalfitObject& fit_object)
{

	out << fit_object.Object_ShortDistanceCoef_globalSystPos; // [nMothers][nColorChannels][nModelSystematicScales]
	out << fit_object.Object_ShortDistanceCoef_globalSystNeg; // [nMothers][nColorChannels][nModelSystematicScales]
	out << fit_object.Object_OctetCompLamth_globalSystPos; // [nMothers][nColorChannels][nModelSystematicScales]
	out << fit_object.Object_OctetCompLamth_globalSystNeg; // [nMothers][nColorChannels][nModelSystematicScales]
	out << fit_object.Object_OctetCompLamph_globalSystPos; // [nMothers][nColorChannels][nModelSystematicScales]
	out << fit_object.Object_OctetCompLamph_globalSystNeg; // [nMothers][nColorChannels][nModelSystematicScales]
	out << fit_object.Object_OctetCompLamtp_globalSystPos; // [nMothers][nColorChannels][nModelSystematicScales]
	out << fit_object.Object_OctetCompLamtp_globalSystNeg; // [nMothers][nColorChannels][nModelSystematicScales]
	out << fit_object.Object_ShortDistanceCoef; // dvector ShortDistanceCoef are the Color Octet Components of one state
	out << fit_object.Object_OctetCompLamth; // empty line signals the end of the matrix  (endl)
	out << fit_object.Object_OctetCompLamph; // empty line signals the end of the matrix  (endl)
	out << fit_object.Object_OctetCompLamtp; // empty line signals the end of the matrix  (endl)
	out << fit_object.Object_ShortDistanceCoef_PointByPointSyst; // [nMothers][nModelSystematicScales]
	out << fit_object.Object_OctetCompLamth_PointByPointSyst; // [nMothers][nModelSystematicScales]
	out << fit_object.Object_OctetCompLamph_PointByPointSyst; // [nMothers][nModelSystematicScales]
	out << fit_object.Object_OctetCompLamtp_PointByPointSyst; // [nMothers][nModelSystematicScales]
	out << fit_object.Object_PolCorrParams << endl;
	out << fit_object.Object_nState << endl; // 1...Jpsi, 2...chic1, 3...chic2, 4...Psi', 5...Y1S, 6...chib1P1, 7...chib1P2, 8...Y2S, 9...chib2P1, 10...chib2P2, 11...Y3S, 12...chib3P1, 13...chib3P2
	out << fit_object.Object_nStateDenom << endl;//State in the denominator of a cross-section ratio or feed-down fraction
	out << fit_object.Object_MeasurementID << endl; //0...Cross section, 1...LamthHX, 2...LamphHX, 3...LamtpHX, 4...Cross section ratio, 5...Feed-down fraction
	out << fit_object.Object_nExperiment << endl;
	out << fit_object.Object_CentralValue << endl;
	out << fit_object.Object_ErrStatPos << " " << fit_object.Object_ErrStatNeg << endl;
	out << fit_object.Object_ErrSystPos << " " << fit_object.Object_ErrSystNeg << endl;
	out << fit_object.Object_ErrTotPos << " " << fit_object.Object_ErrTotNeg  << endl;
	out << fit_object.Object_ErrGlobalPos << " " << fit_object.Object_ErrGlobalNeg << endl;
	out << fit_object.Object_yMin << " " << fit_object.Object_yMax << " " << fit_object.Object_yMean << endl;
	out << fit_object.Object_pTMin << " " << fit_object.Object_pTMax << " " << fit_object.Object_pTMean << endl;
	out << fit_object.Object_ObjectID << endl; // Experiment=XXX Energy=XXXTeV Observable=XXX Unit=XXX State=XXX
	out << fit_object.Object_isDataValid << endl;
	out << fit_object.Object_isModelValid << endl;
	out << fit_object.Object_isAbsRap;

    return out;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                //
//                                              Overloaded >> istream operator for NRQCDglobalfitObject                                           //
//                                                                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

istream& operator>> (istream &in, NRQCDglobalfitObject& fit_object)
{
	in >> fit_object.Object_ShortDistanceCoef_globalSystPos; // [nMothers][nColorChannels][nModelSystematicScales]
	in >> fit_object.Object_ShortDistanceCoef_globalSystNeg; // [nMothers][nColorChannels][nModelSystematicScales]
	in >> fit_object.Object_OctetCompLamth_globalSystPos; // [nMothers][nColorChannels][nModelSystematicScales]
	in >> fit_object.Object_OctetCompLamth_globalSystNeg; // [nMothers][nColorChannels][nModelSystematicScales]
	in >> fit_object.Object_OctetCompLamph_globalSystPos; // [nMothers][nColorChannels][nModelSystematicScales]
	in >> fit_object.Object_OctetCompLamph_globalSystNeg; // [nMothers][nColorChannels][nModelSystematicScales]
	in >> fit_object.Object_OctetCompLamtp_globalSystPos; // [nMothers][nColorChannels][nModelSystematicScales]
	in >> fit_object.Object_OctetCompLamtp_globalSystNeg; // [nMothers][nColorChannels][nModelSystematicScales]
	in >> fit_object.Object_ShortDistanceCoef; // dvector ShortDistanceCoef are the Color Octet Components of one nState
	in >> fit_object.Object_OctetCompLamth; // empty line signals the end of the matrix  (endl)
	in >> fit_object.Object_OctetCompLamph; // empty line signals the end of the matrix  (endl)
	in >> fit_object.Object_OctetCompLamtp; // empty line signals the end of the matrix  (endl)
	in >> fit_object.Object_ShortDistanceCoef_PointByPointSyst; // [nMothers][nModelSystematicScales]
	in >> fit_object.Object_OctetCompLamth_PointByPointSyst; // [nMothers][nModelSystematicScales]
	in >> fit_object.Object_OctetCompLamph_PointByPointSyst; // [nMothers][nModelSystematicScales]
	in >> fit_object.Object_OctetCompLamtp_PointByPointSyst; // [nMothers][nModelSystematicScales]
	in >> fit_object.Object_PolCorrParams;
	in >> fit_object.Object_nState; // 1...Jpsi, 2...chic1, 3...chic2, 4...Psi', 5...Y1S, 6...chib1P1, 7...chib1P2, 8...Y2S, 9...chib2P1, 10...chib2P2, 11...Y3S, 12...chib3P1, 13...chib3P2
	in >> fit_object.Object_nStateDenom;//State in the denominator of a cross-section ratio or feed-down fraction
	in >> fit_object.Object_MeasurementID; //0...Cross section, 1...LamthHX, 2...LamphHX, 3...LamtpHX, 4...Cross section ratio, 5...Feed-down fraction
	in >> fit_object.Object_nExperiment;
	in >> fit_object.Object_CentralValue;
	in >> fit_object.Object_ErrStatPos >>fit_object.Object_ErrStatNeg;
	in >> fit_object.Object_ErrSystPos >>fit_object.Object_ErrSystNeg;
	in >> fit_object.Object_ErrTotPos >>fit_object.Object_ErrTotNeg ;
	in >> fit_object.Object_ErrGlobalPos >>fit_object.Object_ErrGlobalNeg;
	in >> fit_object.Object_yMin >>fit_object.Object_yMax >>fit_object.Object_yMean;
	in >> fit_object.Object_pTMin >>fit_object.Object_pTMax >>fit_object.Object_pTMean;
	in >> fit_object.Object_ObjectID; // Experiment=XXX Energy=XXXTeV Observable=XXX Unit=XXX State=XXX
	in >> fit_object.Object_isDataValid;
	in >> fit_object.Object_isModelValid;
	in >> fit_object.Object_isAbsRap;

    return in;
}

//ClassImp(NRQCDglobalfitObject)

// Creates a vector with nEls initialized to 0

template <class T>
T* newCVector(const int& nEls){
	T *vector;
	vector = new T[nEls];
	if( !vector ){
		cerr << "Error: Unable to allocate buffer space!" << endl;
		exit(1);
	}
	memset(vector, 0, nEls*sizeof(T));
	return vector;
}

// Creates a Float_t matrix with nRows x nCols initialized to 0

template <class T>
T** newCMatrix(const int& nRows, const int& nCols){
	T **matrix = newCVector<T*>(nRows);
	matrix[0] = newCVector<T>(nRows*nCols);
	for(int row=1; row < nRows; ++row) matrix[row]=matrix[row-1]+nCols;
	return matrix;
}

// Deletes a C vector

template <class T>
void deleteCVector(T* vector){
	delete [] vector;
	vector = NULL;
}

// Deletes a C matrix

template <class T>
void deleteCMatrix(T** matrix){
	deleteCVector(matrix[0]); // release buffer space
	delete [] matrix;
	matrix = NULL;
}



//void extractFromProposalPDF( double& lth_candidate,     double& lph_candidate,     double& ltp_candidate,
//                             double& lth,               double& lph,               double& ltp,
//                             double& proposalWidth_lth, double& proposalWidth_lph, double& proposalWidth_ltp ) {
//
//  if ( proposalWidth_lth < 0. || proposalWidth_lph < 0. || proposalWidth_ltp < 0. ) {
//       do {
//           lth_candidate = gRandom->Gaus( lth, proposalWidthBurnIn ); // Gaussian proposal pdf with "large" sigma and positivity constraints
//           lph_candidate = gRandom->Gaus( lph, proposalWidthBurnIn ); // 0.10 all default
//           ltp_candidate = gRandom->Gaus( ltp, proposalWidthBurnIn );
//       }
//
///*    do {
//        lth_candidate = gRandom->Uniform(-1,1); // Gaussian proposal pdf with "large" sigma and positivity constraints
//        lph_candidate = gRandom->Uniform(-1,1);
//        ltp_candidate = gRandom->Uniform(-0.707106781186548,0.707106781186548);
//    }
//*/	while ( TMath::Abs( lph_candidate ) > 0.5*( 1 + lth_candidate ) || lth_candidate*lth_candidate + 2.*ltp_candidate*ltp_candidate > 1
//            || TMath::Abs( ltp_candidate ) > 0.5*( 1 - lph_candidate )
//            || (  (1.+2.*lph_candidate)*(1.+2.*lph_candidate) + 2.*ltp_candidate*ltp_candidate > 1 && lph_candidate < -1./3. ) );
//  }
//  else {
//      lth_candidate = gRandom->Gaus ( lth, proposalWidth_lth ); // Gaussian proposal pdf with possibly smaller sigmas
//      lph_candidate = gRandom->Gaus ( lph, proposalWidth_lph ); // and no positivity constraints
//      ltp_candidate = gRandom->Gaus ( ltp, proposalWidth_ltp );
//  }
//
//}


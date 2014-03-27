#ifndef __GlobalVar_h__
#define __GlobalVar_h__

#include "TLorentzVector.h"
#include "TMath.h"
#include <algorithm>

namespace NRQCDvars{

//'global' variables

//worked for Psi2S only fit:
//const double proposalWidthBurnIn_f=5e-3;
//const double proposalWidthBurnIn_R=2.5e-1;
//const double proposalWidthBurnIn_Np_relToUnc=1e-1;//width proposalWidthBurnIn_Np = Np_uncertainty*proposalWidthBurnIn_Np_relToUnc

const double proposalWidthBurnIn_f=5e-3;
const double proposalWidthBurnIn_R=5.e-2;//1.e-1
const double proposalWidthBurnIn_Np_relToUnc=5e-2;//width proposalWidthBurnIn_Np = Np_uncertainty*proposalWidthBurnIn_Np_relToUnc


bool debug=false;
enum {charm, bottom};
enum {Metropolis, MetropolisHastings}; // modified by Joao: choice of kernel
enum {MH, Minuit}; //choice of Minimizer
enum {PSI_1S, CHIC1_1P, CHIC2_1P, PSI_2S, UPS_1S, CHIB1_1P, CHIB2_1P, UPS_2S,
	  CHIB1_2P, CHIB2_2P, UPS_3S, CHIB1_3P, CHIB2_3P};
const int nStates = (CHIB2_3P-PSI_1S)+1;// modified by Joao: always use last state - first state to set this number
enum {CrossSection, Lamth, Lamph, Lamtp, CrossSectionRatio, FeedDownFraction};
const int nMeasurementIDs = (FeedDownFraction-CrossSection)+1; // modified by Joao: : always use last measur. - first measur. to set this number

//enum {CMS, LHCb, ATLAS, ALICE, CDF};
//const int nExperiments = (CDF-CMS)+1; // modified by Joao: always use last exp. - first exp. to set this number
//Char_t *ExpName[nExperiments] = {"CMS", "LHCb", "ATLAS", "ALICE", "CDF"};
//bool isAbsRapExp[nExperiments] = {true, false, true, false, true};

enum {CMS2010, CMS2011, CMS2012, LHCb2010, LHCb2011, LHCb2012, ATLAS2010, ATLAS2011, ATLAS2012};
const int nExperiments = (ATLAS2012-CMS2010)+1; // modified by Joao: always use last exp. - first exp. to set this number
Char_t *ExpName[nExperiments] = {"CMS2010", "CMS2011", "CMS2012", "LHCb2010", "LHCb2011", "LHCb2012", "ATLAS2010", "ATLAS2011", "ATLAS2012"};
Char_t *ExpNameTex[nExperiments] = {"CMS 2010", "CMS 2011", "CMS 2012", "LHCb 2010", "LHCb 2011", "LHCb 2012", "ATLAS 2010", "ATLAS 2011", "ATLAS 2012"};
bool isAbsRapExp[nExperiments] = {true, true, true, false, false, false, true, true, true};

Char_t *StateName[nStates] = {"PSI_1S", "CHIC1_1P", "CHIC2_1P", "PSI_2S", "UPS_1S", "CHIB1_1P", "CHIB2_1P", "UPS_2S", "CHIB1_2P", "CHIB2_2P", "UPS_3S", "CHIB1_3P", "CHIB2_3P"};
Char_t *StateNameTex[nStates] = {"#psi(1#it{S})", "#chi_{c1}(1#it{P})", "#chi_{c2}(1#it{P})", "#psi(2#it{S})", "#Upsilon(1#it{S})", "#chi_{b1}(1#it{P})", "#chi_{b2}(1#it{P})", "#Upsilon(2#it{S})", "#chi_{b1}(2#it{P})", "#chi_{b2}(2#it{P})", "#Upsilon(3#it{S})", "#chi_{b1}(3#it{P})", "#chi_{b2}(3#it{P})"};
Char_t *MeasurementIDName[nMeasurementIDs] = {"CrossSection", "Lamth", "Lamph", "Lamtp", "CrossSectionRatio", "FeedDownFraction"};
Char_t *MeasurementIDNameTex[nMeasurementIDs] = {"d#sigma / d#it{p}_{T}d#it{y} [nb/GeV]", "#lambda_{#vartheta}^{HX}", "#lambda_{#varphi}^{HX}", "#lambda_{#vartheta#varphi}^{HX}", "CrossSectionRatio", "FeedDownFraction"};
Double_t mass[nStates] = {3.096916, 3.51066, 3.55620, 3.686109, 9.46030, 9.89278, 9.91221, 10.02326, 10.25546,
		                   10.26865, 10.3552, 10.530, 10.530};

enum {quID_S, quID_P1, quID_P2}; //Definition of QuantumID (S, P1, P2 states)
int StateQuantumID[nStates]={quID_S, quID_P1, quID_P2, quID_S, quID_S, quID_P1, quID_P2, quID_S, quID_P1, quID_P2, quID_S, quID_P1, quID_P2};

const int nModelSystematicScales = 0;
const int nDataSystematicScales = nExperiments;

//const int nColorChannels_S=4;//includes CS
//const int nColorChannels_P=2;//includes CS
//Char_t *ColorChannelNameTexS[nColorChannels_S] = {"^{3}#it{S}_{1}^{[1]}", "^{1}#it{S}_{0}^{[8]}", "^{3}#it{S}_{1}^{[8]}", "^{3}#it{P}_{J}^{[8]}"};
//Char_t *ColorChannelNameTexP[nColorChannels_P] = {"^{3}#it{P}_{J}^{[1]}", "^{3}#it{S}_{1}^{[8]}"};

const int nColorChannels_S=3;//includes CS
const int nColorChannels_P=3;//includes CS
Char_t *ColorChannelNameTexS[nColorChannels_S] = {"^{3}#it{S}_{1}^{[1]}", "^{1}#it{S}_{0}^{[8]}", "^{3}#it{S}_{1}^{[8]}"};
Char_t *ColorChannelNameTexP[nColorChannels_P] = {"^{3}#it{P}_{J}^{[1]}", "^{1}#it{S}_{0}^{[8]}", "^{3}#it{S}_{1}^{[8]}"};

//const int nColorChannels_S=2;//includes CS
//const int nColorChannels_P=2;//includes CS
//Char_t *ColorChannelNameTexS[nColorChannels_S] = {"^{3}S#it{S}_{1}^{[1]}", "^{1}#it{S}_{0}^{[8]}"};
//Char_t *ColorChannelNameTexP[nColorChannels_P] = {"^{3}#it{P}_{J}^{[1]}", "^{3}#it{S}_{1}^{[8]}"};

const int nColorChannels = std::max(nColorChannels_S, nColorChannels_P);
//const int nColorChannels4States[nStates]={nColorChannels_S, nColorChannels_P, nColorChannels_P, nColorChannels_S, nColorChannels_S, nColorChannels_P, nColorChannels_P, nColorChannels_S, nColorChannels_P, nColorChannels_P, nColorChannels_S, nColorChannels_P, nColorChannels_P};

double FeedDownBranchingRatio[nStates][nStates]={
		//one row for each state; elements: BR of BR(x->y), [iDaughter][iMother]
		{0,   34.8,   19.8,   60.3,     0,   0,     0,      0,       0,     0,   0,   0,   0},   // 1...Jpsi
		{0,   0,      0,      9.3,      0,   0,     0,      0,       0,     0,   0,   0,   0},   // 2...chic1
		{0,   0,      0,      8.76,     0,   0,     0,      0,       0,     0,   0,   0,   0},   // 3...chic2
		{0,   0,      0,      0,        0,   0,     0,      0,       0,     0,   0,   0,   0},   // 4...Psi'


		{0,   0,      0,      0,        0,   35,   22,   26.7,   10.13,   8.2,   6.6,   0,   0},   // 5...Y1S
		{0,   0,      0,      0,        0,   0,     0,    6.9,       0,     0,     0,   0,   0},   // 6...chib1P1
		{0,   0,      0,      0,        0,   0,     0,   7.15,       0,     0,     0,   0,   0},   // 7...chib1P2
		{0,   0,      0,      0,        0,   0,     0,      0,      21,  16.2,  10.6,   0,   0},   // 8...Y2S
		{0,   0,      0,      0,        0,   0,     0,      0,       0,     0,  12.6,   0,   0},   // 9...chib2P1
		{0,   0,      0,      0,        0,   0,     0,      0,       0,     0,  13.1,   0,   0},   // 10..chib2P2
		{0,   0,      0,      0,        0,   0,     0,      0,       0,     0,     0,   0,   0},   // 11..Y3S
		{0,   0,      0,      0,        0,   0,     0,      0,       0,     0,     0,   0,   0},   // 12..chib3P1
		{0,   0,      0,      0,        0,   0,     0,      0,       0,     0,     0,   0,   0}    // 13..chib3P2
};


double errFeedDownBranchingRatio[nStates][nStates]={
		//one row for each state; elements: BR of BR(x->y)
		{0,   1.5,   0.8,   0.7,        0,   0,   0,   0,   0,   0,   0,   0,   0},   // 1...Jpsi
		{0,   0,     0,     0.4,        0,   0,   0,   0,   0,   0,   0,   0,   0},   // 2...chic1
		{0,   0,     0,     0.34,       0,   0,   0,   0,   0,   0,   0,   0,   0},   // 3...chic2
		{0,   0,     0,     0,          0,   0,   0,   0,   0,   0,   0,   0,   0},   // 4...Psi'


		{0,   0,     0,     0,        0,   8,   4,   0.57,   1.36,  1.06,  0.16,   0,   0},   // 5...Y1S
		{0,   0,     0,     0,        0,   0,   0,    0.4,      0,     0,     0,   0,   0},   // 6...chib1P1
		{0,   0,     0,     0,        0,   0,   0,   0.35,      0,     0,     0,   0,   0},   // 7...chib1P2
		{0,   0,     0,     0,        0,   0,   0,      0,      4,   2.4,   0.8,   0,   0},   // 8...Y2S
		{0,   0,     0,     0,        0,   0,   0,      0,      0,     0,   1.2,   0,   0},   // 9...chib2P1
		{0,   0,     0,     0,        0,   0,   0,      0,      0,     0,   1.6,   0,   0},   // 10..chib2P2
		{0,   0,     0,     0,        0,   0,   0,      0,      0,     0,     0,   0,   0},   // 11..Y3S
		{0,   0,     0,     0,        0,   0,   0,      0,      0,     0,     0,   0,   0},   // 12..chib3P1
		{0,   0,     0,     0,        0,   0,   0,      0,      0,     0,     0,   0,   0}    // 13..chib3P2
};

//Info FeedDownBranchingRatioComment[N_STATES][N_STATES]={
//		//one row for each state; elements: BR of BR(x->y)
//		{0,   "gamma",   "gamma",            0,   0,   0,   0,   0,   0,   0,   0,   0},   // 1...Jpsi
//		{0,   0,   0,    "gamma",            0,   0,   0,   0,   0,   0,   0,   0,   0},   // 2...chic1
//		{0,   0,   0,    "gamma",            0,   0,   0,   0,   0,   0,   0,   0,   0},   // 3...chic2
//		{0,   0,   0, "anything",            0,   0,   0,   0,   0,   0,   0,   0,   0},   // 4...Psi'
//
//
//		{0,   0,   0,          0,            0,   "gamma",   "gamma",   "anything",   "anything",  "anything",                    "hadronic",  "not measured",  "not measured"},   // 5...Y1S
//		{0,   0,   0,          0,            0,         0,         0,      "gamma",            0,           0,            "upper limit only",               0,               0},   // 6...chib1P1
//		{0,   0,   0,          0,            0,         0,         0,      "gamma",            0,           0,            "upper limit only",               0,               0},   // 7...chib1P2
//		{0,   0,   0,          0,            0,         0,         0,            0,   "anything",  "anything",  "anything incl. gamma+gamma",  "not measured",  "not measured"},   // 8...Y2S
//		{0,   0,   0,          0,            0,         0,         0,            0,            0,           0,                       "gamma",               0,               0},   // 9...chib2P1
//		{0,   0,   0,          0,            0,         0,         0,            0,            0,           0,                       "gamma",               0,               0},   // 10..chib2P2
//		{0,   0,   0,          0,            0,         0,         0,            0,            0,           0,                             0,  "not measured",  "not measured"},   // 11..Y3S
//		{0,   0,   0,          0,            0,         0,         0,            0,            0,           0,                             0,               0,               0},   // 12..chib3P1
//		{0,   0,   0,          0,            0,         0,         0,            0,            0,           0,                             0,               0,               0}    // 13..chib3P2
//};

//TODO:
double ColorSingletME[nStates]={1.32, 1e-5,  1e-5, 0.76, 1., 1., 1., 1., 1., 1., 1., 1., 1.}; //ColorSinglet matrix elements taken from literature
double errColorSingletME[nStates]={0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}; //uncertainty on ColorSinglet matrix elements taken from literature
//const double pT_star=35.;//20;
const double pT_star_over_m=6.;
const double rap_star=0.;//0.5;

int randomSeed = 23101987;
const int nMaxCascades = 20;//Used to define the size of an array containing all decay cascades of a certain mother->daughter link
const int nMaxRapBins = 5;
const int nMaxPtBins = 50;

}

#endif

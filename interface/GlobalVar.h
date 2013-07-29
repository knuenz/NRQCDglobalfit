#ifndef __GlobalVar_h__
#define __GlobalVar_h__

#include "TLorentzVector.h"
#include "TMath.h"
#include <algorithm>

namespace NRQCDvars{

//'global' variables
//const int nStates=13;
//const int nExperiments = 5;
//const int nMeasurementIDs = 6;
const int nModelSystematicScales = 0;
const int nDataSystematicScales = 1;
const int nFractionDim = 3;     // modified by Joao: number of elements in the fraction (and width) vector

const double proposalWidthBurnIn_f=5e-3;
const double proposalWidthBurnIn_R=5e-3;//TODO: trial and error

bool debug=false;
enum {charm, bottom};
enum {Metropolis, MetropolisHastings}; // modified by Joao: choice of kernel
enum {PSI_1S, CHIC1_1P, CHIC2_1P, PSI_2S, UPS_1S, CHIB1_1P, CHIB2_1P, UPS_2S,
	  CHIB1_2P, CHIB2_2P, UPS_3S, CHIB1_3P, CHIB2_3P};
const int nStates = (CHIB2_3P-PSI_1S)+1;// modified by Joao: always use last state - first state to set this number
enum {CrossSection, Lamth, Lamph, Lamtp, CrossSectionRatio, FeedDownFraction};
const int nMeasurementIDs = (FeedDownFraction-CrossSection)+1; // modified by Joao: : always use last measur. - first measur. to set this number
enum {CMS, LHCb, ATLAS, ALICE, CDF};
const int nExperiments = (CDF-CMS)+1; // modified by Joao: always use last exp. - first exp. to set this number
Char_t *StateName[nStates] = {"PSI_1S", "CHIC1_1P", "CHIC2_1P", "PSI_2S", "UPS_1S", "CHIB1_1P", "CHIB2_1P", "UPS_2S", "CHIB1_2P", "CHIB2_2P", "UPS_3S", "CHIB1_3P", "CHIB2_3P"};
Char_t *StateNameTex[nStates] = {"#psi(1S)", "#chi_{c1}(1P)", "#chi_{c2}(1P)", "#psi(2S)", "#Upsilon(1S)", "#chi_{b1}(1P)", "#chi_{b2}(1P)", "#Upsilon(2S)", "#chi_{b1}(2P)", "#chi_{b2}(2P)", "#Upsilon(3S)", "#chi_{b1}(3P)", "#chi_{b2}(3P)"};
Char_t *MeasurementIDName[nMeasurementIDs] = {"CrossSection", "Lamth", "Lamph", "Lamtp", "CrossSectionRatio", "FeedDownFraction"};
Char_t *MeasurementIDNameTex[nMeasurementIDs] = {"#sigma [nB]", "#lambda_{#vartheta}^{HX}", "#lambda_{#varphi}^{HX}", "#lambda_{#vartheta#varphi}^{HX}", "CrossSectionRatio", "FeedDownFraction"};
Char_t *ExpName[nExperiments] = {"CMS", "LHCb", "ATLAS", "ALICE", "CDF"};
Double_t mass[nStates] = {3.096916, 3.51066, 3.55620, 3.686109, 9.46030, 9.89278, 9.91221, 10.02326, 10.25546,
		                   10.26865, 10.3552, 10.530, 10.530};

enum {quID_S, quID_P1, quID_P2}; //Definition of QuantumID (S, P1, P2 states)
int StateQuantumID[nStates]={quID_S, quID_P1, quID_P2, quID_S, quID_S, quID_P1, quID_P2, quID_S, quID_P1, quID_P2, quID_S, quID_P1, quID_P2};


const int nColorChannels_S=3;//includes CS
const int nColorChannels_P=2;//includes CS

Char_t *ColorChannelNameTexS[nColorChannels_S] = {"^{3}S_{1}^{(1)}", "^{3}S_{1}^{(8)}", "^{1}S_{0}^{(8)}"};//, "^{3}P_{1}^{(8)}"};
Char_t *ColorChannelNameTexP[nColorChannels_P] = {"^{3}P_{J}^{(1)}", "^{3}S_{1}^{(8)}"};

const int nColorChannels = std::max(nColorChannels_S, nColorChannels_P);
//const int nColorChannels4States[nStates]={nColorChannels_S, nColorChannels_P, nColorChannels_P, nColorChannels_S, nColorChannels_S, nColorChannels_P, nColorChannels_P, nColorChannels_S, nColorChannels_P, nColorChannels_P, nColorChannels_S, nColorChannels_P, nColorChannels_P};

double FeedDownBranchingRatio[nStates][nStates]={
		//one row for each state; elements: BR of BR(x->y), [iDaughter][iMother]
		{0,   34.1,   19.4,   58.7,     0,   0,     0,      0,       0,     0,   0,   0,   0},   // 1...Jpsi
		{0,   0,      0,      9.2,      0,   0,     0,      0,       0,     0,   0,   0,   0},   // 2...chic1
		{0,   0,      0,      8.69,     0,   0,     0,      0,       0,     0,   0,   0,   0},   // 3...chic2
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
		{0,   1.5,   0.8,   0.8,        0,   0,   0,   0,   0,   0,   0,   0,   0},   // 1...Jpsi
		{0,   0,     0,     0.4,        0,   0,   0,   0,   0,   0,   0,   0,   0},   // 2...chic1
		{0,   0,     0,     0.35,       0,   0,   0,   0,   0,   0,   0,   0,   0},   // 3...chic2
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
double ColorSingletME[nStates]={1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.}; //ColorSinglet matrix elements taken from literature
double errColorSingletME[nStates]={0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}; //uncertainty on ColorSinglet matrix elements taken from literature
const double pT_star=20;
const double rap_star=0.5;

int randomSeed = 23101987;
const int nMaxCascades = 20;//Used to define the size of an array containing all decay cascades of a certain mother->daughter link
const int nMaxRapBins = 5;
const int nMaxPtBins = 50;

}

#endif

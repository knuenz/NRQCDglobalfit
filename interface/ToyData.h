/*
 * ToyData.h
 *
 *  Created on: Jul 26, 2013
 *      Author: valentinknuenz
 */

#ifndef TOYDATA_H_
#define TOYDATA_H_


#include "TLorentzVector.h"
#include "TMath.h"
#include <algorithm>
#include "../interface/GlobalVar.h"

namespace toyData{


// TRUTH input

	bool isTRUTH_Ops=false;

	const int nStates_HardCoded=13;
	const int nColorChannels_HardCoded=3;


	double TRUTH_matrix[nStates_HardCoded][nColorChannels_HardCoded]={//sum(f_i_1,f_i_n-1)=1! f_i_0 > -1!
			{1.,     0.,     1.    },   // 1...Jpsi
			{1.,     1.,     0.    },   // 2...chic1
			{1.,     1.,     0.    },   // 3...chic2
			{1.,     1./2.,  1./2. },   // 4...Psi'

			{0.,     0.,     0.    },   // 5...Y1S
			{0.,     0.,     0.    },   // 6...chib1P1
			{0.,     0.,     0.    },   // 7...chib1P2
			{0.,     0.,     0.    },   // 8...Y2S
			{0.,     0.,     0.    },   // 9...chib2P1
			{0.,     0.,     0.    },   // 10..chib2P2
			{0.,     0.,     0.    },   // 11..Y3S
			{0.,     0.,     0.    },   // 12..chib3P1
			{0.,     0.,     0.    }    // 13..chib3P2
	};

	 //Cross section: relative uncertainties
	double Toy_totalUncertaintyPos=0.1;
	double Toy_totalUncertaintyNeg=0.1;
	double Toy_globalUncertaintyPos=0.05;
	double Toy_globalUncertaintyNeg=0.05;
	double Toy_polarizationUncertaintyLongHX=+0.2;
	double Toy_polarizationUncertaintyTransHX=-0.2;

	//Polarization: absolute uncertainties
	double Toy_totalUncertaintyPos_Lamth=0.1;
	double Toy_totalUncertaintyNeg_Lamth=0.1;
	double Toy_totalUncertaintyPos_Lamph=0.05;
	double Toy_totalUncertaintyNeg_Lamph=0.05;
	double Toy_totalUncertaintyPos_Lamtp=0.05;
	double Toy_totalUncertaintyNeg_Lamtp=0.05;



/// Definition of measurements and kinematic ranges
	const int Toy_nStates=4;
	const int Toy_nMeasurementIDs=4;
	const int Toy_nExperiments=1;

	const int Toy_nRapBins=5;
	const int Toy_nPtBins=10;

	double rapRange[Toy_nRapBins+1]={0, 0.5, 1., 1.5, 2., 4.};
	double pTRange[Toy_nRapBins][Toy_nPtBins+1]={{10,15,20,25,30,35,40,45,50,70,100},{10,15,20,25,30,35,40,45,50,70,100},{10,15,20,25,30,35,40,45,50,70,100},{10,15,20,25,30,35,40,45,50,70,100},{10,15,20,25,30,35,40,45,50,70,100}};





}

#endif /* TOYDATA_H_ */

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


	//Fractions: Used to generate toy data SDC and Lamth
	double TRUTH_matrix[nStates_HardCoded][nColorChannels_HardCoded]={//sum(f_i_1,f_i_n-1)=1! f_i_0 > -1!
			{1,     0.5,    0.5   },   // 1...Jpsi
			{1.,     0.5,    0.5   },   // 2...chic1
			{1.,     0.5,    0.5   },   // 3...chic2
			{26.5536, 0.913125, 0.0769133},   // 4...Psi'

			{1.,     0.5,    0.5   },   // 5...Y1S
			{1.,     0.5,    0.5   },   // 6...chib1P1
			{1.,     0.5,    0.5   },   // 7...chib1P2
			{1.,     0.5,    0.5   },   // 8...Y2S
			{1.,     0.5,    0.5   },   // 9...chib2P1
			{1.,     0.5,    0.5   },   // 10..chib2P2
			{20.9611, 0.412066, 0.583231},   // 11..Y3S
			{1.,     0.5,    0.5   },   // 12..chib3P1
			{1.,     0.5,    0.5   }    // 13..chib3P2
	};


	//LDMEs: Only used to calculate MeanPt of generated toy data SDC and Lamth points
	double TRUTH_matrix_O[nStates_HardCoded][nColorChannels_HardCoded]={
			{0., 0., 0.},   // 1...Jpsi
			{0., 0., 0.},   // 2...chic1
			{0., 0., 0.},   // 3...chic2
			{0.76, 0.0446935, 0.000101407},   // 4...Psi'

			{0., 0., 0.},   // 5...Y1S
			{0., 0., 0.},   // 6...chib1P1
			{0., 0., 0.},   // 7...chib1P2
			{0., 0., 0.},   // 8...Y2S
			{0., 0., 0.},   // 9...chib2P1
			{0., 0., 0.},   // 10..chib2P2
			{1, 0.030424, 0.000652217},   // 11..Y3S
			{0., 0., 0.},   // 12..chib3P1
			{0., 0., 0.}    // 13..chib3P2
	};


	 //Cross section: relative uncertainties
	double Toy_totalUncertaintyPos=1e-10;//0.1;
	double Toy_totalUncertaintyNeg=1e-10;//0.1;
	double Toy_globalUncertaintyPos=1e-10;//0.05;
	double Toy_globalUncertaintyNeg=1e-10;//0.05;
	double Toy_polarizationUncertaintyLongHX=1e-10;//+0.2;
	double Toy_polarizationUncertaintyTransHX=1e-10;//-0.2;

	//Polarization: absolute uncertainties
	double Toy_totalUncertaintyPos_Lamth=1e-10;//0.1;
	double Toy_totalUncertaintyNeg_Lamth=1e-10;//0.1;
	double Toy_totalUncertaintyPos_Lamph=1e-10;//0.05;
	double Toy_totalUncertaintyNeg_Lamph=1e-10;//0.05;
	double Toy_totalUncertaintyPos_Lamtp=1e-10;//0.05;
	double Toy_totalUncertaintyNeg_Lamtp=1e-10;//0.05;



/// Definition of measurements and kinematic ranges
	const int Toy_nStates=13;
	bool Toy_StatesToGenerate[Toy_nStates]={false, false, false, true, false, false, false, false, false, false, true, false, false};
	const int Toy_nMeasurementIDs=2;
	const int Toy_nExperiments=1;

	const int Toy_nRapBins=1;
	const int Toy_nPtBins=50;

	double rapRange[Toy_nRapBins+1]={0, 1.2};
	double pTRange[Toy_nStates][Toy_nPtBins+1]={
			{5,9,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85},
			{5,9,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85},
			{5,9,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85},

			{5,9,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85},

			{5,9,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85},
			{5,9,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85},
			{5,9,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85},
			{5,9,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85},
			{5,9,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85},
			{5,9,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85},

			{0,22.5,37.5,51,65,79,93,107,121,135,150,165,180,195,210,225,240},

			{5,9,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85},
			{5,9,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85},
};





}

#endif /* TOYDATA_H_ */

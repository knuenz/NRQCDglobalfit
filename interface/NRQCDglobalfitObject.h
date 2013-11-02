#ifndef __NRQCDglobalfitObject_h__
#define __NRQCDglobalfitObject_h__

#include "../interface/GlobalVar.h"

//#include "TString.h"
#include <string>
#include <cstring>
#include <vector>
#include <fstream>

//#define N_STATES 13

using namespace std;

// define types to make coding easier
typedef vector<vector<vector<double> > > dcube;
typedef vector<vector<vector<float> > > fcube;
typedef vector<vector<vector<int> > > icube;
typedef vector<vector<vector<bool> > > bcube;
typedef vector<vector<double> > dmatrix;
typedef vector<vector<float> > fmatrix;
typedef vector<vector<int> > imatrix;
typedef vector<vector<bool> > bmatrix;
typedef vector<double> dvector;
typedef vector<float> fvector;
typedef vector<int> ivector;
typedef vector<bool> bvector;
//

class NRQCDglobalfitObject {
 public:

	  // constructors
	  NRQCDglobalfitObject(){}
//	  NRQCDglobalfitObject(const int& numStates);
//	  NRQCDglobalfitObject(const NRQCDglobalfitObject &other);

//	  void initFractions(const int& nDim);
//	  NRQCDglobalfitObject(const int nrows, const int ncolumns);
  
	  void Dump(const int Input_nStates, bool DumpData, bool DumpModel); //Dumps all available information about the object

	  void setData(
			  int Input_nState,
			  int Input_nStateDenom,
			  double Input_CentralValue,
			  double Input_ErrStatPos, double Input_ErrStatNeg,
			  double Input_ErrSystPos, double Input_ErrSystNeg,
			  double Input_ErrTotPos, double Input_ErrTotNeg,
			  double Input_ErrGlobalPos, double Input_ErrGlobalNeg,
			  double Input_yMin, double Input_yMax, double Input_yMean,
			  double Input_pTMin, double Input_pTMax, double Input_pTMean,
			  dvector Input_PolCorrParams,
			  int Input_MeasurementID,
			  int Input_nExperiment,
			  string Input_ObjectID,
			  bool Input_isDataValid,
			  bool Input_isAbsRap
			  )
	  {
		  Object_nState = Input_nState;
		  Object_nStateDenom = Input_nStateDenom;
		  Object_CentralValue = Input_CentralValue;
		  Object_ErrStatPos = Input_ErrStatPos; Object_ErrStatNeg = Input_ErrStatNeg;
		  Object_ErrSystPos = Input_ErrSystPos; Object_ErrSystNeg = Input_ErrSystNeg;
		  Object_ErrTotPos = Input_ErrTotPos; Object_ErrTotNeg = Input_ErrTotNeg;
		  Object_ErrGlobalPos = Input_ErrGlobalPos; Object_ErrGlobalNeg = Input_ErrGlobalNeg;
		  Object_yMin = Input_yMin; Object_yMax = Input_yMax; Object_yMean = Input_yMean;
		  Object_pTMin = Input_pTMin; Object_pTMax = Input_pTMax; Object_pTMean = Input_pTMean;
		  Object_PolCorrParams = Input_PolCorrParams;
		  Object_MeasurementID = Input_MeasurementID;
		  Object_nExperiment = Input_nExperiment;
		  Object_ObjectID = Input_ObjectID;
		  Object_isDataValid = Input_isDataValid;
		  Object_isAbsRap = Input_isAbsRap;
	  }

	  void setModel(int Input_nState,
			  dvector Input_ShortDistanceCoef,
			  dvector Input_OctetCompLamth,
			  dvector Input_OctetCompLamph,
			  dvector Input_OctetCompLamtp,
			  dvector Input_ShortDistanceCoef_PointByPointSyst,
			  dvector Input_OctetCompLamth_PointByPointSyst,
			  dvector Input_OctetCompLamph_PointByPointSyst,
			  dvector Input_OctetCompLamtp_PointByPointSyst,
			  dmatrix Input_ShortDistanceCoef_globalSystPos,
			  dmatrix Input_ShortDistanceCoef_globalSystNeg,
			  dmatrix Input_OctetCompLamth_globalSystPos,
			  dmatrix Input_OctetCompLamth_globalSystNeg,
			  dmatrix Input_OctetCompLamph_globalSystPos,
			  dmatrix Input_OctetCompLamph_globalSystNeg,
			  dmatrix Input_OctetCompLamtp_globalSystPos,
			  dmatrix Input_OctetCompLamtp_globalSystNeg,

			  bool Input_isModelValid) //setModel has to be used consecutively for all states starting with nState=0
	  {
		  Object_ShortDistanceCoef.push_back(Input_ShortDistanceCoef);
		  Object_OctetCompLamth.push_back(Input_OctetCompLamth);
		  Object_OctetCompLamph.push_back(Input_OctetCompLamph);
		  Object_OctetCompLamtp.push_back(Input_OctetCompLamtp);
		  Object_ShortDistanceCoef_PointByPointSyst.push_back(Input_ShortDistanceCoef_PointByPointSyst);
		  Object_OctetCompLamth_PointByPointSyst.push_back(Input_OctetCompLamth_PointByPointSyst);
		  Object_OctetCompLamph_PointByPointSyst.push_back(Input_OctetCompLamph_PointByPointSyst);
		  Object_OctetCompLamtp_PointByPointSyst.push_back(Input_OctetCompLamtp_PointByPointSyst);
		  Object_ShortDistanceCoef_globalSystPos.push_back(Input_ShortDistanceCoef_globalSystPos);
		  Object_ShortDistanceCoef_globalSystNeg.push_back(Input_ShortDistanceCoef_globalSystNeg);
		  Object_OctetCompLamth_globalSystPos.push_back(Input_OctetCompLamth_globalSystPos);
		  Object_OctetCompLamth_globalSystNeg.push_back(Input_OctetCompLamth_globalSystNeg);
		  Object_OctetCompLamph_globalSystPos.push_back(Input_OctetCompLamph_globalSystPos);
		  Object_OctetCompLamph_globalSystNeg.push_back(Input_OctetCompLamph_globalSystNeg);
		  Object_OctetCompLamtp_globalSystPos.push_back(Input_OctetCompLamtp_globalSystPos);
		  Object_OctetCompLamtp_globalSystNeg.push_back(Input_OctetCompLamtp_globalSystNeg);
		  Object_isModelValid = Input_isModelValid;
	  }

	  void setInvalidModel(int Input_nStates, int Input_nColorChannels, int Input_nModelSystematicScales);

	  void setisDataValid(bool Input_isDataValid) {Object_isDataValid = Input_isDataValid;}
	  void setisModelValid(bool Input_isModelValid) {Object_isModelValid = Input_isModelValid;}

	  int getState() {return Object_nState;}
	  int getStateDenom() {return Object_nStateDenom;}
	  int getMeasurementID() {return Object_MeasurementID;}
	  int getExperiment() {return Object_nExperiment;}
	  double getCentralValue() {return Object_CentralValue;}
	  double getErrStatPos() {return Object_ErrStatPos;}
	  double getErrStatNeg() {return Object_ErrStatNeg;}
	  double getErrStat() {return (Object_ErrStatPos+Object_ErrStatNeg)/2.;}
	  double getErrSystPos() {return Object_ErrSystPos;}
	  double getErrSystNeg() {return Object_ErrSystNeg;}
	  double getErrSyst() {return (Object_ErrSystPos+Object_ErrSystNeg)/2.;}
	  double getErrTotPos() {return Object_ErrTotPos;}
	  double getErrTotNeg() {return Object_ErrTotNeg;}
	  double getErrTot() {return (fabs(Object_ErrTotPos)+fabs(Object_ErrTotNeg))/2.;}
	  double getErrGlobalPos() {return Object_ErrGlobalPos;}
	  double getErrGlobalNeg() {return Object_ErrGlobalNeg;}
	  double getErrGlobal() {return (Object_ErrGlobalPos+Object_ErrGlobalNeg)/2.;}
	  double getyMin() {return Object_yMin;}
	  double getyMax() {return Object_yMax;}
	  double getyMean() {return Object_yMean;}
	  double getpTMin() {return Object_pTMin;}
	  double getpTMax() {return Object_pTMax;}
	  double getpTMean() {return Object_pTMean;}
	  dvector getPolCorrParams() {return Object_PolCorrParams;}
	  string getObjectID() {return Object_ObjectID;}
	  bool getisDataValid() {return Object_isDataValid;}
	  bool getisModelValid() {return Object_isModelValid;}
	  bool getisAbsRap() {return Object_isAbsRap;}
//	  void getFractionValues();

	  dvector getShortDistanceCoef(int Input_nState) {return Object_ShortDistanceCoef[Input_nState];}
	  dvector getOctetCompLamth(int Input_nState) {return Object_OctetCompLamth[Input_nState];}
	  dvector getOctetCompLamph(int Input_nState) {return Object_OctetCompLamph[Input_nState];}
	  dvector getOctetCompLamtp(int Input_nState) {return Object_OctetCompLamtp[Input_nState];}

	  dvector getShortDistanceCoef_PointByPointSyst(int Input_nState) {return Object_ShortDistanceCoef_PointByPointSyst[Input_nState];} // Use in function getDirectProduction() not yet implemented
	  dvector getOctetCompLamth_PointByPointSyst(int Input_nState) {return Object_OctetCompLamth_PointByPointSyst[Input_nState];}       // Use in function getDirectProduction() not yet implemented
	  dvector getOctetCompLamph_PointByPointSyst(int Input_nState) {return Object_OctetCompLamph_PointByPointSyst[Input_nState];}       // Use in function getDirectProduction() not yet implemented
	  dvector getOctetCompLamtp_PointByPointSyst(int Input_nState) {return Object_OctetCompLamtp_PointByPointSyst[Input_nState];}       // Use in function getDirectProduction() not yet implemented

	  dmatrix getShortDistanceCoef_globalSystPos(int Input_nState) {return Object_ShortDistanceCoef_globalSystPos[Input_nState];}
	  dmatrix getShortDistanceCoef_globalSystNeg(int Input_nState) {return Object_ShortDistanceCoef_globalSystNeg[Input_nState];}
	  dmatrix getOctetCompLamth_globalSystPos(int Input_nState) {return Object_OctetCompLamth_globalSystPos[Input_nState];}
	  dmatrix getOctetCompLamth_globalSystNeg(int Input_nState) {return Object_OctetCompLamth_globalSystNeg[Input_nState];}
	  dmatrix getOctetCompLamph_globalSystPos(int Input_nState) {return Object_OctetCompLamph_globalSystPos[Input_nState];}
	  dmatrix getOctetCompLamph_globalSystNeg(int Input_nState) {return Object_OctetCompLamph_globalSystNeg[Input_nState];}
	  dmatrix getOctetCompLamtp_globalSystPos(int Input_nState) {return Object_OctetCompLamtp_globalSystPos[Input_nState];}
	  dmatrix getOctetCompLamtp_globalSystNeg(int Input_nState) {return Object_OctetCompLamtp_globalSystNeg[Input_nState];}
	  dvector getDirectProduction(int Input_nState, dmatrix &Op, dmatrix &Np_BR, dmatrix &Np_US, bool returnMPDetails, dmatrix &directProductionMatrix);
	  dvector getPromptProduction(int Input_nState, dmatrix &Op, dmatrix &Np_BR, dmatrix &Np_US, bool returnMPDetails, dcube &directProductionCube, dmatrix &promptProductionMatrix);
	  double getCorrPromptCrossSect(dvector &PromptCrossSect, dmatrix &Np_BR, dmatrix &Np_US, bool returnMPDetails, double &polCorrFactor);
	  dvector getObjectLikelihood(dmatrix &Op, dmatrix &Np_BR, dmatrix &Np_US, bool returnMPDetails, dcube &directProductionCube, dmatrix &promptProductionMatrix, double &polCorrFactor);

	  void setState(const int& state) {Object_nState=state;}
	  void setStateDenom(const int& stateDenom) {Object_nStateDenom=stateDenom;}
	  void setMeasurementID(const int& measurementID) {Object_MeasurementID=measurementID;}
	  void setExperiment(const int& nMeasurement) {Object_nExperiment=nMeasurement;}
	  void setCentralValue(const double& value) {Object_CentralValue=value;}
	  void setErrStatPos(const double& ErrStatPos) {Object_ErrStatPos=ErrStatPos;}
	  void setErrStatNeg(const double& ErrStatNeg) {Object_ErrStatNeg=ErrStatNeg;}
	  void setErrSystPos(const double& ErrSystPos) {Object_ErrSystPos=ErrSystPos;}
	  void setErrSystNeg(const double& ErrSystNeg) {Object_ErrSystNeg=ErrSystNeg;}
	  void setErrTotPos(const double& ErrTotPos) {Object_ErrTotPos=ErrTotPos;}
	  void setErrTotNeg(const double& ErrTotNeg) {Object_ErrTotNeg=ErrTotNeg;}
	  void setErrGlobalPos(const double& ErrGlobalPos) {Object_ErrGlobalPos=ErrGlobalPos;}
	  void setErrGlobalNeg(const double& ErrGlobalNeg) {Object_ErrGlobalNeg=ErrGlobalNeg;}
	  void setyMin(const double& yMin) {Object_yMin=yMin;}
	  void setyMax(const double& yMax) {Object_yMax=yMax;}
	  void setyMean(const double& yMean) {Object_yMean=yMean;}
	  void setpTMin(const double& pTMin) {Object_pTMin=pTMin;}
	  void setpTMax(const double& pTMax) {Object_pTMax=pTMax;}
	  void setpTMean(const double& pTMean) {Object_pTMean=pTMean;}
	  void setPolCorrParams(const dvector& PolCorrParams) {Object_PolCorrParams=PolCorrParams;}
	  void setObjectID(const string& ObjectID) {Object_ObjectID=ObjectID;}
	  void setisDataValid(const bool& DataValidFlag) {Object_isDataValid=DataValidFlag;}
	  void setisModelValid(const bool& ModelValidFlag) {Object_isModelValid=ModelValidFlag;}
	  void setisAbsRap(const bool& isAbsRapFlag) {Object_isAbsRap=isAbsRapFlag;}
//	  void setFractionsDim(const int& nDim ) {nFractionDimension=nDim;}
//	  void setFractionsBasis(const int& numOperatorStates);

//	  dvector setShortDistanceCoef(int Input_nState) {return Object_ShortDistanceCoef[Input_nState];}
//	  dvector setOctetCompLamth(int Input_nState) {return Object_OctetCompLamth[Input_nState];}
//	  dvector setOctetCompLamph(int Input_nState) {return Object_OctetCompLamph[Input_nState];}
//	  dvector setOctetCompLamtp(int Input_nState) {return Object_OctetCompLamtp[Input_nState];}
//
//	  dvector setShortDistanceCoef_PointByPointSyst(int Input_nState) {return Object_ShortDistanceCoef_PointByPointSyst[Input_nState];} // Use in function setDirectProduction() not yet implemented
//	  dvector setOctetCompLamth_PointByPointSyst(int Input_nState) {return Object_OctetCompLamth_PointByPointSyst[Input_nState];}       // Use in function setDirectProduction() not yet implemented
//	  dvector setOctetCompLamph_PointByPointSyst(int Input_nState) {return Object_OctetCompLamph_PointByPointSyst[Input_nState];}       // Use in function setDirectProduction() not yet implemented
//	  dvector setOctetCompLamtp_PointByPointSyst(int Input_nState) {return Object_OctetCompLamtp_PointByPointSyst[Input_nState];}       // Use in function setDirectProduction() not yet implemented
//
//	  dmatrix setShortDistanceCoef_globalSystPos(int Input_nState) {return Object_ShortDistanceCoef_globalSystPos[Input_nState];}
//	  dmatrix setShortDistanceCoef_globalSystNeg(int Input_nState) {return Object_ShortDistanceCoef_globalSystNeg[Input_nState];}
//	  dmatrix setOctetCompLamth_globalSystPos(int Input_nState) {return Object_OctetCompLamth_globalSystPos[Input_nState];}
//	  dmatrix setOctetCompLamth_globalSystNeg(int Input_nState) {return Object_OctetCompLamth_globalSystNeg[Input_nState];}
//	  dmatrix setOctetCompLamph_globalSystPos(int Input_nState) {return Object_OctetCompLamph_globalSystPos[Input_nState];}
//	  dmatrix setOctetCompLamph_globalSystNeg(int Input_nState) {return Object_OctetCompLamph_globalSystNeg[Input_nState];}
//	  dmatrix setOctetCompLamtp_globalSystPos(int Input_nState) {return Object_OctetCompLamtp_globalSystPos[Input_nState];}
//	  dmatrix setOctetCompLamtp_globalSystNeg(int Input_nState) {return Object_OctetCompLamtp_globalSystNeg[Input_nState];}


	  friend ostream& operator<< (ostream& out, dvector& my_vector);
	  friend ostream& operator<< (ostream& out, dmatrix& my_matrix);
	  friend ostream& operator<< (ostream& out, dcube& my_cube);
	  friend ostream& operator<< (ostream& out, fvector& my_vector);
	  friend ostream& operator<< (ostream& out, fmatrix& my_matrix);
	  friend ostream& operator<< (ostream& out, fcube& my_cube);
	  friend ostream& operator<< (ostream& out, ivector& my_vector);
	  friend ostream& operator<< (ostream& out, imatrix& my_matrix);
	  friend ostream& operator<< (ostream& out, icube& my_cube);
	  friend istream& operator>> (istream& in, dvector& my_vector);
	  friend istream& operator>> (istream& in, dmatrix& my_matrix);
	  friend istream& operator>> (istream& in, dcube& my_cube);
	  friend istream& operator>> (istream& in, fvector& my_vector);
	  friend istream& operator>> (istream& in, fmatrix& my_matrix);
	  friend istream& operator>> (istream& in, fcube& my_cube);
	  friend istream& operator>> (istream& in, ivector& my_vector);
	  friend istream& operator>> (istream& in, imatrix& my_matrix);
	  friend istream& operator>> (istream& in, icube& my_cube);
	  friend ostream& operator<< (ostream& out, NRQCDglobalfitObject& fit_object);
	  friend istream& operator>> (istream& in, NRQCDglobalfitObject &fit_object);

	  void TestFunction(double &internalReturn){internalReturn=4.;}


 protected:

// Measurement members
  int Object_nState;
  int Object_nStateDenom;//State in the denominator of a cross-section ratio or feed-down fraction
  int Object_MeasurementID; //0...Cross section, 1...LamthHX, 2...LamphHX, 3...LamtpHX, 4...Cross section ratio, 5...Feed-down fraction
  int Object_nExperiment;
  double Object_CentralValue;
  double Object_ErrStatPos, Object_ErrStatNeg; //absolute uncertainties
  double Object_ErrSystPos, Object_ErrSystNeg; //absolute uncertainties
  double Object_ErrTotPos, Object_ErrTotNeg; //absolute uncertainties
  double Object_ErrGlobalPos, Object_ErrGlobalNeg; //absolute uncertainties
  double Object_yMin, Object_yMax, Object_yMean;
  double Object_pTMin, Object_pTMax, Object_pTMean;

  dvector Object_PolCorrParams;//Containing at [0] the cross section for LongHX, at [1] the cross section for TransHX
  string Object_ObjectID; // Experiment=XXX Energy=XXXTeV Observable=XXX Unit=XXX State=XXX
  bool Object_isDataValid;
  bool Object_isAbsRap;

//Model members
  dmatrix Object_ShortDistanceCoef; // [nMothers][nColorChannels], vector<double> ShortDistanceCoef are the Color Octet Components of one nState
  dmatrix Object_OctetCompLamth; // [nMothers][nColorChannels], vector<double> OctetCompLamth are the Color Octet Components of one nState
  dmatrix Object_OctetCompLamph; // [nMothers][nColorChannels], vector<double> OctetCompLamph are the Color Octet Components of one nState
  dmatrix Object_OctetCompLamtp; // [nMothers][nColorChannels], vector<double> OctetCompLamtp are the Color Octet Components of one nState

  dmatrix Object_ShortDistanceCoef_PointByPointSyst; // [nMothers][nColorChannels]
  dmatrix Object_OctetCompLamth_PointByPointSyst; // [nMothers][nColorChannels]
  dmatrix Object_OctetCompLamph_PointByPointSyst; // [nMothers][nColorChannels]
  dmatrix Object_OctetCompLamtp_PointByPointSyst; // [nMothers][nColorChannels]

  dcube Object_ShortDistanceCoef_globalSystPos; // [nMothers][nColorChannels][nModelSystematicScales] // Default + syst
  dcube Object_ShortDistanceCoef_globalSystNeg; // [nMothers][nColorChannels][nModelSystematicScales] // Default - syst
  dcube Object_OctetCompLamth_globalSystPos; // [nMothers][nColorChannels][nModelSystematicScales] // Default + syst
  dcube Object_OctetCompLamth_globalSystNeg; // [nMothers][nColorChannels][nModelSystematicScales] // Default - syst
  dcube Object_OctetCompLamph_globalSystPos; // [nMothers][nColorChannels][nModelSystematicScales] // Default + syst
  dcube Object_OctetCompLamph_globalSystNeg; // [nMothers][nColorChannels][nModelSystematicScales] // Default - syst
  dcube Object_OctetCompLamtp_globalSystPos; // [nMothers][nColorChannels][nModelSystematicScales] // Default + syst
  dcube Object_OctetCompLamtp_globalSystNeg; // [nMothers][nColorChannels][nModelSystematicScales] // Default - syst

  bool Object_isModelValid;


  // auxiliary functions
  void dmatrix_copy(dmatrix& original, dmatrix& copy);
  dmatrix* dmatrix_allocate(const int nrow, const int ncolumn, const int value=0);
  void addRows(dmatrix& matrix, const int nrows = 1, const int value = 0);
  void addColumns(dmatrix& matrix, const int ncols = 1, const int value = 0);
};

// Creates a vector with nEls initialized to 0

template <class T>
T* newCVector(const int& nEls);

// Creates a matrix with nRows x nCols initialized to 0

template <class T>
T** newCMatrix(const int& nRows, const int& nCols);

// Deletes a vector

template <class T>
void deleteCVector(T* vector);

// Deletes a matrix

template <class T>
void deleteCMatrix(T** matrix);

#endif

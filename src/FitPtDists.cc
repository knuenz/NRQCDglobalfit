/*
 * FitPtDists.cc
 *
 *  Created on: Oct 10, 2013
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
#include "TMarker.h"
#include "Riostream.h"
#include "TString.h"
#include "TROOT.h"
#include "TMath.h"
#include "TLine.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TRotation.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TEllipse.h"
#include "TPad.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLatex.h"


using namespace NRQCDvars;

double func_pT_gen(double* x, double* par);
double func_pT_gen_ratio(double* x, double* par);
double exponential(double* x, double* par);

bool scale_pT_over_m=true;
bool scale_p_over_m=false;
bool scale_mT_over_m=false;

int main(int argc, char** argv) {

  	Char_t *DataID = "Default";
	Char_t *storagedir = "Default"; //Storage Directory

  	for( int i=0;i < argc; ++i ) {
		if(std::string(argv[i]).find("DataID") != std::string::npos) {char* DataIDchar = argv[i]; char* DataIDchar2 = strtok (DataIDchar, "="); DataID = DataIDchar2; cout<<"DataID = "<<DataID<<endl;}
		if(std::string(argv[i]).find("storagedir") != std::string::npos) {char* storagedirchar = argv[i]; char* storagedirchar2 = strtok (storagedirchar, "="); storagedir = storagedirchar2; cout<<"storagedir = "<<storagedir<<endl;}
  	}


  	gROOT->SetBatch();

	delete gRandom;
	gRandom = new TRandom3(NRQCDvars::randomSeed);  // better random generator


	// TODO: add loop for nSystematics (now nSystematics will be 0)

	char outname[200];
	char inname[200];
	char datadirname[200];
	char datadirfigname[200];
	sprintf(datadirname,"%s/DataID", storagedir);
	gSystem->mkdir(datadirname);
	sprintf(datadirname,"%s/%s",datadirname,DataID);
	gSystem->mkdir(datadirname);
	sprintf(datadirfigname,"%s/Figures", datadirname);
	gSystem->mkdir(datadirfigname);


	char xTitle[200];
	if(scale_pT_over_m) sprintf(xTitle, "dimuon p_{T}/m");
	if(scale_p_over_m) sprintf(xTitle, "dimuon p/m");
	if(scale_mT_over_m) sprintf(xTitle, "dimuon m_{T}/m");

	char yTitle[200];
	sprintf(yTitle, "d#sigma/(dp_{T}dy) [nb/GeV]");

	TGraphAsymmErrors* PtDataDists[NRQCDvars::nStates][NRQCDvars::nMeasurementIDs][NRQCDvars::nExperiments][NRQCDvars::nMaxRapBins];

	bool doesAnyDataExist_all[NRQCDvars::nStates][NRQCDvars::nMeasurementIDs][NRQCDvars::nExperiments];
	bool doesAnyDataExist_rap[NRQCDvars::nStates][NRQCDvars::nMeasurementIDs][NRQCDvars::nExperiments][NRQCDvars::nMaxRapBins];

	double rapMin[NRQCDvars::nStates][NRQCDvars::nMeasurementIDs][NRQCDvars::nExperiments][NRQCDvars::nMaxRapBins];
	double rapMax[NRQCDvars::nStates][NRQCDvars::nMeasurementIDs][NRQCDvars::nExperiments][NRQCDvars::nMaxRapBins];

	for(int iState=0;iState<NRQCDvars::nStates;iState++){
		for(int iMeasurementID=0;iMeasurementID<NRQCDvars::nMeasurementIDs;iMeasurementID++){
			for(int iExperiment=0;iExperiment<NRQCDvars::nExperiments;iExperiment++){
				for(int iRap = 0; iRap < NRQCDvars::nMaxRapBins; iRap++){

					doesAnyDataExist_rap[iState][iMeasurementID][iExperiment][iRap]=false;

				}
			}
		}
	}



	for(int iState=0;iState<NRQCDvars::nStates;iState++){
		for(int iMeasurementID=0;iMeasurementID<NRQCDvars::nMeasurementIDs;iMeasurementID++){
			for(int iExperiment=0;iExperiment<NRQCDvars::nExperiments;iExperiment++){

				doesAnyDataExist_all[iState][iMeasurementID][iExperiment]=false;
				if(iExperiment>2 || iMeasurementID!=0) continue;

				    	cout<<"build TGraph from data points for state="<<StateName[iState]<<", iMeasurementID="<<MeasurementIDName[iMeasurementID]<<", iExperiment="<<ExpName[iExperiment]<<endl;


				    	NRQCDglobalfitObject *DataObject[NRQCDvars::nMaxRapBins][NRQCDvars::nMaxPtBins];
				    	NRQCDglobalfitObject *ModelObject[NRQCDvars::nMaxRapBins][NRQCDvars::nMaxPtBins];
				    	bool doesDataExist[NRQCDvars::nMaxRapBins][NRQCDvars::nMaxPtBins];
						bool doesAnyDataExist=false;


						//cout<<"Reading in all DataObjects (for all rap/pT bins)"<<endl;
						for(int iRap = 0; iRap < NRQCDvars::nMaxRapBins; iRap++){

							int nPt_var=0;
							vector<double> v_pTmean;
							vector<double> v_model;
							vector<double> v_model_errlow;
							vector<double> v_model_errhigh;

						    for(int iP = 0; iP < NRQCDvars::nMaxPtBins; iP++){

								//get data input file
								sprintf(inname,"%s/ConvertedData_%s_%s_%s_rap%d_pT%d.txt",datadirname, StateName[iState],  MeasurementIDName[iMeasurementID],  ExpName[iExperiment], iRap, iP);
								ifstream in;
								in.open(inname);

								doesDataExist[iRap][iP]=false;
								//read NRQCDglobalfitObject members of data file
								if(in!=NULL){
									doesAnyDataExist_rap[iState][iMeasurementID][iExperiment][iRap]=true;
									nPt_var++;
									doesAnyDataExist=true;
									doesDataExist[iRap][iP]=true;
									DataObject[iRap][iP] = new NRQCDglobalfitObject();
									in >> *DataObject[iRap][iP];
									in.close();

									rapMin[iState][iMeasurementID][iExperiment][iRap]=DataObject[iRap][iP]->getyMin();
									rapMax[iState][iMeasurementID][iExperiment][iRap]=DataObject[iRap][iP]->getyMax();

									v_pTmean.push_back(DataObject[iRap][iP]->getpTMean());
									v_model.push_back(DataObject[iRap][iP]->getCentralValue());
									v_model_errlow.push_back(DataObject[iRap][iP]->getErrTotNeg());
									v_model_errhigh.push_back(DataObject[iRap][iP]->getErrTotPos());

								}


						    }


							doesAnyDataExist_all[iState][iMeasurementID][iExperiment]=doesAnyDataExist;


							if(doesAnyDataExist_rap[iState][iMeasurementID][iExperiment][iRap]){

								int nPt=nPt_var;
								double d_pTmean[nPt];
								double d_model[nPt];
								double d_model_errlow[nPt];
								double d_model_errhigh[nPt];
								double d_zero[nPt];

								for(int m=0;m<nPt;m++){
									d_pTmean[m]=v_pTmean[m];
									d_model[m]=v_model[m];
									d_model_errlow[m]=v_model_errlow[m];
									d_model_errhigh[m]=v_model_errhigh[m];
									d_zero[m]=0;
								}

								PtDataDists[iState][iMeasurementID][iExperiment][iRap] = new TGraphAsymmErrors(nPt,d_pTmean,d_model,d_zero,d_zero,d_model_errlow,d_model_errhigh);

							}


						}
						//cout<<"End reading in all DataObjects (for all rap/pT bins)"<<endl;

						if(!doesAnyDataExist) {cout<<"Measurement does not exist"<<endl; continue;}


			}//iExperiment
		}//iMeasurementID
	}//iState




	double fit_norm[NRQCDvars::nStates][NRQCDvars::nMeasurementIDs][NRQCDvars::nExperiments][NRQCDvars::nMaxRapBins];
	double fit_beta[NRQCDvars::nStates][NRQCDvars::nMeasurementIDs][NRQCDvars::nExperiments][NRQCDvars::nMaxRapBins];
	double fit_pT2[NRQCDvars::nStates][NRQCDvars::nMeasurementIDs][NRQCDvars::nExperiments][NRQCDvars::nMaxRapBins];
	double err_fit_norm[NRQCDvars::nStates][NRQCDvars::nMeasurementIDs][NRQCDvars::nExperiments][NRQCDvars::nMaxRapBins];
	double err_fit_beta[NRQCDvars::nStates][NRQCDvars::nMeasurementIDs][NRQCDvars::nExperiments][NRQCDvars::nMaxRapBins];
	double err_fit_pT2[NRQCDvars::nStates][NRQCDvars::nMeasurementIDs][NRQCDvars::nExperiments][NRQCDvars::nMaxRapBins];
	double fix_pTscale[NRQCDvars::nStates][NRQCDvars::nMeasurementIDs][NRQCDvars::nExperiments][NRQCDvars::nMaxRapBins];

	double red_chi2[NRQCDvars::nStates][NRQCDvars::nMeasurementIDs][NRQCDvars::nExperiments][NRQCDvars::nMaxRapBins];

	double ScaleFromMass=1;

	double pTmin=0;
	double pTmax=100;
	double sigmamin=1e-5;
	double sigmamax=1e3;

	char name[200];

	//FIT pT distributions
	for(int iState=0;iState<NRQCDvars::nStates;iState++){
		for(int iMeasurementID=0;iMeasurementID<NRQCDvars::nMeasurementIDs;iMeasurementID++){
			for(int iExperiment=0;iExperiment<NRQCDvars::nExperiments;iExperiment++){
				for(int iRap = 0; iRap < NRQCDvars::nMaxRapBins; iRap++){

					if(doesAnyDataExist_rap[iState][iMeasurementID][iExperiment][iRap]){


						sprintf(name, "fPtDist");
						TF1* fPtDist = new TF1(name, func_pT_gen, pTmin, pTmax, 6);

						double norm=1.;
						double beta=4.;
						double pT2=NRQCDvars::mass[iState]*NRQCDvars::mass[iState];
						double pTscale=1.;

						//double norm=1.;
						//double beta=3.;
						//double pT2=50.;
						//double pTscale=1.;
                        //
						//double norm=1.;
						//double beta=4.;
						//double pT2=100.;
						//double pTscale=1.;

						fPtDist->SetParameter(0,norm);
						fPtDist->SetParameter(1,beta);
						fPtDist->SetParameter(2,pT2);
						fPtDist->FixParameter(3,pTscale);

						//4...mass, 5...av_rap. In this configuration, x-axis is pT (p = pT):
						fPtDist->FixParameter(4,0);
						fPtDist->FixParameter(5,0);

						PtDataDists[iState][iMeasurementID][iExperiment][iRap]->Fit(fPtDist, "", "", pTmin, pTmax);


						fit_norm[iState][iMeasurementID][iExperiment][iRap]=fPtDist->GetParameter(0);
						fit_beta[iState][iMeasurementID][iExperiment][iRap]=fPtDist->GetParameter(1);
						fit_pT2[iState][iMeasurementID][iExperiment][iRap]=fPtDist->GetParameter(2);
						err_fit_norm[iState][iMeasurementID][iExperiment][iRap]=fPtDist->GetParError(0);
						err_fit_beta[iState][iMeasurementID][iExperiment][iRap]=fPtDist->GetParError(1);
						err_fit_pT2[iState][iMeasurementID][iExperiment][iRap]=fPtDist->GetParError(2);
						fix_pTscale[iState][iMeasurementID][iExperiment][iRap]=NRQCDvars::mass[iState]/ScaleFromMass;;

						red_chi2[iState][iMeasurementID][iExperiment][iRap]=fPtDist->GetChisquare()/fPtDist->GetNDF();

						PtDataDists[iState][iMeasurementID][iExperiment][iRap]->SetMarkerStyle(20);
						PtDataDists[iState][iMeasurementID][iExperiment][iRap]->SetMarkerSize(1.2);
						PtDataDists[iState][iMeasurementID][iExperiment][iRap]->SetMarkerColor(kRed);
						PtDataDists[iState][iMeasurementID][iExperiment][iRap]->SetLineColor(kBlack);

						fPtDist->SetLineColor(kBlack);
						fPtDist->SetLineWidth(2.);


						//TLegend* DMLegend=new TLegend(0.2,0.785,0.6,0.925);
						//DMLegend->SetFillColor(0);
						//DMLegend->SetTextFont(72);
						//DMLegend->SetTextSize(0.0245);
						//DMLegend->SetBorderSize(0);
						//DMLegend->SetMargin(0.135);
						//char DMLegendEntry[200];


							  char DMyTitle[200];
							  sprintf(DMyTitle,yTitle);
							  TCanvas* c2 = new TCanvas("c2","c2",1200,1100);
							  gStyle->SetPalette(1);
						 	  gPad->SetFillColor(kWhite);

						 	  double MargLeft=0.125;//0.175;
						 	  double MargRight=0.025;//0.1;
						 	  double MargTop=0.05;//0.175;

						 	  double MarkerSizeTom=1.75;

						 	  gPad->SetLeftMargin(MargLeft);
						      gPad->SetRightMargin(MargRight);
						      gPad->SetTopMargin(MargTop);


							  TH1F *histAxis = new TH1F;
							  histAxis->SetTitle(0);
							  histAxis->SetStats(0);
							  histAxis->SetMarkerStyle(25);
							  histAxis->SetMarkerSize(MarkerSizeTom);
							  histAxis->SetMarkerColor(kBlack);
							  histAxis->SetLineColor(kBlack);
							  histAxis->SetTitle(0);
							  histAxis->SetStats(0);
							  histAxis->SetMarkerColor(kRed);
							  histAxis->SetMarkerStyle(20);
							  histAxis->SetLineColor(kRed);
							  histAxis->SetMarkerSize(MarkerSizeTom);

							  histAxis = c2->DrawFrame(pTmin, sigmamin, pTmax, sigmamax);


							  histAxis->SetXTitle("dimuon p_{T} [GeV]");
							  histAxis->SetYTitle(DMyTitle);
							  histAxis->GetXaxis()->SetTitleOffset(1.2);
							  histAxis->GetYaxis()->SetTitleOffset(1.35);

							  histAxis->Draw("");

								//sprintf(DMLegendEntry,"2011 Parametrization");
								//DMLegend->AddEntry(fCONT,DMLegendEntry,"l");
								//sprintf(DMLegendEntry,"2012 data");
								//DMLegend->AddEntry(graphMassRap,DMLegendEntry,"lp");
								//sprintf(DMLegendEntry,"2012 Parametrization");
								//DMLegend->AddEntry(fParabola,DMLegendEntry,"l");
                                //
                                //
								//DMLegend->Draw("same");

							  PtDataDists[iState][iMeasurementID][iExperiment][iRap]->Draw("p,same");
							  fPtDist->Draw("l,same");

						      char text[200];
						      sprintf(text,"%s, %s, %1.1f<|y|<%1.1f", NRQCDvars::StateNameTex[iState], NRQCDvars::ExpNameTex[iExperiment], rapMin[iState][iMeasurementID][iExperiment][iRap], rapMax[iState][iMeasurementID][iExperiment][iRap]);
						      TLatex *tex = new TLatex(pTmin+(pTmax-pTmin)*0.35, sigmamin+(sigmamax-sigmamin)*0.2, text);
						      tex->SetTextSize(0.045)                                                                                                                                                                                                                                             ;
						      tex->Draw( "same" )                                                                                                                                                                                                                                                 ;

						      sprintf(text,"#sigma = Norm*p_{T}*[1+1/(beta-2)*p_{T}^{2}/<p_{T}^{2}>]^{-#beta}");
						      tex = new TLatex(pTmin+(pTmax-pTmin)*0.40, sigmamin+(sigmamax-sigmamin)*0.04, text);
						      tex->SetTextSize(0.030)                                                                                                                                                                                                                                             ;
						      tex->Draw( "same" )                                                                                                                                                                                                                                                 ;

						      sprintf(text,"Norm = %1.3f #pm %1.3f nb", fit_norm[iState][iMeasurementID][iExperiment][iRap], err_fit_norm[iState][iMeasurementID][iExperiment][iRap]);
						      tex = new TLatex(pTmin+(pTmax-pTmin)*0.40, sigmamin+(sigmamax-sigmamin)*0.015, text);
						      tex->SetTextSize(0.030)                                                                                                                                                                                                                                             ;
						      tex->Draw( "same" )                                                                                                                                                                                                                                                 ;

						      sprintf(text,"beta = %1.3f #pm %1.3f", fit_beta[iState][iMeasurementID][iExperiment][iRap], err_fit_beta[iState][iMeasurementID][iExperiment][iRap]);
						      tex = new TLatex(pTmin+(pTmax-pTmin)*0.40, sigmamin+(sigmamax-sigmamin)*0.005, text);
						      tex->SetTextSize(0.030)                                                                                                                                                                                                                                             ;
						      tex->Draw( "same" )                                                                                                                                                                                                                                                 ;

						      sprintf(text,"<p_{T}^{2}> = %1.3f #pm %1.3f GeV^{2}", fit_pT2[iState][iMeasurementID][iExperiment][iRap], err_fit_pT2[iState][iMeasurementID][iExperiment][iRap]);
						      tex = new TLatex(pTmin+(pTmax-pTmin)*0.40, sigmamin+(sigmamax-sigmamin)*0.0015, text);
						      tex->SetTextSize(0.030)                                                                                                                                                                                                                                             ;
						      tex->Draw( "same" )                                                                                                                                                                                                                                                 ;

						      sprintf(text,"#chi^{2}/ndf = %1.3f", red_chi2[iState][iMeasurementID][iExperiment][iRap]);
						      tex = new TLatex(pTmin+(pTmax-pTmin)*0.40, sigmamin+(sigmamax-sigmamin)*0.0005, text);
						      tex->SetTextSize(0.030)                                                                                                                                                                                                                                             ;
						      tex->Draw( "same" )                                                                                                                                                                                                                                                 ;



								char savename[200];
								c2->SetFrameBorderMode(0);
							  	sprintf(savename,"%s/FitPtDist_%s_%s_rap%d.pdf", datadirfigname, ExpName[iExperiment], StateName[iState], iRap);
							  	c2->SetLogy(true);
							  	c2->SaveAs(savename);


					}//TGraph!=0

				}//iRap
			}//iExperiment
		}//iMeasurementID
	}//iState










	double pTmin_scaled=3.4;
	double pTmax_scaled=10;
	double sigmamin_scaled=1e-6;
	double sigmamax_scaled=5e2;

	const int nCurves=5;
	int nStateCurves[nCurves]={0,3,4,7,10};
	int nColorCurves[nCurves]={418, 616, 800, 632, 1};
	int nMarkerStyleCurves[nCurves]={20, 21, 22, 23, 33};


	int iCurveRap=0;
	int iCurveExp=0;
	int iCurveMeasID=0;

	bool PlotCurves=true;
	if(PlotCurves){

		  char DMyTitle[200];
		  sprintf(DMyTitle,yTitle);
		  TCanvas* c2 = new TCanvas("c2","c2",1200,1100);
		  gStyle->SetPalette(1);
		  gPad->SetFillColor(kWhite);

		  double MargLeft=0.125;//0.175;
		  double MargRight=0.025;//0.1;
		  double MargTop=0.05;//0.175;

		  double MarkerSizeTom=1.75;

		gPad->SetLeftMargin(MargLeft);
		gPad->SetRightMargin(MargRight);
		gPad->SetTopMargin(MargTop);


		  TH1F *histAxis = new TH1F;
		  histAxis->SetTitle(0);
		  histAxis->SetStats(0);
		  histAxis->SetMarkerStyle(25);
		  histAxis->SetMarkerSize(MarkerSizeTom);
		  histAxis->SetMarkerColor(kBlack);
		  histAxis->SetLineColor(kBlack);
		  histAxis->SetTitle(0);
		  histAxis->SetStats(0);
		  histAxis->SetMarkerColor(kRed);
		  histAxis->SetMarkerStyle(20);
		  histAxis->SetLineColor(kRed);
		  histAxis->SetMarkerSize(MarkerSizeTom);

		  histAxis = c2->DrawFrame(pTmin_scaled, sigmamin_scaled, pTmax_scaled, sigmamax_scaled);


		  histAxis->SetXTitle(xTitle);
		  histAxis->SetYTitle(DMyTitle);
		  histAxis->GetXaxis()->SetTitleOffset(1.2);
		  histAxis->GetYaxis()->SetTitleOffset(1.35);

		  histAxis->Draw("");




			TLegend* DMLegend=new TLegend(0.65,0.685,0.9,0.925);
			DMLegend->SetFillColor(0);
			DMLegend->SetTextFont(72);
			DMLegend->SetTextSize(0.0245);
			DMLegend->SetBorderSize(0);
			DMLegend->SetMargin(0.135);
			char DMLegendEntry[200];



		sprintf(name, "fPtDist_scaled_1");


		TF1* fPtDist_scaled[nCurves];

		for(int iCurve=0;iCurve<nCurves;iCurve++){

			sprintf(name, "fPtDist_scaled_%d",iCurve);
			fPtDist_scaled[iCurve] = new TF1(name, func_pT_gen, pTmin_scaled, pTmax_scaled, 6);
			fPtDist_scaled[iCurve]->FixParameter(0,fit_norm[nStateCurves[iCurve]][iCurveMeasID][iCurveExp][iCurveRap]);
			fPtDist_scaled[iCurve]->FixParameter(1,fit_beta[nStateCurves[iCurve]][iCurveMeasID][iCurveExp][iCurveRap]);
			fPtDist_scaled[iCurve]->FixParameter(2,fit_pT2[nStateCurves[iCurve]][iCurveMeasID][iCurveExp][iCurveRap]);
			fPtDist_scaled[iCurve]->FixParameter(3,fix_pTscale[nStateCurves[iCurve]][iCurveMeasID][iCurveExp][iCurveRap]);

			fPtDist_scaled[iCurve]->FixParameter(4,NRQCDvars::mass[nStateCurves[iCurve]]);
			fPtDist_scaled[iCurve]->FixParameter(5,(rapMin[nStateCurves[iCurve]][iCurveMeasID][iCurveExp][iCurveRap]+rapMax[nStateCurves[iCurve]][iCurveMeasID][iCurveExp][iCurveRap])/2.);

			fPtDist_scaled[iCurve]->SetLineColor(nColorCurves[iCurve]);
			fPtDist_scaled[iCurve]->SetLineWidth(2.);

			fPtDist_scaled[iCurve]->Draw("l,same");

		    sprintf(DMLegendEntry,"%s, %s, %1.1f<|y|<%1.1f", NRQCDvars::StateNameTex[nStateCurves[iCurve]], NRQCDvars::ExpNameTex[iCurveExp], rapMin[nStateCurves[iCurve]][iCurveMeasID][iCurveExp][iCurveRap], rapMax[nStateCurves[iCurve]][iCurveMeasID][iCurveExp][iCurveRap]);
			DMLegend->AddEntry(fPtDist_scaled[iCurve],DMLegendEntry,"l");

		}

		DMLegend->Draw("same");




			char savename[200];
			c2->SetFrameBorderMode(0);
			sprintf(savename,"%s/FitPtDist_scaled_%s.pdf", datadirfigname, ExpName[iCurveExp]);
			c2->SetLogy(true);
			c2->SaveAs(savename);


	}


	double ratio_integral[NRQCDvars::nStates][NRQCDvars::nMeasurementIDs][NRQCDvars::nExperiments][NRQCDvars::nMaxRapBins];
	TGraph* NormScale_vs_mass;
	TGraphAsymmErrors* PtDataDists_scaled[nCurves];

	if(PlotCurves){

		int BuildRatioWithState=0;
		double ratiomin=0.3;
		double ratiomax=1.7;

		  char DMyTitle[200];
		  sprintf(DMyTitle,"d#sigma/(dp_{T}dy) / d#sigma(%s)/(dp_{T}dy) (norm.)",NRQCDvars::StateNameTex[BuildRatioWithState]);
		  TCanvas* c2 = new TCanvas("c2","c2",1200,1100);
		  gStyle->SetPalette(1);
		  gPad->SetFillColor(kWhite);

		  double MargLeft=0.125;//0.175;
		  double MargRight=0.025;//0.1;
		  double MargTop=0.05;//0.175;

		  double MarkerSizeTom=1.75;

		gPad->SetLeftMargin(MargLeft);
		gPad->SetRightMargin(MargRight);
		gPad->SetTopMargin(MargTop);


		  TH1F *histAxis = new TH1F;
		  histAxis->SetTitle(0);
		  histAxis->SetStats(0);
		  histAxis->SetMarkerStyle(25);
		  histAxis->SetMarkerSize(MarkerSizeTom);
		  histAxis->SetMarkerColor(kBlack);
		  histAxis->SetLineColor(kBlack);
		  histAxis->SetTitle(0);
		  histAxis->SetStats(0);
		  histAxis->SetMarkerColor(kRed);
		  histAxis->SetMarkerStyle(20);
		  histAxis->SetLineColor(kRed);
		  histAxis->SetMarkerSize(MarkerSizeTom);

		  histAxis = c2->DrawFrame(pTmin_scaled, ratiomin, pTmax_scaled, ratiomax);

		  histAxis->SetXTitle(xTitle);
		  histAxis->SetYTitle(DMyTitle);
		  histAxis->GetXaxis()->SetTitleOffset(1.2);
		  histAxis->GetYaxis()->SetTitleOffset(1.35);

		  histAxis->Draw("");




			TLegend* DMLegend=new TLegend(0.65,0.685,0.9,0.925);
			DMLegend->SetFillColor(0);
			DMLegend->SetTextFont(72);
			DMLegend->SetTextSize(0.0245);
			DMLegend->SetBorderSize(0);
			DMLegend->SetMargin(0.135);
			char DMLegendEntry[200];



		sprintf(name, "fPtDist_scaled_1");

		TF1* fPtDist_scaled_ratio[nCurves];

		TGraphAsymmErrors* PtDataDists_scaleduncertainty[nCurves];

		for(int iCurve=0;iCurve<nCurves;iCurve++){


			sprintf(name, "fPtDist_scaled_ratio%d",iCurve);
			fPtDist_scaled_ratio[iCurve] = new TF1(name, func_pT_gen_ratio, pTmin_scaled, pTmax_scaled, 12);
			fPtDist_scaled_ratio[iCurve]->FixParameter(0,fit_norm[nStateCurves[iCurve]][iCurveMeasID][iCurveExp][iCurveRap]);
			fPtDist_scaled_ratio[iCurve]->FixParameter(1,fit_beta[nStateCurves[iCurve]][iCurveMeasID][iCurveExp][iCurveRap]);
			fPtDist_scaled_ratio[iCurve]->FixParameter(2,fit_pT2[nStateCurves[iCurve]][iCurveMeasID][iCurveExp][iCurveRap]);
			fPtDist_scaled_ratio[iCurve]->FixParameter(3,fix_pTscale[nStateCurves[iCurve]][iCurveMeasID][iCurveExp][iCurveRap]);

			fPtDist_scaled_ratio[iCurve]->FixParameter(4,NRQCDvars::mass[nStateCurves[iCurve]]);
			fPtDist_scaled_ratio[iCurve]->FixParameter(5,(rapMin[nStateCurves[iCurve]][iCurveMeasID][iCurveExp][iCurveRap]+rapMax[nStateCurves[iCurve]][iCurveMeasID][iCurveExp][iCurveRap])/2.);

			fPtDist_scaled_ratio[iCurve]->FixParameter(6,fit_norm[BuildRatioWithState][iCurveMeasID][iCurveExp][iCurveRap]);
			fPtDist_scaled_ratio[iCurve]->FixParameter(7,fit_beta[BuildRatioWithState][iCurveMeasID][iCurveExp][iCurveRap]);
			fPtDist_scaled_ratio[iCurve]->FixParameter(8,fit_pT2[BuildRatioWithState][iCurveMeasID][iCurveExp][iCurveRap]);
			fPtDist_scaled_ratio[iCurve]->FixParameter(9,fix_pTscale[BuildRatioWithState][iCurveMeasID][iCurveExp][iCurveRap]);

			fPtDist_scaled_ratio[iCurve]->FixParameter(10,NRQCDvars::mass[BuildRatioWithState]);
			fPtDist_scaled_ratio[iCurve]->FixParameter(11,(rapMin[BuildRatioWithState][iCurveMeasID][iCurveExp][iCurveRap]+rapMax[BuildRatioWithState][iCurveMeasID][iCurveExp][iCurveRap])/2.);

			ratio_integral[nStateCurves[iCurve]][iCurveMeasID][iCurveExp][iCurveRap] = fPtDist_scaled_ratio[iCurve]->Integral(pTmin_scaled, pTmax_scaled)/(pTmax_scaled-pTmin_scaled);
			cout<<"Ratio norm from state "<<NRQCDvars::StateNameTex[nStateCurves[iCurve]]<<" = "<<ratio_integral[nStateCurves[iCurve]][iCurveMeasID][iCurveExp][iCurveRap]<<endl;
			fPtDist_scaled_ratio[iCurve]->FixParameter(0,fit_norm[nStateCurves[iCurve]][iCurveMeasID][iCurveExp][iCurveRap]/ratio_integral[nStateCurves[iCurve]][iCurveMeasID][iCurveExp][iCurveRap]);

			//if(nStateCurves[iCurve]==BuildRatioWithState) continue;

			fPtDist_scaled_ratio[iCurve]->SetLineColor(nColorCurves[iCurve]);
			fPtDist_scaled_ratio[iCurve]->SetLineWidth(2.);
			fPtDist_scaled_ratio[iCurve]->Draw("l,same");

		    sprintf(DMLegendEntry,"%s, %s, %1.1f<|y|<%1.1f", NRQCDvars::StateNameTex[nStateCurves[iCurve]], NRQCDvars::ExpNameTex[iCurveExp], rapMin[nStateCurves[iCurve]][iCurveMeasID][iCurveExp][iCurveRap], rapMax[nStateCurves[iCurve]][iCurveMeasID][iCurveExp][iCurveRap]);
			DMLegend->AddEntry(fPtDist_scaled_ratio[iCurve],DMLegendEntry,"l");

			//add data uncertainties to ratio:

			int nPt=PtDataDists[nStateCurves[iCurve]][iCurveMeasID][iCurveExp][iCurveRap]->GetN();
			double d_pTmean[nPt];
			double d_model[nPt];
			double d_model_errlow[nPt];
			double d_model_errhigh[nPt];
			double d_zero[nPt];

			for(int m=0;m<nPt;m++){


				double pTscale=fix_pTscale[nStateCurves[iCurve]][iCurveMeasID][iCurveExp][iCurveRap];
				double mass=NRQCDvars::mass[nStateCurves[iCurve]];
				double av_rap=(rapMin[nStateCurves[iCurve]][iCurveMeasID][iCurveExp][iCurveRap]+rapMax[nStateCurves[iCurve]][iCurveMeasID][iCurveExp][iCurveRap])/2.;
				double pTmean_graph, model_graph;
				PtDataDists[nStateCurves[iCurve]][iCurveMeasID][iCurveExp][iCurveRap]->GetPoint(m,pTmean_graph, model_graph);

				if(scale_pT_over_m)
					d_pTmean[m]=pTmean_graph/pTscale;
				if(scale_p_over_m)
					d_pTmean[m]=TMath::Sqrt((pTmean_graph*pTmean_graph+mass*mass)*TMath::SinH(av_rap)*TMath::SinH(av_rap)+pTmean_graph*pTmean_graph)/pTscale;
				if(scale_mT_over_m)
					d_pTmean[m]=TMath::Sqrt((pTmean_graph*pTmean_graph+mass*mass))/pTscale;


				d_model[m]=fPtDist_scaled_ratio[iCurve]->Eval(d_pTmean[m]);
				d_model_errlow[m]=PtDataDists[nStateCurves[iCurve]][iCurveMeasID][iCurveExp][iCurveRap]->GetErrorYlow(m)*d_model[m]/model_graph;
				d_model_errhigh[m]=PtDataDists[nStateCurves[iCurve]][iCurveMeasID][iCurveExp][iCurveRap]->GetErrorYhigh(m)*d_model[m]/model_graph;
				d_zero[m]=0;
				//cout<<"GetErrorYlow "<<PtDataDists[nStateCurves[iCurve]][iCurveMeasID][iCurveExp][iCurveRap]->GetErrorYlow(m)<<endl;
				//cout<<"GetErrorYlhigh "<<PtDataDists[nStateCurves[iCurve]][iCurveMeasID][iCurveExp][iCurveRap]->GetErrorYhigh(m)<<endl;
				//cout<<"model_graph"<<model_graph<<endl;
				//cout<<"pTmean_graph "<<pTmean_graph<<endl;
			}

			PtDataDists_scaleduncertainty[iCurve] = new TGraphAsymmErrors(nPt,d_pTmean,d_model,d_zero,d_zero,d_model_errlow,d_model_errhigh);

			PtDataDists_scaleduncertainty[iCurve]->SetMarkerStyle(1);
			PtDataDists_scaleduncertainty[iCurve]->SetMarkerColor(nColorCurves[iCurve]);
			PtDataDists_scaleduncertainty[iCurve]->SetLineColor(nColorCurves[iCurve]);

			if(nStateCurves[iCurve]!=BuildRatioWithState) PtDataDists_scaleduncertainty[iCurve]->Draw("p,same");

			for(int m=0;m<nPt;m++){


				double pTscale=fix_pTscale[nStateCurves[iCurve]][iCurveMeasID][iCurveExp][iCurveRap];
				double mass=NRQCDvars::mass[nStateCurves[iCurve]];
				double av_rap=(rapMin[nStateCurves[iCurve]][iCurveMeasID][iCurveExp][iCurveRap]+rapMax[nStateCurves[iCurve]][iCurveMeasID][iCurveExp][iCurveRap])/2.;
				double pTmean_graph, model_graph;
				PtDataDists[nStateCurves[iCurve]][iCurveMeasID][iCurveExp][iCurveRap]->GetPoint(m,pTmean_graph, model_graph);

				if(scale_pT_over_m)
					d_pTmean[m]=pTmean_graph/pTscale;
				if(scale_p_over_m)
					d_pTmean[m]=TMath::Sqrt((pTmean_graph*pTmean_graph+mass*mass)*TMath::SinH(av_rap)*TMath::SinH(av_rap)+pTmean_graph*pTmean_graph)/pTscale;
				if(scale_mT_over_m)
					d_pTmean[m]=TMath::Sqrt((pTmean_graph*pTmean_graph+mass*mass))/pTscale;


				d_model[m]=model_graph;
				d_model_errlow[m]=PtDataDists[nStateCurves[iCurve]][iCurveMeasID][iCurveExp][iCurveRap]->GetErrorYlow(m);
				d_model_errhigh[m]=PtDataDists[nStateCurves[iCurve]][iCurveMeasID][iCurveExp][iCurveRap]->GetErrorYhigh(m);
				d_zero[m]=0;
			}

			PtDataDists_scaled[iCurve] = new TGraphAsymmErrors(nPt,d_pTmean,d_model,d_zero,d_zero,d_model_errlow,d_model_errhigh);

			cout<<"iCurve "<<iCurve<<endl;
			cout<<"PtDataDists_scaled[iCurve]->Print();"<<endl;
			PtDataDists_scaled[iCurve]->Print();

			PtDataDists_scaled[iCurve]->SetMarkerSize(1.5);
			PtDataDists_scaled[iCurve]->SetMarkerStyle(nMarkerStyleCurves[iCurve]);
			PtDataDists_scaled[iCurve]->SetMarkerColor(nColorCurves[iCurve]);
			PtDataDists_scaled[iCurve]->SetLineColor(nColorCurves[iCurve]);

		}



		DMLegend->Draw("same");




			char savename[200];
			c2->SetFrameBorderMode(0);
			sprintf(savename,"%s/FitPtDist_scaled_ratio_%s.pdf", datadirfigname, ExpName[iCurveExp]);
			//c2->SetLogy(true);
			c2->SaveAs(savename);


			double d_mass[nCurves];
			double d_norm[nCurves];

			for(int iCurve=0;iCurve<nCurves;iCurve++){

				d_mass[iCurve]=NRQCDvars::mass[nStateCurves[iCurve]];
				d_norm[iCurve]=ratio_integral[nStateCurves[iCurve]][iCurveMeasID][iCurveExp][iCurveRap];


			}
			NormScale_vs_mass = new TGraph(nCurves,d_mass,d_norm);


	}





	if(PlotCurves){

		int BuildRatioWithState=0;
		double ratiomin=1e-4;
		double ratiomax=5;
		double massmin=1;
		double massmax=13;

		  char DMyTitle[200];
		  sprintf(DMyTitle,"Normalization scale",NRQCDvars::StateNameTex[BuildRatioWithState]);
		  TCanvas* c2 = new TCanvas("c2","c2",1200,1100);
		  gStyle->SetPalette(1);
		  gPad->SetFillColor(kWhite);

		  double MargLeft=0.125;//0.175;
		  double MargRight=0.025;//0.1;
		  double MargTop=0.05;//0.175;

		  double MarkerSizeTom=1.75;

		gPad->SetLeftMargin(MargLeft);
		gPad->SetRightMargin(MargRight);
		gPad->SetTopMargin(MargTop);


		  TH1F *histAxis = new TH1F;
		  histAxis->SetTitle(0);
		  histAxis->SetStats(0);
		  histAxis->SetMarkerStyle(25);
		  histAxis->SetMarkerSize(MarkerSizeTom);
		  histAxis->SetMarkerColor(kBlack);
		  histAxis->SetLineColor(kBlack);
		  histAxis->SetTitle(0);
		  histAxis->SetStats(0);
		  histAxis->SetMarkerColor(kRed);
		  histAxis->SetMarkerStyle(20);
		  histAxis->SetLineColor(kRed);
		  histAxis->SetMarkerSize(MarkerSizeTom);

		  histAxis = c2->DrawFrame(massmin, ratiomin, massmax, ratiomax);

		  histAxis->SetXTitle("dimuon m [GeV]");
		  histAxis->SetYTitle(DMyTitle);
		  histAxis->GetXaxis()->SetTitleOffset(1.2);
		  histAxis->GetYaxis()->SetTitleOffset(1.35);

		  histAxis->Draw("");



			NormScale_vs_mass->SetMarkerStyle(20);
			NormScale_vs_mass->SetMarkerSize(2.);
			NormScale_vs_mass->SetMarkerColor(kRed);
			NormScale_vs_mass->SetLineColor(kBlack);


			NormScale_vs_mass->Draw("p,same");

		//	TLegend* DMLegend=new TLegend(0.65,0.685,0.9,0.925);
		//	DMLegend->SetFillColor(0);
		//	DMLegend->SetTextFont(72);
		//	DMLegend->SetTextSize(0.0245);
		//	DMLegend->SetBorderSize(0);
		//	DMLegend->SetMargin(0.135);
		//	char DMLegendEntry[200];
        //


			int point1=3;
			int point2=10;
			double exponent=log(ratio_integral[point1][iCurveMeasID][iCurveExp][iCurveRap]/ratio_integral[point2][iCurveMeasID][iCurveExp][iCurveRap])/(NRQCDvars::mass[point1]-NRQCDvars::mass[point2]);
			double expnorm=ratio_integral[point1][iCurveMeasID][iCurveExp][iCurveRap]/exp(exponent*NRQCDvars::mass[point1]);
		sprintf(name, "expo_fit");
		TF1* expo_fit  = new TF1(name, exponential, massmin, massmax, 2);
		expo_fit->FixParameter(0,expnorm);
		expo_fit->FixParameter(1,exponent);

		cout<<"expnorm = "<<expnorm<<endl;
		cout<<"exponent = "<<exponent<<endl;

		expo_fit->SetLineColor(kBlack);
		expo_fit->SetLineWidth(1.5);
		expo_fit->SetLineStyle(2);
		expo_fit->Draw("l,same");


		double directFraction[NRQCDvars::nStates];

		directFraction[0] = expo_fit->Eval(NRQCDvars::mass[0])/ratio_integral[0][iCurveMeasID][iCurveExp][iCurveRap];
		directFraction[4] = expo_fit->Eval(NRQCDvars::mass[4])/ratio_integral[4][iCurveMeasID][iCurveExp][iCurveRap];
		directFraction[7] = expo_fit->Eval(NRQCDvars::mass[7])/ratio_integral[7][iCurveMeasID][iCurveExp][iCurveRap];

		cout<<"fraction of directly produced quarkonia: "<<endl;
		cout<<"Psi(1S): "<<directFraction[0]<<endl;
		cout<<"Ups(1S): "<<directFraction[4]<<endl;
		cout<<"Ups(2S): "<<directFraction[7]<<endl;


		//
        //
        //
	    //sprintf(DMLegendEntry,"%s, %s, %1.1f<|y|<%1.1f", NRQCDvars::StateNameTex[nStateCurves[iCurve]], NRQCDvars::ExpName[iCurveExp], rapMin[nStateCurves[iCurve]][iCurveMeasID][iCurveExp][iCurveRap], rapMax[nStateCurves[iCurve]][iCurveMeasID][iCurveExp][iCurveRap]);
		//DMLegend->AddEntry(fPtDist_scaled_ratio[iCurve],DMLegendEntry,"l");
        //
		//DMLegend->Draw("same");
        //
        //


			char savename[200];
			c2->SetFrameBorderMode(0);
			sprintf(savename,"%s/NormScale_vs_mass_%s.pdf", datadirfigname, ExpName[iCurveExp]);
			c2->SetLogy(true);
			c2->SaveAs(savename);



	}



	if(PlotCurves){






	double scaled_fit_norm[NRQCDvars::nStates][NRQCDvars::nMeasurementIDs][NRQCDvars::nExperiments][NRQCDvars::nMaxRapBins];
	double scaled_fit_beta[NRQCDvars::nStates][NRQCDvars::nMeasurementIDs][NRQCDvars::nExperiments][NRQCDvars::nMaxRapBins];
	double scaled_fit_pT2[NRQCDvars::nStates][NRQCDvars::nMeasurementIDs][NRQCDvars::nExperiments][NRQCDvars::nMaxRapBins];
	double scaled_err_fit_norm[NRQCDvars::nStates][NRQCDvars::nMeasurementIDs][NRQCDvars::nExperiments][NRQCDvars::nMaxRapBins];
	double scaled_err_fit_beta[NRQCDvars::nStates][NRQCDvars::nMeasurementIDs][NRQCDvars::nExperiments][NRQCDvars::nMaxRapBins];
	double scaled_err_fit_pT2[NRQCDvars::nStates][NRQCDvars::nMeasurementIDs][NRQCDvars::nExperiments][NRQCDvars::nMaxRapBins];


	//FIT scaled pT distributions

	cout<<"FIT SCALED PT CURVES"<<endl;

	for(int iCurve=0;iCurve<nCurves;iCurve++){


						sprintf(name, "fPtDist_scaled_2");
						TF1* fPtDist_ = new TF1(name, func_pT_gen, pTmin_scaled, pTmax_scaled, 6);

						double norm=1.;
						double beta=3.;
						double pT2=1.;
						double pTscale=1.;

						fPtDist_->SetParameter(0,norm);
						fPtDist_->SetParameter(1,beta);
						fPtDist_->SetParameter(2,pT2);
						fPtDist_->FixParameter(3,pTscale);

						//4...mass, 5...av_rap. In this configuration, x-axis is pT (p = pT):
						fPtDist_->FixParameter(4,0);
						fPtDist_->FixParameter(5,0);

						PtDataDists_scaled[iCurve]->Fit(fPtDist_, "", "", pTmin_scaled, pTmax_scaled);
						//PtDataDists_scaled[iCurve]->Fit(fPtDist_, "", "", 1, 10);


						scaled_fit_norm[nStateCurves[iCurve]][iCurveMeasID][iCurveExp][iCurveRap]=fPtDist_->GetParameter(0);
						scaled_fit_beta[nStateCurves[iCurve]][iCurveMeasID][iCurveExp][iCurveRap]=fPtDist_->GetParameter(1);
						scaled_fit_pT2[nStateCurves[iCurve]][iCurveMeasID][iCurveExp][iCurveRap]=fPtDist_->GetParameter(2);
						scaled_err_fit_norm[nStateCurves[iCurve]][iCurveMeasID][iCurveExp][iCurveRap]=fPtDist_->GetParError(0);
						scaled_err_fit_beta[nStateCurves[iCurve]][iCurveMeasID][iCurveExp][iCurveRap]=fPtDist_->GetParError(1);
						scaled_err_fit_pT2[nStateCurves[iCurve]][iCurveMeasID][iCurveExp][iCurveRap]=fPtDist_->GetParError(2);



	}









	const int meanmodel_nStates=nCurves;
	double d_meanmodel_x[meanmodel_nStates];
	double d_meanmodel_beta[meanmodel_nStates];
	double d_meanmodel_beta_err[meanmodel_nStates];
	double d_meanmodel_pT2[meanmodel_nStates];
	double d_meanmodel_pT2_err[meanmodel_nStates];
	double d_meanmodel_zero[meanmodel_nStates];
	int meanmodelState[meanmodel_nStates]={0,3,4,7,10};

	for(int m=0;m<meanmodel_nStates;m++){
		d_meanmodel_x[m]=double(m);
		d_meanmodel_beta[m]=scaled_fit_beta[nStateCurves[m]][0][0][0];
		d_meanmodel_beta_err[m]=scaled_err_fit_beta[nStateCurves[m]][0][0][0];
		d_meanmodel_pT2[m]=scaled_fit_pT2[nStateCurves[m]][0][0][0];
		d_meanmodel_pT2_err[m]=scaled_err_fit_pT2[nStateCurves[m]][0][0][0];
		d_meanmodel_zero[m]=0;
	}

	TGraphAsymmErrors* meanmodel_beta = new TGraphAsymmErrors(meanmodel_nStates,d_meanmodel_x,d_meanmodel_beta,d_meanmodel_zero,d_meanmodel_zero,d_meanmodel_beta_err,d_meanmodel_beta_err);
	TGraphAsymmErrors* meanmodel_pT2 = new TGraphAsymmErrors(meanmodel_nStates,d_meanmodel_x,d_meanmodel_pT2,d_meanmodel_zero,d_meanmodel_zero,d_meanmodel_pT2_err,d_meanmodel_pT2_err);

	TGraphErrors* meanmodel_beta_pT2 =  new TGraphErrors(meanmodel_nStates,d_meanmodel_beta, d_meanmodel_pT2,d_meanmodel_beta_err,d_meanmodel_pT2_err);

    char FitOptions[200];
    sprintf(FitOptions,"EFNR");

	cout<<"meanmodel_beta->Print();"<<endl;
    meanmodel_beta->Print();
	cout<<"meanmodel_pT2->Print();"<<endl;
	meanmodel_pT2->Print();

    TF1* fConst = new TF1("fConst","pol0",-1,meanmodel_nStates+1);
    meanmodel_beta->Fit("fConst",FitOptions);
    double Const_beta = fConst->GetParameter(0);
    meanmodel_pT2->Fit("fConst",FitOptions);
    double Const_pT2 = fConst->GetParameter(0);


	cout<<"Const_beta "<<Const_beta<<endl;
	cout<<"Const_pT2 "<<Const_pT2<<endl;


		int BuildRatioWithState=0;
		double ratiomin=0.3;
		double ratiomax=1.7;


		for(int iPlot=0;iPlot<2;iPlot++){

		  char DMyTitle[200];
		  sprintf(DMyTitle,"d#sigma/(dp_{T}dy) [nb/GeV]");
		  TCanvas* c2 = new TCanvas("c2","c2",1200,1100);
		  gStyle->SetPalette(1);
		  gPad->SetFillColor(kWhite);

		  double MargLeft=0.125;//0.175;
		  double MargRight=0.025;//0.1;
		  double MargTop=0.05;//0.175;

		  double MarkerSizeTom=1.75;

		gPad->SetLeftMargin(MargLeft);
		gPad->SetRightMargin(MargRight);
		gPad->SetTopMargin(MargTop);


		  TH1F *histAxis = new TH1F;
		  histAxis->SetTitle(0);
		  histAxis->SetStats(0);
		  histAxis->SetMarkerStyle(25);
		  histAxis->SetMarkerSize(MarkerSizeTom);
		  histAxis->SetMarkerColor(kBlack);
		  histAxis->SetLineColor(kBlack);
		  histAxis->SetTitle(0);
		  histAxis->SetStats(0);
		  histAxis->SetMarkerColor(kRed);
		  histAxis->SetMarkerStyle(20);
		  histAxis->SetLineColor(kRed);
		  histAxis->SetMarkerSize(MarkerSizeTom);

		  histAxis = c2->DrawFrame(pTmin_scaled, sigmamin_scaled, pTmax_scaled, sigmamax_scaled);

		  histAxis->SetXTitle(xTitle);
		  histAxis->SetYTitle(DMyTitle);
		  histAxis->GetXaxis()->SetTitleOffset(1.2);
		  histAxis->GetYaxis()->SetTitleOffset(1.35);

		  histAxis->Draw("");




			TLegend* DMLegend=new TLegend(0.65,0.685,0.9,0.925);
			DMLegend->SetFillColor(0);
			DMLegend->SetTextFont(72);
			DMLegend->SetTextSize(0.0245);
			DMLegend->SetBorderSize(0);
			DMLegend->SetMargin(0.135);
			char DMLegendEntry[200];



		sprintf(name, "fPtDist_scaled_1");

		TF1* fPtDist_scaled[nCurves];

		cout<<"FIT NORM OF MEAN-SCALED PT CURVES, PLOT"<<endl;

		for(int iCurve=0;iCurve<nCurves;iCurve++){

			sprintf(name, "fPtDist_scaled%d",iCurve);
			fPtDist_scaled[iCurve] = new TF1(name, func_pT_gen, pTmin_scaled, pTmax_scaled, 6);

			double norm=1.;
			double beta=4.;
			double pT2=NRQCDvars::mass[nStateCurves[iCurve]]*NRQCDvars::mass[nStateCurves[iCurve]];
			double pTscale=1.;

			if(iPlot==0){
			//if you want to plot the 'mean' curve
			fPtDist_scaled[iCurve]->FixParameter(1,Const_beta);
			fPtDist_scaled[iCurve]->FixParameter(2,Const_pT2);
			fPtDist_scaled[iCurve]->FixParameter(3,1.);
			}
			if(iPlot==1){
			//if you want to plot the individually fitted curves
			fPtDist_scaled[iCurve]->SetParameter(0,norm);
			fPtDist_scaled[iCurve]->FixParameter(1,scaled_fit_beta[nStateCurves[iCurve]][0][0][0]);
			fPtDist_scaled[iCurve]->FixParameter(2,scaled_fit_pT2[nStateCurves[iCurve]][0][0][0]);
			fPtDist_scaled[iCurve]->FixParameter(3,1.);
			}


			fPtDist_scaled[iCurve]->FixParameter(4,NRQCDvars::mass[nStateCurves[iCurve]]);
			fPtDist_scaled[iCurve]->FixParameter(5,(rapMin[nStateCurves[iCurve]][iCurveMeasID][iCurveExp][iCurveRap]+rapMax[nStateCurves[iCurve]][iCurveMeasID][iCurveExp][iCurveRap])/2.);

			fPtDist_scaled[iCurve]->SetLineWidth(0.);
			PtDataDists_scaled[iCurve]->Fit(fPtDist_scaled[iCurve], "", "", pTmin_scaled, pTmax_scaled);

			fPtDist_scaled[iCurve]->SetLineColor(nColorCurves[iCurve]);
			fPtDist_scaled[iCurve]->SetLineWidth(2.);
			PtDataDists_scaled[iCurve]->Draw("p,same");
			fPtDist_scaled[iCurve]->Draw("l,same");

		    sprintf(DMLegendEntry,"%s, %s, %1.1f<|y|<%1.1f", NRQCDvars::StateNameTex[nStateCurves[iCurve]], NRQCDvars::ExpNameTex[iCurveExp], rapMin[nStateCurves[iCurve]][iCurveMeasID][iCurveExp][iCurveRap], rapMax[nStateCurves[iCurve]][iCurveMeasID][iCurveExp][iCurveRap]);
			DMLegend->AddEntry(PtDataDists_scaled[iCurve],DMLegendEntry,"lp");

		}


		DMLegend->Draw("same");




			char savename[200];
			c2->SetFrameBorderMode(0);
			if(iPlot==0) sprintf(savename,"%s/FitPtDist_scaledAllSame_withdata_%s.pdf", datadirfigname, ExpName[iCurveExp]);
			if(iPlot==1) sprintf(savename,"%s/FitPtDist_scaledIndividualFits_withdata_%s.pdf", datadirfigname, ExpName[iCurveExp]);
			c2->SetLogy(true);
			c2->SaveAs(savename);



			//// Make Pulls:::

			double pullmin, pullmax;

			pullmin=-5;
			pullmax=5;

			int nPullLines=(pullmax-pullmin-2)/2;

			  sprintf(DMyTitle,"Standard scores #frac{data-fit}{#sigma_{data}}");
			  TCanvas* c2pull = new TCanvas("c2pull","c2pull",1200,500);
			  gStyle->SetPalette(1);
			  gPad->SetFillColor(kWhite);

			  MargLeft=0.125;//0.175;
			  MargRight=0.025;//0.1;
			  MargTop=0.05;//0.175;

			  MarkerSizeTom=1.75;

			gPad->SetLeftMargin(MargLeft);
			gPad->SetRightMargin(MargRight);
			gPad->SetTopMargin(MargTop);


			  TH1F *histAxispull = new TH1F;
			  histAxispull->SetTitle(0);
			  histAxispull->SetStats(0);
			  histAxispull->SetMarkerStyle(25);
			  histAxispull->SetMarkerSize(MarkerSizeTom);
			  histAxispull->SetMarkerColor(kBlack);
			  histAxispull->SetLineColor(kBlack);
			  histAxispull->SetTitle(0);
			  histAxispull->SetStats(0);
			  histAxispull->SetMarkerColor(kRed);
			  histAxispull->SetMarkerStyle(20);
			  histAxispull->SetLineColor(kRed);
			  histAxispull->SetMarkerSize(MarkerSizeTom);

			  histAxispull = c2pull->DrawFrame(pTmin_scaled, pullmin, pTmax_scaled, pullmax);

			  histAxispull->SetXTitle(xTitle);
			  histAxispull->SetYTitle(DMyTitle);
			  histAxispull->GetXaxis()->SetTitleOffset(1.2);
			  histAxispull->GetYaxis()->SetTitleOffset(1.35);

			  histAxispull->Draw("");




				TLine* extreme0 = new TLine( pTmin_scaled, 0, pTmax_scaled, 0);
				extreme0->SetLineWidth( 1 );
				extreme0->SetLineStyle( 1 );
				extreme0->SetLineColor( kBlack );
				extreme0->Draw( "same" );

				for(int iLine=1;iLine<nPullLines+1;iLine++){

					TLine* extremeP = new TLine( pTmin_scaled, iLine, pTmax_scaled, iLine);
					extremeP->SetLineWidth( 1 );
					extremeP->SetLineStyle( 2 );
					extremeP->SetLineColor( kBlack );
					extremeP->Draw( "same" );

					TLine* extremeN = new TLine( pTmin_scaled, -iLine, pTmax_scaled, -iLine);
					extremeN->SetLineWidth( 1 );
					extremeN->SetLineStyle( 2 );
					extremeN->SetLineColor( kBlack );
					extremeN->Draw( "same" );

			  }

			  TGraph* PtPullDists_scaled[nCurves];

				for(int iCurve=0;iCurve<nCurves;iCurve++){

					int nPt=PtDataDists_scaled[iCurve]->GetN();
					double d_pTmean[nPt];
					double d_model[nPt];
					double d_model_errlow[nPt];
					double d_model_errhigh[nPt];
					double d_model_err[nPt];
					double d_zero[nPt];

					double pTmean_graph;
					double model_graph;

					for(int m=0;m<nPt;m++){


						PtDataDists_scaled[iCurve]->GetPoint(m,pTmean_graph, model_graph);
						d_model_errlow[m]=PtDataDists_scaled[iCurve]->GetErrorYlow(m);
						d_model_errhigh[m]=PtDataDists_scaled[iCurve]->GetErrorYhigh(m);
						d_model_err[m]=(d_model_errlow[m]+d_model_errhigh[m])/2.;

						d_pTmean[m]=pTmean_graph;
						d_model[m]=(model_graph-fPtDist_scaled[iCurve]->Eval(pTmean_graph))/d_model_err[m];
						d_zero[m]=0;

					}

					PtPullDists_scaled[iCurve] = new TGraph(nPt,d_pTmean,d_model);

					PtPullDists_scaled[iCurve]->SetMarkerSize(1.5);
					PtPullDists_scaled[iCurve]->SetMarkerStyle(nMarkerStyleCurves[iCurve]);
					PtPullDists_scaled[iCurve]->SetMarkerColor(nColorCurves[iCurve]);
					PtPullDists_scaled[iCurve]->SetLineColor(nColorCurves[iCurve]);

					PtPullDists_scaled[iCurve]->Draw("psame");


				}



				c2pull->SetFrameBorderMode(0);
				if(iPlot==0) sprintf(savename,"%s/FitPtDist_scaledAllSame_pull_withdata_%s.pdf", datadirfigname, ExpName[iCurveExp]);
				if(iPlot==1) sprintf(savename,"%s/FitPtDist_scaledIndividualFits_pull_withdata_%s.pdf", datadirfigname, ExpName[iCurveExp]);
				c2pull->SetLogy(false);
				c2pull->SaveAs(savename);





				//// Make Pulls:::

				double ellmin, ellmax;

				double betamin, pT2min, betamax, pT2max;
				betamin=2;
				betamax=6;
				pT2min=-2;
				pT2max=5;

				char ell_xTitle[200];
				sprintf(ell_xTitle, "#beta");

				  sprintf(DMyTitle,"#gamma");
				  TCanvas* c2ell = new TCanvas("c2ell","c2ell",1200,1100);
				  gStyle->SetPalette(1);
				  gPad->SetFillColor(kWhite);

				  MargLeft=0.125;//0.175;
				  MargRight=0.025;//0.1;
				  MargTop=0.05;//0.175;

				  MarkerSizeTom=1.75;

				gPad->SetLeftMargin(MargLeft);
				gPad->SetRightMargin(MargRight);
				gPad->SetTopMargin(MargTop);


				  TH1F *histAxisell = new TH1F;
				  histAxisell->SetTitle(0);
				  histAxisell->SetStats(0);

				  histAxisell = c2ell->DrawFrame(betamin, pT2min, betamax, pT2max);

				  histAxisell->SetXTitle(ell_xTitle);
				  histAxisell->SetYTitle(DMyTitle);
				  histAxisell->GetXaxis()->SetTitleOffset(1.2);
				  histAxisell->GetYaxis()->SetTitleOffset(1.35);

				  histAxisell->Draw("AXIS");

				  DMLegend->Draw("same");



					for(int iCurve=0;iCurve<nCurves;iCurve++){

						double buffx, buffy, bufferrx, bufferry;

						meanmodel_beta_pT2->GetPoint(iCurve, buffx, buffy);
						bufferrx=meanmodel_beta_pT2->GetErrorX(iCurve);
						bufferry=meanmodel_beta_pT2->GetErrorY(iCurve);

						TEllipse* ellipse = new TEllipse(buffx, buffy, bufferrx, bufferry);
						TMarker* point = new TMarker(buffx, buffy, nMarkerStyleCurves[iCurve]);

						point->SetMarkerSize(1.5);
						point->SetMarkerStyle(nMarkerStyleCurves[iCurve]);
						point->SetMarkerColor(nColorCurves[iCurve]);

						ellipse->SetLineColor(nColorCurves[iCurve]);
						ellipse->SetLineStyle(1);
						ellipse->SetFillStyle(0);
						ellipse->SetLineWidth(1.5);

						ellipse->Draw("lsame");
						point->Draw("psame");


					}



					c2ell->SetFrameBorderMode(0);
					sprintf(savename,"%s/FitPtDist_beta_gamma_ellipse_%s.pdf", datadirfigname, ExpName[iCurveExp]);
					c2ell->SetLogy(false);
					c2ell->SaveAs(savename);


		}

	}







return 0;

}


double func_pT_gen(double* x, double* par) {


	double funcval;

	double norm=par[0];
	double beta=par[1];
	double pT2=par[2];
	double pTscale=par[3];

	double mass=par[4];
	double av_rap=par[5];

	double pTFit;


	if(scale_p_over_m){
		double pFit=x[0]*pTscale;
		pTFit=TMath::Sqrt((pFit*pFit-mass*mass*TMath::SinH(av_rap)*TMath::SinH(av_rap))/(1+TMath::SinH(av_rap)*TMath::SinH(av_rap)));
	}
	if(scale_mT_over_m){
		pTFit=TMath::Sqrt(x[0]*x[0]-1)*pTscale;
	}
	if(scale_pT_over_m || mass<1e-5){
		pTFit=x[0]*pTscale;
	}

	funcval = norm*pTFit*pow( 1+1/(beta-2)*pTFit*pTFit/(pT2), -beta );


	return funcval;

}





double func_pT_gen_ratio(double* x, double* par) {


	double funcval;

	double funcval1;
	double norm1=par[0];
	double beta1=par[1];
	double pT21=par[2];
	double pTscale1=par[3];

	double mass1=par[4];
	double av_rap1=par[5];

	double pTFit1;



	if(scale_p_over_m){
		double pFit=x[0]*pTscale1;
		pTFit1=TMath::Sqrt((pFit*pFit-mass1*mass1*TMath::SinH(av_rap1)*TMath::SinH(av_rap1))/(1+TMath::SinH(av_rap1)*TMath::SinH(av_rap1)));
	}
	if(scale_mT_over_m){
		pTFit1=TMath::Sqrt(x[0]*x[0]-1)*pTscale1;
	}
	if(scale_pT_over_m || mass1<1e-5){
		pTFit1=x[0]*pTscale1;
	}

	funcval1 = norm1*pTFit1*pow( 1+1/(beta1-2)*pTFit1*pTFit1/(pT21), -beta1 );



	double funcval2;
	double norm2=par[6];
	double beta2=par[7];
	double pT22=par[8];
	double pTscale2=par[9];

	double mass2=par[10];
	double av_rap2=par[11];

	double pTFit2;

	if(scale_p_over_m){
		double pFit=x[0]*pTscale2;
		pTFit2=TMath::Sqrt((pFit*pFit-mass2*mass2*TMath::SinH(av_rap2)*TMath::SinH(av_rap2))/(1+TMath::SinH(av_rap2)*TMath::SinH(av_rap2)));
	}
	if(scale_mT_over_m){
		pTFit2=TMath::Sqrt(x[0]*x[0]-1)*pTscale2;
	}
	if(scale_pT_over_m || mass2<1e-5){
		pTFit2=x[0]*pTscale2;
	}

	funcval2 = norm2*pTFit2*pow( 1+1/(beta2-2)*pTFit2*pTFit2/(pT22), -beta2 );

	funcval = funcval1/funcval2;

	return funcval;

}

double exponential(double* x, double* par) {

	double funcval;

	funcval=par[0]*exp(par[1]*x[0]);

	return funcval;

}


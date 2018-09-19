///// Head file ///////

///// C++ libraries

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <ctime>
#include <vector>
#include <string>
#include <sstream>

///// ROOT libraries

#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMatrix.h"
#include "TRandom.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"
#include "TLegend.h"
#include "TSystem.h"
#include "TF1.h"
#include "TMinuit.h"
#include "TFitter.h"
#include "TObjArray.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TRandom3.h"
#include "TMinuit.h"
#include "TMath.h"
#include "TH2.h"
#include "TH1.h"
#include "TAxis.h"

#include "TFractionFitter.h"

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "RooAddPdf.h"
#include "RooCBShape.h"
#include "RooBifurGauss.h"
#include "RooMinuit.h"
#include "RooArgList.h"
#include "RooMsgService.h"


//////////////////////////// Generation Parameters ////////////////////////////

///////////////////////////////////////////////////////////////////////////////

const static UInt_t random_seed = 22;

const static unsigned long int all_sky_LS_events = 1e+6;
const static unsigned long int all_sky_HS_events = 1e+9;

const static unsigned long int data_all_sky_LS_events = 1e+3;
const static unsigned long int data_all_sky_HS_events = 1e+6;

const static time_t time_stamp=time(0);                                   //Setting timestamp for the out files

/////////////////////////////// Dependency Paths //////////////////////////////

///////////////////////////////////////////////////////////////////////////////

// ================= HOME SERVER

//////////////////////////// Outh paths for logs and ROOT files:

const static std::string output_log = "/home/enrico/Documents/DAMPE/MyRepos/MapsFitting/SafePlace/logs/";
const static std::string output_root = "/home/enrico/Documents/DAMPE/MyRepos/MapsFitting/Safelace/results/";

//////////////////////////////// SBI Parameters ///////////////////////////////

///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////// Functions         ;-)

////////////////// Stuff functions //////////////////

extern std::string output_path_creator(const Int_t out_choose);
extern void create_and_initialize_log(std::ofstream &log_file);
extern void log_file_init(std::ofstream &out_file);
extern TH1D* TH2toTH1(TH2D* h2);
extern void create_binning(Int_t n_bin_lat,Double_t lat_bin_min,Double_t lat_bin_max,Double_t* &binning,Bool_t cos_scale);
extern void load_1D_histos(TH1D Templates_LS[],TH1D &DataHisto_I_LS,TH1D &DataHisto_NS_LS,TH1D &DataHisto_EW_LS,TH1D &DataHisto_FB_LS,TH1D &MixedDataHisto_NS_EW_LS,TH1D &MixedDataHisto_NS_FB_LS,TH1D &MixedDataHisto_EW_FB_LS,TH1D &FullMixedDataHisto_LS,TH1D Templates_HS[],TH1D &DataHisto_I_HS,TH1D &DataHisto_NS_HS,TH1D &DataHisto_EW_HS,TH1D &DataHisto_FB_HS,TH1D &MixedDataHisto_NS_EW_HS,TH1D &MixedDataHisto_NS_FB_HS,TH1D &MixedDataHisto_EW_FB_HS,TH1D &FullMixedDataHisto_HS,std::string template_path,std::string data_path,std::ofstream &log_file);
extern void read_from_file(std::string template_path,std::string data_path,std::ofstream &output_log_file);
extern void generate_and_fit(std::ofstream &output_log_file);

////////////////// Analysis functions //////////////////

extern void generate_LS_templates(TH2D* Template_Iso_LS,TH2D* Template_AniNS_LS,TH2D* Template_AniEW_LS,TH2D* Template_AniFB_LS,TH1D &Template_hwNS,TH1D &Template_hwEW,TH1D &Template_hwFB,std::ofstream &log_file,TRandom3 &r_gen);
extern void generate_HS_templates(TH2D* Template_Iso_HS,TH2D* Template_AniNS_HS,TH2D* Template_AniEW_HS,TH2D* Template_AniFB_HS,TH1D &Template_hwNS,TH1D &Template_hwEW,TH1D &Template_hwFB,std::ofstream &log_file,TRandom3 &r_gen);

extern void generate_LS_data(TH2D* Data_Iso_LS,TH2D* Data_AniNS_LS,TH2D* Data_AniEW_LS,TH2D* Data_AniFB_LS,TH2D* MixedData_NS_EW_LS,TH2D* MixedData_NS_FB_LS,TH2D* MixedData_EW_FB_LS,TH2D* FullMixedData_LS,TH1D &Data_hwNS,TH1D &Data_hwEW,TH1D &Data_hwFB,std::ofstream &log_file,TRandom3 &r_gen);
extern void generate_HS_data(TH2D* Data_Iso_HS,TH2D* Data_AniNS_HS,TH2D* Data_AniEW_HS,TH2D* Data_AniFB_HS,TH2D* MixedData_NS_EW_HS,TH2D* MixedData_NS_FB_HS,TH2D* MixedData_EW_FB_HS,TH2D* FullMixedData_HS,TH1D &Data_hwNS,TH1D &Data_hwEW,TH1D &Data_hwFB,std::ofstream &log_file,TRandom3 &r_gen);

extern void normalize_LS_templates(TH2D* Template_Iso_LS,TH2D* Template_AniNS_LS,TH2D* Template_AniEW_LS,TH2D* Template_AniFB_LS,std::ofstream &log_file);
extern void normalize_HS_templates(TH2D* Template_Iso_HS,TH2D* Template_AniNS_HS,TH2D* Template_AniEW_HS,TH2D* Template_AniFB_HS,std::ofstream &log_file);

extern void normalize_LS_data(TH2D* Data_Iso_LS,TH2D* Data_AniNS_LS,TH2D* Data_AniEW_LS,TH2D* Data_AniFB_LS,TH2D* MixedData_NS_EW_LS,TH2D* MixedData_NS_FB_LS,TH2D* MixedData_EW_FB_LS,TH2D* FullMixedData_LS,std::ofstream &log_file);
extern void normalize_HS_data(TH2D* Data_Iso_HS,TH2D* Data_AniNS_HS,TH2D* Data_AniEW_HS,TH2D* Data_AniFB_HS,TH2D* MixedData_NS_EW_HS,TH2D* MixedData_NS_FB_HS,TH2D* MixedData_EW_FB_HS,TH2D* FullMixedData_HS,std::ofstream &log_file);

////////////////// RooFit TemplateFit functions //////////////////

void ZeroRooFitVerbosity();
void RemoveZeroes(TH1* h);
void ResetRooFitVerbosity();
void ZeroRooFitVerbosity();
void fcnchisq(int& npar, double* deriv, double& f, double par[], int flag);
void fcnlike(int& npar, double* deriv, double& f, double par[], int flag);
double func(int npar, double par[], double _dochisqfit);

TH1* TemplateFitRF(TH1* dat, int ncomp, TH1* temp[], double _res[], double _reserr[], double initialguess[], bool fixtoinitialguess[], bool quiet, bool resettemplates, bool kClamping);
TH1* TemplateFitBH(TH1* dat, int ncomp, TH1* temp[], double _res[], double _reserr[], double initialguess[], bool quiet, bool resettemplates, bool kClamping, bool kChiSq,std::ofstream &log_file);





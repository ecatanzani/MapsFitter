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

///////////////////////// RooFit

using namespace RooFit;

static TH1* __dat__;
static TH1** __temp__;
static int NDATA;

//////////////////////////////////


//////////////////////////// Generation Parameters ////////////////////////////

///////////////////////////////////////////////////////////////////////////////

const static UInt_t random_seed = 22;
const static Int_t number_SBI_files = 3;
const static Int_t all_sky_events = 1e+9;
const static Int_t data_all_sky_events = 1e+5;

const static time_t time_stamp=time(0);                                     //Setting timestamp for the out files

/////////////////////////////// Dependency Paths //////////////////////////////

///////////////////////////////////////////////////////////////////////////////


// ================= LOCAL COMPUTER

//////////////////////////// Input path for SBI data files: 

const static TString sbi_path = "/Users/enrico/Documents/Università/Magistrale/Tesi/MyCode/MyRepos/GitHub/MapsFitting/SampleSBI/";

//////////////////////////// Input path for the acceptance final histo:

const static std::string IsoMap_final_plot = "/Users/enrico/Documents/Università/Magistrale/Tesi/MyCode/MyRepos/GitHub/Salomon/results/010_1528749854_maps_result.root";

//////////////////////////// Outh paths for logs and ROOT files:

const static std::string output_log = "/Users/enrico/Documents/Università/Magistrale/Tesi/MyCode/MyRepos/GitHub/MapsFitting/logs/";
const static std::string output_root = "/Users/enrico/Documents/Università/Magistrale/Tesi/MyCode/MyRepos/GitHub/MapsFitting/results/";


// ================= INFN CNAF

//////////////////////////// Input path for SBI data files:

//const static TString sbi_path = "/storage/gpfs_data/dampe/users/ecatanzani/MyRepos/DAMPE/Salomon/SampleSBI/";

//////////////////////////// Input path for the acceptance final histo:

//const static string evDistr_final_plot = "/storage/gpfs_data/dampe/users/ecatanzani/MyRepos/DAMPE/Stuff/DAMPE-GAcceptance/results/1527159147_acceptance_result.root";

//////////////////////////// Outh paths for logs and ROOT files:

//const static string output_log = "/storage/gpfs_data/dampe/users/ecatanzani/MyRepos/DAMPE/Salomon/logs/";
//const static string output_root = "/storage/gpfs_data/dampe/users/ecatanzani/MyRepos/DAMPE/Salomon/results/";

//////////////////////////////// SBI Parameters ///////////////////////////////

///////////////////////////////////////////////////////////////////////////////

const static TString sbi_subsample = "010";
const static std::string string_sbi_subsample = "010";                           //This is usefull into the function that writes log files and output ROOT files

///////////////////////////////////////////////////////////////////////////////////////////// Functions         ;-)

////////////////// Stuff functions //////////////////

extern std::string output_path_creator(const Int_t out_choose);
extern void create_and_initialize_log(std::ofstream &log_file);
extern void log_file_init(std::ofstream &out_file);
extern void create_binning(Int_t n_bin_lat,Double_t lat_bin_min,Double_t lat_bin_max,Double_t* &binning,Bool_t cos_scale);
extern void obtain_IsoMap(TH2D &IsoMap,std::ofstream &log_file);
extern void normalize_map(TH2D &Map);
extern Double_t get_tot_events(TH2D &Map);
extern Double_t get_bin_events(TH2D &IsoMap,Double_t glon,Double_t glat);
extern void BinXBin_difference(TH2D *AniMap,TH2D *IsoMap);
extern void EvaluatePercentage(Double_t &percentage,Int_t entries,Int_t idx,std::ofstream &log_file,bool allsky,bool isomap,bool phishing,bool data);

////////////////// Analysis functions //////////////////

extern void generate_templates(TH2D &TemplateIso,TH2D &TemplateAniNS,TH2D &TemplateAniEW,TH2D &TemplateAniFB,std::ofstream &log_file);
extern void generate_data(TH2D &dataI,TH2D &dataNS,TH2D &dataEW,TH2D &dataFB,std::ofstream &log_file);
extern void normalize_templates(TH2D &TemplateIso,TH2D &TemplateAniNS,TH2D &TemplateAniEW,TH2D &TemplateAniFB,std::ofstream &log_file);
extern void normalize_data(TH2D &dataI,TH2D &dataNS,TH2D &dataEW,TH2D &dataFB,std::ofstream &log_file);

/*
extern void gerenate_ani_map(TH2D &IsoMap,TH2D &AniMap_NS,TH2D &AniMap_EW,TH2D &AniMap_FB,TH1D &h_wNS,TH1D &h_wEW,TH1D &h_wFB,Int_t sky_events,std::ofstream &log_file);
extern void gerenate_ani_map_phishing(TH2D &IsoMap,TH2D &AniMap_NS,TH2D &AniMap_EW,TH2D &AniMap_FB,TH1D &h_wNS,TH1D &h_wEW,TH1D &h_wFB,Int_t sky_events,std::ofstream &log_file);
extern void compute_maps_ratio(TH2D* iso_map[],Int_t maps_number,TH2D* iso_map_ModIris,TH2D* mapsRatio,TH1D* hRatio);
extern void create_ratio_distribution(TH2D* mapsRatio,TH1D* hRatio);
extern void maps_ratio(TH2D &RatioAniMap_NS,TH2D &NAniMap_NS,TH2D &phNAniMap_NS,TH2D &RatioAniMap_EW,TH2D &NAniMap_EW,TH2D &phNAniMap_EW,TH2D &RatioAniMap_FB,TH2D &NAniMap_FB,TH2D &phNAniMap_FB,TH2D &Isolate_NS,TH2D &Isolate_EW,TH2D &Isolate_FB,TH2D &NIsoMap);
extern void fit_map(TH2D &Isolate_NS,TH2D &Isolate_EW,TH2D &Isolate_FB);
extern void fit_NN_map(TH2D &Isolate_NN_NS,TH2D &Isolate_NN_EW,TH2D &Isolate_NN_FB);
extern void isolate_dipoles(TH2D &IsoMap,TH2D &AniMap_NS,TH2D &AniMap_EW,TH2D &AniMap_FB,TH2D &Isolate_NN_NS,TH2D &Isolate_NN_EW,TH2D &Isolate_NN_FB);
extern void generate_allsky_isomap(TH2D &AllSky_IsoMap,std::ofstream &log_file);
extern void generate_allsky_animap(TH2D &AllSky_AniMap_NS,TH2D &AllSky_AniMap_EW,TH2D &AllSky_AniMap_FB,std::ofstream &log_file);
extern void AllSky_MapFit(TH2D &AllSky_Isolate_NS,TH2D &AllSky_Isolate_EW,TH2D &AllSky_Isolate_FB);
extern void subtract_by_bin(TH2D &NAllSky_AniMap_NS,TH2D &NAllSky_IsoMap);
extern void Chi2Fit(TH2D &MapsIsoRaio_NS,TH2D &MapsIsoRaio_EW,TH2D &MapsIsoRaio_FB);
*/


////////////////// RooFit functions //////////////////

void ZeroRooFitVerbosity();
void RemoveZeroes(TH1* h);
void ResetRooFitVerbosity();
void ZeroRooFitVerbosity();
void fcnchisq(int& npar, double* deriv, double& f, double par[], int flag);
void fcnlike(int& npar, double* deriv, double& f, double par[], int flag);
double func(int npar, double par[], double _dochisqfit);
TH1D* TH2toTH1(  TH2D* h2 );

TH1* TemplateFitRF(TH1* dat, int ncomp, TH1* temp[], double _res[], double _reserr[], double initialguess[], bool fixtoinitialguess[], bool quiet, bool resettemplates, bool kClamping);
TH1* TemplateFitBH(TH1* dat, int ncomp, TH1* temp[], double _res[], double _reserr[], double initialguess[], bool quiet, bool resettemplates, bool kClamping, bool kChiSq);



////////////////// TMinuit Fitting functions //////////////////

extern void fcn_NS(int &npar, double *deriv, double &f, double par[], int flag);
extern void fcn_EW(int &npar, double *deriv, double &f, double par[], int flag);
extern void fcn_FB(int &npar, double *deriv, double &f, double par[], int flag);

extern void fcn_NN_NS(int &npar, double *deriv, double &f, double par[], int flag);
extern void fcn_NN_EW(int &npar, double *deriv, double &f, double par[], int flag);
extern void fcn_NN_FB(int &npar, double *deriv, double &f, double par[], int flag);

extern void fcn_AllSky_NS(int &npar, double *deriv, double &f, double par[], int flag);
extern void fcn_AllSky_EW(int &npar, double *deriv, double &f, double par[], int flag);
extern void fcn_AllSky_FB(int &npar, double *deriv, double &f, double par[], int flag);

extern void fcn_Chi2_NS(int &npar, double *deriv, double &f, double par[], int flag);
extern void fcn_Chi2_EW(int &npar, double *deriv, double &f, double par[], int flag);
extern void fcn_Chi2_FB(int &npar, double *deriv, double &f, double par[], int flag);

extern void fit_map(TH2D &Isolate_NS,TH2D &Isolate_EW,TH2D &Isolate_FB);
extern void AllSky_MapFit(TH2D &AllSky_Isolate_NS,TH2D &AllSky_Isolate_EW,TH2D &AllSky_Isolate_FB);
extern void Chi2Fit(TH2D &MapsIsoRaio_NS,TH2D &MapsIsoRaio_EW,TH2D &MapsIsoRaio_FB);






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

#include <boost/progress.hpp>

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
#include "TF2.h"
#include "TAxis.h"
#include "TStyle.h"

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

//const static unsigned long int all_sky_LS_events = 1e+5;
//const static unsigned long int all_sky_HS_events = 1e+7;

const static unsigned long int data_all_sky_LS_events = 1e+3;
const static unsigned long int data_all_sky_HS_events = 1e+9;

const static time_t time_stamp=time(0);                                   //Setting timestamp for the out files

/////////////////////////////// Dependency Paths //////////////////////////////

///////////////////////////////////////////////////////////////////////////////

// ================= HOME SERVER

//////////////////////////// Outh paths for logs and ROOT files:

const static std::string output_log = "/home/enrico/Documents/DAMPE/MyRepos/MapsFitting/AnaliticalTemplates/logs/";
const static std::string output_root = "/home/enrico/Documents/DAMPE/MyRepos/MapsFitting/AnaliticalTemplates/results/";
const static std::string DAMPE_Iso_Map = "/home/enrico/Documents/DAMPE/MyRepos/Salomon/results/FullHistos.root";

//////////////////////////////// SBI Parameters ///////////////////////////////

///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////// Functions         ;-)

////////////////// Stuff functions //////////////////

extern std::string output_path_creator(const Int_t out_choose);
extern void create_and_initialize_log(std::ofstream &log_file);
extern void log_file_init(std::ofstream &out_file);
extern void TH2toTH1(TH1D &Histo1D,TH2D &Histo2D);
extern void create_binning(Int_t n_bin_lat,Double_t lat_bin_min,Double_t lat_bin_max,Double_t* &binning,Bool_t cos_scale);
extern void load_1D_histos(TH1D Templates_LS[],TH1D &DataHisto_I_LS,TH1D &DataHisto_NS_LS,TH1D &DataHisto_EW_LS,TH1D &DataHisto_FB_LS,TH1D &MixedDataHisto_NS_EW_LS,TH1D &MixedDataHisto_NS_FB_LS,TH1D &MixedDataHisto_EW_FB_LS,TH1D &FullMixedDataHisto_LS,TH1D Templates_HS[],TH1D &DataHisto_I_HS,TH1D &DataHisto_NS_HS,TH1D &DataHisto_EW_HS,TH1D &DataHisto_FB_HS,TH1D &MixedDataHisto_NS_EW_HS,TH1D &MixedDataHisto_NS_FB_HS,TH1D &MixedDataHisto_EW_FB_HS,TH1D &FullMixedDataHisto_HS,std::string template_path,std::string data_path,std::ofstream &log_file);
extern void read_from_file(std::string template_path,std::string data_path,std::ofstream &output_log_file);
extern void generate_and_fit(std::ofstream &output_log_file);
extern void generate_templates(std::ofstream &log_file,std::string &templates_path);
extern void generate_data(std::ofstream &log_file,std::string &data_path);
extern void read_DAMPE_FullIso(TH2D &DAMPE_FullIso);
extern void obtain_DAMPE_binning_info(TH2D &DAMPE_FullIso,Int_t &DAMPE_n_bin_lon,Double_t &DAMPE_lon_bin_min,Double_t &DAMPE_lon_bin_max,Int_t &DAMPE_n_bin_lat,Double_t &DAMPE_lat_bin_min,Double_t &DAMPE_lat_bin_max,Double_t* &DAMPE_binning);

////////////////// Analysis functions //////////////////

extern void generate_LS_templates(
                                  TH2D &Template_Iso_LS,
                                  TH2D &Template_AniNS_LS,
                                  TH2D &Template_AniEW_LS,
                                  TH2D &Template_AniFB_LS,
                                  TH1D &Template_hwI,
                                  TH1D &Template_hwNS,
                                  TH1D &Template_hwEW,
                                  TH1D &Template_hwFB,
                                  std::ofstream &log_file,
                                  TCanvas &FCanvas,
                                  TF2 &dI,
                                  TF2 &dNS,
                                  TF2 &dEW,
                                  TF2 &dFB
                                  );

extern void generate_HS_templates(
                                  TH2D &Template_Iso_HS,
                                  TH2D &Template_AniNS_HS,
                                  TH2D &Template_AniEW_HS,
                                  TH2D &Template_AniFB_HS,
                                  TH1D &Template_hwI,
                                  TH1D &Template_hwNS,
                                  TH1D &Template_hwEW,
                                  TH1D &Template_hwFB,
                                  std::ofstream &log_file,
                                  TF2 &dI,
                                  TF2 &dNS,
                                  TF2 &dEW,
                                  TF2 &dFB
                                  );

extern void generate_LS_data(
                             TH2D &Data_Iso_LS,
                             TH2D &Data_AniNS_LS,
                             TH2D &Data_AniEW_LS,
                             TH2D &Data_AniFB_LS,
                             TH2D &MixedData_NS_EW_LS,
                             TH2D &MixedData_NS_FB_LS,
                             TH2D &MixedData_EW_FB_LS,
                             TH2D &FullMixedData_LS,
                             TH1D &Data_hwI,
                             TH1D &Data_hwNS,
                             TH1D &Data_hwEW,
                             TH1D &Data_hwFB,
                             TF2 &dI,
                             TF2 &dNS,
                             TF2 &dEW,
                             TF2 &dFB,
                             std::ofstream &log_file,
                             TRandom3 &r_gen
                             );

extern void generate_HS_data(
                             TH2D &Data_Iso_HS,
                             TH2D &Data_AniNS_HS,
                             TH2D &Data_AniEW_HS,
                             TH2D &Data_AniFB_HS,
                             TH2D &MixedData_NS_EW_HS,
                             TH2D &MixedData_NS_FB_HS,
                             TH2D &MixedData_EW_FB_HS,
                             TH2D &FullMixedData_HS,
                             TH1D &Data_hwI,
                             TH1D &Data_hwNS,
                             TH1D &Data_hwEW,
                             TH1D &Data_hwFB,
                             TF2 &dI,
                             TF2 &dNS,
                             TF2 &dEW,
                             TF2 &dFB,
                             std::ofstream &log_file,
                             TRandom3 &r_gen
                             );

extern Double_t get_TF2_max(TF2 &function,Double_t theta);
extern Double_t get_XY_TF2_max(TF2 &function,Double_t theta,Double_t phi);

extern void fill_iso_background(TH2D &Data_Iso,TRandom3 &r_gen,std::ofstream &log_file,bool high_statistic);

extern void get_cleaned_LS_histos(
                                  TH2D &Data_AniNS_LS,
                                  TH2D &Data_AniEW_LS,
                                  TH2D &Data_AniFB_LS,
                                  TH2D &CData_AniNS_LS,
                                  TH2D &CData_AniEW_LS,
                                  TH2D &CData_AniFB_LS,
                                  TH2D &Data_Iso_LS
                                  );

extern void get_cleaned_HS_histos(
                                  TH2D &Data_AniNS_HS,
                                  TH2D &Data_AniEW_HS,
                                  TH2D &Data_AniFB_HS,
                                  TH2D &CData_AniNS_HS,
                                  TH2D &CData_AniEW_HS,
                                  TH2D &CData_AniFB_HS,
                                  TH2D &Data_Iso_HS
                                  );

extern void clean_bin2bin(TH2D* HDipole,TH2D &Data_Iso);

extern void normalize_LS_templates(TH2D &Template_Iso_LS,TH2D &Template_AniNS_LS,TH2D &Template_AniEW_LS,TH2D &Template_AniFB_LS,std::ofstream &log_file);
extern void normalize_HS_templates(TH2D &Template_Iso_HS,TH2D &Template_AniNS_HS,TH2D &Template_AniEW_HS,TH2D &Template_AniFB_HS,std::ofstream &log_file);

extern void normalize_LS_data(TH2D &Data_Iso_LS,TH2D &Data_AniNS_LS,TH2D &Data_AniEW_LS,TH2D &Data_AniFB_LS,TH2D &MixedData_NS_EW_LS,TH2D &MixedData_NS_FB_LS,TH2D &MixedData_EW_FB_LS,TH2D &FullMixedData_LS,std::ofstream &log_file);
extern void normalize_HS_data(TH2D &Data_Iso_HS,TH2D &Data_AniNS_HS,TH2D &Data_AniEW_HS,TH2D &Data_AniFB_HS,TH2D &MixedData_NS_EW_HS,TH2D &MixedData_NS_FB_HS,TH2D &MixedData_EW_FB_HS,TH2D &FullMixedData_HS,std::ofstream &log_file);

extern void getPull(TH1* Data,TH1* Templates[],Double_t res[],TH1D &hPull);

extern void generate_DAMPE_templates(TH2D &DAMPE_FullIso,TH2D &DAMPE_Template_Iso,TH2D &DAMPE_Template_AniNS,TH2D &DAMPE_Template_AniEW,TH2D &DAMPE_Template_AniFB,TH1D &DAMPE_Template_hwNS,TH1D &DAMPE_Template_hwEW,TH1D &DAMPE_Template_hwFB,std::ofstream &log_file,TF2 &dI,TF2 &dNS,TF2 &dEW,TF2 &dFB,TCanvas &FCanvas);

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

///// Head file ///////

///// C++ libraries

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <ctime>
#include <vector>
#include <string>
#include <sys/stat.h>
#include <unistd.h>

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
#include "TRandom2.h"
#include "TRandom3.h"
#include "TMinuit.h"
#include "TMath.h"
#include "TH2.h"
#include "TH1.h"
#include "TF2.h"
#include "TAxis.h"
#include "TStyle.h"
#include "THStack.h"
#include "TMatrixD.h"

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

///////////////////////////////////////////////////////////////////////////////////////////// Functions         ;-)

////////////////// Stuff functions //////////////////

extern void config_reader(
                            ULong64_t &data_LS_events,
                            ULong64_t &data_HS_events,
                            UInt_t &NbinX,
                            UInt_t &NbinY,
                            Bool_t &write_tmp_histos,
                            Bool_t &all_sky_simulation,
                            Bool_t &DAMPE_simulation,
                            Bool_t &DAMPE_relative_simulation,
                            std::vector<Double_t> &NS_anisotropy,
                            std::vector<Double_t> &EW_anisotropy,
                            std::vector<Double_t> &FB_anisotropy,
                            std::string &output_log,
                            std::string &output_root,
                            std::string &DAMPE_Iso_Map,
                            std::string &seeds_path
                          );

extern int numerize(std::string tmp_string);

extern std::string output_path_creator(
                                           std::string output_log,
                                           std::string output_root,
                                           const Int_t out_choose,
                                           time_t time_stamp,
                                           Double_t NS_dipole = 0,
                                           Double_t EW_dipole = 0,
                                           Double_t FB_dipole = 0,
                                           bool DAMPE = false
                                       );

extern void create_and_initialize_log(std::ofstream &log_file,ULong64_t data_LS_events,ULong64_t data_HS_events,time_t time_stamp);
extern void log_file_init(std::ofstream &out_file);

extern void TH2toTH1_obj(TH1D &Histo1D,TH2D &Histo2D);
extern void TH2toTH1_ptr(TH1D &Histo1D,TH2D* Histo2D);

extern void create_binning(Int_t n_bin_lat,Double_t lat_bin_min,Double_t lat_bin_max,Double_t* &binning,Bool_t cos_scale);
extern void load_1D_histos(TH1D Templates_LS[],TH1D &DataHisto_I_LS,TH1D &DataHisto_NS_LS,TH1D &DataHisto_EW_LS,TH1D &DataHisto_FB_LS,TH1D &MixedDataHisto_NS_EW_LS,TH1D &MixedDataHisto_NS_FB_LS,TH1D &MixedDataHisto_EW_FB_LS,TH1D &FullMixedDataHisto_LS,TH1D Templates_HS[],TH1D &DataHisto_I_HS,TH1D &DataHisto_NS_HS,TH1D &DataHisto_EW_HS,TH1D &DataHisto_FB_HS,TH1D &MixedDataHisto_NS_EW_HS,TH1D &MixedDataHisto_NS_FB_HS,TH1D &MixedDataHisto_EW_FB_HS,TH1D &FullMixedDataHisto_HS,std::string template_path,std::string data_path,std::ofstream &log_file);
extern void read_from_file(std::string template_path,std::string data_path,std::ofstream &output_log_file);

extern void generate_and_fit(
                                 ULong64_t data_LS_events,
                                 ULong64_t data_HS_events,
                                 Bool_t write_tmp_histos,
                                 Bool_t all_sky_simulation,
                                 Bool_t DAMPE_simulation,
                                 Bool_t DAMPE_relative_simulation,
                                 std::vector<Double_t> &NS_anisotropy,
                                 std::vector<Double_t> &EW_anisotropy,
                                 std::vector<Double_t> &FB_anisotropy,
                                 std::string output_log,
                                 std::string output_root,
                                 std::string DAMPE_Iso_Map,
                                 std::string seeds_path,
                                 UInt_t ani_values,
                                 std::ofstream &output_log_file,
                                 time_t time_stamp
                             );

extern void generate_templates(
                                   std::ofstream &output_log_file,
                                   ULong64_t data_LS_events,
                                   ULong64_t data_HS_events,
                                   std::string output_log,
                                   std::string output_root,
                                   std::string DAMPE_Iso_Map,
                                   time_t time_stamp
                               );

extern void generate_data_interface(
                                        std::ofstream &output_log_file,
                                        std::string output_log,
                                        std::string output_root,
                                        std::string DAMPE_Iso_Map,
                                        std::string seeds_path,
                                        time_t time_stamp,
                                        std::vector<Double_t> &NS_anisotropy,
                                        std::vector<Double_t> &EW_anisotropy,
                                        std::vector<Double_t> &FB_anisotropy,
                                        ULong64_t data_LS_events,
                                        ULong64_t data_HS_events,
                                        UInt_t ani_values,
                                        Bool_t all_sky_simulation,
                                        Bool_t DAMPE_simulation,
                                        Bool_t DAMPE_relative_simulation
                                    );

extern void read_DAMPE_FullIso(
                                TH2D &DAMPE_ReferenceMap_LS,
                                TH2D &DAMPE_ReferenceMap_HS,
                                ULong64_t data_LS_events,
                                ULong64_t data_HS_events,
                                std::string DAMPE_Iso_Map
                               );

extern void get_DAMPE_templates(
                                    TH2D &DAMPE_Template_Iso_LS,
                                    TH2D &DAMPE_Template_AniNS_LS,
                                    TH2D &DAMPE_Template_AniEW_LS,
                                    TH2D &DAMPE_Template_AniFB_LS,
                                    TH2D &DAMPE_Template_Iso_HS,
                                    TH2D &DAMPE_Template_AniNS_HS,
                                    TH2D &DAMPE_Template_AniEW_HS,
                                    TH2D &DAMPE_Template_AniFB_HS,
                                    TH2D &DAMPE_ReferenceMap_LS,
                                    TH2D &DAMPE_ReferenceMap_HS,
                                    TH2D &Template_Iso_LS,
                                    TH2D &Template_AniNS_LS,
                                    TH2D &Template_AniEW_LS,
                                    TH2D &Template_AniFB_LS,
                                    TH2D &Template_Iso_HS,
                                    TH2D &Template_AniNS_HS,
                                    TH2D &Template_AniEW_HS,
                                    TH2D &Template_AniFB_HS
                                );

extern void read_templates(
                                TH2D &Template_Iso_LS,
                                TH2D &Template_AniNS_LS,
                                TH2D &Template_AniEW_LS,
                                TH2D &Template_AniFB_LS,
                                TH2D &Template_Iso_HS,
                                TH2D &Template_AniNS_HS,
                                TH2D &Template_AniEW_HS,
                                TH2D &Template_AniFB_HS,
                                std::ofstream &output_log_file,
                                bool DAMPE,
                                std::string template_out_path,
                                std::string DAMPE_template_out_path
                           );

extern void read_relative_templates(
                                        TH2D &relative_DAMPE_Template_Iso_LS,
                                        TH2D &relative_DAMPE_Template_AniNS_LS,
                                        TH2D &relative_DAMPE_Template_AniEW_LS,
                                        TH2D &relative_DAMPE_Template_AniFB_LS,
                                        TH2D &relative_DAMPE_Template_Iso_HS,
                                        TH2D &relative_DAMPE_Template_AniNS_HS,
                                        TH2D &relative_DAMPE_Template_AniEW_HS,
                                        TH2D &relative_DAMPE_Template_AniFB_HS,
                                        std::ofstream &output_log_file,
                                        std::string DAMPE_template_out_path
                                    );

extern void normalize_DAMPE_templates(
                                        TH2D* tmp_DAMPE_Template_Iso_LS,
                                        TH2D* tmp_DAMPE_Template_AniNS_LS,
                                        TH2D* tmp_DAMPE_Template_AniEW_LS,
                                        TH2D* tmp_DAMPE_Template_AniFB_LS,
                                        TH2D* tmp_DAMPE_Template_Iso_HS,
                                        TH2D* tmp_DAMPE_Template_AniNS_HS,
                                        TH2D* tmp_DAMPE_Template_AniEW_HS,
                                        TH2D* tmp_DAMPE_Template_AniFB_HS
                                      );

extern Double_t compute_integral(TH2D* histo);

extern void create_simulation_seeds(std::string seeds_path,UInt_t ani_values);
extern inline bool path_exists (const std::string& path);
extern bool repetition_seed(UInt_t tmp_seed,std::vector<UInt_t> &vseeds,Int_t position);

extern void scale_reference_map(
                                    TH2D* DAMPE_ReferenceMap,
                                    TH2D &DAMPE_ReferenceMap_scaled,
                                    const ULong64_t n_events,
                                    bool low_stat
                                );

extern void get_relative_histo(TH2D &relative_DAMPE_histo,TH2D* data_DAMPE_histo,TH2D &reference_DAMPE_histo);

////////////////// Analysis functions //////////////////

extern void compute_templates(
                                std::string output_log,
                                std::string output_root,
                                time_t time_stamp,
                                std::ofstream &output_log_file,
                                TH2D &DAMPE_ReferenceMap_LS,
                                TH2D &DAMPE_ReferenceMap_HS,
                                std::string template_out_path,
                                std::string DAMPE_template_out_path
                              );

extern void generate_LS_templates(
                                      TH2D &Template_Iso_LS,
                                      TH2D &Template_AniNS_LS,
                                      TH2D &Template_AniEW_LS,
                                      TH2D &Template_AniFB_LS,
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
                                      std::ofstream &log_file,
                                      TF2 &dI,
                                      TF2 &dNS,
                                      TF2 &dEW,
                                      TF2 &dFB
                                  );

extern void generate_LS_data(
                                Double_t NS_anisotropy,
                                Double_t EW_anisotropy,
                                Double_t FB_anisotropy,
                                TH2D* Data_Iso_LS,
                                TH2D* Data_AniNS_LS,
                                TH2D* Data_AniEW_LS,
                                TH2D* Data_AniFB_LS,
                                TH2D* MixedData_NS_EW_LS,
                                TH2D* MixedData_NS_FB_LS,
                                TH2D* MixedData_EW_FB_LS,
                                TH2D* FullMixedData_LS,
                                TH2D &Template_Iso_LS,
                                TH2D &Template_AniNS_LS,
                                TH2D &Template_AniEW_LS,
                                TH2D &Template_AniFB_LS,
                                std::ofstream &log_file,
                                UInt_t tmp_seed,
                                ULong64_t data_LS_events
                             );

extern void generate_HS_data(
                                Double_t NS_anisotropy,
                                Double_t EW_anisotropy,
                                Double_t FB_anisotropy,
                                TH2D* Data_Iso_HS,
                                TH2D* Data_AniNS_HS,
                                TH2D* Data_AniEW_HS,
                                TH2D* Data_AniFB_HS,
                                TH2D* MixedData_NS_EW_HS,
                                TH2D* MixedData_NS_FB_HS,
                                TH2D* MixedData_EW_FB_HS,
                                TH2D* FullMixedData_HS,
                                TH2D &Template_Iso_HS,
                                TH2D &Template_AniNS_HS,
                                TH2D &Template_AniEW_HS,
                                TH2D &Template_AniFB_HS,
                                std::ofstream &log_file,
                                UInt_t tmp_seed,
                                ULong64_t data_HS_events
                             );

extern void compute_data_template(
                                    TH2D &h_Template_Data,
                                    TH2D &h_Template_Iso,
                                    TH2D &h_Template_AniNS,
                                    TH2D &h_Template_AniEW,
                                    TH2D &h_Template_AniFB,
                                    Double_t NS_anisotropy,
                                    Double_t EW_anisotropy,
                                    Double_t FB_anisotropy,
                                    std::string histo_name
                                  );

extern void interface_allSky_simulation(
                                            std::ofstream &output_log_file,
                                            std::string output_log,
                                            std::string output_root,
                                            time_t time_stamp,
                                            Double_t NS_anisotropy,
                                            Double_t EW_anisotropy,
                                            Double_t FB_anisotropy,
                                            UInt_t tmp_seed,
                                            ULong64_t data_LS_events,
                                            ULong64_t data_HS_events
                                        );

extern void interface_DAMPE_simulation(
                                           std::ofstream &output_log_file,
                                           std::string output_log,
                                           std::string output_root,
                                           time_t time_stamp,
                                           Double_t NS_anisotropy,
                                           Double_t EW_anisotropy,
                                           Double_t FB_anisotropy,
                                           UInt_t tmp_seed,
                                           ULong64_t data_LS_events,
                                           ULong64_t data_HS_events
                                       );

extern void interface_DAMPE_relative_simulation(
                                                    std::ofstream &output_log_file,
                                                    std::string output_log,
                                                    std::string output_root,
                                                    time_t time_stamp,
                                                    Double_t NS_anisotropy,
                                                    Double_t EW_anisotropy,
                                                    Double_t FB_anisotropy,
                                                    UInt_t tmp_seed,
                                                    ULong64_t data_LS_events,
                                                    ULong64_t data_HS_events,
                                                    TH2D &DAMPE_ReferenceMap_LS,
                                                    TH2D &DAMPE_ReferenceMap_HS
                                                );

extern Double_t get_TF2_max(TF2 &function,Double_t theta);
extern Double_t get_XY_TF2_max(TF2 &function,Double_t theta,Double_t phi);

extern void fill_iso_background(TH2D &Data_Iso,TRandom3 &r_gen,std::ofstream &log_file,bool high_statistic);

extern void get_cleaned_LS_histos(
                                  TH2D* Data_AniNS_LS,
                                  TH2D* Data_AniEW_LS,
                                  TH2D* Data_AniFB_LS,
                                  TH2D &CData_AniNS_LS,
                                  TH2D &CData_AniEW_LS,
                                  TH2D &CData_AniFB_LS,
                                  TH2D* Data_Iso_LS
                                  );

extern void get_cleaned_HS_histos(
                                  TH2D* Data_AniNS_HS,
                                  TH2D* Data_AniEW_HS,
                                  TH2D* Data_AniFB_HS,
                                  TH2D &CData_AniNS_HS,
                                  TH2D &CData_AniEW_HS,
                                  TH2D &CData_AniFB_HS,
                                  TH2D* Data_Iso_HS
                                  );

extern void clean_bin2bin(TH2D* HDipole,TH2D* Data_Iso);

extern void normalize_LS_templates(TH2D &Template_Iso_LS,TH2D &Template_AniNS_LS,TH2D &Template_AniEW_LS,TH2D &Template_AniFB_LS,std::ofstream &log_file);
extern void normalize_HS_templates(TH2D &Template_Iso_HS,TH2D &Template_AniNS_HS,TH2D &Template_AniEW_HS,TH2D &Template_AniFB_HS,std::ofstream &log_file);

extern void normalize_LS_data(TH2D &Data_Iso_LS,TH2D &Data_AniNS_LS,TH2D &Data_AniEW_LS,TH2D &Data_AniFB_LS,TH2D &MixedData_NS_EW_LS,TH2D &MixedData_NS_FB_LS,TH2D &MixedData_EW_FB_LS,TH2D &FullMixedData_LS,std::ofstream &log_file);
extern void normalize_HS_data(TH2D &Data_Iso_HS,TH2D &Data_AniNS_HS,TH2D &Data_AniEW_HS,TH2D &Data_AniFB_HS,TH2D &MixedData_NS_EW_HS,TH2D &MixedData_NS_FB_HS,TH2D &MixedData_EW_FB_HS,TH2D &FullMixedData_HS,std::ofstream &log_file);

extern void getPull(TH1* Data,TH1* Templates[],Double_t res[],TH1D &hPull,bool relative_fit=false,bool ex_fit=false);

extern void generate_DAMPE_templates(TH2D &DAMPE_FullIso,TH2D &DAMPE_Template_Iso,TH2D &DAMPE_Template_AniNS,TH2D &DAMPE_Template_AniEW,TH2D &DAMPE_Template_AniFB,TH1D &DAMPE_Template_hwNS,TH1D &DAMPE_Template_hwEW,TH1D &DAMPE_Template_hwFB,std::ofstream &log_file,TF2 &dI,TF2 &dNS,TF2 &dEW,TF2 &dFB,TCanvas &FCanvas);

extern void components_analysis(
                                    std::vector<TH1D*> &TemplatesProjections,
                                    std::vector<TH1D*> &DataProjections,
                                    double resfullFitResults_HS[][8],
                                    std::ofstream &log_file,
                                    Double_t NS_anisotropy,
                                    Double_t EW_anisotropy,
                                    Double_t FB_anisotropy,
                                    std::string output_log,
                                    std::string output_root,
                                    time_t time_stamp,
                                    Bool_t write_tmp_histos
                                );

extern double compute_ani_level(double covMatrix[][4],double parameters[],double par_err[],std::ofstream &log_file);

extern void allSky_singleTry_fit(
                                    std::ofstream &output_log_file,
                                    std::string output_log,
                                    std::string output_root,
                                    time_t time_stamp,
                                    UInt_t tmp_seed,
                                    Double_t NS_anisotropy,
                                    Double_t EW_anisotropy,
                                    Double_t FB_anisotropy,
                                    std::string template_out_path,
                                    std::string DAMPE_template_out_path,
                                    Bool_t write_tmp_histos,
                                    ULong64_t data_LS_events,
                                    ULong64_t data_HS_events
                                 );

extern void DAMPE_singleTry_fit(
                                    std::ofstream &output_log_file,
                                    std::string output_log,
                                    std::string output_root,
                                    time_t time_stamp,
                                    UInt_t tmp_seed,
                                    Double_t NS_anisotropy,
                                    Double_t EW_anisotropy,
                                    Double_t FB_anisotropy,
                                    std::string template_out_path,
                                    std::string DAMPE_template_out_path,
                                    ULong64_t data_LS_events,
                                    ULong64_t data_HS_events,
                                    Bool_t write_tmp_histos
                                );

extern void generate_DAMPE_LS_data(
                                       Double_t NS_anisotropy,
                                       Double_t EW_anisotropy,
                                       Double_t FB_anisotropy,
                                       TH2D* DAMPE_Data_Iso_LS,
                                       TH2D* DAMPE_Data_AniNS_LS,
                                       TH2D* DAMPE_Data_AniEW_LS,
                                       TH2D* DAMPE_Data_AniFB_LS,
                                       TH2D* DAMPE_MixedData_NS_EW_LS,
                                       TH2D* DAMPE_MixedData_NS_FB_LS,
                                       TH2D* DAMPE_MixedData_EW_FB_LS,
                                       TH2D* DAMPE_FullMixedData_LS,
                                       TH2D &DAMPE_Template_Iso_LS,
                                       TH2D &DAMPE_Template_AniNS_LS,
                                       TH2D &DAMPE_Template_AniEW_LS,
                                       TH2D &DAMPE_Template_AniFB_LS,
                                       std::ofstream &output_log_file,
                                       UInt_t tmp_seed,
                                       ULong64_t data_LS_events
                                   );

extern void generate_DAMPE_HS_data(
                                       Double_t NS_anisotropy,
                                       Double_t EW_anisotropy,
                                       Double_t FB_anisotropy,
                                       TH2D* DAMPE_Data_Iso_HS,
                                       TH2D* DAMPE_Data_AniNS_HS,
                                       TH2D* DAMPE_Data_AniEW_HS,
                                       TH2D* DAMPE_Data_AniFB_HS,
                                       TH2D* DAMPE_MixedData_NS_EW_HS,
                                       TH2D* DAMPE_MixedData_NS_FB_HS,
                                       TH2D* DAMPE_MixedData_EW_FB_HS,
                                       TH2D* DAMPE_FullMixedData_HS,
                                       TH2D &DAMPE_Template_Iso_HS,
                                       TH2D &DAMPE_Template_AniNS_HS,
                                       TH2D &DAMPE_Template_AniEW_HS,
                                       TH2D &DAMPE_Template_AniFB_HS,
                                       std::ofstream &output_log_file,
                                       UInt_t tmp_seed,
                                       ULong64_t data_HS_events
                                   );

extern void DAMPE_relative_singleTry_fit(
                                            std::ofstream &output_log_file,
                                            std::string output_log,
                                            std::string output_root,
                                            time_t time_stamp,
                                            UInt_t tmp_seed,
                                            Double_t NS_anisotropy,
                                            Double_t EW_anisotropy,
                                            Double_t FB_anisotropy,
                                            std::string template_out_path,
                                            std::string DAMPE_template_out_path,
                                            TH2D &DAMPE_ReferenceMap_LS,
                                            TH2D &DAMPE_ReferenceMap_HS,
                                            ULong64_t data_LS_events,
                                            ULong64_t data_HS_events,
                                            Bool_t write_tmp_histos
                                         );

////////////////// RooFit TemplateFit functions //////////////////

void ZeroRooFitVerbosity();
void RemoveZeroes(TH1* h);
void ResetRooFitVerbosity();
void ZeroRooFitVerbosity();
void fcnchisq(int& npar, double* deriv, double& f, double par[], int flag);
void fcnlike(int& npar, double* deriv, double& f, double par[], int flag);
double func(int npar, double par[], double _dochisqfit);

TH1* TemplateFitRF(
                        TH1* dat,
                        int ncomp,
                        TH1* temp[],
                        double _res[],
                        double _reserr[],
                        double initialguess[],
                        bool fixtoinitialguess[],
                        bool quiet,
                        bool resettemplates,
                        bool kClamping
                   );

TH1* TemplateFitBH(
                       TH1* dat,
                       int ncomp,
                       TH1* temp[],
                       double _res[],
                       double _reserr[],
                       double initialguess[],
                       bool quiet,
                       bool resettemplates,
                       bool kClamping,
                       bool kChiSq,
                       std::ofstream &log_file
                   );

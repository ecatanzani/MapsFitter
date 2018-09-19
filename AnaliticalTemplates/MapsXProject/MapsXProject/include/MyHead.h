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


//////////////////////////// Generation Parameters ////////////////////////////

///////////////////////////////////////////////////////////////////////////////

const static unsigned long int data_all_sky_LS_events = 3e+6;
const static unsigned long int data_all_sky_HS_events = 6e+6;

const static time_t time_stamp=time(0);                                   //Setting timestamp for the out files

//const static Int_t ntry = 1e+2;
const static Int_t ntry = 1;

const static bool write_tmp_histos = false;

const static Int_t ani_values = 3;
static Double_t NS_anisotropy[ani_values] = {0.1,0.01,0.001};
static Double_t EW_anisotropy[ani_values] = {0.1,0.01,0.001};
static Double_t FB_anisotropy[ani_values] = {0.1,0.01,0.001};

/////////////////////////////// Dependency Paths //////////////////////////////

///////////////////////////////////////////////////////////////////////////////


// ================= macOS

//////////////////////////// Outh paths for logs and ROOT files:

const static std::string output_log = "/Users/enrico/Documents/Università/Magistrale/Tesi/MyCode/MyRepos/GitHub/MapsFitting/AnaliticalTemplates/logs/";
const static std::string output_root = "/Users/enrico/Documents/Università/Magistrale/Tesi/MyCode/MyRepos/GitHub/MapsFitting/AnaliticalTemplates/results/";
const static std::string DAMPE_Iso_Map = "/Users/enrico/Documents/Università/Magistrale/Tesi/MyCode/MyRepos/GitHub/Salomon/results/fullHistos.root";
const static std::string seeds_path = "/Users/enrico/Documents/Università/Magistrale/Tesi/MyCode/MyRepos/GitHub/MapsFitting/AnaliticalTemplates/seeds.txt";

//////////////////////////////// SBI Parameters ///////////////////////////////

///////////////////////////////////////////////////////////////////////////////

class fitResult
{
    public:
    
    Double_t inputAni[3];
    
    Double_t chi2_LS[8];
    Double_t ndf_LS[8];
    Double_t fit_par_LS[8][4];
    Double_t fit_err_LS[8][4];
    Double_t delta_LS[8];
    Double_t sum_par_LS[8];
    
    Double_t CMatrix_Iso_LS[4][4];
    Double_t CMatrix_NS_LS[4][4];
    Double_t CMatrix_EW_LS[4][4];
    Double_t CMatrix_FB_LS[4][4];
    Double_t CMatrix_NS_EW_LS[4][4];
    Double_t CMatrix_NS_FB_LS[4][4];
    Double_t CMatrix_EW_FB_LS[4][4];
    Double_t CMatrix_Full_LS[4][4];
    
    Int_t theta_binHistos_LS;
    Int_t phi_binHistos_LS;
    
    ULong64_t events_LS;
    
    Double_t chi2_HS[8];
    Double_t ndf_HS[8];
    Double_t fit_par_HS[8][4];
    Double_t fit_err_HS[8][4];
    Double_t delta_HS[8];
    Double_t sum_par_HS[8];
    
    Double_t CMatrix_Iso_HS[4][4];
    Double_t CMatrix_NS_HS[4][4];
    Double_t CMatrix_EW_HS[4][4];
    Double_t CMatrix_FB_HS[4][4];
    Double_t CMatrix_NS_EW_HS[4][4];
    Double_t CMatrix_NS_FB_HS[4][4];
    Double_t CMatrix_EW_FB_HS[4][4];
    Double_t CMatrix_Full_HS[4][4];
    
    Int_t theta_binHistos_HS;
    Int_t phi_binHistos_HS;
    
    ULong64_t events_HS;
    
    
    fitResult()
    {
        
        for(Int_t idx = 0; idx < 8; ++idx)
        {
            chi2_LS[idx] = -1;
            ndf_LS[idx] = -1;
            delta_LS[idx] = -1;
            sum_par_LS[idx] = -999;
            
            chi2_HS[idx] = -1;
            ndf_HS[idx] = -1;
            delta_HS[idx] = -1;
            sum_par_HS[idx] = -999;
            
            for(Int_t k = 0; k < 4; ++k)
            {
                fit_par_LS[idx][k] = -1;
                fit_err_LS[idx][k] = -1;
                
                fit_par_HS[idx][k] = -1;
                fit_err_HS[idx][k] = -1;
                
                for(int j = 0; j < 4; ++j)
                {
                    CMatrix_Iso_LS[k][j] = -1;
                    CMatrix_NS_LS[k][j] = -1;
                    CMatrix_EW_LS[k][j] = -1;
                    CMatrix_FB_LS[k][j] = -1;
                    CMatrix_NS_EW_LS[k][j] = -1;
                    CMatrix_NS_FB_LS[k][j] = -1;
                    CMatrix_EW_FB_LS[k][j] = -1;
                    CMatrix_Full_LS[k][j] = -1;
                    
                    CMatrix_Iso_HS[k][j] = -1;
                    CMatrix_NS_HS[k][j] = -1;
                    CMatrix_EW_HS[k][j] = -1;
                    CMatrix_FB_HS[k][j] = -1;
                    CMatrix_NS_EW_HS[k][j] = -1;
                    CMatrix_NS_FB_HS[k][j] = -1;
                    CMatrix_EW_FB_HS[k][j] = -1;
                    CMatrix_Full_HS[k][j] = -1;
                    
                }
            }
        }
        
        theta_binHistos_LS = 0;
        phi_binHistos_LS = 0;
        theta_binHistos_HS = 0;
        phi_binHistos_HS = 0;
        
        for(Int_t idx = 0; idx < 3; ++idx)
            inputAni[idx] = -1;
        
        events_LS = 0;
        events_HS = 0;
        
    }
    
    ~fitResult() { }
    
};


///////////////////////////////////////////////////////////////////////////////////////////// Functions         ;-)

////////////////// Stuff functions //////////////////

extern std::string output_path_creator(const Int_t out_choose);
extern void create_and_initialize_log(std::ofstream &log_file);
extern void log_file_init(std::ofstream &out_file);

extern void TH2toTH1_obj(TH1D &Histo1D,TH2D &Histo2D);
extern void TH2toTH1_ptr(TH1D &Histo1D,TH2D* Histo2D);

extern void create_binning(Int_t n_bin_lat,Double_t lat_bin_min,Double_t lat_bin_max,Double_t* &binning,Bool_t cos_scale);
extern void load_1D_histos(TH1D Templates_LS[],TH1D &DataHisto_I_LS,TH1D &DataHisto_NS_LS,TH1D &DataHisto_EW_LS,TH1D &DataHisto_FB_LS,TH1D &MixedDataHisto_NS_EW_LS,TH1D &MixedDataHisto_NS_FB_LS,TH1D &MixedDataHisto_EW_FB_LS,TH1D &FullMixedDataHisto_LS,TH1D Templates_HS[],TH1D &DataHisto_I_HS,TH1D &DataHisto_NS_HS,TH1D &DataHisto_EW_HS,TH1D &DataHisto_FB_HS,TH1D &MixedDataHisto_NS_EW_HS,TH1D &MixedDataHisto_NS_FB_HS,TH1D &MixedDataHisto_EW_FB_HS,TH1D &FullMixedDataHisto_HS,std::string template_path,std::string data_path,std::ofstream &log_file);
extern void read_from_file(std::string template_path,std::string data_path,std::ofstream &output_log_file);
extern void generate_and_fit(std::ofstream &output_log_file);
extern void generate_templates(std::ofstream &log_file,std::string &templates_path);
extern void generate_data(std::ofstream &log_file,std::string &data_path);
extern void read_DAMPE_FullIso(TH2D &DAMPE_FullIso);

extern void get_DAMPE_templates(
                                    TH2D &DAMPE_Template_Iso_LS,
                                    TH2D &DAMPE_Template_AniNS_LS,
                                    TH2D &DAMPE_Template_AniEW_LS,
                                    TH2D &DAMPE_Template_AniFB_LS,
                                    TH2D &DAMPE_Template_Iso_HS,
                                    TH2D &DAMPE_Template_AniNS_HS,
                                    TH2D &DAMPE_Template_AniEW_HS,
                                    TH2D &DAMPE_Template_AniFB_HS,
                                    TH2D &DAMPE_ReferenceMap,
                                    TH2D &Template_Iso_LS,
                                    TH2D &Template_AniNS_LS,
                                    TH2D &Template_AniEW_LS,
                                    TH2D &Template_AniFB_LS,
                                    TH2D &Template_Iso_HS,
                                    TH2D &Template_AniNS_HS,
                                    TH2D &Template_AniEW_HS,
                                    TH2D &Template_AniFB_HS
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

extern void create_simulation_seeds();
extern inline bool path_exists (const std::string& path);
extern bool repetition_seed(UInt_t tmp_seed,std::vector<UInt_t> &vseeds,Int_t position);

extern void scale_reference_map(TH2D &DAMPE_ReferenceMap,TH2D* histo,bool LS);
extern void get_relative_histo(TH2D &relative_DAMPE_histo,TH2D* data_DAMPE_histo,TH2D* reference_DAMPE_histo);

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
                             TH1D &Data_hwI,
                             TH1D &Data_hwNS,
                             TH1D &Data_hwEW,
                             TH1D &Data_hwFB,
                             TF2 &dI,
                             TF2 &dNS,
                             TF2 &dEW,
                             TF2 &dFB,
                             TH2D &Template_Iso_LS,
                             TH2D &Template_AniNS_LS,
                             TH2D &Template_AniEW_LS,
                             TH2D &Template_AniFB_LS,
                             std::ofstream &log_file,
                             TRandom3 &r_gen
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
                             TH1D &Data_hwI,
                             TH1D &Data_hwNS,
                             TH1D &Data_hwEW,
                             TH1D &Data_hwFB,
                             TF2 &dI,
                             TF2 &dNS,
                             TF2 &dEW,
                             TF2 &dFB,
                             TH2D &Template_Iso_HS,
                             TH2D &Template_AniNS_HS,
                             TH2D &Template_AniEW_HS,
                             TH2D &Template_AniFB_HS,
                             std::ofstream &log_file,
                             TRandom3 &r_gen
                             );


extern void MC_generate_LS_data(
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

extern void MC_generate_HS_data(
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
                                    TRandom3 &r_gen
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
                                    TRandom3 &r_gen
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

extern void getPull(TH1* Data,TH1* Templates[],Double_t res[],TH1D &hPull);

extern void generate_DAMPE_templates(TH2D &DAMPE_FullIso,TH2D &DAMPE_Template_Iso,TH2D &DAMPE_Template_AniNS,TH2D &DAMPE_Template_AniEW,TH2D &DAMPE_Template_AniFB,TH1D &DAMPE_Template_hwNS,TH1D &DAMPE_Template_hwEW,TH1D &DAMPE_Template_hwFB,std::ofstream &log_file,TF2 &dI,TF2 &dNS,TF2 &dEW,TF2 &dFB,TCanvas &FCanvas);

extern void components_analysis(
                                std::vector<TH1D*> &TemplatesProjections,
                                std::vector<TH1D*> &DataProjections,
                                double resfullFitResults_HS[][8],
                                std::ofstream &log_file
                                );

extern double compute_ani_level(double covMatrix[][4],double parameters[],double par_err[],std::ofstream &log_file);

extern void allSky_singleTry_fit(
                                    std::ofstream &output_log_file,
                                    std::string &template_out_path,
                                    std::string &data_out_path,
                                    std::string &pools_out_path,
                                    std::string &projections_out_path,
                                    TH2D &Template_Iso_LS,
                                    TH2D &Template_AniNS_LS,
                                    TH2D &Template_AniEW_LS,
                                    TH2D &Template_AniFB_LS,
                                    TH2D &Template_Iso_HS,
                                    TH2D &Template_AniNS_HS,
                                    TH2D &Template_AniEW_HS,
                                    TH2D &Template_AniFB_HS,
                                    TH1D &Template_hwI,
                                    TH1D &Template_hwNS,
                                    TH1D &Template_hwEW,
                                    TH1D &Template_hwFB,
                                    TF2 &dI,
                                    TF2 &dNS,
                                    TF2 &dEW,
                                    TF2 &dFB,
                                    TTree &fiTree,
                                    fitResult &tmp_fit,
                                    Double_t NS_anisotropy,
                                    Double_t EW_anisotropy,
                                    Double_t FB_anisotropy,
                                    TCanvas &LS_Templates_Canvas,
                                    TCanvas &HS_Templates_Canvas,
                                    UInt_t tmp_seed
                                 );

extern void DAMPE_singleTry_fit(
                                    std::ofstream &output_log_file,
                                    std::string &data_out_path,
                                    std::string &pools_out_path,
                                    TH2D &DAMPE_Template_Iso_LS,
                                    TH2D &DAMPE_Template_AniNS_LS,
                                    TH2D &DAMPE_Template_AniEW_LS,
                                    TH2D &DAMPE_Template_AniFB_LS,
                                    TH2D &DAMPE_Template_Iso_HS,
                                    TH2D &DAMPE_Template_AniNS_HS,
                                    TH2D &DAMPE_Template_AniEW_HS,
                                    TH2D &DAMPE_Template_AniFB_HS,
                                    Double_t NS_anisotropy,
                                    Double_t EW_anisotropy,
                                    Double_t FB_anisotropy,
                                    fitResult &tmp_fit,
                                    UInt_t tmp_seed,
                                    TTree &fiTree
                                );

extern void DAMPE_relative_singleTry_fit(
                                            std::ofstream &output_log_file,
                                            std::string &data_out_path,
                                            std::string &pools_out_path,
                                            TH2D &Template_Iso_LS,
                                            TH2D &Template_AniNS_LS,
                                            TH2D &Template_AniEW_LS,
                                            TH2D &Template_AniFB_LS,
                                            TH2D &Template_Iso_HS,
                                            TH2D &Template_AniNS_HS,
                                            TH2D &Template_AniEW_HS,
                                            TH2D &Template_AniFB_HS,
                                            TH2D &DAMPE_ReferenceMap,
                                            Double_t NS_anisotropy,
                                            Double_t EW_anisotropy,
                                            Double_t FB_anisotropy,
                                            fitResult &tmp_fit,
                                            UInt_t tmp_seed,
                                            TTree &fiTree
                                         );

////////////////// RooFit TemplateFit functions //////////////////

void ZeroRooFitVerbosity();
void RemoveZeroes(TH1* h);
void ResetRooFitVerbosity();
void ZeroRooFitVerbosity();
void fcnchisq(int& npar, double* deriv, double& f, double par[], int flag);
void fcnlike(int& npar, double* deriv, double& f, double par[], int flag);
double func(int npar, double par[], double _dochisqfit);

TH1* TemplateFitRF(TH1* dat, int ncomp, TH1* temp[], double _res[], double _reserr[], double initialguess[], bool fixtoinitialguess[], bool quiet, bool resettemplates, bool kClamping);
TH1* TemplateFitBH(TH1* dat, int ncomp, TH1* temp[], double _res[], double _reserr[], double initialguess[], bool quiet, bool resettemplates, bool kClamping, bool kChiSq,std::ofstream &log_file,fitResult &tmp_fit,Int_t fit_element,bool HS);

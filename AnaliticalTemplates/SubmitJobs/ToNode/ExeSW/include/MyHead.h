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
#include <sys/stat.h>
#include <unistd.h>
#include <boost/progress.hpp>
#include <limits>

///// ROOT libraries

#include "TCanvas.h"
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

const static unsigned long int data_all_sky_LS_events = 6928;
const static unsigned long int data_all_sky_HS_events = 13856;

const static time_t time_stamp=time(0);                                   //Setting timestamp for the out files

const static bool write_tmp_histos = true;

const static Int_t ani_values = 3;
static Double_t NS_anisotropy[ani_values] = {0.1,0.01,0.001};
static Double_t EW_anisotropy[ani_values] = {0.1,0.01,0.001};
static Double_t FB_anisotropy[ani_values] = {0.1,0.01,0.001};

const static bool all_sky_simulation = false;
const static bool DAMPE_simulation = true;
const static bool DAMPE_relative_simulation = false;
const static bool fitDistribution_test = false;

const static Int_t Nbin = 9*18;
// const static Nbin = 18*36;

/////////////////////////////// Dependency Paths //////////////////////////////

///////////////////////////////////////////////////////////////////////////////

// ================= INFN CNAF

//////////////////////////// Outh paths for logs and ROOT files:

const static std::string DAMPE_Iso_Map = "/storage/gpfs_data/dampe/users/ecatanzani/MyRepos/DAMPE/Salomon/results/FullHistos.root";
const static std::string DAMPE_Iso_scaled_Maps = "/storage/gpfs_data/dampe/users/ecatanzani/MyRepos/DAMPE/MapsFitter/AnaliticalTemplates/SubmitJobs/ToNode/ExeSW/assets/scalingSoftware/scaled_reference_Isotropic_histos.root";
const static std::string templates_path = "/storage/gpfs_data/dampe/users/ecatanzani/MyRepos/DAMPE/MapsFitter/AnaliticalTemplates/SubmitJobs/ToNode/ExeSW/assets/computeTemplates/results/AllSkyTemplates.root";
const static std::string DAMPE_templates_path = "/storage/gpfs_data/dampe/users/ecatanzani/MyRepos/DAMPE/MapsFitter/AnaliticalTemplates/SubmitJobs/ToNode/ExeSW/assets/computeTemplates/results/DAMPETemplates.root";
const static std::string seeds_path = "/storage/gpfs_data/dampe/users/ecatanzani/MyRepos/DAMPE/MapsFitter/AnaliticalTemplates/SubmitJobs/ToNode/ExeSW/assets/produceSeeds/seeds.txt";


//////////////////////////////// SBI Parameters ///////////////////////////////

///////////////////////////////////////////////////////////////////////////////

class fitResult
{
public:
    
    Double_t inputAni[3];
    ULong64_t seed;
    ULong64_t seed_list_line;
    
    Double_t chi2_LS[8];
    Double_t ndf_LS[8];
    Double_t chi2_r_LS[8];
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
    
    Int_t entries_Iso_Map_LS[ani_values][Nbin];
    Int_t entries_NS_Map_LS[ani_values][Nbin];
    Int_t entries_EW_Map_LS[ani_values][Nbin];
    Int_t entries_FB_Map_LS[ani_values][Nbin];
    Int_t entries_NS_EW_Map_LS[ani_values][Nbin];
    Int_t entries_NS_FB_Map_LS[ani_values][Nbin];
    Int_t entries_EW_FB_Map_LS[ani_values][Nbin];
    Int_t entries_Full_Map_LS[ani_values][Nbin];
    
    ULong64_t events_LS;
    
    Double_t chi2_HS[8];
    Double_t ndf_HS[8];
    Double_t chi2_r_HS[8];
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
    
    Int_t entries_Iso_Map_HS[ani_values][Nbin];
    Int_t entries_NS_Map_HS[ani_values][Nbin];
    Int_t entries_EW_Map_HS[ani_values][Nbin];
    Int_t entries_FB_Map_HS[ani_values][Nbin];
    Int_t entries_NS_EW_Map_HS[ani_values][Nbin];
    Int_t entries_NS_FB_Map_HS[ani_values][Nbin];
    Int_t entries_EW_FB_Map_HS[ani_values][Nbin];
    Int_t entries_Full_Map_HS[ani_values][Nbin];
    
    ULong64_t events_HS;
    
    
    fitResult()
    {
        
        for(Int_t idx = 0; idx < 8; ++idx)
        {
            chi2_LS[idx] = -1;
            ndf_LS[idx] = -1;
            chi2_r_LS[idx] = -1;
            delta_LS[idx] = -1;
            sum_par_LS[idx] = -999;
            
            chi2_HS[idx] = -1;
            ndf_HS[idx] = -1;
            chi2_r_HS[idx] = -1;
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
        
        for(Int_t a_idx=0; a_idx < ani_values; ++a_idx)
            for(Int_t idx = 0; idx < Nbin; ++idx)
            {
            
                entries_Iso_Map_LS[a_idx][idx] = -999;
                entries_NS_Map_LS[a_idx][idx] = -999;
                entries_EW_Map_LS[a_idx][idx] = -999;
                entries_FB_Map_LS[a_idx][idx] = -999;
                entries_NS_EW_Map_LS[a_idx][idx] = -999;
                entries_NS_FB_Map_LS[a_idx][idx] = -999;
                entries_EW_FB_Map_LS[a_idx][idx] = -999;
                entries_Full_Map_LS[a_idx][idx] = -999;
            
                entries_Iso_Map_HS[a_idx][idx] = -999;
                entries_NS_Map_HS[a_idx][idx] = -999;
                entries_EW_Map_HS[a_idx][idx] = -999;
                entries_FB_Map_HS[a_idx][idx] = -999;
                entries_NS_EW_Map_HS[a_idx][idx] = -999;
                entries_NS_FB_Map_HS[a_idx][idx] = -999;
                entries_EW_FB_Map_HS[a_idx][idx] = -999;
                entries_Full_Map_HS[a_idx][idx] = -999;
            
            }
        
        events_LS = 0;
        events_HS = 0;
        
        seed = 0;
        seed_list_line = 0;
        
    }
    
    ~fitResult() { }
    
};

class relative_fitResult
{
public:
    
    Double_t inputAni[3];
    ULong64_t seed;
    ULong64_t seed_list_line;
    
    Double_t chi2_LS[7];
    Double_t ndf_LS[7];
    Double_t chi2_r_LS[7];
    Double_t fit_par_LS[7][3];
    Double_t fit_err_LS[7][3];
    Double_t delta_LS[7];
    Double_t sum_par_LS[7];
    
    Double_t CMatrix_NS_LS[3][3];
    Double_t CMatrix_EW_LS[3][3];
    Double_t CMatrix_FB_LS[3][3];
    Double_t CMatrix_NS_EW_LS[3][3];
    Double_t CMatrix_NS_FB_LS[3][3];
    Double_t CMatrix_EW_FB_LS[3][3];
    Double_t CMatrix_Full_LS[3][3];
    
    Int_t entries_NS_Map_LS[ani_values][Nbin];
    Int_t entries_EW_Map_LS[ani_values][Nbin];
    Int_t entries_FB_Map_LS[ani_values][Nbin];
    Int_t entries_NS_EW_Map_LS[ani_values][Nbin];
    Int_t entries_NS_FB_Map_LS[ani_values][Nbin];
    Int_t entries_EW_FB_Map_LS[ani_values][Nbin];
    Int_t entries_Full_Map_LS[ani_values][Nbin];
    
    Int_t theta_binHistos_LS;
    Int_t phi_binHistos_LS;
    
    ULong64_t events_LS;
    
    Double_t chi2_HS[7];
    Double_t ndf_HS[7];
    Double_t chi2_r_HS[7];
    Double_t fit_par_HS[7][3];
    Double_t fit_err_HS[7][3];
    Double_t delta_HS[7];
    Double_t sum_par_HS[7];
    
    Double_t CMatrix_NS_HS[3][3];
    Double_t CMatrix_EW_HS[3][3];
    Double_t CMatrix_FB_HS[3][3];
    Double_t CMatrix_NS_EW_HS[3][3];
    Double_t CMatrix_NS_FB_HS[3][3];
    Double_t CMatrix_EW_FB_HS[3][3];
    Double_t CMatrix_Full_HS[3][3];
    
    Int_t entries_NS_Map_HS[ani_values][Nbin];
    Int_t entries_EW_Map_HS[ani_values][Nbin];
    Int_t entries_FB_Map_HS[ani_values][Nbin];
    Int_t entries_NS_EW_Map_HS[ani_values][Nbin];
    Int_t entries_NS_FB_Map_HS[ani_values][Nbin];
    Int_t entries_EW_FB_Map_HS[ani_values][Nbin];
    Int_t entries_Full_Map_HS[ani_values][Nbin];
    
    Int_t theta_binHistos_HS;
    Int_t phi_binHistos_HS;
    
    ULong64_t events_HS;
    
    
    relative_fitResult()
    {
        
        for(Int_t idx = 0; idx < 7; ++idx)
        {
            chi2_LS[idx] = -1;
            ndf_LS[idx] = -1;
            chi2_r_LS[idx] = -1;
            delta_LS[idx] = -1;
            sum_par_LS[idx] = -999;
            
            chi2_HS[idx] = -1;
            ndf_HS[idx] = -1;
            chi2_r_HS[idx] = -1;
            delta_HS[idx] = -1;
            sum_par_HS[idx] = -999;
            
            for(Int_t k = 0; k < 3; ++k)
            {
                fit_par_LS[idx][k] = -1;
                fit_err_LS[idx][k] = -1;
                
                fit_par_HS[idx][k] = -1;
                fit_err_HS[idx][k] = -1;
                
                for(int j = 0; j < 3; ++j)
                {
                    CMatrix_NS_LS[k][j] = -1;
                    CMatrix_EW_LS[k][j] = -1;
                    CMatrix_FB_LS[k][j] = -1;
                    CMatrix_NS_EW_LS[k][j] = -1;
                    CMatrix_NS_FB_LS[k][j] = -1;
                    CMatrix_EW_FB_LS[k][j] = -1;
                    CMatrix_Full_LS[k][j] = -1;
                
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
        
        for(Int_t a_idx=0; a_idx < ani_values; ++a_idx)
            for(Int_t idx = 0; idx < 648; ++idx)
            {
            
                entries_NS_Map_LS[a_idx][idx] = -999;
                entries_EW_Map_LS[a_idx][idx] = -999;
                entries_FB_Map_LS[a_idx][idx] = -999;
                entries_NS_EW_Map_LS[a_idx][idx] = -999;
                entries_NS_FB_Map_LS[a_idx][idx] = -999;
                entries_EW_FB_Map_LS[a_idx][idx] = -999;
                entries_Full_Map_LS[a_idx][idx] = -999;
            
                entries_NS_Map_HS[a_idx][idx] = -999;
                entries_EW_Map_HS[a_idx][idx] = -999;
                entries_FB_Map_HS[a_idx][idx] = -999;
                entries_NS_EW_Map_HS[a_idx][idx] = -999;
                entries_NS_FB_Map_HS[a_idx][idx] = -999;
                entries_EW_FB_Map_HS[a_idx][idx] = -999;
                entries_Full_Map_HS[a_idx][idx] = -999;
            
            }
        
        events_LS = 0;
        events_HS = 0;
        
        seed = 0;
        seed_list_line = 0;
        
    }
    
    ~relative_fitResult() { }
    
};

class test_fitResult
{
public:
    
    Double_t inputAni[3];
    
    ULong64_t seed;
    ULong64_t seed_list_line;
    
    Double_t chi2_LS[4];
    Double_t ndf_LS[4];
    Double_t chi2_r_LS[4];
    Double_t fit_par_LS[4];
    Double_t fit_err_LS[4];
    Double_t delta_LS[4];
    
    Int_t theta_binHistos_LS;
    Int_t phi_binHistos_LS;
    
    ULong64_t events_LS;
    
    Double_t chi2_HS[4];
    Double_t ndf_HS[4];
    Double_t chi2_r_HS[4];
    Double_t fit_par_HS[4];
    Double_t fit_err_HS[4];
    Double_t delta_HS[4];
    
    Int_t theta_binHistos_HS;
    Int_t phi_binHistos_HS;
    
    ULong64_t events_HS;
    
    
    test_fitResult()
    {
        
        for(Int_t idx=0; idx<4; ++idx)
        {
            chi2_LS[idx] = -1;
            ndf_LS[idx] = -1;
            chi2_r_LS[idx] = -1;
            fit_par_LS[idx] = -1;
            fit_err_LS[idx] = -1;
            delta_LS[idx] = -1;
        
            chi2_HS[idx] = -1;
            ndf_HS[idx] = -1;
            chi2_r_HS[idx] = -1;
            fit_par_HS[idx] = -1;
            fit_err_HS[idx] = -1;
            delta_HS[idx] = -1;
        }
            
        theta_binHistos_LS = 0;
        phi_binHistos_LS = 0;
        theta_binHistos_HS = 0;
        phi_binHistos_HS = 0;
        
        for(Int_t idx = 0; idx < 3; ++idx)
            inputAni[idx] = -1;
        
        events_LS = 0;
        events_HS = 0;
        
        seed = 0;
        seed_list_line = 0;
        
    }
    
    ~test_fitResult() { }
    
};

///////////////////////////////////////////////////////////////////////////////////////////// Functions         ;-)

////////////////// Stuff functions //////////////////

extern std::string output_path_creator(
                                        const Int_t try_idx,
                                        const Int_t out_choose,
                                        Double_t NS_dipole = 0,
                                        Double_t EW_dipole = 0,
                                        Double_t FB_dipole = 0,
                                        bool DAMPE = false,
                                        bool test_fit = false
                                       );

extern void create_and_initialize_log(
                                        std::ofstream &log_file,
                                        Int_t s_idx,
                                        Int_t s_batch,
                                        Int_t n_try
                                      );

extern void log_file_init(
                            std::ofstream &out_file,
                            Int_t s_idx,
                            Int_t s_batch,
                            Int_t n_try
                          );

extern void TH2toTH1_obj(TH1D &Histo1D,TH2D &Histo2D);
extern void TH2toTH1_ptr(TH1D &Histo1D,TH2D* Histo2D);
extern void create_binning(Int_t n_bin_lat,Double_t lat_bin_min,Double_t lat_bin_max,Double_t* &binning,Bool_t cos_scale);
extern void generate_and_fit(std::ofstream &output_log_file,Int_t s_idx,Int_t s_batch);
extern void generate_and_fit_relative(std::ofstream &output_log_file,Int_t s_idx,Int_t s_batch);
extern void generate_and_fit_test(std::ofstream &output_log_file,Int_t s_idx,Int_t s_batch);
extern void create_tree_branches(TTree &fiTree,fitResult &tmp_fit);
extern void create_rel_tree_branches(TTree &fiTree,relative_fitResult &tmp_fit);
extern void create_test_tree_branches(TTree &fiTree,test_fitResult &tmp_fit);
extern void read_DAMPE_FullIso(TH2D &DAMPE_ReferenceMap);
extern void get_scaled_isotropic_DAMPE_maps(TH2D &DAMPE_ReferenceMap_LS,TH2D &DAMPE_ReferenceMap_HS);

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
                            bool DAMPE = false
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
                                        std::ofstream &output_log_file
                                    );

extern Double_t compute_integral(TH2D* histo);
extern void get_relative_histo(TH2D &relative_DAMPE_histo,TH2D* data_DAMPE_histo,TH2D &reference_DAMPE_histo);

////////////////// Analysis functions //////////////////

extern void templates_computation(
                                  std::ofstream &output_log_file,
                                  TH2D &DAMPE_ReferenceMap,
                                  std::string &templates_path,
                                  std::string &DAMPE_templates_path
                                  );

extern void deploy_simulation(
                                std::ofstream &output_log_file,
                                std::ifstream &inSeed,
                                Int_t s_idx,
                                Int_t s_batch,
                                std::string tmp_seed_str,
                                UInt_t tmp_seed,
                                TTree &fiTree,
                                fitResult &tmp_fit
                              );

extern void deploy_relative_simulation(
                                        std::ofstream &output_log_file,
                                        std::ifstream &inSeed,
                                        Int_t s_idx,
                                        Int_t s_batch,
                                        std::string tmp_seed_str,
                                        UInt_t tmp_seed,
                                        TTree &fiTree,
                                        relative_fitResult &tmp_fit,
                                        TH2D &DAMPE_ReferenceMap_LS,
                                        TH2D &DAMPE_ReferenceMap_HS
                                    );

extern void deploy_test_simulation(
                                    std::ofstream &output_log_file,
                                    std::ifstream &inSeed,
                                    Int_t s_idx,
                                    Int_t s_batch,
                                    std::string tmp_seed_str,
                                    UInt_t tmp_seed,
                                    TTree &fiTree,
                                    test_fitResult &tmp_fit
                                   );

extern void generate_LS_templates(
                                  TH2D &Template_Iso_LS,
                                  TH2D &Template_AniNS_LS,
                                  TH2D &Template_AniEW_LS,
                                  TH2D &Template_AniFB_LS,
                                  std::ofstream &output_log_file,
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
                                  std::ofstream &output_log_file,
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
                             UInt_t tmp_seed
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
                             UInt_t tmp_seed
                             );


extern void generate_LS_example_data(
                                        Double_t NS_anisotropy,
                                        Double_t EW_anisotropy,
                                        Double_t FB_anisotropy,
                                        TH2D* Data_Iso_LS,
                                        TH2D* Data_AniNS_LS,
                                        TH2D* Data_AniEW_LS,
                                        TH2D* Data_AniFB_LS,
                                        TH2D &Template_Iso_LS,
                                        TH2D &Template_AniNS_LS,
                                        TH2D &Template_AniEW_LS,
                                        TH2D &Template_AniFB_LS,
                                        std::ofstream &output_log_file,
                                        UInt_t tmp_seed
                                     );

extern void generate_HS_example_data(
                                        Double_t NS_anisotropy,
                                        Double_t EW_anisotropy,
                                        Double_t FB_anisotropy,
                                        TH2D* Data_Iso_HS,
                                        TH2D* Data_AniNS_HS,
                                        TH2D* Data_AniEW_HS,
                                        TH2D* Data_AniFB_HS,
                                        TH2D &Template_Iso_HS,
                                        TH2D &Template_AniNS_HS,
                                        TH2D &Template_AniEW_HS,
                                        TH2D &Template_AniFB_HS,
                                        std::ofstream &output_log_file,
                                        UInt_t tmp_seed
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
                                   UInt_t tmp_seed
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
                                   UInt_t tmp_seed
                                   );

extern void getPull(TH1* Data,TH1* Templates[],Double_t res[],TH1D &hPull,bool relative_fit=false,bool ex_fit=false);

extern void generate_DAMPE_templates(
                                     TH2D &DAMPE_FullIso,
                                     TH2D &DAMPE_Template_Iso,
                                     TH2D &DAMPE_Template_AniNS,
                                     TH2D &DAMPE_Template_AniEW,
                                     TH2D &DAMPE_Template_AniFB,
                                     TH1D &DAMPE_Template_hwNS,
                                     TH1D &DAMPE_Template_hwEW,
                                     TH1D &DAMPE_Template_hwFB,
                                     std::ofstream &log_file,
                                     TF2 &dI,
                                     TF2 &dNS,
                                     TF2 &dEW,
                                     TF2 &dFB,
                                     TCanvas &FCanvas
                                     );

extern void components_analysis(
                                    std::vector<TH1D*> &TemplatesProjections,
                                    std::vector<TH1D*> &DataProjections,
                                    double resfullFitResults_HS[][8],
                                    std::ofstream &log_file,
                                    Double_t NS_anisotropy,
                                    Double_t EW_anisotropy,
                                    Double_t FB_anisotropy,
                                    Int_t try_idx
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

extern double compute_ani_level(double covMatrix[][4],double parameters[],double par_err[],std::ofstream &log_file);
extern double compute_relative_ani_level(double covMatrix[][3],double parameters[],double par_err[],std::ofstream &log_file);

extern void allSky_singleTry_fit(
                                    std::ofstream &output_log_file,
                                    TTree &fiTree,
                                    fitResult &tmp_fit,
                                    UInt_t tmp_seed,
                                    Int_t itry,
                                    Double_t NS_anisotropy,
                                    Double_t EW_anisotropy,
                                    Double_t FB_anisotropy,
                                    Int_t idx_ani,
                                    Int_t try_idx,
                                    UInt_t seed_line
                                 );

extern void DAMPE_singleTry_fit(
                                    std::ofstream &output_log_file,
                                    TTree &fiTree,
                                    fitResult &tmp_fit,
                                    UInt_t tmp_seed,
                                    Int_t itry,
                                    Double_t NS_anisotropy,
                                    Double_t EW_anisotropy,
                                    Double_t FB_anisotropy,
                                    Int_t idx_ani,
                                    Int_t try_idx,
                                    UInt_t seed_line
                                );

extern void DAMPE_relative_singleTry_fit(
                                            std::ofstream &output_log_file,
                                            TTree &fiTree,
                                            relative_fitResult &tmp_fit,
                                            UInt_t tmp_seed,
                                            Int_t itry,
                                            Double_t NS_anisotropy,
                                            Double_t EW_anisotropy,
                                            Double_t FB_anisotropy,
                                            Int_t idx_ani,
                                            TH2D &DAMPE_ReferenceMap_LS,
                                            TH2D &DAMPE_ReferenceMap_HS,
                                            Int_t try_idx,
                                            UInt_t seed_line
                                         );

extern void AllSky_test_simulation(
                                    std::ofstream &output_log_file,
                                    TTree &fiTree,
                                    test_fitResult &tmp_fit,
                                    UInt_t tmp_seed,
                                    Int_t itry,
                                    Double_t NS_anisotropy,
                                    Double_t EW_anisotropy,
                                    Double_t FB_anisotropy,
                                    Int_t idx_ani,
                                    Int_t s_idx,
                                    UInt_t seed_line
                                   );

////////////////// RooFit TemplateFit functions //////////////////

extern void ZeroRooFitVerbosity();
extern void RemoveZeroes(TH1* h);
extern void ResetRooFitVerbosity();
extern void ZeroRooFitVerbosity();
extern void fcnchisq(int& npar, double* deriv, double& f, double par[], int flag);
extern void fcnlike(int& npar, double* deriv, double& f, double par[], int flag);
extern double func(int npar, double par[], double _dochisqfit);

extern TH1* TemplateFitRF(
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

extern TH1* TemplateFitBH(
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
                            std::ofstream &log_file,
                            fitResult &tmp_fit,
                            Int_t fit_element,
                            bool HS
                          );

extern TH1* TemplateFitBH_rel(
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
                                std::ofstream &log_file,
                                relative_fitResult &tmp_fit,
                                Int_t fit_element,
                                bool HS
                              );

extern TH1* TemplateFitBH_test(
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
                                std::ofstream &log_file,
                                test_fitResult &tmp_fit,
                                Int_t fit_element,
                                bool HS
                               );


#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <fstream>

#include "TFile.h"
#include "TH1D.h"
#include "TMath.h"
#include "TAxis.h"
#include "TH2D.h"
#include "TF2.h"

const static time_t time_stamp=time(0);

const static std::string path_location = "/storage/gpfs_data/dampe/users/ecatanzani/MyRepos/DAMPE/MapsFitting/AnaliticalTemplates/SubmitJobs/ToNode/ExeSW/assets/computeTemplates/logs/";
const static std::string templates_path = "/storage/gpfs_data/dampe/users/ecatanzani/MyRepos/DAMPE/MapsFitting/AnaliticalTemplates/SubmitJobs/ToNode/ExeSW/assets/computeTemplates/results/AllSkyTemplates.root";
const static std::string DAMPE_templates_path = "/storage/gpfs_data/dampe/users/ecatanzani/MyRepos/DAMPE/MapsFitting/AnaliticalTemplates/SubmitJobs/ToNode/ExeSW/assets/computeTemplates/results/DAMPETemplates.root";

//////////////////////// Custom functions

extern void read_DAMPE_FullIso(TH2D &DAMPE_FullIso);
extern void create_binning(Int_t n_bin_lat,Double_t lat_bin_min,Double_t lat_bin_max,Double_t* &binning,Bool_t cos_scale);

extern void get_scaled_isotropic_DAMPE_maps(TH2D &DAMPE_ReferenceMap_LS,TH2D &DAMPE_ReferenceMap_HS,const std::string DAMPE_Iso_scaled_Maps);

extern void templates_computation(
                                    std::ofstream &output_log_file,
                                    TH2D &DAMPE_ReferenceMap_LS,
                                    TH2D &DAMPE_ReferenceMap_HS
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



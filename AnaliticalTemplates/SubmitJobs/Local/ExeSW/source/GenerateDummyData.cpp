
#include "MyHead.h"

void generate_LS_example_data(
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
                              )

{
    
    ///////////// Building data templates
    
    TH2D* pTemplate_Iso_LS = &Template_Iso_LS;
    TH2D* tmpTemplate_AniNS_LS = (TH2D*) Template_AniNS_LS.Clone("tmpTemplate_AniNS_LS");
    TH2D* tmpTemplate_AniEW_LS = (TH2D*) Template_AniEW_LS.Clone("tmpTemplate_AniEW_LS");
    TH2D* tmpTemplate_AniFB_LS = (TH2D*) Template_AniFB_LS.Clone("tmpTemplate_AniFB_LS");
    
    tmpTemplate_AniNS_LS->Scale(NS_anisotropy);
    tmpTemplate_AniEW_LS->Scale(EW_anisotropy);
    tmpTemplate_AniFB_LS->Scale(FB_anisotropy);
    
    ///////////// Filling maps...
    
    gRandom->SetSeed(tmp_seed);
    
    Data_Iso_LS->FillRandom(pTemplate_Iso_LS,gRandom->Poisson(data_all_sky_LS_events));
    Data_AniNS_LS->FillRandom(tmpTemplate_AniNS_LS,gRandom->Poisson(data_all_sky_LS_events));
    Data_AniEW_LS->FillRandom(tmpTemplate_AniEW_LS,gRandom->Poisson(data_all_sky_LS_events));
    Data_AniFB_LS->FillRandom(tmpTemplate_AniFB_LS,gRandom->Poisson(data_all_sky_LS_events));
    
}

void generate_HS_example_data(
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
                              )

{
    
    ///////////// Building data templates
    
    TH2D* pTemplate_Iso_HS = &Template_Iso_HS;
    TH2D* tmpTemplate_AniNS_HS = (TH2D*) Template_AniNS_HS.Clone("tmpTemplate_AniNS_HS");
    TH2D* tmpTemplate_AniEW_HS = (TH2D*) Template_AniEW_HS.Clone("tmpTemplate_AniEW_HS");
    TH2D* tmpTemplate_AniFB_HS = (TH2D*) Template_AniFB_HS.Clone("tmpTemplate_AniFB_HS");
    
    tmpTemplate_AniNS_HS->Scale(NS_anisotropy);
    tmpTemplate_AniEW_HS->Scale(EW_anisotropy);
    tmpTemplate_AniFB_HS->Scale(FB_anisotropy);
    
    ///////////// Filling maps...
    
    gRandom->SetSeed(tmp_seed);
    
    Data_Iso_HS->FillRandom(pTemplate_Iso_HS,gRandom->Poisson(data_all_sky_HS_events));
    Data_AniNS_HS->FillRandom(tmpTemplate_AniNS_HS,gRandom->Poisson(data_all_sky_HS_events));
    Data_AniEW_HS->FillRandom(tmpTemplate_AniEW_HS,gRandom->Poisson(data_all_sky_HS_events));
    Data_AniFB_HS->FillRandom(tmpTemplate_AniFB_HS,gRandom->Poisson(data_all_sky_HS_events));
    
}

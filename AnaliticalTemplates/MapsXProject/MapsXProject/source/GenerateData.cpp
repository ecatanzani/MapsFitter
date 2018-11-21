
#include "MyHead.h"

void generate_LS_data(
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
                      )

{
    
    std::cout << "\n\n Building LS Data maps...";
    
    ///////////// Building data templates
    
    TH2D Template_Data_AniNS_LS;
    TH2D Template_Data_AniEW_LS;
    TH2D Template_Data_AniFB_LS;
    TH2D Template_MixedData_NS_EW_LS;
    TH2D Template_MixedData_NS_FB_LS;
    TH2D Template_MixedData_EW_FB_LS;
    TH2D Template_FullMixedData_LS;
    
    compute_data_template(
                            Template_Data_AniNS_LS,
                            Template_Iso_LS,
                            Template_AniNS_LS,
                            Template_AniEW_LS,
                            Template_AniFB_LS,
                            NS_anisotropy,
                            0,
                            0,
                            "Template_Data_AniNS_LS"
                          );
    
    compute_data_template(
                            Template_Data_AniEW_LS,
                            Template_Iso_LS,
                            Template_AniNS_LS,
                            Template_AniEW_LS,
                            Template_AniFB_LS,
                            0,
                            EW_anisotropy,
                            0,
                            "Template_Data_AniEW_LS"
                          );
    
    compute_data_template(
                            Template_Data_AniFB_LS,
                            Template_Iso_LS,
                            Template_AniNS_LS,
                            Template_AniEW_LS,
                            Template_AniFB_LS,
                            0,
                            0,
                            FB_anisotropy,
                            "Template_Data_AniFB_LS"
                          );
    
    compute_data_template(
                            Template_MixedData_NS_EW_LS,
                            Template_Iso_LS,
                            Template_AniNS_LS,
                            Template_AniEW_LS,
                            Template_AniFB_LS,
                            NS_anisotropy,
                            EW_anisotropy,
                            0,
                            "Template_MixedData_NS_EW_LS"
                          );
    
    compute_data_template(
                            Template_MixedData_NS_FB_LS,
                            Template_Iso_LS,
                            Template_AniNS_LS,
                            Template_AniEW_LS,
                            Template_AniFB_LS,
                            NS_anisotropy,
                            0,
                            FB_anisotropy,
                            "Template_MixedData_NS_FB_LS"
                          );
    
    compute_data_template(
                            Template_MixedData_EW_FB_LS,
                            Template_Iso_LS,
                            Template_AniNS_LS,
                            Template_AniEW_LS,
                            Template_AniFB_LS,
                            0,
                            EW_anisotropy,
                            FB_anisotropy,
                            "Template_MixedData_EW_FB_LS"
                          );
    
    compute_data_template(
                            Template_FullMixedData_LS,
                            Template_Iso_LS,
                            Template_AniNS_LS,
                            Template_AniEW_LS,
                            Template_AniFB_LS,
                            NS_anisotropy,
                            EW_anisotropy,
                            FB_anisotropy,
                            "Template_FullMixedData_LS"
                          );
    
    
    ///////////// Linking to the templates
    
    TH2D* pTemplate_Data_Iso_LS = &Template_Iso_LS;
    TH2D* pTemplate_Data_AniNS_LS = &Template_Data_AniNS_LS;
    TH2D* pTemplate_Data_AniEW_LS = &Template_Data_AniEW_LS;
    TH2D* pTemplate_Data_AniFB_LS = &Template_Data_AniFB_LS;
    TH2D* pTemplate_MixedData_NS_EW_LS = &Template_MixedData_NS_EW_LS;
    TH2D* pTemplate_MixedData_NS_FB_LS = &Template_MixedData_NS_FB_LS;
    TH2D* pTemplate_MixedData_EW_FB_LS = &Template_MixedData_EW_FB_LS;
    TH2D* pTemplate_FullMixedData_LS = &Template_FullMixedData_LS;
    
    ///////////// Filling maps...
    
    gRandom->SetSeed(tmp_seed);
    
    Data_Iso_LS->FillRandom(pTemplate_Data_Iso_LS,gRandom->Poisson(data_LS_events));
    Data_AniNS_LS->FillRandom(pTemplate_Data_AniNS_LS,gRandom->Poisson(data_LS_events));
    Data_AniEW_LS->FillRandom(pTemplate_Data_AniEW_LS,gRandom->Poisson(data_LS_events));
    Data_AniFB_LS->FillRandom(pTemplate_Data_AniFB_LS,gRandom->Poisson(data_LS_events));
    MixedData_NS_EW_LS->FillRandom(pTemplate_MixedData_NS_EW_LS,gRandom->Poisson(data_LS_events));
    MixedData_NS_FB_LS->FillRandom(pTemplate_MixedData_NS_FB_LS,gRandom->Poisson(data_LS_events));
    MixedData_EW_FB_LS->FillRandom(pTemplate_MixedData_EW_FB_LS,gRandom->Poisson(data_LS_events));
    FullMixedData_LS->FillRandom(pTemplate_FullMixedData_LS,gRandom->Poisson(data_LS_events));
    
}


void generate_HS_data(
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
                      )

{
    
    std::cout << "\n\n Building LS Data maps...";
    
    ///////////// Building data templates
    
    TH2D Template_Data_AniNS_HS;
    TH2D Template_Data_AniEW_HS;
    TH2D Template_Data_AniFB_HS;
    TH2D Template_MixedData_NS_EW_HS;
    TH2D Template_MixedData_NS_FB_HS;
    TH2D Template_MixedData_EW_FB_HS;
    TH2D Template_FullMixedData_HS;
    
    compute_data_template(
                            Template_Data_AniNS_HS,
                            Template_Iso_HS,
                            Template_AniNS_HS,
                            Template_AniEW_HS,
                            Template_AniFB_HS,
                            NS_anisotropy,
                            0,
                            0,
                            "Template_Data_AniNS_HS"
                          );
    
    compute_data_template(
                            Template_Data_AniEW_HS,
                            Template_Iso_HS,
                            Template_AniNS_HS,
                            Template_AniEW_HS,
                            Template_AniFB_HS,
                            0,
                            EW_anisotropy,
                            0,
                            "Template_Data_AniEW_HS"
                          );
    
    compute_data_template(
                            Template_Data_AniFB_HS,
                            Template_Iso_HS,
                            Template_AniNS_HS,
                            Template_AniEW_HS,
                            Template_AniFB_HS,
                            0,
                            0,
                            FB_anisotropy,
                            "Template_Data_AniFB_HS"
                          );
    
    compute_data_template(
                            Template_MixedData_NS_EW_HS,
                            Template_Iso_HS,
                            Template_AniNS_HS,
                            Template_AniEW_HS,
                            Template_AniFB_HS,
                            NS_anisotropy,
                            EW_anisotropy,
                            0,
                            "Template_MixedData_NS_EW_HS"
                          );
    
    compute_data_template(
                            Template_MixedData_NS_FB_HS,
                            Template_Iso_HS,
                            Template_AniNS_HS,
                            Template_AniEW_HS,
                            Template_AniFB_HS,
                            NS_anisotropy,
                            0,
                            FB_anisotropy,
                            "Template_MixedData_NS_FB_HS"
                          );
    
    compute_data_template(
                            Template_MixedData_EW_FB_HS,
                            Template_Iso_HS,
                            Template_AniNS_HS,
                            Template_AniEW_HS,
                            Template_AniFB_HS,
                            0,
                            EW_anisotropy,
                            FB_anisotropy,
                            "Template_MixedData_EW_FB_HS"
                          );
    
    compute_data_template(
                            Template_FullMixedData_HS,
                            Template_Iso_HS,
                            Template_AniNS_HS,
                            Template_AniEW_HS,
                            Template_AniFB_HS,
                            NS_anisotropy,
                            EW_anisotropy,
                            FB_anisotropy,
                            "Template_FullMixedData_HS"
                          );
    
    
    ///////////// Linking to the templates
    
    TH2D* pTemplate_Data_Iso_HS = &Template_Iso_HS;
    TH2D* pTemplate_Data_AniNS_HS = &Template_Data_AniNS_HS;
    TH2D* pTemplate_Data_AniEW_HS = &Template_Data_AniEW_HS;
    TH2D* pTemplate_Data_AniFB_HS = &Template_Data_AniFB_HS;
    TH2D* pTemplate_MixedData_NS_EW_HS = &Template_MixedData_NS_EW_HS;
    TH2D* pTemplate_MixedData_NS_FB_HS = &Template_MixedData_NS_FB_HS;
    TH2D* pTemplate_MixedData_EW_FB_HS = &Template_MixedData_EW_FB_HS;
    TH2D* pTemplate_FullMixedData_HS = &Template_FullMixedData_HS;
    
    ///////////// Filling maps...
    
    gRandom->SetSeed(tmp_seed);
    
    Data_Iso_HS->FillRandom(pTemplate_Data_Iso_HS,gRandom->Poisson(data_HS_events));
    Data_AniNS_HS->FillRandom(pTemplate_Data_AniNS_HS,gRandom->Poisson(data_HS_events));
    Data_AniEW_HS->FillRandom(pTemplate_Data_AniEW_HS,gRandom->Poisson(data_HS_events));
    Data_AniFB_HS->FillRandom(pTemplate_Data_AniFB_HS,gRandom->Poisson(data_HS_events));
    MixedData_NS_EW_HS->FillRandom(pTemplate_MixedData_NS_EW_HS,gRandom->Poisson(data_HS_events));
    MixedData_NS_FB_HS->FillRandom(pTemplate_MixedData_NS_FB_HS,gRandom->Poisson(data_HS_events));
    MixedData_EW_FB_HS->FillRandom(pTemplate_MixedData_EW_FB_HS,gRandom->Poisson(data_HS_events));
    FullMixedData_HS->FillRandom(pTemplate_FullMixedData_HS,gRandom->Poisson(data_HS_events));
    
    
}

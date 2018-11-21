
#include "MyHead.h"

void generate_DAMPE_LS_data(
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
                            )

{
    
    std::cout << "\n\n Building LS Data maps...\n";
    
    ///////////// Building data templates
    
    TH2D DAMPE_Template_Data_AniNS_LS;
    TH2D DAMPE_Template_Data_AniEW_LS;
    TH2D DAMPE_Template_Data_AniFB_LS;
    TH2D DAMPE_Template_MixedData_NS_EW_LS;
    TH2D DAMPE_Template_MixedData_NS_FB_LS;
    TH2D DAMPE_Template_MixedData_EW_FB_LS;
    TH2D DAMPE_Template_FullMixedData_LS;
    
    compute_data_template(
                              DAMPE_Template_Data_AniNS_LS,
                              DAMPE_Template_Iso_LS,
                              DAMPE_Template_AniNS_LS,
                              DAMPE_Template_AniEW_LS,
                              DAMPE_Template_AniFB_LS,
                              NS_anisotropy,
                              0,
                              0,
                              "DAMPE_Template_Data_AniNS_LS"
                          );
    
    compute_data_template(
                              DAMPE_Template_Data_AniEW_LS,
                              DAMPE_Template_Iso_LS,
                              DAMPE_Template_AniNS_LS,
                              DAMPE_Template_AniEW_LS,
                              DAMPE_Template_AniFB_LS,
                              0,
                              EW_anisotropy,
                              0,
                              "DAMPE_Template_Data_AniEW_LS"
                          );
    
    compute_data_template(
                              DAMPE_Template_Data_AniFB_LS,
                              DAMPE_Template_Iso_LS,
                              DAMPE_Template_AniNS_LS,
                              DAMPE_Template_AniEW_LS,
                              DAMPE_Template_AniFB_LS,
                              0,
                              0,
                              FB_anisotropy,
                              "DAMPE_Template_Data_AniFB_LS"
                          );
    
    compute_data_template(
                              DAMPE_Template_MixedData_NS_EW_LS,
                              DAMPE_Template_Iso_LS,
                              DAMPE_Template_AniNS_LS,
                              DAMPE_Template_AniEW_LS,
                              DAMPE_Template_AniFB_LS,
                              NS_anisotropy,
                              EW_anisotropy,
                              0,
                              "DAMPE_Template_MixedData_NS_EW_LS"
                          );
    
    compute_data_template(
                              DAMPE_Template_MixedData_NS_FB_LS,
                              DAMPE_Template_Iso_LS,
                              DAMPE_Template_AniNS_LS,
                              DAMPE_Template_AniEW_LS,
                              DAMPE_Template_AniFB_LS,
                              NS_anisotropy,
                              0,
                              FB_anisotropy,
                              "DAMPE_Template_MixedData_NS_FB_LS"
                          );
    
    compute_data_template(
                              DAMPE_Template_MixedData_EW_FB_LS,
                              DAMPE_Template_Iso_LS,
                              DAMPE_Template_AniNS_LS,
                              DAMPE_Template_AniEW_LS,
                              DAMPE_Template_AniFB_LS,
                              0,
                              EW_anisotropy,
                              FB_anisotropy,
                              "DAMPE_Template_MixedData_EW_FB_LS"
                          );
    
    compute_data_template(
                              DAMPE_Template_FullMixedData_LS,
                              DAMPE_Template_Iso_LS,
                              DAMPE_Template_AniNS_LS,
                              DAMPE_Template_AniEW_LS,
                              DAMPE_Template_AniFB_LS,
                              NS_anisotropy,
                              EW_anisotropy,
                              FB_anisotropy,
                              "DAMPE_Template_FullMixedData_LS"
                          );
    
    
    ///////////// Linking to the templates
    
    TH2D* pDAMPE_Template_Data_Iso_LS = &DAMPE_Template_Iso_LS;
    TH2D* pDAMPE_Template_Data_AniNS_LS = &DAMPE_Template_Data_AniNS_LS;
    TH2D* pDAMPE_Template_Data_AniEW_LS = &DAMPE_Template_Data_AniEW_LS;
    TH2D* pDAMPE_Template_Data_AniFB_LS = &DAMPE_Template_Data_AniFB_LS;
    TH2D* pDAMPE_Template_MixedData_NS_EW_LS = &DAMPE_Template_MixedData_NS_EW_LS;
    TH2D* pDAMPE_Template_MixedData_NS_FB_LS = &DAMPE_Template_MixedData_NS_FB_LS;
    TH2D* pDAMPE_Template_MixedData_EW_FB_LS = &DAMPE_Template_MixedData_EW_FB_LS;
    TH2D* pDAMPE_Template_FullMixedData_LS = &DAMPE_Template_FullMixedData_LS;
    
    ///////////// Filling maps...
    
    gRandom->SetSeed(tmp_seed);
    
    DAMPE_Data_Iso_LS->FillRandom(pDAMPE_Template_Data_Iso_LS,gRandom->Poisson(data_LS_events));
    DAMPE_Data_AniNS_LS->FillRandom(pDAMPE_Template_Data_AniNS_LS,gRandom->Poisson(data_LS_events));
    DAMPE_Data_AniEW_LS->FillRandom(pDAMPE_Template_Data_AniEW_LS,gRandom->Poisson(data_LS_events));
    DAMPE_Data_AniFB_LS->FillRandom(pDAMPE_Template_Data_AniFB_LS,gRandom->Poisson(data_LS_events));
    DAMPE_MixedData_NS_EW_LS->FillRandom(pDAMPE_Template_MixedData_NS_EW_LS,gRandom->Poisson(data_LS_events));
    DAMPE_MixedData_NS_FB_LS->FillRandom(pDAMPE_Template_MixedData_NS_FB_LS,gRandom->Poisson(data_LS_events));
    DAMPE_MixedData_EW_FB_LS->FillRandom(pDAMPE_Template_MixedData_EW_FB_LS,gRandom->Poisson(data_LS_events));
    DAMPE_FullMixedData_LS->FillRandom(pDAMPE_Template_FullMixedData_LS,gRandom->Poisson(data_LS_events));
    
    
    /*
     for(Int_t bX=1; bX <= DAMPE_Data_Iso_LS->GetNbinsX(); ++bX)
     {
     for(Int_t bY=1; bY <= DAMPE_Data_Iso_LS->GetNbinsY(); ++bY)
     {
     DAMPE_Data_Iso_LS->SetBinError(bX,bY,DAMPE_Data_Iso_LS->GetBinContent(bX,bY)*(1/TMath::Sqrt(pDAMPE_Template_Data_Iso_LS->GetBinContent(bX,bY))));
     DAMPE_Data_AniNS_LS->SetBinError(bX,bY,DAMPE_Data_AniNS_LS->GetBinContent(bX,bY)*(1/TMath::Sqrt(pDAMPE_Template_Data_AniNS_LS->GetBinContent(bX,bY))));
     DAMPE_Data_AniEW_LS->SetBinError(bX,bY,DAMPE_Data_AniEW_LS->GetBinContent(bX,bY)*(1/TMath::Sqrt(pDAMPE_Template_Data_AniEW_LS->GetBinContent(bX,bY))));
     DAMPE_Data_AniFB_LS->SetBinError(bX,bY,DAMPE_Data_AniFB_LS->GetBinContent(bX,bY)*(1/TMath::Sqrt(pDAMPE_Template_Data_AniFB_LS->GetBinContent(bX,bY))));
     DAMPE_MixedData_NS_EW_LS->SetBinError(bX,bY,DAMPE_MixedData_NS_EW_LS->GetBinContent(bX,bY)*(1/TMath::Sqrt(pDAMPE_Template_MixedData_NS_EW_LS->GetBinContent(bX,bY))));
     DAMPE_MixedData_NS_FB_LS->SetBinError(bX,bY,DAMPE_MixedData_NS_FB_LS->GetBinContent(bX,bY)*(1/TMath::Sqrt(pDAMPE_Template_MixedData_NS_FB_LS->GetBinContent(bX,bY))));
     DAMPE_MixedData_EW_FB_LS->SetBinError(bX,bY,DAMPE_MixedData_EW_FB_LS->GetBinContent(bX,bY)*(1/TMath::Sqrt(pDAMPE_Template_MixedData_EW_FB_LS->GetBinContent(bX,bY))));
     DAMPE_FullMixedData_LS->SetBinError(bX,bY,DAMPE_FullMixedData_LS->GetBinContent(bX,bY)*(1/TMath::Sqrt(pDAMPE_Template_FullMixedData_LS->GetBinContent(bX,bY))));
     }
     }
     */
    
}


void generate_DAMPE_HS_data(
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
                            )

{
    std::cout << "\n Building HS Data maps... \n\n";
    
    ///////////// Building data templates
    
    TH2D DAMPE_Template_Data_AniNS_HS;
    TH2D DAMPE_Template_Data_AniEW_HS;
    TH2D DAMPE_Template_Data_AniFB_HS;
    TH2D DAMPE_Template_MixedData_NS_EW_HS;
    TH2D DAMPE_Template_MixedData_NS_FB_HS;
    TH2D DAMPE_Template_MixedData_EW_FB_HS;
    TH2D DAMPE_Template_FullMixedData_HS;
    
    compute_data_template(
                              DAMPE_Template_Data_AniNS_HS,
                              DAMPE_Template_Iso_HS,
                              DAMPE_Template_AniNS_HS,
                              DAMPE_Template_AniEW_HS,
                              DAMPE_Template_AniFB_HS,
                              NS_anisotropy,
                              0,
                              0,
                              "DAMPE_Template_Data_AniNS_HS"
                          );
    
    compute_data_template(
                              DAMPE_Template_Data_AniEW_HS,
                              DAMPE_Template_Iso_HS,
                              DAMPE_Template_AniNS_HS,
                              DAMPE_Template_AniEW_HS,
                              DAMPE_Template_AniFB_HS,
                              0,
                              EW_anisotropy,
                              0,
                              "DAMPE_Template_Data_AniEW_HS"
                          );
    
    compute_data_template(
                              DAMPE_Template_Data_AniFB_HS,
                              DAMPE_Template_Iso_HS,
                              DAMPE_Template_AniNS_HS,
                              DAMPE_Template_AniEW_HS,
                              DAMPE_Template_AniFB_HS,
                              0,
                              0,
                              FB_anisotropy,
                              "DAMPE_Template_Data_AniFB_HS"
                          );
    
    compute_data_template(
                              DAMPE_Template_MixedData_NS_EW_HS,
                              DAMPE_Template_Iso_HS,
                              DAMPE_Template_AniNS_HS,
                              DAMPE_Template_AniEW_HS,
                              DAMPE_Template_AniFB_HS,
                              NS_anisotropy,
                              EW_anisotropy,
                              0,
                              "DAMPE_Template_MixedData_NS_EW_HS"
                          );
    
    compute_data_template(
                              DAMPE_Template_MixedData_NS_FB_HS,
                              DAMPE_Template_Iso_HS,
                              DAMPE_Template_AniNS_HS,
                              DAMPE_Template_AniEW_HS,
                              DAMPE_Template_AniFB_HS,
                              NS_anisotropy,
                              0,
                              FB_anisotropy,
                              "DAMPE_Template_MixedData_NS_FB_HS"
                          );
    
    compute_data_template(
                              DAMPE_Template_MixedData_EW_FB_HS,
                              DAMPE_Template_Iso_HS,
                              DAMPE_Template_AniNS_HS,
                              DAMPE_Template_AniEW_HS,
                              DAMPE_Template_AniFB_HS,
                              0,
                              EW_anisotropy,
                              FB_anisotropy,
                              "DAMPE_Template_MixedData_EW_FB_HS"
                          );
    
    compute_data_template(
                              DAMPE_Template_FullMixedData_HS,
                              DAMPE_Template_Iso_HS,
                              DAMPE_Template_AniNS_HS,
                              DAMPE_Template_AniEW_HS,
                              DAMPE_Template_AniFB_HS,
                              NS_anisotropy,
                              EW_anisotropy,
                              FB_anisotropy,
                              "DAMPE_Template_FullMixedData_HS"
                          );
    
    
    ///////////// Linking to the templates
    
    TH2D* pDAMPE_Template_Data_Iso_HS = &DAMPE_Template_Iso_HS;
    TH2D* pDAMPE_Template_Data_AniNS_HS = &DAMPE_Template_Data_AniNS_HS;
    TH2D* pDAMPE_Template_Data_AniEW_HS = &DAMPE_Template_Data_AniEW_HS;
    TH2D* pDAMPE_Template_Data_AniFB_HS = &DAMPE_Template_Data_AniFB_HS;
    TH2D* pDAMPE_Template_MixedData_NS_EW_HS = &DAMPE_Template_MixedData_NS_EW_HS;
    TH2D* pDAMPE_Template_MixedData_NS_FB_HS = &DAMPE_Template_MixedData_NS_FB_HS;
    TH2D* pDAMPE_Template_MixedData_EW_FB_HS = &DAMPE_Template_MixedData_EW_FB_HS;
    TH2D* pDAMPE_Template_FullMixedData_HS = &DAMPE_Template_FullMixedData_HS;
    
    ///////////// Filling maps...
    
    gRandom->SetSeed(tmp_seed);
    
    DAMPE_Data_Iso_HS->FillRandom(pDAMPE_Template_Data_Iso_HS,gRandom->Poisson(data_HS_events));
    DAMPE_Data_AniNS_HS->FillRandom(pDAMPE_Template_Data_AniNS_HS,gRandom->Poisson(data_HS_events));
    DAMPE_Data_AniEW_HS->FillRandom(pDAMPE_Template_Data_AniEW_HS,gRandom->Poisson(data_HS_events));
    DAMPE_Data_AniFB_HS->FillRandom(pDAMPE_Template_Data_AniFB_HS,gRandom->Poisson(data_HS_events));
    DAMPE_MixedData_NS_EW_HS->FillRandom(pDAMPE_Template_MixedData_NS_EW_HS,gRandom->Poisson(data_HS_events));
    DAMPE_MixedData_NS_FB_HS->FillRandom(pDAMPE_Template_MixedData_NS_FB_HS,gRandom->Poisson(data_HS_events));
    DAMPE_MixedData_EW_FB_HS->FillRandom(pDAMPE_Template_MixedData_EW_FB_HS,gRandom->Poisson(data_HS_events));
    DAMPE_FullMixedData_HS->FillRandom(pDAMPE_Template_FullMixedData_HS,gRandom->Poisson(data_HS_events));
    
    /*
     for(Int_t bX=1; bX <= DAMPE_Data_Iso_HS->GetNbinsX(); ++bX)
     {
     for(Int_t bY=1; bY <= DAMPE_Data_Iso_HS->GetNbinsY(); ++bY)
     {
     DAMPE_Data_Iso_HS->SetBinError(bX,bY,DAMPE_Data_Iso_HS->GetBinContent(bX,bY)*(1/TMath::Sqrt(pDAMPE_Template_Data_Iso_HS->GetBinContent(bX,bY))));
     DAMPE_Data_AniNS_HS->SetBinError(bX,bY,DAMPE_Data_AniNS_HS->GetBinContent(bX,bY)*(1/TMath::Sqrt(pDAMPE_Template_Data_AniNS_HS->GetBinContent(bX,bY))));
     DAMPE_Data_AniEW_HS->SetBinError(bX,bY,DAMPE_Data_AniEW_HS->GetBinContent(bX,bY)*(1/TMath::Sqrt(pDAMPE_Template_Data_AniEW_HS->GetBinContent(bX,bY))));
     DAMPE_Data_AniFB_HS->SetBinError(bX,bY,DAMPE_Data_AniFB_HS->GetBinContent(bX,bY)*(1/TMath::Sqrt(pDAMPE_Template_Data_AniFB_HS->GetBinContent(bX,bY))));
     DAMPE_MixedData_NS_EW_HS->SetBinError(bX,bY,DAMPE_MixedData_NS_EW_HS->GetBinContent(bX,bY)*(1/TMath::Sqrt(pDAMPE_Template_MixedData_NS_EW_HS->GetBinContent(bX,bY))));
     DAMPE_MixedData_NS_FB_HS->SetBinError(bX,bY,DAMPE_MixedData_NS_FB_HS->GetBinContent(bX,bY)*(1/TMath::Sqrt(pDAMPE_Template_MixedData_NS_FB_HS->GetBinContent(bX,bY))));
     DAMPE_MixedData_EW_FB_HS->SetBinError(bX,bY,DAMPE_MixedData_EW_FB_HS->GetBinContent(bX,bY)*(1/TMath::Sqrt(pDAMPE_Template_MixedData_EW_FB_HS->GetBinContent(bX,bY))));
     DAMPE_FullMixedData_HS->SetBinError(bX,bY,DAMPE_FullMixedData_HS->GetBinContent(bX,bY)*(1/TMath::Sqrt(pDAMPE_Template_FullMixedData_HS->GetBinContent(bX,bY))));
     }
     }
     */
    
}

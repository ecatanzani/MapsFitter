
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
                                TRandom3 &r_gen
                            )

{
        
    std::cout << "\n\n Building LS Data maps\n\n";
        
    ///////////// Pointers to template histos
        
    TH2D* pDAMPE_Template_Data_Iso_LS = &DAMPE_Template_Iso_LS;
    TH2D* DAMPE_Template_Data_AniNS_LS = (TH2D*)DAMPE_Template_Iso_LS.Clone("DAMPE_Template_Data_AniNS_LS");
    TH2D* DAMPE_Template_Data_AniEW_LS = (TH2D*)DAMPE_Template_Iso_LS.Clone("DAMPE_Template_Data_AniEW_LS");
    TH2D* DAMPE_Template_Data_AniFB_LS = (TH2D*)DAMPE_Template_Iso_LS.Clone("DAMPE_Template_Data_AniFB_LS");
    TH2D* DAMPE_Template_MixedData_NS_EW_LS = (TH2D*)DAMPE_Template_Iso_LS.Clone("DAMPE_Template_MixedData_NS_EW_LS");
    TH2D* DAMPE_Template_MixedData_NS_FB_LS = (TH2D*)DAMPE_Template_Iso_LS.Clone("DAMPE_Template_MixedData_NS_FB_LS");
    TH2D* DAMPE_Template_MixedData_EW_FB_LS = (TH2D*)DAMPE_Template_Iso_LS.Clone("DAMPE_Template_MixedData_EW_FB_LS");
    TH2D* DAMPE_Template_FullMixedData_LS = (TH2D*)DAMPE_Template_Iso_LS.Clone("DAMPE_Template_FullMixedData_LS");
        
    DAMPE_Template_Data_AniNS_LS->Reset();
    DAMPE_Template_Data_AniEW_LS->Reset();
    DAMPE_Template_Data_AniFB_LS->Reset();
    DAMPE_Template_MixedData_NS_EW_LS->Reset();
    DAMPE_Template_MixedData_NS_FB_LS->Reset();
    DAMPE_Template_FullMixedData_LS->Reset();
        
    for(Int_t bX = 1; bX <= DAMPE_Template_Data_AniNS_LS->GetNbinsX(); ++bX)
    {
        for(Int_t bY = 1; bY <= DAMPE_Template_Data_AniNS_LS->GetNbinsY(); ++bY)
        {
            DAMPE_Template_Data_AniNS_LS->SetBinContent(bX,bY, (1-NS_anisotropy)*DAMPE_Template_Iso_LS.GetBinContent(bX,bY) + NS_anisotropy*DAMPE_Template_AniNS_LS.GetBinContent(bX,bY) );
            DAMPE_Template_Data_AniEW_LS->SetBinContent(bX,bY, (1-EW_anisotropy)*DAMPE_Template_Iso_LS.GetBinContent(bX,bY) + EW_anisotropy*DAMPE_Template_AniEW_LS.GetBinContent(bX,bY) );
            DAMPE_Template_Data_AniFB_LS->SetBinContent(bX,bY, (1-FB_anisotropy)*DAMPE_Template_Iso_LS.GetBinContent(bX,bY) + FB_anisotropy*DAMPE_Template_AniFB_LS.GetBinContent(bX,bY) );
            DAMPE_Template_MixedData_NS_EW_LS->SetBinContent(bX,bY, (1- NS_anisotropy - EW_anisotropy)*DAMPE_Template_Iso_LS.GetBinContent(bX,bY) + NS_anisotropy*DAMPE_Template_AniNS_LS.GetBinContent(bX,bY) + EW_anisotropy*DAMPE_Template_AniEW_LS.GetBinContent(bX,bY) );
            DAMPE_Template_MixedData_NS_FB_LS->SetBinContent(bX,bY, (1- NS_anisotropy - FB_anisotropy)*DAMPE_Template_Iso_LS.GetBinContent(bX,bY) + NS_anisotropy*DAMPE_Template_AniNS_LS.GetBinContent(bX,bY) + FB_anisotropy*DAMPE_Template_AniFB_LS.GetBinContent(bX,bY) );
            DAMPE_Template_MixedData_EW_FB_LS->SetBinContent(bX,bY, (1- EW_anisotropy - FB_anisotropy)*DAMPE_Template_Iso_LS.GetBinContent(bX,bY) + EW_anisotropy*DAMPE_Template_AniEW_LS.GetBinContent(bX,bY) + FB_anisotropy*DAMPE_Template_AniFB_LS.GetBinContent(bX,bY) );
            DAMPE_Template_FullMixedData_LS->SetBinContent(bX,bY, (1 - NS_anisotropy - EW_anisotropy - FB_anisotropy)*DAMPE_Template_Iso_LS.GetBinContent(bX,bY) + NS_anisotropy*DAMPE_Template_AniNS_LS.GetBinContent(bX,bY) + EW_anisotropy*DAMPE_Template_AniEW_LS.GetBinContent(bX,bY) + FB_anisotropy*DAMPE_Template_AniFB_LS.GetBinContent(bX,bY) );
        }
    }
        
    ///////////// Filling maps...
        
    DAMPE_Data_Iso_LS->FillRandom(pDAMPE_Template_Data_Iso_LS,data_all_sky_LS_events);
    DAMPE_Data_AniNS_LS->FillRandom(DAMPE_Template_Data_AniNS_LS,data_all_sky_LS_events);
    DAMPE_Data_AniEW_LS->FillRandom(DAMPE_Template_Data_AniEW_LS,data_all_sky_LS_events);
    DAMPE_Data_AniFB_LS->FillRandom(DAMPE_Template_Data_AniFB_LS,data_all_sky_LS_events);
    DAMPE_MixedData_NS_EW_LS->FillRandom(DAMPE_Template_MixedData_NS_EW_LS,data_all_sky_LS_events);
    DAMPE_MixedData_NS_FB_LS->FillRandom(DAMPE_Template_MixedData_NS_FB_LS,data_all_sky_LS_events);
    DAMPE_MixedData_EW_FB_LS->FillRandom(DAMPE_Template_MixedData_EW_FB_LS,data_all_sky_LS_events);
    DAMPE_FullMixedData_LS->FillRandom(DAMPE_Template_FullMixedData_LS,data_all_sky_LS_events);

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
                                TRandom3 &r_gen
                            )

{
    std::cout << "\n\n Building HS Data maps\n\n";
    
    ///////////// Pointers to template histos
    
    TH2D* pDAMPE_Template_Data_Iso_HS = &DAMPE_Template_Iso_HS;
    TH2D* DAMPE_Template_Data_AniNS_HS = (TH2D*)DAMPE_Template_Iso_HS.Clone("DAMPE_Template_Data_AniNS_HS");
    TH2D* DAMPE_Template_Data_AniEW_HS = (TH2D*)DAMPE_Template_Iso_HS.Clone("DAMPE_Template_Data_AniEW_HS");
    TH2D* DAMPE_Template_Data_AniFB_HS = (TH2D*)DAMPE_Template_Iso_HS.Clone("DAMPE_Template_Data_AniFB_HS");
    TH2D* DAMPE_Template_MixedData_NS_EW_HS = (TH2D*)DAMPE_Template_Iso_HS.Clone("DAMPE_Template_MixedData_NS_EW_HS");
    TH2D* DAMPE_Template_MixedData_NS_FB_HS = (TH2D*)DAMPE_Template_Iso_HS.Clone("DAMPE_Template_MixedData_NS_FB_HS");
    TH2D* DAMPE_Template_MixedData_EW_FB_HS = (TH2D*)DAMPE_Template_Iso_HS.Clone("DAMPE_Template_MixedData_EW_FB_HS");
    TH2D* DAMPE_Template_FullMixedData_HS = (TH2D*)DAMPE_Template_Iso_HS.Clone("DAMPE_Template_FullMixedData_HS");
    
    DAMPE_Template_Data_AniNS_HS->Reset();
    DAMPE_Template_Data_AniEW_HS->Reset();
    DAMPE_Template_Data_AniFB_HS->Reset();
    DAMPE_Template_MixedData_NS_EW_HS->Reset();
    DAMPE_Template_MixedData_NS_FB_HS->Reset();
    DAMPE_Template_FullMixedData_HS->Reset();
    
    for(Int_t bX = 1; bX <= DAMPE_Template_Data_AniNS_HS->GetNbinsX(); ++bX)
    {
        for(Int_t bY = 1; bY <= DAMPE_Template_Data_AniNS_HS->GetNbinsY(); ++bY)
        {
            DAMPE_Template_Data_AniNS_HS->SetBinContent(bX,bY, (1 - NS_anisotropy)*DAMPE_Template_Iso_HS.GetBinContent(bX,bY) + NS_anisotropy*DAMPE_Template_AniNS_HS.GetBinContent(bX,bY) );
            DAMPE_Template_Data_AniEW_HS->SetBinContent(bX,bY, (1 - EW_anisotropy)*DAMPE_Template_Iso_HS.GetBinContent(bX,bY) + EW_anisotropy*DAMPE_Template_AniEW_HS.GetBinContent(bX,bY) );
            DAMPE_Template_Data_AniFB_HS->SetBinContent(bX,bY, (1 - FB_anisotropy)*DAMPE_Template_Iso_HS.GetBinContent(bX,bY) + FB_anisotropy*DAMPE_Template_AniFB_HS.GetBinContent(bX,bY) );
            DAMPE_Template_MixedData_NS_EW_HS->SetBinContent(bX,bY, (1 - NS_anisotropy - EW_anisotropy)*DAMPE_Template_Iso_HS.GetBinContent(bX,bY) + NS_anisotropy*DAMPE_Template_AniNS_HS.GetBinContent(bX,bY) + EW_anisotropy*DAMPE_Template_AniEW_HS.GetBinContent(bX,bY) );
            DAMPE_Template_MixedData_NS_FB_HS->SetBinContent(bX,bY, (1 - NS_anisotropy - FB_anisotropy)*DAMPE_Template_Iso_HS.GetBinContent(bX,bY) + NS_anisotropy*DAMPE_Template_AniNS_HS.GetBinContent(bX,bY) + FB_anisotropy*DAMPE_Template_AniFB_HS.GetBinContent(bX,bY) );
            DAMPE_Template_MixedData_EW_FB_HS->SetBinContent(bX,bY, (1 - EW_anisotropy - FB_anisotropy)*DAMPE_Template_Iso_HS.GetBinContent(bX,bY) + EW_anisotropy*DAMPE_Template_AniEW_HS.GetBinContent(bX,bY) + FB_anisotropy*DAMPE_Template_AniFB_HS.GetBinContent(bX,bY) );
            DAMPE_Template_FullMixedData_HS->SetBinContent(bX,bY, (1 - NS_anisotropy - EW_anisotropy - FB_anisotropy)*DAMPE_Template_Iso_HS.GetBinContent(bX,bY) + NS_anisotropy*DAMPE_Template_AniNS_HS.GetBinContent(bX,bY) + EW_anisotropy*DAMPE_Template_AniEW_HS.GetBinContent(bX,bY) + FB_anisotropy*DAMPE_Template_AniFB_HS.GetBinContent(bX,bY) );
        }
    }
    
    
    ///////////// Filling maps...
    
    DAMPE_Data_Iso_HS->FillRandom(pDAMPE_Template_Data_Iso_HS,data_all_sky_HS_events);
    DAMPE_Data_AniNS_HS->FillRandom(DAMPE_Template_Data_AniNS_HS,data_all_sky_HS_events);
    DAMPE_Data_AniEW_HS->FillRandom(DAMPE_Template_Data_AniEW_HS,data_all_sky_HS_events);
    DAMPE_Data_AniFB_HS->FillRandom(DAMPE_Template_Data_AniFB_HS,data_all_sky_HS_events);
    DAMPE_MixedData_NS_EW_HS->FillRandom(DAMPE_Template_MixedData_NS_EW_HS,data_all_sky_HS_events);
    DAMPE_MixedData_NS_FB_HS->FillRandom(DAMPE_Template_MixedData_NS_FB_HS,data_all_sky_HS_events);
    DAMPE_MixedData_EW_FB_HS->FillRandom(DAMPE_Template_MixedData_EW_FB_HS,data_all_sky_HS_events);
    DAMPE_FullMixedData_HS->FillRandom(DAMPE_Template_FullMixedData_HS,data_all_sky_HS_events);
    
}

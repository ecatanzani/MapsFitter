
#include "MyHead.h"

void DAMPE_singleTry_fit(
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
                    )

{
    
    std::string data_out_path = output_path_creator(try_idx,2,NS_anisotropy,EW_anisotropy,FB_anisotropy,true);
    std::string pools_out_path = output_path_creator(try_idx,3,NS_anisotropy,EW_anisotropy,FB_anisotropy,true);
    
    tmp_fit.inputAni[0] = NS_anisotropy;
    tmp_fit.inputAni[1] = EW_anisotropy;
    tmp_fit.inputAni[2] = FB_anisotropy;
    
    tmp_fit.seed = tmp_seed;
    tmp_fit.seed_list_line = seed_line;
    
    //////////// TemplateFit variables
    
    Double_t res[4],res_err[4];
    Double_t initialValues[4]={1,0,0,0};
    
    //////////// Open templates
    
    TH2D DAMPE_Template_Iso_LS;
    TH2D DAMPE_Template_AniNS_LS;
    TH2D DAMPE_Template_AniEW_LS;
    TH2D DAMPE_Template_AniFB_LS;
    
    TH2D DAMPE_Template_Iso_HS;
    TH2D DAMPE_Template_AniNS_HS;
    TH2D DAMPE_Template_AniEW_HS;
    TH2D DAMPE_Template_AniFB_HS;
    
    read_templates(
                    DAMPE_Template_Iso_LS,
                    DAMPE_Template_AniNS_LS,
                    DAMPE_Template_AniEW_LS,
                    DAMPE_Template_AniFB_LS,
                    DAMPE_Template_Iso_HS,
                    DAMPE_Template_AniNS_HS,
                    DAMPE_Template_AniEW_HS,
                    DAMPE_Template_AniFB_HS,
                    output_log_file,
                    true
                   );

    
    //////////// Low statistics data
    
    TH2D* DAMPE_Data_Iso_LS = (TH2D*) DAMPE_Template_Iso_LS.Clone("DAMPE_Data_Iso_LS");
    TH2D* DAMPE_Data_AniNS_LS = (TH2D*) DAMPE_Template_Iso_LS.Clone("DAMPE_Data_AniNS_LS");
    TH2D* DAMPE_Data_AniEW_LS = (TH2D*) DAMPE_Template_Iso_LS.Clone("DAMPE_Data_AniEW_LS");
    TH2D* DAMPE_Data_AniFB_LS = (TH2D*) DAMPE_Template_Iso_LS.Clone("DAMPE_Data_AniFB_LS");
    TH2D* DAMPE_MixedData_NS_EW_LS = (TH2D*) DAMPE_Template_Iso_LS.Clone("DAMPE_MixedData_NS_EW_LS");
    TH2D* DAMPE_MixedData_NS_FB_LS = (TH2D*) DAMPE_Template_Iso_LS.Clone("DAMPE_MixedData_NS_FB_LS");
    TH2D* DAMPE_MixedData_EW_FB_LS = (TH2D*) DAMPE_Template_Iso_LS.Clone("DAMPE_MixedData_EW_FB_LS");
    TH2D* DAMPE_FullMixedData_LS = (TH2D*) DAMPE_Template_Iso_LS.Clone("DAMPE_FullMixedData_LS");
    
    DAMPE_Data_Iso_LS->SetTitle("Isotropic LS DAMPE (data) Map; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries");
    DAMPE_Data_AniNS_LS->SetTitle("Anisotropic LS DAMPE (data) Map (NS); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries");
    DAMPE_Data_AniEW_LS->SetTitle("Anisotropic LS DAMPE (data) Map (EW); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries");
    DAMPE_Data_AniFB_LS->SetTitle("Anisotropic LS DAMPE (data) Map (FB); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries");
    DAMPE_MixedData_NS_EW_LS->SetTitle("Anisotropic LS DAMPE (NS-EW data) Map; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries");
    DAMPE_MixedData_NS_FB_LS->SetTitle("Anisotropic LS DAMPE (NS-FB data) Map; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries");
    DAMPE_MixedData_EW_FB_LS->SetTitle("Anisotropic LS DAMPE (EW-FB data) Map; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries");
    DAMPE_FullMixedData_LS->SetTitle("Anisotropic LS DAMPE (NS-EW-FB data) Map; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries");
    
    DAMPE_Data_Iso_LS->Reset();
    DAMPE_Data_AniNS_LS->Reset();
    DAMPE_Data_AniEW_LS->Reset();
    DAMPE_Data_AniFB_LS->Reset();
    DAMPE_MixedData_NS_EW_LS->Reset();
    DAMPE_MixedData_NS_FB_LS->Reset();
    DAMPE_MixedData_EW_FB_LS->Reset();
    DAMPE_FullMixedData_LS->Reset();
    
    //////////// High statistics data
    
    TH2D* DAMPE_Data_Iso_HS = (TH2D*) DAMPE_Template_Iso_HS.Clone("DAMPE_Data_Iso_HS");
    TH2D* DAMPE_Data_AniNS_HS = (TH2D*) DAMPE_Template_Iso_HS.Clone("DAMPE_Data_AniNS_HS");
    TH2D* DAMPE_Data_AniEW_HS = (TH2D*) DAMPE_Template_Iso_HS.Clone("DAMPE_Data_AniEW_HS");
    TH2D* DAMPE_Data_AniFB_HS = (TH2D*) DAMPE_Template_Iso_HS.Clone("DAMPE_Data_AniFB_HS");
    TH2D* DAMPE_MixedData_NS_EW_HS = (TH2D*) DAMPE_Template_Iso_HS.Clone("DAMPE_MixedData_NS_EW_HS");
    TH2D* DAMPE_MixedData_NS_FB_HS = (TH2D*) DAMPE_Template_Iso_HS.Clone("DAMPE_MixedData_NS_FB_HS");
    TH2D* DAMPE_MixedData_EW_FB_HS = (TH2D*) DAMPE_Template_Iso_HS.Clone("DAMPE_MixedData_EW_FB_HS");
    TH2D* DAMPE_FullMixedData_HS = (TH2D*) DAMPE_Template_Iso_HS.Clone("DAMPE_FullMixedData_HS");
    
    DAMPE_Data_Iso_HS->SetTitle("Isotropic HS DAMPE (data) Map; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries");
    DAMPE_Data_AniNS_HS->SetTitle("Anisotropic HS DAMPE (data) Map (NS); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries");
    DAMPE_Data_AniEW_HS->SetTitle("Anisotropic HS DAMPE (data) Map (EW); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries");
    DAMPE_Data_AniFB_HS->SetTitle("Anisotropic HS DAMPE (data) Map (FB); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries");
    DAMPE_MixedData_NS_EW_HS->SetTitle("Anisotropic HS DAMPE (NS-EW data) Map; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries");
    DAMPE_MixedData_NS_FB_HS->SetTitle("Anisotropic HS DAMPE (NS-FB data) Map; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries");
    DAMPE_MixedData_EW_FB_HS->SetTitle("Anisotropic HS DAMPE (EW-FB data) Map; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries");
    DAMPE_FullMixedData_HS->SetTitle("Anisotropic HS DAMPE (NS-EW-FB data) Map; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries");
    
    DAMPE_Data_Iso_HS->Reset();
    DAMPE_Data_AniNS_HS->Reset();
    DAMPE_Data_AniEW_HS->Reset();
    DAMPE_Data_AniFB_HS->Reset();
    DAMPE_MixedData_NS_EW_HS->Reset();
    DAMPE_MixedData_NS_FB_HS->Reset();
    DAMPE_MixedData_EW_FB_HS->Reset();
    DAMPE_FullMixedData_HS->Reset();
    
    //////////// Pool's histos
    
    //////////// Low statistics pools
    
    TH1D DAMPE_hPull_Iso_LS("DAMPE_hPull_Iso_LS","Pull Isotropic DAMPE Sky Map LS",50,-5,5);
    TH1D DAMPE_hPull_AniNS_LS("DAMPE_hPull_AniNS_LS","Pull Anisotropic DAMPE Sky Map (NS) LS",50,-5,5);
    TH1D DAMPE_hPull_AniEW_LS("DAMPE_hPull_AniEW_LS","Pull Anisotropic DAMPE Sky Map (EW) LS",50,-5,5);
    TH1D DAMPE_hPull_AniFB_LS("DAMPE_hPull_AniFB_LS","Pull Anisotropic DAMPE Sky Map (FB) LS",50,-5,5);
    TH1D DAMPE_hPull_Mixed_NS_EW_LS("DAMPE_hPull_Mixed_NS_EW_LS","Pull Anisotropic DAMPE Sky Map (NS-EW) LS",50,-5,5);
    TH1D DAMPE_hPull_Mixed_NS_FB_LS("DAMPE_hPull_Mixed_NS_FB_LS","Pull Anisotropic DAMPE Sky Map (NS-FB) LS",50,-5,5);
    TH1D DAMPE_hPull_Mixed_EW_FB_LS("DAMPE_hPull_Mixed_EW_FB_LS","Pull Anisotropic DAMPE Sky Map (EW-FB) LS",50,-5,5);
    TH1D DAMPE_hPull_FullMixed_LS("DAMPE_hPull_FullMixed_LS","Pull Full Mixed Anisotropic DAMPE Sky Map LS",50,-5,5);
    
    //////////// High statistics pools
    
    TH1D DAMPE_hPull_Iso_HS("DAMPE_hPull_Iso_HS","Pull Isotropic DAMPE Sky Map HS",100,-100,100);
    TH1D DAMPE_hPull_AniNS_HS("DAMPE_hPull_AniNS_HS","Pull Anisotropic DAMPE Sky Map (NS) HS",100,-100,100);
    TH1D DAMPE_hPull_AniEW_HS("DAMPE_hPull_AniEW_HS","Pull Anisotropic DAMPE Sky Map (EW) HS",100,-100,100);
    TH1D DAMPE_hPull_AniFB_HS("DAMPE_hPull_AniFB_HS","Pull Anisotropic DAMPE Sky Map (FB) HS",100,-100,100);
    TH1D DAMPE_hPull_Mixed_NS_EW_HS("DAMPE_hPull_Mixed_NS_EW_HS","Pull Anisotropic DAMPE Sky Map (NS-EW) HS",100,-100,100);
    TH1D DAMPE_hPull_Mixed_NS_FB_HS("DAMPE_hPull_Mixed_NS_FB_HS","Pull Anisotropic DAMPE Sky Map (NS-FB) HS",100,-100,100);
    TH1D DAMPE_hPull_Mixed_EW_FB_HS("DAMPE_hPull_Mixed_EW_FB_HS","Pull Anisotropic DAMPE Sky Map (EW-FB) HS",100,-100,100);
    TH1D DAMPE_hPull_FullMixed_HS("DAMPE_hPull_FullMixed_HS","Pull Full Mixed Anisotropic DAMPE Sky Map HS",100,-100,100);
    
    
    //////////// 1D Histos
    
    TH1D DAMPE_Templates_LS[4];
    TH1D DAMPE_Templates_HS[4];
    
    TH1D DAMPE_DataHisto_I_LS;
    TH1D DAMPE_DataHisto_NS_LS;
    TH1D DAMPE_DataHisto_EW_LS;
    TH1D DAMPE_DataHisto_FB_LS;
    TH1D DAMPE_MixedDataHisto_NS_EW_LS;
    TH1D DAMPE_MixedDataHisto_NS_FB_LS;
    TH1D DAMPE_MixedDataHisto_EW_FB_LS;
    TH1D DAMPE_FullMixedDataHisto_LS;
    
    TH1D DAMPE_DataHisto_I_HS;
    TH1D DAMPE_DataHisto_NS_HS;
    TH1D DAMPE_DataHisto_EW_HS;
    TH1D DAMPE_DataHisto_FB_HS;
    TH1D DAMPE_MixedDataHisto_NS_EW_HS;
    TH1D DAMPE_MixedDataHisto_NS_FB_HS;
    TH1D DAMPE_MixedDataHisto_EW_FB_HS;
    TH1D DAMPE_FullMixedDataHisto_HS;
    
    TH1* DAMPE_Templates_1D_LS[4] = {nullptr,nullptr,nullptr,nullptr};
    TH1* DAMPE_Templates_1D_HS[4] = {nullptr,nullptr,nullptr,nullptr};
    
    TH1D* DAMPE_DataHisto_1D_I_LS = nullptr;
    TH1D* DAMPE_DataHisto_1D_AniNS_LS = nullptr;
    TH1D* DAMPE_DataHisto_1D_AniEW_LS = nullptr;
    TH1D* DAMPE_DataHisto_1D_AniFB_LS = nullptr;
    TH1D* DAMPE_MixedHisto_1D_NS_EW_LS = nullptr;
    TH1D* DAMPE_MixedHisto_1D_NS_FB_LS = nullptr;
    TH1D* DAMPE_MixedHisto_1D_EW_FB_LS = nullptr;
    TH1D* DAMPE_FullMixedHisto_1D_LS = nullptr;
    
    TH1D* DAMPE_DataHisto_1D_I_HS = nullptr;
    TH1D* DAMPE_DataHisto_1D_AniNS_HS = nullptr;
    TH1D* DAMPE_DataHisto_1D_AniEW_HS = nullptr;
    TH1D* DAMPE_DataHisto_1D_AniFB_HS = nullptr;
    TH1D* DAMPE_MixedHisto_1D_NS_EW_HS = nullptr;
    TH1D* DAMPE_MixedHisto_1D_NS_FB_HS = nullptr;
    TH1D* DAMPE_MixedHisto_1D_EW_FB_HS = nullptr;
    TH1D* DAMPE_FullMixedHisto_1D_HS = nullptr;
    
    tmp_fit.theta_binHistos_LS = DAMPE_Data_Iso_LS->GetNbinsY();
    tmp_fit.phi_binHistos_LS = DAMPE_Data_Iso_LS->GetNbinsX();
    tmp_fit.theta_binHistos_HS = DAMPE_Data_Iso_HS->GetNbinsY();
    tmp_fit.phi_binHistos_HS = DAMPE_Data_Iso_HS->GetNbinsX();
    
    tmp_fit.events_LS = data_all_sky_LS_events;
    tmp_fit.events_HS = data_all_sky_HS_events;
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    generate_DAMPE_LS_data(
                            NS_anisotropy,
                            EW_anisotropy,
                            FB_anisotropy,
                            DAMPE_Data_Iso_LS,
                            DAMPE_Data_AniNS_LS,
                            DAMPE_Data_AniEW_LS,
                            DAMPE_Data_AniFB_LS,
                            DAMPE_MixedData_NS_EW_LS,
                            DAMPE_MixedData_NS_FB_LS,
                            DAMPE_MixedData_EW_FB_LS,
                            DAMPE_FullMixedData_LS,
                            DAMPE_Template_Iso_LS,
                            DAMPE_Template_AniNS_LS,
                            DAMPE_Template_AniEW_LS,
                            DAMPE_Template_AniFB_LS,
                            output_log_file,
                            tmp_seed
                           );
    
    generate_DAMPE_HS_data(
                            NS_anisotropy,
                            EW_anisotropy,
                            FB_anisotropy,
                            DAMPE_Data_Iso_HS,
                            DAMPE_Data_AniNS_HS,
                            DAMPE_Data_AniEW_HS,
                            DAMPE_Data_AniFB_HS,
                            DAMPE_MixedData_NS_EW_HS,
                            DAMPE_MixedData_NS_FB_HS,
                            DAMPE_MixedData_EW_FB_HS,
                            DAMPE_FullMixedData_HS,
                            DAMPE_Template_Iso_HS,
                            DAMPE_Template_AniNS_HS,
                            DAMPE_Template_AniEW_HS,
                            DAMPE_Template_AniFB_HS,
                            output_log_file,
                            tmp_seed
                           );
    
    
    ///////////////////////// Converting Templates Maps from TH2 to TH1
    
    TH2toTH1_obj(DAMPE_Templates_LS[0],DAMPE_Template_Iso_LS);
    TH2toTH1_obj(DAMPE_Templates_LS[1],DAMPE_Template_AniNS_LS);
    TH2toTH1_obj(DAMPE_Templates_LS[2],DAMPE_Template_AniEW_LS);
    TH2toTH1_obj(DAMPE_Templates_LS[3],DAMPE_Template_AniFB_LS);
    
    TH2toTH1_obj(DAMPE_Templates_HS[0],DAMPE_Template_Iso_HS);
    TH2toTH1_obj(DAMPE_Templates_HS[1],DAMPE_Template_AniNS_HS);
    TH2toTH1_obj(DAMPE_Templates_HS[2],DAMPE_Template_AniEW_HS);
    TH2toTH1_obj(DAMPE_Templates_HS[3],DAMPE_Template_AniFB_HS);
    
    ///////////////////////// Converting Data Maps from TH2 to TH1
    
    TH2toTH1_ptr(DAMPE_DataHisto_I_LS,DAMPE_Data_Iso_LS);
    TH2toTH1_ptr(DAMPE_DataHisto_NS_LS,DAMPE_Data_AniNS_LS);
    TH2toTH1_ptr(DAMPE_DataHisto_EW_LS,DAMPE_Data_AniEW_LS);
    TH2toTH1_ptr(DAMPE_DataHisto_FB_LS,DAMPE_Data_AniFB_LS);
    
    TH2toTH1_ptr(DAMPE_MixedDataHisto_NS_EW_LS,DAMPE_MixedData_NS_EW_LS);
    TH2toTH1_ptr(DAMPE_MixedDataHisto_NS_FB_LS,DAMPE_MixedData_NS_FB_LS);
    TH2toTH1_ptr(DAMPE_MixedDataHisto_EW_FB_LS,DAMPE_MixedData_EW_FB_LS);
    TH2toTH1_ptr(DAMPE_FullMixedDataHisto_LS,DAMPE_FullMixedData_LS);
    
    TH2toTH1_ptr(DAMPE_DataHisto_I_HS,DAMPE_Data_Iso_HS);
    TH2toTH1_ptr(DAMPE_DataHisto_NS_HS,DAMPE_Data_AniNS_HS);
    TH2toTH1_ptr(DAMPE_DataHisto_EW_HS,DAMPE_Data_AniEW_HS);
    TH2toTH1_ptr(DAMPE_DataHisto_FB_HS,DAMPE_Data_AniFB_HS);
    
    TH2toTH1_ptr(DAMPE_MixedDataHisto_NS_EW_HS,DAMPE_MixedData_NS_EW_HS);
    TH2toTH1_ptr(DAMPE_MixedDataHisto_NS_FB_HS,DAMPE_MixedData_NS_FB_HS);
    TH2toTH1_ptr(DAMPE_MixedDataHisto_EW_FB_HS,DAMPE_MixedData_EW_FB_HS);
    TH2toTH1_ptr(DAMPE_FullMixedDataHisto_HS,DAMPE_FullMixedData_HS);
    
    //////////// Saving bin contents into the Tree
    
    for(Int_t idx = 0; idx < DAMPE_DataHisto_I_LS.GetNbinsX(); ++idx)
    {
        
        tmp_fit.entries_Iso_Map_LS[idx_ani][idx] = DAMPE_DataHisto_I_LS.GetBinContent(idx+1);
        tmp_fit.entries_NS_Map_LS[idx_ani][idx] = DAMPE_DataHisto_NS_LS.GetBinContent(idx+1);
        tmp_fit.entries_EW_Map_LS[idx_ani][idx] = DAMPE_DataHisto_EW_LS.GetBinContent(idx+1);
        tmp_fit.entries_FB_Map_LS[idx_ani][idx] = DAMPE_DataHisto_FB_LS.GetBinContent(idx+1);
        tmp_fit.entries_NS_EW_Map_LS[idx_ani][idx] = DAMPE_MixedDataHisto_NS_EW_LS.GetBinContent(idx+1);
        tmp_fit.entries_NS_FB_Map_LS[idx_ani][idx] = DAMPE_MixedDataHisto_NS_FB_LS.GetBinContent(idx+1);
        tmp_fit.entries_EW_FB_Map_LS[idx_ani][idx] = DAMPE_MixedDataHisto_EW_FB_LS.GetBinContent(idx+1);
        tmp_fit.entries_Full_Map_LS[idx_ani][idx] = DAMPE_FullMixedDataHisto_LS.GetBinContent(idx+1);
        
    }
    
    for(Int_t idx = 0; idx < DAMPE_DataHisto_I_HS.GetNbinsX(); ++idx)
    {
        
        tmp_fit.entries_Iso_Map_HS[idx_ani][idx] = DAMPE_DataHisto_I_HS.GetBinContent(idx+1);
        tmp_fit.entries_NS_Map_HS[idx_ani][idx] = DAMPE_DataHisto_NS_HS.GetBinContent(idx+1);
        tmp_fit.entries_EW_Map_HS[idx_ani][idx] = DAMPE_DataHisto_EW_HS.GetBinContent(idx+1);
        tmp_fit.entries_FB_Map_HS[idx_ani][idx] = DAMPE_DataHisto_FB_HS.GetBinContent(idx+1);
        tmp_fit.entries_NS_EW_Map_HS[idx_ani][idx] = DAMPE_MixedDataHisto_NS_EW_HS.GetBinContent(idx+1);
        tmp_fit.entries_NS_FB_Map_HS[idx_ani][idx] = DAMPE_MixedDataHisto_NS_FB_HS.GetBinContent(idx+1);
        tmp_fit.entries_EW_FB_Map_HS[idx_ani][idx] = DAMPE_MixedDataHisto_EW_FB_HS.GetBinContent(idx+1);
        tmp_fit.entries_Full_Map_HS[idx_ani][idx] = DAMPE_FullMixedDataHisto_HS.GetBinContent(idx+1);
        
    }
    
    ///////////////////////// Linking variables to the pointers, prepearing for the fit
    
    for(Int_t idx_t = 0; idx_t < 4; idx_t++) {
        DAMPE_Templates_1D_LS[idx_t] = &DAMPE_Templates_LS[idx_t];
        DAMPE_Templates_1D_HS[idx_t] = &DAMPE_Templates_HS[idx_t];
    }
    
    DAMPE_DataHisto_1D_I_LS = &DAMPE_DataHisto_I_LS;
    DAMPE_DataHisto_1D_AniNS_LS = &DAMPE_DataHisto_NS_LS;
    DAMPE_DataHisto_1D_AniEW_LS = &DAMPE_DataHisto_EW_LS;
    DAMPE_DataHisto_1D_AniFB_LS = &DAMPE_DataHisto_FB_LS;
    
    DAMPE_MixedHisto_1D_NS_EW_LS = &DAMPE_MixedDataHisto_NS_EW_LS;
    DAMPE_MixedHisto_1D_NS_FB_LS = &DAMPE_MixedDataHisto_NS_FB_LS;
    DAMPE_MixedHisto_1D_EW_FB_LS = &DAMPE_MixedDataHisto_EW_FB_LS;
    DAMPE_FullMixedHisto_1D_LS = &DAMPE_FullMixedDataHisto_LS;
    
    DAMPE_DataHisto_1D_I_HS = &DAMPE_DataHisto_I_HS;
    DAMPE_DataHisto_1D_AniNS_HS = &DAMPE_DataHisto_NS_HS;
    DAMPE_DataHisto_1D_AniEW_HS = &DAMPE_DataHisto_EW_HS;
    DAMPE_DataHisto_1D_AniFB_HS = &DAMPE_DataHisto_FB_HS;
    
    DAMPE_MixedHisto_1D_NS_EW_HS = &DAMPE_MixedDataHisto_NS_EW_HS;
    DAMPE_MixedHisto_1D_NS_FB_HS = &DAMPE_MixedDataHisto_NS_FB_HS;
    DAMPE_MixedHisto_1D_EW_FB_HS = &DAMPE_MixedDataHisto_EW_FB_HS;
    DAMPE_FullMixedHisto_1D_HS = &DAMPE_FullMixedDataHisto_HS;
    
    /////////////////////////////////////////////////////////////////////////////// Fitting !!!
    
    std::cout << "\n\n //////////////////////////////////////// Low Statistics DAMPE fits //////////////////////////////////////// \n\n";
    output_log_file << "\n\n //////////////////////////////////////// Low Statistics DAMPE fits //////////////////////////////////////// \n\n";
    
    std::cout << "\n\n -------------- > Isotropic Monopole \n\n";
    output_log_file << "\n\n -------------- > Isotropic Monopole \n\n";
    
    TemplateFitBH(DAMPE_DataHisto_1D_I_LS,4,DAMPE_Templates_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,0,false);
    getPull(DAMPE_DataHisto_1D_I_LS,DAMPE_Templates_1D_LS,res,DAMPE_hPull_Iso_LS);
    
    std::cout << "\n\n -------------- > Anisotropic Monopole (NS) \n\n";
    output_log_file << "\n\n -------------- > Anisotropic Monopole (NS) \n\n";
    
    TemplateFitBH(DAMPE_DataHisto_1D_AniNS_LS,4,DAMPE_Templates_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,1,false);
    getPull(DAMPE_DataHisto_1D_AniNS_LS,DAMPE_Templates_1D_LS,res,DAMPE_hPull_AniNS_LS);
    
    std::cout << "\n\n -------------- > Anisotropic Monopole (EW) \n\n";
    output_log_file << "\n\n -------------- > Anisotropic Monopole (EW) \n\n";
    
    TemplateFitBH(DAMPE_DataHisto_1D_AniEW_LS,4,DAMPE_Templates_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,2,false);
    getPull(DAMPE_DataHisto_1D_AniEW_LS,DAMPE_Templates_1D_LS,res,DAMPE_hPull_AniEW_LS);
    
    std::cout << "\n\n -------------- > Anisotropic Monopole (FB) \n\n";
    output_log_file << "\n\n -------------- > Anisotropic Monopole (FB) \n\n";
    
    TemplateFitBH(DAMPE_DataHisto_1D_AniFB_LS,4,DAMPE_Templates_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,3,false);
    getPull(DAMPE_DataHisto_1D_AniFB_LS,DAMPE_Templates_1D_LS,res,DAMPE_hPull_AniFB_LS);
    
    std::cout << "\n\n -------------- > NS + EW linear combination \n\n";
    output_log_file << "\n\n -------------- > NS + EW linear combination \n\n";
    
    TemplateFitBH(DAMPE_MixedHisto_1D_NS_EW_LS,4,DAMPE_Templates_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,4,false);
    getPull(DAMPE_MixedHisto_1D_NS_EW_LS,DAMPE_Templates_1D_LS,res,DAMPE_hPull_Mixed_NS_EW_LS);
    
    std::cout << "\n\n -------------- > NS + FB linear combination \n\n";
    output_log_file << "\n\n -------------- > NS + FB linear combination \n\n";
    
    TemplateFitBH(DAMPE_MixedHisto_1D_NS_FB_LS,4,DAMPE_Templates_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,5,false);
    getPull(DAMPE_MixedHisto_1D_NS_FB_LS,DAMPE_Templates_1D_LS,res,DAMPE_hPull_Mixed_NS_FB_LS);
    
    std::cout << "\n\n -------------- > EW + FB linear combination \n\n";
    output_log_file << "\n\n -------------- > EW + FB linear combination \n\n";
    
    TemplateFitBH(DAMPE_MixedHisto_1D_EW_FB_LS,4,DAMPE_Templates_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,6,false);
    getPull(DAMPE_MixedHisto_1D_EW_FB_LS,DAMPE_Templates_1D_LS,res,DAMPE_hPull_Mixed_EW_FB_LS);
    
    std::cout << "\n\n -------------- > FullMixed linear combination \n\n";
    output_log_file << "\n\n -------------- > FullMixed linear combination \n\n";
    
    TemplateFitBH(DAMPE_FullMixedHisto_1D_LS,4,DAMPE_Templates_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,7,false);
    getPull(DAMPE_FullMixedHisto_1D_LS,DAMPE_Templates_1D_LS,res,DAMPE_hPull_FullMixed_LS);
    
    
    std::cout << "\n\n //////////////////////////////////////// High Statistics DAMPE fits //////////////////////////////////////// \n\n";
    output_log_file << "\n\n //////////////////////////////////////// High Statistics DAMPE fits //////////////////////////////////////// \n\n";
    
    std::cout << "\n\n -------------- > Isotropic Monopole \n\n";
    output_log_file << "\n\n -------------- > Isotropic Monopole \n\n";
    
    TemplateFitBH(DAMPE_DataHisto_1D_I_HS,4,DAMPE_Templates_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,0,true);
    getPull(DAMPE_DataHisto_1D_I_HS,DAMPE_Templates_1D_HS,res,DAMPE_hPull_Iso_HS);
    
    std::cout << "\n\n -------------- > Anisotropic Monopole (NS) \n\n";
    output_log_file << "\n\n -------------- > Anisotropic Monopole (NS) \n\n";
    
    TemplateFitBH(DAMPE_DataHisto_1D_AniNS_HS,4,DAMPE_Templates_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,1,true);
    getPull(DAMPE_DataHisto_1D_AniNS_HS,DAMPE_Templates_1D_HS,res,DAMPE_hPull_AniNS_HS);
    
    std::cout << "\n\n -------------- > Anisotropic Monopole (EW) \n\n";
    output_log_file << "\n\n -------------- > Anisotropic Monopole (EW) \n\n";
    
    TemplateFitBH(DAMPE_DataHisto_1D_AniEW_HS,4,DAMPE_Templates_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,2,true);
    getPull(DAMPE_DataHisto_1D_AniEW_HS,DAMPE_Templates_1D_HS,res,DAMPE_hPull_AniEW_HS);
    
    std::cout << "\n\n -------------- > Anisotropic Monopole (FB) \n\n";
    output_log_file << "\n\n -------------- > Anisotropic Monopole (FB) \n\n";
    
    TemplateFitBH(DAMPE_DataHisto_1D_AniFB_HS,4,DAMPE_Templates_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,3,true);
    getPull(DAMPE_DataHisto_1D_AniFB_HS,DAMPE_Templates_1D_HS,res,DAMPE_hPull_AniFB_HS);
    
    std::cout << "\n\n -------------- > NS + EW linear combination \n\n";
    output_log_file << "\n\n -------------- > NS + EW linear combination \n\n";
    
    TemplateFitBH(DAMPE_MixedHisto_1D_NS_EW_HS,4,DAMPE_Templates_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,4,true);
    getPull(DAMPE_MixedHisto_1D_NS_EW_HS,DAMPE_Templates_1D_HS,res,DAMPE_hPull_Mixed_NS_EW_HS);
    
    std::cout << "\n\n -------------- > NS + FB linear combination \n\n";
    output_log_file << "\n\n -------------- > NS + FB linear combination \n\n";
    
    TemplateFitBH(DAMPE_MixedHisto_1D_NS_FB_HS,4,DAMPE_Templates_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,5,true);
    getPull(DAMPE_MixedHisto_1D_NS_FB_HS,DAMPE_Templates_1D_HS,res,DAMPE_hPull_Mixed_NS_FB_HS);
    
    std::cout << "\n\n -------------- > EW + FB linear combination \n\n";
    output_log_file << "\n\n -------------- > EW + FB linear combination \n\n";
    
    TemplateFitBH(DAMPE_MixedHisto_1D_EW_FB_HS,4,DAMPE_Templates_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,6,true);
    getPull(DAMPE_MixedHisto_1D_EW_FB_HS,DAMPE_Templates_1D_HS,res,DAMPE_hPull_Mixed_EW_FB_HS);
    
    std::cout << "\n\n -------------- > FullMixed linear combination \n\n";
    output_log_file << "\n\n -------------- > FullMixed linear combination \n\n";
    
    TemplateFitBH(DAMPE_FullMixedHisto_1D_HS,4,DAMPE_Templates_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,7,true);
    getPull(DAMPE_FullMixedHisto_1D_HS,DAMPE_Templates_1D_HS,res,DAMPE_hPull_FullMixed_HS);
    
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    std::cout<<"\n\nSimulation Completed !\n\n";
    output_log_file << "\n\nSimulation Completed !\n\n";
    
    
    if(write_tmp_histos && itry==0)
    {
        
        //////////////////////////////////// Creating Data out file
        
        TFile DAMPE_data_file(data_out_path.c_str(),"RECREATE");
        if(DAMPE_data_file.IsZombie()) {
            std::cout << "\n\nError writing ROOT DAMPE Data TFile. Prorgram finished \n\n";
            output_log_file << "\n\nError writing ROOT DAMPE Data TFile. Prorgram finished \n\n";
            exit(100);
        }
        
        
        //////////////////////////////////// Writing Data
        
        DAMPE_Data_Iso_LS->Write();
        DAMPE_Data_AniNS_LS->Write();
        DAMPE_Data_AniEW_LS->Write();
        DAMPE_Data_AniFB_LS->Write();
        
        DAMPE_MixedData_NS_EW_LS->Write();
        DAMPE_MixedData_NS_FB_LS->Write();
        DAMPE_MixedData_EW_FB_LS->Write();
        DAMPE_FullMixedData_LS->Write();
        
        DAMPE_Data_Iso_HS->Write();
        DAMPE_Data_AniNS_HS->Write();
        DAMPE_Data_AniEW_HS->Write();
        DAMPE_Data_AniFB_HS->Write();
        
        DAMPE_MixedData_NS_EW_HS->Write();
        DAMPE_MixedData_NS_FB_HS->Write();
        DAMPE_MixedData_EW_FB_HS->Write();
        DAMPE_FullMixedData_HS->Write();
        
        DAMPE_data_file.Write();
        DAMPE_data_file.Close();
        
        
        //////////////////////////////////// Creating Pools out file
        
        TFile DAMPE_pool_file(pools_out_path.c_str(),"RECREATE");
        if(DAMPE_pool_file.IsZombie()) {
            std::cout << "\n\nError writing ROOT DAMPE Pools TFile. Prorgram finished \n\n";
            output_log_file << "\n\nError writing ROOT DAMPE Pools TFile. Prorgram finished \n\n";
            exit(100);
        }
        
        //////////////////////////////////// Gaussian Fit LS hPools
        
        DAMPE_hPull_Iso_LS.Fit("gaus","L,E,M");
        DAMPE_hPull_AniNS_LS.Fit("gaus","L,E,M");
        DAMPE_hPull_AniEW_LS.Fit("gaus","L,E,M");
        DAMPE_hPull_AniFB_LS.Fit("gaus","L,E,M");
        DAMPE_hPull_Mixed_NS_EW_LS.Fit("gaus","L,E,M");
        DAMPE_hPull_Mixed_NS_FB_LS.Fit("gaus","L,E,M");
        DAMPE_hPull_Mixed_EW_FB_LS.Fit("gaus","L,E,M");
        DAMPE_hPull_FullMixed_LS.Fit("gaus","L,E,M");
        
        //////////////////////////////////// Gaussian Fit HS hPools
        
        DAMPE_hPull_Iso_HS.Fit("gaus","L,E,M");
        DAMPE_hPull_AniNS_HS.Fit("gaus","L,E,M");
        DAMPE_hPull_AniEW_HS.Fit("gaus","L,E,M");
        DAMPE_hPull_AniFB_HS.Fit("gaus","L,E,M");
        DAMPE_hPull_Mixed_NS_EW_HS.Fit("gaus","L,E,M");
        DAMPE_hPull_Mixed_NS_FB_HS.Fit("gaus","L,E,M");
        DAMPE_hPull_Mixed_EW_FB_HS.Fit("gaus","L,E,M");
        DAMPE_hPull_FullMixed_HS.Fit("gaus","L,E,M");
        
        //////////////////////////////////// Writing Pools
        
        DAMPE_hPull_Iso_LS.Write();
        DAMPE_hPull_AniNS_LS.Write();
        DAMPE_hPull_AniEW_LS.Write();
        DAMPE_hPull_AniFB_LS.Write();
        DAMPE_hPull_Mixed_NS_EW_LS.Write();
        DAMPE_hPull_Mixed_NS_FB_LS.Write();
        DAMPE_hPull_Mixed_EW_FB_LS.Write();
        DAMPE_hPull_FullMixed_LS.Write();
        
        DAMPE_hPull_Iso_HS.Write();
        DAMPE_hPull_AniNS_HS.Write();
        DAMPE_hPull_AniEW_HS.Write();
        DAMPE_hPull_AniFB_HS.Write();
        DAMPE_hPull_Mixed_NS_EW_HS.Write();
        DAMPE_hPull_Mixed_NS_FB_HS.Write();
        DAMPE_hPull_Mixed_EW_FB_HS.Write();
        DAMPE_hPull_FullMixed_HS.Write();
        
        DAMPE_pool_file.Write();
        DAMPE_pool_file.Close();
        
    }
    
    fiTree.Fill();
    
}

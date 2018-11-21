
#include "MyHead.h"

void DAMPE_relative_singleTry_fit(
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
                                  )

{
    
    std::string data_out_path = output_path_creator(
                                                        output_log,
                                                        output_root,
                                                        2,
                                                        time_stamp,
                                                        NS_anisotropy,
                                                        EW_anisotropy,
                                                        FB_anisotropy,
                                                        true
                                                    );
    
    std::string pools_out_path = output_path_creator(
                                                         output_log,
                                                         output_root,
                                                         3,
                                                         time_stamp,
                                                         NS_anisotropy,
                                                         EW_anisotropy,
                                                         FB_anisotropy,
                                                         true
                                                     );
    
    //////////// TemplateFit variables
    
    Double_t res[3],res_err[3];
    Double_t initialValues[3]={0,0,0};
    
    //////////// Open templates
    
    TH2D DAMPE_Template_Iso_LS;
    TH2D DAMPE_Template_AniNS_LS;
    TH2D DAMPE_Template_AniEW_LS;
    TH2D DAMPE_Template_AniFB_LS;
    
    TH2D DAMPE_Template_Iso_HS;
    TH2D DAMPE_Template_AniNS_HS;
    TH2D DAMPE_Template_AniEW_HS;
    TH2D DAMPE_Template_AniFB_HS;
    
    TH2D Template_Iso_LS;
    TH2D Template_AniNS_LS;
    TH2D Template_AniEW_LS;
    TH2D Template_AniFB_LS;
    
    TH2D Template_Iso_HS;
    TH2D Template_AniNS_HS;
    TH2D Template_AniEW_HS;
    TH2D Template_AniFB_HS;
    
    
    read_templates(
                       Template_Iso_LS,
                       Template_AniNS_LS,
                       Template_AniEW_LS,
                       Template_AniFB_LS,
                       Template_Iso_HS,
                       Template_AniNS_HS,
                       Template_AniEW_HS,
                       Template_AniFB_HS,
                       output_log_file,
                       false,
                       template_out_path,
                       DAMPE_template_out_path
                   );
    
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
                       true,
                       template_out_path,
                       DAMPE_template_out_path
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
    
    TH2D relative_DAMPE_Data_AniNS_LS;
    TH2D relative_DAMPE_Data_AniEW_LS;
    TH2D relative_DAMPE_Data_AniFB_LS;
    TH2D relative_DAMPE_MixedData_NS_EW_LS;
    TH2D relative_DAMPE_MixedData_NS_FB_LS;
    TH2D relative_DAMPE_MixedData_EW_FB_LS;
    TH2D relative_DAMPE_FullMixedData_LS;
    
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
    
    TH2D relative_DAMPE_Data_AniNS_HS;
    TH2D relative_DAMPE_Data_AniEW_HS;
    TH2D relative_DAMPE_Data_AniFB_HS;
    TH2D relative_DAMPE_MixedData_NS_EW_HS;
    TH2D relative_DAMPE_MixedData_NS_FB_HS;
    TH2D relative_DAMPE_MixedData_EW_FB_HS;
    TH2D relative_DAMPE_FullMixedData_HS;
    
    //////////// Pool's histos
    
    //////////// Low statistics pools
    
    TH1D relative_DAMPE_hPull_AniNS_LS("relative_DAMPE_hPull_AniNS_LS","Pull Anisotropic relative DAMPE Sky Map (NS) LS",100,-5,5);
    TH1D relative_DAMPE_hPull_AniEW_LS("relative_DAMPE_hPull_AniEW_LS","Pull Anisotropic relative DAMPE Sky Map (EW) LS",100,-5,5);
    TH1D relative_DAMPE_hPull_AniFB_LS("relative_DAMPE_hPull_AniFB_LS","Pull Anisotropic relative DAMPE Sky Map (FB) LS",100,-5,5);
    TH1D relative_DAMPE_hPull_Mixed_NS_EW_LS("relative_DAMPE_hPull_Mixed_NS_EW_LS","Pull Anisotropic relative DAMPE Sky Map (NS-EW) LS",100,-5,5);
    TH1D relative_DAMPE_hPull_Mixed_NS_FB_LS("relative_DAMPE_hPull_Mixed_NS_FB_LS","Pull Anisotropic relative DAMPE Sky Map (NS-FB) LS",100,-5,5);
    TH1D relative_DAMPE_hPull_Mixed_EW_FB_LS("relative_DAMPE_hPull_Mixed_EW_FB_LS","Pull Anisotropic relative DAMPE Sky Map (EW-FB) LS",100,-5,5);
    TH1D relative_DAMPE_hPull_FullMixed_LS("relative_DAMPE_hPull_FullMixed_LS","Pull Full Mixed Anisotropic relative DAMPE Sky Map LS",100,-5,5);
    
    //////////// High statistics pools
    
    TH1D relative_DAMPE_hPull_AniNS_HS("relative_DAMPE_hPull_AniNS_HS","Pull Anisotropic relative DAMPE Sky Map (NS) HS",50,-10,10);
    TH1D relative_DAMPE_hPull_AniEW_HS("relative_DAMPE_hPull_AniEW_HS","Pull Anisotropic relative DAMPE Sky Map (EW) HS",50,-10,10);
    TH1D relative_DAMPE_hPull_AniFB_HS("relative_DAMPE_hPull_AniFB_HS","Pull Anisotropic relative DAMPE Sky Map (FB) HS",50,-10,10);
    TH1D relative_DAMPE_hPull_Mixed_NS_EW_HS("relative_DAMPE_hPull_Mixed_NS_EW_HS","Pull Anisotropic relative DAMPE Sky Map (NS-EW) HS",50,-10,10);
    TH1D relative_DAMPE_hPull_Mixed_NS_FB_HS("relative_DAMPE_hPull_Mixed_NS_FB_HS","Pull Anisotropic relative DAMPE Sky Map (NS-FB) HS",50,-10,10);
    TH1D relative_DAMPE_hPull_Mixed_EW_FB_HS("relative_DAMPE_hPull_Mixed_EW_FB_HS","Pull Anisotropic relative DAMPE Sky Map (EW-FB) HS",50,-10,10);
    TH1D relative_DAMPE_hPull_FullMixed_HS("relative_DAMPE_hPull_FullMixed_HS","Pull Full Mixed Anisotropic relative DAMPE Sky Map HS",50,-10,10);
    
    
    //////////// 1D Histos
    
    TH1D Templates_LS[3];
    TH1D Templates_HS[3];
    
    TH1D relative_DAMPE_DataHisto_NS_LS;
    TH1D relative_DAMPE_DataHisto_EW_LS;
    TH1D relative_DAMPE_DataHisto_FB_LS;
    TH1D relative_DAMPE_MixedDataHisto_NS_EW_LS;
    TH1D relative_DAMPE_MixedDataHisto_NS_FB_LS;
    TH1D relative_DAMPE_MixedDataHisto_EW_FB_LS;
    TH1D relative_DAMPE_FullMixedDataHisto_LS;
    
    TH1D relative_DAMPE_DataHisto_NS_HS;
    TH1D relative_DAMPE_DataHisto_EW_HS;
    TH1D relative_DAMPE_DataHisto_FB_HS;
    TH1D relative_DAMPE_MixedDataHisto_NS_EW_HS;
    TH1D relative_DAMPE_MixedDataHisto_NS_FB_HS;
    TH1D relative_DAMPE_MixedDataHisto_EW_FB_HS;
    TH1D relative_DAMPE_FullMixedDataHisto_HS;
    
    TH1* Templates_1D_LS[3] = {nullptr,nullptr,nullptr};
    TH1* Templates_1D_HS[3] = {nullptr,nullptr,nullptr};
    
    TH1D* relative_DAMPE_DataHisto_1D_AniNS_LS = nullptr;
    TH1D* relative_DAMPE_DataHisto_1D_AniEW_LS = nullptr;
    TH1D* relative_DAMPE_DataHisto_1D_AniFB_LS = nullptr;
    TH1D* relative_DAMPE_MixedHisto_1D_NS_EW_LS = nullptr;
    TH1D* relative_DAMPE_MixedHisto_1D_NS_FB_LS = nullptr;
    TH1D* relative_DAMPE_MixedHisto_1D_EW_FB_LS = nullptr;
    TH1D* relative_DAMPE_FullMixedHisto_1D_LS = nullptr;
    
    TH1D* relative_DAMPE_DataHisto_1D_AniNS_HS = nullptr;
    TH1D* relative_DAMPE_DataHisto_1D_AniEW_HS = nullptr;
    TH1D* relative_DAMPE_DataHisto_1D_AniFB_HS = nullptr;
    TH1D* relative_DAMPE_MixedHisto_1D_NS_EW_HS = nullptr;
    TH1D* relative_DAMPE_MixedHisto_1D_NS_FB_HS = nullptr;
    TH1D* relative_DAMPE_MixedHisto_1D_EW_FB_HS = nullptr;
    TH1D* relative_DAMPE_FullMixedHisto_1D_HS = nullptr;
    
    
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
                               tmp_seed,
                               data_LS_events
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
                               tmp_seed,
                               data_HS_events
                           );
    
    
    ///////////////////////// Computing relative LS maps
    
    get_relative_histo(relative_DAMPE_Data_AniNS_LS,DAMPE_Data_AniNS_LS,DAMPE_ReferenceMap_LS);
    get_relative_histo(relative_DAMPE_Data_AniEW_LS,DAMPE_Data_AniEW_LS,DAMPE_ReferenceMap_LS);
    get_relative_histo(relative_DAMPE_Data_AniFB_LS,DAMPE_Data_AniFB_LS,DAMPE_ReferenceMap_LS);
    get_relative_histo(relative_DAMPE_MixedData_NS_EW_LS,DAMPE_MixedData_NS_EW_LS,DAMPE_ReferenceMap_LS);
    get_relative_histo(relative_DAMPE_MixedData_NS_FB_LS,DAMPE_MixedData_NS_FB_LS,DAMPE_ReferenceMap_LS);
    get_relative_histo(relative_DAMPE_MixedData_EW_FB_LS,DAMPE_MixedData_EW_FB_LS,DAMPE_ReferenceMap_LS);
    get_relative_histo(relative_DAMPE_FullMixedData_LS,DAMPE_FullMixedData_LS,DAMPE_ReferenceMap_LS);
    
    ///////////////////////// Computing relative HS maps
    
    get_relative_histo(relative_DAMPE_Data_AniNS_HS,DAMPE_Data_AniNS_HS,DAMPE_ReferenceMap_HS);
    get_relative_histo(relative_DAMPE_Data_AniEW_HS,DAMPE_Data_AniEW_HS,DAMPE_ReferenceMap_HS);
    get_relative_histo(relative_DAMPE_Data_AniFB_HS,DAMPE_Data_AniFB_HS,DAMPE_ReferenceMap_HS);
    get_relative_histo(relative_DAMPE_MixedData_NS_EW_HS,DAMPE_MixedData_NS_EW_HS,DAMPE_ReferenceMap_HS);
    get_relative_histo(relative_DAMPE_MixedData_NS_FB_HS,DAMPE_MixedData_NS_FB_HS,DAMPE_ReferenceMap_HS);
    get_relative_histo(relative_DAMPE_MixedData_EW_FB_HS,DAMPE_MixedData_EW_FB_HS,DAMPE_ReferenceMap_HS);
    get_relative_histo(relative_DAMPE_FullMixedData_HS,DAMPE_FullMixedData_HS,DAMPE_ReferenceMap_HS);
    
    ///////////////////////// Converting Templates Maps from TH2 to TH1
    
    TH2toTH1_obj(Templates_LS[0],Template_AniNS_LS);
    TH2toTH1_obj(Templates_LS[1],Template_AniEW_LS);
    TH2toTH1_obj(Templates_LS[2],Template_AniFB_LS);
    
    TH2toTH1_obj(Templates_HS[0],Template_AniNS_HS);
    TH2toTH1_obj(Templates_HS[1],Template_AniEW_HS);
    TH2toTH1_obj(Templates_HS[2],Template_AniFB_HS);
    
    ///////////////////////// Converting Data Maps from TH2 to TH1
    
    TH2toTH1_obj(relative_DAMPE_DataHisto_NS_LS,relative_DAMPE_Data_AniNS_LS);
    TH2toTH1_obj(relative_DAMPE_DataHisto_EW_LS,relative_DAMPE_Data_AniEW_LS);
    TH2toTH1_obj(relative_DAMPE_DataHisto_FB_LS,relative_DAMPE_Data_AniFB_LS);
    
    TH2toTH1_obj(relative_DAMPE_MixedDataHisto_NS_EW_LS,relative_DAMPE_MixedData_NS_EW_LS);
    TH2toTH1_obj(relative_DAMPE_MixedDataHisto_NS_FB_LS,relative_DAMPE_MixedData_NS_FB_LS);
    TH2toTH1_obj(relative_DAMPE_MixedDataHisto_EW_FB_LS,relative_DAMPE_MixedData_EW_FB_LS);
    TH2toTH1_obj(relative_DAMPE_FullMixedDataHisto_LS,relative_DAMPE_FullMixedData_LS);
    
    TH2toTH1_obj(relative_DAMPE_DataHisto_NS_HS,relative_DAMPE_Data_AniNS_HS);
    TH2toTH1_obj(relative_DAMPE_DataHisto_EW_HS,relative_DAMPE_Data_AniEW_HS);
    TH2toTH1_obj(relative_DAMPE_DataHisto_FB_HS,relative_DAMPE_Data_AniFB_HS);
    
    TH2toTH1_obj(relative_DAMPE_MixedDataHisto_NS_EW_HS,relative_DAMPE_MixedData_NS_EW_HS);
    TH2toTH1_obj(relative_DAMPE_MixedDataHisto_NS_FB_HS,relative_DAMPE_MixedData_NS_FB_HS);
    TH2toTH1_obj(relative_DAMPE_MixedDataHisto_EW_FB_HS,relative_DAMPE_MixedData_EW_FB_HS);
    TH2toTH1_obj(relative_DAMPE_FullMixedDataHisto_HS,relative_DAMPE_FullMixedData_HS);
    
    ///////////////////////// Linking variables to the pointers, prepearing for the fit
    
    for(Int_t idx_t = 0; idx_t < 3; idx_t++) {
        Templates_1D_LS[idx_t] = &Templates_LS[idx_t];
        Templates_1D_HS[idx_t] = &Templates_HS[idx_t];
    }
    
    relative_DAMPE_DataHisto_1D_AniNS_LS = &relative_DAMPE_DataHisto_NS_LS;
    relative_DAMPE_DataHisto_1D_AniEW_LS = &relative_DAMPE_DataHisto_EW_LS;
    relative_DAMPE_DataHisto_1D_AniFB_LS = &relative_DAMPE_DataHisto_FB_LS;
    
    relative_DAMPE_MixedHisto_1D_NS_EW_LS = &relative_DAMPE_MixedDataHisto_NS_EW_LS;
    relative_DAMPE_MixedHisto_1D_NS_FB_LS = &relative_DAMPE_MixedDataHisto_NS_FB_LS;
    relative_DAMPE_MixedHisto_1D_EW_FB_LS = &relative_DAMPE_MixedDataHisto_EW_FB_LS;
    relative_DAMPE_FullMixedHisto_1D_LS = &relative_DAMPE_FullMixedDataHisto_LS;
    
    relative_DAMPE_DataHisto_1D_AniNS_HS = &relative_DAMPE_DataHisto_NS_HS;
    relative_DAMPE_DataHisto_1D_AniEW_HS = &relative_DAMPE_DataHisto_EW_HS;
    relative_DAMPE_DataHisto_1D_AniFB_HS = &relative_DAMPE_DataHisto_FB_HS;
    
    relative_DAMPE_MixedHisto_1D_NS_EW_HS = &relative_DAMPE_MixedDataHisto_NS_EW_HS;
    relative_DAMPE_MixedHisto_1D_NS_FB_HS = &relative_DAMPE_MixedDataHisto_NS_FB_HS;
    relative_DAMPE_MixedHisto_1D_EW_FB_HS = &relative_DAMPE_MixedDataHisto_EW_FB_HS;
    relative_DAMPE_FullMixedHisto_1D_HS = &relative_DAMPE_FullMixedDataHisto_HS;
    
    /////////////////////////////////////////////////////////////////////////////// Fitting !!!
    
    std::cout << "\n\n //////////////////////////////////////// Low Statistics DAMPE fits //////////////////////////////////////// \n\n";
    output_log_file << "\n\n //////////////////////////////////////// Low Statistics DAMPE fits //////////////////////////////////////// \n\n";
    
    std::cout << "\n\n -------------- > Anisotropic Monopole (NS) \n\n";
    output_log_file << "\n\n -------------- > Anisotropic Monopole (NS) \n\n";
    
    TemplateFitBH(relative_DAMPE_DataHisto_1D_AniNS_LS,3,Templates_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file);
    getPull(relative_DAMPE_DataHisto_1D_AniNS_LS,Templates_1D_LS,res,relative_DAMPE_hPull_AniNS_LS,true);
    
    std::cout << "\n\n -------------- > Anisotropic Monopole (EW) \n\n";
    output_log_file << "\n\n -------------- > Anisotropic Monopole (EW) \n\n";
    
    TemplateFitBH(relative_DAMPE_DataHisto_1D_AniEW_LS,3,Templates_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file);
    getPull(relative_DAMPE_DataHisto_1D_AniEW_LS,Templates_1D_LS,res,relative_DAMPE_hPull_AniEW_LS,true);
    
    std::cout << "\n\n -------------- > Anisotropic Monopole (FB) \n\n";
    output_log_file << "\n\n -------------- > Anisotropic Monopole (FB) \n\n";
    
    TemplateFitBH(relative_DAMPE_DataHisto_1D_AniFB_LS,3,Templates_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file);
    getPull(relative_DAMPE_DataHisto_1D_AniFB_LS,Templates_1D_LS,res,relative_DAMPE_hPull_AniFB_LS,true);
    
    std::cout << "\n\n -------------- > NS + EW linear combination \n\n";
    output_log_file << "\n\n -------------- > NS + EW linear combination \n\n";
    
    TemplateFitBH(relative_DAMPE_MixedHisto_1D_NS_EW_LS,3,Templates_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file);
    getPull(relative_DAMPE_MixedHisto_1D_NS_EW_LS,Templates_1D_LS,res,relative_DAMPE_hPull_Mixed_NS_EW_LS,true);
    
    std::cout << "\n\n -------------- > NS + FB linear combination \n\n";
    output_log_file << "\n\n -------------- > NS + FB linear combination \n\n";
    
    TemplateFitBH(relative_DAMPE_MixedHisto_1D_NS_FB_LS,3,Templates_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file);
    getPull(relative_DAMPE_MixedHisto_1D_NS_FB_LS,Templates_1D_LS,res,relative_DAMPE_hPull_Mixed_NS_FB_LS,true);
    
    std::cout << "\n\n -------------- > EW + FB linear combination \n\n";
    output_log_file << "\n\n -------------- > EW + FB linear combination \n\n";
    
    TemplateFitBH(relative_DAMPE_MixedHisto_1D_EW_FB_LS,3,Templates_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file);
    getPull(relative_DAMPE_MixedHisto_1D_EW_FB_LS,Templates_1D_LS,res,relative_DAMPE_hPull_Mixed_EW_FB_LS,true);
    
    std::cout << "\n\n -------------- > FullMixed linear combination \n\n";
    output_log_file << "\n\n -------------- > FullMixed linear combination \n\n";
    
    TemplateFitBH(relative_DAMPE_FullMixedHisto_1D_LS,3,Templates_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file);
    getPull(relative_DAMPE_FullMixedHisto_1D_LS,Templates_1D_LS,res,relative_DAMPE_hPull_FullMixed_LS,true);
    
    
    std::cout << "\n\n //////////////////////////////////////// High Statistics DAMPE fits //////////////////////////////////////// \n\n";
    output_log_file << "\n\n //////////////////////////////////////// High Statistics DAMPE fits //////////////////////////////////////// \n\n";
    
    std::cout << "\n\n -------------- > Anisotropic Monopole (NS) \n\n";
    output_log_file << "\n\n -------------- > Anisotropic Monopole (NS) \n\n";
    
    TemplateFitBH(relative_DAMPE_DataHisto_1D_AniNS_HS,3,Templates_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file);
    getPull(relative_DAMPE_DataHisto_1D_AniNS_HS,Templates_1D_HS,res,relative_DAMPE_hPull_AniNS_HS,true);
    
    std::cout << "\n\n -------------- > Anisotropic Monopole (EW) \n\n";
    output_log_file << "\n\n -------------- > Anisotropic Monopole (EW) \n\n";
    
    TemplateFitBH(relative_DAMPE_DataHisto_1D_AniEW_HS,3,Templates_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file);
    getPull(relative_DAMPE_DataHisto_1D_AniEW_HS,Templates_1D_HS,res,relative_DAMPE_hPull_AniEW_HS,true);
    
    std::cout << "\n\n -------------- > Anisotropic Monopole (FB) \n\n";
    output_log_file << "\n\n -------------- > Anisotropic Monopole (FB) \n\n";
    
    TemplateFitBH(relative_DAMPE_DataHisto_1D_AniFB_HS,3,Templates_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file);
    getPull(relative_DAMPE_DataHisto_1D_AniFB_HS,Templates_1D_HS,res,relative_DAMPE_hPull_AniFB_HS,true);
    
    std::cout << "\n\n -------------- > NS + EW linear combination \n\n";
    output_log_file << "\n\n -------------- > NS + EW linear combination \n\n";
    
    TemplateFitBH(relative_DAMPE_MixedHisto_1D_NS_EW_HS,3,Templates_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file);
    getPull(relative_DAMPE_MixedHisto_1D_NS_EW_HS,Templates_1D_HS,res,relative_DAMPE_hPull_Mixed_NS_EW_HS,true);
    
    std::cout << "\n\n -------------- > NS + FB linear combination \n\n";
    output_log_file << "\n\n -------------- > NS + FB linear combination \n\n";
    
    TemplateFitBH(relative_DAMPE_MixedHisto_1D_NS_FB_HS,3,Templates_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file);
    getPull(relative_DAMPE_MixedHisto_1D_NS_FB_HS,Templates_1D_HS,res,relative_DAMPE_hPull_Mixed_NS_FB_HS,true);
    
    std::cout << "\n\n -------------- > EW + FB linear combination \n\n";
    output_log_file << "\n\n -------------- > EW + FB linear combination \n\n";
    
    TemplateFitBH(relative_DAMPE_MixedHisto_1D_EW_FB_HS,3,Templates_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file);
    getPull(relative_DAMPE_MixedHisto_1D_EW_FB_HS,Templates_1D_HS,res,relative_DAMPE_hPull_Mixed_EW_FB_HS,true);
    
    std::cout << "\n\n -------------- > FullMixed linear combination \n\n";
    output_log_file << "\n\n -------------- > FullMixed linear combination \n\n";
    
    TemplateFitBH(relative_DAMPE_FullMixedHisto_1D_HS,3,Templates_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file);
    getPull(relative_DAMPE_FullMixedHisto_1D_HS,Templates_1D_HS,res,relative_DAMPE_hPull_FullMixed_HS,true);
    
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    std::cout<<"\n\nSimulation Completed !\n\n";
    output_log_file << "\n\nSimulation Completed !\n\n";
    
    
    if(write_tmp_histos)
    {
        
        //////////////////////////////////// Creating Data out file
        
        TFile relative_DAMPE_data_file(data_out_path.c_str(),"RECREATE");
        if(relative_DAMPE_data_file.IsZombie()) {
            std::cout << "\n\nError writing ROOT DAMPE Data TFile. Prorgram finished \n\n";
            output_log_file << "\n\nError writing ROOT DAMPE Data TFile. Prorgram finished \n\n";
            exit(-2);
        }
        
        
        //////////////////////////////////// Writing Data
        
        relative_DAMPE_Data_AniNS_LS.Write();
        relative_DAMPE_Data_AniEW_LS.Write();
        relative_DAMPE_Data_AniFB_LS.Write();
        
        relative_DAMPE_MixedData_NS_EW_LS.Write();
        relative_DAMPE_MixedData_NS_FB_LS.Write();
        relative_DAMPE_MixedData_EW_FB_LS.Write();
        relative_DAMPE_FullMixedData_LS.Write();
        
        relative_DAMPE_Data_AniNS_HS.Write();
        relative_DAMPE_Data_AniEW_HS.Write();
        relative_DAMPE_Data_AniFB_HS.Write();
        
        relative_DAMPE_MixedData_NS_EW_HS.Write();
        relative_DAMPE_MixedData_NS_FB_HS.Write();
        relative_DAMPE_MixedData_EW_FB_HS.Write();
        relative_DAMPE_FullMixedData_HS.Write();
        
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
        
        relative_DAMPE_data_file.Write();
        relative_DAMPE_data_file.Close();
        
        
        //////////////////////////////////// Creating Pools out file
        
        TFile relative_DAMPE_pool_file(pools_out_path.c_str(),"RECREATE");
        if(relative_DAMPE_pool_file.IsZombie()) {
            std::cout << "\n\nError writing ROOT DAMPE Pools TFile. Prorgram finished \n\n";
            output_log_file << "\n\nError writing ROOT DAMPE Pools TFile. Prorgram finished \n\n";
            exit(100);
        }
        
        //////////////////////////////////// Gaussian Fit HS hPools
        
        relative_DAMPE_hPull_AniNS_HS.Fit("gaus","L,E,M");
        relative_DAMPE_hPull_AniEW_HS.Fit("gaus","L,E,M");
        relative_DAMPE_hPull_AniFB_HS.Fit("gaus","L,E,M");
        relative_DAMPE_hPull_Mixed_NS_EW_HS.Fit("gaus","L,E,M");
        relative_DAMPE_hPull_Mixed_NS_FB_HS.Fit("gaus","L,E,M");
        relative_DAMPE_hPull_Mixed_EW_FB_HS.Fit("gaus","L,E,M");
        relative_DAMPE_hPull_FullMixed_HS.Fit("gaus","L,E,M");
        
        //////////////////////////////////// Writing Pools
        
        relative_DAMPE_hPull_AniNS_LS.Write();
        relative_DAMPE_hPull_AniEW_LS.Write();
        relative_DAMPE_hPull_AniFB_LS.Write();
        relative_DAMPE_hPull_Mixed_NS_EW_LS.Write();
        relative_DAMPE_hPull_Mixed_NS_FB_LS.Write();
        relative_DAMPE_hPull_Mixed_EW_FB_LS.Write();
        relative_DAMPE_hPull_FullMixed_LS.Write();
        
        relative_DAMPE_hPull_AniNS_HS.Write();
        relative_DAMPE_hPull_AniEW_HS.Write();
        relative_DAMPE_hPull_AniFB_HS.Write();
        relative_DAMPE_hPull_Mixed_NS_EW_HS.Write();
        relative_DAMPE_hPull_Mixed_NS_FB_HS.Write();
        relative_DAMPE_hPull_Mixed_EW_FB_HS.Write();
        relative_DAMPE_hPull_FullMixed_HS.Write();
        
        relative_DAMPE_pool_file.Write();
        relative_DAMPE_pool_file.Close();
        
    }
    
}

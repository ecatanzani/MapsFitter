
#include "MyHead.h"

void interface_allSky_simulation(
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
                                 )
{
    
    std::string template_out_path = output_path_creator(output_log,output_root,1,time_stamp);
    std::string DAMPE_template_out_path = output_path_creator(output_log,output_root,1,time_stamp,0,0,0,true);
    
    std::string data_out_path = output_path_creator(
                                                        output_log,
                                                        output_root,
                                                        2,
                                                        time_stamp,
                                                        NS_anisotropy,
                                                        EW_anisotropy,
                                                        FB_anisotropy,
                                                        false
                                                    );
    
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
    
    //////////// Low statistics data
    
    TH2D* Data_Iso_LS = (TH2D*) Template_Iso_LS.Clone("Data_Iso_LS");
    TH2D* Data_AniNS_LS = (TH2D*) Template_Iso_LS.Clone("Data_AniNS_LS");
    TH2D* Data_AniEW_LS = (TH2D*) Template_Iso_LS.Clone("Data_AniEW_LS");
    TH2D* Data_AniFB_LS = (TH2D*) Template_Iso_LS.Clone("Data_AniFB_LS");
    TH2D* MixedData_NS_EW_LS = (TH2D*) Template_Iso_LS.Clone("MixedData_NS_EW_LS");
    TH2D* MixedData_NS_FB_LS = (TH2D*) Template_Iso_LS.Clone("MixedData_NS_FB_LS");
    TH2D* MixedData_EW_FB_LS = (TH2D*) Template_Iso_LS.Clone("MixedData_EW_FB_LS");
    TH2D* FullMixedData_LS = (TH2D*) Template_Iso_LS.Clone("FullMixedData_LS");
    
    Data_Iso_LS->SetTitle("Isotropic LS All Sky (data) Map; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries");
    Data_AniNS_LS->SetTitle("Anisotropic LS All Sky (data) Map (NS); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries");
    Data_AniEW_LS->SetTitle("Anisotropic LS All Sky (data) Map (EW); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries");
    Data_AniFB_LS->SetTitle("Anisotropic LS All Sky (data) Map (FB); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries");
    MixedData_NS_EW_LS->SetTitle("Anisotropic LS All Sky (NS-EW data) Map; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries");
    MixedData_NS_FB_LS->SetTitle("Anisotropic LS All Sky (NS-FB data) Map; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries");
    MixedData_EW_FB_LS->SetTitle("Anisotropic LS All Sky (EW-FB data) Map; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries");
    FullMixedData_LS->SetTitle("Anisotropic LS All Sky (NS-EW-FB data) Map; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries");
    
    Data_Iso_LS->Reset();
    Data_AniNS_LS->Reset();
    Data_AniEW_LS->Reset();
    Data_AniFB_LS->Reset();
    MixedData_NS_EW_LS->Reset();
    MixedData_NS_FB_LS->Reset();
    MixedData_EW_FB_LS->Reset();
    FullMixedData_LS->Reset();
    
    //////////// High statistics data
    
    TH2D* Data_Iso_HS = (TH2D*) Template_Iso_HS.Clone("Data_Iso_HS");
    TH2D* Data_AniNS_HS = (TH2D*) Template_Iso_HS.Clone("Data_AniNS_HS");
    TH2D* Data_AniEW_HS = (TH2D*) Template_Iso_HS.Clone("Data_AniEW_HS");
    TH2D* Data_AniFB_HS = (TH2D*) Template_Iso_HS.Clone("Data_AniFB_HS");
    TH2D* MixedData_NS_EW_HS = (TH2D*) Template_Iso_HS.Clone("MixedData_NS_EW_HS");
    TH2D* MixedData_NS_FB_HS = (TH2D*) Template_Iso_HS.Clone("MixedData_NS_FB_HS");
    TH2D* MixedData_EW_FB_HS = (TH2D*) Template_Iso_HS.Clone("MixedData_EW_FB_HS");
    TH2D* FullMixedData_HS = (TH2D*) Template_Iso_HS.Clone("FullMixedData_HS");
    
    Data_Iso_HS->SetTitle("Isotropic HS All Sky (data) Map; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries");
    Data_AniNS_HS->SetTitle("Anisotropic HS All Sky (data) Map (NS); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries");
    Data_AniEW_HS->SetTitle("Anisotropic HS All Sky (data) Map (EW); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries");
    Data_AniFB_HS->SetTitle("Anisotropic HS All Sky (data) Map (FB); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries");
    MixedData_NS_EW_HS->SetTitle("Anisotropic HS All Sky (NS-EW data) Map; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries");
    MixedData_NS_FB_HS->SetTitle("Anisotropic HS All Sky (NS-FB data) Map; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries");
    MixedData_EW_FB_HS->SetTitle("Anisotropic HS All Sky (EW-FB data) Map; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries");
    FullMixedData_HS->SetTitle("Anisotropic HS All Sky (NS-EW-FB data) Map; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries");
    
    Data_Iso_HS->Reset();
    Data_AniNS_HS->Reset();
    Data_AniEW_HS->Reset();
    Data_AniFB_HS->Reset();
    MixedData_NS_EW_HS->Reset();
    MixedData_NS_FB_HS->Reset();
    MixedData_EW_FB_HS->Reset();
    FullMixedData_HS->Reset();
    
    
    ///////////// All sky maps
    
    generate_LS_data(
                         NS_anisotropy,
                         EW_anisotropy,
                         FB_anisotropy,
                         Data_Iso_LS,
                         Data_AniNS_LS,
                         Data_AniEW_LS,
                         Data_AniFB_LS,
                         MixedData_NS_EW_LS,
                         MixedData_NS_FB_LS,
                         MixedData_EW_FB_LS,
                         FullMixedData_LS,
                         Template_Iso_LS,
                         Template_AniNS_LS,
                         Template_AniEW_LS,
                         Template_AniFB_LS,
                         output_log_file,
                         tmp_seed,
                         data_LS_events
                     );
    
    generate_HS_data(
                         NS_anisotropy,
                         EW_anisotropy,
                         FB_anisotropy,
                         Data_Iso_HS,
                         Data_AniNS_HS,
                         Data_AniEW_HS,
                         Data_AniFB_HS,
                         MixedData_NS_EW_HS,
                         MixedData_NS_FB_HS,
                         MixedData_EW_FB_HS,
                         FullMixedData_HS,
                         Template_Iso_HS,
                         Template_AniNS_HS,
                         Template_AniEW_HS,
                         Template_AniFB_HS,
                         output_log_file,
                         tmp_seed,
                         data_HS_events
                     );
 
    //////////////////////////////////// Creating Data out file
    
    TFile data_file(data_out_path.c_str(),"RECREATE");
    if(data_file.IsZombie()) {
        std::cerr << "\n\nError writing ROOT Data TFile. Prorgram finished \n\n";
        output_log_file << "\n\nError writing ROOT Data TFile. Prorgram finished \n\n";
        exit(100);
    }
    
    
    //////////////////////////////////// Writing Data
    
    Data_Iso_LS->Write();
    Data_AniNS_LS->Write();
    Data_AniEW_LS->Write();
    Data_AniFB_LS->Write();
    
    MixedData_NS_EW_LS->Write();
    MixedData_NS_FB_LS->Write();
    MixedData_EW_FB_LS->Write();
    FullMixedData_LS->Write();
    
    Data_Iso_HS->Write();
    Data_AniNS_HS->Write();
    Data_AniEW_HS->Write();
    Data_AniFB_HS->Write();
    
    MixedData_NS_EW_HS->Write();
    MixedData_NS_FB_HS->Write();
    MixedData_EW_FB_HS->Write();
    FullMixedData_HS->Write();
    
    data_file.Write();
    data_file.Close();
    
}

void interface_DAMPE_simulation(
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
                                 )
{
    
    std::string template_out_path = output_path_creator(output_log,output_root,1,time_stamp);
    std::string DAMPE_template_out_path = output_path_creator(output_log,output_root,1,time_stamp,0,0,0,true);
    
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
    
    
    ///////////// All sky maps
    
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
    
    //////////////////////////////////// Creating Data out file
    
    TFile DAMPE_data_file(data_out_path.c_str(),"RECREATE");
    if(DAMPE_data_file.IsZombie()) {
        std::cerr << "\n\nError writing ROOT DAMPE Data TFile. Prorgram finished \n\n";
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
    
}

void interface_DAMPE_relative_simulation(
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
                                         )
{
    
    std::string template_out_path = output_path_creator(output_log,output_root,1,time_stamp);
    std::string DAMPE_template_out_path = output_path_creator(output_log,output_root,1,time_stamp,0,0,0,true);
    
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
    
}

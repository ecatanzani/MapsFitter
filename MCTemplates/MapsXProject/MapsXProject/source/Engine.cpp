
#include "MyHead.h"

void read_from_file(std::string template_path,std::string data_path,std::ofstream &output_log_file) {
    /*
    std::string pools_out_path = output_path_creator(3);
    
    //////////// TemplateFit variables
    
    Double_t res[4],res_err[4];
    Double_t initialValues[4]={1,0,0,0};
    
    TH1D Templates_LS[4];
    TH1D Templates_HS[4];
    
    TH1D DataHisto_I_LS;
    TH1D DataHisto_NS_LS;
    TH1D DataHisto_EW_LS;
    TH1D DataHisto_FB_LS;
    TH1D MixedDataHisto_NS_EW_LS;
    TH1D MixedDataHisto_NS_FB_LS;
    TH1D MixedDataHisto_EW_FB_LS;
    TH1D FullMixedDataHisto_LS;
    
    TH1D DataHisto_I_HS;
    TH1D DataHisto_NS_HS;
    TH1D DataHisto_EW_HS;
    TH1D DataHisto_FB_HS;
    TH1D MixedDataHisto_NS_EW_HS;
    TH1D MixedDataHisto_NS_FB_HS;
    TH1D MixedDataHisto_EW_FB_HS;
    TH1D FullMixedDataHisto_HS;
    
    TH1* Templates_1D_LS[4] = {nullptr,nullptr,nullptr,nullptr};
    TH1* Templates_1D_HS[4] = {nullptr,nullptr,nullptr,nullptr};
    
    TH1D* DataHisto_1D_I_LS = nullptr;
    TH1D* DataHisto_1D_AniNS_LS = nullptr;
    TH1D* DataHisto_1D_AniEW_LS = nullptr;
    TH1D* DataHisto_1D_AniFB_LS = nullptr;
    TH1D* MixedHisto_1D_NS_EW_LS = nullptr;
    TH1D* MixedHisto_1D_NS_FB_LS = nullptr;
    TH1D* MixedHisto_1D_EW_FB_LS = nullptr;
    TH1D* FullMixedHisto_1D_LS = nullptr;
    
    TH1D* DataHisto_1D_I_HS = nullptr;
    TH1D* DataHisto_1D_AniNS_HS = nullptr;
    TH1D* DataHisto_1D_AniEW_HS = nullptr;
    TH1D* DataHisto_1D_AniFB_HS = nullptr;
    TH1D* MixedHisto_1D_NS_EW_HS = nullptr;
    TH1D* MixedHisto_1D_NS_FB_HS = nullptr;
    TH1D* MixedHisto_1D_EW_FB_HS = nullptr;
    TH1D* FullMixedHisto_1D_HS = nullptr;
    
    /////////////////////////////////////////////////////////////////////////////////////////
    
    load_1D_histos(Templates_LS,DataHisto_I_LS,DataHisto_NS_LS,DataHisto_EW_LS,DataHisto_FB_LS,MixedDataHisto_NS_EW_LS,MixedDataHisto_NS_FB_LS,MixedDataHisto_EW_FB_LS,FullMixedDataHisto_LS,Templates_HS,DataHisto_I_HS,DataHisto_NS_HS,DataHisto_EW_HS,DataHisto_FB_HS,MixedDataHisto_NS_EW_HS,MixedDataHisto_NS_FB_HS,MixedDataHisto_EW_FB_HS,FullMixedDataHisto_HS,template_path,data_path,output_log_file);
    
    /////////////////////////////////////////////////////////////////////////////////////////
        
    ///////////////////////// Linking variables to the pointers, prepearing for the fit
    
    for(Int_t idx_t = 0; idx_t < 4; idx_t++) {
        Templates_1D_LS[idx_t] = &Templates_LS[idx_t];
        Templates_1D_HS[idx_t] = &Templates_HS[idx_t];
    }
    
    DataHisto_1D_I_LS = &DataHisto_I_LS;
    DataHisto_1D_AniNS_LS = &DataHisto_NS_LS;
    DataHisto_1D_AniEW_LS = &DataHisto_EW_LS;
    DataHisto_1D_AniFB_LS = &DataHisto_FB_LS;
    MixedHisto_1D_NS_EW_LS = &MixedDataHisto_NS_EW_LS;
    MixedHisto_1D_NS_FB_LS = &MixedDataHisto_NS_FB_LS;
    MixedHisto_1D_EW_FB_LS = &MixedDataHisto_EW_FB_LS;
    FullMixedHisto_1D_LS = &FullMixedDataHisto_LS;
    
    DataHisto_1D_I_HS = &DataHisto_I_HS;
    DataHisto_1D_AniNS_HS = &DataHisto_NS_HS;
    DataHisto_1D_AniEW_HS = &DataHisto_EW_HS;
    DataHisto_1D_AniFB_HS = &DataHisto_FB_HS;
    MixedHisto_1D_NS_EW_HS = &MixedDataHisto_NS_EW_HS;
    MixedHisto_1D_NS_FB_HS = &MixedDataHisto_NS_FB_HS;
    MixedHisto_1D_EW_FB_HS = &MixedDataHisto_EW_FB_HS;
    FullMixedHisto_1D_HS = &FullMixedDataHisto_HS;
    
    /////////////////////////////////////////////////////////////////////////////// Fitting !!!
    
    
    std::cout << "\n\n //////////////////////////////////////// Low Statistics fits //////////////////////////////////////// \n\n";
    output_log_file << "\n\n //////////////////////////////////////// Low Statistics fits //////////////////////////////////////// \n\n";
    
    std::cout << "\n\n -------------- > Isotropic Monopole \n\n";
    output_log_file << "\n\n -------------- > Isotropic Monopole \n\n";
    
    TemplateFitBH(DataHisto_1D_I_LS,4,Templates_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file);
    getPull(DataHisto_1D_I_LS,Templates_1D_LS,res,hPull_Iso_LS);
    
    std::cout << "\n\n -------------- > Anisotropic Monopole (NS) \n\n";
    output_log_file << "\n\n -------------- > Anisotropic Monopole (NS) \n\n";
    
    TemplateFitBH(DataHisto_1D_AniNS_LS,4,Templates_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file);
    getPull(DataHisto_1D_AniNS_LS,Templates_1D_LS,res,hPull_AniNS_LS);
    
    std::cout << "\n\n -------------- > Anisotropic Monopole (EW) \n\n";
    output_log_file << "\n\n -------------- > Anisotropic Monopole (EW) \n\n";
    
    TemplateFitBH(DataHisto_1D_AniEW_LS,4,Templates_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file);
    getPull(DataHisto_1D_AniEW_LS,Templates_1D_LS,res,hPull_AniEW_LS);
    
    std::cout << "\n\n -------------- > Anisotropic Monopole (FB) \n\n";
    output_log_file << "\n\n -------------- > Anisotropic Monopole (FB) \n\n";
    
    TemplateFitBH(DataHisto_1D_AniFB_LS,4,Templates_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file);
    getPull(DataHisto_1D_AniFB_LS,Templates_1D_LS,res,hPull_AniFB_LS);
    
    std::cout << "\n\n -------------- > NS + EW linear combination \n\n";
    output_log_file << "\n\n -------------- > NS + EW linear combination \n\n";
    
    TemplateFitBH(MixedHisto_1D_NS_EW_LS,4,Templates_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file);
    getPull(MixedHisto_1D_NS_EW_LS,Templates_1D_LS,res,hPull_Mixed_NS_EW_LS);
    
    std::cout << "\n\n -------------- > NS + FB linear combination \n\n";
    output_log_file << "\n\n -------------- > NS + FB linear combination \n\n";
    
    TemplateFitBH(MixedHisto_1D_NS_FB_LS,4,Templates_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file);
    getPull(MixedHisto_1D_NS_FB_LS,Templates_1D_LS,res,hPull_Mixed_NS_FB_LS);
    
    std::cout << "\n\n -------------- > EW + FB linear combination \n\n";
    output_log_file << "\n\n -------------- > EW + FB linear combination \n\n";
    
    TemplateFitBH(MixedHisto_1D_EW_FB_LS,4,Templates_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file);
    getPull(MixedHisto_1D_EW_FB_LS,Templates_1D_LS,res,hPull_Mixed_EW_FB_LS);
    
    std::cout << "\n\n -------------- > FullMixed linear combination \n\n";
    output_log_file << "\n\n -------------- > FullMixed linear combination \n\n";
    
    TemplateFitBH(FullMixedHisto_1D_LS,4,Templates_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file);
    getPull(FullMixedHisto_1D_LS,Templates_1D_LS,res,hPull_FullMixed_LS);
    
    
    std::cout << "\n\n //////////////////////////////////////// High Statistics fits //////////////////////////////////////// \n\n";
    output_log_file << "\n\n //////////////////////////////////////// High Statistics fits //////////////////////////////////////// \n\n";
    
    std::cout << "\n\n -------------- > Isotropic Monopole \n\n";
    output_log_file << "\n\n -------------- > Isotropic Monopole \n\n";
    
    TemplateFitBH(DataHisto_1D_I_HS,4,Templates_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file);
    getPull(DataHisto_1D_I_HS,Templates_1D_HS,res,hPull_Iso_HS);
    
    std::cout << "\n\n -------------- > Anisotropic Monopole (NS) \n\n";
    output_log_file << "\n\n -------------- > Anisotropic Monopole (NS) \n\n";
    
    TemplateFitBH(DataHisto_1D_AniNS_HS,4,Templates_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file);
    getPull(DataHisto_1D_AniNS_HS,Templates_1D_HS,res,hPull_AniNS_HS);
    
    std::cout << "\n\n -------------- > Anisotropic Monopole (EW) \n\n";
    output_log_file << "\n\n -------------- > Anisotropic Monopole (EW) \n\n";
    
    TemplateFitBH(DataHisto_1D_AniEW_HS,4,Templates_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file);
    getPull(DataHisto_1D_AniEW_HS,Templates_1D_HS,res,hPull_AniEW_HS);
    
    std::cout << "\n\n -------------- > Anisotropic Monopole (FB) \n\n";
    output_log_file << "\n\n -------------- > Anisotropic Monopole (FB) \n\n";
    
    TemplateFitBH(DataHisto_1D_AniFB_HS,4,Templates_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file);
    getPull(DataHisto_1D_AniFB_HS,Templates_1D_HS,res,hPull_AniFB_HS);
    
    std::cout << "\n\n -------------- > NS + EW linear combination \n\n";
    output_log_file << "\n\n -------------- > NS + EW linear combination \n\n";
    
    TemplateFitBH(MixedHisto_1D_NS_EW_HS,4,Templates_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file);
    getPull(MixedHisto_1D_NS_EW_HS,Templates_1D_HS,res,hPull_Mixed_NS_EW_HS);
    
    std::cout << "\n\n -------------- > NS + FB linear combination \n\n";
    output_log_file << "\n\n -------------- > NS + FB linear combination \n\n";
    
    TemplateFitBH(MixedHisto_1D_NS_FB_HS,4,Templates_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file);
    getPull(MixedHisto_1D_NS_FB_HS,Templates_1D_HS,res,hPull_Mixed_NS_FB_HS);
    
    std::cout << "\n\n -------------- > EW + FB linear combination \n\n";
    output_log_file << "\n\n -------------- > EW + FB linear combination \n\n";
    
    TemplateFitBH(MixedHisto_1D_EW_FB_HS,4,Templates_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file);
    getPull(MixedHisto_1D_EW_FB_HS,Templates_1D_HS,res,hPull_Mixed_EW_FB_HS);
    
    std::cout << "\n\n -------------- > FullMixed linear combination \n\n";
    output_log_file << "\n\n -------------- > FullMixed linear combination \n\n";
    
    TemplateFitBH(FullMixedHisto_1D_HS,4,Templates_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file);
    getPull(FullMixedHisto_1D_HS,Templates_1D_HS,res,hPull_FullMixed_HS);
    
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    std::cout << "\n\nSimulation Completed !\n\n";
    output_log_file << "\n\nSimulation Completed !\n\n";
 
    //////////////////////////////////// Creating Pools out file
    
    TFile pool_file(pools_out_path.c_str(),"RECREATE");
    if(pool_file.IsZombie()) {
        std::cout << "\n\nError writing ROOT Pools TFile. Prorgram finished \n\n";
        output_log_file << "\n\nError writing ROOT Pools TFile. Prorgram finished \n\n";
        exit(-2);
    }
    
    
    //////////////////////////////////// Writing Pools
    
    hPull_Iso_LS.Write()
    hPull_AniNS_LS.Write();
    hPull_AniEW_LS.Write();
    hPull_AniFB_LS.Write();
    hPull_Mixed_NS_EW_LS.Write();
    hPull_Mixed_NS_FB_LS.Write();
    hPull_Mixed_EW_FB_LS.Write();
    hPull_FullMixed_LS.Write();
    
    hPull_Iso_HS.Write()
    hPull_AniNS_HS.Write();
    hPull_AniEW_HS.Write();
    hPull_AniFB_HS.Write();
    hPull_Mixed_HS_EW_LS.Write();
    hPull_Mixed_HS_FB_LS.Write();
    hPull_Mixed_HW_FB_LS.Write();
    hPull_FullMixed_HS.Write();
    
    pool_file.Write();
    pool_file.Close();
     */
}

void generate_and_fit(std::ofstream &output_log_file) {
    
    ////////////////////////////////////////////////////////// Variables declaration ///////////////////////////////////////////////////
    
    std::string template_out_path = output_path_creator(1);
    std::string data_out_path = output_path_creator(2);
    std::string pools_out_path = output_path_creator(3);
    TRandom3 r_gen(random_seed);
    
    //////////// Costheta flat binning variables
    
    Int_t n_bin_lon_LS = 36;
    Double_t lon_bin_min_LS = -180;
    Double_t lon_bin_max_LS = 180;
    
    Int_t n_bin_lat_LS = 18;
    Double_t lat_bin_min_LS = -90;
    Double_t lat_bin_max_LS = 90;
    
    Double_t* binning_LS = nullptr;
    
    Int_t n_bin_lon_HS = 36;
    Double_t lon_bin_min_HS = -180;
    Double_t lon_bin_max_HS = 180;
    
    Int_t n_bin_lat_HS = 18;
    Double_t lat_bin_min_HS = -90;
    Double_t lat_bin_max_HS = 90;
    
    Double_t* binning_HS = nullptr;
    
    //////////// TemplateFit variables
    
    Double_t res[4],res_err[4];
    Double_t initialValues[4]={1,0,0,0};
    
    TH1D Templates_LS[4];
    TH1D Templates_HS[4];
    
    TH1D DataHisto_I_LS;
    TH1D DataHisto_NS_LS;
    TH1D DataHisto_EW_LS;
    TH1D DataHisto_FB_LS;
    TH1D MixedDataHisto_NS_EW_LS;
    TH1D MixedDataHisto_NS_FB_LS;
    TH1D MixedDataHisto_EW_FB_LS;
    TH1D FullMixedDataHisto_LS;
    
    TH1D DataHisto_I_HS;
    TH1D DataHisto_NS_HS;
    TH1D DataHisto_EW_HS;
    TH1D DataHisto_FB_HS;
    TH1D MixedDataHisto_NS_EW_HS;
    TH1D MixedDataHisto_NS_FB_HS;
    TH1D MixedDataHisto_EW_FB_HS;
    TH1D FullMixedDataHisto_HS;
    
    TH1* Templates_1D_LS[4] = {nullptr,nullptr,nullptr,nullptr};
    TH1* Templates_1D_HS[4] = {nullptr,nullptr,nullptr,nullptr};
    
    TH1D* DataHisto_1D_I_LS = nullptr;
    TH1D* DataHisto_1D_AniNS_LS = nullptr;
    TH1D* DataHisto_1D_AniEW_LS = nullptr;
    TH1D* DataHisto_1D_AniFB_LS = nullptr;
    TH1D* MixedHisto_1D_NS_EW_LS = nullptr;
    TH1D* MixedHisto_1D_NS_FB_LS = nullptr;
    TH1D* MixedHisto_1D_EW_FB_LS = nullptr;
    TH1D* FullMixedHisto_1D_LS = nullptr;
    
    TH1D* DataHisto_1D_I_HS = nullptr;
    TH1D* DataHisto_1D_AniNS_HS = nullptr;
    TH1D* DataHisto_1D_AniEW_HS = nullptr;
    TH1D* DataHisto_1D_AniFB_HS = nullptr;
    TH1D* MixedHisto_1D_NS_EW_HS = nullptr;
    TH1D* MixedHisto_1D_NS_FB_HS = nullptr;
    TH1D* MixedHisto_1D_EW_FB_HS = nullptr;
    TH1D* FullMixedHisto_1D_HS = nullptr;
     
    TCanvas LS_Templates_Canvas("LS_Templates_Canvas","Low Statistics Templates");
    TCanvas HS_Templates_Canvas("HS_Templates_Canvas","High Statistics Templates");
    
    LS_Templates_Canvas.Divide(2,2);
    HS_Templates_Canvas.Divide(2,2);
    
    
    //////////// Histos
    
    create_binning(n_bin_lat_LS,lat_bin_min_LS,lat_bin_max_LS,binning_LS,true);
    create_binning(n_bin_lat_HS,lat_bin_min_HS,lat_bin_max_HS,binning_HS,true);
    
    
    ///////////////////////////////////////// All Sky Histos
    
    /////////// Low statistics template
    
    TH2D Template_Iso_LS("Template_Iso_LS","Isotropic LS Template All Sky Map; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon_LS,lon_bin_min_LS,lon_bin_max_LS,n_bin_lat_LS,binning_LS);
    TH2D Template_AniNS_LS("Template_AniNS_LS","Anisotropic LS Template All Sky Map (NS); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon_LS,lon_bin_min_LS,lon_bin_max_LS,n_bin_lat_LS,binning_LS);
    TH2D Template_AniEW_LS("Template_AniEW_LS","Anisotropic LS Template All Sky Map (EW); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon_LS,lon_bin_min_LS,lon_bin_max_LS,n_bin_lat_LS,binning_LS);
    TH2D Template_AniFB_LS("Template_AniFB_LS","Anisotropic LS Template All Sky Map (FB); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon_LS,lon_bin_min_LS,lon_bin_max_LS,n_bin_lat_LS,binning_LS);
    
    /////////// High statistics template
    
    TH2D Template_Iso_HS("Template_Iso_HS","Isotropic HS Template All Sky Map; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon_HS,lon_bin_min_HS,lon_bin_max_HS,n_bin_lat_HS,binning_HS);
    TH2D Template_AniNS_HS("Template_AniNS_HS","Anisotropic HS Template All Sky Map (NS); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon_HS,lon_bin_min_HS,lon_bin_max_HS,n_bin_lat_HS,binning_HS);
    TH2D Template_AniEW_HS("Template_AniEW_HS","Anisotropic HS Template All Sky Map (EW); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon_HS,lon_bin_min_HS,lon_bin_max_HS,n_bin_lat_HS,binning_HS);
    TH2D Template_AniFB_HS("Template_AniFB_HS","Anisotropic HS Template All Sky Map (FB); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon_HS,lon_bin_min_HS,lon_bin_max_HS,n_bin_lat_HS,binning_HS);
    
    /////////// Low statistics data
    
    TH2D Data_Iso_LS("Data_Iso_LS","Isotropic LS All Sky (data) Map; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon_LS,lon_bin_min_LS,lon_bin_max_LS,n_bin_lat_LS,binning_LS);
    TH2D Data_AniNS_LS("Data_AniNS_LS","Anisotropic LS All Sky (data) Map (NS); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon_LS,lon_bin_min_LS,lon_bin_max_LS,n_bin_lat_LS,binning_LS);
    TH2D Data_AniEW_LS("Data_AniEW_LS","Anisotropic LS All Sky (data) Map (EW); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon_LS,lon_bin_min_LS,lon_bin_max_LS,n_bin_lat_LS,binning_LS);
    TH2D Data_AniFB_LS("Data_AniFB_LS","Anisotropic LS All Sky (data) Map (FB); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon_LS,lon_bin_min_LS,lon_bin_max_LS,n_bin_lat_LS,binning_LS);
    
    TH2D MixedData_NS_EW_LS("MixedData_NS_EW_LS","Anisotropic LS All Sky (NS-EW data) Map; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon_LS,lon_bin_min_LS,lon_bin_max_LS,n_bin_lat_LS,binning_LS);
    TH2D MixedData_NS_FB_LS("MixedData_NS_FB_LS","Anisotropic LS All Sky (NS-FB data) Map; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon_LS,lon_bin_min_LS,lon_bin_max_LS,n_bin_lat_LS,binning_LS);
    TH2D MixedData_EW_FB_LS("MixedData_EW_FB_LS","Anisotropic LS All Sky (EW-FB data) Map; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon_LS,lon_bin_min_LS,lon_bin_max_LS,n_bin_lat_LS,binning_LS);
    TH2D FullMixedData_LS("FullMixedData_LS","Anisotropic LS All Sky (NS-EW-FB data) Map; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon_LS,lon_bin_min_LS,lon_bin_max_LS,n_bin_lat_LS,binning_LS);
    
    /////////// High statistics data
    
    TH2D Data_Iso_HS("Data_Iso_HS","Isotropic HS All Sky (data) Map; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon_HS,lon_bin_min_HS,lon_bin_max_HS,n_bin_lat_HS,binning_HS);
    TH2D Data_AniNS_HS("Data_AniNS_HS","Anisotropic HS All Sky (data) Map (NS); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon_HS,lon_bin_min_HS,lon_bin_max_HS,n_bin_lat_HS,binning_HS);
    TH2D Data_AniEW_HS("Data_AniEW_HS","Anisotropic HS All Sky (data) Map (EW); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon_HS,lon_bin_min_HS,lon_bin_max_HS,n_bin_lat_HS,binning_HS);
    TH2D Data_AniFB_HS("Data_AniFB_HS","Anisotropic HS All Sky (data) Map (FB); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon_HS,lon_bin_min_HS,lon_bin_max_HS,n_bin_lat_HS,binning_HS);
    
    TH2D MixedData_NS_EW_HS("MixedData_NS_EW_HS","Anisotropic HS All Sky (NS-EW data) Map; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon_HS,lon_bin_min_HS,lon_bin_max_HS,n_bin_lat_HS,binning_HS);
    TH2D MixedData_NS_FB_HS("MixedData_NS_FB_HS","Anisotropic HS All Sky (NS-FB data) Map; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon_HS,lon_bin_min_HS,lon_bin_max_HS,n_bin_lat_HS,binning_HS);
    TH2D MixedData_EW_FB_HS("MixedData_EW_FB_HS","Anisotropic HS All Sky (EW-FB data) Map; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon_HS,lon_bin_min_HS,lon_bin_max_HS,n_bin_lat_HS,binning_HS);
    TH2D FullMixedData_HS("FullMixedData_HS","Anisotropic HS All Sky (NS-EW-FB data) Map; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon_HS,lon_bin_min_HS,lon_bin_max_HS,n_bin_lat_HS,binning_HS);
    
    //////////// Stellite's Histos
    
    
    //////////// Weight Histos
    
    TH1D Template_hwNS("template_hWNS","NS Weight",100,-1,1);
    TH1D Template_hwEW("template_hWEW","EW Weight",100,-1,1);
    TH1D Template_hwFB("template_hWFB","FB Weight",100,-1,1);
    
    TH1D Data_hwNS("data_hWNS","NS Weight",100,-2,2);
    TH1D Data_hwEW("data_hWEW","EW Weight",100,-2,2);
    TH1D Data_hwFB("data_hWFB","FB Weight",100,-2,2);
    
    
    ///////////////////////////////////////// Pool's histos
    
    /////////// Low statistics pools
    
    TH1D hPull_Iso_LS("hPull_Iso_LS","Pull Isotropic Sky Map LS",100,-5,5);
    TH1D hPull_AniNS_LS("hPull_AniNS_LS","Pull Anisotropic Sky Map (NS) LS",100,-5,5);
    TH1D hPull_AniEW_LS("hPull_AniEW_LS","Pull Anisotropic Sky Map (EW) LS",100,-5,5);
    TH1D hPull_AniFB_LS("hPull_AniFB_LS","Pull Anisotropic Sky Map (FB) LS",100,-5,5);
    TH1D hPull_Mixed_NS_EW_LS("hPull_Mixed_NS_EW_LS","Pull Anisotropic Sky Map (NS-EW) LS",100,-5,5);
    TH1D hPull_Mixed_NS_FB_LS("hPull_Mixed_NS_FB_LS","Pull Anisotropic Sky Map (NS-FB) LS",100,-5,5);
    TH1D hPull_Mixed_EW_FB_LS("hPull_Mixed_EW_FB_LS","Pull Anisotropic Sky Map (EW-FB) LS",100,-5,5);
    TH1D hPull_FullMixed_LS("hPull_FullMixed_LS","Pull Full Mixed Anisotropic Sky Map LS",100,-5,5);
    
    /////////// High statistics pools
    
    TH1D hPull_Iso_HS("hPull_Iso_HS","Pull Isotropic Sky Map HS",50,-0.2,0.2);
    TH1D hPull_AniNS_HS("hPull_AniNS_HS","Pull Anisotropic Sky Map (NS) HS",50,-0.2,0.2);
    TH1D hPull_AniEW_HS("hPull_AniEW_HS","Pull Anisotropic Sky Map (EW) HS",50,-0.2,0.2);
    TH1D hPull_AniFB_HS("hPull_AniFB_HS","Pull Anisotropic Sky Map (FB) HS",50,-0.2,0.2);
    TH1D hPull_Mixed_NS_EW_HS("hPull_Mixed_NS_EW_HS","Pull Anisotropic Sky Map (NS-EW) HS",50,-0.2,0.2);
    TH1D hPull_Mixed_NS_FB_HS("hPull_Mixed_NS_FB_HS","Pull Anisotropic Sky Map (NS-FB) HS",50,-0.2,0.2);
    TH1D hPull_Mixed_EW_FB_HS("hPull_Mixed_EW_FB_HS","Pull Anisotropic Sky Map (EW-FB) HS",50,-0.2,0.2);
    TH1D hPull_FullMixed_HS("hPull_FullMixed_HS","Pull Full Mixed Anisotropic Sky Map HS",50,-0.2,0.2);
    
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    create_and_initialize_log(output_log_file);
    
    ///////////// All sky maps
    
    generate_LS_templates(Template_Iso_LS,Template_AniNS_LS,Template_AniEW_LS,Template_AniFB_LS,Template_hwNS,Template_hwEW,Template_hwFB,output_log_file,r_gen);
    generate_HS_templates(Template_Iso_HS,Template_AniNS_HS,Template_AniEW_HS,Template_AniFB_HS,Template_hwNS,Template_hwEW,Template_hwFB,output_log_file,r_gen);
    
    generate_LS_data(Data_Iso_LS,Data_AniNS_LS,Data_AniEW_LS,Data_AniFB_LS,MixedData_NS_EW_LS,MixedData_NS_FB_LS,MixedData_EW_FB_LS,FullMixedData_LS,Data_hwNS,Data_hwEW,Data_hwFB,output_log_file,r_gen);
    generate_HS_data(Data_Iso_HS,Data_AniNS_HS,Data_AniEW_HS,Data_AniFB_HS,MixedData_NS_EW_HS,MixedData_NS_FB_HS,MixedData_EW_FB_HS,FullMixedData_HS,Data_hwNS,Data_hwEW,Data_hwFB,output_log_file,r_gen);
    
    //normalize_LS_templates(Template_Iso_LS,Template_AniNS_LS,Template_AniEW_LS,Template_AniFB_LS,output_log_file);
    //normalize_HS_templates(Template_Iso_HS,Template_AniNS_HS,Template_AniEW_HS,Template_AniFB_HS,output_log_file);
    
    //normalize_LS_data(Data_Iso_LS,Data_AniNS_LS,Data_AniEW_LS,Data_AniFB_LS,MixedData_NS_EW_LS,MixedData_NS_FB_LS,MixedData_EW_FB_LS,FullMixedData_LS,output_log_file);
    //normalize_HS_data(Data_Iso_HS,Data_AniNS_HS,Data_AniEW_HS,Data_AniFB_HS,MixedData_NS_EW_HS,MixedData_NS_FB_HS,MixedData_EW_FB_HS,FullMixedData_HS,output_log_file);
    
    ///////////////////////// Converting Templates Maps from TH2 to TH1
    
    TH2toTH1(Templates_LS[0],Template_Iso_LS);
    TH2toTH1(Templates_LS[1],Template_AniNS_LS);
    TH2toTH1(Templates_LS[2],Template_AniEW_LS);
    TH2toTH1(Templates_LS[3],Template_AniFB_LS);
    
    TH2toTH1(Templates_HS[0],Template_Iso_HS);
    TH2toTH1(Templates_HS[1],Template_AniNS_HS);
    TH2toTH1(Templates_HS[2],Template_AniEW_HS);
    TH2toTH1(Templates_HS[3],Template_AniFB_HS);
    
    ///////////////////////// Converting Data Maps from TH2 to TH1
    
    TH2toTH1(DataHisto_I_LS,Data_Iso_LS);
    TH2toTH1(DataHisto_NS_LS,Data_AniNS_LS);
    TH2toTH1(DataHisto_EW_LS,Data_AniEW_LS);
    TH2toTH1(DataHisto_FB_LS,Data_AniFB_LS);
    TH2toTH1(MixedDataHisto_NS_EW_LS,MixedData_NS_EW_LS);
    TH2toTH1(MixedDataHisto_NS_FB_LS,MixedData_NS_FB_LS);
    TH2toTH1(MixedDataHisto_EW_FB_LS,MixedData_EW_FB_LS);
    TH2toTH1(FullMixedDataHisto_LS,FullMixedData_LS);
    
    TH2toTH1(DataHisto_I_HS,Data_Iso_HS);
    TH2toTH1(DataHisto_NS_HS,Data_AniNS_HS);
    TH2toTH1(DataHisto_EW_HS,Data_AniEW_HS);
    TH2toTH1(DataHisto_FB_HS,Data_AniFB_HS);
    TH2toTH1(MixedDataHisto_NS_EW_HS,MixedData_NS_EW_HS);
    TH2toTH1(MixedDataHisto_NS_FB_HS,MixedData_NS_FB_HS);
    TH2toTH1(MixedDataHisto_EW_FB_HS,MixedData_EW_FB_HS);
    TH2toTH1(FullMixedDataHisto_HS,FullMixedData_HS);
     
    ///////////////////////// Drawing Templates Maps to canvas
    
    for(Int_t idx_comp = 0; idx_comp < 4; idx_comp++) {
        LS_Templates_Canvas.cd(idx_comp+1);
        Templates_LS[idx_comp].Draw();
        HS_Templates_Canvas.cd(idx_comp+1);
        Templates_HS[idx_comp].Draw();
    }
    
    ///////////////////////// Linking variables to the pointers, prepearing for the fit
    
    for(Int_t idx_t = 0; idx_t < 4; idx_t++) {
        Templates_1D_LS[idx_t] = &Templates_LS[idx_t];
        Templates_1D_HS[idx_t] = &Templates_HS[idx_t];
    }
    
    DataHisto_1D_I_LS = &DataHisto_I_LS;
    DataHisto_1D_AniNS_LS = &DataHisto_NS_LS;
    DataHisto_1D_AniEW_LS = &DataHisto_EW_LS;
    DataHisto_1D_AniFB_LS = &DataHisto_FB_LS;
    MixedHisto_1D_NS_EW_LS = &MixedDataHisto_NS_EW_LS;
    MixedHisto_1D_NS_FB_LS = &MixedDataHisto_NS_FB_LS;
    MixedHisto_1D_EW_FB_LS = &MixedDataHisto_EW_FB_LS;
    FullMixedHisto_1D_LS = &FullMixedDataHisto_LS;
    
    DataHisto_1D_I_HS = &DataHisto_I_HS;
    DataHisto_1D_AniNS_HS = &DataHisto_NS_HS;
    DataHisto_1D_AniEW_HS = &DataHisto_EW_HS;
    DataHisto_1D_AniFB_HS = &DataHisto_FB_HS;
    MixedHisto_1D_NS_EW_HS = &MixedDataHisto_NS_EW_HS;
    MixedHisto_1D_NS_FB_HS = &MixedDataHisto_NS_FB_HS;
    MixedHisto_1D_EW_FB_HS = &MixedDataHisto_EW_FB_HS;
    FullMixedHisto_1D_HS = &FullMixedDataHisto_HS;
    
    /////////////////////////////////////////////////////////////////////////////// Fitting !!!
    
    
    std::cout << "\n\n //////////////////////////////////////// Low Statistics fits //////////////////////////////////////// \n\n";
    output_log_file << "\n\n //////////////////////////////////////// Low Statistics fits //////////////////////////////////////// \n\n";
    
    std::cout << "\n\n -------------- > Isotropic Monopole \n\n";
    output_log_file << "\n\n -------------- > Isotropic Monopole \n\n";
    
    TemplateFitBH(DataHisto_1D_I_LS,4,Templates_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file);
    getPull(DataHisto_1D_I_LS,Templates_1D_LS,res,hPull_Iso_LS);
    
    std::cout << "\n\n -------------- > Anisotropic Monopole (NS) \n\n";
    output_log_file << "\n\n -------------- > Anisotropic Monopole (NS) \n\n";
    
    TemplateFitBH(DataHisto_1D_AniNS_LS,4,Templates_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file);
    getPull(DataHisto_1D_AniNS_LS,Templates_1D_LS,res,hPull_AniNS_LS);
    
    std::cout << "\n\n -------------- > Anisotropic Monopole (EW) \n\n";
    output_log_file << "\n\n -------------- > Anisotropic Monopole (EW) \n\n";
    
    TemplateFitBH(DataHisto_1D_AniEW_LS,4,Templates_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file);
    getPull(DataHisto_1D_AniEW_LS,Templates_1D_LS,res,hPull_AniEW_LS);
    
    std::cout << "\n\n -------------- > Anisotropic Monopole (FB) \n\n";
    output_log_file << "\n\n -------------- > Anisotropic Monopole (FB) \n\n";
    
    TemplateFitBH(DataHisto_1D_AniFB_LS,4,Templates_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file);
    getPull(DataHisto_1D_AniFB_LS,Templates_1D_LS,res,hPull_AniFB_LS);
    
    std::cout << "\n\n -------------- > NS + EW linear combination \n\n";
    output_log_file << "\n\n -------------- > NS + EW linear combination \n\n";
    
    TemplateFitBH(MixedHisto_1D_NS_EW_LS,4,Templates_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file);
    getPull(MixedHisto_1D_NS_EW_LS,Templates_1D_LS,res,hPull_Mixed_NS_EW_LS);
    
    std::cout << "\n\n -------------- > NS + FB linear combination \n\n";
    output_log_file << "\n\n -------------- > NS + FB linear combination \n\n";
    
    TemplateFitBH(MixedHisto_1D_NS_FB_LS,4,Templates_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file);
    getPull(MixedHisto_1D_NS_FB_LS,Templates_1D_LS,res,hPull_Mixed_NS_FB_LS);
    
    std::cout << "\n\n -------------- > EW + FB linear combination \n\n";
    output_log_file << "\n\n -------------- > EW + FB linear combination \n\n";
    
    TemplateFitBH(MixedHisto_1D_EW_FB_LS,4,Templates_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file);
    getPull(MixedHisto_1D_EW_FB_LS,Templates_1D_LS,res,hPull_Mixed_EW_FB_LS);
    
    std::cout << "\n\n -------------- > FullMixed linear combination \n\n";
    output_log_file << "\n\n -------------- > FullMixed linear combination \n\n";
    
    TemplateFitBH(FullMixedHisto_1D_LS,4,Templates_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file);
    getPull(FullMixedHisto_1D_LS,Templates_1D_LS,res,hPull_FullMixed_LS);
    
    
    std::cout << "\n\n //////////////////////////////////////// High Statistics fits //////////////////////////////////////// \n\n";
    output_log_file << "\n\n //////////////////////////////////////// High Statistics fits //////////////////////////////////////// \n\n";
    
    std::cout << "\n\n -------------- > Isotropic Monopole \n\n";
    output_log_file << "\n\n -------------- > Isotropic Monopole \n\n";
    
    TemplateFitBH(DataHisto_1D_I_HS,4,Templates_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file);
    getPull(DataHisto_1D_I_HS,Templates_1D_HS,res,hPull_Iso_HS);
    
    std::cout << "\n\n -------------- > Anisotropic Monopole (NS) \n\n";
    output_log_file << "\n\n -------------- > Anisotropic Monopole (NS) \n\n";
    
    TemplateFitBH(DataHisto_1D_AniNS_HS,4,Templates_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file);
    getPull(DataHisto_1D_AniNS_HS,Templates_1D_HS,res,hPull_AniNS_HS);
    
    std::cout << "\n\n -------------- > Anisotropic Monopole (EW) \n\n";
    output_log_file << "\n\n -------------- > Anisotropic Monopole (EW) \n\n";
    
    TemplateFitBH(DataHisto_1D_AniEW_HS,4,Templates_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file);
    getPull(DataHisto_1D_AniEW_HS,Templates_1D_HS,res,hPull_AniEW_HS);
    
    std::cout << "\n\n -------------- > Anisotropic Monopole (FB) \n\n";
    output_log_file << "\n\n -------------- > Anisotropic Monopole (FB) \n\n";
    
    TemplateFitBH(DataHisto_1D_AniFB_HS,4,Templates_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file);
    getPull(DataHisto_1D_AniFB_HS,Templates_1D_HS,res,hPull_AniFB_HS);
    
    std::cout << "\n\n -------------- > NS + EW linear combination \n\n";
    output_log_file << "\n\n -------------- > NS + EW linear combination \n\n";
    
    TemplateFitBH(MixedHisto_1D_NS_EW_HS,4,Templates_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file);
    getPull(MixedHisto_1D_NS_EW_HS,Templates_1D_HS,res,hPull_Mixed_NS_EW_HS);
    
    std::cout << "\n\n -------------- > NS + FB linear combination \n\n";
    output_log_file << "\n\n -------------- > NS + FB linear combination \n\n";
    
    TemplateFitBH(MixedHisto_1D_NS_FB_HS,4,Templates_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file);
    getPull(MixedHisto_1D_NS_FB_HS,Templates_1D_HS,res,hPull_Mixed_NS_FB_HS);
    
    std::cout << "\n\n -------------- > EW + FB linear combination \n\n";
    output_log_file << "\n\n -------------- > EW + FB linear combination \n\n";
    
    TemplateFitBH(MixedHisto_1D_EW_FB_HS,4,Templates_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file);
    getPull(MixedHisto_1D_EW_FB_HS,Templates_1D_HS,res,hPull_Mixed_EW_FB_HS);
    
    std::cout << "\n\n -------------- > FullMixed linear combination \n\n";
    output_log_file << "\n\n -------------- > FullMixed linear combination \n\n";
    
    TemplateFitBH(FullMixedHisto_1D_HS,4,Templates_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file);
    getPull(FullMixedHisto_1D_HS,Templates_1D_HS,res,hPull_FullMixed_HS);
    
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    std::cout<<"\n\nSimulation Completed !\n\n";
    output_log_file << "\n\nSimulation Completed !\n\n";
    
    //////////////////////////////////// Creating Template out file
    
    TFile template_file(template_out_path.c_str(),"RECREATE");
    if(template_file.IsZombie()) {
        std::cout << "\n\nError writing ROOT Templates TFile. Prorgram finished \n\n";
        output_log_file << "\n\nError writing ROOT Templates TFile. Prorgram finished \n\n";
        exit(-2);
    }
    
    //////////////////////////////////// Writing Templates
    
    Template_Iso_LS.Write();
    Template_AniNS_LS.Write();
    Template_AniEW_LS.Write();
    Template_AniFB_LS.Write();
    
    Template_Iso_HS.Write();
    Template_AniNS_HS.Write();
    Template_AniEW_HS.Write();
    Template_AniFB_HS.Write();
    
    Template_hwNS.Write();
    Template_hwEW.Write();
    Template_hwFB.Write();
    
    LS_Templates_Canvas.Write();
    HS_Templates_Canvas.Write();
    
    template_file.Write();
    template_file.Close();
    
    
    
    //////////////////////////////////// Creating Data out file
    
    TFile data_file(data_out_path.c_str(),"RECREATE");
    if(data_file.IsZombie()) {
        std::cout << "\n\nError writing ROOT Data TFile. Prorgram finished \n\n";
        output_log_file << "\n\nError writing ROOT Data TFile. Prorgram finished \n\n";
        exit(-2);
    }
    
    
    //////////////////////////////////// Writing Data
    
    Data_Iso_LS.Write();
    Data_AniNS_LS.Write();
    Data_AniEW_LS.Write();
    Data_AniFB_LS.Write();
    
    MixedData_NS_EW_LS.Write();
    MixedData_NS_FB_LS.Write();
    MixedData_EW_FB_LS.Write();
    FullMixedData_LS.Write();
    
    Data_Iso_HS.Write();
    Data_AniNS_HS.Write();
    Data_AniEW_HS.Write();
    Data_AniFB_HS.Write();
    
    MixedData_NS_EW_HS.Write();
    MixedData_NS_FB_HS.Write();
    MixedData_EW_FB_HS.Write();
    FullMixedData_HS.Write();
    
    Data_hwNS.Write();
    Data_hwEW.Write();
    Data_hwFB.Write();
    
    data_file.Write();
    data_file.Close();
    
    
    //////////////////////////////////// Creating Pools out file
    
    TFile pool_file(pools_out_path.c_str(),"RECREATE");
    if(pool_file.IsZombie()) {
        std::cout << "\n\nError writing ROOT Pools TFile. Prorgram finished \n\n";
        output_log_file << "\n\nError writing ROOT Pools TFile. Prorgram finished \n\n";
        exit(-2);
    }
    
    //////////////////////////////////// Gaussian Fit HS hPools
    
    hPull_Iso_HS.Fit("gaus","L,E,M");
    hPull_AniNS_HS.Fit("gaus","L,E,M");
    hPull_AniEW_HS.Fit("gaus","L,E,M");
    hPull_AniFB_HS.Fit("gaus","L,E,M");
    hPull_Mixed_NS_EW_HS.Fit("gaus","L,E,M");
    hPull_Mixed_NS_FB_HS.Fit("gaus","L,E,M");
    hPull_Mixed_EW_FB_HS.Fit("gaus","L,E,M");
    hPull_FullMixed_HS.Fit("gaus","L,E,M");
    
    //////////////////////////////////// Writing Pools
    
    hPull_Iso_LS.Write();
    hPull_AniNS_LS.Write();
    hPull_AniEW_LS.Write();
    hPull_AniFB_LS.Write();
    hPull_Mixed_NS_EW_LS.Write();
    hPull_Mixed_NS_FB_LS.Write();
    hPull_Mixed_EW_FB_LS.Write();
    hPull_FullMixed_LS.Write();
    
    hPull_Iso_HS.Write();
    hPull_AniNS_HS.Write();
    hPull_AniEW_HS.Write();
    hPull_AniFB_HS.Write();
    hPull_Mixed_NS_EW_HS.Write();
    hPull_Mixed_NS_FB_HS.Write();
    hPull_Mixed_EW_FB_HS.Write();
    hPull_FullMixed_HS.Write();
    
    pool_file.Write();
    pool_file.Close();
    
}

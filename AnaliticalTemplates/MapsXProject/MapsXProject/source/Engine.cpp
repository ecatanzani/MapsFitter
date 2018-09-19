
#include "MyHead.h"


//////////////////// TF2 Template Functions ////////////////////


Double_t IsoMonopole(Double_t *val, Double_t *par) {
    Double_t f = 1/TMath::Sqrt(4*TMath::Pi());
    return f;
}

Double_t NSMonopole(Double_t *val, Double_t *par) {
    Double_t b = val[1];
    Double_t f =  - 0.5*TMath::Sqrt(3/TMath::Pi())*TMath::Sin(b);
    return f;
}

Double_t EWMonopole(Double_t *val, Double_t *par) {
    Double_t l = val[0];
    Double_t b = val[1];
    Double_t f =  - (1/TMath::Sqrt(2))*TMath::Sqrt(3/(2*TMath::Pi()))*TMath::Cos(b)*TMath::Sin(l);
    return f;
}

Double_t FBMonopole(Double_t *val, Double_t *par) {
    Double_t l = val[0];
    Double_t b = val[1];
    Double_t f =  (1/TMath::Sqrt(2))*TMath::Sqrt(3/(2*TMath::Pi()))*TMath::Cos(b)*TMath::Cos(l);
    return f;
}

///////////////////////////////////////////////////////////////////////////////////////////

#if 0
void generate_data(std::ofstream &log_file,std::string &data_path) {
    
    data_path = output_path_creator(3);
    
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
    
    create_binning(n_bin_lat_LS,lat_bin_min_LS,lat_bin_max_LS,binning_LS,true);
    create_binning(n_bin_lat_HS,lat_bin_min_HS,lat_bin_max_HS,binning_HS,true);
    
    //////////// TF2
    
    static TF2 dI("dI",IsoMonopole,-TMath::Pi(),TMath::Pi(),-TMath::Pi()/2,TMath::Pi()/2);
    static TF2 dNS("dNS",NSMonopole,-TMath::Pi(),TMath::Pi(),-TMath::Pi()/2,TMath::Pi()/2);
    static TF2 dEW("dEW",EWMonopole,-TMath::Pi(),TMath::Pi(),-TMath::Pi()/2,TMath::Pi()/2);
    static TF2 dFB("dFB",FBMonopole,-TMath::Pi(),TMath::Pi(),-TMath::Pi()/2,TMath::Pi()/2);
    
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
    
       //////////// Weight Histos
    
    TH1D Data_hwI("data_hWI","Iso Weight",100,0,1);
    TH1D Data_hwNS("data_hWNS","NS Weight",100,-2,2);
    TH1D Data_hwEW("data_hWEW","EW Weight",100,-2,2);
    TH1D Data_hwFB("data_hWFB","FB Weight",100,-2,2);
    
   //////////////////////////////////////////////////////////////////
    
    TRandom3 r_gen(random_seed);

    /*
    generate_LS_data(
                     Data_Iso_LS,
                     Data_AniNS_LS,
                     Data_AniEW_LS,
                     Data_AniFB_LS,
                     MixedData_NS_EW_LS,
                     MixedData_NS_FB_LS,
                     MixedData_EW_FB_LS,
                     FullMixedData_LS,
                     Data_hwI,
                     Data_hwNS,
                     Data_hwEW,
                     Data_hwFB,
                     dI,
                     dNS,
                     dEW,
                     dFB,
                     Template_Iso_LS,
                     Template_AniNS_LS,
                     Template_AniEW_LS,
                     Template_AniFB_LS,
                     output_log_file,
                     r_gen
                     );
    
    generate_HS_data(
                     Data_Iso_HS,
                     Data_AniNS_HS,
                     Data_AniEW_HS,
                     Data_AniFB_HS,
                     MixedData_NS_EW_HS,
                     MixedData_NS_FB_HS,
                     MixedData_EW_FB_HS,
                     FullMixedData_HS,
                     Data_hwI,
                     Data_hwNS,
                     Data_hwEW,
                     Data_hwFB,
                     dI,
                     dNS,
                     dEW,
                     dFB,
                     Template_Iso_HS,
                     Template_AniNS_HS,
                     Template_AniEW_HS,
                     Template_AniFB_HS,
                     output_log_file,
                     r_gen
                     );
    
     */
     
    //normalize_LS_data(Data_Iso_LS,Data_AniNS_LS,Data_AniEW_LS,Data_AniFB_LS,MixedData_NS_EW_LS,MixedData_NS_FB_LS,MixedData_EW_FB_LS,FullMixedData_LS,log_file);
    //normalize_HS_data(Data_Iso_HS,Data_AniNS_HS,Data_AniEW_HS,Data_AniFB_HS,MixedData_NS_EW_HS,MixedData_NS_FB_HS,MixedData_EW_FB_HS,FullMixedData_HS,log_file);
    
    ////////////////////////////////////// Writing results
    
    //////////////////////////////////// Creating Data out file
    
    TFile data_file(data_path.c_str(),"RECREATE");
    if(data_file.IsZombie()) {
        std::cout << "\n\nError writing ROOT Data TFile. Prorgram finished \n\n";
        log_file << "\n\nError writing ROOT Data TFile. Prorgram finished \n\n";
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
    
    Data_hwI.Write();
    Data_hwNS.Write();
    Data_hwEW.Write();
    Data_hwFB.Write();
    
    data_file.Write();
    data_file.Close();
    
}

void generate_templates(std::ofstream &log_file,std::string &templates_path) {
    
    templates_path = output_path_creator(1);
    
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
    
    create_binning(n_bin_lat_LS,lat_bin_min_LS,lat_bin_max_LS,binning_LS,true);
    create_binning(n_bin_lat_HS,lat_bin_min_HS,lat_bin_max_HS,binning_HS,true);
    
    TH1D Templates_LS[4];
    TH1D Templates_HS[4];
    
    TCanvas LS_Templates_Canvas("LS_Templates_Canvas","Low Statistics Templates");
    TCanvas HS_Templates_Canvas("HS_Templates_Canvas","High Statistics Templates");
    TCanvas TemplateFunctions("TemplateFulctions","Templates Functions");
    
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
    
    
    //////////// Weight Histos
    
    TH1D Template_hwI("template_hWI","Iso Weight",100,0,1);
    TH1D Template_hwNS("template_hWNS","NS Weight",100,-1,1);
    TH1D Template_hwEW("template_hWEW","EW Weight",100,-1,1);
    TH1D Template_hwFB("template_hWFB","FB Weight",100,-1,1);
    
    //////////// TF2
    
    static TF2 dI("dI",IsoMonopole,-TMath::Pi(),TMath::Pi(),-TMath::Pi()/2,TMath::Pi()/2);
    static TF2 dNS("dNS",NSMonopole,-TMath::Pi(),TMath::Pi(),-TMath::Pi()/2,TMath::Pi()/2);
    static TF2 dEW("dEW",EWMonopole,-TMath::Pi(),TMath::Pi(),-TMath::Pi()/2,TMath::Pi()/2);
    static TF2 dFB("dFB",FBMonopole,-TMath::Pi(),TMath::Pi(),-TMath::Pi()/2,TMath::Pi()/2);
    
    ////////////////////////////////////////////////////////////////
    
    generate_LS_templates(
                            Template_Iso_LS,
                            Template_AniNS_LS,
                            Template_AniEW_LS,
                            Template_AniFB_LS,
                            Template_hwI,
                            Template_hwNS,
                            Template_hwEW,
                            Template_hwFB,
                            log_file,
                            TemplateFunctions,
                            dI,
                            dNS,
                            dEW,
                            dFB
                          );
    
    generate_HS_templates(
                            Template_Iso_HS,
                            Template_AniNS_HS,
                            Template_AniEW_HS,
                            Template_AniFB_HS,
                            Template_hwI,
                            Template_hwNS,
                            Template_hwEW,
                            Template_hwFB,
                            log_file,
                            dI,
                            dNS,
                            dEW,
                            dFB
                          );
    
    ///////////////////////// Converting Templates Maps from TH2 to TH1
    
    TH2toTH1_obj(Templates_LS[0],Template_Iso_LS);
    TH2toTH1_obj(Templates_LS[1],Template_AniNS_LS);
    TH2toTH1_obj(Templates_LS[2],Template_AniEW_LS);
    TH2toTH1_obj(Templates_LS[3],Template_AniFB_LS);
    
    TH2toTH1_obj(Templates_HS[0],Template_Iso_HS);
    TH2toTH1_obj(Templates_HS[1],Template_AniNS_HS);
    TH2toTH1_obj(Templates_HS[2],Template_AniEW_HS);
    TH2toTH1_obj(Templates_HS[3],Template_AniFB_HS);
    
    ///////////////////////// Drawing Templates Maps to canvas
    
    for(Int_t idx_comp = 0; idx_comp < 4; idx_comp++) {
        LS_Templates_Canvas.cd(idx_comp+1);
        Templates_LS[idx_comp].Draw();
        HS_Templates_Canvas.cd(idx_comp+1);
        Templates_HS[idx_comp].Draw();
    }
    
    
    //////////////////////////////////// Creating Template out file
    
    TFile template_file(templates_path.c_str(),"RECREATE");
    if(template_file.IsZombie()) {
        std::cout << "\n\nError writing ROOT Templates TFile. Prorgram finished \n\n";
        log_file << "\n\nError writing ROOT Templates TFile. Prorgram finished \n\n";
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
    
    Template_hwI.Write();
    Template_hwNS.Write();
    Template_hwEW.Write();
    Template_hwFB.Write();
    
    LS_Templates_Canvas.Write();
    HS_Templates_Canvas.Write();
    TemplateFunctions.Write();
    
    template_file.Write();
    template_file.Close();
    
    
}

void read_from_file(std::string template_path,std::string data_path,std::ofstream &output_log_file) {
    
    std::string pools_out_path = output_path_creator(4);
    
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
    
    TH1D hPull_Iso_HS("hPull_Iso_HS","Pull Isotropic Sky Map HS",50,-10,10);
    TH1D hPull_AniNS_HS("hPull_AniNS_HS","Pull Anisotropic Sky Map (NS) HS",50,-10,10);
    TH1D hPull_AniEW_HS("hPull_AniEW_HS","Pull Anisotropic Sky Map (EW) HS",50,-10,10);
    TH1D hPull_AniFB_HS("hPull_AniFB_HS","Pull Anisotropic Sky Map (FB) HS",50,-10,10);
    TH1D hPull_Mixed_NS_EW_HS("hPull_Mixed_NS_EW_HS","Pull Anisotropic Sky Map (NS-EW) HS",50,-10,10);
    TH1D hPull_Mixed_NS_FB_HS("hPull_Mixed_NS_FB_HS","Pull Anisotropic Sky Map (NS-FB) HS",50,-10,10);
    TH1D hPull_Mixed_EW_FB_HS("hPull_Mixed_EW_FB_HS","Pull Anisotropic Sky Map (EW-FB) HS",50,-10,10);
    TH1D hPull_FullMixed_HS("hPull_FullMixed_HS","Pull Full Mixed Anisotropic Sky Map HS",50,-10,10);
    
    
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
    
    
    ////////////////////////////////////
    
    TTree fiTree("fiTree","TemplateFit results TTree");
    fitResult tmp_fit;
    
    fiTree.Branch("chi2_LS",tmp_fit.chi2_LS,"chi2_LS[8]/D");
    fiTree.Branch("ndf_LS",tmp_fit.ndf_LS,"ndf_LS[8]/D");
    fiTree.Branch("fit_par_LS",tmp_fit.fit_par_LS,"fit_par_LS[8][4]/D");
    fiTree.Branch("fit_err_LS",tmp_fit.fit_err_LS,"fit_err_LS[8][4]/D");
    fiTree.Branch("CMatrix_Iso_LS",tmp_fit.CMatrix_Iso_LS,"CMatrix_Iso_LS[4][4]/D");
    fiTree.Branch("CMatrix_NS_LS",tmp_fit.CMatrix_NS_LS,"CMatrix_NS_LS[4][4]/D");
    fiTree.Branch("CMatrix_EW_LS",tmp_fit.CMatrix_EW_LS,"CMatrix_EW_LS[4][4]/D");
    fiTree.Branch("CMatrix_FB_LS",tmp_fit.CMatrix_FB_LS,"CMatrix_FB_LS[4][4]/D");
    fiTree.Branch("CMatrix_NS_EW_LS",tmp_fit.CMatrix_NS_EW_LS,"CMatrix_NS_EW_LS[4][4]/D");
    fiTree.Branch("CMatrix_NS_FB_LS",tmp_fit.CMatrix_NS_FB_LS,"CMatrix_NS_FB_LS[4][4]/D");
    fiTree.Branch("CMatrix_EW_FB_LS",tmp_fit.CMatrix_EW_FB_LS,"CMatrix_EW_FB_LS[4][4]/D");
    fiTree.Branch("CMatrix_Full_LS",tmp_fit.CMatrix_Full_LS,"CMatrix_Full_LS[4][4]/D");
    
    fiTree.Branch("theta_binHistos_LS",&tmp_fit.theta_binHistos_LS,"theta_binHistos_LS/I");
    fiTree.Branch("phi_binhistos_LS",&tmp_fit.phi_binHistos_LS,"phi_binHistos_LS/I");
    
    fiTree.Branch("events_LS",&tmp_fit.events_LS,"events_LS/l");
    
    fiTree.Branch("chi2_HS",tmp_fit.chi2_HS,"chi2_HS[8]/D");
    fiTree.Branch("ndf_HS",tmp_fit.ndf_HS,"ndf_HS[8]/D");
    fiTree.Branch("fit_par_HS",tmp_fit.fit_par_HS,"fit_par_HS[8][4]/D");
    fiTree.Branch("fit_err_HS",tmp_fit.fit_err_HS,"fit_err_HS[8][4]/D");
    fiTree.Branch("CMatrix_Iso_HS",tmp_fit.CMatrix_Iso_HS,"CMatrix_Iso_HS[4][4]/D");
    fiTree.Branch("CMatrix_NS_HS",tmp_fit.CMatrix_NS_HS,"CMatrix_NS_HS[4][4]/D");
    fiTree.Branch("CMatrix_EW_HS",tmp_fit.CMatrix_EW_HS,"CMatrix_EW_HS[4][4]/D");
    fiTree.Branch("CMatrix_FB_HS",tmp_fit.CMatrix_FB_HS,"CMatrix_FB_HS[4][4]/D");
    fiTree.Branch("CMatrix_NS_EW_HS",tmp_fit.CMatrix_NS_EW_HS,"CMatrix_NS_EW_HS[4][4]/D");
    fiTree.Branch("CMatrix_NS_FB_HS",tmp_fit.CMatrix_NS_FB_HS,"CMatrix_NS_FB_HS[4][4]/D");
    fiTree.Branch("CMatrix_EW_FB_HS",tmp_fit.CMatrix_EW_FB_HS,"CMatrix_EW_FB_HS[4][4]/D");
    fiTree.Branch("CMatrix_Full_HS",tmp_fit.CMatrix_Full_HS,"CMatrix_Full_HS[4][4]/D");
    
    fiTree.Branch("theta_binHistos_HS",&tmp_fit.theta_binHistos_HS,"theta_binHistos_HS/I");
    fiTree.Branch("phi_binhistos_HS",&tmp_fit.phi_binHistos_HS,"phi_binHistos_HS/I");
    
    fiTree.Branch("events_HS",&tmp_fit.events_HS,"events_HS/l");
    
    
    fiTree.Branch("inputAni",tmp_fit.inputAni,"inputAni[3]/D");
    
    
    /////////////////////////////////////////////////////////////////////////////// Fitting !!!
    
    std::cout << "\n\n //////////////////////////////////////// Low Statistics fits //////////////////////////////////////// \n\n";
    output_log_file << "\n\n //////////////////////////////////////// Low Statistics fits //////////////////////////////////////// \n\n";
    
    std::cout << "\n\n -------------- > Isotropic Monopole \n\n";
    output_log_file << "\n\n -------------- > Isotropic Monopole \n\n";
    
    TemplateFitBH(DataHisto_1D_I_LS,4,Templates_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,0,false);
    getPull(DataHisto_1D_I_LS,Templates_1D_LS,res,hPull_Iso_LS);
    
    std::cout << "\n\n -------------- > Anisotropic Monopole (NS) \n\n";
    output_log_file << "\n\n -------------- > Anisotropic Monopole (NS) \n\n";
    
    TemplateFitBH(DataHisto_1D_AniNS_LS,4,Templates_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,1,false);
    getPull(DataHisto_1D_AniNS_LS,Templates_1D_LS,res,hPull_AniNS_LS);
    
    std::cout << "\n\n -------------- > Anisotropic Monopole (EW) \n\n";
    output_log_file << "\n\n -------------- > Anisotropic Monopole (EW) \n\n";
    
    TemplateFitBH(DataHisto_1D_AniEW_LS,4,Templates_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,2,false);
    getPull(DataHisto_1D_AniEW_LS,Templates_1D_LS,res,hPull_AniEW_LS);
    
    std::cout << "\n\n -------------- > Anisotropic Monopole (FB) \n\n";
    output_log_file << "\n\n -------------- > Anisotropic Monopole (FB) \n\n";
    
    TemplateFitBH(DataHisto_1D_AniFB_LS,4,Templates_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,3,false);
    getPull(DataHisto_1D_AniFB_LS,Templates_1D_LS,res,hPull_AniFB_LS);
    
    std::cout << "\n\n -------------- > NS + EW linear combination \n\n";
    output_log_file << "\n\n -------------- > NS + EW linear combination \n\n";
    
    TemplateFitBH(MixedHisto_1D_NS_EW_LS,4,Templates_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,4,false);
    getPull(MixedHisto_1D_NS_EW_LS,Templates_1D_LS,res,hPull_Mixed_NS_EW_LS);
    
    std::cout << "\n\n -------------- > NS + FB linear combination \n\n";
    output_log_file << "\n\n -------------- > NS + FB linear combination \n\n";
    
    TemplateFitBH(MixedHisto_1D_NS_FB_LS,4,Templates_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,5,false);
    getPull(MixedHisto_1D_NS_FB_LS,Templates_1D_LS,res,hPull_Mixed_NS_FB_LS);
    
    std::cout << "\n\n -------------- > EW + FB linear combination \n\n";
    output_log_file << "\n\n -------------- > EW + FB linear combination \n\n";
    
    TemplateFitBH(MixedHisto_1D_EW_FB_LS,4,Templates_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,6,false);
    getPull(MixedHisto_1D_EW_FB_LS,Templates_1D_LS,res,hPull_Mixed_EW_FB_LS);
    
    std::cout << "\n\n -------------- > FullMixed linear combination \n\n";
    output_log_file << "\n\n -------------- > FullMixed linear combination \n\n";
    
    TemplateFitBH(FullMixedHisto_1D_LS,4,Templates_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,7,false);
    getPull(FullMixedHisto_1D_LS,Templates_1D_LS,res,hPull_FullMixed_LS);
    
    
    std::cout << "\n\n //////////////////////////////////////// High Statistics fits //////////////////////////////////////// \n\n";
    output_log_file << "\n\n //////////////////////////////////////// High Statistics fits //////////////////////////////////////// \n\n";
    
    std::cout << "\n\n -------------- > Isotropic Monopole \n\n";
    output_log_file << "\n\n -------------- > Isotropic Monopole \n\n";
    
    TemplateFitBH(DataHisto_1D_I_HS,4,Templates_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,0,true);
    getPull(DataHisto_1D_I_HS,Templates_1D_HS,res,hPull_Iso_HS);
    
    std::cout << "\n\n -------------- > Anisotropic Monopole (NS) \n\n";
    output_log_file << "\n\n -------------- > Anisotropic Monopole (NS) \n\n";
    
    TemplateFitBH(DataHisto_1D_AniNS_HS,4,Templates_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,1,true);
    getPull(DataHisto_1D_AniNS_HS,Templates_1D_HS,res,hPull_AniNS_HS);
    
    std::cout << "\n\n -------------- > Anisotropic Monopole (EW) \n\n";
    output_log_file << "\n\n -------------- > Anisotropic Monopole (EW) \n\n";
    
    TemplateFitBH(DataHisto_1D_AniEW_HS,4,Templates_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,2,true);
    getPull(DataHisto_1D_AniEW_HS,Templates_1D_HS,res,hPull_AniEW_HS);
    
    std::cout << "\n\n -------------- > Anisotropic Monopole (FB) \n\n";
    output_log_file << "\n\n -------------- > Anisotropic Monopole (FB) \n\n";
    
    TemplateFitBH(DataHisto_1D_AniFB_HS,4,Templates_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,3,true);
    getPull(DataHisto_1D_AniFB_HS,Templates_1D_HS,res,hPull_AniFB_HS);
    
    std::cout << "\n\n -------------- > NS + EW linear combination \n\n";
    output_log_file << "\n\n -------------- > NS + EW linear combination \n\n";
    
    TemplateFitBH(MixedHisto_1D_NS_EW_HS,4,Templates_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,4,true);
    getPull(MixedHisto_1D_NS_EW_HS,Templates_1D_HS,res,hPull_Mixed_NS_EW_HS);
    
    std::cout << "\n\n -------------- > NS + FB linear combination \n\n";
    output_log_file << "\n\n -------------- > NS + FB linear combination \n\n";
    
    TemplateFitBH(MixedHisto_1D_NS_FB_HS,4,Templates_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,5,true);
    getPull(MixedHisto_1D_NS_FB_HS,Templates_1D_HS,res,hPull_Mixed_NS_FB_HS);
    
    std::cout << "\n\n -------------- > EW + FB linear combination \n\n";
    output_log_file << "\n\n -------------- > EW + FB linear combination \n\n";
    
    TemplateFitBH(MixedHisto_1D_EW_FB_HS,4,Templates_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,6,true);
    getPull(MixedHisto_1D_EW_FB_HS,Templates_1D_HS,res,hPull_Mixed_EW_FB_HS);
    
    std::cout << "\n\n -------------- > FullMixed linear combination \n\n";
    output_log_file << "\n\n -------------- > FullMixed linear combination \n\n";
    
    TemplateFitBH(FullMixedHisto_1D_HS,4,Templates_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,7,true);
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

#else

void generate_and_fit(std::ofstream &output_log_file)
{
    
    ////////////////////////////////////////////////////////// Variables declaration ///////////////////////////////////////////////////
    
    std::string template_out_path = output_path_creator(1);
    std::string DAMPE_template_out_path = output_path_creator(2);
    std::string DAMPE_data_out_path = output_path_creator(8);
    std::string data_out_path = output_path_creator(3);
    std::string pools_out_path = output_path_creator(4);
    std::string projections_out_path = output_path_creator(5);
    std::string tree_out_path = output_path_creator(7);
    
    //////////////////////////////////////////////////////////// Histos declaration ////////////////////////////////////////////////////
    
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
    
    create_binning(n_bin_lat_LS,lat_bin_min_LS,lat_bin_max_LS,binning_LS,true);
    create_binning(n_bin_lat_HS,lat_bin_min_HS,lat_bin_max_HS,binning_HS,true);
    
    
    
    //////////// Canvases

    TCanvas TemplateFunctions("TemplateFulctions","Templates Functions");
    
    TCanvas LS_Templates_Canvas("LS_Templates_Canvas","Low Statistics Templates");
    TCanvas HS_Templates_Canvas("HS_Templates_Canvas","High Statistics Templates");
    
    LS_Templates_Canvas.Divide(2,2);
    HS_Templates_Canvas.Divide(2,2);
    
    //////////// All Sky Histos
    
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
    
    //////////// Weight Histos
    
    TH1D Template_hwI("template_hWI","Iso Weight",100,0,1);
    TH1D Template_hwNS("template_hWNS","NS Weight",100,-1,1);
    TH1D Template_hwEW("template_hWEW","EW Weight",100,-1,1);
    TH1D Template_hwFB("template_hWFB","FB Weight",100,-1,1);
    
    
    ///////////////////////////////////////// DAMPE's Histos
    
    TH2D DAMPE_ReferenceMap;
    
    /////////// Templates
    
    TH2D DAMPE_Template_Iso_LS;
    TH2D DAMPE_Template_AniNS_LS;
    TH2D DAMPE_Template_AniEW_LS;
    TH2D DAMPE_Template_AniFB_LS;
    
    TH2D DAMPE_Template_Iso_HS;
    TH2D DAMPE_Template_AniNS_HS;
    TH2D DAMPE_Template_AniEW_HS;
    TH2D DAMPE_Template_AniFB_HS;
    
    //////////////////////////////////////////////////////// TF2s declaration //////////////////////////////////////////////////////////
    
    static TF2 dI("dI",IsoMonopole,-TMath::Pi(),TMath::Pi(),-TMath::Pi()/2,TMath::Pi()/2);
    static TF2 dNS("dNS",NSMonopole,-TMath::Pi(),TMath::Pi(),-TMath::Pi()/2,TMath::Pi()/2);
    static TF2 dEW("dEW",EWMonopole,-TMath::Pi(),TMath::Pi(),-TMath::Pi()/2,TMath::Pi()/2);
    static TF2 dFB("dFB",FBMonopole,-TMath::Pi(),TMath::Pi(),-TMath::Pi()/2,TMath::Pi()/2);
    
    //////////// Simulation data parameters
    
    create_simulation_seeds();
    
    std::string tmp_seed_str;
    UInt_t tmp_seed;
    std::ifstream inSeed(seeds_path);
    if(!inSeed.is_open()) {
        std::cout << "\n\nCannot open input seeds file! Program finished !" << std::endl;
        exit(100);
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    //////////// Generating All Sky templates
    
    generate_LS_templates(
                            Template_Iso_LS,
                            Template_AniNS_LS,
                            Template_AniEW_LS,
                            Template_AniFB_LS,
                            Template_hwI,
                            Template_hwNS,
                            Template_hwEW,
                            Template_hwFB,
                            output_log_file,
                            TemplateFunctions,
                            dI,
                            dNS,
                            dEW,
                            dFB
                          );
    
    generate_HS_templates(
                            Template_Iso_HS,
                            Template_AniNS_HS,
                            Template_AniEW_HS,
                            Template_AniFB_HS,
                            Template_hwI,
                            Template_hwNS,
                            Template_hwEW,
                            Template_hwFB,
                            output_log_file,
                            dI,
                            dNS,
                            dEW,
                            dFB
                          );
    
    //////////// Generating DAMPE templates
    
    read_DAMPE_FullIso(DAMPE_ReferenceMap);
    
    get_DAMPE_templates(
                            DAMPE_Template_Iso_LS,
                            DAMPE_Template_AniNS_LS,
                            DAMPE_Template_AniEW_LS,
                            DAMPE_Template_AniFB_LS,
                            DAMPE_Template_Iso_HS,
                            DAMPE_Template_AniNS_HS,
                            DAMPE_Template_AniEW_HS,
                            DAMPE_Template_AniFB_HS,
                            DAMPE_ReferenceMap,
                            Template_Iso_LS,
                            Template_AniNS_LS,
                            Template_AniEW_LS,
                            Template_AniFB_LS,
                            Template_Iso_HS,
                            Template_AniNS_HS,
                            Template_AniEW_HS,
                            Template_AniFB_HS
                        );
    
    
    
    //////////// Generating TTree and linking class variables
    
    TTree fiTree("fiTree","TemplateFit results TTree");
    fitResult tmp_fit;
    
    fiTree.Branch("chi2_LS",tmp_fit.chi2_LS,"chi2_LS[8]/D");
    fiTree.Branch("ndf_LS",tmp_fit.ndf_LS,"ndf_LS[8]/D");
    fiTree.Branch("fit_par_LS",tmp_fit.fit_par_LS,"fit_par_LS[8][4]/D");
    fiTree.Branch("fit_err_LS",tmp_fit.fit_err_LS,"fit_err_LS[8][4]/D");
    fiTree.Branch("delta_LS",tmp_fit.delta_LS,"delta_LS[8]/D");
    fiTree.Branch("sum_par_LS",tmp_fit.sum_par_LS,"sum_par_LS[8]/D");
    fiTree.Branch("CMatrix_Iso_LS",tmp_fit.CMatrix_Iso_LS,"CMatrix_Iso_LS[4][4]/D");
    fiTree.Branch("CMatrix_NS_LS",tmp_fit.CMatrix_NS_LS,"CMatrix_NS_LS[4][4]/D");
    fiTree.Branch("CMatrix_EW_LS",tmp_fit.CMatrix_EW_LS,"CMatrix_EW_LS[4][4]/D");
    fiTree.Branch("CMatrix_FB_LS",tmp_fit.CMatrix_FB_LS,"CMatrix_FB_LS[4][4]/D");
    fiTree.Branch("CMatrix_NS_EW_LS",tmp_fit.CMatrix_NS_EW_LS,"CMatrix_NS_EW_LS[4][4]/D");
    fiTree.Branch("CMatrix_NS_FB_LS",tmp_fit.CMatrix_NS_FB_LS,"CMatrix_NS_FB_LS[4][4]/D");
    fiTree.Branch("CMatrix_EW_FB_LS",tmp_fit.CMatrix_EW_FB_LS,"CMatrix_EW_FB_LS[4][4]/D");
    fiTree.Branch("CMatrix_Full_LS",tmp_fit.CMatrix_Full_LS,"CMatrix_Full_LS[4][4]/D");
    
    fiTree.Branch("theta_binHistos_LS",&tmp_fit.theta_binHistos_LS,"theta_binHistos_LS/I");
    fiTree.Branch("phi_binhistos_LS",&tmp_fit.phi_binHistos_LS,"phi_binHistos_LS/I");
    
    fiTree.Branch("events_LS",&tmp_fit.events_LS,"events_LS/l");
    
    fiTree.Branch("chi2_HS",tmp_fit.chi2_HS,"chi2_HS[8]/D");
    fiTree.Branch("ndf_HS",tmp_fit.ndf_HS,"ndf_HS[8]/D");
    fiTree.Branch("fit_par_HS",tmp_fit.fit_par_HS,"fit_par_HS[8][4]/D");
    fiTree.Branch("fit_err_HS",tmp_fit.fit_err_HS,"fit_err_HS[8][4]/D");
    fiTree.Branch("delta_HS",tmp_fit.delta_HS,"delta_HS[8]/D");
    fiTree.Branch("sum_par_HS",tmp_fit.sum_par_HS,"sum_par_HS[8]/D");
    fiTree.Branch("CMatrix_Iso_HS",tmp_fit.CMatrix_Iso_HS,"CMatrix_Iso_HS[4][4]/D");
    fiTree.Branch("CMatrix_NS_HS",tmp_fit.CMatrix_NS_HS,"CMatrix_NS_HS[4][4]/D");
    fiTree.Branch("CMatrix_EW_HS",tmp_fit.CMatrix_EW_HS,"CMatrix_EW_HS[4][4]/D");
    fiTree.Branch("CMatrix_FB_HS",tmp_fit.CMatrix_FB_HS,"CMatrix_FB_HS[4][4]/D");
    fiTree.Branch("CMatrix_NS_EW_HS",tmp_fit.CMatrix_NS_EW_HS,"CMatrix_NS_EW_HS[4][4]/D");
    fiTree.Branch("CMatrix_NS_FB_HS",tmp_fit.CMatrix_NS_FB_HS,"CMatrix_NS_FB_HS[4][4]/D");
    fiTree.Branch("CMatrix_EW_FB_HS",tmp_fit.CMatrix_EW_FB_HS,"CMatrix_EW_FB_HS[4][4]/D");
    fiTree.Branch("CMatrix_Full_HS",tmp_fit.CMatrix_Full_HS,"CMatrix_Full_HS[4][4]/D");
    
    fiTree.Branch("theta_binHistos_HS",&tmp_fit.theta_binHistos_HS,"theta_binHistos_HS/I");
    fiTree.Branch("phi_binhistos_HS",&tmp_fit.phi_binHistos_HS,"phi_binHistos_HS/I");
    
    fiTree.Branch("events_HS",&tmp_fit.events_HS,"events_HS/l");
    
    fiTree.Branch("inputAni",tmp_fit.inputAni,"inputAni[3]/D");
    
    
    
    for(Int_t itry = 0; itry < ntry; ++itry)
    {
        
        std::cout << "\n\n ////////////////////////////// Try " << itry << " ////////////////////////////// \n\n";
        output_log_file << "\n\n ////////////////////////////// Try " << itry << " ////////////////////////////// \n\n";
        
        for(Int_t idx_ani = 0; idx_ani < ani_values; ++idx_ani)
        {
            std::getline(inSeed,tmp_seed_str);
            tmp_seed = (UInt_t)std::stoul(tmp_seed_str,nullptr,10);
            
            /*
            allSky_singleTry_fit(
                                    output_log_file,
                                    template_out_path,
                                    data_out_path,
                                    pools_out_path,
                                    projections_out_path,
                                    Template_Iso_LS,
                                    Template_AniNS_LS,
                                    Template_AniEW_LS,
                                    Template_AniFB_LS,
                                    Template_Iso_HS,
                                    Template_AniNS_HS,
                                    Template_AniEW_HS,
                                    Template_AniFB_HS,
                                    Template_hwI,
                                    Template_hwNS,
                                    Template_hwEW,
                                    Template_hwFB,
                                    dI,
                                    dNS,
                                    dEW,
                                    dFB,
                                    fiTree,
                                    tmp_fit,
                                    NS_anisotropy[idx_ani],
                                    EW_anisotropy[idx_ani],
                                    FB_anisotropy[idx_ani],
                                    LS_Templates_Canvas,
                                    HS_Templates_Canvas,
                                    tmp_seed
                                 );
            
            
             DAMPE_singleTry_fit(
                                    output_log_file,
                                    data_out_path,
                                    pools_out_path,
                                    DAMPE_Template_Iso_LS,
                                    DAMPE_Template_AniNS_LS,
                                    DAMPE_Template_AniEW_LS,
                                    DAMPE_Template_AniFB_LS,
                                    DAMPE_Template_Iso_HS,
                                    DAMPE_Template_AniNS_HS,
                                    DAMPE_Template_AniEW_HS,
                                    DAMPE_Template_AniFB_HS,
                                    NS_anisotropy[idx_ani],
                                    EW_anisotropy[idx_ani],
                                    FB_anisotropy[idx_ani],
                                    tmp_fit,
                                    tmp_seed,
                                    fiTree
                                );
            
            */
            
            DAMPE_relative_singleTry_fit(
                                            output_log_file,
                                            data_out_path,
                                            pools_out_path,
                                            Template_Iso_LS,
                                            Template_AniNS_LS,
                                            Template_AniEW_LS,
                                            Template_AniFB_LS,
                                            Template_Iso_HS,
                                            Template_AniNS_HS,
                                            Template_AniEW_HS,
                                            Template_AniFB_HS,
                                            DAMPE_ReferenceMap,
                                            NS_anisotropy[idx_ani],
                                            EW_anisotropy[idx_ani],
                                            FB_anisotropy[idx_ani],
                                            tmp_fit,
                                            tmp_seed,
                                            fiTree
                                         );
            
        }
        
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    //////////////////////////////////// Creating Template out file
    
    TFile template_file(template_out_path.c_str(),"RECREATE");
    if(template_file.IsZombie()) {
        std::cout << "\n\nError writing ROOT Templates TFile. Prorgram finished \n\n";
        output_log_file << "\n\nError writing ROOT Templates TFile. Prorgram finished \n\n";
        exit(100);
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
    
    Template_hwI.Write();
    Template_hwNS.Write();
    Template_hwEW.Write();
    Template_hwFB.Write();
    
    LS_Templates_Canvas.Write();
    HS_Templates_Canvas.Write();
    TemplateFunctions.Write();
    
    template_file.Write();
    template_file.Close();
    
    //////////////////////////////////// Creating DAMPE Template out file
    
    TFile DAMPE_template_file(DAMPE_template_out_path.c_str(),"RECREATE");
    if(DAMPE_template_file.IsZombie()) {
        std::cout << "\n\nError writing ROOT DAMPE Templates TFile. Prorgram finished \n\n";
        output_log_file << "\n\nError writing ROOT DAMPE Templates TFile. Prorgram finished \n\n";
        exit(100);
    }
    
    //////////////////////////////////// Writing Templates
    
    DAMPE_Template_Iso_LS.Write();
    DAMPE_Template_AniNS_LS.Write();
    DAMPE_Template_AniEW_LS.Write();
    DAMPE_Template_AniFB_LS.Write();
    
    DAMPE_Template_Iso_HS.Write();
    DAMPE_Template_AniNS_HS.Write();
    DAMPE_Template_AniEW_HS.Write();
    DAMPE_Template_AniFB_HS.Write();
    
    DAMPE_template_file.Write();
    DAMPE_template_file.Close();
    
    
    //////////////////////////////////// Creating TTree out file
    
    TFile tree_file(tree_out_path.c_str(),"RECREATE");
    if(template_file.IsZombie()) {
        std::cout << "\n\nError writing ROOT TTree TFile. Prorgram finished \n\n";
        output_log_file << "\n\nError writing ROOT TTree TFile. Prorgram finished \n\n";
        exit(100);
    }
    
    fiTree.Write();
    
    tree_file.Write();
    tree_file.Close();
    
}



void allSky_singleTry_fit(
                            std::ofstream &output_log_file,
                            std::string &template_out_path,
                            std::string &data_out_path,
                            std::string &pools_out_path,
                            std::string &projections_out_path,
                            TH2D &Template_Iso_LS,
                            TH2D &Template_AniNS_LS,
                            TH2D &Template_AniEW_LS,
                            TH2D &Template_AniFB_LS,
                            TH2D &Template_Iso_HS,
                            TH2D &Template_AniNS_HS,
                            TH2D &Template_AniEW_HS,
                            TH2D &Template_AniFB_HS,
                            TH1D &Template_hwI,
                            TH1D &Template_hwNS,
                            TH1D &Template_hwEW,
                            TH1D &Template_hwFB,
                            TF2 &dI,
                            TF2 &dNS,
                            TF2 &dEW,
                            TF2 &dFB,
                            TTree &fiTree,
                            fitResult &tmp_fit,
                            Double_t NS_anisotropy,
                            Double_t EW_anisotropy,
                            Double_t FB_anisotropy,
                            TCanvas &LS_Templates_Canvas,
                            TCanvas &HS_Templates_Canvas,
                            UInt_t tmp_seed
                          )
{
    
    TRandom3 r_gen(tmp_seed);
    
    tmp_fit.inputAni[0] = NS_anisotropy;
    tmp_fit.inputAni[1] = EW_anisotropy;
    tmp_fit.inputAni[2] = FB_anisotropy;
    
    //////////// TemplateFit variables
    
    Double_t res[4],res_err[4];
    Double_t initialValues[4]={1,0,0,0};
    
    double fullFitResults_HS[4][8];
    
    std::vector<TH1D*> DataProjections;
    DataProjections.resize(16);

    std::vector<TH1D*> TemplatesProjections;
    TemplatesProjections.resize(8);
    
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
    
    //////////// Cleaned Low Statistics data
    
    TH2D CData_AniNS_LS;
    TH2D CData_AniEW_LS;
    TH2D CData_AniFB_LS;
    
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
    
    //////////// Cleaned High Statistics data
    
    TH2D CData_AniNS_HS;
    TH2D CData_AniEW_HS;
    TH2D CData_AniFB_HS;
    
    //////////// Weight Histos
    
    TH1D Data_hwI("data_hWI","Iso Weight",100,0,1);
    TH1D Data_hwNS("data_hWNS","NS Weight",100,-2,2);
    TH1D Data_hwEW("data_hWEW","EW Weight",100,-2,2);
    TH1D Data_hwFB("data_hWFB","FB Weight",100,-2,2);
    
    
    //////////// Pool's histos
    
    //////////// Low statistics pools
    
    TH1D hPull_Iso_LS("hPull_Iso_LS","Pull Isotropic Sky Map LS",100,-5,5);
    TH1D hPull_AniNS_LS("hPull_AniNS_LS","Pull Anisotropic Sky Map (NS) LS",100,-5,5);
    TH1D hPull_AniEW_LS("hPull_AniEW_LS","Pull Anisotropic Sky Map (EW) LS",100,-5,5);
    TH1D hPull_AniFB_LS("hPull_AniFB_LS","Pull Anisotropic Sky Map (FB) LS",100,-5,5);
    TH1D hPull_Mixed_NS_EW_LS("hPull_Mixed_NS_EW_LS","Pull Anisotropic Sky Map (NS-EW) LS",100,-5,5);
    TH1D hPull_Mixed_NS_FB_LS("hPull_Mixed_NS_FB_LS","Pull Anisotropic Sky Map (NS-FB) LS",100,-5,5);
    TH1D hPull_Mixed_EW_FB_LS("hPull_Mixed_EW_FB_LS","Pull Anisotropic Sky Map (EW-FB) LS",100,-5,5);
    TH1D hPull_FullMixed_LS("hPull_FullMixed_LS","Pull Full Mixed Anisotropic Sky Map LS",100,-5,5);
    
    //////////// High statistics pools
    
    TH1D hPull_Iso_HS("hPull_Iso_HS","Pull Isotropic Sky Map HS",50,-10,10);
    TH1D hPull_AniNS_HS("hPull_AniNS_HS","Pull Anisotropic Sky Map (NS) HS",50,-10,10);
    TH1D hPull_AniEW_HS("hPull_AniEW_HS","Pull Anisotropic Sky Map (EW) HS",50,-10,10);
    TH1D hPull_AniFB_HS("hPull_AniFB_HS","Pull Anisotropic Sky Map (FB) HS",50,-10,10);
    TH1D hPull_Mixed_NS_EW_HS("hPull_Mixed_NS_EW_HS","Pull Anisotropic Sky Map (NS-EW) HS",50,-10,10);
    TH1D hPull_Mixed_NS_FB_HS("hPull_Mixed_NS_FB_HS","Pull Anisotropic Sky Map (NS-FB) HS",50,-10,10);
    TH1D hPull_Mixed_EW_FB_HS("hPull_Mixed_EW_FB_HS","Pull Anisotropic Sky Map (EW-FB) HS",50,-10,10);
    TH1D hPull_FullMixed_HS("hPull_FullMixed_HS","Pull Full Mixed Anisotropic Sky Map HS",50,-10,10);
    
    
    //////////// 1D Histos
    
    TH1D Templates_LS[4];
    TH1D Templates_HS[4];
    
    TH1D DataHisto_I_LS;
    TH1D DataHisto_NS_LS;
    TH1D DataHisto_EW_LS;
    TH1D DataHisto_FB_LS;
    TH1D CDataHisto_NS_LS;
    TH1D CDataHisto_EW_LS;
    TH1D CDataHisto_FB_LS;
    TH1D MixedDataHisto_NS_EW_LS;
    TH1D MixedDataHisto_NS_FB_LS;
    TH1D MixedDataHisto_EW_FB_LS;
    TH1D FullMixedDataHisto_LS;
    
    TH1D DataHisto_I_HS;
    TH1D DataHisto_NS_HS;
    TH1D DataHisto_EW_HS;
    TH1D DataHisto_FB_HS;
    TH1D CDataHisto_NS_HS;
    TH1D CDataHisto_EW_HS;
    TH1D CDataHisto_FB_HS;
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
    
    TH1D* CDataHisto_1D_AniNS_LS = nullptr;
    TH1D* CDataHisto_1D_AniEW_LS = nullptr;
    TH1D* CDataHisto_1D_AniFB_LS = nullptr;
    
    TH1D* CDataHisto_1D_AniNS_HS = nullptr;
    TH1D* CDataHisto_1D_AniEW_HS = nullptr;
    TH1D* CDataHisto_1D_AniFB_HS = nullptr;
    
    
    TH1D* Template_Iso_HS_ProjX = nullptr;
    TH1D* Template_Iso_HS_ProjY = nullptr;
    
    TH1D* Template_AniNS_HS_ProjX = nullptr;
    TH1D* Template_AniNS_HS_ProjY = nullptr;
    
    TH1D* Template_AniEW_HS_ProjX = nullptr;
    TH1D* Template_AniEW_HS_ProjY = nullptr;
    
    TH1D* Template_AniFB_HS_ProjX = nullptr;
    TH1D* Template_AniFB_HS_ProjY = nullptr;
    
    
    TH1D* Data_Iso_HS_ProjX = nullptr;
    TH1D* Data_Iso_HS_ProjY = nullptr;
    
    TH1D* Data_AniNS_HS_ProjX = nullptr;
    TH1D* Data_AniNS_HS_ProjY = nullptr;
    
    TH1D* Data_AniEW_HS_ProjX = nullptr;
    TH1D* Data_AniEW_HS_ProjY = nullptr;
    
    TH1D* Data_AniFB_HS_ProjX = nullptr;
    TH1D* Data_AniFB_HS_ProjY = nullptr;
    
    TH1D* MixedData_NS_EW_HS_ProjX = nullptr;
    TH1D* MixedData_NS_EW_HS_ProjY = nullptr;
    
    TH1D* MixedData_NS_FB_HS_ProjX = nullptr;
    TH1D* MixedData_NS_FB_HS_ProjY = nullptr;
    
    TH1D* MixedData_EW_FB_HS_ProjX = nullptr;
    TH1D* MixedData_EW_FB_HS_ProjY = nullptr;
    
    TH1D* FullMixedData_HS_ProjX = nullptr;
    TH1D* FullMixedData_HS_ProjY = nullptr;
    
    
    tmp_fit.theta_binHistos_LS = Data_Iso_LS->GetNbinsY();
    tmp_fit.phi_binHistos_LS = Data_Iso_LS->GetNbinsX();
    tmp_fit.theta_binHistos_HS = Data_Iso_HS->GetNbinsY();
    tmp_fit.phi_binHistos_HS = Data_Iso_HS->GetNbinsX();
    
    tmp_fit.events_LS = data_all_sky_LS_events;
    tmp_fit.events_HS = data_all_sky_HS_events;
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
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
                        Data_hwI,
                        Data_hwNS,
                        Data_hwEW,
                        Data_hwFB,
                        dI,
                        dNS,
                        dEW,
                        dFB,
                        Template_Iso_LS,
                        Template_AniNS_LS,
                        Template_AniEW_LS,
                        Template_AniFB_LS,
                        output_log_file,
                        r_gen
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
                        Data_hwI,
                        Data_hwNS,
                        Data_hwEW,
                        Data_hwFB,
                        dI,
                        dNS,
                        dEW,
                        dFB,
                        Template_Iso_HS,
                        Template_AniNS_HS,
                        Template_AniEW_HS,
                        Template_AniFB_HS,
                        output_log_file,
                        r_gen
                     );
    
    
    /*
     
     MC_generate_LS_data(
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
                            Data_hwI,
                            Data_hwNS,
                            Data_hwEW,
                            Data_hwFB,
                            dI,
                            dNS,
                            dEW,
                            dFB,
                            output_log_file,
                            r_gen
                        );
     
     MC_generate_HS_data(
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
                            Data_hwI,
                            Data_hwNS,
                            Data_hwEW,
                            Data_hwFB,
                            dI,
                            dNS,
                            dEW,
                            dFB,
                            output_log_file,
                            r_gen
                        );
     
    */
    
    
    get_cleaned_LS_histos(
                            Data_AniNS_LS,
                            Data_AniEW_LS,
                            Data_AniFB_LS,
                            CData_AniNS_LS,
                            CData_AniEW_LS,
                            CData_AniFB_LS,
                            Data_Iso_LS
                          );
    
    get_cleaned_HS_histos(
                            Data_AniNS_HS,
                            Data_AniEW_HS,
                            Data_AniFB_HS,
                            CData_AniNS_HS,
                            CData_AniEW_HS,
                            CData_AniFB_HS,
                            Data_Iso_HS
                          );
    
    //normalize_LS_templates(Template_Iso_LS,Template_AniNS_LS,Template_AniEW_LS,Template_AniFB_LS,output_log_file);
    //normalize_HS_templates(Template_Iso_HS,Template_AniNS_HS,Template_AniEW_HS,Template_AniFB_HS,output_log_file);
    
    //normalize_LS_data(Data_Iso_LS,Data_AniNS_LS,Data_AniEW_LS,Data_AniFB_LS,MixedData_NS_EW_LS,MixedData_NS_FB_LS,MixedData_EW_FB_LS,FullMixedData_LS,output_log_file);
    //normalize_HS_data(Data_Iso_HS,Data_AniNS_HS,Data_AniEW_HS,Data_AniFB_HS,MixedData_NS_EW_HS,MixedData_NS_FB_HS,MixedData_EW_FB_HS,FullMixedData_HS,output_log_file);
    
    
    ///////////////////////// Obtaining projections
    
    
    //////////// Templates
    
    Template_Iso_HS_ProjX = (TH1D*) Template_Iso_HS.ProjectionX();
    Template_Iso_HS_ProjX->SetName("Template_Iso_HS_ProjX");
    Template_Iso_HS_ProjX->SetTitle("X Projection Iso Template HS");
    
    Template_Iso_HS_ProjY = (TH1D*) Template_Iso_HS.ProjectionY();
    Template_Iso_HS_ProjY->SetName("Template_Iso_HS_ProjY");
    Template_Iso_HS_ProjY->SetTitle("Y Projection Iso Template HS");
    
    Template_AniNS_HS_ProjX = (TH1D*) Template_AniNS_HS.ProjectionX();
    Template_AniNS_HS_ProjX->SetName("Template_AniNS_HS_ProjX");
    Template_AniNS_HS_ProjX->SetTitle("X Projection NS Template HS");
    
    Template_AniNS_HS_ProjY = (TH1D*) Template_AniNS_HS.ProjectionY();
    Template_AniNS_HS_ProjY->SetName("Template_AniNS_HS_ProjY");
    Template_AniNS_HS_ProjY->SetTitle("Y Projection NS Template HS");
    
    Template_AniEW_HS_ProjX = (TH1D*) Template_AniEW_HS.ProjectionX();
    Template_AniEW_HS_ProjX->SetName("Template_AniEW_HS_ProjX");
    Template_AniEW_HS_ProjX->SetTitle("X Projection EW Template HS");
    
    Template_AniEW_HS_ProjY = (TH1D*) Template_AniEW_HS.ProjectionY();
    Template_AniEW_HS_ProjY->SetName("Template_AniEW_HS_ProjY");
    Template_AniEW_HS_ProjY->SetTitle("Y Projection EW Template HS");
    
    Template_AniFB_HS_ProjX = (TH1D*) Template_AniFB_HS.ProjectionX();
    Template_AniFB_HS_ProjX->SetName("Template_AniFB_HS_ProjX");
    Template_AniFB_HS_ProjX->SetTitle("X Projection FB Template HS");
    
    Template_AniFB_HS_ProjY = (TH1D*) Template_AniFB_HS.ProjectionY();
    Template_AniFB_HS_ProjY->SetName("Template_AniFB_HS_ProjY");
    Template_AniFB_HS_ProjY->SetTitle("Y Projection FB Template HS");

    
    TemplatesProjections[0] = Template_Iso_HS_ProjX;
    TemplatesProjections[1] = Template_Iso_HS_ProjY;
    TemplatesProjections[2] = Template_AniNS_HS_ProjX;
    TemplatesProjections[3] = Template_AniNS_HS_ProjY;
    TemplatesProjections[4] = Template_AniEW_HS_ProjX;
    TemplatesProjections[5] = Template_AniEW_HS_ProjY;
    TemplatesProjections[6] = Template_AniFB_HS_ProjX;
    TemplatesProjections[7] = Template_AniFB_HS_ProjY;
    
    //////////// Data
    
    Data_Iso_HS_ProjX = (TH1D*) Data_Iso_HS->ProjectionX();
    Data_Iso_HS_ProjX->SetName("Data_Iso_HS_ProjX");
    Data_Iso_HS_ProjX->SetTitle("X Projection Iso Data");
    
    Data_Iso_HS_ProjY = (TH1D*) Data_Iso_HS->ProjectionY();
    Data_Iso_HS_ProjY->SetName("Data_Iso_HS_ProjY");
    Data_Iso_HS_ProjY->SetTitle("Y Projection Iso Data");
    
    Data_AniNS_HS_ProjX = (TH1D*) Data_AniNS_HS->ProjectionX();
    Data_AniNS_HS_ProjX->SetName("Data_AniNS_HS_ProjX");
    Data_AniNS_HS_ProjX->SetTitle("X Projection NS Data");
    
    Data_AniNS_HS_ProjY = (TH1D*) Data_AniNS_HS->ProjectionY();
    Data_AniNS_HS_ProjY->SetName("Data_AniNS_HS_ProjY");
    Data_AniNS_HS_ProjY->SetTitle("Y Projection NS Data");
    
    Data_AniEW_HS_ProjX = (TH1D*) Data_AniEW_HS->ProjectionX();
    Data_AniEW_HS_ProjX->SetName("Data_AniEW_HS_ProjX");
    Data_AniEW_HS_ProjX->SetTitle("X Projection EW Data");
    
    Data_AniEW_HS_ProjY = (TH1D*) Data_AniEW_HS->ProjectionY();
    Data_AniEW_HS_ProjY->SetName("Data_AniEW_HS_ProjY");
    Data_AniEW_HS_ProjY->SetTitle("Y Projection EW Data");
    
    Data_AniFB_HS_ProjX = (TH1D*) Data_AniFB_HS->ProjectionX();
    Data_AniFB_HS_ProjX->SetName("Data_AniFB_HS_ProjX");
    Data_AniFB_HS_ProjX->SetTitle("X Projection FB Data");
    
    Data_AniFB_HS_ProjY = (TH1D*) Data_AniFB_HS->ProjectionY();
    Data_AniFB_HS_ProjY->SetName("Data_AniFB_HS_ProjY");
    Data_AniFB_HS_ProjY->SetTitle("Y Projection FB Data");
    
    MixedData_NS_EW_HS_ProjX = (TH1D*) MixedData_NS_EW_HS->ProjectionX();
    MixedData_NS_EW_HS_ProjX->SetName("MixedData_NS_EW_HS_ProjX");
    MixedData_NS_EW_HS_ProjX->SetTitle("X Projection Mixed NS-EW Data");
    
    MixedData_NS_EW_HS_ProjY = (TH1D*) MixedData_NS_EW_HS->ProjectionY();
    MixedData_NS_EW_HS_ProjY->SetName("MixedData_NS_EW_HS_ProjY");
    MixedData_NS_EW_HS_ProjY->SetTitle("Y Projection Mixed NS-EW Data");
    
    MixedData_NS_FB_HS_ProjX = (TH1D*) MixedData_NS_FB_HS->ProjectionX();
    MixedData_NS_FB_HS_ProjX->SetName("MixedData_NS_FB_HS_ProjX");
    MixedData_NS_FB_HS_ProjX->SetTitle("X Projection Mixed NS-FB Data");
    
    MixedData_NS_FB_HS_ProjY = (TH1D*) MixedData_NS_FB_HS->ProjectionY();
    MixedData_NS_FB_HS_ProjY->SetName("MixedData_NS_FB_HS_ProjY");
    MixedData_NS_FB_HS_ProjY->SetTitle("Y Projection Mixed NS-FB Data");
    
    MixedData_EW_FB_HS_ProjX = (TH1D*) MixedData_EW_FB_HS->ProjectionX();
    MixedData_EW_FB_HS_ProjX->SetName("MixedData_EW_FB_HS_ProjX");
    MixedData_EW_FB_HS_ProjX->SetTitle("X Projection Mixed EW-FB Data");
    
    MixedData_EW_FB_HS_ProjY = (TH1D*) MixedData_EW_FB_HS->ProjectionY();
    MixedData_EW_FB_HS_ProjY->SetName("MixedData_EW_FB_HS_ProjY");
    MixedData_EW_FB_HS_ProjY->SetTitle("Y Projection Mixed EW-FB Data");
    
    FullMixedData_HS_ProjX = (TH1D*) FullMixedData_HS->ProjectionX();
    FullMixedData_HS_ProjX->SetName("FullMixedData_HS_ProjX");
    FullMixedData_HS_ProjX->SetTitle("X Projection Full Mixed Data");
    
    FullMixedData_HS_ProjY = (TH1D*) FullMixedData_HS->ProjectionY();
    FullMixedData_HS_ProjY->SetName("FullMixedData_HS_ProjY");
    FullMixedData_HS_ProjY->SetTitle("Y Projection Full Mixed Data");
    
    
    DataProjections[0] = Data_Iso_HS_ProjX;
    DataProjections[1] = Data_Iso_HS_ProjY;
    DataProjections[2] = Data_AniNS_HS_ProjX;
    DataProjections[3] = Data_AniNS_HS_ProjY;
    DataProjections[4] = Data_AniEW_HS_ProjX;
    DataProjections[5] = Data_AniEW_HS_ProjY;
    DataProjections[6] = Data_AniFB_HS_ProjX;
    DataProjections[7] = Data_AniFB_HS_ProjY;
    DataProjections[8] = MixedData_NS_EW_HS_ProjX;
    DataProjections[9] = MixedData_NS_EW_HS_ProjY;
    DataProjections[10] = MixedData_NS_FB_HS_ProjX;
    DataProjections[11] = MixedData_NS_FB_HS_ProjY;
    DataProjections[12] = MixedData_EW_FB_HS_ProjX;
    DataProjections[13] = MixedData_EW_FB_HS_ProjY;
    DataProjections[14] = FullMixedData_HS_ProjX;
    DataProjections[15] = FullMixedData_HS_ProjY;
    
    
    ///////////////////////// Converting Templates Maps from TH2 to TH1
    
    TH2toTH1_obj(Templates_LS[0],Template_Iso_LS);
    TH2toTH1_obj(Templates_LS[1],Template_AniNS_LS);
    TH2toTH1_obj(Templates_LS[2],Template_AniEW_LS);
    TH2toTH1_obj(Templates_LS[3],Template_AniFB_LS);
    
    TH2toTH1_obj(Templates_HS[0],Template_Iso_HS);
    TH2toTH1_obj(Templates_HS[1],Template_AniNS_HS);
    TH2toTH1_obj(Templates_HS[2],Template_AniEW_HS);
    TH2toTH1_obj(Templates_HS[3],Template_AniFB_HS);
    
    ///////////////////////// Converting Data Maps from TH2 to TH1
    
    TH2toTH1_ptr(DataHisto_I_LS,Data_Iso_LS);
    TH2toTH1_ptr(DataHisto_NS_LS,Data_AniNS_LS);
    TH2toTH1_ptr(DataHisto_EW_LS,Data_AniEW_LS);
    TH2toTH1_ptr(DataHisto_FB_LS,Data_AniFB_LS);
    
    TH2toTH1_obj(CDataHisto_NS_LS,CData_AniNS_LS);
    TH2toTH1_obj(CDataHisto_EW_LS,CData_AniEW_LS);
    TH2toTH1_obj(CDataHisto_FB_LS,CData_AniFB_LS);
    
    TH2toTH1_ptr(MixedDataHisto_NS_EW_LS,MixedData_NS_EW_LS);
    TH2toTH1_ptr(MixedDataHisto_NS_FB_LS,MixedData_NS_FB_LS);
    TH2toTH1_ptr(MixedDataHisto_EW_FB_LS,MixedData_EW_FB_LS);
    TH2toTH1_ptr(FullMixedDataHisto_LS,FullMixedData_LS);
    
    TH2toTH1_ptr(DataHisto_I_HS,Data_Iso_HS);
    TH2toTH1_ptr(DataHisto_NS_HS,Data_AniNS_HS);
    TH2toTH1_ptr(DataHisto_EW_HS,Data_AniEW_HS);
    TH2toTH1_ptr(DataHisto_FB_HS,Data_AniFB_HS);
    
    TH2toTH1_ptr(CDataHisto_NS_HS,Data_AniNS_HS);
    TH2toTH1_ptr(CDataHisto_EW_HS,Data_AniEW_HS);
    TH2toTH1_ptr(CDataHisto_FB_HS,Data_AniFB_HS);
    
    TH2toTH1_ptr(MixedDataHisto_NS_EW_HS,MixedData_NS_EW_HS);
    TH2toTH1_ptr(MixedDataHisto_NS_FB_HS,MixedData_NS_FB_HS);
    TH2toTH1_ptr(MixedDataHisto_EW_FB_HS,MixedData_EW_FB_HS);
    TH2toTH1_ptr(FullMixedDataHisto_HS,FullMixedData_HS);
    
    
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
    
    CDataHisto_1D_AniNS_LS = &CDataHisto_NS_LS;
    CDataHisto_1D_AniEW_LS = &CDataHisto_EW_LS;
    CDataHisto_1D_AniFB_LS = &CDataHisto_FB_LS;
    
    MixedHisto_1D_NS_EW_LS = &MixedDataHisto_NS_EW_LS;
    MixedHisto_1D_NS_FB_LS = &MixedDataHisto_NS_FB_LS;
    MixedHisto_1D_EW_FB_LS = &MixedDataHisto_EW_FB_LS;
    FullMixedHisto_1D_LS = &FullMixedDataHisto_LS;
    
    DataHisto_1D_I_HS = &DataHisto_I_HS;
    DataHisto_1D_AniNS_HS = &DataHisto_NS_HS;
    DataHisto_1D_AniEW_HS = &DataHisto_EW_HS;
    DataHisto_1D_AniFB_HS = &DataHisto_FB_HS;
    
    CDataHisto_1D_AniNS_HS = &CDataHisto_NS_HS;
    CDataHisto_1D_AniEW_HS = &CDataHisto_EW_HS;
    CDataHisto_1D_AniFB_HS = &CDataHisto_FB_HS;
    
    MixedHisto_1D_NS_EW_HS = &MixedDataHisto_NS_EW_HS;
    MixedHisto_1D_NS_FB_HS = &MixedDataHisto_NS_FB_HS;
    MixedHisto_1D_EW_FB_HS = &MixedDataHisto_EW_FB_HS;
    FullMixedHisto_1D_HS = &FullMixedDataHisto_HS;
    
    /////////////////////////////////////////////////////////////////////////////// Fitting !!!
    
    std::cout << "\n\n //////////////////////////////////////// Low Statistics fits //////////////////////////////////////// \n\n";
    output_log_file << "\n\n //////////////////////////////////////// Low Statistics fits //////////////////////////////////////// \n\n";
    
    std::cout << "\n\n -------------- > Isotropic Monopole \n\n";
    output_log_file << "\n\n -------------- > Isotropic Monopole \n\n";
    
    TemplateFitBH(DataHisto_1D_I_LS,4,Templates_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,0,false);
    getPull(DataHisto_1D_I_LS,Templates_1D_LS,res,hPull_Iso_LS);
    
    std::cout << "\n\n -------------- > Anisotropic Monopole (NS) \n\n";
    output_log_file << "\n\n -------------- > Anisotropic Monopole (NS) \n\n";
    
    TemplateFitBH(DataHisto_1D_AniNS_LS,4,Templates_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,1,false);
    getPull(DataHisto_1D_AniNS_LS,Templates_1D_LS,res,hPull_AniNS_LS);
    
    std::cout << "\n\n -------------- > Anisotropic Monopole (EW) \n\n";
    output_log_file << "\n\n -------------- > Anisotropic Monopole (EW) \n\n";
    
    TemplateFitBH(DataHisto_1D_AniEW_LS,4,Templates_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,2,false);
    getPull(DataHisto_1D_AniEW_LS,Templates_1D_LS,res,hPull_AniEW_LS);
    
    std::cout << "\n\n -------------- > Anisotropic Monopole (FB) \n\n";
    output_log_file << "\n\n -------------- > Anisotropic Monopole (FB) \n\n";
    
    TemplateFitBH(DataHisto_1D_AniFB_LS,4,Templates_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,3,false);
    getPull(DataHisto_1D_AniFB_LS,Templates_1D_LS,res,hPull_AniFB_LS);
    
    std::cout << "\n\n -------------- > NS + EW linear combination \n\n";
    output_log_file << "\n\n -------------- > NS + EW linear combination \n\n";
    
    TemplateFitBH(MixedHisto_1D_NS_EW_LS,4,Templates_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,4,false);
    getPull(MixedHisto_1D_NS_EW_LS,Templates_1D_LS,res,hPull_Mixed_NS_EW_LS);
    
    std::cout << "\n\n -------------- > NS + FB linear combination \n\n";
    output_log_file << "\n\n -------------- > NS + FB linear combination \n\n";
    
    TemplateFitBH(MixedHisto_1D_NS_FB_LS,4,Templates_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,5,false);
    getPull(MixedHisto_1D_NS_FB_LS,Templates_1D_LS,res,hPull_Mixed_NS_FB_LS);
    
    std::cout << "\n\n -------------- > EW + FB linear combination \n\n";
    output_log_file << "\n\n -------------- > EW + FB linear combination \n\n";
    
    TemplateFitBH(MixedHisto_1D_EW_FB_LS,4,Templates_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,6,false);
    getPull(MixedHisto_1D_EW_FB_LS,Templates_1D_LS,res,hPull_Mixed_EW_FB_LS);
    
    std::cout << "\n\n -------------- > FullMixed linear combination \n\n";
    output_log_file << "\n\n -------------- > FullMixed linear combination \n\n";
    
    TemplateFitBH(FullMixedHisto_1D_LS,4,Templates_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,7,false);
    getPull(FullMixedHisto_1D_LS,Templates_1D_LS,res,hPull_FullMixed_LS);
    
    
    std::cout << "\n\n //////////////////////////////////////// High Statistics fits //////////////////////////////////////// \n\n";
    output_log_file << "\n\n //////////////////////////////////////// High Statistics fits //////////////////////////////////////// \n\n";
    
    std::cout << "\n\n -------------- > Isotropic Monopole \n\n";
    output_log_file << "\n\n -------------- > Isotropic Monopole \n\n";
    
    TemplateFitBH(DataHisto_1D_I_HS,4,Templates_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,0,true);
    getPull(DataHisto_1D_I_HS,Templates_1D_HS,res,hPull_Iso_HS);
    
    for(Int_t idxr=0; idxr < 4; ++idxr)
        fullFitResults_HS[idxr][0] = res[idxr];
    
    
    std::cout << "\n\n -------------- > Anisotropic Monopole (NS) \n\n";
    output_log_file << "\n\n -------------- > Anisotropic Monopole (NS) \n\n";
    
    TemplateFitBH(DataHisto_1D_AniNS_HS,4,Templates_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,1,true);
    getPull(DataHisto_1D_AniNS_HS,Templates_1D_HS,res,hPull_AniNS_HS);
    
    for(Int_t idxr=0; idxr < 4; ++idxr)
        fullFitResults_HS[idxr][1] = res[idxr];
    
    std::cout << "\n\n -------------- > Anisotropic Monopole (EW) \n\n";
    output_log_file << "\n\n -------------- > Anisotropic Monopole (EW) \n\n";
    
    TemplateFitBH(DataHisto_1D_AniEW_HS,4,Templates_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,2,true);
    getPull(DataHisto_1D_AniEW_HS,Templates_1D_HS,res,hPull_AniEW_HS);
    
    for(Int_t idxr=0; idxr < 4; ++idxr)
        fullFitResults_HS[idxr][2] = res[idxr];
    
    std::cout << "\n\n -------------- > Anisotropic Monopole (FB) \n\n";
    output_log_file << "\n\n -------------- > Anisotropic Monopole (FB) \n\n";
    
    TemplateFitBH(DataHisto_1D_AniFB_HS,4,Templates_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,3,true);
    getPull(DataHisto_1D_AniFB_HS,Templates_1D_HS,res,hPull_AniFB_HS);
    
    for(Int_t idxr=0; idxr < 4; ++idxr)
        fullFitResults_HS[idxr][3] = res[idxr];
    
    std::cout << "\n\n -------------- > NS + EW linear combination \n\n";
    output_log_file << "\n\n -------------- > NS + EW linear combination \n\n";
    
    TemplateFitBH(MixedHisto_1D_NS_EW_HS,4,Templates_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,4,true);
    getPull(MixedHisto_1D_NS_EW_HS,Templates_1D_HS,res,hPull_Mixed_NS_EW_HS);
    
    for(Int_t idxr=0; idxr < 4; ++idxr)
        fullFitResults_HS[idxr][4] = res[idxr];
    
    std::cout << "\n\n -------------- > NS + FB linear combination \n\n";
    output_log_file << "\n\n -------------- > NS + FB linear combination \n\n";
    
    TemplateFitBH(MixedHisto_1D_NS_FB_HS,4,Templates_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,5,true);
    getPull(MixedHisto_1D_NS_FB_HS,Templates_1D_HS,res,hPull_Mixed_NS_FB_HS);
    
    for(Int_t idxr=0; idxr < 4; ++idxr)
        fullFitResults_HS[idxr][5] = res[idxr];
    
    std::cout << "\n\n -------------- > EW + FB linear combination \n\n";
    output_log_file << "\n\n -------------- > EW + FB linear combination \n\n";
    
    TemplateFitBH(MixedHisto_1D_EW_FB_HS,4,Templates_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,6,true);
    getPull(MixedHisto_1D_EW_FB_HS,Templates_1D_HS,res,hPull_Mixed_EW_FB_HS);
    
    for(Int_t idxr=0; idxr < 4; ++idxr)
        fullFitResults_HS[idxr][6] = res[idxr];
    
    std::cout << "\n\n -------------- > FullMixed linear combination \n\n";
    output_log_file << "\n\n -------------- > FullMixed linear combination \n\n";
    
    TemplateFitBH(FullMixedHisto_1D_HS,4,Templates_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,7,true);
    getPull(FullMixedHisto_1D_HS,Templates_1D_HS,res,hPull_FullMixed_HS);
    
    for(Int_t idxr=0; idxr < 4; ++idxr)
        fullFitResults_HS[idxr][7] = res[idxr];
    
    
    components_analysis(
                            TemplatesProjections,
                            DataProjections,
                            fullFitResults_HS,
                            output_log_file
                        );
    
    /*
    
    std::cout << "\n\n //////////////////////////////////////// Low Statistics Cleaned fits //////////////////////////////////////// \n\n";
    output_log_file << "\n\n //////////////////////////////////////// Low Statistics Cleaned fits //////////////////////////////////////// \n\n";
    
    std::cout << "\n\n -------------- > Anisotropic Monopole (NS) \n\n";
    output_log_file << "\n\n -------------- > Anisotropic Monopole (NS) \n\n";
    
    TemplateFitBH(CDataHisto_1D_AniNS_LS,4,Templates_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file);
    //getPull(DataHisto_1D_AniNS_LS,Templates_1D_LS,res,hPull_AniNS_LS);
    
    std::cout << "\n\n -------------- > Anisotropic Monopole (EW) \n\n";
    output_log_file << "\n\n -------------- > Anisotropic Monopole (EW) \n\n";
    
    TemplateFitBH(CDataHisto_1D_AniEW_LS,4,Templates_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file);
    //getPull(DataHisto_1D_AniEW_LS,Templates_1D_LS,res,hPull_AniEW_LS);
    
    std::cout << "\n\n -------------- > Anisotropic Monopole (FB) \n\n";
    output_log_file << "\n\n -------------- > Anisotropic Monopole (FB) \n\n";
    
    TemplateFitBH(CDataHisto_1D_AniFB_LS,4,Templates_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file);
    //getPull(DataHisto_1D_AniFB_LS,Templates_1D_LS,res,hPull_AniFB_LS);
    
    
    
    
    std::cout << "\n\n //////////////////////////////////////// High Statistics Cleaned fits //////////////////////////////////////// \n\n";
    output_log_file << "\n\n //////////////////////////////////////// High Cleaned Statistics fits //////////////////////////////////////// \n\n";
    
    std::cout << "\n\n -------------- > Anisotropic Monopole (NS) \n\n";
    output_log_file << "\n\n -------------- > Anisotropic Monopole (NS) \n\n";
    
    TemplateFitBH(CDataHisto_1D_AniNS_HS,4,Templates_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file);
    //getPull(DataHisto_1D_AniNS_HS,Templates_1D_HS,res,hPull_AniNS_HS);
    
    std::cout << "\n\n -------------- > Anisotropic Monopole (EW) \n\n";
    output_log_file << "\n\n -------------- > Anisotropic Monopole (EW) \n\n";
    
    TemplateFitBH(CDataHisto_1D_AniEW_HS,4,Templates_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file);
    //getPull(DataHisto_1D_AniEW_HS,Templates_1D_HS,res,hPull_AniEW_HS);
    
    std::cout << "\n\n -------------- > Anisotropic Monopole (FB) \n\n";
    output_log_file << "\n\n -------------- > Anisotropic Monopole (FB) \n\n";
    
    TemplateFitBH(CDataHisto_1D_AniFB_HS,4,Templates_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file);
    //getPull(DataHisto_1D_AniFB_HS,Templates_1D_HS,res,hPull_AniFB_HS);
    
    */
    
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    std::cout<<"\n\nSimulation Completed !\n\n";
    output_log_file << "\n\nSimulation Completed !\n\n";
    
    
    if(write_tmp_histos)
    {
    
        //////////////////////////////////// Creating Data out file
    
        TFile data_file(data_out_path.c_str(),"RECREATE");
        if(data_file.IsZombie()) {
            std::cout << "\n\nError writing ROOT Data TFile. Prorgram finished \n\n";
            output_log_file << "\n\nError writing ROOT Data TFile. Prorgram finished \n\n";
            exit(-2);
        }
    
    
        //////////////////////////////////// Writing Data
    
        Data_Iso_LS->Write();
        Data_AniNS_LS->Write();
        Data_AniEW_LS->Write();
        Data_AniFB_LS->Write();
    
        CData_AniNS_LS.Write();
        CData_AniEW_LS.Write();
        CData_AniFB_LS.Write();
    
        MixedData_NS_EW_LS->Write();
        MixedData_NS_FB_LS->Write();
        MixedData_EW_FB_LS->Write();
        FullMixedData_LS->Write();
    
        Data_Iso_HS->Write();
        Data_AniNS_HS->Write();
        Data_AniEW_HS->Write();
        Data_AniFB_HS->Write();
    
        CData_AniNS_HS.Write();
        CData_AniEW_HS.Write();
        CData_AniFB_HS.Write();
    
        MixedData_NS_EW_HS->Write();
        MixedData_NS_FB_HS->Write();
        MixedData_EW_FB_HS->Write();
        FullMixedData_HS->Write();
    
        Data_hwI.Write();
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
            exit(100);
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
    
    
        //////////////////////////////////// Creating Projections out file
        
        TFile projection_file(projections_out_path.c_str(),"RECREATE");
        if(projection_file.IsZombie()) {
            std::cout << "\n\nError writing ROOT Projection TFile. Prorgram finished \n\n";
            output_log_file << "\n\nError writing ROOT Projection TFile. Prorgram finished \n\n";
            exit(100);
        }
    
        //////////////////////////////////// Writing Projections
    
        Template_Iso_HS_ProjX->Write();
        Template_Iso_HS_ProjY->Write();
    
        Template_AniNS_HS_ProjX->Write();
        Template_AniNS_HS_ProjY->Write();
    
        Template_AniEW_HS_ProjX->Write();
        Template_AniEW_HS_ProjY->Write();
    
        Template_AniFB_HS_ProjX->Write();
        Template_AniFB_HS_ProjY->Write();
    
        Data_Iso_HS_ProjX->Write();
        Data_Iso_HS_ProjY->Write();
    
        Data_AniNS_HS_ProjX->Write();
        Data_AniNS_HS_ProjY->Write();
    
        Data_AniEW_HS_ProjX->Write();
        Data_AniEW_HS_ProjY->Write();
    
        Data_AniFB_HS_ProjX->Write();
        Data_AniFB_HS_ProjY->Write();
    
        MixedData_NS_EW_HS_ProjX->Write();
        MixedData_NS_EW_HS_ProjY->Write();
    
        MixedData_NS_FB_HS_ProjX->Write();
        MixedData_NS_FB_HS_ProjY->Write();
    
        MixedData_EW_FB_HS_ProjX->Write();
        MixedData_EW_FB_HS_ProjY->Write();
    
        FullMixedData_HS_ProjX->Write();
        FullMixedData_HS_ProjY->Write();
    
        projection_file.Close();
    
    }

    fiTree.Fill();
    
}


void DAMPE_singleTry_fit(
                            std::ofstream &output_log_file,
                            std::string &data_out_path,
                            std::string &pools_out_path,
                            TH2D &DAMPE_Template_Iso_LS,
                            TH2D &DAMPE_Template_AniNS_LS,
                            TH2D &DAMPE_Template_AniEW_LS,
                            TH2D &DAMPE_Template_AniFB_LS,
                            TH2D &DAMPE_Template_Iso_HS,
                            TH2D &DAMPE_Template_AniNS_HS,
                            TH2D &DAMPE_Template_AniEW_HS,
                            TH2D &DAMPE_Template_AniFB_HS,
                            Double_t NS_anisotropy,
                            Double_t EW_anisotropy,
                            Double_t FB_anisotropy,
                            fitResult &tmp_fit,
                            UInt_t tmp_seed,
                            TTree &fiTree
                         )

{
    TRandom3 r_gen(tmp_seed);
    
    tmp_fit.inputAni[0] = NS_anisotropy;
    tmp_fit.inputAni[1] = EW_anisotropy;
    tmp_fit.inputAni[2] = FB_anisotropy;
    
    //////////// TemplateFit variables
    
    Double_t res[4],res_err[4];
    Double_t initialValues[4]={1,0,0,0};
    
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
    
    TH1D DAMPE_hPull_Iso_LS("DAMPE_hPull_Iso_LS","Pull Isotropic DAMPE Sky Map LS",100,-5,5);
    TH1D DAMPE_hPull_AniNS_LS("DAMPE_hPull_AniNS_LS","Pull Anisotropic DAMPE Sky Map (NS) LS",100,-5,5);
    TH1D DAMPE_hPull_AniEW_LS("DAMPE_hPull_AniEW_LS","Pull Anisotropic DAMPE Sky Map (EW) LS",100,-5,5);
    TH1D DAMPE_hPull_AniFB_LS("DAMPE_hPull_AniFB_LS","Pull Anisotropic DAMPE Sky Map (FB) LS",100,-5,5);
    TH1D DAMPE_hPull_Mixed_NS_EW_LS("DAMPE_hPull_Mixed_NS_EW_LS","Pull Anisotropic DAMPE Sky Map (NS-EW) LS",100,-5,5);
    TH1D DAMPE_hPull_Mixed_NS_FB_LS("DAMPE_hPull_Mixed_NS_FB_LS","Pull Anisotropic DAMPE Sky Map (NS-FB) LS",100,-5,5);
    TH1D DAMPE_hPull_Mixed_EW_FB_LS("DAMPE_hPull_Mixed_EW_FB_LS","Pull Anisotropic DAMPE Sky Map (EW-FB) LS",100,-5,5);
    TH1D DAMPE_hPull_FullMixed_LS("DAMPE_hPull_FullMixed_LS","Pull Full Mixed Anisotropic DAMPE Sky Map LS",100,-5,5);
    
    //////////// High statistics pools
    
    TH1D DAMPE_hPull_Iso_HS("DAMPE_hPull_Iso_HS","Pull Isotropic DAMPE Sky Map HS",50,-10,10);
    TH1D DAMPE_hPull_AniNS_HS("DAMPE_hPull_AniNS_HS","Pull Anisotropic DAMPE Sky Map (NS) HS",50,-10,10);
    TH1D DAMPE_hPull_AniEW_HS("DAMPE_hPull_AniEW_HS","Pull Anisotropic DAMPE Sky Map (EW) HS",50,-10,10);
    TH1D DAMPE_hPull_AniFB_HS("DAMPE_hPull_AniFB_HS","Pull Anisotropic DAMPE Sky Map (FB) HS",50,-10,10);
    TH1D DAMPE_hPull_Mixed_NS_EW_HS("DAMPE_hPull_Mixed_NS_EW_HS","Pull Anisotropic DAMPE Sky Map (NS-EW) HS",50,-10,10);
    TH1D DAMPE_hPull_Mixed_NS_FB_HS("DAMPE_hPull_Mixed_NS_FB_HS","Pull Anisotropic DAMPE Sky Map (NS-FB) HS",50,-10,10);
    TH1D DAMPE_hPull_Mixed_EW_FB_HS("DAMPE_hPull_Mixed_EW_FB_HS","Pull Anisotropic DAMPE Sky Map (EW-FB) HS",50,-10,10);
    TH1D DAMPE_hPull_FullMixed_HS("DAMPE_hPull_FullMixed_HS","Pull Full Mixed Anisotropic DAMPE Sky Map HS",50,-10,10);
    
    
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
                            r_gen
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
                            r_gen
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
    
    
    if(write_tmp_histos)
    {
        
        //////////////////////////////////// Creating Data out file
        
        TFile DAMPE_data_file(data_out_path.c_str(),"RECREATE");
        if(DAMPE_data_file.IsZombie()) {
            std::cout << "\n\nError writing ROOT DAMPE Data TFile. Prorgram finished \n\n";
            output_log_file << "\n\nError writing ROOT DAMPE Data TFile. Prorgram finished \n\n";
            exit(-2);
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
        
        //////////////////////////////////// Gaussian Fit HS hPools
        
        DAMPE_hPull_Iso_HS.Fit("DAMPE_gaus","L,E,M");
        DAMPE_hPull_AniNS_HS.Fit("DAMPE_gaus","L,E,M");
        DAMPE_hPull_AniEW_HS.Fit("DAMPE_gaus","L,E,M");
        DAMPE_hPull_AniFB_HS.Fit("DAMPE_gaus","L,E,M");
        DAMPE_hPull_Mixed_NS_EW_HS.Fit("DAMPE_gaus","L,E,M");
        DAMPE_hPull_Mixed_NS_FB_HS.Fit("DAMPE_gaus","L,E,M");
        DAMPE_hPull_Mixed_EW_FB_HS.Fit("DAMPE_gaus","L,E,M");
        DAMPE_hPull_FullMixed_HS.Fit("DAMPE_gaus","L,E,M");
        
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


void DAMPE_relative_singleTry_fit(
                                    std::ofstream &output_log_file,
                                    std::string &data_out_path,
                                    std::string &pools_out_path,
                                    TH2D &Template_Iso_LS,
                                    TH2D &Template_AniNS_LS,
                                    TH2D &Template_AniEW_LS,
                                    TH2D &Template_AniFB_LS,
                                    TH2D &Template_Iso_HS,
                                    TH2D &Template_AniNS_HS,
                                    TH2D &Template_AniEW_HS,
                                    TH2D &Template_AniFB_HS,
                                    TH2D &DAMPE_ReferenceMap,
                                    Double_t NS_anisotropy,
                                    Double_t EW_anisotropy,
                                    Double_t FB_anisotropy,
                                    fitResult &tmp_fit,
                                    UInt_t tmp_seed,
                                    TTree &fiTree
                                  )

{
    TRandom3 r_gen(tmp_seed);
    
    tmp_fit.inputAni[0] = NS_anisotropy;
    tmp_fit.inputAni[1] = EW_anisotropy;
    tmp_fit.inputAni[2] = FB_anisotropy;
    
    //////////// TemplateFit variables
    
    Double_t res[4],res_err[4];
    Double_t initialValues[4]={1,0,0,0};
    
    //////////// Low statistics reference data
    
    TH2D* pDAMPE_ReferenceMap_LS = (TH2D*) DAMPE_ReferenceMap.Clone("pDAMPE_ReferenceMap_LS");
    TH2D* pDAMPE_ReferenceMap_HS = (TH2D*) DAMPE_ReferenceMap.Clone("pDAMPE_ReferenceMap_HS");
    
    pDAMPE_ReferenceMap_LS->Reset();
    pDAMPE_ReferenceMap_HS->Reset();
    
    scale_reference_map(DAMPE_ReferenceMap,pDAMPE_ReferenceMap_LS,true);
    scale_reference_map(DAMPE_ReferenceMap,pDAMPE_ReferenceMap_HS,false);
    
    pDAMPE_ReferenceMap_LS->Sumw2();
    pDAMPE_ReferenceMap_HS->Sumw2();
    
    //////////// Low statistics data
    
    TH2D* DAMPE_Data_Iso_LS = (TH2D*) Template_Iso_LS.Clone("DAMPE_Data_Iso_LS");
    TH2D* DAMPE_Data_AniNS_LS = (TH2D*) Template_Iso_LS.Clone("DAMPE_Data_AniNS_LS");
    TH2D* DAMPE_Data_AniEW_LS = (TH2D*) Template_Iso_LS.Clone("DAMPE_Data_AniEW_LS");
    TH2D* DAMPE_Data_AniFB_LS = (TH2D*) Template_Iso_LS.Clone("DAMPE_Data_AniFB_LS");
    TH2D* DAMPE_MixedData_NS_EW_LS = (TH2D*) Template_Iso_LS.Clone("DAMPE_MixedData_NS_EW_LS");
    TH2D* DAMPE_MixedData_NS_FB_LS = (TH2D*) Template_Iso_LS.Clone("DAMPE_MixedData_NS_FB_LS");
    TH2D* DAMPE_MixedData_EW_FB_LS = (TH2D*) Template_Iso_LS.Clone("DAMPE_MixedData_EW_FB_LS");
    TH2D* DAMPE_FullMixedData_LS = (TH2D*) Template_Iso_LS.Clone("DAMPE_FullMixedData_LS");
    
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
    
    DAMPE_Data_Iso_LS->Sumw2();
    DAMPE_Data_AniNS_LS->Sumw2();
    DAMPE_Data_AniEW_LS->Sumw2();
    DAMPE_Data_AniFB_LS->Sumw2();
    DAMPE_MixedData_NS_EW_LS->Sumw2();
    DAMPE_MixedData_NS_FB_LS->Sumw2();
    DAMPE_MixedData_EW_FB_LS->Sumw2();
    DAMPE_FullMixedData_LS->Sumw2();
    
    TH2D relative_DAMPE_Data_Iso_LS;
    TH2D relative_DAMPE_Data_AniNS_LS;
    TH2D relative_DAMPE_Data_AniEW_LS;
    TH2D relative_DAMPE_Data_AniFB_LS;
    TH2D relative_DAMPE_MixedData_NS_EW_LS;
    TH2D relative_DAMPE_MixedData_NS_FB_LS;
    TH2D relative_DAMPE_MixedData_EW_FB_LS;
    TH2D relative_DAMPE_FullMixedData_LS;
    
    //////////// High statistics data
    
    TH2D* DAMPE_Data_Iso_HS = (TH2D*) Template_Iso_HS.Clone("DAMPE_Data_Iso_HS");
    TH2D* DAMPE_Data_AniNS_HS = (TH2D*) Template_Iso_HS.Clone("DAMPE_Data_AniNS_HS");
    TH2D* DAMPE_Data_AniEW_HS = (TH2D*) Template_Iso_HS.Clone("DAMPE_Data_AniEW_HS");
    TH2D* DAMPE_Data_AniFB_HS = (TH2D*) Template_Iso_HS.Clone("DAMPE_Data_AniFB_HS");
    TH2D* DAMPE_MixedData_NS_EW_HS = (TH2D*) Template_Iso_HS.Clone("DAMPE_MixedData_NS_EW_HS");
    TH2D* DAMPE_MixedData_NS_FB_HS = (TH2D*) Template_Iso_HS.Clone("DAMPE_MixedData_NS_FB_HS");
    TH2D* DAMPE_MixedData_EW_FB_HS = (TH2D*) Template_Iso_HS.Clone("DAMPE_MixedData_EW_FB_HS");
    TH2D* DAMPE_FullMixedData_HS = (TH2D*) Template_Iso_HS.Clone("DAMPE_FullMixedData_HS");
    
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
    
    DAMPE_Data_Iso_HS->Sumw2();
    DAMPE_Data_AniNS_HS->Sumw2();
    DAMPE_Data_AniEW_HS->Sumw2();
    DAMPE_Data_AniFB_HS->Sumw2();
    DAMPE_MixedData_NS_EW_HS->Sumw2();
    DAMPE_MixedData_NS_FB_HS->Sumw2();
    DAMPE_MixedData_EW_FB_HS->Sumw2();
    DAMPE_FullMixedData_HS->Sumw2();
    
    TH2D relative_DAMPE_Data_Iso_HS;
    TH2D relative_DAMPE_Data_AniNS_HS;
    TH2D relative_DAMPE_Data_AniEW_HS;
    TH2D relative_DAMPE_Data_AniFB_HS;
    TH2D relative_DAMPE_MixedData_NS_EW_HS;
    TH2D relative_DAMPE_MixedData_NS_FB_HS;
    TH2D relative_DAMPE_MixedData_EW_FB_HS;
    TH2D relative_DAMPE_FullMixedData_HS;
    
    //////////// Pool's histos
    
    //////////// Low statistics pools
    
    TH1D relative_DAMPE_hPull_Iso_LS("relative_DAMPE_hPull_Iso_LS","Pull Isotropic relative DAMPE Sky Map LS",100,-5,5);
    TH1D relative_DAMPE_hPull_AniNS_LS("relative_DAMPE_hPull_AniNS_LS","Pull Anisotropic relative DAMPE Sky Map (NS) LS",100,-5,5);
    TH1D relative_DAMPE_hPull_AniEW_LS("relative_DAMPE_hPull_AniEW_LS","Pull Anisotropic relative DAMPE Sky Map (EW) LS",100,-5,5);
    TH1D relative_DAMPE_hPull_AniFB_LS("relative_DAMPE_hPull_AniFB_LS","Pull Anisotropic relative DAMPE Sky Map (FB) LS",100,-5,5);
    TH1D relative_DAMPE_hPull_Mixed_NS_EW_LS("relative_DAMPE_hPull_Mixed_NS_EW_LS","Pull Anisotropic relative DAMPE Sky Map (NS-EW) LS",100,-5,5);
    TH1D relative_DAMPE_hPull_Mixed_NS_FB_LS("relative_DAMPE_hPull_Mixed_NS_FB_LS","Pull Anisotropic relative DAMPE Sky Map (NS-FB) LS",100,-5,5);
    TH1D relative_DAMPE_hPull_Mixed_EW_FB_LS("relative_DAMPE_hPull_Mixed_EW_FB_LS","Pull Anisotropic relative DAMPE Sky Map (EW-FB) LS",100,-5,5);
    TH1D relative_DAMPE_hPull_FullMixed_LS("relative_DAMPE_hPull_FullMixed_LS","Pull Full Mixed Anisotropic relative DAMPE Sky Map LS",100,-5,5);
    
    //////////// High statistics pools
    
    TH1D relative_DAMPE_hPull_Iso_HS("relative_DAMPE_hPull_Iso_HS","Pull Isotropic relative DAMPE Sky Map HS",50,-10,10);
    TH1D relative_DAMPE_hPull_AniNS_HS("relative_DAMPE_hPull_AniNS_HS","Pull Anisotropic relative DAMPE Sky Map (NS) HS",50,-10,10);
    TH1D relative_DAMPE_hPull_AniEW_HS("relative_DAMPE_hPull_AniEW_HS","Pull Anisotropic relative DAMPE Sky Map (EW) HS",50,-10,10);
    TH1D relative_DAMPE_hPull_AniFB_HS("relative_DAMPE_hPull_AniFB_HS","Pull Anisotropic relative DAMPE Sky Map (FB) HS",50,-10,10);
    TH1D relative_DAMPE_hPull_Mixed_NS_EW_HS("relative_DAMPE_hPull_Mixed_NS_EW_HS","Pull Anisotropic relative DAMPE Sky Map (NS-EW) HS",50,-10,10);
    TH1D relative_DAMPE_hPull_Mixed_NS_FB_HS("relative_DAMPE_hPull_Mixed_NS_FB_HS","Pull Anisotropic relative DAMPE Sky Map (NS-FB) HS",50,-10,10);
    TH1D relative_DAMPE_hPull_Mixed_EW_FB_HS("relative_DAMPE_hPull_Mixed_EW_FB_HS","Pull Anisotropic relative DAMPE Sky Map (EW-FB) HS",50,-10,10);
    TH1D relative_DAMPE_hPull_FullMixed_HS("relative_DAMPE_hPull_FullMixed_HS","Pull Full Mixed Anisotropic relative DAMPE Sky Map HS",50,-10,10);
    
    
    //////////// 1D Histos
    
    TH1D Templates_LS[4];
    TH1D Templates_HS[4];
    
    TH1D relative_DAMPE_DataHisto_I_LS;
    TH1D relative_DAMPE_DataHisto_NS_LS;
    TH1D relative_DAMPE_DataHisto_EW_LS;
    TH1D relative_DAMPE_DataHisto_FB_LS;
    TH1D relative_DAMPE_MixedDataHisto_NS_EW_LS;
    TH1D relative_DAMPE_MixedDataHisto_NS_FB_LS;
    TH1D relative_DAMPE_MixedDataHisto_EW_FB_LS;
    TH1D relative_DAMPE_FullMixedDataHisto_LS;
    
    TH1D relative_DAMPE_DataHisto_I_HS;
    TH1D relative_DAMPE_DataHisto_NS_HS;
    TH1D relative_DAMPE_DataHisto_EW_HS;
    TH1D relative_DAMPE_DataHisto_FB_HS;
    TH1D relative_DAMPE_MixedDataHisto_NS_EW_HS;
    TH1D relative_DAMPE_MixedDataHisto_NS_FB_HS;
    TH1D relative_DAMPE_MixedDataHisto_EW_FB_HS;
    TH1D relative_DAMPE_FullMixedDataHisto_HS;
    
    TH1* Templates_1D_LS[4] = {nullptr,nullptr,nullptr,nullptr};
    TH1* Templates_1D_HS[4] = {nullptr,nullptr,nullptr,nullptr};
    
    TH1D* relative_DAMPE_DataHisto_1D_I_LS = nullptr;
    TH1D* relative_DAMPE_DataHisto_1D_AniNS_LS = nullptr;
    TH1D* relative_DAMPE_DataHisto_1D_AniEW_LS = nullptr;
    TH1D* relative_DAMPE_DataHisto_1D_AniFB_LS = nullptr;
    TH1D* relative_DAMPE_MixedHisto_1D_NS_EW_LS = nullptr;
    TH1D* relative_DAMPE_MixedHisto_1D_NS_FB_LS = nullptr;
    TH1D* relative_DAMPE_MixedHisto_1D_EW_FB_LS = nullptr;
    TH1D* relative_DAMPE_FullMixedHisto_1D_LS = nullptr;
    
    TH1D* relative_DAMPE_DataHisto_1D_I_HS = nullptr;
    TH1D* relative_DAMPE_DataHisto_1D_AniNS_HS = nullptr;
    TH1D* relative_DAMPE_DataHisto_1D_AniEW_HS = nullptr;
    TH1D* relative_DAMPE_DataHisto_1D_AniFB_HS = nullptr;
    TH1D* relative_DAMPE_MixedHisto_1D_NS_EW_HS = nullptr;
    TH1D* relative_DAMPE_MixedHisto_1D_NS_FB_HS = nullptr;
    TH1D* relative_DAMPE_MixedHisto_1D_EW_FB_HS = nullptr;
    TH1D* relative_DAMPE_FullMixedHisto_1D_HS = nullptr;
    
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
                            Template_Iso_LS,
                            Template_AniNS_LS,
                            Template_AniEW_LS,
                            Template_AniFB_LS,
                            output_log_file,
                            r_gen
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
                            Template_Iso_HS,
                            Template_AniNS_HS,
                            Template_AniEW_HS,
                            Template_AniFB_HS,
                            output_log_file,
                            r_gen
                           );
    
    
    ///////////////////////// Computing relative LS maps
    
    get_relative_histo(relative_DAMPE_Data_Iso_LS,DAMPE_Data_Iso_LS,pDAMPE_ReferenceMap_LS);
    get_relative_histo(relative_DAMPE_Data_AniNS_LS,DAMPE_Data_AniNS_LS,pDAMPE_ReferenceMap_LS);
    get_relative_histo(relative_DAMPE_Data_AniEW_LS,DAMPE_Data_AniEW_LS,pDAMPE_ReferenceMap_LS);
    get_relative_histo(relative_DAMPE_Data_AniFB_LS,DAMPE_Data_AniFB_LS,pDAMPE_ReferenceMap_LS);
    get_relative_histo(relative_DAMPE_MixedData_NS_EW_LS,DAMPE_MixedData_NS_EW_LS,pDAMPE_ReferenceMap_LS);
    get_relative_histo(relative_DAMPE_MixedData_NS_FB_LS,DAMPE_MixedData_NS_FB_LS,pDAMPE_ReferenceMap_LS);
    get_relative_histo(relative_DAMPE_MixedData_EW_FB_LS,DAMPE_MixedData_EW_FB_LS,pDAMPE_ReferenceMap_LS);
    get_relative_histo(relative_DAMPE_FullMixedData_LS,DAMPE_FullMixedData_LS,pDAMPE_ReferenceMap_LS);
    
    ///////////////////////// Computing relative HS maps
    
    get_relative_histo(relative_DAMPE_Data_Iso_HS,DAMPE_Data_Iso_HS,pDAMPE_ReferenceMap_HS);
    get_relative_histo(relative_DAMPE_Data_AniNS_HS,DAMPE_Data_AniNS_HS,pDAMPE_ReferenceMap_HS);
    get_relative_histo(relative_DAMPE_Data_AniEW_HS,DAMPE_Data_AniEW_HS,pDAMPE_ReferenceMap_HS);
    get_relative_histo(relative_DAMPE_Data_AniFB_HS,DAMPE_Data_AniFB_HS,pDAMPE_ReferenceMap_HS);
    get_relative_histo(relative_DAMPE_MixedData_NS_EW_HS,DAMPE_MixedData_NS_EW_HS,pDAMPE_ReferenceMap_HS);
    get_relative_histo(relative_DAMPE_MixedData_NS_FB_HS,DAMPE_MixedData_NS_FB_HS,pDAMPE_ReferenceMap_HS);
    get_relative_histo(relative_DAMPE_MixedData_EW_FB_HS,DAMPE_MixedData_EW_FB_HS,pDAMPE_ReferenceMap_HS);
    get_relative_histo(relative_DAMPE_FullMixedData_HS,DAMPE_FullMixedData_HS,pDAMPE_ReferenceMap_HS);
    
    ///////////////////////// Converting Templates Maps from TH2 to TH1
    
    TH2toTH1_obj(Templates_LS[0],Template_Iso_LS);
    TH2toTH1_obj(Templates_LS[1],Template_AniNS_LS);
    TH2toTH1_obj(Templates_LS[2],Template_AniEW_LS);
    TH2toTH1_obj(Templates_LS[3],Template_AniFB_LS);
    
    TH2toTH1_obj(Templates_HS[0],Template_Iso_HS);
    TH2toTH1_obj(Templates_HS[1],Template_AniNS_HS);
    TH2toTH1_obj(Templates_HS[2],Template_AniEW_HS);
    TH2toTH1_obj(Templates_HS[3],Template_AniFB_HS);
    
    ///////////////////////// Converting Data Maps from TH2 to TH1
    
    TH2toTH1_obj(relative_DAMPE_DataHisto_I_LS,relative_DAMPE_Data_Iso_LS);
    TH2toTH1_obj(relative_DAMPE_DataHisto_NS_LS,relative_DAMPE_Data_AniNS_LS);
    TH2toTH1_obj(relative_DAMPE_DataHisto_EW_LS,relative_DAMPE_Data_AniEW_LS);
    TH2toTH1_obj(relative_DAMPE_DataHisto_FB_LS,relative_DAMPE_Data_AniFB_LS);
    
    TH2toTH1_obj(relative_DAMPE_MixedDataHisto_NS_EW_LS,relative_DAMPE_MixedData_NS_EW_LS);
    TH2toTH1_obj(relative_DAMPE_MixedDataHisto_NS_FB_LS,relative_DAMPE_MixedData_NS_FB_LS);
    TH2toTH1_obj(relative_DAMPE_MixedDataHisto_EW_FB_LS,relative_DAMPE_MixedData_EW_FB_LS);
    TH2toTH1_obj(relative_DAMPE_FullMixedDataHisto_LS,relative_DAMPE_FullMixedData_LS);
    
    TH2toTH1_obj(relative_DAMPE_DataHisto_I_HS,relative_DAMPE_Data_Iso_HS);
    TH2toTH1_obj(relative_DAMPE_DataHisto_NS_HS,relative_DAMPE_Data_AniNS_HS);
    TH2toTH1_obj(relative_DAMPE_DataHisto_EW_HS,relative_DAMPE_Data_AniEW_HS);
    TH2toTH1_obj(relative_DAMPE_DataHisto_FB_HS,relative_DAMPE_Data_AniFB_HS);
    
    TH2toTH1_obj(relative_DAMPE_MixedDataHisto_NS_EW_HS,relative_DAMPE_MixedData_NS_EW_HS);
    TH2toTH1_obj(relative_DAMPE_MixedDataHisto_NS_FB_HS,relative_DAMPE_MixedData_NS_FB_HS);
    TH2toTH1_obj(relative_DAMPE_MixedDataHisto_EW_FB_HS,relative_DAMPE_MixedData_EW_FB_HS);
    TH2toTH1_obj(relative_DAMPE_FullMixedDataHisto_HS,relative_DAMPE_FullMixedData_HS);
    
    ///////////////////////// Linking variables to the pointers, prepearing for the fit
    
    for(Int_t idx_t = 0; idx_t < 4; idx_t++) {
        Templates_1D_LS[idx_t] = &Templates_LS[idx_t];
        Templates_1D_HS[idx_t] = &Templates_HS[idx_t];
    }
    
    relative_DAMPE_DataHisto_1D_I_LS = &relative_DAMPE_DataHisto_I_LS;
    relative_DAMPE_DataHisto_1D_AniNS_LS = &relative_DAMPE_DataHisto_NS_LS;
    relative_DAMPE_DataHisto_1D_AniEW_LS = &relative_DAMPE_DataHisto_EW_LS;
    relative_DAMPE_DataHisto_1D_AniFB_LS = &relative_DAMPE_DataHisto_FB_LS;
    
    relative_DAMPE_MixedHisto_1D_NS_EW_LS = &relative_DAMPE_MixedDataHisto_NS_EW_LS;
    relative_DAMPE_MixedHisto_1D_NS_FB_LS = &relative_DAMPE_MixedDataHisto_NS_FB_LS;
    relative_DAMPE_MixedHisto_1D_EW_FB_LS = &relative_DAMPE_MixedDataHisto_EW_FB_LS;
    relative_DAMPE_FullMixedHisto_1D_LS = &relative_DAMPE_FullMixedDataHisto_LS;
    
    relative_DAMPE_DataHisto_1D_I_HS = &relative_DAMPE_DataHisto_I_HS;
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
    
    std::cout << "\n\n -------------- > Isotropic Monopole \n\n";
    output_log_file << "\n\n -------------- > Isotropic Monopole \n\n";
    
    TemplateFitBH(relative_DAMPE_DataHisto_1D_I_LS,4,Templates_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,0,false);
    getPull(relative_DAMPE_DataHisto_1D_I_LS,Templates_1D_LS,res,relative_DAMPE_hPull_Iso_LS);
    
    std::cout << "\n\n -------------- > Anisotropic Monopole (NS) \n\n";
    output_log_file << "\n\n -------------- > Anisotropic Monopole (NS) \n\n";
    
    TemplateFitBH(relative_DAMPE_DataHisto_1D_AniNS_LS,4,Templates_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,1,false);
    getPull(relative_DAMPE_DataHisto_1D_AniNS_LS,Templates_1D_LS,res,relative_DAMPE_hPull_AniNS_LS);
    
    std::cout << "\n\n -------------- > Anisotropic Monopole (EW) \n\n";
    output_log_file << "\n\n -------------- > Anisotropic Monopole (EW) \n\n";
    
    TemplateFitBH(relative_DAMPE_DataHisto_1D_AniEW_LS,4,Templates_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,2,false);
    getPull(relative_DAMPE_DataHisto_1D_AniEW_LS,Templates_1D_LS,res,relative_DAMPE_hPull_AniEW_LS);
    
    std::cout << "\n\n -------------- > Anisotropic Monopole (FB) \n\n";
    output_log_file << "\n\n -------------- > Anisotropic Monopole (FB) \n\n";
    
    TemplateFitBH(relative_DAMPE_DataHisto_1D_AniFB_LS,4,Templates_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,3,false);
    getPull(relative_DAMPE_DataHisto_1D_AniFB_LS,Templates_1D_LS,res,relative_DAMPE_hPull_AniFB_LS);
    
    std::cout << "\n\n -------------- > NS + EW linear combination \n\n";
    output_log_file << "\n\n -------------- > NS + EW linear combination \n\n";
    
    TemplateFitBH(relative_DAMPE_MixedHisto_1D_NS_EW_LS,4,Templates_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,4,false);
    getPull(relative_DAMPE_MixedHisto_1D_NS_EW_LS,Templates_1D_LS,res,relative_DAMPE_hPull_Mixed_NS_EW_LS);
    
    std::cout << "\n\n -------------- > NS + FB linear combination \n\n";
    output_log_file << "\n\n -------------- > NS + FB linear combination \n\n";
    
    TemplateFitBH(relative_DAMPE_MixedHisto_1D_NS_FB_LS,4,Templates_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,5,false);
    getPull(relative_DAMPE_MixedHisto_1D_NS_FB_LS,Templates_1D_LS,res,relative_DAMPE_hPull_Mixed_NS_FB_LS);
    
    std::cout << "\n\n -------------- > EW + FB linear combination \n\n";
    output_log_file << "\n\n -------------- > EW + FB linear combination \n\n";
    
    TemplateFitBH(relative_DAMPE_MixedHisto_1D_EW_FB_LS,4,Templates_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,6,false);
    getPull(relative_DAMPE_MixedHisto_1D_EW_FB_LS,Templates_1D_LS,res,relative_DAMPE_hPull_Mixed_EW_FB_LS);
    
    std::cout << "\n\n -------------- > FullMixed linear combination \n\n";
    output_log_file << "\n\n -------------- > FullMixed linear combination \n\n";
    
    TemplateFitBH(relative_DAMPE_FullMixedHisto_1D_LS,4,Templates_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,7,false);
    getPull(relative_DAMPE_FullMixedHisto_1D_LS,Templates_1D_LS,res,relative_DAMPE_hPull_FullMixed_LS);
    
    
    std::cout << "\n\n //////////////////////////////////////// High Statistics DAMPE fits //////////////////////////////////////// \n\n";
    output_log_file << "\n\n //////////////////////////////////////// High Statistics DAMPE fits //////////////////////////////////////// \n\n";
    
    std::cout << "\n\n -------------- > Isotropic Monopole \n\n";
    output_log_file << "\n\n -------------- > Isotropic Monopole \n\n";
    
    TemplateFitBH(relative_DAMPE_DataHisto_1D_I_HS,4,Templates_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,0,true);
    getPull(relative_DAMPE_DataHisto_1D_I_HS,Templates_1D_HS,res,relative_DAMPE_hPull_Iso_HS);
    
    std::cout << "\n\n -------------- > Anisotropic Monopole (NS) \n\n";
    output_log_file << "\n\n -------------- > Anisotropic Monopole (NS) \n\n";
    
    TemplateFitBH(relative_DAMPE_DataHisto_1D_AniNS_HS,4,Templates_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,1,true);
    getPull(relative_DAMPE_DataHisto_1D_AniNS_HS,Templates_1D_HS,res,relative_DAMPE_hPull_AniNS_HS);
    
    std::cout << "\n\n -------------- > Anisotropic Monopole (EW) \n\n";
    output_log_file << "\n\n -------------- > Anisotropic Monopole (EW) \n\n";
    
    TemplateFitBH(relative_DAMPE_DataHisto_1D_AniEW_HS,4,Templates_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,2,true);
    getPull(relative_DAMPE_DataHisto_1D_AniEW_HS,Templates_1D_HS,res,relative_DAMPE_hPull_AniEW_HS);
    
    std::cout << "\n\n -------------- > Anisotropic Monopole (FB) \n\n";
    output_log_file << "\n\n -------------- > Anisotropic Monopole (FB) \n\n";
    
    TemplateFitBH(relative_DAMPE_DataHisto_1D_AniFB_HS,4,Templates_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,3,true);
    getPull(relative_DAMPE_DataHisto_1D_AniFB_HS,Templates_1D_HS,res,relative_DAMPE_hPull_AniFB_HS);
    
    std::cout << "\n\n -------------- > NS + EW linear combination \n\n";
    output_log_file << "\n\n -------------- > NS + EW linear combination \n\n";
    
    TemplateFitBH(relative_DAMPE_MixedHisto_1D_NS_EW_HS,4,Templates_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,4,true);
    getPull(relative_DAMPE_MixedHisto_1D_NS_EW_HS,Templates_1D_HS,res,relative_DAMPE_hPull_Mixed_NS_EW_HS);
    
    std::cout << "\n\n -------------- > NS + FB linear combination \n\n";
    output_log_file << "\n\n -------------- > NS + FB linear combination \n\n";
    
    TemplateFitBH(relative_DAMPE_MixedHisto_1D_NS_FB_HS,4,Templates_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,5,true);
    getPull(relative_DAMPE_MixedHisto_1D_NS_FB_HS,Templates_1D_HS,res,relative_DAMPE_hPull_Mixed_NS_FB_HS);
    
    std::cout << "\n\n -------------- > EW + FB linear combination \n\n";
    output_log_file << "\n\n -------------- > EW + FB linear combination \n\n";
    
    TemplateFitBH(relative_DAMPE_MixedHisto_1D_EW_FB_HS,4,Templates_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,6,true);
    getPull(relative_DAMPE_MixedHisto_1D_EW_FB_HS,Templates_1D_HS,res,relative_DAMPE_hPull_Mixed_EW_FB_HS);
    
    std::cout << "\n\n -------------- > FullMixed linear combination \n\n";
    output_log_file << "\n\n -------------- > FullMixed linear combination \n\n";
    
    TemplateFitBH(relative_DAMPE_FullMixedHisto_1D_HS,4,Templates_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,7,true);
    getPull(relative_DAMPE_FullMixedHisto_1D_HS,Templates_1D_HS,res,relative_DAMPE_hPull_FullMixed_HS);
    
    
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
        
        relative_DAMPE_Data_Iso_LS.Write();
        relative_DAMPE_Data_AniNS_LS.Write();
        relative_DAMPE_Data_AniEW_LS.Write();
        relative_DAMPE_Data_AniFB_LS.Write();
        
        relative_DAMPE_MixedData_NS_EW_LS.Write();
        relative_DAMPE_MixedData_NS_FB_LS.Write();
        relative_DAMPE_MixedData_EW_FB_LS.Write();
        relative_DAMPE_FullMixedData_LS.Write();
        
        relative_DAMPE_Data_Iso_HS.Write();
        relative_DAMPE_Data_AniNS_HS.Write();
        relative_DAMPE_Data_AniEW_HS.Write();
        relative_DAMPE_Data_AniFB_HS.Write();
        
        relative_DAMPE_MixedData_NS_EW_HS.Write();
        relative_DAMPE_MixedData_NS_FB_HS.Write();
        relative_DAMPE_MixedData_EW_FB_HS.Write();
        relative_DAMPE_FullMixedData_HS.Write();
        
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
        
        relative_DAMPE_hPull_Iso_HS.Fit("relative_DAMPE_gaus","L,E,M");
        relative_DAMPE_hPull_AniNS_HS.Fit("relative_DAMPE_gaus","L,E,M");
        relative_DAMPE_hPull_AniEW_HS.Fit("relative_DAMPE_gaus","L,E,M");
        relative_DAMPE_hPull_AniFB_HS.Fit("relative_DAMPE_gaus","L,E,M");
        relative_DAMPE_hPull_Mixed_NS_EW_HS.Fit("relative_DAMPE_gaus","L,E,M");
        relative_DAMPE_hPull_Mixed_NS_FB_HS.Fit("relative_DAMPE_gaus","L,E,M");
        relative_DAMPE_hPull_Mixed_EW_FB_HS.Fit("relative_DAMPE_gaus","L,E,M");
        relative_DAMPE_hPull_FullMixed_HS.Fit("relative_DAMPE_gaus","L,E,M");
        
        //////////////////////////////////// Writing Pools
        
        relative_DAMPE_hPull_Iso_LS.Write();
        relative_DAMPE_hPull_AniNS_LS.Write();
        relative_DAMPE_hPull_AniEW_LS.Write();
        relative_DAMPE_hPull_AniFB_LS.Write();
        relative_DAMPE_hPull_Mixed_NS_EW_LS.Write();
        relative_DAMPE_hPull_Mixed_NS_FB_LS.Write();
        relative_DAMPE_hPull_Mixed_EW_FB_LS.Write();
        relative_DAMPE_hPull_FullMixed_LS.Write();
        
        relative_DAMPE_hPull_Iso_HS.Write();
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
    
    fiTree.Fill();
    
}


#endif

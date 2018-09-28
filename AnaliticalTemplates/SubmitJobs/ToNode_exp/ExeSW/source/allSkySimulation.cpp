
#include "MyHead.h"

void allSky_singleTry_fit(
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
    
    std::string data_out_path = output_path_creator(try_idx,2,NS_anisotropy,EW_anisotropy,FB_anisotropy);
    std::string pools_out_path = output_path_creator(try_idx,3,NS_anisotropy,EW_anisotropy,FB_anisotropy);
    std::string projections_out_path = output_path_creator(try_idx,4,NS_anisotropy,EW_anisotropy,FB_anisotropy);
    
    tmp_fit.inputAni[0] = NS_anisotropy;
    tmp_fit.inputAni[1] = EW_anisotropy;
    tmp_fit.inputAni[2] = FB_anisotropy;
    
    tmp_fit.seed = tmp_seed;
    tmp_fit.seed_list_line = seed_line;
    
    //////////// TemplateFit variables
    
    Double_t res[4],res_err[4];
    Double_t initialValues[4]={1,0,0,0};
    
    double fullFitResults_HS[4][8];
    
    std::vector<TH1D*> DataProjections;
    DataProjections.resize(16);
    
    std::vector<TH1D*> TemplatesProjections;
    TemplatesProjections.resize(8);
    
    //////////// Open templates
    
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
                    output_log_file
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
    
    //////////// Pool's histos
   
    //////////// Low statistics pools
    
    TH1D hPull_Iso_LS("hPull_Iso_LS","Pull Isotropic Sky Map LS",50,-10,10);
    TH1D hPull_AniNS_LS("hPull_AniNS_LS","Pull Anisotropic Sky Map (NS) LS",50,-10,10);
    TH1D hPull_AniEW_LS("hPull_AniEW_LS","Pull Anisotropic Sky Map (EW) LS",50,-10,10);
    TH1D hPull_AniFB_LS("hPull_AniFB_LS","Pull Anisotropic Sky Map (FB) LS",50,-10,10);
    TH1D hPull_Mixed_NS_EW_LS("hPull_Mixed_NS_EW_LS","Pull Anisotropic Sky Map (NS-EW) LS",50,-10,10);
    TH1D hPull_Mixed_NS_FB_LS("hPull_Mixed_NS_FB_LS","Pull Anisotropic Sky Map (NS-FB) LS",50,-10,10);
    TH1D hPull_Mixed_EW_FB_LS("hPull_Mixed_EW_FB_LS","Pull Anisotropic Sky Map (EW-FB) LS",50,-10,10);
    TH1D hPull_FullMixed_LS("hPull_FullMixed_LS","Pull Full Mixed Anisotropic Sky Map LS",50,-10,10);
    
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
                        Template_Iso_LS,
                        Template_AniNS_LS,
                        Template_AniEW_LS,
                        Template_AniFB_LS,
                        output_log_file,
                        tmp_seed
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
                        tmp_seed
                     );
    
    /*
    //////////// Set fixed poissonian errors on isotropic maps
    
    for(Int_t bX=1; bX <= Data_Iso_LS->GetNbinsX(); ++bX)
        for(Int_t bY=1; bY <= Data_Iso_LS->GetNbinsX(); ++bY)
            Data_Iso_LS->SetBinError(bX,bY,TMath::Sqrt(data_all_sky_LS_events/(Data_Iso_LS->GetNbinsX()*Data_Iso_LS->GetNbinsY())));
    
    for(Int_t bX=1; bX <= Data_Iso_HS->GetNbinsX(); ++bX)
        for(Int_t bY=1; bY <= Data_Iso_HS->GetNbinsX(); ++bY)
            Data_Iso_HS->SetBinError(bX,bY,TMath::Sqrt(data_all_sky_HS_events/(Data_Iso_HS->GetNbinsX()*Data_Iso_HS->GetNbinsY())));
    */
    
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
    
    TH2toTH1_ptr(MixedDataHisto_NS_EW_LS,MixedData_NS_EW_LS);
    TH2toTH1_ptr(MixedDataHisto_NS_FB_LS,MixedData_NS_FB_LS);
    TH2toTH1_ptr(MixedDataHisto_EW_FB_LS,MixedData_EW_FB_LS);
    TH2toTH1_ptr(FullMixedDataHisto_LS,FullMixedData_LS);
    
    TH2toTH1_ptr(DataHisto_I_HS,Data_Iso_HS);
    TH2toTH1_ptr(DataHisto_NS_HS,Data_AniNS_HS);
    TH2toTH1_ptr(DataHisto_EW_HS,Data_AniEW_HS);
    TH2toTH1_ptr(DataHisto_FB_HS,Data_AniFB_HS);
    
    TH2toTH1_ptr(MixedDataHisto_NS_EW_HS,MixedData_NS_EW_HS);
    TH2toTH1_ptr(MixedDataHisto_NS_FB_HS,MixedData_NS_FB_HS);
    TH2toTH1_ptr(MixedDataHisto_EW_FB_HS,MixedData_EW_FB_HS);
    TH2toTH1_ptr(FullMixedDataHisto_HS,FullMixedData_HS);
    
    //////////// Saving bin contents into the Tree
    
    for(Int_t idx = 0; idx < DataHisto_I_LS.GetNbinsX(); ++idx)
    {
        
        tmp_fit.entries_Iso_Map_LS[idx_ani][idx] = DataHisto_I_LS.GetBinContent(idx+1);
        tmp_fit.entries_NS_Map_LS[idx_ani][idx] = DataHisto_NS_LS.GetBinContent(idx+1);
        tmp_fit.entries_EW_Map_LS[idx_ani][idx] = DataHisto_EW_LS.GetBinContent(idx+1);
        tmp_fit.entries_FB_Map_LS[idx_ani][idx] = DataHisto_FB_LS.GetBinContent(idx+1);
        tmp_fit.entries_NS_EW_Map_LS[idx_ani][idx] = MixedDataHisto_NS_EW_LS.GetBinContent(idx+1);
        tmp_fit.entries_NS_FB_Map_LS[idx_ani][idx] = MixedDataHisto_NS_FB_LS.GetBinContent(idx+1);
        tmp_fit.entries_EW_FB_Map_LS[idx_ani][idx] = MixedDataHisto_EW_FB_LS.GetBinContent(idx+1);
        tmp_fit.entries_Full_Map_LS[idx_ani][idx] = FullMixedDataHisto_LS.GetBinContent(idx+1);
        
    }
    
    for(Int_t idx = 0; idx < DataHisto_I_HS.GetNbinsX(); ++idx)
    {
        
        tmp_fit.entries_Iso_Map_HS[idx_ani][idx] = DataHisto_I_HS.GetBinContent(idx+1);
        tmp_fit.entries_NS_Map_HS[idx_ani][idx] = DataHisto_NS_HS.GetBinContent(idx+1);
        tmp_fit.entries_EW_Map_HS[idx_ani][idx] = DataHisto_EW_HS.GetBinContent(idx+1);
        tmp_fit.entries_FB_Map_HS[idx_ani][idx] = DataHisto_FB_HS.GetBinContent(idx+1);
        tmp_fit.entries_NS_EW_Map_HS[idx_ani][idx] = MixedDataHisto_NS_EW_HS.GetBinContent(idx+1);
        tmp_fit.entries_NS_FB_Map_HS[idx_ani][idx] = MixedDataHisto_NS_FB_HS.GetBinContent(idx+1);
        tmp_fit.entries_EW_FB_Map_HS[idx_ani][idx] = MixedDataHisto_EW_FB_HS.GetBinContent(idx+1);
        tmp_fit.entries_Full_Map_HS[idx_ani][idx] = FullMixedDataHisto_HS.GetBinContent(idx+1);
        
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
                            output_log_file,
                            NS_anisotropy,
                            EW_anisotropy,
                            FB_anisotropy,
                            try_idx
                        );
    
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    std::cout<<"\n\nSimulation Completed !\n\n";
    output_log_file << "\n\nSimulation Completed !\n\n";
    
    
    if(write_tmp_histos && itry==0)
    {
        
        //////////////////////////////////// Creating Data out file
        
        TFile data_file(data_out_path.c_str(),"RECREATE");
        if(data_file.IsZombie()) {
            std::cout << "\n\nError writing ROOT Data TFile. Prorgram finished \n\n";
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
        
        
        //////////////////////////////////// Creating Pools out file
        
        TFile pool_file(pools_out_path.c_str(),"RECREATE");
        if(pool_file.IsZombie()) {
            std::cout << "\n\nError writing ROOT Pools TFile. Prorgram finished \n\n";
            output_log_file << "\n\nError writing ROOT Pools TFile. Prorgram finished \n\n";
            exit(100);
        }
        
        //////////////////////////////////// Gaussian Fit LS hPools
        
        hPull_Iso_LS.Fit("gaus","L,E,M");
        hPull_AniNS_LS.Fit("gaus","L,E,M");
        hPull_AniEW_LS.Fit("gaus","L,E,M");
        hPull_AniFB_LS.Fit("gaus","L,E,M");
        hPull_Mixed_NS_EW_LS.Fit("gaus","L,E,M");
        hPull_Mixed_NS_FB_LS.Fit("gaus","L,E,M");
        hPull_Mixed_EW_FB_LS.Fit("gaus","L,E,M");
        hPull_FullMixed_LS.Fit("gaus","L,E,M");
        
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

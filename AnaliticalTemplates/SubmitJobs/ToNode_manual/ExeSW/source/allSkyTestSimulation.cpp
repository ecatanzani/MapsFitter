
#include "MyHead.h"

void AllSky_test_simulation(
                                std::ofstream &output_log_file,
                                TTree &fiTree,
                                test_fitResult &tmp_fit,
                                UInt_t tmp_seed,
                                Int_t itry,
                                Double_t NS_anisotropy,
                                Double_t EW_anisotropy,
                                Double_t FB_anisotropy,
                                Int_t idx_ani,
                                Int_t s_idx,
                                UInt_t seed_line
                            )

{
    
    std::string data_out_path = output_path_creator(s_idx,2,NS_anisotropy,EW_anisotropy,FB_anisotropy,false,true);
    std::string pools_out_path = output_path_creator(s_idx,3,NS_anisotropy,EW_anisotropy,FB_anisotropy,false,true);
    
    tmp_fit.inputAni[0] = NS_anisotropy;
    tmp_fit.inputAni[1] = EW_anisotropy;
    tmp_fit.inputAni[2] = FB_anisotropy;
    
    tmp_fit.seed = tmp_seed;
    tmp_fit.seed_list_line = seed_line;
    
    //////////// TemplateFit variables
    
    Double_t res[1],res_err[1],initialValues[1];
    
    res[0] = 0;
    res_err[0] = 0;
    initialValues[0] = 1;

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
    
    Data_Iso_LS->SetTitle("Isotropic (test LS) All Sky (data) Map; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries");
    Data_AniNS_LS->SetTitle("Anisotropic (test LS) All Sky (data) Map (NS); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries");
    Data_AniEW_LS->SetTitle("Anisotropic (test LS) All Sky (data) Map (EW); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries");
    Data_AniFB_LS->SetTitle("Anisotropic (test LS) All Sky (data) Map (FB); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries");
    
    Data_Iso_LS->Reset();
    Data_AniNS_LS->Reset();
    Data_AniEW_LS->Reset();
    Data_AniFB_LS->Reset();
    
    //////////// High statistics data
    
    TH2D* Data_Iso_HS = (TH2D*) Template_Iso_HS.Clone("Data_Iso_HS");
    TH2D* Data_AniNS_HS = (TH2D*) Template_Iso_HS.Clone("Data_AniNS_HS");
    TH2D* Data_AniEW_HS = (TH2D*) Template_Iso_HS.Clone("Data_AniEW_HS");
    TH2D* Data_AniFB_HS = (TH2D*) Template_Iso_HS.Clone("Data_AniFB_HS");
    
    Data_Iso_HS->SetTitle("Isotropic (test HS) All Sky (data) Map; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries");
    Data_AniNS_HS->SetTitle("Anisotropic (test HS) All Sky (data) Map (NS); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries");
    Data_AniEW_HS->SetTitle("Anisotropic (test HS) All Sky (data) Map (EW); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries");
    Data_AniFB_HS->SetTitle("Anisotropic (test HS) All Sky (data) Map (FB); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries");
    
    Data_Iso_HS->Reset();
    Data_AniNS_HS->Reset();
    Data_AniEW_HS->Reset();
    Data_AniFB_HS->Reset();
    
    
    //////////// Pool's histos
    
    //////////// Low statistics pools
    
    TH1D hPull_Iso_LS("hPull_Iso_LS","Pull Isotropic Sky Map (test LS)",50,-10,10);
    TH1D hPull_AniNS_LS("hPull_AniNS_LS","Pull Anisotropic Sky Map (NS) (test LS)",50,-10,10);
    TH1D hPull_AniEW_LS("hPull_AniEW_LS","Pull Anisotropic Sky Map (EW) (test LS)",50,-10,10);
    TH1D hPull_AniFB_LS("hPull_AniFB_LS","Pull Anisotropic Sky Map (FB) (test LS)",50,-10,10);
    
    //////////// High statistics pools
    
    TH1D hPull_Iso_HS("hPull_Iso_HS","Pull Isotropic Sky Map (test HS)",50,-10,10);
    TH1D hPull_AniNS_HS("hPull_AniNS_HS","Pull Anisotropic Sky Map (NS) (test HS)",50,-10,10);
    TH1D hPull_AniEW_HS("hPull_AniEW_HS","Pull Anisotropic Sky Map (EW) (test HS)",50,-10,10);
    TH1D hPull_AniFB_HS("hPull_AniFB_HS","Pull Anisotropic Sky Map (FB) (test HS)",50,-10,10);
    
    //////////// 1D Histos
    
    TH1D Templates_LS[4];
    TH1D Templates_HS[4];
    
    TH1D DataHisto_I_LS;
    TH1D DataHisto_NS_LS;
    TH1D DataHisto_EW_LS;
    TH1D DataHisto_FB_LS;
    
    TH1D DataHisto_I_HS;
    TH1D DataHisto_NS_HS;
    TH1D DataHisto_EW_HS;
    TH1D DataHisto_FB_HS;
    
    TH1* Iso_Template_1D_LS[1] = {nullptr};
    TH1* NS_Template_1D_LS[1] = {nullptr};
    TH1* EW_Template_1D_LS[1] = {nullptr};
    TH1* FB_Template_1D_LS[1] = {nullptr};
    
    TH1* Iso_Template_1D_HS[1] = {nullptr};
    TH1* NS_Template_1D_HS[1] = {nullptr};
    TH1* EW_Template_1D_HS[1] = {nullptr};
    TH1* FB_Template_1D_HS[1] = {nullptr};
    
    TH1D* DataHisto_1D_I_LS = nullptr;
    TH1D* DataHisto_1D_AniNS_LS = nullptr;
    TH1D* DataHisto_1D_AniEW_LS = nullptr;
    TH1D* DataHisto_1D_AniFB_LS = nullptr;
    
    TH1D* DataHisto_1D_I_HS = nullptr;
    TH1D* DataHisto_1D_AniNS_HS = nullptr;
    TH1D* DataHisto_1D_AniEW_HS = nullptr;
    TH1D* DataHisto_1D_AniFB_HS = nullptr;
    
    
    tmp_fit.theta_binHistos_LS = Data_Iso_LS->GetNbinsY();
    tmp_fit.phi_binHistos_LS = Data_Iso_LS->GetNbinsX();
    tmp_fit.theta_binHistos_HS = Data_Iso_HS->GetNbinsY();
    tmp_fit.phi_binHistos_HS = Data_Iso_HS->GetNbinsX();
    
    tmp_fit.events_LS = data_all_sky_LS_events;
    tmp_fit.events_HS = data_all_sky_HS_events;
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    ///////////// All sky maps
    
    generate_LS_example_data(
                                NS_anisotropy,
                                EW_anisotropy,
                                FB_anisotropy,
                                Data_Iso_LS,
                                Data_AniNS_LS,
                                Data_AniEW_LS,
                                Data_AniFB_LS,
                                Template_Iso_LS,
                                Template_AniNS_LS,
                                Template_AniEW_LS,
                                Template_AniFB_LS,
                                output_log_file,
                                tmp_seed
                            );
    
    generate_HS_example_data(
                                NS_anisotropy,
                                EW_anisotropy,
                                FB_anisotropy,
                                Data_Iso_HS,
                                Data_AniNS_HS,
                                Data_AniEW_HS,
                                Data_AniFB_HS,
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
    
    TH2toTH1_ptr(DataHisto_I_HS,Data_Iso_HS);
    TH2toTH1_ptr(DataHisto_NS_HS,Data_AniNS_HS);
    TH2toTH1_ptr(DataHisto_EW_HS,Data_AniEW_HS);
    TH2toTH1_ptr(DataHisto_FB_HS,Data_AniFB_HS);
    
    ///////////////////////// Linking variables to the pointers, prepearing for the fit
    
    Iso_Template_1D_LS[0] = &Templates_LS[0];
    NS_Template_1D_LS[0] = &Templates_LS[1];
    EW_Template_1D_LS[0] = &Templates_LS[2];
    FB_Template_1D_LS[0] = &Templates_LS[3];
    
    Iso_Template_1D_HS[0] = &Templates_HS[0];
    NS_Template_1D_HS[0] = &Templates_HS[1];
    EW_Template_1D_HS[0] = &Templates_HS[2];
    FB_Template_1D_HS[0] = &Templates_HS[3];
    
    DataHisto_1D_I_LS = &DataHisto_I_LS;
    DataHisto_1D_AniNS_LS = &DataHisto_NS_LS;
    DataHisto_1D_AniEW_LS = &DataHisto_EW_LS;
    DataHisto_1D_AniFB_LS = &DataHisto_FB_LS;
    
    DataHisto_1D_I_HS = &DataHisto_I_HS;
    DataHisto_1D_AniNS_HS = &DataHisto_NS_HS;
    DataHisto_1D_AniEW_HS = &DataHisto_EW_HS;
    DataHisto_1D_AniFB_HS = &DataHisto_FB_HS;
   
    /////////////////////////////////////////////////////////////////////////////// Fitting !!!
    
    std::cout << "\n\n //////////////////////////////////////// Low Statistics fits //////////////////////////////////////// \n\n";
    output_log_file << "\n\n //////////////////////////////////////// Low Statistics fits //////////////////////////////////////// \n\n";
    
    std::cout << "\n\n -------------- > Isotropic Monopole \n\n";
    output_log_file << "\n\n -------------- > Isotropic Monopole \n\n";
    
    TemplateFitBH_test(DataHisto_1D_I_LS,1,Iso_Template_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,0,false);
    getPull(DataHisto_1D_I_LS,Iso_Template_1D_LS,res,hPull_Iso_LS,false,true);
    
    std::cout << "\n\n -------------- > Anisotropic Monopole (NS) \n\n";
    output_log_file << "\n\n -------------- > Anisotropic Monopole (NS) \n\n";
    
    TemplateFitBH_test(DataHisto_1D_AniNS_LS,1,NS_Template_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,1,false);
    getPull(DataHisto_1D_AniNS_LS,NS_Template_1D_LS,res,hPull_AniNS_LS,false,true);
    
    std::cout << "\n\n -------------- > Anisotropic Monopole (EW) \n\n";
    output_log_file << "\n\n -------------- > Anisotropic Monopole (EW) \n\n";
    
    TemplateFitBH_test(DataHisto_1D_AniEW_LS,1,EW_Template_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,2,false);
    getPull(DataHisto_1D_AniEW_LS,EW_Template_1D_LS,res,hPull_AniEW_LS,false,true);
    
    std::cout << "\n\n -------------- > Anisotropic Monopole (FB) \n\n";
    output_log_file << "\n\n -------------- > Anisotropic Monopole (FB) \n\n";
    
    TemplateFitBH_test(DataHisto_1D_AniFB_LS,1,FB_Template_1D_LS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,3,false);
    getPull(DataHisto_1D_AniFB_LS,FB_Template_1D_LS,res,hPull_AniFB_LS,false,true);
    
    std::cout << "\n\n //////////////////////////////////////// High Statistics fits //////////////////////////////////////// \n\n";
    output_log_file << "\n\n //////////////////////////////////////// High Statistics fits //////////////////////////////////////// \n\n";
    
    std::cout << "\n\n -------------- > Isotropic Monopole \n\n";
    output_log_file << "\n\n -------------- > Isotropic Monopole \n\n";
    
    TemplateFitBH_test(DataHisto_1D_I_HS,1,Iso_Template_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,0,true);
    getPull(DataHisto_1D_I_HS,Iso_Template_1D_HS,res,hPull_Iso_HS,false,true);
    
    std::cout << "\n\n -------------- > Anisotropic Monopole (NS) \n\n";
    output_log_file << "\n\n -------------- > Anisotropic Monopole (NS) \n\n";
    
    TemplateFitBH_test(DataHisto_1D_AniNS_HS,1,NS_Template_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,1,true);
    getPull(DataHisto_1D_AniNS_HS,NS_Template_1D_HS,res,hPull_AniNS_HS,false,true);
    
    std::cout << "\n\n -------------- > Anisotropic Monopole (EW) \n\n";
    output_log_file << "\n\n -------------- > Anisotropic Monopole (EW) \n\n";
    
    TemplateFitBH_test(DataHisto_1D_AniEW_HS,1,EW_Template_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,2,true);
    getPull(DataHisto_1D_AniEW_HS,EW_Template_1D_HS,res,hPull_AniEW_HS,false,true);
    
    std::cout << "\n\n -------------- > Anisotropic Monopole (FB) \n\n";
    output_log_file << "\n\n -------------- > Anisotropic Monopole (FB) \n\n";
    
    TemplateFitBH_test(DataHisto_1D_AniFB_HS,1,FB_Template_1D_HS,res,res_err,initialValues,false,false,false,true,output_log_file,tmp_fit,3,true);
    getPull(DataHisto_1D_AniFB_HS,FB_Template_1D_HS,res,hPull_AniFB_HS,false,true);
    
    
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
        
        Data_Iso_HS->Write();
        Data_AniNS_HS->Write();
        Data_AniEW_HS->Write();
        Data_AniFB_HS->Write();
    
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
        
        //////////////////////////////////// Gaussian Fit HS hPools
        
        hPull_Iso_HS.Fit("gaus","L,E,M");
        hPull_AniNS_HS.Fit("gaus","L,E,M");
        hPull_AniEW_HS.Fit("gaus","L,E,M");
        hPull_AniFB_HS.Fit("gaus","L,E,M");
        
        //////////////////////////////////// Writing Pools
        
        hPull_Iso_LS.Write();
        hPull_AniNS_LS.Write();
        hPull_AniEW_LS.Write();
        hPull_AniFB_LS.Write();
        
        hPull_Iso_HS.Write();
        hPull_AniNS_HS.Write();
        hPull_AniEW_HS.Write();
        hPull_AniFB_HS.Write();
        
        pool_file.Write();
        pool_file.Close();
        
    }
    
        fiTree.Fill();
        
}

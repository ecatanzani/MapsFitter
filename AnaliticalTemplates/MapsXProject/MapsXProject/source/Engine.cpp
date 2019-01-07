
#include "MyHead.h"


void generate_data_interface(
                                std::ofstream &output_log_file,
                                std::string output_log,
                                std::string output_root,
                                std::string DAMPE_Iso_Map,
                                std::string seeds_path,
                                time_t time_stamp,
                                std::vector<Double_t> &NS_anisotropy,
                                std::vector<Double_t> &EW_anisotropy,
                                std::vector<Double_t> &FB_anisotropy,
                                ULong64_t data_LS_events,
                                ULong64_t data_HS_events,
                                UInt_t ani_values,
                                Bool_t all_sky_simulation,
                                Bool_t DAMPE_simulation,
                                Bool_t DAMPE_relative_simulation
                             )
{
    
    //////////// Seed stuff
    
    std::string tmp_seed_str;
    UInt_t tmp_seed = 0;
    std::ifstream inSeed(seeds_path);
    if(!inSeed.is_open()) {
        std::cerr << "\n\nCannot open input seeds file! Program finished !" << std::endl;
        exit(100);
    }
    
    
    //////////// Compute scaled reference maps
    
    TH2D DAMPE_ReferenceMap_LS,DAMPE_ReferenceMap_HS;
    
    read_DAMPE_FullIso(
                           DAMPE_ReferenceMap_LS,
                           DAMPE_ReferenceMap_HS,
                           data_LS_events,
                           data_HS_events,
                           DAMPE_Iso_Map
                       );
    
    
    for(Int_t idx_ani = 0; idx_ani < ani_values; ++idx_ani)
    {
        
        std::cout << "\n\nNS Anisotropy: " << NS_anisotropy[idx_ani] << "\n";
        std::cout << "EW Anisotropy: " << EW_anisotropy[idx_ani] << "\n";
        std::cout << "FB Anisotropy: " << FB_anisotropy[idx_ani] << "\n\n";
        
        output_log_file << "\n\nNS Anisotropy: " << NS_anisotropy[idx_ani] << "\n";
        output_log_file << "EW Anisotropy: " << EW_anisotropy[idx_ani] << "\n";
        output_log_file << "FB Anisotropy: " << FB_anisotropy[idx_ani] << "\n\n";
        
        inSeed >> tmp_seed_str;
        tmp_seed = (UInt_t)std::stoul(tmp_seed_str,nullptr,10);
    
        if(all_sky_simulation)
            interface_allSky_simulation(
                                            output_log_file,
                                            output_log,
                                            output_root,
                                            time_stamp,
                                            NS_anisotropy[idx_ani],
                                            EW_anisotropy[idx_ani],
                                            FB_anisotropy[idx_ani],
                                            tmp_seed,
                                            data_LS_events,
                                            data_HS_events
                                        );
 
        
    
        if(DAMPE_simulation)
            interface_DAMPE_simulation(
                                            output_log_file,
                                            output_log,
                                            output_root,
                                            time_stamp,
                                            NS_anisotropy[idx_ani],
                                            EW_anisotropy[idx_ani],
                                            FB_anisotropy[idx_ani],
                                            tmp_seed,
                                            data_LS_events,
                                            data_HS_events
                                        );
        
    
        if(DAMPE_relative_simulation)
            interface_DAMPE_relative_simulation(
                                                    output_log_file,
                                                    output_log,
                                                    output_root,
                                                    time_stamp,
                                                    NS_anisotropy[idx_ani],
                                                    EW_anisotropy[idx_ani],
                                                    FB_anisotropy[idx_ani],
                                                    tmp_seed,
                                                    data_LS_events,
                                                    data_HS_events,
                                                    DAMPE_ReferenceMap_LS,
                                                    DAMPE_ReferenceMap_HS
                                                );
    }
        
}


void generate_templates(
                            std::ofstream &output_log_file,
                            ULong64_t data_LS_events,
                            ULong64_t data_HS_events,
                            std::string output_log,
                            std::string output_root,
                            std::string DAMPE_Iso_Map,
                            time_t time_stamp
                        )
{
    //////////// Compute scaled reference maps
    
    TH2D DAMPE_ReferenceMap_LS,DAMPE_ReferenceMap_HS;
    
    read_DAMPE_FullIso(
                           DAMPE_ReferenceMap_LS,
                           DAMPE_ReferenceMap_HS,
                           data_LS_events,
                           data_HS_events,
                           DAMPE_Iso_Map
                       );
    
    //////////// Compute templates...
    
    std::string template_out_path = output_path_creator(output_log,output_root,1,time_stamp);
    std::string DAMPE_template_out_path = output_path_creator(output_log,output_root,1,time_stamp,0,0,0,true);
    
    compute_templates(
                          output_log,
                          output_root,
                          time_stamp,
                          output_log_file,
                          DAMPE_ReferenceMap_LS,
                          DAMPE_ReferenceMap_HS,
                          template_out_path,
                          DAMPE_template_out_path
                      );
    
}

/*

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
*/

void generate_and_fit(
                        ULong64_t data_LS_events,
                        ULong64_t data_HS_events,
                        Bool_t write_tmp_histos,
                        Bool_t all_sky_simulation,
                        Bool_t DAMPE_simulation,
                        Bool_t DAMPE_relative_simulation,
                        std::vector<Double_t> &NS_anisotropy,
                        std::vector<Double_t> &EW_anisotropy,
                        std::vector<Double_t> &FB_anisotropy,
                        std::string output_log,
                        std::string output_root,
                        std::string DAMPE_Iso_Map,
                        std::string seeds_path,
                        UInt_t ani_values,
                        std::ofstream &output_log_file,
                        time_t time_stamp
                      )
{

    //////////// Seed stuff
    
    std::string tmp_seed_str;
    UInt_t tmp_seed = 0;
    std::ifstream inSeed(seeds_path);
    if(!inSeed.is_open()) {
        std::cerr << "\n\nCannot open input seeds file! Program finished !" << std::endl;
        exit(100);
    }
    
    //////////// Compute scaled reference maps
    
    TH2D DAMPE_ReferenceMap_LS,DAMPE_ReferenceMap_HS;
    
    read_DAMPE_FullIso(
                           DAMPE_ReferenceMap_LS,
                           DAMPE_ReferenceMap_HS,
                           data_LS_events,
                           data_HS_events,
                           DAMPE_Iso_Map
                       );
    
    //////////// Compute templates...
    
    std::string template_out_path = output_path_creator(output_log,output_root,1,time_stamp);
    std::string DAMPE_template_out_path = output_path_creator(output_log,output_root,1,time_stamp,0,0,0,true);
    
    compute_templates(
                        output_log,
                        output_root,
                        time_stamp,
                        output_log_file,
                        DAMPE_ReferenceMap_LS,
                        DAMPE_ReferenceMap_HS,
                        template_out_path,
                        DAMPE_template_out_path
                      );
    
    ////////////////////////////////////////////////////////////////////// Start Simulation ////////////////////////////////
    
    for(Int_t idx_ani = 0; idx_ani < ani_values; ++idx_ani)
    {
        
        std::cout << "\n\nNS Anisotropy: " << NS_anisotropy[idx_ani] << "\n";
        std::cout << "EW Anisotropy: " << EW_anisotropy[idx_ani] << "\n";
        std::cout << "FB Anisotropy: " << FB_anisotropy[idx_ani] << "\n\n";
        
        output_log_file << "\n\nNS Anisotropy: " << NS_anisotropy[idx_ani] << "\n";
        output_log_file << "EW Anisotropy: " << EW_anisotropy[idx_ani] << "\n";
        output_log_file << "FB Anisotropy: " << FB_anisotropy[idx_ani] << "\n\n";
        
        inSeed >> tmp_seed_str;
        tmp_seed = (UInt_t)std::stoul(tmp_seed_str,nullptr,10);
        
        if(all_sky_simulation)
            allSky_singleTry_fit(
                                    output_log_file,
                                    output_log,
                                    output_root,
                                    time_stamp,
                                    tmp_seed,
                                    NS_anisotropy[idx_ani],
                                    EW_anisotropy[idx_ani],
                                    FB_anisotropy[idx_ani],
                                    template_out_path,
                                    DAMPE_template_out_path,
                                    write_tmp_histos,
                                    data_LS_events,
                                    data_HS_events
                                 );
        
        if(DAMPE_simulation)
            DAMPE_singleTry_fit(
                                    output_log_file,
                                    output_log,
                                    output_root,
                                    time_stamp,
                                    tmp_seed,
                                    NS_anisotropy[idx_ani],
                                    EW_anisotropy[idx_ani],
                                    FB_anisotropy[idx_ani],
                                    template_out_path,
                                    DAMPE_template_out_path,
                                    data_LS_events,
                                    data_HS_events,
                                    write_tmp_histos
                                );
        
        
        if(DAMPE_relative_simulation)
            DAMPE_relative_singleTry_fit(
                                            output_log_file,
                                            output_log,
                                            output_root,
                                            time_stamp,
                                            tmp_seed,
                                            NS_anisotropy[idx_ani],
                                            EW_anisotropy[idx_ani],
                                            FB_anisotropy[idx_ani],
                                            template_out_path,
                                            DAMPE_template_out_path,
                                            DAMPE_ReferenceMap_LS,
                                            DAMPE_ReferenceMap_HS,
                                            data_LS_events,
                                            data_HS_events,
                                            write_tmp_histos
                                         );
        
        
    }
    
    inSeed.close();
    
}

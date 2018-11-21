
#include "MyHead.h"

void generate_and_fit(std::ofstream &output_log_file,Int_t s_idx,Int_t s_batch)
{
    
    //////////// Generating TTree

    TTree fiTree("fiTree","TemplateFit results TTree");
    fitResult tmp_fit;
    std::string tree_out_path = output_path_creator(s_idx,6);

    //////////// Simulation data parameters
    
    std::string tmp_seed_str;
    UInt_t tmp_seed = 0;
    std::ifstream inSeed(seeds_path);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    create_tree_branches(fiTree,tmp_fit);
    
    deploy_simulation(
                        output_log_file,
                        inSeed,
                        s_idx,
                        s_batch,
                        tmp_seed_str,
                        tmp_seed,
                        fiTree,
                        tmp_fit
                      );
    
    //////////////////////////////////// Creating TTree out file
    
    TFile tree_file(tree_out_path.c_str(),"RECREATE");
    if(tree_file.IsZombie()) {
        std::cerra << "\n\nError writing ROOT TTree TFile. Prorgram finished \n\n";
        output_log_file << "\n\nError writing ROOT TTree TFile. Prorgram finished \n\n";
        exit(100);
    }
    
    fiTree.Write();
    
    tree_file.Write();
    tree_file.Close();
    
}

void generate_and_fit_relative(std::ofstream &output_log_file,Int_t s_idx,Int_t s_batch)
{
    
    //////////// Generating TTree
    
    TTree fiTree("fiTree","TemplateFit results TTree");
    relative_fitResult tmp_fit;
    std::string tree_out_path = output_path_creator(s_idx,6);
    
    //////////// Simulation data parameters
    
    std::string tmp_seed_str;
    UInt_t tmp_seed = 0;
    std::ifstream inSeed(seeds_path);
    
    ///////////////////////////////////////// DAMPE's Reference Histos
    
    //TH2D DAMPE_ReferenceMap;
    TH2D DAMPE_ReferenceMap_LS;
    TH2D DAMPE_ReferenceMap_HS;
    
    //read_DAMPE_FullIso(DAMPE_ReferenceMap);
    get_scaled_isotropic_DAMPE_maps(DAMPE_ReferenceMap_LS,DAMPE_ReferenceMap_HS);
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    create_rel_tree_branches(fiTree,tmp_fit);
    
    deploy_relative_simulation(
                                output_log_file,
                                inSeed,
                                s_idx,
                                s_batch,
                                tmp_seed_str,
                                tmp_seed,
                                fiTree,
                                tmp_fit,
                                DAMPE_ReferenceMap_LS,
                                DAMPE_ReferenceMap_HS
                               );
    
    //////////////////////////////////// Creating TTree out file
    
    TFile tree_file(tree_out_path.c_str(),"RECREATE");
    if(tree_file.IsZombie()) {
        std::cout << "\n\nError writing ROOT TTree TFile. Prorgram finished \n\n";
        output_log_file << "\n\nError writing ROOT TTree TFile. Prorgram finished \n\n";
        exit(100);
    }
    
    fiTree.Write();
    
    tree_file.Write();
    tree_file.Close();
    
}


void generate_and_fit_test(std::ofstream &output_log_file,Int_t s_idx,Int_t s_batch)
{
    
    //////////// Generating TTree
    
    TTree fiTree("fiTree","TemplateFit results TTree");
    test_fitResult tmp_fit;
    std::string tree_out_path = output_path_creator(s_idx,6);
    
    //////////// Simulation data parameters
    
    std::string tmp_seed_str;
    UInt_t tmp_seed = 0;
    std::ifstream inSeed(seeds_path);
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    create_test_tree_branches(fiTree,tmp_fit);
    
    deploy_test_simulation(
                            output_log_file,
                            inSeed,
                            s_idx,
                            s_batch,
                            tmp_seed_str,
                            tmp_seed,
                            fiTree,
                            tmp_fit
                           );
    
    //////////////////////////////////// Creating TTree out file
    
    TFile tree_file(tree_out_path.c_str(),"RECREATE");
    if(tree_file.IsZombie()) {
        std::cout << "\n\nError writing ROOT TTree TFile. Prorgram finished \n\n";
        output_log_file << "\n\nError writing ROOT TTree TFile. Prorgram finished \n\n";
        exit(100);
    }
    
    fiTree.Write();
    
    tree_file.Write();
    tree_file.Close();
    
    
}


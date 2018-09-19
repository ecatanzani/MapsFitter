
#include <iostream>
#include <stdlib.h>
#include <stdio.h>

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"

#define ani_values 3

class relative_fitResult
{
public:
    
    Double_t inputAni[3];
    ULong64_t seed;
    ULong64_t seed_list_line;
    
    Double_t chi2_LS[7];
    Double_t ndf_LS[7];
    Double_t chi2_r_LS[7];
    Double_t fit_par_LS[7][3];
    Double_t fit_err_LS[7][3];
    Double_t delta_LS[7];
    Double_t sum_par_LS[7];
    
    Double_t CMatrix_NS_LS[3][3];
    Double_t CMatrix_EW_LS[3][3];
    Double_t CMatrix_FB_LS[3][3];
    Double_t CMatrix_NS_EW_LS[3][3];
    Double_t CMatrix_NS_FB_LS[3][3];
    Double_t CMatrix_EW_FB_LS[3][3];
    Double_t CMatrix_Full_LS[3][3];
    
    Int_t entries_NS_Map_LS[ani_values][648];
    Int_t entries_EW_Map_LS[ani_values][648];
    Int_t entries_FB_Map_LS[ani_values][648];
    Int_t entries_NS_EW_Map_LS[ani_values][648];
    Int_t entries_NS_FB_Map_LS[ani_values][648];
    Int_t entries_EW_FB_Map_LS[ani_values][648];
    Int_t entries_Full_Map_LS[ani_values][648];
    
    Int_t theta_binHistos_LS;
    Int_t phi_binHistos_LS;
    
    ULong64_t events_LS;
    
    Double_t chi2_HS[7];
    Double_t ndf_HS[7];
    Double_t chi2_r_HS[7];
    Double_t fit_par_HS[7][3];
    Double_t fit_err_HS[7][3];
    Double_t delta_HS[7];
    Double_t sum_par_HS[7];
    
    Double_t CMatrix_NS_HS[3][3];
    Double_t CMatrix_EW_HS[3][3];
    Double_t CMatrix_FB_HS[3][3];
    Double_t CMatrix_NS_EW_HS[3][3];
    Double_t CMatrix_NS_FB_HS[3][3];
    Double_t CMatrix_EW_FB_HS[3][3];
    Double_t CMatrix_Full_HS[3][3];
    
    Int_t entries_NS_Map_HS[ani_values][648];
    Int_t entries_EW_Map_HS[ani_values][648];
    Int_t entries_FB_Map_HS[ani_values][648];
    Int_t entries_NS_EW_Map_HS[ani_values][648];
    Int_t entries_NS_FB_Map_HS[ani_values][648];
    Int_t entries_EW_FB_Map_HS[ani_values][648];
    Int_t entries_Full_Map_HS[ani_values][648];
    
    Int_t theta_binHistos_HS;
    Int_t phi_binHistos_HS;
    
    ULong64_t events_HS;
    
    
    relative_fitResult()
    {
        
        for(Int_t idx = 0; idx < 7; ++idx)
        {
            chi2_LS[idx] = -1;
            ndf_LS[idx] = -1;
            chi2_r_LS[idx] = -1;
            delta_LS[idx] = -1;
            sum_par_LS[idx] = -999;
            
            chi2_HS[idx] = -1;
            ndf_HS[idx] = -1;
            chi2_r_HS[idx] = -1;
            delta_HS[idx] = -1;
            sum_par_HS[idx] = -999;
            
            for(Int_t k = 0; k < 3; ++k)
            {
                fit_par_LS[idx][k] = -1;
                fit_err_LS[idx][k] = -1;
                
                fit_par_HS[idx][k] = -1;
                fit_err_HS[idx][k] = -1;
                
                for(int j = 0; j < 3; ++j)
                {
                    CMatrix_NS_LS[k][j] = -1;
                    CMatrix_EW_LS[k][j] = -1;
                    CMatrix_FB_LS[k][j] = -1;
                    CMatrix_NS_EW_LS[k][j] = -1;
                    CMatrix_NS_FB_LS[k][j] = -1;
                    CMatrix_EW_FB_LS[k][j] = -1;
                    CMatrix_Full_LS[k][j] = -1;
                    
                    CMatrix_NS_HS[k][j] = -1;
                    CMatrix_EW_HS[k][j] = -1;
                    CMatrix_FB_HS[k][j] = -1;
                    CMatrix_NS_EW_HS[k][j] = -1;
                    CMatrix_NS_FB_HS[k][j] = -1;
                    CMatrix_EW_FB_HS[k][j] = -1;
                    CMatrix_Full_HS[k][j] = -1;
                    
                }
            }
        }
        
        theta_binHistos_LS = 0;
        phi_binHistos_LS = 0;
        theta_binHistos_HS = 0;
        phi_binHistos_HS = 0;
        
        for(Int_t idx = 0; idx < 3; ++idx)
            inputAni[idx] = -1;
        
        for(Int_t a_idx=0; a_idx < ani_values; ++a_idx)
            for(Int_t idx = 0; idx < 648; ++idx)
            {
                
                entries_NS_Map_LS[a_idx][idx] = -999;
                entries_EW_Map_LS[a_idx][idx] = -999;
                entries_FB_Map_LS[a_idx][idx] = -999;
                entries_NS_EW_Map_LS[a_idx][idx] = -999;
                entries_NS_FB_Map_LS[a_idx][idx] = -999;
                entries_EW_FB_Map_LS[a_idx][idx] = -999;
                entries_Full_Map_LS[a_idx][idx] = -999;
                
                entries_NS_Map_HS[a_idx][idx] = -999;
                entries_EW_Map_HS[a_idx][idx] = -999;
                entries_FB_Map_HS[a_idx][idx] = -999;
                entries_NS_EW_Map_HS[a_idx][idx] = -999;
                entries_NS_FB_Map_HS[a_idx][idx] = -999;
                entries_EW_FB_Map_HS[a_idx][idx] = -999;
                entries_Full_Map_HS[a_idx][idx] = -999;
                
            }
        
        events_LS = 0;
        events_HS = 0;
        
        seed = 0;
        seed_list_line = 0;
        
    }
    
    ~relative_fitResult() { }
    
};


void link_tree_branches(TChain &tree,relative_fitResult &tmp_fit,std::string tree_dir_path);
void create_tree_branches(TTree &fiTree,relative_fitResult &tmp_fit);

void adder(std::string tree_dir_path,std::string full_tree_out_path)
{
    TChain tree("fiTree");
    TTree myFinalTree("fiTree","TemplateFit (full) results TTree");
    relative_fitResult tmp_fit;
    
    full_tree_out_path.append("full_MapsFitTree.root");
    
    link_tree_branches(tree,tmp_fit,tree_dir_path);
    create_tree_branches(myFinalTree,tmp_fit);
    
    for(UInt_t tree_idx=0; tree_idx<tree.GetEntries(); ++tree_idx)
    {
        tree.GetEntry(tree_idx);
        myFinalTree.Fill();
    }
 
    TFile outTree(full_tree_out_path.c_str(),"RECREATE");
    if(outTree.IsZombie())
    {
        std::cerr << "\n\n Error writing final TTree \n\n";
        exit(100);
    }
    
    myFinalTree.Write();
    
    outTree.Write();
    outTree.Close();
    
}

void link_tree_branches(TChain &tree,relative_fitResult &tmp_fit,std::string tree_dir_path)
{
    tree.Add(Form("%s/*_tree.root", tree_dir_path.c_str()));
    
    tree.SetBranchAddress("chi2_LS",tmp_fit.chi2_LS);
    tree.SetBranchAddress("chi2_r_LS",tmp_fit.chi2_r_LS);
    tree.SetBranchAddress("ndf_LS",tmp_fit.ndf_LS);
    tree.SetBranchAddress("fit_par_LS",tmp_fit.fit_par_LS);
    tree.SetBranchAddress("fit_err_LS",tmp_fit.fit_err_LS);
    tree.SetBranchAddress("delta_LS",tmp_fit.delta_LS);
    tree.SetBranchAddress("sum_par_LS",tmp_fit.sum_par_LS);
    
    tree.SetBranchAddress("CMatrix_NS_LS",tmp_fit.CMatrix_NS_LS);
    tree.SetBranchAddress("CMatrix_EW_LS",tmp_fit.CMatrix_EW_LS);
    tree.SetBranchAddress("CMatrix_FB_LS",tmp_fit.CMatrix_FB_LS);
    tree.SetBranchAddress("CMatrix_NS_EW_LS",tmp_fit.CMatrix_NS_EW_LS);
    tree.SetBranchAddress("CMatrix_NS_FB_LS",tmp_fit.CMatrix_NS_FB_LS);
    tree.SetBranchAddress("CMatrix_EW_FB_LS",tmp_fit.CMatrix_EW_FB_LS);
    tree.SetBranchAddress("CMatrix_Full_LS",tmp_fit.CMatrix_Full_LS);
    
    tree.SetBranchAddress("entries_NS_Map_LS",tmp_fit.entries_NS_Map_LS);
    tree.SetBranchAddress("entries_EW_Map_LS",tmp_fit.entries_EW_Map_LS);
    tree.SetBranchAddress("entries_FB_Map_LS",tmp_fit.entries_FB_Map_LS);
    tree.SetBranchAddress("entries_NS_EW_Map_LS",tmp_fit.entries_NS_EW_Map_LS);
    tree.SetBranchAddress("entries_NS_FB_Map_LS",tmp_fit.entries_NS_FB_Map_LS);
    tree.SetBranchAddress("entries_EW_FB_Map_LS",tmp_fit.entries_EW_FB_Map_LS);
    tree.SetBranchAddress("entries_Full_Map_LS",tmp_fit.entries_Full_Map_LS);
    
    tree.SetBranchAddress("theta_binHistos_LS",&tmp_fit.theta_binHistos_LS);
    tree.SetBranchAddress("phi_binhistos_LS",&tmp_fit.phi_binHistos_LS);
    tree.SetBranchAddress("events_LS",&tmp_fit.events_LS);
    
    tree.SetBranchAddress("chi2_HS",tmp_fit.chi2_HS);
    tree.SetBranchAddress("chi2_r_HS",tmp_fit.chi2_r_HS);
    tree.SetBranchAddress("ndf_HS",tmp_fit.ndf_HS);
    tree.SetBranchAddress("fit_par_HS",tmp_fit.fit_par_HS);
    tree.SetBranchAddress("fit_err_HS",tmp_fit.fit_err_HS);
    tree.SetBranchAddress("delta_HS",tmp_fit.delta_HS);
    tree.SetBranchAddress("sum_par_HS",tmp_fit.sum_par_HS);
    
    tree.SetBranchAddress("CMatrix_NS_HS",tmp_fit.CMatrix_NS_HS);
    tree.SetBranchAddress("CMatrix_EW_HS",tmp_fit.CMatrix_EW_HS);
    tree.SetBranchAddress("CMatrix_FB_HS",tmp_fit.CMatrix_FB_HS);
    tree.SetBranchAddress("CMatrix_NS_EW_HS",tmp_fit.CMatrix_NS_EW_HS);
    tree.SetBranchAddress("CMatrix_NS_FB_HS",tmp_fit.CMatrix_NS_FB_HS);
    tree.SetBranchAddress("CMatrix_EW_FB_HS",tmp_fit.CMatrix_EW_FB_HS);
    tree.SetBranchAddress("CMatrix_Full_HS",tmp_fit.CMatrix_Full_HS);
    
    tree.SetBranchAddress("entries_NS_Map_HS",tmp_fit.entries_NS_Map_HS);
    tree.SetBranchAddress("entries_EW_Map_HS",tmp_fit.entries_EW_Map_HS);
    tree.SetBranchAddress("entries_FB_Map_HS",tmp_fit.entries_FB_Map_HS);
    tree.SetBranchAddress("entries_NS_EW_Map_HS",tmp_fit.entries_NS_EW_Map_HS);
    tree.SetBranchAddress("entries_NS_FB_Map_HS",tmp_fit.entries_NS_FB_Map_HS);
    tree.SetBranchAddress("entries_EW_FB_Map_HS",tmp_fit.entries_EW_FB_Map_HS);
    tree.SetBranchAddress("entries_Full_Map_HS",tmp_fit.entries_Full_Map_HS);
    
    tree.SetBranchAddress("theta_binHistos_HS",&tmp_fit.theta_binHistos_HS);
    tree.SetBranchAddress("phi_binhistos_HS",&tmp_fit.phi_binHistos_HS);
    tree.SetBranchAddress("events_HS",&tmp_fit.events_HS);
    
    tree.SetBranchAddress("inputAni",tmp_fit.inputAni);
    
    tree.SetBranchAddress("seed",&tmp_fit.seed);
    tree.SetBranchAddress("seed_list_line",&tmp_fit.seed_list_line);
    
}

void create_tree_branches(TTree &fiTree,relative_fitResult &tmp_fit)
{
    
    fiTree.Branch("chi2_LS",tmp_fit.chi2_LS,"chi2_LS[7]/D");
    fiTree.Branch("chi2_r_LS",tmp_fit.chi2_r_LS,"chi2_r_LS[7]/D");
    fiTree.Branch("ndf_LS",tmp_fit.ndf_LS,"ndf_LS[7]/D");
    fiTree.Branch("fit_par_LS",tmp_fit.fit_par_LS,"fit_par_LS[7][3]/D");
    fiTree.Branch("fit_err_LS",tmp_fit.fit_err_LS,"fit_err_LS[7][3]/D");
    fiTree.Branch("delta_LS",tmp_fit.delta_LS,"delta_LS[7]/D");
    fiTree.Branch("sum_par_LS",tmp_fit.sum_par_LS,"sum_par_LS[7]/D");
    
    fiTree.Branch("CMatrix_NS_LS",tmp_fit.CMatrix_NS_LS,"CMatrix_NS_LS[3][3]/D");
    fiTree.Branch("CMatrix_EW_LS",tmp_fit.CMatrix_EW_LS,"CMatrix_EW_LS[3][3]/D");
    fiTree.Branch("CMatrix_FB_LS",tmp_fit.CMatrix_FB_LS,"CMatrix_FB_LS[3][3]/D");
    fiTree.Branch("CMatrix_NS_EW_LS",tmp_fit.CMatrix_NS_EW_LS,"CMatrix_NS_EW_LS[3][3]/D");
    fiTree.Branch("CMatrix_NS_FB_LS",tmp_fit.CMatrix_NS_FB_LS,"CMatrix_NS_FB_LS[3][3]/D");
    fiTree.Branch("CMatrix_EW_FB_LS",tmp_fit.CMatrix_EW_FB_LS,"CMatrix_EW_FB_LS[3][3]/D");
    fiTree.Branch("CMatrix_Full_LS",tmp_fit.CMatrix_Full_LS,"CMatrix_Full_LS[3][3]/D");
    
    fiTree.Branch("entries_NS_Map_LS",tmp_fit.entries_NS_Map_LS,"entries_NS_Map_LS[648]/I");
    fiTree.Branch("entries_EW_Map_LS",tmp_fit.entries_EW_Map_LS,"entries_EW_Map_LS[648]/I");
    fiTree.Branch("entries_FB_Map_LS",tmp_fit.entries_FB_Map_LS,"entries_FB_Map_LS[648]/I");
    fiTree.Branch("entries_NS_EW_Map_LS",tmp_fit.entries_NS_EW_Map_LS,"entries_NS_EW_Map_LS[648]/I");
    fiTree.Branch("entries_NS_FB_Map_LS",tmp_fit.entries_NS_FB_Map_LS,"entries_NS_FB_Map_LS[648]/I");
    fiTree.Branch("entries_EW_FB_Map_LS",tmp_fit.entries_EW_FB_Map_LS,"entries_EW_FB_Map_LS[648]/I");
    fiTree.Branch("entries_Full_Map_LS",tmp_fit.entries_Full_Map_LS,"entries_Full_Map_LS[648]/I");
    
    fiTree.Branch("theta_binHistos_LS",&tmp_fit.theta_binHistos_LS,"theta_binHistos_LS/I");
    fiTree.Branch("phi_binhistos_LS",&tmp_fit.phi_binHistos_LS,"phi_binHistos_LS/I");
    
    fiTree.Branch("seed",&tmp_fit.seed,"seed/l");
    fiTree.Branch("seed_list_line",&tmp_fit.seed_list_line,"seed_list_line/l");
    
    fiTree.Branch("events_LS",&tmp_fit.events_LS,"events_LS/l");
    
    fiTree.Branch("chi2_HS",tmp_fit.chi2_HS,"chi2_HS[7]/D");
    fiTree.Branch("chi2_r_HS",tmp_fit.chi2_r_HS,"chi2_r_HS[7]/D");
    fiTree.Branch("ndf_HS",tmp_fit.ndf_HS,"ndf_HS[7]/D");
    fiTree.Branch("fit_par_HS",tmp_fit.fit_par_HS,"fit_par_HS[7][3]/D");
    fiTree.Branch("fit_err_HS",tmp_fit.fit_err_HS,"fit_err_HS[7][3]/D");
    fiTree.Branch("delta_HS",tmp_fit.delta_HS,"delta_HS[7]/D");
    fiTree.Branch("sum_par_HS",tmp_fit.sum_par_HS,"sum_par_HS[7]/D");
    
    fiTree.Branch("CMatrix_NS_HS",tmp_fit.CMatrix_NS_HS,"CMatrix_NS_HS[3][3]/D");
    fiTree.Branch("CMatrix_EW_HS",tmp_fit.CMatrix_EW_HS,"CMatrix_EW_HS[3][3]/D");
    fiTree.Branch("CMatrix_FB_HS",tmp_fit.CMatrix_FB_HS,"CMatrix_FB_HS[3][3]/D");
    fiTree.Branch("CMatrix_NS_EW_HS",tmp_fit.CMatrix_NS_EW_HS,"CMatrix_NS_EW_HS[3][3]/D");
    fiTree.Branch("CMatrix_NS_FB_HS",tmp_fit.CMatrix_NS_FB_HS,"CMatrix_NS_FB_HS[3][3]/D");
    fiTree.Branch("CMatrix_EW_FB_HS",tmp_fit.CMatrix_EW_FB_HS,"CMatrix_EW_FB_HS[3][3]/D");
    fiTree.Branch("CMatrix_Full_HS",tmp_fit.CMatrix_Full_HS,"CMatrix_Full_HS[3][3]/D");
    
    fiTree.Branch("entries_NS_Map_HS",tmp_fit.entries_NS_Map_HS,"entries_NS_Map_HS[648]/I");
    fiTree.Branch("entries_EW_Map_HS",tmp_fit.entries_EW_Map_HS,"entries_EW_Map_HS[648]/I");
    fiTree.Branch("entries_FB_Map_HS",tmp_fit.entries_FB_Map_HS,"entries_FB_Map_HS[648]/I");
    fiTree.Branch("entries_NS_EW_Map_HS",tmp_fit.entries_NS_EW_Map_HS,"entries_NS_EW_Map_HS[648]/I");
    fiTree.Branch("entries_NS_FB_Map_HS",tmp_fit.entries_NS_FB_Map_HS,"entries_NS_FB_Map_HS[648]/I");
    fiTree.Branch("entries_EW_FB_Map_HS",tmp_fit.entries_EW_FB_Map_HS,"entries_EW_FB_Map_HS[648]/I");
    fiTree.Branch("entries_Full_Map_HS",tmp_fit.entries_Full_Map_HS,"entries_Full_Map_HS[648]/I");
    
    fiTree.Branch("theta_binHistos_HS",&tmp_fit.theta_binHistos_HS,"theta_binHistos_HS/I");
    fiTree.Branch("phi_binhistos_HS",&tmp_fit.phi_binHistos_HS,"phi_binHistos_HS/I");
    
    fiTree.Branch("events_HS",&tmp_fit.events_HS,"events_HS/l");
    
    fiTree.Branch("inputAni",tmp_fit.inputAni,"inputAni[3]/D");
    
}

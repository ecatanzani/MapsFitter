#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"

#define ani_values 3

class fitResult
{
public:
    
    Double_t inputAni[3];
    ULong64_t seed;
    ULong64_t seed_list_line;
    
    Double_t chi2_LS[8];
    Double_t ndf_LS[8];
    Double_t chi2_r_LS[8];
    Double_t fit_par_LS[8][4];
    Double_t fit_err_LS[8][4];
    Double_t delta_LS[8];
    Double_t sum_par_LS[8];
    
    
    Double_t CMatrix_Iso_LS[4][4];
    Double_t CMatrix_NS_LS[4][4];
    Double_t CMatrix_EW_LS[4][4];
    Double_t CMatrix_FB_LS[4][4];
    Double_t CMatrix_NS_EW_LS[4][4];
    Double_t CMatrix_NS_FB_LS[4][4];
    Double_t CMatrix_EW_FB_LS[4][4];
    Double_t CMatrix_Full_LS[4][4];
    
    Int_t theta_binHistos_LS;
    Int_t phi_binHistos_LS;
    
    Int_t entries_Iso_Map_LS[ani_values][648];
    Int_t entries_NS_Map_LS[ani_values][648];
    Int_t entries_EW_Map_LS[ani_values][648];
    Int_t entries_FB_Map_LS[ani_values][648];
    Int_t entries_NS_EW_Map_LS[ani_values][648];
    Int_t entries_NS_FB_Map_LS[ani_values][648];
    Int_t entries_EW_FB_Map_LS[ani_values][648];
    Int_t entries_Full_Map_LS[ani_values][648];
    
    ULong64_t events_LS;
    
    Double_t chi2_HS[8];
    Double_t ndf_HS[8];
    Double_t chi2_r_HS[8];
    Double_t fit_par_HS[8][4];
    Double_t fit_err_HS[8][4];
    Double_t delta_HS[8];
    Double_t sum_par_HS[8];
    
    Double_t CMatrix_Iso_HS[4][4];
    Double_t CMatrix_NS_HS[4][4];
    Double_t CMatrix_EW_HS[4][4];
    Double_t CMatrix_FB_HS[4][4];
    Double_t CMatrix_NS_EW_HS[4][4];
    Double_t CMatrix_NS_FB_HS[4][4];
    Double_t CMatrix_EW_FB_HS[4][4];
    Double_t CMatrix_Full_HS[4][4];
    
    Int_t theta_binHistos_HS;
    Int_t phi_binHistos_HS;
    
    Int_t entries_Iso_Map_HS[ani_values][648];
    Int_t entries_NS_Map_HS[ani_values][648];
    Int_t entries_EW_Map_HS[ani_values][648];
    Int_t entries_FB_Map_HS[ani_values][648];
    Int_t entries_NS_EW_Map_HS[ani_values][648];
    Int_t entries_NS_FB_Map_HS[ani_values][648];
    Int_t entries_EW_FB_Map_HS[ani_values][648];
    Int_t entries_Full_Map_HS[ani_values][648];
    
    ULong64_t events_HS;
    
    
    fitResult()
    {
        
        for(Int_t idx = 0; idx < 8; ++idx)
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
            
            for(Int_t k = 0; k < 4; ++k)
            {
                fit_par_LS[idx][k] = -1;
                fit_err_LS[idx][k] = -1;
                
                fit_par_HS[idx][k] = -1;
                fit_err_HS[idx][k] = -1;
                
                for(int j = 0; j < 4; ++j)
                {
                    CMatrix_Iso_LS[k][j] = -1;
                    CMatrix_NS_LS[k][j] = -1;
                    CMatrix_EW_LS[k][j] = -1;
                    CMatrix_FB_LS[k][j] = -1;
                    CMatrix_NS_EW_LS[k][j] = -1;
                    CMatrix_NS_FB_LS[k][j] = -1;
                    CMatrix_EW_FB_LS[k][j] = -1;
                    CMatrix_Full_LS[k][j] = -1;
                    
                    CMatrix_Iso_HS[k][j] = -1;
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
                
                entries_Iso_Map_LS[a_idx][idx] = -999;
                entries_NS_Map_LS[a_idx][idx] = -999;
                entries_EW_Map_LS[a_idx][idx] = -999;
                entries_FB_Map_LS[a_idx][idx] = -999;
                entries_NS_EW_Map_LS[a_idx][idx] = -999;
                entries_NS_FB_Map_LS[a_idx][idx] = -999;
                entries_EW_FB_Map_LS[a_idx][idx] = -999;
                entries_Full_Map_LS[a_idx][idx] = -999;
                
                entries_Iso_Map_HS[a_idx][idx] = -999;
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
    
    ~fitResult() { }
    
};

void tree_linking(TTree* myTree,fitResult &tmp_result,UInt_t &tree_entries);

void print_data_distribution(std::string input_tree,std::string output_file,UInt_t tree_event,Int_t ani_idx)
{
    
    UInt_t tree_entries = 0;
    
    ////////////////////// Data event distribution
    
    TH1D entries_Iso_LS("entries_Iso_LS","events distribution LS (Isotropic Map)",10,0,10);
    TH1D entries_NS_LS("entries_NS_Map_LS","events distribution LS (NS Map)",10,0,10);
    TH1D entries_EW_LS("entries_EW_Map_LS","events distribution LS (EW Map)",10,0,10);
    TH1D entries_FB_LS("entries_FB_Map_LS","events distribution LS (FB Map)",10,0,10);
    TH1D entries_NS_EW_LS("entries_NS_EW_Map_LS","events distribution LS (NS-EW Map)",10,0,10);
    TH1D entries_NS_FB_LS("entries_NS_FB_Map_LS","events distribution LS (NS-FB Map)",10,0,10);
    TH1D entries_EW_FB_LS("entries_EW_FB_Map_LS","events distribution LS (EW-FB Map)",10,0,10);
    TH1D entries_Full_LS("entries_Full_Map_LS","events distribution LS (Full Map)",10,0,10);
    
    TH1D entries_Iso_HS("entries_Iso_HS","events distribution HS (Isotropic Map)",100,2000,7000);
    TH1D entries_NS_HS("entries_NS_Map_HS","events distribution HS (NS Map)",100,2000,7000);
    TH1D entries_EW_HS("entries_EW_Map_HS","events distribution HS (EW Map)",100,2000,7000);
    TH1D entries_FB_HS("entries_FB_Map_HS","events distribution HS (FB Map)",100,2000,7000);
    TH1D entries_NS_EW_HS("entries_NS_EW_Map_HS","events distribution HS (NS-EW Map)",100,2000,7000);
    TH1D entries_NS_FB_HS("entries_NS_FB_Map_HS","events distribution HS (NS-FB Map)",100,2000,7000);
    TH1D entries_EW_FB_HS("entries_EW_FB_Map_HS","events distribution HS (EW-FB Map)",100,2000,7000);
    TH1D entries_Full_HS("entries_Full_Map_HS","events distribution HS (Full Map)",100,2000,7000);
    
    ///////////////////////////////////////////////////
    
    TFile inFile(input_tree.c_str(),"READ");
    if(inFile.IsZombie())
    {
        std::cout << "\n\nError opening TTree ROOT file \n\n";
        exit(100);
    }
    
    TTree* myTree = (TTree*)inFile.Get("fiTree");
    fitResult tmp_result;
    
    tree_linking(myTree,tmp_result,tree_entries);
    
    ///////////////////////////////////////////////////
    
    myTree->GetEntry(tree_event);
    
    for(Int_t idx=0; idx<648; ++idx)
    {
        entries_Iso_LS.Fill(tmp_result.entries_Iso_Map_LS[ani_idx][idx]);
        entries_NS_LS.Fill(tmp_result.entries_NS_Map_LS[ani_idx][idx]);
        entries_EW_LS.Fill(tmp_result.entries_EW_Map_LS[ani_idx][idx]);
        entries_FB_LS.Fill(tmp_result.entries_FB_Map_LS[ani_idx][idx]);
        entries_NS_EW_LS.Fill(tmp_result.entries_NS_EW_Map_LS[ani_idx][idx]);
        entries_NS_FB_LS.Fill(tmp_result.entries_NS_FB_Map_LS[ani_idx][idx]);
        entries_EW_FB_LS.Fill(tmp_result.entries_EW_FB_Map_LS[ani_idx][idx]);
        entries_Full_LS.Fill(tmp_result.entries_Full_Map_LS[ani_idx][idx]);
        
        entries_Iso_HS.Fill(tmp_result.entries_Iso_Map_HS[ani_idx][idx]);
        entries_NS_HS.Fill(tmp_result.entries_NS_Map_HS[ani_idx][idx]);
        entries_EW_HS.Fill(tmp_result.entries_EW_Map_HS[ani_idx][idx]);
        entries_FB_HS.Fill(tmp_result.entries_FB_Map_HS[ani_idx][idx]);
        entries_NS_EW_HS.Fill(tmp_result.entries_NS_EW_Map_HS[ani_idx][idx]);
        entries_NS_FB_HS.Fill(tmp_result.entries_NS_FB_Map_HS[ani_idx][idx]);
        entries_EW_FB_HS.Fill(tmp_result.entries_EW_FB_Map_HS[ani_idx][idx]);
        entries_Full_HS.Fill(tmp_result.entries_Full_Map_HS[ani_idx][idx]);
    }
    
    TFile outFile(output_file.c_str(),"RECREATE");
    if(outFile.IsZombie())
    {
        std::cout << "\n\nError opening TTree ROOT out file \n\n";
        exit(100);
    }
    
    entries_Iso_LS.Write();
    entries_NS_LS.Write();
    entries_EW_LS.Write();
    entries_FB_LS.Write();
    entries_NS_EW_LS.Write();
    entries_NS_FB_LS.Write();
    entries_EW_FB_LS.Write();
    entries_Full_LS.Write();
    
    entries_Iso_HS.Write();
    entries_NS_HS.Write();
    entries_EW_HS.Write();
    entries_FB_HS.Write();
    entries_NS_EW_HS.Write();
    entries_NS_FB_HS.Write();
    entries_EW_FB_HS.Write();
    entries_Full_HS.Write();
    
    outFile.Write();
    outFile.Close();
    
}

void tree_linking(TTree* myTree,fitResult &tmp_result,UInt_t &tree_entries)
{
    myTree->SetBranchAddress("chi2_LS",tmp_result.chi2_LS);
    myTree->SetBranchAddress("chi2_r_LS",tmp_result.chi2_r_LS);
    myTree->SetBranchAddress("ndf_LS",tmp_result.ndf_LS);
    myTree->SetBranchAddress("fit_par_LS",tmp_result.fit_par_LS);
    myTree->SetBranchAddress("fit_err_LS",tmp_result.fit_err_LS);
    myTree->SetBranchAddress("delta_LS",tmp_result.delta_LS);
    myTree->SetBranchAddress("sum_par_LS",tmp_result.sum_par_LS);
    
    myTree->SetBranchAddress("CMatrix_Iso_LS",tmp_result.CMatrix_Iso_LS);
    myTree->SetBranchAddress("CMatrix_NS_LS",tmp_result.CMatrix_NS_LS);
    myTree->SetBranchAddress("CMatrix_EW_LS",tmp_result.CMatrix_EW_LS);
    myTree->SetBranchAddress("CMatrix_FB_LS",tmp_result.CMatrix_FB_LS);
    myTree->SetBranchAddress("CMatrix_NS_EW_LS",tmp_result.CMatrix_NS_EW_LS);
    myTree->SetBranchAddress("CMatrix_NS_FB_LS",tmp_result.CMatrix_NS_FB_LS);
    myTree->SetBranchAddress("CMatrix_EW_FB_LS",tmp_result.CMatrix_EW_FB_LS);
    myTree->SetBranchAddress("CMatrix_Full_LS",tmp_result.CMatrix_Full_LS);
    
    myTree->SetBranchAddress("entries_Iso_Map_LS",tmp_result.entries_Iso_Map_LS);
    myTree->SetBranchAddress("entries_NS_Map_LS",tmp_result.entries_NS_Map_LS);
    myTree->SetBranchAddress("entries_EW_Map_LS",tmp_result.entries_EW_Map_LS);
    myTree->SetBranchAddress("entries_FB_Map_LS",tmp_result.entries_FB_Map_LS);
    myTree->SetBranchAddress("entries_NS_EW_Map_LS",tmp_result.entries_NS_EW_Map_LS);
    myTree->SetBranchAddress("entries_NS_FB_Map_LS",tmp_result.entries_NS_FB_Map_LS);
    myTree->SetBranchAddress("entries_EW_FB_Map_LS",tmp_result.entries_EW_FB_Map_LS);
    myTree->SetBranchAddress("entries_Full_Map_LS",tmp_result.entries_Full_Map_LS);
    
    myTree->SetBranchAddress("theta_binHistos_LS",&tmp_result.theta_binHistos_LS);
    myTree->SetBranchAddress("phi_binhistos_LS",&tmp_result.phi_binHistos_LS);
    myTree->SetBranchAddress("events_LS",&tmp_result.events_LS);
    
    myTree->SetBranchAddress("chi2_HS",tmp_result.chi2_HS);
    myTree->SetBranchAddress("chi2_r_HS",tmp_result.chi2_r_HS);
    myTree->SetBranchAddress("ndf_HS",tmp_result.ndf_HS);
    myTree->SetBranchAddress("fit_par_HS",tmp_result.fit_par_HS);
    myTree->SetBranchAddress("fit_err_HS",tmp_result.fit_err_HS);
    myTree->SetBranchAddress("delta_HS",tmp_result.delta_HS);
    myTree->SetBranchAddress("sum_par_HS",tmp_result.sum_par_HS);
    
    myTree->SetBranchAddress("CMatrix_Iso_HS",tmp_result.CMatrix_Iso_HS);
    myTree->SetBranchAddress("CMatrix_NS_HS",tmp_result.CMatrix_NS_HS);
    myTree->SetBranchAddress("CMatrix_EW_HS",tmp_result.CMatrix_EW_HS);
    myTree->SetBranchAddress("CMatrix_FB_HS",tmp_result.CMatrix_FB_HS);
    myTree->SetBranchAddress("CMatrix_NS_EW_HS",tmp_result.CMatrix_NS_EW_HS);
    myTree->SetBranchAddress("CMatrix_NS_FB_HS",tmp_result.CMatrix_NS_FB_HS);
    myTree->SetBranchAddress("CMatrix_EW_FB_HS",tmp_result.CMatrix_EW_FB_HS);
    myTree->SetBranchAddress("CMatrix_Full_HS",tmp_result.CMatrix_Full_HS);
    
    myTree->SetBranchAddress("entries_Iso_Map_HS",tmp_result.entries_Iso_Map_HS);
    myTree->SetBranchAddress("entries_NS_Map_HS",tmp_result.entries_NS_Map_HS);
    myTree->SetBranchAddress("entries_EW_Map_HS",tmp_result.entries_EW_Map_HS);
    myTree->SetBranchAddress("entries_FB_Map_HS",tmp_result.entries_FB_Map_HS);
    myTree->SetBranchAddress("entries_NS_EW_Map_HS",tmp_result.entries_NS_EW_Map_HS);
    myTree->SetBranchAddress("entries_NS_FB_Map_HS",tmp_result.entries_NS_FB_Map_HS);
    myTree->SetBranchAddress("entries_EW_FB_Map_HS",tmp_result.entries_EW_FB_Map_HS);
    myTree->SetBranchAddress("entries_Full_Map_HS",tmp_result.entries_Full_Map_HS);
    
    myTree->SetBranchAddress("theta_binHistos_HS",&tmp_result.theta_binHistos_HS);
    myTree->SetBranchAddress("phi_binhistos_HS",&tmp_result.phi_binHistos_HS);
    myTree->SetBranchAddress("events_HS",&tmp_result.events_HS);
    
    myTree->SetBranchAddress("inputAni",tmp_result.inputAni);
    
    myTree->SetBranchAddress("seed",&tmp_result.seed);
    myTree->SetBranchAddress("seed_list_line",&tmp_result.seed_list_line);
    
    tree_entries = myTree->GetEntries();
    
}

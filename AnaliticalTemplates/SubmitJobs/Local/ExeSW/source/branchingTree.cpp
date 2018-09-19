
#include "MyHead.h"

void create_tree_branches(TTree &fiTree,fitResult &tmp_fit)
{
    
    fiTree.Branch("chi2_LS",tmp_fit.chi2_LS,"chi2_LS[8]/D");
    fiTree.Branch("chi2_r_LS",tmp_fit.chi2_r_LS,"chi2_r_LS[8]/D");
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
    
    fiTree.Branch("entries_Iso_Map_LS",tmp_fit.entries_Iso_Map_LS,"entries_Iso_Map_LS[3][648]/I");
    fiTree.Branch("entries_NS_Map_LS",tmp_fit.entries_NS_Map_LS,"entries_NS_Map_LS[3][648]/I");
    fiTree.Branch("entries_EW_Map_LS",tmp_fit.entries_EW_Map_LS,"entries_EW_Map_LS[3][648]/I");
    fiTree.Branch("entries_FB_Map_LS",tmp_fit.entries_FB_Map_LS,"entries_FB_Map_LS[3][648]/I");
    fiTree.Branch("entries_NS_EW_Map_LS",tmp_fit.entries_NS_EW_Map_LS,"entries_NS_EW_Map_LS[3][648]/I");
    fiTree.Branch("entries_NS_FB_Map_LS",tmp_fit.entries_NS_FB_Map_LS,"entries_NS_FB_Map_LS[3][648]/I");
    fiTree.Branch("entries_EW_FB_Map_LS",tmp_fit.entries_EW_FB_Map_LS,"entries_EW_FB_Map_LS[3][648]/I");
    fiTree.Branch("entries_Full_Map_LS",tmp_fit.entries_Full_Map_LS,"entries_Full_Map_LS[3][648]/I");
    
    fiTree.Branch("theta_binHistos_LS",&tmp_fit.theta_binHistos_LS,"theta_binHistos_LS/I");
    fiTree.Branch("phi_binhistos_LS",&tmp_fit.phi_binHistos_LS,"phi_binHistos_LS/I");
    
    fiTree.Branch("seed",&tmp_fit.seed,"seed/l");
    fiTree.Branch("seed_list_line",&tmp_fit.seed_list_line,"seed_list_line/l");
    
    fiTree.Branch("events_LS",&tmp_fit.events_LS,"events_LS/l");
    
    fiTree.Branch("chi2_HS",tmp_fit.chi2_HS,"chi2_HS[8]/D");
    fiTree.Branch("chi2_r_HS",tmp_fit.chi2_r_HS,"chi2_r_HS[8]/D");
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
    
    fiTree.Branch("entries_Iso_Map_HS",tmp_fit.entries_Iso_Map_HS,"entries_Iso_Map_HS[3][648]/I");
    fiTree.Branch("entries_NS_Map_HS",tmp_fit.entries_NS_Map_HS,"entries_NS_Map_HS[3][648]/I");
    fiTree.Branch("entries_EW_Map_HS",tmp_fit.entries_EW_Map_HS,"entries_EW_Map_HS[3][648]/I");
    fiTree.Branch("entries_FB_Map_HS",tmp_fit.entries_FB_Map_HS,"entries_FB_Map_HS[3][648]/I");
    fiTree.Branch("entries_NS_EW_Map_HS",tmp_fit.entries_NS_EW_Map_HS,"entries_NS_EW_Map_HS[3][648]/I");
    fiTree.Branch("entries_NS_FB_Map_HS",tmp_fit.entries_NS_FB_Map_HS,"entries_NS_FB_Map_HS[3][648]/I");
    fiTree.Branch("entries_EW_FB_Map_HS",tmp_fit.entries_EW_FB_Map_HS,"entries_EW_FB_Map_HS[3][648]/I");
    fiTree.Branch("entries_Full_Map_HS",tmp_fit.entries_Full_Map_HS,"entries_Full_Map_HS[3][648]/I");
    
    fiTree.Branch("theta_binHistos_HS",&tmp_fit.theta_binHistos_HS,"theta_binHistos_HS/I");
    fiTree.Branch("phi_binhistos_HS",&tmp_fit.phi_binHistos_HS,"phi_binHistos_HS/I");
    
    fiTree.Branch("events_HS",&tmp_fit.events_HS,"events_HS/l");
    
    fiTree.Branch("inputAni",tmp_fit.inputAni,"inputAni[3]/D");
}


void create_rel_tree_branches(TTree &fiTree,relative_fitResult &tmp_fit)
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
    
    fiTree.Branch("entries_NS_Map_LS",tmp_fit.entries_NS_Map_LS,"entries_NS_Map_LS[3][648]/I");
    fiTree.Branch("entries_EW_Map_LS",tmp_fit.entries_EW_Map_LS,"entries_EW_Map_LS[3][648]/I");
    fiTree.Branch("entries_FB_Map_LS",tmp_fit.entries_FB_Map_LS,"entries_FB_Map_LS[3][648]/I");
    fiTree.Branch("entries_NS_EW_Map_LS",tmp_fit.entries_NS_EW_Map_LS,"entries_NS_EW_Map_LS[3][648]/I");
    fiTree.Branch("entries_NS_FB_Map_LS",tmp_fit.entries_NS_FB_Map_LS,"entries_NS_FB_Map_LS[3][648]/I");
    fiTree.Branch("entries_EW_FB_Map_LS",tmp_fit.entries_EW_FB_Map_LS,"entries_EW_FB_Map_LS[3][648]/I");
    fiTree.Branch("entries_Full_Map_LS",tmp_fit.entries_Full_Map_LS,"entries_Full_Map_LS[3][648]/I");
    
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
    
    fiTree.Branch("entries_NS_Map_HS",tmp_fit.entries_NS_Map_HS,"entries_NS_Map_HS[3][648]/I");
    fiTree.Branch("entries_EW_Map_HS",tmp_fit.entries_EW_Map_HS,"entries_EW_Map_HS[3][648]/I");
    fiTree.Branch("entries_FB_Map_HS",tmp_fit.entries_FB_Map_HS,"entries_FB_Map_HS[3][648]/I");
    fiTree.Branch("entries_NS_EW_Map_HS",tmp_fit.entries_NS_EW_Map_HS,"entries_NS_EW_Map_HS[3][648]/I");
    fiTree.Branch("entries_NS_FB_Map_HS",tmp_fit.entries_NS_FB_Map_HS,"entries_NS_FB_Map_HS[3][648]/I");
    fiTree.Branch("entries_EW_FB_Map_HS",tmp_fit.entries_EW_FB_Map_HS,"entries_EW_FB_Map_HS[3][648]/I");
    fiTree.Branch("entries_Full_Map_HS",tmp_fit.entries_Full_Map_HS,"entries_Full_Map_HS[3][648]/I");
    
    fiTree.Branch("theta_binHistos_HS",&tmp_fit.theta_binHistos_HS,"theta_binHistos_HS/I");
    fiTree.Branch("phi_binhistos_HS",&tmp_fit.phi_binHistos_HS,"phi_binHistos_HS/I");
        
    fiTree.Branch("events_HS",&tmp_fit.events_HS,"events_HS/l");
    
    fiTree.Branch("inputAni",tmp_fit.inputAni,"inputAni[3]/D");
    
}

void create_test_tree_branches(TTree &fiTree,test_fitResult &tmp_fit)
{
    
    fiTree.Branch("chi2_LS",tmp_fit.chi2_LS,"chi2_LS[4]/D");
    fiTree.Branch("chi2_r_LS",tmp_fit.chi2_r_LS,"chi2_r_LS[4]/D");
    fiTree.Branch("ndf_LS",tmp_fit.ndf_LS,"ndf_LS[4]/D");
    fiTree.Branch("fit_par_LS",tmp_fit.fit_par_LS,"fit_par_LS[4]/D");
    fiTree.Branch("fit_err_LS",tmp_fit.fit_err_LS,"fit_err_LS[4]/D");
    fiTree.Branch("delta_LS",tmp_fit.delta_LS,"delta_LS[4]/D");
    
    fiTree.Branch("theta_binHistos_LS",&tmp_fit.theta_binHistos_LS,"theta_binHistos_LS/I");
    fiTree.Branch("phi_binhistos_LS",&tmp_fit.phi_binHistos_LS,"phi_binHistos_LS/I");
    
    fiTree.Branch("seed",&tmp_fit.seed,"seed/l");
    fiTree.Branch("seed_list_line",&tmp_fit.seed_list_line,"seed_list_line/l");
    
    fiTree.Branch("events_LS",&tmp_fit.events_LS,"events_LS/l");
    
    fiTree.Branch("chi2_HS",tmp_fit.chi2_HS,"chi2_HS[4]/D");
    fiTree.Branch("chi2_r_HS",tmp_fit.chi2_r_HS,"chi2_r_HS[4]/D");
    fiTree.Branch("ndf_HS",tmp_fit.ndf_HS,"ndf_HS[4]/D");
    fiTree.Branch("fit_par_HS",tmp_fit.fit_par_HS,"fit_par_HS[4]/D");
    fiTree.Branch("fit_err_HS",tmp_fit.fit_err_HS,"fit_err_HS[4]/D");
    fiTree.Branch("delta_HS",tmp_fit.delta_HS,"delta_HS[4]/D");
    
    fiTree.Branch("theta_binHistos_HS",&tmp_fit.theta_binHistos_HS,"theta_binHistos_HS/I");
    fiTree.Branch("phi_binhistos_HS",&tmp_fit.phi_binHistos_HS,"phi_binHistos_HS/I");
    
    fiTree.Branch("events_HS",&tmp_fit.events_HS,"events_HS/l");
    
    fiTree.Branch("inputAni",tmp_fit.inputAni,"inputAni[3]/D");
    
}

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string>

#include "TTree.h"
#include "TH1.h"
#include "TFile.h"

class fitResult
{
public:
    
    Double_t inputAni[3];
    
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
        
        events_LS = 0;
        events_HS = 0;
        
    }
    
    ~fitResult() { }
    
};


void tree_linking(TTree* myTree,fitResult &tmp_result,UInt_t &tree_entries);
std::string get_base_path(std::string tree_path,std::string output_dir);
std::string get_final_res_path(std::string base_output_path,std::string ani_level,bool low_statistic);

void extractor(std::string tree_path,std::string output_dir)
{
    
    std::string base_output_path = get_base_path(tree_path,output_dir);
    
    ///////////////// Creating output path strings
    
    std::string results_path_10_LS = get_final_res_path(base_output_path,"10",true);
    std::string results_path_1_LS = get_final_res_path(base_output_path,"1",true);
    std::string results_path_01_LS = get_final_res_path(base_output_path,"01",true);
    
    std::string results_path_10_HS = get_final_res_path(base_output_path,"10",false);
    std::string results_path_1_HS = get_final_res_path(base_output_path,"1",false);
    std::string results_path_01_HS = get_final_res_path(base_output_path,"01",false);
    
    ////////////////////////////////////////
    
    UInt_t tree_entries = 0;
    TFile inFile(tree_path.c_str(),"READ");
    if(inFile.IsZombie())
    {
        std::cout << "\n\nError opening TTree ROOT file \n\n";
        exit(100);
    }
    
    TTree* myTree = (TTree*)inFile.Get("fiTree");
    fitResult tmp_result;
    
    tree_linking(myTree,tmp_result,tree_entries);
    
    ////////////////////////////// Histos //////////////////////////////
    
    ////////////////////////////////////////////////////////// LS Histos
    
    ///////////////////////////// 10% anisotropy
    
    ///////// HS Isotropic Sky
    
    TH1D par_Iso_Iso_LS_10("par_Iso_Iso_LS_10","Iso Parameter (LS Full Isotropic Sky)",100,16398,16420);
    TH1D par_NS_Iso_LS_10("par_NS_Iso_LS_10","NS Parameter (LS Full Isotropic Sky)",100,-100,100);
    TH1D par_EW_Iso_LS_10("par_EW_Iso_LS_10","EW Parameter (LS Full Isotropic Sky)",100,-100,100);
    TH1D par_FB_Iso_LS_10("par_FB_Iso_LS_10","FB Parameter (LS Full Isotropic Sky)",100,-100,100);
    
    TH1D parerr_Iso_Iso_LS_10("parerr_Iso_Iso_LS_10","Iso Parameter Error (LS Full Isotropic Sky)",100,8,12);
    TH1D parerr_NS_Iso_LS_10("parerr_NS_Iso_LS_10","NS Parameter Error (LS Full Isotropic Sky)",100,8,12);
    TH1D parerr_EW_Iso_LS_10("parerr_EW_Iso_LS_10","EW Parameter Error (LS Full Isotropic Sky)",100,8,12);
    TH1D parerr_FB_Iso_LS_10("parerr_FB_Iso_LS_10","FB Parameter Error (LS Full Isotropic Sky)",100,8,12);
    
    ///////// HS NS Dipole Sky
    
    TH1D par_Iso_NS_LS_10("par_Iso_NS_LS_10","Iso Parameter (LS NS Sky)",100,16300,16500);
    TH1D par_NS_NS_LS_10("par_NS_NS_LS_10","NS Parameter (LS NS Sky)",100,1500,2000);
    TH1D par_EW_NS_LS_10("par_EW_NS_LS_10","EW Parameter (LS NS Sky)",100,-100,100);
    TH1D par_FB_NS_LS_10("par_FB_NS_LS_10","FB Parameter (LS NS Sky)",100,-100,100);
    
    TH1D parerr_Iso_NS_LS_10("parerr_Iso_NS_LS_10","Iso Parameter Error (LS NS Sky)",100,8,12);
    TH1D parerr_NS_NS_LS_10("parerr_NS_NS_LS_10","NS Parameter Error (LS NS Sky)",100,8,12);
    TH1D parerr_EW_NS_LS_10("parerr_EW_NS_LS_10","EW Parameter Error (LS NS Sky)",100,8,12);
    TH1D parerr_FB_NS_LS_10("parerr_FB_NS_LS_10","FB Parameter Error (LS NS Sky)",100,8,12);
    
    ///////// HS EW Dipole Sky
    
    TH1D par_Iso_EW_LS_10("par_Iso_EW_LS_10","Iso Parameter (LS EW Sky)",100,16300,16500);
    TH1D par_NS_EW_LS_10("par_NS_EW_LS_10","NS Parameter (LS EW Sky)",100,-100,100);
    TH1D par_EW_EW_LS_10("par_EW_EW_LS_10","EW Parameter (LS EW Sky)",100,1500,2000);
    TH1D par_FB_EW_LS_10("par_FB_EW_LS_10","FB Parameter (LS EW Sky)",100,-100,100);
    
    TH1D parerr_Iso_EW_LS_10("parerr_Iso_EW_LS_10","Iso Parameter Error (LS EW Sky)",100,8,12);
    TH1D parerr_NS_EW_LS_10("parerr_NS_EW_LS_10","NS Parameter Error (LS EW Sky)",100,8,12);
    TH1D parerr_EW_EW_LS_10("parerr_EW_EW_LS_10","EW Parameter Error (LS EW Sky)",100,8,12);
    TH1D parerr_FB_EW_LS_10("parerr_FB_EW_LS_10","FB Parameter Error (LS EW Sky)",100,8,12);
    
    ///////// HS FB Dipole Sky
    
    TH1D par_Iso_FB_LS_10("par_Iso_FB_LS_10","Iso Parameter (LS FB Sky)",100,16300,16500);
    TH1D par_NS_FB_LS_10("par_NS_FB_LS_10","NS Parameter (LS FB Sky)",100,-100,100);
    TH1D par_EW_FB_LS_10("par_EW_FB_LS_10","EW Parameter (LS FB Sky)",100,-100,100);
    TH1D par_FB_FB_LS_10("par_FB_FB_LS_10","FB Parameter (LS FB Sky)",100,1500,2000);
    
    TH1D parerr_Iso_FB_LS_10("parerr_Iso_FB_LS_10","Iso Parameter Error (LS FB Sky)",100,8,12);
    TH1D parerr_NS_FB_LS_10("parerr_NS_FB_LS_10","NS Parameter Error (LS FB Sky)",100,8,12);
    TH1D parerr_EW_FB_LS_10("parerr_EW_FB_LS_10","EW Parameter Error (LS FB Sky)",100,8,12);
    TH1D parerr_FB_FB_LS_10("parerr_FB_FB_LS_10","FB Parameter Error (LS FB Sky)",100,8,12);
    
    ///////// HS NS+EW Dipole Sky
    
    TH1D par_Iso_NSEW_LS_10("par_Iso_NSEW_LS_10","Iso Parameter (LS NSEW Sky)",100,16300,16500);
    TH1D par_NS_NSEW_LS_10("par_NS_NSEW_LS_10","NS Parameter (LS NSEW Sky)",100,1800,2200);
    TH1D par_EW_NSEW_LS_10("par_EW_NSEW_LS_10","EW Parameter (LS NSEW Sky)",100,1800,2200);
    TH1D par_FB_NSEW_LS_10("par_FB_NSEW_LS_10","FB Parameter (LS NSEW Sky)",100,-100,100);
    
    TH1D parerr_Iso_NSEW_LS_10("parerr_Iso_NSEW_LS_10","Iso Parameter Error (LS NSEW Sky)",100,8,12);
    TH1D parerr_NS_NSEW_LS_10("parerr_NS_NSEW_LS_10","NS Parameter Error (LS NSEW Sky)",100,8,12);
    TH1D parerr_EW_NSEW_LS_10("parerr_EW_NSEW_LS_10","EW Parameter Error (LS NSEW Sky)",100,8,12);
    TH1D parerr_FB_NSEW_LS_10("parerr_FB_NSEW_LS_10","FB Parameter Error (LS NSEW Sky)",100,8,12);
    
    ///////// HS NS+FB Dipole Sky
    
    TH1D par_Iso_NSFB_LS_10("par_Iso_NSFB_LS_10","Iso Parameter (LS NSFB Sky)",100,16300,16500);
    TH1D par_NS_NSFB_LS_10("par_NS_NSFB_LS_10","NS Parameter (LS NSFB Sky)",100,1800,2200);
    TH1D par_EW_NSFB_LS_10("par_EW_NSFB_LS_10","EW Parameter (LS NSFB Sky)",100,-100,100);
    TH1D par_FB_NSFB_LS_10("par_FB_NSFB_LS_10","FB Parameter (LS NSFB Sky)",100,1800,2200);
    
    TH1D parerr_Iso_NSFB_LS_10("parerr_Iso_NSFB_LS_10","Iso Parameter Error (LS NSFB Sky)",100,8,12);
    TH1D parerr_NS_NSFB_LS_10("parerr_NS_NSFB_LS_10","NS Parameter Error (LS NSFB Sky)",100,8,12);
    TH1D parerr_EW_NSFB_LS_10("parerr_EW_NSFB_LS_10","EW Parameter Error (LS NSFB Sky)",100,8,12);
    TH1D parerr_FB_NSFB_LS_10("parerr_FB_NSFB_LS_10","FB Parameter Error (LS NSFB Sky)",100,8,12);
    
    ///////// HS EW+FB Dipole Sky
    
    TH1D par_Iso_EWFB_LS_10("par_Iso_EWFB_LS_10","Iso Parameter (LS EWFB Sky)",100,16300,16500);
    TH1D par_NS_EWFB_LS_10("par_NS_EWFB_LS_10","NS Parameter (LS EWFB Sky)",100,-100,100);
    TH1D par_EW_EWFB_LS_10("par_EW_EWFB_LS_10","EW Parameter (LS EWFB Sky)",100,1800,2200);
    TH1D par_FB_EWFB_LS_10("par_FB_EWFB_LS_10","FB Parameter (LS EWFB Sky)",100,1800,2200);
    
    TH1D parerr_Iso_EWFB_LS_10("parerr_Iso_EWFB_LS_10","Iso Parameter Error (LS EWFB Sky)",100,8,12);
    TH1D parerr_NS_EWFB_LS_10("parerr_NS_EWFB_LS_10","NS Parameter Error (LS EWFB Sky)",100,8,12);
    TH1D parerr_EW_EWFB_LS_10("parerr_EW_EWFB_LS_10","EW Parameter Error (LS EWFB Sky)",100,8,12);
    TH1D parerr_FB_EWFB_LS_10("parerr_FB_EWFB_LS_10","FB Parameter Error (LS EWFB Sky)",100,8,12);
    
    ///////// HS Full Dipole Sky
    
    TH1D par_Iso_full_LS_10("par_Iso_full_LS_10","Iso Parameter (LS full Sky)",100,16300,16500);
    TH1D par_NS_full_LS_10("par_NS_full_LS_10","NS Parameter (LS full Sky)",100,1800,2200);
    TH1D par_EW_full_LS_10("par_EW_full_LS_10","EW Parameter (LS full Sky)",100,1800,2200);
    TH1D par_FB_full_LS_10("par_FB_full_LS_10","FB Parameter (LS full Sky)",100,1800,2200);
    
    TH1D parerr_Iso_full_LS_10("parerr_Iso_full_LS_10","Iso Parameter Error (LS full Sky)",100,8,12);
    TH1D parerr_NS_full_LS_10("parerr_NS_full_LS_10","NS Parameter Error (LS full Sky)",100,8,12);
    TH1D parerr_EW_full_LS_10("parerr_EW_full_LS_10","EW Parameter Error (LS full Sky)",100,8,12);
    TH1D parerr_FB_full_LS_10("parerr_FB_full_LS_10","FB Parameter Error (LS full Sky)",100,8,12);
    
    /////////////////////////////////
    
    TH1D chi2_Iso_LS_10("chi2_Iso_LS_10","#chi^2 Iso LS",1000,500,1000);
    TH1D chi2_NS_LS_10("chi2_NS_LS_10","#chi^2 NS LS",1000,500,1000);
    TH1D chi2_EW_LS_10("chi2_EW_LS_10","#chi^2 EW LS",1000,500,1000);
    TH1D chi2_FB_LS_10("chi2_FB_LS_10","#chi^2 FB LS",1000,500,1000);
    TH1D chi2_NSEW_LS_10("chi2_NSEW_LS_10","#chi^2 NS+EW LS",1000,500,1000);
    TH1D chi2_NSFB_LS_10("chi2_NSFB_LS_10","#chi^2 NS+FB LS",1000,500,1000);
    TH1D chi2_EWFB_LS_10("chi2_EWFB_LS_10","#chi^2 EW+FB LS",1000,500,1000);
    TH1D chi2_full_LS_10("chi2_full_LS_10","#chi^2 LS",1000,500,1000);
    
    TH1D ndf_Iso_LS_10("ndf_Iso_LS_10","ndf Iso LS",1000,500,1000);
    TH1D ndf_NS_LS_10("ndf_NS_LS_10","ndf NS LS",1000,500,1000);
    TH1D ndf_EW_LS_10("ndf_EW_LS_10","ndf EW LS",1000,500,1000);
    TH1D ndf_FB_LS_10("ndf_FB_LS_10","ndf FB LS",1000,500,1000);
    TH1D ndf_NSEW_LS_10("ndf_NSEW_LS_10","ndf NS+EW LS",1000,500,1000);
    TH1D ndf_NSFB_LS_10("ndf_NSFB_LS_10","ndf NS+FB LS",1000,500,1000);
    TH1D ndf_EWFB_LS_10("ndf_EWFB_LS_10","ndf EW+FB LS",1000,500,1000);
    TH1D ndf_full_LS_10("ndf_full_LS_10","ndf LS",1000,500,1000);
    
    ///////////////////////////// 1% anisotropy
    
    ///////// HS Isotropic Sky
    
    TH1D par_Iso_Iso_LS_1("par_Iso_Iso_LS_1","Iso Parameter (LS Full Isotropic Sky)",100,16300,16500);
    TH1D par_NS_Iso_LS_1("par_NS_Iso_LS_1","NS Parameter (LS Full Isotropic Sky)",100,-100,100);
    TH1D par_EW_Iso_LS_1("par_EW_Iso_LS_1","EW Parameter (LS Full Isotropic Sky)",100,-100,100);
    TH1D par_FB_Iso_LS_1("par_FB_Iso_LS_1","FB Parameter (LS Full Isotropic Sky)",100,-100,100);
    
    TH1D parerr_Iso_Iso_LS_1("parerr_Iso_Iso_LS_1","Iso Parameter Error (LS Full Isotropic Sky)",100,-100,100);
    TH1D parerr_NS_Iso_LS_1("parerr_NS_Iso_LS_1","NS Parameter Error (LS Full Isotropic Sky)",100,-100,100);
    TH1D parerr_EW_Iso_LS_1("parerr_EW_Iso_LS_1","EW Parameter Error (LS Full Isotropic Sky)",100,-100,100);
    TH1D parerr_FB_Iso_LS_1("parerr_FB_Iso_LS_1","FB Parameter Error (LS Full Isotropic Sky)",100,-100,100);
    
    ///////// HS NS Dipole Sky
    
    TH1D par_Iso_NS_LS_1("par_Iso_NS_LS_1","Iso Parameter (LS NS Sky)",100,16300,16500);
    TH1D par_NS_NS_LS_1("par_NS_NS_LS_1","NS Parameter (LS NS Sky)",100,1500,2000);
    TH1D par_EW_NS_LS_1("par_EW_NS_LS_1","EW Parameter (LS NS Sky)",100,-100,100);
    TH1D par_FB_NS_LS_1("par_FB_NS_LS_1","FB Parameter (LS NS Sky)",100,-100,100);
    
    TH1D parerr_Iso_NS_LS_1("parerr_Iso_NS_LS_1","Iso Parameter Error (LS NS Sky)",100,-100,100);
    TH1D parerr_NS_NS_LS_1("parerr_NS_NS_LS_1","NS Parameter Error (LS NS Sky)",100,-100,100);
    TH1D parerr_EW_NS_LS_1("parerr_EW_NS_LS_1","EW Parameter Error (LS NS Sky)",100,-100,100);
    TH1D parerr_FB_NS_LS_1("parerr_FB_NS_LS_1","FB Parameter Error (LS NS Sky)",100,-100,100);
    
    ///////// HS EW Dipole Sky
    
    TH1D par_Iso_EW_LS_1("par_Iso_EW_LS_1","Iso Parameter (LS EW Sky)",100,16300,16500);
    TH1D par_NS_EW_LS_1("par_NS_EW_LS_1","NS Parameter (LS EW Sky)",100,-100,100);
    TH1D par_EW_EW_LS_1("par_EW_EW_LS_1","EW Parameter (LS EW Sky)",100,1500,2000);
    TH1D par_FB_EW_LS_1("par_FB_EW_LS_1","FB Parameter (LS EW Sky)",100,-100,100);
    
    TH1D parerr_Iso_EW_LS_1("parerr_Iso_EW_LS_1","Iso Parameter Error (LS EW Sky)",100,-100,100);
    TH1D parerr_NS_EW_LS_1("parerr_NS_EW_LS_1","NS Parameter Error (LS EW Sky)",100,-100,100);
    TH1D parerr_EW_EW_LS_1("parerr_EW_EW_LS_1","EW Parameter Error (LS EW Sky)",100,-100,100);
    TH1D parerr_FB_EW_LS_1("parerr_FB_EW_LS_1","FB Parameter Error (LS EW Sky)",100,-100,100);
    
    ///////// HS FB Dipole Sky
    
    TH1D par_Iso_FB_LS_1("par_Iso_FB_LS_1","Iso Parameter (LS FB Sky)",100,16300,16500);
    TH1D par_NS_FB_LS_1("par_NS_FB_LS_1","NS Parameter (LS FB Sky)",100,-100,100);
    TH1D par_EW_FB_LS_1("par_EW_FB_LS_1","EW Parameter (LS FB Sky)",100,-100,100);
    TH1D par_FB_FB_LS_1("par_FB_FB_LS_1","FB Parameter (LS FB Sky)",100,1500,2000);
    
    TH1D parerr_Iso_FB_LS_1("parerr_Iso_FB_LS_1","Iso Parameter Error (LS FB Sky)",100,-100,100);
    TH1D parerr_NS_FB_LS_1("parerr_NS_FB_LS_1","NS Parameter Error (LS FB Sky)",100,-100,100);
    TH1D parerr_EW_FB_LS_1("parerr_EW_FB_LS_1","EW Parameter Error (LS FB Sky)",100,-100,100);
    TH1D parerr_FB_FB_LS_1("parerr_FB_FB_LS_1","FB Parameter Error (LS FB Sky)",100,-100,100);
    
    ///////// HS NS+EW Dipole Sky
    
    TH1D par_Iso_NSEW_LS_1("par_Iso_NSEW_LS_1","Iso Parameter (LS NSEW Sky)",100,16300,16500);
    TH1D par_NS_NSEW_LS_1("par_NS_NSEW_LS_1","NS Parameter (LS NSEW Sky)",100,1800,2200);
    TH1D par_EW_NSEW_LS_1("par_EW_NSEW_LS_1","EW Parameter (LS NSEW Sky)",100,1800,2200);
    TH1D par_FB_NSEW_LS_1("par_FB_NSEW_LS_1","FB Parameter (LS NSEW Sky)",100,-100,100);
    
    TH1D parerr_Iso_NSEW_LS_1("parerr_Iso_NSEW_LS_1","Iso Parameter Error (LS NSEW Sky)",100,-100,100);
    TH1D parerr_NS_NSEW_LS_1("parerr_NS_NSEW_LS_1","NS Parameter Error (LS NSEW Sky)",100,-100,100);
    TH1D parerr_EW_NSEW_LS_1("parerr_EW_NSEW_LS_1","EW Parameter Error (LS NSEW Sky)",100,-100,100);
    TH1D parerr_FB_NSEW_LS_1("parerr_FB_NSEW_LS_1","FB Parameter Error (LS NSEW Sky)",100,-100,100);
    
    ///////// HS NS+FB Dipole Sky
    
    TH1D par_Iso_NSFB_LS_1("par_Iso_NSFB_LS_1","Iso Parameter (LS NSFB Sky)",100,16300,16500);
    TH1D par_NS_NSFB_LS_1("par_NS_NSFB_LS_1","NS Parameter (LS NSFB Sky)",100,1800,2200);
    TH1D par_EW_NSFB_LS_1("par_EW_NSFB_LS_1","EW Parameter (LS NSFB Sky)",100,-100,100);
    TH1D par_FB_NSFB_LS_1("par_FB_NSFB_LS_1","FB Parameter (LS NSFB Sky)",100,1800,2200);
    
    TH1D parerr_Iso_NSFB_LS_1("parerr_Iso_NSFB_LS_1","Iso Parameter Error (LS NSFB Sky)",100,-100,100);
    TH1D parerr_NS_NSFB_LS_1("parerr_NS_NSFB_LS_1","NS Parameter Error (LS NSFB Sky)",100,-100,100);
    TH1D parerr_EW_NSFB_LS_1("parerr_EW_NSFB_LS_1","EW Parameter Error (LS NSFB Sky)",100,-100,100);
    TH1D parerr_FB_NSFB_LS_1("parerr_FB_NSFB_LS_1","FB Parameter Error (LS NSFB Sky)",100,-100,100);
    
    ///////// HS EW+FB Dipole Sky
    
    TH1D par_Iso_EWFB_LS_1("par_Iso_EWFB_LS_1","Iso Parameter (LS EWFB Sky)",100,16300,16500);
    TH1D par_NS_EWFB_LS_1("par_NS_EWFB_LS_1","NS Parameter (LS EWFB Sky)",100,-100,100);
    TH1D par_EW_EWFB_LS_1("par_EW_EWFB_LS_1","EW Parameter (LS EWFB Sky)",100,1800,2200);
    TH1D par_FB_EWFB_LS_1("par_FB_EWFB_LS_1","FB Parameter (LS EWFB Sky)",100,1800,2200);
    
    TH1D parerr_Iso_EWFB_LS_1("parerr_Iso_EWFB_LS_1","Iso Parameter Error (LS EWFB Sky)",100,-100,100);
    TH1D parerr_NS_EWFB_LS_1("parerr_NS_EWFB_LS_1","NS Parameter Error (LS EWFB Sky)",100,-100,100);
    TH1D parerr_EW_EWFB_LS_1("parerr_EW_EWFB_LS_1","EW Parameter Error (LS EWFB Sky)",100,-100,100);
    TH1D parerr_FB_EWFB_LS_1("parerr_FB_EWFB_LS_1","FB Parameter Error (LS EWFB Sky)",100,-100,100);
    
    ///////// HS Full Dipole Sky
    
    TH1D par_Iso_full_LS_1("par_Iso_full_LS_1","Iso Parameter (LS full Sky)",100,16300,16500);
    TH1D par_NS_full_LS_1("par_NS_full_LS_1","NS Parameter (LS full Sky)",100,1800,2200);
    TH1D par_EW_full_LS_1("par_EW_full_LS_1","EW Parameter (LS full Sky)",100,1800,2200);
    TH1D par_FB_full_LS_1("par_FB_full_LS_1","FB Parameter (LS full Sky)",100,1800,2200);
    
    TH1D parerr_Iso_full_LS_1("parerr_Iso_full_LS_1","Iso Parameter Error (LS full Sky)",100,-100,100);
    TH1D parerr_NS_full_LS_1("parerr_NS_full_LS_1","NS Parameter Error (LS full Sky)",100,-100,100);
    TH1D parerr_EW_full_LS_1("parerr_EW_full_LS_1","EW Parameter Error (LS full Sky)",100,-100,100);
    TH1D parerr_FB_full_LS_1("parerr_FB_full_LS_1","FB Parameter Error (LS full Sky)",100,-100,100);
    
    /////////////////////////////////
    
    TH1D chi2_Iso_LS_1("chi2_Iso_LS_1","#chi^2 Iso LS",1000,500,1000);
    TH1D chi2_NS_LS_1("chi2_NS_LS_1","#chi^2 NS LS",1000,500,1000);
    TH1D chi2_EW_LS_1("chi2_EW_LS_1","#chi^2 EW LS",1000,500,1000);
    TH1D chi2_FB_LS_1("chi2_FB_LS_1","#chi^2 FB LS",1000,500,1000);
    TH1D chi2_NSEW_LS_1("chi2_NSEW_LS_1","#chi^2 NS+EW LS",1000,500,1000);
    TH1D chi2_NSFB_LS_1("chi2_NSFB_LS_1","#chi^2 NS+FB LS",1000,500,1000);
    TH1D chi2_EWFB_LS_1("chi2_EWFB_LS_1","#chi^2 EW+FB LS",1000,500,1000);
    TH1D chi2_full_LS_1("chi2_full_LS_1","#chi^2 LS",1000,500,1000);
    
    TH1D ndf_Iso_LS_1("ndf_Iso_LS_1","ndf Iso LS",1000,500,1000);
    TH1D ndf_NS_LS_1("ndf_NS_LS_1","ndf NS LS",1000,500,1000);
    TH1D ndf_EW_LS_1("ndf_EW_LS_1","ndf EW LS",1000,500,1000);
    TH1D ndf_FB_LS_1("ndf_FB_LS_1","ndf FB LS",1000,500,1000);
    TH1D ndf_NSEW_LS_1("ndf_NSEW_LS_1","ndf NS+EW LS",1000,500,1000);
    TH1D ndf_NSFB_LS_1("ndf_NSFB_LS_1","ndf NS+FB LS",1000,500,1000);
    TH1D ndf_EWFB_LS_1("ndf_EWFB_LS_1","ndf EW+FB LS",1000,500,1000);
    TH1D ndf_full_LS_1("ndf_full_LS_1","ndf LS",1000,500,1000);
    
    
    ///////////////////////////// 0.1% anisotropy
    
    ///////// HS Isotropic Sky
    
    TH1D par_Iso_Iso_LS_01("par_Iso_Iso_LS_01","Iso Parameter (LS Full Isotropic Sky)",100,16300,16500);
    TH1D par_NS_Iso_LS_01("par_NS_Iso_LS_01","NS Parameter (LS Full Isotropic Sky)",100,-100,100);
    TH1D par_EW_Iso_LS_01("par_EW_Iso_LS_01","EW Parameter (LS Full Isotropic Sky)",100,-100,100);
    TH1D par_FB_Iso_LS_01("par_FB_Iso_LS_01","FB Parameter (LS Full Isotropic Sky)",100,-100,100);
    
    TH1D parerr_Iso_Iso_LS_01("parerr_Iso_Iso_LS_01","Iso Parameter Error (LS Full Isotropic Sky)",100,-100,100);
    TH1D parerr_NS_Iso_LS_01("parerr_NS_Iso_LS_01","NS Parameter Error (LS Full Isotropic Sky)",100,-100,100);
    TH1D parerr_EW_Iso_LS_01("parerr_EW_Iso_LS_01","EW Parameter Error (LS Full Isotropic Sky)",100,-100,100);
    TH1D parerr_FB_Iso_LS_01("parerr_FB_Iso_LS_01","FB Parameter Error (LS Full Isotropic Sky)",100,-100,100);
    
    ///////// HS NS Dipole Sky
    
    TH1D par_Iso_NS_LS_01("par_Iso_NS_LS_01","Iso Parameter (LS NS Sky)",100,16300,16500);
    TH1D par_NS_NS_LS_01("par_NS_NS_LS_01","NS Parameter (LS NS Sky)",100,-100,100);
    TH1D par_EW_NS_LS_01("par_EW_NS_LS_01","EW Parameter (LS NS Sky)",100,-100,100);
    TH1D par_FB_NS_LS_01("par_FB_NS_LS_01","FB Parameter (LS NS Sky)",100,-100,100);
    
    TH1D parerr_Iso_NS_LS_01("parerr_Iso_NS_LS_01","Iso Parameter (LS NS Sky)",100,-100,100);
    TH1D parerr_NS_NS_LS_01("parerr_NS_NS_LS_01","NS Parameter (LS NS Sky)",100,-100,100);
    TH1D parerr_EW_NS_LS_01("parerr_EW_NS_LS_01","EW Parameter (LS NS Sky)",100,-100,100);
    TH1D parerr_FB_NS_LS_01("parerr_FB_NS_LS_01","FB Parameter (LS NS Sky)",100,-100,100);
    
    ///////// HS EW Dipole Sky
    
    TH1D par_Iso_EW_LS_01("par_Iso_EW_LS_01","Iso Parameter (LS EW Sky)",100,16300,16500);
    TH1D par_NS_EW_LS_01("par_NS_EW_LS_01","NS Parameter (LS EW Sky)",100,-100,100);
    TH1D par_EW_EW_LS_01("par_EW_EW_LS_01","EW Parameter (LS EW Sky)",100,-100,100);
    TH1D par_FB_EW_LS_01("par_FB_EW_LS_01","FB Parameter (LS EW Sky)",100,-100,100);
    
    TH1D parerr_Iso_EW_LS_01("parerr_Iso_EW_LS_01","Iso Parameter Error (LS EW Sky)",100,-100,100);
    TH1D parerr_NS_EW_LS_01("parerr_NS_EW_LS_01","NS Parameter Error (LS EW Sky)",100,-100,100);
    TH1D parerr_EW_EW_LS_01("parerr_EW_EW_LS_01","EW Parameter Error (LS EW Sky)",100,-100,100);
    TH1D parerr_FB_EW_LS_01("parerr_FB_EW_LS_01","FB Parameter Error (LS EW Sky)",100,-100,100);
    
    ///////// HS FB Dipole Sky
    
    TH1D par_Iso_FB_LS_01("par_Iso_FB_LS_01","Iso Parameter (LS FB Sky)",100,16300,16500);
    TH1D par_NS_FB_LS_01("par_NS_FB_LS_01","NS Parameter (LS FB Sky)",100,-100,100);
    TH1D par_EW_FB_LS_01("par_EW_FB_LS_01","EW Parameter (LS FB Sky)",100,-100,100);
    TH1D par_FB_FB_LS_01("par_FB_FB_LS_01","FB Parameter (LS FB Sky)",100,-100,100);
    
    TH1D parerr_Iso_FB_LS_01("parerr_Iso_FB_LS_01","Iso Parameter Error (LS FB Sky)",100,-100,100);
    TH1D parerr_NS_FB_LS_01("parerr_NS_FB_LS_01","NS Parameter Error (LS FB Sky)",100,-100,100);
    TH1D parerr_EW_FB_LS_01("parerr_EW_FB_LS_01","EW Parameter Error (LS FB Sky)",100,-100,100);
    TH1D parerr_FB_FB_LS_01("parerr_FB_FB_LS_01","FB Parameter Error (LS FB Sky)",100,-100,100);
    
    ///////// HS NS+EW Dipole Sky
    
    TH1D par_Iso_NSEW_LS_01("par_Iso_NSEW_LS_01","Iso Parameter (LS NSEW Sky)",100,16300,16500);
    TH1D par_NS_NSEW_LS_01("par_NS_NSEW_LS_01","NS Parameter (LS NSEW Sky)",100,-100,100);
    TH1D par_EW_NSEW_LS_01("par_EW_NSEW_LS_01","EW Parameter (LS NSEW Sky)",100,-100,100);
    TH1D par_FB_NSEW_LS_01("par_FB_NSEW_LS_01","FB Parameter (LS NSEW Sky)",100,-100,100);
    
    TH1D parerr_Iso_NSEW_LS_01("parerr_Iso_NSEW_LS_01","Iso Parameter Error (LS NSEW Sky)",100,-100,100);
    TH1D parerr_NS_NSEW_LS_01("parerr_NS_NSEW_LS_01","NS Parameter Error (LS NSEW Sky)",100,-100,100);
    TH1D parerr_EW_NSEW_LS_01("parerr_EW_NSEW_LS_01","EW Parameter Error (LS NSEW Sky)",100,-100,100);
    TH1D parerr_FB_NSEW_LS_01("parerr_FB_NSEW_LS_01","FB Parameter Error (LS NSEW Sky)",100,-100,100);
    
    ///////// HS NS+FB Dipole Sky
    
    TH1D par_Iso_NSFB_LS_01("par_Iso_NSFB_LS_01","Iso Parameter (LS NSFB Sky)",100,16300,16500);
    TH1D par_NS_NSFB_LS_01("par_NS_NSFB_LS_01","NS Parameter (LS NSFB Sky)",100,-100,100);
    TH1D par_EW_NSFB_LS_01("par_EW_NSFB_LS_01","EW Parameter (LS NSFB Sky)",100,-100,100);
    TH1D par_FB_NSFB_LS_01("par_FB_NSFB_LS_01","FB Parameter (LS NSFB Sky)",100,-100,100);
    
    TH1D parerr_Iso_NSFB_LS_01("parerr_Iso_NSFB_LS_01","Iso Parameter Error (LS NSFB Sky)",100,-100,100);
    TH1D parerr_NS_NSFB_LS_01("parerr_NS_NSFB_LS_01","NS Parameter Error (LS NSFB Sky)",100,-100,100);
    TH1D parerr_EW_NSFB_LS_01("parerr_EW_NSFB_LS_01","EW Parameter Error (LS NSFB Sky)",100,-100,100);
    TH1D parerr_FB_NSFB_LS_01("parerr_FB_NSFB_LS_01","FB Parameter Error (LS NSFB Sky)",100,-100,100);
    
    ///////// HS EW+FB Dipole Sky
    
    TH1D par_Iso_EWFB_LS_01("par_Iso_EWFB_LS_01","Iso Parameter (LS EWFB Sky)",100,16300,16500);
    TH1D par_NS_EWFB_LS_01("par_NS_EWFB_LS_01","NS Parameter (LS EWFB Sky)",100,-100,100);
    TH1D par_EW_EWFB_LS_01("par_EW_EWFB_LS_01","EW Parameter (LS EWFB Sky)",100,-100,100);
    TH1D par_FB_EWFB_LS_01("par_FB_EWFB_LS_01","FB Parameter (LS EWFB Sky)",100,-100,100);
    
    TH1D parerr_Iso_EWFB_LS_01("parerr_Iso_EWFB_LS_01","Iso Parameter Error (LS EWFB Sky)",100,-100,100);
    TH1D parerr_NS_EWFB_LS_01("parerr_NS_EWFB_LS_01","NS Parameter Error (LS EWFB Sky)",100,-100,100);
    TH1D parerr_EW_EWFB_LS_01("parerr_EW_EWFB_LS_01","EW Parameter Error (LS EWFB Sky)",100,-100,100);
    TH1D parerr_FB_EWFB_LS_01("parerr_FB_EWFB_LS_01","FB Parameter Error (LS EWFB Sky)",100,-100,100);
    
    ///////// HS Full Dipole Sky
    
    TH1D par_Iso_full_LS_01("par_Iso_full_LS_01","Iso Parameter (LS full Sky)",100,16300,16500);
    TH1D par_NS_full_LS_01("par_NS_full_LS_01","NS Parameter (LS full Sky)",100,-100,100);
    TH1D par_EW_full_LS_01("par_EW_full_LS_01","EW Parameter (LS full Sky)",100,-100,100);
    TH1D par_FB_full_LS_01("par_FB_full_LS_01","FB Parameter (LS full Sky)",100,-100,100);
    
    TH1D parerr_Iso_full_LS_01("parerr_Iso_full_LS_01","Iso Parameter Error (LS full Sky)",100,-100,100);
    TH1D parerr_NS_full_LS_01("parerr_NS_full_LS_01","NS Parameter Error (LS full Sky)",100,-100,100);
    TH1D parerr_EW_full_LS_01("parerr_EW_full_LS_01","EW Parameter Error (LS full Sky)",100,-100,100);
    TH1D parerr_FB_full_LS_01("parerr_FB_full_LS_01","FB Parameter Error (LS full Sky)",100,-100,100);
    
    /////////////////////////////////
    
    TH1D chi2_Iso_LS_01("chi2_Iso_LS_01","#chi^2 Iso LS",1000,500,1000);
    TH1D chi2_NS_LS_01("chi2_NS_LS_01","#chi^2 NS LS",1000,500,1000);
    TH1D chi2_EW_LS_01("chi2_EW_LS_01","#chi^2 EW LS",1000,500,1000);
    TH1D chi2_FB_LS_01("chi2_FB_LS_01","#chi^2 FB LS",1000,500,1000);
    TH1D chi2_NSEW_LS_01("chi2_NSEW_LS_01","#chi^2 NS+EW LS",1000,500,1000);
    TH1D chi2_NSFB_LS_01("chi2_NSFB_LS_01","#chi^2 NS+FB LS",1000,500,1000);
    TH1D chi2_EWFB_LS_01("chi2_EWFB_LS_01","#chi^2 EW+FB LS",1000,500,1000);
    TH1D chi2_full_LS_01("chi2_full_LS_01","#chi^2 LS",1000,500,1000);
    
    TH1D ndf_Iso_LS_01("ndf_Iso_LS_01","ndf Iso LS",1000,500,1000);
    TH1D ndf_NS_LS_01("ndf_NS_LS_01","ndf NS LS",1000,500,1000);
    TH1D ndf_EW_LS_01("ndf_EW_LS_01","ndf EW LS",1000,500,1000);
    TH1D ndf_FB_LS_01("ndf_FB_LS_01","ndf FB LS",1000,500,1000);
    TH1D ndf_NSEW_LS_01("ndf_NSEW_LS_01","ndf NS+EW LS",1000,500,1000);
    TH1D ndf_NSFB_LS_01("ndf_NSFB_LS_01","ndf NS+FB LS",1000,500,1000);
    TH1D ndf_EWFB_LS_01("ndf_EWFB_LS_01","ndf EW+FB LS",1000,500,1000);
    TH1D ndf_full_LS_01("ndf_full_LS_01","ndf LS",1000,500,1000);
    
    
    ////////////////////////////////////////////////////////////////////
    
    ////////////////////////////////////////////////////////// HS Histos
    
    ///////////////////////////// 10% anisotropy
    
    ///////// HS Isotropic Sky
    
    TH1D par_Iso_Iso_HS_10("par_Iso_Iso_HS_10","Iso Parameter (HS Full Isotropic Sky)",100,32000,33200);
    TH1D par_NS_Iso_HS_10("par_NS_Iso_HS_10","NS Parameter (HS Full Isotropic Sky)",100,-100,100);
    TH1D par_EW_Iso_HS_10("par_EW_Iso_HS_10","EW Parameter (HS Full Isotropic Sky)",100,-100,100);
    TH1D par_FB_Iso_HS_10("par_FB_Iso_HS_10","FB Parameter (HS Full Isotropic Sky)",100,-100,100);
    
    TH1D parerr_Iso_Iso_HS_10("parerr_Iso_Iso_HS_10","Iso Parameter Error (HS Full Isotropic Sky)",100,-100,100);
    TH1D parerr_NS_Iso_HS_10("parerr_NS_Iso_HS_10","NS Parameter Error (HS Full Isotropic Sky)",100,-100,100);
    TH1D parerr_EW_Iso_HS_10("parerr_EW_Iso_HS_10","EW Parameter Error (HS Full Isotropic Sky)",100,-100,100);
    TH1D parerr_FB_Iso_HS_10("parerr_FB_Iso_HS_10","FB Parameter Error (HS Full Isotropic Sky)",100,-100,100);
    
    ///////// HS NS Dipole Sky
    
    TH1D par_Iso_NS_HS_10("par_Iso_NS_HS_10","Iso Parameter (HS NS Sky)",100,32000,33200);
    TH1D par_NS_NS_HS_10("par_NS_NS_HS_10","NS Parameter (HS NS Sky)",100,3000,5000);
    TH1D par_EW_NS_HS_10("par_EW_NS_HS_10","EW Parameter (HS NS Sky)",100,-100,100);
    TH1D par_FB_NS_HS_10("par_FB_NS_HS_10","FB Parameter (HS NS Sky)",100,-100,100);
    
    TH1D parerr_Iso_NS_HS_10("parerr_Iso_NS_HS_10","Iso Parameter Error (HS NS Sky)",100,-100,100);
    TH1D parerr_NS_NS_HS_10("parerr_NS_NS_HS_10","NS Parameter Error (HS NS Sky)",100,-100,100);
    TH1D parerr_EW_NS_HS_10("parerr_EW_NS_HS_10","EW Parameter Error (HS NS Sky)",100,-100,100);
    TH1D parerr_FB_NS_HS_10("parerr_FB_NS_HS_10","FB Parameter Error (HS NS Sky)",100,-100,100);
    
    ///////// HS EW Dipole Sky
    
    TH1D par_Iso_EW_HS_10("par_Iso_EW_HS_10","Iso Parameter (HS EW Sky)",100,32000,33200);
    TH1D par_NS_EW_HS_10("par_NS_EW_HS_10","NS Parameter (HS EW Sky)",100,-100,100);
    TH1D par_EW_EW_HS_10("par_EW_EW_HS_10","EW Parameter (HS EW Sky)",100,3000,5000);
    TH1D par_FB_EW_HS_10("par_FB_EW_HS_10","FB Parameter (HS EW Sky)",100,-100,100);
    
    TH1D parerr_Iso_EW_HS_10("parerr_Iso_EW_HS_10","Iso Parameter Error (HS EW Sky)",100,-100,100);
    TH1D parerr_NS_EW_HS_10("parerr_NS_EW_HS_10","NS Parameter Error (HS EW Sky)",100,-100,100);
    TH1D parerr_EW_EW_HS_10("parerr_EW_EW_HS_10","EW Parameter Error (HS EW Sky)",100,-100,100);
    TH1D parerr_FB_EW_HS_10("parerr_FB_EW_HS_10","FB Parameter Error (HS EW Sky)",100,-100,100);
    
    ///////// HS FB Dipole Sky
    
    TH1D par_Iso_FB_HS_10("par_Iso_FB_HS_10","Iso Parameter (HS FB Sky)",100,32000,33200);
    TH1D par_NS_FB_HS_10("par_NS_FB_HS_10","NS Parameter (HS FB Sky)",100,-100,100);
    TH1D par_EW_FB_HS_10("par_EW_FB_HS_10","EW Parameter (HS FB Sky)",100,-100,100);
    TH1D par_FB_FB_HS_10("par_FB_FB_HS_10","FB Parameter (HS FB Sky)",100,3000,5000);
    
    TH1D parerr_Iso_FB_HS_10("parerr_Iso_FB_HS_10","Iso Parameter Error (HS FB Sky)",100,-100,100);
    TH1D parerr_NS_FB_HS_10("parerr_NS_FB_HS_10","NS Parameter Error (HS FB Sky)",100,-100,100);
    TH1D parerr_EW_FB_HS_10("parerr_EW_FB_HS_10","EW Parameter Error (HS FB Sky)",100,-100,100);
    TH1D parerr_FB_FB_HS_10("parerr_FB_FB_HS_10","FB Parameter Error (HS FB Sky)",100,-100,100);
    
    ///////// HS NS+EW Dipole Sky
    
    TH1D par_Iso_NSEW_HS_10("par_Iso_NSEW_HS_10","Iso Parameter (HS NSEW Sky)",100,32000,33200);
    TH1D par_NS_NSEW_HS_10("par_NS_NSEW_HS_10","NS Parameter (HS NSEW Sky)",100,3000,5000);
    TH1D par_EW_NSEW_HS_10("par_EW_NSEW_HS_10","EW Parameter (HS NSEW Sky)",100,3000,5000);
    TH1D par_FB_NSEW_HS_10("par_FB_NSEW_HS_10","FB Parameter (HS NSEW Sky)",100,-100,100);
    
    TH1D parerr_Iso_NSEW_HS_10("parerr_Iso_NSEW_HS_10","Iso Parameter Error (HS NSEW Sky)",100,-100,100);
    TH1D parerr_NS_NSEW_HS_10("parerr_NS_NSEW_HS_10","NS Parameter Error (HS NSEW Sky)",100,-100,100);
    TH1D parerr_EW_NSEW_HS_10("parerr_EW_NSEW_HS_10","EW Parameter Error (HS NSEW Sky)",100,-100,100);
    TH1D parerr_FB_NSEW_HS_10("parerr_FB_NSEW_HS_10","FB Parameter Error (HS NSEW Sky)",100,-100,100);
    
    ///////// HS NS+FB Dipole Sky
    
    TH1D par_Iso_NSFB_HS_10("par_Iso_NSFB_HS_10","Iso Parameter (HS NSFB Sky)",100,32000,33200);
    TH1D par_NS_NSFB_HS_10("par_NS_NSFB_HS_10","NS Parameter (HS NSFB Sky)",100,3000,5000);
    TH1D par_EW_NSFB_HS_10("par_EW_NSFB_HS_10","EW Parameter (HS NSFB Sky)",100,-100,100);
    TH1D par_FB_NSFB_HS_10("par_FB_NSFB_HS_10","FB Parameter (HS NSFB Sky)",100,3000,5000);
    
    TH1D parerr_Iso_NSFB_HS_10("parerr_Iso_NSFB_HS_10","Iso Parameter Error (HS NSFB Sky)",100,-100,100);
    TH1D parerr_NS_NSFB_HS_10("parerr_NS_NSFB_HS_10","NS Parameter Error (HS NSFB Sky)",100,-100,100);
    TH1D parerr_EW_NSFB_HS_10("parerr_EW_NSFB_HS_10","EW Parameter Error (HS NSFB Sky)",100,-100,100);
    TH1D parerr_FB_NSFB_HS_10("parerr_FB_NSFB_HS_10","FB Parameter Error (HS NSFB Sky)",100,-100,100);
    
    ///////// HS EW+FB Dipole Sky
    
    TH1D par_Iso_EWFB_HS_10("par_Iso_EWFB_HS_10","Iso Parameter (HS EWFB Sky)",100,32000,33200);
    TH1D par_NS_EWFB_HS_10("par_NS_EWFB_HS_10","NS Parameter (HS EWFB Sky)",100,-100,100);
    TH1D par_EW_EWFB_HS_10("par_EW_EWFB_HS_10","EW Parameter (HS EWFB Sky)",100,3000,5000);
    TH1D par_FB_EWFB_HS_10("par_FB_EWFB_HS_10","FB Parameter (HS EWFB Sky)",100,3000,5000);
    
    TH1D parerr_Iso_EWFB_HS_10("parerr_Iso_EWFB_HS_10","Iso Parameter Error (HS EWFB Sky)",100,-100,100);
    TH1D parerr_NS_EWFB_HS_10("parerr_NS_EWFB_HS_10","NS Parameter Error (HS EWFB Sky)",100,-100,100);
    TH1D parerr_EW_EWFB_HS_10("parerr_EW_EWFB_HS_10","EW Parameter Error (HS EWFB Sky)",100,-100,100);
    TH1D parerr_FB_EWFB_HS_10("parerr_FB_EWFB_HS_10","FB Parameter Error (HS EWFB Sky)",100,-100,100);
    
    ///////// HS Full Dipole Sky
    
    TH1D par_Iso_full_HS_10("par_Iso_full_HS_10","Iso Parameter (HS full Sky)",100,32000,33200);
    TH1D par_NS_full_HS_10("par_NS_full_HS_10","NS Parameter (HS full Sky)",100,3000,5000);
    TH1D par_EW_full_HS_10("par_EW_full_HS_10","EW Parameter (HS full Sky)",100,3000,5000);
    TH1D par_FB_full_HS_10("par_FB_full_HS_10","FB Parameter (HS full Sky)",100,3000,5000);
    
    TH1D parerr_Iso_full_HS_10("parerr_Iso_full_HS_10","Iso Parameter Error (HS full Sky)",100,-100,100);
    TH1D parerr_NS_full_HS_10("parerr_NS_full_HS_10","NS Parameter Error (HS full Sky)",100,-100,100);
    TH1D parerr_EW_full_HS_10("parerr_EW_full_HS_10","EW Parameter Error (HS full Sky)",100,-100,100);
    TH1D parerr_FB_full_HS_10("parerr_FB_full_HS_10","FB Parameter Error (HS full Sky)",100,-100,100);
    
    /////////////////////////////////
    
    TH1D chi2_Iso_HS_10("chi2_Iso_HS_10","#chi^2 Iso LS",1000,500,1000);
    TH1D chi2_NS_HS_10("chi2_NS_HS_10","#chi^2 NS LS",1000,500,1000);
    TH1D chi2_EW_HS_10("chi2_EW_HS_10","#chi^2 EW LS",1000,500,1000);
    TH1D chi2_FB_HS_10("chi2_FB_HS_10","#chi^2 FB LS",1000,500,1000);
    TH1D chi2_NSEW_HS_10("chi2_NSEW_HS_10","#chi^2 NS+EW LS",1000,500,1000);
    TH1D chi2_NSFB_HS_10("chi2_NSFB_HS_10","#chi^2 NS+FB LS",1000,500,1000);
    TH1D chi2_EWFB_HS_10("chi2_EWFB_HS_10","#chi^2 EW+FB LS",1000,500,1000);
    TH1D chi2_full_HS_10("chi2_full_HS_10","#chi^2 LS",1000,500,1000);
    
    TH1D ndf_Iso_HS_10("ndf_Iso_HS_10","ndf Iso LS",1000,500,1000);
    TH1D ndf_NS_HS_10("ndf_NS_HS_10","ndf NS LS",1000,500,1000);
    TH1D ndf_EW_HS_10("ndf_EW_HS_10","ndf EW LS",1000,500,1000);
    TH1D ndf_FB_HS_10("ndf_FB_HS_10","ndf FB LS",1000,500,1000);
    TH1D ndf_NSEW_HS_10("ndf_NSEW_HS_10","ndf NS+EW LS",1000,500,1000);
    TH1D ndf_NSFB_HS_10("ndf_NSFB_HS_10","ndf NS+FB LS",1000,500,1000);
    TH1D ndf_EWFB_HS_10("ndf_EWFB_HS_10","ndf EW+FB LS",1000,500,1000);
    TH1D ndf_full_HS_10("ndf_full_HS_10","ndf LS",1000,500,1000);
    
    ///////////////////////////// 1% anisotropy
    
    ///////// HS Isotropic Sky
    
    TH1D par_Iso_Iso_HS_1("par_Iso_Iso_HS_1","Iso Parameter (HS Full Isotropic Sky)",100,16300,16500);
    TH1D par_NS_Iso_HS_1("par_NS_Iso_HS_1","NS Parameter (HS Full Isotropic Sky)",100,-100,100);
    TH1D par_EW_Iso_HS_1("par_EW_Iso_HS_1","EW Parameter (HS Full Isotropic Sky)",100,-100,100);
    TH1D par_FB_Iso_HS_1("par_FB_Iso_HS_1","FB Parameter (HS Full Isotropic Sky)",100,-100,100);
    
    TH1D parerr_Iso_Iso_HS_1("parerr_Iso_Iso_HS_1","Iso Parameter Error (HS Full Isotropic Sky)",100,-100,100);
    TH1D parerr_NS_Iso_HS_1("parerr_NS_Iso_HS_1","NS Parameter Error (HS Full Isotropic Sky)",100,-100,100);
    TH1D parerr_EW_Iso_HS_1("parerr_EW_Iso_HS_1","EW Parameter Error (HS Full Isotropic Sky)",100,-100,100);
    TH1D parerr_FB_Iso_HS_1("parerr_FB_Iso_HS_1","FB Parameter Error (HS Full Isotropic Sky)",100,-100,100);
    
    ///////// HS NS Dipole Sky
    
    TH1D par_Iso_NS_HS_1("par_Iso_NS_HS_1","Iso Parameter (HS NS Sky)",100,16300,16500);
    TH1D par_NS_NS_HS_1("par_NS_NS_HS_1","NS Parameter (HS NS Sky)",100,1500,2000);
    TH1D par_EW_NS_HS_1("par_EW_NS_HS_1","EW Parameter (HS NS Sky)",100,-100,100);
    TH1D par_FB_NS_HS_1("par_FB_NS_HS_1","FB Parameter (HS NS Sky)",100,-100,100);
    
    TH1D parerr_Iso_NS_HS_1("parerr_Iso_NS_HS_1","Iso Parameter Error (HS NS Sky)",100,-100,100);
    TH1D parerr_NS_NS_HS_1("parerr_NS_NS_HS_1","NS Parameter Error (HS NS Sky)",100,-100,100);
    TH1D parerr_EW_NS_HS_1("parerr_EW_NS_HS_1","EW Parameter Error (HS NS Sky)",100,-100,100);
    TH1D parerr_FB_NS_HS_1("parerr_FB_NS_HS_1","FB Parameter Error (HS NS Sky)",100,-100,100);
    
    ///////// HS EW Dipole Sky
    
    TH1D par_Iso_EW_HS_1("par_Iso_EW_HS_1","Iso Parameter (HS EW Sky)",100,16300,16500);
    TH1D par_NS_EW_HS_1("par_NS_EW_HS_1","NS Parameter (HS EW Sky)",100,-100,100);
    TH1D par_EW_EW_HS_1("par_EW_EW_HS_1","EW Parameter (HS EW Sky)",100,1500,2000);
    TH1D par_FB_EW_HS_1("par_FB_EW_HS_1","FB Parameter (HS EW Sky)",100,-100,100);
    
    TH1D parerr_Iso_EW_HS_1("parerr_Iso_EW_HS_1","Iso Parameter Error (HS EW Sky)",100,-100,100);
    TH1D parerr_NS_EW_HS_1("parerr_NS_EW_HS_1","NS Parameter Error (HS EW Sky)",100,-100,100);
    TH1D parerr_EW_EW_HS_1("parerr_EW_EW_HS_1","EW Parameter Error (HS EW Sky)",100,-100,100);
    TH1D parerr_FB_EW_HS_1("parerr_FB_EW_HS_1","FB Parameter Error (HS EW Sky)",100,-100,100);
    
    ///////// HS FB Dipole Sky
    
    TH1D par_Iso_FB_HS_1("par_Iso_FB_HS_1","Iso Parameter (HS FB Sky)",100,16300,16500);
    TH1D par_NS_FB_HS_1("par_NS_FB_HS_1","NS Parameter (HS FB Sky)",100,-100,100);
    TH1D par_EW_FB_HS_1("par_EW_FB_HS_1","EW Parameter (HS FB Sky)",100,-100,100);
    TH1D par_FB_FB_HS_1("par_FB_FB_HS_1","FB Parameter (HS FB Sky)",100,1500,2000);
    
    TH1D parerr_Iso_FB_HS_1("parerr_Iso_FB_HS_1","Iso Parameter Error (HS FB Sky)",100,-100,100);
    TH1D parerr_NS_FB_HS_1("parerr_NS_FB_HS_1","NS Parameter Error (HS FB Sky)",100,-100,100);
    TH1D parerr_EW_FB_HS_1("parerr_EW_FB_HS_1","EW Parameter Error (HS FB Sky)",100,-100,100);
    TH1D parerr_FB_FB_HS_1("parerr_FB_FB_HS_1","FB Parameter Error (HS FB Sky)",100,-100,100);
    
    ///////// HS NS+EW Dipole Sky
    
    TH1D par_Iso_NSEW_HS_1("par_Iso_NSEW_HS_1","Iso Parameter (HS NSEW Sky)",100,16300,16500);
    TH1D par_NS_NSEW_HS_1("par_NS_NSEW_HS_1","NS Parameter (HS NSEW Sky)",100,1800,2200);
    TH1D par_EW_NSEW_HS_1("par_EW_NSEW_HS_1","EW Parameter (HS NSEW Sky)",100,1800,2200);
    TH1D par_FB_NSEW_HS_1("par_FB_NSEW_HS_1","FB Parameter (HS NSEW Sky)",100,-100,100);
    
    TH1D parerr_Iso_NSEW_HS_1("parerr_Iso_NSEW_HS_1","Iso Parameter Error (HS NSEW Sky)",100,-100,100);
    TH1D parerr_NS_NSEW_HS_1("parerr_NS_NSEW_HS_1","NS Parameter Error (HS NSEW Sky)",100,-100,100);
    TH1D parerr_EW_NSEW_HS_1("parerr_EW_NSEW_HS_1","EW Parameter Error (HS NSEW Sky)",100,-100,100);
    TH1D parerr_FB_NSEW_HS_1("parerr_FB_NSEW_HS_1","FB Parameter Error (HS NSEW Sky)",100,-100,100);
    
    ///////// HS NS+FB Dipole Sky
    
    TH1D par_Iso_NSFB_HS_1("par_Iso_NSFB_HS_1","Iso Parameter (HS NSFB Sky)",100,16300,16500);
    TH1D par_NS_NSFB_HS_1("par_NS_NSFB_HS_1","NS Parameter (HS NSFB Sky)",100,1800,2200);
    TH1D par_EW_NSFB_HS_1("par_EW_NSFB_HS_1","EW Parameter (HS NSFB Sky)",100,-100,100);
    TH1D par_FB_NSFB_HS_1("par_FB_NSFB_HS_1","FB Parameter (HS NSFB Sky)",100,1800,2200);
    
    TH1D parerr_Iso_NSFB_HS_1("parerr_Iso_NSFB_HS_1","Iso Parameter Error (HS NSFB Sky)",100,-100,100);
    TH1D parerr_NS_NSFB_HS_1("parerr_NS_NSFB_HS_1","NS Parameter Error (HS NSFB Sky)",100,-100,100);
    TH1D parerr_EW_NSFB_HS_1("parerr_EW_NSFB_HS_1","EW Parameter Error (HS NSFB Sky)",100,-100,100);
    TH1D parerr_FB_NSFB_HS_1("parerr_FB_NSFB_HS_1","FB Parameter Error (HS NSFB Sky)",100,-100,100);
    
    ///////// HS EW+FB Dipole Sky
    
    TH1D par_Iso_EWFB_HS_1("par_Iso_EWFB_HS_1","Iso Parameter (HS EWFB Sky)",100,16300,16500);
    TH1D par_NS_EWFB_HS_1("par_NS_EWFB_HS_1","NS Parameter (HS EWFB Sky)",100,-100,100);
    TH1D par_EW_EWFB_HS_1("par_EW_EWFB_HS_1","EW Parameter (HS EWFB Sky)",100,1800,2200);
    TH1D par_FB_EWFB_HS_1("par_FB_EWFB_HS_1","FB Parameter (HS EWFB Sky)",100,1800,2200);
    
    TH1D parerr_Iso_EWFB_HS_1("parerr_Iso_EWFB_HS_1","Iso Parameter Error (HS EWFB Sky)",100,-100,100);
    TH1D parerr_NS_EWFB_HS_1("parerr_NS_EWFB_HS_1","NS Parameter Error (HS EWFB Sky)",100,-100,100);
    TH1D parerr_EW_EWFB_HS_1("parerr_EW_EWFB_HS_1","EW Parameter Error (HS EWFB Sky)",100,-100,100);
    TH1D parerr_FB_EWFB_HS_1("parerr_FB_EWFB_HS_1","FB Parameter Error (HS EWFB Sky)",100,-100,100);
    
    ///////// HS Full Dipole Sky
    
    TH1D par_Iso_full_HS_1("par_Iso_full_HS_1","Iso Parameter (HS full Sky)",100,16300,16500);
    TH1D par_NS_full_HS_1("par_NS_full_HS_1","NS Parameter (HS full Sky)",100,1800,2200);
    TH1D par_EW_full_HS_1("par_EW_full_HS_1","EW Parameter (HS full Sky)",100,1800,2200);
    TH1D par_FB_full_HS_1("par_FB_full_HS_1","FB Parameter (HS full Sky)",100,1800,2200);
    
    TH1D parerr_Iso_full_HS_1("parerr_Iso_full_HS_1","Iso Parameter Error (HS full Sky)",100,-100,100);
    TH1D parerr_NS_full_HS_1("parerr_NS_full_HS_1","NS Parameter Error (HS full Sky)",100,-100,100);
    TH1D parerr_EW_full_HS_1("parerr_EW_full_HS_1","EW Parameter Error (HS full Sky)",100,-100,100);
    TH1D parerr_FB_full_HS_1("parerr_FB_full_HS_1","FB Parameter Error (HS full Sky)",100,-100,100);
    
    /////////////////////////////////
    
    TH1D chi2_Iso_HS_1("chi2_Iso_HS_1","#chi^2 Iso LS",1000,500,1000);
    TH1D chi2_NS_HS_1("chi2_NS_HS_1","#chi^2 NS LS",1000,500,1000);
    TH1D chi2_EW_HS_1("chi2_EW_HS_1","#chi^2 EW LS",1000,500,1000);
    TH1D chi2_FB_HS_1("chi2_FB_HS_1","#chi^2 FB LS",1000,500,1000);
    TH1D chi2_NSEW_HS_1("chi2_NSEW_HS_1","#chi^2 NS+EW LS",1000,500,1000);
    TH1D chi2_NSFB_HS_1("chi2_NSFB_HS_1","#chi^2 NS+FB LS",1000,500,1000);
    TH1D chi2_EWFB_HS_1("chi2_EWFB_HS_1","#chi^2 EW+FB LS",1000,500,1000);
    TH1D chi2_full_HS_1("chi2_full_HS_1","#chi^2 LS",1000,500,1000);
    
    TH1D ndf_Iso_HS_1("ndf_Iso_HS_1","ndf Iso LS",1000,500,1000);
    TH1D ndf_NS_HS_1("ndf_NS_HS_1","ndf NS LS",1000,500,1000);
    TH1D ndf_EW_HS_1("ndf_EW_HS_1","ndf EW LS",1000,500,1000);
    TH1D ndf_FB_HS_1("ndf_FB_HS_1","ndf FB LS",1000,500,1000);
    TH1D ndf_NSEW_HS_1("ndf_NSEW_HS_1","ndf NS+EW LS",1000,500,1000);
    TH1D ndf_NSFB_HS_1("ndf_NSFB_HS_1","ndf NS+FB LS",1000,500,1000);
    TH1D ndf_EWFB_HS_1("ndf_EWFB_HS_1","ndf EW+FB LS",1000,500,1000);
    TH1D ndf_full_HS_1("ndf_full_HS_1","ndf LS",1000,500,1000);
    
    
    ///////////////////////////// 0.1% anisotropy
    
    ///////// HS Isotropic Sky
    
    TH1D par_Iso_Iso_HS_01("par_Iso_Iso_HS_01","Iso Parameter (HS Full Isotropic Sky)",100,16300,16500);
    TH1D par_NS_Iso_HS_01("par_NS_Iso_HS_01","NS Parameter (HS Full Isotropic Sky)",100,-100,100);
    TH1D par_EW_Iso_HS_01("par_EW_Iso_HS_01","EW Parameter (HS Full Isotropic Sky)",100,-100,100);
    TH1D par_FB_Iso_HS_01("par_FB_Iso_HS_01","FB Parameter (HS Full Isotropic Sky)",100,-100,100);
    
    TH1D parerr_Iso_Iso_HS_01("parerr_Iso_Iso_HS_01","Iso Parameter Error (HS Full Isotropic Sky)",100,-100,100);
    TH1D parerr_NS_Iso_HS_01("parerr_NS_Iso_HS_01","NS Parameter Error (HS Full Isotropic Sky)",100,-100,100);
    TH1D parerr_EW_Iso_HS_01("parerr_EW_Iso_HS_01","EW Parameter Error (HS Full Isotropic Sky)",100,-100,100);
    TH1D parerr_FB_Iso_HS_01("parerr_FB_Iso_HS_01","FB Parameter Error (HS Full Isotropic Sky)",100,-100,100);
    
    ///////// HS NS Dipole Sky
    
    TH1D par_Iso_NS_HS_01("par_Iso_NS_HS_01","Iso Parameter (HS NS Sky)",100,16300,16500);
    TH1D par_NS_NS_HS_01("par_NS_NS_HS_01","NS Parameter (HS NS Sky)",100,-100,100);
    TH1D par_EW_NS_HS_01("par_EW_NS_HS_01","EW Parameter (HS NS Sky)",100,-100,100);
    TH1D par_FB_NS_HS_01("par_FB_NS_HS_01","FB Parameter (HS NS Sky)",100,-100,100);
    
    TH1D parerr_Iso_NS_HS_01("parerr_Iso_NS_HS_01","Iso Parameter (HS NS Sky)",100,-100,100);
    TH1D parerr_NS_NS_HS_01("parerr_NS_NS_HS_01","NS Parameter (HS NS Sky)",100,-100,100);
    TH1D parerr_EW_NS_HS_01("parerr_EW_NS_HS_01","EW Parameter (HS NS Sky)",100,-100,100);
    TH1D parerr_FB_NS_HS_01("parerr_FB_NS_HS_01","FB Parameter (HS NS Sky)",100,-100,100);
    
    ///////// HS EW Dipole Sky
    
    TH1D par_Iso_EW_HS_01("par_Iso_EW_HS_01","Iso Parameter (HS EW Sky)",100,16300,16500);
    TH1D par_NS_EW_HS_01("par_NS_EW_HS_01","NS Parameter (HS EW Sky)",100,-100,100);
    TH1D par_EW_EW_HS_01("par_EW_EW_HS_01","EW Parameter (HS EW Sky)",100,-100,100);
    TH1D par_FB_EW_HS_01("par_FB_EW_HS_01","FB Parameter (HS EW Sky)",100,-100,100);
    
    TH1D parerr_Iso_EW_HS_01("parerr_Iso_EW_HS_01","Iso Parameter Error (HS EW Sky)",100,-100,100);
    TH1D parerr_NS_EW_HS_01("parerr_NS_EW_HS_01","NS Parameter Error (HS EW Sky)",100,-100,100);
    TH1D parerr_EW_EW_HS_01("parerr_EW_EW_HS_01","EW Parameter Error (HS EW Sky)",100,-100,100);
    TH1D parerr_FB_EW_HS_01("parerr_FB_EW_HS_01","FB Parameter Error (HS EW Sky)",100,-100,100);
    
    ///////// HS FB Dipole Sky
    
    TH1D par_Iso_FB_HS_01("par_Iso_FB_HS_01","Iso Parameter (HS FB Sky)",100,16300,16500);
    TH1D par_NS_FB_HS_01("par_NS_FB_HS_01","NS Parameter (HS FB Sky)",100,-100,100);
    TH1D par_EW_FB_HS_01("par_EW_FB_HS_01","EW Parameter (HS FB Sky)",100,-100,100);
    TH1D par_FB_FB_HS_01("par_FB_FB_HS_01","FB Parameter (HS FB Sky)",100,-100,100);
    
    TH1D parerr_Iso_FB_HS_01("parerr_Iso_FB_HS_01","Iso Parameter Error (HS FB Sky)",100,-100,100);
    TH1D parerr_NS_FB_HS_01("parerr_NS_FB_HS_01","NS Parameter Error (HS FB Sky)",100,-100,100);
    TH1D parerr_EW_FB_HS_01("parerr_EW_FB_HS_01","EW Parameter Error (HS FB Sky)",100,-100,100);
    TH1D parerr_FB_FB_HS_01("parerr_FB_FB_HS_01","FB Parameter Error (HS FB Sky)",100,-100,100);
    
    ///////// HS NS+EW Dipole Sky
    
    TH1D par_Iso_NSEW_HS_01("par_Iso_NSEW_HS_01","Iso Parameter (HS NSEW Sky)",100,16300,16500);
    TH1D par_NS_NSEW_HS_01("par_NS_NSEW_HS_01","NS Parameter (HS NSEW Sky)",100,-100,100);
    TH1D par_EW_NSEW_HS_01("par_EW_NSEW_HS_01","EW Parameter (HS NSEW Sky)",100,-100,100);
    TH1D par_FB_NSEW_HS_01("par_FB_NSEW_HS_01","FB Parameter (HS NSEW Sky)",100,-100,100);
    
    TH1D parerr_Iso_NSEW_HS_01("parerr_Iso_NSEW_HS_01","Iso Parameter Error (HS NSEW Sky)",100,-100,100);
    TH1D parerr_NS_NSEW_HS_01("parerr_NS_NSEW_HS_01","NS Parameter Error (HS NSEW Sky)",100,-100,100);
    TH1D parerr_EW_NSEW_HS_01("parerr_EW_NSEW_HS_01","EW Parameter Error (HS NSEW Sky)",100,-100,100);
    TH1D parerr_FB_NSEW_HS_01("parerr_FB_NSEW_HS_01","FB Parameter Error (HS NSEW Sky)",100,-100,100);
    
    ///////// HS NS+FB Dipole Sky
    
    TH1D par_Iso_NSFB_HS_01("par_Iso_NSFB_HS_01","Iso Parameter (HS NSFB Sky)",100,16300,16500);
    TH1D par_NS_NSFB_HS_01("par_NS_NSFB_HS_01","NS Parameter (HS NSFB Sky)",100,-100,100);
    TH1D par_EW_NSFB_HS_01("par_EW_NSFB_HS_01","EW Parameter (HS NSFB Sky)",100,-100,100);
    TH1D par_FB_NSFB_HS_01("par_FB_NSFB_HS_01","FB Parameter (HS NSFB Sky)",100,-100,100);
    
    TH1D parerr_Iso_NSFB_HS_01("parerr_Iso_NSFB_HS_01","Iso Parameter Error (HS NSFB Sky)",100,-100,100);
    TH1D parerr_NS_NSFB_HS_01("parerr_NS_NSFB_HS_01","NS Parameter Error (HS NSFB Sky)",100,-100,100);
    TH1D parerr_EW_NSFB_HS_01("parerr_EW_NSFB_HS_01","EW Parameter Error (HS NSFB Sky)",100,-100,100);
    TH1D parerr_FB_NSFB_HS_01("parerr_FB_NSFB_HS_01","FB Parameter Error (HS NSFB Sky)",100,-100,100);
    
    ///////// HS EW+FB Dipole Sky
    
    TH1D par_Iso_EWFB_HS_01("par_Iso_EWFB_HS_01","Iso Parameter (HS EWFB Sky)",100,16300,16500);
    TH1D par_NS_EWFB_HS_01("par_NS_EWFB_HS_01","NS Parameter (HS EWFB Sky)",100,-100,100);
    TH1D par_EW_EWFB_HS_01("par_EW_EWFB_HS_01","EW Parameter (HS EWFB Sky)",100,-100,100);
    TH1D par_FB_EWFB_HS_01("par_FB_EWFB_HS_01","FB Parameter (HS EWFB Sky)",100,-100,100);
    
    TH1D parerr_Iso_EWFB_HS_01("parerr_Iso_EWFB_HS_01","Iso Parameter Error (HS EWFB Sky)",100,-100,100);
    TH1D parerr_NS_EWFB_HS_01("parerr_NS_EWFB_HS_01","NS Parameter Error (HS EWFB Sky)",100,-100,100);
    TH1D parerr_EW_EWFB_HS_01("parerr_EW_EWFB_HS_01","EW Parameter Error (HS EWFB Sky)",100,-100,100);
    TH1D parerr_FB_EWFB_HS_01("parerr_FB_EWFB_HS_01","FB Parameter Error (HS EWFB Sky)",100,-100,100);
    
    ///////// HS Full Dipole Sky
    
    TH1D par_Iso_full_HS_01("par_Iso_full_HS_01","Iso Parameter (HS full Sky)",100,16300,16500);
    TH1D par_NS_full_HS_01("par_NS_full_HS_01","NS Parameter (HS full Sky)",100,-100,100);
    TH1D par_EW_full_HS_01("par_EW_full_HS_01","EW Parameter (HS full Sky)",100,-100,100);
    TH1D par_FB_full_HS_01("par_FB_full_HS_01","FB Parameter (HS full Sky)",100,-100,100);
    
    TH1D parerr_Iso_full_HS_01("parerr_Iso_full_HS_01","Iso Parameter Error (HS full Sky)",100,-100,100);
    TH1D parerr_NS_full_HS_01("parerr_NS_full_HS_01","NS Parameter Error (HS full Sky)",100,-100,100);
    TH1D parerr_EW_full_HS_01("parerr_EW_full_HS_01","EW Parameter Error (HS full Sky)",100,-100,100);
    TH1D parerr_FB_full_HS_01("parerr_FB_full_HS_01","FB Parameter Error (HS full Sky)",100,-100,100);
    
    /////////////////////////////////
    
    TH1D chi2_Iso_HS_01("chi2_Iso_HS_01","#chi^2 Iso LS",1000,500,1000);
    TH1D chi2_NS_HS_01("chi2_NS_HS_01","#chi^2 NS LS",1000,500,1000);
    TH1D chi2_EW_HS_01("chi2_EW_HS_01","#chi^2 EW LS",1000,500,1000);
    TH1D chi2_FB_HS_01("chi2_FB_HS_01","#chi^2 FB LS",1000,500,1000);
    TH1D chi2_NSEW_HS_01("chi2_NSEW_HS_01","#chi^2 NS+EW LS",1000,500,1000);
    TH1D chi2_NSFB_HS_01("chi2_NSFB_HS_01","#chi^2 NS+FB LS",1000,500,1000);
    TH1D chi2_EWFB_HS_01("chi2_EWFB_HS_01","#chi^2 EW+FB LS",1000,500,1000);
    TH1D chi2_full_HS_01("chi2_full_HS_01","#chi^2 LS",1000,500,1000);
    
    TH1D ndf_Iso_HS_01("ndf_Iso_HS_01","ndf Iso LS",1000,500,1000);
    TH1D ndf_NS_HS_01("ndf_NS_HS_01","ndf NS LS",1000,500,1000);
    TH1D ndf_EW_HS_01("ndf_EW_HS_01","ndf EW LS",1000,500,1000);
    TH1D ndf_FB_HS_01("ndf_FB_HS_01","ndf FB LS",1000,500,1000);
    TH1D ndf_NSEW_HS_01("ndf_NSEW_HS_01","ndf NS+EW LS",1000,500,1000);
    TH1D ndf_NSFB_HS_01("ndf_NSFB_HS_01","ndf NS+FB LS",1000,500,1000);
    TH1D ndf_EWFB_HS_01("ndf_EWFB_HS_01","ndf EW+FB LS",1000,500,1000);
    TH1D ndf_full_HS_01("ndf_full_HS_01","ndf LS",1000,500,1000);
    
    
    ////////////////////////////////////////////////////////////////////

    
    /////////////////////////////////////////////////////////////////////////////////////////////
    
    for(UInt_t tree_idx = 0; tree_idx < tree_entries; ++tree_idx)
    {
        myTree->GetEntry(tree_idx);
        
        if(tmp_result.inputAni[0] == 0.1)
        {
            
            ////////////// Filling parameter LS histos
            
            par_Iso_Iso_LS_10.Fill(tmp_result.fit_par_LS[0][0]);
            par_NS_Iso_LS_10.Fill(tmp_result.fit_par_LS[0][1]);
            par_EW_Iso_LS_10.Fill(tmp_result.fit_par_LS[0][2]);
            par_FB_Iso_LS_10.Fill(tmp_result.fit_par_LS[0][3]);
            
            par_Iso_Iso_LS_1.Fill(tmp_result.fit_par_LS[0][0]);
            par_NS_Iso_LS_1.Fill(tmp_result.fit_par_LS[0][1]);
            par_EW_Iso_LS_1.Fill(tmp_result.fit_par_LS[0][2]);
            par_FB_Iso_LS_1.Fill(tmp_result.fit_par_LS[0][3]);
            
            par_Iso_Iso_LS_01.Fill(tmp_result.fit_par_LS[0][0]);
            par_NS_Iso_LS_01.Fill(tmp_result.fit_par_LS[0][1]);
            par_EW_Iso_LS_01.Fill(tmp_result.fit_par_LS[0][2]);
            par_FB_Iso_LS_01.Fill(tmp_result.fit_par_LS[0][3]);
            
            par_Iso_NS_LS_10.Fill(tmp_result.fit_par_LS[1][0]);
            par_NS_NS_LS_10.Fill(tmp_result.fit_par_LS[1][1]);
            par_EW_NS_LS_10.Fill(tmp_result.fit_par_LS[1][2]);
            par_FB_NS_LS_10.Fill(tmp_result.fit_par_LS[1][3]);
            
            par_Iso_EW_LS_10.Fill(tmp_result.fit_par_LS[2][0]);
            par_NS_EW_LS_10.Fill(tmp_result.fit_par_LS[2][1]);
            par_EW_EW_LS_10.Fill(tmp_result.fit_par_LS[2][2]);
            par_FB_EW_LS_10.Fill(tmp_result.fit_par_LS[2][3]);
            
            par_Iso_FB_LS_10.Fill(tmp_result.fit_par_LS[3][0]);
            par_NS_FB_LS_10.Fill(tmp_result.fit_par_LS[3][1]);
            par_EW_FB_LS_10.Fill(tmp_result.fit_par_LS[3][2]);
            par_FB_FB_LS_10.Fill(tmp_result.fit_par_LS[3][3]);
            
            par_Iso_NSEW_LS_10.Fill(tmp_result.fit_par_LS[4][0]);
            par_NS_NSEW_LS_10.Fill(tmp_result.fit_par_LS[4][1]);
            par_EW_NSEW_LS_10.Fill(tmp_result.fit_par_LS[4][2]);
            par_FB_NSEW_LS_10.Fill(tmp_result.fit_par_LS[4][3]);
            
            par_Iso_NSFB_LS_10.Fill(tmp_result.fit_par_LS[5][0]);
            par_NS_NSFB_LS_10.Fill(tmp_result.fit_par_LS[5][1]);
            par_EW_NSFB_LS_10.Fill(tmp_result.fit_par_LS[5][2]);
            par_FB_NSFB_LS_10.Fill(tmp_result.fit_par_LS[5][3]);
            
            par_Iso_EWFB_LS_10.Fill(tmp_result.fit_par_LS[6][0]);
            par_NS_EWFB_LS_10.Fill(tmp_result.fit_par_LS[6][1]);
            par_EW_EWFB_LS_10.Fill(tmp_result.fit_par_LS[6][2]);
            par_FB_EWFB_LS_10.Fill(tmp_result.fit_par_LS[6][3]);
            
            par_Iso_full_LS_10.Fill(tmp_result.fit_par_LS[7][0]);
            par_NS_full_LS_10.Fill(tmp_result.fit_par_LS[7][1]);
            par_EW_full_LS_10.Fill(tmp_result.fit_par_LS[7][2]);
            par_FB_full_LS_10.Fill(tmp_result.fit_par_LS[7][3]);
            
            ////////////// Filling error LS histos
            
            parerr_Iso_Iso_LS_10.Fill(tmp_result.fit_err_LS[0][0]);
            parerr_NS_Iso_LS_10.Fill(tmp_result.fit_err_LS[0][1]);
            parerr_EW_Iso_LS_10.Fill(tmp_result.fit_err_LS[0][2]);
            parerr_FB_Iso_LS_10.Fill(tmp_result.fit_err_LS[0][3]);
            
            parerr_Iso_Iso_LS_1.Fill(tmp_result.fit_err_LS[0][0]);
            parerr_NS_Iso_LS_1.Fill(tmp_result.fit_err_LS[0][1]);
            parerr_EW_Iso_LS_1.Fill(tmp_result.fit_err_LS[0][2]);
            parerr_FB_Iso_LS_1.Fill(tmp_result.fit_err_LS[0][3]);
            
            parerr_Iso_Iso_LS_01.Fill(tmp_result.fit_err_LS[0][0]);
            parerr_NS_Iso_LS_01.Fill(tmp_result.fit_err_LS[0][1]);
            parerr_EW_Iso_LS_01.Fill(tmp_result.fit_err_LS[0][2]);
            parerr_FB_Iso_LS_01.Fill(tmp_result.fit_err_LS[0][3]);
            
            parerr_Iso_NS_LS_10.Fill(tmp_result.fit_err_LS[1][0]);
            parerr_NS_NS_LS_10.Fill(tmp_result.fit_err_LS[1][1]);
            parerr_EW_NS_LS_10.Fill(tmp_result.fit_err_LS[1][2]);
            parerr_FB_NS_LS_10.Fill(tmp_result.fit_err_LS[1][3]);
            
            parerr_Iso_EW_LS_10.Fill(tmp_result.fit_err_LS[2][0]);
            parerr_NS_EW_LS_10.Fill(tmp_result.fit_err_LS[2][1]);
            parerr_EW_EW_LS_10.Fill(tmp_result.fit_err_LS[2][2]);
            parerr_FB_EW_LS_10.Fill(tmp_result.fit_err_LS[2][3]);
            
            parerr_Iso_FB_LS_10.Fill(tmp_result.fit_err_LS[3][0]);
            parerr_NS_FB_LS_10.Fill(tmp_result.fit_err_LS[3][1]);
            parerr_EW_FB_LS_10.Fill(tmp_result.fit_err_LS[3][2]);
            parerr_FB_FB_LS_10.Fill(tmp_result.fit_err_LS[3][3]);
            
            parerr_Iso_NSEW_LS_10.Fill(tmp_result.fit_err_LS[4][0]);
            parerr_NS_NSEW_LS_10.Fill(tmp_result.fit_err_LS[4][1]);
            parerr_EW_NSEW_LS_10.Fill(tmp_result.fit_err_LS[4][2]);
            parerr_FB_NSEW_LS_10.Fill(tmp_result.fit_err_LS[4][3]);
            
            parerr_Iso_NSFB_LS_10.Fill(tmp_result.fit_err_LS[5][0]);
            parerr_NS_NSFB_LS_10.Fill(tmp_result.fit_err_LS[5][1]);
            parerr_EW_NSFB_LS_10.Fill(tmp_result.fit_err_LS[5][2]);
            parerr_FB_NSFB_LS_10.Fill(tmp_result.fit_err_LS[5][3]);
            
            parerr_Iso_EWFB_LS_10.Fill(tmp_result.fit_err_LS[6][0]);
            parerr_NS_EWFB_LS_10.Fill(tmp_result.fit_err_LS[6][1]);
            parerr_EW_EWFB_LS_10.Fill(tmp_result.fit_err_LS[6][2]);
            parerr_FB_EWFB_LS_10.Fill(tmp_result.fit_err_LS[6][3]);
            
            parerr_Iso_full_LS_10.Fill(tmp_result.fit_err_LS[7][0]);
            parerr_NS_full_LS_10.Fill(tmp_result.fit_err_LS[7][1]);
            parerr_EW_full_LS_10.Fill(tmp_result.fit_err_LS[7][2]);
            parerr_FB_full_LS_10.Fill(tmp_result.fit_err_LS[7][3]);
            
            ////////////// Filling chi2 LS histos
            
            chi2_Iso_LS_10.Fill(tmp_result.chi2_LS[0]);
            chi2_NS_LS_10.Fill(tmp_result.chi2_LS[1]);
            chi2_EW_LS_10.Fill(tmp_result.chi2_LS[2]);
            chi2_FB_LS_10.Fill(tmp_result.chi2_LS[3]);
            chi2_NSEW_LS_10.Fill(tmp_result.chi2_LS[4]);
            chi2_NSFB_LS_10.Fill(tmp_result.chi2_LS[5]);
            chi2_EWFB_LS_10.Fill(tmp_result.chi2_LS[6]);
            chi2_full_LS_10.Fill(tmp_result.chi2_LS[7]);
            
            ////////////// Filling ndf LS histos
            
            ndf_Iso_LS_10.Fill(tmp_result.ndf_LS[0]);
            ndf_NS_LS_10.Fill(tmp_result.ndf_LS[1]);
            ndf_EW_LS_10.Fill(tmp_result.ndf_LS[2]);
            ndf_FB_LS_10.Fill(tmp_result.ndf_LS[3]);
            ndf_NSEW_LS_10.Fill(tmp_result.ndf_LS[4]);
            ndf_NSFB_LS_10.Fill(tmp_result.ndf_LS[5]);
            ndf_EWFB_LS_10.Fill(tmp_result.ndf_LS[6]);
            ndf_full_LS_10.Fill(tmp_result.ndf_LS[7]);
            
            
            
            
            ////////////// Filling parameter HS histos
            
            par_Iso_Iso_HS_10.Fill(tmp_result.fit_par_HS[0][0]);
            par_NS_Iso_HS_10.Fill(tmp_result.fit_par_HS[0][1]);
            par_EW_Iso_HS_10.Fill(tmp_result.fit_par_HS[0][2]);
            par_FB_Iso_HS_10.Fill(tmp_result.fit_par_HS[0][3]);
            
            par_Iso_Iso_HS_1.Fill(tmp_result.fit_par_HS[0][0]);
            par_NS_Iso_HS_1.Fill(tmp_result.fit_par_HS[0][1]);
            par_EW_Iso_HS_1.Fill(tmp_result.fit_par_HS[0][2]);
            par_FB_Iso_HS_1.Fill(tmp_result.fit_par_HS[0][3]);
            
            par_Iso_Iso_HS_01.Fill(tmp_result.fit_par_HS[0][0]);
            par_NS_Iso_HS_01.Fill(tmp_result.fit_par_HS[0][1]);
            par_EW_Iso_HS_01.Fill(tmp_result.fit_par_HS[0][2]);
            par_FB_Iso_HS_01.Fill(tmp_result.fit_par_HS[0][3]);
  
            par_Iso_NS_HS_10.Fill(tmp_result.fit_par_HS[1][0]);
            par_NS_NS_HS_10.Fill(tmp_result.fit_par_HS[1][1]);
            par_EW_NS_HS_10.Fill(tmp_result.fit_par_HS[1][2]);
            par_FB_NS_HS_10.Fill(tmp_result.fit_par_HS[1][3]);
            
            par_Iso_EW_HS_10.Fill(tmp_result.fit_par_HS[2][0]);
            par_NS_EW_HS_10.Fill(tmp_result.fit_par_HS[2][1]);
            par_EW_EW_HS_10.Fill(tmp_result.fit_par_HS[2][2]);
            par_FB_EW_HS_10.Fill(tmp_result.fit_par_HS[2][3]);
            
            par_Iso_FB_HS_10.Fill(tmp_result.fit_par_HS[3][0]);
            par_NS_FB_HS_10.Fill(tmp_result.fit_par_HS[3][1]);
            par_EW_FB_HS_10.Fill(tmp_result.fit_par_HS[3][2]);
            par_FB_FB_HS_10.Fill(tmp_result.fit_par_HS[3][3]);
            
            par_Iso_NSEW_HS_10.Fill(tmp_result.fit_par_HS[4][0]);
            par_NS_NSEW_HS_10.Fill(tmp_result.fit_par_HS[4][1]);
            par_EW_NSEW_HS_10.Fill(tmp_result.fit_par_HS[4][2]);
            par_FB_NSEW_HS_10.Fill(tmp_result.fit_par_HS[4][3]);
            
            par_Iso_NSFB_HS_10.Fill(tmp_result.fit_par_HS[5][0]);
            par_NS_NSFB_HS_10.Fill(tmp_result.fit_par_HS[5][1]);
            par_EW_NSFB_HS_10.Fill(tmp_result.fit_par_HS[5][2]);
            par_FB_NSFB_HS_10.Fill(tmp_result.fit_par_HS[5][3]);
            
            par_Iso_EWFB_HS_10.Fill(tmp_result.fit_par_HS[6][0]);
            par_NS_EWFB_HS_10.Fill(tmp_result.fit_par_HS[6][1]);
            par_EW_EWFB_HS_10.Fill(tmp_result.fit_par_HS[6][2]);
            par_FB_EWFB_HS_10.Fill(tmp_result.fit_par_HS[6][3]);
            
            par_Iso_full_HS_10.Fill(tmp_result.fit_par_HS[7][0]);
            par_NS_full_HS_10.Fill(tmp_result.fit_par_HS[7][1]);
            par_EW_full_HS_10.Fill(tmp_result.fit_par_HS[7][2]);
            par_FB_full_HS_10.Fill(tmp_result.fit_par_HS[7][3]);
            
            ////////////// Filling error HS histos
            
            parerr_Iso_Iso_HS_10.Fill(tmp_result.fit_err_HS[0][0]);
            parerr_NS_Iso_HS_10.Fill(tmp_result.fit_err_HS[0][1]);
            parerr_EW_Iso_HS_10.Fill(tmp_result.fit_err_HS[0][2]);
            parerr_FB_Iso_HS_10.Fill(tmp_result.fit_err_HS[0][3]);
            
            parerr_Iso_Iso_HS_1.Fill(tmp_result.fit_err_HS[0][0]);
            parerr_NS_Iso_HS_1.Fill(tmp_result.fit_err_HS[0][1]);
            parerr_EW_Iso_HS_1.Fill(tmp_result.fit_err_HS[0][2]);
            parerr_FB_Iso_HS_1.Fill(tmp_result.fit_err_HS[0][3]);
            
            parerr_Iso_Iso_HS_01.Fill(tmp_result.fit_err_HS[0][0]);
            parerr_NS_Iso_HS_01.Fill(tmp_result.fit_err_HS[0][1]);
            parerr_EW_Iso_HS_01.Fill(tmp_result.fit_err_HS[0][2]);
            parerr_FB_Iso_HS_01.Fill(tmp_result.fit_err_HS[0][3]);
            
            parerr_Iso_NS_HS_10.Fill(tmp_result.fit_err_HS[1][0]);
            parerr_NS_NS_HS_10.Fill(tmp_result.fit_err_HS[1][1]);
            parerr_EW_NS_HS_10.Fill(tmp_result.fit_err_HS[1][2]);
            parerr_FB_NS_HS_10.Fill(tmp_result.fit_err_HS[1][3]);
            
            parerr_Iso_EW_HS_10.Fill(tmp_result.fit_err_HS[2][0]);
            parerr_NS_EW_HS_10.Fill(tmp_result.fit_err_HS[2][1]);
            parerr_EW_EW_HS_10.Fill(tmp_result.fit_err_HS[2][2]);
            parerr_FB_EW_HS_10.Fill(tmp_result.fit_err_HS[2][3]);
            
            parerr_Iso_FB_HS_10.Fill(tmp_result.fit_err_HS[3][0]);
            parerr_NS_FB_HS_10.Fill(tmp_result.fit_err_HS[3][1]);
            parerr_EW_FB_HS_10.Fill(tmp_result.fit_err_HS[3][2]);
            parerr_FB_FB_HS_10.Fill(tmp_result.fit_err_HS[3][3]);
            
            parerr_Iso_NSEW_HS_10.Fill(tmp_result.fit_err_HS[4][0]);
            parerr_NS_NSEW_HS_10.Fill(tmp_result.fit_err_HS[4][1]);
            parerr_EW_NSEW_HS_10.Fill(tmp_result.fit_err_HS[4][2]);
            parerr_FB_NSEW_HS_10.Fill(tmp_result.fit_err_HS[4][3]);
            
            parerr_Iso_NSFB_HS_10.Fill(tmp_result.fit_err_HS[5][0]);
            parerr_NS_NSFB_HS_10.Fill(tmp_result.fit_err_HS[5][1]);
            parerr_EW_NSFB_HS_10.Fill(tmp_result.fit_err_HS[5][2]);
            parerr_FB_NSFB_HS_10.Fill(tmp_result.fit_err_HS[5][3]);
            
            parerr_Iso_EWFB_HS_10.Fill(tmp_result.fit_err_HS[6][0]);
            parerr_NS_EWFB_HS_10.Fill(tmp_result.fit_err_HS[6][1]);
            parerr_EW_EWFB_HS_10.Fill(tmp_result.fit_err_HS[6][2]);
            parerr_FB_EWFB_HS_10.Fill(tmp_result.fit_err_HS[6][3]);
            
            parerr_Iso_full_HS_10.Fill(tmp_result.fit_err_HS[7][0]);
            parerr_NS_full_HS_10.Fill(tmp_result.fit_err_HS[7][1]);
            parerr_EW_full_HS_10.Fill(tmp_result.fit_err_HS[7][2]);
            parerr_FB_full_HS_10.Fill(tmp_result.fit_err_HS[7][3]);
            
            ////////////// Filling chi2 HS histos
            
            chi2_Iso_HS_10.Fill(tmp_result.chi2_HS[0]);
            chi2_NS_HS_10.Fill(tmp_result.chi2_HS[1]);
            chi2_EW_HS_10.Fill(tmp_result.chi2_HS[2]);
            chi2_FB_HS_10.Fill(tmp_result.chi2_HS[3]);
            chi2_NSEW_HS_10.Fill(tmp_result.chi2_HS[4]);
            chi2_NSFB_HS_10.Fill(tmp_result.chi2_HS[5]);
            chi2_EWFB_HS_10.Fill(tmp_result.chi2_HS[6]);
            chi2_full_HS_10.Fill(tmp_result.chi2_HS[7]);
            
            ////////////// Filling ndf HS histos
            
            ndf_Iso_HS_10.Fill(tmp_result.ndf_HS[0]);
            ndf_NS_HS_10.Fill(tmp_result.ndf_HS[1]);
            ndf_EW_HS_10.Fill(tmp_result.ndf_HS[2]);
            ndf_FB_HS_10.Fill(tmp_result.ndf_HS[3]);
            ndf_NSEW_HS_10.Fill(tmp_result.ndf_HS[4]);
            ndf_NSFB_HS_10.Fill(tmp_result.ndf_HS[5]);
            ndf_EWFB_HS_10.Fill(tmp_result.ndf_HS[6]);
            ndf_full_HS_10.Fill(tmp_result.ndf_HS[7]);

            
        }
        
        else if(tmp_result.inputAni[0] == 0.01)
        {
            
            ////////////// Filling parameter LS histos
            
            par_Iso_Iso_LS_10.Fill(tmp_result.fit_par_LS[0][0]);
            par_NS_Iso_LS_10.Fill(tmp_result.fit_par_LS[0][1]);
            par_EW_Iso_LS_10.Fill(tmp_result.fit_par_LS[0][2]);
            par_FB_Iso_LS_10.Fill(tmp_result.fit_par_LS[0][3]);
            
            par_Iso_Iso_LS_1.Fill(tmp_result.fit_par_LS[0][0]);
            par_NS_Iso_LS_1.Fill(tmp_result.fit_par_LS[0][1]);
            par_EW_Iso_LS_1.Fill(tmp_result.fit_par_LS[0][2]);
            par_FB_Iso_LS_1.Fill(tmp_result.fit_par_LS[0][3]);
            
            par_Iso_Iso_LS_01.Fill(tmp_result.fit_par_LS[0][0]);
            par_NS_Iso_LS_01.Fill(tmp_result.fit_par_LS[0][1]);
            par_EW_Iso_LS_01.Fill(tmp_result.fit_par_LS[0][2]);
            par_FB_Iso_LS_01.Fill(tmp_result.fit_par_LS[0][3]);
            
            par_Iso_NS_LS_1.Fill(tmp_result.fit_par_LS[1][0]);
            par_NS_NS_LS_1.Fill(tmp_result.fit_par_LS[1][1]);
            par_EW_NS_LS_1.Fill(tmp_result.fit_par_LS[1][2]);
            par_FB_NS_LS_1.Fill(tmp_result.fit_par_LS[1][3]);
            
            par_Iso_EW_LS_1.Fill(tmp_result.fit_par_LS[2][0]);
            par_NS_EW_LS_1.Fill(tmp_result.fit_par_LS[2][1]);
            par_EW_EW_LS_1.Fill(tmp_result.fit_par_LS[2][2]);
            par_FB_EW_LS_1.Fill(tmp_result.fit_par_LS[2][3]);
            
            par_Iso_FB_LS_1.Fill(tmp_result.fit_par_LS[3][0]);
            par_NS_FB_LS_1.Fill(tmp_result.fit_par_LS[3][1]);
            par_EW_FB_LS_1.Fill(tmp_result.fit_par_LS[3][2]);
            par_FB_FB_LS_1.Fill(tmp_result.fit_par_LS[3][3]);
            
            par_Iso_NSEW_LS_1.Fill(tmp_result.fit_par_LS[4][0]);
            par_NS_NSEW_LS_1.Fill(tmp_result.fit_par_LS[4][1]);
            par_EW_NSEW_LS_1.Fill(tmp_result.fit_par_LS[4][2]);
            par_FB_NSEW_LS_1.Fill(tmp_result.fit_par_LS[4][3]);
            
            par_Iso_NSFB_LS_1.Fill(tmp_result.fit_par_LS[5][0]);
            par_NS_NSFB_LS_1.Fill(tmp_result.fit_par_LS[5][1]);
            par_EW_NSFB_LS_1.Fill(tmp_result.fit_par_LS[5][2]);
            par_FB_NSFB_LS_1.Fill(tmp_result.fit_par_LS[5][3]);
            
            par_Iso_EWFB_LS_1.Fill(tmp_result.fit_par_LS[6][0]);
            par_NS_EWFB_LS_1.Fill(tmp_result.fit_par_LS[6][1]);
            par_EW_EWFB_LS_1.Fill(tmp_result.fit_par_LS[6][2]);
            par_FB_EWFB_LS_1.Fill(tmp_result.fit_par_LS[6][3]);
            
            par_Iso_full_LS_1.Fill(tmp_result.fit_par_LS[7][0]);
            par_NS_full_LS_1.Fill(tmp_result.fit_par_LS[7][1]);
            par_EW_full_LS_1.Fill(tmp_result.fit_par_LS[7][2]);
            par_FB_full_LS_1.Fill(tmp_result.fit_par_LS[7][3]);
            
            ////////////// Filling error LS histos
            
            parerr_Iso_Iso_LS_10.Fill(tmp_result.fit_err_LS[0][0]);
            parerr_NS_Iso_LS_10.Fill(tmp_result.fit_err_LS[0][1]);
            parerr_EW_Iso_LS_10.Fill(tmp_result.fit_err_LS[0][2]);
            parerr_FB_Iso_LS_10.Fill(tmp_result.fit_err_LS[0][3]);
            
            parerr_Iso_Iso_LS_1.Fill(tmp_result.fit_err_LS[0][0]);
            parerr_NS_Iso_LS_1.Fill(tmp_result.fit_err_LS[0][1]);
            parerr_EW_Iso_LS_1.Fill(tmp_result.fit_err_LS[0][2]);
            parerr_FB_Iso_LS_1.Fill(tmp_result.fit_err_LS[0][3]);
            
            parerr_Iso_Iso_LS_01.Fill(tmp_result.fit_err_LS[0][0]);
            parerr_NS_Iso_LS_01.Fill(tmp_result.fit_err_LS[0][1]);
            parerr_EW_Iso_LS_01.Fill(tmp_result.fit_err_LS[0][2]);
            parerr_FB_Iso_LS_01.Fill(tmp_result.fit_err_LS[0][3]);
            
            parerr_Iso_NS_LS_1.Fill(tmp_result.fit_err_LS[1][0]);
            parerr_NS_NS_LS_1.Fill(tmp_result.fit_err_LS[1][1]);
            parerr_EW_NS_LS_1.Fill(tmp_result.fit_err_LS[1][2]);
            parerr_FB_NS_LS_1.Fill(tmp_result.fit_err_LS[1][3]);
            
            parerr_Iso_EW_LS_1.Fill(tmp_result.fit_err_LS[2][0]);
            parerr_NS_EW_LS_1.Fill(tmp_result.fit_err_LS[2][1]);
            parerr_EW_EW_LS_1.Fill(tmp_result.fit_err_LS[2][2]);
            parerr_FB_EW_LS_1.Fill(tmp_result.fit_err_LS[2][3]);
            
            parerr_Iso_FB_LS_1.Fill(tmp_result.fit_err_LS[3][0]);
            parerr_NS_FB_LS_1.Fill(tmp_result.fit_err_LS[3][1]);
            parerr_EW_FB_LS_1.Fill(tmp_result.fit_err_LS[3][2]);
            parerr_FB_FB_LS_1.Fill(tmp_result.fit_err_LS[3][3]);
            
            parerr_Iso_NSEW_LS_1.Fill(tmp_result.fit_err_LS[4][0]);
            parerr_NS_NSEW_LS_1.Fill(tmp_result.fit_err_LS[4][1]);
            parerr_EW_NSEW_LS_1.Fill(tmp_result.fit_err_LS[4][2]);
            parerr_FB_NSEW_LS_1.Fill(tmp_result.fit_err_LS[4][3]);
            
            parerr_Iso_NSFB_LS_1.Fill(tmp_result.fit_err_LS[5][0]);
            parerr_NS_NSFB_LS_1.Fill(tmp_result.fit_err_LS[5][1]);
            parerr_EW_NSFB_LS_1.Fill(tmp_result.fit_err_LS[5][2]);
            parerr_FB_NSFB_LS_1.Fill(tmp_result.fit_err_LS[5][3]);
            
            parerr_Iso_EWFB_LS_1.Fill(tmp_result.fit_err_LS[6][0]);
            parerr_NS_EWFB_LS_1.Fill(tmp_result.fit_err_LS[6][1]);
            parerr_EW_EWFB_LS_1.Fill(tmp_result.fit_err_LS[6][2]);
            parerr_FB_EWFB_LS_1.Fill(tmp_result.fit_err_LS[6][3]);
            
            parerr_Iso_full_LS_1.Fill(tmp_result.fit_err_LS[7][0]);
            parerr_NS_full_LS_1.Fill(tmp_result.fit_err_LS[7][1]);
            parerr_EW_full_LS_1.Fill(tmp_result.fit_err_LS[7][2]);
            parerr_FB_full_LS_1.Fill(tmp_result.fit_err_LS[7][3]);

            ////////////// Filling chi2 LS histos
            
            chi2_Iso_LS_1.Fill(tmp_result.chi2_LS[0]);
            chi2_NS_LS_1.Fill(tmp_result.chi2_LS[1]);
            chi2_EW_LS_1.Fill(tmp_result.chi2_LS[2]);
            chi2_FB_LS_1.Fill(tmp_result.chi2_LS[3]);
            chi2_NSEW_LS_1.Fill(tmp_result.chi2_LS[4]);
            chi2_NSFB_LS_1.Fill(tmp_result.chi2_LS[5]);
            chi2_EWFB_LS_1.Fill(tmp_result.chi2_LS[6]);
            chi2_full_LS_1.Fill(tmp_result.chi2_LS[7]);
            
            ////////////// Filling ndf LS histos
            
            ndf_Iso_LS_1.Fill(tmp_result.ndf_LS[0]);
            ndf_NS_LS_1.Fill(tmp_result.ndf_LS[1]);
            ndf_EW_LS_1.Fill(tmp_result.ndf_LS[2]);
            ndf_FB_LS_1.Fill(tmp_result.ndf_LS[3]);
            ndf_NSEW_LS_1.Fill(tmp_result.ndf_LS[4]);
            ndf_NSFB_LS_1.Fill(tmp_result.ndf_LS[5]);
            ndf_EWFB_LS_1.Fill(tmp_result.ndf_LS[6]);
            ndf_full_LS_1.Fill(tmp_result.ndf_LS[7]);
            
            
            
            
            
            
            ////////////// Filling parameter HS histos
            
            par_Iso_Iso_HS_10.Fill(tmp_result.fit_par_HS[0][0]);
            par_NS_Iso_HS_10.Fill(tmp_result.fit_par_HS[0][1]);
            par_EW_Iso_HS_10.Fill(tmp_result.fit_par_HS[0][2]);
            par_FB_Iso_HS_10.Fill(tmp_result.fit_par_HS[0][3]);
            
            par_Iso_Iso_HS_1.Fill(tmp_result.fit_par_HS[0][0]);
            par_NS_Iso_HS_1.Fill(tmp_result.fit_par_HS[0][1]);
            par_EW_Iso_HS_1.Fill(tmp_result.fit_par_HS[0][2]);
            par_FB_Iso_HS_1.Fill(tmp_result.fit_par_HS[0][3]);
            
            par_Iso_Iso_HS_01.Fill(tmp_result.fit_par_HS[0][0]);
            par_NS_Iso_HS_01.Fill(tmp_result.fit_par_HS[0][1]);
            par_EW_Iso_HS_01.Fill(tmp_result.fit_par_HS[0][2]);
            par_FB_Iso_HS_01.Fill(tmp_result.fit_par_HS[0][3]);
            
            par_Iso_NS_HS_1.Fill(tmp_result.fit_par_HS[1][0]);
            par_NS_NS_HS_1.Fill(tmp_result.fit_par_HS[1][1]);
            par_EW_NS_HS_1.Fill(tmp_result.fit_par_HS[1][2]);
            par_FB_NS_HS_1.Fill(tmp_result.fit_par_HS[1][3]);
            
            par_Iso_EW_HS_1.Fill(tmp_result.fit_par_HS[2][0]);
            par_NS_EW_HS_1.Fill(tmp_result.fit_par_HS[2][1]);
            par_EW_EW_HS_1.Fill(tmp_result.fit_par_HS[2][2]);
            par_FB_EW_HS_1.Fill(tmp_result.fit_par_HS[2][3]);
            
            par_Iso_FB_HS_1.Fill(tmp_result.fit_par_HS[3][0]);
            par_NS_FB_HS_1.Fill(tmp_result.fit_par_HS[3][1]);
            par_EW_FB_HS_1.Fill(tmp_result.fit_par_HS[3][2]);
            par_FB_FB_HS_1.Fill(tmp_result.fit_par_HS[3][3]);
            
            par_Iso_NSEW_HS_1.Fill(tmp_result.fit_par_HS[4][0]);
            par_NS_NSEW_HS_1.Fill(tmp_result.fit_par_HS[4][1]);
            par_EW_NSEW_HS_1.Fill(tmp_result.fit_par_HS[4][2]);
            par_FB_NSEW_HS_1.Fill(tmp_result.fit_par_HS[4][3]);
            
            par_Iso_NSFB_HS_1.Fill(tmp_result.fit_par_HS[5][0]);
            par_NS_NSFB_HS_1.Fill(tmp_result.fit_par_HS[5][1]);
            par_EW_NSFB_HS_1.Fill(tmp_result.fit_par_HS[5][2]);
            par_FB_NSFB_HS_1.Fill(tmp_result.fit_par_HS[5][3]);
            
            par_Iso_EWFB_HS_1.Fill(tmp_result.fit_par_HS[6][0]);
            par_NS_EWFB_HS_1.Fill(tmp_result.fit_par_HS[6][1]);
            par_EW_EWFB_HS_1.Fill(tmp_result.fit_par_HS[6][2]);
            par_FB_EWFB_HS_1.Fill(tmp_result.fit_par_HS[6][3]);
            
            par_Iso_full_HS_1.Fill(tmp_result.fit_par_HS[7][0]);
            par_NS_full_HS_1.Fill(tmp_result.fit_par_HS[7][1]);
            par_EW_full_HS_1.Fill(tmp_result.fit_par_HS[7][2]);
            par_FB_full_HS_1.Fill(tmp_result.fit_par_HS[7][3]);
            
            ////////////// Filling error HS histos
            
            parerr_Iso_Iso_HS_10.Fill(tmp_result.fit_err_HS[0][0]);
            parerr_NS_Iso_HS_10.Fill(tmp_result.fit_err_HS[0][1]);
            parerr_EW_Iso_HS_10.Fill(tmp_result.fit_err_HS[0][2]);
            parerr_FB_Iso_HS_10.Fill(tmp_result.fit_err_HS[0][3]);
            
            parerr_Iso_Iso_HS_1.Fill(tmp_result.fit_err_HS[0][0]);
            parerr_NS_Iso_HS_1.Fill(tmp_result.fit_err_HS[0][1]);
            parerr_EW_Iso_HS_1.Fill(tmp_result.fit_err_HS[0][2]);
            parerr_FB_Iso_HS_1.Fill(tmp_result.fit_err_HS[0][3]);
            
            parerr_Iso_Iso_HS_01.Fill(tmp_result.fit_err_HS[0][0]);
            parerr_NS_Iso_HS_01.Fill(tmp_result.fit_err_HS[0][1]);
            parerr_EW_Iso_HS_01.Fill(tmp_result.fit_err_HS[0][2]);
            parerr_FB_Iso_HS_01.Fill(tmp_result.fit_err_HS[0][3]);
            
            parerr_Iso_NS_HS_1.Fill(tmp_result.fit_err_HS[1][0]);
            parerr_NS_NS_HS_1.Fill(tmp_result.fit_err_HS[1][1]);
            parerr_EW_NS_HS_1.Fill(tmp_result.fit_err_HS[1][2]);
            parerr_FB_NS_HS_1.Fill(tmp_result.fit_err_HS[1][3]);
            
            parerr_Iso_EW_HS_1.Fill(tmp_result.fit_err_HS[2][0]);
            parerr_NS_EW_HS_1.Fill(tmp_result.fit_err_HS[2][1]);
            parerr_EW_EW_HS_1.Fill(tmp_result.fit_err_HS[2][2]);
            parerr_FB_EW_HS_1.Fill(tmp_result.fit_err_HS[2][3]);
            
            parerr_Iso_FB_HS_1.Fill(tmp_result.fit_err_HS[3][0]);
            parerr_NS_FB_HS_1.Fill(tmp_result.fit_err_HS[3][1]);
            parerr_EW_FB_HS_1.Fill(tmp_result.fit_err_HS[3][2]);
            parerr_FB_FB_HS_1.Fill(tmp_result.fit_err_HS[3][3]);
            
            parerr_Iso_NSEW_HS_1.Fill(tmp_result.fit_err_HS[4][0]);
            parerr_NS_NSEW_HS_1.Fill(tmp_result.fit_err_HS[4][1]);
            parerr_EW_NSEW_HS_1.Fill(tmp_result.fit_err_HS[4][2]);
            parerr_FB_NSEW_HS_1.Fill(tmp_result.fit_err_HS[4][3]);
            
            parerr_Iso_NSFB_HS_1.Fill(tmp_result.fit_err_HS[5][0]);
            parerr_NS_NSFB_HS_1.Fill(tmp_result.fit_err_HS[5][1]);
            parerr_EW_NSFB_HS_1.Fill(tmp_result.fit_err_HS[5][2]);
            parerr_FB_NSFB_HS_1.Fill(tmp_result.fit_err_HS[5][3]);
            
            parerr_Iso_EWFB_HS_1.Fill(tmp_result.fit_err_HS[6][0]);
            parerr_NS_EWFB_HS_1.Fill(tmp_result.fit_err_HS[6][1]);
            parerr_EW_EWFB_HS_1.Fill(tmp_result.fit_err_HS[6][2]);
            parerr_FB_EWFB_HS_1.Fill(tmp_result.fit_err_HS[6][3]);
            
            parerr_Iso_full_HS_1.Fill(tmp_result.fit_err_HS[7][0]);
            parerr_NS_full_HS_1.Fill(tmp_result.fit_err_HS[7][1]);
            parerr_EW_full_HS_1.Fill(tmp_result.fit_err_HS[7][2]);
            parerr_FB_full_HS_1.Fill(tmp_result.fit_err_HS[7][3]);

            ////////////// Filling chi2 HS histos
            
            chi2_Iso_HS_1.Fill(tmp_result.chi2_HS[0]);
            chi2_NS_HS_1.Fill(tmp_result.chi2_HS[1]);
            chi2_EW_HS_1.Fill(tmp_result.chi2_HS[2]);
            chi2_FB_HS_1.Fill(tmp_result.chi2_HS[3]);
            chi2_NSEW_HS_1.Fill(tmp_result.chi2_HS[4]);
            chi2_NSFB_HS_1.Fill(tmp_result.chi2_HS[5]);
            chi2_EWFB_HS_1.Fill(tmp_result.chi2_HS[6]);
            chi2_full_HS_1.Fill(tmp_result.chi2_HS[7]);
            
            ////////////// Filling ndf HS histos
            
            ndf_Iso_HS_1.Fill(tmp_result.ndf_HS[0]);
            ndf_NS_HS_1.Fill(tmp_result.ndf_HS[1]);
            ndf_EW_HS_1.Fill(tmp_result.ndf_HS[2]);
            ndf_FB_HS_1.Fill(tmp_result.ndf_HS[3]);
            ndf_NSEW_HS_1.Fill(tmp_result.ndf_HS[4]);
            ndf_NSFB_HS_1.Fill(tmp_result.ndf_HS[5]);
            ndf_EWFB_HS_1.Fill(tmp_result.ndf_HS[6]);
            ndf_full_HS_1.Fill(tmp_result.ndf_HS[7]);

            
        }
        
        else
        {
            
            ////////////// Filling parameter LS histos
            
            par_Iso_Iso_LS_10.Fill(tmp_result.fit_par_LS[0][0]);
            par_NS_Iso_LS_10.Fill(tmp_result.fit_par_LS[0][1]);
            par_EW_Iso_LS_10.Fill(tmp_result.fit_par_LS[0][2]);
            par_FB_Iso_LS_10.Fill(tmp_result.fit_par_LS[0][3]);
            
            par_Iso_Iso_LS_1.Fill(tmp_result.fit_par_LS[0][0]);
            par_NS_Iso_LS_1.Fill(tmp_result.fit_par_LS[0][1]);
            par_EW_Iso_LS_1.Fill(tmp_result.fit_par_LS[0][2]);
            par_FB_Iso_LS_1.Fill(tmp_result.fit_par_LS[0][3]);
            
            par_Iso_Iso_LS_01.Fill(tmp_result.fit_par_LS[0][0]);
            par_NS_Iso_LS_01.Fill(tmp_result.fit_par_LS[0][1]);
            par_EW_Iso_LS_01.Fill(tmp_result.fit_par_LS[0][2]);
            par_FB_Iso_LS_01.Fill(tmp_result.fit_par_LS[0][3]);
            
            par_Iso_NS_LS_01.Fill(tmp_result.fit_par_LS[1][0]);
            par_NS_NS_LS_01.Fill(tmp_result.fit_par_LS[1][1]);
            par_EW_NS_LS_01.Fill(tmp_result.fit_par_LS[1][2]);
            par_FB_NS_LS_01.Fill(tmp_result.fit_par_LS[1][3]);
            
            par_Iso_EW_LS_01.Fill(tmp_result.fit_par_LS[2][0]);
            par_NS_EW_LS_01.Fill(tmp_result.fit_par_LS[2][1]);
            par_EW_EW_LS_01.Fill(tmp_result.fit_par_LS[2][2]);
            par_FB_EW_LS_01.Fill(tmp_result.fit_par_LS[2][3]);
            
            par_Iso_FB_LS_01.Fill(tmp_result.fit_par_LS[3][0]);
            par_NS_FB_LS_01.Fill(tmp_result.fit_par_LS[3][1]);
            par_EW_FB_LS_01.Fill(tmp_result.fit_par_LS[3][2]);
            par_FB_FB_LS_01.Fill(tmp_result.fit_par_LS[3][3]);
            
            par_Iso_NSEW_LS_01.Fill(tmp_result.fit_par_LS[4][0]);
            par_NS_NSEW_LS_01.Fill(tmp_result.fit_par_LS[4][1]);
            par_EW_NSEW_LS_01.Fill(tmp_result.fit_par_LS[4][2]);
            par_FB_NSEW_LS_01.Fill(tmp_result.fit_par_LS[4][3]);
            
            par_Iso_NSFB_LS_01.Fill(tmp_result.fit_par_LS[5][0]);
            par_NS_NSFB_LS_01.Fill(tmp_result.fit_par_LS[5][1]);
            par_EW_NSFB_LS_01.Fill(tmp_result.fit_par_LS[5][2]);
            par_FB_NSFB_LS_01.Fill(tmp_result.fit_par_LS[5][3]);
            
            par_Iso_EWFB_LS_01.Fill(tmp_result.fit_par_LS[6][0]);
            par_NS_EWFB_LS_01.Fill(tmp_result.fit_par_LS[6][1]);
            par_EW_EWFB_LS_01.Fill(tmp_result.fit_par_LS[6][2]);
            par_FB_EWFB_LS_01.Fill(tmp_result.fit_par_LS[6][3]);
            
            par_Iso_full_LS_01.Fill(tmp_result.fit_par_LS[7][0]);
            par_NS_full_LS_01.Fill(tmp_result.fit_par_LS[7][1]);
            par_EW_full_LS_01.Fill(tmp_result.fit_par_LS[7][2]);
            par_FB_full_LS_01.Fill(tmp_result.fit_par_LS[7][3]);
            
            ////////////// Filling error LS histos
            
            parerr_Iso_Iso_LS_10.Fill(tmp_result.fit_err_LS[0][0]);
            parerr_NS_Iso_LS_10.Fill(tmp_result.fit_err_LS[0][1]);
            parerr_EW_Iso_LS_10.Fill(tmp_result.fit_err_LS[0][2]);
            parerr_FB_Iso_LS_10.Fill(tmp_result.fit_err_LS[0][3]);
            
            parerr_Iso_Iso_LS_1.Fill(tmp_result.fit_err_LS[0][0]);
            parerr_NS_Iso_LS_1.Fill(tmp_result.fit_err_LS[0][1]);
            parerr_EW_Iso_LS_1.Fill(tmp_result.fit_err_LS[0][2]);
            parerr_FB_Iso_LS_1.Fill(tmp_result.fit_err_LS[0][3]);
            
            parerr_Iso_Iso_LS_01.Fill(tmp_result.fit_err_LS[0][0]);
            parerr_NS_Iso_LS_01.Fill(tmp_result.fit_err_LS[0][1]);
            parerr_EW_Iso_LS_01.Fill(tmp_result.fit_err_LS[0][2]);
            parerr_FB_Iso_LS_01.Fill(tmp_result.fit_err_LS[0][3]);
            
            parerr_Iso_NS_LS_01.Fill(tmp_result.fit_err_LS[1][0]);
            parerr_NS_NS_LS_01.Fill(tmp_result.fit_err_LS[1][1]);
            parerr_EW_NS_LS_01.Fill(tmp_result.fit_err_LS[1][2]);
            parerr_FB_NS_LS_01.Fill(tmp_result.fit_err_LS[1][3]);
            
            parerr_Iso_EW_LS_01.Fill(tmp_result.fit_err_LS[2][0]);
            parerr_NS_EW_LS_01.Fill(tmp_result.fit_err_LS[2][1]);
            parerr_EW_EW_LS_01.Fill(tmp_result.fit_err_LS[2][2]);
            parerr_FB_EW_LS_01.Fill(tmp_result.fit_err_LS[2][3]);
            
            parerr_Iso_FB_LS_01.Fill(tmp_result.fit_err_LS[3][0]);
            parerr_NS_FB_LS_01.Fill(tmp_result.fit_err_LS[3][1]);
            parerr_EW_FB_LS_01.Fill(tmp_result.fit_err_LS[3][2]);
            parerr_FB_FB_LS_01.Fill(tmp_result.fit_err_LS[3][3]);
            
            parerr_Iso_NSEW_LS_01.Fill(tmp_result.fit_err_LS[4][0]);
            parerr_NS_NSEW_LS_01.Fill(tmp_result.fit_err_LS[4][1]);
            parerr_EW_NSEW_LS_01.Fill(tmp_result.fit_err_LS[4][2]);
            parerr_FB_NSEW_LS_01.Fill(tmp_result.fit_err_LS[4][3]);
            
            parerr_Iso_NSFB_LS_01.Fill(tmp_result.fit_err_LS[5][0]);
            parerr_NS_NSFB_LS_01.Fill(tmp_result.fit_err_LS[5][1]);
            parerr_EW_NSFB_LS_01.Fill(tmp_result.fit_err_LS[5][2]);
            parerr_FB_NSFB_LS_01.Fill(tmp_result.fit_err_LS[5][3]);
            
            parerr_Iso_EWFB_LS_01.Fill(tmp_result.fit_err_LS[6][0]);
            parerr_NS_EWFB_LS_01.Fill(tmp_result.fit_err_LS[6][1]);
            parerr_EW_EWFB_LS_01.Fill(tmp_result.fit_err_LS[6][2]);
            parerr_FB_EWFB_LS_01.Fill(tmp_result.fit_err_LS[6][3]);
            
            parerr_Iso_full_LS_01.Fill(tmp_result.fit_err_LS[7][0]);
            parerr_NS_full_LS_01.Fill(tmp_result.fit_err_LS[7][1]);
            parerr_EW_full_LS_01.Fill(tmp_result.fit_err_LS[7][2]);
            parerr_FB_full_LS_01.Fill(tmp_result.fit_err_LS[7][3]);

            ////////////// Filling chi2 LS histos
            
            chi2_Iso_LS_01.Fill(tmp_result.chi2_LS[0]);
            chi2_NS_LS_01.Fill(tmp_result.chi2_LS[1]);
            chi2_EW_LS_01.Fill(tmp_result.chi2_LS[2]);
            chi2_FB_LS_01.Fill(tmp_result.chi2_LS[3]);
            chi2_NSEW_LS_01.Fill(tmp_result.chi2_LS[4]);
            chi2_NSFB_LS_01.Fill(tmp_result.chi2_LS[5]);
            chi2_EWFB_LS_01.Fill(tmp_result.chi2_LS[6]);
            chi2_full_LS_01.Fill(tmp_result.chi2_LS[7]);
            
            ////////////// Filling ndf LS histos
            
            ndf_Iso_LS_01.Fill(tmp_result.ndf_LS[0]);
            ndf_NS_LS_01.Fill(tmp_result.ndf_LS[1]);
            ndf_EW_LS_01.Fill(tmp_result.ndf_LS[2]);
            ndf_FB_LS_01.Fill(tmp_result.ndf_LS[3]);
            ndf_NSEW_LS_01.Fill(tmp_result.ndf_LS[4]);
            ndf_NSFB_LS_01.Fill(tmp_result.ndf_LS[5]);
            ndf_EWFB_LS_01.Fill(tmp_result.ndf_LS[6]);
            ndf_full_LS_01.Fill(tmp_result.ndf_LS[7]);
            
            
            
            
            
            
            ////////////// Filling parameter HS histos
            
            par_Iso_Iso_HS_10.Fill(tmp_result.fit_par_HS[0][0]);
            par_NS_Iso_HS_10.Fill(tmp_result.fit_par_HS[0][1]);
            par_EW_Iso_HS_10.Fill(tmp_result.fit_par_HS[0][2]);
            par_FB_Iso_HS_10.Fill(tmp_result.fit_par_HS[0][3]);
            
            par_Iso_Iso_HS_1.Fill(tmp_result.fit_par_HS[0][0]);
            par_NS_Iso_HS_1.Fill(tmp_result.fit_par_HS[0][1]);
            par_EW_Iso_HS_1.Fill(tmp_result.fit_par_HS[0][2]);
            par_FB_Iso_HS_1.Fill(tmp_result.fit_par_HS[0][3]);
            
            par_Iso_Iso_HS_01.Fill(tmp_result.fit_par_HS[0][0]);
            par_NS_Iso_HS_01.Fill(tmp_result.fit_par_HS[0][1]);
            par_EW_Iso_HS_01.Fill(tmp_result.fit_par_HS[0][2]);
            par_FB_Iso_HS_01.Fill(tmp_result.fit_par_HS[0][3]);
            
            par_Iso_NS_HS_01.Fill(tmp_result.fit_par_HS[1][0]);
            par_NS_NS_HS_01.Fill(tmp_result.fit_par_HS[1][1]);
            par_EW_NS_HS_01.Fill(tmp_result.fit_par_HS[1][2]);
            par_FB_NS_HS_01.Fill(tmp_result.fit_par_HS[1][3]);
            
            par_Iso_EW_HS_01.Fill(tmp_result.fit_par_HS[2][0]);
            par_NS_EW_HS_01.Fill(tmp_result.fit_par_HS[2][1]);
            par_EW_EW_HS_01.Fill(tmp_result.fit_par_HS[2][2]);
            par_FB_EW_HS_01.Fill(tmp_result.fit_par_HS[2][3]);
            
            par_Iso_FB_HS_01.Fill(tmp_result.fit_par_HS[3][0]);
            par_NS_FB_HS_01.Fill(tmp_result.fit_par_HS[3][1]);
            par_EW_FB_HS_01.Fill(tmp_result.fit_par_HS[3][2]);
            par_FB_FB_HS_01.Fill(tmp_result.fit_par_HS[3][3]);
            
            par_Iso_NSEW_HS_01.Fill(tmp_result.fit_par_HS[4][0]);
            par_NS_NSEW_HS_01.Fill(tmp_result.fit_par_HS[4][1]);
            par_EW_NSEW_HS_01.Fill(tmp_result.fit_par_HS[4][2]);
            par_FB_NSEW_HS_01.Fill(tmp_result.fit_par_HS[4][3]);
            
            par_Iso_NSFB_HS_01.Fill(tmp_result.fit_par_HS[5][0]);
            par_NS_NSFB_HS_01.Fill(tmp_result.fit_par_HS[5][1]);
            par_EW_NSFB_HS_01.Fill(tmp_result.fit_par_HS[5][2]);
            par_FB_NSFB_HS_01.Fill(tmp_result.fit_par_HS[5][3]);
            
            par_Iso_EWFB_HS_01.Fill(tmp_result.fit_par_HS[6][0]);
            par_NS_EWFB_HS_01.Fill(tmp_result.fit_par_HS[6][1]);
            par_EW_EWFB_HS_01.Fill(tmp_result.fit_par_HS[6][2]);
            par_FB_EWFB_HS_01.Fill(tmp_result.fit_par_HS[6][3]);
            
            par_Iso_full_HS_01.Fill(tmp_result.fit_par_HS[7][0]);
            par_NS_full_HS_01.Fill(tmp_result.fit_par_HS[7][1]);
            par_EW_full_HS_01.Fill(tmp_result.fit_par_HS[7][2]);
            par_FB_full_HS_01.Fill(tmp_result.fit_par_HS[7][3]);
            
            ////////////// Filling error HS histos
            
            parerr_Iso_Iso_HS_10.Fill(tmp_result.fit_err_HS[0][0]);
            parerr_NS_Iso_HS_10.Fill(tmp_result.fit_err_HS[0][1]);
            parerr_EW_Iso_HS_10.Fill(tmp_result.fit_err_HS[0][2]);
            parerr_FB_Iso_HS_10.Fill(tmp_result.fit_err_HS[0][3]);
            
            parerr_Iso_Iso_HS_1.Fill(tmp_result.fit_err_HS[0][0]);
            parerr_NS_Iso_HS_1.Fill(tmp_result.fit_err_HS[0][1]);
            parerr_EW_Iso_HS_1.Fill(tmp_result.fit_err_HS[0][2]);
            parerr_FB_Iso_HS_1.Fill(tmp_result.fit_err_HS[0][3]);
            
            parerr_Iso_Iso_HS_01.Fill(tmp_result.fit_err_HS[0][0]);
            parerr_NS_Iso_HS_01.Fill(tmp_result.fit_err_HS[0][1]);
            parerr_EW_Iso_HS_01.Fill(tmp_result.fit_err_HS[0][2]);
            parerr_FB_Iso_HS_01.Fill(tmp_result.fit_err_HS[0][3]);
            
            parerr_Iso_NS_HS_01.Fill(tmp_result.fit_err_HS[1][0]);
            parerr_NS_NS_HS_01.Fill(tmp_result.fit_err_HS[1][1]);
            parerr_EW_NS_HS_01.Fill(tmp_result.fit_err_HS[1][2]);
            parerr_FB_NS_HS_01.Fill(tmp_result.fit_err_HS[1][3]);
            
            parerr_Iso_EW_HS_01.Fill(tmp_result.fit_err_HS[2][0]);
            parerr_NS_EW_HS_01.Fill(tmp_result.fit_err_HS[2][1]);
            parerr_EW_EW_HS_01.Fill(tmp_result.fit_err_HS[2][2]);
            parerr_FB_EW_HS_01.Fill(tmp_result.fit_err_HS[2][3]);
            
            parerr_Iso_FB_HS_01.Fill(tmp_result.fit_err_HS[3][0]);
            parerr_NS_FB_HS_01.Fill(tmp_result.fit_err_HS[3][1]);
            parerr_EW_FB_HS_01.Fill(tmp_result.fit_err_HS[3][2]);
            parerr_FB_FB_HS_01.Fill(tmp_result.fit_err_HS[3][3]);
            
            parerr_Iso_NSEW_HS_01.Fill(tmp_result.fit_err_HS[4][0]);
            parerr_NS_NSEW_HS_01.Fill(tmp_result.fit_err_HS[4][1]);
            parerr_EW_NSEW_HS_01.Fill(tmp_result.fit_err_HS[4][2]);
            parerr_FB_NSEW_HS_01.Fill(tmp_result.fit_err_HS[4][3]);
            
            parerr_Iso_NSFB_HS_01.Fill(tmp_result.fit_err_HS[5][0]);
            parerr_NS_NSFB_HS_01.Fill(tmp_result.fit_err_HS[5][1]);
            parerr_EW_NSFB_HS_01.Fill(tmp_result.fit_err_HS[5][2]);
            parerr_FB_NSFB_HS_01.Fill(tmp_result.fit_err_HS[5][3]);
            
            parerr_Iso_EWFB_HS_01.Fill(tmp_result.fit_err_HS[6][0]);
            parerr_NS_EWFB_HS_01.Fill(tmp_result.fit_err_HS[6][1]);
            parerr_EW_EWFB_HS_01.Fill(tmp_result.fit_err_HS[6][2]);
            parerr_FB_EWFB_HS_01.Fill(tmp_result.fit_err_HS[6][3]);
            
            parerr_Iso_full_HS_01.Fill(tmp_result.fit_err_HS[7][0]);
            parerr_NS_full_HS_01.Fill(tmp_result.fit_err_HS[7][1]);
            parerr_EW_full_HS_01.Fill(tmp_result.fit_err_HS[7][2]);
            parerr_FB_full_HS_01.Fill(tmp_result.fit_err_HS[7][3]);

            ////////////// Filling chi2 HS histos
            
            chi2_Iso_HS_01.Fill(tmp_result.chi2_HS[0]);
            chi2_NS_HS_01.Fill(tmp_result.chi2_HS[1]);
            chi2_EW_HS_01.Fill(tmp_result.chi2_HS[2]);
            chi2_FB_HS_01.Fill(tmp_result.chi2_HS[3]);
            chi2_NSEW_HS_01.Fill(tmp_result.chi2_HS[4]);
            chi2_NSFB_HS_01.Fill(tmp_result.chi2_HS[5]);
            chi2_EWFB_HS_01.Fill(tmp_result.chi2_HS[6]);
            chi2_full_HS_01.Fill(tmp_result.chi2_HS[7]);
            
            ////////////// Filling ndf HS histos
            
            ndf_Iso_HS_01.Fill(tmp_result.ndf_HS[0]);
            ndf_NS_HS_01.Fill(tmp_result.ndf_HS[1]);
            ndf_EW_HS_01.Fill(tmp_result.ndf_HS[2]);
            ndf_FB_HS_01.Fill(tmp_result.ndf_HS[3]);
            ndf_NSEW_HS_01.Fill(tmp_result.ndf_HS[4]);
            ndf_NSFB_HS_01.Fill(tmp_result.ndf_HS[5]);
            ndf_EWFB_HS_01.Fill(tmp_result.ndf_HS[6]);
            ndf_full_HS_01.Fill(tmp_result.ndf_HS[7]);
            
        }
        
    }
 
    
    
    /////////////////////////////// Writing results to file
    
    TFile outFile_LS_10(results_path_10_LS.c_str(),"RECREATE");
    if(outFile_LS_10.IsZombie())
    {
        std::cout << "\n\nError writing ROOT output file: " << results_path_10_LS << "\n\n";
        exit(100);
    }
    
    par_Iso_Iso_LS_10.Write();
    par_NS_Iso_LS_10.Write();
    par_EW_Iso_LS_10.Write();
    par_FB_Iso_LS_10.Write();
    
    par_Iso_NS_LS_10.Write();
    par_NS_NS_LS_10.Write();
    par_EW_NS_LS_10.Write();
    par_FB_NS_LS_10.Write();
    
    par_Iso_EW_LS_10.Write();
    par_NS_EW_LS_10.Write();
    par_EW_EW_LS_10.Write();
    par_FB_EW_LS_10.Write();
    
    par_Iso_FB_LS_10.Write();
    par_NS_FB_LS_10.Write();
    par_EW_FB_LS_10.Write();
    par_FB_FB_LS_10.Write();
    
    par_Iso_NSEW_LS_10.Write();
    par_NS_NSEW_LS_10.Write();
    par_EW_NSEW_LS_10.Write();
    par_FB_NSEW_LS_10.Write();
    
    par_Iso_NSFB_LS_10.Write();
    par_NS_NSFB_LS_10.Write();
    par_EW_NSFB_LS_10.Write();
    par_FB_NSFB_LS_10.Write();
    
    par_Iso_EWFB_LS_10.Write();
    par_NS_EWFB_LS_10.Write();
    par_EW_EWFB_LS_10.Write();
    par_FB_EWFB_LS_10.Write();
    
    par_Iso_full_LS_10.Write();
    par_NS_full_LS_10.Write();
    par_EW_full_LS_10.Write();
    par_FB_full_LS_10.Write();
    
    parerr_Iso_Iso_LS_10.Write();
    parerr_NS_Iso_LS_10.Write();
    parerr_EW_Iso_LS_10.Write();
    parerr_FB_Iso_LS_10.Write();
    
    parerr_Iso_NS_LS_10.Write();
    parerr_NS_NS_LS_10.Write();
    parerr_EW_NS_LS_10.Write();
    parerr_FB_NS_LS_10.Write();
    
    parerr_Iso_EW_LS_10.Write();
    parerr_NS_EW_LS_10.Write();
    parerr_EW_EW_LS_10.Write();
    parerr_FB_EW_LS_10.Write();
    
    parerr_Iso_FB_LS_10.Write();
    parerr_NS_FB_LS_10.Write();
    parerr_EW_FB_LS_10.Write();
    parerr_FB_FB_LS_10.Write();
    
    parerr_Iso_NSEW_LS_10.Write();
    parerr_NS_NSEW_LS_10.Write();
    parerr_EW_NSEW_LS_10.Write();
    parerr_FB_NSEW_LS_10.Write();
    
    parerr_Iso_NSFB_LS_10.Write();
    parerr_NS_NSFB_LS_10.Write();
    parerr_EW_NSFB_LS_10.Write();
    parerr_FB_NSFB_LS_10.Write();
    
    parerr_Iso_EWFB_LS_10.Fill(tmp_result.fit_err_LS[6][0]);
    parerr_NS_EWFB_LS_10.Fill(tmp_result.fit_err_LS[6][1]);
    parerr_EW_EWFB_LS_10.Fill(tmp_result.fit_err_LS[6][2]);
    parerr_FB_EWFB_LS_10.Fill(tmp_result.fit_err_LS[6][3]);
    
    parerr_Iso_full_LS_10.Fill(tmp_result.fit_err_LS[7][0]);
    parerr_NS_full_LS_10.Fill(tmp_result.fit_err_LS[7][1]);
    parerr_EW_full_LS_10.Fill(tmp_result.fit_err_LS[7][2]);
    parerr_FB_full_LS_10.Fill(tmp_result.fit_err_LS[7][3]);
    
    chi2_Iso_LS_10.Write();
    chi2_NS_LS_10.Write();
    chi2_EW_LS_10.Write();
    chi2_FB_LS_10.Write();
    chi2_NSEW_LS_10.Write();
    chi2_NSFB_LS_10.Write();
    chi2_EWFB_LS_10.Write();
    chi2_full_LS_10.Write();
    
    ndf_Iso_LS_10.Write();
    ndf_NS_LS_10.Write();
    ndf_EW_LS_10.Write();
    ndf_FB_LS_10.Write();
    ndf_NSEW_LS_10.Write();
    ndf_NSFB_LS_10.Write();
    ndf_EWFB_LS_10.Write();
    ndf_full_LS_10.Write();
    
    outFile_LS_10.Write();
    outFile_LS_10.Close();
    
    //////////////////////////////
    
    TFile outFile_LS_1(results_path_1_LS.c_str(),"RECREATE");
    if(outFile_LS_1.IsZombie())
    {
        std::cout << "\n\nError writing ROOT output file: " << results_path_1_LS << "\n\n";
        exit(100);
    }
    
    par_Iso_Iso_LS_1.Write();
    par_NS_Iso_LS_1.Write();
    par_EW_Iso_LS_1.Write();
    par_FB_Iso_LS_1.Write();
    
    par_Iso_NS_LS_1.Write();
    par_NS_NS_LS_1.Write();
    par_EW_NS_LS_1.Write();
    par_FB_NS_LS_1.Write();
    
    par_Iso_EW_LS_1.Write();
    par_NS_EW_LS_1.Write();
    par_EW_EW_LS_1.Write();
    par_FB_EW_LS_1.Write();
    
    par_Iso_FB_LS_1.Write();
    par_NS_FB_LS_1.Write();
    par_EW_FB_LS_1.Write();
    par_FB_FB_LS_1.Write();
    
    par_Iso_NSEW_LS_1.Write();
    par_NS_NSEW_LS_1.Write();
    par_EW_NSEW_LS_1.Write();
    par_FB_NSEW_LS_1.Write();
    
    par_Iso_NSFB_LS_1.Write();
    par_NS_NSFB_LS_1.Write();
    par_EW_NSFB_LS_1.Write();
    par_FB_NSFB_LS_1.Write();
    
    par_Iso_EWFB_LS_1.Write();
    par_NS_EWFB_LS_1.Write();
    par_EW_EWFB_LS_1.Write();
    par_FB_EWFB_LS_1.Write();
    
    par_Iso_full_LS_1.Write();
    par_NS_full_LS_1.Write();
    par_EW_full_LS_1.Write();
    par_FB_full_LS_1.Write();
    
    parerr_Iso_Iso_LS_1.Write();
    parerr_NS_Iso_LS_1.Write();
    parerr_EW_Iso_LS_1.Write();
    parerr_FB_Iso_LS_1.Write();
    
    parerr_Iso_NS_LS_1.Write();
    parerr_NS_NS_LS_1.Write();
    parerr_EW_NS_LS_1.Write();
    parerr_FB_NS_LS_1.Write();
    
    parerr_Iso_EW_LS_1.Write();
    parerr_NS_EW_LS_1.Write();
    parerr_EW_EW_LS_1.Write();
    parerr_FB_EW_LS_1.Write();
    
    parerr_Iso_FB_LS_1.Write();
    parerr_NS_FB_LS_1.Write();
    parerr_EW_FB_LS_1.Write();
    parerr_FB_FB_LS_1.Write();
    
    parerr_Iso_NSEW_LS_1.Write();
    parerr_NS_NSEW_LS_1.Write();
    parerr_EW_NSEW_LS_1.Write();
    parerr_FB_NSEW_LS_1.Write();
    
    parerr_Iso_NSFB_LS_1.Write();
    parerr_NS_NSFB_LS_1.Write();
    parerr_EW_NSFB_LS_1.Write();
    parerr_FB_NSFB_LS_1.Write();
    
    parerr_Iso_EWFB_LS_1.Write();
    parerr_NS_EWFB_LS_1.Write();
    parerr_EW_EWFB_LS_1.Write();
    parerr_FB_EWFB_LS_1.Write();
    
    parerr_Iso_full_LS_1.Write();
    parerr_NS_full_LS_1.Write();
    parerr_EW_full_LS_1.Write();
    parerr_FB_full_LS_1.Write();
    
    chi2_Iso_LS_1.Write();
    chi2_NS_LS_1.Write();
    chi2_EW_LS_1.Write();
    chi2_FB_LS_1.Write();
    chi2_NSEW_LS_1.Write();
    chi2_NSFB_LS_1.Write();
    chi2_EWFB_LS_1.Write();
    chi2_full_LS_1.Write();
    
    ndf_Iso_LS_1.Write();
    ndf_NS_LS_1.Write();
    ndf_EW_LS_1.Write();
    ndf_FB_LS_1.Write();
    ndf_NSEW_LS_1.Write();
    ndf_NSFB_LS_1.Write();
    ndf_EWFB_LS_1.Write();
    ndf_full_LS_1.Write();
    
    outFile_LS_1.Write();
    outFile_LS_1.Close();
    
    //////////////////////////////
    
    TFile outFile_LS_01(results_path_01_LS.c_str(),"RECREATE");
    if(outFile_LS_01.IsZombie())
    {
        std::cout << "\n\nError writing ROOT output file: " << results_path_01_LS << "\n\n";
        exit(100);
    }
    
    par_Iso_Iso_LS_01.Write();
    par_NS_Iso_LS_01.Write();
    par_EW_Iso_LS_01.Write();
    par_FB_Iso_LS_01.Write();
    
    par_Iso_NS_LS_01.Write();
    par_NS_NS_LS_01.Write();
    par_EW_NS_LS_01.Write();
    par_FB_NS_LS_01.Write();
    
    par_Iso_EW_LS_01.Write();
    par_NS_EW_LS_01.Write();
    par_EW_EW_LS_01.Write();
    par_FB_EW_LS_01.Write();
    
    par_Iso_FB_LS_01.Write();
    par_NS_FB_LS_01.Write();
    par_EW_FB_LS_01.Write();
    par_FB_FB_LS_01.Write();
    
    par_Iso_NSEW_LS_01.Write();
    par_NS_NSEW_LS_01.Write();
    par_EW_NSEW_LS_01.Write();
    par_FB_NSEW_LS_01.Write();
    
    par_Iso_NSFB_LS_01.Write();
    par_NS_NSFB_LS_01.Write();
    par_EW_NSFB_LS_01.Write();
    par_FB_NSFB_LS_01.Write();
    
    par_Iso_EWFB_LS_01.Write();
    par_NS_EWFB_LS_01.Write();
    par_EW_EWFB_LS_01.Write();
    par_FB_EWFB_LS_01.Write();
    
    par_Iso_full_LS_01.Write();
    par_NS_full_LS_01.Write();
    par_EW_full_LS_01.Write();
    par_FB_full_LS_01.Write();

    parerr_Iso_Iso_LS_01.Write();
    parerr_NS_Iso_LS_01.Write();
    parerr_EW_Iso_LS_01.Write();
    parerr_FB_Iso_LS_01.Write();
    
    parerr_Iso_NS_LS_01.Write();
    parerr_NS_NS_LS_01.Write();
    parerr_EW_NS_LS_01.Write();
    parerr_FB_NS_LS_01.Write();
    
    parerr_Iso_EW_LS_01.Write();
    parerr_NS_EW_LS_01.Write();
    parerr_EW_EW_LS_01.Write();
    parerr_FB_EW_LS_01.Write();
    
    parerr_Iso_FB_LS_01.Write();
    parerr_NS_FB_LS_01.Write();
    parerr_EW_FB_LS_01.Write();
    parerr_FB_FB_LS_01.Write();
    
    parerr_Iso_NSEW_LS_01.Write();
    parerr_NS_NSEW_LS_01.Write();
    parerr_EW_NSEW_LS_01.Write();
    parerr_FB_NSEW_LS_01.Write();
    
    parerr_Iso_NSFB_LS_01.Write();
    parerr_NS_NSFB_LS_01.Write();
    parerr_EW_NSFB_LS_01.Write();
    parerr_FB_NSFB_LS_01.Write();
    
    parerr_Iso_EWFB_LS_01.Write();
    parerr_NS_EWFB_LS_01.Write();
    parerr_EW_EWFB_LS_01.Write();
    parerr_FB_EWFB_LS_01.Write();
    
    parerr_Iso_full_LS_01.Write();
    parerr_NS_full_LS_01.Write();
    parerr_EW_full_LS_01.Write();
    parerr_FB_full_LS_01.Write();
    
    chi2_Iso_LS_01.Write();
    chi2_NS_LS_01.Write();
    chi2_EW_LS_01.Write();
    chi2_FB_LS_01.Write();
    chi2_NSEW_LS_01.Write();
    chi2_NSFB_LS_01.Write();
    chi2_EWFB_LS_01.Write();
    chi2_full_LS_01.Write();
    
    ndf_Iso_LS_01.Write();
    ndf_NS_LS_01.Write();
    ndf_EW_LS_01.Write();
    ndf_FB_LS_01.Write();
    ndf_NSEW_LS_01.Write();
    ndf_NSFB_LS_01.Write();
    ndf_EWFB_LS_01.Write();
    ndf_full_LS_01.Write();
    
    outFile_LS_01.Write();
    outFile_LS_01.Close();
    

    //////////////////////////////
    
    TFile outFile_HS_10(results_path_10_HS.c_str(),"RECREATE");
    if(outFile_HS_10.IsZombie())
    {
        std::cout << "\n\nError writing ROOT output file: " << results_path_10_HS << "\n\n";
        exit(100);
    }
    
    par_Iso_Iso_HS_10.Write();
    par_NS_Iso_HS_10.Write();
    par_EW_Iso_HS_10.Write();
    par_FB_Iso_HS_10.Write();
    
    par_Iso_NS_HS_10.Write();
    par_NS_NS_HS_10.Write();
    par_EW_NS_HS_10.Write();
    par_FB_NS_HS_10.Write();
    
    par_Iso_EW_HS_10.Write();
    par_NS_EW_HS_10.Write();
    par_EW_EW_HS_10.Write();
    par_FB_EW_HS_10.Write();
    
    par_Iso_FB_HS_10.Write();
    par_NS_FB_HS_10.Write();
    par_EW_FB_HS_10.Write();
    par_FB_FB_HS_10.Write();
    
    par_Iso_NSEW_HS_10.Write();
    par_NS_NSEW_HS_10.Write();
    par_EW_NSEW_HS_10.Write();
    par_FB_NSEW_HS_10.Write();
    
    par_Iso_NSFB_HS_10.Write();
    par_NS_NSFB_HS_10.Write();
    par_EW_NSFB_HS_10.Write();
    par_FB_NSFB_HS_10.Write();
    
    par_Iso_EWFB_HS_10.Write();
    par_NS_EWFB_HS_10.Write();
    par_EW_EWFB_HS_10.Write();
    par_FB_EWFB_HS_10.Write();
    
    par_Iso_full_HS_10.Write();
    par_NS_full_HS_10.Write();
    par_EW_full_HS_10.Write();
    par_FB_full_HS_10.Write();
    
    parerr_Iso_Iso_HS_10.Write();
    parerr_NS_Iso_HS_10.Write();
    parerr_EW_Iso_HS_10.Write();
    parerr_FB_Iso_HS_10.Write();
    
    parerr_Iso_NS_HS_10.Write();
    parerr_NS_NS_HS_10.Write();
    parerr_EW_NS_HS_10.Write();
    parerr_FB_NS_HS_10.Write();
    
    parerr_Iso_EW_HS_10.Write();
    parerr_NS_EW_HS_10.Write();
    parerr_EW_EW_HS_10.Write();
    parerr_FB_EW_HS_10.Write();
    
    parerr_Iso_FB_HS_10.Write();
    parerr_NS_FB_HS_10.Write();
    parerr_EW_FB_HS_10.Write();
    parerr_FB_FB_HS_10.Write();
    
    parerr_Iso_NSEW_HS_10.Write();
    parerr_NS_NSEW_HS_10.Write();
    parerr_EW_NSEW_HS_10.Write();
    parerr_FB_NSEW_HS_10.Write();
    
    parerr_Iso_NSFB_HS_10.Write();
    parerr_NS_NSFB_HS_10.Write();
    parerr_EW_NSFB_HS_10.Write();
    parerr_FB_NSFB_HS_10.Write();
    
    parerr_Iso_EWFB_HS_10.Write();
    parerr_NS_EWFB_HS_10.Write();
    parerr_EW_EWFB_HS_10.Write();
    parerr_FB_EWFB_HS_10.Write();
    
    parerr_Iso_full_HS_10.Write();
    parerr_NS_full_HS_10.Write();
    parerr_EW_full_HS_10.Write();
    parerr_FB_full_HS_10.Write();
    
    chi2_Iso_HS_10.Write();
    chi2_NS_HS_10.Write();
    chi2_EW_HS_10.Write();
    chi2_FB_HS_10.Write();
    chi2_NSEW_HS_10.Write();
    chi2_NSFB_HS_10.Write();
    chi2_EWFB_HS_10.Write();
    chi2_full_HS_10.Write();
    
    ndf_Iso_HS_10.Write();
    ndf_NS_HS_10.Write();
    ndf_EW_HS_10.Write();
    ndf_FB_HS_10.Write();
    ndf_NSEW_HS_10.Write();
    ndf_NSFB_HS_10.Write();
    ndf_EWFB_HS_10.Write();
    ndf_full_HS_10.Write();
    
    outFile_HS_10.Write();
    outFile_HS_10.Close();
    
    //////////////////////////////
    
    TFile outFile_HS_1(results_path_1_HS.c_str(),"RECREATE");
    if(outFile_HS_1.IsZombie())
    {
        std::cout << "\n\nError writing ROOT output file: " << results_path_1_HS << "\n\n";
        exit(100);
    }
    
    par_Iso_Iso_HS_1.Write();
    par_NS_Iso_HS_1.Write();
    par_EW_Iso_HS_1.Write();
    par_FB_Iso_HS_1.Write();
    
    par_Iso_NS_HS_1.Write();
    par_NS_NS_HS_1.Write();
    par_EW_NS_HS_1.Write();
    par_FB_NS_HS_1.Write();
    
    par_Iso_EW_HS_1.Write();
    par_NS_EW_HS_1.Write();
    par_EW_EW_HS_1.Write();
    par_FB_EW_HS_1.Write();
    
    par_Iso_FB_HS_1.Write();
    par_NS_FB_HS_1.Write();
    par_EW_FB_HS_1.Write();
    par_FB_FB_HS_1.Write();
    
    par_Iso_NSEW_HS_1.Write();
    par_NS_NSEW_HS_1.Write();
    par_EW_NSEW_HS_1.Write();
    par_FB_NSEW_HS_1.Write();
    
    par_Iso_NSFB_HS_1.Write();
    par_NS_NSFB_HS_1.Write();
    par_EW_NSFB_HS_1.Write();
    par_FB_NSFB_HS_1.Write();
    
    par_Iso_EWFB_HS_1.Write();
    par_NS_EWFB_HS_1.Write();
    par_EW_EWFB_HS_1.Write();
    par_FB_EWFB_HS_1.Write();
    
    par_Iso_full_HS_1.Write();
    par_NS_full_HS_1.Write();
    par_EW_full_HS_1.Write();
    par_FB_full_HS_1.Write();
    
    parerr_Iso_Iso_HS_1.Write();
    parerr_NS_Iso_HS_1.Write();
    parerr_EW_Iso_HS_1.Write();
    parerr_FB_Iso_HS_1.Write();

    parerr_Iso_NS_HS_1.Write();
    parerr_NS_NS_HS_1.Write();
    parerr_EW_NS_HS_1.Write();
    parerr_FB_NS_HS_1.Write();
    
    parerr_Iso_EW_HS_1.Write();
    parerr_NS_EW_HS_1.Write();
    parerr_EW_EW_HS_1.Write();
    parerr_FB_EW_HS_1.Write();
    
    parerr_Iso_FB_HS_1.Write();
    parerr_NS_FB_HS_1.Write();
    parerr_EW_FB_HS_1.Write();
    parerr_FB_FB_HS_1.Write();
    
    parerr_Iso_NSEW_HS_1.Write();
    parerr_NS_NSEW_HS_1.Write();
    parerr_EW_NSEW_HS_1.Write();
    parerr_FB_NSEW_HS_1.Write();
    
    parerr_Iso_NSFB_HS_1.Write();
    parerr_NS_NSFB_HS_1.Write();
    parerr_EW_NSFB_HS_1.Write();
    parerr_FB_NSFB_HS_1.Write();
    
    parerr_Iso_EWFB_HS_1.Write();
    parerr_NS_EWFB_HS_1.Write();
    parerr_EW_EWFB_HS_1.Write();
    parerr_FB_EWFB_HS_1.Write();
    
    parerr_Iso_full_HS_1.Write();
    parerr_NS_full_HS_1.Write();
    parerr_EW_full_HS_1.Write();
    parerr_FB_full_HS_1.Write();
    
    chi2_Iso_HS_1.Write();
    chi2_NS_HS_1.Write();
    chi2_EW_HS_1.Write();
    chi2_FB_HS_1.Write();
    chi2_NSEW_HS_1.Write();
    chi2_NSFB_HS_1.Write();
    chi2_EWFB_HS_1.Write();
    chi2_full_HS_1.Write();
    
    ndf_Iso_HS_1.Write();
    ndf_NS_HS_1.Write();
    ndf_EW_HS_1.Write();
    ndf_FB_HS_1.Write();
    ndf_NSEW_HS_1.Write();
    ndf_NSFB_HS_1.Write();
    ndf_EWFB_HS_1.Write();
    ndf_full_HS_1.Write();
    
    outFile_HS_1.Write();
    outFile_HS_1.Close();
    
    //////////////////////////////
    
    TFile outFile_HS_01(results_path_01_HS.c_str(),"RECREATE");
    if(outFile_HS_01.IsZombie())
    {
        std::cout << "\n\nError writing ROOT output file: " << results_path_01_HS << "\n\n";
        exit(100);
    }
    
    par_Iso_Iso_HS_01.Write();
    par_NS_Iso_HS_01.Write();
    par_EW_Iso_HS_01.Write();
    par_FB_Iso_HS_01.Write();
    
    par_Iso_NS_HS_01.Write();
    par_NS_NS_HS_01.Write();
    par_EW_NS_HS_01.Write();
    par_FB_NS_HS_01.Write();
    
    par_Iso_EW_HS_01.Write();
    par_NS_EW_HS_01.Write();
    par_EW_EW_HS_01.Write();
    par_FB_EW_HS_01.Write();
    
    par_Iso_FB_HS_01.Write();
    par_NS_FB_HS_01.Write();
    par_EW_FB_HS_01.Write();
    par_FB_FB_HS_01.Write();
    
    par_Iso_NSEW_HS_01.Write();
    par_NS_NSEW_HS_01.Write();
    par_EW_NSEW_HS_01.Write();
    par_FB_NSEW_HS_01.Write();
    
    par_Iso_NSFB_HS_01.Write();
    par_NS_NSFB_HS_01.Write();
    par_EW_NSFB_HS_01.Write();
    par_FB_NSFB_HS_01.Write();
    
    par_Iso_EWFB_HS_01.Write();
    par_NS_EWFB_HS_01.Write();
    par_EW_EWFB_HS_01.Write();
    par_FB_EWFB_HS_01.Write();
    
    par_Iso_full_HS_01.Write();
    par_NS_full_HS_01.Write();
    par_EW_full_HS_01.Write();
    par_FB_full_HS_01.Write();
    
    parerr_Iso_Iso_HS_01.Write();
    parerr_NS_Iso_HS_01.Write();
    parerr_EW_Iso_HS_01.Write();
    parerr_FB_Iso_HS_01.Write();
    
    parerr_Iso_NS_HS_01.Write();
    parerr_NS_NS_HS_01.Write();
    parerr_EW_NS_HS_01.Write();
    parerr_FB_NS_HS_01.Write();
    
    parerr_Iso_EW_HS_01.Write();
    parerr_NS_EW_HS_01.Write();
    parerr_EW_EW_HS_01.Write();
    parerr_FB_EW_HS_01.Write();
    
    parerr_Iso_FB_HS_01.Write();
    parerr_NS_FB_HS_01.Write();
    parerr_EW_FB_HS_01.Write();
    parerr_FB_FB_HS_01.Write();
    
    parerr_Iso_NSEW_HS_01.Write();
    parerr_NS_NSEW_HS_01.Write();
    parerr_EW_NSEW_HS_01.Write();
    parerr_FB_NSEW_HS_01.Write();
    
    parerr_Iso_NSFB_HS_01.Write();
    parerr_NS_NSFB_HS_01.Write();
    parerr_EW_NSFB_HS_01.Write();
    parerr_FB_NSFB_HS_01.Write();
    
    parerr_Iso_EWFB_HS_01.Write();
    parerr_NS_EWFB_HS_01.Write();
    parerr_EW_EWFB_HS_01.Write();
    parerr_FB_EWFB_HS_01.Write();
    
    parerr_Iso_full_HS_01.Write();
    parerr_NS_full_HS_01.Write();
    parerr_EW_full_HS_01.Write();
    parerr_FB_full_HS_01.Write();
    
    chi2_Iso_HS_01.Write();
    chi2_NS_HS_01.Write();
    chi2_EW_HS_01.Write();
    chi2_FB_HS_01.Write();
    chi2_NSEW_HS_01.Write();
    chi2_NSFB_HS_01.Write();
    chi2_EWFB_HS_01.Write();
    chi2_full_HS_01.Write();
    
    ndf_Iso_HS_01.Write();
    ndf_NS_HS_01.Write();
    ndf_EW_HS_01.Write();
    ndf_FB_HS_01.Write();
    ndf_NSEW_HS_01.Write();
    ndf_NSFB_HS_01.Write();
    ndf_EWFB_HS_01.Write();
    ndf_full_HS_01.Write();
    
    outFile_HS_01.Write();
    outFile_HS_01.Close();
    

    
}


std::string get_base_path(std::string tree_path,std::string output_dir)
{
    size_t s_pos = tree_path.find_last_of("/");
    size_t f_pos = tree_path.find_last_of("_");
    std::string out_string = output_dir;
    out_string.append("/");
    out_string.append(tree_path.substr(s_pos,(f_pos-s_pos)));
    return out_string;
}

std::string get_final_res_path(std::string base_output_path,std::string ani_level,bool low_statistic)
{
    std::string result_path = base_output_path;
    
    result_path.append("_extracted_results_");
    result_path.append(ani_level);
    if(low_statistic)
        result_path.append("_LS.root");
    else
        result_path.append("_HS.root");
    return result_path;
}

void tree_linking(TTree* myTree,fitResult &tmp_result,UInt_t &tree_entries)
{
    myTree->SetBranchAddress("chi2_LS",&tmp_result.chi2_LS);
    myTree->SetBranchAddress("chi2_r_LS",&tmp_result.chi2_r_LS);
    myTree->SetBranchAddress("ndf_LS",&tmp_result.ndf_LS);
    myTree->SetBranchAddress("fit_par_LS",&tmp_result.fit_par_LS);
    myTree->SetBranchAddress("fit_err_LS",&tmp_result.fit_err_LS);
    myTree->SetBranchAddress("CMatrix_Iso_LS",&tmp_result.CMatrix_Iso_LS);
    myTree->SetBranchAddress("CMatrix_NS_LS",&tmp_result.CMatrix_NS_LS);
    myTree->SetBranchAddress("CMatrix_EW_LS",&tmp_result.CMatrix_EW_LS);
    myTree->SetBranchAddress("CMatrix_FB_LS",&tmp_result.CMatrix_FB_LS);
    myTree->SetBranchAddress("CMatrix_NS_EW_LS",&tmp_result.CMatrix_NS_EW_LS);
    myTree->SetBranchAddress("CMatrix_NS_FB_LS",&tmp_result.CMatrix_NS_FB_LS);
    myTree->SetBranchAddress("CMatrix_EW_FB_LS",&tmp_result.CMatrix_EW_FB_LS);
    myTree->SetBranchAddress("CMatrix_Full_LS",&tmp_result.CMatrix_Full_LS);
    
    myTree->SetBranchAddress("theta_binHistos_LS",&tmp_result.theta_binHistos_LS);
    myTree->SetBranchAddress("phi_binhistos_LS",&tmp_result.phi_binHistos_LS);
    
    myTree->SetBranchAddress("events_LS",&tmp_result.events_LS);
    
    myTree->SetBranchAddress("chi2_HS",&tmp_result.chi2_HS);
    myTree->SetBranchAddress("chi2_r_HS",&tmp_result.chi2_r_HS);
    myTree->SetBranchAddress("ndf_HS",&tmp_result.ndf_HS);
    myTree->SetBranchAddress("fit_par_HS",&tmp_result.fit_par_HS);
    myTree->SetBranchAddress("fit_err_HS",&tmp_result.fit_err_HS);
    myTree->SetBranchAddress("CMatrix_Iso_HS",&tmp_result.CMatrix_Iso_HS);
    myTree->SetBranchAddress("CMatrix_NS_HS",&tmp_result.CMatrix_NS_HS);
    myTree->SetBranchAddress("CMatrix_EW_HS",&tmp_result.CMatrix_EW_HS);
    myTree->SetBranchAddress("CMatrix_FB_HS",&tmp_result.CMatrix_FB_HS);
    myTree->SetBranchAddress("CMatrix_NS_EW_HS",&tmp_result.CMatrix_NS_EW_HS);
    myTree->SetBranchAddress("CMatrix_NS_FB_HS",&tmp_result.CMatrix_NS_FB_HS);
    myTree->SetBranchAddress("CMatrix_EW_FB_HS",&tmp_result.CMatrix_EW_FB_HS);
    myTree->SetBranchAddress("CMatrix_Full_HS",&tmp_result.CMatrix_Full_HS);
    
    myTree->SetBranchAddress("theta_binHistos_HS",&tmp_result.theta_binHistos_HS);
    myTree->SetBranchAddress("phi_binhistos_HS",&tmp_result.phi_binHistos_HS);
    
    myTree->SetBranchAddress("events_HS",&tmp_result.events_HS);
    
    myTree->SetBranchAddress("inputAni",&tmp_result.inputAni);
 
    tree_entries = myTree->GetEntries();
    
}
    
    
    
    


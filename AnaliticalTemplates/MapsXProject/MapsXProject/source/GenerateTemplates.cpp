
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



void compute_templates(
                        std::string output_log,
                        std::string output_root,
                        time_t time_stamp,
                        std::ofstream &output_log_file,
                        TH2D &DAMPE_ReferenceMap_LS,
                        TH2D &DAMPE_ReferenceMap_HS,
                        std::string template_out_path,
                        std::string DAMPE_template_out_path
                       )
{
    
    //////////// Histos stuff
    
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
    
    //////////// Templates Canvas
    
    TCanvas TemplateFunctions("TemplateFulctions","Templates Functions");
    
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
    
    /////////// DAMPE Templates
    
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
    
    
    //////////// Generating All Sky and DAMPE templates
    
    generate_LS_templates(
                              Template_Iso_LS,
                              Template_AniNS_LS,
                              Template_AniEW_LS,
                              Template_AniFB_LS,
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
                              output_log_file,
                              dI,
                              dNS,
                              dEW,
                              dFB
                          );
    
    //////////// Generating DAMPE templates
    
    get_DAMPE_templates(
                            DAMPE_Template_Iso_LS,
                            DAMPE_Template_AniNS_LS,
                            DAMPE_Template_AniEW_LS,
                            DAMPE_Template_AniFB_LS,
                            DAMPE_Template_Iso_HS,
                            DAMPE_Template_AniNS_HS,
                            DAMPE_Template_AniEW_HS,
                            DAMPE_Template_AniFB_HS,
                            DAMPE_ReferenceMap_LS,
                            DAMPE_ReferenceMap_HS,
                            Template_Iso_LS,
                            Template_AniNS_LS,
                            Template_AniEW_LS,
                            Template_AniFB_LS,
                            Template_Iso_HS,
                            Template_AniNS_HS,
                            Template_AniEW_HS,
                            Template_AniFB_HS
                        );

    
    //////////////////////////////////// Creating Template out file
    
    TFile template_file(template_out_path.c_str(),"RECREATE");
    if(template_file.IsZombie()) {
        std::cerr << "\n\nError writing ROOT Templates TFile. Prorgram finished \n\n";
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
    
    TemplateFunctions.Write();
    
    template_file.Write();
    template_file.Close();
    
    //////////////////////////////////// Creating DAMPE Template out file
    
    TFile DAMPE_template_file(DAMPE_template_out_path.c_str(),"RECREATE");
    if(DAMPE_template_file.IsZombie()) {
        std::cerr << "\n\nError writing ROOT DAMPE Templates TFile. Prorgram finished \n\n";
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
    
}

void generate_LS_templates(
                            TH2D &Template_Iso_LS,
                            TH2D &Template_AniNS_LS,
                            TH2D &Template_AniEW_LS,
                            TH2D &Template_AniFB_LS,
                            std::ofstream &log_file,
                            TCanvas &FCanvas,
                            TF2 &dI,
                            TF2 &dNS,
                            TF2 &dEW,
                            TF2 &dFB
                           )

{
    
    Double_t l = 0,b = 0;
    Double_t w_NS = 0, w_EW = 0, w_FB = 0,w_I = 0;
    
    TAxis* theta_axis = nullptr;
    TAxis* phi_axis = nullptr;
    FCanvas.Divide(2,2);
    gStyle->SetPalette(60);
    
    ////////////////// Writing Functions Canvas
    
    FCanvas.cd(1);
    dI.Draw("surf2");
    FCanvas.cd(2);
    dNS.Draw("surf2");
    FCanvas.cd(3);
    dEW.Draw("surf2");
    FCanvas.cd(4);
    dFB.Draw("surf2");
    
    ///////////////////////////////////////////
    
    theta_axis = Template_Iso_LS.GetYaxis();
    phi_axis = Template_Iso_LS.GetXaxis();
    
    for(Int_t Yidx = 1; Yidx <= Template_Iso_LS.GetNbinsY(); Yidx++)
    {
        b = theta_axis->GetBinCenter(Yidx)*TMath::DegToRad();
        for(Int_t Xidx = 1; Xidx <= Template_Iso_LS.GetNbinsX(); Xidx++)
        {
            l = phi_axis->GetBinCenter(Xidx)*TMath::DegToRad();
        
            w_I = dI.Eval(l,b,0);
            w_NS = dNS.Eval(l,b,0);
            w_EW = dEW.Eval(l,b,0);
            w_FB = dFB.Eval(l,b,0);
            
            Template_Iso_LS.Fill(l*TMath::RadToDeg(),b*TMath::RadToDeg(),w_I);
            Template_AniNS_LS.Fill(l*TMath::RadToDeg(),b*TMath::RadToDeg(),w_NS);
            Template_AniEW_LS.Fill(l*TMath::RadToDeg(),b*TMath::RadToDeg(),w_EW);
            Template_AniFB_LS.Fill(l*TMath::RadToDeg(),b*TMath::RadToDeg(),w_FB);
            
        }
    }
    
    for(Int_t bX=1; bX <= Template_Iso_LS.GetNbinsX(); ++bX)
    {
        for(Int_t bY=1; bY <= Template_Iso_LS.GetNbinsY(); ++bY)
        {
            Template_Iso_LS.SetBinError(bX,bY,0);
            Template_AniNS_LS.SetBinError(bX,bY,0);
            Template_AniEW_LS.SetBinError(bX,bY,0);
            Template_AniFB_LS.SetBinError(bX,bY,0);
        }
    }
    
}


void generate_HS_templates(
                            TH2D &Template_Iso_HS,
                            TH2D &Template_AniNS_HS,
                            TH2D &Template_AniEW_HS,
                            TH2D &Template_AniFB_HS,
                            std::ofstream &log_file,
                            TF2 &dI,
                            TF2 &dNS,
                            TF2 &dEW,
                            TF2 &dFB
                           )

{
    
    Double_t l = 0,b = 0;
    Double_t w_NS = 0, w_EW = 0, w_FB = 0,w_I = 0;
    
    TAxis* theta_axis = nullptr;
    TAxis* phi_axis = nullptr;
    
    theta_axis = Template_Iso_HS.GetYaxis();
    phi_axis = Template_Iso_HS.GetXaxis();
    
    for(Int_t Yidx = 1; Yidx <= Template_Iso_HS.GetNbinsY(); Yidx++)
    {
        b = theta_axis->GetBinCenter(Yidx)*TMath::DegToRad();
        for(Int_t Xidx = 1; Xidx <= Template_Iso_HS.GetNbinsX(); Xidx++)
        {
            l = phi_axis->GetBinCenter(Xidx)*TMath::DegToRad();
            
            w_I = dI.Eval(l,b,0);
            w_NS = dNS.Eval(l,b,0);
            w_EW = dEW.Eval(l,b,0);
            w_FB = dFB.Eval(l,b,0);
            
            Template_Iso_HS.Fill(l*TMath::RadToDeg(),b*TMath::RadToDeg(),w_I);
            Template_AniNS_HS.Fill(l*TMath::RadToDeg(),b*TMath::RadToDeg(),w_NS);
            Template_AniEW_HS.Fill(l*TMath::RadToDeg(),b*TMath::RadToDeg(),w_EW);
            Template_AniFB_HS.Fill(l*TMath::RadToDeg(),b*TMath::RadToDeg(),w_FB);
            
        }
    }
    
    for(Int_t bX=1; bX <= Template_Iso_HS.GetNbinsX(); ++bX)
    {
        for(Int_t bY=1; bY <= Template_Iso_HS.GetNbinsY(); ++bY)
        {
            Template_Iso_HS.SetBinError(bX,bY,0);
            Template_AniNS_HS.SetBinError(bX,bY,0);
            Template_AniEW_HS.SetBinError(bX,bY,0);
            Template_AniFB_HS.SetBinError(bX,bY,0);
        }
    }
    
}

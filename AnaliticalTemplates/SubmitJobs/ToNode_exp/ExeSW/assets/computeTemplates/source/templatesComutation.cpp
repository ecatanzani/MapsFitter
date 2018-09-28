
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

void templates_computation(
                            std::ofstream &output_log_file,
                            TH2D &DAMPE_ReferenceMap_LS,
                            TH2D &DAMPE_ReferenceMap_HS
                           )

{
    
    //////////// Costheta flat binning variables
    
    Int_t n_bin_lon_LS = 18;
    Double_t lon_bin_min_LS = -180;
    Double_t lon_bin_max_LS = 180;
    
    Int_t n_bin_lat_LS = 9;
    Double_t lat_bin_min_LS = -90;
    Double_t lat_bin_max_LS = 90;
    
    Double_t* binning_LS = nullptr;
    
    Int_t n_bin_lon_HS = 18;
    Double_t lon_bin_min_HS = -180;
    Double_t lon_bin_max_HS = 180;
    
    Int_t n_bin_lat_HS = 9;
    Double_t lat_bin_min_HS = -90;
    Double_t lat_bin_max_HS = 90;
    
    Double_t* binning_HS = nullptr;
    
    create_binning(n_bin_lat_LS,lat_bin_min_LS,lat_bin_max_LS,binning_LS,true);
    create_binning(n_bin_lat_HS,lat_bin_min_HS,lat_bin_max_HS,binning_HS,true);
    
    //////////////////////////////////////////////////////////// Histos declaration ////////////////////////////////////////////////////
    
    //////////// All Sky Histos
    
    /////////// Low statistics template
    
    TH2D Template_Iso_LS("Template_Iso_LS","Isotropic LS Template All Sky Map; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon_LS,lon_bin_min_LS,lon_bin_max_LS,n_bin_lat_LS,binning_LS);
    TH2D Template_AniNS_LS("Template_AniNS_LS","Anisotropic LS Template All Sky Map (NS); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon_LS,lon_bin_min_LS,lon_bin_max_LS,n_bin_lat_LS,binning_LS);
    TH2D Template_AniEW_LS("Template_AniEW_LS","Anisotropic LS Template All Sky Map (EW); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon_LS,lon_bin_min_LS,lon_bin_max_LS,n_bin_lat_LS,binning_LS);
    TH2D Template_AniFB_LS("Template_AniFB_LS","Anisotropic LS Template All Sky Map (FB); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon_LS,lon_bin_min_LS,lon_bin_max_LS,n_bin_lat_LS,binning_LS);
    
    Template_Iso_LS.Sumw2();
    Template_AniNS_LS.Sumw2();
    Template_AniEW_LS.Sumw2();
    Template_AniFB_LS.Sumw2();
    
    /////////// High statistics template
    
    TH2D Template_Iso_HS("Template_Iso_HS","Isotropic HS Template All Sky Map; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon_HS,lon_bin_min_HS,lon_bin_max_HS,n_bin_lat_HS,binning_HS);
    TH2D Template_AniNS_HS("Template_AniNS_HS","Anisotropic HS Template All Sky Map (NS); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon_HS,lon_bin_min_HS,lon_bin_max_HS,n_bin_lat_HS,binning_HS);
    TH2D Template_AniEW_HS("Template_AniEW_HS","Anisotropic HS Template All Sky Map (EW); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon_HS,lon_bin_min_HS,lon_bin_max_HS,n_bin_lat_HS,binning_HS);
    TH2D Template_AniFB_HS("Template_AniFB_HS","Anisotropic HS Template All Sky Map (FB); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon_HS,lon_bin_min_HS,lon_bin_max_HS,n_bin_lat_HS,binning_HS);
    
    Template_Iso_HS.Sumw2();
    Template_AniNS_HS.Sumw2();
    Template_AniEW_HS.Sumw2();
    Template_AniFB_HS.Sumw2();
    
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
    
    ///////////////////////// Generating templates with a certain anisotropy level
    
    generate_LS_templates(
                            Template_Iso_LS,
                            Template_AniNS_LS,
                            Template_AniEW_LS,
                            Template_AniFB_LS,
                            output_log_file,
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
    
    
    
    TFile template_results(templates_path.c_str(),"RECREATE");
    
    if(template_results.IsZombie())
    {
        std::cerr << "\n\nError writing output ROOT All Sky template out file \n\n";
        exit(100);
    }
    
    Template_Iso_LS.Write();
    Template_AniNS_LS.Write();
    Template_AniEW_LS.Write();
    Template_AniFB_LS.Write();
    
    Template_Iso_HS.Write();
    Template_AniNS_HS.Write();
    Template_AniEW_HS.Write();
    Template_AniFB_HS.Write();
    
    dI.Write();
    dNS.Write();
    dEW.Write();
    dFB.Write();
    
    template_results.Write();
    template_results.Close();
    
    
    TFile DAMPE_template_results(DAMPE_templates_path.c_str(),"RECREATE");
    
    if(DAMPE_template_results.IsZombie())
    {
        std::cerr << "\n\nError writing output ROOT DAMPE template out file \n\n";
        exit(100);
    }
    
    DAMPE_Template_Iso_LS.Write();
    DAMPE_Template_AniNS_LS.Write();
    DAMPE_Template_AniEW_LS.Write();
    DAMPE_Template_AniFB_LS.Write();
    
    DAMPE_Template_Iso_HS.Write();
    DAMPE_Template_AniNS_HS.Write();
    DAMPE_Template_AniEW_HS.Write();
    DAMPE_Template_AniFB_HS.Write();
    
    DAMPE_template_results.Write();
    DAMPE_template_results.Close();
    
}





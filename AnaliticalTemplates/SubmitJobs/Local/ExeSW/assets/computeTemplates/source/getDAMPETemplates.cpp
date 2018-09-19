
#include "MyHead.h"

void get_DAMPE_templates(
                            TH2D &DAMPE_Template_Iso_LS,
                            TH2D &DAMPE_Template_AniNS_LS,
                            TH2D &DAMPE_Template_AniEW_LS,
                            TH2D &DAMPE_Template_AniFB_LS,
                            TH2D &DAMPE_Template_Iso_HS,
                            TH2D &DAMPE_Template_AniNS_HS,
                            TH2D &DAMPE_Template_AniEW_HS,
                            TH2D &DAMPE_Template_AniFB_HS,
                            TH2D &DAMPE_ReferenceMap_LS,
                            TH2D &DAMPE_ReferenceMap_HS,
                            TH2D &Template_Iso_LS,
                            TH2D &Template_AniNS_LS,
                            TH2D &Template_AniEW_LS,
                            TH2D &Template_AniFB_LS,
                            TH2D &Template_Iso_HS,
                            TH2D &Template_AniNS_HS,
                            TH2D &Template_AniEW_HS,
                            TH2D &Template_AniFB_HS
                        )
{
    
    TH2D* pTemplate_Iso_LS = &Template_Iso_LS;
    TH2D* pTemplate_AniNS_LS = &Template_AniNS_LS;
    TH2D* pTemplate_AniEW_LS = &Template_AniEW_LS;
    TH2D* pTemplate_AniFB_LS = &Template_AniFB_LS;
    
    TH2D* pTemplate_Iso_HS = &Template_Iso_HS;
    TH2D* pTemplate_AniNS_HS = &Template_AniNS_HS;
    TH2D* pTemplate_AniEW_HS = &Template_AniEW_HS;
    TH2D* pTemplate_AniFB_HS = &Template_AniFB_HS;
    
    TH2D* pDAMPE_ReferenceMap_LS = &DAMPE_ReferenceMap_LS;
    TH2D* pDAMPE_ReferenceMap_HS = &DAMPE_ReferenceMap_HS;
    
    TH2D* tmp_DAMPE_Template_Iso_LS = (TH2D*) pDAMPE_ReferenceMap_LS->Clone("tmp_DAMPE_Template_Iso_LS");
    TH2D* tmp_DAMPE_Template_AniNS_LS = (TH2D*) pDAMPE_ReferenceMap_LS->Clone("tmp_DAMPE_Template_AniNS_LS");
    TH2D* tmp_DAMPE_Template_AniEW_LS = (TH2D*) pDAMPE_ReferenceMap_LS->Clone("tmp_DAMPE_Template_AniEW_LS");
    TH2D* tmp_DAMPE_Template_AniFB_LS = (TH2D*) pDAMPE_ReferenceMap_LS->Clone("tmp_DAMPE_Template_AniFB_LS");
    
    TH2D* tmp_DAMPE_Template_Iso_HS = (TH2D*) pDAMPE_ReferenceMap_HS->Clone("tmp_DAMPE_Template_Iso_HS");
    TH2D* tmp_DAMPE_Template_AniNS_HS = (TH2D*) pDAMPE_ReferenceMap_HS->Clone("tmp_DAMPE_Template_AniNS_HS");
    TH2D* tmp_DAMPE_Template_AniEW_HS = (TH2D*) pDAMPE_ReferenceMap_HS->Clone("tmp_DAMPE_Template_AniEW_HS");
    TH2D* tmp_DAMPE_Template_AniFB_HS = (TH2D*) pDAMPE_ReferenceMap_HS->Clone("tmp_DAMPE_Template_AniFB_HS");
    
    tmp_DAMPE_Template_Iso_LS->SetTitle("Isotropic Template LS DAMPE Sky Map; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries");
    tmp_DAMPE_Template_AniNS_LS->SetTitle("Anisotropic Template LS DAMPE Sky Map (NS); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries");
    tmp_DAMPE_Template_AniEW_LS->SetTitle("Anisotropic Template LS DAMPE Sky Map (EW); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries");
    tmp_DAMPE_Template_AniFB_LS->SetTitle("Anisotropic Template LS DAMPE Sky Map (FB); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries");
    
    tmp_DAMPE_Template_Iso_HS->SetTitle("Isotropic Template HS DAMPE Sky Map; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries");
    tmp_DAMPE_Template_AniNS_HS->SetTitle("Anisotropic Template HS DAMPE Sky Map (NS); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries");
    tmp_DAMPE_Template_AniEW_HS->SetTitle("Anisotropic Template HS DAMPE Sky Map (EW); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries");
    tmp_DAMPE_Template_AniFB_HS->SetTitle("Anisotropic Template HS DAMPE Sky Map (FB); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries");
    
    tmp_DAMPE_Template_Iso_LS->Reset();
    tmp_DAMPE_Template_AniNS_LS->Reset();
    tmp_DAMPE_Template_AniEW_LS->Reset();
    tmp_DAMPE_Template_AniFB_LS->Reset();
    
    tmp_DAMPE_Template_Iso_HS->Reset();
    tmp_DAMPE_Template_AniNS_HS->Reset();
    tmp_DAMPE_Template_AniEW_HS->Reset();
    tmp_DAMPE_Template_AniFB_HS->Reset();
    
    tmp_DAMPE_Template_Iso_LS->Multiply(pTemplate_Iso_LS,pDAMPE_ReferenceMap_LS);
    tmp_DAMPE_Template_AniNS_LS->Multiply(pTemplate_AniNS_LS,pDAMPE_ReferenceMap_LS);
    tmp_DAMPE_Template_AniEW_LS->Multiply(pTemplate_AniEW_LS,pDAMPE_ReferenceMap_LS);
    tmp_DAMPE_Template_AniFB_LS->Multiply(pTemplate_AniFB_LS,pDAMPE_ReferenceMap_LS);
    
    tmp_DAMPE_Template_Iso_HS->Multiply(pTemplate_Iso_HS,pDAMPE_ReferenceMap_HS);
    tmp_DAMPE_Template_AniNS_HS->Multiply(pTemplate_AniNS_HS,pDAMPE_ReferenceMap_HS);
    tmp_DAMPE_Template_AniEW_HS->Multiply(pTemplate_AniEW_HS,pDAMPE_ReferenceMap_HS);
    tmp_DAMPE_Template_AniFB_HS->Multiply(pTemplate_AniFB_HS,pDAMPE_ReferenceMap_HS);
    
    /*
    ////////////////// Scaling DAMPE templates for sqrt(4*pi) to obtain comparable fit parameters
    
    tmp_DAMPE_Template_Iso_LS->Scale(TMath::Sqrt(4*TMath::Pi()));
    tmp_DAMPE_Template_AniNS_LS->Scale(TMath::Sqrt(4*TMath::Pi()));
    tmp_DAMPE_Template_AniEW_LS->Scale(TMath::Sqrt(4*TMath::Pi()));
    tmp_DAMPE_Template_AniFB_LS->Scale(TMath::Sqrt(4*TMath::Pi()));
    
    tmp_DAMPE_Template_Iso_HS->Scale(TMath::Sqrt(4*TMath::Pi()));
    tmp_DAMPE_Template_AniNS_HS->Scale(TMath::Sqrt(4*TMath::Pi()));
    tmp_DAMPE_Template_AniEW_HS->Scale(TMath::Sqrt(4*TMath::Pi()));
    tmp_DAMPE_Template_AniFB_HS->Scale(TMath::Sqrt(4*TMath::Pi()));
    */
     
    ////////////////// Linking with main histos
    
    new (&DAMPE_Template_Iso_LS) (TH2D) (*(TH2D*)tmp_DAMPE_Template_Iso_LS->Clone("DAMPE_Template_Iso_LS"));
    new (&DAMPE_Template_AniNS_LS) (TH2D) (*(TH2D*)tmp_DAMPE_Template_AniNS_LS->Clone("DAMPE_Template_AniNS_LS"));
    new (&DAMPE_Template_AniEW_LS) (TH2D) (*(TH2D*)tmp_DAMPE_Template_AniEW_LS->Clone("DAMPE_Template_AniEW_LS"));
    new (&DAMPE_Template_AniFB_LS) (TH2D) (*(TH2D*)tmp_DAMPE_Template_AniFB_LS->Clone("DAMPE_Template_AniFB_LS"));
    
    new (&DAMPE_Template_Iso_HS) (TH2D) (*(TH2D*)tmp_DAMPE_Template_Iso_HS->Clone("DAMPE_Template_Iso_HS"));
    new (&DAMPE_Template_AniNS_HS) (TH2D) (*(TH2D*)tmp_DAMPE_Template_AniNS_HS->Clone("DAMPE_Template_AniNS_HS"));
    new (&DAMPE_Template_AniEW_HS) (TH2D) (*(TH2D*)tmp_DAMPE_Template_AniEW_HS->Clone("DAMPE_Template_AniEW_HS"));
    new (&DAMPE_Template_AniFB_HS) (TH2D) (*(TH2D*)tmp_DAMPE_Template_AniFB_HS->Clone("DAMPE_Template_AniFB_HS"));
    
}


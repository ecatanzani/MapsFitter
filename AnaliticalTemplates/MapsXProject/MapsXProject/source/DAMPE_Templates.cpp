
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
                            TH2D &DAMPE_ReferenceMap,
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
    TH2D* pDAMPE_ReferenceMap = &DAMPE_ReferenceMap;
    
    TH2D* tmp_DAMPE_Template_Iso_LS = (TH2D*) pDAMPE_ReferenceMap->Clone("tmp_DAMPE_Template_Iso_LS");
    TH2D* tmp_DAMPE_Template_AniNS_LS = (TH2D*) pDAMPE_ReferenceMap->Clone("tmp_DAMPE_Template_AniNS_LS");
    TH2D* tmp_DAMPE_Template_AniEW_LS = (TH2D*) pDAMPE_ReferenceMap->Clone("tmp_DAMPE_Template_AniEW_LS");
    TH2D* tmp_DAMPE_Template_AniFB_LS = (TH2D*) pDAMPE_ReferenceMap->Clone("tmp_DAMPE_Template_AniFB_LS");
    
    TH2D* tmp_DAMPE_Template_Iso_HS = (TH2D*) pDAMPE_ReferenceMap->Clone("tmp_DAMPE_Template_Iso_HS");
    TH2D* tmp_DAMPE_Template_AniNS_HS = (TH2D*) pDAMPE_ReferenceMap->Clone("tmp_DAMPE_Template_AniNS_HS");
    TH2D* tmp_DAMPE_Template_AniEW_HS = (TH2D*) pDAMPE_ReferenceMap->Clone("tmp_DAMPE_Template_AniEW_HS");
    TH2D* tmp_DAMPE_Template_AniFB_HS = (TH2D*) pDAMPE_ReferenceMap->Clone("tmp_DAMPE_Template_AniFB_HS");
    
    tmp_DAMPE_Template_Iso_LS->Reset();
    tmp_DAMPE_Template_AniNS_LS->Reset();
    tmp_DAMPE_Template_AniEW_LS->Reset();
    tmp_DAMPE_Template_AniFB_LS->Reset();
    
    tmp_DAMPE_Template_Iso_HS->Reset();
    tmp_DAMPE_Template_AniNS_HS->Reset();
    tmp_DAMPE_Template_AniEW_HS->Reset();
    tmp_DAMPE_Template_AniFB_HS->Reset();
    
    tmp_DAMPE_Template_Iso_LS->SetTitle("Isotropic Template LS DAMPE Sky Map; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries");
    tmp_DAMPE_Template_AniNS_LS->SetTitle("Anisotropic Template LS DAMPE Sky Map (NS); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries");
    tmp_DAMPE_Template_AniEW_LS->SetTitle("Anisotropic Template LS DAMPE Sky Map (EW); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries");
    tmp_DAMPE_Template_AniFB_LS->SetTitle("Anisotropic Template LS DAMPE Sky Map (FB); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries");
    
    tmp_DAMPE_Template_Iso_HS->SetTitle("Isotropic Template HS DAMPE Sky Map; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries");
    tmp_DAMPE_Template_AniNS_HS->SetTitle("Anisotropic Template HS DAMPE Sky Map (NS); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries");
    tmp_DAMPE_Template_AniEW_HS->SetTitle("Anisotropic Template HS DAMPE Sky Map (EW); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries");
    tmp_DAMPE_Template_AniFB_HS->SetTitle("Anisotropic Template HS DAMPE Sky Map (FB); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries");
    
    
    
    for(Int_t bX = 1; bX <= pDAMPE_ReferenceMap->GetNbinsX(); ++bX)
        for(Int_t bY = 1; bY <= pDAMPE_ReferenceMap->GetNbinsY(); ++bY)
        {
            tmp_DAMPE_Template_Iso_LS->SetBinContent(bX,bY,pDAMPE_ReferenceMap->GetBinContent(bX,bY)*Template_Iso_LS.GetBinContent(bX,bY));
            tmp_DAMPE_Template_AniNS_LS->SetBinContent(bX,bY,pDAMPE_ReferenceMap->GetBinContent(bX,bY)*Template_AniNS_LS.GetBinContent(bX,bY));
            tmp_DAMPE_Template_AniEW_LS->SetBinContent(bX,bY,pDAMPE_ReferenceMap->GetBinContent(bX,bY)*Template_AniEW_LS.GetBinContent(bX,bY));
            tmp_DAMPE_Template_AniFB_LS->SetBinContent(bX,bY,pDAMPE_ReferenceMap->GetBinContent(bX,bY)*Template_AniFB_LS.GetBinContent(bX,bY));
            
            tmp_DAMPE_Template_Iso_HS->SetBinContent(bX,bY,pDAMPE_ReferenceMap->GetBinContent(bX,bY)*Template_Iso_HS.GetBinContent(bX,bY));
            tmp_DAMPE_Template_AniNS_HS->SetBinContent(bX,bY,pDAMPE_ReferenceMap->GetBinContent(bX,bY)*Template_AniNS_HS.GetBinContent(bX,bY));
            tmp_DAMPE_Template_AniEW_HS->SetBinContent(bX,bY,pDAMPE_ReferenceMap->GetBinContent(bX,bY)*Template_AniEW_HS.GetBinContent(bX,bY));
            tmp_DAMPE_Template_AniFB_HS->SetBinContent(bX,bY,pDAMPE_ReferenceMap->GetBinContent(bX,bY)*Template_AniFB_HS.GetBinContent(bX,bY));
        }
    
    normalize_DAMPE_templates(
                                tmp_DAMPE_Template_Iso_LS,
                                tmp_DAMPE_Template_AniNS_LS,
                                tmp_DAMPE_Template_AniEW_LS,
                                tmp_DAMPE_Template_AniFB_LS,
                                tmp_DAMPE_Template_Iso_HS,
                                tmp_DAMPE_Template_AniNS_HS,
                                tmp_DAMPE_Template_AniEW_HS,
                                tmp_DAMPE_Template_AniFB_HS
                              );
    
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

void normalize_DAMPE_templates(
                                TH2D* tmp_DAMPE_Template_Iso_LS,
                                TH2D* tmp_DAMPE_Template_AniNS_LS,
                                TH2D* tmp_DAMPE_Template_AniEW_LS,
                                TH2D* tmp_DAMPE_Template_AniFB_LS,
                                TH2D* tmp_DAMPE_Template_Iso_HS,
                                TH2D* tmp_DAMPE_Template_AniNS_HS,
                                TH2D* tmp_DAMPE_Template_AniEW_HS,
                                TH2D* tmp_DAMPE_Template_AniFB_HS
                               )

{
    
    for(Int_t bX = 1; bX <= tmp_DAMPE_Template_Iso_LS->GetNbinsX(); ++bX)
        for(Int_t bY = 1; bY <= tmp_DAMPE_Template_Iso_LS->GetNbinsY(); ++bY)
        {
            tmp_DAMPE_Template_Iso_LS->SetBinContent(bX,bY,tmp_DAMPE_Template_Iso_LS->GetBinContent(bX,bY)/data_all_sky_LS_events);
            tmp_DAMPE_Template_AniNS_LS->SetBinContent(bX,bY,tmp_DAMPE_Template_AniNS_LS->GetBinContent(bX,bY)/data_all_sky_LS_events);
            tmp_DAMPE_Template_AniEW_LS->SetBinContent(bX,bY,tmp_DAMPE_Template_AniEW_LS->GetBinContent(bX,bY)/data_all_sky_LS_events);
            tmp_DAMPE_Template_AniFB_LS->SetBinContent(bX,bY,tmp_DAMPE_Template_AniFB_LS->GetBinContent(bX,bY)/data_all_sky_LS_events);
            
            tmp_DAMPE_Template_Iso_HS->SetBinContent(bX,bY,tmp_DAMPE_Template_Iso_HS->GetBinContent(bX,bY)/data_all_sky_HS_events);
            tmp_DAMPE_Template_AniNS_HS->SetBinContent(bX,bY,tmp_DAMPE_Template_AniNS_HS->GetBinContent(bX,bY)/data_all_sky_HS_events);
            tmp_DAMPE_Template_AniEW_HS->SetBinContent(bX,bY,tmp_DAMPE_Template_AniEW_HS->GetBinContent(bX,bY)/data_all_sky_HS_events);
            tmp_DAMPE_Template_AniFB_HS->SetBinContent(bX,bY,tmp_DAMPE_Template_AniFB_HS->GetBinContent(bX,bY)/data_all_sky_HS_events);
        }
    
    
    
    
    /*          WRONG WAY TO NORMALIZE -> All the templates have to be normalized to the same value, otherwise the proportions between fitted values are not still valid.
    
    Double_t normIso = 1/TMath::Sqrt(compute_integral(tmp_DAMPE_Template_Iso));
    Double_t normNS = 1/TMath::Sqrt(compute_integral(tmp_DAMPE_Template_AniNS));
    Double_t normEW = 1/TMath::Sqrt(compute_integral(tmp_DAMPE_Template_AniEW));
    Double_t normFB = 1/TMath::Sqrt(compute_integral(tmp_DAMPE_Template_AniFB));
    
     for(Int_t bX = 1; bX <= tmp_DAMPE_Template_Iso->GetNbinsX(); ++bX)
        for(Int_t bY = 1; bY <= tmp_DAMPE_Template_Iso->GetNbinsY(); ++bY)
        {
            tmp_DAMPE_Template_Iso->SetBinContent(bX,bY,tmp_DAMPE_Template_Iso->GetBinContent(bX,bY)*normIso);
            tmp_DAMPE_Template_AniNS->SetBinContent(bX,bY,tmp_DAMPE_Template_AniNS->GetBinContent(bX,bY)*normNS);
            tmp_DAMPE_Template_AniEW->SetBinContent(bX,bY,tmp_DAMPE_Template_AniEW->GetBinContent(bX,bY)*normEW);
            tmp_DAMPE_Template_AniFB->SetBinContent(bX,bY,tmp_DAMPE_Template_AniFB->GetBinContent(bX,bY)*normFB);
        }
    
    std::cout << "\n\nIntegral Iso DAMPE template: " << compute_integral(tmp_DAMPE_Template_Iso);
    std::cout << "\nIntegral NS DAMPE template: " << compute_integral(tmp_DAMPE_Template_AniNS);
    std::cout << "\nIntegral EW DAMPE template: " << compute_integral(tmp_DAMPE_Template_AniEW);
    std::cout << "\nIntegral NS DAMPE template: " << compute_integral(tmp_DAMPE_Template_AniFB) << "\n\n";
    
    */
     
}

Double_t compute_integral(TH2D* histo)
{
    TAxis* theta_axis = histo->GetYaxis();
    TAxis* phi_axis = histo->GetXaxis();
    
    Double_t sum = 0;
    Double_t theta = 0, phi = 0;
    Double_t jacob = 0;
    
    for(Int_t bX = 1; bX <= histo->GetNbinsX(); ++bX)
    {
        phi = (phi_axis->GetBinWidth(bX))*TMath::DegToRad();
        
        for(Int_t bY = 1; bY <= histo->GetNbinsY(); ++bY)
        {
            theta = (theta_axis->GetBinWidth(bY))*TMath::DegToRad();
            jacob = TMath::Cos(theta_axis->GetBinCenter(bY)*TMath::DegToRad());
            sum += jacob*theta*phi*TMath::Power(histo->GetBinContent(bX,bY),2);
        }
    }
    
    return sum;
    
}


#include "MyHead.h"

void generate_LS_templates(
                            TH2D &Template_Iso_LS,
                            TH2D &Template_AniNS_LS,
                            TH2D &Template_AniEW_LS,
                            TH2D &Template_AniFB_LS,
                            std::ofstream &output_log_file,
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
    
    ///////////////////////////////////////////
    
    theta_axis = Template_Iso_LS.GetYaxis();
    phi_axis = Template_Iso_LS.GetXaxis();
    
    for(Int_t Yidx = 1; Yidx <= Template_Iso_LS.GetNbinsY(); Yidx++) {
        b = theta_axis->GetBinCenter(Yidx)*TMath::DegToRad();
        for(Int_t Xidx = 1; Xidx <= Template_Iso_LS.GetNbinsX(); Xidx++) {
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
                            std::ofstream &output_log_file,
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
    
    for(Int_t Yidx = 1; Yidx <= Template_Iso_HS.GetNbinsY(); Yidx++) {
        b = theta_axis->GetBinCenter(Yidx)*TMath::DegToRad();
        for(Int_t Xidx = 1; Xidx <= Template_Iso_HS.GetNbinsX(); Xidx++) {
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


#include "MyHead.h"

void generate_LS_templates(
                            TH2D &Template_Iso_LS,
                            TH2D &Template_AniNS_LS,
                            TH2D &Template_AniEW_LS,
                            TH2D &Template_AniFB_LS,
                            TH1D &Template_hwI,
                            TH1D &Template_hwNS,
                            TH1D &Template_hwEW,
                            TH1D &Template_hwFB,
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
    
    std::cout<<"\n\n";
    
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
    
    for(Int_t Yidx = 1; Yidx <= Template_Iso_LS.GetNbinsY(); Yidx++) {
        b = theta_axis->GetBinCenter(Yidx)*TMath::DegToRad();
        for(Int_t Xidx = 1; Xidx <= Template_Iso_LS.GetNbinsX(); Xidx++) {
            l = phi_axis->GetBinCenter(Xidx)*TMath::DegToRad();
        
            w_I = dI.Eval(l,b,0);
            w_NS = dNS.Eval(l,b,0);
            w_EW = dEW.Eval(l,b,0);
            w_FB = dFB.Eval(l,b,0);
        
            Template_hwI.Fill(dI.Eval(l,b,0));
            Template_hwNS.Fill(dNS.Eval(l,b,0));
            Template_hwEW.Fill(dEW.Eval(l,b,0));
            Template_hwFB.Fill(dFB.Eval(l,b,0));
            
            Template_Iso_LS.Fill(l*TMath::RadToDeg(),b*TMath::RadToDeg(),w_I);
            Template_AniNS_LS.Fill(l*TMath::RadToDeg(),b*TMath::RadToDeg(),w_NS);
            Template_AniEW_LS.Fill(l*TMath::RadToDeg(),b*TMath::RadToDeg(),w_EW);
            Template_AniFB_LS.Fill(l*TMath::RadToDeg(),b*TMath::RadToDeg(),w_FB);
            
        }
    }
}


void generate_HS_templates(
                            TH2D &Template_Iso_HS,
                            TH2D &Template_AniNS_HS,
                            TH2D &Template_AniEW_HS,
                            TH2D &Template_AniFB_HS,
                            TH1D &Template_hwI,
                            TH1D &Template_hwNS,
                            TH1D &Template_hwEW,
                            TH1D &Template_hwFB,
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
    std::cout<<"\n\n";
    
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
            
            Template_hwI.Fill(dI.Eval(l,b,0));
            Template_hwNS.Fill(dNS.Eval(l,b,0));
            Template_hwEW.Fill(dEW.Eval(l,b,0));
            Template_hwFB.Fill(dFB.Eval(l,b,0));
            
            Template_Iso_HS.Fill(l*TMath::RadToDeg(),b*TMath::RadToDeg(),w_I);
            Template_AniNS_HS.Fill(l*TMath::RadToDeg(),b*TMath::RadToDeg(),w_NS);
            Template_AniEW_HS.Fill(l*TMath::RadToDeg(),b*TMath::RadToDeg(),w_EW);
            Template_AniFB_HS.Fill(l*TMath::RadToDeg(),b*TMath::RadToDeg(),w_FB);
            
        }
    }
}

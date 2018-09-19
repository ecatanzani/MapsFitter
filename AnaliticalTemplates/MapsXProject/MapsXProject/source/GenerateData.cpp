
#include "MyHead.h"

void generate_LS_data(
                        Double_t NS_anisotropy,
                        Double_t EW_anisotropy,
                        Double_t FB_anisotropy,
                        TH2D* Data_Iso_LS,
                        TH2D* Data_AniNS_LS,
                        TH2D* Data_AniEW_LS,
                        TH2D* Data_AniFB_LS,
                        TH2D* MixedData_NS_EW_LS,
                        TH2D* MixedData_NS_FB_LS,
                        TH2D* MixedData_EW_FB_LS,
                        TH2D* FullMixedData_LS,
                        TH1D &Data_hwI,
                        TH1D &Data_hwNS,
                        TH1D &Data_hwEW,
                        TH1D &Data_hwFB,
                        TF2 &dI,
                        TF2 &dNS,
                        TF2 &dEW,
                        TF2 &dFB,
                        TH2D &Template_Iso_LS,
                        TH2D &Template_AniNS_LS,
                        TH2D &Template_AniEW_LS,
                        TH2D &Template_AniFB_LS,
                        std::ofstream &log_file,
                        TRandom3 &r_gen
                      )

{
    
    std::cout << "\n\n Building LS Data maps\n\n";
        
    ///////////// Pointers to template histos
        
    TH2D* pTemplate_Data_Iso_LS = &Template_Iso_LS;
    TH2D* Template_Data_AniNS_LS = (TH2D*)Template_Iso_LS.Clone("Template_Data_AniNS_LS");
    TH2D* Template_Data_AniEW_LS = (TH2D*)Template_Iso_LS.Clone("Template_Data_AniEW_LS");
    TH2D* Template_Data_AniFB_LS = (TH2D*)Template_Iso_LS.Clone("Template_Data_AniFB_LS");
    TH2D* Template_MixedData_NS_EW_LS = (TH2D*)Template_Iso_LS.Clone("Template_MixedData_NS_EW_LS");
    TH2D* Template_MixedData_NS_FB_LS = (TH2D*)Template_Iso_LS.Clone("Template_MixedData_NS_FB_LS");
    TH2D* Template_MixedData_EW_FB_LS = (TH2D*)Template_Iso_LS.Clone("Template_MixedData_EW_FB_LS");
    TH2D* Template_FullMixedData_LS = (TH2D*)Template_Iso_LS.Clone("Template_FullMixedData_LS");
    
    Template_Data_AniNS_LS->Reset();
    Template_Data_AniEW_LS->Reset();
    Template_Data_AniFB_LS->Reset();
    Template_MixedData_NS_EW_LS->Reset();
    Template_MixedData_NS_FB_LS->Reset();
    Template_FullMixedData_LS->Reset();
    
    for(Int_t bX = 1; bX <= Template_Data_AniNS_LS->GetNbinsX(); ++bX)
    {
        for(Int_t bY = 1; bY <= Template_Data_AniNS_LS->GetNbinsY(); ++bY)
        {
            Template_Data_AniNS_LS->SetBinContent(bX,bY, (1-NS_anisotropy)*Template_Iso_LS.GetBinContent(bX,bY) + NS_anisotropy*Template_AniNS_LS.GetBinContent(bX,bY) );
            Template_Data_AniEW_LS->SetBinContent(bX,bY, (1-EW_anisotropy)*Template_Iso_LS.GetBinContent(bX,bY) + EW_anisotropy*Template_AniEW_LS.GetBinContent(bX,bY) );
            Template_Data_AniFB_LS->SetBinContent(bX,bY, (1-FB_anisotropy)*Template_Iso_LS.GetBinContent(bX,bY) + FB_anisotropy*Template_AniFB_LS.GetBinContent(bX,bY) );
            Template_MixedData_NS_EW_LS->SetBinContent(bX,bY, (1- NS_anisotropy - EW_anisotropy)*Template_Iso_LS.GetBinContent(bX,bY) + NS_anisotropy*Template_AniNS_LS.GetBinContent(bX,bY) + EW_anisotropy*Template_AniEW_LS.GetBinContent(bX,bY) );
            Template_MixedData_NS_FB_LS->SetBinContent(bX,bY, (1- NS_anisotropy - FB_anisotropy)*Template_Iso_LS.GetBinContent(bX,bY) + NS_anisotropy*Template_AniNS_LS.GetBinContent(bX,bY) + FB_anisotropy*Template_AniFB_LS.GetBinContent(bX,bY) );
            Template_MixedData_EW_FB_LS->SetBinContent(bX,bY, (1- EW_anisotropy - FB_anisotropy)*Template_Iso_LS.GetBinContent(bX,bY) + EW_anisotropy*Template_AniEW_LS.GetBinContent(bX,bY) + FB_anisotropy*Template_AniFB_LS.GetBinContent(bX,bY) );
            Template_FullMixedData_LS->SetBinContent(bX,bY, (1 - NS_anisotropy - EW_anisotropy - FB_anisotropy)*Template_Iso_LS.GetBinContent(bX,bY) + NS_anisotropy*Template_AniNS_LS.GetBinContent(bX,bY) + EW_anisotropy*Template_AniEW_LS.GetBinContent(bX,bY) + FB_anisotropy*Template_AniFB_LS.GetBinContent(bX,bY) );
        }
    }
        
    ///////////// Filling maps...
        
    Data_Iso_LS->FillRandom(pTemplate_Data_Iso_LS,data_all_sky_LS_events);
    Data_AniNS_LS->FillRandom(Template_Data_AniNS_LS,data_all_sky_LS_events);
    Data_AniEW_LS->FillRandom(Template_Data_AniEW_LS,data_all_sky_LS_events);
    Data_AniFB_LS->FillRandom(Template_Data_AniFB_LS,data_all_sky_LS_events);
    MixedData_NS_EW_LS->FillRandom(Template_MixedData_NS_EW_LS,data_all_sky_LS_events);
    MixedData_NS_FB_LS->FillRandom(Template_MixedData_NS_FB_LS,data_all_sky_LS_events);
    MixedData_EW_FB_LS->FillRandom(Template_MixedData_EW_FB_LS,data_all_sky_LS_events);
    FullMixedData_LS->FillRandom(Template_FullMixedData_LS,data_all_sky_LS_events);
   
}


/*
 
void generate_LS_data(
                      TH2D* Data_Iso_LS,
                      TH2D* Data_AniNS_LS,
                      TH2D* Data_AniEW_LS,
                      TH2D* Data_AniFB_LS,
                      TH2D* MixedData_NS_EW_LS,
                      TH2D* MixedData_NS_FB_LS,
                      TH2D* MixedData_EW_FB_LS,
                      TH2D* FullMixedData_LS,
                      TH1D &Data_hwI,
                      TH1D &Data_hwNS,
                      TH1D &Data_hwEW,
                      TH1D &Data_hwFB,
                      TF2 &dI,
                      TF2 &dNS,
                      TF2 &dEW,
                      TF2 &dFB,
                      TH2D &Template_Iso_LS,
                      TH2D &Template_AniNS_LS,
                      TH2D &Template_AniEW_LS,
                      TH2D &Template_AniFB_LS,
                      std::ofstream &log_file,
                      TRandom3 &r_gen
                      )

{

    Double_t coin = 0;
    Double_t l = 0,b = 0, cosb = 0;

    /////////////// Template probabilities...

    Double_t probI = 0;
    Double_t probNS = 0;
    Double_t probEW = 0;
    Double_t probFB = 0;

    std::cout << "\n\nBuilding LS data maps...";
    boost::progress_display progress(data_all_sky_LS_events);

    for(unsigned long int ev_idx = 0; ev_idx < data_all_sky_LS_events; ++ev_idx)
    {
    
        l = r_gen.Uniform(-180,180);
        cosb = r_gen.Uniform(-1,1);
        b = TMath::RadToDeg()*TMath::ACos(cosb)-90;
        coin = r_gen.Uniform(0,1);
    
        l *= TMath::DegToRad();
        b *= TMath::DegToRad();
    
        probI = dI.Eval(l,b,0);
        probNS = dNS.Eval(l,b,0);
        probEW = dEW.Eval(l,b,0);
        probFB = dFB.Eval(l,b,0);
    
        Data_hwI.Fill(probI);
        Data_hwNS.Fill(probNS);
        Data_hwEW.Fill(probEW);
        Data_hwFB.Fill(probFB);
    
        l *= TMath::RadToDeg();
        b *= TMath::RadToDeg();
    
        /////////////// Filling histos (Hit & Miss method)
    
    
        if( TMath::Power(probI,2) > coin)
            Data_Iso_LS->Fill(l,b,1);
        if( (TMath::Power(TMath::Sqrt(0.9)*probI + TMath::Sqrt(0.1)*probNS,2)) > coin )
            Data_AniNS_LS->Fill(l,b,1);
        if( (TMath::Power(TMath::Sqrt(0.9)*probI + TMath::Sqrt(0.1)*probEW,2)) > coin )
            Data_AniEW_LS->Fill(l,b,1);
        if( (TMath::Power(TMath::Sqrt(0.9)*probI + TMath::Sqrt(0.1)*probFB,2)) > coin )
            Data_AniFB_LS->Fill(l,b,1);
        if( (TMath::Power(TMath::Sqrt(0.8)*probI + TMath::Sqrt(0.1)*probNS + TMath::Sqrt(0.1)*probEW,2)) > coin )
            MixedData_NS_EW_LS->Fill(l,b,1);
        if( (TMath::Power(TMath::Sqrt(0.8)*probI + TMath::Sqrt(0.1)*probNS + TMath::Sqrt(0.1)*probFB,2)) > coin )
            MixedData_NS_FB_LS->Fill(l,b,1);
        if( (TMath::Power(TMath::Sqrt(0.8)*probI + TMath::Sqrt(0.1)*probEW + TMath::Sqrt(0.1)*probFB,2)) > coin )
            MixedData_EW_FB_LS->Fill(l,b,1);
        if( (TMath::Power(TMath::Sqrt(0.7)*probI + TMath::Sqrt(0.1)*probNS + TMath::Sqrt(0.1)*probEW + TMath::Sqrt(0.1)*probFB,2)) > coin )
            FullMixedData_LS->Fill(l,b,1);
 
     
        if( TMath::Power(probI,2) > coin)
            Data_Iso_LS.Fill(l,b,1);
        if( (TMath::Power(0.9*probI + 0.1*probNS,2)) > coin )
            Data_AniNS_LS.Fill(l,b,1);
        if( (TMath::Power(0.9*probI + 0.1*probEW,2)) > coin )
            Data_AniEW_LS.Fill(l,b,1);
        if( (TMath::Power(0.9*probI + 0.1*probFB,2)) > coin )
            Data_AniFB_LS.Fill(l,b,1);
        if( (TMath::Power(0.8*probI + 0.1*probNS + 0.1*probEW,2)) > coin )
            MixedData_NS_EW_LS.Fill(l,b,1);
        if( (TMath::Power(0.8*probI + 0.1*probNS + 0.1*probFB,2)) > coin )
            MixedData_NS_FB_LS.Fill(l,b,1);
        if( (TMath::Power(0.8*probI + 0.1*probEW + 0.1*probFB,2)) > coin )
            MixedData_EW_FB_LS.Fill(l,b,1);
        if( (TMath::Power(0.7*probI + 0.1*probNS + 0.1*probEW + 0.1*probFB,2)) > coin )
            FullMixedData_LS.Fill(l,b,1);

        ++progress;
    }
}

*/

void MC_generate_LS_data(
                         TH2D* Data_Iso_LS,
                         TH2D* Data_AniNS_LS,
                         TH2D* Data_AniEW_LS,
                         TH2D* Data_AniFB_LS,
                         TH2D* MixedData_NS_EW_LS,
                         TH2D* MixedData_NS_FB_LS,
                         TH2D* MixedData_EW_FB_LS,
                         TH2D* FullMixedData_LS,
                         TH1D &Data_hwI,
                         TH1D &Data_hwNS,
                         TH1D &Data_hwEW,
                         TH1D &Data_hwFB,
                         TF2 &dI,
                         TF2 &dNS,
                         TF2 &dEW,
                         TF2 &dFB,
                         std::ofstream &log_file,
                         TRandom3 &r_gen
                         )

{
    Double_t coin = 0;
    Double_t l = 0,b = 0, cosb = 0;
    
    /////////////// Template probabilities...
    
    Double_t probI = 0;
    Double_t probNS = 0;
    Double_t probEW = 0;
    Double_t probFB = 0;
    
    
    std::cout << "\n\nBuilding LS MC data maps...";
    boost::progress_display progress(data_all_sky_LS_events);
    
    for(unsigned long int ev_idx = 0; ev_idx < data_all_sky_LS_events; ++ev_idx)
    {
        
        l = r_gen.Uniform(-180,180);
        cosb = r_gen.Uniform(-1,1);
        b = TMath::RadToDeg()*TMath::ACos(cosb)-90;
        coin = r_gen.Uniform(0,1);
        
        l *= TMath::DegToRad();
        b *= TMath::DegToRad();
        
        probI = dI.Eval(l,b,0);
        probNS = dNS.Eval(l,b,0);
        probEW = dEW.Eval(l,b,0);
        probFB = dFB.Eval(l,b,0);
        
        Data_hwI.Fill(probI);
        Data_hwNS.Fill(probNS);
        Data_hwEW.Fill(probEW);
        Data_hwFB.Fill(probFB);
        
        l *= TMath::RadToDeg();
        b *= TMath::RadToDeg();
    
        
        Data_Iso_LS->Fill(l,r_gen.Poisson(probI));
        Data_AniNS_LS->Fill(l,b,r_gen.Poisson(TMath::Sqrt(0.9)*probI + TMath::Sqrt(0.1)*probNS));
        Data_AniEW_LS->Fill(l,b,r_gen.Poisson(TMath::Sqrt(0.9)*probI + TMath::Sqrt(0.1)*probEW));
        Data_AniFB_LS->Fill(l,b,r_gen.Poisson(TMath::Sqrt(0.9)*probI + TMath::Sqrt(0.1)*probFB));
        MixedData_NS_EW_LS->Fill(l,b,r_gen.Poisson(TMath::Sqrt(0.8)*probI + TMath::Sqrt(0.1)*probNS + TMath::Sqrt(0.1)*probEW));
        MixedData_NS_FB_LS->Fill(l,b,r_gen.Poisson(TMath::Sqrt(0.8)*probI + TMath::Sqrt(0.1)*probNS + TMath::Sqrt(0.1)*probFB));
        MixedData_EW_FB_LS->Fill(l,b,r_gen.Poisson(TMath::Sqrt(0.8)*probI + TMath::Sqrt(0.1)*probEW + TMath::Sqrt(0.1)*probFB));
        FullMixedData_LS->Fill(l,b,r_gen.Poisson(TMath::Sqrt(0.7)*probI + TMath::Sqrt(0.1)*probNS + TMath::Sqrt(0.1)*probEW + TMath::Sqrt(0.1)*probFB));
    
        ++progress;
    
    }
       
}

void generate_HS_data(
                        Double_t NS_anisotropy,
                        Double_t EW_anisotropy,
                        Double_t FB_anisotropy,
                        TH2D* Data_Iso_HS,
                        TH2D* Data_AniNS_HS,
                        TH2D* Data_AniEW_HS,
                        TH2D* Data_AniFB_HS,
                        TH2D* MixedData_NS_EW_HS,
                        TH2D* MixedData_NS_FB_HS,
                        TH2D* MixedData_EW_FB_HS,
                        TH2D* FullMixedData_HS,
                        TH1D &Data_hwI,
                        TH1D &Data_hwNS,
                        TH1D &Data_hwEW,
                        TH1D &Data_hwFB,
                        TF2 &dI,
                        TF2 &dNS,
                        TF2 &dEW,
                        TF2 &dFB,
                        TH2D &Template_Iso_HS,
                        TH2D &Template_AniNS_HS,
                        TH2D &Template_AniEW_HS,
                        TH2D &Template_AniFB_HS,
                        std::ofstream &log_file,
                        TRandom3 &r_gen
                      )

{

    std::cout << "\n\n Building HS Data maps\n\n";
    
    ///////////// Pointers to template histos
    
    TH2D* pTemplate_Data_Iso_HS = &Template_Iso_HS;
    TH2D* Template_Data_AniNS_HS = (TH2D*)Template_Iso_HS.Clone("Template_Data_AniNS_HS");
    TH2D* Template_Data_AniEW_HS = (TH2D*)Template_Iso_HS.Clone("Template_Data_AniEW_HS");
    TH2D* Template_Data_AniFB_HS = (TH2D*)Template_Iso_HS.Clone("Template_Data_AniFB_HS");
    TH2D* Template_MixedData_NS_EW_HS = (TH2D*)Template_Iso_HS.Clone("Template_MixedData_NS_EW_HS");
    TH2D* Template_MixedData_NS_FB_HS = (TH2D*)Template_Iso_HS.Clone("Template_MixedData_NS_FB_HS");
    TH2D* Template_MixedData_EW_FB_HS = (TH2D*)Template_Iso_HS.Clone("Template_MixedData_EW_FB_HS");
    TH2D* Template_FullMixedData_HS = (TH2D*)Template_Iso_HS.Clone("Template_FullMixedData_HS");
    
    Template_Data_AniNS_HS->Reset();
    Template_Data_AniEW_HS->Reset();
    Template_Data_AniFB_HS->Reset();
    Template_MixedData_NS_EW_HS->Reset();
    Template_MixedData_NS_FB_HS->Reset();
    Template_FullMixedData_HS->Reset();
    
    for(Int_t bX = 1; bX <= Template_Data_AniNS_HS->GetNbinsX(); ++bX)
    {
        for(Int_t bY = 1; bY <= Template_Data_AniNS_HS->GetNbinsY(); ++bY)
        {
            Template_Data_AniNS_HS->SetBinContent(bX,bY, (1 - NS_anisotropy)*Template_Iso_HS.GetBinContent(bX,bY) + NS_anisotropy*Template_AniNS_HS.GetBinContent(bX,bY) );
            Template_Data_AniEW_HS->SetBinContent(bX,bY, (1 - EW_anisotropy)*Template_Iso_HS.GetBinContent(bX,bY) + EW_anisotropy*Template_AniEW_HS.GetBinContent(bX,bY) );
            Template_Data_AniFB_HS->SetBinContent(bX,bY, (1 - FB_anisotropy)*Template_Iso_HS.GetBinContent(bX,bY) + FB_anisotropy*Template_AniFB_HS.GetBinContent(bX,bY) );
            Template_MixedData_NS_EW_HS->SetBinContent(bX,bY, (1 - NS_anisotropy - EW_anisotropy)*Template_Iso_HS.GetBinContent(bX,bY) + NS_anisotropy*Template_AniNS_HS.GetBinContent(bX,bY) + EW_anisotropy*Template_AniEW_HS.GetBinContent(bX,bY) );
            Template_MixedData_NS_FB_HS->SetBinContent(bX,bY, (1 - NS_anisotropy - FB_anisotropy)*Template_Iso_HS.GetBinContent(bX,bY) + NS_anisotropy*Template_AniNS_HS.GetBinContent(bX,bY) + FB_anisotropy*Template_AniFB_HS.GetBinContent(bX,bY) );
            Template_MixedData_EW_FB_HS->SetBinContent(bX,bY, (1 - EW_anisotropy - FB_anisotropy)*Template_Iso_HS.GetBinContent(bX,bY) + EW_anisotropy*Template_AniEW_HS.GetBinContent(bX,bY) + FB_anisotropy*Template_AniFB_HS.GetBinContent(bX,bY) );
            Template_FullMixedData_HS->SetBinContent(bX,bY, (1 - NS_anisotropy - EW_anisotropy - FB_anisotropy)*Template_Iso_HS.GetBinContent(bX,bY) + NS_anisotropy*Template_AniNS_HS.GetBinContent(bX,bY) + EW_anisotropy*Template_AniEW_HS.GetBinContent(bX,bY) + FB_anisotropy*Template_AniFB_HS.GetBinContent(bX,bY) );
        }
    }
    
    
    ///////////// Filling maps...
    
    Data_Iso_HS->FillRandom(pTemplate_Data_Iso_HS,data_all_sky_HS_events);
    Data_AniNS_HS->FillRandom(Template_Data_AniNS_HS,data_all_sky_HS_events);
    Data_AniEW_HS->FillRandom(Template_Data_AniEW_HS,data_all_sky_HS_events);
    Data_AniFB_HS->FillRandom(Template_Data_AniFB_HS,data_all_sky_HS_events);
    MixedData_NS_EW_HS->FillRandom(Template_MixedData_NS_EW_HS,data_all_sky_HS_events);
    MixedData_NS_FB_HS->FillRandom(Template_MixedData_NS_FB_HS,data_all_sky_HS_events);
    MixedData_EW_FB_HS->FillRandom(Template_MixedData_EW_FB_HS,data_all_sky_HS_events);
    FullMixedData_HS->FillRandom(Template_FullMixedData_HS,data_all_sky_HS_events);

     
}


/*

void generate_HS_data(
                        TH2D* Data_Iso_HS,
                        TH2D* Data_AniNS_HS,
                        TH2D* Data_AniEW_HS,
                        TH2D* Data_AniFB_HS,
                        TH2D* MixedData_NS_EW_HS,
                        TH2D* MixedData_NS_FB_HS,
                        TH2D* MixedData_EW_FB_HS,
                        TH2D* FullMixedData_HS,
                        TH1D &Data_hwI,
                        TH1D &Data_hwNS,
                        TH1D &Data_hwEW,
                        TH1D &Data_hwFB,
                        TF2 &dI,
                        TF2 &dNS,
                        TF2 &dEW,
                        TF2 &dFB,
                        TH2D &Template_Iso_HS,
                        TH2D &Template_AniNS_HS,
                        TH2D &Template_AniEW_HS,
                        TH2D &Template_AniFB_HS,
                        std::ofstream &log_file,
                        TRandom3 &r_gen
                      )

{
 
    Double_t coin = 0;
    Double_t l = 0,b = 0, cosb = 0;
 
    /////////////// Template probabilities...
 
    Double_t probI = 0;
    Double_t probNS = 0;
    Double_t probEW = 0;
    Double_t probFB = 0;
 
 
    std::cout << "\n\nBuilding HS data maps...";
    boost::progress_display progress(data_all_sky_HS_events);
 
    for(unsigned long int ev_idx = 0; ev_idx < data_all_sky_HS_events; ++ev_idx)
    {
 
        l = r_gen.Uniform(-180,180);
        cosb = r_gen.Uniform(-1,1);
        b = TMath::RadToDeg()*TMath::ACos(cosb)-90;
        coin = r_gen.Uniform(0,1);
 
        l *= TMath::DegToRad();
        b *= TMath::DegToRad();
 
        probI = dI.Eval(l,b,0);
        probNS = dNS.Eval(l,b,0);
        probEW = dEW.Eval(l,b,0);
        probFB = dFB.Eval(l,b,0);
 
        Data_hwI.Fill(probI);
        Data_hwNS.Fill(probNS);
        Data_hwEW.Fill(probEW);
        Data_hwFB.Fill(probFB);
 
        l *= TMath::RadToDeg();
        b *= TMath::RadToDeg();
 
 
        /////////////// Filling histos (Hit & Miss method)
 
        if( TMath::Power(probI,2) > coin)
            Data_Iso_HS.Fill(l,b,1);
        if( (TMath::Power(TMath::Sqrt(0.9)*probI + TMath::Sqrt(0.1)*probNS,2)) > coin )
            Data_AniNS_HS.Fill(l,b,1);
        if( (TMath::Power(TMath::Sqrt(0.9)*probI + TMath::Sqrt(0.1)*probEW,2)) > coin )
            Data_AniEW_HS.Fill(l,b,1);
        if( (TMath::Power(TMath::Sqrt(0.9)*probI + TMath::Sqrt(0.1)*probFB,2)) > coin )
            Data_AniFB_HS.Fill(l,b,1);
        if( (TMath::Power(TMath::Sqrt(0.8)*probI + TMath::Sqrt(0.1)*probNS + TMath::Sqrt(0.1)*probEW,2)) > coin )
            MixedData_NS_EW_HS.Fill(l,b,1);
        if( (TMath::Power(TMath::Sqrt(0.8)*probI + TMath::Sqrt(0.1)*probNS + TMath::Sqrt(0.1)*probFB,2)) > coin )
            MixedData_NS_FB_HS.Fill(l,b,1);
        if( (TMath::Power(TMath::Sqrt(0.8)*probI + TMath::Sqrt(0.1)*probEW + TMath::Sqrt(0.1)*probFB,2)) > coin )
            MixedData_EW_FB_HS.Fill(l,b,1);
        if( (TMath::Power(TMath::Sqrt(0.7)*probI + TMath::Sqrt(0.1)*probNS + TMath::Sqrt(0.1)*probEW + TMath::Sqrt(0.1)*probFB,2)) > coin )
            FullMixedData_HS.Fill(l,b,1);
 
 
        if( TMath::Power(probI,2) > coin)
            Data_Iso_HS.Fill(l,b,1);
        if( (TMath::Power(0.9*probI + 0.1*probNS,2)) > coin )
            Data_AniNS_HS.Fill(l,b,1);
        if( (TMath::Power(0.9*probI + 0.1*probEW,2)) > coin )
            Data_AniEW_HS.Fill(l,b,1);
        if( (TMath::Power(0.9*probI + 0.1*probFB,2)) > coin )
            Data_AniFB_HS.Fill(l,b,1);
        if( (TMath::Power(0.8*probI + 0.1*probNS + 0.1*probEW,2)) > coin )
            MixedData_NS_EW_HS.Fill(l,b,1);
        if( (TMath::Power(0.8*probI + 0.1*probNS + 0.1*probFB,2)) > coin )
            MixedData_NS_FB_HS.Fill(l,b,1);
        if( (TMath::Power(0.8*probI + 0.1*probEW + 0.1*probFB,2)) > coin )
            MixedData_EW_FB_HS.Fill(l,b,1);
        if( (TMath::Power(0.7*probI + 0.1*probNS + 0.1*probEW + 0.1*probFB,2)) > coin )
            FullMixedData_HS.Fill(l,b,1);
 
        ++progress;
 
    }
}
 
*/


void MC_generate_HS_data(
                            TH2D* Data_Iso_HS,
                            TH2D* Data_AniNS_HS,
                            TH2D* Data_AniEW_HS,
                            TH2D* Data_AniFB_HS,
                            TH2D* MixedData_NS_EW_HS,
                            TH2D* MixedData_NS_FB_HS,
                            TH2D* MixedData_EW_FB_HS,
                            TH2D* FullMixedData_HS,
                            TH1D &Data_hwI,
                            TH1D &Data_hwNS,
                            TH1D &Data_hwEW,
                            TH1D &Data_hwFB,
                            TF2 &dI,
                            TF2 &dNS,
                            TF2 &dEW,
                            TF2 &dFB,
                            std::ofstream &log_file,
                            TRandom3 &r_gen
                        )
    
{
    
    Double_t coin = 0;
    Double_t l = 0,b = 0, cosb = 0;
    
    /////////////// Template probabilities...
        
    Double_t probI = 0;
    Double_t probNS = 0;
    Double_t probEW = 0;
    Double_t probFB = 0;
        
        
    std::cout << "\n\nBuilding HS MC data maps...";
    boost::progress_display progress(data_all_sky_HS_events);
        
    for(unsigned long int ev_idx = 0; ev_idx < data_all_sky_HS_events; ++ev_idx)
    {
            
        l = r_gen.Uniform(-180,180);
        cosb = r_gen.Uniform(-1,1);
        b = TMath::RadToDeg()*TMath::ACos(cosb)-90;
        coin = r_gen.Uniform(0,1);
        
        l *= TMath::DegToRad();
        b *= TMath::DegToRad();
            
        probI = dI.Eval(l,b,0);
        probNS = dNS.Eval(l,b,0);
        probEW = dEW.Eval(l,b,0);
        probFB = dFB.Eval(l,b,0);
            
        Data_hwI.Fill(probI);
        Data_hwNS.Fill(probNS);
        Data_hwEW.Fill(probEW);
        Data_hwFB.Fill(probFB);
            
        l *= TMath::RadToDeg();
        b *= TMath::RadToDeg();
    
            
        Data_Iso_HS->Fill(l,r_gen.Poisson(probI));
        Data_AniNS_HS->Fill(l,b,r_gen.Poisson(TMath::Sqrt(0.9)*probI + TMath::Sqrt(0.1)*probNS));
        Data_AniEW_HS->Fill(l,b,r_gen.Poisson(TMath::Sqrt(0.9)*probI + TMath::Sqrt(0.1)*probEW));
        Data_AniFB_HS->Fill(l,b,r_gen.Poisson(TMath::Sqrt(0.9)*probI + TMath::Sqrt(0.1)*probFB));
        MixedData_NS_EW_HS->Fill(l,b,r_gen.Poisson(TMath::Sqrt(0.8)*probI + TMath::Sqrt(0.1)*probNS + TMath::Sqrt(0.1)*probEW));
        MixedData_NS_FB_HS->Fill(l,b,r_gen.Poisson(TMath::Sqrt(0.8)*probI + TMath::Sqrt(0.1)*probNS + TMath::Sqrt(0.1)*probFB));
        MixedData_EW_FB_HS->Fill(l,b,r_gen.Poisson(TMath::Sqrt(0.8)*probI + TMath::Sqrt(0.1)*probEW + TMath::Sqrt(0.1)*probFB));
        FullMixedData_HS->Fill(l,b,r_gen.Poisson(TMath::Sqrt(0.7)*probI + TMath::Sqrt(0.1)*probNS + TMath::Sqrt(0.1)*probEW + TMath::Sqrt(0.1)*probFB));
        
        ++progress;
            
    }
    
}
    
    
    
    

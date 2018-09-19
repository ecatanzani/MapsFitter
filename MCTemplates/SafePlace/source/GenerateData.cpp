
#include "MyHead.h"

void generate_LS_data(TH2D* Data_Iso_LS,TH2D* Data_AniNS_LS,TH2D* Data_AniEW_LS,TH2D* Data_AniFB_LS,TH2D* MixedData_NS_EW_LS,TH2D* MixedData_NS_FB_LS,TH2D* MixedData_EW_FB_LS,TH2D* FullMixedData_LS,TH1D &Data_hwNS,TH1D &Data_hwEW,TH1D &Data_hwFB,std::ofstream &log_file,TRandom3 &r_gen) {
    
    Double_t perc = 0;
    Double_t l = 0,b = 0,cosb = 0;
    Double_t w_NS = 0, w_EW = 0, w_FB = 0;
    
    std::cout<<"\n\n";
    
    for(unsigned long int idx_ev = 0; idx_ev < data_all_sky_LS_events; idx_ev++) {
        
        if(((Double_t)idx_ev/data_all_sky_LS_events)>(perc*0.01)) {
            std::cout << "-> Generating LS Data Map: [ " << perc << " % ]" << std::endl;
            log_file << "-> Generating LS Data Map: [ " << perc << " % ]" << std::endl;
            perc++;
        }
        
        /*
        l = r_gen.Uniform(0,360);
        cosb = r_gen.Uniform(-1., 1.);
        b = TMath::RadToDeg()*TMath::ACos(cosb);
        
        w_NS = TMath::Sqrt(3./(4*TMath::Pi()))*TMath::Cos(b*TMath::DegToRad());
        w_EW = TMath::Sqrt(2.)*TMath::Sqrt(3./(4*TMath::Pi()))*TMath::Sin(b*TMath::DegToRad())*TMath::Sin(l*TMath::DegToRad());
        w_FB = - TMath::Sqrt(2)*TMath::Sqrt(3./(4*TMath::Pi()))*TMath::Sin(b*TMath::DegToRad())*TMath::Cos(l*TMath::DegToRad());
        
        l-=180;
        b-=90;
        */
        
        l = r_gen.Uniform(-180,180);
        cosb = r_gen.Uniform(-1,1);
        b = TMath::RadToDeg()*TMath::ACos(cosb)-90;
        
        w_NS = - 0.5*TMath::Sqrt(3/TMath::Pi())*TMath::Sin(b*TMath::DegToRad());
        w_EW = - TMath::Sqrt(3/(2*TMath::Pi()))*TMath::Cos(b*TMath::DegToRad())*TMath::Sin(l*TMath::DegToRad());
        w_FB = TMath::Sqrt(3/(2*TMath::Pi()))*TMath::Cos(b*TMath::DegToRad())*TMath::Cos(l*TMath::DegToRad());
        
        Data_hwNS.Fill(w_NS);
        Data_hwEW.Fill(w_EW);
        Data_hwFB.Fill(w_FB);
        
        Data_Iso_LS->Fill(l,b,1);
        Data_AniNS_LS->Fill(l,b,0.9 + 0.1*w_NS);
        Data_AniEW_LS->Fill(l,b,0.9 + 0.1*w_EW);
        Data_AniFB_LS->Fill(l,b,0.9 + 0.1*w_FB);
        
        MixedData_NS_EW_LS->Fill(l,b,0.8 + 0.1*w_NS + 0.1*w_EW);
        MixedData_NS_FB_LS->Fill(l,b,0.8 + 0.1*w_NS + 0.1*w_FB);
        MixedData_EW_FB_LS->Fill(l,b,0.8 + 0.1*w_EW + 0.1*w_FB);
        
        FullMixedData_LS->Fill(l,b,0.7 + 0.1*w_NS + 0.1*w_EW + 0.1*w_FB);
        
    }
    
}


void generate_HS_data(TH2D* Data_Iso_HS,TH2D* Data_AniNS_HS,TH2D* Data_AniEW_HS,TH2D* Data_AniFB_HS,TH2D* MixedData_NS_EW_HS,TH2D* MixedData_NS_FB_HS,TH2D* MixedData_EW_FB_HS,TH2D* FullMixedData_HS,TH1D &Data_hwNS,TH1D &Data_hwEW,TH1D &Data_hwFB,std::ofstream &log_file,TRandom3 &r_gen) {
    
    Double_t perc = 0;
    Double_t l = 0,b = 0,cosb = 0;
    Double_t w_NS = 0, w_EW = 0, w_FB = 0;
    
    std::cout<<"\n\n";
    
    for(unsigned long int idx_ev = 0; idx_ev < data_all_sky_HS_events; idx_ev++) {
        
        if(((Double_t)idx_ev/data_all_sky_HS_events)>(perc*0.01)) {
            std::cout << "-> Generating HS Data Map: [ " << perc << " % ]" << std::endl;
            log_file << "-> Generating HS Data Map: [ " << perc << " % ]" << std::endl;
            perc++;
        }
        
        /*
         l = r_gen.Uniform(0,360);
         cosb = r_gen.Uniform(-1., 1.);
         b = TMath::RadToDeg()*TMath::ACos(cosb);
         
         w_NS = TMath::Sqrt(3./(4*TMath::Pi()))*TMath::Cos(b*TMath::DegToRad());
         w_EW = TMath::Sqrt(2.)*TMath::Sqrt(3./(4*TMath::Pi()))*TMath::Sin(b*TMath::DegToRad())*TMath::Sin(l*TMath::DegToRad());
         w_FB = - TMath::Sqrt(2)*TMath::Sqrt(3./(4*TMath::Pi()))*TMath::Sin(b*TMath::DegToRad())*TMath::Cos(l*TMath::DegToRad());
         
         l-=180;
         b-=90;
         */
        
        l = r_gen.Uniform(-180,180);
        cosb = r_gen.Uniform(-1,1);
        b = TMath::RadToDeg()*TMath::ACos(cosb)-90;
        
        w_NS = - 0.5*TMath::Sqrt(3/TMath::Pi())*TMath::Sin(b*TMath::DegToRad());
        w_EW = - TMath::Sqrt(3/(2*TMath::Pi()))*TMath::Cos(b*TMath::DegToRad())*TMath::Sin(l*TMath::DegToRad());
        w_FB = TMath::Sqrt(3/(2*TMath::Pi()))*TMath::Cos(b*TMath::DegToRad())*TMath::Cos(l*TMath::DegToRad());
        
        Data_hwNS.Fill(w_NS);
        Data_hwEW.Fill(w_EW);
        Data_hwFB.Fill(w_FB);
        
        Data_Iso_HS->Fill(l,b,1);
        Data_AniNS_HS->Fill(l,b,0.9 + 0.1*w_NS);
        Data_AniEW_HS->Fill(l,b,0.9 + 0.1*w_EW);
        Data_AniFB_HS->Fill(l,b,0.9 + 0.1*w_FB);
        
        MixedData_NS_EW_HS->Fill(l,b,0.8 + 0.1*w_NS + 0.1*w_EW);
        MixedData_NS_FB_HS->Fill(l,b,0.8 + 0.1*w_NS + 0.1*w_FB);
        MixedData_EW_FB_HS->Fill(l,b,0.8 + 0.1*w_EW + 0.1*w_FB);
        
        FullMixedData_HS->Fill(l,b,0.7 + 0.1*w_NS + 0.1*w_EW + 0.1*w_FB);
        
    }
    
}



#include "MyHead.h"

void generate_LS_templates(TH2D* Template_Iso_LS,TH2D* Template_AniNS_LS,TH2D* Template_AniEW_LS,TH2D* Template_AniFB_LS,TH1D &Template_hwNS,TH1D &Template_hwEW,TH1D &Template_hwFB,std::ofstream &log_file,TRandom3 &r_gen) {
    
    Double_t perc = 0;
    Double_t l = 0,b = 0,cosb = 0;
    Double_t w_NS = 0, w_EW = 0, w_FB = 0;
    
    std::cout<<"\n\n";
    
    for(unsigned long int idx_ev = 0; idx_ev < all_sky_LS_events; idx_ev++) {
        
        if(((Double_t)idx_ev/all_sky_LS_events)>(perc*0.01)) {
            std::cout << "-> Generating LS Templates Map: [ " << perc << " % ]" << std::endl;
            log_file << "-> Generating LS Templates Map: [ " << perc << " % ]" << std::endl;
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
        
        Template_hwNS.Fill(w_NS);
        Template_hwEW.Fill(w_EW);
        Template_hwFB.Fill(w_FB);
        
        Template_Iso_LS->Fill(l,b,1);
        Template_AniNS_LS->Fill(l,b,w_NS);
        Template_AniEW_LS->Fill(l,b,w_EW);
        Template_AniFB_LS->Fill(l,b,w_FB);
        
    }
    
}


void generate_HS_templates(TH2D* Template_Iso_HS,TH2D* Template_AniNS_HS,TH2D* Template_AniEW_HS,TH2D* Template_AniFB_HS,TH1D &Template_hwNS,TH1D &Template_hwEW,TH1D &Template_hwFB,std::ofstream &log_file,TRandom3 &r_gen) {
    
    Double_t perc = 0;
    Double_t l = 0,b = 0,cosb = 0;
    Double_t w_NS = 0, w_EW = 0, w_FB = 0;
    
    std::cout<<"\n\n";
    
    for(unsigned long int idx_ev = 0; idx_ev < all_sky_HS_events; idx_ev++) {
        
        if(((Double_t)idx_ev/all_sky_HS_events)>(perc*0.01)) {
            std::cout << "-> Generating HS Templates Map: [ " << perc << " % ]" << std::endl;
            log_file << "-> Generating HS Templates Map: [ " << perc << " % ]" << std::endl;
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
        
        Template_hwNS.Fill(w_NS);
        Template_hwEW.Fill(w_EW);
        Template_hwFB.Fill(w_FB);
        
        Template_Iso_HS->Fill(l,b,1);
        Template_AniNS_HS->Fill(l,b,w_NS);
        Template_AniEW_HS->Fill(l,b,w_EW);
        Template_AniFB_HS->Fill(l,b,w_FB);
        
    }
    
}

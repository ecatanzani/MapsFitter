
#include "MyHead.h"

void gerenate_ani_map(TH2D &IsoMap,TH2D &AniMap_NS,TH2D &AniMap_EW,TH2D &AniMap_FB,TH1D &h_wNS,TH1D &h_wEW,TH1D &h_wFB,Int_t sky_events,std::ofstream &log_file) {
    
    TRandom3 r_gen(random_seed);
    Double_t perc = 0;
    Double_t l = 0,b = 0,cosb = 0;
    Double_t w_NS = 0, w_EW = 0, w_FB = 0;
    
    std::cout<<"\n\n";
    
    for(Int_t idx_ev = 0; idx_ev < sky_events; idx_ev++) {
        
        EvaluatePercentage(perc,sky_events,idx_ev,log_file,false,false,false,false);
        
        l = r_gen.Uniform(-180.,180.);
        cosb = r_gen.Uniform(-1., 1.);
        b = TMath::RadToDeg()*TMath::ACos(cosb) - 90.;
        
        /*
        w_NS = 1. + TMath::Sin(l*TMath::DegToRad())*TMath::Cos(b*TMath::DegToRad());
        w_EW = 1. + TMath::Cos(l*TMath::DegToRad())*TMath::Cos(b*TMath::DegToRad());
        w_FB = 1. + TMath::Cos(b*TMath::DegToRad());
        */
        
        w_NS = 1 + TMath::Sqrt(3./(4*TMath::Pi()))*TMath::Cos((b+90)*TMath::DegToRad());
        w_EW = 1 + TMath::Sqrt(2.)*TMath::Sqrt(3./(4*TMath::Pi()))*TMath::Sin((b+90)*TMath::DegToRad())*TMath::Sin((l+180)*TMath::DegToRad());
        w_FB = 1 - TMath::Sqrt(2)*TMath::Sqrt(3./(4*TMath::Pi()))*TMath::Sin((b+90)*TMath::DegToRad())*TMath::Cos((l+180)*TMath::DegToRad());
        
        
        h_wNS.Fill(w_NS);
        h_wEW.Fill(w_EW);
        h_wFB.Fill(w_FB);
        
        AniMap_NS.Fill(l,b,w_NS*get_bin_events(IsoMap,l,b));
        AniMap_EW.Fill(l,b,w_EW*get_bin_events(IsoMap,l,b));
        AniMap_FB.Fill(l,b,w_FB*get_bin_events(IsoMap,l,b));
        
    }
    
}

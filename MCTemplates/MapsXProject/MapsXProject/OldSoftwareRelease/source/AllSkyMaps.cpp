
#include "MyHead.h"

void generate_allsky_isomap(TH2D &AllSky_IsoMap,std::ofstream &log_file) {
    
    TRandom3 r_gen(random_seed);
    Double_t perc = 0;
    Double_t l = 0,b = 0,cosb = 0;
    
    std::cout<<"\n\n";
    
    for(Int_t idx_ev = 0; idx_ev < all_sky_events; idx_ev++) {
        
        EvaluatePercentage(perc,all_sky_events,idx_ev,log_file,true,true,false,false);
        
        l = r_gen.Uniform(-180.,180.);
        cosb = r_gen.Uniform(-1., 1.);
        b = TMath::RadToDeg()*TMath::ACos(cosb) - 90.;
        
        AllSky_IsoMap.Fill(l,b,1);
        
    }
    
}

void generate_allsky_animap(TH2D &AllSky_AniMap_NS,TH2D &AllSky_AniMap_EW,TH2D &AllSky_AniMap_FB,std::ofstream &log_file) {
    
    TRandom3 r_gen(random_seed);
    Double_t perc = 0;
    Double_t l = 0,b = 0,cosb = 0;
    Double_t w_NS = 0, w_EW = 0, w_FB = 0;
    
    std::cout<<"\n\n";
    
    for(Int_t idx_ev = 0; idx_ev < all_sky_events; idx_ev++) {
        
        EvaluatePercentage(perc,all_sky_events,idx_ev,log_file,true,false,false,false);
        
        l = r_gen.Uniform(-180.,180.);
        cosb = r_gen.Uniform(-1., 1.);
        b = TMath::RadToDeg()*TMath::ACos(cosb) - 90.;
        
        w_NS = 1 + TMath::Sqrt(3./(4*TMath::Pi()))*TMath::Cos((b+90)*TMath::DegToRad());
        w_EW = 1 + TMath::Sqrt(2.)*TMath::Sqrt(3./(4*TMath::Pi()))*TMath::Sin((b+90)*TMath::DegToRad())*TMath::Sin((l+180)*TMath::DegToRad());
        w_FB = 1 - TMath::Sqrt(2)*TMath::Sqrt(3./(4*TMath::Pi()))*TMath::Sin((b+90)*TMath::DegToRad())*TMath::Cos((l+180)*TMath::DegToRad());
        
        AllSky_AniMap_NS.Fill(l,b,w_NS);
        AllSky_AniMap_EW.Fill(l,b,w_EW);
        AllSky_AniMap_FB.Fill(l,b,w_FB);
        
    }
    
    
}

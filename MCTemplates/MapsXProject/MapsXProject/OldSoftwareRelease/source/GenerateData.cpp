
#include "MyHead.h"

void generate_data(TH2D &dataI,TH2D &dataNS,TH2D &dataEW,TH2D &dataFB,std::ofstream &log_file) {
    
    TRandom3 r_gen(random_seed);
    Double_t perc = 0;
    Double_t l = 0,b = 0,cosb = 0;
    Double_t w_NS = 0, w_EW = 0, w_FB = 0;
    
    std::cout<<"\n\n";
    
    for(Int_t idx_ev = 0; idx_ev < data_all_sky_events; idx_ev++) {
        
        EvaluatePercentage(perc,data_all_sky_events,idx_ev,log_file,true,false,false,true);
        
        l = r_gen.Uniform(0,360);
        cosb = r_gen.Uniform(-1., 1.);
        b = TMath::RadToDeg()*TMath::ACos(cosb);
        
        w_NS = 1 + TMath::Sqrt(3./(4*TMath::Pi()))*TMath::Cos(b*TMath::DegToRad());
        w_EW = 1 + TMath::Sqrt(2.)*TMath::Sqrt(3./(4*TMath::Pi()))*TMath::Sin(b*TMath::DegToRad())*TMath::Sin(l*TMath::DegToRad());
        w_FB = 1 - TMath::Sqrt(2)*TMath::Sqrt(3./(4*TMath::Pi()))*TMath::Sin(b*TMath::DegToRad())*TMath::Cos(l*TMath::DegToRad());
        
        dataI.Fill(l,b,1);
        dataNS.Fill(l,b,w_NS);
        dataEW.Fill(l,b,w_EW);
        dataFB.Fill(l,b,w_FB);
        
    }
    
}


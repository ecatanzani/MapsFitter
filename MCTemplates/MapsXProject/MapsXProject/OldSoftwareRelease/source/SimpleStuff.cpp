
#include "MyHead.h"

void normalize_map(TH2D &Map) {
    
    Double_t tot_events=get_tot_events(Map);
    
    for(Int_t bX=1; bX <= Map.GetNbinsX(); bX++)
        for(Int_t bY=1; bY <= Map.GetNbinsY(); bY++)
            Map.SetBinContent(bX, bY, Map.GetBinContent(bX,bY)/tot_events);
    
}

Double_t get_tot_events(TH2D &Map) {
    
    Double_t tot=0;
    
    for(Int_t bX=1; bX <= Map.GetNbinsX(); bX++)
        for(Int_t bY=1; bY <= Map.GetNbinsY(); bY++)
            tot+=Map.GetBinContent(bX,bY);
    
    return tot;
}

Double_t get_bin_events(TH2D &IsoMap,Double_t glon,Double_t glat) {
    
    Int_t Xbin=0,Ybin=0;
    
    TAxis *Yaxis = IsoMap.GetYaxis();
    TAxis *Xaxis = IsoMap.GetXaxis();
    Ybin=Yaxis->FindBin(glat);
    Xbin=Xaxis->FindBin(glon);
    
    return (Double_t)IsoMap.GetBinContent(Xbin,Ybin);
    
}

void EvaluatePercentage(Double_t &percentage,Int_t entries,Int_t idx,std::ofstream &log_file,bool allsky,bool isomap,bool phishing,bool data) {
    if(((Double_t)idx/entries)>(percentage*0.01)) {
        if(!allsky){
            if(!phishing) {
                std::cout << "-> Generating Anisotropic Map: [ " << percentage << " % ]" << std::endl;
                log_file << "-> Generating Anisotropic Map: [ " << percentage << " % ]" << std::endl;
            }
            else {
                std::cout << "-> Generating Phishing Anisotropic Map: [ " << percentage << " % ]" << std::endl;
                log_file << "-> Generating Phishing Anisotropic Map: [ " << percentage << " % ]" << std::endl;
            }
        }
        else {
            if(!data) {
                if(isomap) {
                    std::cout << "-> Generating All Sky Isotropic Map: [ " << percentage << " % ]" << std::endl;
                    log_file << "-> Generating All Sky Isotropic Map: [ " << percentage << " % ]" << std::endl;
                }
                else {
                    std::cout << "-> Generating All Sky Anisotropic Map: [ " << percentage << " % ]" << std::endl;
                    log_file << "-> Generating All Sky Anisotropic Map: [ " << percentage << " % ]" << std::endl;
                }
            }
            else {
                if(isomap) {
                    std::cout << "-> Generating All Sky Isotropic Data Map: [ " << percentage << " % ]" << std::endl;
                    log_file << "-> Generating All Sky Isotropic Data Map: [ " << percentage << " % ]" << std::endl;
                }
                else {
                    std::cout << "-> Generating All Sky Anisotropic Data Map: [ " << percentage << " % ]" << std::endl;
                    log_file << "-> Generating All Sky Anisotropic Data Map: [ " << percentage << " % ]" << std::endl;
                }
            }
        }
        
         percentage++;
    }
}

void subtract_by_bin(TH2D &NAllSky_AniMap_NS,TH2D &NAllSky_IsoMap) {
    
    for(Int_t bX=1; bX<=NAllSky_IsoMap.GetNbinsX(); bX++)
        for(Int_t bY=1; bY<=NAllSky_IsoMap.GetNbinsY(); bY++)
            NAllSky_AniMap_NS.SetBinContent(bX, bY,(Double_t)NAllSky_AniMap_NS.GetBinContent(bX,bY)-(Double_t)NAllSky_IsoMap.GetBinContent(bX,bY));

}

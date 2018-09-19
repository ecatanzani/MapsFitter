
#include "MyHead.h"

void getPull(TH1* Data,TH1* Templates[],Double_t res[],TH1D &hPull) {
    
    Double_t prob[4] = {0,0,0,0};
    Double_t bW = 0;
    Double_t nevs = 0,nexp = 0,s_exp = 0;
    
    for(Int_t bX = 1; bX <= Data->GetNbinsX(); bX++) {
    
        for(Int_t idx_comp = 0; idx_comp < 4; idx_comp++) {
            if(bX>1)
                prob[idx_comp] = 0;
            
            prob[idx_comp] = Templates[idx_comp]->GetBinContent(bX);
        
        }
        
        bW = Data->GetXaxis()->GetBinWidth(bX);
        nevs = Data->GetBinContent(bX);
        
        if(bX>1)
            nexp = 0;
        
        for(Int_t idx_comp = 0; idx_comp < 4; idx_comp++)
            nexp += res[idx_comp]*prob[idx_comp]*bW;
        
        s_exp = nexp;
        if(s_exp>0)
            //hPull.Fill( TMath::Power((nevs - nexp),2)/s_exp ); that's chi2 !
            hPull.Fill( (nevs - nexp)/s_exp );
    }
    
}

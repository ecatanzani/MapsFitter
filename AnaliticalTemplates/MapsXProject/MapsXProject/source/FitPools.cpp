
#include "MyHead.h"

void getPull(TH1* Data,TH1* Templates[],Double_t res[],TH1D &hPull,bool relative_fit,bool ex_fit) {
    
    Int_t n_pos = 4;
    if(relative_fit)
    n_pos = 3;
    if(ex_fit)
    n_pos = 1;
    
    Double_t prob[n_pos];
    Double_t bW = 0;
    Double_t nevs = 0,nexp = 0,s_exp = 0;
    
    for(Int_t idx=0; idx < n_pos; ++idx)
    prob[idx] = 0;
    
    
    for(Int_t bX = 1; bX <= Data->GetNbinsX(); ++bX) {
        
        for(Int_t idx_comp = 0; idx_comp < n_pos; ++idx_comp) {
            if(bX>1)
            prob[idx_comp] = 0;
            
            prob[idx_comp] = Templates[idx_comp]->GetBinContent(bX);
            
        }
        
        bW = Data->GetXaxis()->GetBinWidth(bX);
        nevs = Data->GetBinContent(bX);
        
        if(bX>1)
        nexp = 0;
        
        for(Int_t idx_comp = 0; idx_comp < n_pos; ++idx_comp)
        nexp += res[idx_comp]*prob[idx_comp]*bW;
        
        s_exp = TMath::Sqrt(nexp);
        if(s_exp>0)
        //hPull.Fill( TMath::Power((nevs - nexp),2)/s_exp ); that's Pearson's chi2 !
        hPull.Fill( (nevs - nexp)/s_exp );
    }
    
}

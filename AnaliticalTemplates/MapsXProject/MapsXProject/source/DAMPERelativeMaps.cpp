
#include "MyHead.h"

void get_relative_histo(TH2D &relative_DAMPE_histo,TH2D* data_DAMPE_histo,TH2D* reference_DAMPE_histo)
{
    TH2D* tmp_relative_DAMPE_histo = (TH2D*) data_DAMPE_histo->Clone("tmp_relative_DAMPE_histo");
    
    for(Int_t bX = 1; bX <= tmp_relative_DAMPE_histo->GetNbinsX(); ++bX)
        for(Int_t bY = 1; bY <= tmp_relative_DAMPE_histo->GetNbinsY(); ++bY)
            tmp_relative_DAMPE_histo->SetBinContent(bX,bY,(tmp_relative_DAMPE_histo->GetBinContent(bX,bY)-reference_DAMPE_histo->GetBinContent(bX,bY))/reference_DAMPE_histo->GetBinContent(bX,bY));
    
    new (&relative_DAMPE_histo) TH2D (*(TH2D*)tmp_relative_DAMPE_histo->Clone("relative_DAMPE_histo"));
    
}

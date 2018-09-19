
#include "MyHead.h"

void scale_reference_map(TH2D &DAMPE_ReferenceMap,TH2D* histo,bool LS)
{
    
    TH1D* pDAMPE_ReferenceMap = (TH1D*)DAMPE_ReferenceMap.Clone("pDAMPE_ReferenceMap");
    
    for(Int_t bX = 1; bX <= histo->GetNbinsX(); ++bX)
    {
        for(Int_t bY = 1; bY <= histo->GetNbinsY(); ++bY)
        {
            if(LS)
                histo->FillRandom(pDAMPE_ReferenceMap,data_all_sky_LS_events);
            else
                histo->FillRandom(pDAMPE_ReferenceMap,data_all_sky_HS_events);
        }
    }
        
}

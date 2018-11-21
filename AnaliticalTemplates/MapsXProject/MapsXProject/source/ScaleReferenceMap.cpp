
#include "MyHead.h"

void scale_reference_map(
                            TH2D* DAMPE_ReferenceMap,
                            TH2D &DAMPE_ReferenceMap_scaled,
                            const ULong64_t n_events,
                            bool low_stat
                         )
{
    
    DAMPE_ReferenceMap->Scale(n_events/(Double_t)DAMPE_ReferenceMap->Integral());
    
    if(low_stat)
        new (&DAMPE_ReferenceMap_scaled) (TH2D) (*(TH2D*)DAMPE_ReferenceMap->Clone("DAMPE_ReferenceMap_LS"));
    else
        new (&DAMPE_ReferenceMap_scaled) (TH2D) (*(TH2D*)DAMPE_ReferenceMap->Clone("DAMPE_ReferenceMap_HS"));
    
}

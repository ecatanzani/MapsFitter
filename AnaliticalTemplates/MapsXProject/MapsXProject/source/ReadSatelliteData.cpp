
#include "MyHead.h"

void read_DAMPE_FullIso(
                            TH2D &DAMPE_ReferenceMap_LS,
                            TH2D &DAMPE_ReferenceMap_HS,
                            ULong64_t data_LS_events,
                            ULong64_t data_HS_events,
                            std::string DAMPE_Iso_Map
                        )
{
    
    TFile inputFile(DAMPE_Iso_Map.c_str(),"READ");
    if(inputFile.IsZombie())
    {
        std::cerr << "\n\n Error loading DAMPE Isotropic whole map: " << DAMPE_Iso_Map << std::endl;
        exit(100);
    }
    
    TH2D* tmp_reference = (TH2D*)inputFile.Get("Iso_SkyMap_1000");
    tmp_reference->Rebin2D(10,10);
    
    scale_reference_map(
                            tmp_reference,
                            DAMPE_ReferenceMap_LS,
                            data_LS_events,
                            true
                        );
    
    scale_reference_map(
                            tmp_reference,
                            DAMPE_ReferenceMap_HS,
                            data_HS_events,
                            false
                        );
    
    inputFile.Close();
    
}

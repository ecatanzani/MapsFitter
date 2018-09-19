
#include "MyHead.h"

void read_DAMPE_FullIso(TH2D &DAMPE_FullIso) {
    
    TFile inputFile(DAMPE_Iso_Map.c_str(),"READ");
    if(inputFile.IsZombie()) {
        std::cerr << "\n\n Error loading DAMPE Isotropic whole map: " << DAMPE_Iso_Map << std::endl;
        exit(100);
    }
    
    TH2D* tmp_reference = (TH2D*)inputFile.Get("Iso_SkyMap_1000");
    tmp_reference->Rebin2D(10,10);
    
    new (&DAMPE_FullIso) (TH2D)(*(TH2D*)tmp_reference->Clone("DAMPE_ReferenceMap"));
    
    inputFile.Close();
    
}

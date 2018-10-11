
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


void get_scaled_isotropic_DAMPE_maps(TH2D &DAMPE_ReferenceMap_LS,TH2D &DAMPE_ReferenceMap_HS)
{
    TFile inputFile(DAMPE_Iso_scaled_Maps.c_str(),"READ");
    if(inputFile.IsZombie())
    {
        std::cerr << "\n\n Error loading scaled DAMPE Isotropic whole map: " << DAMPE_Iso_scaled_Maps << std::endl;
        exit(100);
    }
    
    /*
    TH2D* tmp_scaled_reference_LS = (TH2D*)inputFile.Get("DAMPE_ReferenceMap_LS")->Clone("tmpDAMPE_ReferenceMap_LS");
    TH2D* tmp_scaled_reference_HS = (TH2D*)inputFile.Get("DAMPE_ReferenceMap_HS")->Clone("tmpDAMPE_ReferenceMap_HS");
    
    new (&DAMPE_ReferenceMap_LS) (TH2D)(*(TH2D*)tmp_scaled_reference_LS->Clone("DAMPE_ReferenceMap_LS"));
    new (&DAMPE_ReferenceMap_HS) (TH2D)(*(TH2D*)tmp_scaled_reference_HS->Clone("DAMPE_ReferenceMap_HS"));
    */
    
    new (&DAMPE_ReferenceMap_LS) (TH2D)(*(TH2D*)inputFile.Get("DAMPE_ReferenceMap_LS"));
    new (&DAMPE_ReferenceMap_HS) (TH2D)(*(TH2D*)inputFile.Get("DAMPE_ReferenceMap_HS"));
    
    inputFile.Close();
}

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string>

#include "TFile.h"
#include "TH2D.h"
#include "TMath.h"

void read_DAMPE_FullIso(TH2D &DAMPE_FullIso,std::string input_path);
void scale_reference_map(TH2D &DAMPE_ReferenceMap,TH2D &DAMPE_ReferenceMap_scaled,const unsigned long int n_events,bool low_stat);

int main(int argc,char* argv[])
{
    
    ///////////////// Maps statistics variables:
    
    const static unsigned long int data_all_sky_LS_events = 10e+3;
    const static unsigned long int data_all_sky_HS_events = 10e+4;
    
    ///////////////// Dependency paths:
    
    std::string full_reference_path(argv[1]);
    std::string scaled_reference_outpath(argv[2]);
    
    ///////////////// Histos declaration:
    
    static TH2D DAMPE_ReferenceMap;
    TH2D DAMPE_ReferenceMap_LS;
    TH2D DAMPE_ReferenceMap_HS;
    
    ///////////////////////////////////////////////////////////////
    
    read_DAMPE_FullIso(DAMPE_ReferenceMap,full_reference_path);
        
    scale_reference_map(DAMPE_ReferenceMap,DAMPE_ReferenceMap_LS,data_all_sky_LS_events,true);
    scale_reference_map(DAMPE_ReferenceMap,DAMPE_ReferenceMap_HS,data_all_sky_HS_events,false);
    
    ///////////////// Setting the correct titles to the histos
    
    DAMPE_ReferenceMap_LS.SetTitle("DAMPE reference Map Scaled (LS = 16e+3)");
    DAMPE_ReferenceMap_HS.SetTitle("DAMPE reference Map Scaled (HS = 32e+6)");
    
    ///////////////// Writing final results
    
    TFile scaledResults(scaled_reference_outpath.c_str(),"RECREATE");
    if(scaledResults.IsZombie())
    {
        std::cout << "\n\n Error writing output Isotropic Sky Map \n\n";
        exit(100);
    }
    
    DAMPE_ReferenceMap_LS.Write();
    DAMPE_ReferenceMap_HS.Write();
    DAMPE_ReferenceMap.Write();
    
    scaledResults.Write();
    scaledResults.Close();
    
    return 1;
    
}

void read_DAMPE_FullIso(TH2D &DAMPE_FullIso,std::string input_path) {
    
    TFile inputFile(input_path.c_str(),"READ");
    if(inputFile.IsZombie())
    {
        std::cerr << "\n\n Error loading DAMPE Isotropic whole map: " << input_path << std::endl;
        exit(100);
    }
    else
    {
        TH2D* tmp_reference = (TH2D*)inputFile.Get("Iso_SkyMap_1000");
        tmp_reference->Rebin2D(10,10);
        
        new (&DAMPE_FullIso) (TH2D)(*(TH2D*)tmp_reference->Clone("DAMPE_ReferenceMap"));
    }
    
    inputFile.Close();
    
}

void scale_reference_map(TH2D &DAMPE_ReferenceMap,TH2D &DAMPE_ReferenceMap_scaled,const unsigned long int n_events,bool low_stat)
{
    
    TH2D* pDAMPE_ReferenceMap = (TH2D*)DAMPE_ReferenceMap.Clone("pDAMPE_ReferenceMap");
    
    pDAMPE_ReferenceMap->Scale(n_events/(Double_t)pDAMPE_ReferenceMap->Integral());
    
    if(low_stat)
        new (&DAMPE_ReferenceMap_scaled) (TH2D) (*(TH2D*)pDAMPE_ReferenceMap->Clone("DAMPE_ReferenceMap_LS"));
    else
        new (&DAMPE_ReferenceMap_scaled) (TH2D) (*(TH2D*)pDAMPE_ReferenceMap->Clone("DAMPE_ReferenceMap_HS"));
     
}



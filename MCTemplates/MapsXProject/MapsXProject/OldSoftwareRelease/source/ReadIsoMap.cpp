
#include "MyHead.h"

void obtain_IsoMap(TH2D &IsoMap,std::ofstream &log_file) {
    
    /////////////////////////////// Opening DAMPE acceptance 2D-histo
    TFile infile(IsoMap_final_plot.c_str());
    if(infile.IsZombie()) {
        std::cout << "\n\nError opening DAMPE Events Distribution TFile. Prorgram finished \n\n";
        log_file << "\n\nError opening DAMPE Events Distribution TFile. Prorgram finished \n\n";
        exit(-2);
    }
    
    new (&IsoMap) TH2D(*(TH2D*)infile.Get("Iso_SkyMap_1000"));
    
    IsoMap.SetName("IsoMap");
    
    infile.Close();
}

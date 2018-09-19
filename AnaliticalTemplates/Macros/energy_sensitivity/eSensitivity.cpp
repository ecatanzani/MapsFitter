#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string>

#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TFile.h"


Double_t wtsydp( Double_t minene, Double_t maxene, Double_t index)
{
    Double_t dene = maxene-minene;
    return pow( fabs( (pow(maxene,index+1)-pow(minene,index+1)) / ((index+1)*(dene)) ),  1./(index) );
}

void energy_plot(std::string path1,std::string path2,std::string path3, std::string path4,std::string path5, std::string path6)
{
    
    Double_t s_energy[6];
    Double_t e_value[6];
    
    Double_t errXlow[6];
    Double_t errXhigh[6];
    
    Double_t errYlow[6] = {0,0,0,0,0,0};
    Double_t errYhigh[6] = {0,0,0,0,0,0};
    
    Double_t min[6] = {40,55,72,95,125,166};
    Double_t max[6] = {55,72,95,125,166,218};
    
    std::string vPath[6] = {path1,path2,path3,path4,path5,path6};
    
    
    for(Int_t idx_plot=0; idx_plot < 6; ++idx_plot)
    {
        TFile tmpFile(vPath[idx_plot].c_str(),"READ");
        Double_t n_im;
        TGraph *tmpGraph = (TGraph*)tmpFile.Get("sensitivity_95_CL");
        tmpGraph->GetPoint(0,n_im,s_energy[idx_plot]);
        e_value[idx_plot] = wtsydp(min[idx_plot],max[idx_plot],-3);
        errXlow[idx_plot] = e_value[idx_plot] - min[idx_plot];
        errXhigh[idx_plot] = max[idx_plot] - e_value[idx_plot];
    }
    
    TFile outFile("energy_sensitivity.root","RECREATE");
    
    TGraphAsymmErrors outGraph(6,e_value,s_energy,errXlow,errXhigh,errYlow,errYhigh);
    outGraph.Write();
    
    outFile.Close();
    
}

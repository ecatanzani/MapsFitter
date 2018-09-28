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

void energy_plot(
                    std::string path1,
                    std::string path2,
                    std::string path3,
                    std::string path4,
                    std::string path5,
                    std::string path6,
                    std::string path7,
                    std::string path8,
                    std::string path9,
                    std::string path10
                 )
{
    
    Double_t s_energy[10];
    Double_t e_value[10];
    
    Double_t errXlow[10];
    Double_t errXhigh[10];
    
    Double_t errYlow[10] = {0,0,0,0,0,0,0,0,0,0};
    Double_t errYhigh[10] = {0,0,0,0,0,0,0,0,0,0};
    
    Double_t min[10] = {40,55,72,95,125,166,218,288,380,758};
    Double_t max[10] = {55,72,95,125,166,218,288,380,758,1995};
    
    std::string vPath[10] = {path1,path2,path3,path4,path5,path6,path7,path8,path9,path10};
    
    
    for(Int_t idx_plot=0; idx_plot < 10; ++idx_plot)
    {
        TFile tmpFile(vPath[idx_plot].c_str(),"READ");
        Double_t n_im;
        TGraph *tmpGraph = (TGraph*)tmpFile.Get("sensitivity_95_CL");
        tmpGraph->GetPoint(3,n_im,s_energy[idx_plot]);
        e_value[idx_plot] = wtsydp(min[idx_plot],max[idx_plot],-3);
        errXlow[idx_plot] = e_value[idx_plot] - min[idx_plot];
        errXhigh[idx_plot] = max[idx_plot] - e_value[idx_plot];
    }
    
    TFile outFile("energy_sensitivity.root","RECREATE");
    
    TGraphAsymmErrors outGraph(10,e_value,s_energy,errXlow,errXhigh,errYlow,errYhigh);
    outGraph.Write();
    
    outFile.Close();
    
}

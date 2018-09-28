#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>

#include "TGraphAsymmErrors.h"
#include "TFile.h"

void mplots_to_graph(std::string data_file1,std::string data_file2,std::string data_file3,std::string data_file4)
{
    
    std::string paths[4] = {data_file1,data_file2,data_file3,data_file4};
    
    Double_t emin[9] = {42,56,75,100,133,178,237,316,562};
    Double_t emax[9] = {56,75,100,133,178,237,316,562, 2000};
    
    Double_t xmin[9];
    Double_t xmax[9];
    
    Double_t ymin[9] = {0,0,0,0,0,0,0,0,0};
    Double_t ymax[9] = {0,0,0,0,0,0,0,0,0};
    
    Double_t ani_value[9];
    Double_t energy[9];

    for(int i=0; i<9; ++i)
    {
        ani_value[i] = 0;
        energy[i] = 0;
    }
    
    Int_t idxp=0;
    
    std::string tmp_value;

    for(int idxf=0; idxf<4; ++idxf) {
        
        idxp = 0;
        
        std::ifstream input_file(paths[idxf].c_str());
        if(!input_file.is_open()) {
            std::cerr << "\nERROR 100! File not open. The path is:\n" << paths[idxf] << "\n\n";
            exit(100);
        }
        
        while(getline(input_file,tmp_value,','))
        {
           
            energy[idxp] +=  stod(tmp_value);
            getline(input_file,tmp_value,'\n');
            ani_value[idxp] +=  stod(tmp_value);
            ++idxp;
        }
        
    }
    
    for(int i=0; i<9; ++i)
    {
        ani_value[i] = ani_value[i]/4;
        energy[i] = energy[i]/4;
        
        xmin[i] = energy[i] - emin[i];
        xmax[i] = emax[i] - energy[i];
        
    }
    
    
    
    TGraphAsymmErrors tmpGraph(9,energy,ani_value,xmin,xmax,ymin,ymax);
    
    TFile outFile("importedMeanGraph.root","RECREATE");
    tmpGraph.Write();
    outFile.Write();
    outFile.Close();
   
    
    
}

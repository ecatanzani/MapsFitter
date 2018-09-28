#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>

#include "TGraph.h"
#include "TFile.h"

void transform_to_graph(std::string data_file)
{
    

    std::vector<Double_t> xplot;
    std::vector<Double_t> yplot;
    
    std::string tmp_value;
    
    std::ifstream input_file(data_file.c_str());
    if(!input_file.is_open()) {
        std::cerr << "\nERROR 100! File not open. The path is:\n" << data_file << "\n\n";
        exit(100);
    }
   
    while(getline(input_file,tmp_value,','))
    {
        xplot.push_back(stod(tmp_value));
        getline(input_file,tmp_value,'\n');
        yplot.push_back(stod(tmp_value));
    }
    
    
    TGraph tmpGraph(xplot.size(),&xplot[0],&yplot[0]);
    
    TFile outFile("importedGraph.root","RECREATE");
    tmpGraph.Write();
    outFile.Write();
    outFile.Close();
    
    
    
    
}


#include "MyHead.h"

double compute_ani_level(double covMatrix[][4],double parameters[],double par_err[],std::ofstream &log_file)
{
    double level = TMath::Sqrt(TMath::Power(parameters[1],2) + TMath::Power(parameters[2],2) + TMath::Power(parameters[3],2))/parameters[0];
    double delta_derivative[4];
    double var_delta = 0;
    
    delta_derivative[0] = - TMath::Sqrt(TMath::Power(parameters[1],2) + TMath::Power(parameters[2],2) + TMath::Power(parameters[3],2))/TMath::Power(parameters[0],2);
    delta_derivative[1] = parameters[1]/(parameters[0]*TMath::Sqrt(TMath::Power(parameters[1],2) + TMath::Power(parameters[2],2) + TMath::Power(parameters[3],2)));
    delta_derivative[2] = parameters[2]/(parameters[0]*TMath::Sqrt(TMath::Power(parameters[1],2) + TMath::Power(parameters[2],2) + TMath::Power(parameters[3],2)));
    delta_derivative[3] = parameters[3]/(parameters[0]*TMath::Sqrt(TMath::Power(parameters[1],2) + TMath::Power(parameters[2],2) + TMath::Power(parameters[3],2)));
    
    for(Int_t idx1=0; idx1 < 4; ++idx1)
        for(Int_t idx2=0; idx2 < 4; ++idx2)
            var_delta += delta_derivative[idx1]*delta_derivative[idx2]*covMatrix[idx1][idx2];
    
    
    std::cout << "\n\nAnisotropy level: " << level << " +/- " << TMath::Sqrt(var_delta) << std::endl << std::endl;
    log_file << "\n\nAnisotropy level: " << level << " +/- " << TMath::Sqrt(var_delta) << std::endl << std::endl;
    
    return level;
    
}

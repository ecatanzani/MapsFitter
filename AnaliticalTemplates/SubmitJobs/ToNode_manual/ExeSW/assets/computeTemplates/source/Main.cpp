
#include "MyHead.h"

int main(int argc,char* argv[])
{
    
    ////////////////// Dependencies:
    
    std::string log_path = path_location;
    const std::string DAMPE_Iso_Map = argv[1];
    const std::string DAMPE_Iso_scaled_Maps = argv[2];
    
    ////////////////////////////////////////////
    
    TH2D DAMPE_ReferenceMap_LS,DAMPE_ReferenceMap_HS;
    
    log_path += std::to_string((long long)time_stamp);
    log_path += ".log";
    
    std::ofstream output_log_file(log_path);
    
    if(!output_log_file.is_open())
    {
        std::cerr << "\n\nError writign output log file\n\n";
        exit(100);
    }
    
    get_scaled_isotropic_DAMPE_maps(DAMPE_ReferenceMap_LS,DAMPE_ReferenceMap_HS,DAMPE_Iso_scaled_Maps);
    templates_computation(
                            output_log_file,
                            DAMPE_ReferenceMap_LS,
                            DAMPE_ReferenceMap_HS
                          );
    
    return 1;
    
}

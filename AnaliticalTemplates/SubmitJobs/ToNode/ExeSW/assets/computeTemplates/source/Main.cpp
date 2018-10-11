
#include "MyHead.h"

int main(int argc,char* argv[])
{
    
    ////////////////// Dependencies:
    
    const std::string DAMPE_Iso_scaled_Maps (argv[1]);
    const std::string templates_path (argv[2]);
    const std::string DAMPE_templates_path (argv[3]);
    std::string log_path (argv[4]);
    
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
    
    get_scaled_isotropic_DAMPE_maps(
                                        DAMPE_ReferenceMap_LS,
                                        DAMPE_ReferenceMap_HS,
                                        DAMPE_Iso_scaled_Maps
                                    );
    templates_computation(
                            output_log_file,
                            templates_path,
                            DAMPE_templates_path,
                            DAMPE_ReferenceMap_LS,
                            DAMPE_ReferenceMap_HS
                          );
    
    return 0;
    
}

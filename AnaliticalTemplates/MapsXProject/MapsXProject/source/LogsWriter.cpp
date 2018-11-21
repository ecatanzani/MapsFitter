
#include "MyHead.h"

void create_and_initialize_log(
                                std::ofstream &log_file,
                                ULong64_t data_LS_events,
                                ULong64_t data_HS_events,
                                time_t time_stamp
                               )
{
    
    log_file << "********************* Automatic Log File Generator *******************" << std::endl << std::endl;
    
    log_file << "////////////////////////// Simulation Parameters //////////////////////////"<<std::endl << std::endl;
    log_file << "Simulation timestamp: "<< time_stamp << std::endl;
    log_file << "All Sky LS events: " << data_LS_events << std::endl;
    log_file << "All Sky HS events: " << data_HS_events << std::endl;
    
    log_file << "*\n*\n*\n" << std::endl;
}


std::string output_path_creator(
                                    std::string output_log,
                                    std::string output_root,
                                    const Int_t out_choose,
                                    time_t time_stamp,
                                    Double_t NS_dipole,
                                    Double_t EW_dipole,
                                    Double_t FB_dipole,
                                    bool DAMPE
                                )
{
    
    std::string output;
    std::string subfix = "_";
    
    subfix += std::to_string(NS_dipole);
    subfix += "_";
    subfix += std::to_string(EW_dipole);
    subfix += "_";
    subfix += std::to_string(FB_dipole);
    
    switch(out_choose) {

        case 0:
            output = output_log;
            output += std::to_string((long long)time_stamp);
            output += ".txt";
            std::cout << "\nSelected log file: -> \t " << output << std::endl;
            break;
        
        case 1:     //Templates
            output = output_root;
            output += std::to_string((long long)time_stamp);
            if(DAMPE)
                output += "_DAMPE_templates.root";
            else
                output += "_templates.root";
            std::cout << "\nSelected ROOT file: -> \t " << output << std::endl;
            break;
            
        case 2:     //Data
            output = output_root;
            output += std::to_string((long long)time_stamp);
            output += subfix;
            if(DAMPE)
                output += "_DAMPE_data.root";
            else
                output += "_data.root";
            std::cout << "\nSelected ROOT file: -> \t " << output << std::endl;
            break;
            
        case 3:     //Pools
            output = output_root;
            output += std::to_string((long long)time_stamp);
            output += subfix;
            output += "_pools.root";
            std::cout << "\nSelected Pools ROOT file: -> \t " << output << std::endl;
            break;
            
        case 4:     //Projections
            output = output_root;
            output += std::to_string((long long)time_stamp);
            output += subfix;
            output += "_projections.root";
            std::cout << "\nSelected Projection ROOT file: -> \t " << output << std::endl;
            break;
            
        case 5:     //Stacks
            output = output_root;
            output += std::to_string((long long)time_stamp);
            output += subfix;
            output += "_stacks.root";
            std::cout << "\nSelected Stacks ROOT file: -> \t " << output << std::endl;
            break;
                
    }
  
    return output;
}

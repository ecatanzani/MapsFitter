
#include "MyHead.h"

void create_and_initialize_log(std::ofstream &log_file) {
    if(!log_file.is_open()) {
        std::cout<<"\n\nCannot create output log file! Program finished !"<<std::endl;
        exit(-2);
    }
    log_file_init(log_file);
}

std::string output_path_creator(const Int_t out_choose) {
    // out_choose == 0 means we are creating the path for a log file
    // out_choose == 0 means we are creating the path for a ROOT output file

    std::string output;
  
    switch(out_choose) {

    case 0:
            output = output_log;
            output += std::to_string((long long)time_stamp);
            output += ".txt";
            std::cout << "\nWritten log file: -> \t " << output << std::endl;
            break;

        case 1:     //Templates
            output = output_root;
            output += std::to_string((long long)time_stamp);
            output += "_templates.root";
            std::cout << "\nWritten ROOT file: -> \t " << output << std::endl;
            break;
            
        case 2:     //Data
            output = output_root;
            output += std::to_string((long long)time_stamp);
            output += "_data.root";
            std::cout << "\nWritten ROOT file: -> \t " << output << std::endl;
            break;
            
        case 3:     //Pools
            output = output_root;
            output += std::to_string((long long)time_stamp);
            output += "_pools.root";
            std::cout << "\nWritten Pools ROOT file: -> \t " << output << std::endl;
            break;
            
    }
  
    return output;
}

void log_file_init(std::ofstream &out_file) {
    
    out_file << "********************* Automatic Log File Generator *******************" << std::endl << std::endl;
    
    out_file << "////////////////////////// Simulation Parameters //////////////////////////"<<std::endl << std::endl;
    out_file << "Simulation timestamp: "<< time_stamp << std::endl;
    out_file << "Simulation TRandom3 seed: "<< random_seed << std::endl;
    out_file << "All Sky Template LS events: " << all_sky_LS_events << std::endl;
    out_file << "All Sky Data LS events: " << data_all_sky_LS_events << std::endl;
    out_file << "All Sky Template HS events: " << all_sky_HS_events << std::endl;
    out_file << "All Sky Data HS events: " << data_all_sky_HS_events <<std::endl;
    
    out_file << "*\n*\n*\n"<<std::endl;
}

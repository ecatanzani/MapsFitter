
#include "MyHead.h"

void create_and_initialize_log(
                                std::ofstream &log_file,
                                Int_t s_idx,
                                Int_t s_batch,
                                Int_t n_try
                               )
{
    
    if(!log_file.is_open()) {
        std::cout<<"\n\nCannot create output log file! Program finished !"<<std::endl;
        exit(100);
    }
    
    log_file_init(
                    log_file,
                    s_idx,
                    s_batch,
                    n_try
                  );

}

std::string output_path_creator(
                                    const Int_t try_idx,
                                    const Int_t out_choose,
                                    Double_t NS_dipole,
                                    Double_t EW_dipole,
                                    Double_t FB_dipole,
                                    bool DAMPE,
                                    bool test_fit
                                )

{
    // out_choose == 0 means we are creating the path for a log file
    // out_choose == 0 means we are creating the path for a ROOT output file
    
    std::string output;
    std::string str_idx = "_";
    std::string subfix = "_";
    std::string ani_NS,ani_EW,ani_FB;
    std::ostringstream tmp_ani;
    
    str_idx += std::to_string(try_idx);
    str_idx += "_";
    
    tmp_ani << NS_dipole;
    ani_NS = tmp_ani.str();
    tmp_ani.clear();
    tmp_ani.str("");
    tmp_ani << EW_dipole;
    ani_EW = tmp_ani.str();
    tmp_ani.clear();
    tmp_ani.str("");
    tmp_ani << FB_dipole;
    ani_FB = tmp_ani.str();
    
    subfix += ani_NS;
    subfix += "_";
    subfix += ani_EW;
    subfix += "_";
    subfix += ani_FB;
    
    switch(out_choose) {

        case 0:
            output = std::to_string((long long)time_stamp);
            output += str_idx;
            output += ".txt";
            std::cout << "\nSelected log file: -> \t " << output << std::endl;
            break;
    
        case 2:     //Data
            output = std::to_string((long long)time_stamp);
            output += str_idx;
            output += subfix;
            if(DAMPE)
                output += "_DAMPE_data.root";
            else if(test_fit)
                output += "_test_data.root";
            else
                output += "_data.root";
            std::cout << "\nSelected ROOT file: -> \t " << output << std::endl;
            break;
            
        case 3:     //Pools
            output = std::to_string((long long)time_stamp);
            output += str_idx;
            output += subfix;
            if(DAMPE)
                output += "_DAMPE_pools.root";
            else if(test_fit)
                output += "_test_pools.root";
            else
                output += "_pools.root";
            std::cout << "\nSelected Pools ROOT file: -> \t " << output << std::endl;
            break;
            
        case 4:     //Projections
            output = std::to_string((long long)time_stamp);
            output += str_idx;
            output += subfix;
            if(DAMPE)
                output += "_DAMPE_projections.root";
            else if(test_fit)
                output += "_test_projections.root";
            else
                output += "_projections.root";
            std::cout << "\nSelected Projection ROOT file: -> \t " << output << std::endl;
            break;
            
        case 5:     //Stacks
            output = std::to_string((long long)time_stamp);
            output += str_idx;
            output += subfix;
            output += "_stacks.root";
            std::cout << "\nSelected Stacks ROOT file: -> \t " << output << std::endl;
            break;
            
        case 6:     //Tree
            output = std::to_string((long long)time_stamp);
            output += str_idx;
            output += "_tree.root";
            std::cout << "\nSelected TTree ROOT file: -> \t " << output << std::endl;
            break;
            
    }
  
    return output;
}

void log_file_init(
                    std::ofstream &out_file,
                    Int_t s_idx,
                    Int_t s_batch,
                    Int_t n_try
                   )

{
    
    out_file << "********************* Automatic Log File Generator *******************" << std::endl << std::endl;
    
    out_file << "////////////////////////// Simulation Parameters //////////////////////////"<<std::endl << std::endl;
    out_file << "Simulation timestamp: "<< time_stamp << std::endl;
    out_file << "All Sky LS events: " << data_all_sky_LS_events << std::endl;
    out_file << "All Sky HS events: " << data_all_sky_HS_events << std::endl;
    out_file << "Seed starting reading point: " << s_idx << std::endl;
    out_file << "Number of simulation per job: " << s_batch << std::endl;
    out_file << "Number of total simulations: " << n_try << std::endl;
    
    out_file << "*\n*\n*\n"<<std::endl;
}

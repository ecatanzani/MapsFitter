
#include "MyHead.h"

int main(int argc,char *argv[]) {
    
    
    //////////// Extracting simultion parameters
    
    std::stringstream tmp_par;
    Int_t s_idx = 0,s_batch = 0,n_try = 0;
    
    tmp_par << argv[1];
    tmp_par >> s_idx;
    
    tmp_par.clear();
    tmp_par.str("");
    
    tmp_par << argv[2];
    tmp_par >> s_batch;
    
    tmp_par.clear();
    tmp_par.str("");
    
    tmp_par << argv[3];
    tmp_par >> n_try;
    
    /////////////////////////
    
    //////////// Stuff variables
    
    std::string log_path = output_path_creator(s_idx,0);
    std::ofstream output_log_file(log_path);
    
    create_and_initialize_log(
                                output_log_file,
                                s_idx,
                                s_batch,
                                n_try
                              );
    
    //////////////////////////////////////////////////////////////////////////
    
    if(fitDistribution_test)
    {
        std::cout << "\n\n ///////////////////// TEST MODE ///////////////////// ";
        generate_and_fit_test(output_log_file,s_idx,s_batch);
    }
    else
    {
        if(DAMPE_relative_simulation)
            generate_and_fit_relative(output_log_file,s_idx,s_batch);
        else
            generate_and_fit(output_log_file,s_idx,s_batch);
        
    }
        

    
    output_log_file.close();
    
    return 0;
  
}

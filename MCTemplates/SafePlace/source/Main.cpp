
#include "MyHead.h"

int main(int argc,char *argv[]) {
    
    //////////// Stuff variables
    
    std::string log_path = output_path_creator(0);
    std::ofstream output_log_file(log_path);

    
    //////////////////////////////////////////////////////////////////////////
    
    if(argc>1) {
        std::string template_path(argv[1]),data_path(argv[2]);
        read_from_file(template_path,data_path,output_log_file);
    }
    else
        generate_and_fit(output_log_file);
    
    output_log_file.close();
    
    return 0;
  
}

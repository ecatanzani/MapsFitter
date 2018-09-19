
#include "MyHead.h"

int main(int argc,char *argv[]) {
    
    //////////// Stuff variables
    
    std::string log_path = output_path_creator(0);
    std::string template_path,data_path;
    std::ofstream output_log_file(log_path);
    std::stringstream sp_vsr;
    
    create_and_initialize_log(output_log_file);
    
    //////////////////////////////////////////////////////////////////////////
    
    if(argc>1) {
        
        if(argc==3) {
            sp_vsr << argv[1];
            sp_vsr >> template_path;
            sp_vsr.clear();
            sp_vsr << argv[2];
            sp_vsr >> data_path;
        }
        if(argc==2) {
            std::string tmp_path(argv[1]);
            if(tmp_path.find("template")!=std::string::npos) {
                sp_vsr << argv[1];
                sp_vsr >> template_path;
 //               generate_data(output_log_file,data_path);
            }
            else if(tmp_path.find("data")!=std::string::npos) {
                sp_vsr << argv[1];
                sp_vsr >> data_path;
 //               generate_templates(output_log_file,template_path);
            }
        }
        
 //       read_from_file(template_path,data_path,output_log_file);
    
    }
    else
        generate_and_fit(output_log_file);
    
    output_log_file.close();
    
    return 0;
  
}

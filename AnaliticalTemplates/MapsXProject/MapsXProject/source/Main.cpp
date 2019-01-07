
#include "MyHead.h"

int main(int argc,char *argv[]) {
    
    //////////// Stuff variables
    
    ULong64_t data_LS_events = 0,data_HS_events = 0;
    UInt_t ani_values = 0;
    UInt_t NbinX=0,NbinY=0;
    Bool_t write_tmp_histos = true;
    Bool_t all_sky_simulation = false;
    Bool_t DAMPE_simulation = false;
    Bool_t DAMPE_relative_simulation = false;
    std::vector<Double_t> NS_anisotropy;
    std::vector<Double_t> EW_anisotropy;
    std::vector<Double_t> FB_anisotropy;
    std::string output_log,output_root,DAMPE_Iso_Map,seeds_path;
    
    config_reader(
                    data_LS_events,
                    data_HS_events,
                    NbinX,
                    NbinY,
                    write_tmp_histos,
                    all_sky_simulation,
                    DAMPE_simulation,
                    DAMPE_relative_simulation,
                    NS_anisotropy,
                    EW_anisotropy,
                    FB_anisotropy,
                    output_log,
                    output_root,
                    DAMPE_Iso_Map,
                    seeds_path
                  );
    
    ani_values = (UInt_t)NS_anisotropy.size();
    
    const static time_t time_stamp=time(0);
    
    std::string log_path = output_path_creator(output_log,output_root,0,time_stamp);
    std::string template_path,data_path;
    std::ofstream output_log_file(log_path);
    std::stringstream sp_vsr;
    
    create_and_initialize_log(output_log_file,data_LS_events,data_HS_events,time_stamp);
    
    //////////////////////////////////////////////////////////////////////////
    
    if(argc>1)
    {
        
        if(argc==3)
        {
            sp_vsr << argv[1];
            sp_vsr >> template_path;
            sp_vsr.clear();
            sp_vsr << argv[2];
            sp_vsr >> data_path;
        }
        if(argc==2)
        {
            std::string tmp_path(argv[1]);
            if(tmp_path.find("template")!=std::string::npos)
            {
                sp_vsr << argv[1];
                sp_vsr >> template_path;
                generate_data_interface(
                                            output_log_file,
                                            output_log,
                                            output_root,
                                            DAMPE_Iso_Map,
                                            seeds_path,
                                            time_stamp,
                                            NS_anisotropy,
                                            EW_anisotropy,
                                            FB_anisotropy,
                                            data_LS_events,
                                            data_HS_events,
                                            ani_values,
                                            all_sky_simulation,
                                            DAMPE_simulation,
                                            DAMPE_relative_simulation
                                        );
            }
            else if(tmp_path.find("data")!=std::string::npos)
            {
                sp_vsr << argv[1];
                sp_vsr >> data_path;
                generate_templates(
                                    output_log_file,
                                    data_LS_events,
                                    data_HS_events,
                                    output_log,
                                    output_root,
                                    DAMPE_Iso_Map,
                                    time_stamp
                                   );
            }
        }
        
        //read_from_file(template_path,data_path,output_log_file);
    
    }
    else
        generate_and_fit(
                            data_LS_events,
                            data_HS_events,
                            write_tmp_histos,
                            all_sky_simulation,
                            DAMPE_simulation,
                            DAMPE_relative_simulation,
                            NS_anisotropy,
                            EW_anisotropy,
                            FB_anisotropy,
                            output_log,
                            output_root,
                            DAMPE_Iso_Map,
                            seeds_path,
                            ani_values,
                            output_log_file,
                            time_stamp
                         );
    
    output_log_file.close();
    
    return 0;
  
}

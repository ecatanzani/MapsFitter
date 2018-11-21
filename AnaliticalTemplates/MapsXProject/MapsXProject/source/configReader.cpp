
#include "MyHead.h"

void config_reader(
                    ULong64_t &data_LS_events,
                    ULong64_t &data_HS_events,
                    UInt_t &NbinX,
                    UInt_t &NbinY,
                    Bool_t &write_tmp_histos,
                    Bool_t &all_sky_simulation,
                    Bool_t &DAMPE_simulation,
                    Bool_t &DAMPE_relative_simulation,
                    std::vector<Double_t> &NS_anisotropy,
                    std::vector<Double_t> &EW_anisotropy,
                    std::vector<Double_t> &FB_anisotropy,
                    std::string &output_log,
                    std::string &output_root,
                    std::string &DAMPE_Iso_Map,
                    std::string &seeds_path
                   )
{
    std::string config_path = "/Users/enrico/myRepos/MapsFitter/AnaliticalTemplates/simulation.conf";
    std::ifstream input_file(config_path.c_str());
    if(!input_file.is_open()) {
        std::cerr << "\nERROR 100! File not open. The path is:\n" << config_path << "\n\n";
        exit(100);
    }
    std::string input_string((std::istreambuf_iterator< char >(input_file)), (std::istreambuf_iterator< char >()));
    std::string tmp_str;
    input_file.close();
    std::istringstream input_stream(input_string);
    while(input_stream >> tmp_str)
    {
        switch(numerize(tmp_str))
        {
            case 0:
                input_stream >> tmp_str;
                data_LS_events = (UInt_t)std::stol(tmp_str,nullptr,10);
                break;
            
            case 1:
                input_stream >> tmp_str;
                data_HS_events = (UInt_t)std::stol(tmp_str,nullptr,10);
                break;
            
            case 2:
                input_stream >> tmp_str;
                write_tmp_histos = (tmp_str == "1");
                break;
            
            case 3:
                while(tmp_str != "EW_anisotropy")
                {
                    input_stream >> tmp_str;
                    NS_anisotropy.push_back((Double_t)std::stod(tmp_str,nullptr));
                }
            
                NS_anisotropy.resize(NS_anisotropy.size());
                EW_anisotropy.resize(NS_anisotropy.size());
                FB_anisotropy.resize(NS_anisotropy.size());
            
                for(int idx=0; idx <NS_anisotropy.size(); ++idx)
                {
                    input_stream >> tmp_str;
                    EW_anisotropy[idx] = (Double_t)std::stod(tmp_str,nullptr);
                }
            
                input_stream >> tmp_str;
            
                for(int idx=0; idx <NS_anisotropy.size(); ++idx)
                {
                    input_stream >> tmp_str;
                    FB_anisotropy[idx] = (Double_t)std::stod(tmp_str,nullptr);
                }
                break;
            
            case 4:
                input_stream >> tmp_str;
                all_sky_simulation = (tmp_str == "1");
                break;
            
            case 5:
                input_stream >> tmp_str;
                DAMPE_simulation = (tmp_str == "1");
                break;
            
            case 6:
                input_stream >> tmp_str;
                DAMPE_relative_simulation = (tmp_str == "1");
                break;
            
            case 7:
                input_stream >> tmp_str;
                NbinX = (UInt_t)std::stoi(tmp_str,nullptr,10);
                break;
            
            case 8:
                input_stream >> tmp_str;
                NbinY = (UInt_t)std::stoi(tmp_str,nullptr,10);
                break;
            
            case 9:
                input_stream >> tmp_str;
                output_log = tmp_str;
                break;
            
            case 10:
                input_stream >> tmp_str;
                output_root = tmp_str;
                break;
            
            case 11:
                input_stream >> tmp_str;
                DAMPE_Iso_Map = tmp_str;
                break;
            
            case 12:
                input_stream >> tmp_str;
                seeds_path = tmp_str;
                break;
        }
    }
    
}

int numerize(std::string tmp_string)
{
    int string_id = -1;
    
    if(tmp_string == "data_LS_events")
        string_id = 0;
    else if(tmp_string == "data_HS_events")
        string_id = 1;
    else if(tmp_string == "write_tmp_histos")
        string_id = 2;
    else if(tmp_string == "NS_anisotropy")
        string_id = 3;
    else if(tmp_string == "all_sky_simulation")
        string_id = 4;
    else if(tmp_string == "DAMPE_simulation")
        string_id = 5;
    else if(tmp_string == "DAMPE_relative_simulation")
        string_id = 6;
    else if(tmp_string == "NbinX")
        string_id = 7;
    else if(tmp_string == "NbinY")
        string_id = 8;
    else if(tmp_string == "output_log")
        string_id = 9;
    else if(tmp_string == "output_root")
        string_id = 10;
    else if(tmp_string == "DAMPE_Iso_Map")
        string_id = 11;
    else if(tmp_string == "seeds_path")
        string_id = 12;
    
    return string_id;
}


#include "MyHead.h"

void deploy_simulation(
                        std::ofstream &output_log_file,
                        std::ifstream &inSeed,
                        Int_t s_idx,
                        Int_t s_batch,
                        std::string tmp_seed_str,
                        UInt_t tmp_seed,
                        TTree &fiTree,
                        fitResult &tmp_fit
                       )
{
    
    UInt_t seed_line = 0;
    
    if(!inSeed.is_open())
    {
        std::cout << "\n\nCannot open input seeds file! Program finished !" << std::endl;
        output_log_file << "\n\nCannot open input seeds file! Program finished !" << std::endl;
        exit(100);
    }
    
    /////// Ignore used seeds values in previous simulations (if any)
    
    for(Int_t i_idx=0; i_idx<s_idx*s_batch*ani_values; ++i_idx)
        inSeed.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    
    for(Int_t itry = 0; itry < s_batch; ++itry)
    {
        
        std::cout << "\n\n ////////////////////////////// Try " << itry << " ////////////////////////////// \n\n";
        output_log_file << "\n\n ////////////////////////////// Try " << itry << " ////////////////////////////// \n\n";
        
        /*
        /////// Ignore used seeds values in previous simulations in this job (if any)
        
        for(Int_t i_idx=0; i_idx<itry*ani_values; ++i_idx)
            inSeed.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        */
        
        for(Int_t idx_ani = 0; idx_ani < ani_values; ++idx_ani)
        {
            
            std::cout << "\n\nNS Anisotropy: " << NS_anisotropy[idx_ani] << "\n";
            std::cout << "EW Anisotropy: " << EW_anisotropy[idx_ani] << "\n";
            std::cout << "FB Anisotropy: " << FB_anisotropy[idx_ani] << "\n\n";
            
            output_log_file << "\n\nNS Anisotropy: " << NS_anisotropy[idx_ani] << "\n";
            output_log_file << "EW Anisotropy: " << EW_anisotropy[idx_ani] << "\n";
            output_log_file << "FB Anisotropy: " << FB_anisotropy[idx_ani] << "\n\n";
            
            inSeed >> tmp_seed_str;
            tmp_seed = (UInt_t)std::stoul(tmp_seed_str,nullptr,10);
            
            seed_line = (s_idx*s_batch*ani_values - 1) + itry*ani_values + (idx_ani + 1);
            
            if(all_sky_simulation)
                allSky_singleTry_fit(
                                        output_log_file,
                                        fiTree,
                                        tmp_fit,
                                        tmp_seed,
                                        itry,
                                        NS_anisotropy[idx_ani],
                                        EW_anisotropy[idx_ani],
                                        FB_anisotropy[idx_ani],
                                        idx_ani,
                                        s_idx,
                                        seed_line
                                     );
            
            if(DAMPE_simulation)
                DAMPE_singleTry_fit(
                                        output_log_file,
                                        fiTree,
                                        tmp_fit,
                                        tmp_seed,
                                        itry,
                                        NS_anisotropy[idx_ani],
                                        EW_anisotropy[idx_ani],
                                        FB_anisotropy[idx_ani],
                                        idx_ani,
                                        s_idx,
                                        seed_line
                                    );
            
        }
    }
}


void deploy_relative_simulation(
                                    std::ofstream &output_log_file,
                                    std::ifstream &inSeed,
                                    Int_t s_idx,
                                    Int_t s_batch,
                                    std::string tmp_seed_str,
                                    UInt_t tmp_seed,
                                    TTree &fiTree,
                                    relative_fitResult &tmp_fit,
                                    TH2D &DAMPE_ReferenceMap_LS,
                                    TH2D &DAMPE_ReferenceMap_HS
                                )
{
    
    UInt_t seed_line = 0;
    
    if(!inSeed.is_open())
    {
        std::cout << "\n\nCannot open input seeds file! Program finished !" << std::endl;
        output_log_file << "\n\nCannot open input seeds file! Program finished !" << std::endl;
        exit(100);
    }
    
    /////// Ignore used seeds values in previous simulations (if any)
    
    for(Int_t i_idx=0; i_idx<s_idx*s_batch*ani_values; ++i_idx)
        inSeed.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    
    for(Int_t itry = 0; itry < s_batch; ++itry)
    {
        
        std::cout << "\n\n ////////////////////////////// Try " << itry << " ////////////////////////////// \n\n";
        output_log_file << "\n\n ////////////////////////////// Try " << itry << " ////////////////////////////// \n\n";
        
        /*
        /////// Ignore used seeds values in previous simulations in this job (if any)
        
        for(Int_t i_idx=0; i_idx<itry*ani_values; ++i_idx)
            inSeed.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        */
        
        for(Int_t idx_ani = 0; idx_ani < ani_values; ++idx_ani)
        {
            
            std::cout << "\n\nNS Anisotropy: " << NS_anisotropy[idx_ani] << "\n";
            std::cout << "EW Anisotropy: " << EW_anisotropy[idx_ani] << "\n";
            std::cout << "FB Anisotropy: " << FB_anisotropy[idx_ani] << "\n\n";
            
            output_log_file << "\n\nNS Anisotropy: " << NS_anisotropy[idx_ani] << "\n";
            output_log_file << "EW Anisotropy: " << EW_anisotropy[idx_ani] << "\n";
            output_log_file << "FB Anisotropy: " << FB_anisotropy[idx_ani] << "\n\n";
            
            inSeed >> tmp_seed_str;
            tmp_seed = (UInt_t)std::stoul(tmp_seed_str,nullptr,10);
            
            seed_line = (s_idx*s_batch*ani_values - 1) + itry*ani_values + (idx_ani + 1);
           
            DAMPE_relative_singleTry_fit(
                                            output_log_file,
                                            fiTree,
                                            tmp_fit,
                                            tmp_seed,
                                            itry,
                                            NS_anisotropy[idx_ani],
                                            EW_anisotropy[idx_ani],
                                            FB_anisotropy[idx_ani],
                                            idx_ani,
                                            DAMPE_ReferenceMap_LS,
                                            DAMPE_ReferenceMap_HS,
                                            s_idx,
                                            seed_line
                                        );
            
        }
    }
}

void deploy_test_simulation(
                                std::ofstream &output_log_file,
                                std::ifstream &inSeed,
                                Int_t s_idx,
                                Int_t s_batch,
                                std::string tmp_seed_str,
                                UInt_t tmp_seed,
                                TTree &fiTree,
                                test_fitResult &tmp_fit
                            )

{
    
    UInt_t seed_line = 0;
    
    if(!inSeed.is_open())
    {
        std::cout << "\n\nCannot open input seeds file! Program finished !" << std::endl;
        output_log_file << "\n\nCannot open input seeds file! Program finished !" << std::endl;
        exit(100);
    }
    
    /////// Ignore used seeds values in previous simulations (if any)
    
    for(Int_t i_idx=0; i_idx<s_idx*s_batch*ani_values; ++i_idx)
        inSeed.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    
    for(Int_t itry = 0; itry < s_batch; ++itry)
    {
        
        std::cout << "\n\n ////////////////////////////// Try " << itry << " ////////////////////////////// \n\n";
        output_log_file << "\n\n ////////////////////////////// Try " << itry << " ////////////////////////////// \n\n";
        
        /*
        /////// Ignore used seeds values in previous simulations in this job (if any)
        
        for(Int_t i_idx=0; i_idx<itry*ani_values; ++i_idx)
            inSeed.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        */
        
        for(Int_t idx_ani = 0; idx_ani < ani_values; ++idx_ani)
        {
            
            std::cout << "\n\nNS Anisotropy: " << NS_anisotropy[idx_ani] << "\n";
            std::cout << "EW Anisotropy: " << EW_anisotropy[idx_ani] << "\n";
            std::cout << "FB Anisotropy: " << FB_anisotropy[idx_ani] << "\n\n";
            
            output_log_file << "\n\nNS Anisotropy: " << NS_anisotropy[idx_ani] << "\n";
            output_log_file << "EW Anisotropy: " << EW_anisotropy[idx_ani] << "\n";
            output_log_file << "FB Anisotropy: " << FB_anisotropy[idx_ani] << "\n\n";
            
            inSeed >> tmp_seed_str;
            tmp_seed = (UInt_t)std::stoul(tmp_seed_str,nullptr,10);
            
            seed_line = (s_idx*s_batch*ani_values - 1) + itry*ani_values + (idx_ani + 1);
            
            AllSky_test_simulation(
                                    output_log_file,
                                    fiTree,
                                    tmp_fit,
                                    tmp_seed,
                                    itry,
                                    NS_anisotropy[idx_ani],
                                    EW_anisotropy[idx_ani],
                                    FB_anisotropy[idx_ani],
                                    idx_ani,
                                    s_idx,
                                    seed_line
                                   );
            
        }
    }
}

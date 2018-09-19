
#include "MyHead.h"

void create_simulation_seeds()
{
    if(!path_exists(seeds_path))
    {
        Int_t n_seeds = ani_values*ntry;
        Int_t tmp_seed;
        std::vector<UInt_t> vseeds;
        vseeds.resize(n_seeds);
        TRandom2 r_seed((UInt_t)time(0));
        std::ofstream seeds_file(seeds_path);
        if(!seeds_file.is_open()) {
            std::cout << "\n\nCannot create output seeds file! Program finished !" << std::endl;
            exit(100);
        }
        for(Int_t idx_s = 0; idx_s < n_seeds; ++idx_s)
        {
            do
            {
                tmp_seed = r_seed.Integer((UInt_t)(n_seeds*1e+3));
            } while(repetition_seed(tmp_seed,vseeds,idx_s)==true);
            vseeds[idx_s] = tmp_seed;
            seeds_file << tmp_seed << std::endl;
        }
        seeds_file.close();
    }
}


inline bool path_exists (const std::string& path) {
    struct stat buffer;
    return (stat (path.c_str(), &buffer) == 0);
}

bool repetition_seed(UInt_t tmp_seed,std::vector<UInt_t> &vseeds,Int_t position)
{
    bool found = false;
    for(Int_t idx = 0; idx < position; ++idx)
        if(vseeds[idx]==tmp_seed)
        {
            found = true;
            break;
        }
    return found;
}

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <string>
#include <fstream>

#include "TRandom3.h"

bool just_used(UInt_t tmp_seed,std::vector<UInt_t> &seed,Int_t s_idx);

int main(int argc,char* argv[])
{
    Int_t ani_values = 3;
    UInt_t max_generation_value = 1e+6;
    std::string seeds_path(argv[1]);
    Int_t n_try = atoi(argv[2]);
    Int_t n_batch_try = atoi(argv[3]);
    Int_t n_seeds = 3*n_try*n_batch_try;
    UInt_t tmp_seed = 0;
    TRandom3 r_gen(time(0));
    std::ofstream out_seeds(seeds_path.c_str());
    std::vector<UInt_t> seeds;
    
    seeds.resize(n_seeds);
    
    for(Int_t s_idx=0; s_idx < n_seeds; ++s_idx)
    {
        tmp_seed = r_gen.Integer(max_generation_value);
        if(just_used(tmp_seed,seeds,s_idx))
        {
            --s_idx;
            continue;
        }
        else
        {
            seeds[s_idx] = tmp_seed;
            out_seeds << tmp_seed << std::endl;
        }
    }
    
    out_seeds.close();
    
    return 1;
    
}

bool just_used(UInt_t tmp_seed,std::vector<UInt_t> &seed,Int_t s_idx)
{
    bool found = false;
    
    for(Int_t idx=0; idx < s_idx; ++idx)
        if(seed[idx]==tmp_seed)
        {
            found = true;
            break;
        }
    
    return found;
}


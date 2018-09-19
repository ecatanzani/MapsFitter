
#include "MyHead.h"

Double_t get_TF2_max(TF2 &function,Double_t theta)
{
    Double_t max_value;
    max_value = function.GetMaximum(theta);
    return max_value;
}

Double_t get_XY_TF2_max(TF2 &function,Double_t theta,Double_t phi)
{
    Double_t max_value;
    //max_value = function.GetMaximumXY(phi,theta);
    return max_value;
}

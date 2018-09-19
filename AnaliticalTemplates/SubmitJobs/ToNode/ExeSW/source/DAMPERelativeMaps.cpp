
#include "MyHead.h"

void get_relative_histo(TH2D &relative_DAMPE_histo,TH2D* data_DAMPE_histo,TH2D &reference_DAMPE_histo)
{
    std::string htmp_name = "relative_";
    std::string htmp_title = "relative ";
    
    htmp_name.append(data_DAMPE_histo->GetName());
    htmp_title.append(data_DAMPE_histo->GetTitle());
    
    relative_DAMPE_histo.SetTitle(htmp_title.c_str());
    
    TH2D* preference_DAMPE_histo = &reference_DAMPE_histo;
    TH2D* tmp_relative_DAMPE_histo = (TH2D*) data_DAMPE_histo->Clone("tmp_relative_DAMPE_histo");
    
    tmp_relative_DAMPE_histo->Add(preference_DAMPE_histo,-1);
    tmp_relative_DAMPE_histo->Divide(preference_DAMPE_histo);
    
    new (&relative_DAMPE_histo) TH2D (*(TH2D*)tmp_relative_DAMPE_histo->Clone(htmp_name.c_str()));
    
}

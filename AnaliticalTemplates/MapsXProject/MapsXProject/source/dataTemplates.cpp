
#include "MyHead.h"

void compute_data_template(
                           TH2D &h_Template_Data,
                           TH2D &h_Template_Iso,
                           TH2D &h_Template_AniNS,
                           TH2D &h_Template_AniEW,
                           TH2D &h_Template_AniFB,
                           Double_t NS_anisotropy,
                           Double_t EW_anisotropy,
                           Double_t FB_anisotropy,
                           std::string histo_name
                           )

{
    
    TH2D* tmp_Template_Iso = (TH2D*) h_Template_Iso.Clone("tmp_Template_Iso");
    TH2D* tmp_Template_AniNS = (TH2D*) h_Template_AniNS.Clone("tmp_Template_AniNS");
    TH2D* tmp_Template_AniEW = (TH2D*) h_Template_AniEW.Clone("tmp_Template_AniEW");
    TH2D* tmp_Template_AniFB = (TH2D*) h_Template_AniFB.Clone("tmp_Template_AniFB");
    
    TH2D* data_template = (TH2D*) tmp_Template_Iso->Clone("data_template");
    data_template->Reset();
    
    tmp_Template_Iso->Scale(1 - NS_anisotropy - EW_anisotropy - FB_anisotropy);
    tmp_Template_AniNS->Scale(NS_anisotropy);
    tmp_Template_AniEW->Scale(EW_anisotropy);
    tmp_Template_AniFB->Scale(FB_anisotropy);
    
    data_template->Add(tmp_Template_Iso);
    data_template->Add(tmp_Template_AniNS);
    data_template->Add(tmp_Template_AniEW);
    data_template->Add(tmp_Template_AniFB);
    
    new (&h_Template_Data) (TH2D) (*(TH2D*)data_template->Clone(histo_name.c_str()));
    
}



#include "MyHead.h"

void get_cleaned_LS_histos(
                            TH2D* Data_AniNS_LS,
                            TH2D* Data_AniEW_LS,
                            TH2D* Data_AniFB_LS,
                            TH2D &CData_AniNS_LS,
                            TH2D &CData_AniEW_LS,
                            TH2D &CData_AniFB_LS,
                            TH2D* Data_Iso_LS
                           )
{
    TH2D* pAniNS_LS = (TH2D*)Data_AniNS_LS->Clone("CData_AniNS_LS");
    TH2D* pAniEW_LS = (TH2D*)Data_AniEW_LS->Clone("CData_AniEW_LS");
    TH2D* pAniFB_LS = (TH2D*)Data_AniFB_LS->Clone("CData_AniFB_LS");
    
    clean_bin2bin(pAniNS_LS,Data_Iso_LS);
    clean_bin2bin(pAniEW_LS,Data_Iso_LS);
    clean_bin2bin(pAniFB_LS,Data_Iso_LS);
    
    new (&CData_AniNS_LS) (TH2D)(*(TH2D*)pAniNS_LS->Clone("ClonedData_AniNS_LS"));
    new (&CData_AniEW_LS) (TH2D)(*(TH2D*)pAniEW_LS->Clone("ClonedData_AniEW_LS"));
    new (&CData_AniFB_LS) (TH2D)(*(TH2D*)pAniFB_LS->Clone("ClonedData_AniFB_LS"));

}

void get_cleaned_HS_histos(
                            TH2D* Data_AniNS_HS,
                            TH2D* Data_AniEW_HS,
                            TH2D* Data_AniFB_HS,
                            TH2D &CData_AniNS_HS,
                            TH2D &CData_AniEW_HS,
                            TH2D &CData_AniFB_HS,
                            TH2D* Data_Iso_HS
                           )
{
    TH2D* pAniNS_HS = (TH2D*)Data_AniNS_HS->Clone("CData_AniNS_HS");
    TH2D* pAniEW_HS = (TH2D*)Data_AniEW_HS->Clone("CData_AniEW_HS");
    TH2D* pAniFB_HS = (TH2D*)Data_AniFB_HS->Clone("CData_AniFB_HS");
    
    clean_bin2bin(pAniNS_HS,Data_Iso_HS);
    clean_bin2bin(pAniEW_HS,Data_Iso_HS);
    clean_bin2bin(pAniFB_HS,Data_Iso_HS);
    
    new (&CData_AniNS_HS) (TH2D)(*(TH2D*)pAniNS_HS->Clone("ClonedData_AniNS_HS"));
    new (&CData_AniEW_HS) (TH2D)(*(TH2D*)pAniEW_HS->Clone("ClonedData_AniEW_HS"));
    new (&CData_AniFB_HS) (TH2D)(*(TH2D*)pAniFB_HS->Clone("ClonedData_AniFB_HS"));
    
}




void clean_bin2bin(TH2D* HDipole,TH2D* Data_Iso)
{
    for(Int_t bY = 1; bY <= HDipole->GetNbinsY(); ++bY)
        for(Int_t bX = 1; bX <=HDipole->GetNbinsX(); ++bX)
            HDipole->SetBinContent(bX,bY,(Double_t)(HDipole->GetBinContent(bX,bY)-Data_Iso->GetBinContent(bX,bY))/Data_Iso->GetBinContent(bX,bY));
}


#include "MyHead.h"

void normalize_LS_data(TH2D &Data_Iso_LS,TH2D &Data_AniNS_LS,TH2D &Data_AniEW_LS,TH2D &Data_AniFB_LS,TH2D &MixedData_NS_EW_LS,TH2D &MixedData_NS_FB_LS,TH2D &MixedData_EW_FB_LS,TH2D &FullMixedData_LS,std::ofstream &log_file) {
    
    Double_t l_bX = ((Double_t)360/Data_Iso_LS.GetNbinsX())*TMath::DegToRad();
    Double_t l_bY = ((Double_t)180/Data_Iso_LS.GetNbinsX())*TMath::DegToRad();
    
    Double_t integral_I = 0, sf_I = 1;
    Double_t integral_NS = 0, sf_NS = 1;
    Double_t integral_EW = 0, sf_EW = 1;
    Double_t integral_FB = 0, sf_FB = 1;
    Double_t integral_data_NS_EW = 0, sf_data_NS_EW = 1;
    Double_t integral_data_NS_FB = 0, sf_data_NS_FB = 1;
    Double_t integral_data_EW_FB = 0, sf_data_EW_FB = 1;
    Double_t integral_data = 0, sf_data = 1;
    
    /////////////////////////////////////// Isotrope Data normalization
    
    for(Int_t bY=1; bY<=Data_Iso_LS.GetNbinsY(); bY++)
        for(Int_t bX=1; bX<=Data_Iso_LS.GetNbinsX(); bX++)
            integral_I += l_bX*l_bY*Data_Iso_LS.GetBinContent(bX,bY);
    
    sf_I = integral_I/(4*TMath::Pi());
    
    for(Int_t bY=1; bY<=Data_Iso_LS.GetNbinsY(); bY++)
        for(Int_t bX=1; bX<=Data_Iso_LS.GetNbinsX(); bX++)
            Data_Iso_LS.SetBinContent(bX,bY,(Double_t)Data_Iso_LS.GetBinContent(bX,bY)/sf_I);
    
    
    /////////////////////////////////////// NS Data normalization
    
    for(Int_t bY=1; bY<=Data_AniNS_LS.GetNbinsY(); bY++)
        for(Int_t bX=1; bX<=Data_AniNS_LS.GetNbinsX(); bX++)
            integral_NS += l_bX*l_bY*Data_AniNS_LS.GetBinContent(bX,bY);
    
    sf_NS = integral_NS/((18/5.)*TMath::Pi());
    
    for(Int_t bY=1; bY<=Data_AniNS_LS.GetNbinsY(); bY++)
        for(Int_t bX=1; bX<=Data_AniNS_LS.GetNbinsX(); bX++)
            Data_AniNS_LS.SetBinContent(bX,bY,(Double_t)Data_AniNS_LS.GetBinContent(bX,bY)/sf_NS);
    
    
    /////////////////////////////////////// EW Data normalization
    
    for(Int_t bY=1; bY<=Data_AniEW_LS.GetNbinsY(); bY++)
        for(Int_t bX=1; bX<=Data_AniEW_LS.GetNbinsX(); bX++)
            integral_EW += l_bX*l_bY*Data_AniEW_LS.GetBinContent(bX,bY);
    
    sf_EW = integral_EW/((18/5.)*TMath::Pi());
    
    for(Int_t bY=1; bY<=Data_AniEW_LS.GetNbinsY(); bY++)
        for(Int_t bX=1; bX<=Data_AniEW_LS.GetNbinsX(); bX++)
            Data_AniEW_LS.SetBinContent(bX,bY,(Double_t)Data_AniEW_LS.GetBinContent(bX,bY)/sf_EW);
    
    
    /////////////////////////////////////// FB Data normalization
    
    for(Int_t bY=1; bY<=Data_AniFB_LS.GetNbinsY(); bY++)
        for(Int_t bX=1; bX<=Data_AniFB_LS.GetNbinsX(); bX++)
            integral_FB += l_bX*l_bY*Data_AniFB_LS.GetBinContent(bX,bY);
    
    sf_FB = integral_FB/((18/5.)*TMath::Pi());
    
    for(Int_t bY=1; bY<=Data_AniFB_LS.GetNbinsY(); bY++)
        for(Int_t bX=1; bX<=Data_AniFB_LS.GetNbinsX(); bX++)
            Data_AniFB_LS.SetBinContent(bX,bY,(Double_t)Data_AniFB_LS.GetBinContent(bX,bY)/sf_FB);
    
    
    /////////////////////////////////////// 0.8 NS + 0.8 EW Data normalization
    
    for(Int_t bY=1; bY<=MixedData_NS_EW_LS.GetNbinsY(); bY++)
        for(Int_t bX=1; bX<=MixedData_NS_EW_LS.GetNbinsX(); bX++)
            integral_data_NS_EW += l_bX*l_bY*MixedData_NS_EW_LS.GetBinContent(bX,bY);
    
    sf_data_NS_EW = integral_data_NS_EW/((16/5.)*TMath::Pi());
    
    for(Int_t bY=1; bY<=MixedData_NS_EW_LS.GetNbinsY(); bY++)
        for(Int_t bX=1; bX<=MixedData_NS_EW_LS.GetNbinsX(); bX++)
            MixedData_NS_EW_LS.SetBinContent(bX,bY,(Double_t)MixedData_NS_EW_LS.GetBinContent(bX,bY)/sf_data_NS_EW);
    
    
    /////////////////////////////////////// 0.8 NS + 0.8 FB Data normalization
    
    for(Int_t bY=1; bY<=MixedData_NS_FB_LS.GetNbinsY(); bY++)
        for(Int_t bX=1; bX<=MixedData_NS_FB_LS.GetNbinsX(); bX++)
            integral_data_NS_FB += l_bX*l_bY*MixedData_NS_FB_LS.GetBinContent(bX,bY);
    
    sf_data_NS_FB = integral_data_NS_FB/((16/5.)*TMath::Pi());
    
    for(Int_t bY=1; bY<=MixedData_NS_FB_LS.GetNbinsY(); bY++)
        for(Int_t bX=1; bX<=MixedData_NS_FB_LS.GetNbinsX(); bX++)
            MixedData_NS_FB_LS.SetBinContent(bX,bY,(Double_t)MixedData_NS_FB_LS.GetBinContent(bX,bY)/sf_data_NS_FB);
    
    
    /////////////////////////////////////// 0.8 EW + 0.8 FB Data normalization
    
    for(Int_t bY=1; bY<=MixedData_EW_FB_LS.GetNbinsY(); bY++)
        for(Int_t bX=1; bX<=MixedData_EW_FB_LS.GetNbinsX(); bX++)
            integral_data_EW_FB += l_bX*l_bY*MixedData_EW_FB_LS.GetBinContent(bX,bY);
    
    sf_data_EW_FB = integral_data_EW_FB/((16/5.)*TMath::Pi());
    
    for(Int_t bY=1; bY<=MixedData_EW_FB_LS.GetNbinsY(); bY++)
        for(Int_t bX=1; bX<=MixedData_EW_FB_LS.GetNbinsX(); bX++)
            MixedData_EW_FB_LS.SetBinContent(bX,bY,(Double_t)MixedData_EW_FB_LS.GetBinContent(bX,bY)/sf_data_EW_FB);
    
    
    /////////////////////////////////////// Data normalization
    
    for(Int_t bY=1; bY<=FullMixedData_LS.GetNbinsY(); bY++)
        for(Int_t bX=1; bX<=FullMixedData_LS.GetNbinsX(); bX++)
            integral_data += l_bX*l_bY*FullMixedData_LS.GetBinContent(bX,bY);
    
    sf_data = integral_data/((14/5.)*TMath::Pi());
    
    for(Int_t bY=1; bY<=FullMixedData_LS.GetNbinsY(); bY++)
        for(Int_t bX=1; bX<=FullMixedData_LS.GetNbinsX(); bX++)
            FullMixedData_LS.SetBinContent(bX,bY,(Double_t)FullMixedData_LS.GetBinContent(bX,bY)/sf_data);
    
    ////////////////////////////////////////////////////////////////////////////////////////////
    
    std::cout << "/////////////////////////////////////////// Data Normalization ///////////////////////////////////////////" << std::endl;
    std::cout << "//////////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl << std::endl;
    
    std::cout << "\n\nData maps have been normalized !\n";
    std::cout << " --> Isotropic \t integral = " << integral_I << "\t scale factor = " << sf_I << std::endl;
    std::cout << " --> Anisotropic NS \t integral = " << integral_NS << "\t scale factor = " << sf_NS << std::endl;
    std::cout << " --> Anisotropic EW \t integral = " << integral_EW << "\t scale factor = " << sf_EW << std::endl;
    std::cout << " --> Anisotropic FB \t integral = " << integral_FB << "\t scale factor = " << sf_FB << std::endl;
    std::cout << " --> Anisotropic NS+EW \t integral = " << integral_data_NS_EW << "\t scale factor = " << sf_data_NS_EW << std::endl;
    std::cout << " --> Anisotropic NS+FB \t integral = " << integral_data_NS_FB << "\t scale factor = " << sf_data_NS_FB << std::endl;
    std::cout << " --> Anisotropic EW+FB \t integral = " << integral_data_EW_FB << "\t scale factor = " << sf_data_EW_FB << std::endl;
    std::cout << " --> Anisotropic FullMixed \t integral = " << integral_data << "\t scale factor = " << sf_data << std::endl;
    
    std::cout << "\n\n //////////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl << std::endl;
    
    log_file << "/////////////////////////////////////////// Data Normalization ///////////////////////////////////////////" << std::endl;
    log_file << "//////////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl << std::endl;
    
    log_file << "\n\nData maps have been normalized !\n";
    log_file << " --> Isotropic \t integral = " << integral_I << "\t scale factor = " << sf_I << std::endl;
    log_file << " --> Anisotropic NS \t integral = " << integral_NS << "\t scale factor = " << sf_NS << std::endl;
    log_file << " --> Anisotropic EW \t integral = " << integral_EW << "\t scale factor = " << sf_EW << std::endl;
    log_file << " --> Anisotropic FB \t integral = " << integral_FB << "\t scale factor = " << sf_FB << std::endl;
    log_file << " --> Anisotropic NS+EW \t integral = " << integral_data_NS_EW << "\t scale factor = " << sf_data_NS_EW << std::endl;
    log_file << " --> Anisotropic NS+FB \t integral = " << integral_data_NS_FB << "\t scale factor = " << sf_data_NS_FB << std::endl;
    log_file << " --> Anisotropic EW+FB \t integral = " << integral_data_EW_FB << "\t scale factor = " << sf_data_EW_FB << std::endl;
    log_file << " --> Anisotropic FullMixed \t integral = " << integral_data << "\t scale factor = " << sf_data << std::endl;
    
    log_file << "\n\n //////////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl << std::endl;
    
    
}


void normalize_HS_data(TH2D &Data_Iso_HS,TH2D &Data_AniNS_HS,TH2D &Data_AniEW_HS,TH2D &Data_AniFB_HS,TH2D &MixedData_NS_EW_HS,TH2D &MixedData_NS_FB_HS,TH2D &MixedData_EW_FB_HS,TH2D &FullMixedData_HS,std::ofstream &log_file) {
    
    Double_t l_bX = ((Double_t)360/Data_Iso_HS.GetNbinsX())*TMath::DegToRad();
    Double_t l_bY = ((Double_t)180/Data_Iso_HS.GetNbinsX())*TMath::DegToRad();
    
    Double_t integral_I = 0, sf_I = 1;
    Double_t integral_NS = 0, sf_NS = 1;
    Double_t integral_EW = 0, sf_EW = 1;
    Double_t integral_FB = 0, sf_FB = 1;
    Double_t integral_data_NS_EW = 0, sf_data_NS_EW = 1;
    Double_t integral_data_NS_FB = 0, sf_data_NS_FB = 1;
    Double_t integral_data_EW_FB = 0, sf_data_EW_FB = 1;
    Double_t integral_data = 0, sf_data = 1;
    
    /////////////////////////////////////// Isotrope Data normalization
    
    for(Int_t bY=1; bY<=Data_Iso_HS.GetNbinsY(); bY++)
        for(Int_t bX=1; bX<=Data_Iso_HS.GetNbinsX(); bX++)
            integral_I += l_bX*l_bY*Data_Iso_HS.GetBinContent(bX,bY);
    
    sf_I = integral_I/(4*TMath::Pi());
    
    for(Int_t bY=1; bY<=Data_Iso_HS.GetNbinsY(); bY++)
        for(Int_t bX=1; bX<=Data_Iso_HS.GetNbinsX(); bX++)
            Data_Iso_HS.SetBinContent(bX,bY,(Double_t)Data_Iso_HS.GetBinContent(bX,bY)/sf_I);
    
    
    /////////////////////////////////////// NS Data normalization
    
    for(Int_t bY=1; bY<=Data_AniNS_HS.GetNbinsY(); bY++)
        for(Int_t bX=1; bX<=Data_AniNS_HS.GetNbinsX(); bX++)
            integral_NS += l_bX*l_bY*Data_AniNS_HS.GetBinContent(bX,bY);
    
    sf_NS = integral_NS/((18/5.)*TMath::Pi());
    
    for(Int_t bY=1; bY<=Data_AniNS_HS.GetNbinsY(); bY++)
        for(Int_t bX=1; bX<=Data_AniNS_HS.GetNbinsX(); bX++)
            Data_AniNS_HS.SetBinContent(bX,bY,(Double_t)Data_AniNS_HS.GetBinContent(bX,bY)/sf_NS);
    
    
    /////////////////////////////////////// EW Data normalization
    
    for(Int_t bY=1; bY<=Data_AniEW_HS.GetNbinsY(); bY++)
        for(Int_t bX=1; bX<=Data_AniEW_HS.GetNbinsX(); bX++)
            integral_EW += l_bX*l_bY*Data_AniEW_HS.GetBinContent(bX,bY);
    
    sf_EW = integral_EW/((18/5.)*TMath::Pi());
    
    for(Int_t bY=1; bY<=Data_AniEW_HS.GetNbinsY(); bY++)
        for(Int_t bX=1; bX<=Data_AniEW_HS.GetNbinsX(); bX++)
            Data_AniEW_HS.SetBinContent(bX,bY,(Double_t)Data_AniEW_HS.GetBinContent(bX,bY)/sf_EW);
    
    
    /////////////////////////////////////// FB Data normalization
    
    for(Int_t bY=1; bY<=Data_AniFB_HS.GetNbinsY(); bY++)
        for(Int_t bX=1; bX<=Data_AniFB_HS.GetNbinsX(); bX++)
            integral_FB += l_bX*l_bY*Data_AniFB_HS.GetBinContent(bX,bY);
    
    sf_FB = integral_FB/((18/5.)*TMath::Pi());
    
    for(Int_t bY=1; bY<=Data_AniFB_HS.GetNbinsY(); bY++)
        for(Int_t bX=1; bX<=Data_AniFB_HS.GetNbinsX(); bX++)
            Data_AniFB_HS.SetBinContent(bX,bY,(Double_t)Data_AniFB_HS.GetBinContent(bX,bY)/sf_FB);
    
    
    /////////////////////////////////////// 0.8 NS + 0.8 EW Data normalization
    
    for(Int_t bY=1; bY<=MixedData_NS_EW_HS.GetNbinsY(); bY++)
        for(Int_t bX=1; bX<=MixedData_NS_EW_HS.GetNbinsX(); bX++)
            integral_data_NS_EW += l_bX*l_bY*MixedData_NS_EW_HS.GetBinContent(bX,bY);
    
    sf_data_NS_EW = integral_data_NS_EW/((16/5.)*TMath::Pi());
    
    for(Int_t bY=1; bY<=MixedData_NS_EW_HS.GetNbinsY(); bY++)
        for(Int_t bX=1; bX<=MixedData_NS_EW_HS.GetNbinsX(); bX++)
            MixedData_NS_EW_HS.SetBinContent(bX,bY,(Double_t)MixedData_NS_EW_HS.GetBinContent(bX,bY)/sf_data_NS_EW);
    
    
    /////////////////////////////////////// 0.8 NS + 0.8 FB Data normalization
    
    for(Int_t bY=1; bY<=MixedData_NS_FB_HS.GetNbinsY(); bY++)
        for(Int_t bX=1; bX<=MixedData_NS_FB_HS.GetNbinsX(); bX++)
            integral_data_NS_FB += l_bX*l_bY*MixedData_NS_FB_HS.GetBinContent(bX,bY);
    
    sf_data_NS_FB = integral_data_NS_FB/((16/5.)*TMath::Pi());
    
    for(Int_t bY=1; bY<=MixedData_NS_FB_HS.GetNbinsY(); bY++)
        for(Int_t bX=1; bX<=MixedData_NS_FB_HS.GetNbinsX(); bX++)
            MixedData_NS_FB_HS.SetBinContent(bX,bY,(Double_t)MixedData_NS_FB_HS.GetBinContent(bX,bY)/sf_data_NS_FB);
    
    
    /////////////////////////////////////// 0.8 EW + 0.8 FB Data normalization
    
    for(Int_t bY=1; bY<=MixedData_EW_FB_HS.GetNbinsY(); bY++)
        for(Int_t bX=1; bX<=MixedData_EW_FB_HS.GetNbinsX(); bX++)
            integral_data_EW_FB += l_bX*l_bY*MixedData_EW_FB_HS.GetBinContent(bX,bY);
    
    sf_data_EW_FB = integral_data_EW_FB/((16/5.)*TMath::Pi());
    
    for(Int_t bY=1; bY<=MixedData_EW_FB_HS.GetNbinsY(); bY++)
        for(Int_t bX=1; bX<=MixedData_EW_FB_HS.GetNbinsX(); bX++)
            MixedData_EW_FB_HS.SetBinContent(bX,bY,(Double_t)MixedData_EW_FB_HS.GetBinContent(bX,bY)/sf_data_EW_FB);
    
    
    /////////////////////////////////////// Data normalization
    
    for(Int_t bY=1; bY<=FullMixedData_HS.GetNbinsY(); bY++)
        for(Int_t bX=1; bX<=FullMixedData_HS.GetNbinsX(); bX++)
            integral_data += l_bX*l_bY*FullMixedData_HS.GetBinContent(bX,bY);
    
    sf_data = integral_data/((14/5.)*TMath::Pi());
    
    for(Int_t bY=1; bY<=FullMixedData_HS.GetNbinsY(); bY++)
        for(Int_t bX=1; bX<=FullMixedData_HS.GetNbinsX(); bX++)
            FullMixedData_HS.SetBinContent(bX,bY,(Double_t)FullMixedData_HS.GetBinContent(bX,bY)/sf_data);
    
    ////////////////////////////////////////////////////////////////////////////////////////////
    
    std::cout << "/////////////////////////////////////////// Data Normalization ///////////////////////////////////////////" << std::endl;
    std::cout << "//////////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl << std::endl;
    
    std::cout << "\n\nData maps have been normalized !\n";
    std::cout << " --> Isotropic \t integral = " << integral_I << "\t scale factor = " << sf_I << std::endl;
    std::cout << " --> Anisotropic NS \t integral = " << integral_NS << "\t scale factor = " << sf_NS << std::endl;
    std::cout << " --> Anisotropic EW \t integral = " << integral_EW << "\t scale factor = " << sf_EW << std::endl;
    std::cout << " --> Anisotropic FB \t integral = " << integral_FB << "\t scale factor = " << sf_FB << std::endl;
    std::cout << " --> Anisotropic NS+EW \t integral = " << integral_data_NS_EW << "\t scale factor = " << sf_data_NS_EW << std::endl;
    std::cout << " --> Anisotropic NS+FB \t integral = " << integral_data_NS_FB << "\t scale factor = " << sf_data_NS_FB << std::endl;
    std::cout << " --> Anisotropic EW+FB \t integral = " << integral_data_EW_FB << "\t scale factor = " << sf_data_EW_FB << std::endl;
    std::cout << " --> Anisotropic FullMixed \t integral = " << integral_data << "\t scale factor = " << sf_data << std::endl;
    
    std::cout << "\n\n //////////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl << std::endl;
    
    log_file << "/////////////////////////////////////////// Data Normalization ///////////////////////////////////////////" << std::endl;
    log_file << "//////////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl << std::endl;
    
    log_file << "\n\nData maps have been normalized !\n";
    log_file << " --> Isotropic \t integral = " << integral_I << "\t scale factor = " << sf_I << std::endl;
    log_file << " --> Anisotropic NS \t integral = " << integral_NS << "\t scale factor = " << sf_NS << std::endl;
    log_file << " --> Anisotropic EW \t integral = " << integral_EW << "\t scale factor = " << sf_EW << std::endl;
    log_file << " --> Anisotropic FB \t integral = " << integral_FB << "\t scale factor = " << sf_FB << std::endl;
    log_file << " --> Anisotropic NS+EW \t integral = " << integral_data_NS_EW << "\t scale factor = " << sf_data_NS_EW << std::endl;
    log_file << " --> Anisotropic NS+FB \t integral = " << integral_data_NS_FB << "\t scale factor = " << sf_data_NS_FB << std::endl;
    log_file << " --> Anisotropic EW+FB \t integral = " << integral_data_EW_FB << "\t scale factor = " << sf_data_EW_FB << std::endl;
    log_file << " --> Anisotropic FullMixed \t integral = " << integral_data << "\t scale factor = " << sf_data << std::endl;
    
    log_file << "\n\n //////////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl << std::endl;
    
    
}


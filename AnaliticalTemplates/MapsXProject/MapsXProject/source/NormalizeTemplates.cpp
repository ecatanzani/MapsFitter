
#include "MyHead.h"

void normalize_LS_templates(TH2D &Template_Iso_LS,TH2D &Template_AniNS_LS,TH2D &Template_AniEW_LS,TH2D &Template_AniFB_LS,std::ofstream &log_file) {
    
    Double_t l_bX = ((Double_t)360/Template_Iso_LS.GetNbinsX())*TMath::DegToRad();
    Double_t l_bY = ((Double_t)180/Template_Iso_LS.GetNbinsX())*TMath::DegToRad();
    
    Double_t integral_I = 0, sf_I = 1;
    Double_t integral_NS = 0, sf_NS = 1;
    Double_t integral_EW = 0, sf_EW = 1;
    Double_t integral_FB = 0, sf_FB = 1;
    TAxis* theta_axis = nullptr;
    Double_t Jacobian = 0;
    
    /////////////////////////////////////// Isotrope Template normalization
    
    for(Int_t bY=1; bY<=Template_Iso_LS.GetNbinsY(); bY++) {
        theta_axis = Template_Iso_LS.GetYaxis();
        Jacobian = TMath::Sin(theta_axis->GetBinCenter(bY)*TMath::DegToRad());
        for(Int_t bX=1; bX<=Template_Iso_LS.GetNbinsX(); bX++)
            integral_I += l_bX*l_bY*Template_Iso_LS.GetBinContent(bX,bY);
    }
    sf_I = (integral_I*TMath::DegToRad())/(4*TMath::Pi());
    
    for(Int_t bY=1; bY<=Template_Iso_LS.GetNbinsY(); bY++)
        for(Int_t bX=1; bX<=Template_Iso_LS.GetNbinsX(); bX++)
            Template_Iso_LS.SetBinContent(bX,bY,(Double_t)Template_Iso_LS.GetBinContent(bX,bY)/sf_I);
    
    /////////////////////////////////////// NS Template normalization
    
    for(Int_t bY=1; bY<=9; bY++) {
        theta_axis = Template_AniNS_LS.GetYaxis();
        Jacobian = TMath::Sin(theta_axis->GetBinCenter(bY)*TMath::DegToRad());
        for(Int_t bX=1; bX<=Template_AniNS_LS.GetNbinsX(); bX++)
            integral_NS += l_bX*l_bY*Template_AniNS_LS.GetBinContent(bX,bY);
    }
    sf_NS = (integral_NS*TMath::DegToRad())/TMath::Sqrt((3/4.)*TMath::Pi());
    
    
    for(Int_t bY=1; bY<=Template_AniNS_LS.GetNbinsY(); bY++)
        for(Int_t bX=1; bX<=Template_AniNS_LS.GetNbinsX(); bX++)
            Template_AniNS_LS.SetBinContent(bX,bY,(Double_t)Template_AniNS_LS.GetBinContent(bX,bY)/sf_NS);
    
    
    /////////////////////////////////////// EW Template normalization
   
    for(Int_t bY=1; bY<=Template_AniEW_LS.GetNbinsY(); bY++) {
        theta_axis = Template_AniEW_LS.GetYaxis();
        Jacobian = TMath::Sin(theta_axis->GetBinCenter(bY)*TMath::DegToRad());
        for(Int_t bX=1; bX<=18; bX++)
            integral_EW += l_bX*l_bY*Template_AniEW_LS.GetBinContent(bX,bY);
    }
    sf_EW = (integral_EW*TMath::DegToRad())/TMath::Sqrt((3/2.)*TMath::Pi());
    
    for(Int_t bY=1; bY<=Template_AniEW_LS.GetNbinsY(); bY++)
        for(Int_t bX=1; bX<=Template_AniEW_LS.GetNbinsX(); bX++)
            Template_AniEW_LS.SetBinContent(bX,bY,(Double_t)Template_AniEW_LS.GetBinContent(bX,bY)/sf_EW);
    
    /////////////////////////////////////// FB Template normalization
    
    for(Int_t bY=1; bY<=Template_AniFB_LS.GetNbinsY(); bY++) {
        theta_axis = Template_AniFB_LS.GetYaxis();
        Jacobian = TMath::Sin(theta_axis->GetBinCenter(bY)*TMath::DegToRad());
        for(Int_t bX=10; bX<=27; bX++)
            integral_FB += l_bX*l_bY*Template_AniFB_LS.GetBinContent(bX,bY);
    }
    sf_FB = (integral_FB*TMath::DegToRad())/TMath::Sqrt((3/2.));
    
     
    for(Int_t bY=1; bY<=Template_AniFB_LS.GetNbinsY(); bY++)
        for(Int_t bX=1; bX<=Template_AniFB_LS.GetNbinsX(); bX++)
            Template_AniFB_LS.SetBinContent(bX,bY,(Double_t)Template_AniFB_LS.GetBinContent(bX,bY)/sf_FB);

    
    ////////////////////////////////////////////////////////////////////////////////////////////
    
    std::cout << "\n\n";
    log_file << "\n\n";
    
    std::cout << "///////////////////////////////////////// Template Normalization /////////////////////////////////////////" << std::endl;
    std::cout << "//////////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl << std::endl;
    
    std::cout << "\n\nTemplate maps have been normalized !\n";
    std::cout << " --> Isotropic \t integral = " <<integral_I << "\t scale factor = " << sf_I << std::endl;
    std::cout << " --> Anisotropic NS \t integral = " <<integral_NS << "\t scale factor = " << sf_NS << std::endl;
    std::cout << " --> Anisotropic EW \t integral = " <<integral_EW << "\t scale factor = " << sf_EW << std::endl;
    std::cout << " --> Anisotropic FB \t integral = " <<integral_FB << "\t scale factor = " << sf_FB << std::endl;
    
    std::cout << "\n\n //////////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl << std::endl;
    
    log_file << "///////////////////////////////////////// Template Normalization /////////////////////////////////////////" << std::endl;
    log_file << "//////////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl << std::endl;
    
    log_file << "\n\nTemplate maps have been normalized !\n";
    log_file << " --> Isotropic \t integral = " <<integral_I << "\t scale factor = " << sf_I << std::endl;
    log_file << " --> Anisotropic NS \t integral = " <<integral_NS << "\t scale factor = " << sf_NS << std::endl;
    log_file << " --> Anisotropic EW \t integral = " <<integral_EW << "\t scale factor = " << sf_EW << std::endl;
    log_file << " --> Anisotropic FB \t integral = " <<integral_FB << "\t scale factor = " << sf_FB << std::endl;
    
    log_file << "\n\n //////////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl << std::endl;
    
}


void normalize_HS_templates(TH2D &Template_Iso_HS,TH2D &Template_AniNS_HS,TH2D &Template_AniEW_HS,TH2D &Template_AniFB_HS,std::ofstream &log_file) {
    
    //Double_t l_bX = ((Double_t)360/Template_Iso_HS.GetNbinsX())*TMath::DegToRad();
    //Double_t l_bY = ((Double_t)180/Template_Iso_HS.GetNbinsX())*TMath::DegToRad();
    Double_t l_bX = ((Double_t)360/Template_Iso_HS.GetNbinsX());
    Double_t l_bY = ((Double_t)180/Template_Iso_HS.GetNbinsX());
    
    Double_t integral_I = 0, sf_I = 1;
    Double_t integral_NS = 0, sf_NS = 1;
    Double_t integral_EW = 0, sf_EW = 1;
    Double_t integral_FB = 0, sf_FB = 1;
    TAxis* theta_axis = nullptr;
    Double_t Jacobian = 0;
    
    /////////////////////////////////////// Isotrope Template normalization
    
    for(Int_t bY=1; bY<=Template_Iso_HS.GetNbinsY(); bY++) {
        theta_axis = Template_Iso_HS.GetYaxis();
        Jacobian = TMath::Cos(theta_axis->GetBinCenter(bY)*TMath::DegToRad());
        for(Int_t bX=1; bX<=Template_Iso_HS.GetNbinsX(); bX++)
            integral_I += (l_bX)*(l_bY)*Jacobian*Template_Iso_HS.GetBinContent(bX,bY);
    }
    //sf_I = (integral_I)/(4*TMath::Pi());
    sf_I = (integral_I)/(4*TMath::Pi()*2);
    
    for(Int_t bY=1; bY<=Template_Iso_HS.GetNbinsY(); bY++)
        for(Int_t bX=1; bX<=Template_Iso_HS.GetNbinsX(); bX++)
            Template_Iso_HS.SetBinContent(bX,bY,(Double_t)Template_Iso_HS.GetBinContent(bX,bY)/sf_I);
    
    
    /////////////////////////////////////// NS Template normalization
    
    for(Int_t bY=1; bY<=9; bY++) {
        theta_axis = Template_AniNS_HS.GetYaxis();
        Jacobian = TMath::Cos(theta_axis->GetBinCenter(bY)*TMath::DegToRad());
        for(Int_t bX=1; bX<=Template_AniNS_HS.GetNbinsX(); bX++)
            integral_NS += l_bX*l_bY*Jacobian*Template_AniNS_HS.GetBinContent(bX,bY);
    }
    sf_NS = (integral_NS*TMath::DegToRad())/TMath::Sqrt((3/4.)*TMath::Pi());
    
    
    for(Int_t bY=1; bY<=Template_AniNS_HS.GetNbinsY(); bY++)
        for(Int_t bX=1; bX<=Template_AniNS_HS.GetNbinsX(); bX++)
            Template_AniNS_HS.SetBinContent(bX,bY,(Double_t)Template_AniNS_HS.GetBinContent(bX,bY)/sf_NS);
    
    
    /////////////////////////////////////// EW Template normalization
    
    for(Int_t bY=1; bY<=Template_AniEW_HS.GetNbinsY(); bY++) {
        theta_axis = Template_AniEW_HS.GetYaxis();
        Jacobian = TMath::Sin(theta_axis->GetBinCenter(bY)*TMath::DegToRad());
        for(Int_t bX=1; bX<=18; bX++)
            integral_EW += l_bX*l_bY*Template_AniEW_HS.GetBinContent(bX,bY);
    }
    sf_EW = (integral_EW*TMath::DegToRad())/TMath::Sqrt((3/2.)*TMath::Pi());
    
    for(Int_t bY=1; bY<=Template_AniEW_HS.GetNbinsY(); bY++)
        for(Int_t bX=1; bX<=Template_AniEW_HS.GetNbinsX(); bX++)
            Template_AniEW_HS.SetBinContent(bX,bY,(Double_t)Template_AniEW_HS.GetBinContent(bX,bY)/sf_EW);
    
    
    /////////////////////////////////////// FB Template normalization
    
    for(Int_t bY=1; bY<=Template_AniFB_HS.GetNbinsY(); bY++) {
        theta_axis = Template_AniFB_HS.GetYaxis();
        Jacobian = TMath::Sin(theta_axis->GetBinCenter(bY)*TMath::DegToRad());
        for(Int_t bX=10; bX<=27; bX++)
            integral_FB += l_bX*l_bY*Template_AniFB_HS.GetBinContent(bX,bY);
    }
    sf_FB = (integral_FB*TMath::DegToRad())/TMath::Sqrt((3/2.)*TMath::Pi());
    
    for(Int_t bY=1; bY<=Template_AniFB_HS.GetNbinsY(); bY++)
        for(Int_t bX=1; bX<=Template_AniFB_HS.GetNbinsX(); bX++)
            Template_AniFB_HS.SetBinContent(bX,bY,(Double_t)Template_AniFB_HS.GetBinContent(bX,bY)/sf_FB);
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////
    
    std::cout << "\n\n";
    log_file << "\n\n";
    
    std::cout << "///////////////////////////////////////// Template Normalization /////////////////////////////////////////" << std::endl;
    std::cout << "//////////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl << std::endl;
    
    std::cout << "\n\nTemplate maps have been normalized !\n";
    std::cout << " --> Isotropic \t integral = " <<integral_I << "\t scale factor = " << sf_I << std::endl;
    std::cout << " --> Anisotropic NS \t integral = " <<integral_NS << "\t scale factor = " << sf_NS << std::endl;
    std::cout << " --> Anisotropic EW \t integral = " <<integral_EW << "\t scale factor = " << sf_EW << std::endl;
    std::cout << " --> Anisotropic FB \t integral = " <<integral_FB << "\t scale factor = " << sf_FB << std::endl;
    
    std::cout << "\n\n //////////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl << std::endl;
    
    log_file << "///////////////////////////////////////// Template Normalization /////////////////////////////////////////" << std::endl;
    log_file << "//////////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl << std::endl;
    
    log_file << "\n\nTemplate maps have been normalized !\n";
    log_file << " --> Isotropic \t integral = " <<integral_I << "\t scale factor = " << sf_I << std::endl;
    log_file << " --> Anisotropic NS \t integral = " <<integral_NS << "\t scale factor = " << sf_NS << std::endl;
    log_file << " --> Anisotropic EW \t integral = " <<integral_EW << "\t scale factor = " << sf_EW << std::endl;
    log_file << " --> Anisotropic FB \t integral = " <<integral_FB << "\t scale factor = " << sf_FB << std::endl;
    
    log_file << "\n\n //////////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl << std::endl;
    
}



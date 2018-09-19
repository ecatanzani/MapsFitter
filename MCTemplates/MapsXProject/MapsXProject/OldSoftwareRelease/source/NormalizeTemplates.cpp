
#include "MyHead.h"

void normalize_templates(TH2D &TemplateIso,TH2D &TemplateAniNS,TH2D &TemplateAniEW,TH2D &TemplateAniFB,std::ofstream &log_file) {
    
    Double_t l_bX = ((Double_t)360/TemplateIso.GetNbinsX())*TMath::DegToRad();
    Double_t l_bY = ((Double_t)180/TemplateIso.GetNbinsX())*TMath::DegToRad();
    
    Double_t integral_I = 0, sf_I = 1;
    Double_t integral_NS = 0, sf_NS = 1;
    Double_t integral_EW = 0, sf_EW = 1;
    Double_t integral_FB = 0, sf_FB = 1;
    
    
    /////////////////////////////////////// Isotrope Template normalization
    
    for(Int_t bY=1; bY<=TemplateIso.GetNbinsY(); bY++)
        for(Int_t bX=1; bX<=TemplateIso.GetNbinsX(); bX++)
            integral_I += l_bX*l_bY*TemplateIso.GetBinContent(bX,bY);
    
    sf_I = integral_I;
    
    for(Int_t bY=1; bY<=TemplateIso.GetNbinsY(); bY++)
        for(Int_t bX=1; bX<=TemplateIso.GetNbinsX(); bX++)
            TemplateIso.SetBinContent(bX,bY,(Double_t)TemplateIso.GetBinContent(bX,bY)/sf_I);
    
    
    /////////////////////////////////////// NS Template normalization
    
    for(Int_t bY=1; bY<=9; bY++)
        for(Int_t bX=1; bX<=TemplateAniNS.GetNbinsX(); bX++)
            integral_NS += l_bX*l_bY*TemplateAniNS.GetBinContent(bX,bY);
    
    sf_NS = integral_NS/TMath::Sqrt((3/4.)*TMath::Pi());
    
    for(Int_t bY=1; bY<=TemplateAniNS.GetNbinsY(); bY++)
        for(Int_t bX=1; bX<=TemplateAniNS.GetNbinsX(); bX++)
            TemplateAniNS.SetBinContent(bX,bY,(Double_t)TemplateAniNS.GetBinContent(bX,bY)/sf_NS);
    
    
    /////////////////////////////////////// EW Template normalization
    
    for(Int_t bY=1; bY<=TemplateAniEW.GetNbinsY(); bY++)
        for(Int_t bX=1; bX<=18; bX++)
            integral_EW += l_bX*l_bY*TemplateAniEW.GetBinContent(bX,bY);
    
    sf_EW = integral_EW/TMath::Sqrt(6*TMath::Pi());

    for(Int_t bY=1; bY<=TemplateAniEW.GetNbinsY(); bY++)
        for(Int_t bX=1; bX<=TemplateAniEW.GetNbinsX(); bX++)
            TemplateAniEW.SetBinContent(bX,bY,(Double_t)TemplateAniEW.GetBinContent(bX,bY)/sf_EW);
    
    
    /////////////////////////////////////// FB Template normalization
    
    for(Int_t bY=1; bY<=TemplateAniFB.GetNbinsY(); bY++)
        for(Int_t bX=10; bX<=27; bX++)
            integral_FB += l_bX*l_bY*TemplateAniFB.GetBinContent(bX,bY);
    
    sf_FB = integral_FB/TMath::Sqrt((3/2.)*TMath::Pi());

    for(Int_t bY=1; bY<=TemplateAniFB.GetNbinsY(); bY++)
        for(Int_t bX=1; bX<=TemplateAniFB.GetNbinsX(); bX++)
            TemplateAniFB.SetBinContent(bX,bY,(Double_t)TemplateAniFB.GetBinContent(bX,bY)/sf_FB);

    
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

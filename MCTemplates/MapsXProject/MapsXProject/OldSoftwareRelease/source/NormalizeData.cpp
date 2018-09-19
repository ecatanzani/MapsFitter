
#include "MyHead.h"

void normalize_data(TH2D &dataI,TH2D &dataNS,TH2D &dataEW,TH2D &dataFB,std::ofstream &log_file) {
    
    Double_t l_bX = ((Double_t)360/dataNS.GetNbinsX())*TMath::DegToRad();
    Double_t l_bY = ((Double_t)180/dataNS.GetNbinsX())*TMath::DegToRad();
    
    Double_t integral_I = 0, sf_I = 1;
    Double_t integral_NS = 0, sf_NS = 1;
    Double_t integral_EW = 0, sf_EW = 1;
    Double_t integral_FB = 0, sf_FB = 1;
    
    /////////////////////////////////////// Isotrope Data normalization
    
    for(Int_t bY=1; bY<=dataI.GetNbinsY(); bY++)
        for(Int_t bX=1; bX<=dataI.GetNbinsX(); bX++)
            integral_I += l_bX*l_bY*dataI.GetBinContent(bX,bY);
    
    sf_I = integral_I;
    
    for(Int_t bY=1; bY<=dataI.GetNbinsY(); bY++)
        for(Int_t bX=1; bX<=dataI.GetNbinsX(); bX++)
            dataI.SetBinContent(bX,bY,(Double_t)dataI.GetBinContent(bX,bY)/sf_I);
    
    
    /////////////////////////////////////// NS Data normalization
    
    for(Int_t bY=1; bY<=9; bY++)
        for(Int_t bX=1; bX<=dataNS.GetNbinsX(); bX++)
            integral_NS += l_bX*l_bY*dataNS.GetBinContent(bX,bY);
    
    sf_NS = integral_NS/TMath::Sqrt((3/4.)*TMath::Pi());
    
    for(Int_t bY=1; bY<=dataNS.GetNbinsY(); bY++)
        for(Int_t bX=1; bX<=dataNS.GetNbinsX(); bX++)
            dataNS.SetBinContent(bX,bY,(Double_t)dataNS.GetBinContent(bX,bY)/sf_NS);
    
    
    /////////////////////////////////////// EW Data normalization
    
    for(Int_t bY=1; bY<=dataEW.GetNbinsY(); bY++)
        for(Int_t bX=1; bX<=18; bX++)
            integral_EW += l_bX*l_bY*dataEW.GetBinContent(bX,bY);
    
    sf_EW = integral_EW/TMath::Sqrt(6*TMath::Pi());
    
    for(Int_t bY=1; bY<=dataEW.GetNbinsY(); bY++)
        for(Int_t bX=1; bX<=dataEW.GetNbinsX(); bX++)
            dataEW.SetBinContent(bX,bY,(Double_t)dataEW.GetBinContent(bX,bY)/sf_EW);
    
    
    /////////////////////////////////////// FB Data normalization
    
    for(Int_t bY=1; bY<=dataFB.GetNbinsY(); bY++)
        for(Int_t bX=10; bX<=27; bX++)
            integral_FB += l_bX*l_bY*dataFB.GetBinContent(bX,bY);
    
    sf_FB = integral_FB/TMath::Sqrt((3/2.)*TMath::Pi());
    
    for(Int_t bY=1; bY<=dataFB.GetNbinsY(); bY++)
        for(Int_t bX=1; bX<=dataFB.GetNbinsX(); bX++)
            dataFB.SetBinContent(bX,bY,(Double_t)dataFB.GetBinContent(bX,bY)/sf_FB);
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////
    
    std::cout << "/////////////////////////////////////////// Data Normalization ///////////////////////////////////////////" << std::endl;
    std::cout << "//////////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl << std::endl;
    
    std::cout << "\n\nData maps have been normalized !\n";
    std::cout << " --> Isotropic \t integral = " <<integral_I << "\t scale factor = " << sf_I << std::endl;
    std::cout << " --> Anisotropic NS \t integral = " <<integral_NS << "\t scale factor = " << sf_NS << std::endl;
    std::cout << " --> Anisotropic EW \t integral = " <<integral_EW << "\t scale factor = " << sf_EW << std::endl;
    std::cout << " --> Anisotropic FB \t integral = " <<integral_FB << "\t scale factor = " << sf_FB << std::endl;
    
    std::cout << "\n\n //////////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl << std::endl;
    
    log_file << "/////////////////////////////////////////// Data Normalization ///////////////////////////////////////////" << std::endl;
    log_file << "//////////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl << std::endl;
    
    log_file << "\n\nData maps have been normalized !\n";
    log_file << " --> Isotropic \t integral = " <<integral_I << "\t scale factor = " << sf_I << std::endl;
    log_file << " --> Anisotropic NS \t integral = " <<integral_NS << "\t scale factor = " << sf_NS << std::endl;
    log_file << " --> Anisotropic EW \t integral = " <<integral_EW << "\t scale factor = " << sf_EW << std::endl;
    log_file << " --> Anisotropic FB \t integral = " <<integral_FB << "\t scale factor = " << sf_FB << std::endl;
    
    log_file << "\n\n //////////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl << std::endl;
    
    
}

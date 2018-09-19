
#include "MyHead.h"

/*
 
TH1D* TH2toTH1(  TH2D* h2 ){
    Int_t nbins = (Int_t)h2->GetNbinsX()*(Int_t)h2->GetNbinsY();
    TH1D *h1 = new TH1D( Form("%s_TH1",h2->GetName()), h2->GetTitle(), nbins, 0, nbins);
    for(int ibinx=1; ibinx<h2->GetNbinsX(); ibinx++){
        for(int ibiny=1; ibiny<h2->GetNbinsY(); ibiny++){
            Int_t ibin = ibinx + h2->GetNbinsX()*(ibiny-1);
            h1->SetBinContent( ibin, h2->GetBinContent(ibinx,ibiny) );
            h1->SetBinError( ibin, h2->GetBinError(ibinx,ibiny) );
        }
    }
    return h1;
}

TH1D* TH2toTH1(  TH2D  &h2 ){
    Int_t nbins = (Int_t)h2.GetNbinsX()*(Int_t)h2.GetNbinsY();
    TH1D *h1 = new TH1D( Form("%s_TH1",h2.GetName()), h2.GetTitle(), nbins, 0, nbins);
    for(int ibinx=1; ibinx<h2.GetNbinsX(); ibinx++){
        for(int ibiny=1; ibiny<h2.GetNbinsY(); ibiny++){
            Int_t ibin = ibinx + h2.GetNbinsX()*(ibiny-1);
            h1->SetBinContent( ibin, h2.GetBinContent(ibinx,ibiny) );
            h1->SetBinError( ibin, h2.GetBinError(ibinx,ibiny) );
        }
    }
    return h1;
}

*/

void TH2toTH1_ptr(TH1D &Histo1D,TH2D *Histo2D) {
    Int_t nbins = (Int_t)Histo2D->GetNbinsX()*(Int_t)Histo2D->GetNbinsY();
    TH1D tmpHisto1D(Form("%s_TH1",Histo2D->GetName()),Histo2D->GetTitle(),nbins,0,nbins);
    for(int ibinx=1; ibinx<=Histo2D->GetNbinsX(); ibinx++){
        for(int ibiny=1; ibiny<=Histo2D->GetNbinsY(); ibiny++){
            Int_t ibin = ibinx + Histo2D->GetNbinsX()*(ibiny-1);
            tmpHisto1D.SetBinContent(ibin,Histo2D->GetBinContent(ibinx,ibiny));
            tmpHisto1D.SetBinError(ibin,Histo2D->GetBinError(ibinx,ibiny));
        }
    }
    new (&Histo1D) (TH1D)(*(TH1D*)tmpHisto1D.Clone(Form("%s_1D",Histo2D->GetName())));
}

void TH2toTH1_obj(TH1D &Histo1D,TH2D &Histo2D) {
    Int_t nbins = (Int_t)Histo2D.GetNbinsX()*(Int_t)Histo2D.GetNbinsY();
    TH1D tmpHisto1D(Form("%s_TH1",Histo2D.GetName()),Histo2D.GetTitle(),nbins,0,nbins);
    for(int ibinx=1; ibinx<=Histo2D.GetNbinsX(); ibinx++){
        for(int ibiny=1; ibiny<=Histo2D.GetNbinsY(); ibiny++){
            Int_t ibin = ibinx + Histo2D.GetNbinsX()*(ibiny-1);
            tmpHisto1D.SetBinContent(ibin,Histo2D.GetBinContent(ibinx,ibiny));
            tmpHisto1D.SetBinError(ibin,Histo2D.GetBinError(ibinx,ibiny));
        }
    }
    new (&Histo1D) (TH1D)(*(TH1D*)tmpHisto1D.Clone(Form("%s_1D",Histo2D.GetName())));
}

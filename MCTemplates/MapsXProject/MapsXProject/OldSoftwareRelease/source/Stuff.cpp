
#include "MyHead.h"

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

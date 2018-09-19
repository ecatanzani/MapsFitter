
#include "MyHead.h"

void maps_ratio(TH2D &RatioAniMap_NS,TH2D &NAniMap_NS,TH2D &phNAniMap_NS,TH2D &RatioAniMap_EW,TH2D &NAniMap_EW,TH2D &phNAniMap_EW,TH2D &RatioAniMap_FB,TH2D &NAniMap_FB,TH2D &phNAniMap_FB,TH2D &Isolate_NS,TH2D &Isolate_EW,TH2D &Isolate_FB,TH2D &NIsoMap) {
    
    TH2D* pNIsoMap = (TH2D*)NIsoMap.Clone("pNIsoMap");
    
    TH2D* tmpR_NS = (TH2D*)NAniMap_NS.Clone("tmpR_NS");
    TH2D* tmpR_EW = (TH2D*)NAniMap_EW.Clone("tmpR_EW");
    TH2D* tmpR_FB = (TH2D*)NAniMap_FB.Clone("tmpR_FB");
    
    TH2D* pIsolate_NS = (TH2D*)NAniMap_NS.Clone("pIsolate_NS");
    TH2D* pIsolate_EW = (TH2D*)NAniMap_EW.Clone("pIsolate_EW");
    TH2D* pIsolate_FB = (TH2D*)NAniMap_FB.Clone("pIsolate_FB");
    
    TH2D* p_phNAniMap_NS = &phNAniMap_NS;
    TH2D* p_phNAniMap_EW = &phNAniMap_EW;
    TH2D* p_phNAniMap_FB = &phNAniMap_FB;
    
    tmpR_NS->Divide(p_phNAniMap_NS);
    tmpR_EW->Divide(p_phNAniMap_EW);
    tmpR_FB->Divide(p_phNAniMap_FB);
    
    new (&RatioAniMap_NS) TH2D(*(TH2D*)tmpR_NS->Clone("RatioAniMap_NS"));
    new (&RatioAniMap_EW) TH2D(*(TH2D*)tmpR_EW->Clone("RatioAniMap_EW"));
    new (&RatioAniMap_FB) TH2D(*(TH2D*)tmpR_FB->Clone("RatioAniMap_FB"));
    
    ////////////////////////// To show NS anisotropy //////////////////////////
    
    BinXBin_difference(pIsolate_NS, pNIsoMap);
    pIsolate_NS->Divide(pNIsoMap);
    new (&Isolate_NS) (TH2D)(*(TH2D*)pIsolate_NS->Clone("Isolate_NS"));
    
    ////////////////////////// To show EW anisotropy //////////////////////////
    
    BinXBin_difference(pIsolate_EW, pNIsoMap);
    pIsolate_EW->Divide(pNIsoMap);
    new (&Isolate_EW) (TH2D)(*(TH2D*)pIsolate_EW->Clone("Isolate_EW"));
    
    ////////////////////////// To show FB anisotropy //////////////////////////
    
    BinXBin_difference(pIsolate_FB, pNIsoMap);
    pIsolate_FB->Divide(pNIsoMap);
    new (&Isolate_FB) (TH2D)(*(TH2D*)pIsolate_FB->Clone("Isolate_FB"));
    
}

void BinXBin_difference(TH2D *AniMap,TH2D *IsoMap) {
    
    for(Int_t bX=1; bX<=AniMap->GetNbinsX(); bX++)
        for(Int_t bY=1; bY<=AniMap->GetNbinsY(); bY++)
            AniMap->SetBinContent(bX, bY,(Double_t)AniMap->GetBinContent(bX,bY)-(Double_t)IsoMap->GetBinContent(bX,bY));
    
}

void isolate_dipoles(TH2D &IsoMap,TH2D &AniMap_NS,TH2D &AniMap_EW,TH2D &AniMap_FB,TH2D &Isolate_NN_NS,TH2D &Isolate_NN_EW,TH2D &Isolate_NN_FB) {
    
    TH2D* pIsoMap = (TH2D*)IsoMap.Clone("pIsoMap");
    
    TH2D* pIsolate_NS = (TH2D*)AniMap_NS.Clone("pIsolate_NS");
    TH2D* pIsolate_EW = (TH2D*)AniMap_EW.Clone("pIsolate_EW");
    TH2D* pIsolate_FB = (TH2D*)AniMap_FB.Clone("pIsolate_FB");
    
    ////////////////////////// To show NS anisotropy //////////////////////////
    
    BinXBin_difference(pIsolate_NS,pIsoMap);
    pIsolate_NS->Divide(pIsoMap);
    new (&Isolate_NN_NS) (TH2D)(*(TH2D*)pIsolate_NS->Clone("Isolate_NN_NS"));

    ////////////////////////// To show NS anisotropy //////////////////////////
    
    BinXBin_difference(pIsolate_EW,pIsoMap);
    pIsolate_EW->Divide(pIsoMap);
    new (&Isolate_NN_EW) (TH2D)(*(TH2D*)pIsolate_EW->Clone("Isolate_NN_EW"));
    
    ////////////////////////// To show NS anisotropy //////////////////////////
    
    BinXBin_difference(pIsolate_FB,pIsoMap);
    pIsolate_FB->Divide(pIsoMap);
    new (&Isolate_NN_FB) (TH2D)(*(TH2D*)pIsolate_FB->Clone("Isolate_NN_FB"));

}

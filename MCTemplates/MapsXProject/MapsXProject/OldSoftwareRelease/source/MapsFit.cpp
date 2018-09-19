
#include "MyHead.h"

/////////////////////////////// Define vectors for TMinuit Minimization

std::vector<std::vector<Double_t> > binEntries_NS = std::vector<std::vector<Double_t> >(360, std::vector<Double_t>(180));
std::vector<std::vector<Double_t> > binEntries_EW = std::vector<std::vector<Double_t> >(360, std::vector<Double_t>(180));
std::vector<std::vector<Double_t> > binEntries_FB = std::vector<std::vector<Double_t> >(360, std::vector<Double_t>(180));

///////////////////////////////////////////////////////////////////////

void fcn_NS(int &npar, double *deriv, double &f, double par[], int flag) {
    
    Float_t lnL = 0.0;
    Float_t f1 = 0;
    Double_t tmp_entries = 0;
    Double_t tot_entries = 0;
    
    Double_t l_min = 0;
    Double_t l_max = 0;
    Double_t b_min = 0;
    Double_t b_max = 0;
    
    Double_t I = 0;
    Double_t D_NS = 0;
    Double_t D_EW = 0;
    Double_t D_FB = 0;
    
//    Double_t l = 0,b = 0;
    
    for(Int_t bX=1; bX<=360; bX++)
        for(Int_t bY=1; bY<=180; bY++)
            tot_entries += binEntries_NS[bX-1][bY-1];
            
    
    for(Int_t bX=1; bX<=360; bX++) {
        for(Int_t bY=1; bY<=180; bY++) {
            tmp_entries = binEntries_NS[bX-1][bY-1];
            if(tmp_entries==0)
                continue;
            
            /*
            l = -180 + .5*((bX-1)+bX);
            b = -90 + .5*((bY-1)+bY);
            */
             
            l_min = -180 + (bX-1);
            l_max = -180 + bX;
            
            b_min = -90 + (bY-1);
            b_max = -90 + bY;
            
            /*
            D_NS = TMath::Sqrt(3./(4*TMath::Pi()))*TMath::Cos(b*TMath::DegToRad());
            D_EW = TMath::Sqrt(2.)*TMath::Sqrt(3./(4*TMath::Pi()))*TMath::Sin(b*TMath::DegToRad())*TMath::Sin(l*TMath::DegToRad());
            D_FB = (-1)*TMath::Sqrt(2)*TMath::Sqrt(3./(4*TMath::Pi()))*TMath::Sin(b*TMath::DegToRad())*TMath::Cos(l*TMath::DegToRad());
            
            D_NS = TMath::Sin(l*TMath::DegToRad())*TMath::Cos(b*TMath::DegToRad());
            D_EW = TMath::Cos(l*TMath::DegToRad())*TMath::Cos(b*TMath::DegToRad());
            D_FB = TMath::Cos(b*TMath::DegToRad());
            */
            
            I = (l_max*TMath::DegToRad() - l_min*TMath::DegToRad())*(TMath::Cos(b_min*TMath::DegToRad())-TMath::Cos(b_max*TMath::DegToRad()));
            D_NS = .5*(TMath::Power(TMath::Sin(b_max*TMath::DegToRad()),2)-TMath::Power(TMath::Sin(b_min*TMath::DegToRad()),2))*(TMath::Cos(l_min*TMath::DegToRad())-TMath::Cos(l_max*TMath::DegToRad()));
            D_EW = .5*(TMath::Power(TMath::Sin(b_max*TMath::DegToRad()),2)-TMath::Power(TMath::Sin(b_min*TMath::DegToRad()),2))*(TMath::Sin(l_max*TMath::DegToRad())-TMath::Sin(l_min*TMath::DegToRad()));
            D_FB = .5*(TMath::Power(TMath::Sin(b_max*TMath::DegToRad()),2)-TMath::Power(TMath::Sin(b_min*TMath::DegToRad()),2))*(l_max*TMath::DegToRad()-l_min*TMath::DegToRad());
            
            f1 = abs(tot_entries*(par[0]*I + par[1]*D_NS + par[2]*D_EW + par[3]*D_FB));
            lnL += tmp_entries*log(f1);
        }
    }
    
    f=-lnL;

}

void fcn_EW(int &npar, double *deriv, double &f, double par[], int flag) {
    
    Float_t lnL = 0.0;
    Float_t f1 = 0;
    Double_t tmp_entries = 0;
    Double_t tot_entries = 0;
    
    Double_t l_min = 0;
    Double_t l_max = 0;
    Double_t b_min = 0;
    Double_t b_max = 0;
    
    Double_t I = 0;
    Double_t D_NS = 0;
    Double_t D_EW = 0;
    Double_t D_FB = 0;
    
 //   Double_t l = 0,b = 0;
    
    for(Int_t bX=1; bX<=360; bX++)
        for(Int_t bY=1; bY<=180; bY++)
            tot_entries += binEntries_EW[bX-1][bY-1];
    
    
    for(Int_t bX=1; bX<=360; bX++) {
        for(Int_t bY=1; bY<=180; bY++) {
            tmp_entries = binEntries_EW[bX-1][bY-1];
            if(tmp_entries==0)
                continue;
            
            /*
            l = -180 + .5*((bX-1)+bX);
            b = -90 + .5*((bY-1)+bY);
            */
             
            l_min = -180 + (bX-1);
            l_max = -180 + bX;
            
            b_min = -90 + (bY-1);
            b_max = -90 + bY;
            
            /*
             D_NS = TMath::Sqrt(3./(4*TMath::Pi()))*TMath::Cos(b*TMath::DegToRad());
             D_EW = TMath::Sqrt(2.)*TMath::Sqrt(3./(4*TMath::Pi()))*TMath::Sin(b*TMath::DegToRad())*TMath::Sin(l*TMath::DegToRad());
             D_FB = (-1)*TMath::Sqrt(2)*TMath::Sqrt(3./(4*TMath::Pi()))*TMath::Sin(b*TMath::DegToRad())*TMath::Cos(l*TMath::DegToRad());
             
             D_NS = TMath::Sin(l*TMath::DegToRad())*TMath::Cos(b*TMath::DegToRad());
             D_EW = TMath::Cos(l*TMath::DegToRad())*TMath::Cos(b*TMath::DegToRad());
             D_FB = TMath::Cos(b*TMath::DegToRad());
             */
            
            I = (l_max*TMath::DegToRad() - l_min*TMath::DegToRad())*(TMath::Cos(b_min*TMath::DegToRad())-TMath::Cos(b_max*TMath::DegToRad()));
            D_NS = .5*(TMath::Power(TMath::Sin(b_max*TMath::DegToRad()),2)-TMath::Power(TMath::Sin(b_min*TMath::DegToRad()),2))*(TMath::Cos(l_min*TMath::DegToRad())-TMath::Cos(l_max*TMath::DegToRad()));
            D_EW = .5*(TMath::Power(TMath::Sin(b_max*TMath::DegToRad()),2)-TMath::Power(TMath::Sin(b_min*TMath::DegToRad()),2))*(TMath::Sin(l_max*TMath::DegToRad())-TMath::Sin(l_min*TMath::DegToRad()));
            D_FB = .5*(TMath::Power(TMath::Sin(b_max*TMath::DegToRad()),2)-TMath::Power(TMath::Sin(b_min*TMath::DegToRad()),2))*(l_max*TMath::DegToRad()-l_min*TMath::DegToRad());
            
            f1 = abs(tot_entries*(par[0]*I + par[1]*D_NS + par[2]*D_EW + par[3]*D_FB));
            lnL += tmp_entries*log(f1);
        }
    }
    
    f=-lnL;
    
}

void fcn_FB(int &npar, double *deriv, double &f, double par[], int flag) {
    
    Float_t lnL = 0.0;
    Float_t f1 = 0;
    Double_t tmp_entries = 0;
    Double_t tot_entries = 0;
    
    Double_t l_min = 0;
    Double_t l_max = 0;
    Double_t b_min = 0;
    Double_t b_max = 0;
    
    Double_t I = 0;
    Double_t D_NS = 0;
    Double_t D_EW = 0;
    Double_t D_FB = 0;
    
  //  Double_t l = 0,b = 0;
    
    for(Int_t bX=1; bX<=360; bX++)
        for(Int_t bY=1; bY<=180; bY++)
            tot_entries += binEntries_FB[bX-1][bY-1];
    
    
    for(Int_t bX=1; bX<=360; bX++) {
        for(Int_t bY=1; bY<=180; bY++) {
            tmp_entries = binEntries_FB[bX-1][bY-1];
            if(tmp_entries==0)
                continue;
            
            /*
            l = -180 + .5*((bX-1)+bX);
            b = -90 + .5*((bY-1)+bY);
            */
             
            l_min = -180 + (bX-1);
            l_max = -180 + bX;
            
            b_min = -90 + (bY-1);
            b_max = -90 + bY;
            
            /*
             D_NS = TMath::Sqrt(3./(4*TMath::Pi()))*TMath::Cos(b*TMath::DegToRad());
             D_EW = TMath::Sqrt(2.)*TMath::Sqrt(3./(4*TMath::Pi()))*TMath::Sin(b*TMath::DegToRad())*TMath::Sin(l*TMath::DegToRad());
             D_FB = (-1)*TMath::Sqrt(2)*TMath::Sqrt(3./(4*TMath::Pi()))*TMath::Sin(b*TMath::DegToRad())*TMath::Cos(l*TMath::DegToRad());
             
             D_NS = TMath::Sin(l*TMath::DegToRad())*TMath::Cos(b*TMath::DegToRad());
             D_EW = TMath::Cos(l*TMath::DegToRad())*TMath::Cos(b*TMath::DegToRad());
             D_FB = TMath::Cos(b*TMath::DegToRad());
             */
            
            I = (l_max*TMath::DegToRad() - l_min*TMath::DegToRad())*(TMath::Cos(b_min*TMath::DegToRad())-TMath::Cos(b_max*TMath::DegToRad()));
            D_NS = .5*(TMath::Power(TMath::Sin(b_max*TMath::DegToRad()),2)-TMath::Power(TMath::Sin(b_min*TMath::DegToRad()),2))*(TMath::Cos(l_min*TMath::DegToRad())-TMath::Cos(l_max*TMath::DegToRad()));
            D_EW = .5*(TMath::Power(TMath::Sin(b_max*TMath::DegToRad()),2)-TMath::Power(TMath::Sin(b_min*TMath::DegToRad()),2))*(TMath::Sin(l_max*TMath::DegToRad())-TMath::Sin(l_min*TMath::DegToRad()));
            D_FB = .5*(TMath::Power(TMath::Sin(b_max*TMath::DegToRad()),2)-TMath::Power(TMath::Sin(b_min*TMath::DegToRad()),2))*(l_max*TMath::DegToRad()-l_min*TMath::DegToRad());
            
            f1 = abs(tot_entries*(par[0]*I + par[1]*D_NS + par[2]*D_EW + par[3]*D_FB));
            lnL += tmp_entries*log(f1);
        }
    }
    
    f=-lnL;
    
}


////////////////////////////////////////////////////// Not normalized maps


void fcn_NN_NS(int &npar, double *deriv, double &f, double par[], int flag) {
    
    Float_t lnL = 0.0;
    Float_t f1 = 0;
    Double_t tmp_entries = 0;
    Double_t tot_entries = 0;
    
    Double_t l_min = 0;
    Double_t l_max = 0;
    Double_t b_min = 0;
    Double_t b_max = 0;
    
    Double_t I = 0;
    Double_t D_NS = 0;
    Double_t D_EW = 0;
    Double_t D_FB = 0;
    
    //    Double_t l = 0,b = 0;
    
    for(Int_t bX=1; bX<=360; bX++)
        for(Int_t bY=1; bY<=180; bY++)
            tot_entries += binEntries_NS[bX-1][bY-1];
    
    
    for(Int_t bX=1; bX<=360; bX++) {
        for(Int_t bY=1; bY<=180; bY++) {
            tmp_entries = binEntries_NS[bX-1][bY-1];
            if(tmp_entries==0)
                continue;
            
            /*
             l = -180 + .5*((bX-1)+bX);
             b = -90 + .5*((bY-1)+bY);
             */
            
            l_min = -180 + (bX-1);
            l_max = -180 + bX;
            
            b_min = -90 + (bY-1);
            b_max = -90 + bY;
            
            /*
             D_NS = TMath::Sqrt(3./(4*TMath::Pi()))*TMath::Cos(b*TMath::DegToRad());
             D_EW = TMath::Sqrt(2.)*TMath::Sqrt(3./(4*TMath::Pi()))*TMath::Sin(b*TMath::DegToRad())*TMath::Sin(l*TMath::DegToRad());
             D_FB = (-1)*TMath::Sqrt(2)*TMath::Sqrt(3./(4*TMath::Pi()))*TMath::Sin(b*TMath::DegToRad())*TMath::Cos(l*TMath::DegToRad());
             
             D_NS = TMath::Sin(l*TMath::DegToRad())*TMath::Cos(b*TMath::DegToRad());
             D_EW = TMath::Cos(l*TMath::DegToRad())*TMath::Cos(b*TMath::DegToRad());
             D_FB = TMath::Cos(b*TMath::DegToRad());
             */
            
            I = (l_max*TMath::DegToRad() - l_min*TMath::DegToRad())*(TMath::Cos(b_min*TMath::DegToRad())-TMath::Cos(b_max*TMath::DegToRad()));
            D_NS = .5*(TMath::Power(TMath::Sin(b_max*TMath::DegToRad()),2)-TMath::Power(TMath::Sin(b_min*TMath::DegToRad()),2))*(TMath::Cos(l_min*TMath::DegToRad())-TMath::Cos(l_max*TMath::DegToRad()));
            D_EW = .5*(TMath::Power(TMath::Sin(b_max*TMath::DegToRad()),2)-TMath::Power(TMath::Sin(b_min*TMath::DegToRad()),2))*(TMath::Sin(l_max*TMath::DegToRad())-TMath::Sin(l_min*TMath::DegToRad()));
            D_FB = .5*(TMath::Power(TMath::Sin(b_max*TMath::DegToRad()),2)-TMath::Power(TMath::Sin(b_min*TMath::DegToRad()),2))*(l_max*TMath::DegToRad()-l_min*TMath::DegToRad());
            
            f1 = abs(tot_entries*(par[0]*I + par[1]*D_NS + par[2]*D_EW + par[3]*D_FB));
            lnL += tmp_entries*log(f1);
        }
    }
    
    f=-lnL;
    
}

void fcn_NN_EW(int &npar, double *deriv, double &f, double par[], int flag) {
    
    Float_t lnL = 0.0;
    Float_t f1 = 0;
    Double_t tmp_entries = 0;
    Double_t tot_entries = 0;
    
    Double_t l_min = 0;
    Double_t l_max = 0;
    Double_t b_min = 0;
    Double_t b_max = 0;
    
    Double_t I = 0;
    Double_t D_NS = 0;
    Double_t D_EW = 0;
    Double_t D_FB = 0;
    
    //   Double_t l = 0,b = 0;
    
    for(Int_t bX=1; bX<=360; bX++)
        for(Int_t bY=1; bY<=180; bY++)
            tot_entries += binEntries_EW[bX-1][bY-1];
    
    
    for(Int_t bX=1; bX<=360; bX++) {
        for(Int_t bY=1; bY<=180; bY++) {
            tmp_entries = binEntries_EW[bX-1][bY-1];
            if(tmp_entries==0)
                continue;
            
            /*
             l = -180 + .5*((bX-1)+bX);
             b = -90 + .5*((bY-1)+bY);
             */
            
            l_min = -180 + (bX-1);
            l_max = -180 + bX;
            
            b_min = -90 + (bY-1);
            b_max = -90 + bY;
            
            /*
             D_NS = TMath::Sqrt(3./(4*TMath::Pi()))*TMath::Cos(b*TMath::DegToRad());
             D_EW = TMath::Sqrt(2.)*TMath::Sqrt(3./(4*TMath::Pi()))*TMath::Sin(b*TMath::DegToRad())*TMath::Sin(l*TMath::DegToRad());
             D_FB = (-1)*TMath::Sqrt(2)*TMath::Sqrt(3./(4*TMath::Pi()))*TMath::Sin(b*TMath::DegToRad())*TMath::Cos(l*TMath::DegToRad());
             
             D_NS = TMath::Sin(l*TMath::DegToRad())*TMath::Cos(b*TMath::DegToRad());
             D_EW = TMath::Cos(l*TMath::DegToRad())*TMath::Cos(b*TMath::DegToRad());
             D_FB = TMath::Cos(b*TMath::DegToRad());
             */
            
            I = (l_max*TMath::DegToRad() - l_min*TMath::DegToRad())*(TMath::Cos(b_min*TMath::DegToRad())-TMath::Cos(b_max*TMath::DegToRad()));
            D_NS = .5*(TMath::Power(TMath::Sin(b_max*TMath::DegToRad()),2)-TMath::Power(TMath::Sin(b_min*TMath::DegToRad()),2))*(TMath::Cos(l_min*TMath::DegToRad())-TMath::Cos(l_max*TMath::DegToRad()));
            D_EW = .5*(TMath::Power(TMath::Sin(b_max*TMath::DegToRad()),2)-TMath::Power(TMath::Sin(b_min*TMath::DegToRad()),2))*(TMath::Sin(l_max*TMath::DegToRad())-TMath::Sin(l_min*TMath::DegToRad()));
            D_FB = .5*(TMath::Power(TMath::Sin(b_max*TMath::DegToRad()),2)-TMath::Power(TMath::Sin(b_min*TMath::DegToRad()),2))*(l_max*TMath::DegToRad()-l_min*TMath::DegToRad());
            
            f1 = abs(tot_entries*(par[0]*I + par[1]*D_NS + par[2]*D_EW + par[3]*D_FB));
            lnL += tmp_entries*log(f1);
        }
    }
    
    f=-lnL;
    
}

void fcn_NN_FB(int &npar, double *deriv, double &f, double par[], int flag) {
    
    Float_t lnL = 0.0;
    Float_t f1 = 0;
    Double_t tmp_entries = 0;
    Double_t tot_entries = 0;
    
    Double_t l_min = 0;
    Double_t l_max = 0;
    Double_t b_min = 0;
    Double_t b_max = 0;
    
    Double_t I = 0;
    Double_t D_NS = 0;
    Double_t D_EW = 0;
    Double_t D_FB = 0;
    
    //  Double_t l = 0,b = 0;
    
    for(Int_t bX=1; bX<=360; bX++)
        for(Int_t bY=1; bY<=180; bY++)
            tot_entries += binEntries_FB[bX-1][bY-1];
    
    
    for(Int_t bX=1; bX<=360; bX++) {
        for(Int_t bY=1; bY<=180; bY++) {
            tmp_entries = binEntries_FB[bX-1][bY-1];
            if(tmp_entries==0)
                continue;
            
            /*
             l = -180 + .5*((bX-1)+bX);
             b = -90 + .5*((bY-1)+bY);
             */
            
            l_min = -180 + (bX-1);
            l_max = -180 + bX;
            
            b_min = -90 + (bY-1);
            b_max = -90 + bY;
            
            /*
             D_NS = TMath::Sqrt(3./(4*TMath::Pi()))*TMath::Cos(b*TMath::DegToRad());
             D_EW = TMath::Sqrt(2.)*TMath::Sqrt(3./(4*TMath::Pi()))*TMath::Sin(b*TMath::DegToRad())*TMath::Sin(l*TMath::DegToRad());
             D_FB = (-1)*TMath::Sqrt(2)*TMath::Sqrt(3./(4*TMath::Pi()))*TMath::Sin(b*TMath::DegToRad())*TMath::Cos(l*TMath::DegToRad());
             
             D_NS = TMath::Sin(l*TMath::DegToRad())*TMath::Cos(b*TMath::DegToRad());
             D_EW = TMath::Cos(l*TMath::DegToRad())*TMath::Cos(b*TMath::DegToRad());
             D_FB = TMath::Cos(b*TMath::DegToRad());
             */
            
            I = (l_max*TMath::DegToRad() - l_min*TMath::DegToRad())*(TMath::Cos(b_min*TMath::DegToRad())-TMath::Cos(b_max*TMath::DegToRad()));
            D_NS = .5*(TMath::Power(TMath::Sin(b_max*TMath::DegToRad()),2)-TMath::Power(TMath::Sin(b_min*TMath::DegToRad()),2))*(TMath::Cos(l_min*TMath::DegToRad())-TMath::Cos(l_max*TMath::DegToRad()));
            D_EW = .5*(TMath::Power(TMath::Sin(b_max*TMath::DegToRad()),2)-TMath::Power(TMath::Sin(b_min*TMath::DegToRad()),2))*(TMath::Sin(l_max*TMath::DegToRad())-TMath::Sin(l_min*TMath::DegToRad()));
            D_FB = .5*(TMath::Power(TMath::Sin(b_max*TMath::DegToRad()),2)-TMath::Power(TMath::Sin(b_min*TMath::DegToRad()),2))*(l_max*TMath::DegToRad()-l_min*TMath::DegToRad());
            
            f1 = abs(tot_entries*(par[0]*I + par[1]*D_NS + par[2]*D_EW + par[3]*D_FB));
            lnL += tmp_entries*log(f1);
        }
    }
    
    f=-lnL;
    
}

/////////////////////////

void fcn_AllSky_NS(int &npar, double *deriv, double &f, double par[], int flag) {
    
    Float_t lnL = 0.0;
    Float_t f1 = 0;
    Double_t tmp_entries = 0;
    Double_t tot_entries = 0;
    
    Double_t l_min = 0;
    Double_t l_max = 0;
    Double_t b_min = 0;
    Double_t b_max = 0;
    
    Double_t I = 0;
    Double_t D_NS = 0;
    Double_t D_EW = 0;
    Double_t D_FB = 0;
    
    //    Double_t l = 0,b = 0;
    
    for(Int_t bX=1; bX<=360; bX++)
        for(Int_t bY=1; bY<=180; bY++)
            tot_entries += binEntries_NS[bX-1][bY-1];
    
    
    for(Int_t bX=1; bX<=360; bX++) {
        for(Int_t bY=1; bY<=180; bY++) {
            tmp_entries = binEntries_NS[bX-1][bY-1];
            if(tmp_entries==0)
                continue;
            
            /*
             l = -180 + .5*((bX-1)+bX);
             b = -90 + .5*((bY-1)+bY);
             */
            
            l_min = -180 + (bX-1);
            l_max = -180 + bX;
            
            b_min = -90 + (bY-1);
            b_max = -90 + bY;
            
            l_min+=180;
            l_max+=180;
            
            b_min+=90;
            b_max+=90;
            
            /*
             D_NS = TMath::Sqrt(3./(4*TMath::Pi()))*TMath::Cos(b*TMath::DegToRad());
             D_EW = TMath::Sqrt(2.)*TMath::Sqrt(3./(4*TMath::Pi()))*TMath::Sin(b*TMath::DegToRad())*TMath::Sin(l*TMath::DegToRad());
             D_FB = (-1)*TMath::Sqrt(2)*TMath::Sqrt(3./(4*TMath::Pi()))*TMath::Sin(b*TMath::DegToRad())*TMath::Cos(l*TMath::DegToRad());
             
             D_NS = TMath::Sin(l*TMath::DegToRad())*TMath::Cos(b*TMath::DegToRad());
             D_EW = TMath::Cos(l*TMath::DegToRad())*TMath::Cos(b*TMath::DegToRad());
             D_FB = TMath::Cos(b*TMath::DegToRad());
             */
            
            I = (l_max*TMath::DegToRad() - l_min*TMath::DegToRad())*(TMath::Cos(b_min*TMath::DegToRad())-TMath::Cos(b_max*TMath::DegToRad()));
            
            
            D_NS = TMath::Sqrt(3/(16*TMath::Pi()))*(l_max*TMath::DegToRad() - l_min*TMath::DegToRad())*(TMath::Power(TMath::Sin(b_max*TMath::DegToRad()),2)-TMath::Power(TMath::Sin(b_min*TMath::DegToRad()),2));
            
            
            D_EW = TMath::Sqrt(3/(8*TMath::Pi()))*(TMath::Sin(b_min*TMath::DegToRad())*TMath::Cos(b_min*TMath::DegToRad()) - TMath::Sin(b_max*TMath::DegToRad())*TMath::Cos(b_max*TMath::DegToRad()) + b_max*TMath::DegToRad()-b_min*TMath::DegToRad())*( TMath::Cos(l_min*TMath::DegToRad()) - TMath::Cos(l_max*TMath::DegToRad()));
            
            
            D_FB = -TMath::Sqrt(3/(8*TMath::Pi()))*(TMath::Sin(b_min*TMath::DegToRad())*TMath::Cos(b_min*TMath::DegToRad()) - TMath::Sin(b_max*TMath::DegToRad())*TMath::Cos(b_max*TMath::DegToRad()) + b_max*TMath::DegToRad()-b_min*TMath::DegToRad())*( TMath::Sin(l_max*TMath::DegToRad()) - TMath::Sin(l_min*TMath::DegToRad()));
            
            
            f1 = abs(tot_entries*(par[0]*I + par[1]*D_NS + par[2]*D_EW + par[3]*D_FB));
            lnL += tmp_entries*log(f1);
        }
    }
    
    f=-lnL;
    
}

void fcn_AllSky_EW(int &npar, double *deriv, double &f, double par[], int flag) {
    
    Float_t lnL = 0.0;
    Float_t f1 = 0;
    Double_t tmp_entries = 0;
    Double_t tot_entries = 0;
    
    Double_t l_min = 0;
    Double_t l_max = 0;
    Double_t b_min = 0;
    Double_t b_max = 0;
    
    Double_t I = 0;
    Double_t D_NS = 0;
    Double_t D_EW = 0;
    Double_t D_FB = 0;
    
    //   Double_t l = 0,b = 0;
    
    for(Int_t bX=1; bX<=360; bX++)
        for(Int_t bY=1; bY<=180; bY++)
            tot_entries += binEntries_EW[bX-1][bY-1];
    
    
    for(Int_t bX=1; bX<=360; bX++) {
        for(Int_t bY=1; bY<=180; bY++) {
            tmp_entries = binEntries_EW[bX-1][bY-1];
            if(tmp_entries==0)
                continue;
            
            l_min = -180 + (bX-1);
            l_max = -180 + bX;
            
            b_min = -90 + (bY-1);
            b_max = -90 + bY;
            
            l_min+=180;
            l_max+=180;
            
            b_min+=90;
            b_max+=90;
            
            /*
             D_NS = TMath::Sqrt(3./(4*TMath::Pi()))*TMath::Cos(b*TMath::DegToRad());
             D_EW = TMath::Sqrt(2.)*TMath::Sqrt(3./(4*TMath::Pi()))*TMath::Sin(b*TMath::DegToRad())*TMath::Sin(l*TMath::DegToRad());
             D_FB = (-1)*TMath::Sqrt(2)*TMath::Sqrt(3./(4*TMath::Pi()))*TMath::Sin(b*TMath::DegToRad())*TMath::Cos(l*TMath::DegToRad());
             
             D_NS = TMath::Sin(l*TMath::DegToRad())*TMath::Cos(b*TMath::DegToRad());
             D_EW = TMath::Cos(l*TMath::DegToRad())*TMath::Cos(b*TMath::DegToRad());
             D_FB = TMath::Cos(b*TMath::DegToRad());
             */
            
            I = (l_max*TMath::DegToRad() - l_min*TMath::DegToRad())*(TMath::Cos(b_min*TMath::DegToRad())-TMath::Cos(b_max*TMath::DegToRad()));
            
            
            D_NS = TMath::Sqrt(3/(16*TMath::Pi()))*(l_max*TMath::DegToRad() - l_min*TMath::DegToRad())*(TMath::Power(TMath::Sin(b_max*TMath::DegToRad()),2)-TMath::Power(TMath::Sin(b_min*TMath::DegToRad()),2));
            
        
            D_EW = TMath::Sqrt(3/(8*TMath::Pi()))*(TMath::Sin(b_min*TMath::DegToRad())*TMath::Cos(b_min*TMath::DegToRad()) - TMath::Sin(b_max*TMath::DegToRad())*TMath::Cos(b_max*TMath::DegToRad()) + b_max*TMath::DegToRad()-b_min*TMath::DegToRad())*( TMath::Cos(l_min*TMath::DegToRad()) - TMath::Cos(l_max*TMath::DegToRad()));
            
            
            D_FB = -TMath::Sqrt(3/(8*TMath::Pi()))*(TMath::Sin(b_min*TMath::DegToRad())*TMath::Cos(b_min*TMath::DegToRad()) - TMath::Sin(b_max*TMath::DegToRad())*TMath::Cos(b_max*TMath::DegToRad()) + b_max*TMath::DegToRad()-b_min*TMath::DegToRad())*( TMath::Sin(l_max*TMath::DegToRad()) - TMath::Sin(l_min*TMath::DegToRad()));
            
            f1 = abs(tot_entries*(par[0]*I + par[1]*D_NS + par[2]*D_EW + par[3]*D_FB));
            lnL += tmp_entries*log(f1);
        }
    }
    
    f=-lnL;
    
}

void fcn_AllSky_FB(int &npar, double *deriv, double &f, double par[], int flag) {
    
    Float_t lnL = 0.0;
    Float_t f1 = 0;
    Double_t tmp_entries = 0;
    Double_t tot_entries = 0;
    
    Double_t l_min = 0;
    Double_t l_max = 0;
    Double_t b_min = 0;
    Double_t b_max = 0;
    
    Double_t I = 0;
    Double_t D_NS = 0;
    Double_t D_EW = 0;
    Double_t D_FB = 0;
    
    //  Double_t l = 0,b = 0;
    
    for(Int_t bX=1; bX<=360; bX++)
        for(Int_t bY=1; bY<=180; bY++)
            tot_entries += binEntries_FB[bX-1][bY-1];
    
    for(Int_t bX=1; bX<=360; bX++) {
        for(Int_t bY=1; bY<=180; bY++) {
            tmp_entries = binEntries_FB[bX-1][bY-1];
            if(tmp_entries==0)
                continue;
            
            /*
             l = -180 + .5*((bX-1)+bX);
             b = -90 + .5*((bY-1)+bY);
             */
            
            l_min = -180 + (bX-1);
            l_max = -180 + bX;
            
            b_min = -90 + (bY-1);
            b_max = -90 + bY;
            
            l_min+=180;
            l_max+=180;
            
            b_min+=90;
            b_max+=90;
            
            /*
             D_NS = TMath::Sqrt(3./(4*TMath::Pi()))*TMath::Cos(b*TMath::DegToRad());
             D_EW = TMath::Sqrt(2.)*TMath::Sqrt(3./(4*TMath::Pi()))*TMath::Sin(b*TMath::DegToRad())*TMath::Sin(l*TMath::DegToRad());
             D_FB = (-1)*TMath::Sqrt(2)*TMath::Sqrt(3./(4*TMath::Pi()))*TMath::Sin(b*TMath::DegToRad())*TMath::Cos(l*TMath::DegToRad());
             
             D_NS = TMath::Sin(l*TMath::DegToRad())*TMath::Cos(b*TMath::DegToRad());
             D_EW = TMath::Cos(l*TMath::DegToRad())*TMath::Cos(b*TMath::DegToRad());
             D_FB = TMath::Cos(b*TMath::DegToRad());
             */
            
            I = (l_max*TMath::DegToRad() - l_min*TMath::DegToRad())*(TMath::Cos(b_min*TMath::DegToRad())-TMath::Cos(b_max*TMath::DegToRad()));
            
            
            D_NS = TMath::Sqrt(3/(16*TMath::Pi()))*(l_max*TMath::DegToRad() - l_min*TMath::DegToRad())*(TMath::Power(TMath::Sin(b_max*TMath::DegToRad()),2)-TMath::Power(TMath::Sin(b_min*TMath::DegToRad()),2));
            
            
            D_EW = TMath::Sqrt(3/(8*TMath::Pi()))*(TMath::Sin(b_min*TMath::DegToRad())*TMath::Cos(b_min*TMath::DegToRad()) - TMath::Sin(b_max*TMath::DegToRad())*TMath::Cos(b_max*TMath::DegToRad()) + b_max*TMath::DegToRad()-b_min*TMath::DegToRad())*( TMath::Cos(l_min*TMath::DegToRad()) - TMath::Cos(l_max*TMath::DegToRad()));
            
            
            D_FB = -TMath::Sqrt(3/(8*TMath::Pi()))*(TMath::Sin(b_min*TMath::DegToRad())*TMath::Cos(b_min*TMath::DegToRad()) - TMath::Sin(b_max*TMath::DegToRad())*TMath::Cos(b_max*TMath::DegToRad()) + b_max*TMath::DegToRad()-b_min*TMath::DegToRad())*( TMath::Sin(l_max*TMath::DegToRad()) - TMath::Sin(l_min*TMath::DegToRad()));
            
            
            f1 = abs(tot_entries*(par[0]*I + par[1]*D_NS + par[2]*D_EW + par[3]*D_FB));
            lnL += tmp_entries*log(f1);
        }
    }
    
    f=-lnL;
    
}






/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void fit_map(TH2D &Isolate_NS,TH2D &Isolate_EW,TH2D &Isolate_FB) {
    
    for(Int_t bX=1; bX<=360; bX++) {
        for(Int_t bY=1; bY<=180; bY++) {
            binEntries_NS[bX-1][bY-1] = (Double_t)Isolate_NS.GetBinContent(bX,bY);
            binEntries_EW[bX-1][bY-1] = (Double_t)Isolate_EW.GetBinContent(bX,bY);
            binEntries_FB[bX-1][bY-1] = (Double_t)Isolate_FB.GetBinContent(bX,bY);
        }
    }
    
    
    //////////////////////////////////// TMinuit Fit
    
    //////////////////// General variables ////////////////////
    
    int err = 0;
    double set_str_list[1];
    set_str_list[0] = 2;
    
    std::string parName[4];
    double par[4],stepSize[4],minVal[4],maxVal[4];
    
    ///////////////////////////////////////////////////////////
    
    
    
    ///////////////////// NS Minimization /////////////////////
    
    std::cout<<"\n\n \\ -------------------------  Fitting NS Dipole  ------------------------- \\ \n\n";
    
    TMinuit minimization_NS(4);
    minimization_NS.SetFCN(fcn_NS);
    
    par[0] = 0.1;
    stepSize[0] = 0.001;
    minVal[0] = 0;
    maxVal[0] = 1;
    parName[0] = "I";
    
    par[1] = 0.1;
    stepSize[1] = 0.001;
    minVal[1] = 0;
    maxVal[1] = 1;
    parName[1] = "D_NS";
    
    par[2] = 0.1;
    stepSize[2] = 0.001;
    minVal[2] = 0;
    maxVal[2] = 1;
    parName[2] = "D_EW";
    
    par[3] = 0.1;
    stepSize[3] = 0.001;
    minVal[3] = 0;
    maxVal[3] = 1;
    parName[3] = "D_FB";
    
    minimization_NS.DefineParameter(0,parName[0].c_str(),par[0],stepSize[0],minVal[0],maxVal[0]);
    minimization_NS.DefineParameter(1,parName[1].c_str(),par[1],stepSize[1],minVal[1],maxVal[1]);
    minimization_NS.DefineParameter(2,parName[2].c_str(),par[2],stepSize[2],minVal[2],maxVal[2]);
    minimization_NS.DefineParameter(3,parName[3].c_str(),par[3],stepSize[3],minVal[3],maxVal[3]);
    
    minimization_NS.SetErrorDef(0.5);
    minimization_NS.mnexcm("SET STR",set_str_list,1,err);

    minimization_NS.Migrad();
    
    
    ///////////////////// EW Minimization /////////////////////
    
    std::cout<<"\n\n \\ -------------------------  Fitting EW Dipole  ------------------------- \\ \n\n";
    
    TMinuit minimization_EW(4);
    minimization_EW.SetFCN(fcn_EW);
    
    par[0] = 0.1;
    stepSize[0] = 0.001;
    minVal[0] = 0;
    maxVal[0] = 1;
    parName[0] = "I";
    
    par[1] = 0.1;
    stepSize[1] = 0.001;
    minVal[1] = 0;
    maxVal[1] = 1;
    parName[1] = "D_NS";
    
    par[2] = 0.1;
    stepSize[2] = 0.001;
    minVal[2] = 0;
    maxVal[2] = 1;
    parName[2] = "D_EW";
    
    par[3] = 0.1;
    stepSize[3] = 0.001;
    minVal[3] = 0;
    maxVal[3] = 1;
    parName[3] = "D_FB";
    
    minimization_EW.DefineParameter(0,parName[0].c_str(),par[0],stepSize[0],minVal[0],maxVal[0]);
    minimization_EW.DefineParameter(1,parName[1].c_str(),par[1],stepSize[1],minVal[1],maxVal[1]);
    minimization_EW.DefineParameter(2,parName[2].c_str(),par[2],stepSize[2],minVal[2],maxVal[2]);
    minimization_EW.DefineParameter(3,parName[3].c_str(),par[3],stepSize[3],minVal[3],maxVal[3]);
    
    minimization_EW.SetErrorDef(0.5);
    minimization_EW.mnexcm("SET STR",set_str_list,1,err);
    
    minimization_EW.Migrad();
    
    
    ///////////////////// FB Minimization /////////////////////
    
    std::cout<<"\n\n \\ -------------------------  Fitting FB Dipole  ------------------------- \\ \n\n";
    
    TMinuit minimization_FB(4);
    minimization_FB.SetFCN(fcn_FB);
    
    par[0] = 0.1;
    stepSize[0] = 0.001;
    minVal[0] = 0;
    maxVal[0] = 1;
    parName[0] = "I";
    
    par[1] = 0.1;
    stepSize[1] = 0.001;
    minVal[1] = 0;
    maxVal[1] = 1;
    parName[1] = "D_NS";
    
    par[2] = 0.1;
    stepSize[2] = 0.001;
    minVal[2] = 0;
    maxVal[2] = 1;
    parName[2] = "D_EW";
    
    par[3] = 0.1;
    stepSize[3] = 0.001;
    minVal[3] = 0;
    maxVal[3] = 1;
    parName[3] = "D_FB";
    
    minimization_FB.DefineParameter(0,parName[0].c_str(),par[0],stepSize[0],minVal[0],maxVal[0]);
    minimization_FB.DefineParameter(1,parName[1].c_str(),par[1],stepSize[1],minVal[1],maxVal[1]);
    minimization_FB.DefineParameter(2,parName[2].c_str(),par[2],stepSize[2],minVal[2],maxVal[2]);
    minimization_FB.DefineParameter(3,parName[3].c_str(),par[3],stepSize[3],minVal[3],maxVal[3]);
    
    minimization_FB.SetErrorDef(0.5);
    minimization_FB.mnexcm("SET STR",set_str_list,1,err);
    
    minimization_FB.Migrad();
    

}


void fit_NN_map(TH2D &Isolate_NN_NS,TH2D &Isolate_NN_EW,TH2D &Isolate_NN_FB) {
    
    for(Int_t bX=1; bX<=360; bX++) {
        for(Int_t bY=1; bY<=180; bY++) {
            binEntries_NS[bX-1][bY-1] = (Double_t)Isolate_NN_NS.GetBinContent(bX,bY);
            binEntries_EW[bX-1][bY-1] = (Double_t)Isolate_NN_EW.GetBinContent(bX,bY);
            binEntries_FB[bX-1][bY-1] = (Double_t)Isolate_NN_FB.GetBinContent(bX,bY);
        }
    }
    
    
    //////////////////////////////////// TMinuit Fit
    
    //////////////////// General variables ////////////////////
    
    int err = 0;
    double set_str_list[1];
    set_str_list[0] = 2;
    
    std::string parName[4];
    double par[4],stepSize[4],minVal[4],maxVal[4];
    
    ///////////////////////////////////////////////////////////
    
    ///////////////////// NS Minimization /////////////////////
    
    std::cout<<"\n\n \\ -------------------------  Fitting NS Dipole  ------------------------- \\ \n\n";
    
    TMinuit minimization_NS(4);
    minimization_NS.SetFCN(fcn_NN_NS);
    
    par[0] = 0.5;
    stepSize[0] = 0.001;
    minVal[0] = 0;
    maxVal[0] = 1;
    parName[0] = "I";
    
    par[1] = 0.5;
    stepSize[1] = 0.001;
    minVal[1] = 0;
    maxVal[1] = 1;
    parName[1] = "D_NS";
    
    par[2] = 0.5;
    stepSize[2] = 0.001;
    minVal[2] = 0;
    maxVal[2] = 1;
    parName[2] = "D_EW";
    
    par[3] = 0.5;
    stepSize[3] = 0.001;
    minVal[3] = 0;
    maxVal[3] = 1;
    parName[3] = "D_FB";
    
    minimization_NS.DefineParameter(0,parName[0].c_str(),par[0],stepSize[0],minVal[0],maxVal[0]);
    minimization_NS.DefineParameter(1,parName[1].c_str(),par[1],stepSize[1],minVal[1],maxVal[1]);
    minimization_NS.DefineParameter(2,parName[2].c_str(),par[2],stepSize[2],minVal[2],maxVal[2]);
    minimization_NS.DefineParameter(3,parName[3].c_str(),par[3],stepSize[3],minVal[3],maxVal[3]);
    
    minimization_NS.SetErrorDef(0.5);
    minimization_NS.mnexcm("SET STR",set_str_list,1,err);
    
    minimization_NS.Migrad();
    
    
    ///////////////////// EW Minimization /////////////////////
    
    std::cout<<"\n\n \\ -------------------------  Fitting EW Dipole  ------------------------- \\ \n\n";
    
    TMinuit minimization_EW(4);
    minimization_EW.SetFCN(fcn_NN_EW);
    
    par[0] = 0.5;
    stepSize[0] = 0.001;
    minVal[0] = 0;
    maxVal[0] = 1;
    parName[0] = "I";
    
    par[1] = 0.5;
    stepSize[1] = 0.001;
    minVal[1] = 0;
    maxVal[1] = 1;
    parName[1] = "D_NS";
    
    par[2] = 0.5;
    stepSize[2] = 0.001;
    minVal[2] = 0;
    maxVal[2] = 1;
    parName[2] = "D_EW";
    
    par[3] = 0.5;
    stepSize[3] = 0.001;
    minVal[3] = 0;
    maxVal[3] = 1;
    parName[3] = "D_FB";
    
    minimization_EW.DefineParameter(0,parName[0].c_str(),par[0],stepSize[0],minVal[0],maxVal[0]);
    minimization_EW.DefineParameter(1,parName[1].c_str(),par[1],stepSize[1],minVal[1],maxVal[1]);
    minimization_EW.DefineParameter(2,parName[2].c_str(),par[2],stepSize[2],minVal[2],maxVal[2]);
    minimization_EW.DefineParameter(3,parName[3].c_str(),par[3],stepSize[3],minVal[3],maxVal[3]);
    
    minimization_EW.SetErrorDef(0.5);
    minimization_EW.mnexcm("SET STR",set_str_list,1,err);
    
    minimization_EW.Migrad();
    
    
    ///////////////////// FB Minimization /////////////////////
    
    std::cout<<"\n\n \\ -------------------------  Fitting FB Dipole  ------------------------- \\ \n\n";
    
    TMinuit minimization_FB(4);
    minimization_FB.SetFCN(fcn_NN_FB);
    
    par[0] = 0.5;
    stepSize[0] = 0.001;
    minVal[0] = 0;
    maxVal[0] = 1;
    parName[0] = "I";
    
    par[1] = 0.5;
    stepSize[1] = 0.001;
    minVal[1] = 0;
    maxVal[1] = 1;
    parName[1] = "D_NS";
    
    par[2] = 0.5;
    stepSize[2] = 0.001;
    minVal[2] = 0;
    maxVal[2] = 1;
    parName[2] = "D_EW";
    
    par[3] = 0.5;
    stepSize[3] = 0.001;
    minVal[3] = 0;
    maxVal[3] = 1;
    parName[3] = "D_FB";
    
    minimization_FB.DefineParameter(0,parName[0].c_str(),par[0],stepSize[0],minVal[0],maxVal[0]);
    minimization_FB.DefineParameter(1,parName[1].c_str(),par[1],stepSize[1],minVal[1],maxVal[1]);
    minimization_FB.DefineParameter(2,parName[2].c_str(),par[2],stepSize[2],minVal[2],maxVal[2]);
    minimization_FB.DefineParameter(3,parName[3].c_str(),par[3],stepSize[3],minVal[3],maxVal[3]);
    
    minimization_FB.SetErrorDef(0.5);
    minimization_FB.mnexcm("SET STR",set_str_list,1,err);
    
    minimization_FB.Migrad();
    
}

void AllSky_MapFit(TH2D &AllSky_Isolate_NS,TH2D &AllSky_Isolate_EW,TH2D &AllSky_Isolate_FB) {
    
    for(Int_t bX=1; bX<=360; bX++) {
        for(Int_t bY=1; bY<=180; bY++) {
            binEntries_NS[bX-1][bY-1] = (Double_t)AllSky_Isolate_NS.GetBinContent(bX,bY);
            binEntries_EW[bX-1][bY-1] = (Double_t)AllSky_Isolate_EW.GetBinContent(bX,bY);
            binEntries_FB[bX-1][bY-1] = (Double_t)AllSky_Isolate_FB.GetBinContent(bX,bY);
        }
    }
    
    Double_t tot_entries=0;
    
    for(Int_t bX=1; bX<=360; bX++)
        for(Int_t bY=1; bY<=180; bY++)
            tot_entries += binEntries_FB[bX-1][bY-1];
    
    std::cout<<"\n\n\n\nTot entries: "<<tot_entries<<"\n\n\n\n";
    
    //////////////////////////////////// TMinuit Fit
    
    //////////////////// General variables ////////////////////
    
    int err = 0;
    double set_str_list[1];
    set_str_list[0] = 2;
    
    std::string parName[4];
    double par[4],stepSize[4],minVal[4],maxVal[4];
    
    ///////////////////////////////////////////////////////////
    
    ///////////////////// NS Minimization /////////////////////
    
    std::cout<<"\n\n \\ -------------------------  Fitting NS Dipole  ------------------------- \\ \n\n";
    
    TMinuit minimization_NS(4);
    minimization_NS.SetFCN(fcn_AllSky_NS);
    
    par[0] = 0.5;
    stepSize[0] = 0.001;
    minVal[0] = 0;
    maxVal[0] = 1;
    parName[0] = "I";
    
    par[1] = 0.5;
    stepSize[1] = 0.001;
    minVal[1] = 0;
    maxVal[1] = 1;
    parName[1] = "D_NS";
    
    par[2] = 0.5;
    stepSize[2] = 0.001;
    minVal[2] = 0;
    maxVal[2] = 1;
    parName[2] = "D_EW";
    
    par[3] = 0.5;
    stepSize[3] = 0.001;
    minVal[3] = 0;
    maxVal[3] = 1;
    parName[3] = "D_FB";
    
    minimization_NS.DefineParameter(0,parName[0].c_str(),par[0],stepSize[0],minVal[0],maxVal[0]);
    minimization_NS.DefineParameter(1,parName[1].c_str(),par[1],stepSize[1],minVal[1],maxVal[1]);
    minimization_NS.DefineParameter(2,parName[2].c_str(),par[2],stepSize[2],minVal[2],maxVal[2]);
    minimization_NS.DefineParameter(3,parName[3].c_str(),par[3],stepSize[3],minVal[3],maxVal[3]);
    
    minimization_NS.SetErrorDef(0.5);
    minimization_NS.mnexcm("SET STR",set_str_list,1,err);
    
    minimization_NS.Migrad();
    
    
    ///////////////////// EW Minimization /////////////////////
    
    std::cout<<"\n\n \\ -------------------------  Fitting EW Dipole  ------------------------- \\ \n\n";
    
    TMinuit minimization_EW(4);
    minimization_EW.SetFCN(fcn_NN_EW);
    
    par[0] = 0.5;
    stepSize[0] = 0.001;
    minVal[0] = 0;
    maxVal[0] = 1;
    parName[0] = "I";
    
    par[1] = 0.5;
    stepSize[1] = 0.001;
    minVal[1] = 0;
    maxVal[1] = 1;
    parName[1] = "D_NS";
    
    par[2] = 0.5;
    stepSize[2] = 0.001;
    minVal[2] = 0;
    maxVal[2] = 1;
    parName[2] = "D_EW";
    
    par[3] = 0.5;
    stepSize[3] = 0.001;
    minVal[3] = 0;
    maxVal[3] = 1;
    parName[3] = "D_FB";
    
    minimization_EW.DefineParameter(0,parName[0].c_str(),par[0],stepSize[0],minVal[0],maxVal[0]);
    minimization_EW.DefineParameter(1,parName[1].c_str(),par[1],stepSize[1],minVal[1],maxVal[1]);
    minimization_EW.DefineParameter(2,parName[2].c_str(),par[2],stepSize[2],minVal[2],maxVal[2]);
    minimization_EW.DefineParameter(3,parName[3].c_str(),par[3],stepSize[3],minVal[3],maxVal[3]);
    
    minimization_EW.SetErrorDef(0.5);
    minimization_EW.mnexcm("SET STR",set_str_list,1,err);
    
    minimization_EW.Migrad();
    
    
    ///////////////////// FB Minimization /////////////////////
    
    std::cout<<"\n\n \\ -------------------------  Fitting FB Dipole  ------------------------- \\ \n\n";
    
    TMinuit minimization_FB(4);
    minimization_FB.SetFCN(fcn_NN_FB);
    
    par[0] = 0.5;
    stepSize[0] = 0.001;
    minVal[0] = 0;
    maxVal[0] = 1;
    parName[0] = "I";
    
    par[1] = 0.5;
    stepSize[1] = 0.001;
    minVal[1] = 0;
    maxVal[1] = 1;
    parName[1] = "D_NS";
    
    par[2] = 0.5;
    stepSize[2] = 0.001;
    minVal[2] = 0;
    maxVal[2] = 1;
    parName[2] = "D_EW";
    
    par[3] = 0.5;
    stepSize[3] = 0.001;
    minVal[3] = 0;
    maxVal[3] = 1;
    parName[3] = "D_FB";
    
    minimization_FB.DefineParameter(0,parName[0].c_str(),par[0],stepSize[0],minVal[0],maxVal[0]);
    minimization_FB.DefineParameter(1,parName[1].c_str(),par[1],stepSize[1],minVal[1],maxVal[1]);
    minimization_FB.DefineParameter(2,parName[2].c_str(),par[2],stepSize[2],minVal[2],maxVal[2]);
    minimization_FB.DefineParameter(3,parName[3].c_str(),par[3],stepSize[3],minVal[3],maxVal[3]);
    
    minimization_FB.SetErrorDef(0.5);
    minimization_FB.mnexcm("SET STR",set_str_list,1,err);
    
    minimization_FB.Migrad();
    
    
}

void fcn_Chi2_NS(int &npar, double *deriv, double &f, double par[], int flag) {
    
    Float_t f1 = 0;
    Double_t tmp_entries = 0;
    Double_t tot_entries = 0;
    
    Double_t l_min = 0;
    Double_t l_max = 0;
    Double_t b_min = 0;
    Double_t b_max = 0;
    
    Double_t I = 0;
    Double_t D_NS = 0;
    Double_t D_EW = 0;
    Double_t D_FB = 0;
    
    //  Double_t l = 0,b = 0;
    
    for(Int_t bX=1; bX<=360; bX++)
        for(Int_t bY=1; bY<=180; bY++)
            tot_entries += binEntries_NS[bX-1][bY-1];
    
    for(Int_t bX=1; bX<=360; bX++) {
        for(Int_t bY=1; bY<=180; bY++) {
            tmp_entries = binEntries_NS[bX-1][bY-1];
            if(tmp_entries==0)
                continue;
            
            /*
             l = -180 + .5*((bX-1)+bX);
             b = -90 + .5*((bY-1)+bY);
             */
            
            l_min = -180 + (bX-1);
            l_max = -180 + bX;
            
            b_min = -90 + (bY-1);
            b_max = -90 + bY;
            
            l_min+=180;
            l_max+=180;
            
            b_min+=90;
            b_max+=90;
            
            I = (l_max*TMath::DegToRad() - l_min*TMath::DegToRad())*(TMath::Cos(b_min*TMath::DegToRad())-TMath::Cos(b_max*TMath::DegToRad()));
            
            
            D_NS = TMath::Sqrt(3/(16*TMath::Pi()))*(l_max*TMath::DegToRad() - l_min*TMath::DegToRad())*(TMath::Power(TMath::Sin(b_max*TMath::DegToRad()),2)-TMath::Power(TMath::Sin(b_min*TMath::DegToRad()),2));
            
            
            D_EW = TMath::Sqrt(3/(8*TMath::Pi()))*(TMath::Sin(b_min*TMath::DegToRad())*TMath::Cos(b_min*TMath::DegToRad()) - TMath::Sin(b_max*TMath::DegToRad())*TMath::Cos(b_max*TMath::DegToRad()) + b_max*TMath::DegToRad()-b_min*TMath::DegToRad())*( TMath::Cos(l_min*TMath::DegToRad()) - TMath::Cos(l_max*TMath::DegToRad()));
            
            
            D_FB = -TMath::Sqrt(3/(8*TMath::Pi()))*(TMath::Sin(b_min*TMath::DegToRad())*TMath::Cos(b_min*TMath::DegToRad()) - TMath::Sin(b_max*TMath::DegToRad())*TMath::Cos(b_max*TMath::DegToRad()) + b_max*TMath::DegToRad()-b_min*TMath::DegToRad())*( TMath::Sin(l_max*TMath::DegToRad()) - TMath::Sin(l_min*TMath::DegToRad()));
            
            
            f1 += TMath::Power( tmp_entries -  tot_entries*(par[0]*I + par[1]*D_NS + par[2]*D_EW + par[3]*D_FB),2)/tot_entries*(par[0]*I + par[1]*D_NS + par[2]*D_EW + par[3]*D_FB);
            
        }
    }
    
    f=f1;
    
}





void fcn_Chi2_EW(int &npar, double *deriv, double &f, double par[], int flag) {
    
    Float_t f1 = 0;
    Double_t tmp_entries = 0;
    Double_t tot_entries = 0;
    
    Double_t l_min = 0;
    Double_t l_max = 0;
    Double_t b_min = 0;
    Double_t b_max = 0;
    
    Double_t I = 0;
    Double_t D_NS = 0;
    Double_t D_EW = 0;
    Double_t D_FB = 0;
    
    //  Double_t l = 0,b = 0;
    
    for(Int_t bX=1; bX<=360; bX++)
        for(Int_t bY=1; bY<=180; bY++)
            tot_entries += binEntries_EW[bX-1][bY-1];
    
    for(Int_t bX=1; bX<=360; bX++) {
        for(Int_t bY=1; bY<=180; bY++) {
            tmp_entries = binEntries_EW[bX-1][bY-1];
            if(tmp_entries==0)
                continue;
            
            /*
             l = -180 + .5*((bX-1)+bX);
             b = -90 + .5*((bY-1)+bY);
             */
            
            l_min = -180 + (bX-1);
            l_max = -180 + bX;
            
            b_min = -90 + (bY-1);
            b_max = -90 + bY;
            
            l_min+=180;
            l_max+=180;
            
            b_min+=90;
            b_max+=90;
            
            I = (l_max*TMath::DegToRad() - l_min*TMath::DegToRad())*(TMath::Cos(b_min*TMath::DegToRad())-TMath::Cos(b_max*TMath::DegToRad()));
            
            
            D_NS = TMath::Sqrt(3/(16*TMath::Pi()))*(l_max*TMath::DegToRad() - l_min*TMath::DegToRad())*(TMath::Power(TMath::Sin(b_max*TMath::DegToRad()),2)-TMath::Power(TMath::Sin(b_min*TMath::DegToRad()),2));
            
            
            D_EW = TMath::Sqrt(3/(8*TMath::Pi()))*(TMath::Sin(b_min*TMath::DegToRad())*TMath::Cos(b_min*TMath::DegToRad()) - TMath::Sin(b_max*TMath::DegToRad())*TMath::Cos(b_max*TMath::DegToRad()) + b_max*TMath::DegToRad()-b_min*TMath::DegToRad())*( TMath::Cos(l_min*TMath::DegToRad()) - TMath::Cos(l_max*TMath::DegToRad()));
            
            
            D_FB = -TMath::Sqrt(3/(8*TMath::Pi()))*(TMath::Sin(b_min*TMath::DegToRad())*TMath::Cos(b_min*TMath::DegToRad()) - TMath::Sin(b_max*TMath::DegToRad())*TMath::Cos(b_max*TMath::DegToRad()) + b_max*TMath::DegToRad()-b_min*TMath::DegToRad())*( TMath::Sin(l_max*TMath::DegToRad()) - TMath::Sin(l_min*TMath::DegToRad()));
            
            
            f1 += TMath::Power( tmp_entries -  tot_entries*(par[0]*I + par[1]*D_NS + par[2]*D_EW + par[3]*D_FB),2)/tot_entries*(par[0]*I + par[1]*D_NS + par[2]*D_EW + par[3]*D_FB);
            
        }
    }
    
    f=f1;
    
}




void fcn_Chi2_FB(int &npar, double *deriv, double &f, double par[], int flag) {
    
    Float_t f1 = 0;
    Double_t tmp_entries = 0;
    Double_t tot_entries = 0;
    
    Double_t l_min = 0;
    Double_t l_max = 0;
    Double_t b_min = 0;
    Double_t b_max = 0;
    
    Double_t I = 0;
    Double_t D_NS = 0;
    Double_t D_EW = 0;
    Double_t D_FB = 0;
    
    //  Double_t l = 0,b = 0;
    
    for(Int_t bX=1; bX<=360; bX++)
        for(Int_t bY=1; bY<=180; bY++)
            tot_entries += binEntries_FB[bX-1][bY-1];
    
    for(Int_t bX=1; bX<=360; bX++) {
        for(Int_t bY=1; bY<=180; bY++) {
            tmp_entries = binEntries_FB[bX-1][bY-1];
            if(tmp_entries==0)
                continue;
            
            /*
             l = -180 + .5*((bX-1)+bX);
             b = -90 + .5*((bY-1)+bY);
             */
            
            l_min = -180 + (bX-1);
            l_max = -180 + bX;
            
            b_min = -90 + (bY-1);
            b_max = -90 + bY;
            
            l_min+=180;
            l_max+=180;
            
            b_min+=90;
            b_max+=90;
            
            I = (l_max*TMath::DegToRad() - l_min*TMath::DegToRad())*(TMath::Cos(b_min*TMath::DegToRad())-TMath::Cos(b_max*TMath::DegToRad()));
            
            
            D_NS = TMath::Sqrt(3/(16*TMath::Pi()))*(l_max*TMath::DegToRad() - l_min*TMath::DegToRad())*(TMath::Power(TMath::Sin(b_max*TMath::DegToRad()),2)-TMath::Power(TMath::Sin(b_min*TMath::DegToRad()),2));
            
            
            D_EW = TMath::Sqrt(3/(8*TMath::Pi()))*(TMath::Sin(b_min*TMath::DegToRad())*TMath::Cos(b_min*TMath::DegToRad()) - TMath::Sin(b_max*TMath::DegToRad())*TMath::Cos(b_max*TMath::DegToRad()) + b_max*TMath::DegToRad()-b_min*TMath::DegToRad())*( TMath::Cos(l_min*TMath::DegToRad()) - TMath::Cos(l_max*TMath::DegToRad()));
            
            
            D_FB = -TMath::Sqrt(3/(8*TMath::Pi()))*(TMath::Sin(b_min*TMath::DegToRad())*TMath::Cos(b_min*TMath::DegToRad()) - TMath::Sin(b_max*TMath::DegToRad())*TMath::Cos(b_max*TMath::DegToRad()) + b_max*TMath::DegToRad()-b_min*TMath::DegToRad())*( TMath::Sin(l_max*TMath::DegToRad()) - TMath::Sin(l_min*TMath::DegToRad()));
            
            
            f1 += TMath::Power( tmp_entries -  tot_entries*(par[0]*I + par[1]*D_NS + par[2]*D_EW + par[3]*D_FB),2)/tot_entries*(par[0]*I + par[1]*D_NS + par[2]*D_EW + par[3]*D_FB);
            
        }
    }
    
    f=f1;
    
}




void Chi2Fit(TH2D &MapsIsoRaio_NS,TH2D &MapsIsoRaio_EW,TH2D &MapsIsoRaio_FB) {
    
    for(Int_t bX=1; bX<=360; bX++) {
        for(Int_t bY=1; bY<=180; bY++) {
            binEntries_NS[bX-1][bY-1] = (Double_t)MapsIsoRaio_NS.GetBinContent(bX,bY);
            binEntries_EW[bX-1][bY-1] = (Double_t)MapsIsoRaio_EW.GetBinContent(bX,bY);
            binEntries_FB[bX-1][bY-1] = (Double_t)MapsIsoRaio_FB.GetBinContent(bX,bY);
        }
    }
    
    //////////////////////////////////// TMinuit Fit
    
    //////////////////// General variables ////////////////////
    
    int err = 0;
    double set_str_list[1];
    set_str_list[0] = 2;
    
    std::string parName[4];
    double par[4],stepSize[4],minVal[4],maxVal[4];
    
    ///////////////////////////////////////////////////////////
    
    ///////////////////// NS Minimization /////////////////////
    
    std::cout<<"\n\n \\ -------------------------  Fitting NS Dipole  ------------------------- \\ \n\n";
    
    TMinuit minimization_NS(4);
    minimization_NS.SetFCN(fcn_Chi2_NS);
    
    par[0] = 1;
    stepSize[0] = 0.001;
    minVal[0] = 0;
    maxVal[0] = 0;
    parName[0] = "I";
    
    par[1] = 0;
    stepSize[1] = 0.001;
    minVal[1] = 0;
    maxVal[1] = 0;
    parName[1] = "D_NS";
    
    par[2] = 0;
    stepSize[2] = 0.001;
    minVal[2] = 0;
    maxVal[2] = 0;
    parName[2] = "D_EW";
    
    par[3] = 0;
    stepSize[3] = 0.001;
    minVal[3] = 0;
    maxVal[3] = 0;
    parName[3] = "D_FB";
    
    minimization_NS.DefineParameter(0,parName[0].c_str(),par[0],stepSize[0],minVal[0],maxVal[0]);
    minimization_NS.DefineParameter(1,parName[1].c_str(),par[1],stepSize[1],minVal[1],maxVal[1]);
    minimization_NS.DefineParameter(2,parName[2].c_str(),par[2],stepSize[2],minVal[2],maxVal[2]);
    minimization_NS.DefineParameter(3,parName[3].c_str(),par[3],stepSize[3],minVal[3],maxVal[3]);
    
    minimization_NS.SetErrorDef(0.5);
    minimization_NS.mnexcm("SET STR",set_str_list,1,err);
    
    minimization_NS.Migrad();
    
    
    ///////////////////// EW Minimization /////////////////////
    
    std::cout<<"\n\n \\ -------------------------  Fitting EW Dipole  ------------------------- \\ \n\n";
    
    TMinuit minimization_EW(4);
    minimization_EW.SetFCN(fcn_Chi2_EW);
    
    par[0] = 1;
    stepSize[0] = 0.001;
    minVal[0] = 0;
    maxVal[0] = 0;
    parName[0] = "I";
    
    par[1] = 0;
    stepSize[1] = 0.001;
    minVal[1] = 0;
    maxVal[1] = 0;
    parName[1] = "D_NS";
    
    par[2] = 0;
    stepSize[2] = 0.001;
    minVal[2] = 0;
    maxVal[2] = 0;
    parName[2] = "D_EW";
    
    par[3] = 0;
    stepSize[3] = 0.001;
    minVal[3] = 0;
    maxVal[3] = 0;
    parName[3] = "D_FB";
    
    minimization_EW.DefineParameter(0,parName[0].c_str(),par[0],stepSize[0],minVal[0],maxVal[0]);
    minimization_EW.DefineParameter(1,parName[1].c_str(),par[1],stepSize[1],minVal[1],maxVal[1]);
    minimization_EW.DefineParameter(2,parName[2].c_str(),par[2],stepSize[2],minVal[2],maxVal[2]);
    minimization_EW.DefineParameter(3,parName[3].c_str(),par[3],stepSize[3],minVal[3],maxVal[3]);
    
    minimization_EW.SetErrorDef(0.5);
    minimization_EW.mnexcm("SET STR",set_str_list,1,err);
    
    minimization_EW.Migrad();
    
    
    ///////////////////// FB Minimization /////////////////////
    
    std::cout<<"\n\n \\ -------------------------  Fitting FB Dipole  ------------------------- \\ \n\n";
    
    TMinuit minimization_FB(4);
    minimization_FB.SetFCN(fcn_Chi2_FB);
    
    par[0] = 1;
    stepSize[0] = 0.001;
    minVal[0] = 0;
    maxVal[0] = 0;
    parName[0] = "I";
    
    par[1] = 0;
    stepSize[1] = 0.001;
    minVal[1] = 0;
    maxVal[1] = 0;
    parName[1] = "D_NS";
    
    par[2] = 0;
    stepSize[2] = 0.001;
    minVal[2] = 0;
    maxVal[2] = 0;
    parName[2] = "D_EW";
    
    par[3] = 0;
    stepSize[3] = 0.001;
    minVal[3] = 0;
    maxVal[3] = 0;
    parName[3] = "D_FB";
    
    minimization_FB.DefineParameter(0,parName[0].c_str(),par[0],stepSize[0],minVal[0],maxVal[0]);
    minimization_FB.DefineParameter(1,parName[1].c_str(),par[1],stepSize[1],minVal[1],maxVal[1]);
    minimization_FB.DefineParameter(2,parName[2].c_str(),par[2],stepSize[2],minVal[2],maxVal[2]);
    minimization_FB.DefineParameter(3,parName[3].c_str(),par[3],stepSize[3],minVal[3],maxVal[3]);
    
    minimization_FB.SetErrorDef(0.5);
    minimization_FB.mnexcm("SET STR",set_str_list,1,err);
    
    minimization_FB.Migrad();
    
    
}

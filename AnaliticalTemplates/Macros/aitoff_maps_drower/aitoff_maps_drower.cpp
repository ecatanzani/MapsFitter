
///////////////////////////////////// Macro description:







////////////////////////////////////////////////////////////

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string>

#include "TFile.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TH2.h"
#include "TH2D.h"
#include "TColor.h"
#include "TMath.h"
#include "TString.h"
#include "TLatex.h"
#include "TText.h"
#include "TVirtualPad.h"
#include "TPaletteAxis.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TRandom3.h"

TVirtualPad* draw_single_map(TH2* h_to_draw,Bool_t logz);
void draw_grid(TVirtualPad* pad,Double_t canv_x_size,Double_t canv_y_size,Bool_t draw_label=false,Int_t n_coo=6,Int_t lat_label=12);

void draw_maps(const std::string resultPath,bool HS=false,bool templates=false) {
    Bool_t logz = false;
    
    ////////////////// Template Maps
    
    TH2* Template_Iso = nullptr;
    TH2* Template_AniNS = nullptr;
    TH2* Template_AniEW = nullptr;
    TH2* Template_AniFB = nullptr;
    
    ////////////////// Data Maps
    
    TH2* Data_Iso_LS = nullptr;
    TH2* Data_AniNS_LS = nullptr;
    TH2* Data_AniEW_LS = nullptr;
    TH2* Data_AniFB_LS = nullptr;
    TH2* MixedData_NS_EW_LS = nullptr;
    TH2* MixedData_NS_FB_LS = nullptr;
    TH2* MixedData_EW_FB_LS = nullptr;
    TH2* FullMixedData_LS = nullptr;
    
    TFile inputFile(resultPath.c_str(),"READ");
    if(inputFile.IsZombie()) {
        std::cerr << "\n\nError opening results TFile " << resultPath.c_str() << "--> Macro finished\n\n";
        exit(100);
    }
    
    if(!templates)
    {
    
        /////////////////////////////////////////////////////////////// Data Maps
    
        if(HS)
        {
        
            Data_Iso_LS = (TH2D*)inputFile.Get("DAMPE_Data_Iso_HS");
            Data_Iso_LS->SetDirectory(0);
    
            Data_AniNS_LS = (TH2D*)inputFile.Get("DAMPE_Data_AniNS_HS");
            Data_AniNS_LS->SetDirectory(0);
    
            Data_AniEW_LS = (TH2D*)inputFile.Get("DAMPE_Data_AniEW_HS");
            Data_AniEW_LS->SetDirectory(0);
    
            Data_AniFB_LS = (TH2D*)inputFile.Get("DAMPE_Data_AniFB_HS");
            Data_AniFB_LS->SetDirectory(0);
    
            MixedData_NS_EW_LS = (TH2D*)inputFile.Get("DAMPE_MixedData_NS_EW_HS");
            MixedData_NS_EW_LS->SetDirectory(0);
    
            MixedData_NS_FB_LS = (TH2D*)inputFile.Get("DAMPE_MixedData_NS_FB_HS");
            MixedData_NS_FB_LS->SetDirectory(0);
    
            MixedData_EW_FB_LS = (TH2D*)inputFile.Get("DAMPE_MixedData_EW_FB_HS");
            MixedData_EW_FB_LS->SetDirectory(0);
    
            FullMixedData_LS = (TH2D*)inputFile.Get("DAMPE_FullMixedData_HS");
            FullMixedData_LS->SetDirectory(0);
    
        }
        else
        {
            
            Data_Iso_LS = (TH2D*)inputFile.Get("DAMPE_Data_Iso_LS");
            Data_Iso_LS->SetDirectory(0);
        
            Data_AniNS_LS = (TH2D*)inputFile.Get("DAMPE_Data_AniNS_LS");
            Data_AniNS_LS->SetDirectory(0);
        
            Data_AniEW_LS = (TH2D*)inputFile.Get("DAMPE_Data_AniEW_LS");
            Data_AniEW_LS->SetDirectory(0);
        
            Data_AniFB_LS = (TH2D*)inputFile.Get("DAMPE_Data_AniFB_LS");
            Data_AniFB_LS->SetDirectory(0);
        
            MixedData_NS_EW_LS = (TH2D*)inputFile.Get("DAMPE_MixedData_NS_EW_LS");
            MixedData_NS_EW_LS->SetDirectory(0);
        
            MixedData_NS_FB_LS = (TH2D*)inputFile.Get("DAMPE_MixedData_NS_FB_LS");
            MixedData_NS_FB_LS->SetDirectory(0);
        
            MixedData_EW_FB_LS = (TH2D*)inputFile.Get("DAMPE_MixedData_EW_FB_LS");
            MixedData_EW_FB_LS->SetDirectory(0);
        
            FullMixedData_LS = (TH2D*)inputFile.Get("DAMPE_FullMixedData_LS");
            FullMixedData_LS->SetDirectory(0);
        
        }
    
        inputFile.Close();
    
        draw_single_map(Data_Iso_LS,logz);
        draw_single_map(Data_AniNS_LS,logz);
        draw_single_map(Data_AniEW_LS,logz);
        draw_single_map(Data_AniFB_LS,logz);
        draw_single_map(MixedData_NS_EW_LS,logz);
        draw_single_map(MixedData_NS_FB_LS,logz);
        draw_single_map(MixedData_EW_FB_LS,logz);
        draw_single_map(FullMixedData_LS,logz);
    
    }
    else
    {
    
        ////////////////// Templates
    
    
        /////////////////////////////////////////////////////////////// Template Maps
    
        //Template_Iso = (TH2D*)inputFile.Get("Template_Iso_HS");
        Template_Iso = (TH2D*)inputFile.Get("DAMPE_Template_Iso_HS");
        Template_Iso->SetDirectory(0);
    
        //Template_AniNS = (TH2D*)inputFile.Get("Template_AniNS_HS");
        Template_AniNS = (TH2D*)inputFile.Get("DAMPE_Template_AniNS_HS");
        Template_AniNS->SetDirectory(0);
    
        //Template_AniEW = (TH2D*)inputFile.Get("Template_AniEW_HS");
        Template_AniEW = (TH2D*)inputFile.Get("DAMPE_Template_AniEW_HS");
        Template_AniEW->SetDirectory(0);
    
        //Template_AniFB = (TH2D*)inputFile.Get("Template_AniFB_HS");
        Template_AniFB = (TH2D*)inputFile.Get("DAMPE_Template_AniFB_HS");
        Template_AniFB->SetDirectory(0);
    
        inputFile.Close();
    
        draw_single_map(Template_Iso,logz);
        draw_single_map(Template_AniNS,logz);
        draw_single_map(Template_AniEW,logz);
        draw_single_map(Template_AniFB,logz);
    
    }
     
}


    //////////////////////////////////////////////////////////////////////////// Drawing functions.....




TVirtualPad* draw_single_map(TH2* h_to_draw,Bool_t logz) {
    
    Double_t canv_x_size=1024;
    Double_t canv_y_size=768;
    
    // Add a little style
    
    gStyle->SetOptStat(0);
    gStyle->SetNumberContours(8);
    gStyle->SetPalette(51);
    //gStyle->SetPalette(60);
    
    TH2* h_to_draw_clone = (TH2D*)(h_to_draw->Clone(Form("%s_clone", h_to_draw->GetName())));
    
    h_to_draw_clone->GetXaxis()->SetLabelSize(0.0);
    h_to_draw_clone->GetXaxis()->SetTickSize(0.0);
    h_to_draw_clone->GetXaxis()->SetTitleOffset(0.5);
    h_to_draw_clone->GetYaxis()->SetLabelSize(0.0);
    h_to_draw_clone->GetYaxis()->SetTickSize(0.0);
    h_to_draw_clone->GetYaxis()->SetTitleOffset(0.5);
    
    TCanvas* c_aitoff = new TCanvas(Form("Aitoff_%s",h_to_draw->GetName()),h_to_draw->GetTitle(),1.1*canv_x_size,canv_y_size);
    c_aitoff->cd();
    
    TPad* pad = new TPad("pad","", 0, 0, 1.0/1.1-0.075, 1.0, 0, 0);
    pad->SetNumber(1);
    pad->Draw();
    c_aitoff->cd(1);
    gPad->SetRightMargin(0.01);
    
    h_to_draw_clone->Draw("aitoff");                  //Draw the histo with the aitoff projection !! You can also use the mercatore projection, if you want !
    
    c_aitoff->cd();
    TPaletteAxis* pal = new TPaletteAxis(1.0/1.1-0.075,0.1,1.0/1.1-0.075+0.05,0.90,h_to_draw_clone);
    pal->Draw();
    
    draw_grid(c_aitoff->cd(1),canv_x_size,canv_y_size);
    
    TVirtualPad* pad_2_return = c_aitoff->cd(1);
    pad_2_return->SetLogz(logz);
    
    c_aitoff->cd();
    c_aitoff->SetLogz(logz);
    
    return pad_2_return;
    
}

void draw_grid(TVirtualPad* pad,Double_t canv_x_size,Double_t canv_y_size,Bool_t draw_label,Int_t n_coo,Int_t lat_label) {
    
    Double_t convx=((1.0/(canv_x_size*1.02))*(10.0/3.0));
    Double_t convy=((1.0/(canv_y_size*0.85*0.8))*(10.0/3.0));
    Bool_t d_grid = true;
    
    Float_t la, lo, _x, _y, z;
    
    const Int_t Nl = 19;    // Number of drawn latitudes
    const Int_t NL = 19;    // Number of drawn longitudes
    Int_t       M  = 90;
    
    Double_t _M;
    
    TGraph  *latitudes[Nl];
    TGraph  *longitudes[NL];
    
    for (int j=0;j<Nl;++j) {
        latitudes[j]=new TGraph();
        latitudes[j]->SetMarkerSize(0.3);
        latitudes[j]->SetMarkerColor(0);      // -> Set grid latitude marker color
        la = -90+180/(Nl-1)*j;
        _M = M;
        
        // The number of point do draw shold be set according to the latitude. Not all the latitudes need the same number of points
        
        if (fabs(la)>50.0) {
            _M = M*3.0/5.0;
        }
        if (fabs(la)>60.0) {
            _M = M*2.0/5.0;
        }
        if (fabs(la)>70.0) {
            _M = M*1.0/5.0;
        }
        for (int i=0;i<_M+1;++i) {
            lo = -180+360/_M*i;
            z  = sqrt(1+cos(la*TMath::DegToRad())*cos(lo*TMath::DegToRad()/2));
            _x  = 180*cos(la*TMath::DegToRad())*sin(lo*TMath::DegToRad()/2)/z;
            _y  = 90*sin(la*TMath::DegToRad())/z;
            latitudes[j]->SetPoint(i,_x*convx,_y*convy);
        }
    }
    
    
    for (int j=0;j<NL;++j) {
        longitudes[j]=new TGraph();
        longitudes[j]->SetMarkerSize(0.3);
        longitudes[j]->SetMarkerColor(0);    //Set grid longitude marker color
        lo = -180+360/(NL-1)*j;
        for (int i=0;i<M+1;++i) {
            la = -90+180/M*i;
            if (fabs(la)<85.0) {
                z  = sqrt(1+cos(la*TMath::DegToRad())*cos(lo*TMath::DegToRad()/2));
                _x  = 180*cos(la*TMath::DegToRad())*sin(lo*TMath::DegToRad()/2)/z;
                _y  = 90*sin(la*TMath::DegToRad())/z;
                longitudes[j]->SetPoint(i,_x*convx,_y*convy);
            }
        }
    }
    
    if (d_grid) {
        
        // Draw the grid. That is done drawing each single latitude and longitude point.
        
        pad->cd();
        for (int j=0;j<Nl;++j) {
            if (j!=0 && j!=(Nl-1)) latitudes[j]->Draw("P");
        }
        for (int j=0;j<NL;++j) {
            longitudes[j]->Draw("P");
        }
    }
    
    // draw label on map
    
    if (draw_label==true){
        char coo_name[255];
        double* vec_xx = (double*)latitudes[lat_label]->GetX();
        int Npoints = (int)latitudes[lat_label]->GetN();
        double delta_coo = (double)Npoints/(double)n_coo;
        double delta_coo_for_label = 360./(double)n_coo;
        double coo_for_label = -180.;
        int num_coo = 0;
        for (int nn=0; nn<=n_coo; nn++){
            double y = latitudes[lat_label]->Eval(vec_xx[num_coo]);
            double x = vec_xx[num_coo];
            sprintf(coo_name,"%.0f.",coo_for_label);
            TText* tex = new TLatex(x,y,coo_name);
            tex->SetTextFont(40);
            tex->SetTextSize(0.03);
            tex->Draw();
            num_coo = num_coo+delta_coo;
            coo_for_label = coo_for_label+delta_coo_for_label;
        }
    }
}

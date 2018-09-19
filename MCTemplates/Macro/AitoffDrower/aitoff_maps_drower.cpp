
///////////////////////////////////// Macro description:







////////////////////////////////////////////////////////////

#include <iostream>
#include <stdlib.h>
#include <stdio.h>

#include "TFile.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TH2.h"
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
void draw_pull(TH2* h_ratio,TH2* h_num,TH2* h_den,TH1D* h_pull[],Int_t idx_pool);
void filling_entries(TH1I* h_entries,TH2* h_map);
void pull_correlation_MC(TString out_path);
void obtain_correlation_map(TH2* h_map1,TH2* h_map2,TH2I* h_correl);

using namespace std;

void draw_maps(const TString results_alias) {
    Bool_t logz = false;
  
    ////////////////////////////// Set input and output path
  
    TString input_results_path = "../../results/";
    input_results_path+=results_alias;
    input_results_path+="_maps_result.root";

    ///////////////////////////////////////////////////////

    gStyle->SetOptStat(0);
    gStyle->SetNumberContours(8);

    TFile *input_results_file = new TFile(input_results_path.Data());
    if(input_results_file->IsZombie()) {
        cout<<"\n\nError opening file. Macro finished\n\n";
        exit(-1);
    }
    
    /////////////////////////////////////////////////////////////// Isotropic and anisotropic maps
    
#if 0
    
    ///////// Draw isotropic map
    
    TH2* Iso_SkyMap = (TH2D*)(input_results_file->Get("IsoMap"));
    TH2* NIso_SkyMap = (TH2D*)(input_results_file->Get("NIsoMap"));
    
    draw_single_map(Iso_SkyMap,logz);
    draw_single_map(NIso_SkyMap,logz);
    
    ///////// Draw anisotropic map
    
    TH2* AniMap_NS = (TH2D*)(input_results_file->Get("AniMap_NS"));
    TH2* AniMap_EW = (TH2D*)(input_results_file->Get("AniMap_EW"));
    TH2* AniMap_FB = (TH2D*)(input_results_file->Get("AniMap_FB"));
    
    TH2* NAniMap_NS = (TH2D*)(input_results_file->Get("NAniMap_NS"));
    TH2* NAniMap_EW = (TH2D*)(input_results_file->Get("NAniMap_EW"));
    TH2* NAniMap_FB = (TH2D*)(input_results_file->Get("NAniMap_FB"));
    
    TH2* ph_AniMap_NS = (TH2D*)(input_results_file->Get("ph_AniMap_NS"));
    TH2* ph_AniMap_EW = (TH2D*)(input_results_file->Get("ph_AniMap_EW"));
    TH2* ph_AniMap_FB = (TH2D*)(input_results_file->Get("ph_AniMap_FB"));
    
    TH2* ph_NAniMap_NS = (TH2D*)(input_results_file->Get("phNAniMap_NS"));
    TH2* ph_NAniMap_EW = (TH2D*)(input_results_file->Get("phNAniMap_EW"));
    TH2* ph_NAniMap_FB = (TH2D*)(input_results_file->Get("phNAniMap_FB"));
    
    TH2* RatioAniMap_NS = (TH2D*)(input_results_file->Get("RatioAniMap_NS"));
    TH2* RatioAniMap_EW = (TH2D*)(input_results_file->Get("RatioAniMap_EW"));
    TH2* RatioAniMap_FB = (TH2D*)(input_results_file->Get("RatioAniMap_FB"));
    
    draw_single_map(AniMap_NS,logz);
    draw_single_map(AniMap_EW,logz);
    draw_single_map(AniMap_FB,logz);
    
    draw_single_map(NAniMap_NS,logz);
    draw_single_map(NAniMap_EW,logz);
    draw_single_map(NAniMap_FB,logz);
    
    draw_single_map(ph_AniMap_NS,logz);
    draw_single_map(ph_AniMap_EW,logz);
    draw_single_map(ph_AniMap_FB,logz);
    
    draw_single_map(ph_NAniMap_NS,logz);
    draw_single_map(ph_NAniMap_EW,logz);
    draw_single_map(ph_NAniMap_FB,logz);
    
    draw_single_map(RatioAniMap_NS,logz);
    draw_single_map(RatioAniMap_EW,logz);
    draw_single_map(RatioAniMap_FB,logz);
    
    ///////// Evidence Simple Dipoles
    
    TH2* Isolate_NS = (TH2D*)(input_results_file->Get("Isolate_NS"));
    TH2* Isolate_EW = (TH2D*)(input_results_file->Get("Isolate_EW"));
    TH2* Isolate_FB = (TH2D*)(input_results_file->Get("Isolate_FB"));
    
    draw_single_map(Isolate_NS,logz);
    draw_single_map(Isolate_EW,logz);
    draw_single_map(Isolate_FB,logz);

#else
    /*
    TH2* Template_HS_Iso = (TH2D*)(input_results_file->Get("Template_Iso_HS"));
    TH2* Template_HS_AniNS = (TH2D*)(input_results_file->Get("Template_AniNS_HS"));
    TH2* Template_HS_AniEW = (TH2D*)(input_results_file->Get("Template_AniEW_HS"));
    TH2* Template_HS_AniFB = (TH2D*)(input_results_file->Get("Template_AniFB_HS"));
    
    TH2* Data_Iso_HS = (TH2D*)(input_results_file->Get("Data_Iso_HS"));
    TH2* Data_AniNS_HS = (TH2D*)(input_results_file->Get("Data_AniNS_HS"));
    TH2* Data_AniEW_HS = (TH2D*)(input_results_file->Get("Data_AniEW_HS"));
    TH2* Data_AniFB_HS = (TH2D*)(input_results_file->Get("Data_AniFB_HS"));
    */
    TH2* MixedData_NS_EW_HS = (TH2D*)(input_results_file->Get("MixedData_NS_EW_HS"));
    TH2* MixedData_NS_FB_HS = (TH2D*)(input_results_file->Get("MixedData_NS_FB_HS"));
    TH2* MixedData_EW_FB_HS = (TH2D*)(input_results_file->Get("MixedData_EW_FB_HS"));
    TH2* FullMixedData_HS = (TH2D*)(input_results_file->Get("FullMixedData_HS"));
    /*
    draw_single_map(Template_HS_Iso,logz);
    draw_single_map(Template_HS_AniNS,logz);
    draw_single_map(Template_HS_AniEW,logz);
    draw_single_map(Template_HS_AniFB,logz);
    
    draw_single_map(Data_Iso_HS,logz);
    draw_single_map(Data_AniNS_HS,logz);
    draw_single_map(Data_AniEW_HS,logz);
    draw_single_map(Data_AniFB_HS,logz);
    */
    draw_single_map(MixedData_NS_EW_HS,logz);
    draw_single_map(MixedData_NS_FB_HS,logz);
    draw_single_map(MixedData_EW_FB_HS,logz);
    draw_single_map(FullMixedData_HS,logz);
    
    /*
    
    TH2* AllSky_IsoMap = (TH2D*)(input_results_file->Get("AllSky_IsoMap"));
    
    TH2* AllSky_AniMap_NS = (TH2D*)(input_results_file->Get("AllSky_AniMap_NS"));
    TH2* AllSky_AniMap_EW = (TH2D*)(input_results_file->Get("AllSky_AniMap_EW"));
    TH2* AllSky_AniMap_FB = (TH2D*)(input_results_file->Get("AllSky_AniMap_FB"));
    
    TH2* AllSky_Isolate_NS = (TH2D*)(input_results_file->Get("AllSky_Isolate_NS"));
    TH2* AllSky_Isolate_EW = (TH2D*)(input_results_file->Get("AllSky_Isolate_EW"));
    TH2* AllSky_Isolate_FB = (TH2D*)(input_results_file->Get("AllSky_Isolate_FB"));
    
    draw_single_map(AllSky_IsoMap,logz);
    
    draw_single_map(AllSky_AniMap_NS,logz);
    draw_single_map(AllSky_AniMap_EW,logz);
    draw_single_map(AllSky_AniMap_FB,logz);
    
    draw_single_map(AllSky_Isolate_NS,logz);
    draw_single_map(AllSky_Isolate_EW,logz);
    draw_single_map(AllSky_Isolate_FB,logz);
    
     */
     
#endif
    
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


    void draw_pull(TH2* h_ratio,TH2* h_num,TH2* h_den,TH1D* h_pull[],Int_t idx_pool) {

      Double_t pull_min = 0;
      Double_t pull_max = 0;

      Double_t y_rec,y_exp,sigma;
      Double_t pull;

      
      /*
       
        This first two loops are just needed because I wanto to obtain max and min value for the pull histo.
        These two values will be needed to proprly set the pool's histos max and min value for the X axis.
        The second loop just fills poll's histo ! So obviously I need a couple of loops.

       */

      for(Int_t xx=1; xx<=h_ratio->GetNbinsX(); xx++) {
        for(Int_t yy=1; yy<=h_ratio->GetNbinsY(); yy++) {
          y_rec = h_num->GetBinContent(xx,yy);
          y_exp = h_den->GetBinContent(xx,yy);
          if(y_rec<10 || y_exp<10) {
        //cout<<"\nLow single bin statistic !!";
        continue;
          }
          //sigma = sqrt(y_exp+y_rec);
          sigma = sqrt(y_exp)+sqrt(y_rec);
          //sigma = sqrt(y_exp);
          if (sigma>0) {
            pull = (y_rec-y_exp)/sigma;
            if (pull<pull_min)
          pull_min=pull;
            if (pull>pull_max)
          pull_max=pull;
          } //end if
        } //end for yy
      } //end for xx

      //h_pull[idx_pool] = new TH1D(Form("%s_pull", h_ratio->GetName()), Form("Pull of %s", h_ratio->GetTitle()), 4*3*5*7, 1.2*pull_min, 1.2*pull_max);
      h_pull[idx_pool] = new TH1D(Form("%s_pull", h_ratio->GetName()), Form("Pull of %s", h_ratio->GetTitle()),150, 1.2*pull_min, 1.2*pull_max);
      h_pull[idx_pool]->GetXaxis()->SetTitle("(y_{rec}-y_{ref})/(#sqrt{#sigma^{2}_{rec}+#sigma^{2}_{ref}})");
      h_pull[idx_pool]->GetYaxis()->SetTitle("Entries");

      for(Int_t xx=1; xx<=h_ratio->GetNbinsX(); xx++) {
        for(Int_t yy=1; yy<=h_ratio->GetNbinsY(); yy++) {
          y_rec = h_num->GetBinContent(xx, yy);
          y_exp = h_den->GetBinContent(xx, yy);
          if(y_rec<10 || y_exp<10) {
        //cout<<"\nLow single bin statistic !!";
        continue;
          }
          sigma = sqrt(y_exp)+sqrt(y_rec);
          //sigma = sqrt(y_exp);
          if (y_rec<0 || y_exp<0)
        printf("%f %f\n", y_rec, y_exp);
          if (sigma>0) {
            pull = (y_rec-y_exp)/sigma;
            h_pull[idx_pool]->Fill(pull);
          }
        }
      }
    }

void filling_entries(TH1I* h_entries,TH2* h_map) {
    for(Int_t idx_bx=1; idx_bx<=h_map->GetNbinsX(); idx_bx++)
        for(Int_t idx_by=1; idx_by<=h_map->GetNbinsY(); idx_by++)
            if(h_map->GetBinContent(idx_bx,idx_by)>0)
                h_entries->Fill(h_map->GetBinContent(idx_bx,idx_by));
    
}

void pull_correlation_MC(TString out_path) {
    UInt_t seed=22;
    Double_t n_events=1e+9;
    Int_t ent_1,ent_2;
    Double_t sigma,pull;
    TRandom3 *rand = new TRandom3(seed);
    Double_t tmp_bin1,tmp_bin2;
    
    TFile *results = new TFile(out_path.Data(),"RECREATE");
    if(results->IsZombie()) {
        cout<<"\n\nError writing pool correlation study output file\n\n";
        exit(-1);
    }
    
    // Two independent uniform distributions
    
    gStyle->SetOptStat(1);
    TH1D *h_uniform1 = new TH1D("h_uniform1","Unifrom Distribution",1e+5,0,1);
    TH1D *h_uniform2 = new TH1D("h_uniform2","Uniform Distribution",1e+5,0,1);
    TH1D *h_pull_uniform = new TH1D("h_pool_uniform","Pool from Uniform Distributions",1000,-5,5);
    TH1I *h_uniform1_entries = new TH1I("h_uniform1_entries","Entries Distribution Uniform 1",1000,9e+3,11e+3);
    TH1I *h_uniform2_entries = new TH1I("h_uniform2_entries","Entries Distribution Uniform 2",1000,9e+3,11e+3);
    TH2D *uniform_correl = new TH2D("uniform_correl","Correlation Uniform Distribution",1000,9e+3,11e+3,1000,9e+3,11e+3);
    
    //TH2I *uniform_correl = new TH2I("uniform_correl","Uniform Maps Correlation",1000)
    
    for(Int_t idx_l=0; idx_l<n_events; idx_l++) {
        h_uniform1->Fill(rand->Uniform(0,1));
        h_uniform2->Fill(rand->Uniform(0,1));
    }
    
    for(Int_t idx_bx=1; idx_bx<=h_uniform1->GetNbinsX(); idx_bx++) {
        h_uniform1_entries->Fill(h_uniform1->GetBinContent(idx_bx));
        h_uniform2_entries->Fill(h_uniform2->GetBinContent(idx_bx));
    }
    
    for(Int_t idx_bx=1; idx_bx<=h_uniform1->GetNbinsX(); idx_bx++) {
        tmp_bin1=h_uniform1->GetBinContent(idx_bx);
        tmp_bin2=h_uniform2->GetBinContent(idx_bx);
        uniform_correl->Fill(tmp_bin1,tmp_bin2);
    }
        
    for(Int_t idx_b=1; idx_b<=h_uniform1->GetNbinsX(); idx_b++) {
        ent_1=h_uniform1->GetBinContent(idx_b);
        ent_2=h_uniform2->GetBinContent(idx_b);
        if(ent_1<10 || ent_2<10)
            continue;
        //sigma=sqrt(ent_1+ent_2);
        sigma=sqrt(ent_1+ent_2+2*(uniform_correl->GetCorrelationFactor())*sqrt(ent_1)*sqrt(ent_2));
        if(sigma>0) {
            pull=(ent_1-ent_2)/sigma;
            h_pull_uniform->Fill(pull);
        }
    }
    
    // Two independent gaussian distributions
    
    gStyle->SetOptStat(1);
    TH1D *h_gaus1 = new TH1D("h_gaus1","Gaussian Distribution",1e+5,-5,5);
    TH1D *h_gaus2 = new TH1D("h_gaus2","Gaussian Distribution",1e+5,-5,5);
    TH1D *h_pull_gaus = new TH1D("h_pool_gaus","Pool from Gaussian Distributions",1000,-5,5);
    TH1I *h_gaus1_entries = new TH1I("h_gaus1_entries","Entries Distribution Uniform 1",1000,0,45e+3);
    TH1I *h_gaus2_entries = new TH1I("h_gaus2_entries","Entries Distribution Uniform 2",1000,0,45e+3);
    TH2D *gaus_correl = new TH2D("gaus_correl","Correlation Gaussian Distribution",1000,0,45e+3,1000,0,45e+3);
    
    for(Int_t idx_l=0; idx_l<n_events; idx_l++) {
        h_gaus1->Fill(rand->Gaus(0,1));
        h_gaus2->Fill(rand->Gaus(0,1));
    }
    
    for(Int_t idx_bx=1; idx_bx<=h_gaus1->GetNbinsX(); idx_bx++) {
        h_gaus1_entries->Fill(h_gaus1->GetBinContent(idx_bx));
        h_gaus2_entries->Fill(h_gaus2->GetBinContent(idx_bx));
    }
    
    for(Int_t idx_bx=1; idx_bx<=h_gaus1->GetNbinsX(); idx_bx++) {
        tmp_bin1=h_gaus1->GetBinContent(idx_bx);
        tmp_bin2=h_gaus2->GetBinContent(idx_bx);
        gaus_correl->Fill(tmp_bin1,tmp_bin2);
    }
    
    for(Int_t idx_b=1; idx_b<=h_gaus1->GetNbinsX(); idx_b++) {
        ent_1=h_gaus1->GetBinContent(idx_b);
        ent_2=h_gaus2->GetBinContent(idx_b);
        if(ent_1<10 || ent_2<10)
            continue;
        sigma=sqrt(ent_1+ent_2+2*(gaus_correl->GetCorrelationFactor())*sqrt(ent_1)*sqrt(ent_2));
        //sigma=sqrt(ent_1)+sqrt(ent_2);
        //sigma=sqrt(ent_1+ent_2);
        if(sigma>0) {
            pull=(ent_1-ent_2)/sigma;
            h_pull_gaus->Fill(pull);
        }
    }
    
    
    results->Write();
    results->Close();
    
}

void obtain_correlation_map(TH2* h_map1,TH2* h_map2,TH2I* h_correl) {
    Int_t tmp_bin1,tmp_bin2;
    for(Int_t idx_bx=1; idx_bx<=h_map1->GetNbinsX(); idx_bx++)
        for(Int_t idx_by=1; idx_by<=h_map1->GetNbinsY(); idx_by++) {
            tmp_bin1=h_map1->GetBinContent(idx_bx,idx_by);
            tmp_bin2=h_map2->GetBinContent(idx_bx,idx_by);
            h_correl->Fill(tmp_bin1,tmp_bin2);
        }
}

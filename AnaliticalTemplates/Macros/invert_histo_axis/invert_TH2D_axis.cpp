
///////////////////////////////////// Macro description:







////////////////////////////////////////////////////////////

#include <string>

#include "TGaxis.h"
#include "TH2D.h"
#include "TCanvas.h"

void ReverseXAxis(TH2D &histo);
void ReverseYAxis(TH2D &histo);
void draw_X_label(TH1 *h);
void draw_Y_label(TH1 *h);
void uniform_axis(TH1 *h,bool Xaxis=false,bool Yaxis=false);

void invert_TH2F_axis(TH2D* hToInvert,bool Xaxis=true,bool Yaxis=false,std::string draw_options="")
{
    
    // Reverse the histo axis
    
    if(Xaxis)
        ReverseXAxis(*hToInvert);
    if(Yaxis)
        ReverseYAxis(*hToInvert);
    
   // Redraw correctly the labels
    
    hToInvert->Draw(draw_options.c_str());
    
    if(Xaxis)
    {
        draw_X_label(hToInvert);
        uniform_axis(hToInvert,false,true);
    }
    if(Yaxis)
    {
        draw_Y_label(hToInvert);
        uniform_axis(hToInvert,true,false);
    }
    
}

void ReverseXAxis (TH2D &histo)
{
    TH2D *tmpHisto = (TH2D*)histo.Clone("tmpHisto");
    tmpHisto->Reset();
    UInt_t nbinX = histo.GetNbinsX();
    
    for(UInt_t bY=1; bY <= (UInt_t)histo.GetNbinsY(); ++bY)
        for(UInt_t bX=1; bX <= nbinX; ++bX)
            tmpHisto->SetBinContent(nbinX - (bX -1),bY,histo.GetBinContent(bX,bY));
    
    new (&histo) (TH2D) (*(TH2D*) tmpHisto->Clone(histo.GetName()));
}

void ReverseYAxis (TH2D &histo)
{
    TH2D *tmpHisto = (TH2D*)histo.Clone("tmpHisto");
    tmpHisto->Reset();
    UInt_t nbinY = histo.GetNbinsY();
    
    for(UInt_t bX=1; bX <= (UInt_t)histo.GetNbinsX(); ++bX)
        for(UInt_t bY=1; bY <= nbinY; ++bY)
            tmpHisto->SetBinContent(bX,nbinY - (bY -1),histo.GetBinContent(bX,bY));
    
    new (&histo) (TH2D) (*(TH2D*) tmpHisto->Clone(histo.GetName()));
}

void draw_X_label(TH1 *h)
{
    // Remove the current axis
    
    h->GetXaxis()->SetLabelOffset(999);
    h->GetXaxis()->SetTickLength(0);
    
    gPad->Update();
    
    TGaxis *newaxis = new TGaxis(
                                    gPad->GetUxmax(),
                                    gPad->GetUymin(),
                                    gPad->GetUxmin(),
                                    gPad->GetUymin(),
                                    h->GetXaxis()->GetXmin(),
                                    h->GetXaxis()->GetXmax(),
                                    510,
                                    "-"
                                 );
    
    newaxis->SetLabelOffset(-0.03);
    
    newaxis->Draw();

}

void draw_Y_label(TH1 *h)
{
    // Remove the current axis
    
    h->GetYaxis()->SetLabelOffset(999);
    h->GetYaxis()->SetTickLength(0);
    
    // Redraw the new axis
    
    gPad->Update();
    
    TGaxis *newaxis = new TGaxis(
                                    gPad->GetUxmin(),
                                    gPad->GetUymax(),
                                    gPad->GetUxmin()-0.001,
                                    gPad->GetUymin(),
                                    h->GetYaxis()->GetXmin(),
                                    h->GetYaxis()->GetXmax(),
                                    510,
                                    "+"
                                 );
    
    newaxis->SetLabelOffset(-0.03);
    newaxis->Draw();
}

void uniform_axis(TH1 *h,bool Xaxis,bool Yaxis)
{
    
    if(Xaxis)
    {
        // Remove the current axis
    
        h->GetXaxis()->SetLabelOffset(999);
        h->GetXaxis()->SetTickLength(0);
    
        // Redraw the new axis
    
        gPad->Update();
    
        TGaxis *newaxis = new TGaxis(
                                        gPad->GetUxmin(),
                                        gPad->GetUymin(),
                                        gPad->GetUxmax(),
                                        gPad->GetUymin(),
                                        h->GetXaxis()->GetXmin(),
                                        h->GetXaxis()->GetXmax(),
                                        510,
                                        "+"
                                    );
    
        newaxis->SetLabelOffset(0.03);
        newaxis->Draw();
    }
    
    if(Yaxis)
    {
        // Remove the current axis
        
        h->GetYaxis()->SetLabelOffset(999);
        h->GetYaxis()->SetTickLength(0);
        
        // Redraw the new axis
        
        gPad->Update();
        
        TGaxis *newaxis = new TGaxis(
                                         gPad->GetUxmin(),
                                         gPad->GetUymin(),
                                         gPad->GetUxmin()-0.001,
                                         gPad->GetUymax(),
                                         h->GetYaxis()->GetXmin(),
                                         h->GetYaxis()->GetXmax(),
                                         510,
                                         "-"
                                     );
        
        newaxis->SetLabelOffset(0.03);
        newaxis->Draw();
        
    }
    
}

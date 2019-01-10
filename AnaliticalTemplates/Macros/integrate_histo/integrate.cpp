#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TMath.h"
#include "TGraph.h"
#include "TAttMarker.h"

class reduced_fitResult
{
public:
    
    Double_t delta_LS[8];
    Double_t delta_HS[8];
    
    
    reduced_fitResult()
    {
        for(Int_t idx=0; idx<8; ++idx)
        {
            delta_LS[idx]=-999;
            delta_HS[idx]=-999;
        }
    }
    
    ~reduced_fitResult() { }
    
};

#define n_pts 4

void normalize_histos(std::vector<TH1D> &hdelta);
Double_t compute_integral(TH1D histo);
void get_68_CL(std::vector<TH1D> &hdelta,Double_t integrals_68[]);
void get_95_CL(std::vector<TH1D> &hdelta,Double_t integrals_95[]);
void get_delta_histos(std::string paths[],std::vector<TH1D> &hdelta);
void tree_linking(TTree* myTree,reduced_fitResult &tmp_result,UInt_t &tree_entries);

Double_t c_fitfunc(Double_t *x, Double_t *par);

int main(int argc,char* argv[])
{
    std::vector<TH1D> hdelta;
    hdelta.resize(n_pts);
    
    Double_t times[n_pts] = {1.5,3,6,10};
    Double_t integrals_68[n_pts];
    Double_t integrals_95[n_pts];
    
    std::string inputPath[2] = {argv[1],argv[2]};
    
    std::string outputPath = argv[3];
    
    ///////////////////// Opening delta histos

    get_delta_histos(inputPath,hdelta);
    
    ///////////////////// Computing integrals
    
    normalize_histos(hdelta);
    get_68_CL(hdelta,integrals_68);
    get_95_CL(hdelta,integrals_95);
    
    TGraph sensitivity_68(hdelta.size(),times,integrals_68);
    sensitivity_68.SetName("sensitivity_68_CL");
    sensitivity_68.SetTitle("Sensitivity simulation (68% CL)");
    sensitivity_68.GetXaxis()->SetTitle("acquisition time (years)");
    sensitivity_68.GetYaxis()->SetTitle("sensitivity");
    sensitivity_68.SetMarkerStyle(22);
    
    
    TGraph sensitivity_95(hdelta.size(),times,integrals_95);
    sensitivity_95.SetName("sensitivity_95_CL");
    sensitivity_95.SetTitle("Sensitivity simulation (95% CL)");
    sensitivity_95.GetXaxis()->SetTitle("acquisition time (years)");
    sensitivity_95.GetYaxis()->SetTitle("sensitivity");
    sensitivity_95.SetMarkerStyle(22);
    
    ///////////////////// Fitting TGraphs
    
    Double_t par_0_68 = 0;
    Double_t par_1_68 = 0;
    Double_t par_2_68 = 0;
    Double_t par_0_95 = 0;
    Double_t par_1_95 = 0;
    Double_t par_2_95 = 0;
    
    Double_t max_fit_vale = 10;
    
    TF1 fitFcn("fitFcn",c_fitfunc,0,max_fit_vale,3);
    
    sensitivity_68.Fit("fitFcn");
    
    par_0_68 = fitFcn.GetParameter(0);
    par_1_68 = fitFcn.GetParameter(1);
    par_2_68 = fitFcn.GetParameter(2);
    
    sensitivity_95.Fit("fitFcn");
    
    par_0_95 = fitFcn.GetParameter(0);
    par_1_95 = fitFcn.GetParameter(1);
    par_2_95 = fitFcn.GetParameter(2);
    
    TF1 cl68_fitfunc("cl68_fitfunc",c_fitfunc,0,max_fit_vale,3);
    cl68_fitfunc.SetParameter(0,par_0_68);
    cl68_fitfunc.SetParameter(1,par_1_68);
    cl68_fitfunc.SetParameter(2,par_2_68);
    
    TF1 cl95_fitfunc("cl95_fitfunc",c_fitfunc,0,max_fit_vale,3);
    cl95_fitfunc.SetParameter(0,par_0_95);
    cl95_fitfunc.SetParameter(1,par_1_95);
    cl95_fitfunc.SetParameter(2,par_2_95);
    
    ///////////////////// Writing final result
    
    TFile outFile(outputPath.c_str(),"RECREATE");
    if(outFile.IsZombie())
    {
        std::cout << "\n\n Error writing ROOT output file \n\n";
        exit(100);
    }
    
    sensitivity_68.Write();
    sensitivity_95.Write();
    
    cl68_fitfunc.Write();
    cl95_fitfunc.Write();
    
    outFile.Write();
    outFile.Close();
    
}

void normalize_histos(std::vector<TH1D> &hdelta)
{
    for(Int_t idx=0; idx<hdelta.size(); ++idx)
        hdelta[idx].Scale(1/compute_integral(hdelta[idx]));
    
}

Double_t compute_integral(TH1D histo)
{
    Double_t sum = 0;
    for(Int_t idx=1; idx<=histo.GetNbinsX(); ++idx)
        sum += histo.GetBinContent(idx)*histo.GetBinWidth(idx);
    
    return sum;
}

void get_68_CL(std::vector<TH1D> &hdelta,Double_t integrals_68[])
{
    for(Int_t idx=0; idx<hdelta.size(); ++idx)
    {
        Double_t tmp_sum = 0;
        for(Int_t bX=1; bX <= hdelta[idx].GetNbinsX(); ++bX)
        {
            tmp_sum += hdelta[idx].GetBinContent(bX)*hdelta[idx].GetBinWidth(bX);
            if(tmp_sum >= 0.68)
            {
                integrals_68[idx] = hdelta[idx].GetBinCenter(bX);
                break;
            }
        }
        
    }
    
}

void get_95_CL(std::vector<TH1D> &hdelta,Double_t integrals_95[])
{
    for(Int_t idx=0; idx<hdelta.size(); ++idx)
    {
        Double_t tmp_sum = 0;
        for(Int_t bX=1; bX <= hdelta[idx].GetNbinsX(); ++bX)
        {
            tmp_sum += hdelta[idx].GetBinContent(bX)*hdelta[idx].GetBinWidth(bX);
            if(tmp_sum >= 0.95)
            {
                integrals_95[idx] = hdelta[idx].GetBinCenter(bX);
                break;
            }
        }
        
    }
    
}


void get_delta_histos(std::string paths[],std::vector<TH1D> &hdelta)
{
    UInt_t tree_entries = 0;
    
    TH1D delta_16("delta_16","delta distribution (16000)",100,0,0.1);
    TH1D delta_32("delta_32","delta distribution (32000)",100,0,0.06);
    TH1D delta_64("delta_64","delta distribution (64000)",100,0,0.04);
    TH1D delta_110("delta_110","delta distribution (110000)",100,0,0.035);
    
    for(Int_t idx=0; idx<2; ++idx)
    {
        TFile inFile(paths[idx].c_str(),"READ");
        if(inFile.IsZombie())
        {
            std::cerr << "\n\nError opening TTree ROOT file \n\n";
            exit(100);
        }
    
        TTree* myTree = (TTree*)inFile.Get("fiTree");
        reduced_fitResult tmp_result;
    
        tree_linking(myTree,tmp_result,tree_entries);
        
        for(Int_t tree_idx=0; tree_idx<tree_entries; ++tree_idx)
        {
            myTree->GetEntry(tree_idx);
            if(idx==0)
            {
                delta_16.Fill(tmp_result.delta_LS[0]);
                delta_32.Fill(tmp_result.delta_HS[0]);
            }
            else
            {
                delta_64.Fill(tmp_result.delta_LS[0]);
                delta_110.Fill(tmp_result.delta_HS[0]);
            }
        }
        
        inFile.Close();
        
    }
    
    new (&hdelta[0]) (TH1D) (*(TH1D*)delta_16.Clone("ddelta_16"));
    new (&hdelta[1]) (TH1D) (*(TH1D*)delta_32.Clone("ddelta_32"));
    new (&hdelta[2]) (TH1D) (*(TH1D*)delta_64.Clone("ddelta_64"));
    new (&hdelta[3]) (TH1D) (*(TH1D*)delta_110.Clone("ddelta_110"));
    
}


void tree_linking(TTree* myTree,reduced_fitResult &tmp_result,UInt_t &tree_entries)
{
    myTree->SetBranchAddress("delta_LS",tmp_result.delta_LS);
    myTree->SetBranchAddress("delta_HS",tmp_result.delta_HS);
    
    tree_entries = myTree->GetEntries();
}


Double_t c_fitfunc(Double_t *x, Double_t *par)
{
    
    return par[0] + par[1]/TMath::Sqrt(x[0]-par[2]);
    
}
        
    


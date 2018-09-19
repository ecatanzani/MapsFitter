
#include "MyHead.h"

void components_analysis(
                            std::vector<TH1D*> &TemplatesProjections,
                            std::vector<TH1D*> &DataProjections,
                            double resfullFitResults_HS[][8],
                            std::ofstream &log_file
                         )

{
    std::string stacks_out_path = output_path_creator(6);
    std::vector<THStack> HStack_Data;
    HStack_Data.resize(16);
    
    std::vector<std::vector<TH1D*>> tmpTemplatesProjections;
    tmpTemplatesProjections.resize(8);
    
    std::for_each( tmpTemplatesProjections.begin(), tmpTemplatesProjections.end(),
                  std::bind2nd( std::mem_fun_ref( &std::vector<TH1D*>::resize ), 8 ) );
    
    TString templateName[8];
    
    for(Int_t f_comp=0; f_comp<8; ++f_comp)
    {
        
        templateName[0] = "cTemplate_Iso_HS_projX_";
        templateName[1] = "cTemplate_Iso_HS_projY_";
        templateName[2] = "cTemplate_AniNS_HS_projX_";
        templateName[3] = "cTemplate_AniNS_HS_projY_";
        templateName[4] = "cTemplate_AniEW_HS_projX_";
        templateName[5] = "cTemplate_AniEW_HS_projY_";
        templateName[6] = "cTemplate_AniFB_HS_projX_";
        templateName[7] = "cTemplate_AniFB_HS_projY_";
        
        for(Int_t idxt=0; idxt<8; ++idxt)
            templateName[idxt] += f_comp;
        
        //////////////// Cloning templates' projections
    
        tmpTemplatesProjections[0][f_comp] = (TH1D*) TemplatesProjections[0]->Clone(templateName[0]);
        tmpTemplatesProjections[1][f_comp] = (TH1D*) TemplatesProjections[1]->Clone(templateName[1]);
    
        tmpTemplatesProjections[2][f_comp] = (TH1D*) TemplatesProjections[2]->Clone(templateName[2]);
        tmpTemplatesProjections[3][f_comp] = (TH1D*) TemplatesProjections[3]->Clone(templateName[3]);
    
        tmpTemplatesProjections[4][f_comp] = (TH1D*) TemplatesProjections[4]->Clone(templateName[4]);
        tmpTemplatesProjections[5][f_comp] = (TH1D*) TemplatesProjections[5]->Clone(templateName[5]);
    
        tmpTemplatesProjections[6][f_comp] = (TH1D*) TemplatesProjections[6]->Clone(templateName[6]);
        tmpTemplatesProjections[7][f_comp] = (TH1D*) TemplatesProjections[7]->Clone(templateName[7]);

        //////////////// Defining new tmp THStacks
    
        TString XStackName = "XStack_";
        TString YStackName = "YStack_";
    
        switch(f_comp)
        {
            case 0:
                XStackName += "Iso_HS";
                YStackName += "Iso_HS";
                break;
            case 1:
                XStackName += "AniNS_HS";
                YStackName += "AniNS_HS";
                break;
            case 2:
                XStackName += "AniEW_HS";
                YStackName += "AniEW_HS";
                break;
            case 3:
                XStackName += "AniFB_HS";
                YStackName += "AniFB_HS";
                break;
            case 4:
                XStackName += "Mixed_NS_EW_HS";
                YStackName += "Mixed_NS_EW_HS";
                break;
            case 5:
                XStackName += "Mixed_NS_FB_HS";
                YStackName += "Mixed_NS_FB_HS";
                break;
            case 6:
                XStackName += "Mixed_EW_FB_HS";
                YStackName += "Mixed_EW_FB_HS";
                break;
            case 7:
                XStackName += "FullMixedData_HS";
                YStackName += "FullMixedData_HS";
                break;
        }
    
        HStack_Data[2*f_comp].SetName(XStackName);
        HStack_Data[2*f_comp].SetTitle(XStackName);
        HStack_Data[2*f_comp+1].SetName(YStackName);
        HStack_Data[2*f_comp+1].SetTitle(YStackName);
        
        ///////////////////////////////////////////////////
    
        for(Int_t idx=1; idx<= tmpTemplatesProjections[0][f_comp]->GetNbinsX(); ++idx)
        {
            tmpTemplatesProjections[0][f_comp]->SetBinContent(idx,tmpTemplatesProjections[0][f_comp]->GetBinContent(idx)*resfullFitResults_HS[0][f_comp]);
            tmpTemplatesProjections[1][f_comp]->SetBinContent(idx,tmpTemplatesProjections[1][f_comp]->GetBinContent(idx)*resfullFitResults_HS[0][f_comp]);
        
            tmpTemplatesProjections[2][f_comp]->SetBinContent(idx,tmpTemplatesProjections[2][f_comp]->GetBinContent(idx)*resfullFitResults_HS[1][f_comp]);
            tmpTemplatesProjections[3][f_comp]->SetBinContent(idx,tmpTemplatesProjections[3][f_comp]->GetBinContent(idx)*resfullFitResults_HS[1][f_comp]);
        
            tmpTemplatesProjections[4][f_comp]->SetBinContent(idx,tmpTemplatesProjections[4][f_comp]->GetBinContent(idx)*resfullFitResults_HS[2][f_comp]);
            tmpTemplatesProjections[5][f_comp]->SetBinContent(idx,tmpTemplatesProjections[5][f_comp]->GetBinContent(idx)*resfullFitResults_HS[2][f_comp]);
        
            tmpTemplatesProjections[6][f_comp]->SetBinContent(idx,tmpTemplatesProjections[6][f_comp]->GetBinContent(idx)*resfullFitResults_HS[3][f_comp]);
            tmpTemplatesProjections[7][f_comp]->SetBinContent(idx,tmpTemplatesProjections[7][f_comp]->GetBinContent(idx)*resfullFitResults_HS[3][f_comp]);
        
        }
    
        tmpTemplatesProjections[0][f_comp]->SetFillColor(kRed);
        tmpTemplatesProjections[1][f_comp]->SetFillColor(kRed);
    
        tmpTemplatesProjections[2][f_comp]->SetFillColor(kGreen);
        tmpTemplatesProjections[3][f_comp]->SetFillColor(kGreen);
    
        tmpTemplatesProjections[4][f_comp]->SetFillColor(kViolet);
        tmpTemplatesProjections[5][f_comp]->SetFillColor(kViolet);
    
        tmpTemplatesProjections[6][f_comp]->SetFillColor(kYellow);
        tmpTemplatesProjections[7][f_comp]->SetFillColor(kYellow);
    
        HStack_Data[2*f_comp].Add(tmpTemplatesProjections[0][f_comp]);
        HStack_Data[2*f_comp].Add(tmpTemplatesProjections[2][f_comp]);
        HStack_Data[2*f_comp].Add(tmpTemplatesProjections[4][f_comp]);
        HStack_Data[2*f_comp].Add(tmpTemplatesProjections[6][f_comp]);
    
        HStack_Data[2*f_comp+1].Add(tmpTemplatesProjections[1][f_comp]);
        HStack_Data[2*f_comp+1].Add(tmpTemplatesProjections[3][f_comp]);
        HStack_Data[2*f_comp+1].Add(tmpTemplatesProjections[5][f_comp]);
        HStack_Data[2*f_comp+1].Add(tmpTemplatesProjections[7][f_comp]);
        
    }
    
    if(write_tmp_histos)
    {
    
        //////////////////////////////////// Creating Projections out file
        
        TFile stack_file(stacks_out_path.c_str(),"RECREATE");
        if(stack_file.IsZombie()) {
            std::cout << "\n\nError writing ROOT Stack TFile. Prorgram finished \n\n";
            log_file << "\n\nError writing ROOT Stack TFile. Prorgram finished \n\n";
            exit(100);
        }
        
        //////////////////////////////////// Writing Projections
    
        for(Int_t idx=0; idx < HStack_Data.size(); ++idx)
            HStack_Data[idx].Write();
    
        stack_file.Close();
    
    }
    
}

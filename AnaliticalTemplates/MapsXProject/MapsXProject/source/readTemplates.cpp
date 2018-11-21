
#include "MyHead.h"

void read_templates(
                        TH2D &Template_Iso_LS,
                        TH2D &Template_AniNS_LS,
                        TH2D &Template_AniEW_LS,
                        TH2D &Template_AniFB_LS,
                        TH2D &Template_Iso_HS,
                        TH2D &Template_AniNS_HS,
                        TH2D &Template_AniEW_HS,
                        TH2D &Template_AniFB_HS,
                        std::ofstream &output_log_file,
                        bool DAMPE,
                        std::string template_out_path,
                        std::string DAMPE_template_out_path
                    )

{
    
    std::string path_to_templates;
    if(!DAMPE)
        path_to_templates = template_out_path;
    else
        path_to_templates = DAMPE_template_out_path;
    
    TH2D* tmpTemplate_Iso_LS = nullptr;
    TH2D* tmpTemplate_AniNS_LS = nullptr;
    TH2D* tmpTemplate_AniEW_LS = nullptr;
    TH2D* tmpTemplate_AniFB_LS = nullptr;
    TH2D* tmpTemplate_Iso_HS = nullptr;
    TH2D* tmpTemplate_AniNS_HS = nullptr;
    TH2D* tmpTemplate_AniEW_HS = nullptr;
    TH2D* tmpTemplate_AniFB_HS = nullptr;
    
    TFile inTemplates(path_to_templates.c_str(),"READ");
    if(inTemplates.IsZombie())
    {
        std::cerr << "\n\nError opening templates ROOT file \n\n";
        output_log_file << "\n\nError opening templates ROOT file \n\n";
        exit(100);
    }
    
    if(DAMPE)
    {
        inTemplates.GetObject("DAMPE_Template_Iso_LS",tmpTemplate_Iso_LS);
        inTemplates.GetObject("DAMPE_Template_AniNS_LS",tmpTemplate_AniNS_LS);
        inTemplates.GetObject("DAMPE_Template_AniEW_LS",tmpTemplate_AniEW_LS);
        inTemplates.GetObject("DAMPE_Template_AniFB_LS",tmpTemplate_AniFB_LS);
        
        inTemplates.GetObject("DAMPE_Template_Iso_HS",tmpTemplate_Iso_HS);
        inTemplates.GetObject("DAMPE_Template_AniNS_HS",tmpTemplate_AniNS_HS);
        inTemplates.GetObject("DAMPE_Template_AniEW_HS",tmpTemplate_AniEW_HS);
        inTemplates.GetObject("DAMPE_Template_AniFB_HS",tmpTemplate_AniFB_HS);
    }
    else
    {
        
        inTemplates.GetObject("Template_Iso_LS",tmpTemplate_Iso_LS);
        inTemplates.GetObject("Template_AniNS_LS",tmpTemplate_AniNS_LS);
        inTemplates.GetObject("Template_AniEW_LS",tmpTemplate_AniEW_LS);
        inTemplates.GetObject("Template_AniFB_LS",tmpTemplate_AniFB_LS);
        
        inTemplates.GetObject("Template_Iso_HS",tmpTemplate_Iso_HS);
        inTemplates.GetObject("Template_AniNS_HS",tmpTemplate_AniNS_HS);
        inTemplates.GetObject("Template_AniEW_HS",tmpTemplate_AniEW_HS);
        inTemplates.GetObject("Template_AniFB_HS",tmpTemplate_AniFB_HS);
        
    }
    
    /*
     tmpTemplate_Iso_LS = (TH2D*) inTemplates.Get("Template_Iso_LS");
     tmpTemplate_AniNS_LS = (TH2D*) inTemplates.Get("Template_AniNS_LS");
     tmpTemplate_AniEW_LS = (TH2D*) inTemplates.Get("Template_AniEW_LS");
     tmpTemplate_AniFB_LS = (TH2D*) inTemplates.Get("Template_AniFB_LS");
     
     tmpTemplate_Iso_HS = (TH2D*) inTemplates.Get("Template_Iso_HS");
     tmpTemplate_AniNS_HS = (TH2D*) inTemplates.Get("Template_AniNS_HS");
     tmpTemplate_AniEW_HS = (TH2D*) inTemplates.Get("Template_AniEW_HS");
     tmpTemplate_AniFB_HS = (TH2D*) inTemplates.Get("Template_AniFB_HS");
     */
    
    if(!DAMPE)
    {
        new (&Template_Iso_LS) (TH2D) (*(TH2D*)tmpTemplate_Iso_LS->Clone("Template_Iso_LS"));
        new (&Template_AniNS_LS) (TH2D) (*(TH2D*)tmpTemplate_AniNS_LS->Clone("Template_AniNS_LS"));
        new (&Template_AniEW_LS) (TH2D) (*(TH2D*)tmpTemplate_AniEW_LS->Clone("Template_AniEW_LS"));
        new (&Template_AniFB_LS) (TH2D) (*(TH2D*)tmpTemplate_AniFB_LS->Clone("Template_AniFB_LS"));
        
        new (&Template_Iso_HS) (TH2D) (*(TH2D*)tmpTemplate_Iso_HS->Clone("Template_Iso_HS"));
        new (&Template_AniNS_HS) (TH2D) (*(TH2D*)tmpTemplate_AniNS_HS->Clone("Template_AniNS_HS"));
        new (&Template_AniEW_HS) (TH2D) (*(TH2D*)tmpTemplate_AniEW_HS->Clone("Template_AniEW_HS"));
        new (&Template_AniFB_HS) (TH2D) (*(TH2D*)tmpTemplate_AniFB_HS->Clone("Template_AniFB_HS"));
    }
    else
    {
        new (&Template_Iso_LS) (TH2D) (*(TH2D*)tmpTemplate_Iso_LS->Clone("DAMPE_Template_Iso_LS"));
        new (&Template_AniNS_LS) (TH2D) (*(TH2D*)tmpTemplate_AniNS_LS->Clone("DAMPE_Template_AniNS_LS"));
        new (&Template_AniEW_LS) (TH2D) (*(TH2D*)tmpTemplate_AniEW_LS->Clone("DAMPE_Template_AniEW_LS"));
        new (&Template_AniFB_LS) (TH2D) (*(TH2D*)tmpTemplate_AniFB_LS->Clone("DAMPE_Template_AniFB_LS"));
        
        new (&Template_Iso_HS) (TH2D) (*(TH2D*)tmpTemplate_Iso_HS->Clone("DAMPE_Template_Iso_HS"));
        new (&Template_AniNS_HS) (TH2D) (*(TH2D*)tmpTemplate_AniNS_HS->Clone("DAMPE_Template_AniNS_HS"));
        new (&Template_AniEW_HS) (TH2D) (*(TH2D*)tmpTemplate_AniEW_HS->Clone("DAMPE_Template_AniEW_HS"));
        new (&Template_AniFB_HS) (TH2D) (*(TH2D*)tmpTemplate_AniFB_HS->Clone("DAMPE_Template_AniFB_HS"));
    }
    
    inTemplates.Close();
    
}


void read_relative_templates(
                                TH2D &relative_DAMPE_Template_Iso_LS,
                                TH2D &relative_DAMPE_Template_AniNS_LS,
                                TH2D &relative_DAMPE_Template_AniEW_LS,
                                TH2D &relative_DAMPE_Template_AniFB_LS,
                                TH2D &relative_DAMPE_Template_Iso_HS,
                                TH2D &relative_DAMPE_Template_AniNS_HS,
                                TH2D &relative_DAMPE_Template_AniEW_HS,
                                TH2D &relative_DAMPE_Template_AniFB_HS,
                                std::ofstream &output_log_file,
                                std::string DAMPE_template_out_path
                             )

{
    
    TH2D* tmpTemplate_Iso_LS = nullptr;
    TH2D* tmpTemplate_AniNS_LS = nullptr;
    TH2D* tmpTemplate_AniEW_LS = nullptr;
    TH2D* tmpTemplate_AniFB_LS = nullptr;
    TH2D* tmpTemplate_Iso_HS = nullptr;
    TH2D* tmpTemplate_AniNS_HS = nullptr;
    TH2D* tmpTemplate_AniEW_HS = nullptr;
    TH2D* tmpTemplate_AniFB_HS = nullptr;
    
    
    TFile inTemplates(DAMPE_template_out_path.c_str(),"READ");
    if(inTemplates.IsZombie())
    {
        std::cerr << "\n\nError opening relative templates ROOT file \n\n";
        output_log_file << "\n\nError opening relative templates ROOT file \n\n";
        exit(100);
    }
    
    inTemplates.GetObject("relative_DAMPE_Template_Iso_LS",tmpTemplate_Iso_LS);
    inTemplates.GetObject("relative_DAMPE_Template_AniNS_LS",tmpTemplate_AniNS_LS);
    inTemplates.GetObject("relative_DAMPE_Template_AniEW_LS",tmpTemplate_AniEW_LS);
    inTemplates.GetObject("relative_DAMPE_Template_AniFB_LS",tmpTemplate_AniFB_LS);
    
    inTemplates.GetObject("relative_DAMPE_Template_Iso_HS",tmpTemplate_Iso_HS);
    inTemplates.GetObject("relative_DAMPE_Template_AniNS_HS",tmpTemplate_AniNS_HS);
    inTemplates.GetObject("relative_DAMPE_Template_AniEW_HS",tmpTemplate_AniEW_HS);
    inTemplates.GetObject("relative_DAMPE_Template_AniFB_HS",tmpTemplate_AniFB_HS);
    
    
    new (&relative_DAMPE_Template_Iso_LS) (TH2D) (*(TH2D*)tmpTemplate_Iso_LS->Clone("relative_DAMPE_Template_Iso_LS"));
    new (&relative_DAMPE_Template_AniNS_LS) (TH2D) (*(TH2D*)tmpTemplate_AniNS_LS->Clone("relative_DAMPE_Template_AniNS_LS"));
    new (&relative_DAMPE_Template_AniEW_LS) (TH2D) (*(TH2D*)tmpTemplate_AniEW_LS->Clone("relative_DAMPE_Template_AniEW_LS"));
    new (&relative_DAMPE_Template_AniFB_LS) (TH2D) (*(TH2D*)tmpTemplate_AniFB_LS->Clone("relative_DAMPE_Template_AniFB_LS"));
    
    new (&relative_DAMPE_Template_Iso_HS) (TH2D) (*(TH2D*)tmpTemplate_Iso_HS->Clone("relative_DAMPE_Template_Iso_HS"));
    new (&relative_DAMPE_Template_AniNS_HS) (TH2D) (*(TH2D*)tmpTemplate_AniNS_HS->Clone("relative_DAMPE_Template_AniNS_HS"));
    new (&relative_DAMPE_Template_AniEW_HS) (TH2D) (*(TH2D*)tmpTemplate_AniEW_HS->Clone("relative_DAMPE_Template_AniEW_HS"));
    new (&relative_DAMPE_Template_AniFB_HS) (TH2D) (*(TH2D*)tmpTemplate_AniFB_HS->Clone("relative_DAMPE_Template_AniFB_HS"));
    
    inTemplates.Close();
    
}


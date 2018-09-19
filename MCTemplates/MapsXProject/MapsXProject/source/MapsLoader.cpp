
#include "MyHead.h"

void load_1D_histos(TH1D Templates_LS[],TH1D &DataHisto_I_LS,TH1D &DataHisto_NS_LS,TH1D &DataHisto_EW_LS,TH1D &DataHisto_FB_LS,TH1D &MixedDataHisto_NS_EW_LS,TH1D &MixedDataHisto_NS_FB_LS,TH1D &MixedDataHisto_EW_FB_LS,TH1D &FullMixedDataHisto_LS,TH1D Templates_HS[],TH1D &DataHisto_I_HS,TH1D &DataHisto_NS_HS,TH1D &DataHisto_EW_HS,TH1D &DataHisto_FB_HS,TH1D &MixedDataHisto_NS_EW_HS,TH1D &MixedDataHisto_NS_FB_HS,TH1D &MixedDataHisto_EW_FB_HS,TH1D &FullMixedDataHisto_HS,std::string template_path,std::string data_path,std::ofstream &log_file) {
    
    TH1D tmpTemplates_LS[4];
    TH1D tmpTemplates_HS[4];
    
    TH1D tmpDataHisto_I_LS;
    TH1D tmpDataHisto_NS_LS;
    TH1D tmpDataHisto_EW_LS;
    TH1D tmpDataHisto_FB_LS;
    TH1D tmpMixedDataHisto_NS_EW_LS;
    TH1D tmpMixedDataHisto_NS_FB_LS;
    TH1D tmpMixedDataHisto_EW_FB_LS;
    TH1D tmpFullMixedDataHisto_LS;
    
    TH1D tmpDataHisto_I_HS;
    TH1D tmpDataHisto_NS_HS;
    TH1D tmpDataHisto_EW_HS;
    TH1D tmpDataHisto_FB_HS;
    TH1D tmpMixedDataHisto_NS_EW_HS;
    TH1D tmpMixedDataHisto_NS_FB_HS;
    TH1D tmpMixedDataHisto_EW_FB_HS;
    TH1D tmpFullMixedDataHisto_HS;
    
    TH2D Template_Iso_LS;
    TH2D Template_AniNS_LS;
    TH2D Template_AniEW_LS;
    TH2D Template_AniFB_LS;
    
    TH2D Template_Iso_HS;
    TH2D Template_AniNS_HS;
    TH2D Template_AniEW_HS;
    TH2D Template_AniFB_HS;
    
    TH2D Data_Iso_LS;
    TH2D Data_AniNS_LS;
    TH2D Data_AniEW_LS;
    TH2D Data_AniFB_LS;
    TH2D MixedData_NS_EW_LS;
    TH2D MixedData_NS_FB_LS;
    TH2D MixedData_EW_FB_LS;
    TH2D FullMixedData_LS;
    
    TH2D Data_Iso_HS;
    TH2D Data_AniNS_HS;
    TH2D Data_AniEW_HS;
    TH2D Data_AniFB_HS;
    TH2D MixedData_NS_EW_HS;
    TH2D MixedData_NS_FB_HS;
    TH2D MixedData_EW_FB_HS;
    TH2D FullMixedData_HS;
    
    
    ////////////// Read histos from file
    
    TFile TemplateFile(template_path.c_str(),"READ");
    if(TemplateFile.IsZombie()) {
        std::cout << "\n\nError reading input template file! \n\n";
        log_file << "\n\nError reading input template file! \n\n";
        exit(-2);
    }
    
    new (&Template_Iso_LS) (TH2D)(*(TH2D*)TemplateFile.Get("Template_Iso_LS"));
    new (&Template_AniNS_LS) (TH2D)(*(TH2D*)TemplateFile.Get("Template_AniNS_LS"));
    new (&Template_AniEW_LS) (TH2D)(*(TH2D*)TemplateFile.Get("Template_AniEW_LS"));
    new (&Template_AniFB_LS) (TH2D)(*(TH2D*)TemplateFile.Get("Template_AniFB_LS"));
    
    new (&Template_Iso_HS) (TH2D)(*(TH2D*)TemplateFile.Get("Template_Iso_HS"));
    new (&Template_AniNS_HS) (TH2D)(*(TH2D*)TemplateFile.Get("Template_AniNS_HS"));
    new (&Template_AniEW_HS) (TH2D)(*(TH2D*)TemplateFile.Get("Template_AniEW_HS"));
    new (&Template_AniFB_HS) (TH2D)(*(TH2D*)TemplateFile.Get("Template_AniFB_HS"));

    TemplateFile.Close();
    
    TFile DataFile(data_path.c_str(),"READ");
    if(DataFile.IsZombie()) {
        std::cout << "\n\nError reading input data file! \n\n";
        log_file << "\n\nError reading input data file! \n\n";
        exit(-2);
    }
    
    new (&Data_Iso_LS) (TH2D)(*(TH2D*)DataFile.Get("Data_Iso_LS"));
    new (&Data_AniNS_LS) (TH2D)(*(TH2D*)DataFile.Get("Data_AniNS_LS"));
    new (&Data_AniEW_LS) (TH2D)(*(TH2D*)DataFile.Get("Data_AniEW_LS"));
    new (&Data_AniFB_LS) (TH2D)(*(TH2D*)DataFile.Get("Data_AniFB_LS"));
    new (&MixedData_NS_EW_LS) (TH2D)(*(TH2D*)DataFile.Get("MixedData_NS_EW_LS"));
    new (&MixedData_NS_FB_LS) (TH2D)(*(TH2D*)DataFile.Get("MixedData_NS_FB_LS"));
    new (&MixedData_EW_FB_LS) (TH2D)(*(TH2D*)DataFile.Get("MixedData_EW_FB_LS"));
    new (&FullMixedData_LS) (TH2D)(*(TH2D*)DataFile.Get("FullMixedData_LS"));
    
    new (&Data_Iso_HS) (TH2D)(*(TH2D*)DataFile.Get("Data_Iso_HS"));
    new (&Data_AniNS_HS) (TH2D)(*(TH2D*)DataFile.Get("Data_AniNS_HS"));
    new (&Data_AniEW_HS) (TH2D)(*(TH2D*)DataFile.Get("Data_AniEW_HS"));
    new (&Data_AniFB_HS) (TH2D)(*(TH2D*)DataFile.Get("Data_AniFB_HS"));
    new (&MixedData_NS_EW_HS) (TH2D)(*(TH2D*)DataFile.Get("MixedData_NS_EW_HS"));
    new (&MixedData_NS_FB_HS) (TH2D)(*(TH2D*)DataFile.Get("MixedData_NS_FB_HS"));
    new (&MixedData_EW_FB_HS) (TH2D)(*(TH2D*)DataFile.Get("MixedData_EW_FB_HS"));
    new (&FullMixedData_HS) (TH2D)(*(TH2D*)DataFile.Get("FullMixedData_HS"));
    
    DataFile.Close();
    
    
#if 0
    ///////////////////////// Converting Templates Maps from TH2 to TH1
    
    TH2toTH1(tmpTemplates_LS[0],Template_Iso_LS);
    TH2toTH1(tmpTemplates_LS[1],Template_AniNS_LS);
    TH2toTH1(tmpTemplates_LS[2],Template_AniEW_LS);
    TH2toTH1(tmpTemplates_LS[3],Template_AniFB_LS);
    
    TH2toTH1(tmpTemplates_HS[0],Template_Iso_HS);
    TH2toTH1(tmpTemplates_HS[1],Template_AniNS_HS);
    TH2toTH1(tmpTemplates_HS[2],Template_AniEW_HS);
    TH2toTH1(tmpTemplates_HS[3],Template_AniFB_HS);
    
    ///////////////////////// Converting Data Maps from TH2 to TH1
    
    TH2toTH1(tmpDataHisto_I_LS,Data_Iso_LS);
    TH2toTH1(tmpDataHisto_NS_LS,Data_AniNS_LS);
    TH2toTH1(tmpDataHisto_EW_LS,Data_AniEW_LS);
    TH2toTH1(tmpDataHisto_FB_LS,Data_AniFB_LS);
    TH2toTH1(tmpMixedDataHisto_NS_EW_LS,MixedData_NS_EW_LS);
    TH2toTH1(tmpMixedDataHisto_NS_FB_LS,MixedData_NS_FB_LS);
    TH2toTH1(tmpMixedDataHisto_EW_FB_LS,MixedData_EW_FB_LS);
    TH2toTH1(tmpFullMixedDataHisto_LS,FullMixedData_LS);
    
    TH2toTH1(tmpDataHisto_I_HS,Data_Iso_HS);
    TH2toTH1(tmpDataHisto_NS_HS,Data_AniNS_HS);
    TH2toTH1(tmpDataHisto_EW_HS,Data_AniEW_HS);
    TH2toTH1(tmpDataHisto_FB_HS,Data_AniFB_HS);
    TH2toTH1(tmpMixedDataHisto_NS_EW_HS,MixedData_NS_EW_HS);
    TH2toTH1(tmpMixedDataHisto_NS_FB_HS,MixedData_NS_FB_HS);
    TH2toTH1(tmpMixedDataHisto_EW_FB_HS,MixedData_EW_FB_HS);
    TH2toTH1(tmpFullMixedDataHisto_HS,FullMixedData_HS);

    /////////////////////////Linking to "Main TH1D's
    
    for(Int_t t_idx = 0; t_idx < 4; t_idx++) {
        new (&Templates_LS[t_idx]) (TH1D)(*(TH1D*)tmpTemplates_LS[0].Clone());
        new (&Templates_HS[t_idx]) (TH1D)(*(TH1D*)tmpTemplates_HS[0].Clone());
    }
    
    new (&DataHisto_I_LS) (TH1D)(*(TH1D*)tmpDataHisto_I_LS.Clone());
    new (&DataHisto_NS_LS) (TH1D)(*(TH1D*)tmpDataHisto_NS_LS.Clone());
    new (&DataHisto_EW_LS) (TH1D)(*(TH1D*)tmpDataHisto_EW_LS.Clone());
    new (&DataHisto_FB_LS) (TH1D)(*(TH1D*)tmpDataHisto_FB_LS.Clone());
    new (&MixedDataHisto_NS_EW_LS) (TH1D)(*(TH1D*)tmpMixedDataHisto_NS_EW_LS.Clone());
    new (&MixedDataHisto_NS_FB_LS) (TH1D)(*(TH1D*)tmpMixedDataHisto_NS_FB_LS.Clone());
    new (&MixedDataHisto_EW_FB_LS) (TH1D)(*(TH1D*)tmpMixedDataHisto_EW_FB_LS.Clone());
    new (&FullMixedDataHisto_LS) (TH1D)(*(TH1D*)tmpFullMixedDataHisto_LS.Clone());
    
    new (&DataHisto_I_HS) (TH1D)(*(TH1D*)tmpDataHisto_I_HS.Clone());
    new (&DataHisto_NS_HS) (TH1D)(*(TH1D*)tmpDataHisto_NS_HS.Clone());
    new (&DataHisto_EW_HS) (TH1D)(*(TH1D*)tmpDataHisto_EW_HS.Clone());
    new (&DataHisto_FB_HS) (TH1D)(*(TH1D*)tmpDataHisto_FB_HS.Clone());
    new (&MixedDataHisto_NS_EW_HS) (TH1D)(*(TH1D*)tmpMixedDataHisto_NS_EW_HS.Clone());
    new (&MixedDataHisto_NS_FB_HS) (TH1D)(*(TH1D*)tmpMixedDataHisto_NS_FB_HS.Clone());
    new (&MixedDataHisto_EW_FB_HS) (TH1D)(*(TH1D*)tmpMixedDataHisto_EW_FB_HS.Clone());
    new (&FullMixedDataHisto_HS) (TH1D)(*(TH1D*)tmpFullMixedDataHisto_HS.Clone());
    
    /////////////////////////////
    
#else
    
    ///////////////////////// Converting Templates Maps from TH2 to TH1
    
    TH2toTH1(Templates_LS[0],Template_Iso_LS);
    TH2toTH1(Templates_LS[1],Template_AniNS_LS);
    TH2toTH1(Templates_LS[2],Template_AniEW_LS);
    TH2toTH1(Templates_LS[3],Template_AniFB_LS);
    
    TH2toTH1(Templates_HS[0],Template_Iso_HS);
    TH2toTH1(Templates_HS[1],Template_AniNS_HS);
    TH2toTH1(Templates_HS[2],Template_AniEW_HS);
    TH2toTH1(Templates_HS[3],Template_AniFB_HS);
    
    ///////////////////////// Converting Data Maps from TH2 to TH1
    
    TH2toTH1(DataHisto_I_LS,Data_Iso_LS);
    TH2toTH1(DataHisto_NS_LS,Data_AniNS_LS);
    TH2toTH1(DataHisto_EW_LS,Data_AniEW_LS);
    TH2toTH1(DataHisto_FB_LS,Data_AniFB_LS);
    TH2toTH1(MixedDataHisto_NS_EW_LS,MixedData_NS_EW_LS);
    TH2toTH1(MixedDataHisto_NS_FB_LS,MixedData_NS_FB_LS);
    TH2toTH1(MixedDataHisto_EW_FB_LS,MixedData_EW_FB_LS);
    TH2toTH1(FullMixedDataHisto_LS,FullMixedData_LS);
    
    TH2toTH1(DataHisto_I_HS,Data_Iso_HS);
    TH2toTH1(DataHisto_NS_HS,Data_AniNS_HS);
    TH2toTH1(DataHisto_EW_HS,Data_AniEW_HS);
    TH2toTH1(DataHisto_FB_HS,Data_AniFB_HS);
    TH2toTH1(MixedDataHisto_NS_EW_HS,MixedData_NS_EW_HS);
    TH2toTH1(MixedDataHisto_NS_FB_HS,MixedData_NS_FB_HS);
    TH2toTH1(MixedDataHisto_EW_FB_HS,MixedData_EW_FB_HS);
    TH2toTH1(FullMixedDataHisto_HS,FullMixedData_HS);

#endif
}

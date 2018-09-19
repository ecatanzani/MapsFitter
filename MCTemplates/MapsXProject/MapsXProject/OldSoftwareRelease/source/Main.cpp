
#include "MyHead.h"

int main(int argc,char *argv[]) {

    if(argc!=1) {
        
        std::stringstream ss;
        std::string input_path;
        
        ss << argv[1];
        ss >> input_path;
        
        //////////// Histos declaration
        
        //// All Sky Maps
        
        TH2D AllSky_IsoMap,AllSky_AniMap_NS,AllSky_AniMap_EW,AllSky_AniMap_FB;
        TH2D NAllSky_IsoMap,NAllSky_AniMap_NS,NAllSky_AniMap_EW,NAllSky_AniMap_FB;
        TH2D AllSky_Isolate_NS,AllSky_Isolate_EW,AllSky_Isolate_FB;
        TH2D MapsIsoRaio_NS,MapsIsoRaio_EW,MapsIsoRaio_FB;
        
        //// Satellite's Maps
        
        TH2D IsoMap;
        TH2D AniMap_NS,AniMap_EW,AniMap_FB;
        TH2D phAniMap_NS,phAniMap_EW,phAniMap_FB;
        TH2D NIsoMap,NAniMap_NS,NAniMap_EW,NAniMap_FB;
        TH2D phNAniMap_NS,phNAniMap_EW,phNAniMap_FB;
        TH2D Isolate_NS,Isolate_EW,Isolate_FB;
        TH2D Isolate_NN_NS,Isolate_NN_EW,Isolate_NN_FB;
        TH2D RatioAniMap_NS,RatioAniMap_EW,RatioAniMap_FB;
    
        TH1D h_wNS,h_wEW,h_wFB;
        TH1D ph_h_wNS,ph_h_wEW,ph_h_wFB;
        
        //////////////////////////////////// Importing histos....
        
        TFile inFile(input_path.c_str(),"READ");
        if(inFile.IsZombie()) {
            std::cout<<"\n\nERROR ! Could not open input TFIle\n\n";
            exit(-2);
        }
        
        new (&AllSky_IsoMap) TH2D(*(TH2D*)inFile.Get("AllSky_IsoMap"));
        new (&AllSky_AniMap_NS) TH2D(*(TH2D*)inFile.Get("AllSky_AniMap_NS"));
        new (&AllSky_AniMap_EW) TH2D(*(TH2D*)inFile.Get("AllSky_AniMap_EW"));
        new (&AllSky_AniMap_FB) TH2D(*(TH2D*)inFile.Get("AllSky_AniMap_FB"));
        new (&NAllSky_IsoMap) TH2D(*(TH2D*)inFile.Get("NAllSky_IsoMap"));
        new (&NAllSky_AniMap_NS) TH2D(*(TH2D*)inFile.Get("NAllSky_AniMap_NS"));
        new (&NAllSky_AniMap_EW) TH2D(*(TH2D*)inFile.Get("NAllSky_AniMap_EW"));
        new (&NAllSky_AniMap_FB) TH2D(*(TH2D*)inFile.Get("NAllSky_AniMap_FB"));
        new (&AllSky_Isolate_NS) TH2D(*(TH2D*)inFile.Get("AllSky_Isolate_NS"));
        new (&AllSky_Isolate_EW) TH2D(*(TH2D*)inFile.Get("AllSky_Isolate_EW"));
        new (&AllSky_Isolate_FB) TH2D(*(TH2D*)inFile.Get("AllSky_Isolate_FB"));
        new (&MapsIsoRaio_NS) TH2D(*(TH2D*)inFile.Get("MapsIsoRaio_NS"));
        new (&MapsIsoRaio_EW) TH2D(*(TH2D*)inFile.Get("MapsIsoRaio_EW"));
        new (&MapsIsoRaio_FB) TH2D(*(TH2D*)inFile.Get("MapsIsoRaio_FB"));
        
        new (&IsoMap) TH2D(*(TH2D*)inFile.Get("IsoMap"));
        new (&AniMap_NS) TH2D(*(TH2D*)inFile.Get("AniMap_NS"));
        new (&AniMap_EW) TH2D(*(TH2D*)inFile.Get("AniMap_EW"));
        new (&AniMap_FB) TH2D(*(TH2D*)inFile.Get("AniMap_FB"));
        new (&phAniMap_NS) TH2D(*(TH2D*)inFile.Get("phAniMap_NS"));
        new (&phAniMap_EW) TH2D(*(TH2D*)inFile.Get("phAniMap_EW"));
        new (&phAniMap_FB) TH2D(*(TH2D*)inFile.Get("phAniMap_FB"));
        new (&NIsoMap) TH2D(*(TH2D*)inFile.Get("NIsoMap"));
        new (&NAniMap_NS) TH2D(*(TH2D*)inFile.Get("NAniMap_NS"));
        new (&NAniMap_EW) TH2D(*(TH2D*)inFile.Get("NAniMap_EW"));
        new (&NAniMap_FB) TH2D(*(TH2D*)inFile.Get("NAniMap_FB"));
        new (&phNAniMap_NS) TH2D(*(TH2D*)inFile.Get("phNAniMap_NS"));
        new (&phNAniMap_EW) TH2D(*(TH2D*)inFile.Get("phNAniMap_EW"));
        new (&phNAniMap_FB) TH2D(*(TH2D*)inFile.Get("phNAniMap_FB"));
        new (&Isolate_NS) TH2D(*(TH2D*)inFile.Get("Isolate_NS"));
        new (&Isolate_EW) TH2D(*(TH2D*)inFile.Get("Isolate_EW"));
        new (&Isolate_FB) TH2D(*(TH2D*)inFile.Get("Isolate_FB"));
        new (&Isolate_NN_NS) TH2D(*(TH2D*)inFile.Get("Isolate_NN_NS"));
        new (&Isolate_NN_EW) TH2D(*(TH2D*)inFile.Get("Isolate_NN_EW"));
        new (&Isolate_NN_FB) TH2D(*(TH2D*)inFile.Get("Isolate_NN_FB"));
        new (&RatioAniMap_NS) TH2D(*(TH2D*)inFile.Get("RatioAniMap_NS"));
        new (&RatioAniMap_EW) TH2D(*(TH2D*)inFile.Get("RatioAniMap_EW"));
        new (&RatioAniMap_FB) TH2D(*(TH2D*)inFile.Get("RatioAniMap_FB"));
        
        new (&h_wNS) (TH1D)(*(TH1D*)inFile.Get("h_wNS"));
        new (&h_wEW) (TH1D)(*(TH1D*)inFile.Get("h_wEW"));
        new (&h_wFB) (TH1D)(*(TH1D*)inFile.Get("h_wFB"));
        new (&ph_h_wNS) (TH1D)(*(TH1D*)inFile.Get("ph_h_wNS"));
        new (&ph_h_wEW) (TH1D)(*(TH1D*)inFile.Get("ph_h_wEW"));
        new (&ph_h_wFB) (TH1D)(*(TH1D*)inFile.Get("ph_h_wFB"));
        
        inFile.Close();
        
        /////////////////////////////////////// Finally making fits !!!
        
        std::cout<< "\n\n --------------------------------------- All Sky Fits --------------------------------------- \n\n";
        
        //Chi2Fit(MapsIsoRaio_NS,MapsIsoRaio_EW,MapsIsoRaio_FB);
        //AllSky_MapFit(AllSky_Isolate_NS,AllSky_Isolate_EW,AllSky_Isolate_FB);
        
        std::cout<< "\n\n -------------------------------------------------------------------------------------------- \n\n";
        
        std::cout<<"\n\n Normalized maps fit";
        std::cout<<"\n -------------------------------------------------- \n\n";
        //fit_map(Isolate_NS,Isolate_EW,Isolate_FB);
        
        std::cout<<"\n\n Maps fit";
        std::cout<<"\n -------------------------------------------------- \n\n";
        //fit_NN_map(Isolate_NN_NS,Isolate_NN_EW,Isolate_NN_FB);
        
        
    }
    else {
        
        ////////////////////////////////////////////////////////// Variables declaration ///////////////////////////////////////////////////
    
        //////////// Costheta flat binning variables
    
        Int_t n_bin_lon = 36;                                                                      // -> Number of bins along longitude axis
        Double_t lon_bin_min = 0;                                                              // -> Set max and min for longitude binning
        Double_t lon_bin_max = 360;
    
        Int_t n_bin_lat = 18;                                                                      // -> Number of bins along latitude axis
        Double_t lat_bin_min = 0;                                                               // -> Set max and min for latitude binning
        Double_t lat_bin_max = 180;
        
        Double_t* binning = nullptr;                                                                // -> Array used to store the custom binning intervals !!!
                                                                                                    // -> Number of events into the isotripic loaded map
    
        //////////// Stuff variables
        
        std::string log_path = output_path_creator(0),root_out_path = output_path_creator(1);       // -> Log end ROOT result file variables
        std::ofstream output_log_file(log_path);                                                    // -> log file creation !
        
        //////////// TemplateFit variables
        
        TH2D* Templates[4] = {nullptr,nullptr,nullptr,nullptr};
        TH2D* DataHisto_I = nullptr;
        TH2D* DataHisto_NS = nullptr;
        TH2D* DataHisto_EW = nullptr;
        TH2D* DataHisto_FB = nullptr;
        
        Double_t res[4],res_err[4];
        Double_t initialValues[4]={1,0,0,0};
        bool fix[4]={false,false,false,false};
        
        
        //////////// Histos
        
        create_binning(n_bin_lat,lat_bin_min,lat_bin_max,binning,true);                             // -> Custom binning for Healpix representation
        
        //////////// All Sky Histos
        
        TH2D TemplateIso("TemplateIso","Isotropic Template All Sky Map; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon,lon_bin_min,lon_bin_max,n_bin_lat,binning);
        TH2D TemplateAniNS("TemplateAniNS","Anisotropic Template All Sky Map (NS); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon,lon_bin_min,lon_bin_max,n_bin_lat,binning);
        TH2D TemplateAniEW("TemplateAniEW","Anisotropic Template All Sky Map (EW); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon,lon_bin_min,lon_bin_max,n_bin_lat,binning);
        TH2D TemplateAniFB("TemplateAniFB","Anisotropic Template All Sky Map (FB); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon,lon_bin_min,lon_bin_max,n_bin_lat,binning);
        
        TH2D dataI("dataI","Isotropic All Sky (data) Map; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon,lon_bin_min,lon_bin_max,n_bin_lat,binning);
        TH2D dataNS("dataNS","Anisotropic All Sky (data) Map (NS); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon,lon_bin_min,lon_bin_max,n_bin_lat,binning);
        TH2D dataEW("dataEW","Anisotropic All Sky (data) Map (EW); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon,lon_bin_min,lon_bin_max,n_bin_lat,binning);
        TH2D dataFB("dataFB","Anisotropic All Sky (data) Map (FB); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon,lon_bin_min,lon_bin_max,n_bin_lat,binning);
        

        /*
        //////////// Stellite's Histos
        
        TH2D IsoMap,NIsoMap;
        TH2D AniMap_NS("AniMap_NS","Anisotropic Sky Map (NS) ; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon,lon_bin_min,lon_bin_max,n_bin_lat,binning);
        TH2D AniMap_EW("AniMap_EW","Anisotropic Sky Map (EW) ; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon,lon_bin_min,lon_bin_max,n_bin_lat,binning);
        TH2D AniMap_FB("AniMap_FB","Anisotropic Sky Map (FB) ; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon,lon_bin_min,lon_bin_max,n_bin_lat,binning);
        
        TH2D phAniMap_NS("ph_AniMap_NS","Anisotropic Sky Map (NS) ; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon,lon_bin_min,lon_bin_max,n_bin_lat,binning);
        TH2D phAniMap_EW("ph_AniMap_EW","Anisotropic Sky Map (EW) ; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon,lon_bin_min,lon_bin_max,n_bin_lat,binning);
        TH2D phAniMap_FB("ph_AniMap_FB","Anisotropic Sky Map (FB) ; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon,lon_bin_min,lon_bin_max,n_bin_lat,binning);
        
        TH2D NAniMap_NS,NAniMap_EW,NAniMap_FB;
        TH2D phNAniMap_NS,phNAniMap_EW,phNAniMap_FB;
        TH2D Isolate_NS,Isolate_EW,Isolate_FB;
        TH2D Isolate_NN_NS,Isolate_NN_EW,Isolate_NN_FB;
        TH2D RatioAniMap_NS,RatioAniMap_EW,RatioAniMap_FB;
        
        TH1D h_wNS("h_wNS","w_{NS} distribution",100,-2,3);
        TH1D h_wEW("h_wEW","w_{EW} distribution",100,-2,3);
        TH1D h_wFB("h_wFB","w_{FB} distribution",100,-2,3);
        
        TH1D ph_h_wNS("ph_wNS","w_{NS} distribution",100,-2,3);
        TH1D ph_h_wEW("ph_wEW","w_{EW} distribution",100,-2,3);
        TH1D ph_h_wFB("ph_wFB","w_{FB} distribution",100,-2,3);
        */
         
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
        create_and_initialize_log(output_log_file);
        
        ////////////////////////////////////////////////////////////// All sky maps
        
        generate_templates(TemplateIso,TemplateAniNS,TemplateAniEW,TemplateAniFB,output_log_file);
        generate_data(dataI,dataNS,dataEW,dataFB,output_log_file);
        normalize_templates(TemplateIso,TemplateAniNS,TemplateAniEW,TemplateAniFB,output_log_file);
        normalize_data(dataI,dataNS,dataEW,dataFB,output_log_file);
        
        Templates[0] = &TemplateIso;
        Templates[1] = &TemplateAniNS;
        Templates[2] = &TemplateAniEW;
        Templates[3] = &TemplateAniFB;
        
        DataHisto_I = &dataI;
        DataHisto_NS = &dataNS;
        DataHisto_EW = &dataEW;
        DataHisto_FB = &dataFB;
        
        TH1 *Template1D[4];
        TH1D *DataHisto1D;
        DataHisto1D=TH2toTH1( DataHisto_I );
        
        for(int icomp=0; icomp<4; icomp++)
            Template1D[icomp] = TH2toTH1( Templates[icomp] );
        
        TCanvas ctemp("ctemp","ctemp");
        ctemp.Divide(2,2);
        
        for (int icomp=0; icomp<4; icomp++) {
            ctemp.cd(icomp+1);
            Template1D[icomp]->Draw();
        }
        
        ////////////////////////////////////////// Finally Fit !!! ;-(
        
        //TemplateFitRF(DataHisto,4,Templates,res,res_err,initialValues,fix,false,false,false);
        TemplateFitBH(DataHisto1D,4,Template1D,res,res_err,initialValues,false,false,false,true); // ---> !!!! Chi2 fit !! Substitute last "true" with "false" for likelihood fit !
        
        
        /*
        new (&NAllSky_IsoMap) TH2D(*(TH2D*)AllSky_IsoMap.Clone("NAllSky_IsoMap"));
        
        new (&NAllSky_AniMap_NS) TH2D(*(TH2D*)AllSky_AniMap_NS.Clone("NAllSky_AniMap_NS"));
        new (&NAllSky_AniMap_EW) TH2D(*(TH2D*)AllSky_AniMap_EW.Clone("NAllSky_AniMap_EW"));
        new (&NAllSky_AniMap_FB) TH2D(*(TH2D*)AllSky_AniMap_FB.Clone("NAllSky_AniMap_FB"));
        
        NAllSky_IsoMap.SetTitle("Isotropic All Sky Normalized Map");
        NAllSky_AniMap_NS.SetTitle("Anisotropic All Sky Normalized Map (NS)");
        NAllSky_AniMap_EW.SetTitle("Anisotropic All Sky Normalized Map (EW)");
        NAllSky_AniMap_FB.SetTitle("Anisotropic All Sky Normalized Map (FB)");
        
        normalize_map(NAllSky_IsoMap);
        normalize_map(NAllSky_AniMap_NS);
        normalize_map(NAllSky_AniMap_EW);
        normalize_map(NAllSky_AniMap_FB);
        
        isolate_dipoles(NAllSky_IsoMap,NAllSky_AniMap_NS,NAllSky_AniMap_EW,NAllSky_AniMap_FB,AllSky_Isolate_NS,AllSky_Isolate_EW,AllSky_Isolate_FB);
        
        AllSky_Isolate_NS.SetName("AllSky_Isolate_NS");
        AllSky_Isolate_EW.SetName("AllSky_Isolate_EW");
        AllSky_Isolate_FB.SetName("AllSky_Isolate_FB");
        
        AllSky_Isolate_NS.SetTitle("AllSky NS Dipole");
        AllSky_Isolate_EW.SetTitle("AllSky EW Dipole");
        AllSky_Isolate_FB.SetTitle("AllSky FB Dipole");
        
        new (&MapsIsoRaio_NS) TH2D(*(TH2D*)NAllSky_AniMap_NS.Clone("MapsIsoRaio_NS"));
        new (&MapsIsoRaio_EW) TH2D(*(TH2D*)NAllSky_AniMap_EW.Clone("MapsIsoRaio_EW"));
        new (&MapsIsoRaio_FB) TH2D(*(TH2D*)NAllSky_AniMap_FB.Clone("MapsIsoRaio_FB"));
        
        MapsIsoRaio_NS.SetTitle("Maps Ratio NS");
        MapsIsoRaio_EW.SetTitle("Maps Ratio EW");
        MapsIsoRaio_FB.SetTitle("Maps Ratio FB");
        
        pMapsIsoRaio_NS->Divide(pNAllSky_IsoMap);
        pMapsIsoRaio_EW->Divide(pNAllSky_IsoMap);
        pMapsIsoRaio_FB->Divide(pNAllSky_IsoMap);
        
        std::cout<< "\n\n --------------------------------------- All Sky Fits --------------------------------------- \n\n";
        
        Chi2Fit(MapsIsoRaio_NS,MapsIsoRaio_EW,MapsIsoRaio_FB);
        
        //AllSky_MapFit(AllSky_Isolate_NS,AllSky_Isolate_EW,AllSky_Isolate_FB);
        
        std::cout<< "\n\n -------------------------------------------------------------------------------------------- \n\n";
        
        ////////////////////////////////////////////////////////////// Satellite's maps
        
        obtain_IsoMap(IsoMap,output_log_file);
        sky_events = (Int_t)IsoMap.GetEntries();
        gerenate_ani_map(IsoMap,AniMap_NS,AniMap_EW,AniMap_FB,h_wNS,h_wEW,h_wFB,sky_events,output_log_file);
        gerenate_ani_map_phishing(IsoMap,phAniMap_NS,phAniMap_EW,phAniMap_FB,ph_h_wNS,ph_h_wEW,ph_h_wFB,sky_events,output_log_file);
        
        new (&NIsoMap) TH2D(*(TH2D*)IsoMap.Clone("NIsoMap"));
        new (&NAniMap_NS) TH2D(*(TH2D*)AniMap_NS.Clone("NAniMap_NS"));
        new (&NAniMap_EW) TH2D(*(TH2D*)AniMap_EW.Clone("NAniMap_EW"));
        new (&NAniMap_FB) TH2D(*(TH2D*)AniMap_FB.Clone("NAniMap_FB"));
        new (&phNAniMap_NS) TH2D(*(TH2D*)phAniMap_NS.Clone("phNAniMap_NS"));
        new (&phNAniMap_EW) TH2D(*(TH2D*)phAniMap_EW.Clone("phNAniMap_EW"));
        new (&phNAniMap_FB) TH2D(*(TH2D*)phAniMap_FB.Clone("phNAniMap_FB"));
        
        normalize_map(NIsoMap);
        normalize_map(NAniMap_NS);
        normalize_map(NAniMap_EW);
        normalize_map(NAniMap_FB);
        normalize_map(phNAniMap_NS);
        normalize_map(phNAniMap_EW);
        normalize_map(phNAniMap_FB);
        
        maps_ratio(RatioAniMap_NS,NAniMap_NS,phNAniMap_NS,RatioAniMap_EW,NAniMap_EW,phNAniMap_EW,RatioAniMap_FB,NAniMap_FB,phNAniMap_FB,Isolate_NS,Isolate_EW,Isolate_FB,NIsoMap);
        
        isolate_dipoles(IsoMap,AniMap_NS,AniMap_EW,AniMap_FB,Isolate_NN_NS,Isolate_NN_EW,Isolate_NN_FB);
        
        std::cout<<"\n\n Normalized maps fit";
        std::cout<<"\n -------------------------------------------------- \n\n";
        fit_map(Isolate_NS,Isolate_EW,Isolate_FB);
        
        std::cout<<"\n\n Maps fit";
        std::cout<<"\n -------------------------------------------------- \n\n";
        fit_NN_map(Isolate_NN_NS,Isolate_NN_EW,Isolate_NN_FB);
        
        */
         
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
        std::cout<<"\n\nSimulation Completed !\n\n";
        output_log_file << "\n\nSimulation Completed !\n\n";
        
        //////////////////////////////////// Creating out file
        
        TFile results_file(root_out_path.c_str(),"RECREATE");
        if(results_file.IsZombie()) {
            std::cout << "\n\nError writing ROOT results TFile. Prorgram finished \n\n";
            output_log_file << "\n\nError ROOT results TFile. Prorgram finished \n\n";
            exit(-2);
        }
        
        //////////////////////////////////// Writing results
        
        /*
         
        IsoMap.Write();
        NIsoMap.Write();
        
        AniMap_NS.Write();
        AniMap_EW.Write();
        AniMap_FB.Write();
        
        NAniMap_NS.Write();
        NAniMap_EW.Write();
        NAniMap_FB.Write();
        
        h_wNS.Write();
        h_wEW.Write();
        h_wFB.Write();
        
        ph_h_wNS.Write();
        ph_h_wEW.Write();
        ph_h_wFB.Write();
        
        phAniMap_NS.Write();
        phAniMap_EW.Write();
        phAniMap_FB.Write();
        
        phNAniMap_NS.Write();
        phNAniMap_EW.Write();
        phNAniMap_FB.Write();
        
        RatioAniMap_NS.Write();
        RatioAniMap_EW.Write();
        RatioAniMap_FB.Write();
        
        Isolate_NS.Write();
        Isolate_EW.Write();
        Isolate_FB.Write();
        
        Isolate_NN_NS.Write();
        Isolate_NN_EW.Write();
        Isolate_NN_FB.Write();
        
        AllSky_IsoMap.Write();
        
        AllSky_AniMap_NS.Write();
        AllSky_AniMap_EW.Write();
        AllSky_AniMap_FB.Write();
        
        NAllSky_IsoMap.Write();
        
        NAllSky_AniMap_NS.Write();
        NAllSky_AniMap_EW.Write();
        NAllSky_AniMap_FB.Write();
        
        AllSky_Isolate_NS.Write();
        AllSky_Isolate_EW.Write();
        AllSky_Isolate_FB.Write();
        
        MapsIsoRaio_NS.Write();
        MapsIsoRaio_EW.Write();
        MapsIsoRaio_FB.Write();
        
        */
        
        TemplateIso.Write();
        TemplateAniNS.Write();
        TemplateAniEW.Write();
        TemplateAniFB.Write();
        
        dataI.Write();
        dataNS.Write();
        dataEW.Write();
        dataFB.Write();
        
        ctemp.Write();
        
        results_file.Write();
        
         
        //////////////////////////////////// Closing files
        
        //results_file.Close();
        
        //output_log_file.close();
        
    }
    
    return 0;
  
}

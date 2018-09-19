//// C++ librariesA

#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMatrix.h"
#include "TRandom.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"
#include "TLegend.h"
#include "TSystem.h"
#include "TF1.h"
#include "TMinuit.h"
#include "TFitter.h"
#include "TObjArray.h"
#include "TProfile.h"
#include "TGraphErrors.h"

#include "TFractionFitter.h"

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "RooAddPdf.h"
#include "RooCBShape.h"
#include "RooBifurGauss.h"
#include "RooMinuit.h"
#include "RooArgList.h"
#include "RooMsgService.h"
#include "TRandom3.h"

using namespace RooFit;

TH1* __dat__;
TH1** __temp__;

int NDATA;

void ZeroRooFitVerbosity();
void RemoveZeroes(TH1* h);
void ResetRooFitVerbosity();
void ZeroRooFitVerbosity();
void fcnchisq(int& npar, double* deriv, double& f, double par[], int flag);
void fcnlike(int& npar, double* deriv, double& f, double par[], int flag);
double func(int npar, double par[], double _dochisqfit);

TH1* TemplateFitRF(TH1* dat, int ncomp, TH1* temp[], double _res[], double _reserr[], double initialguess[], bool fixtoinitialguess[], bool quiet, bool resettemplates, bool kClamping);
TH1* TemplateFitBH(TH1* dat, int ncomp, TH1* temp[], double _res[], double _reserr[], double initialguess[], bool quiet, bool resettemplates, bool kClamping, bool kChiSq);







void create_binning(Int_t n_bin_lat,Double_t lat_bin_min,Double_t lat_bin_max,Double_t* &binning,Bool_t cos_scale);


double func(int npar, double par[], double _dochisqfit){
  
  NDATA=0;

  double ChiSq = 0.0;
  double LogL = 0.0;
  
  for (int bb=0; bb<=__dat__->GetNbinsX()+1; bb++){
    
    double prob[npar];
    for (int icomp=0; icomp<npar; icomp++) {
      prob[icomp] = __temp__[icomp]->GetBinContent(bb);
    }
    
    double nevs = __dat__->GetBinContent(bb);
    double binw = __dat__->GetXaxis()->GetBinWidth(bb);

    // the bindiwth is required since we've to put the pdf
    // (a real pdf: sum of bincontents * binwidth =1) integrated over the bin
    // (or alternatively the bincontent of a histo-'pdf' normalized to have Integral(1))
    double nu[npar];
    double nu_i=0;
    for (int icomp=0; icomp<npar; icomp++) {
      nu[icomp] = par[icomp]*prob[icomp]*binw;
      nu_i += nu[icomp];
    }
    
    if (_dochisqfit) {
      double yexp = nu_i;
      double yobs = nevs;
      double sigmasq = yobs;//MLS
      //     double sigmasq = yexp;//LS (this would also imply to 'throw away' not-empty bins, if the model is 0...) 
      if (sigmasq>0) {
	ChiSq += pow((yobs-yexp), 2.0)/sigmasq;
	NDATA++;
      }
    }
    else {//LIKELIHOOD
      if (nu_i>0) {//this implies to 'throw away' not-empty bins, if the model is 0... 
	LogL += nevs*log(nu_i);
	NDATA++;
      }
    }
  }
  if (!_dochisqfit) { //LIKELIHOOD
    for (int icomp=0; icomp<npar; icomp++) {
      //needed for Extended ML
      LogL -= par[icomp]; 
    }
  }

  if (_dochisqfit) {
    return ChiSq;
  }
  else {//LIKELIHOOD
    return  -2.0*LogL;// '-' to find the MAXimim with MINuit, factor '2' so minuit gets the errors right (without changing ErrDef(1.0) to ErrDef(0.5)
  }
}


void RemoveZeroes(TH1* h) {

  if (!h) {
    printf("RemoveZeroes: the passed pointer is NULL!\n");
    return;
  }

  for (int ii=0; ii<=h->GetNbinsX()+1; ii++) {
    double bc = h->GetBinContent(ii);
    if (bc<1e-30) {
      h->SetBinContent(ii, 1e-30);
    }
  }
  
  return;
}




TH1* TemplateFitRF(TH1* dat, int ncomp, TH1* temp[], double _res[], double _reserr[], double initialguess[], bool fixtoinitialguess[], bool quiet, bool resettemplates, bool kClamping) {

  //  printf("<TemplateFitRF>\n");
  
  static int ncomp_old = ncomp;
  if (ncomp_old != ncomp) resettemplates=true;
  
  static TH1** tempnorm = NULL;
  bool templatesresetted = false;
  if (!tempnorm || resettemplates) {
    if (tempnorm) {
      for (int icomp=0; icomp<ncomp_old; icomp++) {
	if (tempnorm[icomp]) delete tempnorm[icomp];  
      }
      delete[] tempnorm;
    }
    templatesresetted = true;
    tempnorm = new TH1*[ncomp];
    for (int icomp=0; icomp<ncomp; icomp++) {
      tempnorm[icomp] = (TH1*)(temp[icomp]->Clone(Form("TemplateFit_tempnorm%d", icomp)));
      tempnorm[icomp]->Scale(1.0/temp[icomp]->Integral());
    }
  }
  
  // il caso NON EXTENDED non funziona...
#define EXTENDED
  
  if (quiet) ZeroRooFitVerbosity();
  
  //------------ RooFit ---------------
  
  static RooRealVar* x = NULL;
  static double datxmin = dat->GetXaxis()->GetXmin();
  static double datxmax = dat->GetXaxis()->GetXmax();
  double xresetted = false;
  if (!x || datxmin!=(dat->GetXaxis()->GetXmin()) || datxmax!=(dat->GetXaxis()->GetXmax())) { 
    if (x) delete x;
    xresetted = true;
    datxmin = dat->GetXaxis()->GetXmin();
    datxmax = dat->GetXaxis()->GetXmax();
    x = new RooRealVar("x", "x", datxmin, datxmax);
  }
  
  static RooRealVar** n = NULL;
  static double datgetentries = dat->GetEntries();
  bool nresetted = false;
  if (!n || templatesresetted || datgetentries!=dat->GetEntries()) {
    if (n) {
      for(int icomp=0; icomp<ncomp_old; icomp++) {
	delete n[icomp];
      }
      delete[] n;
    }
    nresetted = true;
    datgetentries = dat->GetEntries();
    n = new RooRealVar*[ncomp];
    for(int icomp=0; icomp<ncomp; icomp++) {
      double _initialguess = 0.5*datgetentries;
      if (initialguess) _initialguess = initialguess[icomp];
      double _lowlimit  = 0.0;
      double _highlimit = datgetentries;
      if (kClamping) {
	_lowlimit  = 0.0;
	_highlimit = datgetentries;
      }
      else {
	_lowlimit  = -7.0*sqrt(datgetentries);
	_highlimit = datgetentries+7.0*sqrt(datgetentries);
      }
      if (initialguess) {
	if (fixtoinitialguess) {
	  if (fixtoinitialguess[icomp]) {
	    _lowlimit  = initialguess[icomp];
	    _highlimit = initialguess[icomp];
	  }
	}
      }
      n[icomp] = new RooRealVar(Form("n_%d", icomp), Form("n_%d", icomp), _initialguess, _lowlimit, _highlimit);
    }
  }
  
  //--------- templates and model ------------
  
  if (resettemplates) {
    for(int icomp=0; icomp<ncomp; icomp++) {
      RemoveZeroes(tempnorm[icomp]);
    }
  }
  
  static TH1** hfromfile     = NULL;
  static RooDataHist** pdfdh = NULL;
  static RooHistPdf** hpdf   = NULL;
  bool hpdfresetted = false;
  if (!hfromfile || templatesresetted || xresetted) {
    if (hfromfile) {
      for (int icomp=0; icomp<ncomp_old; icomp++) {
	//	if (hfromfile[icomp]) delete hfromfile[icomp];//no!!! this is only a container for the tempnorm pointers!
	if (pdfdh[icomp]) delete pdfdh[icomp];
	if (hpdf[icomp]) delete hpdf[icomp];
      }
      delete[] hfromfile;
      if (pdfdh) delete[] pdfdh;
      if (hpdf) delete[] hpdf;
    }
    hpdfresetted = true;
    hfromfile = new TH1*[ncomp];
    pdfdh = new RooDataHist*[ncomp];
    hpdf = new RooHistPdf*[ncomp];
    for (int icomp=0; icomp<ncomp; icomp++) {
      hfromfile[icomp] = tempnorm[icomp];
      pdfdh[icomp] = new RooDataHist(Form("pdfdh_%d", icomp), Form("pdfdh_%d", icomp), RooArgList(*x), Import(*hfromfile[icomp]));
      hpdf[icomp] = new RooHistPdf(Form("hpdf_%d", icomp), Form("hpdf_%d", icomp), RooArgSet(*x), *pdfdh[icomp]);
    }
  }
  
  static RooAddPdf* model = NULL;
  bool modelresetted = false;
  if(!model || hpdfresetted || nresetted) {
    if (model) delete model;
    modelresetted = true;
    RooArgList ralh;
    RooArgList raln;
    for (int icomp=0; icomp<ncomp-1; icomp++) {
      ralh.add(*hpdf[icomp]);
      raln.add(*n[icomp]);
    }
    ralh.add(*hpdf[ncomp-1]);
#ifdef EXTENDED
    raln.add(*n[ncomp-1]);
#endif
    model = new RooAddPdf(Form("model"), Form("model"), ralh, raln);
  }

  //---------- dataset -----------------------
  RooDataHist sumhist(Form("sumhist"), Form("sumhist"), *x, Import(*dat));
  //  sumhist.Print();
  //  ovl = OVL(h[0], h[1]);

  //-----------Fit ---------------------------
#ifdef EXTENDED
  bool ext=kTRUE;
#else
  bool ext=kFALSE;
#endif
  model->fitTo(sumhist, Extended(ext), Save(kTRUE), SumW2Error(kFALSE), Verbose(!quiet), PrintLevel(quiet?-1:3), Warnings(!quiet));
  //  model->fitTo(sumhist, Extended(ext), Save(kTRUE), SumW2Error(kFALSE), Verbose(!quiet), PrintLevel(quiet?-1:3));

  double allbutlast = 0.0;
  double allbutlast_err = 0.0;
  for (int icomp=0; icomp<ncomp-1; icomp++) {
    _res[icomp] = n[icomp]->getVal();
    _reserr[icomp] = n[icomp]->getError();
    allbutlast += _res[icomp];
    allbutlast_err += pow(_reserr[icomp], 2.0);
    //FIXME: here a real propagation error (covmat) would be needed. In this way the error is much larger while should be much lower!
  }
  allbutlast_err = sqrt(allbutlast_err);
#ifdef EXTENDED
  _res[ncomp-1] = n[ncomp-1]->getVal();
  _reserr[ncomp-1] = n[ncomp-1]->getError();
#else
  _res[ncomp-1] = dat->Integral() - allbutlast;
  _reserr[ncomp-1] = allbutlast_err;
#endif
  
  TH1* result = (TH1*)(dat->Clone(Form("%s_TemplateFitRF_result", dat->GetName())));
  result->Reset();
  for (int icomp=0; icomp<ncomp; icomp++) {
    result->Add(temp[icomp], _res[icomp]/temp[icomp]->Integral());
  }
  
  if (quiet) ResetRooFitVerbosity();

  ncomp_old = ncomp;

  //  printf("</TemplateFitRF>\n");
  
  return result;
}

void ZeroRooFitVerbosity() {
  
  RooMsgService::instance().getStream(1).removeTopic(Minimization);
  RooMsgService::instance().getStream(1).removeTopic(Plotting);
  RooMsgService::instance().getStream(1).removeTopic(Fitting);
  RooMsgService::instance().getStream(1).removeTopic(Eval);
  RooMsgService::instance().getStream(1).removeTopic(Caching);
  RooMsgService::instance().getStream(1).removeTopic(ObjectHandling);
  RooMsgService::instance().getStream(1).removeTopic(InputArguments);
  RooMsgService::instance().getStream(1).removeTopic(DataHandling);
  RooMsgService::instance().getStream(1).removeTopic(NumIntegration);
  
  RooMsgService::instance().getStream(0).removeTopic(Minimization);
  RooMsgService::instance().getStream(0).removeTopic(Plotting);
  RooMsgService::instance().getStream(0).removeTopic(Fitting);
  RooMsgService::instance().getStream(0).removeTopic(Eval);
  RooMsgService::instance().getStream(0).removeTopic(Caching);
  RooMsgService::instance().getStream(0).removeTopic(ObjectHandling);
  RooMsgService::instance().getStream(0).removeTopic(InputArguments);
  RooMsgService::instance().getStream(0).removeTopic(DataHandling);
  RooMsgService::instance().getStream(0).removeTopic(Generation);
  RooMsgService::instance().getStream(0).removeTopic(Integration);
  RooMsgService::instance().getStream(0).removeTopic(LinkStateMgmt);
  RooMsgService::instance().getStream(0).removeTopic(Optimization);
  RooMsgService::instance().getStream(0).removeTopic(Tracing);
  RooMsgService::instance().getStream(0).removeTopic(Contents);
  RooMsgService::instance().getStream(0).removeTopic(NumIntegration);
  
  // cout << endl ;
  // RooMsgService::instance().Print() ;
  // cout << endl ;
  // sleep(3);

  // RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  // RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

  RooMsgService::instance().setSilentMode(true); 
  
  return;
}

void ResetRooFitVerbosity() {

  RooMsgService::instance().reset(); 

  return;
}


void create_binning(Int_t n_bin_lat,Double_t lat_bin_min,Double_t lat_bin_max,Double_t* &binning,Bool_t cos_scale) {
  binning = new Double_t[n_bin_lat+1];

  if (cos_scale==false) {
    binning[0]=lat_bin_min;
    for(Int_t idx_b=1; idx_b<=n_bin_lat; idx_b++)
      binning[idx_b] = binning[idx_b-1] + (lat_bin_max-lat_bin_min)/n_bin_lat;
  }
  else {
    Double_t binning_shift = (TMath::Cos(TMath::DegToRad()*(lat_bin_max-90.0))-TMath::Cos(TMath::DegToRad()*(lat_bin_min-90.0)))/n_bin_lat;  //I subtrack 90 degrees to refer the costheta to the y axis !!
    Double_t prev_bin;  //Define a new variable for the "previous" bin value
    binning[0]=lat_bin_min;
    for(Int_t idx_b=1; idx_b<=n_bin_lat; idx_b++) {
      prev_bin = TMath::DegToRad()*(binning[idx_b-1]-90.0);
      binning[idx_b] = 90.0-(TMath::RadToDeg()*TMath::ACos(TMath::Cos(prev_bin)+binning_shift));
    }
  }
}


TH1* TemplateFitBH(TH1* dat, int ncomp, TH1* temp[], double _res[], double _reserr[], double initialguess[], bool quiet, bool resettemplates, bool kClamping, bool kChiSq) {
  
  //  printf("<TemplateFitBH>\n");

  //--------- templates and data ------------
  
  static int ncomp_old = ncomp;
  if (ncomp_old != ncomp) resettemplates=true;
  
  static TH1** tempnorm = NULL;
  bool templatesresetted = false;
  if (!tempnorm || resettemplates) {
    if (tempnorm) {
      for (int icomp=0; icomp<ncomp_old; icomp++) {
	if (tempnorm[icomp]) delete tempnorm[icomp];  
      }
      delete[] tempnorm;
    }
    templatesresetted = true;
    tempnorm = new TH1*[ncomp];
    for (int icomp=0; icomp<ncomp; icomp++) {
      tempnorm[icomp] = (TH1*)(temp[icomp]->Clone(Form("TemplateFit_tempnorm%d", icomp)));
      tempnorm[icomp]->Scale(1.0/temp[icomp]->Integral("width"));//"width" makes a real pdf
      //      printf("TemplateFit_tempnorm%d) Integral=%f, Integral(\"width\")=%f\n", icomp, tempnorm[icomp]->Integral(), tempnorm[icomp]->Integral("width"));
    }
  }
  
  if (resettemplates) {
    for(int icomp=0; icomp<ncomp; icomp++) {
      RemoveZeroes(tempnorm[icomp]);
    }
  }

  //----- set the global variables ----
  __dat__ = dat;
  //  printf("%s) Integral=%f, Integral(\"width\")=%f\n", __dat__->GetName(), __dat__->Integral(), __dat__->Integral("width"));
  __temp__ = new TH1*[ncomp];
  for(int icomp=0; icomp<ncomp; icomp++) {
    __temp__[icomp] = tempnorm[icomp];
}

  //--------- minimization ------------
  
  int npar=ncomp;
  
  TMinuit minuit(npar);
  minuit.SetPrintLevel(2*(int)(!quiet)-1);

  if (kChiSq) {
    minuit.SetFCN(fcnchisq);
  }
  else {
    minuit.SetFCN(fcnlike);
  }
  
  double par[npar];
  double stepSize[npar];
  double minVal[npar];
  double maxVal[npar];
  string parName[npar];
  double arglist[10];

  for (int icomp=0; icomp<npar; icomp++) {
    par[icomp] = 0;
    if (initialguess) par[icomp] = initialguess[icomp];
    stepSize[icomp] = 0.0001;
    //minVal[icomp] = 0.0;
    minVal[icomp] = -1.;         // Modified parameter low limit ! Before set at 0 !
    if (kClamping) {
      maxVal[icomp] = dat->Integral();
    }
    else {
      maxVal[icomp] = dat->Integral()+7.0*sqrt(dat->Integral());
    }
    parName[icomp] = Form("N_%d", icomp);
  }
  
  for (int i=0; i<npar; i++){
    minuit.DefineParameter(i, parName[i].c_str(), 
			   par[i], stepSize[i], minVal[i], maxVal[i]);
  }
  
  minuit.SetErrorDef(1.0);//1.0 for chisq or -2*lnL, 0.5 for -lnL
  
  arglist[0]=2;
  int ierflag=0;
  minuit.mnexcm("SET STR",arglist,1,ierflag);
  minuit.Migrad();

  if (!kClamping) {
    //the problem is that if all the parameters go <0 the fcn starts to return very negatives values
    // and so the minimization continues findind the minimum (-inf) when the parameters are -inf
    //the only way I found is to perform a "clamped" (on the bottom) fit
    // and then replicate with the un-clamped one
    // I'M NOT SURE IF THE SUBSEQUENT CALL TO MIGRAD REALLY DOES ANYTHING...
    for (int icomp=0; icomp<npar; icomp++) {
      minVal[icomp] = 0.0-7.0*sqrt(dat->Integral());
      maxVal[icomp] = dat->Integral()+7.0*sqrt(dat->Integral());
      // minVal[icomp] = 0.0;
      // maxVal[icomp] = 0.0;
    }
    minuit.Migrad();
  }
  
  //--------- get the result ------------  
  
  if (!quiet) printf("\n");
  
  for (int icomp=0; icomp<npar; icomp++){
    minuit.GetParameter(icomp, _res[icomp], _reserr[icomp]);
    if (!quiet) {
      cout << " Parameter " << parName[icomp].c_str() << " value " << _res[icomp] << " error " << _reserr[icomp]<<endl;
    }
  }

  double chisq = func(npar, _res, true);
  int freepars = minuit.GetNumFreePars();
  int ndof = NDATA-freepars;
  if (!quiet) printf("chisq: %f/%d = %f\n", chisq, ndof, chisq/ndof);
  
  TH1* result = (TH1*)(dat->Clone(Form("%s_TemplateFitBH_result", dat->GetName())));
  result->Reset();
  for (int icomp=0; icomp<ncomp; icomp++) {
    result->Add(temp[icomp], _res[icomp]/temp[icomp]->Integral());
  }
  
  ncomp_old = ncomp;

  //  printf("</TemplateFitBH>\n");
  
  return result;
}

void fcnchisq(int& npar, double* deriv, double& f, double par[], int flag){
  flag=flag;
  deriv=deriv;
  f=func(npar, par, true);
  return;
}

void fcnlike(int& npar, double* deriv, double& f, double par[], int flag){
  flag=flag;
  deriv=deriv;
  f=func(npar, par, false);
  return;
}


void CreateMySample(const char* data_path = "./ROOTFiles/010_1528984302_maps_result.root", const char* sample_out_path = "./ROOTFiles/SampleTemplateFit.root") {

  //////////////////////////////////////// Variables and Hitsos //////////////////////////////////////
  
  TH2D AniTemplateNS,AniTemplateEW,AniTemplateFB,IsoTemplate;

  Int_t n_bin_lon = 360;                                                                      // -> Number of bins along longitude axis                                                                                                                      
  Double_t lon_bin_min = -180.0;                                                              // -> Set max and min for longitude binning                                                                                                                    
  Double_t lon_bin_max = 180.0;

  Int_t n_bin_lat = 180;                                                                      // -> Number of bins along latitude axis
  Double_t lat_bin_min = -90.0;                                                               // -> Set max and min for latitude binning                                                                                                                     
  Double_t lat_bin_max = 90.0;

  Double_t* binning = nullptr;

  TRandom3 r_gen(90);
  Double_t l,b,cosb,w_NS;
  Int_t data_events = 1e+5;

  create_binning(n_bin_lat,lat_bin_min,lat_bin_max,binning,true);
  TH2D DataHisto("data_histo","Anisotropic All Sky Map (NS); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon,lon_bin_min,lon_bin_max,n_bin_lat,binning);
  
  //////////////////////////////////////// Filling randomly data histo with NS dipole ///////////////////////                                                                                                                                                 

  for(Int_t i=0; i<data_events; i++) {
    l = r_gen.Uniform(-180.,180.);
    cosb = r_gen.Uniform(-1., 1.);
    b = TMath::RadToDeg()*TMath::ACos(cosb) - 90.;

    w_NS = 1 + TMath::Sqrt(3./(4*TMath::Pi()))*TMath::Cos((b+90)*TMath::DegToRad());

    DataHisto.Fill(l,b,w_NS);

  }


  //////////////////////////////////////// Reading templates from file and saving histo /////////////////////                                                                                                                                                 

  TFile infile(data_path,"READ");
  if(infile.IsZombie()) {
    std::cout<<"\n\nError opening input file\n\n";
    exit(-2);
  }
  new (&IsoTemplate) TH2D(*(TH2D*)infile.Get("AllSky_IsoMap"));
  new (&AniTemplateNS) TH2D(*(TH2D*)infile.Get("AllSky_AniMap_NS"));
  new (&AniTemplateEW) TH2D(*(TH2D*)infile.Get("AllSky_AniMap_EW"));
  new (&AniTemplateFB) TH2D(*(TH2D*)infile.Get("AllSky_AniMap_FB"));
  infile.Close();

  IsoTemplate.SetName("IsoTemplate");
  IsoTemplate.SetTitle("Isotropic Map Template");

  AniTemplateNS.SetName("AniTemplateNS");
  AniTemplateNS.SetTitle("Anisotropic (NS) Map Template");

  AniTemplateEW.SetName("AniTemplateEW");
  AniTemplateEW.SetTitle("Anisotropic (EW) Map Template");

  AniTemplateFB.SetName("AniTemplateFB");
  AniTemplateFB.SetTitle("Anisotropic (FB) Map Template");

  //////////////////////////////////////////// Writing sampe !
  
  TFile outfile(sample_out_path,"RECREATE");

  
  IsoTemplate.Write();
  AniTemplateNS.Write();
  AniTemplateEW.Write();
  AniTemplateFB.Write();

  DataHisto.Write();

  outfile.Close();


}

void MyTemplateFit(const char* TemplateAssets = "./ROOTFiles/SampleTemplateFit.root") {


  //////////////////////////////////////// Variables and Hitsos //////////////////////////////////////
  
  TH1* Templates[4] = {nullptr,nullptr,nullptr,nullptr};
  TH2D* DataHisto = nullptr;

  Double_t res[4],res_err[4];
  Double_t initialValues[4]={1,0,0,0};
  bool fix[4]={false,false,false,false};
  
  //////////////////////////////////////// Reading SampleTemplate File ///////// /////////////////////
  
  TFile infile(TemplateAssets,"READ");
  if(infile.IsZombie()) {
    std::cout<<"\n\nError opening input file\n\n";
    exit(-2);
  }
  
  Templates[0] = (TH2D*)infile.Get("IsoTemplate")->Clone("IsoTemplate");
  Templates[1] = (TH2D*)infile.Get("AniTemplateNS")->Clone("AniTemplateNS");
  Templates[2] = (TH2D*)infile.Get("AniTemplateEW")->Clone("AniTemplateEW");
  Templates[3] = (TH2D*)infile.Get("AniTemplateFB")->Clone("AniTemplateFB");
  DataHisto =  (TH2D*)infile.Get("data_histo")->Clone("data_histo");

  Templates[0]->SetDirectory(0);
  Templates[1]->SetDirectory(0);
  Templates[2]->SetDirectory(0);
  Templates[3]->SetDirectory(0);
  DataHisto->SetDirectory(0);
  
  infile.Close();

  ////////////////////////////////////////// Finally Fit !!! ;-(
  
  //TemplateFitRF(DataHisto,4,Templates,res,res_err,initialValues,fix,false,false,false);
  TemplateFitBH(DataHisto,4,Templates,res,res_err,initialValues,false,false,false,true); // ---> !!!! Chi2 fit !! Substitute last "true" with "false" for likelihood fit !
  
}


#include "MyHead.h"

///////////////////////// RooFit

using namespace RooFit;

static TH1* __dat__;
static TH1** __temp__;
static int NDATA;
static int chi2_filter_bc = 5;
static bool relative_fit = false;

//////////////////////////////////


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
            //double sigmasq = yobs;//MLS
            double sigmasq = TMath::Power(__dat__->GetBinError(bb),2);
            //     double sigmasq = yexp;//LS (this would also imply to 'throw away' not-empty bins, if the model is 0...)
            if(relative_fit)
            {
                if(sigmasq!=0)
                {
                    ChiSq += pow((yobs-yexp), 2.0)/sigmasq;
                    NDATA++;
                }
            }
            else
            {
                if(yobs >= chi2_filter_bc && sigmasq!=0) {
                    ChiSq += pow((yobs-yexp), 2.0)/sigmasq;
                    NDATA++;
                }
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
                //    if (hfromfile[icomp]) delete hfromfile[icomp];//no!!! this is only a container for the tempnorm pointers!
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
    //bool modelresetted = false;
    if(!model || hpdfresetted || nresetted) {
        if (model) delete model;
        //modelresetted = true;
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

TH1* TemplateFitBH(
                    TH1* dat,
                    int ncomp,
                    TH1* temp[],
                    double _res[],
                    double _reserr[],
                    double initialguess[],
                    bool quiet,
                    bool resettemplates,
                    bool kClamping,
                    bool kChiSq,
                    std::ofstream &log_file,
                    fitResult &tmp_fit,
                    Int_t fit_element,
                    bool HS
                   )

{
    
    //  printf("<TemplateFitBH>\n");
    
    //--------- templates and data ------------
    
    double ani_level = 0;
    
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
            //tempnorm[icomp]->Scale(1.0/temp[icomp]->Integral("width"));//"width" makes a real pdf
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
    std::string parName[npar];
    double arglist[10];
    
    for (int icomp=0; icomp<npar; icomp++) {
        par[icomp] = 0;
        if (initialguess) par[icomp] = initialguess[icomp];
        stepSize[icomp] = 0.0001;
        //minVal[icomp] = 0.0;
        //minVal[icomp] = -1.;         // Modified parameter low limit ! Before set at 0 !
        minVal[icomp] = -(dat->Integral()+7.0*sqrt(dat->Integral()));
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
    
    
    for (int icomp=0; icomp<npar; icomp++)
    {
        minuit.GetParameter(icomp, _res[icomp], _reserr[icomp]);
        if (!quiet)
        {
            std::cout << "Parameter " << parName[icomp].c_str() << " value " << _res[icomp] << " error " << _reserr[icomp]<<std::endl;
            log_file << "Parameter " << parName[icomp].c_str() << " value " << _res[icomp] << " error " << _reserr[icomp]<<std::endl;
        }
    }

    /*
    for (int icomp=0; icomp<npar; icomp++)
    {
        minuit.GetParameter(icomp, _res[icomp], _reserr[icomp]);
        if (!quiet)
        {
            if(icomp==0)
            {
                std::cout << "Parameter " << parName[icomp].c_str() << " value " << _res[icomp] << " error " << _reserr[icomp]<<std::endl;
                log_file << "Parameter " << parName[icomp].c_str() << " value " << _res[icomp] << " error " << _reserr[icomp]<<std::endl;
            }
            else if(icomp==1)
            {
                std::cout << "Parameter " << parName[icomp].c_str() << " value " << _res[icomp]*TMath::Sqrt((3/(4*TMath::Pi()))) << " error " << _reserr[icomp]*TMath::Sqrt((3/(4*TMath::Pi()))) <<std::endl;
                log_file << "Parameter " << parName[icomp].c_str() << " value " << _res[icomp]*TMath::Sqrt((3/(4*TMath::Pi()))) << " error " << _reserr[icomp]*TMath::Sqrt((3/(4*TMath::Pi()))) <<std::endl;
                //std::cout << "Parameter " << parName[icomp].c_str() << " value " << _res[icomp]*TMath::Sqrt((3/(4*TMath::Pi()))) << " error " << _reserr[icomp] << std::endl;
                //log_file << "Parameter " << parName[icomp].c_str() << " value " << _res[icomp]*TMath::Sqrt((3/(4*TMath::Pi()))) << " error " << _reserr[icomp] << std::endl;
            }
            else
            {
                std::cout << "Parameter " << parName[icomp].c_str() << " value " << _res[icomp]*TMath::Sqrt((6/(4*TMath::Pi()))) << " error " << _reserr[icomp]*TMath::Sqrt((6/(4*TMath::Pi()))) <<std::endl;
                log_file << "Parameter " << parName[icomp].c_str() << " value " << _res[icomp]*TMath::Sqrt((6/(4*TMath::Pi()))) << " error " << _reserr[icomp]*TMath::Sqrt((6/(4*TMath::Pi()))) <<std::endl;
                //std::cout << "Parameter " << parName[icomp].c_str() << " value " << _res[icomp]*TMath::Sqrt((6/(4*TMath::Pi()))) << " error " << _reserr[icomp] << std::endl;
                //log_file << "Parameter " << parName[icomp].c_str() << " value " << _res[icomp]*TMath::Sqrt((6/(4*TMath::Pi()))) << " error " << _reserr[icomp]  << std::endl;
            }
        }
    }
    */
    
    double chisq = func(npar, _res, true);
    int freepars = minuit.GetNumFreePars();
    int ndof = NDATA-freepars;
    if (!quiet)
    {
        printf("chisq: %f/%d = %f\n", chisq, ndof, chisq/ndof);
        log_file << "\nchisq: " << chisq << "ndof: " << ndof << "chisq/ndof: " << chisq/ndof;
    }
    
    
    TH1* result = (TH1*)(dat->Clone(Form("%s_TemplateFitBH_result", dat->GetName())));
    result->Reset();
    for (int icomp=0; icomp<ncomp; icomp++) {
        result->Add(temp[icomp], _res[icomp]/temp[icomp]->Integral());
    }
    
    ncomp_old = ncomp;
    
    //  printf("</TemplateFitBH>\n");
    
    
    ///////////////// Filling class
    
    if(!HS)
    {
            
        double covMatrix_LS[4][4];
        minuit.mnemat(&covMatrix_LS[0][0],4);
        ani_level = compute_ani_level(covMatrix_LS,_res,_reserr,log_file);
            
        for(Int_t par_idx = 0; par_idx < 4; ++par_idx)
        {
            tmp_fit.fit_par_LS[fit_element][par_idx] = _res[par_idx];
            tmp_fit.fit_err_LS[fit_element][par_idx] = _reserr[par_idx];
        }
            
        tmp_fit.chi2_LS[fit_element] = chisq;
        tmp_fit.ndf_LS[fit_element] = ndof;
        tmp_fit.chi2_r_LS[fit_element] = chisq/ndof;
        tmp_fit.delta_LS[fit_element] = ani_level;
        tmp_fit.sum_par_LS[fit_element] = _res[0] + _res[1] + _res[2] + _res[3];
            
        switch(fit_element)
        {
            case 0:
                for(Int_t Ridx = 0; Ridx < 4; ++Ridx)
                    for(Int_t Cidx = 0; Cidx < 4; ++Cidx)
                        tmp_fit.CMatrix_Iso_LS[Ridx][Cidx] = covMatrix_LS[Ridx][Cidx];
                break;
                
            case 1:
                for(Int_t Ridx = 0; Ridx < 4; ++Ridx)
                    for(Int_t Cidx = 0; Cidx < 4; ++Cidx)
                        tmp_fit.CMatrix_NS_LS[Ridx][Cidx] = covMatrix_LS[Ridx][Cidx];
                break;
                    
            case 2:
                for(Int_t Ridx = 0; Ridx < 4; ++Ridx)
                    for(Int_t Cidx = 0; Cidx < 4; ++Cidx)
                        tmp_fit.CMatrix_EW_LS[Ridx][Cidx] = covMatrix_LS[Ridx][Cidx];
                break;
                
            case 3:
                for(Int_t Ridx = 0; Ridx < 4; ++Ridx)
                    for(Int_t Cidx = 0; Cidx < 4; ++Cidx)
                        tmp_fit.CMatrix_FB_LS[Ridx][Cidx] = covMatrix_LS[Ridx][Cidx];
                break;
                    
            case 4:
                for(Int_t Ridx = 0; Ridx < 4; ++Ridx)
                    for(Int_t Cidx = 0; Cidx < 4; ++Cidx)
                        tmp_fit.CMatrix_NS_EW_LS[Ridx][Cidx] = covMatrix_LS[Ridx][Cidx];
                break;
                    
            case 5:
                for(Int_t Ridx = 0; Ridx < 4; ++Ridx)
                    for(Int_t Cidx = 0; Cidx < 4; ++Cidx)
                        tmp_fit.CMatrix_NS_FB_LS[Ridx][Cidx] = covMatrix_LS[Ridx][Cidx];
                break;
                    
            case 6:
                for(Int_t Ridx = 0; Ridx < 4; ++Ridx)
                    for(Int_t Cidx = 0; Cidx < 4; ++Cidx)
                        tmp_fit.CMatrix_EW_FB_LS[Ridx][Cidx] = covMatrix_LS[Ridx][Cidx];
                break;
                    
            case 7:
                for(Int_t Ridx = 0; Ridx < 4; ++Ridx)
                    for(Int_t Cidx = 0; Cidx < 4; ++Cidx)
                        tmp_fit.CMatrix_Full_LS[Ridx][Cidx] = covMatrix_LS[Ridx][Cidx];
                break;
                    
        }
    }
    else
    {
            
        Double_t covMatrix_HS[4][4];
        minuit.mnemat(&covMatrix_HS[0][0],4);
        ani_level = compute_ani_level(covMatrix_HS,_res,_reserr,log_file);
            
        for(Int_t par_idx = 0; par_idx < 4; ++par_idx)
        {
            tmp_fit.fit_par_HS[fit_element][par_idx] = _res[par_idx];
            tmp_fit.fit_err_HS[fit_element][par_idx] = _reserr[par_idx];
        }
            
        tmp_fit.chi2_HS[fit_element] = chisq;
        tmp_fit.ndf_HS[fit_element] = ndof;
        tmp_fit.chi2_r_HS[fit_element] = chisq/ndof;
        tmp_fit.delta_HS[fit_element] = ani_level;
        tmp_fit.sum_par_HS[fit_element] = _res[0] + _res[1] + _res[2] + _res[3];
         
        switch(fit_element)
        {
            case 0:
                for(Int_t Ridx = 0; Ridx < 4; ++Ridx)
                    for(Int_t Cidx = 0; Cidx < 4; ++Cidx)
                        tmp_fit.CMatrix_Iso_LS[Ridx][Cidx] = covMatrix_HS[Ridx][Cidx];
                break;
                    
            case 1:
                for(Int_t Ridx = 0; Ridx < 4; ++Ridx)
                    for(Int_t Cidx = 0; Cidx < 4; ++Cidx)
                        tmp_fit.CMatrix_NS_LS[Ridx][Cidx] = covMatrix_HS[Ridx][Cidx];
                break;
                    
            case 2:
                for(Int_t Ridx = 0; Ridx < 4; ++Ridx)
                    for(Int_t Cidx = 0; Cidx < 4; ++Cidx)
                        tmp_fit.CMatrix_EW_LS[Ridx][Cidx] = covMatrix_HS[Ridx][Cidx];
                break;
                    
            case 3:
                for(Int_t Ridx = 0; Ridx < 4; ++Ridx)
                    for(Int_t Cidx = 0; Cidx < 4; ++Cidx)
                        tmp_fit.CMatrix_FB_LS[Ridx][Cidx] = covMatrix_HS[Ridx][Cidx];
                break;
                    
            case 4:
                for(Int_t Ridx = 0; Ridx < 4; ++Ridx)
                    for(Int_t Cidx = 0; Cidx < 4; ++Cidx)
                        tmp_fit.CMatrix_NS_EW_LS[Ridx][Cidx] = covMatrix_HS[Ridx][Cidx];
                break;
                    
            case 5:
                for(Int_t Ridx = 0; Ridx < 4; ++Ridx)
                    for(Int_t Cidx = 0; Cidx < 4; ++Cidx)
                        tmp_fit.CMatrix_NS_FB_LS[Ridx][Cidx] = covMatrix_HS[Ridx][Cidx];
                break;
                    
            case 6:
                for(Int_t Ridx = 0; Ridx < 4; ++Ridx)
                    for(Int_t Cidx = 0; Cidx < 4; ++Cidx)
                        tmp_fit.CMatrix_EW_FB_LS[Ridx][Cidx] = covMatrix_HS[Ridx][Cidx];
                break;
                    
            case 7:
                for(Int_t Ridx = 0; Ridx < 4; ++Ridx)
                    for(Int_t Cidx = 0; Cidx < 4; ++Cidx)
                        tmp_fit.CMatrix_Full_LS[Ridx][Cidx] = covMatrix_HS[Ridx][Cidx];
                break;
                    
        }
    }
   
    //////////////////////////////////////
    
    return result;
}



TH1* TemplateFitBH_rel(
                        TH1* dat,
                        int ncomp,
                        TH1* temp[],
                        double _res[],
                        double _reserr[],
                        double initialguess[],
                        bool quiet,
                        bool resettemplates,
                        bool kClamping,
                        bool kChiSq,
                        std::ofstream &log_file,
                        relative_fitResult &tmp_fit,
                        Int_t fit_element,
                        bool HS
                    )

{
    
    //  printf("<TemplateFitBH>\n");
    
    //--------- templates and data ------------
    
    double ani_level = 0;
    relative_fit = true;
    
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
            //tempnorm[icomp]->Scale(1.0/temp[icomp]->Integral("width"));//"width" makes a real pdf
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
    std::string parName[npar];
    double arglist[10];
    
    for (int icomp=0; icomp<npar; icomp++) {
        par[icomp] = 0;
        if (initialguess) par[icomp] = initialguess[icomp];
        stepSize[icomp] = 0.0001;
        minVal[icomp] = -100;
        maxVal[icomp] = 100;
        /*
        //minVal[icomp] = 0.0;
        //minVal[icomp] = -1.;         // Modified parameter low limit ! Before set at 0 !
        minVal[icomp] = -(dat->Integral()+7.0*sqrt(dat->Integral()));
        if (kClamping) {
            maxVal[icomp] = dat->Integral();
        }
        else {
            maxVal[icomp] = dat->Integral()+7.0*sqrt(dat->Integral());
        }
         */
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
    
    
    for (int icomp=0; icomp<npar; icomp++)
    {
        minuit.GetParameter(icomp, _res[icomp], _reserr[icomp]);
        if (!quiet)
        {
            std::cout << "Parameter " << parName[icomp].c_str() << " value " << _res[icomp] << " error " << _reserr[icomp]<<std::endl;
            log_file << "Parameter " << parName[icomp].c_str() << " value " << _res[icomp] << " error " << _reserr[icomp]<<std::endl;
        }
    }
    
    /*
     for (int icomp=0; icomp<npar; icomp++)
     {
     minuit.GetParameter(icomp, _res[icomp], _reserr[icomp]);
     if (!quiet)
     {
     if(icomp==0)
     {
     std::cout << "Parameter " << parName[icomp].c_str() << " value " << _res[icomp] << " error " << _reserr[icomp]<<std::endl;
     log_file << "Parameter " << parName[icomp].c_str() << " value " << _res[icomp] << " error " << _reserr[icomp]<<std::endl;
     }
     else if(icomp==1)
     {
     std::cout << "Parameter " << parName[icomp].c_str() << " value " << _res[icomp]*TMath::Sqrt((3/(4*TMath::Pi()))) << " error " << _reserr[icomp]*TMath::Sqrt((3/(4*TMath::Pi()))) <<std::endl;
     log_file << "Parameter " << parName[icomp].c_str() << " value " << _res[icomp]*TMath::Sqrt((3/(4*TMath::Pi()))) << " error " << _reserr[icomp]*TMath::Sqrt((3/(4*TMath::Pi()))) <<std::endl;
     //std::cout << "Parameter " << parName[icomp].c_str() << " value " << _res[icomp]*TMath::Sqrt((3/(4*TMath::Pi()))) << " error " << _reserr[icomp] << std::endl;
     //log_file << "Parameter " << parName[icomp].c_str() << " value " << _res[icomp]*TMath::Sqrt((3/(4*TMath::Pi()))) << " error " << _reserr[icomp] << std::endl;
     }
     else
     {
     std::cout << "Parameter " << parName[icomp].c_str() << " value " << _res[icomp]*TMath::Sqrt((6/(4*TMath::Pi()))) << " error " << _reserr[icomp]*TMath::Sqrt((6/(4*TMath::Pi()))) <<std::endl;
     log_file << "Parameter " << parName[icomp].c_str() << " value " << _res[icomp]*TMath::Sqrt((6/(4*TMath::Pi()))) << " error " << _reserr[icomp]*TMath::Sqrt((6/(4*TMath::Pi()))) <<std::endl;
     //std::cout << "Parameter " << parName[icomp].c_str() << " value " << _res[icomp]*TMath::Sqrt((6/(4*TMath::Pi()))) << " error " << _reserr[icomp] << std::endl;
     //log_file << "Parameter " << parName[icomp].c_str() << " value " << _res[icomp]*TMath::Sqrt((6/(4*TMath::Pi()))) << " error " << _reserr[icomp]  << std::endl;
     }
     }
     }
     */
    
    double chisq = func(npar, _res, true);
    int freepars = minuit.GetNumFreePars();
    int ndof = NDATA-freepars;
    if (!quiet)
    {
        printf("chisq: %f/%d = %f\n", chisq, ndof, chisq/ndof);
        log_file << "\nchisq: " << chisq << "ndof: " << ndof << "chisq/ndof: " << chisq/ndof;
    }
    
    
    TH1* result = (TH1*)(dat->Clone(Form("%s_TemplateFitBH_result", dat->GetName())));
    result->Reset();
    for (int icomp=0; icomp<ncomp; icomp++) {
        result->Add(temp[icomp], _res[icomp]/temp[icomp]->Integral());
    }
    
    ncomp_old = ncomp;
    
    //  printf("</TemplateFitBH>\n");
    
    
    ///////////////// Filling class
    
    if(!HS)
    {
            
        double covMatrix_LS[3][3];
        minuit.mnemat(&covMatrix_LS[0][0],3);
        ani_level = compute_relative_ani_level(covMatrix_LS,_res,_reserr,log_file);
            
        for(Int_t par_idx = 0; par_idx < 3; ++par_idx)
        {
            tmp_fit.fit_par_LS[fit_element][par_idx] = _res[par_idx];
            tmp_fit.fit_err_LS[fit_element][par_idx] = _reserr[par_idx];
        }
            
        tmp_fit.chi2_LS[fit_element] = chisq;
        tmp_fit.ndf_LS[fit_element] = ndof;
        tmp_fit.chi2_r_LS[fit_element] = chisq/ndof;
        tmp_fit.delta_LS[fit_element] = ani_level;
        tmp_fit.sum_par_LS[fit_element] = _res[0] + _res[1] + _res[2];
            
        switch(fit_element)
        {
                    
            case 0:
                for(Int_t Ridx = 0; Ridx < 3; ++Ridx)
                    for(Int_t Cidx = 0; Cidx < 3; ++Cidx)
                        tmp_fit.CMatrix_NS_LS[Ridx][Cidx] = covMatrix_LS[Ridx][Cidx];
                break;
                    
            case 1:
                for(Int_t Ridx = 0; Ridx < 3; ++Ridx)
                    for(Int_t Cidx = 0; Cidx < 3; ++Cidx)
                        tmp_fit.CMatrix_EW_LS[Ridx][Cidx] = covMatrix_LS[Ridx][Cidx];
                break;
                    
            case 2:
                for(Int_t Ridx = 0; Ridx < 3; ++Ridx)
                    for(Int_t Cidx = 0; Cidx < 3; ++Cidx)
                        tmp_fit.CMatrix_FB_LS[Ridx][Cidx] = covMatrix_LS[Ridx][Cidx];
                break;
                    
            case 3:
                for(Int_t Ridx = 0; Ridx < 3; ++Ridx)
                    for(Int_t Cidx = 0; Cidx < 3; ++Cidx)
                        tmp_fit.CMatrix_NS_EW_LS[Ridx][Cidx] = covMatrix_LS[Ridx][Cidx];
                break;
                    
            case 4:
                for(Int_t Ridx = 0; Ridx < 3; ++Ridx)
                    for(Int_t Cidx = 0; Cidx < 3; ++Cidx)
                        tmp_fit.CMatrix_NS_FB_LS[Ridx][Cidx] = covMatrix_LS[Ridx][Cidx];
                break;
                    
            case 5:
                for(Int_t Ridx = 0; Ridx < 3; ++Ridx)
                    for(Int_t Cidx = 0; Cidx < 3; ++Cidx)
                        tmp_fit.CMatrix_EW_FB_LS[Ridx][Cidx] = covMatrix_LS[Ridx][Cidx];
                break;
                    
            case 6:
                for(Int_t Ridx = 0; Ridx < 3; ++Ridx)
                    for(Int_t Cidx = 0; Cidx < 3; ++Cidx)
                        tmp_fit.CMatrix_Full_LS[Ridx][Cidx] = covMatrix_LS[Ridx][Cidx];
                break;
    
        }
    }
    else
    {
            
        Double_t covMatrix_HS[3][3];
        minuit.mnemat(&covMatrix_HS[0][0],3);
        ani_level = compute_relative_ani_level(covMatrix_HS,_res,_reserr,log_file);
            
        for(Int_t par_idx = 0; par_idx < 3; ++par_idx)
        {
            tmp_fit.fit_par_HS[fit_element][par_idx] = _res[par_idx];
            tmp_fit.fit_err_HS[fit_element][par_idx] = _reserr[par_idx];
        }
            
        tmp_fit.chi2_HS[fit_element] = chisq;
        tmp_fit.ndf_HS[fit_element] = ndof;
        tmp_fit.chi2_r_HS[fit_element] = chisq/ndof;
        tmp_fit.delta_HS[fit_element] = ani_level;
        tmp_fit.sum_par_HS[fit_element] = _res[0] + _res[1] + _res[2];
        
        switch(fit_element)
        {
                    
            case 0:
                for(Int_t Ridx = 0; Ridx < 3; ++Ridx)
                    for(Int_t Cidx = 0; Cidx < 3; ++Cidx)
                        tmp_fit.CMatrix_NS_LS[Ridx][Cidx] = covMatrix_HS[Ridx][Cidx];
                break;
                    
            case 1:
                for(Int_t Ridx = 0; Ridx < 3; ++Ridx)
                    for(Int_t Cidx = 0; Cidx < 3; ++Cidx)
                        tmp_fit.CMatrix_EW_LS[Ridx][Cidx] = covMatrix_HS[Ridx][Cidx];
                break;
                    
            case 2:
                for(Int_t Ridx = 0; Ridx < 3; ++Ridx)
                    for(Int_t Cidx = 0; Cidx < 3; ++Cidx)
                        tmp_fit.CMatrix_FB_LS[Ridx][Cidx] = covMatrix_HS[Ridx][Cidx];
                break;
                    
            case 3:
                for(Int_t Ridx = 0; Ridx < 3; ++Ridx)
                    for(Int_t Cidx = 0; Cidx < 3; ++Cidx)
                        tmp_fit.CMatrix_NS_EW_LS[Ridx][Cidx] = covMatrix_HS[Ridx][Cidx];
                break;
                    
            case 4:
                for(Int_t Ridx = 0; Ridx < 3; ++Ridx)
                    for(Int_t Cidx = 0; Cidx < 3; ++Cidx)
                        tmp_fit.CMatrix_NS_FB_LS[Ridx][Cidx] = covMatrix_HS[Ridx][Cidx];
                break;
                    
            case 5:
                for(Int_t Ridx = 0; Ridx < 3; ++Ridx)
                    for(Int_t Cidx = 0; Cidx < 3; ++Cidx)
                        tmp_fit.CMatrix_EW_FB_LS[Ridx][Cidx] = covMatrix_HS[Ridx][Cidx];
                break;
                    
            case 6:
                for(Int_t Ridx = 0; Ridx < 3; ++Ridx)
                    for(Int_t Cidx = 0; Cidx < 3; ++Cidx)
                        tmp_fit.CMatrix_Full_LS[Ridx][Cidx] = covMatrix_HS[Ridx][Cidx];
                break;
                    
        }
    }
    
    //////////////////////////////////////
    
    return result;
}



TH1* TemplateFitBH_test(
                            TH1* dat,
                            int ncomp,
                            TH1* temp[],
                            double _res[],
                            double _reserr[],
                            double initialguess[],
                            bool quiet,
                            bool resettemplates,
                            bool kClamping,
                            bool kChiSq,
                            std::ofstream &log_file,
                            test_fitResult &tmp_fit,
                            Int_t fit_element,
                            bool HS
                        )

{
    
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
            tempnorm[icomp] = (TH1*)(temp[icomp]->Clone("TemplateFit_tempnorm"));
            //tempnorm[icomp]->Scale(1.0/temp[icomp]->Integral("width"));//"width" makes a real pdf
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
    std::string parName[npar];
    double arglist[10];
    
    for (int icomp=0; icomp<npar; icomp++) {
        par[icomp] = 0;
        if (initialguess) par[icomp] = initialguess[icomp];
        stepSize[icomp] = 0.0001;
        
        //minVal[icomp] = -100;
        //maxVal[icomp] = 100;
        
        //minVal[icomp] = 0.0;
        //minVal[icomp] = -1.;         // Modified parameter low limit ! Before set at 0 !
        
        minVal[icomp] = -(dat->Integral()+7.0*sqrt(dat->Integral()));
        if (kClamping)
            maxVal[icomp] = dat->Integral();
        else
            maxVal[icomp] = dat->Integral()+7.0*sqrt(dat->Integral());
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
            
            
    for (int icomp=0; icomp<npar; icomp++)
    {
        minuit.GetParameter(icomp, _res[icomp], _reserr[icomp]);
        if (!quiet)
        {
            std::cout << "Parameter " << parName[icomp].c_str() << " value " << _res[icomp] << " error " << _reserr[icomp]<<std::endl;
            log_file << "Parameter " << parName[icomp].c_str() << " value " << _res[icomp] << " error " << _reserr[icomp]<<std::endl;
        }
    }
            
            /*
             for (int icomp=0; icomp<npar; icomp++)
             {
             minuit.GetParameter(icomp, _res[icomp], _reserr[icomp]);
             if (!quiet)
             {
             if(icomp==0)
             {
             std::cout << "Parameter " << parName[icomp].c_str() << " value " << _res[icomp] << " error " << _reserr[icomp]<<std::endl;
             log_file << "Parameter " << parName[icomp].c_str() << " value " << _res[icomp] << " error " << _reserr[icomp]<<std::endl;
             }
             else if(icomp==1)
             {
             std::cout << "Parameter " << parName[icomp].c_str() << " value " << _res[icomp]*TMath::Sqrt((3/(4*TMath::Pi()))) << " error " << _reserr[icomp]*TMath::Sqrt((3/(4*TMath::Pi()))) <<std::endl;
             log_file << "Parameter " << parName[icomp].c_str() << " value " << _res[icomp]*TMath::Sqrt((3/(4*TMath::Pi()))) << " error " << _reserr[icomp]*TMath::Sqrt((3/(4*TMath::Pi()))) <<std::endl;
             //std::cout << "Parameter " << parName[icomp].c_str() << " value " << _res[icomp]*TMath::Sqrt((3/(4*TMath::Pi()))) << " error " << _reserr[icomp] << std::endl;
             //log_file << "Parameter " << parName[icomp].c_str() << " value " << _res[icomp]*TMath::Sqrt((3/(4*TMath::Pi()))) << " error " << _reserr[icomp] << std::endl;
             }
             else
             {
             std::cout << "Parameter " << parName[icomp].c_str() << " value " << _res[icomp]*TMath::Sqrt((6/(4*TMath::Pi()))) << " error " << _reserr[icomp]*TMath::Sqrt((6/(4*TMath::Pi()))) <<std::endl;
             log_file << "Parameter " << parName[icomp].c_str() << " value " << _res[icomp]*TMath::Sqrt((6/(4*TMath::Pi()))) << " error " << _reserr[icomp]*TMath::Sqrt((6/(4*TMath::Pi()))) <<std::endl;
             //std::cout << "Parameter " << parName[icomp].c_str() << " value " << _res[icomp]*TMath::Sqrt((6/(4*TMath::Pi()))) << " error " << _reserr[icomp] << std::endl;
             //log_file << "Parameter " << parName[icomp].c_str() << " value " << _res[icomp]*TMath::Sqrt((6/(4*TMath::Pi()))) << " error " << _reserr[icomp]  << std::endl;
             }
             }
             }
             */
            
    double chisq = func(npar, _res, true);
    int freepars = minuit.GetNumFreePars();
    int ndof = NDATA-freepars;
    if (!quiet)
    {
        printf("chisq: %f/%d = %f\n", chisq, ndof, chisq/ndof);
        log_file << "\nchisq: " << chisq << "ndof: " << ndof << "chisq/ndof: " << chisq/ndof;
    }
            
            
    TH1* result = (TH1*)(dat->Clone(Form("%s_TemplateFitBH_result", dat->GetName())));
    result->Reset();
    for (int icomp=0; icomp<ncomp; icomp++) {
    result->Add(temp[icomp], _res[icomp]/temp[icomp]->Integral());
    }
            
    ncomp_old = ncomp;
            
    //  printf("</TemplateFitBH>\n");
            
            
    ///////////////// Filling class
            
    if(!HS)
    {
                
    
        tmp_fit.fit_par_LS[fit_element] = _res[0];
        tmp_fit.fit_err_LS[fit_element] = _reserr[0];
        
        tmp_fit.chi2_LS[fit_element] = chisq;
        tmp_fit.ndf_LS[fit_element] = ndof;
        tmp_fit.chi2_r_LS[fit_element] = chisq/ndof;
        tmp_fit.delta_LS[fit_element] = _res[0];
                
    }
    else
    {
                
        tmp_fit.fit_par_HS[fit_element] = _res[0];
        tmp_fit.fit_err_HS[fit_element] = _reserr[0];
                
        tmp_fit.chi2_HS[fit_element] = chisq;
        tmp_fit.ndf_HS[fit_element] = ndof;
        tmp_fit.chi2_r_HS[fit_element] = chisq/ndof;
        tmp_fit.delta_HS[fit_element] = _res[0];
                
    }
            
    //////////////////////////////////////
            
    return result;
    
}



void fcnchisq(int& npar, double* deriv, double& f, double par[], int flag){
    //flag=flag;
    //deriv=deriv;
    f=func(npar, par, true);
    return;
}

void fcnlike(int& npar, double* deriv, double& f, double par[], int flag){
    //flag=flag;
    //deriv=deriv;
    f=func(npar, par, false);
    return;
}


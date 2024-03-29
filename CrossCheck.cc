#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <cassert>

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TRandom.h"

#include "Rtypes.h"

#include "RooMsgService.h"
#include "RooFitResult.h"
#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooPlot.h"
#include "RooGlobalFunc.h"
#include "RooFit.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooGaussian.h"
#include "RooGenericPdf.h"
#include "RooAddPdf.h"
#include "RooExponential.h"
#include "RooPolynomial.h"
#include "RooDataSet.h"
#include "RooWorkspace.h"
#include "RooFitResult.h"
#include "RooChi2Var.h"
#include "RooNLLVar.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooBernstein.h"
#include "RooChebychev.h"
#include "RooRandom.h"

using namespace std;
using namespace RooFit;

const int nFuncs=13;
string funcNames[nFuncs] = {"2pol","3pol","4pol","5pol","1exp","2exp","3exp","1pow","2pow","3pow","2lau","4lau","6lau"};

void checkFunc(string name){
  if (name!="2pol" && name!="3pol" && name!="4pol" && name !="5pol" && name!="1exp" && name!="2exp" && name!="3exp" && name!="1pow" && name!="2pow" && name!="3pow" && name!="2lau" && name!="4lau" && name!="6lau") {
    cout << "Invalid function: " << name << endl;
    cout << "Options are: " << endl;
    for (int f=0; f<nFuncs; f++){
      cout << "   " << funcNames[f].c_str() << endl;
    }
    exit(1);
  }
}

const int getPar(string name){
  if (name=="1exp" || name=="1pow" || name=="2lau") return 1;
  else if (name=="2pol") return 2;
  else if (name=="3pol" || name=="2exp" || name=="2pow" || name=="4lau") return 3;
  else if (name=="4pol") return 4;
  else if (name=="5pol" || name=="3exp" || name=="3pow" || name=="6lau") return 5;
  else exit(1);
}

string getType(string name){
  if (name=="1exp" || name=="2exp" || name=="3exp") return "exp";
  else if (name=="2pol" || name=="3pol" || name=="4pol" || name=="5pol") return "pol";
  else if (name=="1pow" || name=="2pow" || name=="3pow") return "pow";
  else if (name=="2lau" || name=="4lau" || name=="6lau") return "lau";
  else exit(1);
}

int main(int argc, char* argv[]){

  bool help=false;
  bool verbose=false;
  bool doBkgInt=false;
  bool plotGen=false;
  bool saveDataFit=false;
  bool wideRange=false;
  int nToys;
  int nJobs=1;
  int jobNumb=0;
  int toyStep;
  string genName;
  string fitName;

  for (int arg=0; arg<argc; arg++) {
    if (string(argv[arg])=="-h" || string(argv[arg])=="--help") help=true;
    if (string(argv[arg])=="-v") verbose=true;
    if (string(argv[arg])=="-bkg") doBkgInt=true;
    if (string(argv[arg])=="-pG") plotGen=true;
    if (string(argv[arg])=="-sDF") saveDataFit=true;
    if (string(argv[arg])=="-gen") genName=string(argv[arg+1]);
    if (string(argv[arg])=="-fit") fitName=string(argv[arg+1]);
    if (string(argv[arg])=="-t") nToys=atoi(argv[arg+1]);
    if (string(argv[arg])=="-p") toyStep=atoi(argv[arg+1]);
    if (string(argv[arg])=="-nJ") nJobs=atoi(argv[arg+1]);
    if (string(argv[arg])=="-j") jobNumb=atoi(argv[arg+1]);
    if (string(argv[arg])=="-wide") wideRange=true;
  }
  if (!help){
    checkFunc(genName);
    checkFunc(fitName);
  }
  if (argc<9 || help){
    cout << "--- Run with following options: ---" << endl;
    cout << "    -t    nToys " << endl;
    cout << "    -p    plotStep " << endl;
    cout << "    -gen  $i to gen with func $i " << endl;
    cout << "    -fit  $i to fit with func $i " << endl;
    cout << "    -nJ   number of jobs " << endl;
    cout << "    -j    job number " << endl;
    cout << "--- Additional options: ---" << endl;
    cout << "    -v    for diagnostics " << endl;
    cout << "    -bkg  to do bkg int " << endl;
    cout << "    -pG   to plot gen func " << endl;
    cout << "    -sDF  to save data fit " << endl;
    cout << "    -wide to fit from 100-180 " << endl;
    exit(1);
  }

  const int nGenPar=getPar(genName);
  const int nFitPar=getPar(fitName);
  string genType=getType(genName);
  string fitType=getType(fitName);
  
  cout << "--- Running with following options ---" << endl;
  cout << "    nToys:            " << nToys << endl;
  cout << "    plotStep:         " << toyStep << endl;
  cout << "    genFunction:      " << genName << endl;
  cout << "    fitFunction:      " << fitName << endl;
  cout << "    nJobs:            " << nJobs << endl;
  cout << "    jobNum:           " << jobNumb << endl;
  if (verbose)     cout << "    print fit results on " << endl;
  if (doBkgInt)    cout << "    bkg integral on " << endl;
  if (plotGen)     cout << "    plot gen function on " << endl;
  if (saveDataFit) cout << "    save data fit on " << endl;
  if (wideRange)   cout << "    wide range fit on " << endl;
  if (saveDataFit && nToys>0){
    cout << " ERROR: CANNOT SAVE STARTING PARAMS AND FIT TOYS SIMULTANEOUSLY" << endl;
    cout << " --> either run without -sDF option or without -t or with -t 0 for no toys" << endl;
    exit(1);
  }

  int toysPerJob = nToys/nJobs;

  if (!verbose) RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  gROOT->SetStyle("Plain");
  gROOT->ForceStyle();
  
  TFile *inFile;
  if (!wideRange) inFile = new TFile("CMS-HGG_4686pb.root");
  else inFile = new TFile("CMS-HGG_4686pb_pol5_110_150_1.root");
  gROOT->SetStyle("Plain");
  gROOT->ForceStyle();
  // get workspace and mass data
  RooWorkspace *dataWS = (RooWorkspace*)inFile->Get("cms_hgg_workspace");
  RooRealVar *mass = (RooRealVar*)dataWS->var("CMS_hgg_mass");
  if (wideRange) {
    mass->setBins(320);
    mass->setBins(80,"chi2binning");
  }
  else {
    mass->setBins(240);
    mass->setBins(60,"chi2binning");
  }

  system("mkdir FitResults");
  //system("mkdir rm -r plots");
  system("mkdir -p plots/toys");
  system("mkdir -p plots/data");
  system("mkdir -p plots/histos");
  
  // --- declare variables for this code -----
  const int nCats=4;
  const int nMasses=8;
  const int nWinds=2;
  string winName[nWinds] = {"all","10"};
  double startPar[nFitPar][nCats];
  double swStartPar[nMasses][nCats][nFitPar];
  RooDataSet *data[nCats];
  RooRealVar *genPars[5][nCats];
  RooRealVar *fitPars[5][nCats];
  RooExponential *rooExp[6][nCats];
  RooAbsPdf *genFcn[nCats];
  RooAbsPdf *fitFcn[nCats];
  double sigSMEvents[nMasses][nCats];
  RooRealVar *mu = new RooRealVar("mu","mu",0.,-10,10);
  RooHistPdf *sigMCPdf[nMasses][nCats];
  RooRealVar *bkgYield[nMasses][nCats];
  RooFormulaVar *sigYield[nMasses][nCats];
  RooAddPdf *sigAndBkg[nMasses][nCats];

  TH1F *bias[nMasses][nCats];
  TH1F *muHist[nWinds][nMasses];

  TFile *outFile = new TFile(Form("SWResults/biasCheck_g%s_f%s_j%d.root",genName.c_str(),fitName.c_str(),jobNumb),"RECREATE");
  TFile *dataFitFile;
  if (saveDataFit) dataFitFile = new TFile("NewResults/fitsTodata.root","UPDATE");
  else {
    dataFitFile = new TFile("FitResults/fitsTodata.root");
    // --- get starting fit parameters
    /*
    for (int cat=0; cat<nCats; cat++){
      for (int par=0; par<nFitPar; par++){
        RooRealVar *datFitRes = (RooRealVar*)((RooFitResult*)dataFitFile->Get(Form("fitRes_%s_toData_cat%d",fitName.c_str(),cat)))->floatParsFinal().find(Form("%s_p%d_cat%d",fitType.c_str(),par,cat));
        startPar[par][cat] = datFitRes->getVal();
      }
    }
    */
  }
  
  int mIt=0;
  for (int mMC=110; mMC<=150; mMC+=5){
    if (mMC==145) continue;
    for (int win=0; win<nWinds; win++) muHist[win][mIt] = new TH1F(Form("muHist_g%s_f%s_m%d_w%s",genName.c_str(),fitName.c_str(),mMC,winName[win].c_str()),Form("muHist_g%s_f%s_m%d_w%s",genName.c_str(),fitName.c_str(),mMC,winName[win].c_str()),100,-10,10);
    if (doBkgInt) for (int cat=0; cat<nCats; cat++) bias[mIt][cat] = new TH1F(Form("b_g%s_f%s_m%d_c%d",genName.c_str(),fitName.c_str(),mMC,cat),Form("b_g%s_f%s_m%d_c%d",genName.c_str(),fitName.c_str(),mMC,cat),100,-100,100);
    mIt++;
  }
  double toData[nMasses][nCats];

  // ------ declare fit variables ------------
  cout << "Declaring variables for fit" << endl;
  for (int cat=0; cat<nCats; cat++){
    // ------------------------------------------
    // ----- params for initial fit to data -----
    // ------------------------------------------
    // --- polynomials ---
    if (genType=="pol"){
      for (int par=0; par<nGenPar; par++) genPars[par][cat] = new RooRealVar(Form("pol_p%d_cat%d",par,cat),Form("pol_p%d_cat%d",par,cat),0.06,-2.,2.); 
      if (genName=="2pol") genFcn[cat] = new RooChebychev(Form("2pol_cat%d",cat),Form("2pol_cat%d",cat),*mass,RooArgList(*genPars[0][cat],*genPars[1][cat]));
      if (genName=="3pol") genFcn[cat] = new RooChebychev(Form("3pol_cat%d",cat),Form("3pol_cat%d",cat),*mass,RooArgList(*genPars[0][cat],*genPars[1][cat],*genPars[2][cat]));
      if (genName=="4pol") genFcn[cat] = new RooChebychev(Form("4pol_cat%d",cat),Form("4pol_cat%d",cat),*mass,RooArgList(*genPars[0][cat],*genPars[1][cat],*genPars[2][cat],*genPars[3][cat]));
      if (genName=="5pol") genFcn[cat] = new RooChebychev(Form("5pol_cat%d",cat),Form("5pol_cat%d",cat),*mass,RooArgList(*genPars[0][cat],*genPars[1][cat],*genPars[2][cat],*genPars[3][cat],*genPars[4][cat]));
    }
    // --- exponentials ---
    if (genType=="exp"){
      for (int par=0; par<nGenPar; par+=2) genPars[par][cat] = new RooRealVar(Form("exp_p%d_cat%d",par,cat),Form("exp_p%d_cat%d",par,cat),-0.1,-2.,0.);
      for (int par=1; par<nGenPar; par+=2) genPars[par][cat] = new RooRealVar(Form("exp_p%d_cat%d",par,cat),Form("exp_p%d_cat%d",par,cat),0.5,0.,1.);
      if (nGenPar>0) rooExp[0][cat] = new RooExponential(Form("rooexp0_cat%d",cat),Form("rooexp0_cat%d",cat),*mass,*genPars[0][cat]);
      if (nGenPar>2) rooExp[1][cat] = new RooExponential(Form("rooexp1_cat%d",cat),Form("rooexp1_cat%d",cat),*mass,*genPars[2][cat]);
      if (nGenPar>4) rooExp[2][cat] = new RooExponential(Form("rooexp2_cat%d",cat),Form("rooexp2_cat%d",cat),*mass,*genPars[4][cat]);
      if (genName=="1exp") genFcn[cat] = (RooExponential*)rooExp[0][cat]->clone(Form("1exp_cat%d",cat)); 
      if (genName=="2exp") genFcn[cat] = new RooAddPdf(Form("2exp_cat%d",cat),Form("2exp_cat%d",cat),RooArgList(*rooExp[0][cat],*rooExp[1][cat]),RooArgList(*genPars[1][cat]));
      if (genName=="3exp") genFcn[cat] = new RooAddPdf(Form("3exp_cat%d",cat),Form("3exp_cat%d",cat),RooArgList(*rooExp[0][cat],*rooExp[1][cat],*rooExp[2][cat]),RooArgList(*genPars[1][cat],*genPars[3][cat]));
    }
    // --- power laws ---
    if (genType=="pow"){
      for (int par=0; par<nGenPar; par+=2) genPars[par][cat] = new RooRealVar(Form("pow_p%d_cat%d",par,cat),Form("pow_p%d_cat%d",par,cat),-1.,-10.,10.);
      for (int par=1; par<nGenPar; par+=2) genPars[par][cat] = new RooRealVar(Form("pow_p%d_cat%d",par,cat),Form("pow_p%d_cat%d",par,cat),0.8,0.,1.);
      if (genName=="1pow") genFcn[cat] = new RooGenericPdf(Form("1pow_cat%d",cat),Form("1pow_cat%d",cat),"pow(@0,@1)",RooArgList(*mass,*genPars[0][cat]));
      if (genName=="2pow") genFcn[cat] = new RooGenericPdf(Form("2pow_cat%d",cat),Form("2pow_cat%d",cat),"(1-@2)*pow(@0,@1)+@2*pow(@0,@3)",RooArgList(*mass,*genPars[0][cat],*genPars[1][cat],*genPars[2][cat]));
      if (genName=="3pow") genFcn[cat] = new RooGenericPdf(Form("3pow_cat%d",cat),Form("3pow_cat%d",cat),"(1-@2-@4)*pow(@0,@1)+@2*pow(@0,@3)+@4*pow(@0,@5)",RooArgList(*mass,*genPars[0][cat],*genPars[1][cat],*genPars[2][cat],*genPars[3][cat],*genPars[4][cat]));
    }
    // --- laurent series ---
    if (genType=="lau"){
      for (int par=0; par<nGenPar; par++) genPars[par][cat] = new RooRealVar(Form("lau_p%d_cat%d",par,cat),Form("lau_p%d_cat%d",par,cat),0.5,0.,1.);
      if (genName=="2lau") genFcn[cat] = new RooGenericPdf(Form("2lau_cat%d",cat),Form("2lau_cat%d",cat),"(1-@1)*pow(@0,-4.0)+@1*pow(@0,-5.0)",RooArgList(*mass,*genPars[0][cat]));
      if (genName=="4lau") genFcn[cat] = new RooGenericPdf(Form("4lau_cat%d",cat),Form("4lau_cat%d",cat),"(1-@1-@2-@3)*pow(@0,-4.0)+@1*pow(@0,-5.0)+@2*pow(@0,-3.0)+@3*pow(@0,-6.0)",RooArgList(*mass,*genPars[0][cat],*genPars[1][cat],*genPars[2][cat]));
      if (genName=="6lau") genFcn[cat] = new RooGenericPdf(Form("6lau_cat%d",cat),Form("6lau_cat%d",cat),"(1-@1-@2-@3-@4-@5)*pow(@0,-4.0)+@1*pow(@0,-5.0)+@2*pow(@0,-3.0)+@3*pow(@0,-6.0)+@4*pow(@0,-2.0)+@5*pow(@0,-7.0)",RooArgList(*mass,*genPars[0][cat],*genPars[1][cat],*genPars[2][cat],*genPars[3][cat],*genPars[4][cat]));
    }

    // ------------------------------------------
    // --- params for fit to gen data ---
    // ------------------------------------------
    // --- polynomials ---
    if (fitType=="pol"){
      for (int par=0; par<nFitPar; par++) fitPars[par][cat] = new RooRealVar(Form("polF_p%d_cat%d",par,cat),Form("polF_p%d_cat%d",par,cat),0.01,-2.,2.); 
      if (fitName=="2pol") fitFcn[cat] = new RooChebychev(Form("2polF_cat%d",cat),Form("2polF_cat%d",cat),*mass,RooArgList(*fitPars[0][cat],*fitPars[1][cat]));
      if (fitName=="3pol") fitFcn[cat] = new RooChebychev(Form("3polF_cat%d",cat),Form("3polF_cat%d",cat),*mass,RooArgList(*fitPars[0][cat],*fitPars[1][cat],*fitPars[2][cat]));
      if (fitName=="4pol") fitFcn[cat] = new RooChebychev(Form("4polF_cat%d",cat),Form("4polF_cat%d",cat),*mass,RooArgList(*fitPars[0][cat],*fitPars[1][cat],*fitPars[2][cat],*fitPars[3][cat]));
      if (fitName=="5pol") fitFcn[cat] = new RooChebychev(Form("5polF_cat%d",cat),Form("5polF_cat%d",cat),*mass,RooArgList(*fitPars[0][cat],*fitPars[1][cat],*fitPars[2][cat],*fitPars[3][cat],*fitPars[4][cat]));
    }
    // --- exponentials ---
    if (fitType=="exp"){
      for (int par=0; par<nFitPar; par+=2) fitPars[par][cat] = new RooRealVar(Form("expF_p%d_cat%d",par,cat),Form("expF_p%d_cat%d",par,cat),-0.1,-1.,0.);
      for (int par=1; par<nFitPar; par+=2) fitPars[par][cat] = new RooRealVar(Form("expF_p%d_cat%d",par,cat),Form("expF_p%d_cat%d",par,cat),0.1,0.,1.);
      if (nFitPar>0) rooExp[3][cat] = new RooExponential(Form("rooexpF3_cat%d",cat),Form("rooexpF3_cat%d",cat),*mass,*fitPars[0][cat]);
      if (nFitPar>2) rooExp[4][cat] = new RooExponential(Form("rooexpF4_cat%d",cat),Form("rooexpF4_cat%d",cat),*mass,*fitPars[2][cat]);
      if (nFitPar>4) rooExp[5][cat] = new RooExponential(Form("rooexpF5_cat%d",cat),Form("rooexpF5_cat%d",cat),*mass,*fitPars[4][cat]);
      if (fitName=="1exp") fitFcn[cat] = (RooExponential*)rooExp[3][cat]->clone(Form("1expF_cat%d",cat)); 
      if (fitName=="2exp") fitFcn[cat] = new RooAddPdf(Form("2expF_cat%d",cat),Form("2expF_cat%d",cat),RooArgList(*rooExp[3][cat],*rooExp[4][cat]),RooArgList(*fitPars[1][cat]));
      if (fitName=="3exp") fitFcn[cat] = new RooAddPdf(Form("3expF_cat%d",cat),Form("3expF_cat%d",cat),RooArgList(*rooExp[3][cat],*rooExp[4][cat],*rooExp[5][cat]),RooArgList(*fitPars[1][cat],*fitPars[3][cat]));
    }
    // --- power laws ---
    if (fitType=="pow"){
      for (int par=0; par<nFitPar; par+=2) fitPars[par][cat] = new RooRealVar(Form("powF_p%d_cat%d",par,cat),Form("powF_p%d_cat%d",par,cat),-1.,-10.,10.);
      for (int par=1; par<nFitPar; par+=2) fitPars[par][cat] = new RooRealVar(Form("powF_p%d_cat%d",par,cat),Form("powF_p%d_cat%d",par,cat),0.8,0.,1.);
      if (fitName=="1pow") fitFcn[cat] = new RooGenericPdf(Form("1powF_cat%d",cat),Form("1powF_cat%d",cat),"pow(@0,@1)",RooArgList(*mass,*fitPars[0][cat]));
      if (fitName=="2pow") fitFcn[cat] = new RooGenericPdf(Form("2powF_cat%d",cat),Form("2powF_cat%d",cat),"(1-@2)*pow(@0,@1)+@2*pow(@0,@3)",RooArgList(*mass,*fitPars[0][cat],*fitPars[1][cat],*fitPars[2][cat]));
      if (fitName=="3pow") fitFcn[cat] = new RooGenericPdf(Form("3powF_cat%d",cat),Form("3powF_cat%d",cat),"(1-@2-@4)*pow(@0,@1)+@2*pow(@0,@3)+@4*pow(@0,@5)",RooArgList(*mass,*fitPars[0][cat],*fitPars[1][cat],*fitPars[2][cat],*fitPars[3][cat],*fitPars[4][cat]));
    }
    // --- laurent series ---
    if (fitType=="lau"){
      for (int par=0; par<nFitPar; par++) fitPars[par][cat] = new RooRealVar(Form("lauF_p%d_cat%d",par,cat),Form("lauF_p%d_cat%d",par,cat),0.5,0.,1.);
      if (fitName=="2lau") fitFcn[cat] = new RooGenericPdf(Form("2lauF_cat%d",cat),Form("2lauF_cat%d",cat),"(1-@1)*pow(@0,-4.0)+@1*pow(@0,-5.0)",RooArgList(*mass,*fitPars[0][cat]));
      if (fitName=="4lau") fitFcn[cat] = new RooGenericPdf(Form("4lauF_cat%d",cat),Form("4lauF_cat%d",cat),"(1-@1-@2-@3)*pow(@0,-4.0)+@1*pow(@0,-5.0)+@2*pow(@0,-3.0)+@3*pow(@0,-6.0)",RooArgList(*mass,*fitPars[0][cat],*fitPars[1][cat],*fitPars[2][cat]));
      if (fitName=="6lau") fitFcn[cat] = new RooGenericPdf(Form("6lauF_cat%d",cat),Form("6lauF_cat%d",cat),"(1-@1-@2-@3-@4-@5)*pow(@0,-4.0)+@1*pow(@0,-5.0)+@2*pow(@0,-3.0)+@3*pow(@0,-6.0)+@4*pow(@0,-2.0)+@5*pow(@0,-7.0)",RooArgList(*mass,*fitPars[0][cat],*fitPars[1][cat],*fitPars[2][cat],*fitPars[3][cat],*fitPars[4][cat]));
    }
  
  // ---- get data in each category
  data[cat] = (RooDataSet*)dataWS->data(Form("data_mass_cat%d",cat));
  // ---- try and estimate starting pars for SW by fitting bkg only to data ----
    for (int win=1; win<nWinds; win++){
      mIt=0;
      for (int mMC=110; mMC<=150; mMC+=5){
        if (mMC==145) continue;
        mass->setRange("testSW",mMC-10.,mMC+10.);
        // --- get copy of bkg func ---
        RooAbsPdf *tempFitFcn;
        if (fitType=="pol") tempFitFcn = (RooChebychev*)fitFcn[cat]->Clone(Form("swFitToData_%s_w%s_m%d_c%d",fitName.c_str(),winName[win].c_str(),mMC,cat));
        if (fitName=="1exp") tempFitFcn = (RooExponential*)fitFcn[cat]->Clone(Form("swFitToData_%s_w%s_m%d_c%d",fitName.c_str(),winName[win].c_str(),mMC,cat));
        if (fitName=="2exp" || fitName=="3exp") tempFitFcn = (RooAddPdf*)fitFcn[cat]->Clone(Form("swFitToData_%s_w%s_m%d_c%d",fitName.c_str(),winName[win].c_str(),mMC,cat));
        if (fitType=="pow" || fitType=="lau") tempFitFcn = (RooGenericPdf*)fitFcn[cat]->Clone(Form("swFitToData_%s_w%s_m%d_c%d",fitName.c_str(),winName[win].c_str(),mMC,cat));
        
        RooFitResult *swFit = tempFitFcn->fitTo(*data[cat],Range("testSW"),SumCoefRange("testSW"),Save(),PrintLevel(-1),Warnings(false),PrintEvalErrors(-1));
        swFit->floatParsFinal().Print("s");
        for (int par=0; par<nFitPar; par++) swStartPar[mIt][cat][par]=fitPars[par][cat]->getVal();
        // ---- draw SW fit to data
        TCanvas *c2 = new TCanvas();
        RooPlot *swFrame = mass->frame(Title(Form("%s fit to data sw%s m%d cat%d",fitName.c_str(),winName[win].c_str(),mMC,cat)));
        data[cat]->plotOn(swFrame,DataError(RooDataSet::SumW2),Binning("chi2binning"));
        tempFitFcn->plotOn(swFrame);
        swFrame->SetXTitle("m_{#gamma#gamma} (GeV/c^{2})");
        swFrame->Draw();
        double chi2 = swFrame->chiSquare(nFitPar);
        double nBins;
        if (wideRange) nBins=80;
        else nBins=60;
        double prob = TMath::Prob(chi2*(nBins-nFitPar),(nBins-nFitPar));
        TPaveText *text = new TPaveText(0.7,0.6,0.85,0.89,"NDC");
        text->SetFillColor(0);
        text->SetLineColor(0);
        text->AddText(Form("#chi^{2} = %1.2f",chi2));
        text->AddText(Form("prob(#chi^{2}) = %1.2f",prob));
        for (int par=0; par<nGenPar; par++) text->AddText(Form("p%d = %1.2f",par,genPars[par][cat]->getVal()));
        text->SetTextAlign(13);
        text->SetTextSize(0.04);
        text->Draw("same");
        c2->Print(Form("plots/data/fitTodat_%s_sw%s_m%d_c%d.png",fitName.c_str(),winName[win].c_str(),mMC,cat));
        mIt++;
        //delete text;
        //delete c2;
      }
    }

  // ----- fit to data in each category ----- 
    RooFitResult *datFit;
    if (verbose) datFit = genFcn[cat]->fitTo(*data[cat],Save());
    else datFit = genFcn[cat]->fitTo(*data[cat],Save(),PrintLevel(-1),Warnings(false),PrintEvalErrors(-1));
    datFit->floatParsFinal().Print("s");
    datFit->SetName(Form("fitRes_%s_toData_cat%d",genName.c_str(),cat));
    if (saveDataFit){
      dataFitFile->cd();
      datFit->Write();
    }
    outFile->cd();
    datFit->Write();
  // ----- draw fit to data ---------------
    TCanvas *c = new TCanvas();
    RooPlot *tFrame = mass->frame(Title(Form("%s fit to data cat%d",genName.c_str(),cat)));
    data[cat]->plotOn(tFrame,DataError(RooDataSet::SumW2),Binning("chi2binning"));
    genFcn[cat]->plotOn(tFrame);
    tFrame->SetXTitle("m_{#gamma#gamma} (GeV/c^{2})");
    tFrame->Draw();
    double chi2 = tFrame->chiSquare(nGenPar);
    double nBins;
    if (wideRange) nBins=80;
    else nBins=60;
    double prob = TMath::Prob(chi2*(nBins-nGenPar),(nBins-nGenPar));
    TPaveText *text = new TPaveText(0.7,0.6,0.85,0.89,"NDC");
    text->SetFillColor(0);
    text->SetLineColor(0);
    text->AddText(Form("#chi^{2} = %1.2f",chi2));
    text->AddText(Form("prob(#chi^{2}) = %1.2f",prob));
    for (int par=0; par<nGenPar; par++) text->AddText(Form("p%d = %1.2f",par,genPars[par][cat]->getVal()));
    text->SetTextAlign(13);
    text->SetTextSize(0.04);
    text->Draw("same");
    c->Print(Form("plots/data/fitTodat_%s_cat%d.png",genName.c_str(),cat),"png");
    //delete c;
  
  // ---- Get mass distributions  and entries ------
    RooRealVar intRange(*mass);
    mIt=0;
    cout << "--- Event check: ----- " << endl;
    for (int mMC=110; mMC<=150; mMC+=5){
      if (mMC==145) continue;
      double mLow = mMC-5.;
      double mHigh = mMC+5.;
      intRange.setRange("sigWindow",mLow,mHigh);
    // ----- calc integral around mass ---------------
      if (doBkgInt){
        RooAbsReal *onData = genFcn[cat]->createIntegral(*mass,NormSet(*mass),Range("sigWindow"));
        toData[mMC][cat] = onData->getVal()*data[cat]->numEntries();
      }
    // get signal data
      RooDataSet *sigData = (RooDataSet*)dataWS->data(Form("sig_ggh_mass_m%d_cat%d",mMC,cat));
      sigData->append(*((RooDataSet*)dataWS->data(Form("sig_vbf_mass_m%d_cat%d",mMC,cat))));
      sigData->append(*((RooDataSet*)dataWS->data(Form("sig_wzh_mass_m%d_cat%d",mMC,cat))));
      sigData->append(*((RooDataSet*)dataWS->data(Form("sig_tth_mass_m%d_cat%d",mMC,cat))));
    // get expected SM events
      sigSMEvents[mIt][cat] = (((TH1F*)inFile->Get(Form("th1f_sig_ggh_mass_m%d_cat%d",mMC,cat)))->Integral())+((TH1F*)inFile->Get(Form("th1f_sig_vbf_mass_m%d_cat%d",mMC,cat)))->Integral()+((TH1F*)inFile->Get(Form("th1f_sig_wzh_mass_m%d_cat%d",mMC,cat)))->Integral()+((TH1F*)inFile->Get(Form("th1f_sig_tth_mass_m%d_cat%d",mMC,cat)))->Integral();
      cout << "   Mass: " << mMC << " cat: " << cat << " intSM: " << sigSMEvents[mIt][cat] << endl;
    // construct s+b model
      RooDataHist *sigMC = new RooDataHist(Form("sigMC_m%d_cat%d",mMC,cat),Form("sigMC_m%d_cat%d",mMC,cat),*mass,*sigData);
      sigMCPdf[mIt][cat] = new RooHistPdf(Form("sigMCPdf_m%d_cat%d",mMC,cat),Form("sigMC_m%d_cat%d",mMC,cat),*mass,*sigMC);
   // ----- def s and b yields and construct s+b model
      bkgYield[mIt][cat] = new RooRealVar(Form("bkgYield_m%d_cat%d",mMC,cat),Form("bkgYield_m%d_cat%d",mMC,cat),500,200,8000);
      sigYield[mIt][cat] = new RooFormulaVar(Form("sigYield_m%d_cat%d",mMC,cat),Form("sigYield_m%d_cat%d",mMC,cat),"@0*@1",RooArgList(RooConst(sigSMEvents[mIt][cat]),*mu));
      sigAndBkg[mIt][cat] = new RooAddPdf(Form("SandB_m%d_cat%d",mMC,cat),Form("SandB_m%d_cat%d",mMC,cat),RooArgList(*fitFcn[cat],*sigMCPdf[mIt][cat]),RooArgList(*bkgYield[mIt][cat],*sigYield[mIt][cat]));
      mIt++;
    }
  }
  cout << "------------------------------------------------" << endl;
  cout << "--- Data fitted. Mass distributions obtained ---" << endl;
  cout << "------------------------------------------------" << endl;
  cout << "-------- Generating and fitting toys -----------" << endl;

  RooDataSet *genDatUnBin[nCats];
  RooDataHist *genDat[nCats];
  RooRandom::randomGenerator()->SetSeed(0);
  for (int itToy=jobNumb*toysPerJob; itToy<(jobNumb+1)*toysPerJob; itToy++){
    cout << "----------------------------------------------" << endl;
    cout << "------------- TOY: " << itToy << "------------" << endl;
    cout << "----------------------------------------------" << endl;
    for (int cat=0; cat<nCats; cat++){
      genDatUnBin[cat] = genFcn[cat]->generate(*mass,data[cat]->numEntries(),Extended());
      genDat[cat] = new RooDataHist("gen","gen",*mass,*genDatUnBin[cat]);
      // --- set starting to vals to that of gen
      //for (int par=0; par<nFitPar; par++) fitPars[par][cat]->setVal(startPar[par][cat]);
      //for (int mIt=0; mIt<nMasses; mIt++) bkgYield[mIt][cat]->setVal(500);
      mu->setVal(0.);
    }
    // --- set up categories and combine data ---- 
    RooCategory category("category","category");
    category.defineType("cat0");
    category.defineType("cat1");
    category.defineType("cat2");
    category.defineType("cat3");
    RooDataHist combData(Form("combData_toy%d",itToy),Form("combData_toy%d",itToy),*mass,Index(category),Import("cat0",*genDat[0]),Import("cat1",*genDat[1]),Import("cat2",*genDat[2]),Import("cat3",*genDat[3]));
    
    // ---- construct RooSimultaneous from 4 cats -----
    mIt=0;
    for (int mMC=110; mMC<=150; mMC+=5){
      if (mMC==145) continue;
      for (int win=1; win<nWinds; win++){
        if (win==0 && wideRange) mass->setRange("SW",100,180);
        else if (win==0 && !wideRange) mass->setRange("SW",100,160);
        else if (win>0) mass->setRange("SW",mMC-10,mMC+10);
        RooSimultaneous simPdf(Form("simPdf%d_toy%d",mMC,itToy),Form("simPdf%d_toy%d",mMC,itToy),category);
        RooAddPdf *SandB[nCats];
        for (int cat=0; cat<nCats; cat++) SandB[cat] = (RooAddPdf*)sigAndBkg[mIt][cat]->Clone(Form("SandB_m%d_w%s_c%d",mMC,winName[win].c_str(),cat));
        simPdf.addPdf(*SandB[0],"cat0");
        simPdf.addPdf(*SandB[1],"cat1");
        simPdf.addPdf(*SandB[2],"cat2");
        simPdf.addPdf(*SandB[3],"cat3");
        // ---- reset starting vals ----
        mu->setVal(0.);
        for (int cat=0; cat<nCats; cat++){
          if (win==0) {
            bkgYield[mIt][cat]->setVal(5000.);
            for (int p=0; p<nFitPar; p++){
              if (fitType=="pol") fitPars[p][cat]->setVal(0.);
              if (fitType=="exp") {
                if (p%2==0) fitPars[p][cat]->setVal(-0.1);
                else fitPars[p][cat]->setVal(0.5);
              }
              if (fitType=="pow") {
                if (p%2==0) fitPars[p][cat]->setVal(-1.);
                else fitPars[p][cat]->setVal(0.5);
              }
              if (fitType=="lau") fitPars[p][cat]->setVal(0.5);
            }
          }
          else {
            bkgYield[mIt][cat]->setVal(1000.);
            for (int p=0; p<nFitPar; p++) fitPars[p][cat]->setVal(swStartPar[mIt][cat][p]);
          }
        }
        // ---- fit and save output -----
        RooFitResult *res;
        if (verbose) res = simPdf.fitTo(combData,Range("SW"),SumCoefRange("SW"),Save(true));
        else res = simPdf.fitTo(combData,Range("SW"),SumCoefRange("SW"),Save(true),PrintLevel(-1),Warnings(false),PrintEvalErrors(-1));
        cout << "------------------------------------------------" << endl;
        cout << "--- TOY: " << itToy << " Gen: " << genName << " Fit: " << fitName << " mass: " << mMC << " window: " << winName[win] << endl;
        cout << "------------------------------------------------" << endl;
        res->floatParsFinal().Print("s");
        outFile->cd();
        res->SetName(Form("fitRes_g%s_f%s_m%d_w%s_t%d",genName.c_str(),fitName.c_str(),mMC,winName[win].c_str(),itToy));
        //res->Write();
        if (itToy>0) muHist[win][mIt]->Fill(mu->getVal()); 
        // --- find bkg integral -----
        for (int cat=0; cat<nCats; cat++){
          if (doBkgInt){
            RooAbsReal *intInRange = fitFcn[cat]->createIntegral(*mass,NormSet(*mass),Range("sigWindow"));
            double toBkg = intInRange->getVal()*genDat[cat]->numEntries();
            bias[mIt][cat]->Fill(toData[mIt][cat]-toBkg);
          }
        // ---- make plots ----
          if (itToy%toyStep==0 && itToy>0){
            TCanvas *c1 = new TCanvas();
            TLegend *leg = new TLegend(0.65,0.65,0.89,0.89);
            leg->SetLineColor(0);
            leg->SetFillColor(0);
            RooPlot *mFrame = mass->frame(Title(Form("Gen: %s. Fit: %s. Mass %d win %s cat %d toy %d",genName.c_str(),fitName.c_str(),mMC,winName[win].c_str(),cat,itToy)));
            if (wideRange) mass->setBins(160,"coarse");
            else mass->setBins(120,"coarse");
            genDat[cat]->plotOn(mFrame,DataError(RooDataSet::SumW2),Binning("coarse"));
            if (plotGen) genFcn[cat]->plotOn(mFrame,LineColor(kMagenta));
            //fitFcn[cat]->plotOn(mFrame,LineColor(kRed),LineStyle(kDashed));
            SandB[cat]->plotOn(mFrame,LineColor(kRed),LineStyle(kDashed),Components(*(SandB[cat]->pdfList().at(0))),Range(100,160));
            SandB[cat]->plotOn(mFrame);
            mFrame->SetXTitle("m_{#gamma#gamma} (GeV/c^{2})");
            mFrame->Draw();
        // ---- make legend -----
            TPaveText *text = new TPaveText(0.55,0.7,0.65,0.89,"NDC");
            text->SetLineColor(0);
            text->SetFillColor(0);
            text->AddText(Form("#mu = %1.2f",mu->getVal()));
            TH1F *h = new TH1F("h","h",1,0,1);
            h->SetLineColor(kMagenta);
            h->SetLineWidth(3);
            TH1F *h1 = new TH1F("h1","h2",1,0,1);
            h1->SetLineColor(kBlue);
            h1->SetLineWidth(3);
            TH1F *h2 = new TH1F("h2","h2",1,0,1);
            h2->SetLineColor(kRed);
            h2->SetLineStyle(kDashed);
            h2->SetLineWidth(3);
            if (plotGen) leg->AddEntry(h,"Truth","l");
            leg->AddEntry(h2,"bkg part of s+b","l");
            leg->AddEntry(h1,"s+b after fit","l");
            leg->Draw("same");
            text->Draw("same");
            c1->Print(Form("plots/toys/fitTogen_g%s_f%s_m%d_w%s_cat%d_toy%d.png",genName.c_str(),fitName.c_str(),mMC,winName[win].c_str(),cat,itToy),"png");
            
            //delete c1;
            //delete h;
            //delete h1;
            //delete h2;
            //delete text;
            //delete leg;
          }
        // --- set starting to vals to that of gen
          //for (int par=0; par<nFitPar; par++) fitPars[par][cat]->setVal(startPar[par][cat]);
          //for (int mIt=0; mIt<nMasses; mIt++) bkgYield[mIt][cat]->setVal(5000);
        }
        mu->setVal(0.);
      }
      mIt++;
    }
  }

  outFile->cd();
  mIt=0;
  for (int mMC=110; mMC<=150; mMC+=5){
    if (mMC==145) continue;
    for (int win=0; win<nWinds; win++){
      muHist[win][mIt]->Write();
      cout << "m" << mMC << " w" << win << " mean: " << muHist[win][mIt]->GetMean() << endl;
      TCanvas *canv = new TCanvas();
      muHist[win][mIt]->Draw();
      canv->Print(Form("plots/histos/mu_g%s_f%s_m%d_w%s_j%d.png",genName.c_str(),fitName.c_str(),mMC,winName[win].c_str(),jobNumb),"png");
      if (doBkgInt){
        for (int cat=0; cat<nCats; cat++){
          bias[mIt][cat]->Write();
          bias[mIt][cat]->Draw();
          canv->Print(Form("plots/histos/biasCheck_g%s_f%s_m%d_cat%d.png",genName.c_str(),fitName.c_str(),mMC,cat),"png");
          canv->Clear();
        }
      }
      //delete canv;
    }
    mIt++;
  }

  ofstream complete(Form("SWResults/g%s_f%s_j%d.txt",genName.c_str(),fitName.c_str(),jobNumb));
  complete << "Job: " << jobNumb << "/" << nJobs << " gen " << genName << " fit " << fitName << " completed successfully" << endl;
  complete.close();
  cout << "Text file written " << endl;
  outFile->Close();
  inFile->Close();

  return 0;
}

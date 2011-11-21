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

using namespace std;
using namespace RooFit;

int main(int argc, char* argv[]){

  const int nToys=atoi(argv[1]);
  const int toyStep=atoi(argv[2]);

//  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  gROOT->SetStyle("Plain");
  gROOT->ForceStyle();
  
  TFile *inFile = new TFile("CMS-HGG_4686pb.root");
  gROOT->SetStyle("Plain");
  gROOT->ForceStyle();
  // get workspace and mass data
  RooWorkspace *dataWS = (RooWorkspace*)inFile->Get("cms_hgg_workspace");
  RooRealVar *mass = (RooRealVar*)dataWS->var("CMS_hgg_mass");
  mass->setBins(240,"fineBins");

  TFile *outFile = new TFile("biasCheck.root","RECREATE");
  
  system("mkdir plots");
  const int nCats=4;
  const int nMasses=8;
  RooDataSet *data[nCats];
  RooRealVar *pol[6][nCats];
  RooGenericPdf *genFcn[nCats];
  RooGenericPdf *fitFcn[nCats];
  double sigSMEvents[nMasses][nCats];
  RooRealVar *mu = new RooRealVar("mu","mu",0.,-10,10);
  RooHistPdf *sigMCPdf[nMasses][nCats];
  RooRealVar *bkgYield[nMasses][nCats];
  RooFormulaVar *sigYield[nMasses][nCats];
  RooAddPdf *sigAndBkg[nMasses][nCats];

  TH1F *bias[nMasses][nCats];
  TH1F *muHist[nMasses];
  int mIt=0;
  for (int mMC=110; mMC<=150; mMC+=5){
    if (mMC==145) continue;
    muHist[mIt] = new TH1F(Form("muHist_m%d",mMC),Form("muHist_m%d",mMC),100,-10,10);
    for (int cat=0; cat<nCats; cat++) bias[mIt][cat] = new TH1F(Form("b_m%d_c%d",mMC,cat),Form("b_m%d_c%d",mMC,cat),100,-100,100);
    mIt++;
  }
  double toData[nCats];

  // ------ declare fit variables ------------
  cout << "Declaring variables for fit" << endl;
  for (int cat=0; cat<nCats; cat++){
    pol[0][cat] = new RooRealVar(Form("pol_p0_cat%d",cat),Form("pol_p0_cat%d",cat),-10.,10.); 
    pol[1][cat] = new RooRealVar(Form("pol_p1_cat%d",cat),Form("pol_p1_cat%d",cat),-10.,10.); 
    pol[2][cat] = new RooRealVar(Form("pol_p2_cat%d",cat),Form("pol_p2_cat%d",cat),-10.,10.); 
    genFcn[cat] = new RooGenericPdf(Form("3pol_cat%d",cat),Form("3pol_cat%d",cat),"(@1*@0)+(@2*pow(@0,2.0))+(@3*pow(@0,3.0))",RooArgSet(*mass,*pol[0][cat],*pol[1][cat],*pol[2][cat]));
    pol[3][cat] = new RooRealVar(Form("polF_p3_cat%d",cat),Form("polF_p3_cat%d",cat),2.,-10.,10.); 
    pol[4][cat] = new RooRealVar(Form("polF_p4_cat%d",cat),Form("polF_p4_cat%d",cat),0.25,-10.,10.); 
    pol[5][cat] = new RooRealVar(Form("polF_p5_cat%d",cat),Form("polF_p5_cat%d",cat),0.01,-10.,10.); 
    fitFcn[cat] = new RooGenericPdf(Form("3polF_cat%d",cat),Form("3polF_cat%d",cat),"(@1*@0)+(@2*pow(@0,2.0))+(@3*pow(@0,3.0))",RooArgSet(*mass,*pol[3][cat],*pol[4][cat],*pol[5][cat]));
  // ----- fit to data in each category ----- 
    pol[0][cat]->setVal(0.5);
    pol[1][cat]->setVal(0.5);
    pol[2][cat]->setVal(0.5);
    data[cat] = (RooDataSet*)dataWS->data(Form("data_mass_cat%d",cat));
    RooFitResult *datFit = genFcn[cat]->fitTo(*data[cat],Save());//,PrintLevel(-1),Warnings(false),PrintEvalErrors(-1));
    datFit->floatParsFinal().Print("s");
  // ----- calc integral around mass ---------------
    RooRealVar intRange(*mass);
    intRange.setRange("sigWindow",120,130);
    RooAbsReal *onData = genFcn[cat]->createIntegral(*mass,NormSet(*mass),Range("sigWindow"));
    toData[cat] = onData->getVal()*data[cat]->numEntries();
  // ----- draw fit to data ---------------
    TCanvas *c = new TCanvas();
    RooPlot *tFrame = mass->frame(Title(Form("2pol fit to data %d",cat)));
    data[cat]->plotOn(tFrame,DataError(RooDataSet::SumW2));
    genFcn[cat]->plotOn(tFrame);
    tFrame->Draw();
    c->Print(Form("plots/fitTodat_cat%d.png",cat),"png");
    delete c;
  
  // ---- Get mass distributions  and entries ------
    mIt=0;
    for (int mMC=110; mMC<=150; mMC+=5){
      if (mMC==145) continue;
      RooDataSet *sigData = (RooDataSet*)dataWS->data(Form("sig_ggh_mass_m%d_cat%d",mMC,cat));
      sigData->append(*((RooDataSet*)dataWS->data(Form("sig_vbf_mass_m%d_cat%d",mMC,cat))));
      sigData->append(*((RooDataSet*)dataWS->data(Form("sig_wzh_mass_m%d_cat%d",mMC,cat))));
      sigData->append(*((RooDataSet*)dataWS->data(Form("sig_tth_mass_m%d_cat%d",mMC,cat))));
      sigSMEvents[mIt][cat] = (((TH1F*)inFile->Get(Form("th1f_sig_ggh_mass_m%d_cat%d",mMC,cat)))->Integral())+((TH1F*)inFile->Get(Form("th1f_sig_vbf_mass_m%d_cat%d",mMC,cat)))->Integral()+((TH1F*)inFile->Get(Form("th1f_sig_wzh_mass_m%d_cat%d",mMC,cat)))->Integral()+((TH1F*)inFile->Get(Form("th1f_sig_tth_mass_m%d_cat%d",mMC,cat)))->Integral();
      cout << "--- Event check: ----- " << endl;
      cout << "Mass: " << mMC << " cat: " << cat << " intSM: " << sigSMEvents[mIt][cat] << endl;
      RooDataHist *sigMC = new RooDataHist(Form("sigMC_m%d_cat%d",mMC,cat),Form("sigMC_m%d_cat%d",mMC,cat),*mass,*sigData);
      sigMCPdf[mIt][cat] = new RooHistPdf(Form("sigMCPdf_m%d_cat%d",mMC,cat),Form("sigMC_m%d_cat%d",mMC,cat),*mass,*sigMC);
   // ----- def s and b yields and construct s+b model
      bkgYield[mIt][cat] = new RooRealVar(Form("bkgYield_m%d_cat%d",mMC,cat),Form("bkgYield_m%d_cat%d",mMC,cat),5000,3000,8000);
      sigYield[mIt][cat] = new RooFormulaVar(Form("sigYield_m%d_cat%d",mMC,cat),Form("sigYield_m%d_cat%d",mMC,cat),"@0*@1",RooArgList(RooConst(sigSMEvents[mIt][cat]),*mu));
      //RooRealVar *sigYield = new RooRealVar(Form("sigYield_m%d_cat%d",mMC,cat),Form("sigYield_m%d_cat%d",mMC,cat),0,-35,35);
      sigAndBkg[mIt][cat] = new RooAddPdf(Form("SandB_m%d_cat%d",mMC,cat),Form("SandB_m%d_cat%d",mMC,cat),RooArgList(*fitFcn[cat],*sigMCPdf[mIt][cat]),RooArgList(*bkgYield[mIt][cat],*sigYield[mIt][cat]));
      mIt++;
    }
  }
  cout << "------------------------------------------------" << endl;
  cout << "--- Data fitted. Mass distributions obtained ---" << endl;
  cout << "------------------------------------------------" << endl;

  RooDataSet *genDat[nCats];

  for (int itToy=0; itToy<nToys; itToy++){
    cout << "----------------------------------------" << endl;
    cout << "------------- TOY: " << itToy << "------" << endl;
    cout << "----------------------------------------" << endl;
    for (int cat=0; cat<nCats; cat++){
      genDat[cat] = genFcn[cat]->generate(*mass,data[cat]->numEntries(),Extended());
      pol[3][cat]->setVal(pol[0][cat]->getVal());
      pol[4][cat]->setVal(pol[1][cat]->getVal());
      pol[5][cat]->setVal(pol[2][cat]->getVal());
      for (int mIt=0; mIt<nMasses; mIt++) bkgYield[mIt][cat]->setVal(5000);
      mu->setVal(0.);
    }
    // --- set up categories and combine data ---- 
    RooCategory category("category","category");
    category.defineType("cat0");
    category.defineType("cat1");
    category.defineType("cat2");
    category.defineType("cat3");
    RooDataSet combData(Form("combData_toy%d",itToy),Form("combData_toy%d",itToy),*mass,Index(category),Import("cat0",*genDat[0]),Import("cat1",*genDat[1]),Import("cat2",*genDat[2]),Import("cat3",*genDat[3]));
    
    // ---- construct RooSimultaneous from 4 cats -----
    mIt=0;
    for (int mMC=110; mMC<=150; mMC+=5){
      if (mMC==145) continue;
      RooSimultaneous simPdf(Form("simPdf%d_toy%d",mMC,itToy),Form("simPdf%d_toy%d",mMC,itToy),category);
      simPdf.addPdf(*sigAndBkg[mIt][0],"cat0");
      simPdf.addPdf(*sigAndBkg[mIt][1],"cat1");
      simPdf.addPdf(*sigAndBkg[mIt][2],"cat2");
      simPdf.addPdf(*sigAndBkg[mIt][3],"cat3");
      // ---- fit and save output -----
      RooFitResult *res = simPdf.fitTo(combData,Save(true));//,PrintLevel(-1),Warnings(false),PrintEvalErrors(-1));
      res->floatParsFinal().Print("s");
      muHist[mIt]->Fill(mu->getVal()); 
      // --- find bkg integral -----
      for (int cat=0; cat<nCats; cat++){
        RooAbsReal *intInRange = fitFcn[cat]->createIntegral(*mass,NormSet(*mass),Range("sigWindow"));
        double toBkg = intInRange->getVal()*genDat[cat]->numEntries();
        bias[mIt][cat]->Fill(toData[cat]-toBkg);
      // ---- make plots ----
        if (itToy%toyStep==0){
          TCanvas *c1 = new TCanvas();
          RooPlot *mFrame = mass->frame(Title(Form("3pol fit to data %d toy%d",cat,itToy)));
          genDat[cat]->plotOn(mFrame,DataError(RooDataSet::SumW2));
          fitFcn[cat]->plotOn(mFrame,LineColor(kRed),LineStyle(kDashed));
          sigAndBkg[mIt][cat]->plotOn(mFrame);
          mFrame->Draw();
          c1->Print(Form("plots/fitTogen_m%d_cat%d_toy%d.png",mMC,cat,itToy),"png");
          delete c1;
        }
        
        pol[3][cat]->setVal(pol[0][cat]->getVal());
        pol[4][cat]->setVal(pol[1][cat]->getVal());
        pol[5][cat]->setVal(pol[2][cat]->getVal());
        bkgYield[mIt][cat]->setVal(5000);
        mu->setVal(0);
      }
      mIt++;
    }
  }

  outFile->cd();
  mIt=0;
  for (int mMC=110; mMC<=150; mMC+=5){
    if (mMC==145) continue;
    muHist[mIt]->Write();
    TCanvas *canv = new TCanvas();
    muHist[mIt]->Draw();
    canv->Print(Form("plots/mu_m%d.png",mMC),"png");
    for (int cat=0; cat<nCats; cat++){
      bias[mIt][cat]->Write();
      bias[mIt][cat]->Draw();
      canv->Print(Form("plots/biasCheck_m%d_cat%d.png",mMC,cat),"png");
      canv->Clear();
    }
    mIt++;
  }

    /*
    for (int cat=0; cat<nCats; cat++){
      RooAbsReal *intInRange = fitFcn[cat]->createIntegral(*mass,NormSet(*mass),Range("sigWindow"));
      double toBkg = intInRange->getVal()*genDat[cat]->numEntries();
      bias[cat]->Fill(toData[cat]-toBkg);
      //if (itToy%20==0){
        TCanvas *c1 = new TCanvas();
        RooPlot *mFrame = mass->frame(Title(Form("3pol fit to data %d toy%d",cat,itToy)));
        genDat[cat]->plotOn(mFrame,DataError(RooDataSet::SumW2));
        fitFcn[cat]->plotOn(mFrame,LineColor(kRed),LineStyle(kDashed));
        sigAndBkg[cat]->plotOn(mFrame);
        mFrame->Draw();
        c1->Print(Form("plots/fitTodat_cat%d_toy%d.png",cat,itToy),"png");
        delete c1;
    //  }
      
      pol[3][cat]->setVal(2.);
      pol[4][cat]->setVal(0.25);
      pol[5][cat]->setVal(0.01);
      bkgYield[cat]->setVal(5000);
    }
    mu->setVal(0);
  }
  muHist->Write();
  TCanvas *canv = new TCanvas();
  muHist->Draw();
  canv->Print(Form("plots/mu.png","png"));
  canv->Clear();
  for (int cat=0; cat<nCats; cat++){
    bias[cat]->Write();
    bias[cat]->Draw();
    canv->Print(Form("plots/biasCheck_cat%d.png",cat),"png");
    canv->Clear();
  }
  */
  inFile->Close();
  outFile->Close();

  return 0;
}

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <cassert>
#include <cmath>
#include <algorithm>

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TF1.h"
#include "TLegend.h"
#include "TGaxis.h"
#include "TPad.h"
#include "TFrame.h"
#include "TPaveText.h"

using namespace std;

void AnalyseResults(){

  //gROOT->SetStyle("Plain");
  //gROOT->ForceStyle();
  gROOT->SetBatch();


  system("mkdir -p plots/Graphs");
  system("mkdir -p plots/TexTables");

  const int nMasses=8;
  
  double mean[nMasses];
  double meanError[nMasses];
  double spread[nMasses];
  double spreadError[nMasses];
  double massP[nMasses];
  double massE[nMasses]={0.0};
  double worstBias=0.0;
  double worstSpread=0.0;
  double minSpread;
  
  TGraphErrors *graphMean;
  TGraphErrors *graphSpread;
  TF1 *line = new TF1("line","0.0",108,152);
  line->SetLineColor(kBlack);
  line->SetLineStyle(2);

  ofstream logFile("anRes.log");
  logFile << setw(6) << "mass " << setw(10) << "mean " << setw(10) << "meanE" << setw(10) << "spread" << setw(10) << "spreadE" << endl; 

  TFile *inFile = new TFile("biasCheck.root");
  // ------ do whole range first ------ HistoFile_allCat_ml120_mh120_win20.root
  minSpread=100.;
  int mIt=0;
  for (int mass=110; mass<=150; mass+=5){
    if (mass==145) continue;
    TH1F *temp = (TH1F*)inFile->Get(Form("muHist_m%d",mass));
    TF1 *fitR = new TF1("fitR","gaus",-100,100);
    temp->Fit(fitR,"q");
    mean[mIt] = fitR->GetParameter(1);
    spread[mIt] = fitR->GetParameter(2);
    meanError[mIt] = fitR->GetParError(1);
    spreadError[mIt] = fitR->GetParError(2);
    logFile << setw(6) << mass << setw(10) << mean[mIt] << setw(10) << meanError[mIt] << setw(10) << spread[mIt] << setw(10) << spreadError[mIt] << endl; 
    worstBias = max(worstBias,fabs(mean[mIt]));
    worstSpread = max(worstSpread,fabs(spread[mIt]));
    minSpread = min(minSpread,spread[mIt]);
    massP[mIt]=double(mass);
    mIt++;
  }
  
  graphMean = new TGraphErrors(nMasses,massP,mean,massE,meanError);
  graphSpread = new TGraphErrors(nMasses,massP,spread,massE,spreadError);
  graphMean->SetLineColor(kBlue);
  graphMean->SetMarkerStyle(20);
  graphMean->SetMarkerColor(kBlue);
  graphSpread->SetLineColor(kRed);
  graphSpread->SetMarkerStyle(20);
  graphSpread->SetMarkerColor(kRed);
  graphMean->GetXaxis()->SetTitle("M_{#gamma#gamma}");
  graphMean->GetYaxis()->SetTitle("Amount of SM fitted");
  graphMean->GetYaxis()->SetRangeUser(-0.2,0.2);
  //graphMean->GetYaxis()->SetRangeUser(-1.1*worstBias,1.1*worstBias);
  graphMean->SetTitle("Gen: 3pol. Fit: 3pol. Window: All");
    
  gROOT->Reset();
  TCanvas *c1 = new TCanvas("c1","",200,10,700,500);
  TPad *pad = new TPad("pad","",0,0,1,1);
  pad->SetFillColor(0);
  pad->SetGrid();
  pad->Draw();
  pad->cd();

  // draw a frame to define the range
  TH1F *hr = c1->DrawFrame(108,-0.25,152,0.25);
  hr->SetXTitle("M_{#gamma#gamma}");
  hr->SetYTitle("Amount of SM fitted");
  pad->GetFrame()->SetFillColor(0);
  pad->GetFrame()->SetBorderSize(0);

  // create first graph
  graphMean->Draw("P");
  line->Draw("same");

  // create second graph
  //create a transparent pad drawn on top of the main pad
  c1->cd();
  TPad *overlay = new TPad("overlay","",0,0,1,1);
  overlay->SetFillStyle(4000);
  overlay->SetFillColor(0);
  overlay->SetFrameFillStyle(4000);
  overlay->Draw();
  overlay->cd();
  Double_t xmin = pad->GetUxmin();
  Double_t ymin = 0.5;
  //Double_t ymin = 0.9*minSpread;
  Double_t xmax = pad->GetUxmax();
  Double_t ymax = 1.5;
  //Double_t ymax = 1.1*worstSpread;
  TH1F *hframe = overlay->DrawFrame(xmin,ymin,xmax,ymax);
  hframe->GetXaxis()->SetLabelOffset(99);
  hframe->GetYaxis()->SetLabelOffset(99);
  hframe->GetYaxis()->SetTicks("-");
  hframe->GetYaxis()->SetDrawOption("Y+");
  graphSpread->Draw("P");


  //Draw an axis on the right side
  TGaxis *axis = new TGaxis(xmax,ymin,xmax, ymax,ymin,ymax,510,"+L");
  axis->SetLineColor(kRed);
  axis->SetLabelColor(kRed);
  axis->SetTitle("1#sigma spread from bias");
  axis->SetTitleColor(kRed);
  axis->SetTitleOffset(1.2);
  axis->Draw();

  TPaveText *text = new TPaveText(0.1,0.92,0.9,1.0,"NDC");
  text->SetFillColor(0);
  text->SetLineColor(0);
  text->AddText("Gen: 3pol. Fit: 3pol. Window: All");
  text->Draw("same");
  
  c1->Print("plots/Graphs/3pol.png","png");
  
  logFile.close();
  inFile->Close();
}

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
  const int nFits=12;
  const int nGens=12;
  string genName[nFits] = {"2pol","3pol","4pol","1exp","2exp","3exp","1pow","2pow","3pow","2lau","4lau","6lau"};
  string fitName[nGens] = {"2pol","3pol","4pol","1exp","2exp","3exp","1pow","2pow","3pow","2lau","4lau","6lau"};

  double mean[nGens][nFits][nMasses];
  double meanError[nGens][nFits][nMasses];
  double spread[nGens][nFits][nMasses];
  double spreadError[nGens][nFits][nMasses];
  double massP[nGens][nFits][nMasses];
  double massE[nGens][nFits][nMasses]={{{0.0}}};
  double worstBias[nGens][nFits]={{0.0}};
  double worstSpread[nGens][nFits]={{0.0}};
  double minSpread[nGens][nFits];
  double yLmax[nGens][nFits]={{0.25}};
  double yRmin[nGens][nFits]={{0.4}};
  double yRmax[nGens][nFits]={{1.6}};
  
  TGraphErrors *graphMean[nGens][nFits];
  TGraphErrors *graphSpread[nGens][nFits];
  TF1 *line = new TF1("line","0.0",108,152);
  line->SetLineColor(kBlack);
  line->SetLineStyle(2);

  ofstream logFile("plots/Graphs/anRes.log");
  logFile << "--- RESULTS FROM FITS ---" << endl;

  // --- cycle over gens and fits ---
  for (int f=0; f<nFits; f++){
    if (f>2) continue;              //only looking at pols
    for (int g=0; g<nGens; g++){
      if (g!=f) continue;           // only looking at same-same
      logFile << "--- FIT: " << fitName[f] << " --- GEN: " << genName[g] << " ---" << endl;
      logFile << setw(6) << "mass " << setw(15) << "mean " << setw(15) << "meanE" << setw(15) << "spread" << setw(15) << "spreadE" << endl; 
      TFile *inFile = new TFile(Form("FitResultsFromVols/biasCheck_g%s_f%s.root",genName[g].c_str(),fitName[f].c_str()));
      minSpread[g][f]=100.;
      int mIt=0;
      for (int mass=110; mass<=150; mass+=5){
        if (mass==145) continue;
        TH1F *temp = (TH1F*)inFile->Get(Form("muHist_g%s_f%s_m%d",genName[g].c_str(),fitName[f].c_str(),mass));
        TF1 *fitR = new TF1("fitR","gaus",-100,100);
        temp->Fit(fitR,"q");
        mean[g][f][mIt] = fitR->GetParameter(1);
        spread[g][f][mIt] = fitR->GetParameter(2);
        meanError[g][f][mIt] = fitR->GetParError(1);
        spreadError[g][f][mIt] = fitR->GetParError(2);
        logFile << setw(6) << mass << setw(15) << mean[g][f][mIt] << setw(15) << meanError[g][f][mIt] << setw(15) << spread[g][f][mIt] << setw(15) << spreadError[g][f][mIt] << endl; 
        worstBias[g][f] = max(worstBias[g][f],fabs(mean[g][f][mIt])+meanError[g][f][mIt]);
        worstSpread[g][f] = max(worstSpread[g][f],fabs(spread[g][f][mIt]));
        minSpread[g][f] = min(minSpread[g][f],spread[g][f][mIt]);
        massP[g][f][mIt]=double(mass);
        mIt++;
      }
     
      // find nearest worst bias to 0.05
      yLmax[g][f] = ceil(worstBias[g][f]*20.0)/20.0;
      yRmin[g][f] = 1.0-yLmax[g][f]*4.0;
      yRmax[g][f] = 1.0+yLmax[g][f]*4.0;

      graphMean[g][f] = new TGraphErrors(nMasses,massP[g][f],mean[g][f],massE[g][f],meanError[g][f]);
      graphSpread[g][f] = new TGraphErrors(nMasses,massP[g][f],spread[g][f],massE[g][f],spreadError[g][f]);
      graphMean[g][f]->SetLineColor(kBlue);
      graphMean[g][f]->SetMarkerStyle(20);
      graphMean[g][f]->SetMarkerColor(kBlue);
      graphSpread[g][f]->SetLineColor(kRed);
      graphSpread[g][f]->SetMarkerStyle(20);
      graphSpread[g][f]->SetMarkerColor(kRed);
      graphMean[g][f]->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV)");
      graphMean[g][f]->GetYaxis()->SetTitle("Amount of SM fitted");
      graphMean[g][f]->GetYaxis()->SetRangeUser(-0.2,0.2);
      //graphMean->GetYaxis()->SetRangeUser(-1.1*worstBias,1.1*worstBias);
      graphMean[g][f]->SetTitle(Form("Gen: %s. Fit: %s. Window: All",genName[g].c_str(),fitName[f].c_str()));
        
      gROOT->Reset();
      TCanvas *c1 = new TCanvas("c1","",200,10,700,500);
      TPad *pad = new TPad("pad","",0,0,1,1);
      pad->SetFillColor(0);
      pad->SetGrid();
      pad->Draw();
      pad->cd();

      // draw a frame to define the range
      TH1F *hr = c1->DrawFrame(108,-1*yLmax[g][f],152,yLmax[g][f]);
      hr->SetXTitle("m_{#gamma#gamma} (GeV)");
      hr->SetYTitle("Amount of SM fitted");
      pad->GetFrame()->SetFillColor(0);
      pad->GetFrame()->SetBorderSize(0);

      // create first graph
      graphMean[g][f]->Draw("P");
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
      Double_t ymin = yRmin[g][f];
      //Double_t ymin = 0.9*minSpread;
      Double_t xmax = pad->GetUxmax();
      Double_t ymax = yRmax[g][f];
      //Double_t ymax = 1.1*worstSpread;
      TH1F *hframe = overlay->DrawFrame(xmin,ymin,xmax,ymax);
      hframe->GetXaxis()->SetLabelOffset(99);
      hframe->GetYaxis()->SetLabelOffset(99);
      hframe->GetYaxis()->SetTicks("-");
      hframe->GetYaxis()->SetDrawOption("Y+");
      graphSpread[g][f]->Draw("P");


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
      text->AddText(Form("Gen: %s. Fit: %s. Window: All. B_{max}=%1.2f. s_{max}=%1.2f",genName[g].c_str(),fitName[f].c_str(),worstBias[g][f],worstSpread[g][f]));
      text->Draw("same");
      
      c1->Print(Form("plots/Graphs/g%s_f%s.png",genName[g].c_str(),fitName[f].c_str()),"png");
      inFile->Close();
    }
  }
  system("python publishToWeb.py");
  system("cp -r plots/Graphs ~/public_html/h2g/BackgroundSystematic/");
  logFile.close();
}

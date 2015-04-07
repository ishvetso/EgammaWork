#include <TFile.h>
#include <TTree.h>
#include <TTreeFormula.h>
#include <TH1.h>
#include <TH2.h>
#include <THStack.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TGraph.h>
#include <iostream>
#include <iostream>
#include <sstream>
#include <map>
#include <algorithm> 
#include </afs/cern.ch/work/i/ishvetso/GitHub/IvanShvetsov/CMS_stylistics/tdrstyle.C>
#include </afs/cern.ch/work/i/ishvetso/GitHub/IvanShvetsov/CMS_stylistics/CMS_lumi.cpp>

using namespace std;

int count(string cut, TTree *tree)
{
  int N = tree -> GetEntries(cut.c_str() );
  return N;
  
}

TGraph * graph_ROC(string var)
{
  TFile file_bkg("/afs/cern.ch/work/i/ishvetso/EgammaWork/ElectronIsolationPUPPI/CMSSW_7_3_0/src/EgammaWork/ElectronNtupler/test/crab_projects/crab_Electron-Isolation_CITK_validation_ttbar_miniAOD_PUPPI/results/ttbar.root");
  TFile file_sig("/afs/cern.ch/work/i/ishvetso/EgammaWork/ElectronIsolationPUPPI/CMSSW_7_3_0/src/EgammaWork/ElectronNtupler/test/crab_projects/crab_Electron-Isolation_CITK_validation_DY_miniAOD_PUPPI/results/DY.root");
  
  TTree *tree_bkg = (TTree*)file_bkg.Get("ntupler/ElectronTree");
  TTree *tree_sig = (TTree*)file_sig.Get("ntupler/ElectronTree");
  
  string SigSelection = "pt > 20 && isTrueElectron == 1";
  string BkgSelection = "pt > 20 && ( isTrueElectron == 0 || isTrueElectron == 3 )";

   ostringstream ss;
   ss.precision(3);
   
   vector <float> sigEff, bkgEff;
   
  for (float iCut = 2.; iCut > 0.01; iCut -= 0.02)
  {
    ss << iCut;
    string myCut = ss.str();
    float effSig = (float) count(SigSelection + " && "+ var + " < " + myCut, tree_sig)/(count(SigSelection, tree_sig));
    float effBkg = (float) count(BkgSelection + " && "+ var + " < " + myCut, tree_bkg)/(count(BkgSelection, tree_bkg));
    sigEff.push_back(effSig);
    bkgEff.push_back(effBkg);
    ss.str("");
  }
  
  TGraph *gr = new TGraph( sigEff.size(), bkgEff.data(), sigEff.data());
  
  return gr;
   
}

void draw_ROC()
{
  gStyle->SetCanvasColor(kWhite);
//gStyle->SetTitleFillColor(kWhite);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
//gPad->SetGrid();
  
  TGraph * gr1 =  graph_ROC("relIsoWithEA");
  TGraph * gr2 =  graph_ROC("reliso_PUPPI");
  gr1 -> GetXaxis() -> SetTitle("bkg eff");
  gr1 -> GetYaxis() -> SetTitle("sig eff");
  gr1 -> SetLineWidth(2.);
  gr2 -> SetLineWidth(2.);
  
  gr1 -> SetLineColor(kBlue);
  gr2 -> SetLineColor(kRed);
  
  setTDRStyle();   
  TCanvas *c1= new TCanvas("c1","canvas",1200,800);

  TLegend *leg = new TLegend(0.6,0.2,0.8,0.5);
  leg -> SetFillColor(kWhite);  
  
  leg->AddEntry(gr1, "relIsoWithEA","l");
  leg->AddEntry(gr2, "reliso_PUPPI","l");
  
 
  c1 -> cd();
  gr1 -> Draw();
  gr2 -> Draw("SAME");
  leg -> Draw("SAME");
  CMS_lumi( c1, 4, 0 );
  c1 -> SaveAs("ROC.png");
}
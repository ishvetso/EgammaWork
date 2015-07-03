#include "ROC_drawer.hpp"
 ROC_Drawer::ROC_Drawer(){}
    
  TGraph * ROC_Drawer::graph_ROC(string var, int Nbins, float xmin, float xmax)
{
  TFile file_bkg("/afs/cern.ch/work/i/ishvetso/EgammaWork/ElectronIsolationPUPPI_745_30June/CMSSW_7_4_5/src/EgammaWork/ElectronNtupler/test/crab_projects/crab_Electron-Isolation_CITK_validation_ttbar_miniAOD_PUPPI_with_NoLeptons_update30June2015/results/ttbar.root");
  TFile file_sig("/afs/cern.ch/work/i/ishvetso/EgammaWork/ElectronIsolationPUPPI_745_30June/CMSSW_7_4_5/src/EgammaWork/ElectronNtupler/test/crab_projects/crab_Electron-Isolation_CITK_validation_DY_miniAOD_PUPPI_with_NoLeptons_update30June2015/results/DY.root");
  
  TTree *tree_bkg = (TTree*)file_bkg.Get("ntupler/ElectronTree");
  TTree *tree_sig = (TTree*)file_sig.Get("ntupler/ElectronTree");
  
  string SigSelection = "pt > 20 && isTrueElectron == 1";
  string BkgSelection = "pt > 20 && ( isTrueElectron == 0 || isTrueElectron == 3 )";

   ostringstream ss;
   ss.precision(3);
   
   vector <float> sigEff, bkgEff;
   
   TH1F *hist_sig = new TH1F("sig", "sig", Nbins, xmin, xmax);
   TH1F *hist_bkg = new TH1F("bkg", "bkg", Nbins, xmin, xmax);
   
   tree_sig -> Project("sig", var.c_str(), (SigSelection + "&&" + addSelection) .c_str());
   tree_bkg -> Project("bkg", var.c_str(), (BkgSelection + "&&" + addSelection) .c_str());
   
   
  for (int iBin = 2; iBin < Nbins ; iBin += 1)
  {
    float effSig = (float) (hist_sig -> Integral(1, iBin))/(hist_sig -> Integral());
    float effBkg = (float) (hist_bkg -> Integral(1, iBin))/(hist_bkg -> Integral());
    
    sigEff.push_back(effSig);
    bkgEff.push_back(1.0 - effBkg);
  }
  
  TGraph *gr = new TGraph( sigEff.size(), sigEff.data(), bkgEff.data());
  
  return gr;   
}

void ROC_Drawer::draw_ROC()
{
  gStyle->SetCanvasColor(kWhite);
//gStyle->SetTitleFillColor(kWhite);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
//gPad->SetGrid();
  
  TGraph * gr1 =  graph_ROC("relIsoWithEA", 10000, 0., 10.);
  TGraph * gr2 =  graph_ROC("relIsoWithDBeta", 10000, 0., 10.);
  TGraph * gr3 =  graph_ROC("reliso_PUPPI", 10000, 0., 10.);
  TGraph * gr4 =  graph_ROC("reliso_PUPPI_NoLeptons", 10000, 0., 10.);
  
  
  gr1 -> GetYaxis() -> SetTitle("bkg rej");
  gr1 -> GetXaxis() -> SetTitle("sig eff");
  
  gr1 -> SetLineWidth(2.);
  gr2 -> SetLineWidth(2.);
  gr3 -> SetLineWidth(2.);
  gr4 -> SetLineWidth(2.);
  
  gr1 -> SetLineColor(kBlue);
  gr2 -> SetLineColor(kGreen+2);
  gr3 -> SetLineColor(kRed);
  gr4 -> SetLineColor(kYellow + 1);
  
  setTDRStyle();   
  TCanvas *c1= new TCanvas("c1","canvas",1200,800);

  TLegend *leg = new TLegend(0.2,0.3,0.5,0.7);
  leg -> SetFillColor(kWhite);  
  
  leg->AddEntry(gr1, "relIsoWithEA","l");
  leg->AddEntry(gr2, "relIsoWithDBeta","l");
  leg->AddEntry(gr3, "reliso_PUPPI","l");
  leg->AddEntry(gr4, "reliso_PUPPI_NoLeptons","l");
  
 
  c1 -> cd();
  gr1 -> Draw();
  gr2 -> Draw("SAME");
  gr3 -> Draw("SAME");
  gr4 -> Draw("SAME");
  
  leg -> Draw("SAME");
  CMS_lumi( c1, 4, 0 );
  c1 -> SaveAs( (name +  ".png").c_str());
  delete c1;
}

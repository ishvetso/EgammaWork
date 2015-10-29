#include "ROC_drawer.hpp"
 ROC_Drawer::ROC_Drawer(){}
    
  TGraph * ROC_Drawer::graph_ROC(string var, int Nbins, float xmin, float xmax)
{
  TFile file(FileName.c_str());
  
  TTree *tree = (TTree*)file.Get("ntupler/PhotonTree");
 

   ostringstream ss;
   ss.precision(3);
   
   vector <float> sigEff, bkgEff;
   
   TH1F *hist_sig = new TH1F("sig", "sig", Nbins, xmin, xmax);
   TH1F *hist_bkg = new TH1F("bkg", "bkg", Nbins, xmin, xmax);

   TH1F *hist_sig_abs = new TH1F("sig_abs", "sig_abs", Nbins, 0., 10.);
   TH1F *hist_bkg_abs = new TH1F("bkg_abs", "bkg_abs", Nbins, 0., 10.);
   
   tree -> Project("sig", var.c_str(), (SigSelection + "&&" + addSelection) .c_str());
   tree -> Project("bkg", var.c_str(), (BkgSelection + "&&" + addSelection) .c_str());
   
   tree -> Project("sig_abs", var.c_str(), (SigSelection + "&&" + addSelection) .c_str());
   tree -> Project("bkg_abs", var.c_str(), (BkgSelection + "&&" + addSelection) .c_str());
   
   
  for (int iBin = 1; iBin <= Nbins ; iBin += 1)
  {
    float effSig = (float) (hist_sig -> Integral(1, iBin))/(hist_sig_abs -> Integral());
    float effBkg = (float) (hist_bkg -> Integral(1, iBin))/(hist_bkg_abs -> Integral());
    
    sigEff.push_back(effSig);
    bkgEff.push_back(1.0 - effBkg);
  }
  
  if (sigEff.size() != bkgEff.size()) cerr << "Oops..."<< endl;

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
  

  TGraph * gr1 =  graph_ROC("relisoWithEA_CITK_", 10000, 0., 10.);
  TGraph * gr2 =  graph_ROC("relisoWithEA_PUPPI", 10000, 0., 10.);
  TGraph * gr3 =  graph_ROC("relisoWithEA_pf", 10000, 0., 10.);
  
  
  gr1 -> GetYaxis() -> SetTitle("bkg rejection");
  gr1 -> GetXaxis() -> SetTitle("sig eff");
  
  gr1 -> SetLineWidth(2.);
  gr2 -> SetLineWidth(2.);
  gr3 -> SetLineWidth(2.);
  
  gr1 -> SetLineColor(kBlue);
  gr2 -> SetLineColor(kGreen);
  gr3 -> SetLineColor(kRed);
  
  setTDRStyle();   
  TCanvas *c1= new TCanvas("c1","canvas",1200,800);

  TLegend *leg = new TLegend(0.2,0.3,0.5,0.7);
  leg -> SetFillColor(kWhite);  
  
 
  leg->AddEntry(gr1, "relisoWithEA_CITK_","l");
  leg->AddEntry(gr2, "relisoWithEA_PUPPI","l");
  leg->AddEntry(gr3, "relisoWithEA_pf","l");
  
 
  c1 -> cd();
  gr1 -> Draw();
  gr2 -> Draw("SAME");
  gr3 -> Draw("SAME");
  
  
  leg -> Draw("SAME");
  CMS_lumi( c1, 4, 0 );
  c1 -> SaveAs( (name +  ".png").c_str());
  delete c1;
}

#include "ROC_drawer.hpp"
 ROC_Drawer::ROC_Drawer(){}
    
  TGraph * ROC_Drawer::graph_ROC(string var, int Nbins, float xmin, float xmax)

{
  std::cout << var << std::endl;
  TFile file_bkg(bkgFileName.c_str());
  TFile file_sig(sigFileName.c_str());
  
  TTree *tree_bkg = (TTree*)file_bkg.Get("ntupler/ElectronTree");
  TTree *tree_sig = (TTree*)file_sig.Get("ntupler/ElectronTree");

   ostringstream ss;
   ss.precision(3);
   
   vector <float> sigEff, bkgEff;
   
   TH1F *hist_sig = new TH1F("sig", "sig", Nbins, xmin, xmax);
   TH1F *hist_bkg = new TH1F("bkg", "bkg", Nbins, xmin, xmax);

   TH1F *hist_sig_abs = new TH1F("sig_abs", "sig_abs", Nbins, 0., 1000.);
   TH1F *hist_bkg_abs = new TH1F("bkg_abs", "bkg_abs", Nbins, 0., 1000.);

   hist_sig -> Sumw2();
   hist_bkg -> Sumw2();
   hist_sig_abs -> Sumw2();
   hist_bkg_abs -> Sumw2();
   
   tree_sig -> Project("sig", var.c_str(), (SigSelection + " && " + addSelection) .c_str());
   tree_bkg -> Project("bkg", var.c_str(), (BkgSelection + " && " + addSelection) .c_str());
   
   tree_sig -> Project("sig_abs", var.c_str(), (SigSelection + " && " + addSelection) .c_str());
   tree_bkg -> Project("bkg_abs", var.c_str(), (BkgSelection + " && " + addSelection) .c_str());
   
   
  for (int iBin = 1; iBin < Nbins ; iBin += 1)
  {
    float effSig = (float) (hist_sig -> Integral(1, iBin))/(hist_sig_abs -> Integral());
    float effBkg = (float) (hist_bkg -> Integral(1, iBin))/(hist_bkg_abs -> Integral());
    std::cout << effSig << "  " << effBkg << std::endl;
    
    if (effSig < 0.99 && effSig >  0.8)
    {
      sigEff.push_back(effSig);
      bkgEff.push_back(1.0 - effBkg);
    }
  }
  
  TGraph *gr = new TGraph( sigEff.size(), sigEff.data(), bkgEff.data());
  
  return gr;   
}

void ROC_Drawer::draw_ROC()
{
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  
  TGraph * gr1 =  graph_ROC("reliso_ConeVeto_EA", 10000, 0., 15.);
  TGraph * gr2 =  graph_ROC("reliso_MapBasedVeto_EA", 10000, 0., 15.);
  TGraph * gr3 =  graph_ROC("relIsoWithDBeta", 10000, 0., 15.);
  TGraph * gr4 =  graph_ROC("reliso_PUPPI_ConeVeto_average", 10000, 0., 15.);
  TGraph * gr5 =  graph_ROC("reliso_PUPPI_MapBasedVeto_average", 10000, 0., 15.);
  TGraph * gr6 =  graph_ROC("reliso_ConeVeto_raw", 10000, 0., 15.);
  TGraph * gr7 =  graph_ROC("reliso_MapBasedVeto_raw", 10000, 0., 15.);

  gr1 -> GetYaxis() -> SetTitle("bkg rejection");
  gr1 -> GetXaxis() -> SetTitle("sig eff");
  
  gr1 -> SetLineColor(kGreen);
  gr2 -> SetLineColor(kGreen);
  gr2 -> SetLineStyle(2);
  gr3 -> SetLineColor(kOrange);
  gr4 -> SetLineColor(kRed);
  gr5 -> SetLineColor(kRed);
  gr5 -> SetLineStyle(2);
  gr6 -> SetLineColor(kCyan);
  gr7 -> SetLineColor(kCyan);
  gr7 -> SetLineStyle(2);

  gr1 -> SetLineWidth(2.);
  gr2 -> SetLineWidth(2.);
  gr3 -> SetLineWidth(2.);
  gr4 -> SetLineWidth(2.);
  gr5 -> SetLineWidth(2.);
  gr6 -> SetLineWidth(2.);
  gr7 -> SetLineWidth(2.);
 
  
  setTDRStyle();   
  TCanvas *c1= new TCanvas("c1","canvas",1200,800);

  TLegend *leg = new TLegend(0.2,0.3,0.5,0.7);
  leg -> SetFillColor(kWhite);  
  
  leg->AddEntry(gr1, "effective area, cone veto","l");
  leg->AddEntry(gr2, "effective area, map based veto","l");
  leg->AddEntry(gr3, "#delta#beta-corrected","l");
  leg->AddEntry(gr4, "PUPPI cone veto combined","l");
  leg->AddEntry(gr5, "PUPPI map based veto combined","l");
  leg->AddEntry(gr6, "raw cone veto","l");
  leg->AddEntry(gr7, "raw map based veto","l");

  c1 -> cd();
  gr1 -> Draw();
  gr2 -> Draw("SAME");
  gr3 -> Draw("SAME");
  gr4 -> Draw("SAME");
  gr5 -> Draw("SAME");
  gr6 -> Draw("SAME");
  gr7 -> Draw("SAME");
  leg -> Draw("SAME");

  CMS_lumi( c1, 4, 0 );

  c1 -> SaveAs( (name +  ".pdf").c_str());
  delete c1;
}

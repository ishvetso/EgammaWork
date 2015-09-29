#include "ROC_drawer.hpp"
 ROC_Drawer::ROC_Drawer(){}
    
  TGraph * ROC_Drawer::graph_ROC(string var, int Nbins, float xmin, float xmax)
{
  TFile file_bkg(bkgFileName.c_str());
  TFile file_sig(sigFileName.c_str());
  
  TTree *tree_bkg = (TTree*)file_bkg.Get("ntupler/ElectronTree");
  TTree *tree_sig = (TTree*)file_sig.Get("ntupler/ElectronTree");

   ostringstream ss;
   ss.precision(3);
   
   vector <float> sigEff, bkgEff;
   
   TH1F *hist_sig = new TH1F("sig", "sig", Nbins, xmin, xmax);
   TH1F *hist_bkg = new TH1F("bkg", "bkg", Nbins, xmin, xmax);

   TH1F *hist_sig_abs = new TH1F("sig_abs", "sig_abs", Nbins, 0., 15.);
   TH1F *hist_bkg_abs = new TH1F("bkg_abs", "bkg_abs", Nbins, 0., 15.);

   hist_sig -> Sumw2();
   hist_bkg -> Sumw2();
   hist_sig_abs -> Sumw2();
   hist_bkg_abs -> Sumw2();
   
   tree_sig -> Project("sig", var.c_str(), (SigSelection + "&&" + addSelection) .c_str());
   tree_bkg -> Project("bkg", var.c_str(), (BkgSelection + "&&" + addSelection) .c_str());
   
   tree_sig -> Project("sig_abs", var.c_str(), (SigSelection + "&&" + addSelection) .c_str());
   tree_bkg -> Project("bkg_abs", var.c_str(), (BkgSelection + "&&" + addSelection) .c_str());
   
   
  for (int iBin = 2; iBin < Nbins ; iBin += 1)
  {
    float effSig = (float) (hist_sig -> Integral(1, iBin))/(hist_sig_abs -> Integral());
    float effBkg = (float) (hist_bkg -> Integral(1, iBin))/(hist_bkg_abs -> Integral());
    
    if (effSig < 0.99 && effSig >  0.8)
    {
      sigEff.push_back(effSig);
      bkgEff.push_back(1.0 - effBkg);
    }

  /*  if (effSig < 0.52 && effSig > 0.48)
      { 
        std::cout << " var "  << var << std::endl;
        std::cout << "signal effeciency " << effSig << std::endl;
        std::cout << "cut " << hist_sig -> GetBinLowEdge(iBin + 1) << std::endl;
      }*/

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
  
  TGraph * gr1 =  graph_ROC("relIsoWithEA", 10000, 0., 15.);
  TGraph * gr2 =  graph_ROC("relIsoWithDBeta", 10000, 0., 15.);
  TGraph * gr3 =  graph_ROC("reliso_PUPPI", 10000, 0., 15.);
  TGraph * gr4 =  graph_ROC("reliso_PUPPI_NoLeptons", 10000, 0., 15.);
  TGraph * gr5 =  graph_ROC("reliso_PUPPI_average", 10000, 0., 15.);
  
  
  gr1 -> GetYaxis() -> SetTitle("bkg rejection");
  gr1 -> GetXaxis() -> SetTitle("sig eff");
  
  gr1 -> SetLineWidth(2.);
  gr2 -> SetLineWidth(2.);
  gr3 -> SetLineWidth(2.);
  gr4 -> SetLineWidth(2.);
  gr5 -> SetLineWidth(2.);
  
  gr1 -> SetLineColor(kCyan);
  gr2 -> SetLineColor(kBlack);
  gr3 -> SetLineColor(kBlue);
  gr4 -> SetLineColor(kRed);
  gr5 -> SetLineColor(kGreen);
 
  
  setTDRStyle();   
  TCanvas *c1= new TCanvas("c1","canvas",1200,800);

  TLegend *leg = new TLegend(0.2,0.3,0.5,0.7);
  leg -> SetFillColor(kWhite);  
  
  leg->AddEntry(gr1, "effective area","l");
  leg->AddEntry(gr2, "#delta#beta-corrected","l");
  leg->AddEntry(gr3, "PUPPI","l");
  leg->AddEntry(gr4, "PUPPINoLeptons","l");
  leg->AddEntry(gr5, "PUPPI combined","l");
 
 
  c1 -> cd();
  gr1 -> Draw();
  gr2 -> Draw("SAME");
  gr3 -> Draw("SAME");
  gr4 -> Draw("SAME");
  gr5 -> Draw("SAME");
  leg -> Draw("SAME");

  CMS_lumi( c1, 4, 0 );

  c1 -> SaveAs( (name +  ".png").c_str());
  delete c1;
}

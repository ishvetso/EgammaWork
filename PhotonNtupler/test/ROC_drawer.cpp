#include "ROC_drawer.hpp"
 ROC_Drawer::ROC_Drawer(){
 }
    
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

TGraph * ROC_Drawer::addPoint(std::string nominator, std::string denominator, std::string _SigSelection, std::string _BkgSelection){

  TFile file(FileName.c_str());
  TTree *tree = (TTree*)file.Get("ntupler/PhotonTree");

  TH1F *hist_nom_sig = new TH1F("nom_sig", "nom_sig", 1, -10., 10.);
  TH1F *hist_denom_sig = new TH1F("denom_sig", "denom_sig", 1, -10., 10.);

  TH1F *hist_nom_bkg = new TH1F("nom_bkg", "nom_bkg", 1, -10., 10.);
  TH1F *hist_denom_bkg = new TH1F("denom_bkg", "denom_bkg", 1, -10., 10.);

  tree -> Project("nom_sig", "phi", (nominator + " && " + _SigSelection  + " && " + addSelection).c_str());
  tree -> Project("denom_sig", "phi", (denominator + " && "+ _SigSelection  + " && " + addSelection).c_str());

  tree -> Project("nom_bkg", "phi", (nominator + " && " + _BkgSelection + " && " + addSelection ).c_str());
  tree -> Project("denom_bkg", "phi", (denominator + " && "+ _BkgSelection + " && " + addSelection ).c_str());

  vector <float> sigEff, bkgRej;
  sigEff.push_back(hist_nom_sig -> Integral() / hist_denom_sig -> Integral());
  bkgRej.push_back(1. - (hist_nom_bkg -> Integral() / hist_denom_bkg -> Integral()));

  TGraph *gr = new TGraph( sigEff.size(), sigEff.data(),  bkgRej.data());

  std::cout << sigEff[0] << "  " << bkgRej[0] << std::endl;

  return gr;
}

void ROC_Drawer::draw_ROC()
{
  gStyle->SetCanvasColor(kWhite);
//gStyle->SetTitleFillColor(kWhite);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
//gPad->SetGrid();
  

  TGraph * gr1 =  graph_ROC("relisoWithEA_CITK", 10000, 0., 10.);
  TGraph * gr2 =  graph_ROC("relisoWithEA_PUPPI", 10000, 0., 10.);
  TGraph * gr3 =  graph_ROC("relisoWithEA_pf", 10000, 0., 10.);

  std::string HoverEAndsigmaIetaIetaCutsTight = "( (isEB == 1 && hOverE < 0.05 && full5x5_sigmaIetaIeta < 0.01) || (isEB == 0 && hOverE < 0.05 && full5x5_sigmaIetaIeta < 0.0268)  )";
  std::string TightCutsBarrel = "isoChargedHadrons_CITK < 0.76 && isoNeutralHadronsWithEA_CITK < (0.97 + 0.014*pt + 0.000019*pt*pt) && isoPhotonsWithEA_CITK < ( 0.08 + 0.0053*pt )";
  std::string TightCutsEndcap = "isoChargedHadrons_CITK < 0.56 && isoNeutralHadronsWithEA_CITK < (2.09 +0.0139*pt+0.000025*pt*pt) && isoPhotonsWithEA_CITK < ( 0.16 + 0.0034*pt )";
  std::string TightCuts = "((isEB == 1 && " + TightCutsBarrel + ") || (isEB == 0 && " + TightCutsEndcap + ") )";
  
  std::string HoverEAndsigmaIetaIetaCutsMedium = "( (isEB == 1 && hOverE < 0.05 && full5x5_sigmaIetaIeta < 0.0102) || (isEB == 0 && hOverE < 0.05 && full5x5_sigmaIetaIeta < 0.0268)  )";
  std::string MediumCutsBarrel = "isoChargedHadrons_CITK < 1.37 && isoNeutralHadronsWithEA_CITK < (1.06 + 0.014*pt + 0.000019*pt*pt) && isoPhotonsWithEA_CITK < ( 0.28 + 0.0053*pt)";
  std::string MediumCutsEndcap = "isoChargedHadrons_CITK < 1.10 && isoNeutralHadronsWithEA_CITK < (2.69 + 0.0139*pt+0.000025*pt*pt ) && isoPhotonsWithEA_CITK < ( 0.39 + 0.0034*pt )";
  std::string MediumCuts = "((isEB == 1 && " + MediumCutsBarrel + ") || (isEB == 0 && " + MediumCutsEndcap + ") )";

  std::string HoverEAndsigmaIetaIetaCutsLoose = "( (isEB == 1 && hOverE < 0.05 && full5x5_sigmaIetaIeta < 0.0102) || (isEB == 0 && hOverE < 0.05 && full5x5_sigmaIetaIeta < 0.0274)  )";
  std::string LooseCutsBarrel = "isoChargedHadrons_CITK < 3.32 && isoNeutralHadronsWithEA_CITK < (1.92 + 0.014*pt + 0.000019*pt*pt ) && isoPhotonsWithEA_CITK < ( 0.81 + 0.0053*pt)";
  std::string LooseCutsEndcap = "isoChargedHadrons_CITK < 1.97 && isoNeutralHadronsWithEA_CITK < (11.86 + 0.0139*pt+0.000025*pt*pt ) && isoPhotonsWithEA_CITK < ( 0.83 + 0.0034*pt )";
  std::string LooseCuts = "((isEB == 1 && " + LooseCutsBarrel + ") || (isEB == 0 && " + LooseCutsEndcap + ") )";
  
  TGraph *gr_Tight = addPoint(HoverEAndsigmaIetaIetaCutsTight + " && " + TightCuts, HoverEAndsigmaIetaIetaCutsTight  , "( pt > 20 && isTrue == 1 )", "( pt > 20 && ( isTrue != 1  ) )");
  TGraph *gr_Medium = addPoint(HoverEAndsigmaIetaIetaCutsMedium + " && " + MediumCuts, HoverEAndsigmaIetaIetaCutsMedium  , "( pt > 20 && isTrue == 1 )", "( pt > 20 && ( isTrue != 1  ) )");
  TGraph *gr_Loose = addPoint(HoverEAndsigmaIetaIetaCutsLoose + " && " + LooseCuts, HoverEAndsigmaIetaIetaCutsLoose  , "( pt > 20 && isTrue == 1 )", "( pt > 20 && ( isTrue != 1  ) )");

  
  gr1 -> GetYaxis() -> SetTitle("bkg rejection");
  gr1 -> GetXaxis() -> SetTitle("sig eff");

  gr1 -> GetYaxis() -> SetRangeUser(0., 1.1);
  gr1 -> GetXaxis() -> SetRangeUser(0.5, 1.1);

  
  gr1 -> SetLineWidth(2.);
  gr2 -> SetLineWidth(2.);
  gr3 -> SetLineWidth(2.);
  gr_Tight -> SetLineWidth(5.);
  
  gr1 -> SetLineColor(kBlue);
  gr2 -> SetLineColor(kGreen);
  gr3 -> SetLineColor(kRed);
  gr_Tight -> SetMarkerColor(kBlack);
  gr_Tight ->SetMarkerSize(3.5);
  gr_Medium -> SetMarkerColor(kRed);
  gr_Medium ->SetMarkerSize(3.5);
  gr_Loose -> SetMarkerColor(kBlue);
  gr_Loose ->SetMarkerSize(3.5);
  
  setTDRStyle();   
  TCanvas *c1= new TCanvas("c1","canvas",1200,800);

  TLegend *leg = new TLegend(0.2,0.65,0.5,0.95);
  leg -> SetFillColor(kWhite);  
  leg->SetFillStyle(0);
  leg -> SetBorderSize(0);
  
 
  leg->AddEntry(gr1, "effective area","l");
  leg->AddEntry(gr2, "PUPPI","l");
  leg->AddEntry(gr3, "effective area, pfIsoVariables","l");
  leg-> AddEntry(gr_Tight, "Tight WP", "p");
  leg-> AddEntry(gr_Medium, "Medium WP", "p");
  leg-> AddEntry(gr_Loose, "Loose WP", "p");
  
 
  c1 -> cd();
  gr1 -> Draw();
  gr2 -> Draw("SAME");
  gr3 -> Draw("SAME");
  gr_Tight -> Draw("P*SAME");
  gr_Medium -> Draw("P*SAME");
  gr_Loose -> Draw("P*SAME");
  
  
  leg -> Draw("SAME");
  CMS_lumi( c1, 4, 0 );
  c1 -> SaveAs( (name +  ".png").c_str());
  delete c1;
}

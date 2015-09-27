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

void draw_var(string var, double xmin, double xmax)
{
	TFile f_miniAOD("/afs/cern.ch/work/i/ishvetso/EgammaWork/PhotonIsolation_final_checks/CMSSW_7_4_7/src/EgammaWork/PhotonNtupler/test/crab_projects/crab_PhotonIsolation_checks_17September2015_miniAOD_fixing_keys_in_love_with_crab/results/GJets.root");
	TFile f_AOD("/afs/cern.ch/work/i/ishvetso/EgammaWork/PhotonIsolation_final_checks/CMSSW_7_4_7/src/EgammaWork/PhotonNtupler/test/crab_projects/crab_PhotonIsolation_checks_17September2015_AOD_fixing_keys_in_love_with_crab/results/tree_AOD.root");

	TTree *tree_miniAOD = (TTree*)f_miniAOD.Get("ntupler/PhotonTree");
	TTree *tree_AOD = (TTree*)f_AOD.Get("ntupler/PhotonTree");

	TH1F *hist_miniAOD = new TH1F("miniAOD", "miniAOD", 40, xmin, xmax);
	TH1F *hist_AOD = new TH1F("AOD", "AOD", 40, xmin, xmax);
  hist_AOD -> Sumw2();
  hist_miniAOD -> Sumw2();
	
	vector <float> * x_AOD = new std::vector<float>();
  vector <float> * x_miniAOD = new std::vector<float>();
  
  vector <float> * hOverE = new std::vector<float>();
  

  tree_miniAOD -> SetBranchAddress(var.c_str(), &x_miniAOD);
  tree_AOD -> SetBranchAddress(var.c_str(), &x_AOD);
  tree_AOD -> SetBranchAddress("hOverE", &hOverE);
  tree_miniAOD -> SetBranchAddress("hOverE", &hOverE);

  for (unsigned int iEntry = 0; iEntry < tree_miniAOD -> GetEntries(); iEntry ++)
  {
    tree_miniAOD -> GetEntry(iEntry);
    for (unsigned int i = 0; i < x_miniAOD -> size(); i++)
    {
      if ( hOverE -> at(i) < 0.15) hist_miniAOD -> Fill(x_miniAOD -> at(i));
      //std::cout << x_miniAOD -> at(i) << std::endl;
    }
    x_miniAOD -> clear();
    hOverE -> clear();
  }

   for (unsigned int iEntry = 0; iEntry < tree_miniAOD -> GetEntries(); iEntry ++)
   {
      tree_AOD -> GetEntry(iEntry);
      for (unsigned int i = 0; i < x_AOD -> size(); i++)
      {
        if ( hOverE -> at(i) < 0.15) hist_AOD -> Fill(x_AOD -> at(i));
        //std::cout << x_AOD -> at(i) << std::endl;
        
      }
       x_AOD -> clear();
        hOverE -> clear();
   }




    setTDRStyle();   
  	TCanvas *c1= new TCanvas("c1","canvas",1200,800);

  	hist_AOD -> SetLineColor(kRed);
  	hist_miniAOD -> SetLineColor(kBlue);
  	hist_AOD -> SetLineWidth(2.);
  	hist_miniAOD -> SetLineWidth(2.);

 	TLegend *leg = new TLegend(0.2,0.6,0.5,0.8);
 	leg -> SetFillStyle(0);
  	leg -> SetFillColor(kWhite);  
  	leg -> SetHeader(var.c_str());
  	leg -> SetBorderSize(0);
  
 	leg->AddEntry(hist_AOD, "AOD/miniAOD","l");
  	//leg->AddEntry(hist_miniAOD, "miniAOD","l");

  	hist_AOD -> Divide(hist_miniAOD);
  	hist_AOD -> Draw("E1");
  	//hist_AOD -> Draw("SAME");
  	
  	leg -> Draw("SAME");

  	c1 -> SaveAs(("miniAOD_vs_AOD_" + var + ".png").c_str());

  	delete c1;
}

void test()
{
	draw_var("isoPhotons_CITK", 0., 40.);
	draw_var("isoNeutralHadrons_CITK", 0., 40.);
	draw_var("isoChargedHadrons_CITK", 0., 40.);
}







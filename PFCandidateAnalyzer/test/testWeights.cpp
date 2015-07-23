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


void testWeights()
{

	TFile file_("/afs/cern.ch/work/i/ishvetso/EgammaWork/ElectronIsolationPUPPI_746_update6July2015/CMSSW_7_4_6_patch2/src/EgammaWork/tree_ele_pfNoLeptons.root");

	TTree * tree = (TTree *)file_.Get("ElectronTree/PFCandidateAnalyzerTree");

	Float_t deltaR, puppiWeight, pdgId, charge;
	Int_t fromPV;
	 vector <float> deltaR_vec, puppiWeight_vec;

	tree -> SetBranchAddress("deltaR", &deltaR);
	tree -> SetBranchAddress("puppiWeight", &puppiWeight);
	tree -> SetBranchAddress("fromPV", &fromPV);
	tree -> SetBranchAddress("pdgId", &pdgId);

	TH2F * hist = new TH2F("hist", "hist", 50, 0.9, 1.0, 50, 0., 1.);

	for (unsigned int iEntry = 0; iEntry < tree -> GetEntries(); iEntry ++)
	{
		tree -> GetEntry(iEntry);
		//std::cout << pdgId << std::endl;
		 if((fromPV == 1 || fromPV == 2 || fromPV == 3))hist -> Fill(puppiWeight, deltaR);
		
		 if((fromPV == 1 || fromPV == 2 || fromPV == 3) && puppiWeight == 1 &&  deltaR < 0.08 && deltaR > 0.02) std::cout << pdgId << " puppiWeight " << puppiWeight << " deltaR " << deltaR << std::endl;
	}

	

	 setTDRStyle(); 
	 gStyle->SetPalette(1);  
  	TCanvas *c1= new TCanvas("c1","canvas",1200,800);
  	//c1->SetLogx(1);

  	c1 -> cd();
  	hist -> GetXaxis() -> SetTitle("puppiWeight");
  	hist -> GetYaxis() -> SetTitle("deltaR");
  	hist -> Draw("COL");

  	c1 -> SaveAs("deltaR_puppiWeight_electrons_pfNoLeptons.png");
  	delete c1;


}
#include <TFile.h>
#include <TTree.h>
#include <TTreeFormula.h>
#include <TGraphAsymmErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <THStack.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TAxis.h>
#include <TGraph.h>
#include <iostream>
#include <iostream>
#include <sstream>
#include <map>
#include <algorithm> 
#include </afs/cern.ch/work/i/ishvetso/GitHub/IvanShvetsov/CMS_stylistics/tdrstyle.C>
#include </afs/cern.ch/work/i/ishvetso/GitHub/IvanShvetsov/CMS_stylistics/CMS_lumi.cpp>



TGraphAsymmErrors * get_graph_bkg(std::string var)
{	
	std::cout << var << std::endl;
	TFile file_("/afs/cern.ch/work/i/ishvetso/EgammaWork/photon_isolations/CMSSW_7_6_4/src/EgammaWork/PhotonNtupler/test/crab_projects/crab_PhotonIsolations/results/GJets.root");
	TTree * tree = (TTree *) file_.Get("ntupler/PhotonTree");

	TH1D * hist_sig = new TH1D(( var + "_sig").c_str(), (var + "_sig").c_str(), 10000., 0., 10.);

	hist_sig -> Sumw2();
	double cutValue;
	
	tree -> Project(( var + "_sig" ).c_str(), var.c_str(),"( pt > 20 && ( isTrue == 1  ) ) && ( (isEB == 1 && hOverE < 0.05 && full5x5_sigmaIetaIeta < 0.01) || (isEB == 0 && hOverE < 0.05 && full5x5_sigmaIetaIeta < 0.0268)  )");


	for (unsigned int iBin = 1; iBin <= hist_sig -> GetNbinsX(); iBin ++){
	
			double _effSig = hist_sig -> Integral(1, iBin)/(hist_sig -> Integral());
			if (_effSig > 0.75) {
				cutValue = (double) hist_sig -> GetBinCenter(iBin);
				std::cout << "eff : " << _effSig << std::endl;
				break;
			}
	}
	//std::string std::to_string(cutValue);
	std::cout << cutValue << std::endl;
	std::cout << "string :" << std::to_string(cutValue) << std::endl;

	TH1D * hist_bkg_nPV_total = new TH1D(( var + "_bkg_nPV_total").c_str(), (var + "_bkg_nPV_total").c_str(), 13, 4.5, 30.5);
	TH1D * hist_bkg_nPV_passed = new TH1D(( var + "_bkg_nPV_passed").c_str(), (var + "_bkg_nPV_passed").c_str(), 13, 4.5, 30.5);

	tree -> Project(( var + "_bkg_nPV_total").c_str(), "nPV","( pt > 20 && isTrue != 1 ) && ( (isEB == 1 && hOverE < 0.05 && full5x5_sigmaIetaIeta < 0.01) || (isEB == 0 && hOverE < 0.05 && full5x5_sigmaIetaIeta < 0.0268)  )");
	tree -> Project(( var + "_bkg_nPV_passed").c_str(), "nPV",(" ( ( pt > 20 && isTrue != 1 ) && ( (isEB == 1 && hOverE < 0.05 && full5x5_sigmaIetaIeta < 0.01) || (isEB == 0 && hOverE < 0.05 && full5x5_sigmaIetaIeta < 0.0268)  ) && " + var + " < " + std::to_string(cutValue) + " )").c_str());

	TGraphAsymmErrors * graph = new TGraphAsymmErrors (hist_bkg_nPV_passed, hist_bkg_nPV_total);

	return graph;


}


TGraphAsymmErrors * get_graph_sig(std::string var)
{	
	std::cout << var << std::endl;
	TFile file_("/afs/cern.ch/work/i/ishvetso/EgammaWork/photon_isolations/CMSSW_7_6_4/src/EgammaWork/PhotonNtupler/test/crab_projects/crab_PhotonIsolations/results/GJets.root");
	TTree * tree = (TTree *) file_.Get("ntupler/PhotonTree");
	TH1D * hist_bkg = new TH1D(( var + "_bkg").c_str(), (var + "_bkg").c_str(), 10000., 0., 10.);

	hist_bkg -> Sumw2();
	double cutValue;
	
	tree -> Project(( var + "_bkg" ).c_str(), var.c_str(),"( pt > 20 && ( isTrue != 1  ) ) && ( (isEB == 1 && hOverE < 0.05 && full5x5_sigmaIetaIeta < 0.01) || (isEB == 0 && hOverE < 0.05 && full5x5_sigmaIetaIeta < 0.0268)  )");

	for (unsigned int iBin = 1; iBin <= hist_bkg -> GetNbinsX(); iBin ++){
	
			double _effBkg = hist_bkg -> Integral(1, iBin)/(hist_bkg -> Integral());
			if (_effBkg > 0.5) {
				cutValue = (double) hist_bkg -> GetBinCenter(iBin);
				std::cout << "eff : " << _effBkg << std::endl;
				break;
			}
	}

	std::cout << cutValue << std::endl;
	std::cout << "string :" << std::to_string(cutValue) << std::endl;

	TH1D * hist_sig_nPV_total = new TH1D(( var + "_sig_nPV_total").c_str(), (var + "_sig_nPV_total").c_str(), 13, 4.5, 30.5);
	TH1D * hist_sig_nPV_passed = new TH1D(( var + "_sig_nPV_passed").c_str(), (var + "_sig_nPV_passed").c_str(), 13, 4.5, 30.5);

	tree -> Project(( var + "_sig_nPV_total").c_str(), "nPV","( pt > 20 && isTrue == 1 ) && ( (isEB == 1 && hOverE < 0.05 && full5x5_sigmaIetaIeta < 0.01) || (isEB == 0 && hOverE < 0.05 && full5x5_sigmaIetaIeta < 0.0268)  )");
	tree -> Project(( var + "_sig_nPV_passed").c_str(), "nPV",(" ( ( pt > 20 && isTrue == 1 ) && ( (isEB == 1 && hOverE < 0.05 && full5x5_sigmaIetaIeta < 0.01) || (isEB == 0 && hOverE < 0.05 && full5x5_sigmaIetaIeta < 0.0268)  ) && " + var + " < " + std::to_string(cutValue) + " )").c_str());
	TGraphAsymmErrors * graph = new TGraphAsymmErrors (hist_sig_nPV_passed, hist_sig_nPV_total);

	return graph;


}


void eff_drawer()
{
	gStyle->SetCanvasColor(kWhite);
	//gStyle->SetTitleFillColor(kWhite);
  	gStyle->SetOptStat(0);
  	gStyle->SetOptTitle(0);
  	setTDRStyle();  
	TCanvas *c1= new TCanvas("c1","canvas",1200,800);	
	
	TGraphAsymmErrors *graph1 = get_graph_bkg("relisoWithEA_CITK");
	TGraphAsymmErrors *graph2 = get_graph_bkg("relisoWithEA_pf");
	TGraphAsymmErrors *graph3 = get_graph_bkg("relisoWithEA_PUPPI");
	TGraphAsymmErrors *graph4 = get_graph_bkg("reliso_raw");
	

	graph1 -> GetYaxis() -> SetTitle("bkg efficiency");
	graph1 -> GetYaxis() -> SetRangeUser(0., 0.8);
  	graph1 -> GetXaxis() -> SetTitle("n_{PV}");

	graph1 -> SetLineColor(kBlue);
	graph2 -> SetLineColor(kGreen);
	graph3 -> SetLineColor(kRed);
	graph4 -> SetLineColor(kBlack);
	

	graph1 -> SetLineWidth(2.);
	graph2 -> SetLineWidth(2.);
	graph3 -> SetLineWidth(2.);
	graph4 -> SetLineWidth(2.);

	TLegend *leg = new TLegend(0.2,0.2,0.5,0.4);
  	leg -> SetFillColor(kWhite);  
  
  	leg->AddEntry(graph1, "effective area, map based veto","l");
  	leg->AddEntry(graph2, "effective area, cone veto","l");
  	leg->AddEntry(graph3, "PUPPI map based veto","l");
  	leg->AddEntry(graph4, "raw, map based veto","l");
  	

	graph1 -> Draw();
	graph2 -> Draw("SAME");
	graph3 -> Draw("SAME");
	graph4 -> Draw("SAME");


	leg -> Draw("SAME");

	CMS_lumi( c1, 4, 0 );
	c1 -> SaveAs("eff_vs_nPV.png");

	
}


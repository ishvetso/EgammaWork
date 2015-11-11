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
	TFile file_sig("/afs/cern.ch/work/i/ishvetso/EgammaWork/ElectronIsolationPUPPI_747_ID_reweighting_PFSelection/CMSSW_7_4_7/src/EgammaWork/ElectronNtupler/test/crab_projects/crab_Electron-Isolation_CITK_PUPPI_DY_matching_by_ref_v6/results/DY.root");
	TFile file_bkg("/afs/cern.ch/work/i/ishvetso/EgammaWork/ElectronIsolationPUPPI_747_ID_reweighting_PFSelection/CMSSW_7_4_7/src/EgammaWork/ElectronNtupler/test/crab_projects/crab_Electron-Isolation_CITK_PUPPI_ttbar_AOD_matching_by_ref_new_splitting/results/ttbar.root");
	TTree * tree_sig = (TTree *) file_sig.Get("ntupler/ElectronTree");
	TTree * tree_bkg = (TTree *) file_bkg.Get("ntupler/ElectronTree");
	TH1D * hist_sig = new TH1D(( var + "_signal").c_str(), (var + "_signal").c_str(), 10000., 0., 10.);
	TH1D * hist_bkg = new TH1D(( var + "_bkg").c_str(), (var + "_bkg").c_str(), 10000., 0., 10.);

	hist_sig -> Sumw2();
	hist_bkg -> Sumw2();
	double cutValue;
	
	tree_sig -> Project(( var + "_signal" ).c_str(), var.c_str(),"(pt > 20 && isTrueElectron == 1  && mvaIDBit_w90 == 1 && PF_ID == 1  )*genWeight");
	//tree_bkg -> Project(( var + "_bkg" ).c_str(), var.c_str(),"pt > 20 && ( isTrueElectron == 0 || isTrueElectron == 3 ) && mvaIDBit_w90 == 1 && PF_ID == 1");

	for (unsigned int iBin = 1; iBin <= hist_sig -> GetNbinsX(); iBin ++){
			double _effSig = hist_sig -> Integral(1, iBin)/(hist_sig -> Integral());
			//double _effBkg = hist_bkg -> Integral(1, iBin)/(hist_bkg -> Integral());
			//std::cout << _effSig << std::endl;
			//std::cout << hist_sig -> GetBinCenter(iBin) << std::endl;
			if (_effSig > 0.95) {
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

	tree_bkg -> Project(( var + "_bkg_nPV_total").c_str(), "nPV","pt > 20. && ( isTrueElectron == 0 || isTrueElectron == 3 ) && mvaIDBit_w90 == 1 && PF_ID == 1");
	tree_bkg -> Project(( var + "_bkg_nPV_passed").c_str(), "nPV",("(pt > 20. && ( isTrueElectron == 0 || isTrueElectron == 3 ) && mvaIDBit_w90 == 1 && PF_ID == 1 && " + var + " < " + std::to_string(cutValue) + " )").c_str());
	TGraphAsymmErrors * graph = new TGraphAsymmErrors (hist_bkg_nPV_passed, hist_bkg_nPV_total);

	return graph;


}


TGraphAsymmErrors * get_graph_sig(std::string var)
{	
	std::cout << var << std::endl;
	TFile file_sig("/afs/cern.ch/work/i/ishvetso/EgammaWork/ElectronIsolationPUPPI_747_ID_reweighting_PFSelection/CMSSW_7_4_7/src/EgammaWork/ElectronNtupler/test/crab_projects/crab_Electron-Isolation_CITK_PUPPI_DY_matching_by_ref_v6/results/DY.root");
	TFile file_bkg("/afs/cern.ch/work/i/ishvetso/EgammaWork/ElectronIsolationPUPPI_747_ID_reweighting_PFSelection/CMSSW_7_4_7/src/EgammaWork/ElectronNtupler/test/crab_projects/crab_Electron-Isolation_CITK_PUPPI_ttbar_AOD_matching_by_ref_new_splitting/results/ttbar.root");
	TTree * tree_sig = (TTree *) file_sig.Get("ntupler/ElectronTree");
	TTree * tree_bkg = (TTree *) file_bkg.Get("ntupler/ElectronTree");
	TH1D * hist_bkg = new TH1D(( var + "_bkg").c_str(), (var + "_bkg").c_str(), 10000., 0., 10.);

	hist_bkg -> Sumw2();
	double cutValue;
	
	tree_bkg -> Project(( var + "_bkg" ).c_str(), var.c_str(),"pt > 20 && ( isTrueElectron == 0 || isTrueElectron == 3 ) && mvaIDBit_w90 == 1 && PF_ID == 1");

	for (unsigned int iBin = 1; iBin <= hist_bkg -> GetNbinsX(); iBin ++){
	
			double _effBkg = hist_bkg -> Integral(1, iBin)/(hist_bkg -> Integral());
			//std::cout << _effSig << std::endl;
			//std::cout << hist_sig -> GetBinCenter(iBin) << std::endl;
			if (_effBkg > 0.05) {
				cutValue = (double) hist_bkg -> GetBinCenter(iBin);
				std::cout << "eff : " << _effBkg << std::endl;
				break;
			}
	}

	std::cout << cutValue << std::endl;
	std::cout << "string :" << std::to_string(cutValue) << std::endl;

	TH1D * hist_sig_nPV_total = new TH1D(( var + "_sig_nPV_total").c_str(), (var + "_sig_nPV_total").c_str(), 13, 4.5, 30.5);
	TH1D * hist_sig_nPV_passed = new TH1D(( var + "_sig_nPV_passed").c_str(), (var + "_sig_nPV_passed").c_str(), 13, 4.5, 30.5);

	tree_sig -> Project(( var + "_sig_nPV_total").c_str(), "nPV","(pt > 20 && isTrueElectron == 1  && mvaIDBit_w90 == 1 && PF_ID == 1  )*genWeight");
	tree_sig -> Project(( var + "_sig_nPV_passed").c_str(), "nPV",("(pt > 20 && isTrueElectron == 1  && mvaIDBit_w90 == 1 && PF_ID == 1  && " + var + " < " + std::to_string(cutValue) + " )*genWeight").c_str());
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
	
	TGraphAsymmErrors *graph1 = get_graph_sig("relIsoWithEA");
	TGraphAsymmErrors *graph2 = get_graph_sig("relIsoWithDBeta");
	TGraphAsymmErrors *graph3 = get_graph_sig("reliso_PUPPI");
	TGraphAsymmErrors *graph4 = get_graph_sig("reliso_PUPPI_NoLeptons");
	TGraphAsymmErrors *graph5 = get_graph_sig("reliso_PUPPI_average");

	graph1 -> GetYaxis() -> SetTitle("signal efficiency");
	graph1 -> GetYaxis() -> SetRangeUser(0.6, 1.1);
  	graph1 -> GetXaxis() -> SetTitle("n_{PV}");

	graph1 -> SetLineColor(kCyan);
	graph2 -> SetLineColor(kBlack);
	graph3 -> SetLineColor(kBlue);
	graph4 -> SetLineColor(kRed);
	graph5 -> SetLineColor(kGreen);

	graph1 -> SetLineWidth(2.);
	graph2 -> SetLineWidth(2.);
	graph3 -> SetLineWidth(2.);
	graph4 -> SetLineWidth(2.);
	graph5 -> SetLineWidth(2.);

	TLegend *leg = new TLegend(0.2,0.7,0.5,0.92);
  	leg -> SetFillColor(kWhite);  
  
  	leg->AddEntry(graph1, "effective area","l");
  	leg->AddEntry(graph2, "#delta#beta-corrected","l");
  	leg->AddEntry(graph3, "PUPPI","l");
  	leg->AddEntry(graph4, "PUPPINoLeptons","l");
  	leg->AddEntry(graph5, "PUPPI combined","l");

	graph1 -> Draw();
	graph2 -> Draw("SAME");
	graph3 -> Draw("SAME");
	graph4 -> Draw("SAME");
	graph5 -> Draw("SAME");

	leg -> Draw("SAME");

	CMS_lumi( c1, 4, 0 );
	c1 -> SaveAs("eff_vs_nPV.png");

	
}


#include <TFile.h>
#include <TTree.h>
#include <TTreeFormula.h>
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



TH1D * get_hist(std::string var)
{	
	std::cout << var << std::endl;
	TFile file("/afs/cern.ch/work/i/ishvetso/EgammaWork/PhotonIsolation_matching/CMSSW_7_4_7/src/EgammaWork/PhotonNtupler/test/crab_projects/crab_PhotonIsolation_PUPPI_AOD_v2/results/GJets.root");
	TTree * tree = (TTree *) file.Get("ntupler/PhotonTree");
	TH2D * hist_sig = new TH2D(( var + "_signal").c_str(), (var + "_signal").c_str(), 26, 4.5, 30.5, 100000., 0., 10.);
	TH2D * hist_bkg = new TH2D(( var + "_bkg").c_str(), (var + "_bkg").c_str(), 26, 4.5, 30.5, 100000., 0., 10.);
	
	tree -> Project(( var + "_signal" ).c_str(), (var + ":nPV").c_str(),"isTrue == 1 && PF_ID == 1 && pt > 20");
	tree -> Project(( var + "_bkg" ).c_str(), (var + ":nPV").c_str(),"isTrue != 1 && PF_ID == 1 && pt > 20");

	TH1D * eff_hist = new TH1D("eff", "eff", 26, 4.5, 30.5);
	eff_hist -> SetDirectory(0);
	double eff = 0.;
	double unc = 0;
	bool EffCalculated = false;

	for (unsigned iPUBin = 1; iPUBin <= 26 ; iPUBin ++) {
		TH1D * _hist_sig = (TH1D *) hist_sig -> ProjectionY("sig", iPUBin  , iPUBin );
		TH1D * _hist_bkg = (TH1D *) hist_bkg -> ProjectionY("bkg", iPUBin  , iPUBin );

		for (unsigned int iBin = 1; iBin < _hist_sig -> GetNbinsX(); iBin ++){
			double _effSig = _hist_sig -> Integral(1, iBin)/(_hist_sig -> Integral());
			double _effBkg = _hist_bkg -> Integral(1, iBin)/(_hist_bkg -> Integral());
			
			
			
			if (_effSig > 0.7) {
				eff = _effBkg;
				double _effBkg_unc_up =  std::abs((_hist_bkg -> Integral(1, iBin) + sqrt(_hist_bkg -> Integral(1, iBin)))/(_hist_bkg -> Integral() + sqrt(_hist_bkg -> Integral())) - _effBkg);
				double _effBkg_unc_down = std::abs((_hist_bkg -> Integral(1, iBin) - sqrt(_hist_bkg -> Integral(1, iBin)))/(_hist_bkg -> Integral() - sqrt(_hist_bkg -> Integral())) - _effBkg );
				unc = std::max(_effBkg_unc_down, _effBkg_unc_up);
				std::cout << "iBin " << iPUBin <<  " sig: " << _effSig << " bkg : " << _effBkg  << " unc " << unc << std::endl;
				EffCalculated = true;
				break;
			}
			 
		}
		if (!EffCalculated) {
			std::cerr << "effeciency wasn't calculated for bin " << iPUBin << std::endl;
			exit(0); 
		}
		eff_hist -> SetBinContent(iPUBin, eff);
		eff_hist -> SetBinError(iPUBin, unc);
	}
	
	
	return eff_hist;


}


void eff_drawer()
{
	gStyle->SetCanvasColor(kWhite);
	//gStyle->SetTitleFillColor(kWhite);
  	gStyle->SetOptStat(0);
  	gStyle->SetOptTitle(0);
  	setTDRStyle();  
	TCanvas *c1= new TCanvas("c1","canvas",1200,800);	
	TH1D * hist1 = (TH1D*) get_hist("relisoWithEA_PUPPI");
	TH1D * hist2 =  get_hist("relisoWithEA_pf");
	TH1D * hist3 =  get_hist("relisoWithEA_CITK_");

	
	hist1 -> GetYaxis() -> SetRangeUser(0.15, 0.85);

	hist1 -> GetYaxis() -> SetTitle("eff_{bkg}(eff_{sig} = 0.7 )");
  	hist1 -> GetXaxis() -> SetTitle("n_{PV}");

	hist1 -> SetLineWidth(2.);
  	hist2 -> SetLineWidth(2.);
  	hist3 -> SetLineWidth(2.);

  	hist1 -> SetMarkerStyle(22);
  	hist2 -> SetMarkerStyle(22);
  	hist3 -> SetMarkerStyle(22);

  	hist1 -> SetMarkerColor(kRed);
  	hist2 -> SetMarkerColor(kBlue);
  	hist3 -> SetMarkerColor(kGreen);

	hist1 -> SetLineColor(kRed);
	hist2 -> SetLineColor(kBlue);
	hist3 -> SetLineColor(kGreen);

	hist1 -> Draw("LPE");
	hist2 -> Draw("LPESAME");
	hist3 -> Draw("LPESAME");

	TLegend *leg = new TLegend(0.2,0.7,0.4,0.9);
  	leg -> SetFillColor(kWhite);  
  
 	leg->AddEntry(hist1, "PUPPI","l");
  	leg->AddEntry(hist2, "PF","l");
  	leg->AddEntry(hist3, "CITK","l");

  	leg -> Draw("SAME");


	CMS_lumi( c1, 4, 0 );
	c1 -> SaveAs("iso_vs_nPV.png");
}


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
#include <iostream>
#include <iostream>
#include <sstream>
#include <map>
#include <algorithm> 
#include </afs/cern.ch/work/i/ishvetso/GitHub/IvanShvetsov/CMS_stylistics/tdrstyle.C>
#include </afs/cern.ch/work/i/ishvetso/GitHub/IvanShvetsov/CMS_stylistics/CMS_lumi.cpp>
void simple_count()
{
  TFile file("/afs/cern.ch/work/i/ishvetso/EgammaWork/Validation_ElectronIsolation/CMSSW_7_3_3/src/EgammaWork/ElectronNtupler/test/crab_projects/crab_Electron-Isolation_CITK_validation_ttbar_AOD_2June2015/results/ttbar.root", "READ");

  TTree * tree = (TTree * )file.Get("ntupler/ElectronTree");
  
  float isoCH, isoCH_CITK;
  Char_t isEB;
  int count_positive_barrel = 0, count_negative_barrel = 0, count_positive_endcap = 0, count_negative_endcap = 0;;
  
  tree -> SetBranchAddress("isoChargedHadrons", &isoCH);
  tree -> SetBranchAddress("sumChargedHadronPt_CITK", &isoCH_CITK);
  tree -> SetBranchAddress("isEB", &isEB);
  
  for (unsigned int iEntry = 0; iEntry < tree -> GetEntries(); iEntry ++)
  {
    tree -> GetEntry(iEntry);
    
    if (((isoCH - isoCH_CITK) < -0.01) && isEB) count_negative_barrel ++;
    if (((isoCH - isoCH_CITK) > 0.01) && isEB) count_positive_barrel ++;
    if (((isoCH - isoCH_CITK) < -0.01) && !isEB) count_negative_endcap ++;
    if (((isoCH - isoCH_CITK) > 0.01) && !isEB) count_positive_endcap ++;
    
  }
  
  std::cout << "positive barrel " << count_positive_barrel << std::endl;
  std::cout << "negative barrel " << count_negative_barrel << std::endl;
  std::cout << "positive endcap " << count_positive_endcap << std::endl;
  std::cout << "negative endcap " << count_negative_endcap << std::endl;
  std::cout << "total " << tree -> GetEntries() << std::endl;
  
  
}
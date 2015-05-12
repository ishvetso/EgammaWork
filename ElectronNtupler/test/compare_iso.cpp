#include "draw.hpp"

void compare_iso()
{
  Sample ttbar, DY;

  
  DY.SetParameters("Drell-Yan");
  DY.SetFileNames("/afs/cern.ch/work/i/ishvetso/EgammaWork/my_puppi_test/CMSSW_7_3_3/src/EgammaWork/ElectronNtupler/test/crab_projects/crab_Electron-Isolation_CITK_validation_DY_miniAOD_PUPPI_with_NoLeptons_update11May2015/results/DY.root");
  
  ttbar.SetParameters("ttbar");
  ttbar.SetFileNames("/afs/cern.ch/work/i/ishvetso/EgammaWork/my_puppi_test/CMSSW_7_3_3/src/EgammaWork/ElectronNtupler/test/crab_projects/crab_Electron-Isolation_CITK_validation_ttbar_miniAOD_PUPPI_with_NoLeptons_update11May2015/results/ttbar.root");
    
  vector <Var> Vars_CH, Vars_NH, Vars_Gamma;
  Var var;
  
  //charged hadrons
  var.VarName = "sumChargedHadronPt_CITK";
  var.color = kBlue;
  var.SetRange(0., 1.5, 30);
  Vars_CH.push_back(var);
  
  var.VarName = "sumChargedHadronPt_PUPPI";
  var.color = kGreen;
  var.SetRange(0., 1.5, 30);
  Vars_CH.push_back(var);
  
  var.VarName = "sumChargedHadronPt_PUPPI_NoLeptons";
  var.color = kRed;
  var.SetRange(0., 1.5, 30);
  Vars_CH.push_back(var);
  
  //neutral hadrons
  var.VarName = "sumNeutralHadronPt_CITK";
  var.color = kBlue;
  var.SetRange(0., 1.5, 30);
  Vars_NH.push_back(var);
  
  var.VarName = "sumNeutralHadronPt_PUPPI";
  var.color = kGreen;
  var.SetRange(0., 1.5, 30);
  Vars_NH.push_back(var);
  
  var.VarName = "sumNeutralHadronPt_PUPPI_NoLeptons";
  var.color = kRed;
  var.SetRange(0., 1.5, 30);
  Vars_NH.push_back(var);
  
  //photons
  var.VarName = "sumPhotonPt_CITK";
  var.color = kBlue;
  var.SetRange(0., 1.5, 30);
  Vars_Gamma.push_back(var);
  
  var.VarName = "sumPhotonPt_PUPPI";
  var.color = kGreen;
  var.SetRange(0., 1.5, 30);
  Vars_Gamma.push_back(var);
  
  var.VarName = "sumPhotonPt_PUPPI_NoLeptons";
  var.color = kRed;
  var.SetRange(0., 1.5, 30);
  Vars_Gamma.push_back(var);

  
  setTDRStyle();
  draw(Vars_CH, ttbar, "ch", "barrel",  "isEB == 1", "iso_matching");
  draw(Vars_CH, ttbar, "ch", "endcap",  "isEB != 1", "iso_matching");
  draw(Vars_CH, DY, "ch", "endcap",  "isEB != 1", "iso_matching");
  draw(Vars_CH, DY, "ch", "barrel",  "isEB == 1", "iso_matching");
  
  draw(Vars_NH, ttbar, "nh", "barrel",  "isEB == 1", "iso_matching");
  draw(Vars_NH, ttbar, "nh", "endcap",  "isEB != 1", "iso_matching");
  draw(Vars_NH, DY, "nh", "endcap",  "isEB != 1", "iso_matching");
  draw(Vars_NH, DY, "nh", "barrel",  "isEB == 1", "iso_matching");
  
  draw(Vars_Gamma  , ttbar, "gamma", "barrel",  "isEB == 1", "iso_matching");
  draw(Vars_Gamma, ttbar, "gamma", "endcap",  "isEB != 1", "iso_matching");
  draw(Vars_Gamma, DY, "gamma", "endcap",  "isEB != 1", "iso_matching");
  draw(Vars_Gamma, DY, "gamma", "barrel",  "isEB == 1", "iso_matching");
 
 
  
   
}


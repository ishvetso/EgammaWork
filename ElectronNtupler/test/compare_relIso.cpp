#include "draw.hpp"

void compare_relIso()
{
  Sample ttbar, DY;

  
  DY.SetParameters("Drell-Yan");
  DY.SetFileNames("/afs/cern.ch/work/i/ishvetso/EgammaWork/my_puppi_test/CMSSW_7_3_3/src/EgammaWork/ElectronNtupler/test/crab_projects/crab_Electron-Isolation_CITK_validation_DY_miniAOD_PUPPI_with_NoLeptons_update11May2015/results/DY.root");
  
  ttbar.SetParameters("ttbar");
  ttbar.SetFileNames("/afs/cern.ch/work/i/ishvetso/EgammaWork/my_puppi_test/CMSSW_7_3_3/src/EgammaWork/ElectronNtupler/test/crab_projects/crab_Electron-Isolation_CITK_validation_ttbar_miniAOD_PUPPI_with_NoLeptons_update11May2015/results/ttbar.root");
    
  vector <Var> relIsoVars_DY, relIsoVars_ttbar;
  
  Var var;
  
  var.VarName = "relIsoWithEA";
  var.color = kBlue;
  var.SetRange(0., 4., 30);
  relIsoVars_ttbar.push_back(var);
  var.SetRange(0., 1.5, 30);
  relIsoVars_DY.push_back(var);
  
  var.VarName = "relIsoWithDBeta";
  var.color = kGreen;
  var.SetRange(0., 4., 30);
  relIsoVars_ttbar.push_back(var);
  var.SetRange(0., 1.5, 30);
  relIsoVars_DY.push_back(var);
  
  var.VarName = "reliso_PUPPI";
  var.color = kRed;
  var.SetRange(0., 4., 30);
  relIsoVars_ttbar.push_back(var);
  var.SetRange(0., 1.5, 30);
  relIsoVars_DY.push_back(var);
  
  var.VarName = "reliso_PUPPI_NoLeptons";
  var.color = kOrange;
  var.SetRange(0., 4., 30);
  relIsoVars_ttbar.push_back(var);
  var.SetRange(0., 1.5, 30);
  relIsoVars_DY.push_back(var);
  
  
  setTDRStyle();
  draw(relIsoVars_ttbar, ttbar, "relIso", "barrel",  "isEB == 1 && pt > 20 && ( isTrueElectron == 0 || isTrueElectron == 3 )", "relIso_distributions_matching");
  draw(relIsoVars_DY, DY, "relIso", "barrel", "isEB == 1 && pt > 20 && isTrueElectron == 1 ", "relIso_distributions_matching");
  draw(relIsoVars_ttbar, ttbar, "relIso", "endcap", "isEB == 0 && pt > 20 && ( isTrueElectron == 0 || isTrueElectron == 3)", "relIso_distributions_matching");
  draw(relIsoVars_DY, DY, "relIso", "endcap", "isEB == 0 && pt > 20 && isTrueElectron == 1", "relIso_distributions_matching");
   
}
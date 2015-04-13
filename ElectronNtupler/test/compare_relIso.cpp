#include "draw.hpp"

void compare_relIso()
{
  Sample ttbar, DY;

  
  DY.SetParameters("Drell-Yan");
  DY.SetFileNames("/afs/cern.ch/work/i/ishvetso/EgammaWork/ElectronIsolationPUPPI/CMSSW_7_3_0/src/EgammaWork/ElectronNtupler/test/crab_projects/crab_Electron-Isolation_CITK_validation_DY_miniAOD_PUPPI/results/DY.root");
  
  ttbar.SetParameters("ttbar");
  ttbar.SetFileNames("/afs/cern.ch/work/i/ishvetso/EgammaWork/ElectronIsolationPUPPI/CMSSW_7_3_0/src/EgammaWork/ElectronNtupler/test/crab_projects/crab_Electron-Isolation_CITK_validation_ttbar_miniAOD_PUPPI/results/ttbar.root");
    
  vector <Var> relIsoVars;
  Var var;
  
  var.VarName = "relIsoWithEA";
  var.color = kBlue;
  var.SetRange(0., 4., 30);
  relIsoVars.push_back(var);
  
  var.VarName = "relIsoWithDBeta";
  var.color = kGreen;
  var.SetRange(0., 4., 30);
  relIsoVars.push_back(var);
  
  var.VarName = "reliso_PUPPI";
  var.color = kRed;
  var.SetRange(0., 4., 30);
  relIsoVars.push_back(var);
  
  
  setTDRStyle();
  draw(relIsoVars, ttbar, "relIso", "barrel",  "isEB == 1 && pt > 20 && ( isTrueElectron == 0 || isTrueElectron == 3 )", "relIso_distributions_matching");
  draw(relIsoVars, DY, "relIso", "barrel", "isEB == 1 && pt > 20 && isTrueElectron == 1 ", "relIso_distributions_matching");
  draw(relIsoVars, ttbar, "relIso", "endcap", "isEB == 0 && pt > 20 && ( isTrueElectron == 0 || isTrueElectron == 3)", "relIso_distributions_matching");
  draw(relIsoVars, DY, "relIso", "endcap", "isEB == 0 && pt > 20 && isTrueElectron == 1", "relIso_distributions_matching");
   
}
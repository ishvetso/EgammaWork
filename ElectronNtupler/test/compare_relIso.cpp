#include "draw.hpp"

void compare_relIso()
{
  Sample ttbar, DY;

  
  DY.SetParameters("Drell-Yan");
  DY.SetFileNames("/afs/cern.ch/work/i/ishvetso/EgammaWork/ElectronIsolationPUPPI_746_update6July2015/CMSSW_7_4_6_patch2/src/EgammaWork/ElectronNtupler/test/crab_projects/crab_Electron-Isolation_CITK_validation_DY_miniAOD_PUPPI_with_NoLeptons_update13July2015_new_selection/results/DY.root");
  
  ttbar.SetParameters("ttbar");
  ttbar.SetFileNames("/afs/cern.ch/work/i/ishvetso/EgammaWork/ElectronIsolationPUPPI_746_update6July2015/CMSSW_7_4_6_patch2/src/EgammaWork/ElectronNtupler/test/crab_projects/crab_Electron-Isolation_CITK_validation_ttbar_miniAOD_PUPPI_with_NoLeptons_update13July2015_new_selection/results/ttbar.root");
    
  vector <Var> relIsoVars_DY, relIsoVars_ttbar;
  
  Var var;
  

  var.VarName = "reliso_raw";
  var.color = kBlack;
  var.SetRange(0., 4., 30);
  relIsoVars_ttbar.push_back(var);
  var.SetRange(0., 0.5, 30);
  relIsoVars_DY.push_back(var);

  var.VarName = "relIsoWithEA";
  var.color = kBlue;
  var.SetRange(0., 4., 30);
  relIsoVars_ttbar.push_back(var);
  var.SetRange(0., 0.5, 30);
  relIsoVars_DY.push_back(var);
  
  var.VarName = "relIsoWithDBeta";
  var.color = kGreen;
  var.SetRange(0., 4., 30);
  relIsoVars_ttbar.push_back(var);
  var.SetRange(0., 0.5, 30);
  relIsoVars_DY.push_back(var);


  var.VarName = "relIsoWithEA_MapBasedVeto";
  var.color = kMagenta;
  var.SetRange(0., 4., 30);
  relIsoVars_ttbar.push_back(var);
  var.SetRange(0., 0.5, 30);
  relIsoVars_DY.push_back(var);
  
  var.VarName = "reliso_PUPPI";
  var.color = kRed;
  var.SetRange(0., 4., 30);
  relIsoVars_ttbar.push_back(var);
  var.SetRange(0., 0.5, 30);
  relIsoVars_DY.push_back(var);
  
  var.VarName = "reliso_PUPPI_NoLeptons";
  var.color = kOrange;
  var.SetRange(0., 4., 30);
  relIsoVars_ttbar.push_back(var);
  var.SetRange(0., 0.5, 30);
  relIsoVars_DY.push_back(var);
  
  
  setTDRStyle();
  draw(relIsoVars_ttbar, ttbar, "relIso", "barrel",  "isEB == 1 && pt > 20 && ( isTrueElectron == 0 || isTrueElectron == 3 )", "relIso_distributions_matching_new_selection");
  draw(relIsoVars_DY, DY, "relIso", "barrel", "isEB == 1 && pt > 20 && isTrueElectron == 1 ", "relIso_distributions_matching_new_selection");
  draw(relIsoVars_ttbar, ttbar, "relIso", "endcap", "isEB == 0 && pt > 20 && ( isTrueElectron == 0 || isTrueElectron == 3)", "relIso_distributions_matching_new_selection");
  draw(relIsoVars_DY, DY, "relIso", "endcap", "isEB == 0 && pt > 20 && isTrueElectron == 1", "relIso_distributions_matching_new_selection");
   
}
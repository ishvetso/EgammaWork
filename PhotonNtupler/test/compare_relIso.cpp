#include "draw.hpp"

void compare_relIso()
{
  Sample GJets;

  
  GJets.SetParameters("GJets");
  GJets.SetFileNames("/afs/cern.ch/work/i/ishvetso/EgammaWork/PhotonIsolation_PUPPI_746_8July/CMSSW_7_4_6_patch2/src/EgammaWork/PhotonNtupler/test/crab_projects/crab_PhotonIsolation_PUPPI_miniAOD_8July2015/results/GJets.root");
  
  
  vector <Var> relIsoVars_true, relIsoVars_fakes;
  
  Var var;
  
  var.VarName = "relisoWithEA";
  var.color = kBlue;
  var.SetRange(0., 2.5, 30);
  relIsoVars_fakes.push_back(var);
  var.SetRange(0., 0.4, 30);
  relIsoVars_true.push_back(var);

  var.VarName = "relisoWithEA_CITK_";
  var.color = kGreen;
  var.SetRange(0., 2.5, 30);
  relIsoVars_fakes.push_back(var);
  var.SetRange(0., 0.4, 30);
  relIsoVars_true.push_back(var);

  
  var.VarName = "relisoWithEA_PUPPI";
  var.color = kRed;
  var.SetRange(0., 2.5, 30);
  relIsoVars_fakes.push_back(var);
  var.SetRange(0., 0.4, 30);
  relIsoVars_true.push_back(var);
  
  var.VarName = "relisoWithEA_pf";
  var.color = kOrange;
  var.SetRange(0., 2.5, 30);
  relIsoVars_fakes.push_back(var);
  var.SetRange(0., 0.4, 30);
  relIsoVars_true.push_back(var);
  
  setTDRStyle();
  draw(relIsoVars_true, GJets, "relIso", "true",  "pt > 20 && ( isTrue == 1 )", "relIso_distributions_matching");
  draw(relIsoVars_fakes, GJets, "relIso", "fakes",  "pt > 20 && ( isTrue == 0 )", "relIso_distributions_matching");
}
  
  
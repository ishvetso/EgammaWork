#include "draw.hpp"

void compare_iso()
{
  Sample GJets; 
  GJets.SetParameters("GJets");
  GJets.SetFileNames("/afs/cern.ch/work/i/ishvetso/EgammaWork/PhotonIsolation_PUPPI_746_8July/CMSSW_7_4_6_patch2/src/EgammaWork/PhotonNtupler/test/crab_projects/crab_PhotonIsolation_PUPPI_miniAOD_8July2015/results/GJets.root");

  vector <Var> Vars_CH, Vars_NH, Vars_Gamma;
  Var var;
  
  //charged hadrons
  var.VarName = "isoChargedHadrons";
  var.color = kBlue;
  var.SetRange(0., 1.5, 30);
  Vars_CH.push_back(var);
  
  var.VarName = "isoChargedHadrons_CITK";
  var.color = kGreen;
  var.SetRange(0., 1.5, 30);
  Vars_CH.push_back(var);
  
  //neutral hadrons
  var.VarName = "isoNeutralHadrons";
  var.color = kBlue;
  var.SetRange(0., 1.5, 30);
  Vars_NH.push_back(var);
  
  var.VarName = "isoNeutralHadrons_CITK";
  var.color = kGreen;
  var.SetRange(0., 1.5, 30);
  Vars_NH.push_back(var);
  
  //photons
  var.VarName = "isoPhotons";
  var.color = kBlue;
  var.SetRange(0., 1.5, 30);
  Vars_Gamma.push_back(var);
  
  var.VarName = "isoPhotons_CITK";
  var.color = kGreen;
  var.SetRange(0., 1.5, 30);
  Vars_Gamma.push_back(var);


  
  setTDRStyle();
 draw(Vars_CH, GJets, "ch", "charged_hadrons",  "1", "isolation_comparison");
 draw(Vars_NH, GJets, "nh", "neutral_hadrons",  "1", "isolation_comparison");
 draw(Vars_Gamma, GJets, "gamma", "gamma",  "1", "isolation_comparison");

  /*draw_difference("isoPhotons","isoPhotons_CITK", GJets, 11, -0.02, 0.02, "Photons", "diff_CITK_PhotonIDValueMapProducer");
  draw_difference("isoChargedHadrons","isoChargedHadrons_CITK", GJets, 11, -0.02, 0.02, "CH", "diff_CITK_PhotonIDValueMapProducer");
  draw_difference("isoNeutralHadrons","isoNeutralHadrons_CITK", GJets, 11, -0.02, 0.02, "NH", "diff_CITK_PhotonIDValueMapProducer"); 

  draw_difference2D("isoNeutralHadrons","isoNeutralHadrons_CITK","pt", GJets, 11, -250.0, 250.0,40, 20., 300., "NeutralHadrons", "CITK_diff_2D");
  draw_difference2D("isoChargedHadrons","isoChargedHadrons_CITK","pt", GJets, 11, -250.0, 250.0,40, 20., 300., "ChargedHadrons", "CITK_diff_2D");
  draw_difference2D("isoPhotons","isoPhotons_CITK","pt", GJets, 11, -250.0, 250.0,40, 20., 300., "Photons", "CITK_diff_2D");*/
   
}


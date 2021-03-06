#include "draw.hpp"

void draw_difference()
{
   Sample ttbar, DY;

  
  DY.SetParameters("Drell-Yan");
  DY.SetFileNames("/afs/cern.ch/work/i/ishvetso/EgammaWork/Validation_ElectronIsolation/CMSSW_7_3_3/src/EgammaWork/ElectronNtupler/test/crab_projects/crab_Electron-Isolation_CITK_validation_DY_AOD_2June2015/results/DY.root");
  
  ttbar.SetParameters("ttbar");
  ttbar.SetFileNames("/afs/cern.ch/work/i/ishvetso/EgammaWork/Validation_ElectronIsolation/CMSSW_7_3_3/src/EgammaWork/ElectronNtupler/test/crab_projects/crab_Electron-Isolation_CITK_validation_ttbar_AOD_2June2015/results/ttbar.root");
  
  setTDRStyle();
   
  draw_difference2D("isoPhotons","sumPhotonPt_CITK","pt", ttbar, 11, -100.002, 100.002,20, 20., 300., "Photons", "CITK_diff_2D");
  draw_difference2D("isoPhotons","sumPhotonPt_CITK","etaSC", ttbar, 11, -100.002, 100.002,20, -3., 3., "Photons", "CITK_diff_2D");
  
  draw_difference2D("isoChargedHadrons","sumChargedHadronPt_CITK","pt", ttbar, 11, -100.002, 100.002,20, 20., 300., "ChargedHadrons", "CITK_diff_2D");
  draw_difference2D("isoChargedHadrons","sumChargedHadronPt_CITK","etaSC", ttbar, 11, -100.002, 100.002,20, -3., 3., "ChargedHadrons", "CITK_diff_2D");
  
  draw_difference2D("isoNeutralHadrons","sumNeutralHadronPt_CITK","pt", ttbar, 11, -100.002, 100.002,20, 20., 300., "NeutralHadrons", "CITK_diff_2D");
  draw_difference2D("isoNeutralHadrons","sumNeutralHadronPt_CITK","etaSC", ttbar, 11, -100.002, 100.002,20, -3., 3., "NeutralHadrons", "CITK_diff_2D");
  
  draw_difference("isoPhotons","sumPhotonPt_CITK", DY, 11, -0.02, 0.02, "Photons", "CITK_diff_miniAOD_1dim");
  draw_difference("isoChargedHadrons","sumChargedHadronPt_CITK", DY, 11, -0.02, 0.02, "ChargedHadrons", "CITK_diff_miniAOD_1dim");
  draw_difference("isoNeutralHadrons","sumNeutralHadronPt_CITK", DY, 11, -0.02, 0.02, "NeutralHadrons", "CITK_diff_miniAOD_1dim");
 
  
  draw_difference("isoPhotons","sumPhotonPt_CITK", ttbar, 11, -0.02, 0.02, "Photons", "CITK_diff_miniAOD_1dim");
  draw_difference("isoChargedHadrons","sumChargedHadronPt_CITK", ttbar, 11, -0.02, 0.02, "ChargedHadrons", "CITK_diff_miniAOD_1dim");
  draw_difference("isoNeutralHadrons","sumNeutralHadronPt_CITK", ttbar, 11, -0.02, 0.02, "NeutralHadrons", "CITK_diff_miniAOD_1dim");
  
  draw_difference2D("isoChargedHadrons","sumChargedHadronPt_CITK","pt", DY, 11, -250.0, 250.0,20, 20., 300., "ChargedHadrons", "CITK_diff_2D_AOD");
  draw_difference2D("isoChargedHadrons","sumChargedHadronPt_CITK","etaSC", DY, 11, -250.0, 250.0,20, -3., 3., "ChargedHadrons", "CITK_diff_2D_AOD");
  
  draw_difference2D("isoNeutralHadrons","sumNeutralHadronPt_CITK","pt", DY, 11, -250.0, 250.0,20, 20., 300., "NeutralHadrons", "CITK_diff_2D_AOD");
  draw_difference2D("isoNeutralHadrons","sumNeutralHadronPt_CITK","etaSC", DY, 11, -250.0, 250.0,20, -3., 3., "NeutralHadrons", "CITK_diff_2D_AOD");
  
  draw_difference2D("isoPhotons","sumPhotonPt_CITK","pt", DY, 11, -250.0, 250.0,20, 20., 300., "Photons", "CITK_diff_2D_AOD");
  draw_difference2D("isoPhotons","sumPhotonPt_CITK","etaSC", DY, 11, -250.0, 250.0,20, -3., 3., "Photons", "CITK_diff_2D_AOD");
  
  draw_difference2D("isoChargedHadrons","sumChargedHadronPt_CITK","pt", ttbar, 11, -250.0, 250.0,20, 20., 300., "ChargedHadrons", "CITK_diff_2D_AOD");
  draw_difference2D("isoChargedHadrons","sumChargedHadronPt_CITK","etaSC", ttbar, 11, -250.0, 250.0,20, -3., 3., "ChargedHadrons", "CITK_diff_2D_AOD");
  
  draw_difference2D("isoNeutralHadrons","sumNeutralHadronPt_CITK","pt", ttbar, 11, -250.0, 250.0,20, 20., 300., "NeutralHadrons", "CITK_diff_2D_AOD");
  draw_difference2D("isoNeutralHadrons","sumNeutralHadronPt_CITK","etaSC", ttbar, 11, -250.0, 250.0,20, -3., 3., "NeutralHadrons", "CITK_diff_2D_AOD");
  
  draw_difference2D("isoPhotons","sumPhotonPt_CITK","pt", ttbar, 11, -250.0, 250.0,20, 20., 300., "Photons", "CITK_diff_2D_AOD");
  draw_difference2D("isoPhotons","sumPhotonPt_CITK","etaSC", ttbar, 11, -250.0, 250.0,20, -3., 3., "Photons", "CITK_diff_2D_AOD");
 
}
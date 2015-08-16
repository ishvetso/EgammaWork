#include "ROC_drawer.cpp"


void draw_ROC()
{
  ROC_Drawer MyROCDrawer;
  MyROCDrawer.bkgFileName = "/afs/cern.ch/work/i/ishvetso/EgammaWork/ElectronIsolationPUPPI_747_ID_reweighting_PFSelection/CMSSW_7_4_7/src/EgammaWork/ElectronNtupler/test/crab_projects/crab_Electron-Isolation_CITK_validation_ttbar_miniAOD_PUPPI_with_NoLeptons_ids_selection_track_matching_16August2015/results/ttbar.root";
  MyROCDrawer.sigFileName = "/afs/cern.ch/work/i/ishvetso/EgammaWork/ElectronIsolationPUPPI_747_ID_reweighting_PFSelection/CMSSW_7_4_7/src/EgammaWork/ElectronNtupler/test/crab_projects/crab_Electron-Isolation_CITK_validation_DY_miniAOD_PUPPI_with_NoLeptons_ids_selection_track_matching_16August2015/results/DY.root";
  MyROCDrawer.SigSelection = "(pt > 20 && isTrueElectron == 1 && mvaIDBit_w90 == 1)*genWeight";
  MyROCDrawer.BkgSelection = "pt > 20 && ( isTrueElectron == 0 || isTrueElectron == 3 ) && mvaIDBit_w90 == 1";
  MyROCDrawer.name = "ROC";
  MyROCDrawer.addSelection = "1";
  MyROCDrawer.draw_ROC();
  
}

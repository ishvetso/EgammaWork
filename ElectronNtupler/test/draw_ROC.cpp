#include "ROC_drawer.cpp"


void draw_ROC()
{
  ROC_Drawer MyROCDrawer;
  MyROCDrawer.bkgFileName = "/afs/cern.ch/work/i/ishvetso/EgammaWork/ElectronIsolationPUPPI_746_update6July2015/CMSSW_7_4_6_patch2/src/EgammaWork/ElectronNtupler/test/crab_projects/crab_Electron-Isolation_CITK_validation_ttbar_miniAOD_PUPPI_with_NoLeptons_update6July2015/results/ttbar.root";
  MyROCDrawer.sigFileName = "/afs/cern.ch/work/i/ishvetso/EgammaWork/ElectronIsolationPUPPI_746_update6July2015/CMSSW_7_4_6_patch2/src/EgammaWork/ElectronNtupler/test/crab_projects/crab_Electron-Isolation_CITK_validation_DY_miniAOD_PUPPI_with_NoLeptons_update6July2015/results/DY.root";
  MyROCDrawer.SigSelection = "pt > 20 && isTrueElectron == 1";
  MyROCDrawer.BkgSelection = "pt > 20 && ( isTrueElectron == 0 || isTrueElectron == 3 )";
  MyROCDrawer.name = "ROC_endcap";
  MyROCDrawer.addSelection = " isEB != 1";
  MyROCDrawer.draw_ROC();
  MyROCDrawer.name = "ROC_barrel";
  MyROCDrawer.addSelection = " isEB == 1";
  MyROCDrawer.draw_ROC();
}
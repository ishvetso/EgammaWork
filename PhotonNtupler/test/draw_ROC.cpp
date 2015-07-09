#include "ROC_drawer.cpp"


void draw_ROC()
{
  ROC_Drawer MyROCDrawer;
  MyROCDrawer.FileName = "/afs/cern.ch/work/i/ishvetso/EgammaWork/PhotonIsolation_PUPPI_746_8July/CMSSW_7_4_6_patch2/src/EgammaWork/PhotonNtupler/test/crab_projects/crab_PhotonIsolation_PUPPI_miniAOD_8July2015/results/GJets.root";
  MyROCDrawer.SigSelection = "pt > 20 && isTrue == 1";
  MyROCDrawer.BkgSelection = "pt > 20 && ( isTrue == 0  )";
  MyROCDrawer.name = "ROC";
  MyROCDrawer.addSelection = "1";
  MyROCDrawer.draw_ROC();
}

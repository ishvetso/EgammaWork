#include "ROC_drawer.cpp"


void draw_ROC()
{
  ROC_Drawer MyROCDrawer;
  MyROCDrawer.bkgFileName = "/afs/cern.ch/work/i/ishvetso/EgammaWork/MuonIsolation/CMSSW_7_4_6_patch2/src/EgammaWork/MuonNtupler/test/crab_projects/crab_MuonIsolation_First/results/DY.root";
  MyROCDrawer.sigFileName = "/afs/cern.ch/work/i/ishvetso/EgammaWork/MuonIsolation/CMSSW_7_4_6_patch2/src/EgammaWork/MuonNtupler/test/crab_projects/crab_MuonIsolation_First_QCD/results/QCD.root";
  MyROCDrawer.SigSelection = "pt > 20. && fabs(eta) < 2.4";
  MyROCDrawer.BkgSelection = "pt > 20. && fabs(eta) < 2.4 ";
  MyROCDrawer.name = "ROC";
  MyROCDrawer.addSelection = "1";
  MyROCDrawer.draw_ROC();
}

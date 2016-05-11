#include "ROC_drawer.cpp"


void draw_ROC()
{
  ROC_Drawer MyROCDrawer;
  MyROCDrawer.FileName = "/afs/cern.ch/work/i/ishvetso/EgammaWork/photon_isolations/CMSSW_7_6_4/src/EgammaWork/PhotonNtupler/test/crab_projects/crab_PhotonIsolations/results/GJets.root";
  MyROCDrawer.SigSelection = "( pt > 20 && isTrue == 1 ) && ( (isEB == 1 && hOverE < 0.05 && full5x5_sigmaIetaIeta < 0.01) || (isEB == 0 && hOverE < 0.05 && full5x5_sigmaIetaIeta < 0.0268)  )";
  MyROCDrawer.BkgSelection = "( pt > 20 && ( isTrue != 1  ) ) && ( (isEB == 1 && hOverE < 0.05 && full5x5_sigmaIetaIeta < 0.01) || (isEB == 0 && hOverE < 0.05 && full5x5_sigmaIetaIeta < 0.0268)  )";
  MyROCDrawer.name = "ROC";
  MyROCDrawer.addSelection = "1";
  MyROCDrawer.draw_ROC();
}

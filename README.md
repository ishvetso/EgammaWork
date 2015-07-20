This repository contains a setup for running electron isolation with CITK (Common Isolation Toolkit) [CITK twiki]. 
It also includes PUPPI-based electron isolation.

<b>IMPORTANT: PHYS14 samples cannot be run in 74X, 73X should be used. To run 74X, use Spring15 samples.</b>

In order to run <b>only electron isolation</b> with CITK (pfIsolationVariables are also there), you should do following (CMSSW verstion is the one I used):

cmsrel CMSSW_7_4_6_patch2 <br />
cd CMSSW_7_4_6_patch2/src <br />
cmsenv  <br />
git clone -b ElectronBranch git@github.com:ishvetso/EgammaWork.git  <br />
scram b -j10  <br />
cmsRun EgammaWork/electron_isolation_CITK.py (for miniAOD)  <br />
cmsRun EgammaWork/electron_isolation_CITK_AOD.py (for AOD)  <br />

In order to run only <b>electron isolation</b> with CITK <b>including PUPPI-based electron isolation</b>, you should following (CMSSW verstion is the one I used): 

1.cmsrel CMSSW_7_4_6_patch2  <br />
  cd CMSSW_7_4_6_patch2/src <br />
  cmsenv <br />
  git clone -b ElectronBranch git@github.com:ishvetso/EgammaWork.git <br />
2. do what is stated on the [PUPPI twiki] <br/>
3. Then: 
  scram b -j10  <br />
  # this is to run electron isolation for puppi with electrons with cone footprint removal <br / >
  cmsRun EgammaWork/electron_isolation_CITK_PUPPI.py (miniAOD) <br />
  # this is  for cross-checks, added isolation with map-based footprint removal (for standard PFCandidates, not for puppi) <br />
  cmsRun EgammaWork/electron_isolation_CITK_PUPPI_working.py <br />

[CITK twiki]:https://twiki.cern.ch/twiki/bin/viewauth/CMS/CommonIDAndIsolationFW
[PUPPI twiki]:https://twiki.cern.ch/twiki/bin/viewauth/CMS/PUPPI#Validation_framework_in_CMSSW_73

This repository contains a setup for running electron isolation with CITK (Common Isolation Toolkit) [CITK twiki]. 
It also includes PUPPI-based electron isolation. 

In order to run *only electron isolation* with CITK (pfIsolationVariables are also there), you should do following (CMSSW verstion is the one I used):

cmsrel CMSSW_7_3_3 <br />
cd CMSSW_7_3_3/src <br />
cmsenv  <br />
git clone -b ElectronBranch git@github.com:ishvetso/EgammaWork.git  <br />
scram b -j10  <br />
cmsRun EgammaWork/electron_isolation_CITK.py (for miniAOD)  <br />
cmsRun EgammaWork/electron_isolation_CITK_AOD.py (for AOD)  <br />

In order to run only *electron isolation* with CITK *including PUPPI-based electron isolation*, you should following (CMSSW verstion is the one I used): 

1.cmsrel CMSSW_7_3_3  <br />
  cd CMSSW_7_3_3/src <br />
  cmsenv <br />
  git clone -b ElectronBranch git@github.com:ishvetso/EgammaWork.git <br />
2. do what is stated on the [PUPPI twiki] 
3. Then: 
  scram b -j10  <br />
  cmsRun EgammaWork/electron_isolation_CITK_PUPPI.py (miniAOD) <br />

[CITK twiki]:https://twiki.cern.ch/twiki/bin/viewauth/CMS/CommonIDAndIsolationFW
[PUPPI twiki]:https://twiki.cern.ch/twiki/bin/viewauth/CMS/PUPPI#Validation_framework_in_CMSSW_73

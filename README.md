This repository contains a setup for running electron isolation with CITK (Common Isolation Toolkit) [CITK twiki]. 
It also includes PUPPI-based electron isolation.

In order to run <b>only electron isolation</b> with CITK (pfIsolationVariables are also there), you should do following (CMSSW verstion is the one I used):

cmsrel CMSSW_7_4_6_patch2 <br />
cd CMSSW_7_4_6_patch2/src <br />
cmsenv  <br />
git clone -b ElectronBranch git@github.com:ishvetso/EgammaWork.git  <br />
scram b -j10  <br />
cmsRun EgammaWork/electron_isolation_CITK_miniAOD.py (for miniAOD)  <br />
cmsRun EgammaWork/electron_isolation_CITK_AOD.py (for AOD)  <br />

In order to run only <b>electron isolation</b> with CITK <b>including PUPPI-based electron isolation</b>, you should do following (CMSSW verstion is the one I used):<br/> 
(<b>CMSSW_7_4_7 should be used or later</b>(mva electrons ids are not available in earlier releases))<br/>
1.cmsrel CMSSW_7_4_7  <br />
  cd CMSSW_7_4_7/src <br />
  cmsenv <br />
  git clone -b ElectronBranch git@github.com:ishvetso/EgammaWork.git <br />
2. Since CMSSW_8_0_20, "v10" of PUPPI is in CMSSW, no action should be taken. However, this might change in the future (check [PUPPI twiki] <br/>
3. Then: 
  scram b -j10  <br />
  # this is to run electron isolation for puppi with electrons with cone footprint removal <br / >
  cmsRun EgammaWork/electron_isolation_CITK_PUPPI_IDs.py (miniAOD) <br />

[CITK twiki]:https://twiki.cern.ch/twiki/bin/viewauth/CMS/CommonIDAndIsolationFW
[PUPPI twiki]:https://twiki.cern.ch/twiki/bin/viewauth/CMS/PUPPI#Validation_framework_in_CMSSW_73

This repository contains a setup for running electron isolation with CITK (Common Isolation Toolkit) [CITK twiki]. 
It also includes PUPPI-based electron isolation. 

In order to run only electron isolation with CITK, you should following (CMSSW verstion is the one I used):

cmsrel CMSSW_7_3_3
cd CMSSW_7_3_3/src
cmsenv
git clone git@github.com:ishvetso/EgammaWork.git
scram b -j10
cmsRun EgammaWork/electron_isolation_CITK.py (for miniAOD)
cmsRun EgammaWork/electron_isolation_CITK_AOD.py (for AOD)
 
In order to run only electron isolation with CITK, you should following (CMSSW verstion is the one I used): 

1. do what is stated on the PUPPI-twiki [2] 
2. Then: 
cmsrel CMSSW_7_3_3
cd CMSSW_7_3_3/src
cmsenv
git clone git@github.com:ishvetso/EgammaWork.git
scram b -j10
cmsRun EgammaWork/electron_isolation_CITK_PUPPI.py (miniAOD)

[CITK]:https://twiki.cern.ch/twiki/bin/viewauth/CMS/CommonIDAndIsolationFW
[PUPPI]:https://twiki.cern.ch/twiki/bin/viewauth/CMS/PUPPI#Validation_framework_in_CMSSW_73

This repository contains a setup for running electron isolation with CITK (Common Isolation Toolkit) [CITK twiki]. 


In order to run <b>only electron isolation</b> with CITK (pfIsolationVariables are also there), you should do following (CMSSW verstion is the one I used):

cmsrel CMSSW_7_6_4 <br />
cd CMSSW_7_6_4/src <br />
cmsenv  <br />
git clone -b CITK_demos_electrons git@github.com:ishvetso/EgammaWork.git  <br />
scram b -  <br />
cmsRun EgammaWork/electron_isolation_CITK_miniAOD.py (for miniAOD)  <br />
cmsRun EgammaWork/electron_isolation_CITK_AOD.py (for AOD)  <br />

[CITK twiki]:https://twiki.cern.ch/twiki/bin/viewauth/CMS/CommonIDAndIsolationFW

This repository contains a setup for running photon isolation with CITK (Common Isolation Toolkit) [CITK twiki]. 


In order to run only photon isolation with CITK (pfIsolationVariables are also there), you should do following (CMSSW verstion is the one I used):



1. cmsrel CMSSW_7_4_7 <br />
cd CMSSW_7_4_7/src <br />
cmsenv <br />
2. git cms-merge-topic ishvetso:PhotonPFIsolationWithMapBasedVeto
3.  git cms-merge-topic ikrav:egm_id_747_v2
4. git clone -b PhotonBranch_v2 git@github.com:ishvetso/EgammaWork.git 
5. scram b -j 20
6. cmsRun EgammaWork/photonIsolation_miniAOD.py


[CITK twiki]:https://twiki.cern.ch/twiki/bin/viewauth/CMS/CommonIDAndIsolationFW

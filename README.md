This repository contains a setup for running photon isolation with CITK (Common Isolation Toolkit) [CITK twiki]. 

In order to run only photon isolation with CITK (pfIsolationVariables are also there), you should do following (CMSSW verstion is the one I used):



1. cmsrel CMSSW_7_4_7 <br />
cd CMSSW_7_4_7/src <br />
cmsenv <br />
2. git cms-merge-topic ikrav:egm_id_74X_v2
3. git cms-merge-topic ishvetso:PhotonPFIsolationWithMapBasedVeto # get module for map-based vetoing
5. git clone -b PhotonBranch git@github.com:ishvetso/EgammaWork.git 
6. scram b -j 20
7. cmsRun EgammaWork/photonIsolation_miniAOD.py


[CITK twiki]:https://twiki.cern.ch/twiki/bin/viewauth/CMS/CommonIDAndIsolationFW


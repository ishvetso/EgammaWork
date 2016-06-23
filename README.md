<b> This branch contains a setup for running photon isolation with CITK (Common Isolation Toolkit):</b> [CITK twiki]. 

In order to run only photon isolation with CITK, you should do following (CMSSW version is the one recipe is produced with):

1. cmsrel CMSSW_7_6_4 <br />
cd CMSSW_7_6_4/src <br />
cmsenv <br />
2. git clone -b  photon_isolations  git@github.com:ishvetso/EgammaWork.git 
3. scram b 
4. cmsRun EgammaWork/photonIsolation_miniAOD.py # miniAOD case <br/>
   cmsRun EgammaWork/photonIsolation_AOD.py # AOD case


[CITK twiki]:https://twiki.cern.ch/twiki/bin/viewauth/CMS/CommonIDAndIsolationFW


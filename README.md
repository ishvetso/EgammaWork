<b> This branch contains a setup for running photon isolation with CITK (Common Isolation Toolkit):</b> [CITK twiki]. 

In order to run only photon isolation with CITK, you should do following (CMSSW version is the one recipe is produced with):

1. cmsrel CMSSW_7_4_7 <br />
cd CMSSW_7_4_7/src <br />
cmsenv <br />
2. git cms-merge-topic ishvetso:PhotonIsolationCITK_747_recipe # get CITK module for photons
3. git clone -b PhotonIsolationCITK git@github.com:ishvetso/EgammaWork.git 
4. scram b -j 10
5. cmsRun EgammaWork/photonIsolation_miniAOD.py # miniAOD case <br/>
   cmsRun EgammaWork/photonIsolation_AOD.py # AOD case


[CITK twiki]:https://twiki.cern.ch/twiki/bin/viewauth/CMS/CommonIDAndIsolationFW


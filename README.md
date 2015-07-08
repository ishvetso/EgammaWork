This repository contains a setup for running photon isolation with CITK (Common Isolation Toolkit) [CITK twiki]. 
It also includes PUPPI-based photon isolation. 

In order to run only photon isolation with CITK (pfIsolationVariables are also there), you should do following (CMSSW verstion is the one I used):

cmsrel CMSSW_7_4_6_patch2 <br />
cd CMSSW_7_4_6_patch2/src <br />

1. do what is stated on the [PUPPI twiki]; <br />
for example: <br />
git cms-addpkg CommonTools/PileupAlgos <br />
git cms-merge-topic nhanvtran:puppi-etadep-746p2-v8 (<b>this might change!</b>) <br />
2. git cms-merge-topic ikrav:egm_id_74X_v2
3. git cms-merge-topic ishvetso:PhotonPFIsolationWithMapBasedVeto # get module for map-based vetoing
4. git cms-merge-topic ishvetso:CITK_PUPPI # module for PUPPI to reweight PFCandidates on the fly using ValueMaps from puppi
5. git clone -b PhotonBranch git@github.com:ishvetso/EgammaWork.git 
6. scram b -j 20


[CITK twiki]:https://twiki.cern.ch/twiki/bin/viewauth/CMS/CommonIDAndIsolationFW
[PUPPI twiki]:https://twiki.cern.ch/twiki/bin/viewauth/CMS/PUPPI#Validation_framework_in_CMSSW_73

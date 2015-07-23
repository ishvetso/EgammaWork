// -*- C++ -*-
//
// Package:    EgammaWork/PFCandidateAnalyzer
// Class:      PFCandidateAnalyzer
// 
/**\class PFCandidateAnalyzer PFCandidateAnalyzer.cc EgammaWork/PFCandidateAnalyzer/plugins/PFCandidateAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Ivan Shvetsov
//         Created:  Thu, 23 Jul 2015 08:21:58 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "TTree.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/ValueMap.h"     
//
// class declaration
//

class PFCandidateAnalyzer : public edm::EDAnalyzer {
   public:
      explicit PFCandidateAnalyzer(const edm::ParameterSet&);
      ~PFCandidateAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      edm::EDGetTokenT<edm::View<reco::Candidate> > leptons_;
      edm::EDGetTokenT<edm::View<pat::PackedCandidate> > cands_;
      edm::EDGetTokenT<edm::ValueMap<float> > puppiValueMapToken_;//for puppiValueMap

      TTree * tree;

      Float_t deltaR_;
      Float_t puppiWeight;
      Float_t pdgId_;
      Float_t charge_;
      Int_t fromPV_;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
PFCandidateAnalyzer::PFCandidateAnalyzer(const edm::ParameterSet& iConfig):
  leptons_(consumes<edm::View<reco::Candidate> > (iConfig.getParameter<edm::InputTag>( "leptons" ) ) ),
  cands_(consumes<edm::View<pat::PackedCandidate> > (iConfig.getParameter<edm::InputTag>( "cand_src" ) ) ),
  puppiValueMapToken_ (consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("puppiValueMap")) )

{
   //now do what ever initialization is needed
  edm::Service<TFileService> fs;
  tree = fs->make<TTree> ("PFCandidateAnalyzerTree", "puppiWeights");
  tree -> Branch("deltaR" , &deltaR_ , "deltaR/F");
  tree -> Branch("puppiWeight" , &puppiWeight , "puppiWeight/F");
  tree -> Branch("pdgId" , &pdgId_ , "pdgId/F");
  tree -> Branch("charge" , &charge_ , "charge/F");
   tree -> Branch("fromPV" , &fromPV_ , "fromPV/I");

}


PFCandidateAnalyzer::~PFCandidateAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
PFCandidateAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  Handle<edm::View<reco::Candidate> > leptons;
  Handle<edm::View<pat::PackedCandidate> > cands;
  Handle<edm::ValueMap<float> > puppiValueMap;

  iEvent.getByToken(leptons_, leptons);
  iEvent.getByToken(cands_, cands);
  iEvent.getByToken(puppiValueMapToken_, puppiValueMap);

  for (unsigned int iLepton = 0; iLepton <  leptons -> size() ; iLepton ++)
  {
    for(unsigned int iCand = 0; iCand < cands -> size(); iCand ++)
    {
      puppiWeight = (*puppiValueMap)[cands -> ptrAt(iCand)];
      deltaR_ = deltaR(leptons -> at(iLepton), cands -> at(iCand));
      pdgId_ = (cands -> at(iCand)).pdgId();
      charge_ = (cands -> at(iCand)).charge();
      fromPV_ = (cands -> at(iCand)).fromPV();
      tree -> Fill();
    } 
  }


  

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
PFCandidateAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PFCandidateAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
PFCandidateAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
PFCandidateAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
PFCandidateAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
PFCandidateAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PFCandidateAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please §e this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PFCandidateAnalyzer);

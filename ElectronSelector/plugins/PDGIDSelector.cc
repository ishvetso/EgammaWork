// -*- C++ -*-
//
// Package:    EgammaWork/PDGIDSelector
// Class:      PDGIDSelector
// 
/**\class PDGIDSelector PDGIDSelector.cc EgammaWork/PDGIDSelector/plugins/PDGIDSelector.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Ivan Shvetsov
//         Created:  Mon, 13 Jul 2015 13:27:48 GMT
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
#include "DataFormats/Common/interface/Ref.h"


//
// class declaration
//

class PDGIDSelector : public edm::EDProducer {
   public:
      explicit PDGIDSelector(const edm::ParameterSet&);
      ~PDGIDSelector();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      edm::EDGetTokenT<std::vector<reco::PFCandidate>> pfCollectionToken;
      std::vector<int> PDGIDs;
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
PDGIDSelector::PDGIDSelector(const edm::ParameterSet& iConfig):
  pfCollectionToken(consumes<std::vector<reco::PFCandidate>> (iConfig.getParameter<edm::InputTag>( "src" ) ) )  
 
{
   //now do what ever initialization is needed
  PDGIDs= iConfig.getParameter<std::vector<int>>("PDGIDs");
   produces<std::vector<reco::PFCandidate>>();

}


PDGIDSelector::~PDGIDSelector()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
void
PDGIDSelector::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
  Handle<std::vector<reco::PFCandidate> > PF;
  iEvent.getByToken(pfCollectionToken, PF);
  std::auto_ptr<std::vector<reco::PFCandidate> > CandidatesSelected(new std::vector<reco::PFCandidate>);

  for (unsigned int iCand = 0; iCand < PF -> size(); iCand ++)
   {
       for (unsigned int iPDGID = 0; iPDGID < PDGIDs.size(); iPDGID++){
        if ( PF -> at(iCand).pdgId() == PDGIDs.at(iPDGID)) CandidatesSelected -> push_back(PF-> at(iCand));

       }
      
   }

   
   
   iEvent.put(CandidatesSelected);
   


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
PDGIDSelector::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PDGIDSelector::endJob() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
PDGIDSelector::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
PDGIDSelector::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
PDGIDSelector::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
PDGIDSelector::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PDGIDSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(PDGIDSelector);

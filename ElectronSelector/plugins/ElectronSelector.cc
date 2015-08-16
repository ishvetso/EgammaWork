// -*- C++ -*-
//
// Package:    EgammaWork/ElectronSelector
// Class:      ElectronSelector
// 
/**\class ElectronSelector ElectronSelector.cc EgammaWork/ElectronSelector/plugins/ElectronSelector.cc

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

class ElectronSelector : public edm::EDFilter {
   public:
      explicit ElectronSelector(const edm::ParameterSet&);
      ~ElectronSelector();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      edm::EDGetTokenT<edm::View<pat::Electron> > Electrons_;
      edm::EDGetTokenT<edm::View<pat::PackedCandidate> > cands_;
      
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
ElectronSelector::ElectronSelector(const edm::ParameterSet& iConfig):
  Electrons_(consumes<edm::View<pat::Electron> > (iConfig.getParameter<edm::InputTag>( "electron_src" ) ) ),
  cands_(consumes<edm::View<pat::PackedCandidate> > (iConfig.getParameter<edm::InputTag>( "cand_src" ) ) )
{
   //now do what ever initialization is needed
   produces<edm::RefVector<std::vector<pat::Electron>>>();

}


ElectronSelector::~ElectronSelector()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
ElectronSelector::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
  Handle<edm::View<pat::Electron> > Electrons;
  Handle<edm::View<pat::PackedCandidate> > cands;

  bool pass = false;

    std::auto_ptr<edm::RefVector<std::vector<pat::Electron> > > ElectronsSelected(new edm::RefVector<std::vector<pat::Electron>>);

    iEvent.getByToken(Electrons_, Electrons);
     iEvent.getByToken(cands_, cands);

  for (unsigned int iElectron = 0; iElectron < Electrons -> size(); iElectron ++)
   {
       pass = false;
       reco::GsfTrackRef trackref_ele = (Electrons -> at(iElectron) ).gsfTrack() ;
       float eta_ele = trackref_ele -> eta();
       float phi_ele = trackref_ele -> phi();
      
       for (unsigned iCand = 0; iCand < cands -> size(); iCand ++)
       {
          reco::Track cand_track = (cands -> at(iCand)).pseudoTrack();
          float eta_cand = cand_track.eta();
          float phi_cand = cand_track.phi();
          if (  fabs((cands -> at(iCand).pdgId()) == 11 ) ) {
            float deltaR_ = deltaR(eta_ele, phi_ele, eta_cand, phi_cand); 
            if (deltaR_ < 0.01 ) pass = true;
           }
       }
       if (pass) ElectronsSelected -> push_back( Electrons -> refAt(iElectron).castTo<edm::Ref<std::vector<pat::Electron> > >());
   }

   bool passEvent = false;
   if (ElectronsSelected -> size() > 0 ) passEvent = true;
   iEvent.put(ElectronsSelected);
    std::cout << passEvent << std::endl;


#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
   return passEvent;
}

// ------------ method called once each job just before starting event loop  ------------
void 
ElectronSelector::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ElectronSelector::endJob() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
ElectronSelector::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
ElectronSelector::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
ElectronSelector::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
ElectronSelector::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ElectronSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(ElectronSelector);

// -*- C++ -*-
//
// Package:    EgammaWork/MyPdgIDPackedCandidateSelector
// Class:      MyPdgIDPackedCandidateSelector
// 
/**\class MyPdgIDPackedCandidateSelector MyPdgIDPackedCandidateSelector.cc EgammaWork/MyPdgIDPackedCandidateSelector/plugins/MyPdgIDPackedCandidateSelector.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Ivan Shvetsov
//         Created:  Wed, 15 Apr 2015 11:57:42 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

//
// class declaration
//

class MyPdgIDPackedCandidateSelector : public edm::EDProducer {
   public:
      explicit MyPdgIDPackedCandidateSelector(const edm::ParameterSet&);
      ~MyPdgIDPackedCandidateSelector();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      
      edm::EDGetTokenT<edm::View<pat::PackedCandidate> > candidates_;
      std::vector<int> pdgIDs;
      
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
MyPdgIDPackedCandidateSelector::MyPdgIDPackedCandidateSelector(const edm::ParameterSet& iConfig):
  candidates_(consumes<edm::View<pat::PackedCandidate> > (iConfig.getParameter<edm::InputTag>( "src" ) ) ),
  pdgIDs(iConfig.getParameter<std::vector<int>>( "pdgId" )) 
{
   //register your products
   produces<std::vector<pat::PackedCandidate>>();
/* Examples
   produces<ExampleData2>();

   //if do put with a label
   produces<ExampleData2>("label");
 
   //if you want to put into the Run
   produces<ExampleData2,InRun>();
*/
   //now do what ever other initialization is needed
  
}


MyPdgIDPackedCandidateSelector::~MyPdgIDPackedCandidateSelector()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
MyPdgIDPackedCandidateSelector::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
  
   Handle <edm::View <pat::PackedCandidate> > candidates;
   std::auto_ptr<std::vector<pat::PackedCandidate> > OutCollection(new std::vector<pat::PackedCandidate>);
     
   iEvent.getByToken(candidates_, candidates);
   
      
   for (unsigned int iCand = 0; iCand < candidates -> size(); iCand ++ )
   {
    for (unsigned iPDG = 0; iPDG < pdgIDs.size(); iPDG ++)
    {
      if ((candidates -> at(iCand)).pdgId() == pdgIDs.at(iPDG)) OutCollection -> push_back(candidates -> at(iCand));
     }  
   }
   iEvent.put(OutCollection);
/* This is an event example
   //Read 'ExampleData' from the Event
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);

   //Use the ExampleData to create an ExampleData2 which 
   // is put into the Event
   std::unique_ptr<ExampleData2> pOut(new ExampleData2(*pIn));
   iEvent.put(std::move(pOut));
*/

/* this is an EventSetup example
   //Read SetupData from the SetupRecord in the EventSetup
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
*/
 
}

// ------------ method called once each job just before starting event loop  ------------
void 
MyPdgIDPackedCandidateSelector::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MyPdgIDPackedCandidateSelector::endJob() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
MyPdgIDPackedCandidateSelector::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
MyPdgIDPackedCandidateSelector::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
MyPdgIDPackedCandidateSelector::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
MyPdgIDPackedCandidateSelector::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MyPdgIDPackedCandidateSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MyPdgIDPackedCandidateSelector);

// -*- C++ -*-
//
// Package:    ElectronWork/ElectronNtupler
// Class:      PhotonValidator
// 
/**\class PhotonValidator PhotonValidator.cc ElectronWork/PhotonValidator/plugins/PhotonValidator.cc
 Description: [one line class summary]
 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Ilya Kravchenko
//         Created:  Thu, 10 Jul 2014 09:54:13 GMT
//
//


// system include files
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"

#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"
#include "Math/VectorUtil.h"



//
// class declaration
//

class PhotonValidator : public edm::EDAnalyzer {
   public:
      explicit PhotonValidator(const edm::ParameterSet&);
       ~PhotonValidator();
       static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      edm::EDGetTokenT<edm::View<pat::Photon> > photonCollectionToken_;
      edm::EDGetTokenT<edm::ValueMap<float> > phoChargedIsolationToken_; 
      edm::EDGetTokenT<edm::ValueMap<float> > phoNeutralHadronIsolationToken_; 
      edm::EDGetTokenT<edm::ValueMap<float> > phoPhotonIsolationToken_; 
      
      //CITK
      edm::EDGetTokenT<edm::ValueMap<float> > phoChargedIsolation_CITK_Token_;
      edm::EDGetTokenT<edm::ValueMap<float> > phoNeutralHadronIsolation_CITK_Token_;
      edm::EDGetTokenT<edm::ValueMap<float> > phoPhotonIsolation_CITK_Token_;

      TTree *photonTree_;

      Float_t isoChargedHadrons_;
      Float_t isoNeutralHadrons_;
      Float_t isoPhotons_;
      
      Float_t isoChargedHadrons_CITK_;
      Float_t isoNeutralHadrons_CITK_;
      Float_t isoPhotons_CITK_;


};



//
// static data member definitions
//

//
// constructors and destructor
//
PhotonValidator::PhotonValidator(const edm::ParameterSet& iConfig):

  photonCollectionToken_(consumes<edm::View<pat::Photon> >(iConfig.getParameter<edm::InputTag>("photons"))),
  phoChargedIsolationToken_(consumes <edm::ValueMap<float> >
			    (iConfig.getParameter<edm::InputTag>("phoChargedIsolation"))),
  phoNeutralHadronIsolationToken_(consumes <edm::ValueMap<float> >
				  (iConfig.getParameter<edm::InputTag>("phoNeutralHadronIsolation"))),
  phoPhotonIsolationToken_(consumes <edm::ValueMap<float> >
			   (iConfig.getParameter<edm::InputTag>("phoPhotonIsolation"))),
  //CITK
  phoChargedIsolation_CITK_Token_(consumes <edm::ValueMap<float> >
			    (iConfig.getParameter<edm::InputTag>("phoChargedIsolation_CITK"))),
  phoNeutralHadronIsolation_CITK_Token_(consumes <edm::ValueMap<float> >
				  (iConfig.getParameter<edm::InputTag>("phoNeutralHadronIsolation_CITK"))),
  phoPhotonIsolation_CITK_Token_(consumes <edm::ValueMap<float> >
			   (iConfig.getParameter<edm::InputTag>("phoPhotonIsolation_CITK")))
{

  edm::Service<TFileService> fs;
  photonTree_ = fs->make<TTree> ("PhotonTree", "Photon data");
  
  photonTree_->Branch("isoChargedHadrons"      , &isoChargedHadrons_, "isoChargedHadrons/F");
  photonTree_->Branch("isoNeutralHadrons"      , &isoNeutralHadrons_, "isoNeutralHadrons/F");
  photonTree_->Branch("isoPhotons"             , &isoPhotons_, "isoPhotons/F");
  
  //CITK
  photonTree_->Branch("isoChargedHadrons_CITK"      , &isoChargedHadrons_CITK_, "isoChargedHadrons_CITK/F");
  photonTree_->Branch("isoNeutralHadrons_CITK"      , &isoNeutralHadrons_CITK_, "isoNeutralHadrons_CITK/F");
  photonTree_->Branch("isoPhotons_CITK"             , &isoPhotons_CITK_, "isoPhotons_CITK/F");

 
}


PhotonValidator::~PhotonValidator()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
PhotonValidator::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace std;
  using namespace edm;
  // using namespace reco;
  
  // Get photon collection
  edm::Handle<edm::View<pat::Photon> > collection;
  iEvent.getByToken(photonCollectionToken_, collection);


  // Get the isolation maps
  edm::Handle<edm::ValueMap<float> > phoChargedIsolationMap;
  iEvent.getByToken(phoChargedIsolationToken_, phoChargedIsolationMap);
  edm::Handle<edm::ValueMap<float> > phoNeutralHadronIsolationMap;
  iEvent.getByToken(phoNeutralHadronIsolationToken_, phoNeutralHadronIsolationMap);
  edm::Handle<edm::ValueMap<float> > phoPhotonIsolationMap;
  iEvent.getByToken(phoPhotonIsolationToken_, phoPhotonIsolationMap);
  
  // Get the isolation maps CITK
  edm::Handle<edm::ValueMap<float> > phoChargedIsolationMap_CITK;
  iEvent.getByToken(phoChargedIsolation_CITK_Token_, phoChargedIsolationMap_CITK);
  edm::Handle<edm::ValueMap<float> > phoNeutralHadronIsolationMap_CITK;
  iEvent.getByToken(phoNeutralHadronIsolation_CITK_Token_, phoNeutralHadronIsolationMap_CITK);
  edm::Handle<edm::ValueMap<float> > phoPhotonIsolationMap_CITK;
  iEvent.getByToken(phoPhotonIsolation_CITK_Token_, phoPhotonIsolationMap_CITK);

  for( View<pat::Photon>::const_iterator pho = collection->begin();
       pho != collection->end(); pho++){
    
    // Kinematics (nobody uses photons below 15 GeV, and they are not stored in miniAOD anyways)
    if( pho->pt() < 15 ) 
      continue;
    
    const edm::Ptr<pat::Photon> phoPtr( collection, pho - collection->begin() );
   
    isoChargedHadrons_ = (*phoChargedIsolationMap)[phoPtr] ;
    isoNeutralHadrons_ =  (*phoNeutralHadronIsolationMap)[phoPtr] ;
    isoPhotons_        = (*phoPhotonIsolationMap)[phoPtr] ;
    
    //CITK
    isoChargedHadrons_CITK_ = (*phoChargedIsolationMap_CITK)[phoPtr] ;
    isoNeutralHadrons_CITK_ =  (*phoNeutralHadronIsolationMap_CITK)[phoPtr] ;
    isoPhotons_CITK_        = (*phoPhotonIsolationMap_CITK)[phoPtr] ;

    // Save the info
    photonTree_->Fill();
   }
   



}


// ------------ method called once each job just before starting event loop  ------------
void 
PhotonValidator::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PhotonValidator::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
PhotonValidator::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
PhotonValidator::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
PhotonValidator::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
PhotonValidator::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------

void
PhotonValidator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PhotonValidator);
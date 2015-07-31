// -*- C++ -*-
//
// Package:    ElectronWork/MuonNtupler
// Class:      MuonNtupler
// 
/**\class MuonNtupler MuonNtupler.cc ElectronWork/MuonNtupler/plugins/MuonNtupler.cc

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

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/Muon.h"

#include "TTree.h"
#include "Math/VectorUtil.h"



//
// class declaration
//

class MuonNtupler : public edm::EDAnalyzer {
public:
  explicit MuonNtupler(const edm::ParameterSet&);
  ~MuonNtupler();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  

  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  

  
  // ----------member data ---------------------------
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  edm::EDGetTokenT<edm::View<reco::Muon>> muonToken_;
 
  
  //CITK
  edm::EDGetTokenT<edm::ValueMap<float> > ValueMaps_ChargedHadrons_;
  edm::EDGetTokenT<edm::ValueMap<float> > ValueMaps_NeutralHadrons_;
  edm::EDGetTokenT<edm::ValueMap<float> > ValueMaps_Photons_;
  
  TTree *muonTree_;
  
  //event info
  int nevent, run, lumi;
  
  Int_t nPV_;        // number of reconsrtucted primary vertices
 
  // all electron variables
  Float_t pt_;
  Float_t eta_;
 
  
  Float_t isoChargedHadrons_;
  Float_t isoNeutralHadrons_;
  Float_t isoPhotons_;
  Float_t isoChargedFromPU_;
  
  Float_t relisoChargedHadrons_;
  Float_t relisoNeutralHadrons_;
  Float_t relisoPhotons_;

  Float_t relIsoWithDBeta_;
  
  //CITK
  Float_t sumChargedHadronPt_CITK;
  Float_t sumNeutralHadronPt_CITK;
  Float_t sumPhotonPt_CITK;
  

  
  Float_t relisoChargedHadronPt_CITK;
  Float_t relisoNeutralHadronPt_CITK;
  Float_t relisoPhotonPt_CITK;
  
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
MuonNtupler::MuonNtupler(const edm::ParameterSet& iConfig):
  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  
  muonToken_(consumes<edm::View<reco::Muon> >(iConfig.getParameter<edm::InputTag>("muons"))),
  //CITK
  ValueMaps_ChargedHadrons_(consumes<edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>( "ValueMaps_ChargedHadrons_src" ) ) ),
  ValueMaps_NeutralHadrons_(consumes<edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>( "ValueMaps_NeutralHadrons_src" ) ) ),
  ValueMaps_Photons_(consumes<edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>( "ValueMaps_Photons_src" ) ) )

{

  edm::Service<TFileService> fs;
  muonTree_ = fs->make<TTree> ("ElectronTree", "Electron data");
  
  //event info
  muonTree_->Branch("event",	      &nevent,    	  "event/I"           );
  muonTree_->Branch("lumi", 	      &lumi,   		  "lumi/I"  		);
  muonTree_->Branch("run",	      &run,		  "run/I"  	       );
  
  muonTree_->Branch("nPV"        ,  &nPV_     , "nPV/I");
  
  muonTree_->Branch("pt"    ,  &pt_    , "pt/F");			    
  muonTree_->Branch("eta" ,  &eta_ , "eta/F");
  
  muonTree_->Branch("isoChargedFromPU"       , &isoChargedFromPU_);
  //muonTree_->Branch("relIsoWithEA"           , &relIsoWithEA_, "relIsoWithEA/F");
  muonTree_->Branch("relIsoWithDBeta"      , &relIsoWithDBeta_, "relIsoWithDBeta/F");
  
  
 
  muonTree_->Branch("isoChargedHadrons"      , &isoChargedHadrons_ , "isoChargedHadrons/F");
  muonTree_->Branch("isoNeutralHadrons"      , &isoNeutralHadrons_, "isoNeutralHadrons/F");
  muonTree_->Branch("isoPhotons"             , &isoPhotons_, "isoPhotons/F");
  
  muonTree_->Branch("relisoChargedHadrons"      , &relisoChargedHadrons_ , "relisoChargedHadrons/F");
  muonTree_->Branch("relisoNeutralHadrons"      , &relisoNeutralHadrons_, "relisoNeutralHadrons/F");
  muonTree_->Branch("relisoPhotons"             , &relisoPhotons_, "relisoPhotons/F");
  
  //CITK
  muonTree_ -> Branch("sumChargedHadronPt_CITK", &sumChargedHadronPt_CITK, "sumChargedHadronPt_CITK/F");
  muonTree_ -> Branch("sumNeutralHadronPt_CITK", &sumNeutralHadronPt_CITK, "sumNeutralHadronPt_CITK/F");
  muonTree_ -> Branch("sumPhotonPt_CITK", &sumPhotonPt_CITK, "sumPhotonPt_CITK/F");
    
  muonTree_ -> Branch("relisoChargedHadronPt_CITK", &relisoChargedHadronPt_CITK, "relisoChargedHadronPt_CITK/F");
  muonTree_ -> Branch("relisoNeutralHadronPt_CITK", &relisoNeutralHadronPt_CITK, "relisoNeutralHadronPt_CITK/F");
  muonTree_ -> Branch("relisoPhotonPt_CITK", &relisoPhotonPt_CITK, "relisoPhotonPt_CITK/F");
  


}


MuonNtupler::~MuonNtupler()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
MuonNtupler::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace std;
  using namespace edm;
  using namespace reco;
  
  //event info
  nevent = iEvent.eventAuxiliary().event();
  run    = iEvent.eventAuxiliary().run();
  lumi   = iEvent.eventAuxiliary().luminosityBlock();
  
 

  // Get PV
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);
  if (vertices->empty()) return; // skip the event if no PV found
  //const reco::Vertex &pv = vertices->front();
  
  nPV_ = vertices -> size();

  
  // Get electron collection
  Handle<edm::View<reco::Muon> > muons;
  iEvent.getByToken(muonToken_, muons);
    
  //CITK
  Handle <edm::ValueMap <float> > ValueMaps_ChargedHadrons, ValueMaps_NeutralHadrons, ValueMaps_Photons;
  Handle <edm::ValueMap <float> > ValueMaps_PUPPI_ChargedHadrons, ValueMaps_PUPPI_NeutralHadrons, ValueMaps_PUPPI_Photons;
  Handle <edm::ValueMap <float> > ValueMaps_PUPPI_NoLeptons_ChargedHadrons, ValueMaps_PUPPI_NoLeptons_NeutralHadrons, ValueMaps_PUPPI_NoLeptons_Photons;
  
  //CITK
  iEvent.getByToken( ValueMaps_ChargedHadrons_ , ValueMaps_ChargedHadrons);
  iEvent.getByToken( ValueMaps_NeutralHadrons_ , ValueMaps_NeutralHadrons);
  iEvent.getByToken( ValueMaps_Photons_ , ValueMaps_Photons);


  //
  // Loop over muons
  //
  // printf("DEBUG: new event\n"); 
  for (unsigned int iMuon = 0; iMuon < muons -> size(); iMuon++) {
    
    
    auto muonPtr = muons -> ptrAt(iMuon);
   
    // Kinematics
    pt_ = muonPtr -> pt();
    
    
    
    eta_ = muonPtr ->eta();
    
        
    // Isolation
    reco::MuonPFIsolation pfIso = muonPtr -> pfIsolationR03();
    isoChargedHadrons_ = pfIso.sumChargedHadronPt ;
    isoNeutralHadrons_ = pfIso.sumNeutralHadronEt ;
    isoPhotons_        = pfIso.sumPhotonEt ;
    isoChargedFromPU_  = pfIso.sumPUPt ;
    
    relisoChargedHadrons_ = pfIso.sumChargedHadronPt/pt_;
    relisoNeutralHadrons_ = pfIso.sumNeutralHadronEt/pt_;
    relisoPhotons_        = pfIso.sumPhotonEt/pt_;
   
    
    // Compute isolation with delta beta correction for PU
    float absiso = pfIso.sumChargedHadronPt + max(0.0 , pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - 0.5 * pfIso.sumPUPt );
    relIsoWithDBeta_ = absiso/pt_;
    
    //CITK
    sumChargedHadronPt_CITK =  (*ValueMaps_ChargedHadrons)[muonPtr];
    sumNeutralHadronPt_CITK =  (*ValueMaps_NeutralHadrons)[muonPtr];
    sumPhotonPt_CITK        =  (*ValueMaps_Photons)[muonPtr];
  
    relisoChargedHadronPt_CITK = sumChargedHadronPt_CITK/pt_;
    relisoNeutralHadronPt_CITK = sumNeutralHadronPt_CITK/pt_;
    relisoPhotonPt_CITK = sumPhotonPt_CITK/pt_;
    
         
    // Save this electron's info
    muonTree_->Fill();
  }
  
}


// ------------ method called once each job just before starting event loop  ------------
void 
MuonNtupler::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuonNtupler::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
MuonNtupler::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
MuonNtupler::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
MuonNtupler::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
MuonNtupler::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuonNtupler::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  
}


//define this as a plug-in
DEFINE_FWK_MODULE(MuonNtupler);

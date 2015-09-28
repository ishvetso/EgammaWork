// -*- C++ -*-
//
// Package:    EgammaWork/PhotonNtupler
// Class:      SimplePhotonNtupler
// 
/**\class SimplePhotonNtupler SimplePhotonNtupler.cc EgammaWork/PhotonNtupler/plugins/SimplePhotonNtupler.cc
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

#include "DataFormats/EgammaCandidates/interface/Photon.h"

#include "DataFormats/PatCandidates/interface/Photon.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "TTree.h"
#include "Math/VectorUtil.h"

//
// class declaration
//

class SimplePhotonNtupler : public edm::EDAnalyzer {
 public:
  explicit SimplePhotonNtupler(const edm::ParameterSet&);
  ~SimplePhotonNtupler();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  enum PhotonMatchType {UNMATCHED = 0, 
			MATCHED_FROM_GUDSCB,
			MATCHED_FROM_PI0,
			MATCHED_FROM_OTHER_SOURCES};
  
 private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  
  int matchToTruth(const reco::Photon &pho, 
		   const edm::Handle<edm::View<reco::GenParticle>>  &genParticles);
  
  void findFirstNonPhotonMother(const reco::Candidate *particle,
				int &ancestorPID, int &ancestorStatus);
  
  // ----------member data ---------------------------

  // Format-independent data members
  edm::EDGetTokenT<double> rhoToken_;
  
  // AOD case data members
  edm::EDGetToken photonsToken_;
  edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesToken_;
  
  // MiniAOD case data members
  edm::EDGetToken photonsMiniAODToken_;
  edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesMiniAODToken_;

  //CITK
  edm::EDGetTokenT<edm::ValueMap<float> > phoChargedIsolationToken_CITK; 
  edm::EDGetTokenT<edm::ValueMap<float> > phoNeutralHadronIsolationToken_CITK; 
  edm::EDGetTokenT<edm::ValueMap<float> > phoPhotonIsolationToken_CITK; 
 

  TTree *photonTree_;
  Float_t rho_;      // the rho variable
  
  // all photon variables
  // CAUTION: whatever vectors are declared here, make sure you clear them before the loop over photons!
  Int_t nPhotons_;

  std::vector<Float_t> pt_;
  std::vector<Float_t> eta_;
  std::vector<Float_t> phi_;

  // Variables typically used for cut based photon ID
  std::vector<Float_t> hOverE_;
  std::vector<Int_t> hasPixelSeed_;

  std::vector<Float_t> isoChargedHadrons_CITK_;
  std::vector<Float_t> isoNeutralHadrons_CITK_;
  std::vector<Float_t> isoPhotons_CITK_;

  std::vector<Float_t> isoChargedHadrons_pf_;
  std::vector<Float_t> isoNeutralHadrons_pf_;
  std::vector<Float_t> isoPhotons_pf_;

  //relative isolation from CITK with map based veto
  std::vector<Float_t> relisoWithEA_CITK_;
  //relative isolation for pf
  std::vector<Float_t> relisoWithEA_pf_;

  std::vector<Float_t> r9;

  std::vector<Int_t> isTrue_;

  std::vector<int> genWeight;

  // Effective area constants for all isolation types
  EffectiveAreas effAreaChHadrons_;
  EffectiveAreas effAreaNeuHadrons_;
  EffectiveAreas effAreaPhotons_;

   edm::EDGetTokenT<GenEventInfoProduct> genInfoToken;

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
SimplePhotonNtupler::SimplePhotonNtupler(const edm::ParameterSet& iConfig):
  rhoToken_(consumes<double> (iConfig.getParameter<edm::InputTag>("rho"))),
  // Isolations from CITK
  phoChargedIsolationToken_CITK(consumes <edm::ValueMap<float> >
			    (iConfig.getParameter<edm::InputTag>("phoChargedIsolation_CITK"))),
  phoNeutralHadronIsolationToken_CITK(consumes <edm::ValueMap<float> >
				  (iConfig.getParameter<edm::InputTag>("phoNeutralHadronIsolation_CITK"))),
  phoPhotonIsolationToken_CITK(consumes <edm::ValueMap<float> >
			   (iConfig.getParameter<edm::InputTag>("phoPhotonIsolation_CITK"))),

  // Objects containing effective area constants
  effAreaChHadrons_( (iConfig.getParameter<edm::FileInPath>("effAreaChHadFile")).fullPath() ),
  effAreaNeuHadrons_( (iConfig.getParameter<edm::FileInPath>("effAreaNeuHadFile")).fullPath() ),
  effAreaPhotons_( (iConfig.getParameter<edm::FileInPath>("effAreaPhoFile")).fullPath() ),
  genInfoToken(consumes<GenEventInfoProduct> (iConfig.getParameter<edm::InputTag>( "genInfo" ) ) )
{

  //
  // Prepare tokens for all input collections and objects
  //
  
  // AOD tokens
  photonsToken_ = mayConsume<edm::View<reco::Photon> >
    (iConfig.getParameter<edm::InputTag>
     ("photons"));
  
  genParticlesToken_ = mayConsume<edm::View<reco::GenParticle> >
    (iConfig.getParameter<edm::InputTag>
     ("genParticles"));
    
  // MiniAOD tokens
  photonsMiniAODToken_ = mayConsume<edm::View<reco::Photon> >
    (iConfig.getParameter<edm::InputTag>
     ("photonsMiniAOD"));

  genParticlesMiniAODToken_ = mayConsume<edm::View<reco::GenParticle> >
    (iConfig.getParameter<edm::InputTag>
     ("genParticlesMiniAOD"));
 
  edm::Service<TFileService> fs;
  photonTree_ = fs->make<TTree> ("PhotonTree", "Photon data");
  
  photonTree_->Branch("rho"        ,  &rho_ , "rho/F");
  photonTree_->Branch("nPho",  &nPhotons_ , "nPho/I");

  // Kinematics
  photonTree_->Branch("pt"  ,  &pt_    );
  photonTree_->Branch("eta" ,  &eta_ );
  photonTree_->Branch("phi" ,  &phi_ );

  // Variables typically used for cut based photon ID
  photonTree_->Branch("hOverE"                 ,  &hOverE_);
  photonTree_->Branch("hasPixelSeed"           ,  &hasPixelSeed_);

  //CITK
  photonTree_->Branch("isoChargedHadrons_CITK"      , &isoChargedHadrons_CITK_);
  photonTree_->Branch("isoNeutralHadrons_CITK"      , &isoNeutralHadrons_CITK_);
  photonTree_->Branch("isoPhotons_CITK"             , &isoPhotons_CITK_);
  
  //pfIsolation variables
  photonTree_->Branch("isoChargedHadrons_pf"      , &isoChargedHadrons_pf_);
  photonTree_->Branch("isoNeutralHadrons_pf"      , &isoNeutralHadrons_pf_);
  photonTree_->Branch("isoPhotons_pf"             , &isoPhotons_pf_);

  //relative isolation
  photonTree_->Branch("relisoWithEA_CITK"                 , &relisoWithEA_CITK_);
  photonTree_->Branch("relisoWithEA_pf"                 , &relisoWithEA_pf_);

  photonTree_->Branch("r9"                 , &r9);

  photonTree_->Branch("isTrue"             , &isTrue_);
  photonTree_ -> Branch("genWeight", &genWeight);
}


SimplePhotonNtupler::~SimplePhotonNtupler()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
SimplePhotonNtupler::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace std;
  using namespace edm;
  using namespace reco;

  // Retrieve the collection of photons from the event.
  // If we fail to retrieve the collection with the standard AOD
  // name, we next look for the one with the stndard miniAOD name. 
  //   We use exactly the same handle for AOD and miniAOD formats
  // since pat::Photon objects can be recast as reco::Photon objects.
  edm::Handle<edm::View<reco::Photon> > photons;
  bool isAOD = true;
  iEvent.getByToken(photonsToken_, photons);
  if( !photons.isValid() ){
    isAOD = false;
    iEvent.getByToken(photonsMiniAODToken_,photons);
  }

  // Get generator level info
  Handle<edm::View<reco::GenParticle> > genParticles;
  if( isAOD )
    iEvent.getByToken(genParticlesToken_,genParticles);
  else
    iEvent.getByToken(genParticlesMiniAODToken_,genParticles);

  // Get rho
  edm::Handle< double > rhoH;
  iEvent.getByToken(rhoToken_,rhoH);
  rho_ = *rhoH;

 // Get the isolation maps for CITK
  edm::Handle<edm::ValueMap<float> > phoChargedIsolationMap_CITK;
  iEvent.getByToken(phoChargedIsolationToken_CITK, phoChargedIsolationMap_CITK);
  edm::Handle<edm::ValueMap<float> > phoNeutralHadronIsolationMap_CITK;
  iEvent.getByToken(phoNeutralHadronIsolationToken_CITK, phoNeutralHadronIsolationMap_CITK);
  edm::Handle<edm::ValueMap<float> > phoPhotonIsolationMap_CITK;
  iEvent.getByToken(phoPhotonIsolationToken_CITK, phoPhotonIsolationMap_CITK);

  //generator info
  Handle <GenEventInfoProduct> genInfo; 
   iEvent.getByToken( genInfoToken , genInfo);

  // Clear vectors
  nPhotons_ = 0;
  pt_.clear();
  eta_.clear();
  phi_.clear();
  //
  hOverE_.clear();
  hasPixelSeed_.clear();
  //
  isoChargedHadrons_CITK_.clear();
  isoNeutralHadrons_CITK_.clear();
  isoPhotons_CITK_.clear();
  //
  isoChargedHadrons_pf_.clear();
  isoNeutralHadrons_pf_.clear();
  isoPhotons_pf_.clear();
  //
  relisoWithEA_CITK_.clear();
  relisoWithEA_pf_.clear();
  r9.clear();
  //
  isTrue_.clear();
  genWeight.clear();

  // Loop over photons
  for (size_t i = 0; i < photons->size(); ++i){
    const auto pho = photons->ptrAt(i);

    // Kinematics
    if( pho->pt() < 15 ) 
      continue;
    
    nPhotons_++;

    //
    // Save photon kinematics
    //
    pt_  .push_back( pho->pt() );
    eta_ .push_back( pho->superCluster()->eta() );
    phi_ .push_back( pho->superCluster()->phi() );

    hOverE_                .push_back( pho->hadTowOverEm() );
    hasPixelSeed_          .push_back( (Int_t)pho->hasPixelSeed() );
    r9.push_back( pho->r9() );

    //isolations from CITK
    float chIso_CITK =  (*phoChargedIsolationMap_CITK)[pho];
    float nhIso_CITK =  (*phoNeutralHadronIsolationMap_CITK)[pho];
    float phIso_CITK = (*phoPhotonIsolationMap_CITK)[pho];

    isoChargedHadrons_CITK_.push_back( chIso_CITK );
    isoNeutralHadrons_CITK_.push_back( nhIso_CITK );
    isoPhotons_CITK_       .push_back( phIso_CITK );

    //isolations with effective area correction
    float abseta = fabs( pho->superCluster()->eta());
    
    //pfIsolation variables
    float chIso_pf = pho -> chargedHadronIso() ;
    float nhIso_pf = pho -> neutralHadronIso();
    float phIso_pf = pho -> photonIso();

    isoChargedHadrons_pf_.push_back( chIso_pf );
    isoNeutralHadrons_pf_.push_back( nhIso_pf );
    isoPhotons_pf_       .push_back( phIso_pf );

    //relative isolations
    float Area = effAreaChHadrons_.getEffectiveArea(abseta) + effAreaNeuHadrons_.getEffectiveArea(abseta) + effAreaPhotons_.getEffectiveArea(abseta);
    relisoWithEA_CITK_.push_back((std::max( (float)0.0, chIso_CITK + nhIso_CITK + phIso_CITK - rho_*Area))/(pho -> pt()) );
    relisoWithEA_pf_.push_back((std::max( (float)0.0, chIso_pf + nhIso_pf + phIso_pf - rho_*Area )) /(pho -> pt()) ); 
    // Save MC truth match
    isTrue_.push_back( matchToTruth(*pho, genParticles) );
    //generator weight
    genWeight.push_back (  (genInfo -> weight()) > 0 ? 1 : -1 ) ;
   }
   
  // Save the info
  photonTree_->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void 
SimplePhotonNtupler::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SimplePhotonNtupler::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
SimplePhotonNtupler::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
SimplePhotonNtupler::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
SimplePhotonNtupler::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
SimplePhotonNtupler::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SimplePhotonNtupler::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

int SimplePhotonNtupler::matchToTruth(const reco::Photon &pho, 
				   const edm::Handle<edm::View<reco::GenParticle>>  
				   &genParticles)
{
  // 
  // Explicit loop and geometric matching method 
  //

  // Find the closest status 1 gen photon to the reco photon
  double dR = 999;
  const reco::Candidate *closestPhoton = 0;
  for(size_t i=0; i<genParticles->size();i++){
    const reco::Candidate *particle = &(*genParticles)[i];
    // Drop everything that is not photon or not status 1
    if( abs(particle->pdgId()) != 22 || particle->status() != 1 )
      continue;
    //
    double dRtmp = ROOT::Math::VectorUtil::DeltaR( pho.p4(), particle->p4() );
    if( dRtmp < dR ){
      dR = dRtmp;
      closestPhoton = particle;
    }
  }
  // See if the closest photon (if it exists) is close enough.
  // If not, no match found.
  if( !(closestPhoton != 0 && dR < 0.1) ) {
    return UNMATCHED;
  }

  // Find ID of the parent of the found generator level photon match
  int ancestorPID = -999; 
  int ancestorStatus = -999;
  findFirstNonPhotonMother(closestPhoton, ancestorPID, ancestorStatus);

  // Allowed parens: quarks pdgId 1-5, or a gluon 21
  std::vector<int> allowedParents { -1, 1, -2, 2, -3, 3, -4, 4, -5, 5, -21, 21 };
  if( !(std::find(allowedParents.begin(), 
		 allowedParents.end(), ancestorPID)
	!= allowedParents.end()) ){
    // So it is not from g, u, d, s, c, b. Check if it is from pi0 or not. 
    if( abs(ancestorPID) == 111 )
      return MATCHED_FROM_PI0;
    else
      return MATCHED_FROM_OTHER_SOURCES;
  }
  return MATCHED_FROM_GUDSCB;
   
}

void SimplePhotonNtupler::findFirstNonPhotonMother(const reco::Candidate *particle,
						int &ancestorPID, int &ancestorStatus){
  
  if( particle == 0 ){
    printf("SimplePhotonNtupler: ERROR! null candidate pointer, this should never happen\n");
    return;
  }

  // Is this the first non-photon parent? If yes, return, otherwise
  // go deeper into recursion
  if( abs(particle->pdgId()) == 22 ){
    findFirstNonPhotonMother(particle->mother(0), ancestorPID, ancestorStatus);
  }else{
    ancestorPID = particle->pdgId();
    ancestorStatus = particle->status();
  }
  
  return;
}


//define this as a plug-in
DEFINE_FWK_MODULE(SimplePhotonNtupler);
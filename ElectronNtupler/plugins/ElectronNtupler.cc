// -*- C++ -*-
//
// Package:    ElectronWork/ElectronNtupler
// Class:      ElectronNtupler
// 
/**\class ElectronNtupler ElectronNtupler.cc ElectronWork/ElectronNtupler/plugins/ElectronNtupler.cc

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
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

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
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefVector.h"
#include "TTree.h"
#include "Math/VectorUtil.h"
#include "DataFormats/Common/interface/RefToPtr.h"



namespace reco {
  typedef edm::Ptr<reco::GsfElectron> GsfElectronPtr;
}

namespace pat {
  typedef edm::Ptr<pat::Electron> PatElectronPtr;
}
//
// class declaration
//

class ElectronNtupler : public edm::EDAnalyzer {
public:
  explicit ElectronNtupler(const edm::ParameterSet&);
  ~ElectronNtupler();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  enum ElectronMatchType {UNMATCHED = 0, 
			  TRUE_PROMPT_ELECTRON, 
			  TRUE_ELECTRON_FROM_TAU,
			  TRUE_NON_PROMPT_ELECTRON}; // The last does not include tau parents
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  
  // MC truth matching utilities
  // The function that uses algorith from Josh Bendavid with 
  // an explicit loop over gen particles. 
  int matchToTruth(const pat::Electron &el, const edm::Handle<edm::View<reco::GenParticle>>  &prunedGenParticles);
  // The function that uses the standard genParticle() matching for electrons.
  int matchToTruthAlternative(const pat::Electron &el);
  
  bool checkAncestor(const reco::Candidate *gen, int ancestorPid);
  void findFirstNonElectronMother(const reco::Candidate *particle, int &ancestorPID, int &ancestorStatus);
  void printAllZeroMothers(const reco::Candidate *particle);
  
  // ----------member data ---------------------------
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  edm::EDGetTokenT<edm::View<PileupSummaryInfo> > pileupToken_;
  edm::EDGetTokenT<edm::View<reco::GsfElectron>> electronToken_;
  edm::EDGetTokenT<edm::View<reco::PFCandidate> > cands_;
  edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenToken_;
  edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > packedGenToken_;
  edm::EDGetTokenT<double> rhoToken_;
  
  //CITK
  edm::EDGetTokenT<edm::ValueMap<float> > ValueMaps_ChargedHadrons_;
  edm::EDGetTokenT<edm::ValueMap<float> > ValueMaps_NeutralHadrons_;
  edm::EDGetTokenT<edm::ValueMap<float> > ValueMaps_Photons_;
  
  //PUPPI
  edm::EDGetTokenT<edm::ValueMap<float> > ValueMaps_PUPPI_ChargedHadrons_;
  edm::EDGetTokenT<edm::ValueMap<float> > ValueMaps_PUPPI_NeutralHadrons_;
  edm::EDGetTokenT<edm::ValueMap<float> > ValueMaps_PUPPI_Photons_;
  
  //PUPPI no leptons
  edm::EDGetTokenT<edm::ValueMap<float> > ValueMaps_PUPPI_NoLeptons_ChargedHadrons_;
  edm::EDGetTokenT<edm::ValueMap<float> > ValueMaps_PUPPI_NoLeptons_NeutralHadrons_;
  edm::EDGetTokenT<edm::ValueMap<float> > ValueMaps_PUPPI_NoLeptons_Photons_;

  edm::EDGetTokenT<edm::ValueMap<bool> > ValueMap_ids_wp80_Token;
  edm::EDGetTokenT<edm::ValueMap<bool> > ValueMap_ids_wp90_Token; 

  edm::EDGetTokenT<GenEventInfoProduct> genInfoToken;

  TTree *electronTree_;
  
  // Vars for pile-up
  Int_t nPUTrue_;    // true pile-up
  Int_t nPU_;        // generated pile-up
  Int_t nPV_;        // number of reconsrtucted primary vertices
  Float_t rho_;      // the rho variable
  
  //event info
  int nevent, run, lumi;

  // all electron variables
  Float_t pt_;
  Float_t etaSC_;
  Float_t dEtaIn_;
  Float_t dPhiIn_;
  Float_t hOverE_;
  // Float_t sigmaIetaIeta_;
  Float_t full5x5_sigmaIetaIeta_;
  Float_t relIsoWithEA_;
  Float_t relIsoWithDBeta_;
  Float_t ooEmooP_;
  Float_t d0_;
  Float_t dz_;
  Int_t   expectedMissingInnerHits_;
  
  Int_t classification;
  
  // I comment this because it is not accessible in AOD
  //Int_t   passConversionVeto_;     
  Int_t   isTrueElectron_;
  
  Float_t isoChargedHadrons_;
  Float_t isoNeutralHadrons_;
  Float_t isoPhotons_;
  Float_t isoChargedFromPU_;
  
  Float_t relisoChargedHadrons_;
  Float_t relisoNeutralHadrons_;
  Float_t relisoPhotons_;
  
  //CITK
  Float_t sumChargedHadronPt_CITK;
  Float_t sumNeutralHadronPt_CITK;
  Float_t sumPhotonPt_CITK;
  
  //PUPPI
  Float_t sumChargedHadronPt_PUPPI;
  Float_t sumNeutralHadronPt_PUPPI;
  Float_t sumPhotonPt_PUPPI;
  
  //PUPPINoLeptons
  Float_t sumChargedHadronPt_PUPPI_NoLeptons;
  Float_t sumNeutralHadronPt_PUPPI_NoLeptons;
  Float_t sumPhotonPt_PUPPI_NoLeptons;
  
  Float_t reliso_PUPPI, reliso_PUPPI_NoLeptons;
  Float_t reliso_PUPPI_average;
  
  Float_t relisoChargedHadronPt_CITK;
  Float_t relisoNeutralHadronPt_CITK;
  Float_t relisoPhotonPt_CITK;
  
  Float_t relisoChargedHadronPt_PUPPI;
  Float_t relisoNeutralHadronPt_PUPPI;
  Float_t relisoPhotonPt_PUPPI;
  
  Float_t relisoChargedHadronPt_PUPPINoLeptons;
  Float_t relisoNeutralHadronPt_PUPPINoLeptons;
  Float_t relisoPhotonPt_PUPPINoLeptons;

  bool mvaIDBit_w80;
  bool mvaIDBit_w90;
  bool PF_ID;
  Int_t genWeight;
  bool isEB;
};

//
// constants, enums and typedefs
//

// Effective areas for electrons from Giovanni P. and Cristina
// distributed as private slides in Jan 2015, derived for PHYS14
namespace EffectiveAreas {
  const int nEtaBins = 5;
  const float etaBinLimits[nEtaBins+1] = {
    0.0, 0.8, 1.3, 2.0, 2.2, 2.5};
  const float effectiveAreaValues[nEtaBins] = {
    0.1013, 0.0988, 0.0572, 0.0842, 0.1530};
}

//
// static data member definitions
//

//
// constructors and destructor
//
ElectronNtupler::ElectronNtupler(const edm::ParameterSet& iConfig):
  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  pileupToken_(consumes<edm::View<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("pileup"))),
  electronToken_(consumes<edm::View<reco::GsfElectron>>(iConfig.getParameter<edm::InputTag>("electrons"))),
  cands_(consumes<edm::View<reco::PFCandidate> > (iConfig.getParameter<edm::InputTag>( "cand_src" ) ) ),
  prunedGenToken_(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("pruned"))),
  rhoToken_(consumes<double> (iConfig.getParameter<edm::InputTag>("rho"))),
  
  //CITK
  ValueMaps_ChargedHadrons_(consumes<edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>( "ValueMaps_ChargedHadrons_src" ) ) ),
  ValueMaps_NeutralHadrons_(consumes<edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>( "ValueMaps_NeutralHadrons_src" ) ) ),
  ValueMaps_Photons_(consumes<edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>( "ValueMaps_Photons_src" ) ) ),
  //PUPPI
  ValueMaps_PUPPI_ChargedHadrons_(consumes<edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>( "ValueMaps_PUPPI_ChargedHadrons_src" ) ) ),
  ValueMaps_PUPPI_NeutralHadrons_(consumes<edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>( "ValueMaps_PUPPI_NeutralHadrons_src" ) ) ),
  ValueMaps_PUPPI_Photons_(consumes<edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>( "ValueMaps_PUPPI_Photons_src" ) ) ),
  //PUPPINoLeptons
  ValueMaps_PUPPI_NoLeptons_ChargedHadrons_(consumes<edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>( "ValueMaps_PUPPI_NoLeptons_ChargedHadrons_src" ) ) ),
  ValueMaps_PUPPI_NoLeptons_NeutralHadrons_(consumes<edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>( "ValueMaps_PUPPI_NoLeptons_NeutralHadrons_src" ) ) ),
  ValueMaps_PUPPI_NoLeptons_Photons_(consumes<edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>( "ValueMaps_PUPPI_NoLeptons_Photons_src" ) ) ),

 ValueMap_ids_wp80_Token(consumes<edm::ValueMap<bool> > (iConfig.getParameter<edm::InputTag>( "mva_idw80_src" ) ) ),
 ValueMap_ids_wp90_Token(consumes<edm::ValueMap<bool> > (iConfig.getParameter<edm::InputTag>( "mva_idw90_src" ) ) ),
 genInfoToken(consumes<GenEventInfoProduct> (iConfig.getParameter<edm::InputTag>( "genInfo" ) ) )
{
  edm::Service<TFileService> fs;
  electronTree_ = fs->make<TTree> ("ElectronTree", "Electron data");
  
  
  //event info
  electronTree_->Branch("event",	      &nevent,    	  "event/I"           );
  electronTree_->Branch("lumi", 	      &lumi,   		  "lumi/I"  		);
  electronTree_->Branch("run",	      &run,		  "run/I"  	       );
  
  electronTree_->Branch("nPV"        ,  &nPV_     , "nPV/I");
  electronTree_->Branch("nPU"        ,  &nPU_     , "nPU/I");
  electronTree_->Branch("nPUTrue"    ,  &nPUTrue_ , "nPUTrue/I");
  electronTree_->Branch("rho"        ,  &rho_ , "rho/F");

  electronTree_->Branch("pt"    ,  &pt_    , "pt/F");			    
  electronTree_->Branch("etaSC" ,  &etaSC_ , "etaSC/F");
  electronTree_->Branch("dEtaIn",  &dEtaIn_, "dEtaIn/F");
  electronTree_->Branch("dPhiIn",  &dPhiIn_, "dPhiIn/F");
  electronTree_->Branch("hOverE",  &hOverE_, "hOverE/F");
  // electronTree_->Branch("sigmaIetaIeta",         &sigmaIetaIeta_, "sigmaIetaIeta/F");
  electronTree_->Branch("full5x5_sigmaIetaIeta", &full5x5_sigmaIetaIeta_, "full5x5_sigmaIetaIeta/F");
 
  electronTree_->Branch("isoChargedFromPU"       , &isoChargedFromPU_);
  electronTree_->Branch("relIsoWithEA"           , &relIsoWithEA_, "relIsoWithEA/F");
  electronTree_->Branch("relIsoWithDBeta"      , &relIsoWithDBeta_, "relIsoWithDBeta/F");
  electronTree_->Branch("ooEmooP", &ooEmooP_, "ooEmooP/F");
  electronTree_->Branch("d0"     , &d0_,      "d0/F");
  electronTree_->Branch("dz"     , &dz_,      "dz/F");
  electronTree_->Branch("expectedMissingInnerHits", &expectedMissingInnerHits_, "expectedMissingInnerHits/I");
  
  electronTree_->Branch("classification"    , &classification,     "classification/I");
  
  //electronTree_->Branch("passConversionVeto", &passConversionVeto_, "passConversionVeto/I");
  electronTree_->Branch("isTrueElectron"    , &isTrueElectron_,     "isTrueElectron/I");
 
  electronTree_->Branch("isoChargedHadrons"      , &isoChargedHadrons_ , "isoChargedHadrons/F");
  electronTree_->Branch("isoNeutralHadrons"      , &isoNeutralHadrons_, "isoNeutralHadrons/F");
  electronTree_->Branch("isoPhotons"             , &isoPhotons_, "isoPhotons/F");
  
  electronTree_->Branch("relisoChargedHadrons"      , &relisoChargedHadrons_ , "relisoChargedHadrons/F");
  electronTree_->Branch("relisoNeutralHadrons"      , &relisoNeutralHadrons_, "relisoNeutralHadrons/F");
  electronTree_->Branch("relisoPhotons"             , &relisoPhotons_, "relisoPhotons/F");
  
  //CITK
  electronTree_ -> Branch("sumChargedHadronPt_CITK", &sumChargedHadronPt_CITK, "sumChargedHadronPt_CITK/F");
  electronTree_ -> Branch("sumNeutralHadronPt_CITK", &sumNeutralHadronPt_CITK, "sumNeutralHadronPt_CITK/F");
  electronTree_ -> Branch("sumPhotonPt_CITK", &sumPhotonPt_CITK, "sumPhotonPt_CITK/F");

  //PUPPI
  electronTree_ -> Branch("sumChargedHadronPt_PUPPI", &sumChargedHadronPt_PUPPI, "sumChargedHadronPt_PUPPI/F");
  electronTree_ -> Branch("sumNeutralHadronPt_PUPPI", &sumNeutralHadronPt_PUPPI, "sumNeutralHadronPt_PUPPI/F");
  electronTree_ -> Branch("sumPhotonPt_PUPPI", &sumPhotonPt_PUPPI, "sumPhotonPt_PUPPI/F");
  
  electronTree_ -> Branch("reliso_PUPPI", &reliso_PUPPI, "reliso_PUPPI/F");
  
    //PUPPI
  electronTree_ -> Branch("sumChargedHadronPt_PUPPI_NoLeptons", &sumChargedHadronPt_PUPPI_NoLeptons, "sumChargedHadronPt_PUPPI_NoLeptons/F");
  electronTree_ -> Branch("sumNeutralHadronPt_PUPPI_NoLeptons", &sumNeutralHadronPt_PUPPI_NoLeptons, "sumNeutralHadronPt_PUPPI_NoLeptons/F");
  electronTree_ -> Branch("sumPhotonPt_PUPPI_NoLeptons", &sumPhotonPt_PUPPI_NoLeptons, "sumPhotonPt_PUPPI_NoLeptons/F");
  
  electronTree_ -> Branch("reliso_PUPPI_NoLeptons", &reliso_PUPPI_NoLeptons, "reliso_PUPPI_NoLeptons/F");

  electronTree_ -> Branch("reliso_PUPPI_average", &reliso_PUPPI_average, "reliso_PUPPI_average/F");
  
    
  electronTree_ -> Branch("relisoChargedHadronPt_CITK", &relisoChargedHadronPt_CITK, "relisoChargedHadronPt_CITK/F");
  electronTree_ -> Branch("relisoNeutralHadronPt_CITK", &relisoNeutralHadronPt_CITK, "relisoNeutralHadronPt_CITK/F");
  electronTree_ -> Branch("relisoPhotonPt_CITK", &relisoPhotonPt_CITK, "relisoPhotonPt_CITK/F");
  
  electronTree_ -> Branch("relisoChargedHadronPt_PUPPI", &relisoChargedHadronPt_PUPPI, "relisoChargedHadronPt_PUPPI/F");
  electronTree_ -> Branch("relisoNeutralHadronPt_PUPPI", &relisoNeutralHadronPt_PUPPI, "relisoNeutralHadronPt_PUPPI/F");
  electronTree_ -> Branch("relisoPhotonPt_PUPPI", &relisoPhotonPt_PUPPI, "relisoPhotonPt_PUPPI/F");
  
  electronTree_ -> Branch("relisoChargedHadronPt_PUPPINoLeptons", &relisoChargedHadronPt_PUPPINoLeptons, "relisoChargedHadronPt_PUPPINoLeptons/F");
  electronTree_ -> Branch("relisoNeutralHadronPt_PUPPINoLeptons", &relisoNeutralHadronPt_PUPPINoLeptons, "relisoNeutralHadronPt_PUPPINoLeptons/F");
  electronTree_ -> Branch("relisoPhotonPt_PUPPINoLeptons", &relisoPhotonPt_PUPPINoLeptons, "relisoPhotonPt_PUPPINoLeptons/F");

  electronTree_ -> Branch("mvaIDBit_w80", &mvaIDBit_w80, "mvaIDBit_w80/B");
  electronTree_ -> Branch("mvaIDBit_w90", &mvaIDBit_w90, "mvaIDBit_w90/B");

  electronTree_ -> Branch("PF_ID", &PF_ID, "PF_ID/B");

  electronTree_ -> Branch("genWeight", &genWeight, "genWeight/I");
  
  electronTree_ -> Branch("isEB", &isEB, "isEB/B");

}


ElectronNtupler::~ElectronNtupler()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
ElectronNtupler::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace std;
  using namespace edm;
  using namespace reco;
  
  
  //event info
  nevent = iEvent.eventAuxiliary().event();
  run    = iEvent.eventAuxiliary().run();
  lumi   = iEvent.eventAuxiliary().luminosityBlock();
  
  // Pruned particles are the one containing "important" stuff
 
  Handle<edm::View<reco::GenParticle> > prunedGenParticles;
  iEvent.getByToken(prunedGenToken_,prunedGenParticles);

  // Get Pileup info
  Handle<edm::View<PileupSummaryInfo> > pileupHandle;
  iEvent.getByToken(pileupToken_, pileupHandle);
  for( auto & puInfoElement : *pileupHandle){
    if( puInfoElement.getBunchCrossing() == 0 ){
      nPU_    = puInfoElement.getPU_NumInteractions();
      nPUTrue_= puInfoElement.getTrueNumInteractions();
    }
  }

  // Get PV
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);
  if (vertices->empty()) return; // skip the event if no PV found
  //const reco::Vertex &pv = vertices->front();
  
  nPV_ = vertices -> size();

  VertexCollection::const_iterator firstGoodVertex = vertices->end();
  int firstGoodVertexIdx = 0;
  for (VertexCollection::const_iterator vtx = vertices->begin(); 
       vtx != vertices->end(); ++vtx, ++firstGoodVertexIdx) {
    // The "good vertex" selection is borrowed from Giovanni Zevi Della Porta
    // Replace isFake() for miniAOD because it requires tracks and miniAOD vertices don't have tracks:
    // Vertex.h: bool isFake() const {return (chi2_==0 && ndof_==0 && tracks_.empty());}
    if (  /*!vtx->isFake() &&*/ 
	!(vtx->chi2()==0 && vtx->ndof()==0) 
	&&  vtx->ndof()>=4. && vtx->position().Rho()<=2.0
	&& std::abs(vtx->position().Z())<=24.0) {
      firstGoodVertex = vtx;
      break;
    }
  }
  
  if ( firstGoodVertex==vertices->end() )
    return; // skip event if there are no good PVs
  
  
  edm::Handle< double > rhoH;
  iEvent.getByToken(rhoToken_,rhoH);
  rho_ = *rhoH;
  
  // Get electron collection
  Handle<edm::View<reco::GsfElectron>> electrons;
  iEvent.getByToken(electronToken_, electrons);
  Handle<edm::View<reco::PFCandidate> > cands;
  iEvent.getByToken(cands_, cands);
  
  //CITK
  Handle <edm::ValueMap <float> > ValueMaps_ChargedHadrons, ValueMaps_NeutralHadrons, ValueMaps_Photons;
  Handle <edm::ValueMap <float> > ValueMaps_PUPPI_ChargedHadrons, ValueMaps_PUPPI_NeutralHadrons, ValueMaps_PUPPI_Photons;
  Handle <edm::ValueMap <float> > ValueMaps_PUPPI_NoLeptons_ChargedHadrons, ValueMaps_PUPPI_NoLeptons_NeutralHadrons, ValueMaps_PUPPI_NoLeptons_Photons;
  
  //CITK
  iEvent.getByToken( ValueMaps_ChargedHadrons_ , ValueMaps_ChargedHadrons);
  iEvent.getByToken( ValueMaps_NeutralHadrons_ , ValueMaps_NeutralHadrons);
  iEvent.getByToken( ValueMaps_Photons_ , ValueMaps_Photons);
  
  //PUPPI 
  iEvent.getByToken( ValueMaps_PUPPI_ChargedHadrons_ , ValueMaps_PUPPI_ChargedHadrons);
  iEvent.getByToken( ValueMaps_PUPPI_NeutralHadrons_ , ValueMaps_PUPPI_NeutralHadrons);
  iEvent.getByToken( ValueMaps_PUPPI_Photons_ , ValueMaps_PUPPI_Photons);
  
  //PUPPI_NoLeptons 
  iEvent.getByToken( ValueMaps_PUPPI_NoLeptons_ChargedHadrons_ , ValueMaps_PUPPI_NoLeptons_ChargedHadrons);
  iEvent.getByToken( ValueMaps_PUPPI_NoLeptons_NeutralHadrons_ , ValueMaps_PUPPI_NoLeptons_NeutralHadrons);
  iEvent.getByToken( ValueMaps_PUPPI_NoLeptons_Photons_ , ValueMaps_PUPPI_NoLeptons_Photons);

  Handle <edm::ValueMap <bool> >  ValueMap_ids_wp80, ValueMap_ids_wp90;

  iEvent.getByToken( ValueMap_ids_wp80_Token , ValueMap_ids_wp80);
  iEvent.getByToken( ValueMap_ids_wp90_Token , ValueMap_ids_wp90);

   Handle <GenEventInfoProduct> genInfo; 
   iEvent.getByToken( genInfoToken , genInfo);
   
  //
  // Loop over electrons
  //
  // printf("DEBUG: new event\n");

  for (unsigned int iElectron = 0; iElectron < electrons -> size(); iElectron++) {
    
    auto elePtr = electrons -> ptrAt(iElectron);
    reco::GsfElectronPtr eleGsfPtr(elePtr);
    // Kinematics
    pt_ = eleGsfPtr -> pt();
   
    // Keep only electrons above 10 GeV.
    // NOTE: miniAOD does not store some of the info for electrons <5 GeV at all!
    if( pt_ < 10 ) 
      continue;
    
    etaSC_ = eleGsfPtr -> superCluster()->eta();
    
    // ID and matching
    dEtaIn_ = eleGsfPtr -> deltaEtaSuperClusterTrackAtVtx();
    dPhiIn_ = eleGsfPtr -> deltaPhiSuperClusterTrackAtVtx();
    hOverE_ = eleGsfPtr -> hcalOverEcal();
    // sigmaIetaIeta_ = el.sigmaIetaIeta();
    full5x5_sigmaIetaIeta_ = eleGsfPtr -> full5x5_sigmaIetaIeta();
    // |1/E-1/p| = |1/E - EoverPinner/E| is computed below
    // The if protects against ecalEnergy == inf or zero (always
    // the case for electrons below 5 GeV in miniAOD)
    if( eleGsfPtr -> ecalEnergy() == 0 ){
      printf("Electron energy is zero!\n");
      ooEmooP_ = 1e30;
    }else if( !std::isfinite(eleGsfPtr -> ecalEnergy())){
      printf("Electron energy is not finite!\n");
      ooEmooP_ = 1e30;
    }else{
      ooEmooP_ = std::abs(1.0/eleGsfPtr -> ecalEnergy() - eleGsfPtr -> eSuperClusterOverP()/eleGsfPtr -> ecalEnergy() );
    }
    
    // Isolation
    GsfElectron::PflowIsolationVariables pfIso = eleGsfPtr -> pfIsolationVariables();
    isoChargedHadrons_ = pfIso.sumChargedHadronPt ;
    isoNeutralHadrons_ = pfIso.sumNeutralHadronEt ;
    isoPhotons_        = pfIso.sumPhotonEt ;
    isoChargedFromPU_  = pfIso.sumPUPt ;
    
    relisoChargedHadrons_ = pfIso.sumChargedHadronPt/pt_;
    relisoNeutralHadrons_ = pfIso.sumNeutralHadronEt/pt_;
    relisoPhotons_        = pfIso.sumPhotonEt/pt_;
    // Compute isolation with effective area correction for PU
    // Find eta bin first. If eta>2.5, the last eta bin is used.
    int etaBin = 0; 
    while ( etaBin < EffectiveAreas::nEtaBins-1 
	    && abs(etaSC_) > EffectiveAreas::etaBinLimits[etaBin+1] )
      { ++etaBin; };
    double area = EffectiveAreas::effectiveAreaValues[etaBin];
    relIsoWithEA_ = ( pfIso.sumChargedHadronPt + max(0.0, pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - rho_ * area ) )/pt_;
    
    // Compute isolation with delta beta correction for PU
    float absiso = pfIso.sumChargedHadronPt + max(0.0 , pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - 0.5 * pfIso.sumPUPt );
    relIsoWithDBeta_ = absiso/pt_;
    
    // Impact parameter
    d0_ = (-1) * eleGsfPtr -> gsfTrack()->dxy(firstGoodVertex->position() );
    dz_ = eleGsfPtr -> gsfTrack()->dz( firstGoodVertex->position() );
    
    // Conversion rejection
    // pre-72X method below is commented out
    //expectedMissingInnerHits_ = el.gsfTrack()->trackerExpectedHitsInner().numberOfLostHits();
    // since 72X, the access of missing hits is this:
    expectedMissingInnerHits_ = eleGsfPtr -> gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);
    // I comment this because it is not accessible in AOD
    //passConversionVeto_ = elePatPtr -> passConversionVeto();
    
     classification = eleGsfPtr -> classification();
    // Match to generator level truth
    
    // 
    // Explicit loop over gen candidates method
    //
    isTrueElectron_ = matchToTruth( eleGsfPtr, prunedGenParticles);
    
     // I comment this because it is not accessible in AOD
    // isTrueElectronAlternative_ = matchToTruthAlternative( el );
    
    // For debug purposes, one can use this utility that prints 
    // the decay history, using standard matching in this case:
    //   printAllZeroMothers( el.genParticle() );

    //CITK
    sumChargedHadronPt_CITK =  (*ValueMaps_ChargedHadrons)[eleGsfPtr];
    sumNeutralHadronPt_CITK =  (*ValueMaps_NeutralHadrons)[eleGsfPtr];
    sumPhotonPt_CITK        =  (*ValueMaps_Photons)[eleGsfPtr];

     
    //PUPPI
    sumChargedHadronPt_PUPPI =  (*ValueMaps_PUPPI_ChargedHadrons)[elePtr];
    sumNeutralHadronPt_PUPPI =  (*ValueMaps_PUPPI_NeutralHadrons)[elePtr];
    sumPhotonPt_PUPPI        =  (*ValueMaps_PUPPI_Photons)[elePtr];
    
    //PUPPINoLeptons
    sumChargedHadronPt_PUPPI_NoLeptons =  (*ValueMaps_PUPPI_NoLeptons_ChargedHadrons)[elePtr];
    sumNeutralHadronPt_PUPPI_NoLeptons =  (*ValueMaps_PUPPI_NoLeptons_NeutralHadrons)[elePtr];
    sumPhotonPt_PUPPI_NoLeptons        =  (*ValueMaps_PUPPI_NoLeptons_Photons)[elePtr];
    
    //CITK
    relisoChargedHadronPt_CITK = sumChargedHadronPt_CITK/pt_;
    relisoNeutralHadronPt_CITK = sumNeutralHadronPt_CITK/pt_;
    relisoPhotonPt_CITK = sumPhotonPt_CITK/pt_;
    
    //PUPPI
    relisoChargedHadronPt_PUPPI = sumChargedHadronPt_PUPPI/pt_;
    relisoNeutralHadronPt_PUPPI = sumNeutralHadronPt_PUPPI/pt_;
    relisoPhotonPt_PUPPI = sumPhotonPt_PUPPI/pt_;
    
    //PUPPI_NoLeptons
    relisoChargedHadronPt_PUPPINoLeptons = sumChargedHadronPt_PUPPI_NoLeptons/pt_;
    relisoNeutralHadronPt_PUPPINoLeptons = sumNeutralHadronPt_PUPPI_NoLeptons/pt_;
    relisoPhotonPt_PUPPINoLeptons = sumPhotonPt_PUPPI_NoLeptons/pt_;
    
    //PUPPI total isolation   
    reliso_PUPPI = (sumChargedHadronPt_PUPPI + sumNeutralHadronPt_PUPPI + sumPhotonPt_PUPPI)/pt_;
    //PUPPINoLeptons total isolation
    reliso_PUPPI_NoLeptons = (sumChargedHadronPt_PUPPI_NoLeptons + sumNeutralHadronPt_PUPPI_NoLeptons + sumPhotonPt_PUPPI_NoLeptons)/pt_;

    //reliso puppi average
    reliso_PUPPI_average = 0.5*(reliso_PUPPI_NoLeptons + reliso_PUPPI);
    
    const reco::CaloClusterPtr& seed = eleGsfPtr -> superCluster()->seed();
    isEB = ( seed->seed().subdetId() == EcalBarrel );

    //saving id bits
    mvaIDBit_w80  =  (*ValueMap_ids_wp80)[elePtr];
    mvaIDBit_w90  =  (*ValueMap_ids_wp90)[elePtr];

    reco::GsfTrackRef trackref_ele = elePtr -> gsfTrack() ;

    PF_ID = false;
      
    for (unsigned iCand = 0; iCand < cands -> size(); iCand ++)
    {
      if (  std::abs((cands -> at(iCand).pdgId()) == 11 ) && trackref_ele  == cands -> at(iCand).gsfTrackRef() )   PF_ID = true;
    }

    genWeight = (genInfo -> weight()) > 0 ? 1 : -1;
         
    // Save this electron's info
    electronTree_->Fill();
  }
  
}


// ------------ method called once each job just before starting event loop  ------------
void 
ElectronNtupler::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ElectronNtupler::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
ElectronNtupler::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
ElectronNtupler::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
ElectronNtupler::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
ElectronNtupler::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ElectronNtupler::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

bool ElectronNtupler::checkAncestor(const reco::Candidate *gen, int ancestorPid){

  // General sanity check
  if( gen == 0 ){
    printf("ElectronNtupler::checkAncestor: ERROR null particle is passed in, ignore it.\n");
    return false;
  }

  // If this is true, we found our target ancestor
  if( abs( gen->pdgId() ) == ancestorPid )
    return true;

  // Go deeper and check all mothers
  for(size_t i=0;i< gen->numberOfMothers();i++) {
    if ( checkAncestor( gen->mother(i), ancestorPid) )
      return true;
  }
  
  return false;
}

// The function that uses algorith from Josh Bendavid with 
// an explicit loop over gen particles. 
int ElectronNtupler::matchToTruth(const pat::Electron &el, 
				  const edm::Handle<edm::View<reco::GenParticle>> &prunedGenParticles){

  // 
  // Explicit loop and geometric matching method (advised by Josh Bendavid)
  //

  // Find the closest status 1 gen electron to the reco electron
  double dR = 999;
  const reco::Candidate *closestElectron = 0;
  for(size_t i=0; i<prunedGenParticles->size();i++){
    const reco::Candidate *particle = &(*prunedGenParticles)[i];
    // Drop everything that is not electron or not status 1
    if( abs(particle->pdgId()) != 11 || particle->status() != 1 )
      continue;
    //
    double dRtmp = ROOT::Math::VectorUtil::DeltaR( el.p4(), particle->p4() );
    if( dRtmp < dR ){
      dR = dRtmp;
      closestElectron = particle;
    }
  }
  // See if the closest electron (if it exists) is close enough.
  // If not, no match found.
  if( !(closestElectron != 0 && dR < 0.1) ) {
    return UNMATCHED;
  }

  // 
  int ancestorPID = -999; 
  int ancestorStatus = -999;
  findFirstNonElectronMother(closestElectron, ancestorPID, ancestorStatus);

  if( ancestorPID == -999 && ancestorStatus == -999 ){
    // No non-electron parent??? This should never happen.
    // Complain.
    printf("ElectronNtupler: ERROR! Electron does not apper to have a non-electron parent\n");
    return UNMATCHED;
  }
  
  if( abs(ancestorPID) > 50 && ancestorStatus == 2 )
    return TRUE_NON_PROMPT_ELECTRON;

  if( abs(ancestorPID) == 15 && ancestorStatus == 2 )
    return TRUE_ELECTRON_FROM_TAU;

  // What remains is true prompt electrons
  return TRUE_PROMPT_ELECTRON;
}

void ElectronNtupler::findFirstNonElectronMother(const reco::Candidate *particle,
						 int &ancestorPID, int &ancestorStatus){

  if( particle == 0 ){
    printf("ElectronNtupler: ERROR! null candidate pointer, this should never happen\n");
    return;
  }

  // Is this the first non-electron parent? If yes, return, otherwise
  // go deeper into recursion
  if( abs(particle->pdgId()) == 11 ){
    findFirstNonElectronMother(particle->mother(0), ancestorPID, ancestorStatus);
  }else{
    ancestorPID = particle->pdgId();
    ancestorStatus = particle->status();
  }

  return;
}

// The function that uses the standard genParticle() matching for electrons.
int ElectronNtupler::matchToTruthAlternative(const pat::Electron &el){

     //
     // genParticle method
     //
     int result = UNMATCHED;

     const reco::GenParticle * gen = el.genParticle();
     if( gen != 0 ){
       int pid = gen->pdgId();
       int status = gen->status();
       bool isFromZ   = checkAncestor(gen, 23);
       bool isFromW   = checkAncestor(gen, 24);
       bool isFromTau = checkAncestor(gen, 15);
       // Check if it is a true prompt electron
       if( abs( pid ) == 11 // this is electron
	   && (status == 1 || status == 22 || status == 23 ) // NOTE: Pythia8 status here 22/23 (for Pythia6 would be 3)
	   && (isFromZ || isFromW ) && !isFromTau // comes from Z or W+-, but not from tau
	   )
	 {
	   result = TRUE_PROMPT_ELECTRON;
	 } else if ( abs( pid ) == 11 
		     && (status == 1 || status == 22 || status == 23 ) 
		     && (isFromTau ) 
		     ) 
	 {
	   // This is a true electron, but it comes from tau
	   result = TRUE_ELECTRON_FROM_TAU;
	 } else if ( abs( pid ) == 11 )
	 {
	   // This is a true electron, but it comes from something else
	   const reco::Candidate *mom = el.mother(0);
	   int momPid = -999;
	   if ( mom != 0 )
	     momPid = mom->pdgId();
	   printf("pid= %d  status= %d isFromZ= %d isFromW= %d  isFromTau= %d  momPid= %d\n", 
		  pid,  status, isFromZ, isFromW, isFromTau, momPid);
	   result = TRUE_NON_PROMPT_ELECTRON;
	 } else {
	 printf("The reco electron has a truth match with pid= %d\n", pid);
       }
     }

     return result;
}

void ElectronNtupler::printAllZeroMothers(const reco::Candidate *particle){
  
  if( particle == 0 ){
    printf("ElectronNtupler::printAllZeroMothers: reached the top of the decay tree\n");
    return;
  }
  
  printf("ElectronNtupler::printAllZeroMothers: ancestor ID= %d, status= %d\n",
	 particle->pdgId(), particle->status() );

  printAllZeroMothers( particle->mother(0) );

  return;
}

//define this as a plug-in
DEFINE_FWK_MODULE(ElectronNtupler);

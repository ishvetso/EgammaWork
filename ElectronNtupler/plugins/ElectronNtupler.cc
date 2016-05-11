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
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"



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
  
  //CITK standard isolations
  edm::EDGetTokenT<edm::ValueMap<float> > ValueMaps_ChargedHadrons_ConeVeto_;
  edm::EDGetTokenT<edm::ValueMap<float> > ValueMaps_NeutralHadrons_ConeVeto_;
  edm::EDGetTokenT<edm::ValueMap<float> > ValueMaps_Photons_ConeVeto_;

  edm::EDGetTokenT<edm::ValueMap<float> > ValueMaps_ChargedHadrons_MapBasedVeto_;
  edm::EDGetTokenT<edm::ValueMap<float> > ValueMaps_NeutralHadrons_MapBasedVeto_;
  edm::EDGetTokenT<edm::ValueMap<float> > ValueMaps_Photons_MapBasedVeto_;
  
  //PUPPI cone based veto
  edm::EDGetTokenT<edm::ValueMap<float> > ValueMaps_PUPPI_ChargedHadrons_ConeVeto_;
  edm::EDGetTokenT<edm::ValueMap<float> > ValueMaps_PUPPI_NeutralHadrons_ConeVeto_;
  edm::EDGetTokenT<edm::ValueMap<float> > ValueMaps_PUPPI_Photons_ConeVeto_;

  //PUPPI map based veto
  edm::EDGetTokenT<edm::ValueMap<float> > ValueMaps_PUPPI_ChargedHadrons_MapBasedVeto_;
  edm::EDGetTokenT<edm::ValueMap<float> > ValueMaps_PUPPI_NeutralHadrons_MapBasedVeto_;
  edm::EDGetTokenT<edm::ValueMap<float> > ValueMaps_PUPPI_Photons_MapBasedVeto_;
  
  //PUPPI no leptons cone veto
  edm::EDGetTokenT<edm::ValueMap<float> > ValueMaps_PUPPI_NoLeptons_ChargedHadrons_ConeVeto_;
  edm::EDGetTokenT<edm::ValueMap<float> > ValueMaps_PUPPI_NoLeptons_NeutralHadrons_ConeVeto_;
  edm::EDGetTokenT<edm::ValueMap<float> > ValueMaps_PUPPI_NoLeptons_Photons_ConeVeto_;

  //PUPPI no leptons map based veto
  edm::EDGetTokenT<edm::ValueMap<float> > ValueMaps_PUPPI_NoLeptons_ChargedHadrons_MapBasedVeto_;
  edm::EDGetTokenT<edm::ValueMap<float> > ValueMaps_PUPPI_NoLeptons_NeutralHadrons_MapBasedVeto_;
  edm::EDGetTokenT<edm::ValueMap<float> > ValueMaps_PUPPI_NoLeptons_Photons_MapBasedVeto_;

  edm::EDGetTokenT<edm::ValueMap<bool> > ValueMap_ids_wp80_Token;
  edm::EDGetTokenT<edm::ValueMap<bool> > ValueMap_ids_wp90_Token; 

  edm::EDGetTokenT<GenEventInfoProduct> genInfoToken;

  EffectiveAreas effArea_;

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
  
  //CITK
  Float_t sumChargedHadronPt_CITK_ConeVeto;
  Float_t sumNeutralHadronPt_CITK_ConeVeto;
  Float_t sumPhotonPt_CITK_ConeVeto;

  //CITK
  Float_t sumChargedHadronPt_CITK_MapBasedVeto;
  Float_t sumNeutralHadronPt_CITK_MapBasedVeto;
  Float_t sumPhotonPt_CITK_MapBasedVeto;
  
  //PUPPI
  Float_t sumChargedHadronPt_PUPPI_ConeVeto;
  Float_t sumNeutralHadronPt_PUPPI_ConeVeto;
  Float_t sumPhotonPt_PUPPI_ConeVeto;

  Float_t sumChargedHadronPt_PUPPI_MapBasedVeto;
  Float_t sumNeutralHadronPt_PUPPI_MapBasedVeto;
  Float_t sumPhotonPt_PUPPI_MapBasedVeto;
  
  //PUPPINoLeptons
  Float_t sumChargedHadronPt_PUPPI_NoLeptons_ConeVeto;
  Float_t sumNeutralHadronPt_PUPPI_NoLeptons_ConeVeto;
  Float_t sumPhotonPt_PUPPI_NoLeptons_ConeVeto;

  Float_t sumChargedHadronPt_PUPPI_NoLeptons_MapBasedVeto;
  Float_t sumNeutralHadronPt_PUPPI_NoLeptons_MapBasedVeto;
  Float_t sumPhotonPt_PUPPI_NoLeptons_MapBasedVeto;
  
  Float_t reliso_PUPPI_ConeVeto, reliso_PUPPI_NoLeptons_ConeVeto, reliso_PUPPI_MapBasedVeto, reliso_PUPPI_NoLeptons_MapBasedVeto;
  Float_t reliso_PUPPI_ConeVeto_average, reliso_PUPPI_MapBasedVeto_average;

  Float_t reliso_ConeVeto_raw, reliso_MapBasedVeto_raw;

  Float_t reliso_ConeVeto_EA, reliso_MapBasedVeto_EA;

  bool mvaIDBit_w80;
  bool mvaIDBit_w90;
  bool PF_ID;
  Int_t genWeight;
  bool isEB;
};

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
  ValueMaps_ChargedHadrons_ConeVeto_(consumes<edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>( "ValueMaps_ChargedHadrons_ConeVeto_src" ) ) ),
  ValueMaps_NeutralHadrons_ConeVeto_(consumes<edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>( "ValueMaps_NeutralHadrons_ConeVeto_src" ) ) ),
  ValueMaps_Photons_ConeVeto_(consumes<edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>( "ValueMaps_Photons_ConeVeto_src" ) ) ),

  ValueMaps_ChargedHadrons_MapBasedVeto_(consumes<edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>( "ValueMaps_ChargedHadrons_MapBasedVeto_src" ) ) ),
  ValueMaps_NeutralHadrons_MapBasedVeto_(consumes<edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>( "ValueMaps_NeutralHadrons_MapBasedVeto_src" ) ) ),
  ValueMaps_Photons_MapBasedVeto_(consumes<edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>( "ValueMaps_Photons_MapBasedVeto_src" ) ) ),
  //PUPPI
  ValueMaps_PUPPI_ChargedHadrons_ConeVeto_(consumes<edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>( "ValueMaps_PUPPI_ChargedHadrons_ConeVeto_src" ) ) ),
  ValueMaps_PUPPI_NeutralHadrons_ConeVeto_(consumes<edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>( "ValueMaps_PUPPI_NeutralHadrons_ConeVeto_src" ) ) ),
  ValueMaps_PUPPI_Photons_ConeVeto_(consumes<edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>( "ValueMaps_PUPPI_Photons_ConeVeto_src" ) ) ),

  ValueMaps_PUPPI_ChargedHadrons_MapBasedVeto_(consumes<edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>( "ValueMaps_PUPPI_ChargedHadrons_MapBasedVeto_src" ) ) ),
  ValueMaps_PUPPI_NeutralHadrons_MapBasedVeto_(consumes<edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>( "ValueMaps_PUPPI_NeutralHadrons_MapBasedVeto_src" ) ) ),
  ValueMaps_PUPPI_Photons_MapBasedVeto_(consumes<edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>( "ValueMaps_PUPPI_Photons_MapBasedVeto_src" ) ) ),
  //PUPPINoLeptons
  ValueMaps_PUPPI_NoLeptons_ChargedHadrons_ConeVeto_(consumes<edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>( "ValueMaps_PUPPI_NoLeptons_ChargedHadrons_ConeVeto_src" ) ) ),
  ValueMaps_PUPPI_NoLeptons_NeutralHadrons_ConeVeto_(consumes<edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>( "ValueMaps_PUPPI_NoLeptons_NeutralHadrons_ConeVeto_src" ) ) ),
  ValueMaps_PUPPI_NoLeptons_Photons_ConeVeto_(consumes<edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>( "ValueMaps_PUPPI_NoLeptons_Photons_ConeVeto_src" ) ) ),

  ValueMaps_PUPPI_NoLeptons_ChargedHadrons_MapBasedVeto_(consumes<edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>( "ValueMaps_PUPPI_NoLeptons_ChargedHadrons_MapBasedVeto_src" ) ) ),
  ValueMaps_PUPPI_NoLeptons_NeutralHadrons_MapBasedVeto_(consumes<edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>( "ValueMaps_PUPPI_NoLeptons_NeutralHadrons_MapBasedVeto_src" ) ) ),
  ValueMaps_PUPPI_NoLeptons_Photons_MapBasedVeto_(consumes<edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>( "ValueMaps_PUPPI_NoLeptons_Photons_MapBasedVeto_src" ) ) ),

 ValueMap_ids_wp80_Token(consumes<edm::ValueMap<bool> > (iConfig.getParameter<edm::InputTag>( "mva_idw80_src" ) ) ),
 ValueMap_ids_wp90_Token(consumes<edm::ValueMap<bool> > (iConfig.getParameter<edm::InputTag>( "mva_idw90_src" ) ) ),
 genInfoToken(consumes<GenEventInfoProduct> (iConfig.getParameter<edm::InputTag>( "genInfo" ) ) ),
 effArea_( (iConfig.getParameter<edm::FileInPath>("effAreaFile")).fullPath() )
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
 
  electronTree_->Branch("isoChargedHadrons"      , &isoChargedHadrons_ , "isoChargedHadrons/F");
  electronTree_->Branch("isoNeutralHadrons"      , &isoNeutralHadrons_, "isoNeutralHadrons/F");
  electronTree_->Branch("isoPhotons"             , &isoPhotons_, "isoPhotons/F");

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
  //CITK
  electronTree_ -> Branch("sumChargedHadronPt_CITK_ConeVeto", &sumChargedHadronPt_CITK_ConeVeto, "sumChargedHadronPt_CITK_ConeVeto/F");
  electronTree_ -> Branch("sumNeutralHadronPt_CITK_ConeVeto", &sumNeutralHadronPt_CITK_ConeVeto, "sumNeutralHadronPt_CITK_ConeVeto/F");
  electronTree_ -> Branch("sumPhotonPt_CITK_ConeVeto", &sumPhotonPt_CITK_ConeVeto, "sumPhotonPt_CITK_ConeVeto/F");

  electronTree_ -> Branch("sumChargedHadronPt_CITK_MapBasedVeto", &sumChargedHadronPt_CITK_MapBasedVeto, "sumChargedHadronPt_CITK_MapBasedVeto/F");
  electronTree_ -> Branch("sumNeutralHadronPt_CITK_MapBasedVeto", &sumNeutralHadronPt_CITK_MapBasedVeto, "sumNeutralHadronPt_CITK_MapBasedVeto/F");
  electronTree_ -> Branch("sumPhotonPt_CITK_MapBasedVeto", &sumPhotonPt_CITK_MapBasedVeto, "sumPhotonPt_CITK_MapBasedVeto/F");

  //PUPPI
  electronTree_ -> Branch("sumChargedHadronPt_PUPPI_ConeVeto", &sumChargedHadronPt_PUPPI_ConeVeto, "sumChargedHadronPt_PUPPI_ConeVeto/F");
  electronTree_ -> Branch("sumNeutralHadronPt_PUPPI_ConeVeto", &sumNeutralHadronPt_PUPPI_ConeVeto, "sumNeutralHadronPt_PUPPI_ConeVeto/F");
  electronTree_ -> Branch("sumPhotonPt_PUPPI_ConeVeto", &sumPhotonPt_PUPPI_ConeVeto, "sumPhotonPt_PUPPI_ConeVeto/F");

  electronTree_ -> Branch("sumChargedHadronPt_PUPPI_MapBasedVeto", &sumChargedHadronPt_PUPPI_MapBasedVeto, "sumChargedHadronPt_PUPPI_MapBasedVeto/F");
  electronTree_ -> Branch("sumNeutralHadronPt_PUPPI_MapBasedVeto", &sumNeutralHadronPt_PUPPI_MapBasedVeto, "sumNeutralHadronPt_PUPPI_MapBasedVeto/F");
  electronTree_ -> Branch("sumPhotonPt_PUPPI_MapBasedVeto", &sumPhotonPt_PUPPI_MapBasedVeto, "sumPhotonPt_PUPPI_MapBasedVeto/F");
  
  //PUPPI No Leptons
  electronTree_ -> Branch("sumChargedHadronPt_PUPPI_NoLeptons_ConeVeto", &sumChargedHadronPt_PUPPI_NoLeptons_ConeVeto, "sumChargedHadronPt_PUPPI_NoLeptons_ConeVeto/F");
  electronTree_ -> Branch("sumNeutralHadronPt_PUPPI_NoLeptons_ConeVeto", &sumNeutralHadronPt_PUPPI_NoLeptons_ConeVeto, "sumNeutralHadronPt_PUPPI_NoLeptons_ConeVeto/F");
  electronTree_ -> Branch("sumPhotonPt_PUPPI_NoLeptons_ConeVeto", &sumPhotonPt_PUPPI_NoLeptons_ConeVeto, "sumPhotonPt_PUPPI_NoLeptons_ConeVeto/F");

  electronTree_ -> Branch("sumChargedHadronPt_PUPPI_NoLeptons_MapBasedVeto", &sumChargedHadronPt_PUPPI_NoLeptons_MapBasedVeto, "sumChargedHadronPt_PUPPI_NoLeptons_MapBasedVeto/F");
  electronTree_ -> Branch("sumNeutralHadronPt_PUPPI_NoLeptons_MapBasedVeto", &sumNeutralHadronPt_PUPPI_NoLeptons_MapBasedVeto, "sumNeutralHadronPt_PUPPI_NoLeptons_MapBasedVeto/F");
  electronTree_ -> Branch("sumPhotonPt_PUPPI_NoLeptons_MapBasedVeto", &sumPhotonPt_PUPPI_NoLeptons_MapBasedVeto, "sumPhotonPt_PUPPI_NoLeptons_MapBasedVeto/F");

  electronTree_ -> Branch("reliso_ConeVeto_raw", &reliso_ConeVeto_raw, "reliso_ConeVeto_raw/F");
  electronTree_ -> Branch("reliso_MapBasedVeto_raw", &reliso_MapBasedVeto_raw, "reliso_MapBasedVeto_raw/F");

  electronTree_ -> Branch("reliso_ConeVeto_EA", &reliso_ConeVeto_EA, "reliso_ConeVeto_EA/F");
  electronTree_ -> Branch("reliso_MapBasedVeto_EA", &reliso_MapBasedVeto_EA, "reliso_MapBasedVeto_EA/F");

  electronTree_ -> Branch("reliso_PUPPI_ConeVeto", &reliso_PUPPI_ConeVeto, "reliso_PUPPI_ConeVeto/F");
  electronTree_ -> Branch("reliso_PUPPI_NoLeptons_ConeVeto", &reliso_PUPPI_NoLeptons_ConeVeto, "reliso_PUPPI_NoLeptons_ConeVeto/F");
  electronTree_ -> Branch("reliso_PUPPI_ConeVeto_average", &reliso_PUPPI_ConeVeto_average, "reliso_PUPPI_ConeVeto_average/F");

  electronTree_ -> Branch("reliso_PUPPI_MapBasedVeto", &reliso_PUPPI_MapBasedVeto, "reliso_PUPPI_MapBasedVeto/F");
  electronTree_ -> Branch("reliso_PUPPI_NoLeptons_MapBasedVeto", &reliso_PUPPI_NoLeptons_MapBasedVeto, "reliso_PUPPI_NoLeptons_MapBasedVeto/F");
  electronTree_ -> Branch("reliso_PUPPI_MapBasedVeto_average", &reliso_PUPPI_MapBasedVeto_average, "reliso_PUPPI_MapBasedVeto_average/F");
  
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
  Handle <edm::ValueMap <float> > ValueMaps_ChargedHadrons_ConeVeto, ValueMaps_NeutralHadrons_ConeVeto, ValueMaps_Photons_ConeVeto;
  Handle <edm::ValueMap <float> > ValueMaps_ChargedHadrons_MapBasedVeto, ValueMaps_NeutralHadrons_MapBasedVeto, ValueMaps_Photons_MapBasedVeto;
  
  Handle <edm::ValueMap <float> > ValueMaps_PUPPI_ChargedHadrons_ConeVeto, ValueMaps_PUPPI_NeutralHadrons_ConeVeto, ValueMaps_PUPPI_Photons_ConeVeto;
  Handle <edm::ValueMap <float> > ValueMaps_PUPPI_ChargedHadrons_MapBasedVeto, ValueMaps_PUPPI_NeutralHadrons_MapBasedVeto, ValueMaps_PUPPI_Photons_MapBasedVeto;

  Handle <edm::ValueMap <float> > ValueMaps_PUPPI_NoLeptons_ChargedHadrons_ConeVeto, ValueMaps_PUPPI_NoLeptons_NeutralHadrons_ConeVeto, ValueMaps_PUPPI_NoLeptons_Photons_ConeVeto;
  Handle <edm::ValueMap <float> > ValueMaps_PUPPI_NoLeptons_ChargedHadrons_MapBasedVeto, ValueMaps_PUPPI_NoLeptons_NeutralHadrons_MapBasedVeto, ValueMaps_PUPPI_NoLeptons_Photons_MapBasedVeto;
  
  //CITK
  iEvent.getByToken( ValueMaps_ChargedHadrons_ConeVeto_ , ValueMaps_ChargedHadrons_ConeVeto);
  iEvent.getByToken( ValueMaps_ChargedHadrons_MapBasedVeto_, ValueMaps_ChargedHadrons_MapBasedVeto);

  iEvent.getByToken( ValueMaps_NeutralHadrons_ConeVeto_ , ValueMaps_NeutralHadrons_ConeVeto);
  iEvent.getByToken( ValueMaps_NeutralHadrons_MapBasedVeto_ , ValueMaps_NeutralHadrons_MapBasedVeto);

  iEvent.getByToken( ValueMaps_Photons_ConeVeto_ , ValueMaps_Photons_ConeVeto);
  iEvent.getByToken( ValueMaps_Photons_MapBasedVeto_ , ValueMaps_Photons_MapBasedVeto);
  
  //PUPPI 
  iEvent.getByToken( ValueMaps_PUPPI_ChargedHadrons_ConeVeto_ , ValueMaps_PUPPI_ChargedHadrons_ConeVeto);
  iEvent.getByToken( ValueMaps_PUPPI_ChargedHadrons_MapBasedVeto_ , ValueMaps_PUPPI_ChargedHadrons_MapBasedVeto);
  
  iEvent.getByToken( ValueMaps_PUPPI_NeutralHadrons_ConeVeto_ , ValueMaps_PUPPI_NeutralHadrons_ConeVeto);
  iEvent.getByToken( ValueMaps_PUPPI_NeutralHadrons_MapBasedVeto_ , ValueMaps_PUPPI_NeutralHadrons_MapBasedVeto);
  
  iEvent.getByToken( ValueMaps_PUPPI_Photons_ConeVeto_ , ValueMaps_PUPPI_Photons_ConeVeto);
  iEvent.getByToken( ValueMaps_PUPPI_Photons_MapBasedVeto_ , ValueMaps_PUPPI_Photons_MapBasedVeto);
  
  //PUPPI_NoLeptons 
  iEvent.getByToken( ValueMaps_PUPPI_NoLeptons_ChargedHadrons_ConeVeto_ , ValueMaps_PUPPI_NoLeptons_ChargedHadrons_ConeVeto);
  iEvent.getByToken( ValueMaps_PUPPI_NoLeptons_ChargedHadrons_MapBasedVeto_ , ValueMaps_PUPPI_NoLeptons_ChargedHadrons_MapBasedVeto);
  
  iEvent.getByToken( ValueMaps_PUPPI_NoLeptons_NeutralHadrons_ConeVeto_ , ValueMaps_PUPPI_NoLeptons_NeutralHadrons_ConeVeto);
  iEvent.getByToken( ValueMaps_PUPPI_NoLeptons_NeutralHadrons_MapBasedVeto_ , ValueMaps_PUPPI_NoLeptons_NeutralHadrons_MapBasedVeto);
  
  iEvent.getByToken( ValueMaps_PUPPI_NoLeptons_Photons_ConeVeto_ , ValueMaps_PUPPI_NoLeptons_Photons_ConeVeto);
  iEvent.getByToken( ValueMaps_PUPPI_NoLeptons_Photons_MapBasedVeto_ , ValueMaps_PUPPI_NoLeptons_Photons_MapBasedVeto);

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

    float Area = effArea_.getEffectiveArea(std::abs(etaSC_));
    
    // Isolation
    GsfElectron::PflowIsolationVariables pfIso = eleGsfPtr -> pfIsolationVariables();
    isoChargedHadrons_ = pfIso.sumChargedHadronPt ;
    isoNeutralHadrons_ = pfIso.sumNeutralHadronEt ;
    isoPhotons_        = pfIso.sumPhotonEt ;
    isoChargedFromPU_  = pfIso.sumPUPt ;
    
    relIsoWithEA_ = ( pfIso.sumChargedHadronPt + std::max(0.0, (double)pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - rho_ * Area ) )/pt_;
    
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
    sumChargedHadronPt_CITK_ConeVeto =  (*ValueMaps_ChargedHadrons_ConeVeto)[eleGsfPtr];
    sumNeutralHadronPt_CITK_ConeVeto =  (*ValueMaps_NeutralHadrons_ConeVeto)[eleGsfPtr];
    sumPhotonPt_CITK_ConeVeto        =  (*ValueMaps_Photons_ConeVeto)[eleGsfPtr];

    sumChargedHadronPt_CITK_MapBasedVeto =  (*ValueMaps_ChargedHadrons_MapBasedVeto)[eleGsfPtr];
    sumNeutralHadronPt_CITK_MapBasedVeto =  (*ValueMaps_NeutralHadrons_MapBasedVeto)[eleGsfPtr];
    sumPhotonPt_CITK_MapBasedVeto        =  (*ValueMaps_Photons_MapBasedVeto)[eleGsfPtr];
     
    //PUPPI
    sumChargedHadronPt_PUPPI_ConeVeto =  (*ValueMaps_PUPPI_ChargedHadrons_ConeVeto)[elePtr];
    sumNeutralHadronPt_PUPPI_ConeVeto =  (*ValueMaps_PUPPI_NeutralHadrons_ConeVeto)[elePtr];
    sumPhotonPt_PUPPI_ConeVeto        =  (*ValueMaps_PUPPI_Photons_ConeVeto)[elePtr];

    sumChargedHadronPt_PUPPI_MapBasedVeto =  (*ValueMaps_PUPPI_ChargedHadrons_MapBasedVeto)[elePtr];
    sumNeutralHadronPt_PUPPI_MapBasedVeto =  (*ValueMaps_PUPPI_NeutralHadrons_MapBasedVeto)[elePtr];
    sumPhotonPt_PUPPI_MapBasedVeto        =  (*ValueMaps_PUPPI_Photons_MapBasedVeto)[elePtr];
    
    //PUPPINoLeptons
    sumChargedHadronPt_PUPPI_NoLeptons_ConeVeto =  (*ValueMaps_PUPPI_NoLeptons_ChargedHadrons_ConeVeto)[elePtr];
    sumNeutralHadronPt_PUPPI_NoLeptons_ConeVeto =  (*ValueMaps_PUPPI_NoLeptons_NeutralHadrons_ConeVeto)[elePtr];
    sumPhotonPt_PUPPI_NoLeptons_ConeVeto        =  (*ValueMaps_PUPPI_NoLeptons_Photons_ConeVeto)[elePtr];

    sumChargedHadronPt_PUPPI_NoLeptons_MapBasedVeto =  (*ValueMaps_PUPPI_NoLeptons_ChargedHadrons_MapBasedVeto)[elePtr];
    sumNeutralHadronPt_PUPPI_NoLeptons_MapBasedVeto =  (*ValueMaps_PUPPI_NoLeptons_NeutralHadrons_MapBasedVeto)[elePtr];
    sumPhotonPt_PUPPI_NoLeptons_MapBasedVeto        =  (*ValueMaps_PUPPI_NoLeptons_Photons_MapBasedVeto)[elePtr];
    
    reliso_ConeVeto_raw = (sumChargedHadronPt_CITK_ConeVeto + sumNeutralHadronPt_CITK_ConeVeto + sumPhotonPt_CITK_ConeVeto)/pt_;
    reliso_MapBasedVeto_raw = (sumChargedHadronPt_CITK_MapBasedVeto + sumNeutralHadronPt_CITK_MapBasedVeto + sumPhotonPt_CITK_MapBasedVeto)/pt_;

    reliso_ConeVeto_EA = (sumChargedHadronPt_CITK_ConeVeto + std::max(0.,(double)sumNeutralHadronPt_CITK_ConeVeto + sumPhotonPt_CITK_ConeVeto - rho_*Area))/pt_;
    reliso_MapBasedVeto_EA = (sumChargedHadronPt_CITK_MapBasedVeto + std::max(0.,(double)sumNeutralHadronPt_CITK_MapBasedVeto + sumPhotonPt_CITK_MapBasedVeto - rho_*Area))/pt_;
    
    //PUPPI total isolation   
    reliso_PUPPI_ConeVeto = (sumChargedHadronPt_PUPPI_ConeVeto + sumNeutralHadronPt_PUPPI_ConeVeto + sumPhotonPt_PUPPI_ConeVeto)/pt_;
    //PUPPINoLeptons total isolation
    reliso_PUPPI_NoLeptons_ConeVeto = (sumChargedHadronPt_PUPPI_NoLeptons_ConeVeto + sumNeutralHadronPt_PUPPI_NoLeptons_ConeVeto + sumPhotonPt_PUPPI_NoLeptons_ConeVeto)/pt_;

    //reliso puppi average
    reliso_PUPPI_ConeVeto_average = 0.5*(reliso_PUPPI_NoLeptons_ConeVeto + reliso_PUPPI_ConeVeto);

    reliso_PUPPI_MapBasedVeto = (sumChargedHadronPt_PUPPI_MapBasedVeto + sumNeutralHadronPt_PUPPI_MapBasedVeto + sumPhotonPt_PUPPI_MapBasedVeto)/pt_;
    //PUPPINoLeptons total isolation
    reliso_PUPPI_NoLeptons_MapBasedVeto = (sumChargedHadronPt_PUPPI_NoLeptons_MapBasedVeto + sumNeutralHadronPt_PUPPI_NoLeptons_MapBasedVeto + sumPhotonPt_PUPPI_NoLeptons_MapBasedVeto)/pt_;

    //reliso puppi average
    reliso_PUPPI_MapBasedVeto_average = 0.5*(reliso_PUPPI_NoLeptons_MapBasedVeto + reliso_PUPPI_MapBasedVeto);
    
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

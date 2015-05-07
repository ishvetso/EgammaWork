// -*- C++ -*-
//
// Package:    ElectronWork/Test_pt
// Class:      Test_pt
// 
/**\class Test_pt Test_pt.cc ElectronWork/Test_pt/plugins/Test_pt.cc

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
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"

#include "TTree.h"
#include "Math/VectorUtil.h"


class Test_pt : public edm::EDAnalyzer {
public:
  explicit Test_pt(const edm::ParameterSet&);
  ~Test_pt();
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  
  
  // ----------member data ---------------------------
  edm::EDGetTokenT<std::vector<reco::PFCandidate>> puppiToken_;
  edm::EDGetTokenT<std::vector<pat::PackedCandidate> > packedCandidatesToken_;
  
  
  TTree *tree;
  
  Float_t pt_puppi, pt_pf,eta_pf, eta_puppi, phi_pf, phi_puppi;
  
};


//
// static data member definitions
//

//
// constructors and destructor
//
Test_pt::Test_pt(const edm::ParameterSet& iConfig):
  puppiToken_(consumes<std::vector<reco::PFCandidate>>(iConfig.getParameter<edm::InputTag>("puppi"))),
  packedCandidatesToken_(consumes<std::vector<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("pf")))
{

  edm::Service<TFileService> fs;
  tree = fs->make<TTree> ("test", "test");
  
  tree -> Branch("pt_puppi"        ,  &pt_puppi     , "pt_puppi/F");
  tree -> Branch("pt_pf"        ,  &pt_pf     , "pt_pf/F");
  tree -> Branch("eta_puppi"        ,  &eta_puppi     , "eta_puppi/F");
  tree -> Branch("eta_pf"        ,  &eta_pf     , "eta_pf/F");
  tree -> Branch("phi_puppi"        ,  &phi_puppi     , "phi_puppi/F");
  tree -> Branch("phi_pf"        ,  &phi_pf     , "phi_pf/F");
  
}


Test_pt::~Test_pt()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
Test_pt::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace std;
  using namespace edm;
  using namespace reco;
  
  // Pruned particles are the one containing "important" stuff
 
  Handle<std::vector<pat::PackedCandidate>> packedCandidates;
  iEvent.getByToken(packedCandidatesToken_,packedCandidates);
  
 
  Handle<std::vector<reco::PFCandidate>> puppi_Candidates;
  iEvent.getByToken(puppiToken_,puppi_Candidates);
    
  for (unsigned int ipacked = 0; ipacked < packedCandidates -> size(); ipacked ++)
  {
    if((packedCandidates -> at(ipacked)).charge() !=0 ) pt_pf = (packedCandidates -> at(ipacked)).pt();
    else pt_pf = -99.;
    if((puppi_Candidates-> at(ipacked)).charge() !=0 ) pt_puppi = (puppi_Candidates -> at(ipacked)).pt();
    else pt_puppi = -99.;
    
    if((packedCandidates -> at(ipacked)).charge() !=0 ) eta_pf = (packedCandidates -> at(ipacked)).eta();
    else eta_pf = -99.;
    if((puppi_Candidates-> at(ipacked)).charge() !=0 ) eta_puppi = (puppi_Candidates -> at(ipacked)).eta();
    else eta_puppi = -99.;
    
    if((packedCandidates -> at(ipacked)).charge() !=0 ) phi_pf = (packedCandidates -> at(ipacked)).phi();
    else phi_pf = -99.;
    if((puppi_Candidates-> at(ipacked)).charge() !=0 ) phi_puppi = (puppi_Candidates -> at(ipacked)).phi();
    else phi_puppi = -99.;
    
    
    tree -> Fill();
  }
  
  
  
}


// ------------ method called once each job just before starting event loop  ------------
void 
Test_pt::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
Test_pt::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
Test_pt::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
Test_pt::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
Test_pt::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
Test_pt::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/


//define this as a plug-in
DEFINE_FWK_MODULE(Test_pt);

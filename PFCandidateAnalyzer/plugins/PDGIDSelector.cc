#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/UtilAlgos/interface/SingleObjectSelector.h"
#include "CommonTools/UtilAlgos/interface/PdgIdSelector.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

typedef SingleObjectSelector <
          edm::View<pat::PackedCandidate>,
          PdgIdSelector
        > PDGIDSelector;

DEFINE_FWK_MODULE( PDGIDSelector );

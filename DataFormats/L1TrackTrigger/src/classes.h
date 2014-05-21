#include "Rtypes.h" 
#include "Math/Cartesian3D.h" 
#include "Math/Polar3D.h" 
#include "Math/CylindricalEta3D.h" 
#include "Math/PxPyPzE4D.h" 
#include <boost/cstdint.hpp> 

#include "DataFormats/L1TrackTrigger/interface/L1TrackPrimaryVertex.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkEtMissParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkEtMissParticleFwd.h"

#include "DataFormats/L1TrackTrigger/interface/L1TkEmParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkEmParticleFwd.h"

#include "DataFormats/L1TrackTrigger/interface/L1TkElectronParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkElectronParticleFwd.h"

#include "DataFormats/L1TrackTrigger/interface/L1TkJetParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkJetParticleFwd.h"

#include "DataFormats/L1TrackTrigger/interface/L1TkHTMissParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkHTMissParticleFwd.h"

#include "DataFormats/L1TrackTrigger/interface/L1TkMuonParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkMuonParticleFwd.h"

#include "DataFormats/L1TrackTrigger/interface/L1TkTauParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkTauParticleFwd.h"

#include "DataFormats/Common/interface/Wrapper.h"
#include "DataFormats/Common/interface/RefVector.h"
#include "DataFormats/Common/interface/RefProd.h"

namespace {
  struct dictionary {


	// L1 Primary Vertex
     L1TrackPrimaryVertex trzvtx;
     edm::Wrapper<L1TrackPrimaryVertexCollection> trzvtxColl;
     edm::Ref< L1TrackPrimaryVertexCollection > trkvtxRef ;


	// L1TkEtMiss... following L1EtMiss...
     l1extra::L1TkEtMissParticle TketMiss ;
     l1extra::L1TkEtMissParticleCollection TketMissColl ;
     edm::Wrapper<l1extra::L1TkEtMissParticle> w_TketMiss;
     edm::Wrapper<l1extra::L1TkEtMissParticleCollection> w_TketMissColl;
     //l1extra::L1TkEtMissParticleRef refTkEtMiss ;
     //l1extra::L1TkEtMissParticleRefVector refTkVecEtMiss ;
     //l1extra::L1TkEtMissParticleVectorRef vecTkRefEtMissColl ;
     //l1extra::L1TkEtMissParticleRefProd refTkProdEtMiss ;
     //edm::reftobase::Holder<reco::Candidate, l1extra::L1TkEtMissParticleRef> rtbTkm1;
     //edm::reftobase::Holder<reco::Candidate, l1extra::L1TkEtMissParticleRefProd> rtbTkm2;

	// L1TkEmParticle
     l1extra::L1TkEmParticleCollection trkemColl ;
     edm::Wrapper<l1extra::L1TkEmParticleCollection> w_trkemColl;
     l1extra::L1TkEmParticleRef reftrkEm ;
     //l1extra::L1TkEmParticleRefVector refVectrkEmColl ;
     //l1extra::L1TkEmParticleVectorRef vecReftrkEmColl ;
     //edm::reftobase::Holder<reco::Candidate, l1extra::L1TkEmParticleRef> rtbtrke;

        // L1TkElectronParticle
     l1extra::L1TkElectronParticleCollection trkeleColl ;
     edm::Wrapper<l1extra::L1TkElectronParticleCollection> w_trkeleColl;
     l1extra::L1TkElectronParticleRef reftrkEle ;

	// L1TkJetParticle
     l1extra::L1TkJetParticleCollection trkjetColl ;
     edm::Wrapper<l1extra::L1TkJetParticleCollection> w_trkjetColl;
     l1extra::L1TkJetParticleRef reftrkJet ;
     l1extra::L1TkJetParticleRefProd refTkProdJet ;

	// L1TkHTMissParticle
     l1extra::L1TkHTMissParticle TkHTMiss ;
     l1extra::L1TkHTMissParticleCollection TkHTMissColl ;
     edm::Wrapper<l1extra::L1TkHTMissParticle> w_TkHTMiss;
     edm::Wrapper<l1extra::L1TkHTMissParticleCollection> w_TkHTMissColl;

        // L1TkMuonParticle
     l1extra::L1TkMuonParticleCollection trkmuColl ;
     edm::Wrapper<l1extra::L1TkMuonParticleCollection> w_trkmuColl;
     l1extra::L1TkMuonParticleRef reftrkMu ;

        // L1TkTauParticle
     l1extra::L1TkTauParticleCollection trktauColl ;
     edm::Wrapper<l1extra::L1TkTauParticleCollection> w_trktauColl;
     l1extra::L1TkTauParticleRef reftrkTau ;


  };
}

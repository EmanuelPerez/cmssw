// -*- C++ -*-
//
// Package:    L1TrackTriggerObjectsAnalyzer
// Class:      L1TrackTriggerObjectsAnalyzer
// 
/**\class L1TrackTriggerObjectsAnalyzer L1TrackTriggerObjectsAnalyzer.cc SLHCUpgradeSimulations/L1TrackTriggerObjectsAnalyzer/src/L1TrackTriggerObjectsAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Emmanuelle Perez,40 1-A28,+41227671915,
//         Created:  Thu Nov 14 11:22:13 CET 2013
// $Id$
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

#include "DataFormats/Common/interface/Handle.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"


// Gen-level stuff:
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/Candidate/interface/Candidate.h"

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


#include "TFile.h"
#include "TH1F.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

using namespace l1extra;



//
// class declaration
//

class L1TrackTriggerObjectsAnalyzer : public edm::EDAnalyzer {
   public:

  typedef L1TkTrack_PixelDigi_                          L1TkTrackType;
  typedef std::vector< L1TkTrackType >                               L1TkTrackCollectionType;

      explicit L1TrackTriggerObjectsAnalyzer(const edm::ParameterSet&);
      ~L1TrackTriggerObjectsAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      //virtual void endRun(edm::Run const&, edm::EventSetup const&);
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------

	// to test the L1TrackPrimaryVertex :
	edm::InputTag L1VtxInputTag;
	//TH1F* h_zgen;
	//TH1F* h_dz1;
	//TH1F* h_dz2;

	// for L1TrackEtmiss:
	edm::InputTag L1TkEtMissInputTag;

	// for L1TkEmParticles
        edm::InputTag L1TkPhotonsInputTag;
	edm::InputTag L1TkElectronsInputTag;

	// for L1TkJetParticles
        edm::InputTag L1TkJetsInputTag;

	// for L1TkHTMParticle
	edm::InputTag L1TkHTMInputTag;

	// for L1TkMuonParticle
	edm::InputTag L1TkMuonsInputTag;
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
L1TrackTriggerObjectsAnalyzer::L1TrackTriggerObjectsAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

/*
  edm::Service<TFileService> fs;
  int nbins = 25;
  float x1=-25.;
  float x2 = 25.;
   h_zgen = fs -> make<TH1F>("h_zgen",";generated z_{vtx} (cm); Evts",nbins, x1,x2);
   nbins=100;
   x1 = -2;
   x2 = 2;
   h_dz1 = fs -> make<TH1F>("h_dz1",";z_{L1} - z_{gen} (cm); Evts",nbins,x1,x2);
   h_dz2 = fs -> make<TH1F>("h_dz2",";z_{L1} - z_{gen} (cm); Evts",nbins, x1, x2);
*/

  L1VtxInputTag = iConfig.getParameter<edm::InputTag>("L1VtxInputTag") ;
  L1TkEtMissInputTag = iConfig.getParameter<edm::InputTag>("L1TkEtMissInputTag");
  L1TkElectronsInputTag = iConfig.getParameter<edm::InputTag>("L1TkElectronsInputTag");
  L1TkPhotonsInputTag = iConfig.getParameter<edm::InputTag>("L1TkPhotonsInputTag");
  L1TkJetsInputTag = iConfig.getParameter<edm::InputTag>("L1TkJetsInputTag");
  L1TkHTMInputTag = iConfig.getParameter<edm::InputTag>("L1TkHTMInputTag");
  L1TkMuonsInputTag = iConfig.getParameter<edm::InputTag>("L1TkMuonsInputTag");

}


L1TrackTriggerObjectsAnalyzer::~L1TrackTriggerObjectsAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
L1TrackTriggerObjectsAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   std::cout << " ----  a new event ----- " << std::endl;

	// First, retrieve the generated primary vertex


  edm::Handle<edm::HepMCProduct> HepMCEvt;
  iEvent.getByLabel("generator",HepMCEvt);

     const HepMC::GenEvent* MCEvt = HepMCEvt->GetEvent();
     const double mm=0.1;

     float zvtx_gen = -999;
     for ( HepMC::GenEvent::vertex_const_iterator ivertex = MCEvt->vertices_begin(); ivertex != MCEvt->vertices_end(); ++ivertex )
         {
             bool hasParentVertex = false;
 
             // Loop over the parents looking to see if they are coming from a production vertex
             for (
                 HepMC::GenVertex::particle_iterator iparent = (*ivertex)->particles_begin(HepMC::parents);
                 iparent != (*ivertex)->particles_end(HepMC::parents);
                 ++iparent
             )
                 if ( (*iparent)->production_vertex() )
                 {
                     hasParentVertex = true;
                     break;
                }
 
             // Reject those vertices with parent vertices
             if (hasParentVertex) continue;
 
             // Get the position of the vertex
             HepMC::FourVector pos = (*ivertex)->position();
	     zvtx_gen = pos.z()*mm; 

	     break;  // there should be one single primary vertex

          }  // end loop over gen vertices

     //h_zgen -> Fill( zvtx_gen );
     std::cout << " Generated zvertex : " << zvtx_gen << std::endl;



	// ----------------------------------------------------------------------
	//
        // retrieve the L1 vertices
        
 edm::Handle<L1TrackPrimaryVertexCollection> L1VertexHandle;
 iEvent.getByLabel(L1VtxInputTag,L1VertexHandle);
 std::vector<L1TrackPrimaryVertex>::const_iterator vtxIter;
 
 if ( L1VertexHandle.isValid() ) {
     std::cout << " -----  L1TrackPrimaryVertex   ----- " << std::endl;
     int ivtx = 0;
	// several algorithms have been run in the L1TrackPrimaryVertexProducer
	// hence there is a collection of L1 primary vertices.
	// the first one is probably the most reliable.

     for (vtxIter = L1VertexHandle->begin(); vtxIter != L1VertexHandle->end(); ++vtxIter) {
        float z = vtxIter -> getZvertex();
        float sum = vtxIter -> getSum();
        std::cout << " a vertex with  z = sum " << z << " " << sum << std::endl;
        ivtx ++;
        //if (ivtx == 1) h_dz1 -> Fill( z - zvtx_gen) ;
        //if (ivtx == 2) h_dz2 -> Fill( z - zvtx_gen);
     }  
 }

	//
        // ----------------------------------------------------------------------
	// retrieve the EtMiss objects
	//

 edm::Handle<L1TkEtMissParticleCollection> L1TkEtMissHandle;
 iEvent.getByLabel(L1TkEtMissInputTag, L1TkEtMissHandle);
 std::vector<L1TkEtMissParticle>::const_iterator etmIter;

 if (L1TkEtMissHandle.isValid() ) {
    std::cout << " -----  L1TkEtMiss   -----  " << std::endl; 
    for (etmIter = L1TkEtMissHandle -> begin(); etmIter != L1TkEtMissHandle->end(); ++etmIter) {
	float etmis = etmIter -> et();
	const edm::Ref< L1TrackPrimaryVertexCollection > vtxRef = etmIter -> getVtxRef();
	float zvtx = vtxRef -> getZvertex();
        float etMissPU = etmIter -> etMissPU();
	std::cout << " ETmiss = " << etmis << " for zvtx = " << zvtx << " and ETmiss from PU = " << etMissPU << std::endl;
    }
 }


        //
        // ----------------------------------------------------------------------
        // retrieve the L1TkJetParticle objects
        //
        
 edm::Handle<L1TkJetParticleCollection> L1TkJetsHandle;
 iEvent.getByLabel(L1TkJetsInputTag, L1TkJetsHandle);
 std::vector<L1TkJetParticle>::const_iterator jetIter ;

 if ( L1TkJetsHandle.isValid() ) {
    std::cout << " -----   L1TkJetParticle  objects -----  " << std::endl;
    for (jetIter = L1TkJetsHandle -> begin(); jetIter != L1TkJetsHandle->end(); ++jetIter) {
        float et = jetIter -> pt();
        float phi = jetIter -> phi();
        float eta = jetIter -> eta();
        int bx = jetIter -> bx() ;
	float jetvtx = jetIter -> getJetVtx();	// this is the z-vertex of the jet
        const edm::Ref< L1JetParticleCollection > Jetref = jetIter -> getJetRef();
        float et_L1Jet = Jetref -> et();

        std::cout << " a Jet candidate ET eta phi zvertex " << et << " " << eta << " " << phi << " " << jetvtx  << std::endl;
        std::cout << "                           Calo  ET " << et_L1Jet << " bx = " << bx << std::endl;

	// retrieve the L1Tracks associated to the jet :
	std::vector< edm::Ptr< L1TkTrackType > > theseTracks = jetIter -> getTrkPtrs();
	std::cout << "    The tracks associated to the jet and within 2 cm of the jet z-vertex : " << std::endl;
        for (int it=0; it<(int)theseTracks.size(); it++) {
        	edm::Ptr< L1TkTrackType > thisTrk = theseTracks.at(it);
                float tmp_trk_z0   = thisTrk->getVertex().z();
		float dz = fabs( jetvtx - tmp_trk_z0 );
		if (dz > 2. ) continue ;	// keep tracks within 2 cm of the jet vertex
        	float tmp_trk_pt   = thisTrk->getMomentum().perp();
        	float tmp_trk_eta  = thisTrk->getMomentum().eta();
        	float tmp_trk_phi  = thisTrk->getMomentum().phi();
        	float tmp_trk_chi2 = thisTrk->getChi2();
        	std::cout << "     trk pt = " << tmp_trk_pt << " eta = " << tmp_trk_eta << " phi = " << tmp_trk_phi << " z0 = " << tmp_trk_z0 << " chi2 = " << tmp_trk_chi2 << std::endl;
      	}

    }
 }

        //
        // ----------------------------------------------------------------------
        // retrieve HT and HTM
	//

 edm::Handle<L1TkHTMissParticleCollection> L1TkHTMHandle;
 iEvent.getByLabel(L1TkHTMInputTag, L1TkHTMHandle);

 if ( L1TkHTMHandle.isValid() ) {
	std::cout << " -----  L1TkHTMissParticle: size (should be 1) = " << L1TkHTMHandle -> size() << std::endl;
	std::vector<L1TkHTMissParticle>::const_iterator HTMIter = L1TkHTMHandle -> begin();
	float HTT = HTMIter -> EtTotal();
	float HTM = HTMIter -> EtMiss();
	//float HTM_the_same = HTMIter -> et();

	// phi of the HTM vector :
	float phi = HTMIter -> phi();
	std::cout << " HTT = " << HTT << " HTM = " << HTM << " " << "phi(HTM) = " << phi << std::endl;
	 
	// access the L1TkJets used to build HT and HTM :
	const edm::RefProd< L1TkJetParticleCollection > jetCollRef = HTMIter -> getjetCollectionRef();
 	std::vector<L1TkJetParticle>::const_iterator jet = jetCollRef -> begin();
	std::cout << " ET of the first L1TkJet = " << jet -> et() << std::endl;
 }
 else {
    std::cout << L1TkHTMInputTag << " is non valid." << std::endl;
 }


        //
        // ----------------------------------------------------------------------
        // retrieve the L1TkEmParticle objects
	//

 edm::Handle<L1TkEmParticleCollection> L1TkPhotonsHandle;
 iEvent.getByLabel(L1TkPhotonsInputTag, L1TkPhotonsHandle);
 std::vector<L1TkEmParticle>::const_iterator phoIter ;

 if ( L1TkPhotonsHandle.isValid() ) {
    std::cout << " -----   L1TkEmParticle  objects -----  " << std::endl;
    for (phoIter = L1TkPhotonsHandle -> begin(); phoIter != L1TkPhotonsHandle->end(); ++phoIter) {
	float et = phoIter -> pt();
	float phi = phoIter -> phi();
	float eta = phoIter -> eta();
	int bx = phoIter -> bx() ;
        float trkisol = phoIter -> getTrkIsol() ;
	const edm::Ref< L1EmParticleCollection > EGref = phoIter -> getEGRef();
	float et_L1Calo = EGref -> et();
	float eta_calo = EGref -> eta();
	float phi_calo = EGref -> phi();

	std::cout << " a photon candidate ET eta phi trkisol " << et << " " << eta << " " << phi << " " << trkisol << std::endl;
	std::cout << "                Calo  ET eta phi " << et_L1Calo << " " << eta_calo << " " << phi_calo << " bx = " << bx << std::endl; 
    }
 }


        //  
        // ----------------------------------------------------------------------
        // retrieve the L1TkMuons
	//

 edm::Handle<L1TkMuonParticleCollection> L1TkMuonsHandle;
 iEvent.getByLabel(L1TkMuonsInputTag, L1TkMuonsHandle);
 std::vector<L1TkMuonParticle>::const_iterator muIter;

 if ( L1TkMuonsHandle.isValid() ) {
    std::cout << " -----   L1TkMuonPaticle objects ---- " << std::endl;
    for (muIter = L1TkMuonsHandle -> begin(); muIter != L1TkMuonsHandle->end(); ++muIter) {
	float pt = muIter -> pt();
	float eta = muIter -> eta();
	float phi = muIter -> phi();
	float zvtx = muIter -> getTrkzVtx();
	unsigned int quality = muIter -> quality();
	// access the quality via the reference :
	const edm::Ref< L1MuonParticleCollection >	MuRef = muIter -> getMuRef();
	unsigned int qualityBis = MuRef -> gmtMuonCand().quality();
	std::cout << " a muon candidate pt eta phi " << pt << " " << eta << " " << phi << " zvertex = " << zvtx << std::endl;
	std::cout << "   quality (the two qual flags are the same by definition) = " << quality << " " << qualityBis << std::endl;
	
    }
 }

        //
        // ----------------------------------------------------------------------
        // retrieve the L1TkElectronParticle objects
        //

 edm::Handle<L1TkElectronParticleCollection> L1TkElectronsHandle;
 iEvent.getByLabel(L1TkElectronsInputTag, L1TkElectronsHandle);
 std::vector<L1TkElectronParticle>::const_iterator eleIter ;

 if ( L1TkElectronsHandle.isValid() ) {
    std::cout << " -----   L1TkElectronParticle  objects -----  " << std::endl;
    for (eleIter = L1TkElectronsHandle -> begin(); eleIter != L1TkElectronsHandle->end(); ++eleIter) {
        float et = eleIter -> pt();
        float phi = eleIter -> phi();
        float eta = eleIter -> eta();
	int bx = eleIter -> bx();
    	float trkisol = eleIter -> getTrkIsol() ;
	float ztr = eleIter -> getTrkzVtx() ;
        std::cout << "an electron candidate ET eta phi bx trkisol ztr " << et << " " << eta << " " << phi << " " << bx << " " << trkisol << " " << ztr << std::endl;

        const edm::Ref< L1EmParticleCollection > EGref = eleIter -> getEGRef();
        if ( EGref.isNonnull() ) {
           float et_L1Calo = EGref -> et();
           float eta_calo = EGref -> eta();
           float phi_calo = EGref -> phi();
           std::cout << "                Calo  ET eta phi " << et_L1Calo << " " << eta_calo << " " << phi_calo << std::endl;
	}
	else {
	    std::cout << " .... edm::Ref to EGamma is unvalid !! ?  " << std::endl;
	}

        const edm::Ptr< L1TkTrackType > TrkRef = eleIter -> getTrkPtr();
	if ( TrkRef.isNonnull() ) {
            float pt_track = TrkRef -> getMomentum().perp();
            float phi_track = TrkRef -> getMomentum().phi();
            float eta_track = TrkRef -> getMomentum().eta();
            float ztrack = TrkRef -> getVertex().z() ;
            std::cout << "                Track PT eta phi ztr " << pt_track << " " << eta_track << " " << phi_track << " " << ztrack << std::endl;
	}
	else {
	    std::cout << " ... edm::Ptr to L1Tracks is unvalid (e.g. electron was matched to stubs) " << std::endl;
	}
    }
 }


}


// ------------ method called once each job just before starting event loop  ------------
void 
L1TrackTriggerObjectsAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
L1TrackTriggerObjectsAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
L1TrackTriggerObjectsAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
L1TrackTriggerObjectsAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
L1TrackTriggerObjectsAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
L1TrackTriggerObjectsAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
L1TrackTriggerObjectsAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1TrackTriggerObjectsAnalyzer);

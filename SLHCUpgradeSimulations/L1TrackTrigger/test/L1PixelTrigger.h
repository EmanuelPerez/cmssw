#ifndef L1PixelTrigger_h
#define L1PixelTrigger_h

// system include files
#include <cmath>

//// system include files
//#include <memory>
//
//// user include files
//#include "FWCore/Framework/interface/Frameworkfwd.h"
//#include "FWCore/Framework/interface/EDAnalyzer.h"
//
//#include "FWCore/Framework/interface/Event.h"
//#include "FWCore/Framework/interface/MakerMacros.h"
//
//#include "FWCore/ParameterSet/interface/ParameterSet.h"


#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// Geometry
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"

// --- L1Trigger dataFormats
#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1EmParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticleFwd.h"

// SiPixelRecHit dataFormat
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"

// Pileup
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

// --- GenParticles
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"



// --- SimTracks & SimVertex
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

// --- Beam Spot
#include "DataFormats/BeamSpot/interface/BeamSpot.h"



// ROOT includes
#include "TTree.h"

class L1PixelTrigger  : public edm::EDAnalyzer {
 public:
      explicit L1PixelTrigger(const edm::ParameterSet&);
      ~L1PixelTrigger();
 private:
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      int                                      run;
      int                                    event;
      int                                  bunch_n;
      std::vector<int>                      pileup;
      int                                genpart_n;
      std::vector<float>                 genpart_e;
      std::vector<float>                genpart_pt;
      std::vector<float>               genpart_eta;
      std::vector<float>               genpart_phi;
      std::vector<int>              genpart_charge;
      std::vector<int>                  genpart_id;

      float beamspot_x0;
      float beamspot_y0;
      float beamspot_z0;
      float beamspot_x0Error;
      float beamspot_y0Error;
      float beamspot_z0Error;
      float beam_widthX;
      float beam_widthY;
      float beam_sigmaZ;
      float beam_widthXError;
      float beam_widthYError;
      float beam_sigmaZError;



      int simtrk_n;
      std::vector<float> simtrk_pt;
      std::vector<float> simtrk_eta;
      std::vector<float> simtrk_phi;
      std::vector<int> simtrk_id;
      std::vector<int> simtrk_type;
      std::vector<float> simtrk_vx;
      std::vector<float> simtrk_vy;
      std::vector<float> simtrk_vz;
      std::vector<float> sim_vx;
      std::vector<float> sim_vy;
      std::vector<float> sim_vz;



/*
      int                                    tau_n;
      std::vector<float>                     tau_e;
      std::vector<float>                    tau_et;
      std::vector<float>                   tau_eta;
      std::vector<float>                   tau_phi;
*/

      int egammaC_n;
      std::vector<float> egammaC_e;
      std::vector<float> egammaC_et;
      std::vector<float> egammaC_eta;
      std::vector<float> egammaC_phi;
      std::vector<float> egammaC_gx;
      std::vector<float> egammaC_gy;
      std::vector<float> egammaC_gz;
      std::vector<int> egammaC_charge;

      int                                 egamma_n;
      std::vector<float>                  egamma_e;
      std::vector<float>                 egamma_et;
      std::vector<float>                egamma_eta;
      std::vector<float>                egamma_phi;

      std::vector<float>                 egamma_gx;
      std::vector<float>                 egamma_gy;
      std::vector<float>                 egamma_gz;


      int                                   bHit_n;
      std::vector<int>                  bHit_layer;
      std::vector<int>                 bHit_ladder;
      std::vector<float>                   bHit_gx;
      std::vector<float>                   bHit_gy;
      std::vector<float>                   bHit_gz;
      std::vector<int>                    bCl_size;
      std::vector<int>                   bCl_sizex;
      std::vector<int>                   bCl_sizey;
      int                                   fHit_n;
      std::vector<int>                   fHit_disk;
      std::vector<int>                  fHit_blade;
      std::vector<int>                   fHit_side;
      std::vector<float>                   fHit_gx;
      std::vector<float>                   fHit_gy;
      std::vector<float>                   fHit_gz;
      std::vector<int>                    fCl_size;
      std::vector<int>                   fCl_sizex;
      std::vector<int>                   fCl_sizey;
      void InitializeVectors();
      TTree* t;
      edm::InputTag egamma_tag;
      edm::InputTag eGammaCrystal_tag;
      edm::InputTag tau_tag;


      GlobalPoint getCalorimeterPosition(double phi, double eta, double e);
};

L1PixelTrigger::L1PixelTrigger(const edm::ParameterSet& iConfig)
{
       //now do what ever initialization is needed
      tau_tag = iConfig.getParameter<edm::InputTag>("tau");
      egamma_tag = iConfig.getParameter<edm::InputTag>("egamma");
      eGammaCrystal_tag = iConfig.getParameter<edm::InputTag>("egcrystal");
     // eGammaCrystal_tag = iConfig.getParameter<edm::InputTag>("egamma");


      edm::Service<TFileService> fs;
      t = fs->make<TTree>("t","t");
      t->Branch("run",               &run);
      t->Branch("event",             &event);
      t->Branch("bunchN",            &bunch_n);
      t->Branch("pileup",            &pileup);
      t->Branch("genPartN",          &genpart_n);
      t->Branch("genPartE",          &genpart_e);
      t->Branch("genPartPt",         &genpart_pt);
      t->Branch("genPartEta",        &genpart_eta);
      t->Branch("genPartPhi",        &genpart_phi);
      t->Branch("genPartCharge",     &genpart_charge);
      t->Branch("genPartId",         &genpart_id);
/*
      t->Branch("tauN",              &tau_n);
      t->Branch("tauE",              &tau_e);
      t->Branch("tauEt",             &tau_et);
      t->Branch("tauEta",            &tau_eta);
      t->Branch("tauPhi",            &tau_phi);
*/
     t->Branch("simTrkN", &simtrk_n);
     t->Branch("simTrkPt", &simtrk_pt);
     t->Branch("simTrkEta", &simtrk_eta);
     t->Branch("simTrkPhi", &simtrk_phi);
     t->Branch("simTrkId", &simtrk_id);
     t->Branch("simTrkType", &simtrk_type);
     t->Branch("simTrkVx", &simtrk_vx);
     t->Branch("simTrkVy", &simtrk_vy);
     t->Branch("simTrkVz", &simtrk_vz);
     t->Branch("simVx", &sim_vx);
     t->Branch("simVy", &sim_vy);
     t->Branch("simVz", &sim_vz);
/////////// Added on 30 Jan 2015 ////////////////////   
     t->Branch("egCrysN", &egammaC_n);
     t->Branch("egCrysE", &egammaC_e);
     t->Branch("egCrysEt", &egammaC_et);
     t->Branch("egCrysEta", &egammaC_eta);
     t->Branch("egCrysPhi", &egammaC_phi);
     t->Branch("egCrysGx", &egammaC_gx);
     t->Branch("egCrysGy", &egammaC_gy);
     t->Branch("egCrysGz", &egammaC_gz);
     t->Branch("egCrysCharge", &egammaC_charge);
/////////////////////////////////////////////////////

      t->Branch("egN",               &egamma_n);
      t->Branch("egE",               &egamma_e);
      t->Branch("egEt",              &egamma_et);
      t->Branch("egEta",             &egamma_eta);
      t->Branch("egPhi",             &egamma_phi);

      t->Branch("egGx",             &egamma_gx);
      t->Branch("egGy",             &egamma_gy);
      t->Branch("egGz",             &egamma_gz);

      t->Branch("fHitN",             &fHit_n);
      t->Branch("fHitDisk",          &fHit_disk);
      t->Branch("fHitBlade",         &fHit_blade);
      t->Branch("fHitSide",          &fHit_side);
      t->Branch("fHitGx",            &fHit_gx);
      t->Branch("fHitGy",            &fHit_gy);
      t->Branch("fHitGz",            &fHit_gz);
      t->Branch("fClSize",           &fCl_size);
      t->Branch("fClSizeX",          &fCl_sizex);
      t->Branch("fClSizeY",          &fCl_sizey);
      t->Branch("bHitN",             &bHit_n);
      t->Branch("bHitLayer",         &bHit_layer);
      t->Branch("bHitLadder",        &bHit_ladder);
      t->Branch("bHitGx",            &bHit_gx);
      t->Branch("bHitGy",            &bHit_gy);
      t->Branch("bHitGz",            &bHit_gz);
      t->Branch("bClSize",           &bCl_size);
      t->Branch("bClSizeX",          &bCl_sizex);
      t->Branch("bClSizeY",          &bCl_sizey);
}

L1PixelTrigger::~L1PixelTrigger()
{
}
#endif


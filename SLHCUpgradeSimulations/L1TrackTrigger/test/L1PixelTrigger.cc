// -*- C++ -*-
//
// Package:    L1PixelTrigger
// Class:      L1PixelTrigger
// 
/**\class L1PixelTrigger L1PixelTrigger.cc SLHCUpgradeSimulations/L1PixelTrigger/plugins/L1PixelTrigger.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Chang-Seong Moon
//         Created:  Thu, 11 Sep 2014 22:03:15 GMT
// $Id$
//
//

#include "/afs/cern.ch/user/e/eperez/work/L1TrackTrigger/SLC6_620_SLHC12/CMSSW_6_2_0_SLHC12_patch1/src/SLHCUpgradeSimulations/L1TrackTrigger/test/L1PixelTrigger.h"

// ------------ method called for each event  ------------
void L1PixelTrigger::analyze(const edm::Event& e, const edm::EventSetup& es)
{
  InitializeVectors();
  run = e.id().run();
  event = e.id().event();
  ///////////////////////////////////////////////////////////
  // Pileup  
  //////////////////////////////////////////////////////////
  edm::Handle< std::vector<PileupSummaryInfo> > puinfo;
  e.getByLabel( "addPileupInfo", puinfo );
  bunch_n = puinfo->size();
  std::vector<PileupSummaryInfo>::const_iterator PVI;
  for ( PVI = puinfo->begin(); PVI != puinfo->end(); ++PVI ) {
      pileup.push_back(PVI->getPU_NumInteractions());
  }
  ///////////////////////////////////////////////////////////
  // Gen Particles  
  //////////////////////////////////////////////////////////
  edm::Handle< reco::GenParticleCollection > genParticles;
  e.getByLabel( "genParticles", genParticles );
  genpart_n = genParticles->size() ;
  for ( int i = 0; i < genpart_n; ++i ){
      const reco::GenParticle & genParticle = genParticles->at(i);
      genpart_e.push_back(      genParticle.energy() );
      genpart_pt.push_back(     genParticle.pt()     );
      genpart_eta.push_back(    genParticle.eta()    );
      genpart_phi.push_back(    genParticle.phi()    );
      genpart_charge.push_back( genParticle.charge() );
      genpart_id.push_back(     genParticle.pdgId()  );
  }

  ///////////////////////////////////////////////////////////
  // SimTracks & SimVertices
  //////////////////////////////////////////////////////////
  edm::Handle< edm::SimTrackContainer >   simTrackHandle;
  edm::Handle< edm::SimVertexContainer >  simVertex;
  e.getByLabel( "g4SimHits", simTrackHandle );
  e.getByLabel( "g4SimHits", simVertex );
  edm::SimTrackContainer::const_iterator iterSimTracks;
  simtrk_n = 0;
  for (iterSimTracks = simTrackHandle->begin(); iterSimTracks != simTrackHandle->end(); ++iterSimTracks) {
    simtrk_n++;
    simtrk_pt.push_back( iterSimTracks->momentum().pt() );
    simtrk_eta.push_back( iterSimTracks->momentum().eta() );
    simtrk_phi.push_back( iterSimTracks->momentum().phi() );
    simtrk_id.push_back( iterSimTracks->trackId() );
    simtrk_type.push_back( iterSimTracks->type() );
    int index = iterSimTracks->vertIndex();
    simtrk_vx.push_back( simVertex->at(index).position().x() );
    simtrk_vy.push_back( simVertex->at(index).position().y() );
    simtrk_vz.push_back( simVertex->at(index).position().z() );
    // index==0 gets the primary vertex;
    if(index==0){
      sim_vx.push_back( simVertex->at(index).position().x() );
      sim_vy.push_back( simVertex->at(index).position().y() );
      sim_vz.push_back( simVertex->at(index).position().z() );
    }
  }

  ///////////////////////////////////////////////////////////
  // Ecal Single Crystal
  //////////////////////////////////////////////////////////
  edm::Handle< l1extra::L1EmParticleCollection > EgammaC;
  e.getByLabel( eGammaCrystal_tag, EgammaC );
  egammaC_n = EgammaC->size();
  for(l1extra::L1EmParticleCollection::const_iterator it = EgammaC->begin(); it!=EgammaC->end(); ++it){
    egammaC_e.push_back(it->energy());
    egammaC_et.push_back(it->et());
    egammaC_eta.push_back(it->eta());
    egammaC_phi.push_back(it->phi());
    egammaC_charge.push_back(it->charge());
    GlobalPoint pos= getCalorimeterPosition(it->phi(), it->eta(), it->energy());
    egammaC_gx.push_back(pos.x());
    egammaC_gy.push_back(pos.y());
    egammaC_gz.push_back(pos.z());
  }

  ///////////////////////////////////////////////////////////
  // Ecal Trigger Primitives
  //////////////////////////////////////////////////////////
  edm::Handle< l1extra::L1EmParticleCollection > egamma;
  e.getByLabel( egamma_tag, egamma );
  egamma_n = egamma->size();
  l1extra::L1EmParticleCollection::const_iterator egIt = egamma->begin();


  for ( ; egIt!=egamma->end(); ++egIt ){

      egamma_e.push_back(   egIt->energy() );
      egamma_et.push_back(  egIt->et()     );
      egamma_eta.push_back( egIt->eta()    );
      egamma_phi.push_back( egIt->phi()    );
      GlobalPoint pos= getCalorimeterPosition(egIt->phi(), egIt->eta(), egIt->energy());
      egamma_gx.push_back(pos.x());
      egamma_gy.push_back(pos.y());
      egamma_gz.push_back(pos.z());

  }
  //////////////////////////////////////////////////////////
  // Geometry record 
  //////////////////////////////////////////////////////////
  edm::ESHandle<TrackerGeometry> geom;
  es.get<TrackerDigiGeometryRecord>().get(geom);
  //////////////////////////////////////////////////////////
  // RecHits 
  //////////////////////////////////////////////////////////
  edm::Handle<SiPixelRecHitCollection> recHits;
  e.getByLabel( "siPixelRecHits", recHits );
  bHit_n = 0;
  fHit_n = 0;
  SiPixelRecHitCollection::const_iterator detUnitIt    = recHits->begin();
  SiPixelRecHitCollection::const_iterator detUnitItEnd = recHits->end();
  for ( ; detUnitIt != detUnitItEnd; detUnitIt++ ) {
      DetId detId = DetId(detUnitIt->detId()); 
      int subid = detId.subdetId();
      SiPixelRecHitCollection::DetSet::const_iterator recHitIt    = detUnitIt->begin();
      SiPixelRecHitCollection::DetSet::const_iterator recHitItEnd = detUnitIt->end();
      for ( ; recHitIt != recHitItEnd; ++recHitIt) {
          LocalPoint  lp = recHitIt->localPosition();
          GlobalPoint gp = ( (geom.product())->idToDet(detId) )->surface().toGlobal(lp);
          SiPixelRecHit::ClusterRef const& Cluster = recHitIt->cluster();
          if ( gp.perp() < 20 ) { // drop outer tracker 
                if ( subid==PixelSubdetector::PixelEndcap ){
                     fHit_n ++;
                     fHit_disk.push_back(  PXFDetId(detId).disk()  );
                     fHit_blade.push_back( PXFDetId(detId).blade() );
                     fHit_side.push_back(  PXFDetId(detId).side()  );
                     fHit_gx.push_back(    gp.x() );
                     fHit_gy.push_back(    gp.y() );
                     fHit_gz.push_back(    gp.z() );
                     fCl_size.push_back(   Cluster->size()  );
                     fCl_sizex.push_back(  Cluster->sizeX() );
                     fCl_sizey.push_back(  Cluster->sizeY() );
                }
                if ( subid==PixelSubdetector::PixelBarrel ){
                     bHit_n ++;
                     bHit_layer.push_back(  PXBDetId(detId).layer() ); 
                     bHit_ladder.push_back( PXBDetId(detId).ladder() ); 
                     bHit_gx.push_back(     gp.x() );
                     bHit_gy.push_back(     gp.y() );
                     bHit_gz.push_back(     gp.z() );
                     bCl_size.push_back(    Cluster->size()  );
                     bCl_sizex.push_back(   Cluster->sizeX() );
                     bCl_sizey.push_back(   Cluster->sizeY() );
                }
          } 
      } // close recHits loop
  } // close detUnits loop

  //////////////////////////////////////////////////////////
  // Fill Tree 
  //////////////////////////////////////////////////////////
  t->Fill();

}// close L1PixelTrigger::analyze

GlobalPoint L1PixelTrigger::getCalorimeterPosition(double phi, double eta, double e) {
  double x = 0;
  double y = 0;
  double z = 0;
  double depth = 0.89*(7.7+ log(e) );
  double theta = 2*atan(exp(-1*eta));
  double r = 0;
  if( fabs(eta) > 1.479 )
    {
      double ecalZ = 315.4*fabs(eta)/eta;

      r = ecalZ / cos( 2*atan( exp( -1*eta ) ) ) + depth;
      x = r * cos( phi ) * sin( theta );
      y = r * sin( phi ) * sin( theta );
      z = r * cos( theta );
    }
  else
    {
      double rperp = 129.0;
      double zface = sqrt( cos( theta ) * cos( theta ) /
                            ( 1 - cos( theta ) * cos( theta ) ) *
                            rperp * rperp ) * fabs( eta ) / eta;
      r = sqrt( rperp * rperp + zface * zface ) + depth;
      x = r * cos( phi ) * sin( theta );
      y = r * sin( phi ) * sin( theta );
      z = r * cos( theta );
    }
  GlobalPoint pos(x,y,z);
  return pos;
}


void L1PixelTrigger::InitializeVectors()
{
             pileup.clear();
          genpart_e.clear();
         genpart_pt.clear();
        genpart_eta.clear();
        genpart_phi.clear();
     genpart_charge.clear();
         genpart_id.clear();
/*
              tau_e.clear();
             tau_et.clear();
            tau_eta.clear();
            tau_phi.clear();
*/

          simtrk_pt.clear();
         simtrk_eta.clear();
         simtrk_phi.clear();
          simtrk_id.clear();
        simtrk_type.clear();
          simtrk_vx.clear();
          simtrk_vy.clear();
          simtrk_vz.clear();
             sim_vx.clear();
             sim_vy.clear();
             sim_vz.clear();

          egammaC_e.clear();
         egammaC_et.clear();
        egammaC_eta.clear();
        egammaC_phi.clear();
         egammaC_gx.clear();
         egammaC_gy.clear();
         egammaC_gz.clear();
     egammaC_charge.clear();



           egamma_e.clear();
          egamma_et.clear();
         egamma_eta.clear();
         egamma_gx.clear();
         egamma_gy.clear();
         egamma_gz.clear();
         egamma_phi.clear();
          fHit_disk.clear();
         fHit_blade.clear();
          fHit_side.clear();
            fHit_gx.clear();
            fHit_gy.clear();
            fHit_gz.clear();
           fCl_size.clear();
          fCl_sizex.clear();
          fCl_sizey.clear();
         bHit_layer.clear();
        bHit_ladder.clear();
            bHit_gx.clear();
            bHit_gy.clear();
            bHit_gz.clear();
           bCl_size.clear();
          bCl_sizex.clear();
          bCl_sizey.clear();
}
//define this as a plug-in
DEFINE_FWK_MODULE(L1PixelTrigger);


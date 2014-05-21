// A.Savin University of Wisconsin, Version 1.0 //
// A.Savin University of Wisconsin, Version 2.0 //

#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1EmParticleFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalBarrelGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalEndcapGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include <iostream>

#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "FastSimulation/CaloGeometryTools/interface/CaloGeometryHelper.h"
#include "SimDataFormats/SLHC/interface/L1EGCrystalCluster.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"

//here
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "SimDataFormats/CaloTest/interface/HcalTestNumbering.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"


#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalRecHit/interface/HcalSourcePositionData.h"


// Numbering scheme
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include <TVector3.h>


class L1EGCrystalClusterProducer : public edm::EDProducer {
public:

  explicit L1EGCrystalClusterProducer(const edm::ParameterSet&);

private:
  virtual void produce(edm::Event&, const edm::EventSetup&);
  unsigned int getMatchedClusterIndex(l1slhc::L1EGCrystalCluster& egxtals, float& eta, float& phi, float& dr_min);
//  virtual void beginRun(const edm::EventSetup & c) ;

  CaloGeometryHelper myGeometry;

  l1slhc::L1EGCrystalCluster ecalhits[100000] ;

  l1slhc::L1EGCrystalCluster hcalhits[100000] ;

  l1slhc::L1EGCrystalCluster clusters ;
 
  bool First = true ;
  bool DEBUG;

  // Concrete collection of output objects (with extra tuning information)
//  typedef struct L1EGCrystCluster TrigCrystalCollection;
//  typedef std::vector<TrigCrystalCollection> TrigCrystalCollectionPtr;

};

namespace objects_ordering {
  class L1EmETComparator {
  public:
    bool operator()(const l1extra::L1EmParticle a, const l1extra::L1EmParticle b) const {
      double et_a = 0.0;
      double et_b = 0.0;
      if (cosh(a.eta()) > 0.0) et_a = a.energy()/cosh(a.eta());
      if (cosh(b.eta()) > 0.0) et_b = b.energy()/cosh(b.eta());

      return et_a > et_b;
    }
  };
  class ClusterETComparator {
  public:
    bool operator()(const l1slhc::L1EGCrystalCluster a, const l1slhc::L1EGCrystalCluster b) const {
      return a.et > b.et;
    }
  };
}


L1EGCrystalClusterProducer::L1EGCrystalClusterProducer(const edm::ParameterSet& iConfig) 
{
  produces<l1slhc::L1EGCrystalClusterCollection>( "EGCrystalCluster" ) ;
  produces < l1extra::L1EmParticleCollection > ( "EGammaCrystal" );

  DEBUG = iConfig.getParameter<bool>("DEBUG");
}



/*
void L1EGCrystalClusterProducer::beginRun(const edm::EventSetup & c){

std::cout << " I am in beginJob " << std::endl ;

  edm::ESHandle<CaloTopology> theCaloTopology;
  c.get<CaloTopologyRecord>().get(theCaloTopology);       
  edm::ESHandle<CaloGeometry> pG;
  c.get<CaloGeometryRecord>().get(pG);     
  // Setup the tools
  double bField000 = 4.;
  myGeometry.setupGeometry(*pG);
  myGeometry.setupTopology(*theCaloTopology);
  myGeometry.initialize(bField000);

std::cout << " I am out of beginJob " << std::endl ;
    
}
*/

void  L1EGCrystalClusterProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

 if (DEBUG)  std::cout << " enter in produce " << std::endl ;

  if (First) {
  edm::ESHandle<CaloTopology> theCaloTopology;
  iSetup.get<CaloTopologyRecord>().get(theCaloTopology);       
  edm::ESHandle<CaloGeometry> pG;
  iSetup.get<CaloGeometryRecord>().get(pG);     
  // Setup the tools
  double bField000 = 4.;
  myGeometry.setupGeometry(*pG);
  myGeometry.setupTopology(*theCaloTopology);
  myGeometry.initialize(bField000);
  First = false ;
  }

  //  std::cout << " I am in producer " << std::endl ;


  std::auto_ptr<l1slhc::L1EGCrystalClusterCollection> trigCrystalClusters (new l1slhc::L1EGCrystalClusterCollection()) ;

/*
  // Lets check Extras
  edm::Handle<l1extra::L1EmParticleCollection> eGammaCollection;
  iEvent.getByLabel("SLHCL1ExtraParticles","EGamma", eGammaCollection);
  l1extra::L1EmParticleCollection egColl = (*eGammaCollection.product());
  std::sort(egColl.begin(), egColl.end(), objects_ordering::L1EmETComparator());
  std::vector<unsigned int> matched_indices;

  for ( size_t i = 0; i < egColl.size(); ++i) {
    float eta = egColl[i].eta();
    float et = 0;
    if (cosh(eta) > 0.0) et = egColl[i].energy()/cosh(eta);
    std::cout << "et, e, eta, phi  of EGamma Object #" << i << " : "
              << et << " "
              << egColl[i].energy() << " "
              << egColl[i].eta() << " "
              << egColl[i].phi() << std::endl;
   }
*/

std::auto_ptr < l1extra::L1EmParticleCollection > l1EGammaCrystal( new l1extra::L1EmParticleCollection );

  // Retrieve the SimHits.
  int nHits = 0 ;
  	// Sasha's FAMOS :
  //edm::Handle<edm::PCaloHitContainer> pcalohits;
  //iEvent.getByLabel("famosSimHits","EcalHitsEB",pcalohits);
	// using RecHits :
  edm::Handle<EcalRecHitCollection> pcalohits;
  iEvent.getByLabel("ecalRecHit","EcalRecHitsEB",pcalohits);
	// Geant's pcaloHits :
////  edm::Handle<edm::PCaloHitContainer> pcalohits;
////  iEvent.getByLabel("g4SimHits","EcalHitsEB",pcalohits);

  std::vector<float> simHitsBarrel;
  simHitsBarrel.resize(62000,0.);
////  edm::PCaloHitContainer::const_iterator it=pcalohits.product()->begin();
////  edm::PCaloHitContainer::const_iterator itend=pcalohits.product()->end();
  EcalRecHitCollection::const_iterator it=pcalohits.product()->begin();
  EcalRecHitCollection::const_iterator itend=pcalohits.product()->end();

  for(;it!=itend;++it)
    {
      simHitsBarrel[EBDetId(it->id()).hashedIndex()]+=it->energy();
      if(it->energy() > 0.2)
	{
          EBDetId EBid = EBDetId(it->id()) ;
	  double theta =  myGeometry.getEcalBarrelGeometry()->getGeometry(it->id())->getPosition().theta() ;
	  double phi =  myGeometry.getEcalBarrelGeometry()->getGeometry(it->id())->getPosition().phi() ;
	  double eta = -1.*log(tan(theta/2.)) ;
	  ecalhits[nHits].et = (it->energy())*sin(theta) ;
	  ecalhits[nHits].eta = eta ;
	  ecalhits[nHits].phi = phi ;
	  ecalhits[nHits].e = (it->energy()) ;
	  ecalhits[nHits].x = myGeometry.getEcalBarrelGeometry()->getGeometry(it->id())->getPosition().x() ;
	  ecalhits[nHits].y = myGeometry.getEcalBarrelGeometry()->getGeometry(it->id())->getPosition().y() ;
	  ecalhits[nHits].z = myGeometry.getEcalBarrelGeometry()->getGeometry(it->id())->getPosition().z() ;
	  nHits++ ;
	  if (DEBUG) std::cout << " EB Hits " <<  it->energy() <<  " phi " << phi << " eta " << eta << " eta1 " << EBid.ieta() << " theta " << theta << std::endl;
	}
    }
  std::vector<float> simHitsEndcap;
  simHitsEndcap.resize(20000,0.);
  //iEvent.getByLabel("famosSimHits","EcalHitsEE",pcalohits);
  iEvent.getByLabel("ecalRecHit","EcalRecHitsEE",pcalohits);
////  iEvent.getByLabel("g4SimHits","EcalHitsEE",pcalohits);
  it=pcalohits.product()->begin();
  itend=pcalohits.product()->end();
  for(;it!=itend;++it)
    {
      simHitsEndcap[EEDetId(it->id()).hashedIndex()]+=it->energy();
      if(it->energy() > 0.2)
	{
	  double theta =  myGeometry.getEcalEndcapGeometry()->getGeometry(it->id())->getPosition().theta() ;
	  double phi =  myGeometry.getEcalEndcapGeometry()->getGeometry(it->id())->getPosition().phi() ;
	  double eta = -1.*log(tan(theta/2.)) ;
	  ecalhits[nHits].et = (it->energy())*sin(theta) ;
	  ecalhits[nHits].eta = eta ;
	  ecalhits[nHits].phi = phi ;
	  ecalhits[nHits].e = (it->energy()) ;
	  ecalhits[nHits].x = myGeometry.getEcalEndcapGeometry()->getGeometry(it->id())->getPosition().x() ;
	  ecalhits[nHits].y = myGeometry.getEcalEndcapGeometry()->getGeometry(it->id())->getPosition().y() ;
	  ecalhits[nHits].z = myGeometry.getEcalEndcapGeometry()->getGeometry(it->id())->getPosition().z() ;
	  nHits++ ;
	  if (DEBUG) std::cout << " EE Hits " << " energy " << it->energy() << " phi " << phi << " eta " << eta <<  std::endl;
	}
    }
  
//  iEvent.getByLabel("g4SimHits","HcalHits",pcalohits);
//  it=pcalohits.product()->begin();
//  itend=pcalohits.product()->end();
  
  edm::ESHandle<CaloGeometry> pG1;
  iSetup.get<CaloGeometryRecord>().get(pG1);
  const CaloGeometry* geometry = pG1.product();
  
  int nHitsHCAL = 0 ;
/*  
  for(;it!=itend;++it)
    {
      HcalDetId myDetId(it->id());
      
      const CaloCellGeometry* cellGeometry =
	geometry->getSubdetectorGeometry(myDetId)->getGeometry(myDetId) ;
      if ( it->energy() > 0.1 )
	{ 
	  hcalhits[nHitsHCAL].e = (it->energy()) ;
	  hcalhits[nHitsHCAL].eta = cellGeometry->getPosition().eta() ;
	  hcalhits[nHitsHCAL].phi = cellGeometry->getPosition().phi() ;
	  //	  if(DEBUG) std::cout << " id " << myDetId << " Energy " << it->energy() ;
	}
      
      nHitsHCAL++ ;
      
    }
*/  

     edm::Handle<HBHERecHitCollection> hbhecoll;
     iEvent.getByLabel("hbheprereco", hbhecoll);

     for (HBHERecHitCollection::const_iterator j=hbhecoll->begin(); j != hbhecoll->end(); j++) {
       HcalDetId cell(j->id());
       const CaloCellGeometry* cellGeometry =
         geometry->getSubdetectorGeometry (cell)->getGeometry (cell) ;
      if ( j->energy() > 0.1 )
      {
      hcalhits[nHitsHCAL].e = (j->energy()) ;      
      hcalhits[nHitsHCAL].eta = cellGeometry->getPosition().eta() ;
      hcalhits[nHitsHCAL].phi = cellGeometry->getPosition().phi() ;
      if(DEBUG && hcalhits[nHitsHCAL].e > 10) std::cout << " id " << cell << " Energy " << j->energy() << " eta " << hcalhits[nHitsHCAL]
.eta << " phi " << hcalhits[nHitsHCAL].phi <<  std::endl ;
      }
      nHitsHCAL++ ;
    }

  
  float Etmax = 0 ;
  while( Etmax != 2.)
    {
      Etmax = 2. ;
      int   centerId = 0 ;
      for(int k=0; k<nHits; k++)
	{
//	  if(ecalhits[k].e > Emax)
	  if(ecalhits[k].et > Etmax)
	    {
	      Etmax = ecalhits[k].et ;
	      centerId = k ;
	    }
	} 

      //  std::cout << " I am here " << std::endl ;


      if( Etmax != 2. )
	{
	  float Total_E = ecalhits[centerId].e ;
//	   std::cout << " E center " << ecalhits[centerId].e << " phi " << ecalhits[centerId].phi << " eta " << ecalhits[centerId].eta <<  std::endl ;
	  float Weightedx = ecalhits[centerId].x*ecalhits[centerId].e ;
	  float Weightedy = ecalhits[centerId].y*ecalhits[centerId].e ;
	  float Weightedz = ecalhits[centerId].z*ecalhits[centerId].e ;
	  for(int k=0; k<nHits; k++)
	    {
	      if ( k != centerId && fabs(ecalhits[centerId].eta-ecalhits[k].eta) < 0.08 && fabs(ecalhits[centerId].phi-ecalhits[k].phi) < 0.1 )
		{
//	  std::cout << " E to add " << ecalhits[k].e << " phi " << ecalhits[k].phi << " eta " << ecalhits[k].eta <<  std::endl ;
		  Total_E = Total_E + ecalhits[k].e ;
		  Weightedx = Weightedx + ecalhits[k].x*ecalhits[k].e;
		  Weightedy = Weightedy + ecalhits[k].y*ecalhits[k].e;
		  Weightedz = Weightedz + ecalhits[k].z*ecalhits[k].e;
		  TVector3 tmp1(ecalhits[k].x,ecalhits[k].y,ecalhits[k].z) ;
		  //   		std::cout << " tmp1 " << " phi " << tmp1.Phi() << " eta " << tmp1.PseudoRapidity() <<  std::endl ;
		  ecalhits[k].et = 0 ;
		}
	    }
	  ecalhits[centerId].et = 0 ;
	  clusters.e = Total_E ;
	  clusters.x = Weightedx/Total_E ;
	  clusters.y = Weightedy/Total_E ;
	  clusters.z = Weightedz/Total_E ;
	  TVector3 tmp(clusters.x,clusters.y,clusters.z) ;
	  clusters.phi = tmp.Phi() ;
	  clusters.eta = tmp.PseudoRapidity() ;
	  clusters.et = Total_E*sin(tmp.Theta()) ;
	  if (DEBUG) {
	    //	    std::cout << " clusters " << Total_E << " et " << clusters.et << " phi " << clusters.phi << " eta " << clusters.eta <<  std::endl ;
	  }

 
           Total_E = 0. ;
           Weightedx = 0 ;
           Weightedy = 0 ;
           Weightedz = 0 ;
           float Total_E1 = 0. ;
           float Weightedx1 = 0. ;
           float Weightedy1 = 0. ;
           float Weightedz1 = 0. ;
 	  for(int k=0; k<nHits; k++)
 	    {
 	      if ( ecalhits[k].et > 0.05 && fabs(ecalhits[centerId].eta-ecalhits[k].eta) < 0.25 && fabs(ecalhits[centerId].phi-ecalhits[k].phi) < 0.25 )
 		{
 		  Total_E1 = Total_E1 + ecalhits[k].e ;
 		  Weightedx1 = Weightedx1 + ecalhits[k].x*ecalhits[k].e;
 		  Weightedy1 = Weightedy1 + ecalhits[k].y*ecalhits[k].e;
 		  Weightedz1 = Weightedz1 + ecalhits[k].z*ecalhits[k].e;
 		}
 	      if ( ecalhits[k].et > 0.05 && ecalhits[k].et < 5. && fabs(ecalhits[centerId].eta-ecalhits[k].eta) < 0.12 && fabs(ecalhits[centerId].phi-ecalhits[k].phi) < 1. )
 		{
 		  Total_E = Total_E + ecalhits[k].e ;
 		  Weightedx = Weightedx + ecalhits[k].x*ecalhits[k].e;
 		  Weightedy = Weightedy + ecalhits[k].y*ecalhits[k].e;
 		  Weightedz = Weightedz + ecalhits[k].z*ecalhits[k].e;
 		}
 	    }
 	  TVector3 tmp2(Weightedx1/Total_E1,Weightedy1/Total_E1,Weightedz1/Total_E1) ;
 	  clusters.ECALiso = Total_E1*sin(tmp2.Theta())/clusters.et ;
 
 	  TVector3 tmp1(Weightedx/Total_E,Weightedy/Total_E,Weightedz/Total_E) ;
 	  clusters.ECALetPUcorr = clusters.et-Total_E*sin(tmp1.Theta())/19. ;
 
	  
	  Total_E = 0 ;
	  
	  for(int k=0; k<nHitsHCAL; k++)
	    {
	      if ( fabs(ecalhits[centerId].eta-hcalhits[k].eta) < 0.15 && fabs(ecalhits[centerId].phi-hcalhits[k].phi) < 0.15 )
		{
		  Total_E = Total_E + hcalhits[k].e ;
		}
	    }


	  clusters.hovere = Total_E/clusters.e ;

	  
	  trigCrystalClusters->push_back(clusters) ;

/* add l1extra 
   match with existing and the next protion just to create new 

    l1extra::L1EmParticle* em_new = new l1extra::L1EmParticle ;
    for ( l1extra::L1EmParticleCollection::const_iterator lIt =egColl.begin(); lIt != egColl.end(); ++lIt ) {
      l1extra::L1EmParticle em( *lIt );
      float eta = em.eta();
      float phi = em.phi();
      float deltaR ;
      float deltaRmin = 0.3 ;
      unsigned int indx = getMatchedClusterIndex(clusters, eta, phi, deltaR);
      if (indx != 9999 && deltaR < deltaRmin) {
        deltaRmin = deltaR ;
//        l1extra::L1EmParticle* em_new = em.clone();
        em_new = em.clone();
        float eta_new = clusters.eta;
        float phi_new = clusters.phi;
        float pt_new  = clusters.et;
        reco::Candidate::PolarLorentzVector p4(pt_new, eta_new, phi_new, em.mass());
        em_new->setP4(p4);
    std::cout << " Old and New " << " old eta " << eta << " new " << em_new->eta() << " old phi " << phi << " new " << em_new->phi() << std::endl ; 
        }
      }
        l1EGammaCrystal->push_back( *em_new );
*/
//        std::cout << " A cluster: et " << clusters.et << " eta " << clusters.eta << " phi " << clusters.phi << " hovere "  << clusters.hovere << " iso  " << clusters.ECALiso <<  std::endl ;
        if (clusters.hovere < 1. && clusters.ECALiso < 2.){
        reco::Candidate::PolarLorentzVector p4(clusters.et, clusters.eta, clusters.phi, 0.);
        edm::Ref< L1GctEmCandCollection > ref();
	l1EGammaCrystal->push_back(l1extra::L1EmParticle( p4,edm::Ref< L1GctEmCandCollection >(),0) );
	//        l1EGammaCrystal->push_back(l1extra::L1EmParticle( p4));
        // std::cout << " Put in a cluster " << std::endl ;
           }
	}
    }
  
  if (DEBUG) {
    //    std::cout << " writing in " << trigCrystalClusters->size() <<   std::endl ;
    for (size_t i = 0; i < trigCrystalClusters->size(); ++i) {
    //  std::cout << "hovere of cluster # " << i << " is " << trigCrystalClusters->at(i).hovere << std::endl;
    }
  }
  
  iEvent.put(trigCrystalClusters,"EGCrystalCluster");
  iEvent.put(l1EGammaCrystal, "EGammaCrystal" );
  
}

unsigned int L1EGCrystalClusterProducer::getMatchedClusterIndex(l1slhc::L1EGCrystalCluster& egxtals, float& eta, float& phi, float& dr_min) {
  dr_min = 999.9;
  size_t index_min = 9999;
  float dr = deltaR(egxtals.eta, egxtals.phi, eta, phi);
  if (dr < dr_min) {
      index_min = 0;
      dr_min = dr;
    }
  return index_min; 
}


DEFINE_FWK_MODULE(L1EGCrystalClusterProducer);

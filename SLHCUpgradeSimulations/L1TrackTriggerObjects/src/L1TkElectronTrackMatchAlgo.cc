// -*- C++ -*-
//
// Package:    L1TrackTriggerObjects
// Class:      L1TkElectronTrackMatchAlgo
// 
/**\class L1TkElectronTrackMatchAlgo 

 Description: Algorithm to match L1EGamma oject with L1Track candidates

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  S. Dutta and A. Modak
//         Created:  Wed Dec 4 12 11:55:35 CET 2013
// $Id$
//
//


// system include files
#include <memory>
#include <cmath>

#include "DataFormats/Math/interface/deltaPhi.h"
#include "SLHCUpgradeSimulations/L1TrackTriggerObjects/interface/L1TkElectronTrackMatchAlgo.h"
namespace L1TkElectronTrackMatchAlgo {
  // ------------ match EGamma and Track
  void doMatch(l1extra::L1EmParticleCollection::const_iterator egIter, L1TkTrackCollectionType::const_iterator trkIter, double& dph, float&  dr, float& deta, float& dphiprime) {
    GlobalPoint egPos = L1TkElectronTrackMatchAlgo::calorimeterPosition(egIter->phi(), egIter->eta(), egIter->energy());
    dph  = L1TkElectronTrackMatchAlgo::deltaPhi(egPos, trkIter);
    dr   = L1TkElectronTrackMatchAlgo::deltaR(egPos, trkIter);
    deta = L1TkElectronTrackMatchAlgo::deltaEta(egPos, trkIter);

    dphiprime = L1TkElectronTrackMatchAlgo::deltaPhiPrime(egPos, trkIter);
    //dphiplus = L1TkElectronTrackMatchAlgo::deltaPhiPlus(egPos, trkIter);
    //dphiminus  = L1TkElectronTrackMatchAlgo::deltaPhiMinus(egPos, trkIter);
  }
  // --------------- calculate deltaR between Track and EGamma object
double deltaPhi(GlobalPoint epos, L1TkTrackCollectionType::const_iterator trkIter){
    double er = epos.perp();

    // Using track fit curvature
    //  double curv = 0.003 * magnetStrength * trk->getCharge()/ trk->getMomentum().perp(); 
    double curv = trkIter->getRInv();
    double x1 = (asin(er*curv/(2.0)));
    double phi1 = reco::deltaPhi(trkIter->getMomentum().phi(), epos.phi());

    double dif1 = phi1 - x1;
    double dif2 = phi1 + x1; 

    if (fabs(dif1) < fabs(dif2)) return dif1;
    else return dif2; 
  
  }

double deltaPhiPrime(GlobalPoint epos, L1TkTrackCollectionType::const_iterator trkIter){
    double er = epos.perp();

    // Using track fit curvature
    //  double curv = 0.003 * magnetStrength * trk->getCharge()/ trk->getMomentum().perp(); 
    double curv = trkIter->getRInv();
    double x1 = (asin(er*curv/(2.0)));
    double phi1 = reco::deltaPhi(trkIter->getMomentum().phi(), epos.phi());

    double dif1 = phi1 - x1;
    double dif2 = phi1 + x1;

    double rinv = trkIter -> getRInv();
    if (rinv < 0) return dif1;
    return dif2;

  }


double deltaPhiPlus(GlobalPoint epos, L1TkTrackCollectionType::const_iterator trkIter){
    double er = epos.perp();

    // Using track fit curvature
    //  double curv = 0.003 * magnetStrength * trk->getCharge()/ trk->getMomentum().perp(); 
    double curv = trkIter->getRInv();
    double x1 = (asin(er*curv/(2.0)));
    double phi1 = reco::deltaPhi(trkIter->getMomentum().phi(), epos.phi());

    double dif2 = phi1 + x1;
    return dif2;
}

double deltaPhiMinus(GlobalPoint epos, L1TkTrackCollectionType::const_iterator trkIter){
    double er = epos.perp();

    // Using track fit curvature
    //  double curv = 0.003 * magnetStrength * trk->getCharge()/ trk->getMomentum().perp(); 
    double curv = trkIter->getRInv();
    double x1 = (asin(er*curv/(2.0)));
    double phi1 = reco::deltaPhi(trkIter->getMomentum().phi(), epos.phi());

    double dif1 = phi1 - x1;
    return dif1;
}


// --------------- calculate deltaPhi between Track and EGamma object                 
float deltaR(GlobalPoint epos, L1TkTrackCollectionType::const_iterator trkIter){
    float dPhi = fabs(reco::deltaPhi(epos.phi(), trkIter->getMomentum().phi()));
    float dEta = (epos.eta() - trkIter->getMomentum().eta());
    return sqrt(dPhi*dPhi + dEta*dEta);
  }
  // --------------- calculate deltaEta between Track and EGamma object                 
float deltaEta(GlobalPoint epos, L1TkTrackCollectionType::const_iterator trkIter){
    float corr_eta = 999.0;
    float er = epos.perp();
    float ez = epos.z();
    //float z0 = trkIter->getMomentum().z();
    float z0 = trkIter->getVertex().z() ;
    float theta = 0.0;
    if (ez >= 0) theta = atan(er/fabs(ez-z0));
    else theta = M_PI - atan(er/fabs(ez-z0));
    corr_eta = -1.0 * log(tan(theta/2.0));
    float deleta = (corr_eta - trkIter->getMomentum().eta());
    return deleta;
  }
  // -------------- get Calorimeter position
  GlobalPoint calorimeterPosition(float phi, float eta, float e) {
    float x = 0; 
    float y = 0;
    float z = 0;
    float depth = 0.89*(7.7+ log(e) );
    float theta = 2*atan(exp(-1*eta));
    float r = 0;
    if( fabs(eta) > 1.479 ) 
      { 
	float ecalZ = 315.4*fabs(eta)/eta;
	
	r = ecalZ / cos( 2*atan( exp( -1*eta ) ) ) + depth;
	x = r * cos( phi ) * sin( theta );
	y = r * sin( phi ) * sin( theta );
	z = r * cos( theta );
      }
    else
      {
	float rperp = 129.0;
	float zface =  sqrt( cos( theta ) * cos( theta ) /
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

}

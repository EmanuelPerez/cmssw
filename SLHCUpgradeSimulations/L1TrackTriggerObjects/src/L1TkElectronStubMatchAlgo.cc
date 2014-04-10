// -*- C++ -*-
//
// Package:    L1TrackTriggerObjects
// Class:      L1TkElectronStubMatchAlgo
// 
/**\class L1TkElectronStubMatchAlgo 

 Description: Algorithm to match L1EGamma oject with L1Track candidates

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  S. Dutta and A. Modak
//         Created:  Fri Feb 14 14:15:38 CET 2014
// $Id$
//
//


// system include files
#include <memory>
#include <cmath>


#include "DataFormats/Math/interface/deltaPhi.h"
#include "SLHCUpgradeSimulations/L1TrackTriggerObjects/interface/L1TkElectronStubMatchAlgo.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/Records/interface/StackedTrackerGeometryRecord.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
 
namespace L1TkElectronStubMatchAlgo {
  const double rcal = 129.0;
  double radii[7] = {0.0, 23.0, 35.7, 50.8, 68.6, 88.8, 108.0};
  // ------------ match EGamma and Track
  unsigned int doMatch(l1extra::L1EmParticleCollection::const_iterator egIter, edm::Handle< L1TkStub_PixelDigi_Collection > stubHandle,
		       const edm::ParameterSet& conf, const edm::EventSetup& setup, GlobalPoint bspot, std::vector<double>& zvals) {
    
    double ptMinCutoff = conf.getParameter<double>("StubMinPt");
    double dPhiCutoff  = conf.getParameter<double>("StubEGammaDeltaPhi");
    double dZCutoff    = conf.getParameter<double>("StubEGammaDeltaZ");
    double phiMissCutoff = conf.getParameter<double>("StubEGammaPhiMiss");
    double zMissCutoff   = conf.getParameter<double>("StubEGammaZMiss");
   
    GlobalPoint egPos = L1TkElectronStubMatchAlgo::calorimeterPosition(egIter->phi(), egIter->eta(), egIter->energy());
    
    // Magnetic Field
    double magnetStrength;
    edm::ESHandle<MagneticField> theMagField;
    setup.get<IdealMagneticFieldRecord>().get(theMagField);
    magnetStrength = theMagField.product()->inTesla(GlobalPoint(0,0,0)).z();

    edm::ESHandle<StackedTrackerGeometry> stackedGeometryHandle;
    setup.get<StackedTrackerGeometryRecord>().get(stackedGeometryHandle);
    
    std::vector<L1TkStubIter> selectedStubs;
    for (L1TkStubIter stubIter = stubHandle->begin(); stubIter != stubHandle->end(); ++stubIter) {
      float stub_pt = stackedGeometryHandle->findRoughPt(magnetStrength,&(*stubIter));
      unsigned int ilayer = getLayerId(stubIter);
      if ((ilayer%10) > 3) continue;
      
      double r   = stackedGeometryHandle->findGlobalPosition(&(*stubIter)).perp();
      double phi = stackedGeometryHandle->findGlobalPosition(&(*stubIter)).phi();
      double z = stackedGeometryHandle->findGlobalPosition(&(*stubIter)).z();
      
      
      double dPhi = getDPhi(egPos, egIter->et(), r, phi, magnetStrength);
      double zIntercept = getZIntercept(egPos, r, z);
      double scaledZInterceptCut;
      double scaledDPhiCut;
      double scaledPtMinCut;
      if (fabs(egIter->eta()) < 1.1) {
	//      if (ilayer < 10) {
	scaledZInterceptCut = getScaledZInterceptCut(ilayer,dZCutoff, 0.75, egPos.eta());
        scaledDPhiCut = dPhiCutoff;
        scaledPtMinCut = ptMinCutoff;
      }	else {
        scaledDPhiCut = 1.6* dPhiCutoff;
        scaledZInterceptCut = dZCutoff;
        scaledPtMinCut = ptMinCutoff*1.0;
      }    
      if (scaledPtMinCut > 0.0 && stub_pt <= scaledPtMinCut) continue;

      if ( (fabs(dPhi) < scaledDPhiCut) && 
	   (fabs(zIntercept)< scaledZInterceptCut) ) selectedStubs.push_back(stubIter);
    }
    sort( selectedStubs.begin(), selectedStubs.end(), compareStubLayer);
    int ncount = 0;  
    // Get two-point stubs
    for (std::vector<L1TkStubIter>::const_iterator istub1 = selectedStubs.begin(); istub1 != selectedStubs.end(); istub1++) {
      for (std::vector<L1TkStubIter>::const_iterator istub2 = istub1+1; istub2 != selectedStubs.end(); istub2++) {
        unsigned layer1 = getLayerId(*istub1);   
        unsigned layer2 = getLayerId(*istub2);   
	if (!selectLayers(egPos.eta(),layer1,layer2)) continue;
        bool barrel = true; 
        if (layer2 > 10) barrel = false;   

        double innerZ = stackedGeometryHandle->findGlobalPosition(&(*(*istub1))).z();
	double outerZ = stackedGeometryHandle->findGlobalPosition(&(*(*istub2))).z();
        double innerR = stackedGeometryHandle->findGlobalPosition(&(*(*istub1))).perp();
	double outerR = stackedGeometryHandle->findGlobalPosition(&(*(*istub2))).perp();
	double innerPhi = stackedGeometryHandle->findGlobalPosition(&(*(*istub1))).phi();
	double outerPhi = stackedGeometryHandle->findGlobalPosition(&(*(*istub2))).phi();
	
	GlobalPoint s1pos(stackedGeometryHandle->findGlobalPosition(&(*(*istub1))).x()-bspot.x(),
			  stackedGeometryHandle->findGlobalPosition(&(*(*istub1))).y()-bspot.y(),
			  stackedGeometryHandle->findGlobalPosition(&(*(*istub1))).z()-bspot.z());
	
	GlobalPoint s2pos(stackedGeometryHandle->findGlobalPosition(&(*(*istub2))).x()-bspot.x(),
			  stackedGeometryHandle->findGlobalPosition(&(*(*istub2))).y()-bspot.y(),
			  stackedGeometryHandle->findGlobalPosition(&(*(*istub2))).z()-bspot.z());
	

	if(!goodTwoPointPhi(innerR, outerR, innerPhi, outerPhi, magnetStrength)) continue;
	if(!goodTwoPointZ(innerR, outerR, innerZ, outerZ)) continue;
        
	double zMiss = getZMiss(egPos, innerR, outerR, innerZ, outerZ, barrel);
	double phiMiss = getPhiMiss(egIter->et(), s1pos, s2pos);
	
        double phiMissScaledCut = phiMissCutoff;
	if (fabs(egPos.eta()) >= 1.1) {
	  if (layer1 <= 3 && layer2 <= 3) phiMissScaledCut *= 1.4;
	  else phiMissScaledCut *= 1.8;
	}
        double zMissScaledCut; 
        if (barrel) {
	  zMissScaledCut = getScaledZMissCut(layer1, layer2,zMissCutoff, 0.04, egPos.eta());
        } else {
          zMissScaledCut = 2.0;
	}
	if(fabs(phiMiss)< phiMissScaledCut && fabs(zMiss) < zMissScaledCut) {
          ncount++;
          zvals.push_back(getCompatibleZPoint(innerR, outerR, innerZ, outerZ));
	}
      }
    }
    return ncount;
  }
  // -------------- get Calorimeter position
  GlobalPoint calorimeterPosition(double phi, double eta, double e) {
    double x = 0.; 
    double y = 0.;
    double z = 0.;
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
	double zface =  sqrt( cos( theta ) * cos( theta ) /
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
  // get Layer of a Stub
  unsigned int getLayerId(const L1TkStubIter& aStub) {
    unsigned int layer = 999999;
    if (StackedTrackerDetId(aStub->getDetId()).isBarrel()) {
      layer = StackedTrackerDetId(aStub->getDetId()).iLayer();
    } else if (StackedTrackerDetId(aStub->getDetId()).isEndcap()) {
      layer = StackedTrackerDetId(aStub->getDetId()).iSide() * 10 + StackedTrackerDetId(aStub->getDetId()).iDisk();
    } else {
      edm::LogError("L1TkElectronStubMatchAlgo")  << "Neither Barrel nor Endcap " << layer << std::endl;
    }
    if (layer > 100) edm::LogError("L1TkElectronStubMatchAlgo")  << " Wrong Layer " << layer << std::endl;
    return layer;
  }

  // Z-compatibility between two stubs
  bool  goodTwoPointZ(double innerR, double outerR, double innerZ, double outerZ ) {
    
    double mIPWidth = 200.0;
    double positiveZBoundary =  (mIPWidth - outerZ) * (outerR - innerR);
    double negativeZBoundary = -(mIPWidth + outerZ) * (outerR - innerR);
    double multipliedLocation = (innerZ - outerZ) * outerR;
    
    
    if( multipliedLocation < positiveZBoundary &&
	multipliedLocation > negativeZBoundary )
      return true;
    return false;
  }
  // Phi-compatibility between two stubs
  bool goodTwoPointPhi(double innerR, double outerR, double innerPhi, double outerPhi, double m_strength) {
  
    // Rebase the angles in terms of 0-2PI, should
    if ( innerPhi < 0.0 ) innerPhi += 2.0 * TMath::Pi();
    if ( outerPhi < 0.0 ) outerPhi += 2.0 * TMath::Pi();
    
    // Check for seed compatibility given a pt cut
    // Threshold computed from radial location of hits
    double mCompatibilityScalingFactor = 
      (100.0 * 2.0e+9 * 2.0) / (TMath::C() * m_strength);
    
    mCompatibilityScalingFactor = 1.0 / mCompatibilityScalingFactor;
    
    double deltaPhiThreshold = 
      (outerR - innerR) * mCompatibilityScalingFactor;  
    
    // Delta phi computed from hit phi locations
    double deltaPhi = outerPhi - innerPhi;
    if(deltaPhi<0) deltaPhi = -deltaPhi;
    
    if(deltaPhi<deltaPhiThreshold) return true;
    else return false;
  }
  // Delta Phi
  double getDPhi(GlobalPoint epos, double eet, double r, double phi, double m_strength) {
    
    double er = epos.perp();
    
    double phiVsRSlope = -3.00e-3 * m_strength / eet / 2.;
    
    // preselecton variable
    double psi = reco::deltaPhi(phi,epos.phi());
    double deltaPsi = psi - (er-r)*phiVsRSlope;
    double antiDeltaPsi = psi - (r-er)*phiVsRSlope;
    double dP;
    if (fabs(deltaPsi)<fabs(antiDeltaPsi)){
      dP = deltaPsi;
    }else{
      dP = antiDeltaPsi;
    }
    return dP;
  }
  // Z-Itercept
  double getZIntercept(GlobalPoint epos, double r, double z) {
	  
    double er = epos.perp();
    double ez = epos.z();
        
    double zint = (er*z - r*ez)/(er-r);
    return zint;
  }
  // PhiMiss
  double getPhiMiss(double eet, GlobalPoint spos1, GlobalPoint spos2) {

    double pT = eet;
    double curv = pT*100*.877;
    if (curv == 0) return 999.9;
    
    
    double r1 = spos1.perp();
    double r2 = spos2.perp();
    
    double phi1 = spos1.phi();
    double phi2 = spos2.phi();
    
    //Predict phi of hit 2
    double a = (r2-r1)/(2*curv);
    double b = reco::deltaPhi(phi2,phi1);
    
    double phiMiss = 0;
    if(fabs(b - a)<fabs(b + a)) phiMiss = b - a;
    if(fabs(b - a)>fabs(b + a)) phiMiss = b + a;
	    
    return phiMiss;

  }

  // ZMiss
  double getZMiss(GlobalPoint epos, double r1, double r2, double z1, double z2, bool bar) {

    double er = epos.perp();
    double ez = epos.z();
    
    
    double missVal ;
    if (bar) {
      missVal = z2 - (r2*(ez-z1)-r1*ez + 
		      er*z1)/(er-r1);
    } else {
      missVal = r2 - er - (er-r1)*(z2-ez)/(ez-z1);
    }
    return missVal;  
  }
  // Z-Intercept Cut
  double getScaledZInterceptCut(unsigned int layer, double cut, double cfac, double eta) {
    double mult  = (rcal-radii[layer])/(rcal-radii[1]);
    return (mult*(cut+cfac*(1.0/(1.0-cos(2*atan(exp(-1.0*fabs(eta)))))))); 
  }
  // Z-issCut
  double getScaledZMissCut(int layer1, int layer2, double cut, double cfac, double eta) {
    double mult = ( (radii[layer2] - radii[layer1])/(rcal-radii[layer1]) )*( (rcal -radii[1])/(radii[2]-radii[1]));
    return  (mult*(cut+cfac*(1.0/(1.0-cos(2*atan(exp(-1.0*fabs(eta))))))));
  }
  bool compareStubLayer( L1TkStubIter s1, L1TkStubIter s2) {
    unsigned int l1 = 0;
    unsigned int l2 = 0;
    if (StackedTrackerDetId(s1->getDetId()).isBarrel()) {
      l1 = StackedTrackerDetId(s1->getDetId()).iLayer();
    } else if (StackedTrackerDetId(s1->getDetId()).isEndcap()) {
      l1 = StackedTrackerDetId(s1->getDetId()).iSide() * 10 + StackedTrackerDetId(s1->getDetId()).iDisk();
    }
    if (StackedTrackerDetId(s2->getDetId()).isBarrel()) {
      l2 = StackedTrackerDetId(s2->getDetId()).iLayer();
    } else if (StackedTrackerDetId(s2->getDetId()).isEndcap()) {
      l2 = StackedTrackerDetId(s2->getDetId()).iSide() * 10 + StackedTrackerDetId(s2->getDetId()).iDisk();
    }
    return l1 < l2;
  }

  bool selectLayers(float eta, int l1, int l2) {
    bool select = false;
    if ((fabs(eta) < 1.3) && !(l1 < 4 && l2 < 4)) return select;
    if ((fabs(eta) >= 1.3 && TMath::Abs(eta) < 1.7) && !(
							 (l1 == 1 && l2 == 2)  ||
							 (l1 == 1 && l2 == 11) ||
							 (l1 == 1 && l2 == 21) ||
							 (l1 == 2 && l2 == 11) ||
							 (l1 == 2 && l2 == 21) ||
							 (l1 == 11 && l2 == 12)||
							 (l1 == 21 && l2 == 22)))  return select;
    if ((fabs(eta) >= 1.7 && TMath::Abs(eta) <= 2.3) && !(
							  (l1 == 1 && l2 == 11) ||
							  (l1 == 1 && l2 == 12) ||
							  (l1 == 1 && l2 == 21) ||
							  (l1 == 1 && l2 == 22) ||
							  (l1 == 11 && l2 == 12) ||
							  (l1 == 21 && l2 == 22) ||
							  (l1 == 12 && l2 == 13) ||
							  (l1 == 22 && l2 == 23))) return select;
    select = true;
    return select;
  }
  double getCompatibleZPoint(double r1, double r2, double z1, double z2) {
    return (z1 - r1*(z2-z1)/(r2-r1));
  }
}

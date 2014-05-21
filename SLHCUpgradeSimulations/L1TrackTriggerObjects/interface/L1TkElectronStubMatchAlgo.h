#ifndef L1TkElectronStubMatchAlgo_HH
#define L1TkElectronStubMatchAlgo_HH

#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1EmParticleFwd.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "SimDataFormats/SLHC/interface/StackedTrackerTypes.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "Geometry/TrackerGeometryBuilder/interface/StackedTrackerGeometry.h"


namespace L1TkElectronStubMatchAlgo {
  typedef L1TkStub_PixelDigi_Collection::const_iterator L1TkStubIter;
  unsigned int doMatch(l1extra::L1EmParticleCollection::const_iterator egIter, edm::Handle< L1TkStub_PixelDigi_Collection > stubHandle,
		       const edm::ParameterSet& conf, const edm::EventSetup& iSetup, GlobalPoint bspot, std::vector<double>& zvals);
  GlobalPoint calorimeterPosition(double phi, double eta, double e);
  unsigned int getLayerId(const L1TkStubIter& aStub);
  bool goodTwoPointZ(double innerR, double outerR, double innerZ, double outerZ );
  bool goodTwoPointPhi(double innerR, double outerR, double innerPhi, double outerPhi, double m_strength);
  double getDPhi(GlobalPoint epos, double eet, double r, double phi, double m_strength);
  double getZIntercept(GlobalPoint epos, double r, double z);
  double getPhiMiss(double eet, GlobalPoint spos1, GlobalPoint spos2);
  double getZMiss(GlobalPoint epos, double r1, double r2, double z1, double z2, bool bar);
  double getScaledZInterceptCut(unsigned int layer, double cut, double cfac, double eta);
  double getScaledZMissCut(int layer1, int layer2, double cut, double cfac, double eta);
  bool compareStubLayer( L1TkStubIter s1, L1TkStubIter s2);
  bool selectLayers(float eta, int l1, int l2); 
  double getCompatibleZPoint(double r1, double r2, double z1, double z2);
}
#endif

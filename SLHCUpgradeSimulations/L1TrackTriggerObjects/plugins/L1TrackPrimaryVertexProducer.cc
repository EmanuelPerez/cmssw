// -*- C++ -*-
//
//
// Original Author:  Emmanuelle Perez,40 1-A28,+41227671915,
//         Created:  Tue Nov 12 17:03:19 CET 2013
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"


#include "SimDataFormats/SLHC/interface/StackedTrackerTypes.h"
#include "SimDataFormats/SLHC/interface/slhcevent.hh"

////////////////////////////
// DETECTOR GEOMETRY HEADERS
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetType.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelTopologyBuilder.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/RectangularPixelTopology.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"

#include "Geometry/Records/interface/StackedTrackerGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/StackedTrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/StackedTrackerDetUnit.h"
#include "DataFormats/SiPixelDetId/interface/StackedTrackerDetId.h"


#include "DataFormats/L1Trigger/interface/L1TrackPrimaryVertex.h"


//
// class declaration
//

class L1TrackPrimaryVertexProducer : public edm::EDProducer {
   public:

  typedef L1TkTrack_PixelDigi_                          L1TkTrackType;
  typedef std::vector< L1TkTrackType >                               L1TkTrackCollectionType;

      explicit L1TrackPrimaryVertexProducer(const edm::ParameterSet&);
      ~L1TrackPrimaryVertexProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


      float MaxPtVertex(const edm::Handle<L1TkTrackCollectionType> & L1TkTrackHandle,
                float& sum,
                int nStubsmin, int nPSmin, float ptmin, int imode) ;

      float SumPtVertex(const edm::Handle<L1TkTrackCollectionType> & L1TkTrackHandle,
                float z, int nStubsmin, int nPSmin, float ptmin, int imode) ;


   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      //virtual void endRun(edm::Run&, edm::EventSetup const&);
      //virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      //virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------

	float ZMAX;	// in cm
	float DeltaZ;	// in cm
	float CHI2MAX;

        const StackedTrackerGeometry*                   theStackedGeometry;

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
L1TrackPrimaryVertexProducer::L1TrackPrimaryVertexProducer(const edm::ParameterSet& iConfig)
{
   //register your products
/* Examples
   produces<ExampleData2>();

   //if do put with a label
   produces<ExampleData2>("label");
 
   //if you want to put into the Run
   produces<ExampleData2,InRun>();
*/
   //now do what ever other initialization is needed
  
  ZMAX = (float)iConfig.getParameter<double>("ZMAX");
  DeltaZ = (float)iConfig.getParameter<double>("DeltaZ");
  CHI2MAX = (float)iConfig.getParameter<double>("CHI2MAX");


  produces<L1TrackPrimaryVertexCollection>();

}


L1TrackPrimaryVertexProducer::~L1TrackPrimaryVertexProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
L1TrackPrimaryVertexProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
/* This is an event example
   //Read 'ExampleData' from the Event
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);

   //Use the ExampleData to create an ExampleData2 which 
   // is put into the Event
   std::auto_ptr<ExampleData2> pOut(new ExampleData2(*pIn));
   iEvent.put(pOut);
*/

/* this is an EventSetup example
   //Read SetupData from the SetupRecord in the EventSetup
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
*/
 
 std::auto_ptr<L1TrackPrimaryVertexCollection> result(new L1TrackPrimaryVertexCollection);


  edm::Handle<L1TkTrackCollectionType> L1TkTrackHandle;
  iEvent.getByLabel("L1Tracks","Level1TkTracks",L1TkTrackHandle);

	// max(Sum PT2), tracks with at least 3 stubs in PS layers
   float sum1 = -999;
   int nStubsmin = 4;
   int nPSmin = 3;
   float ptmin = 2. ;
   int imode = 2;
   float z1 = MaxPtVertex( L1TkTrackHandle, sum1, nStubsmin, nPSmin, ptmin, imode );
   L1TrackPrimaryVertex vtx1( z1, sum1 );

	// max(Sum PT2), no constraint on the number of stubs in PS
   float sum2 = -999;
   nPSmin = 0;
   float z2 = MaxPtVertex( L1TkTrackHandle, sum2, nStubsmin, nPSmin, ptmin, imode );
   L1TrackPrimaryVertex vtx2( z2, sum2 );

	// max(Sum PT), tracks with at least 3 stubs in PS layers
   float sum3 = -999;
   nPSmin = 3;
   imode = 1;
   float z3 = MaxPtVertex( L1TkTrackHandle, sum3, nStubsmin, nPSmin, ptmin, imode );
   L1TrackPrimaryVertex vtx3( z3, sum3 );

 result -> push_back( vtx1 );
 result -> push_back( vtx2 );
 result -> push_back( vtx3 );

 iEvent.put( result);
}


float L1TrackPrimaryVertexProducer::MaxPtVertex(const edm::Handle<L1TkTrackCollectionType> & L1TkTrackHandle,
 		float& sum,
		int nStubsmin, int nPSmin, float ptmin, int imode) {
        // return the zvtx corresponding to the max(SumPT)
        // of tracks with at least nPSmin stubs in PS modules
   
      float sumMax = 0;
      float zvtxmax = -999;
      int nIter = (int)(ZMAX * 10. * 2.) ;
      for (int itest = 0; itest <= nIter; itest ++) {
   
        //float z = -100 + itest;         // z in mm
	float z = -ZMAX * 10 + itest ;  	// z in mm
        z = z/10.  ;   // z in cm
        float sum = SumPtVertex(L1TkTrackHandle, z, nStubsmin, nPSmin, ptmin, imode);
        if (sumMax >0 && sum == sumMax) {
          //cout << " Note: Several vertices have the same sum " << zvtxmax << " " << z << " " << sumMax << endl;
        }
   
        if (sum > sumMax) {
           sumMax = sum;
           zvtxmax = z;
        }
      }
   
 sum = sumMax;
 return zvtxmax;
}  


float L1TrackPrimaryVertexProducer::SumPtVertex(const edm::Handle<L1TkTrackCollectionType> & L1TkTrackHandle,
		float z, int nStubsmin, int nPSmin, float ptmin, int imode) {

        // sumPT of tracks with >= nPSmin stubs in PS modules
        // z in cm
 float sumpt = 0;


  L1TkTrackCollectionType::const_iterator trackIter;

  for (trackIter = L1TkTrackHandle->begin(); trackIter != L1TkTrackHandle->end(); ++trackIter) {

    float pt = trackIter->getMomentum().perp();
    float chi2 = trackIter->getChi2();
    float ztr  = trackIter->getVertex().z();

    if (pt < ptmin) continue;
    if (fabs(ztr) > ZMAX ) continue;
    if (chi2 > CHI2MAX) continue;
    if ( fabs(ztr - z) > DeltaZ) continue;   // eg DeltaZ = 1 mm


	// get the number of stubs and the number of stubs in PS layers
    float nPS = 0.;     // number of stubs in PS modules
    float nstubs = 0;

      // get pointers to stubs associated to the L1 track
      std::vector< edm::Ptr< L1TkStub_PixelDigi_ > > theStubs = trackIter ->getStubPtrs();
      int tmp_trk_nstub = (int) theStubs.size();
      if ( tmp_trk_nstub < 0) continue;

      // loop over the stubs
      for (unsigned int istub=0; istub<(unsigned int)theStubs.size(); istub++) {
        bool genuine = theStubs.at(istub)->isGenuine();
        if (genuine) {
           nstubs ++;
           StackedTrackerDetId detIdStub( theStubs.at(istub)->getDetId() );
           bool isPS = theStackedGeometry -> isPSModule( detIdStub );
           //if (isPS) cout << " this is a stub in a PS module " << endl;
           if (isPS) nPS ++;
	} // endif genuine
       } // end loop over stubs

        if (imode == 1 || imode == 2 ) {
            if (nPS < nPSmin) continue;
        }
	if ( nstubs < nStubsmin) continue;

        if (imode == 2) sumpt += pt*pt;
        if (imode == 1) sumpt += pt;

  } // end loop over the tracks

 return sumpt;

}



// ------------ method called once each job just before starting event loop  ------------
void 
L1TrackPrimaryVertexProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
L1TrackPrimaryVertexProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void
L1TrackPrimaryVertexProducer::beginRun(edm::Run& iRun, edm::EventSetup const& iSetup)
{

  /// Geometry handles etc
  edm::ESHandle<TrackerGeometry>                               geometryHandle;
  //const TrackerGeometry*                                       theGeometry;
  edm::ESHandle<StackedTrackerGeometry>           stackedGeometryHandle;

  /// Geometry setup
  /// Set pointers to Geometry
  iSetup.get<TrackerDigiGeometryRecord>().get(geometryHandle);
  //theGeometry = &(*geometryHandle);
  /// Set pointers to Stacked Modules
  iSetup.get<StackedTrackerGeometryRecord>().get(stackedGeometryHandle);
  theStackedGeometry = stackedGeometryHandle.product(); /// Note this is different 
                                                        /// from the "global" geometry
  if (theStackedGeometry == 0) cout << " theStackedGeometry = 0 " << endl;      // for compil when not used...

}
 
// ------------ method called when ending the processing of a run  ------------
/*
void
L1TrackPrimaryVertexProducer::endRun(edm::Run&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
L1TrackPrimaryVertexProducer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
L1TrackPrimaryVertexProducer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
L1TrackPrimaryVertexProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1TrackPrimaryVertexProducer);

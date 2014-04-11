#ifndef L1TTRACK_DEGRADE_H
#define L1TTRACK_DEGRADE_H

#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
//
// #include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
//
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/EDProduct.h"
#include "SimDataFormats/SLHC/interface/StackedTrackerTypes.h"
//#include "SimDataFormats/SLHC/interface/slhcevent.hh"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

//#include "SimDataFormats/SLHC/interface/L1TBarrel.hh"
//#include "SimDataFormats/SLHC/interface/L1TDisk.hh"
//#include "SimDataFormats/SLHC/interface/L1TStub.hh"

#include "TRandom.h"
#include "TMath.h"


using namespace edm;

class L1TrackDegrader : public edm::EDProducer
{ 
public:

  typedef L1TkTrack_PixelDigi_                          L1TkTrackType;
  typedef std::vector< L1TkTrackType >                               L1TkTrackCollectionType;

  /// Constructor/destructor
  explicit L1TrackDegrader(const edm::ParameterSet& iConfig);
  virtual ~L1TrackDegrader();

protected:

private:
         edm::InputTag L1TrackInputTag;

	 TRandom ran;


  /// ///////////////// ///
  /// MANDATORY METHODS ///
  virtual void beginRun( edm::Run& run, const edm::EventSetup& iSetup );
  virtual void endRun( edm::Run& run, const edm::EventSetup& iSetup );
  virtual void produce( edm::Event& iEvent, const edm::EventSetup& iSetup );
};


L1TrackDegrader::L1TrackDegrader(edm::ParameterSet const& iConfig) // :   config(iConfig)
{

   L1TrackInputTag = iConfig.getParameter<edm::InputTag>("L1TrackInputTag");
   ran.SetSeed(0);

  produces< L1TkTrackCollectionType >( "Level1TkTracks" ).setBranchAlias("Level1TkTracks");

}

L1TrackDegrader::~L1TrackDegrader()
{ 
  /// Insert here what you need to delete
  /// when you close the class instance
}


void L1TrackDegrader::beginRun(edm::Run& run, const edm::EventSetup& iSetup )
{
}


void L1TrackDegrader::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

 edm::Handle<L1TkTrackCollectionType> L1TkTrackHandle;
 iEvent.getByLabel(L1TrackInputTag, L1TkTrackHandle); 
 L1TkTrackCollectionType::const_iterator trackIter;

  std::auto_ptr< L1TkTrackCollectionType > L1TkTracksForOutput( new L1TkTrackCollectionType );

        for (trackIter = L1TkTrackHandle->begin(); trackIter != L1TkTrackHandle->end(); ++trackIter) {

	   float Eta = fabs( trackIter->getMomentum().eta() );
	   //float Pt = trackIter->getMomentum().perp();
	   float z0 = trackIter->getVertex().z();

	   // hard-coding of the z-smearing based on Stefano's plots

	   float sigma= 0;

	   if (Eta >=  0.5 && Eta < 2) sigma = 0.168 * Eta - 0.08 ;
	   if (Eta >= 2 && Eta < 2.3) sigma = 0.256;

	   float deltaZ = ran.Gaus(0.,sigma);
	   float smeared_z = z0 + deltaZ ;

	   float x0 = trackIter->getVertex().x();
	   float y0 = trackIter->getVertex().y();

	   GlobalPoint smearedVertex(x0, y0, smeared_z);
	   
	   L1TkTrackType smearedTrack = *trackIter;
	   smearedTrack.setVertex( smearedVertex );

	   L1TkTracksForOutput->push_back( smearedTrack );

	}

  iEvent.put( L1TkTracksForOutput, "Level1TkTracks");

}


void L1TrackDegrader::endRun(edm::Run& run, const edm::EventSetup& iSetup)
{
  /// Things to be done at the exit of the event Loop 

}

// ///////////////////////////
// // DEFINE THIS AS A PLUG-IN
DEFINE_FWK_MODULE(L1TrackDegrader);

#endif








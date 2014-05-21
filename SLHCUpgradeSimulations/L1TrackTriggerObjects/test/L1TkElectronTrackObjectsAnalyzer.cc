// -*- C++ -*-
//
// Package:    L1TkElectronTrackObjectsAnalyzer
// Class:      L1TkElectronTrackObjectsAnalyzer
// 
/**\class L1TkElectronTrackObjectsAnalyzer L1TkElectronTrackObjectsAnalyzer.cc SLHCUpgradeSimulations/L1TkElectronTrackObjectsAnalyzer/src/L1TkElectronTrackObjectsAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
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
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "DataFormats/L1TrackTrigger/interface/L1TkEmParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkEmParticleFwd.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkElectronParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkElectronParticleFwd.h"


#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"


#include "DataFormats/Math/interface/deltaPhi.h"
#include "SLHCUpgradeSimulations/L1TrackTriggerObjects/interface/L1TkElectronEtComparator.h"
#include "TH1F.h"


using namespace l1extra;



//
// class declaration
//

class L1TkElectronTrackObjectsAnalyzer : public edm::EDAnalyzer {
   public:

  typedef L1TkTrack_PixelDigi_                          L1TkTrackType;
  typedef std::vector< L1TkTrackType >                               L1TkTrackCollectionType;

  explicit L1TkElectronTrackObjectsAnalyzer(const edm::ParameterSet&);
  ~L1TkElectronTrackObjectsAnalyzer();
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  void fillIntegralHistos(TH1F* th, float var);  
  void checkEfficiency();
  void checkRate();
  int matchEGWithGenParticle();
  int matchEGWithSimTrack();
  // for L1TrackEmParticles
  edm::InputTag L1EGammaInputTag;
  edm::InputTag L1TkElectronsInputTag;

  int selectedEGTot_;
  int selectedEGTrkTot_;

  TH1F* etaGen_;
  TH1F* etaEGamma_;
  TH1F* etaEGammaTrk_;

  TH1F* etGen_;
  TH1F* etGenEGamma_;
  TH1F* etGenEGammaTurnOn_;
  TH1F* etGenEGammaTrk_;
  TH1F* etGenEGammaTrkTurnOn_;
  TH1F* etEGamma_;
  TH1F* etEGammaTrk_;
  TH1F* isoEGammaTrk_;


  std::string analysisOption_;
  float etaCutoff_;
  float trkPtCutoff_;

  edm::SimTrackContainer simTracks_;
  l1extra::L1EmParticleCollection eGammaCollection;
  edm::Handle<L1TkElectronParticleCollection> L1TrackElectronsHandle;
  edm::Handle<L1EmParticleCollection> EGammaHandle;
  edm::Handle<reco::GenParticleCollection> genParticleHandle;
  int ievent; 
};

L1TkElectronTrackObjectsAnalyzer::L1TkElectronTrackObjectsAnalyzer(const edm::ParameterSet& iConfig) {

  edm::Service<TFileService> fs;
  L1EGammaInputTag = iConfig.getParameter<edm::InputTag>("L1EGammaInputTag") ;
  L1TkElectronsInputTag = iConfig.getParameter<edm::InputTag>("L1TkElectronsInputTag");

  analysisOption_ = iConfig.getParameter<std::string>("AnalysisOption");
  etaCutoff_ = iConfig.getParameter<double>("EtaCutOff");
  trkPtCutoff_ = iConfig.getParameter<double>("TrackPtCutOff");

  
}
void L1TkElectronTrackObjectsAnalyzer::beginJob() {
  edm::Service<TFileService> fs;

  etaGen_ = fs->make<TH1F>("Eta_Gen","Eta of GenParticle", 31, -3.14, 3.14);
  etaEGamma_ = fs->make<TH1F>("Eta_EGamma","Eta of EGamma", 31, -3.14, 3.14);
  etaEGammaTrk_ = fs->make<TH1F>("Eta_EGammaTrk","Eta of TrkEGamma", 31, -3.14, 3.14);
  isoEGammaTrk_ = fs->make<TH1F>("Isolation_EGammaTrk","Isolation of TrkEGamma", 200, 0.0, 2.0);
    
  if (analysisOption_ == "Efficiency") {
    etGen_ = fs->make<TH1F>("GenEt","Et of GenParticle", 20, 9.5, 59.5);
    etGenEGamma_    = fs->make<TH1F>("GenEt_EGamma","Et of GenParticle (EG > 0)", 20, 9.5, 59.5);
    etGenEGammaTurnOn_ = fs->make<TH1F>("GenEt_EGammaEt","Et of GenParticle (EG > 20)", 20, 9.5, 59.5);
    etGenEGammaTrk_ = fs->make<TH1F>("GenEt_EGammaTrk","Et of GenParticle (EGTrk > 0)", 20, 9.5, 59.5);
    etGenEGammaTrkTurnOn_ = fs->make<TH1F>("GenEt_EGammaTrkEt","Et of GenParticle (EGTrk 20)", 20, 9.5, 59.5);
  } else {
    etEGamma_ = fs->make<TH1F>("EGammaEtThresholdEvt_Ref","Et of EGamma (EventEt threshold)", 90, 4.5, 94.5);
    etEGammaTrk_ = fs->make<TH1F>("EGammaEtThresholdEvt_Track","Et of TrkEGamma( Event Et threshold)", 90, 4.5, 94.5);
  }

  selectedEGTot_ = 0;
  selectedEGTrkTot_ = 0;
  ievent = 0;
}

L1TkElectronTrackObjectsAnalyzer::~L1TkElectronTrackObjectsAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
L1TkElectronTrackObjectsAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  ievent++;  

  iEvent.getByLabel(L1TkElectronsInputTag, L1TrackElectronsHandle);

  iEvent.getByLabel(L1EGammaInputTag,EGammaHandle);
  eGammaCollection = (*EGammaHandle.product());
  
  Handle<edm::SimTrackContainer> simTrackHandle;
  iEvent.getByLabel( "g4SimHits", "", simTrackHandle);
  simTracks_ = (*simTrackHandle.product());
   
  iEvent.getByLabel("genParticles", genParticleHandle);
  if (analysisOption_ == "Efficiency") checkEfficiency();
  else checkRate();
}
void L1TkElectronTrackObjectsAnalyzer::endJob() {
  std::cout << " Selected EGammas " << selectedEGTot_ << std::endl;
  std::cout << " Selected Track EGammas " << selectedEGTrkTot_ << std::endl;
  if (analysisOption_ == "Rate") {
    float scale_fac = 30000.0/ievent;
    etEGamma_->Scale(scale_fac);
    etEGammaTrk_->Scale(scale_fac);
  }     
}
void L1TkElectronTrackObjectsAnalyzer::checkEfficiency() {
  float genPt = (*genParticleHandle)[0].pt();
  float genEta = (*genParticleHandle)[0].eta();
  //  float genPhi = (*genParticleHandle)[0].phi();

  int nSelectedEG = matchEGWithGenParticle();
  if (nSelectedEG == 0 ) return;
  
  std::vector<L1TkElectronParticle>::const_iterator egTrkIter ;
  int nSelectedEGTrk = 0;
  int nSelectedEGTrkEt = 0;
  //  float dRmin = 999.9;
  //  float eta_min;
  for (egTrkIter = L1TrackElectronsHandle -> begin(); egTrkIter != L1TrackElectronsHandle->end(); ++egTrkIter) {

    if (fabs(egTrkIter->eta()) < etaCutoff_ && egTrkIter->pt() > 0) {
      if ( egTrkIter->getTrkPtr().isNonnull() && egTrkIter->getTrkPtr()->getMomentum().perp() <= trkPtCutoff_) continue;
      nSelectedEGTrk++;
      if (egTrkIter->pt() > 20) {
	nSelectedEGTrkEt++;
	//	float dPhi = reco::deltaPhi(egTrkIter->phi(), genPhi);
	//	float dEta = (egTrkIter->eta() - genEta);
	//	float dR =  sqrt(dPhi*dPhi + dEta*dEta);
	//	if (dR < dRmin) {
	//	  dRmin = dR;
	//	  eta_min = egTrkIter->eta(); 
      }
    }
  }
  if (nSelectedEGTrk > 0) {
    std::cout<< "Event # " << ievent << " Selected  EGamma "<< nSelectedEG << " matched EGamma Trk " << nSelectedEGTrk << std::endl;     
    etGenEGammaTrk_->Fill(genPt);
    if (nSelectedEGTrkEt > 0) etGenEGammaTrkTurnOn_->Fill(genPt);
    if (nSelectedEGTrkEt > 0 && genPt > 30) etaEGammaTrk_->Fill(genEta);
  }

  selectedEGTot_ += nSelectedEG;
  selectedEGTrkTot_ += nSelectedEGTrk;
}
void L1TkElectronTrackObjectsAnalyzer::checkRate() {

  int nSelectedEGTrk = 0;
  int nSelectedEG = 0;
  sort(eGammaCollection.begin(), eGammaCollection.end(), L1TkElectron::EtComparator());
  l1extra::L1EmParticleCollection::const_iterator egIter; 
  for (egIter = eGammaCollection.begin();  egIter != eGammaCollection.end(); ++egIter) {
    
    //    int ibx = egIter -> bx();
    //    if (ibx != 0) continue;
    
    float e_ele   = egIter->energy();
    float eta_ele = egIter->eta(); 
    float phi_ele = egIter->phi(); 
    float et_ele = 0;
    if (cosh(eta_ele) > 0.0) et_ele = e_ele/cosh(eta_ele);
    else et_ele = -1.0;
    if (fabs(eta_ele) >= etaCutoff_) continue;
    nSelectedEG++;
    if (nSelectedEG == 1) {
      fillIntegralHistos(etEGamma_, et_ele);
    }

    float et_min; 
    float dRmin = 999.9;
    float iso_min = 999.0; 
    std::vector<L1TkElectronParticle>::const_iterator egTrkIter ;
    for (egTrkIter = L1TrackElectronsHandle -> begin(); egTrkIter != L1TrackElectronsHandle->end(); ++egTrkIter) {
    if (fabs(egTrkIter->eta()) < etaCutoff_ && egTrkIter->pt() > 0) {
      if ( egTrkIter->getTrkPtr().isNonnull() && egTrkIter->getTrkPtr()->getMomentum().perp() <= trkPtCutoff_) continue;
	//      if (fabs(egTrkIter->getEGRef()->eta()) >= etaCutoff_) continue;
	float dPhi = reco::deltaPhi(phi_ele, egTrkIter->getEGRef()->phi());
	float dEta = (eta_ele - egTrkIter->getEGRef()->eta());
	float dR =  sqrt(dPhi*dPhi + dEta*dEta);
	if (dR < dRmin) {
	  dRmin = dR;
          iso_min = egTrkIter->getTrkIsol();
	  et_min = egTrkIter->pt();
	}
      }
    }
    if (dRmin < 0.1) {
      nSelectedEGTrk++;
      if (nSelectedEGTrk == 1) {
        if (et_min> 20.0) isoEGammaTrk_->Fill(iso_min);
	std::cout << "Selected EGammaTrk objet in the event " << et_ele << " with pt "<< egTrkIter->pt() << std::endl;           
	fillIntegralHistos(etEGammaTrk_, et_min);
      }
    }         
  }
  selectedEGTot_ += nSelectedEG;
  selectedEGTrkTot_ += nSelectedEGTrk;
}
void L1TkElectronTrackObjectsAnalyzer::fillIntegralHistos(TH1F* th, float var){
  int nbin = th->FindBin(var); 
  for (int ibin = 1; ibin < nbin+1; ibin++) th->Fill(th->GetBinCenter(ibin));
}
int L1TkElectronTrackObjectsAnalyzer::matchEGWithSimTrack() {
  if ( fabs(simTracks_[0].momentum().eta())> etaCutoff_ || simTracks_[0].momentum().pt() <= 0.0) return -1;
  etGen_->Fill(simTracks_[0].momentum().pt());
  etaGen_->Fill(simTracks_[0].momentum().eta()); 
  int nEG = 0;
  int nEGEt = 0;
  for (unsigned int igam = 0; igam != eGammaCollection.size(); igam++) {
    int ibx = eGammaCollection[igam].bx();
    if (ibx != 0) continue;
    
    float e_ele   = eGammaCollection[igam].energy();
    float eta_ele = eGammaCollection[igam].eta(); 
    //    float phi_ele = eGammaCollection[igam].phi(); 
    float et_ele = 0;
    if (cosh(eta_ele) > 0.0) et_ele = e_ele/cosh(eta_ele);
    else et_ele = -1.0;
    if ( fabs(eta_ele) > etaCutoff_ || et_ele <= 0.0) continue;
    nEG++;
    if (et_ele > 20.0)   nEGEt++;
  } 

  if (nEG > 0) etGenEGamma_->Fill(simTracks_[0].momentum().pt());
  if (nEGEt > 0) etGenEGammaTurnOn_->Fill(simTracks_[0].momentum().pt());
  return nEG;
}
int L1TkElectronTrackObjectsAnalyzer::matchEGWithGenParticle() {
  const reco::Candidate & p = (*genParticleHandle)[0];
  if ( fabs(p.eta()) > etaCutoff_ || p.pt() <= 0.0) return -1;
  etGen_->Fill(p.pt()); 
  if (p.pt() > 30.0) etaGen_->Fill(p.eta()); 
  int nEG = 0;
  int nEGEt = 0;
  //  float dRmin = 999.9;
  //  float eta_min = -5.0;
  for (unsigned int igam = 0; igam != eGammaCollection.size(); igam++) {
    int ibx = eGammaCollection[igam].bx();
    if (ibx != 0) continue;
    
    float eta_ele = eGammaCollection[igam].eta(); 
    float phi_ele = eGammaCollection[igam].phi(); 
    float e_ele   = eGammaCollection[igam].energy();
    float et_ele = 0;
    if (cosh(eta_ele) > 0.0) et_ele = e_ele/cosh(eta_ele);
    else et_ele = -1.0;
    if ( fabs(eta_ele) > etaCutoff_ || et_ele <= 0.0) continue;
    nEG++;
    if (et_ele > 20.0) {
      float dPhi = reco::deltaPhi(p.phi(), phi_ele);
      float dEta = (p.eta() - eta_ele);
      float dR =  sqrt(dPhi*dPhi + dEta*dEta);
      //      if (dR < dRmin) {
      //        dRmin = dR;
      //        eta_min = eta_ele; 
      if (dR < 0.5) nEGEt++;
    }
  }
  if (nEG > 0) etGenEGamma_->Fill(p.pt());
  if (nEGEt > 0) etGenEGammaTurnOn_->Fill(p.pt());
  if (nEGEt > 0 && p.pt() > 30.0) etaEGamma_->Fill(p.eta());
  return nEG;
}
//define this as a plug-in
DEFINE_FWK_MODULE(L1TkElectronTrackObjectsAnalyzer);

#define GCC_VERSION ( 10000 * __GNUC__ + 100 * __GNUC_MINOR__ + __GNUC_PATCHLEVEL__ )
// for ultra-legacy AOD (10_6_4 or higher)
#if GCC_VERSION > 70400
// this one comes on top of the previous
#define CMSSW106plus
#endif
#if GCC_VERSION > 80300
// for Run 3 MC studies
// this one comes on top of the previous
// GCC_VERSION preliminary, might need to be changed/sharpened
#define CMSSW11plus
#endif
#if GCC_VERSION > 90299
// for 2021 pilot data 
// this one comes on top of the previous
#define CMSSW12plus
#endif

// system include files
#include <memory>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include "Math/VectorUtil.h"

// Root
#include "TMath.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TTree.h"

// Utilities
#include "../interface/helperfunc.h"
#include "../interface/MyStruct.h"


//*****************************
// general user include files *
//*****************************
#include "FWCore/Framework/interface/Frameworkfwd.h"
#ifndef CMSSW12plus
// Run 1 and 2
#include "FWCore/Framework/interface/EDAnalyzer.h"
#else
// Run 3
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#endif

// FWCore
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

// ------ EXTRA HEADER FILES--------------------//
#include "FWCore/Framework/interface/EventSetup.h"
#ifndef CMSSW12plus
#include "FWCore/Framework/interface/ESHandle.h"
#endif

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Math/interface/deltaR.h"

//**************************
// for trigger information *
//**************************
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/L1Trigger/interface/BXVector.h"
#include "DataFormats/L1Trigger/interface/EGamma.h"

//***************************
// for tracking information *
//***************************
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

//*************************
// for vertex information *
//*************************
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include <DataFormats/VertexReco/interface/Vertex.h>
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

//***********************
// for electron information *
//***********************
#include "DataFormats/PatCandidates/interface/Electron.h"

//*******************************
// for gen particle information *
//*******************************
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

// set namespaces
using namespace edm;
using namespace reco;
using namespace std;


//***********************************
// main analyzer class (EDAnalyzer) *
//***********************************

#ifndef CMSSW12plus
// Run 1 and 2
class NanoAnalyzerDoubleMuDebug : public edm::EDAnalyzer
#else
// Run 3
class NanoAnalyzerDoubleMuDebug : public edm::one::EDAnalyzer<>
#endif
{
public:
  explicit NanoAnalyzerDoubleMuDebug(const edm::ParameterSet&);
  ~NanoAnalyzerDoubleMuDebug();

  static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);

private:

  virtual void beginJob(const edm::ParameterSet& iConfig);
  virtual void beginRun(const edm::Run &iRun, const edm::EventSetup &iStp);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void endJob();

  void createBranch();

  void reset();

  const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> bFieldToken_;
  EDGetTokenT<reco::VertexCollection> verticeToken_;
  EDGetTokenT<pat::ElectronCollection> electronToken_;
  EDGetTokenT<edm::TriggerResults> triggerToken_;  
  EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerobjectToken_;
  edm::EDGetTokenT<l1t::EGammaBxCollection> l1EG_;      

  Handle<pat::ElectronCollection> electrons_;
  Handle< reco::VertexCollection > vertices_;
  Handle< edm::TriggerResults> HLTtriggers_;
  Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;


/////////////////////////////////////////////////////////////////////////
////////////////////////// declare tree, file, //////////////////////////
/////////////////////////////////////////////////////////////////////////
  
  edm::Service<TFileService> fs;
  TTree* tree_;

  /// original nanos
  UInt_t run;
  ULong64_t event;
  UInt_t luminosityBlock;
  int nvtx;

  float  Jpsi_e1_pt;
  float  Jpsi_e1_eta;
  float  Jpsi_e1_phi;
  float  Jpsi_e1_mass;
  int    Jpsi_e1_q;   
  float  Jpsi_e1_isElectron;   
  float  Jpsi_e1_passMVA;   
  // 
  float  Jpsi_e1_bestL1pt;
  float  Jpsi_e1_bestL1eta;
  float  Jpsi_e1_bestL1phi;
  float  Jpsi_e1_bestL1dR;
  //
  //
  float  Jpsi_e2_pt      ;
  float  Jpsi_e2_eta     ;
  float  Jpsi_e2_phi     ;
  float  Jpsi_e2_mass     ;
  int    Jpsi_e2_q   ;   
  float  Jpsi_e2_isElectron;     
  float  Jpsi_e2_passMVA  ;
  // 
  float  Jpsi_e2_bestL1pt ;
  float  Jpsi_e2_bestL1eta ;
  float  Jpsi_e2_bestL1phi ;
  float  Jpsi_e2_bestL1dR ;
  //
  //
  float  Jpsi_e1_trgobj_pt      ;
  float  Jpsi_e1_trgobj_eta     ;
  float  Jpsi_e1_trgobj_phi     ;
  float  Jpsi_e1_trgobj_mass     ;
  int    Jpsi_e1_trgobj_q      ;   
  float  Jpsi_e1_trgobj_dR      ;
  // 
  float  Jpsi_e2_trgobj_pt ;
  float  Jpsi_e2_trgobj_eta ;
  float  Jpsi_e2_trgobj_phi ;
  float  Jpsi_e2_trgobj_mass ;
  int    Jpsi_e2_trgobj_q ;   
  float  Jpsi_e2_trgobj_dR ;
  //
  //
  int DoubleMu_fired=0;
  std::vector<int> DoubleEle_fired{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  std::vector<float> dieleobj1_pt{-999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999.};
  std::vector<float> dieleobj1_eta{-999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999.};
  std::vector<float> dieleobj1_phi{-999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999.};
  std::vector<float> dieleobj2_pt{-999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999.};
  std::vector<float> dieleobj2_eta{-999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999.};
  std::vector<float> dieleobj2_phi{-999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999.};
  
  std::vector<int> ele1_matchedDiEle{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  std::vector<int> ele2_matchedDiEle{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  std::vector<float> ele1_matchedDiEle_pt{-999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999.};
  std::vector<float> ele2_matchedDiEle_pt{-999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999.};
  std::vector<float> ele1_matchedDiEle_eta{-999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999.};
  std::vector<float> ele2_matchedDiEle_eta{-999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999.};
  std::vector<float> ele1_matchedDiEle_phi{-999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999.};
  std::vector<float> ele2_matchedDiEle_phi{-999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999., -999.};
  
  float ele1_matchedDiEle_dR  = 0;
  float ele2_matchedDiEle_dR  = 0;
  // float ele1_matchedDiEle_pt  = 0;
  // float ele2_matchedDiEle_pt  = 0;
  // float ele1_matchedDiEle_eta = 0;
  // float ele2_matchedDiEle_eta = 0;
  // float ele1_matchedDiEle_phi = 0;
  // float ele2_matchedDiEle_phi = 0;
  //
  //
  float  Jpsi_fit_pt;
  float  Jpsi_nonfit_pt;
  float  Jpsi_fit_eta;
  float  Jpsi_nonfit_eta;
  float  Jpsi_fit_phi;
  float  Jpsi_nonfit_phi;
  float  Jpsi_fit_mass;
  float  Jpsi_nonfit_mass;
  float  Jpsi_fit_vprob;
  float  Jpsi_electronsDr;


  // Not for tree, other variables declaration
  helperfunc aux;
  float chi = 0.;
  float ndf = 0.;

  // float dieleobj1_pt  = -999.;
  // float dieleobj1_eta = -999.;
  // float dieleobj1_phi = -999.;
  // float dieleobj2_pt  = -999.;
  // float dieleobj2_eta = -999.;
  // float dieleobj2_phi = -999.;
  TH1F * hist; 

}; // end of class member


NanoAnalyzerDoubleMuDebug::NanoAnalyzerDoubleMuDebug(const edm::ParameterSet& iConfig): 
  bFieldToken_(esConsumes<MagneticField, IdealMagneticFieldRecord>())
{
  electronToken_           = consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"));

  verticeToken_       = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));

  triggerToken_	      = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("HLT"));
  triggerobjectToken_ = consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("triggerobjects"));

  l1EG_               = consumes<l1t::EGammaBxCollection>(iConfig.getParameter<edm::InputTag>("l1EG"));

  hist = fs->make<TH1F>("cutflow", "cutflow", 10,0,10);
  tree_ = fs->make<TTree>( "tree", "tree" );

  createBranch();

} // end of constructor

NanoAnalyzerDoubleMuDebug::~NanoAnalyzerDoubleMuDebug() { }

void NanoAnalyzerDoubleMuDebug::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  reset();

  const auto& bField = iSetup.getData(bFieldToken_);

  iEvent.getByToken(verticeToken_, vertices_ );
  iEvent.getByToken(triggerToken_, HLTtriggers_);
  iEvent.getByToken(triggerobjectToken_ , triggerObjects); 

  nvtx = vertices_->size();

  run = (iEvent.id()).run();
  event = (iEvent.id()).event();
  luminosityBlock = (iEvent.id()).luminosityBlock();

  // check if the reference trigger and/or the analysis trigger fired
  bool isTriggered = false;
  const edm::TriggerNames& trigNames = iEvent.triggerNames(*HLTtriggers_);
  std::string refTriggerName="";    
  std::string finalTriggerName="";  

  // pT thresholds for di-ele trigger
  std::vector<std::string> pt_thr_v_string{"10", "9p5", "9", "8p5", "8", "7p5", "7", "6p5", "6", "5p5", "5", "4p5", "4"};

  // set of DoubleMu reference triggers
  std::vector<std::string> dmu_triggers{"HLT_DoubleMu4_3_Bs",
				      "HLT_DoubleMu4_3_Jpsi",
				      "HLT_DoubleMu4_3_LowMass",
				      "HLT_DoubleMu4_LowMass_Displaced",
				      "HLT_Mu0_L1DoubleMu",
				      "HLT_Mu4_L1DoubleMu",
				      "HLT_DoubleMu3_Trk_Tau3mu",
				      "HLT_DoubleMu3_TkMu_DsTau3Mu",
				      "HLT_DoubleMu4_MuMuTrk_Displaced",
				      "HLT_DoubleMu4_Jpsi_Displaced",
				      "HLT_DoubleMu4_Jpsi_NoVertexing",
				      "HLT_DoubleMu4_JpsiTrkTrk_Displaced",
				      "HLT_DoubleMu4_JpsiTrk_Bc*"
				      };

  for (unsigned int i = 0, n = HLTtriggers_->size(); i < n; ++i) {

    // std::cout << "Loop over all triggers: i = " << i << ", name = " << trigNames.triggerName(i) << ", fired = " << HLTtriggers_->accept(i) << std::endl;
   
    // Check if the reference trigger path fired
    for(int j=0; j<int(pt_thr_v_string.size()); j++){
      if(trigNames.triggerName(i).find(dmu_triggers[j])!= std::string::npos){             // chiara
	if(HLTtriggers_->accept(i)){
	  isTriggered = true;
	  refTriggerName = trigNames.triggerName(i);  
	  DoubleMu_fired = 1;
	}
      }
    }
    
    // Check if any of the di-ele trigger paths fired
    for(int j=0; j<int(pt_thr_v_string.size()); j++){
      std::string pt_thr_string = pt_thr_v_string[j];
      std::string DoubleEleTrigName = "HLT_DoubleEle"+pt_thr_string+"_eta1p22_mMax6";
      if(trigNames.triggerName(i).find(DoubleEleTrigName)!= std::string::npos){
	      if(HLTtriggers_->accept(i)){
	        DoubleEle_fired[j]=1;
	        //std::cout << "Event = " << event << " : DoubleEle_fired[j] fired for j = " << j << std::endl;
	      }
      }
    }
  }

  /*
  std::cout << std::endl;
  std::cout << "Ref trigger? " << isTriggered << std::endl;
  std::cout << "AN trigger? "  << DoubleEle_fired << std::endl;
  std::cout << "finalTriggerName = " << finalTriggerName << std::endl;
  std::cout << std::endl;
  */

  // Our reference path fired
  if(!isTriggered) return;    
  hist->Fill(0);


  // Take the electrons collection
  iEvent.getByToken(electronToken_ , electrons_    );

  // Output collections
  std::vector<pat::Electron> electroncollection;                          // collection with offline electrons passing minimal selection
  std::vector<pat::TriggerObjectStandAlone> trg_obj_collection;   // collection with HLT candidate with best match with offline electrons
                                                                  // (passing matching criteria) - at most one per electron
  std::vector<int> electronmatched;                                   // integer with the position in the trg_obj_collection of the HLT object matched to electron
  // 
  electroncollection.clear();
  trg_obj_collection.clear();
  electronmatched.clear();


  // Offline electrons
  for(size_t ielectron = 0; ielectron < electrons_->size(); ielectron++){
    const pat::Electron & electron = (*electrons_)[ielectron];

    // std::cout << "Offline: " << electron.pt() << " " << electron.eta() << " " << electron.phi() << std::endl;

    // Offline cuts
    if(electron.pt() < 2.5) continue;
    hist->Fill(1);
    if(fabs(electron.eta()) > 2.4) continue;
    hist->Fill(2);
    if(!(electron.gsfTrack().isNonnull())) continue;
    hist->Fill(3);
    const reco::GsfTrackRef gsfTrk = electron.gsfTrack();
    
    // This is to further process 
    electroncollection.push_back(electron);

  } // Loop over offline electrons

  
  if (electroncollection.size() < 2) return; 
  hist->Fill(4);


  // Prepare offline electron pairs
  float jpsi_max_pt = -1;
  int mcidx_e1 = -1;
  int mcidx_e2 = -1;
  TLorentzVector jpsi_tlv_highest;

  for(int ie = 0; ie < (int)electroncollection.size(); ie++){
    for(int je = ie+1; je < (int)electroncollection.size(); je++){
            

      const pat::Electron ele1 = electroncollection[ie];
      const pat::Electron ele2 = electroncollection[je];

      TLorentzVector tlv_e1;
      TLorentzVector tlv_e2;
      tlv_e1.SetPtEtaPhiM(ele1.pt(), ele1.eta(), ele1.phi(), aux.mass_electron);
      tlv_e2.SetPtEtaPhiM(ele2.pt(), ele2.eta(), ele2.phi(), aux.mass_electron);

      TLorentzVector tlv_jpsi = (tlv_e1 + tlv_e2);
      float jpsi_mass = tlv_jpsi.M();
      float jpsi_pt = tlv_jpsi.Pt();

      std::cout << "jpsi_mass = " << jpsi_mass << ", jpsi_max_pt = " << jpsi_max_pt << ", jpsi_pt = " << jpsi_pt << std::endl;      
      std::cout << "event = " << event << "; In the loop: mcidx_e1 = " << mcidx_e1 << ", mcidx_e2 = " << mcidx_e2 << std::endl;

      if (ele1.charge() + ele2.charge() !=0) continue;
      if (jpsi_mass < 2.0) continue; 
      if (jpsi_mass > 4.0) continue;
      
      if(jpsi_max_pt < jpsi_pt){
        jpsi_max_pt = jpsi_pt;
        mcidx_e1 = ie;
        mcidx_e2 = je;
        jpsi_tlv_highest = tlv_jpsi;
      }
    }
  }

  // At least 1 reco J/psi
  if(jpsi_max_pt == -1) return;  
  hist->Fill(6);


  // Loop over di-ele path, when fired
  for(int j=0; j<int(pt_thr_v_string.size()); j++){

    if (DoubleEle_fired[j]==1) {   

      std::string pt_thr_string = pt_thr_v_string[j];
      std::string DoubleEleTrigName = "HLT_DoubleEle"+pt_thr_string+"_eta1p22_mMax6";
      std::string DoubleEleObjName  = "hltDoubleEle"+pt_thr_string+"eta1p22mMax6ValidHitsFilter"; 
    
      // Loop over trigger objects matching the dielectron path 
      for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
      
	// check if this obj comes from the wanted di-electron path
	obj.unpackPathNames(trigNames);
	obj.unpackFilterLabels(iEvent, *HLTtriggers_);
	std::vector<std::string> pathNamesAll = obj.pathNames(false);
	bool isPathExist = false;
	for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h) {
	  // std::cout << "event " << event << " => path " << h << " => " << pathNamesAll[h] << std::endl;
	  if(pathNamesAll[h].find(DoubleEleTrigName)!= std::string::npos) {
	    isPathExist = true;   
	    // std::cout << "event " << event << " => found" << std::endl;
	  }
	}
	if(!isPathExist) continue;
	
	// std::cout << "ok, this object matches our diele path " << DoubleEleTrigName << std::endl;
	
	// check if the object is from the correct filter
	for (unsigned hh = 0; hh < obj.filterLabels().size(); ++hh){	
	  // std::cout << "event " << event << " => hh = " << hh << ", obj.filterLabels()[hh] = " << obj.filterLabels()[hh] << std::endl;
	  if(obj.filterLabels()[hh].find(DoubleEleObjName) != std::string::npos) {
	    if (dieleobj1_pt[j]<0) {
	      dieleobj1_pt[j]  = obj.pt();
	      dieleobj1_eta[j] = obj.eta();
	      dieleobj1_phi[j] = obj.phi();
	    } else if (dieleobj2_pt[j]<0) {
	      dieleobj2_pt[j]  = obj.pt();
	      dieleobj2_eta[j] = obj.eta();
	      dieleobj2_phi[j] = obj.phi();
	    }
	  }
	}
    }
   }
  } // path fired


  // Kinematic fit to offline electrons
  float reco_e1_pt  = electroncollection[mcidx_e1].pt(); 
  float reco_e1_eta = electroncollection[mcidx_e1].eta(); 
  float reco_e1_phi = electroncollection[mcidx_e1].phi(); 
  float reco_e2_pt  = electroncollection[mcidx_e2].pt(); 
  float reco_e2_eta = electroncollection[mcidx_e2].eta(); 
  float reco_e2_phi = electroncollection[mcidx_e2].phi(); 
  TVector3 ele1TV3, ele2TV3;
  ele1TV3.SetPtEtaPhi( reco_e1_pt, reco_e1_eta, reco_e1_phi );
  ele2TV3.SetPtEtaPhi( reco_e2_pt, reco_e2_eta, reco_e2_phi );
  // std::cout << "Offline1: " << reco_e1_pt << " " << reco_e1_eta << " " << reco_e1_phi << std::endl;
  // std::cout << "Offline2: " << reco_e2_pt << " " << reco_e2_eta << " " << reco_e2_phi << std::endl;

  // Kin fit
  const reco::TransientTrack electron1TT((*(electroncollection[mcidx_e1].bestTrack())),&bField);  
  const reco::TransientTrack electron2TT((*(electroncollection[mcidx_e2].bestTrack())),&bField);
  KinematicParticleFactoryFromTransientTrack pFactory;
  std::vector<RefCountedKinematicParticle> electronParticles;
  electronParticles.push_back(pFactory.particle(electron1TT, aux.electron_mass, chi, ndf, aux.electron_sigma));
  electronParticles.push_back(pFactory.particle(electron2TT, aux.electron_mass, chi, ndf, aux.electron_sigma));

  // Kinematic fit of the two electrons to a common vtx
  RefCountedKinematicParticle jpsi_part;
  RefCountedKinematicVertex jpsi_vertex;
  RefCountedKinematicTree jpTree;
  Bool_t jpsifit_flag;
  std::tie(jpsifit_flag, jpsi_part, jpsi_vertex, jpTree) = aux.KinematicFit(electronParticles, -1, -1);
  
  // Successfull kin fit
  if(!jpsifit_flag) return;
  hist->Fill(7);

  Jpsi_fit_pt    = jpsi_part->currentState().globalMomentum().perp();
  Jpsi_fit_eta   = jpsi_part->currentState().globalMomentum().eta();
  Jpsi_fit_phi   = jpsi_part->currentState().globalMomentum().phi();
  Jpsi_fit_mass  = jpsi_part->currentState().mass();
  Jpsi_fit_vprob = TMath::Prob(jpsi_part->chiSquared(), jpsi_part->degreesOfFreedom());


  // Match analysis HLT / offline
  float ele11_matchedDiEle_dR  = -999.;
  float ele11_matchedDiEle_pt  = -999.;
  float ele11_matchedDiEle_eta = -999.;
  float ele11_matchedDiEle_phi = -999.;
  //
  float ele21_matchedDiEle_dR  =  999.;
  float ele21_matchedDiEle_pt  = -999.;
  float ele21_matchedDiEle_eta = -999.;
  float ele21_matchedDiEle_phi = -999.;
  //
  float ele12_matchedDiEle_dR  =  999.;
  float ele12_matchedDiEle_pt  = -999.;
  float ele12_matchedDiEle_eta = -999.;
  float ele12_matchedDiEle_phi = -999.;
  //
  float ele22_matchedDiEle_dR  =  999.;
  float ele22_matchedDiEle_pt  = -999.;
  float ele22_matchedDiEle_eta = -999.;
  float ele22_matchedDiEle_phi = -999.;

    // Match HLT / offline
  for(int j=0; j<int(pt_thr_v_string.size()); j++){
    
    if (dieleobj1_pt[j]>=0) {
      float obj1_pt  = dieleobj1_pt[j];
      float obj1_eta = dieleobj1_eta[j];
      float obj1_phi = dieleobj1_phi[j];
      TVector3 obj1TV3;
      obj1TV3.SetPtEtaPhi( obj1_pt, obj1_eta, obj1_phi );
      Float_t deltaEta11 = fabs(reco_e1_eta - obj1_eta);
      Float_t deltaPhi11 = fabs(ele1TV3.DeltaPhi(obj1TV3));
      Float_t deltaEta21 = fabs(reco_e2_eta - obj1_eta);
      Float_t deltaPhi21 = fabs(ele2TV3.DeltaPhi(obj1TV3));
      if (deltaEta11<0.07 && deltaPhi11<0.2) {
	ele1_matchedDiEle[j] = 1;      
	ele1_matchedDiEle_pt[j]  = dieleobj1_pt[j];
	ele1_matchedDiEle_eta[j] = dieleobj1_eta[j];
	ele1_matchedDiEle_phi[j] = dieleobj1_phi[j];
	//std::cout << "Ok match offline1 / online1 " << std::endl;
      }
      if (deltaEta21<0.07 && deltaPhi21<0.2) {
	ele2_matchedDiEle[j] = 1; 
	ele2_matchedDiEle_pt[j]  = dieleobj1_pt[j];
	ele2_matchedDiEle_eta[j] = dieleobj1_eta[j];
	ele2_matchedDiEle_phi[j] = dieleobj1_phi[j];
	//std::cout << "Ok match offline2 / online1 " << std::endl;
      }
    }

    if (dieleobj2_pt[j]>=0) {
      float obj2_pt  = dieleobj2_pt[j];
      float obj2_eta = dieleobj2_eta[j];
      float obj2_phi = dieleobj2_phi[j];
      TVector3 obj2TV3;
      obj2TV3.SetPtEtaPhi( obj2_pt, obj2_eta, obj2_phi );
      Float_t deltaEta12 = fabs(reco_e1_eta - obj2_eta);
      Float_t deltaPhi12 = fabs(ele1TV3.DeltaPhi(obj2TV3));
      Float_t deltaEta22 = fabs(reco_e2_eta - obj2_eta);
      Float_t deltaPhi22 = fabs(ele2TV3.DeltaPhi(obj2TV3));
      if (deltaEta12<0.07 && deltaPhi12<0.2) {
	ele1_matchedDiEle[j] = 1; 
	ele1_matchedDiEle_pt[j]  = dieleobj2_pt[j];
	ele1_matchedDiEle_eta[j] = dieleobj2_eta[j];
	ele1_matchedDiEle_phi[j] = dieleobj2_phi[j];
	//std::cout << "Ok match offline1 / online2 " << std::endl;
      }
      if (deltaEta22<0.07 && deltaPhi22<0.2) {
	ele2_matchedDiEle[j] = 1; 
	ele2_matchedDiEle_pt[j]  = dieleobj2_pt[j];
	ele2_matchedDiEle_eta[j] = dieleobj2_eta[j];
	ele2_matchedDiEle_phi[j] = dieleobj2_phi[j];
	//std::cout << "Ok match offline2 / online2 " << std::endl;
      }
    }
  }




  // Match L1 / offline
  float bestMatchE1_eta = -99.;
  float bestMatchE1_phi = -99.;
  float bestMatchE1_pt  = -99.;
  float bestMatchE1_dR  = 999.;
  float bestMatchE2_eta = -99.;
  float bestMatchE2_phi = -99.;
  float bestMatchE2_pt  = -99.;
  float bestMatchE2_dR  = 999.;
  const auto &l1EG = iEvent.get(l1EG_); 
  for (l1t::EGammaBxCollection::const_iterator it = l1EG.begin(0); it != l1EG.end(0); it++) {
    pat::TriggerObjectStandAlone l1obj(it->p4());
    
    TVector3 l1objTV3;
    l1objTV3.SetPtEtaPhi( it->pt(), it->eta(), it->phi() );
    Float_t deltaRE1 = fabs(ele1TV3.DeltaR(l1objTV3));
    Float_t deltaRE2 = fabs(ele2TV3.DeltaR(l1objTV3));
    
    if (deltaRE1<bestMatchE1_dR){
      bestMatchE1_eta  = it->eta(); 
      bestMatchE1_phi  = it->phi(); 
      bestMatchE1_pt   = it->pt();    
      bestMatchE1_dR   = deltaRE1;
    }
    if (deltaRE2<bestMatchE2_dR){
      bestMatchE2_eta  = it->eta(); 
      bestMatchE2_phi  = it->phi(); 
      bestMatchE2_pt   = it->pt(); 
      bestMatchE2_dR   = deltaRE2;
    }
  }

  // Infos about JPsi candidate
  Jpsi_e1_pt   = electroncollection[mcidx_e1].pt();
  Jpsi_e1_eta  = electroncollection[mcidx_e1].eta();
  Jpsi_e1_phi  = electroncollection[mcidx_e1].phi();
  Jpsi_e1_mass = electroncollection[mcidx_e1].mass();
  Jpsi_e1_q    = electroncollection[mcidx_e1].charge();
  Jpsi_e1_isElectron = electroncollection[mcidx_e1].isElectron();
  Jpsi_e1_passMVA = electroncollection[mcidx_e1].electronID("mvaEleID-Fall17-noIso-V2-wpLoose");
  Jpsi_e1_bestL1pt   = bestMatchE1_pt;
  Jpsi_e1_bestL1eta  = bestMatchE1_eta;
  Jpsi_e1_bestL1phi  = bestMatchE1_phi;
  Jpsi_e1_bestL1dR   = bestMatchE1_dR;  

  Jpsi_e2_pt   = electroncollection[mcidx_e2].pt();
  Jpsi_e2_eta  = electroncollection[mcidx_e2].eta();
  Jpsi_e2_phi  = electroncollection[mcidx_e2].phi();
  Jpsi_e2_mass = electroncollection[mcidx_e2].mass();
  Jpsi_e2_q    = electroncollection[mcidx_e2].charge();
  Jpsi_e2_isElectron = electroncollection[mcidx_e2].isElectron();
  Jpsi_e2_passMVA = electroncollection[mcidx_e2].electronID("mvaEleID-Fall17-noIso-V2-wpLoose");
  Jpsi_e2_bestL1pt  = bestMatchE2_pt;
  Jpsi_e2_bestL1eta = bestMatchE2_eta;
  Jpsi_e2_bestL1phi = bestMatchE2_phi;
  Jpsi_e2_bestL1dR  = bestMatchE2_dR;  
  
  Jpsi_nonfit_pt   = jpsi_tlv_highest.Pt();
  Jpsi_nonfit_eta  = jpsi_tlv_highest.Eta();
  Jpsi_nonfit_phi  = jpsi_tlv_highest.Phi();
  Jpsi_nonfit_mass = jpsi_tlv_highest.M(); 

  Jpsi_electronsDr = ele1TV3.DeltaR(ele2TV3); 

  tree_->Fill();

  return;

} //NanoAnalyzer::analyze ends



//**************************************************
//************* additional methods *****************
//**************************************************

void NanoAnalyzerDoubleMuDebug::beginJob(const edm::ParameterSet& iConfig) { }

void NanoAnalyzerDoubleMuDebug::beginRun(const edm::Run &iRun, const edm::EventSetup &iStp) { }

void NanoAnalyzerDoubleMuDebug::fillDescriptions(edm::ConfigurationDescriptions & descriptions) { }

void NanoAnalyzerDoubleMuDebug::endRun(edm::Run const&, edm::EventSetup const&) { }

void NanoAnalyzerDoubleMuDebug::endJob() { }

//define this as a plug-in

// branch title creation
void NanoAnalyzerDoubleMuDebug::createBranch() { 

  tree_->Branch("run", &run, "run/i");
  tree_->Branch("event", &event, "event/l");
  tree_->Branch("luminosityBlock", &luminosityBlock, "luminosityBlock/i");
  tree_->Branch("nvtx", &nvtx, "nvtx/i");

  // Offline
  tree_->Branch("Jpsi_e1_pt", &Jpsi_e1_pt );
  tree_->Branch("Jpsi_e1_eta", &Jpsi_e1_eta );
  tree_->Branch("Jpsi_e1_phi", &Jpsi_e1_phi );
  tree_->Branch("Jpsi_e1_mass", &Jpsi_e1_mass );
  tree_->Branch("Jpsi_e1_q", &Jpsi_e1_q );
  tree_->Branch("Jpsi_e1_isElectron", &Jpsi_e1_isElectron );
  tree_->Branch("Jpsi_e1_passMVA"   , &Jpsi_e1_passMVA    );
  tree_->Branch("Jpsi_e1_bestL1pt"  , &Jpsi_e1_bestL1pt   );
  tree_->Branch("Jpsi_e1_bestL1eta" , &Jpsi_e1_bestL1eta  );
  tree_->Branch("Jpsi_e1_bestL1phi" , &Jpsi_e1_bestL1phi  );
  tree_->Branch("Jpsi_e1_bestL1dR",   &Jpsi_e1_bestL1dR );

  tree_->Branch("Jpsi_e2_pt", &Jpsi_e2_pt );
  tree_->Branch("Jpsi_e2_eta", &Jpsi_e2_eta );
  tree_->Branch("Jpsi_e2_phi", &Jpsi_e2_phi );
  tree_->Branch("Jpsi_e2_mass", &Jpsi_e2_mass );
  tree_->Branch("Jpsi_e2_q", &Jpsi_e2_q );
  tree_->Branch("Jpsi_e2_isElectron"   , &Jpsi_e2_isElectron    );
  tree_->Branch("Jpsi_e2_passMVA"   , &Jpsi_e2_passMVA    );
  tree_->Branch("Jpsi_e2_bestL1pt",   &Jpsi_e2_bestL1pt );
  tree_->Branch("Jpsi_e2_bestL1eta" , &Jpsi_e2_bestL1eta );
  tree_->Branch("Jpsi_e2_bestL1phi" , &Jpsi_e2_bestL1phi );
  tree_->Branch("Jpsi_e2_bestL1dR" ,  &Jpsi_e2_bestL1dR );

  // Analysis double-ele trigger
  tree_->Branch("Jpsi_e1_Diele_dR",     &ele1_matchedDiEle_dR );
  tree_->Branch("Jpsi_e1_Diele_pt",     &ele1_matchedDiEle_pt );
  tree_->Branch("Jpsi_e1_Diele_eta",    &ele1_matchedDiEle_eta );
  tree_->Branch("Jpsi_e1_Diele_phi",    &ele1_matchedDiEle_phi );

  tree_->Branch("Jpsi_e2_Diele_dR",     &ele2_matchedDiEle_dR );
  tree_->Branch("Jpsi_e2_Diele_pt",     &ele2_matchedDiEle_pt );
  tree_->Branch("Jpsi_e2_Diele_eta",    &ele2_matchedDiEle_eta );
  tree_->Branch("Jpsi_e2_Diele_phi",    &ele2_matchedDiEle_phi );

  tree_->Branch("DoubleMu_fired", &DoubleMu_fired );
  tree_->Branch("DoubleEle10_fired", &DoubleEle_fired[0] );
  tree_->Branch("DoubleEle9p5_fired", &DoubleEle_fired[1] );
  tree_->Branch("DoubleEle9_fired", &DoubleEle_fired[2] );
  tree_->Branch("DoubleEle8p5_fired", &DoubleEle_fired[3] );
  tree_->Branch("DoubleEle8_fired", &DoubleEle_fired[4] );
  tree_->Branch("DoubleEle7p5_fired", &DoubleEle_fired[5] );
  tree_->Branch("DoubleEle7_fired", &DoubleEle_fired[6] );
  tree_->Branch("DoubleEle6p5_fired", &DoubleEle_fired[7] );
  tree_->Branch("DoubleEle6_fired", &DoubleEle_fired[8] );
  tree_->Branch("DoubleEle5p5_fired", &DoubleEle_fired[9] );
  tree_->Branch("DoubleEle5_fired", &DoubleEle_fired[10] );
  tree_->Branch("DoubleEle4p5_fired", &DoubleEle_fired[11] );
  tree_->Branch("DoubleEle4_fired", &DoubleEle_fired[12] );
  
  // JPsi and daughters
  tree_->Branch("Jpsi_fit_pt",      &Jpsi_fit_pt );
  tree_->Branch("Jpsi_nonfit_pt",   &Jpsi_nonfit_pt );
  tree_->Branch("Jpsi_fit_eta",     &Jpsi_fit_eta );
  tree_->Branch("Jpsi_nonfit_eta",  &Jpsi_nonfit_eta );
  tree_->Branch("Jpsi_fit_phi",     &Jpsi_fit_phi );
  tree_->Branch("Jpsi_nonfit_phi",  &Jpsi_nonfit_phi);
  tree_->Branch("Jpsi_fit_mass",    &Jpsi_fit_mass );
  tree_->Branch("Jpsi_nonfit_mass", &Jpsi_nonfit_mass );
  tree_->Branch("Jpsi_fit_vprob",   &Jpsi_fit_vprob );
  tree_->Branch("Jpsi_electronsDr",     &Jpsi_electronsDr );
}

void NanoAnalyzerDoubleMuDebug::reset(void){

  run = -1;
  event = -1;
  luminosityBlock = -1;
  nvtx = -1;

  Jpsi_e1_pt = -99;
  Jpsi_e1_eta = -99;
  Jpsi_e1_phi = -99;
  Jpsi_e1_mass = -99;
  Jpsi_e1_q = -99;
  Jpsi_e1_isElectron = -99;
  Jpsi_e1_passMVA = -99;
  Jpsi_e1_bestL1pt = -99;
  Jpsi_e1_bestL1eta = -99;
  Jpsi_e1_bestL1phi = -99;
  Jpsi_e1_bestL1dR = -99;

  Jpsi_e2_pt = -99;
  Jpsi_e2_eta = -99;
  Jpsi_e2_phi = -99;
  Jpsi_e2_mass = -99;
  Jpsi_e2_q = -99;
  Jpsi_e2_isElectron = -99;
  Jpsi_e2_passMVA = -99;
  Jpsi_e2_bestL1pt = -99;
  Jpsi_e2_bestL1eta = -99;
  Jpsi_e2_bestL1phi = -99;
  Jpsi_e2_bestL1dR = -99;

  DoubleMu_fired = 0.;
  for (int ii=0; ii<13; ii++) DoubleEle_fired[ii] = 0.;

  for (int ii=0; ii<13; ii++) dieleobj1_pt[ii]  = -999.;
  for (int ii=0; ii<13; ii++) dieleobj1_eta[ii] = -999.;
  for (int ii=0; ii<13; ii++) dieleobj1_phi[ii] = -999.;
  for (int ii=0; ii<13; ii++) dieleobj2_pt[ii]  = -999.;
  for (int ii=0; ii<13; ii++) dieleobj2_eta[ii] = -999.;
  for (int ii=0; ii<13; ii++) dieleobj2_phi[ii] = -999.;

  for (int ii=0; ii<13; ii++) ele1_matchedDiEle[ii] = 0.;
  for (int ii=0; ii<13; ii++) ele2_matchedDiEle[ii] = 0.;
  for (int ii=0; ii<13; ii++) ele1_matchedDiEle_pt[ii]  = 0.;
  for (int ii=0; ii<13; ii++) ele2_matchedDiEle_pt[ii]  = 0.;
  for (int ii=0; ii<13; ii++) ele1_matchedDiEle_eta[ii] = 0.;
  for (int ii=0; ii<13; ii++) ele2_matchedDiEle_eta[ii] = 0.;
  for (int ii=0; ii<13; ii++) ele1_matchedDiEle_phi[ii] = 0.;
  for (int ii=0; ii<13; ii++) ele2_matchedDiEle_phi[ii] = 0.;

  ele1_matchedDiEle_dR = 0.;
  ele2_matchedDiEle_dR = 0.;

  Jpsi_fit_pt      = -99;
  Jpsi_nonfit_pt   = -99;
  Jpsi_fit_eta     = -99;
  Jpsi_nonfit_eta  = -99;
  Jpsi_fit_phi     = -99;
  Jpsi_nonfit_phi  = -99;
  Jpsi_fit_mass    = -99;
  Jpsi_nonfit_mass = -99;
  Jpsi_fit_vprob   = -99;
  Jpsi_electronsDr     = -99;

  // Not for tree
  // dieleobj1_pt  = -999.;
  // dieleobj1_eta = -999.;
  // dieleobj1_phi = -999.;
  // dieleobj2_pt  = -999.;
  // dieleobj2_eta = -999.;
  // dieleobj2_phi = -999.;
}

DEFINE_FWK_MODULE(NanoAnalyzerDoubleMuDebug);

// v.0.19.SS

#ifndef APAnalyzers_GammaJetAnalyzer_GammaJetAnalyzer_H
#define APAnalyzers_GammaJetAnalyzer_GammaJetAnalyzer_H



//
// Package:    GammaJetAnalyzer
// Class:      GammaJetAnalyzer
// 
// class GammaJetAnalyzer
//       GammaJetAnalyzer.cc 
//       APAnalyzers/GammaJetAnalyzer/src/GammaJetAnalyzer.cc
//
// Original Author:  Andrea Perrotta
// Additional Work:  Stefano Sinigardi
//         Created:  Thu Sep 27 17:00:00 CEST 2007
//




// system include files
#include <memory>

// user include files

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"

#include "TH1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"


class TFile;

namespace HepMC {
  class GenEvent;
  class GenParticle;
}

namespace reco {
  class Photon;
  class CaloJet;
  class GenJet;
  class Track;
}



//
// class declaration
//


class GammaJetAnalyzer : public edm::EDAnalyzer 
{
   public:
      explicit GammaJetAnalyzer(const edm::ParameterSet&);
      ~GammaJetAnalyzer();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      std::string pdgName (int pdgId);

      double deltaR(double eta1, double phi1, double eta2, double phi2);
      double deltaR(HepMC::GenParticle * , double eta2, double phi2);
      double deltaR(reco::CaloJet & , double eta2, double phi2);
      double deltaR(reco::GenJet & , double eta2, double phi2);
      double deltaR(reco::CaloJet & , reco::CaloJet &);
      double deltaR(const reco::Track & , double eta2, double phi2);

      HepMC::GenParticle * getMother(HepMC::GenParticle *) ;
      int getMotherId(HepMC::GenParticle *) ;
      HepMC::GenParticle * mcMatched(const HepMC::GenEvent *, double, double, double, bool);

      void printConeParticles (const HepMC::GenEvent *, double eta, double phi, double dimCone);
      int photonFlag (const HepMC::GenEvent *, HepMC::GenParticle *, double eta, double phi, double dimCone);
      bool detBTag (const HepMC::GenEvent *, double , double, double);

      // ----------member data ---------------------------
   private:
 
      std::string outputFile_;
      TFile*  rootFile_;
 
      std::string mcProducer_;
      std::string photonCollectionProducer_;
      std::string photonCorrCollectionProducer_;
      std::string uncorrectedPhotonCollection_;
      std::string vertexProducer_;
      std::string correctedPhotonCollection_;
      std::string CaloJetAlgorithm_ , GenJetAlgorithm_ ;
      std::string trackProducer_;

      //
      
      bool debug_particle;
      bool debug_photon;
      bool debug_jet;
      bool debug_;

      double theHiggsReferenceMass;
      double minimumPhotonPt;
      double maximumPhotonEta;
      double minimumDRPhotonJet;
      double minimumJetHadEnergyFraction;
      double maximumPhotonDR;
      double minimumParticlePt;
      int maximumChargedInCone;
      double ptThresholdInCone;
      double minimumLJetPt, minimumLJetEnergy, maximumLJetEta , minimumLJetDeltaEta, minimumLJetMass;
      double minimumBJetPt, minimumBJetEnergy, maximumBJetEta;
      double minimumAnyJetPt, minimumAnyJetEnergy, maximumAnyJetEta;
      double minimumJetPhotonDeltaR;

      //
      
      edm::ESHandle<ParticleDataTable> Particle_Data_Table;

      // Counters

      int nEvtTot, nEvtPho, nEvtJet, nEvtSel;

  
      // -- Histograms and root trees (if any...) 

      TH1F*  hPhoton_e;
      TH1F*  hPhoton_pt;
      TH1F*  hPhoton_eta;
      TH1F*  hPhoton_phi;
      TH1F*  hPhoton_r9;
      TH1F*  hMatchedPhoton_recEoverTrueE;

      TH1F*  hJet_pt;
      TH1F*  hJet_eta;
      TH1F*  hJet_phi;
      TH1F*  hNJets;

      TH1F*  hDrPhoton;
      TH1F*  hDrPizero;


      TTree *theTree;

      Int_t run;            // nb: queste due sono riempite all'inizio dell'analyzer
      Int_t evt;
     
      Float_t photon_pt;         
      Float_t photon_eta;
      Float_t photon_phi;
      Int_t photon_pdgId;
      Int_t photon_momId;
      Int_t photon_flag;
      Int_t photon_nTk;
      Float_t photon_r9;

      Float_t photonJet_pt;
      Float_t photonJet_eta;
      Float_t photonJet_phi;
      Float_t photonJet_eecal;
      Float_t photonJet_ehcal;

      Int_t nComb;
      Int_t nJets;
  
      Float_t mbb;
      Float_t mbbgam;

      Float_t jet_pt[4];        
      Float_t jet_eta[4];
      Float_t jet_phi[4];
      Float_t jet_eecal[4];
      Float_t jet_ehcal[4];
      Int_t jet_btag[4];


};

#endif

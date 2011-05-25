// v. 0.19.SS


// class GammaJetAnalyzer 
//       GammaJetAnalyzer.cc 
//       APAnalyzers/GammaJetAnalyzer/src/GammaJetAnalyzer.cc
//
// Original Author:  Andrea Perrotta
// Additional Work:  Stefano Sinigardi
//         Created:  Thu Sep 27 17:00:00 CEST 2007
//

#include "APAnalyzers/GammaJetAnalyzer/interface/GammaJetAnalyzer.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"

#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

// Generator level + CLHEP
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "HepMC/GenEvent.h"
#include "HepMC/GenParticle.h"

#include "CLHEP/Units/PhysicalConstants.h"
#include "TFile.h"

enum {PHOTON, NEUTRAL_HAD, NEUTRAL_JET, CHARGED_JET};


GammaJetAnalyzer::GammaJetAnalyzer(const edm::ParameterSet& ps)
{

  // Event samples

  photonCollectionProducer_ = ps.getParameter<std::string>("phoProducer");
  photonCorrCollectionProducer_ = ps.getParameter<std::string>("corrPhoProducer");
  uncorrectedPhotonCollection_ = ps.getParameter<std::string>("uncorrectedPhotonCollection");
  correctedPhotonCollection_ = ps.getParameter<std::string>("correctedPhotonCollection");

  trackProducer_ = ps.getParameter<std::string>("trackProducer");

  mcProducer_ = ps.getParameter<std::string>("mcProducer");

  vertexProducer_ = ps.getParameter<std::string>("primaryVertexProducer");
 
  CaloJetAlgorithm_ = ps.getParameter<std::string>("CaloJetAlgorithm") ;
  GenJetAlgorithm_ = ps.getParameter<std::string>("GenJetAlgorithm") ;

  // Output files

  outputFile_   = ps.getParameter<std::string>("outputFile");
  rootFile_ = TFile::Open(outputFile_.c_str(),"RECREATE"); // open output file to store histograms

  // Debug options

  debug_particle = ps.getParameter<bool>("debug_particle");
  debug_photon   = ps.getParameter<bool>("debug_photon");
  debug_jet      = ps.getParameter<bool>("debug_jet");
  debug_  =  debug_particle || debug_photon || debug_jet;

  // Qua iniziano i parametri per customizzare l'analisi:

  minimumPhotonPt =  ps.getParameter<double>("minimumPhotonPt");
  maximumPhotonEta = ps.getParameter<double>("maximumPhotonEta");

  minimumDRPhotonJet = ps.getParameter<double>("minimumDRPhotonJet");

  minimumJetHadEnergyFraction = ps.getParameter<double>("minimumHadronicRatioInJets");
  maximumPhotonDR = ps.getParameter<double>("coneAroundPhoton");

  minimumParticlePt = ps.getParameter<double>("minimumParticlePt");
  maximumChargedInCone = ps.getParameter<int>("maximumChargedInCone");
  ptThresholdInCone = ps.getParameter<double>("ptThresholdInCone");

  minimumLJetPt =  ps.getParameter<double>("minimumLJetPt");
  minimumBJetPt =  ps.getParameter<double>("minimumBJetPt");
  minimumLJetEnergy = ps.getParameter<double>("minimumLJetEnergy");
  minimumBJetEnergy = ps.getParameter<double>("minimumBJetEnergy");
  maximumLJetEta = ps.getParameter<double>("maximumLJetEta");
  maximumBJetEta = ps.getParameter<double>("maximumBJetEta");
  minimumLJetDeltaEta = ps.getParameter<double>("minimumLJetDeltaEta");
  minimumLJetMass = ps.getParameter<double>("minimumLJetMass");
  minimumAnyJetPt = (minimumLJetPt<minimumBJetPt?minimumLJetPt:minimumBJetPt);
  minimumAnyJetEnergy = (minimumLJetEnergy<minimumBJetEnergy?minimumLJetEnergy:minimumBJetEnergy);
  maximumAnyJetEta = (maximumLJetEta>maximumBJetEta?maximumLJetEta:maximumBJetEta);

  minimumJetPhotonDeltaR = ps.getParameter<double>("minimumJetPhotonDeltaR");
}




GammaJetAnalyzer::~GammaJetAnalyzer()
{
  delete rootFile_;
}




// -------------------------- ANALYZER ----------------------------

void GammaJetAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;                           // needed for all framework related classes

  run =  iEvent.id().run();
  evt =  iEvent.id().event();
  nEvtTot++;

  if (debug_) std::cout << "*** Analyzing event # " << nEvtTot 
			<< " : run " << run << " , event " << evt << "\n";

//
// =================PHOTONS=================== 
//

  /*
  // Get the uncorrected  photon collection
  Handle<reco::PhotonCollection> uncorrectedPhotonHandle;
  iEvent.getByLabel(photonCollectionProducer_, uncorrectedPhotonCollection_ , uncorrectedPhotonHandle);
  const reco::PhotonCollection phoCollection = *(uncorrectedPhotonHandle.product());
  */


  // Get the corrected  photon collection
  Handle<reco::PhotonCollection> correctedPhotonHandle;
  iEvent.getByLabel(photonCorrCollectionProducer_, correctedPhotonCollection_ , correctedPhotonHandle);
  const reco::PhotonCollection corrPhoCollection = *(correctedPhotonHandle.product());

  // Get the primary event vertex
  Handle<reco::VertexCollection> vertexHandle;
  iEvent.getByLabel(vertexProducer_, vertexHandle);
  reco::VertexCollection vertexCollection = *(vertexHandle.product());
  math::XYZPoint vtx(0.,0.,0.);
  if (vertexCollection.size()>0) vtx = vertexCollection.begin()->position();

  //Get the CaloJet collection
  Handle<reco::CaloJetCollection> caloJets;
  iEvent.getByLabel( CaloJetAlgorithm_, caloJets );


  /*
  // Get the GenJet collection
  Handle<reco::GenJetCollection> genJets;
  iEvent.getByLabel( GenJetAlgorithm_, genJets );
  */


  // Get the MC truth
  Handle< HepMCProduct > hepProd ;
  iEvent.getByLabel( mcProducer_.c_str(),  hepProd ) ;
  const HepMC::GenEvent * myGenEvent = hepProd->GetEvent();  

  if (debug_particle) myGenEvent->print(); 


  int photonNumber=0;


  //
  // Loop over corrected Photon candidates -> qui inizia l'analisi!
  //

  photonNumber=0;
  reco::PhotonCollection vectorCandInterestingPhotons;
  //  reco::CaloJetCollection vectorJetsInterestingPhotons;
  std::vector<reco::CaloJetCollection::const_iterator> vectorJetsInterestingPhotons;
  std::vector<int> vectorFlagInterestingPhotons;

  HepMC::GenParticle *mcMatchedPhoton = 0;
  HepMC::GenParticle *mcMatchedParticle = 0;

  reco::PhotonCollection::const_iterator pho;
  for( reco::PhotonCollection::const_iterator pho = corrPhoCollection.begin(); pho != corrPhoCollection.end(); pho++) 
  {
    // Set event vertex
    reco::Photon localPho(*pho);
    localPho.setVertex(vtx);


    //  Fill histos

    hPhoton_e->Fill( localPho.energy() );
    hPhoton_pt->Fill( localPho.pt() );
    hPhoton_eta->Fill( localPho.eta() );
    hPhoton_phi->Fill( localPho.phi() );
    hPhoton_r9->Fill( localPho.r9() );

    double ptPhoton=localPho.pt();
    double etaPhoton=localPho.eta();
    double phiPhoton=localPho.phi();

    bool onlyPhoton = true;
    mcMatchedPhoton = mcMatched(myGenEvent,etaPhoton,phiPhoton,maximumPhotonDR,onlyPhoton);
    onlyPhoton = false;
    mcMatchedParticle = mcMatched(myGenEvent,etaPhoton,phiPhoton,maximumPhotonDR,onlyPhoton);
 
    int momId = 0;
    if (mcMatchedPhoton) 
    {
      momId = getMotherId(mcMatchedPhoton);
      hMatchedPhoton_recEoverTrueE -> Fill( localPho.energy()/ mcMatchedPhoton->momentum().e() );
    }


//  DEBUG
    if (debug_photon) {
      printConeParticles (myGenEvent, etaPhoton, phiPhoton,maximumPhotonDR);
      if (mcMatchedPhoton)      {
	momId = getMotherId(mcMatchedPhoton);
	double ptgen = mcMatchedPhoton->momentum().px()*mcMatchedPhoton->momentum().px() +
	  mcMatchedPhoton->momentum().py()*mcMatchedPhoton->momentum().py();
	if (ptgen>0) ptgen = sqrt(ptgen); else ptgen =0.;
	std::cout << "RecoPH #" << photonNumber 
		  << " (" << localPho.pt() << ", " << localPho.eta() << ", " << localPho.phi() 
		  <<") matched with: GenPH #" << mcMatchedPhoton->barcode()  
		  << " (" << ptgen << ", " << mcMatchedPhoton->momentum().eta() << ", " 
		  << mcMatchedPhoton->momentum().phi()  << ") - deltaR: " << deltaR(mcMatchedPhoton,localPho.eta(),localPho.phi())
		  << " - mother: " << pdgName(momId) << std::endl;
      }
      else      {
	std::cout << "Unable to match a GenPH  with RecoPH #" << photonNumber << std::endl;
      }
      
      if (mcMatchedParticle)    {
	if (mcMatchedParticle->pdg_id() != 22)   {
	  double ptgen = mcMatchedParticle->momentum().px()*mcMatchedParticle->momentum().px() +
	    mcMatchedParticle->momentum().py()*mcMatchedParticle->momentum().py();
 	  if (ptgen>0) ptgen = sqrt(ptgen); else ptgen =0.;
	  std::cout << "RecoPH #" << photonNumber 
		    << " (" << localPho.pt() << ", " << localPho.eta() << ", " << localPho.phi() 
		    << ") matched with: GenPart #" << mcMatchedParticle->barcode() 
		    << " (" << ptgen << ", " << mcMatchedParticle->momentum().eta() << ", " 
		    << mcMatchedParticle->momentum().phi()  << ") - deltaR: " << deltaR(mcMatchedParticle,localPho.eta(),localPho.phi())
		    << " - pdgId: " << mcMatchedParticle->pdg_id() << " (" << pdgName(mcMatchedParticle->pdg_id())
		    << ")" << std::endl;
	}
      }
      else      {
	std::cout << "\tUnable to match this RecoPH with GenPARTICLE #" << photonNumber << std::endl;
      }
    }
    
// END DEBUG


    photonNumber++;
    if ( ptPhoton > minimumPhotonPt  &&  etaPhoton > -maximumPhotonEta  &&  etaPhoton < maximumPhotonEta ) {

      double drMin = 999.;
      double drTwo = 999.;
      double ptOne = 0.;

      reco::CaloJetCollection::const_iterator  thePhotonJet;

      for(reco::CaloJetCollection::const_iterator cal = caloJets->begin(); cal != caloJets->end(); cal++)  {
	reco::CaloJet localJet(*cal);
        double ptJet = localJet.pt();
	double dr = deltaR(localJet,etaPhoton,phiPhoton);
        if (dr<drMin) {
          if (ptOne>minimumAnyJetPt) drTwo = drMin;
          drMin = dr;
          ptOne = ptJet;
          thePhotonJet = cal;
	}
	else if (dr<drTwo) {
	  if (ptJet>minimumAnyJetPt)  drTwo = dr;
        }
      }

      hDrPhoton -> Fill(drTwo);
      if ( drTwo > minimumDRPhotonJet ) {
	vectorCandInterestingPhotons.push_back(localPho);
        vectorJetsInterestingPhotons.push_back(thePhotonJet);
	// photonFlag riempie anche l'histo con la DR per pi0, eta, etc..
 	int theFlag = photonFlag (myGenEvent, mcMatchedParticle, etaPhoton, phiPhoton, maximumPhotonDR);
	vectorFlagInterestingPhotons.push_back(theFlag);
      }
    }


  }  //  End loop over reconstructed corrected photons

  if (vectorCandInterestingPhotons.size() <= 0)  return;
  nEvtPho++;

  //
  // Decidi quale fotone tenere per l'analisi:
  //


  std::vector<reco::CaloJetCollection::const_iterator>::const_iterator cal=vectorJetsInterestingPhotons.begin();
  std::vector<int>::const_iterator flg=vectorFlagInterestingPhotons.begin();
  int thePhotonFlag = 99;
  double thePhotonPt = 0.;

  reco::Photon thePhoton;
  reco::CaloJetCollection::const_iterator thePhotonJet;

  for (pho = vectorCandInterestingPhotons.begin(); pho != vectorCandInterestingPhotons.end(); pho++) {
    int thisFlag = ( ( (*flg)==PHOTON || (*flg)==NEUTRAL_HAD )? PHOTON : (*flg) );
    if ( thisFlag <= thePhotonFlag  &&  (*pho).pt()>thePhotonPt ) {
      thePhotonPt = (*pho).pt();
      thePhotonFlag = (*flg);
      thePhoton = (*pho);
      thePhotonJet = (*cal);
    }
    cal++;
    flg++;
  }
  double etaPhoton = thePhoton.eta();
  double phiPhoton = thePhoton.phi();

  // Count the number of reconstructed tracks around the selected photon:
  int ntrks = 0;
  Handle<reco::TrackCollection> theTracks;
  iEvent.getByLabel(trackProducer_,theTracks);
  for ( reco::TrackCollection::const_iterator trk=theTracks->begin(); trk!=theTracks->end(); ++trk) {
    if ( (*trk).pt() > minimumParticlePt ) {
      if ( deltaR((*trk),etaPhoton,phiPhoton) < maximumPhotonDR ) {
	ntrks++;
      }
    }
  }
    

//
// =================JETS===================
//

  //
  // Loop on Calo Jets
  //
  int jetInd = 0;
  
  std::vector<reco::CaloJet> vectorCandInterestingCaloJets;
  std::vector<reco::CaloJet> vectorCandInterestingCaloBJets;

  for( reco::CaloJetCollection::const_iterator cal = caloJets->begin(); cal != caloJets->end(); cal++ ) 
  {

    reco::CaloJet localJet(*cal);

    // DEBUG 
    if (debug_jet) {
      std::cout << "CALO JET #" << jetInd << std::endl << localJet.print() << std::endl;
      if (cal == thePhotonJet) 	{
	std::cout << "*** THIS IS THE PHOTON JET ***" << std::endl;
        std::cout << "*** PHOTON pt, eta phi : " << thePhoton.pt() << " , " << thePhoton.eta() << " , " << thePhoton.phi() << std::endl;
      }
      jetInd++;
    }
    // END DEBUG

    if (cal == thePhotonJet) {
      continue;
    }

    // Fill vectors with interesting jets:

    double ptJet = localJet.pt();
    double eJet = localJet.energy();
    double etaJet = localJet.eta();
    if ( localJet.energyFractionHadronic() > minimumJetHadEnergyFraction ) {
      if ( ptJet > minimumAnyJetPt
           && (etaJet > -maximumAnyJetEta) &&  (etaJet < maximumAnyJetEta)
	   && eJet > minimumAnyJetEnergy  )
	vectorCandInterestingCaloJets.push_back(localJet);
    }

    hJet_pt->Fill(ptJet);
    hJet_eta->Fill(etaJet);
    hJet_phi->Fill( localJet.phi() );
  }



  /*
  // 
  // Loop on GenJets
  //

  jetInd = 0;

  for( reco::GenJetCollection::const_iterator gen = genJets->begin(); gen != genJets->end(); gen++ ) 
  {

    reco::GenJet localJet(*gen);

    // DEBUG
    if (debug_jet) {
      std::cout << "GEN JET #" << jetInd << std::endl << localJet.print() << std::endl;
      jetInd++;
    }
    // END DEBUG

    //    h_ptGen_ ->Fill( localJet.pt() );   
    //    h_etaGen_->Fill( localJet.eta() );
    //    h_phiGen_->Fill( localJet.phi() );
    

  }

  */



  //
  // E qua utilizziamo i fotoni e i jets che abbiamo selezionato
  //
  
  int nJetsAccepted = vectorCandInterestingCaloJets.size();
  hNJets -> Fill ((float) nJetsAccepted);

  if ( nJetsAccepted < 4 )  return;
  nEvtJet++;

  int interestingEvent=0;
  
  std::vector<reco::CaloJet>::const_iterator cal1, cal2, bcal1, bcal2;

  for (cal1 = vectorCandInterestingCaloJets.begin(); cal1 != vectorCandInterestingCaloJets.end(); cal1++)  {
    reco::CaloJet localLJet1(*cal1);
    if ( localLJet1.pt() < minimumLJetPt
	 ||  localLJet1.eta() < -maximumLJetEta  ||  localLJet1.eta() > maximumLJetEta
	 ||  localLJet1.energy() < minimumLJetEnergy  )  continue;  

    for (cal2 = cal1+1; cal2 != vectorCandInterestingCaloJets.end(); cal2++)   {
      reco::CaloJet localLJet2(*cal2);
      if ( localLJet2.pt() < minimumLJetPt
	   ||  localLJet2.eta() < -maximumLJetEta  ||  localLJet2.eta() > maximumLJetEta
	   ||  localLJet2.energy() < minimumLJetEnergy  )  continue;  
      double deltaLJetEta = localLJet1.eta()-localLJet2.eta();
      if (deltaLJetEta < 0) deltaLJetEta = -deltaLJetEta;
      double invMassLJet = (localLJet1.p4()+localLJet2.p4()).mass();
      bool okTagJets = deltaLJetEta > minimumLJetDeltaEta &&
                       invMassLJet > minimumLJetMass;
      if (okTagJets)    {

	for (bcal1 = vectorCandInterestingCaloJets.begin(); bcal1 != vectorCandInterestingCaloJets.end(); bcal1++)  { 
          if (bcal1==cal1 || bcal1==cal2) continue;
	  reco::CaloJet localBJet1(*bcal1);
	  if ( localBJet1.pt() < minimumBJetPt
	       ||  localBJet1.eta() < -maximumBJetEta  ||  localBJet1.eta() > maximumBJetEta
	       ||  localBJet1.energy() < minimumBJetEnergy  )  continue;  

	  for (bcal2 = bcal1+1; bcal2 != vectorCandInterestingCaloJets.end(); bcal2++)  {
	    if (bcal2==cal1 || bcal2==cal2) continue;
	    reco::CaloJet localBJet2(*bcal2);
	    if ( localBJet2.pt() < minimumBJetPt
		 ||  localBJet2.eta() < -maximumBJetEta  ||  localBJet2.eta() > maximumBJetEta
		 ||  localBJet2.energy() < minimumBJetEnergy  )  continue;  

	    //	    double etaPhoton = thePhoton.eta();
	    //	    double phiPhoton = thePhoton.phi();

	    double deltaR_Pho_BJet1 = deltaR(localBJet1,etaPhoton,phiPhoton);
	    double deltaR_Pho_BJet2 = deltaR(localBJet2,etaPhoton,phiPhoton);
	    double deltaR_Pho_LJet1 = deltaR(localLJet1,etaPhoton,phiPhoton);
	    double deltaR_Pho_LJet2 = deltaR(localLJet2,etaPhoton,phiPhoton);
 
            bool okEvent = deltaR_Pho_BJet1 > minimumJetPhotonDeltaR && deltaR_Pho_BJet2 > minimumJetPhotonDeltaR &&
	                   deltaR_Pho_LJet1 > minimumJetPhotonDeltaR && deltaR_Pho_LJet2 > minimumJetPhotonDeltaR    ;
	    if (okEvent) { 
	      
	      interestingEvent++;

	      photon_pt = thePhoton.pt();
              photon_eta = thePhoton.eta();
              photon_phi = thePhoton.phi();
	      if (mcMatchedParticle) photon_pdgId = mcMatchedParticle->pdg_id();
              else photon_pdgId = 0;
	      if (mcMatchedPhoton) photon_momId = (getMother(mcMatchedPhoton))->pdg_id();
              else photon_momId = 0;
              photon_flag = thePhotonFlag;
	      photon_nTk = ntrks;
              photon_r9 = thePhoton.r9();

              photonJet_pt = thePhotonJet->pt();
              photonJet_eta = thePhotonJet->eta();
              photonJet_phi = thePhotonJet->phi();
              photonJet_eecal = thePhotonJet->emEnergyFraction();
              photonJet_ehcal = thePhotonJet->energyFractionHadronic();

	      nComb = interestingEvent;
              nJet = vectorCandInterestingCaloJets.size();

              jet_pt[0] = localLJet1.pt();
              jet_eta[0] = localLJet1.eta();
              jet_phi[0] = localLJet1.phi();
              jet_eecal[0] = localLJet1.emEnergyFraction();
              jet_ehcal[0] = localLJet1.energyFractionHadronic();
              jet_btag[0] = detBTag(myGenEvent, localLJet1.eta(), localLJet1.phi(), maximumPhotonDR);

              jet_pt[1] = localLJet2.pt();
              jet_eta[1] = localLJet2.eta();
              jet_phi[1] = localLJet2.phi();
              jet_eecal[1] = localLJet2.emEnergyFraction();
              jet_ehcal[1] = localLJet2.energyFractionHadronic();
              jet_btag[1] = detBTag(myGenEvent, localLJet2.eta(), localLJet2.phi(), maximumPhotonDR);
	      
               jet_pt[2] = localBJet1.pt();
              jet_eta[2] = localBJet1.eta();
              jet_phi[2] = localBJet1.phi();
              jet_eecal[2] = localBJet1.emEnergyFraction();
              jet_ehcal[2] = localBJet1.energyFractionHadronic();
              jet_btag[2] = detBTag(myGenEvent, localBJet1.eta(), localBJet1.phi(), maximumPhotonDR);

              jet_pt[3] = localBJet2.pt();
              jet_eta[3] = localBJet2.eta();
              jet_phi[3] = localBJet2.phi();
              jet_eecal[3] = localBJet2.emEnergyFraction();
              jet_ehcal[3] = localBJet2.energyFractionHadronic();
              jet_btag[3] = detBTag(myGenEvent, localBJet2.eta(), localBJet2.phi(), maximumPhotonDR);

	      // Cosi' riempie tutte le combinazioni possibili
              // Poi, nell'analisi col root tree dobbiamo, ad esempio, scegliere solo gli eventi con
              //       "interestingEvent == 1" per prendere solo la prima combinazione.
	      theTree->Fill();

	    }
	  }
	}
      }
    }
  }

  if (interestingEvent>0) {
    // Se c'e' piu' di una combinazione gli facciamo riempire solo l'ultima?
    nEvtSel++;
    //    theTree->Fill();
  }



  // DEBUG
  if (debug_) {
    std::cout << "\nThere are " << caloJets->size() << " calo Jets " << std::endl;
    //    std::cout << "There are " << genJets->size() << " generator level Jets" << std::endl;

    if (vectorCandInterestingCaloJets.size() != 1)    {
      std::cout << "There are " << vectorCandInterestingCaloJets.size() << " interesting hadronic jets" << std::endl;
    }
    else    {
      std::cout << "There is " << vectorCandInterestingCaloJets.size() << " interesting hadronic jet" << std::endl;
    }

    if (vectorCandInterestingPhotons.size() != 1)  {
      std::cout << "There are " << vectorCandInterestingPhotons.size() << " interesting photons" << std::endl;
    }
    else  {
      std::cout << "There is " << vectorCandInterestingPhotons.size() << " interesting photon" << std::endl;
    }
 
    if (interestingEvent != 0)  {
      std::cout << "This event has " << interestingEvent << " combinations interesting for the analysis" << std::endl;
    }
    else if (vectorCandInterestingCaloJets.size() != 0 && vectorCandInterestingPhotons.size() != 0)   {
      std::cout << "This event has not passed cuts even if it had some good starting points" << std::endl;
    }

    std::cout << "\n\n" << std::endl;

  }
  // END DEBUG


  vectorCandInterestingCaloJets.clear();
  vectorCandInterestingCaloBJets.clear();
  vectorCandInterestingPhotons.clear();
  vectorJetsInterestingPhotons.clear();
  vectorFlagInterestingPhotons.clear();
}








HepMC::GenParticle*  GammaJetAnalyzer::getMother(HepMC::GenParticle * gp)
{
// 
// Trova la mamma del mcMatchedPhoton navigando in GenEvent
//
    HepMC::GenParticle* mom = 0;
    int PartId = gp->pdg_id();
    int MomId = 0;

    if (  gp->production_vertex() ) {
      if ( gp->production_vertex()->particles_begin(HepMC::parents) !=
	   gp->production_vertex()->particles_end(HepMC::parents)  )	{
	mom = *( gp->production_vertex()->particles_begin(HepMC::parents) ) ;
	MomId = mom->pdg_id();
      }
      // Let try once more if MomId == PartId...
      if (MomId == PartId && mom->production_vertex() ) {
	if ( mom->production_vertex()->particles_begin(HepMC::parents) !=
	     mom->production_vertex()->particles_end(HepMC::parents) )	{
	  mom = *( mom->production_vertex()->particles_begin(HepMC::parents) ) ;
	  MomId = mom->pdg_id();
	}
      }
    }
    return mom;
}

int GammaJetAnalyzer::getMotherId(HepMC::GenParticle * gp)
{
    int MomId = -1;
    HepMC::GenParticle* mom = getMother(gp);
    if (mom)  MomId = mom->pdg_id();
    return MomId;
}


std::string GammaJetAnalyzer::pdgName (int pdgId) 
{
  std::string partname="NoName";

  const ParticleData* pdt = Particle_Data_Table->particle(pdgId);
  if (pdt)
    partname = pdt->name();
  else
    throw edm::Exception( edm::errors::InvalidReference ) 
      << std::endl << "*** WARNING: particle PDG Id code = " << pdgId << " has no particle data" << std::endl;

  return partname;
}

double GammaJetAnalyzer::deltaR(double eta1, double phi1, double eta2, double phi2)
{
  double deltaEta = eta1-eta2;
  deltaEta=pow(deltaEta,2);
  double deltaPhi = phi1-phi2;
  if ( deltaPhi > pi )  deltaPhi -= twopi;
  if ( deltaPhi < -pi) deltaPhi += twopi;
  deltaPhi=pow(deltaPhi,2);
  double dr = sqrt( deltaPhi+deltaEta);
  return dr;
}

double GammaJetAnalyzer::deltaR(HepMC::GenParticle * mcp, double eta2, double phi2)
{
  double eta1 = mcp->momentum().eta();
  double phi1 = mcp->momentum().phi();
  return deltaR(eta1,phi1,eta2,phi2);
}

double GammaJetAnalyzer::deltaR(reco::CaloJet & jet, double eta2, double phi2)
{
  double eta1 = jet.momentum().eta();
  double phi1 = jet.momentum().phi();
  return deltaR(eta1,phi1,eta2,phi2);
}

double GammaJetAnalyzer::deltaR(reco::CaloJet & jet1, reco::CaloJet & jet2)
{
  double eta1 = jet1.momentum().eta();
  double phi1 = jet1.momentum().phi();
  double eta2 = jet1.momentum().eta();
  double phi2 = jet1.momentum().phi();
  return deltaR(eta1,phi1,eta2,phi2);
}

double GammaJetAnalyzer::deltaR(reco::GenJet & jet, double eta2, double phi2)
{
  double eta1 = jet.momentum().eta();
  double phi1 = jet.momentum().phi();
  return deltaR(eta1,phi1,eta2,phi2);
}

double GammaJetAnalyzer::deltaR(const reco::Track & trk, double eta2, double phi2)
{
  double eta1 = trk.eta();
  double phi1 = trk.phi();
  return deltaR(eta1,phi1,eta2,phi2);
}



HepMC::GenParticle *  GammaJetAnalyzer::mcMatched(const HepMC::GenEvent *myGE,
						  double eta, double phi, double dimCone, bool onlyPhoton)
{
  double minDeltaPhoton = dimCone;
  HepMC::GenParticle * mcp=0;
  for ( HepMC::GenEvent::particle_const_iterator p = myGE->particles_begin(); p != myGE->particles_end(); p++ )  {
    if ((*p)->status()== 1 && ( !onlyPhoton || (*p)->pdg_id() == 22 ) ) {
      double deltaPhoton =  deltaR((*p),eta,phi);
      if ( deltaPhoton < minDeltaPhoton )  {
	minDeltaPhoton=deltaPhoton;
	mcp = (*p);
      }
    }
  }
  return mcp;
}


void GammaJetAnalyzer::printConeParticles(const HepMC::GenEvent *myGE, double eta, double phi, double dimCone)
{
  for ( HepMC::GenEvent::particle_const_iterator p = myGE->particles_begin(); p != myGE->particles_end(); p++ )    {
    if ( (*p)->status()== 1 )  {
      double deltaParticle = deltaR((*p),eta,phi);
      if (deltaParticle <= dimCone)  {
	std::cout << "GenPart #" << (*p)->barcode() << " : " << (*p)->pdg_id() << " , E = " <<(*p)->momentum().e()
		  << " , deltaR = " << deltaParticle << std::endl;
      }
    }
  }
}


int GammaJetAnalyzer::photonFlag(const HepMC::GenEvent *myGE, HepMC::GenParticle * mcp, double eta, double phi, double dimCone)
{
  int theFlag = -1;
  if (!mcp) return theFlag;

  HepMC::GenParticle * mcMother = 0;
  HepMC::GenParticle * mcSister = 0;
  int pdgId = mcp->pdg_id();
  int momId = -1;
  if (pdgId==22) {
    mcMother = getMother(mcp);
    momId = mcMother->pdg_id();
    for  ( HepMC::GenVertex::particle_iterator p = mcp->production_vertex()->particles_begin(HepMC::descendants);
           p != mcp->production_vertex()->particles_end(HepMC::descendants); p++) {
      if ( (*p)->status()==1 && (*p)->pdg_id()==22 && (*p)!=mcp ) {
	mcSister = (*p);
        hDrPizero -> Fill(deltaR(mcSister,eta,phi));
      }
    }
  }

  int nCha = 0;
  double ptSum = 0;
  //  std::cout << "*** In photonFlag(): " << mcp->barcode() << " , " << pdgId << " ( " << momId << " ) " << std::endl;
  for ( HepMC::GenEvent::particle_const_iterator p = myGE->particles_begin(); p != myGE->particles_end(); p++ )    {
    if ( (*p)->status()== 1 )  {
      double deltaParticle = deltaR((*p),eta,phi);
      if (deltaParticle <= dimCone)  {
        double ptParticle = sqrt( pow((*p)->momentum().px(),2) + pow((*p)->momentum().py(),2) );
        if (ptParticle > minimumParticlePt) {
	  if ( Particle_Data_Table->particle((*p)->pdg_id())->charge() != 0 ) nCha++;
          if ( (*p)!=mcp && (*p)!=mcSister ) ptSum += ptParticle;
	  //	  std::cout << "*** GenPart #" << (*p)->barcode() << " : " << (*p)->pdg_id() << " , pT = " << ptParticle
	  //		    << " , deltaR = " << deltaParticle << " ( ncha = " << nCha << " ) " << std::endl;
	}
      }
    }
  }

  if (nCha > maximumChargedInCone)  theFlag = CHARGED_JET;
  else if (ptSum > ptThresholdInCone) theFlag = NEUTRAL_JET;
  else if (pdgId |= 22 || 
	                  ( momId >= 100 && momId < 1000) ) theFlag = NEUTRAL_HAD;
  else theFlag = PHOTON;
  // in ogni caso salviamo la pdgId e la mamma nel root tree per ri-flaggare, nel caso, sull'output finale.

  //  std::cout << "*** In photonFlag(), nCha, ptSum = " << nCha << " , " << ptSum << " => flag = " << theFlag << std::endl;
  return theFlag;
}


bool GammaJetAnalyzer::detBTag(const HepMC::GenEvent *myGE, double eta, double phi, double dimCone)
{
    bool bTagged = false;

    for ( HepMC::GenEvent::particle_const_iterator p = myGE->particles_begin(); p != myGE->particles_end(); p++ ) 
    {
      double deltaParticle = deltaR((*p), eta, phi);

      if (deltaParticle <= dimCone)
      {
        int pdgB = abs((*p)->pdg_id());
        pdgB = pdgB/100;
        int overBMeson = pdgB % 10;
        bool bMeson = (overBMeson==5?true:false);

        pdgB = pdgB/10;
        int overBBaryon = pdgB % 10;
        bool bBaryon = (overBBaryon==5?true:false);

        if (bMeson || bBaryon) {
	  bTagged = true;
          return bTagged;
	}
      }
    }
    return bTagged;
}


// -------------------------- BEGINJOB -----------------------------

void  GammaJetAnalyzer::beginJob(const edm::EventSetup& iSetup)
{
  // Initialize counters:
  nEvtTot = 0;
  nEvtPho = 0;
  nEvtJet = 0;
  nEvtSel = 0;

  // Initialize the ParticleDataTable
  iSetup.getData(Particle_Data_Table);

  // go to *OUR* rootfile and book histograms
  rootFile_  -> cd();


  theTree = new TTree("T","GammaJetAnalyzer root tree",0);

  theTree->Branch("run",&run,"run/I");
  theTree->Branch("evt",&evt,"evt/I");

  theTree->Branch("photon_pt",&photon_pt,"photon_pt/F");
  theTree->Branch("photon_eta",&photon_eta,"photon_eta/F");
  theTree->Branch("photon_phi",&photon_phi,"photon_phi/F");
  theTree->Branch("photon_pdgId",&photon_pdgId,"photon_pdgId/I");
  theTree->Branch("photon_momId",&photon_momId,"photon_momId/I");
  theTree->Branch("photon_flag",&photon_flag,"photon_flag/I");
  theTree->Branch("photon_nTk",&photon_nTk,"photon_nTk/I");
  theTree->Branch("photon_r9",&photon_r9,"photon_r9/F");

  theTree->Branch("photonJet_pt",&photonJet_pt,"photonJet_pt/F");
  theTree->Branch("photonJet_eta",&photonJet_eta,"photonJet_eta/F");
  theTree->Branch("photonJet_phi",&photonJet_phi,"photonJet_phi/F");
  theTree->Branch("photonJet_eecal",&photonJet_eecal,"photonJet_eecal/F");
  theTree->Branch("photonJet_ehcal",&photonJet_ehcal,"photonJet_ehcal/F");

  theTree->Branch("jet_pt",&jet_pt,"jet_pt[4]/F");
  theTree->Branch("jet_eta",&jet_eta,"jet_eta[4]/F");
  theTree->Branch("jet_phi",&jet_phi,"jet_phi[4]/F");
  theTree->Branch("jet_eecal",&jet_eecal,"jet_eecal[4]/F");
  theTree->Branch("jet_ehcal",&jet_ehcal,"jet_ehcal[4]/F");
  theTree->Branch("jet_btag",&jet_ehcal,"jet_btag[4]/I");


  
  hPhoton_e                      = new TH1F("Photon_e","Photons : Energy ",200,0.,400.);
  hPhoton_pt                     = new TH1F("Photon_pt","Photons:  p_{T} ",200,0.,400.);
  hPhoton_eta                    = new TH1F("Photon_eta","Photons: #eta ",200,-2.5, 2.5);
  hPhoton_phi                    = new TH1F("Photon_phi","Photons: #phi ",200,-pi, pi);
  hPhoton_r9                     = new TH1F("Photon_r9","Photons: R9 ",200, 0., 2.);
  hMatchedPhoton_recEoverTrueE   = new TH1F("MatchedPhoton_recEoverTrueE","Matched Photons: recE / trueE ",200, 0., 10.);

  hJet_pt                        = new TH1F("Jet_pt","Jets: p_{T} ", 200, 0, 1000);
  hJet_eta                       = new TH1F("Jet_eta","Jets: #eta ", 200,-5. , 5.);
  hJet_phi                       = new TH1F("Jet_phi","Jets: #phi ", 200, -pi, pi);
  hNJets                         = new TH1F("nJets","Number of accepted Jets",21,-0.5,20.5);

  hDrPhoton                      = new TH1F("drPhoton","Photon isolation (#Delta R)",200,0.,10.);
  hDrPizero                      = new TH1F("drPizero","#Delta R in #pi",200,0.,10.);

}




// -------------------------- ENDJOB ----------------------------

void GammaJetAnalyzer::endJob() 
{

  std::cout << "*** -----------------------------------------------------" << std::endl;
  std::cout << "*** GAMMAJETANALYZER final report : " << std::endl;
  std::cout << "*** Numero totale di eventi:                             " << nEvtTot << std::endl;
  std::cout << "*** Numero di eventi con almeno un fotone selezionato:   " << nEvtPho << std::endl;
  std::cout << "*** Numero di eventi con un fotone e 4 jets selezionati: " << nEvtJet << std::endl;
  std::cout << "*** Numero di eventi che passano tutti i tagli:          " << nEvtSel << std::endl;
  std::cout << "*** -----------------------------------------------------" << std::endl;

  rootFile_                     -> cd();
 
  theTree                       -> Write();

  hPhoton_e                     -> Write();
  hPhoton_pt                    -> Write();
  hPhoton_eta                   -> Write();
  hPhoton_phi                   -> Write();
  hPhoton_r9                    -> Write();
  hMatchedPhoton_recEoverTrueE  -> Write();

  hJet_pt                       -> Write();
  hJet_eta                      -> Write();
  hJet_phi                      -> Write();
  hNJets                        -> Write();

  hDrPhoton                     -> Write();
  hDrPizero                     -> Write();

  rootFile_                     -> Close();
}




// -------------- Define this as a plug-in ---------------------------------

DEFINE_FWK_MODULE(GammaJetAnalyzer);

//=================================  Introduction  =================================//<
// -*- C++ -*-
//Mohamed Aly and Ahmed Abdelmotteleb
//Majorana Neutrino Sensitivity at ATLAS
//Rivet analysis code that runs over our data [signal + background] and applies cuts

//=================================  Libraries  =================================//
#include <iostream>
#include <fstream>
#include <ostream>
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/SmearedJets.hh"
#include "Rivet/Projections/SmearedParticles.hh"
#include "Rivet/Projections/SmearedMET.hh"

#include <math.h>

//==========================================================================//
//==========================================================================//
//=================================  MAIN  =================================//
//==========================================================================//
//==========================================================================//
namespace Rivet {

  //divided into constructor + 3 loops: init, analyze and finalize

  //==================================================================================//
  //=================================  Constructer ==================================//
  //=================================================================================//
  class MLRSM_mumujj : public Analysis {
  public:

    /// Constructor
    MLRSM_mumujj() : Analysis("MLRSM_mumujj") {   }

    //=================================  Declaring global variables  =================================//

    double sumW = 0;  //sum of weights
    double var_sumW = 0; //variance of sum of weights
    double fid_sumW = 0; //fiducial sum of weights
    double var_fid_sumW = 0; //variance of fiducial sum of weights
    //=================================  Declaring global functions  =================================//

    //=====LEPTONS
    //Lepton centrality
    double calcLC( const Jets& jets, const FourMomentum& lep ){

      //Define important jets
      // .at() checks if we are within bounds of our array number of elements
      // use .at() instead of [] because it will give an error if we exceed elements however [] won't
      FourMomentum J1 = jets.at(0).momentum(); //highest pT jet
      FourMomentum J2 = jets.at(1).momentum(); //second-highest pT jet

      double LCval = -999; //should be negative??
      if (jets.size() < 2 ) {
        //std::cout << "ERROR: should not allow calculation of Lepton Centrality on jets with less than multiplicity two! " << std::endl;
        return LCval;
      }

      else {  //centrality equation
        double LCval = fabs( ( lep.rapidity() - ( (J1.rapidity() + J2.rapidity() ) /2.0  ) ) / ( J1.rapidity() - J2.rapidity() ) );
        return LCval;
      }

    }

    //boolean that checks if your lepton passes the centrality cut (ATLAS detector)
    bool PassesLeptonCentralityCut( const Jets& jets, const Particle& lep ){
      if ( jets.size() < 2 ) return false;

      if ( calcLC(jets,lep) < 0.4 ) {
        return true;
      }

      else {
        return false;
      }
    }

    //=====JETS
    //Calculate Jet Centrality
    double calcJC( const Jets& jets ){
      //Define important jets
      FourMomentum J1 = jets.at(0).momentum();  //highest pT jet
      FourMomentum J2 = jets.at(1).momentum();  //second-highest pT jet
      double JCval = 999; //should be positive??
      if ( jets.size() < 3 ) {
        //std::cout << "ERROR: should not allow calculation of Jet Centrality on jets with less than multiplicity three! " << std::endl;
        return JCval;
      }

      else {
        FourMomentum J = jets.at(2).momentum();  //take the third jet (jets are always ordered in descending pT order)
        //calculating the centrality
        JCval = fabs(( J.rapidity() - ( (J1.rapidity() + J2.rapidity() ) /2.0  ) ) / ( J1.rapidity() - J2.rapidity() ));
        return JCval;
      }
    }

    //boolean that checks if your lepton (Jet? Ben 14/10/2020) passes the centrality cut (ATLAS detector)
    bool PassesJetCentralityCut( const Jets& jets ){
      if ( jets.size() < 3 ) return true;

      double JCval = calcJC(jets);
      if (JCval < 0.4) { 
        return false;
      }

      else {
        return true;
      }
    }

    //function to find jets between two jet rapidities --> returns true or false
    bool _isBetween(const Jet j, const Jet bj1, const Jet bj2) {
      double y_j = j.rapidity(); //probe jet
      double y_bj1 = bj1.rapidity(); //boundary jet 1
      double y_bj2 = bj2.rapidity(); //boundary jet 2

      double y_min = std::min(y_bj1, y_bj2); //boundary jet with lower rapidity
      double y_max = std::max(y_bj1, y_bj2); //boundary jet with higher rapidity

      if (y_j > y_min && y_j < y_max) return true; //check if our probe lies between the boundary jets
      else return false;
    }


    //function to find the number of gap jets
    int _getNumGapJets(const Jets jets, FourMomentum& thirdJet) {
      if (jets.size() < 2) return 0;

      // The vector of jets is already sorted by pT. So the boundary jets will be the first two.
      const Jet bj1 = jets.at(0);
      const Jet bj2 = jets.at(1);

      int n_gap_jets = 0;
      // Start loop at the 3rd hardest pT jet
      for (size_t i = 2; i < jets.size(); i++) {
        const Jet j = jets.at(i);
        // If this jet is between the boundary jets and is hard enough, increment counter
        if (_isBetween(j, bj1, bj2)) {
          if (n_gap_jets == 0) thirdJet = j.momentum();
          n_gap_jets++;
        }
      }
      return n_gap_jets;
    }

    //Event counter
    int event_number = 0;
    int num_events_before_a_cut=0;
    int num_events_after_a_cut=0;
    //====================================================================================//
    //=================================  Initialisation  =================================//
    //====================================================================================//

    //// Book histograms and initialise projections before the run
    void init() {
      // 0 == no object smearing, 1 == ATLAS Run2 smearing, 2 == CMS Run2 smearing
      SMEAR = 0;
      // 0 == no smearing, 1 == smear
      MET_SMEAR = 1;
      JET_SMEAR = 1;
      MUON_SMEAR = 1;

      //std::string outfile = Analysis::getOption<std::string>("outfile", 0);

      //std::string string = Analysis::getOption(string);

      //std::string outfile = getOption("outfile")
      //std::string outfile = "rivet_outfile";
      //if (getOption("OUTFILENAME") == "0") { outfile = outfile + "1.txt"; } 
      //if (getOption("OUTFILENAME") == "1") { outfile = outfile + "2.txt"; } 
      //if (getOption("OUTFILENAME") == "2") { outfile = outfile + "3.txt"; } 

      // Create output text file
      //outFile.open(outfile);
      //outFile.open("tmp_file.txt");
      //outFile << "test" << std::endl; 
      //outFile.close();

      outFile.open("rivet_output.txt");

      //// Initialise and register projections (declaring particles)
      //=================================  General Final States Declaration  =================================//
      ///define final state
      const FinalState fs(Cuts::abseta < 5);

      ///define visible final states (not pT miss)
      VisibleFinalState vfs(Cuts::abseta < 4.5);
      declare(vfs,"vfs");
    //=================================  MET Declaration  =================================//
      // MET, should remove muons/electrons that fail cuts from this as vetoed final state?
      MissingMomentum missing_et(vfs);
      declare(missing_et,"missing_et");

      // define MET smeared with ATLAS Run2
      SmearedMET met_smeared_atlas(missing_et, MET_SMEAR_ATLAS_RUN2);
      declare(met_smeared_atlas,"met_smeared_atlas");

      // define MET smeared with CMS Run2
      SmearedMET met_smeared_cms(missing_et, MET_SMEAR_CMS_RUN2);
      declare(met_smeared_cms,"met_smeared_cms");

    //=================================  Photon and Neutrino Decleration  =================================//
      //define photon
      IdentifiedFinalState photon(fs);
      photon.acceptIdPair(PID::PHOTON); 
      declare(photon,"photon");


      //Identify the Neutrinos
      IdentifiedFinalState neutrinos(fs);
      neutrinos.acceptNeutrinos();
      declare(neutrinos, "neutrinos");

    //=================================  Muons Decleration  =================================//
      //All Muons
      IdentifiedFinalState muon(fs);
      muon.acceptIdPair(PID::MUON);
      declare(muon, "muon");

      //Muons that are visible in a detector 
      VetoedFinalState vfs_muons(vfs);
      vfs_muons.addVetoPair(PID::MUON,Cuts::abseta > 2.7);
      declare(vfs_muons,"vfs_muons");

      MissingMomentum missing_et_muons(vfs_muons);
      declare(missing_et_muons,"missing_et_muons");

      // define MET smeared with ATLAS Run2
      SmearedMET met_smeared_atlas_muons(missing_et_muons, MET_SMEAR_ATLAS_RUN2);
      declare(met_smeared_atlas_muons,"met_smeared_atlas_muons");

      // define MET smeared with CMS Run2
      SmearedMET met_smeared_cms_muons(missing_et_muons, MET_SMEAR_CMS_RUN2);
      declare(met_smeared_cms_muons,"met_smeared_cms_muons");

      //define dressed muons as muons observed with a photon at radius 0.1 at  most from them
      DressedLeptons dressed_muons(photon, muon, 0.1, Cuts::abseta < 2.4 && Cuts::pT > 10*GeV);
      declare(dressed_muons, "dressed_muons");

      SmearedParticles muons_eff_smear_atlas(dressed_muons, MUON_EFF_ATLAS_RUN2, MUON_SMEAR_ATLAS_RUN2);
      declare(muons_eff_smear_atlas, "muons_eff_smear_atlas");

      SmearedParticles muons_eff_smear_cms(dressed_muons, MUON_EFF_CMS_RUN2, MUON_SMEAR_CMS_RUN2);
      declare(muons_eff_smear_cms, "muons_eff_smear_cms");

    //=================================  Electrons Decleration  =================================//
      /*//All Muons
      IdentifiedFinalState elec(fs);
      elec.acceptIdPair(PID::ELECTRON);
      declare(elec, "elec");
      //Muons that are visible in a detector       
      VetoedFinalState vfs_elecs(vfs);
      vfs_elecs.addVetoPair(PID::ELECTRON,Cuts::abseta > 2.5);
      declare(vfs_elecs,"vfs_elecs");

      MissingMomentum missing_et_elecs(vfs_elecs);
      declare(missing_et_elecs,"missing_et_elecs");

      // define MET smeared with ATLAS Run2
      SmearedMET met_smeared_atlas_elecs(missing_et_elecs, MET_SMEAR_ATLAS_RUN2);
      declare(met_smeared_atlas_elecs,"met_smeared_atlas_elecs");

      // define MET smeared with CMS Run2
      SmearedMET met_smeared_cms_elecs(missing_et_elecs, MET_SMEAR_CMS_RUN2);
      declare(met_smeared_cms_elecs,"met_smeared_cms_elecs");

      //define dressed electrons as muons observed with a photon at radius 0.1 at  most from them
      DressedLeptons dressed_elecs(photon, elec, 0.1, Cuts::abseta < 2.4 && Cuts::pT > 10*GeV);
      declare(dressed_elecs, "dressed_elecs");

      SmearedParticles elecs_eff_smear_atlas(dressed_elecs, ELECTRON_EFF_ATLAS_RUN2, ELECTRON_SMEAR_ATLAS_RUN2);
      declare(elecs_eff_smear_atlas, "elecs_eff_smear_atlas");

      SmearedParticles elecs_eff_smear_cms(dressed_elecs, ELECTRON_EFF_CMS_RUN2, ELECTRON_SMEAR_CMS_RUN2);
      declare(elecs_eff_smear_cms, "elecs_eff_smear_cms");*/

    //=================================  Jets Decleration  =================================//
      ///Use FastJet to select Jets 
      FastJets fj(fs, FastJets::ANTIKT, 0.4);
      declare(fj,"jets");

      // define jets smeared with ATLAS Run2
      SmearedJets jets_smeared_atlas(fj,JET_SMEAR_ATLAS_RUN2);
      declare(jets_smeared_atlas,"jets_smeared_atlas");

      // define jets smeared with CMS Run2
      SmearedJets jets_smeared_cms(fj,JET_SMEAR_CMS_RUN2);
      declare(jets_smeared_cms,"jets_smeared_cms");

    //=================================  OutFile -- CSV header  =================================//

      //write to csv file (note different from MSG_DEBUG)
      outFile << "weight,eTmiss,lep1_pT,lep2_pT,dilep_pT,jet1_pT,jet2_pT,dijet_pT,HT,mjj,mll,mlljj,delta_rap_jj,delta_rap_ll,rap1,rap2,delta_phi_jj,delta_phi_ll,delta_R_jj,delta_R_ll,delta_R_lljj,lep1_centrality,lep2_centrality,jets_centrality,dilep_centrality,ngapjets,njets,mt,eTmiss_in,eTmiss_out,rapl1,rapl2,delta_R_l1j1,delta_R_l2j2,delta_phi_met_j1,delta_phi_met_j2,nu_from_hadron,et_nu,delta_phi_met_jj,eTmiss_muons,MET,MET_smeared_atlas,MET_smeared_cms,MET_muons,MET_smeared_atlas_muons,MET_smeared_cms_muons" << std::endl;
    }

    //===============================================================================//
    //=================================  Analysis  =================================//
    //==============================================================================//

    /// Perform the per-event analysis
    //select particles, filter according to cuts, loop over combinations, construct observables, fill histograms.
    // This is where the per-event aspect of the analysis algorithm goes.

    ///Do the event by event analysis here (loop over every event)


    void analyze(const Event& event) {
      
      if(event_number % 1000 == 0){std::cout << "Proccessing event: " << event_number << std::endl;}

    //================================= VISIBLE FINAL STATE PARTICLES IN EVENT  =================================//
      //Get visible final state particles
      Particles vfs_particles = apply<VisibleFinalState>(event, "vfs").particles();
      //Visible final state particles with muon rapidity restriction
      Particles vfs_particles_muons = apply<VetoedFinalState>(event, "vfs_muons").particles();
      //Visible final state particles with electron rapidity restriction
/*      Particles vfs_particles_elecs = apply<VetoedFinalState>(event, "vfs_elecs").particles();
*/
      //array of leptons after applying second cut (first cut is below-->cand_electrons & cand_muons)
      Particles recon_leptons;

      //array of jets after applying first cut
      Jets cand_jets;

      //array of jets after applying second cut
      Jets recon_jets;

      //define a 4-momentum array for missing pT
      FourMomentum pTmiss;
      FourMomentum pTmiss_muons;
      FourMomentum pTmiss_elecs;

      //define gap jets
      FourMomentum gapjet(0., 0., 0., 0.);

      //=================== weight and variance ===================//

      //This will calculate the weight for each of the N events that rivet processes

      //weight of each event
      const double weight = event.weight();

      //summing over all the weights (sum of weights)
      sumW += weight;

      //Sum up the weights^2 --> This happens to be equal to the variance
      var_sumW += std::pow(weight,2);

    //================================= Building Candidate Leptons  -- Dressed =================================//
        //=================== Muons ===================//
      const Particles cand_muons = apply<DressedLeptons>(event, "dressed_muons").particlesByPt( Cuts::abseta < 2.4 && Cuts::pT > 10*GeV);
      const Particles& smeared_muons_atlas = apply<ParticleFinder>(event, "muons_eff_smear_atlas").particlesByPt();
      const Particles& smeared_muons_cms = apply<ParticleFinder>(event, "muons_eff_smear_cms").particlesByPt();       
        //=================== Electrons ===================//
/*      const Particles cand_elecs = apply<DressedLeptons>(event, "dressed_elecs").particlesByPt( Cuts::abseta < 2.4 && Cuts::pT > 10*GeV);
      const Particles& smeared_elecs_atlas = apply<ParticleFinder>(event, "elecs_eff_smear_atlas").particlesByPt();
      const Particles& smeared_elecs_cms = apply<ParticleFinder>(event, "elecs_eff_smear_cms").particlesByPt();       */   
        //=================== All Leptons -- Dressed ===================//
      Particles cand_leptons;
      foreach(const Particle& mu, cand_muons)  cand_leptons.push_back(mu);
/*      foreach(const Particle& e, cand_elecs)  cand_leptons.push_back(e);
*/

      //MSG_DEBUG("Number of dressed muons = " << cand_muons.size());

      //=================== Building Candidate Neutrinos and Hadronic Neutrino Flag  ===================//

    //=====Neutrinos 
    Particles neutrinos = applyProjection<IdentifiedFinalState>(event, "neutrinos").particlesByPt();
    //create a binary that is 0 if neutrino not from hadron, and 1 if neutrino comes form hadron
    int nu_from_hadron=0;
    foreach(Particle nu, neutrinos){
      if(nu.fromHadron()==true){
        nu_from_hadron=1;
      }
    }

    //=================== Building Candidate Jets ===================//

    //Create FastJets array with FastJets projection conditions defined earlier

    const Jets& truth_jets = apply<FastJets>(event, "jets").jetsByPt( (Cuts::abseta < 4.5) && (Cuts::pT > 20*GeV) );
    const Jets& jets_smeared_atlas = apply<JetAlg>(event, "jets_smeared_atlas").jetsByPt( (Cuts::abseta < 4.5) && (Cuts::pT > 20*GeV) );
    const Jets& jets_smeared_cms = apply<JetAlg>(event, "jets_smeared_cms").jetsByPt( (Cuts::abseta < 4.5) && (Cuts::pT > 20*GeV) );

    //Defining candidate Jets based on smearing setting
    if ( SMEAR==0 ){ foreach ( const Jet& jet, truth_jets ) { cand_jets.push_back(jet); } }
    if ( SMEAR==1 && JET_SMEAR ){ foreach ( const Jet& jet, jets_smeared_atlas ) { cand_jets.push_back(jet); } }
    if ( SMEAR==2 && JET_SMEAR ){ foreach ( const Jet& jet, jets_smeared_cms ) { cand_jets.push_back(jet); } }

    //MSG_DEBUG("Jet multiplicity = " << cand_jets.size()); //print for debug
    //=================== Choosing non-hadronic Leptons ===================//
    Particles non_hadronic_leps;
    foreach ( const Particle& lep, cand_leptons) {  //loop over candidate jets
      bool away_from_jet = true;  //define condition
        //deltaR is a distance defined by eta and phi
        if (lep.fromHadron()) {
          away_from_jet = false;
          break;
        }
      if ( away_from_jet ) {
        if ( SMEAR==0 ){
            non_hadronic_leps.push_back(lep);
        }
        else if ( SMEAR==1 && MUON_SMEAR ) {
            non_hadronic_leps.push_back(lep);
        }
        else if ( SMEAR==2 && MUON_SMEAR ) {
            non_hadronic_leps.push_back(lep);
        }
      }
    }    

    //=================== Choosing reconstructed jets ===================//
    //Use non-hadronic leptons to select jets that are far from prompt leptons 
    foreach ( const Jet& jet, cand_jets ) {  //loop over candidate jets
      bool away_from_lep = true;  //define condition
      foreach ( const Particle& lep, non_hadronic_leps) { //loop over reconstructed leptons
        //deltaR is a distance defined by eta and phi
        if (deltaR(lep.momentum(),jet.momentum()) <= 0.3) {
          away_from_lep = false;
          break;
        }
      }
      if ( away_from_lep ) recon_jets.push_back( jet );
    }

    //=================== Choosing and sorting reconstructed Leptons ===================//
      //Use non hadronic leptons to select leptons that are far enough from all other leptons
      foreach(const Particle& lep1, non_hadronic_leps){
        bool isolated_lep = true;
        foreach(const Particle& lep2, non_hadronic_leps){
          if(!isSame(lep1,lep2)){
            if (deltaR(lep1.momentum(),lep2.momentum()) < 0.2) {
              isolated_lep=false; 
              break;
            }
          }
        }
        if(isolated_lep) recon_leptons.push_back(lep1); 
      }  
      //Sort the leptons by pT
      std::sort (recon_leptons.begin(), recon_leptons.end(), [](Particle const &a, Particle const &b) { return a.momentum().pT() > b.momentum().pT(); });

    //=================== Leptons and Jets multiplicity vetoes ===================//
    // Ensure Exactly 2 leptons for each event
    if ( recon_leptons.size() != 2) vetoEvent;

    //Ensure we have at least 2 jets for each event
    if (recon_jets.size() < 2 ) vetoEvent;

    // Veto events with opposite charged leptons
    if(recon_leptons.at(0).pid()*recon_leptons.at(1).pid()<0) vetoEvent;

    //=================== Additional fixes on leptons and jets arrays ===================//

    // Reorder the leptons reconstructed so that higher pT electron is first
    if(recon_leptons.at(0).perp()<recon_leptons.at(1).perp())
    std::swap(recon_leptons.at(0),recon_leptons.at(1));

    //=================== Summing up weight of surviving events ===================//

    //This will calculate the weight for the N_selected events after some cuts. If an event i passes the cuts, its weight will be calculated here
    //and in the weight variable. If it didn't pass, it will only be in the all_weight variable.

    //deduce the weight of the event that made it here
    const double selected_weight = event.weight();

    //sum up the weights
    //fiducial means after applying cuts (our selected events)
    fid_sumW += selected_weight;

    //Calculate the variance
    var_fid_sumW += std::pow(selected_weight,2);

    //=================== Calculating Missing Transverse Energy (MET) ===================//

    const double MET = apply<MissingMomentum>(event, "missing_et").vectorEt().mod();
    const double MET_smeared_atlas = apply<SmearedMET>(event, "met_smeared_atlas").vectorEt().mod();
    const double MET_smeared_cms = apply<SmearedMET>(event, "met_smeared_cms").vectorEt().mod();

    const double MET_muons = apply<MissingMomentum>(event, "missing_et_muons").vectorEt().mod();
    const double MET_smeared_atlas_muons = apply<SmearedMET>(event, "met_smeared_atlas_muons").vectorEt().mod();
    const double MET_smeared_cms_muons = apply<SmearedMET>(event, "met_smeared_cms_muons").vectorEt().mod();


    //calculate MET from VFS, considerign only muons with (|y| < 2.7)
    foreach ( const Particle & p, vfs_particles_muons ) {
      pTmiss_muons -= p.momentum();
    }
    const double eTmiss_muons = pTmiss_muons.pT();
    //calculate MET from VFS, considering all muons with (|y| < 4.5)
    foreach ( const Particle & p, vfs_particles ) {
      pTmiss -= p.momentum();
    }

    const double eTmiss = pTmiss.pT();

    double eTmiss_in=0; // MET vals from MET vectors within |y| < 4.5
    double eTmiss_out=0; // MET vals from MET vectors within |y| > 4.5
    if (fabs( pTmiss.eta() ) <= 4.5 ) {
      eTmiss_in = pTmiss.pT();
    }
    else{
      eTmiss_out=pTmiss.pT();
    }
    //MET of the true Neutrinos 
    FourMomentum pt_nu;
    foreach(Particle& nu, neutrinos){
      pt_nu+=nu.momentum();
    }
    const double et_nu=pt_nu.pT();

//=================== Calculating the kinematic variables ===================//

//=====Defining some final arrays
//highest pT jet
FourMomentum j1 = recon_jets.at(0).momentum();

//second-highest pT jet
FourMomentum j2 = recon_jets.at(1).momentum();

//highest pt lepton
FourMomentum l1 = recon_leptons.at(0).momentum();

//second-highest lepton
FourMomentum l2 = recon_leptons.at(1).momentum();

//di-jet
FourMomentum diJet = (j1+j2);

//di-lepton
FourMomentum diLep = (l1 + l2);

const double lep1_centrality = calcLC(recon_jets,recon_leptons.at(0));  //Lepton 1 Centrality
MSG_DEBUG("Lepton 1 Centrality = " << lep1_centrality);
const double lep2_centrality = calcLC(recon_jets,recon_leptons.at(1));  // Lepton 2 Centrality
MSG_DEBUG("Lepton 2 Centrality = " << lep2_centrality);
const double jets_centrality = calcJC(recon_jets);  // Third Jet centrality
MSG_DEBUG("jets_centrality = " << jets_centrality);
const double dilep_centrality = calcLC(recon_jets,diLep); // Dilepton centrality
MSG_DEBUG("Dileptons_centrality = " << jets_centrality);
const double lep1_pT=recon_leptons.at(0).pT()/1000;  // Lepton 1 pT
MSG_DEBUG("lep1_pT = " << lep1_pT );
const double lep2_pT=recon_leptons.at(1).pT()/1000; // Lepton 2 pT
MSG_DEBUG("lep2_pT = " << lep2_pT);
const double jet1_pT=recon_jets.at(0).momentum().pT()/1000;  // Jet 1 pT
MSG_DEBUG("jet1_pT = " << jet1_pT);
const double jet2_pT=recon_jets.at(1).momentum().pT()/1000; // Jet 2 pT
MSG_DEBUG("jet2_pT = " << jet2_pT);
const double dilep_pT=(recon_leptons.at(0).pT()+recon_leptons.at(1).pT())/1000; // Dilepton pT
MSG_DEBUG("dilep_pT = " << dilep_pT);
const double dijet_pT=(jet1_pT+jet2_pT); // Dijets pT
MSG_DEBUG("dijet_pT = " << dijet_pT);
const double mjj = (j1+j2).mass()/1000; // Invariant mass of Dijets
MSG_DEBUG("mjj = " << mjj);
const double mll = (l1+l2).mass()/1000;  //Invariant Mass of Dileptons
MSG_DEBUG("mll = " << mll);
const double mlljj = (l1+l2+j1+j2).mass()/1000; //Invariant Mass of Final State
MSG_DEBUG("mlljj = " << mlljj);
//change the following back to fabs later
const double delta_rap_jj= j1.rapidity() - j2.rapidity();  // Difference in rapidity between jets
MSG_DEBUG("delta_rap_jj = " << delta_rap_jj);
const double delta_rap_ll= fabs(l1.rapidity() - l2.rapidity()); // Difference in rapidity between leptons
MSG_DEBUG("delta_rap_ll = " << delta_rap_ll);
const double rapj1=j1.rapidity(); // Rapidity of Jet 1
MSG_DEBUG("rap1 = " << rapj1);
const double rapj2=j2.rapidity();  // Rapidity of Jet 2
MSG_DEBUG("rap2 = " << rapj2);
const double rapl1=l1.rapidity(); // Rapidity of lepton 1
MSG_DEBUG("rap1 = " << rapl1);
const double rapl2=l2.rapidity();  // Rapidity of lepton 2
MSG_DEBUG("rap2 = " << rapl2);
const double delta_phi_jj = deltaPhi(j1, j2);  // Change in Phi between jets (maps difference to 0->PI)
MSG_DEBUG("Azimuthal Angular Seperation between the Jets " << delta_phi_jj);
const double delta_phi_ll = deltaPhi(l1, l2);  // Change in Phi between leptons (maps difference to 0->PI)
MSG_DEBUG("Azimuthal Angular Seperation between the leptons " << delta_phi_ll);
const double delta_phi_met_j1=deltaPhi(pTmiss,j1); // azimuthal angle between jet1 and MET
MSG_DEBUG("Azimuthal Angular Seperation between the jet1 and MET " << delta_phi_met_j1);
const double delta_phi_met_j2=deltaPhi(pTmiss,j2); // azimuthal angle between jet1 and MET
MSG_DEBUG("Azimuthal Angular Seperation between the jet2 and MET " << delta_phi_met_j2);
const double HT=fabs(jet1_pT+jet2_pT+ lep1_pT+lep2_pT); // Scalar Sum of Transverse momentum of final state
MSG_DEBUG("HT = " << HT);
const double delta_R_jj = sqrt( sqr(delta_rap_jj) + sqr(delta_phi_jj) ); // Delta R for jj system
MSG_DEBUG("delta_R_jj = " << delta_R_jj);
const double delta_R_ll = sqrt( sqr(delta_rap_ll) + sqr(delta_phi_ll) );  // Delta R for ll system
MSG_DEBUG("delta_R_ll = " << delta_R_ll);
const double delta_rap_l1j1= fabs(j1.rapidity() - l1.rapidity());
const double delta_rap_l2j2= fabs(j2.rapidity() - l2.rapidity());
const double delta_phi_l1j1= deltaPhi(j1,l1);
const double delta_phi_l2j2= deltaPhi(j2,l2);
const double delta_phi_met_jj= deltaPhi(pTmiss,diJet);
const double delta_R_l1j1 = sqrt( sqr(delta_rap_l1j1) + sqr(delta_phi_l1j1) );  // Delta R between l1 and j1
MSG_DEBUG("delta_R_l1j1 = " << delta_R_l1j1);
const double delta_R_l2j2 = sqrt( sqr(delta_rap_l2j2) + sqr(delta_phi_l2j2) );  // Delta R between l2 and j2
MSG_DEBUG("delta_R_l2j2 = " << delta_R_l2j2);

const double delta_R_lljj = deltaR(diJet,diLep);  //Delta R of the total final state
MSG_DEBUG("delta_R_lljj = " << delta_R_lljj);
const int ngapjets = _getNumGapJets(recon_jets, gapjet); // Number of gap Jets
MSG_DEBUG("ngapjets= " << ngapjets);
const int njets= recon_jets.size();  // Number of Jets
MSG_DEBUG("njets= " << njets);
const double mt_sq= pow(mlljj,2)+pow((j1+j2+l1+l2).pT(),0.5); //transverse mass squared
const double mt = pow(mt_sq,0.5);
MSG_DEBUG("Done with Calculations");

event_number++;
MSG_DEBUG("event number: " << event_number);
//Print the values of the kinematic variables
 outFile << selected_weight << "," << eTmiss << "," << lep1_pT << "," << lep2_pT << "," << dilep_pT << "," << jet1_pT  << "," << jet2_pT << "," << dijet_pT << "," << HT << "," << mjj << "," << mll << "," << mlljj << "," << delta_rap_jj << "," << delta_rap_ll << "," << rapj1 << ","<< rapj2 << "," << delta_phi_jj << "," << delta_phi_ll << "," << delta_R_jj << "," << delta_R_ll << "," << delta_R_lljj << "," << lep1_centrality << "," << lep2_centrality << "," << jets_centrality << "," <<  dilep_centrality << "," << ngapjets << "," << njets << "," << mt << "," << eTmiss_in << "," << eTmiss_out << "," << rapl1 << "," << rapl2 << "," << delta_R_l1j1 << "," << delta_R_l2j2 << "," << delta_phi_met_j1 << "," << delta_phi_met_j2 << "," << nu_from_hadron << "," << et_nu << "," << delta_phi_met_jj << ","  << eTmiss_muons << "," <<  MET << "," << MET_smeared_atlas << "," << MET_smeared_cms << "," << MET_muons << "," << MET_smeared_atlas_muons << "," << MET_smeared_cms_muons << std::endl;

}

//================================================================================//
//=================================  Finalising  =================================//
//================================================================================//

/// Normalise histograms, scale/divide histograms etc., after the run
void finalize() {
  //uncertainty on sum of weights
  const double err_sumW=std::pow(var_sumW,0.5);

  //uncertainty on fiducial sum of weights
  const double err_fid_sumW=std::pow(var_fid_sumW,0.5);

  //total cross-section in femtobarns
  const double total_xs = crossSection()/femtobarn;

  //fiducial cross-section in femtobarns
  const double fid_xs = fid_sumW*total_xs/sumW;

  //error on cross-section generated from MadGraph (hard-coded in femtobarns)
  const double err_xs_mg = 0 ;// 0.18; //too big?

  //error on fiducial cross-section
  const double err_fid_xs = pow(pow(err_fid_sumW,2)*pow(total_xs/sumW,2)+pow(err_xs_mg,2)*pow(fid_sumW/sumW,2)+pow(err_sumW,2)*pow(fid_sumW*total_xs/pow(sumW,2),2),0.5);
  std::cout << "SMEAR OPTION: " << SMEAR << ", JET_SMEAR: " << JET_SMEAR << ", MUON_SMEAR: " << MUON_SMEAR << std::endl;
  std::cout << total_xs << std::endl;
  std::cout << "******* IMPORTANT INFORMATION *********" << std::endl;
  std::cout << "Number of Events surviving" << event_number << std::endl;
  std::cout << "Number of Events Before Num Jets cut" << num_events_before_a_cut << std::endl;
  std::cout << "Number of Events After Num Jets cut" << num_events_after_a_cut << std::endl;
  std::cout << " Total XS MadGraph:  " << total_xs << std::endl;
  std::cout << " Total XS error: " << err_xs_mg << std::endl;
  std::cout << " Total Sum of Weights: " << sumW << std::endl;
  std::cout << " Error on Total Sum of Weights: " <<  err_sumW << std::endl;
  std::cout << " Selected Events Sum of Weights: " << fid_sumW << std::endl;
  std::cout << " Error on Selected Events Sum of Weights: " << err_fid_sumW << std::endl;
  std::cout << " Fiducial crossection: " << fid_xs << std::endl;
  std::cout << " Error of fiducial crossection: " << err_fid_xs << std::endl;
  std::cout << "************************" << std::endl;

  outFile.close();
}

private:
  std::ofstream outFile;
  int SMEAR;
  bool MET_SMEAR;
  bool JET_SMEAR;
  bool MUON_SMEAR;

}; //closes the majorana class at the very top

//=================================  Plugin hook  =================================//
// The hook for the plugin system
DECLARE_RIVET_PLUGIN(MLRSM_mumujj);


}  //closes the Rivet namespace

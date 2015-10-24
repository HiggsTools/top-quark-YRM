// -*- C++ -*-
#include <iostream>
#include <sstream>
#include <string>

#include "Rivet/Analysis.hh"

#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/UnstableFinalState.hh"

#include "Rivet/Jet.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/ChargedLeptons.hh"

#include "fastjet/internal/base.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/PseudoJet.hh"

namespace Rivet {
    
    class HIGGSTOOLS_2015 : public Analysis {
    public:
        
        /// Constructor
        HIGGSTOOLS_2015()
        : Analysis("HIGGSTOOLS_2015")
        {
        }
        
        
        /// Book histograms and initialise projections before the run
        void init() {
            FinalState fs;
            addProjection(fs, "FS");
            
            FastJets antikt4fj(FinalState(-5.0,5.0), FastJets::ANTIKT, 0.4);
            addProjection(antikt4fj, "AntiKt4Jets");
            
            IdentifiedFinalState photonfs(Cuts::abseta < 2.37 && Cuts::pT > 16*GeV);
            photonfs.acceptId(PID::PHOTON);
            addProjection(photonfs, "Photon");
            
            addProjection(UnstableFinalState(), "UFS");
            
            // bare leptons
            ChargedLeptons bare_leptons(fs);
            
            // dressed leptons
            Cut cuts = (Cuts::abseta < 2.5) & (Cuts::pT > 25*GeV);
            DressedLeptons leptons(fs, bare_leptons, 0.1, cuts);
            addProjection(leptons, "leptons");
            
            // Don't make pT or eta cuts on the neutrino
            IdentifiedFinalState neutrinos(fs);
            neutrinos.acceptNeutrinos();
            addProjection(neutrinos, "Neutrinos");
            
            _h_Et_photon_0 = bookHisto1D("Et_photon_0",100,0,100);
            _h_Et_photon_1 = bookHisto1D("Et_photon_1",100,0,100);
            _h_Eta_photon_0 = bookHisto1D("Eta_photon_0",100,-3,3);
            _h_Eta_photon_1 = bookHisto1D("Eta_photon_1",100,-3,3);
            _h_diphoton_mass = bookHisto1D("diphoton_mass",200,0,200);
            
            _h_costheta_labF_ll = bookHisto1D("costheta_labF_ll",50,-1,1);
            _h_costheta_labF_bb = bookHisto1D("costheta_labF_bb",50,-1,1);
            _h_deltaeta_labF_ll = bookHisto1D("deltaeta_labF_ll",50,0,3);
            _h_deltaeta_labF_bb = bookHisto1D("deltaeta_labF_bb",50,0,3);
        }
        
        /// Perform the per-event analysis
        void analyze(const Event& event) {
            const double weight = event.weight();
            
            FourMomentum top;
            FourMomentum antiTop;
            
            bool foundTop=false;
            bool foundAntiTop=false;
            bool foundW_TauX=false;
            
            Particles neutrinos;
           
            foreach (const GenParticle* p, Rivet::particles(event.genEvent())) {
                
                // first try to find Z->ee
                // skip if already found Z->ee
                if( abs(p->pdg_id())==6 && (!foundTop || !foundAntiTop) ) {
                    
                    const GenVertex* pv = p->production_vertex();
                    if (pv) {
                        HepMC::GenVertex::particle_iterator firstMother, thisMother, lastMother;
                        const HepMC::GenParticle* mother;
                        
                        firstMother = p->production_vertex()->particles_begin(HepMC::parents);
                        lastMother = p->production_vertex()->particles_end(HepMC::parents);
                        //Main loop:  iterate over all mothers, call method recursively.
                        for( thisMother = firstMother; thisMother != lastMother++; ++thisMother)   {
                            mother = (*thisMother);
                            int mother_id = abs(mother->pdg_id());
                            if( mother_id==1 || mother_id==2 || mother_id==3 || mother_id==4 || mother_id==21 ){
                                if( p->pdg_id()==6 ) {
                                    foundTop=true;
                                    top = p->momentum();
                                }
                                else {
                                    foundAntiTop=true;
                                    antiTop = p->momentum();
                                }
                            }// end conditional over electron/positrom from ME
                        }// end mother loop
                    }// end mother loop
                }// loop over parents only to search for electron/positron from ME              
                if( abs(p->pdg_id())==24) {
                    // loop over childs
                    const GenVertex* dv = p->end_vertex();
                    if (dv) {
                        for (GenVertex::particles_out_const_iterator pp = dv->particles_out_const_begin() ;
                             pp != dv->particles_out_const_end() ; ++pp) {// look for Hbb decays
                            if( abs((*pp)->pdg_id())==15 ) {foundW_TauX=true;}
                            if( abs((*pp)->pdg_id())==12 || abs((*pp)->pdg_id())==14 || abs((*pp)->pdg_id())==16){
                                Particle dummy((*pp)->pdg_id(),p->momentum());
                                neutrinos.push_back(dummy);
                            }
                        }
                    }
                }// end conditional over Wdecay into Tau
            }// end loop over particles in event record

            if(foundW_TauX){vetoEvent;}
            if(neutrinos.size()<2) {vetoEvent;}
            
            Particles photons = applyProjection<IdentifiedFinalState>(event, "Photon").particlesByPt();
            
            // Isolate photons with ET_sum in cone
            Particles isolated_photons;
            Particles fs = applyProjection<FinalState>(event, "FS").particles();
            foreach (const Particle& photon, photons) {
                FourMomentum mom_in_cone;
                double eta_P = photon.eta();
                double phi_P = photon.phi();
                foreach (const Particle& p, fs) {
                    if (deltaR(eta_P, phi_P, p.eta(), p.phi()) < 0.4) {
                        mom_in_cone += p.momentum();
                    }
                }
                if (mom_in_cone.Et()-photon.Et() < 1.0*GeV) {
                    isolated_photons.push_back(photon);
                }
            }
            
            if (isolated_photons.size() != 2) {
                vetoEvent;
            }
            
            const UnstableFinalState &ufs = applyProjection<UnstableFinalState> (event, "UFS");
            
            int bhadron_counter = 0;
            vector<Particle> bhadrons;
            /// Do the event by event analysis here
            foreach( const Particle &p, ufs.particles() )
            {
                const PdgId id = abs(p.pdgId());
                /// Only interested in bottom hadrons
                if( (PID::isHadron(id) && PID::hasBottom(id)) ) {
                    // Only interested in final bottom hadrons
                    bool isFinal = true;
                    foreach( const Particle&pp, p.allDescendants() ) {
                        if( PID::hasBottom( abs(pp.pdgId()) ) ) {isFinal = false;}
                    }
                    if(isFinal){
                        bhadrons.push_back(p);
                        //cout << "  pT " << p.pT() << " eta " << p.eta() << " phi " << p.phi() << endl;
                        bhadron_counter++;
                    }
                }// end b-hadron conditional
            }// end loop over particles
            //cout << " nbhadrons " << bhadron_counter << endl;
            bhadrons = sortByPt( bhadrons );
            
            Jets small_jets;
            foreach ( const Jet& jet, applyProjection<FastJets>(event, "AntiKt4Jets").jetsByPt(25.0*GeV) ){
                if ( fabs( jet.eta() ) < 3.0 ) {
                    small_jets.push_back(jet);
                }
            }// end jet loop
            
            Jets btagged_small_jets;
            vector<int> tagged_jets;
            if(small_jets.size()>0){
                for (int i_hadron=0; i_hadron<(int)bhadrons.size(); i_hadron++){
                    int jet_idx = 0;// start with leading jet
                    float dR = deltaR(bhadrons[i_hadron],small_jets[jet_idx]);
                    for (unsigned int i_jet=0; i_jet < small_jets.size(); i_jet++){
                        float _dR_ = deltaR(bhadrons[i_hadron],small_jets[i_jet]);
                        if(_dR_<dR) {
                            jet_idx=i_jet;
                            dR = _dR_;
                        }
                    }// end loop over jets
                    if(dR < 0.75*0.4) {
                        // tag only once
                        if(find(tagged_jets.begin(),tagged_jets.end(),jet_idx)==tagged_jets.end()){
                            btagged_small_jets.push_back(small_jets[jet_idx]);
                            tagged_jets.push_back(jet_idx);
                        }
                    }
                }// end loop over hadrons
            }// need at least one jet
            
            if (btagged_small_jets.size() != 2) {
                vetoEvent;
            }
            
            const vector<DressedLepton>& leptons = applyProjection<DressedLeptons>(event, "leptons").dressedLeptons();
            if ( leptons.size() < 2 )  vetoEvent;
            
            _h_Et_photon_0->fill(isolated_photons[0].Et(), weight);
            _h_Eta_photon_0->fill(isolated_photons[0].eta(), weight);
            
            _h_Et_photon_1->fill(isolated_photons[1].Et(), weight);
            _h_Eta_photon_1->fill(isolated_photons[1].eta(), weight);
            
            FourMomentum higgs = isolated_photons[0].momentum() + isolated_photons[1].momentum();
            _h_diphoton_mass->fill(higgs.mass(), weight);
            
            _h_costheta_labF_ll->fill(cos(angle(leptons[0].momentum(),leptons[1].momentum())), weight);
            _h_deltaeta_labF_ll->fill(abs(leptons[0].eta()-leptons[1].eta()), weight);
            _h_costheta_labF_bb->fill(cos(angle(btagged_small_jets[0].momentum(),btagged_small_jets[1].momentum())), weight);
            _h_deltaeta_labF_bb->fill(abs(btagged_small_jets[0].eta()-btagged_small_jets[1].eta()), weight);
            
            //std::cout << " found top "<< foundTop << "  found anti-top " << foundAntiTop << std::endl;
            //const Particles& neutrinos = applyProjection<IdentifiedFinalState>(event, "Neutrinos").particlesByPt();
            //std::cout << "  number of neutrinos " << neutrinos.size() << std::endl;
            //std::cout << "  found tau " << foundW_TauX << std::endl;
            
            // find the right neutrino for each lepton
            //FourMomentum W0_0=leptons[0].momentum()+neutrinos[0].momentum();
            //FourMomentum W0_1=leptons[0].momentum()+neutrinos[1].momentum();
            //FourMomentum W1_0=leptons[1].momentum()+neutrinos[0].momentum();
            //FourMomentum W1_1=leptons[1].momentum()+neutrinos[1].momentum();
            
            //std::cout << " mass W0_0 "<<W0_0.mass()<<std::endl;
            //std::cout << " mass W0_1 "<<W0_1.mass()<<std::endl;
            //std::cout << " mass W1_0 "<<W1_0.mass()<<std::endl;
            //std::cout << " mass W1_1 "<<W1_1.mass()<<std::endl;
        }
        
        
        /// Normalise histograms etc., after the run
        void finalize() {
            scale(_h_Et_photon_0, 1./sumOfWeights());
            scale(_h_Eta_photon_0, 1./sumOfWeights());
            
            scale(_h_Et_photon_1, 1./sumOfWeights());
            scale(_h_Eta_photon_1, 1./sumOfWeights());
            
            scale(_h_diphoton_mass, 1./sumOfWeights());
            
            scale(_h_costheta_labF_ll, 1./sumOfWeights());
            scale(_h_costheta_labF_bb, 1./sumOfWeights());
            scale(_h_deltaeta_labF_ll, 1./sumOfWeights());
            scale(_h_deltaeta_labF_bb, 1./sumOfWeights());
        }
        
    private:
        
        Histo1DPtr _h_Et_photon_0;
        Histo1DPtr _h_Eta_photon_0;
        Histo1DPtr _h_Et_photon_1;
        Histo1DPtr _h_Eta_photon_1;
        Histo1DPtr _h_diphoton_mass;
        
        Histo1DPtr _h_costheta_labF_ll;
        Histo1DPtr _h_costheta_labF_bb;
        Histo1DPtr _h_deltaeta_labF_ll;
        Histo1DPtr _h_deltaeta_labF_bb;
    };
    
    // The hook for the plugin system
    DECLARE_RIVET_PLUGIN(HIGGSTOOLS_2015);
    
}

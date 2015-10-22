// -*- C++ -*-
#include <iostream>
#include <sstream>
#include <string>

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Jet.hh"
#include "Rivet/Projections/FastJets.hh"

#include "fastjet/internal/base.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/PseudoJet.hh"

namespace Rivet {


  class HIGGSTOOLS_2015_I1 : public Analysis {
  public:

    /// Constructor
    HIGGSTOOLS_2015_I1()
      : Analysis("HIGGSTOOLS_2015_I1")
    {
    }


    /// Book histograms and initialise projections before the run
    void init() {
      FinalState fs;
      addProjection(fs, "FS");

      IdentifiedFinalState photonfs(Cuts::abseta < 2.37 && Cuts::pT > 16*GeV);
      photonfs.acceptId(PID::PHOTON);
      addProjection(photonfs, "Photon");

      _h_Et_photon_0 = bookHisto1D("Et_photon_0",100,0,100);
      _h_Et_photon_1 = bookHisto1D("Et_photon_1",100,0,100);
      _h_Eta_photon_0 = bookHisto1D("Eta_photon_0",100,-3,3);
      _h_Eta_photon_1 = bookHisto1D("Eta_photon_1",100,-3,3);
      _h_diphoton_mass = bookHisto1D("diphoton_mass",200,0,200);
      
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      Particles photons = applyProjection<IdentifiedFinalState>(event, "Photon").particlesByPt();

      std::cout << " number of photons " << photons.size() << std::endl;

      if (photons.size() < 1) {
        vetoEvent;
      }

      FourMomentum leadingPhoton = photons[0].momentum();

      _h_Et_photon_0->fill(leadingPhoton.Et(), weight);
      _h_Eta_photon_0->fill(leadingPhoton.eta(), weight);
      
      if (photons.size() < 2) {
        vetoEvent;
      }

      _h_Et_photon_1->fill(photons[1].Et(), weight);
      _h_Eta_photon_1->fill(photons[1].eta(), weight);
      
      FourMomentum higgs = photons[0].momentum() + photons[1].momentum();
      _h_diphoton_mass->fill(higgs.mass(), weight);

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_Et_photon_0, crossSection()/sumOfWeights());
      scale(_h_Eta_photon_0, crossSection()/sumOfWeights());

      scale(_h_Et_photon_1, crossSection()/sumOfWeights());
      scale(_h_Eta_photon_1, crossSection()/sumOfWeights());

      scale(_h_diphoton_mass, crossSection()/sumOfWeights());
    }


  private:

    Histo1DPtr _h_Et_photon_0;
    Histo1DPtr _h_Eta_photon_0;
    Histo1DPtr _h_Et_photon_1;
    Histo1DPtr _h_Eta_photon_1;
    Histo1DPtr _h_diphoton_mass;

    //    fastjet::AreaDefinition* _area_def;

    //    std::vector<float> _eta_bins;
    //   std::vector<float> _eta_bins_areaoffset;

    //   std::vector<float> _ptDensity;
    //   std::vector<float> _sigma;
    //   std::vector<float> _Njets;
  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(HIGGSTOOLS_2015_I1);

}

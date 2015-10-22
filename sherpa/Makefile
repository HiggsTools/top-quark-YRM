
# Paths/Versions
SHERPA=Sherpa2.2.0
RIVET=rivet

# Files
SIGNAL=tth_yy.hepmc2g
BACKGROUND=ttyy.hepmc2g
SIGNAL_ANALYSIS=HIGGSTOOLS_2015_I1
SIGNAL_RUNCARD=Run_TTH_YY.dat
BACKGROUND_ANALYSIS=HIGGSTOOLS_2015_I2
BACKGROUND_RUNCARD=Run_TTYY.dat

plots: signal background

signal: $(SIGNAL)
	$(RIVET)-buildplugin Rivet_$(SIGNAL_ANALYSIS).so $(SIGNAL_ANALYSIS).cc
	$(RIVET) -a $(SIGNAL_ANALYSIS) $(SIGNAL) --pwd
	$(RIVET)-mkhtml Rivet.yoda

background: $(BACKGROUND)
	$(RIVET)-buildplugin Rivet_$(BACKGROUND_ANALYSIS).so $(BACKGROUND_ANALYSIS).cc
	$(RIVET) -a $(BACKGROUND_ANALYSIS) $(BACKGROUND) --pwd
	$(RIVET)-mkhtml Rivet.yoda

$(SIGNAL):
	$(SHERPA) -f $(SIGNAL_RUNCARD) 	

$(BACKGROUND):
	$(SHERPA) -f $(BACKGROUND_RUNCARD) 	

.PHONY: clean very-clean

clean:
	rm -f Rivet_$(SIGNAL_ANALYSIS).so
	rm -f Rivet.yoda
	rm -f Sherpa_References.tex

very-clean: clean
	rm -f $(SIGNAL)
	rm -f $(BACKGROUND)
	rm -rf Results/
	rm -rf Process/
	rm -f MIG_P+P+_14000_NNPDF_2.db
	rm -f MPI_Cross_Sections.dat
	rm -f Results.db
	rm -f Results.db.bak
	rm -rf plots/

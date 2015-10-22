
# Paths/Versions
SHERPA=Sherpa2.2.0
RIVET=rivet

# Files
EVENTS=tth_yy.hepmc2g
ANALYSIS=HIGGSTOOLS_2015_I1
RUNCARD=Run_TTH_YY.dat

plots: $(EVENTS)
	$(RIVET)-buildplugin Rivet_$(ANALYSIS).so $(ANALYSIS).cc
	$(RIVET) -a $(ANALYSIS) $(EVENTS) --pwd
	$(RIVET)-mkhtml Rivet.yoda


$(EVENTS):
	$(SHERPA) -f $(RUNCARD) 	

.PHONY: clean very-clean

clean:
	rm -f Rivet_$(ANALYSIS).so
	rm -f Rivet.yoda
	rm -f Sherpa_References.tex

very-clean: clean
	rm -f $(EVENTS)
	rm -rf Results/
	rm -rf Process/
	rm -f MIG_P+P+_14000_NNPDF_2.db
	rm -f MPI_Cross_Sections.dat
	rm -f Results.db
	rm -f Results.db.bak
	rm -rf plots/

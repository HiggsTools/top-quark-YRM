
# Paths/Versions
SHERPA=Sherpa2.2.0
RIVET=rivet

# Files
EVENTS=tth_yy.hepmc2g
ANALYSIS=HIGGSTOOLS_2015_I1
RUNCARD=Run_TTH_YY.dat

$(EVENTS):
	$(SHERPA) -f $(RUNCARD) 	

plots: $(EVENTS)
	$(RIVET)-buildplugin Rivet_$(ANALYSIS).so $(ANALYSIS).cc
	$(RIVET) -a $(ANALYSIS) $(EVENTS) --pwd
	$(RIVET)-mkhtml Rivet.yoda

.PHONY: clean very-clean

clean:
	rm Rivet_$(ANALYSIS).so
	rm Rivet.yoda

very-clean:
	rm -f $(EVENTS)
	rm -rf Results/
	rm -rf Process/

# Paths/Versions
include Makefile.inc

# Files/Folders
PLOTFILE = plot.plot
OUTPUT = plots
SIGNAL=tth
BACKGROUND=ttyy
ANALYSIS = HIGGSTOOLS_2015

plots: Rivet.yoda
	$(RIVET)-mkhtml --mc-errs -c $(PLOTFILE) -o $(OUTPUT) Rivet.yoda

Rivet.yoda: $(SIGNAL).hepmc2g $(BACKGROUND).hepmc2g
	$(RIVET)-buildplugin Rivet_$(ANALYSIS).so $(ANALYSIS).cc
	$(RIVET) -a $(ANALYSIS) $(SIGNAL).hepmc2g $(BACKGROUND).hepmc2g --pwd

%.hepmc2g : %.dat
	$(SHERPA) -f $< -R 2006,2010
	
.PHONY: clean very-clean

clean:
	rm -f *.db
	rm -f *.bak
	rm -f *.so
	rm -f *.yoda
	rm -f MPI_Cross_Sections.dat
	rm -f Sherpa_References.tex
	rm -rf Process/

very-clean: clean
	rm -f *.hepmc2g
	rm -rf plots/

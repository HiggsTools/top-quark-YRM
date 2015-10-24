# Paths/Versions
include Makefile.inc

%.hepmcg2 : %.dat
	$(SHERPA) -f $< -R RANDOM_SEED 2006 2010
	
.PHONY: clean very-clean

clean:

very-clean: clean


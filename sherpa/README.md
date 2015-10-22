# Sherpa Signal & Background Validation

Software
========

Install the following prerequisites:
* rivet
* LHAPDF
* root

install sherpa with:
~~~
$ svn co http://sherpa.hepforge.org/svn/branches/rel-2-2-0
$ cd rel-2-2-0
$ ./configure --enable-analysis --enable-root --enable-lhole --enable-versioning --enable-fastjet=path/to/rivet/fastjet-3.0.6/bin --enable-hepmc2=path/to/rivet/HepMC-2.06.09/bin --enable-rivet --enable-lhapdf=path/to/LHAPDF/LHAPDF-6.1.5/64 --libdir=path/to/Sherpa/rel-2-2-0 --prefix=path/to/Sherpa (--enable-mpi)
$ make
$ make install
~~~

Analysis
========

Our analysis is performed using rivet. We consider the tree level processes:
* pp -> t tb (H -> y y) __(Signal)__
* pp -> t tb y y __(Background)__

To alter the distributions produced edit __HIGGSTOOLS_2015_I1.cc__. 
The default analysis produces the following distributions:
* __TODO__

Plots
=====

First edit the Makefile to ensure the path to Sherpa and Rivet is correct.

Run:
~~~
$ make plots
~~~

Output can be viewed in a browser by navigating to  plots/index.html


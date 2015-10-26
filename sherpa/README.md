# Sherpa Signal & Background Validation

Installation
============

Install the following prerequisites:
* Rivet >2.4.0 (Including fastjet >3.0.6, HepMC >2.06.09, YODA >1.5.5): https://rivet.hepforge.org/
* LHAPDF >6.1.5 (Optional): https://lhapdf.hepforge.org/
* Root >6.04/06 (Optional): https://root.cern.ch/
* OpenMPI/MPI (Optional): various

Install Sherpa >2.2.0 with (e.g.):
~~~
$ svn co http://sherpa.hepforge.org/svn/branches/rel-2-2-0
$ cd rel-2-2-0
$ ./configure --enable-analysis --enable-root --enable-lhole --enable-versioning --enable-fastjet=path/to/rivet/fastjet-3.0.6/bin --enable-hepmc2=path/to/rivet/HepMC-2.06.09/bin --enable-rivet --enable-lhapdf=path/to/LHAPDF/LHAPDF-6.1.5/64 --libdir=path/to/Sherpa/rel-2-2-0 --prefix=path/to/Sherpa (--enable-mpi)
$ make
$ make install
~~~

Some of the above options (e.g. Root, LHAPDF) are not required but recommended. The enable-mpi option is not required but allows Sherpa to use multiple cores/ machines.

Plot Generation
===============

1. Edit Makefile.inc to ensure the path to Sherpa and Rivet is correct.
2. Run:
    ```
    $ make
    ```
3. (Optional) To remove intermediate files run:
    ```
    $ make clean
    ```

Output can be viewed in a browser by navigating to plots/index.html

Analysis
========

Our analysis is performed using Rivet. We consider the tree level processes:
* pp > t t~ (H > y y) __(Signal)__
* pp > t t~ y y __(Background)__

Also included are the top decays and some (rudimentary) detector effects. Any W+/- bosons produced (by the top decay) are assumed to decay leptonically, specifically, events in which they decay to q = u, d, s, c are not simulated.

Files
=====

Description:
* tth.dat (Sherpa run card for signal)
* ttyy.dat (Sherpa run card for background)
* HIGGSTOOLS_2015.cc (Rivet analysis file)
* plot.plot (Rivet plot settings)

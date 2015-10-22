# top-quark-YRM

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

Plots
=====

Run:
~~~
make plots
~~~


# top-quark-YRM

### Goal as agreed on 10.10:
Reproduce polarization distribution plots from
  * http://arxiv.org/pdf/1403.1790v3.pdf
  * http://arxiv.org/pdf/1503.07787v2.pdf


Software
========

Install the following prerequisites:

  * madgraph 5: https://launchpad.net/mg5amcnlo
  * recent version of python2: (2.7 with numpy 1.9). If you don't have, get from here: https://www.continuum.io/downloads
  * ROOT 6: download source from https://root.cern.ch/download/root_v6.05.02.source.tar.gz and compile with `./configure && make`


To check python version
~~~
$ python --version
Python 2.7.9
$ python -c "import numpy; print numpy.__version__"
1.9.2
~~~


Check that you have ROOT working in python
~~~
$ python
>>> import ROOT
>>> ROOT.gROOT.GetVersion()
'6.05/02'
~~~

In madgraph, install pythia with
~~~
./mg5_aMC
MG5_aMC> install pythia-pgs
~~~

Run cards
=========

## runcards/madgraph/tth_haa_madspin
Simulation of MadGraph+MadSpin+Pythia with `p p > t t~ h`, inclusive higgs decays.

To run, call
~~~
$ cd runcards/madgraph
./mg5_aMC
MG5_aMC> import tth_haa_madspin/Cards/proc_card_mg5.dat
MG5_aMC> launch
MG5_aMC> *** manually select MadSpin (4) and pythia (1), type "done"
~~~

**TODO**: figure out hadronization & detector simulation.

The LHE output will be in `runcards/madgraph/tth_haa_madspin//Events/run_01_decayed_1/tag_1_pythia_events.lhe.gz`

Unpack with
~~~
$ gzip -d runcards/madgraph/tth_haa_madspin//Events/run_01_decayed_1/tag_1_pythia_events.lhe.gz
~~~

## runcards/madgraph/ttgammagamma

**TODO**.


LHE tools
=========

## lhetools/lhe2root.py

Run as `python lhe2root.py file.lhe`. Will convert all events and particles to a root tree.
~~~
$ python lhetools/lhe2root.py runcards/madgraph/tth_haa_madspin//Events/run_01_decayed_1/tag_1_pythia_events.lhe
Setup complete
Opened file runcards/madgraph/tth_haa_madspin//Events/run_01_decayed_1/tag_1_pythia_events.lhe
Converting to .root format and outputing to lhe.root
$ du -csh lhe.root
8.7M	lhe.root
~~~

# Analysis tools

## analysis/lheanalysis.py

Run as `python analysis/lheanalysis.py file.root`. Will loop over all particles and print out the decay chain.

**TODO**: add simple reconstruction example.

Here is one event, commented.
~~~
---
#gluon gluon fusion (four-momentum) status mothers
p[g] (-61.71, 50.64, 381.63, 389.89) -1 []
p[g] (128.62, -13.42, -184.32, 225.16) -1 []

#top, antitop and higgs, coming from gluons
p[t] (82.03, 8.04, -45.68, 198.13) 2 ['g']
p[~t] (-25.77, 17.96, 229.70, 290.26) 2 ['g']
p[h] (10.64, 11.22, 13.29, 126.65) 2 ['g']

#top decay
p[b] (39.25, -25.99, 44.96, 65.27) 2 ['t']
p[W+] (42.78, 34.03, -90.64, 132.87) 2 ['t']

#antitop decay
p[~b] (42.48, 33.87, 23.98, 59.57) 2 ['~t']
p[W-] (-68.25, -15.91, 205.72, 230.69) 2 ['~t']

#W+ decay
p[u] (24.60, 22.90, -102.92, 108.26) 2 ['W+']
p[~d] (18.19, 11.13, 12.28, 24.60) 2 ['W+']

#W- decay
p[d] (1.04, -2.32, -6.53, 7.01) 2 ['W-']
p[~u] (-69.29, -13.59, 212.25, 223.68) 2 ['W-']

#higgs decay
p[b] (13.27, 31.94, 63.27, 72.27) 2 ['h']
p[~b] (-2.62, -20.72, -49.98, 54.38) 2 ['h']

#intermediate particles
p[g] (-130.00, 18.94, -2519.77, 2523.23) 1 []
p[b] (88.80, -45.54, 103.13, 145.91) 1 []
p[b] (46.59, 40.22, 21.29, 66.82) 1 []
p[g] (-60.77, 2.53, 189.89, 199.68) 1 []
p[g] (35.25, 5.74, -64.14, 74.70) 1 []
p[g] (-21.46, -18.56, -3.37, 29.95) 1 []
p[g] (7.57, 24.29, -42.94, 50.64) 1 []
~~~

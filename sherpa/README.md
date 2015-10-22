install rivet
install LHAPDF
install root

install sherpa with:
svn co http://sherpa.hepforge.org/svn/branches/rel-2-2-0
cd rel-2-2-0
./configure --enable-analysis --enable-root --enable-lhole --enable-versioning --enable-fastjet=/home/pcl335b/sjones/Documents/rivet/fastjet-3.0.6/bin --enable-hepmc2=/home/pcl335b/sjones/Documents/rivet/HepMC-2.06.09/bin --enable-rivet --enable-lhapdf=/home/pcl335b/sjones/Documents/LHAPDF/LHAPDF-6.1.5/64 --libdir=/home/pcl335b/sjones/Documents/Sherpa/rel-2-1-1 --prefix=/home/pcl335b/sjones/Documents/Sherpa (--enable-mpi)
make
make install

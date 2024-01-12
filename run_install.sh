
GCC="$(which gcc)"
GPLUSPLUS="$(which g++)"
GFORTRAN="$(which gfortran)"
mainDir="/home/swh514/Projects/dependencies/fftw"
export LDFLAGS="-L$mainDir/lib -L$mainDir/lib64"
export CPPFLAGS="-I$mainDir/include"
export CC=$GCC
export CXX=$GPLUSPLUS
export FC=$GFORTRAN

CC=$GCC CXX=$GPLUSPLUS ./configure CXXFLAGS=-std=c++11 --enable-mmdb --enable-minimol --enable-cif --enable-fortran --enable-ccp4 enable-gemmi F77=gfortran --prefix=$mainDir
make
make install

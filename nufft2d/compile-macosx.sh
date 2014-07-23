
PATH=/Applications/MATLAB_R2011b.app/bin/:$PATH

make -f makefile.mwrap -j4 TARGET=matlab-maci64-openmp clean
make -f makefile.mwrap -j4 TARGET=matlab-maci64-openmp  

make -f makefile.mwrap clean distclean


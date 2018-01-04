# SlipModel

A code (Slipfmtstopscmp.f90) and an example (s2012SUMATR01YUEx.fsp from equake-rc.info). 
The code Slipfmtstopscmp.f90 transforms different slip format to be compatible with Dr. Rongjiang Wang's code PSCMP.
The input includes slip models in Caltech format, USGS format, SIV format, GCMT format and more.

Compile it with gfortran. 
Run it like this:
./Slipfmtstopscmp
s2012SUMATR01YUEx.fsp fault_pscmp.dat 'SIV' 0
where, s2012SUMATR01YUEx.fsp is the input file, fault_pscmp.dat is the output file, 'SIV' is the format type of the input, 0 is the ocean thickness in km to be used.


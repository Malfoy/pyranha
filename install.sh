mkdir bin
cd thirdparty/KMC
make clean; make -j2
mv bin/kmc bin/kmc_tools bin/kmc_dump ../../bin; cd ../..
cd tools/piRANHA_mapping/
make clean; make -j2; cd ../kIWI_denovo
make clean; make -j2

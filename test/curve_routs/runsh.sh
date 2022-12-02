cd ..
cd ..
make install
cd test/curve_routs
gfortran -std=legacy test_near_stab_eval.f -o int2-stab -L/home/travis/lib -lfmm2dbie -lfmm2d -llapack -lblas
./int2-stab

cd ..
cd ..
make install
cd test/curve_routs
gfortran -std=legacy test_near_stab_eval.f -L/usr/local/lib -lfmm2dbie -lfmm2d -framework accelerate
./a.out

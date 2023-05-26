rm test-index-*.exe
rm test-index.log

# Run some basic leakage checks
gcc -g -o test-index-gcc-g.exe -Wall test-index.c -lm
valgrind ./test-index-gcc-g.exe 1000 1.0 0.10 5 
clang -fsanitize=address -g -o test-index-clang-g.exe -Wall test-index.c -lm
./test-index-clang-g.exe 1000 1.0 0.10 5 

# Run basic correctness tests
gcc -o test-index-gcc.exe -O2 -Wall test-index.c -lm
clang -o test-index-clang.exe -O2 -Wall test-index.c -lm

N=15000
D=1.2345
qhw=0.141

touch test-index.log
echo "*** CORRECTNESS: START ***" >> test-index.log
echo "using N=$N uniformly random test points..."

for i in 1 2 4 8 16 32
do
  echo "--- --- gcc $i --- ---" >> test-index.log
  ./test-index-gcc.exe $N $D $qhw $i >> test-index.log && echo "gcc OK $i"
  echo "--- --- clang $i --- ---" >> test-index.log
  ./test-index-clang.exe $N $D $qhw $i >> test-index.log && echo "clang OK $i"
done

echo "*** CORRECTNESS: DONE ***" >> test-index.log

# Build Octave executable (used in Octave test codes)
octave --eval "build_mex('qtreehat');"

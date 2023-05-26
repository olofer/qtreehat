rm test-index-gcc.exe
rm test-index-clang.exe
rm test-index.log

gcc -o test-index-gcc.exe -O2 -Wall test-index.c -lm
clang -o test-index-clang.exe -O2 -Wall test-index.c -lm

N=15000
D=1.2345
qhw=0.141

echo "using N=$N uniformly random test points..."

touch test-index.log
for i in 1 2 4 8 16 32
do
  echo "--- --- gcc $i --- ---" >> test-index.log
  ./test-index-gcc.exe $N $D $qhw $i >> test-index.log && echo "gcc OK $i"
  echo "--- --- clang $i --- ---" >> test-index.log
  ./test-index-clang.exe $N $D $qhw $i >> test-index.log && echo "clang OK $i"
done

# octave --eval "build_mex('qtreehat');"

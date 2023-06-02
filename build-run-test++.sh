rm test-index++*.exe
rm test-index++.log

g++ -std=c++14 -Wall -o test-index++g.exe -O2 test-index.cpp
clang++ -std=c++14 -Wall -o test-index++c.exe -O2 test-index.cpp

N=15000
D=1.2345
qhw=0.141

touch test-index++.log
echo "*** CORRECTNESS: START ***" >> test-index++.log

echo "--- --- g++ $i --- ---" >> test-index++.log
./test-index++g.exe $N $D $qhw >> test-index++.log && echo "g++ OK"

echo "--- --- clang++ $i --- ---" >> test-index++.log
./test-index++c.exe $N $D $qhw >> test-index++.log && echo "clang++ OK"

echo "*** CORRECTNESS: DONE ***" >> test-index++.log

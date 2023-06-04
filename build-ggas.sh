rm ggas.wasm
emcc -std=c++14 -Wall -msimd128 -O3 gravity-gas.cpp -s STANDALONE_WASM -fno-exceptions -DNDEBUG --no-entry -o ggas.wasm;

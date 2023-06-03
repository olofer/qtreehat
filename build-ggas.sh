rm ggas.wasm
emcc gravity-gas.cpp -s STANDALONE_WASM -fno-exceptions -DNDEBUG -std=c++14 -Wall -O3 --no-entry -o ggas.wasm;

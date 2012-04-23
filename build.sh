# C_OPTIONS="-ggdb -Wall -march=native -fopenmp" 
C_OPTIONS="-fwhole-program -fopenmp -msse2 -mfpmath=both -Wall -funsafe-math-optimizations -Wunsafe-loop-optimizations -march=native -O3 -DNDEBUG"

mkdir -p bin
mkdir -p images
g++ buddhabrot.cpp  $C_OPTIONS -o bin/buddhabrot
g++ test_random.cpp $C_OPTIONS -o bin/test_random



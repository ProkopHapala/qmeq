#!/bin/bash

# Get ASan library path
ASAN_LIB=$(gcc -print-file-name=libasan.so)
echo "Using ASan library: $ASAN_LIB"

# Set environment variables for ASan
export LD_PRELOAD=$ASAN_LIB
export ASAN_OPTIONS=detect_leaks=0

rm -f ../cpp/pauli_solver_wrapper.o
rm -f ../cpp/pauli_solver_wrapper.so

# Run the Python script with proper ASan preloading
python3 compare_solvers.py | tee compare_solvers.log

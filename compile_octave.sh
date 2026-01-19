#!/bin/bash
# Compile generateEllipseCandidates.cpp for Octave with OpenCV 4.12.0

set -e

# Get OpenCV flags
OPENCV_CFLAGS=$(pkg-config --cflags opencv4)
OPENCV_LIBS=$(pkg-config --libs opencv4)

# Use mkoctfile for Octave mex compilation
echo "Compiling with OpenCV CFLAGS: $OPENCV_CFLAGS"
echo "Linking with OpenCV LIBS: $OPENCV_LIBS"

# Compile with mkoctfile --mex
mkoctfile --mex generateEllipseCandidates.cpp $OPENCV_CFLAGS $OPENCV_LIBS -llapack -lblas 2>&1 | tee compile.log

if [ $? -eq 0 ]; then
    echo "Compilation successful"
    
    # Check if mex file was created
    if [ -f "generateEllipseCandidates.mex" ]; then
        echo "Mex file generated: generateEllipseCandidates.mex"
    else
        echo "Warning: Mex file not generated as expected"
        ls -la *.mex* 2>/dev/null || echo "No mex files found"
    fi
else
    echo "Compilation failed. See compile.log for details."
    exit 1
fi
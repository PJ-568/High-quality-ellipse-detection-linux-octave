#!/bin/bash
# Batch Ellipse Detection Script Runner
# Usage: ./batch_process.sh [input_folder] [output_folder] [Tac] [Tr] [polarity]
#
# Parameters:
#   input_folder:  Folder containing input images (default: 'pics')
#   output_folder: Folder to save results (default: 'results')
#   Tac:           Elliptic angular coverage (default: 165)
#   Tr:            Ratio of support inliers (default: 0.6)
#   polarity:      0: all ellipses, 1: positive, -1: negative (default: 0)
#
# Example:
#   ./batch_process.sh my_images my_results 165 0.6 0

set -e

# Default values
INPUT_FOLDER="pics"
OUTPUT_FOLDER="results"
TAC="165"
TR="0.6"
POLARITY="0"

# Parse command line arguments
if [ $# -ge 1 ]; then
    INPUT_FOLDER="$1"
fi
if [ $# -ge 2 ]; then
    OUTPUT_FOLDER="$2"
fi
if [ $# -ge 3 ]; then
    TAC="$3"
fi
if [ $# -ge 4 ]; then
    TR="$4"
fi
if [ $# -ge 5 ]; then
    POLARITY="$5"
fi

# Check if Octave is available
if ! command -v octave &> /dev/null; then
    echo "Error: Octave is not installed or not in PATH"
    exit 1
fi

# Check if input folder exists
if [ ! -d "$INPUT_FOLDER" ]; then
    echo "Error: Input folder '$INPUT_FOLDER' does not exist"
    exit 1
fi

# Create output folder if it doesn't exist
mkdir -p "$OUTPUT_FOLDER"

echo "========================================"
echo "Batch Ellipse Detection"
echo "========================================"
echo "Input folder:  $INPUT_FOLDER"
echo "Output folder: $OUTPUT_FOLDER"
echo "Parameters:"
echo "  Tac (elliptic angular coverage): $TAC"
echo "  Tr (ratio of support inliers):   $TR"
echo "  Polarity:                        $POLARITY"
echo "========================================"

echo "Starting ellipse detection..."
echo ""

# Run Octave with the batch function
octave --no-gui --eval "
try
    batch_ellipse_detection('$INPUT_FOLDER', '$OUTPUT_FOLDER', $TAC, $TR, $POLARITY);
catch err
    fprintf('Error: %s\\n', err.message);
    exit(1);
end
"

echo ""
echo "========================================"
echo "Batch processing finished!"
echo "Results are in: $OUTPUT_FOLDER"
echo "========================================"
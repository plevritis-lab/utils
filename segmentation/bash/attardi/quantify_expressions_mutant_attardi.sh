#!/bin/bash

SEGMENTATION_METHOD="mesmer"

IMAGE_DIRECTORY="$HOME/Desktop/attardi_analysis/data/mutant"
QUANTIFICATION_SCRIPT="$HOME/Desktop/utils/segmentation/source/quantify_expression.py"

for SAMPLE_DIRECTORY in "$IMAGE_DIRECTORY"/*; do
    SAMPLE_NAME=$(basename "$SAMPLE_DIRECTORY")
    SAMPLE_DATA="$SAMPLE_DIRECTORY/data/${SAMPLE_NAME}.qptiff"
    SAMPLE_SEGMENTATION="$SAMPLE_DIRECTORY/full/$SEGMENTATION_METHOD/E_CADHERIN/image_1_seg.npy"

    if [[ "$SAMPLE_NAME" != "histology" && "$SAMPLE_NAME" != "quantifications" && "$SAMPLE_NAME" != "assignments" ]]; then
        python3 "$QUANTIFICATION_SCRIPT" \
            --image_path "$SAMPLE_DATA" \
            --mask_path "$SAMPLE_SEGMENTATION" \
            --save_path "$IMAGE_DIRECTORY" \
            --apply_$SEGMENTATION_METHOD
    fi
done
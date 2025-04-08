#!/bin/bash

SEGMENTATION_METHOD="<TODO>"

PANEL_PATH="<TODO>"
IMAGE_DIRECTORY="<TODO>"
QUANTIFICATION_SCRIPT="<TODO>"

for SAMPLE_DIRECTORY in "$IMAGE_DIRECTORY"/*; do
    SAMPLE_NAME=$(basename "$SAMPLE_DIRECTORY")
    SAMPLE_DATA="$SAMPLE_DIRECTORY/data/${SAMPLE_NAME}.<TODO>"
    SAMPLE_SEGMENTATION="$SAMPLE_DIRECTORY/full/$SEGMENTATION_METHOD/<TODO>/image_1_seg.npy"

    if [[ "$SAMPLE_NAME" != "histology" && "$SAMPLE_NAME" != "quantifications" && "$SAMPLE_NAME" != "assignments" ]]; then
        python3 "$QUANTIFICATION_SCRIPT" \
            --image_path "$SAMPLE_DATA" \
            --mask_path "$SAMPLE_SEGMENTATION" \
            --panel_path "$PANEL_PATH" \
            --save_path "$IMAGE_DIRECTORY" \
            --apply_$SEGMENTATION_METHOD
    fi
done
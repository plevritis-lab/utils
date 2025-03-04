#!/bin/bash

PANEL_PATH="<TODO>"
IMAGE_DIRECTORY="<TODO>"
SEGMENTATION_SCRIPT="<TODO>"

for SAMPLE_DIRECTORY in "$IMAGE_DIRECTORY"/*; do
    SAMPLE_NAME=$(basename "$SAMPLE_DIRECTORY")
    
    if [[ "$SAMPLE_NAME" != "histology" ]]; then
        SAMPLE_DATA="$SAMPLE_DIRECTORY/data/${SAMPLE_NAME}.tif"

        python3 "$SEGMENTATION_SCRIPT" \
            --image_path "$SAMPLE_DATA" \
            --panel_path $PANEL_PATH \
            --nuclear_channel "<TODO>" \
            --segment_channel "<TODO>" \
            --<TODO>
    fi
done